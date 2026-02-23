#
# Adjacently: Julia Complex Directed Networks Library
# Copyright (C) 2016-2025 Jimmy Dubuisson <jimmy.dubuisson@gmail.com>
#
# This program is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.
#
# This program is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.
#
# RCGE: Reversible Coarsening Graph Encoding

module RCGE

using LightGraphs
using LightGraphs: AbstractGraph, outneighbors, nv, is_directed

using ...IO: BitWriter, BitReader, write_bit, write_bits, write_value, flush_bitwriter,
    read_bit, read_bits, read_value
import ..Compression
using ..Compression: write_encoded_value, write_delta, write_truncated_binary_coding,
    write_hybrid_mix_encoded_list, delta_encode_vector, write_elias_fano,
    read_encoded_value, read_delta, read_elias_fano,
    read_small_count, read_rle_ones_deltas, read_bitmap_rle_ones_deltas,
    read_bitpacked_bitmap, read_intervals_and_residuals,
    read_compressed_graph_data

export RCGEParams, RCGEStats, encode_level, decode_level, load_rcge_graph

"""
    RCGEStats

Lightweight counters (in bits) for encode_level sections.
"""
mutable struct RCGEStats
    bits_membership::Int
    bits_intra::Int
    bits_intra_headers::Int
    bits_intra_copy::Int
    bits_intra_add::Int
    bits_intra_raw::Int
    bits_intra_ref_small_headers::Int
    intra_ref_used::Int
    intra_no_ref::Int
    bits_inter_headers::Int
    bits_inter_degrees::Int
    bits_inter_perms::Int
    bits_inter_lists::Int
    ab_metrics::Vector{Tuple{Int,Int,Int,Int,Int,Int}}
    RCGEStats() = new(0,0,0,0,0,0,0,0,0,0,0,0,0, Tuple{Int,Int,Int,Int,Int,Int}[])
end

@inline function _total_bits(w::BitWriter)
    # Bytes written to underlying IO + buffered bits
    return Int(position(w.io)) * 8 + (w.index - 1)
end

"""
    RCGEParams(; L=32, varint=:fibonacci, count_varint=:fibonacci, gap=:fibonacci, degree=:golomb, perm_strategy=:lehmer, undirected_pairs=true, membership=:delta, inter_strategy=:perm,
                 intra_ref_enabled::Bool=true, intra_ref_window::Int=16, intra_ref_min_overlap::Float64=0.3, intra_ref_rle::Bool=true,
                 intra_block_try::Bool=false)

Parameters for RCGE encoding.
- L: cluster size threshold for bit-matrix vs list encoding
- varint: integer encoding used for sizes/lengths (positive only)
- count_varint: encoding for counts that may be zero (default :golomb supports zero)
- gap: integer encoding for gap-coded sorted lists
- degree: integer encoding for degree vectors
- perm_strategy: strategy for permutation encoding (:lehmer, :raw, :blockpos)
- undirected_pairs: if true and graph is undirected, only encode pairs A<=B
- membership: :delta or :elias_fano for cluster membership lists
 - inter_strategy: :perm (degree vectors + permutation) or :lists (explicit per-u neighbor lists in B)
"""
Base.@kwdef struct RCGEParams
    L::Int = 32
    varint::Symbol = :fibonacci
    count_varint::Symbol = :fibonacci
    gap::Symbol = :fibonacci
    degree::Symbol = :golomb
    perm_strategy::Symbol = :lehmer
    undirected_pairs::Bool = true
    membership::Symbol = :stop        # :stop | :delta | :elias_fano | :implicit_ranges
    inter_strategy::Symbol = :lists
    intra_ref_enabled::Bool = true
    intra_ref_window::Int = 16
    intra_ref_min_overlap::Float64 = 0.3
    intra_ref_rle::Bool = true
    intra_block_try::Bool = false
    positions_mode::Symbol = :delta   # :delta | :bitmap | :auto
    additions_mode::Symbol = :auto    # :auto | :delta | :intervals
    min_cluster_density::Float64 = 0.0
    intra_intervals::Bool = false     # Use MGS-style interval+residual encoding
    intra_mil::Int = 4                # Min interval length for interval detection
    intra_greedy_mil::Bool = false    # Per-vertex greedy mil search over {2,3,4,5}
    intra_mgs::Bool = false           # Use full MGS encoder (reference+interval+recursive) per cluster
    intra_zigzag::Bool = false        # Zigzag relative first-value encoding (offset from local vertex index)
    intra_stop_deltas::Bool = false   # Use STOP-terminated delta lists (eliminates per-vertex count prefix)
    intra_copy_blocks::Bool = false   # WebGraph-style copy-blocks for reference positions (vs RLE bitmap)
    intra_ref_fixwidth::Bool = false  # Fixed-width ref delta encoding (1-bit flag + ceil(log2(window)) bits per ref)
    intra_ref_vlc::Bool = false       # VLC (Fibonacci) for ref delta instead of fixed-width (requires intra_ref_fixwidth=true)
    intra_add_adaptive::Bool = false  # Per-vertex adaptive additions: pick cheaper of STOP-delta vs intervals (requires intra_stop_deltas=true)
    intra_raw_adaptive::Bool = false  # Per-vertex adaptive raw: pick cheaper of STOP-delta vs intervals (requires intra_stop_deltas=true)
    intra_adapt_mil::Int = 2          # MIL for adaptive interval encoding (2=most aggressive)
end

# --------------------------
# Helper encoders
# --------------------------

"""
    write_varint(w, v; encoding=:elias_gamma)

Write an unsigned integer `v` using the configured varint encoding.
"""
write_varint(w::BitWriter, v::T; encoding::Symbol=:elias_gamma) where {T<:Unsigned} = write_encoded_value(w, v, encoding)

"""
    write_gap_coded(w, sorted_vals; encoding=:elias_gamma)

Encode a sorted vector of positive integers via delta/gap coding.
"""
function write_gap_coded(w::BitWriter, sorted_vals::Vector{T}; encoding::Symbol=:elias_gamma) where {T<:Unsigned}
    isempty(sorted_vals) && return
    write_delta(w, sorted_vals, encoding)
end

"""
    write_stop_delta_list(w, sorted_vals; encoding=:fibonacci)

Write a STOP-terminated list with 1/0 prelude and delta-coded values.
For a sorted increasing list `sorted_vals`, emits for each value v:
 - bit '1', then varint-encoded delta (v - prev), with prev initialized to 0;
and finally a single '0' bit to terminate the list.
"""
function write_stop_delta_list(w::BitWriter, sorted_vals::Vector{T}; encoding::Symbol=:fibonacci) where {T<:Unsigned}
    prev = zero(T)
    for v in sorted_vals
        @assert v >= prev
        write_bit(w, true)
        write_encoded_value(w, T(v - prev), encoding)
        prev = v
    end
    write_bit(w, false) # STOP
end

"""
    _write_stop_delta_zigzag(w, sorted_vals, encoding, vertex_id)

STOP-terminated delta list with optional zigzag first value.
Format: for each value: bit '1' + varint(gap); then bit '0'.
Gaps are v[i] - v[i-1] (>= 1 for sorted unique lists, Fibonacci-safe).
When vertex_id is provided, the first value is zigzag-encoded relative to it.
"""
function _write_stop_delta_zigzag(w::BitWriter, sorted_vals::Vector{T}, encoding::Symbol, vertex_id) where {T<:Unsigned}
    if isempty(sorted_vals)
        write_bit(w, false)  # STOP
        return
    end
    # First value
    write_bit(w, true)
    if vertex_id !== nothing
        offset = Int64(sorted_vals[1]) - Int64(vertex_id)
        write_encoded_value(w, T(Compression._rl_zigzag_encode(offset) + 1), encoding)
    else
        write_encoded_value(w, sorted_vals[1], encoding)
    end
    # Remaining: gaps (>= 1 for sorted unique lists)
    for i in 2:length(sorted_vals)
        write_bit(w, true)
        write_encoded_value(w, sorted_vals[i] - sorted_vals[i-1], encoding)
    end
    write_bit(w, false)  # STOP
end

"""
    _read_stop_delta_zigzag(r, encoding, T, vertex_id)

Read a STOP-terminated delta list written by `_write_stop_delta_zigzag`.
"""
function _read_stop_delta_zigzag(r::BitReader, encoding::Symbol, ::Type{T}, vertex_id) where {T<:Unsigned}
    result = T[]
    first = true
    prev = zero(T)
    while read_bit(r)  # '1' = more values, '0' = STOP
        raw = read_encoded_value(r, encoding, T)
        if first && vertex_id !== nothing
            val = T(Int64(vertex_id) + Compression._rl_zigzag_decode(UInt64(raw - 1)))
            first = false
        else
            val = prev + raw
            first = false
        end
        push!(result, val)
        prev = val
    end
    return result
end

"""
    _write_copy_blocks(w, positions, encoding)

Write sorted position indices as WebGraph-style copy-blocks.
Positions are grouped into contiguous runs (copy blocks). Format:
  small_count(num_blocks), then for each block:
  - first block: varint(start), varint(length)
  - subsequent blocks: varint(gap_from_prev_end), varint(length)
All values >= 1 (1-based positions, gaps >= 1, lengths >= 1).
"""
function _write_copy_blocks(w::BitWriter, positions::Vector{Int}, encoding::Symbol)
    if isempty(positions)
        Compression.write_small_count(w, UInt32(0), encoding)
        return
    end
    # Compute copy blocks: contiguous runs of sorted positions
    blocks = Tuple{Int,Int}[]  # (start, length)
    i = 1
    while i <= length(positions)
        start = positions[i]
        len = 1
        while i + len <= length(positions) && positions[i + len] == start + len
            len += 1
        end
        push!(blocks, (start, len))
        i += len
    end

    ncb = length(blocks)
    Compression.write_small_count(w, UInt32(ncb), encoding)
    # First block: start (>= 1) and length (>= 1)
    write_encoded_value(w, UInt32(blocks[1][1]), encoding)
    write_encoded_value(w, UInt32(blocks[1][2]), encoding)
    # Subsequent blocks: gap from previous block end, then length
    prev_end = blocks[1][1] + blocks[1][2]
    for i in 2:ncb
        gap = blocks[i][1] - prev_end  # >= 1 for sorted unique positions
        write_encoded_value(w, UInt32(gap), encoding)
        write_encoded_value(w, UInt32(blocks[i][2]), encoding)
        prev_end = blocks[i][1] + blocks[i][2]
    end
end

"""
    _read_copy_blocks(r, encoding, T) → Vector{Int}

Read copy-blocks written by `_write_copy_blocks`. Returns sorted position indices.
"""
function _read_copy_blocks(r::BitReader, encoding::Symbol, ::Type{T}) where {T<:Unsigned}
    ncb = Int(Compression.read_small_count(r, encoding, T))
    if ncb == 0
        return Int[]
    end
    positions = Int[]
    # First block
    start = Int(read_encoded_value(r, encoding, T))
    len = Int(read_encoded_value(r, encoding, T))
    for j in 0:(len-1)
        push!(positions, start + j)
    end
    prev_end = start + len
    # Subsequent blocks
    for i in 2:ncb
        gap = Int(read_encoded_value(r, encoding, T))
        start = prev_end + gap
        len = Int(read_encoded_value(r, encoding, T))
        for j in 0:(len-1)
            push!(positions, start + j)
        end
        prev_end = start + len
    end
    return positions
end

"""
    _write_fixwidth_ref_deltas(w, use_ref_vec, ref_delta_vec, window; vlc=false)

Per-vertex ref delta encoding. For each vertex in order:
  - 1 bit: has_ref flag
  - If has_ref and vlc=false: (delta - 1) in ceil(log2(window)) fixed bits
  - If has_ref and vlc=true:  Fibonacci(delta) VLC code (saves bits for small deltas)
"""
function _write_fixwidth_ref_deltas(w::BitWriter, use_ref_vec::Vector{Bool},
                                    ref_delta_vec::Vector{UInt32}, window::Int;
                                    vlc::Bool=false)
    if vlc
        for idx in 1:length(use_ref_vec)
            write_bit(w, use_ref_vec[idx])
            if use_ref_vec[idx]
                write_encoded_value(w, ref_delta_vec[idx], :fibonacci)  # delta ≥ 1
            end
        end
    else
        nbits = max(1, ceil(Int, log2(window)))
        for idx in 1:length(use_ref_vec)
            if use_ref_vec[idx]
                write_bit(w, true)
                d = Int(ref_delta_vec[idx]) - 1  # delta ∈ [1, window] → [0, window-1]
                for b in (nbits-1):-1:0
                    write_bit(w, ((d >> b) & 1) == 1)
                end
            else
                write_bit(w, false)
            end
        end
    end
end

"""
    _read_fixwidth_ref_deltas(r, s, window; vlc=false) → (use_ref_vec, ref_delta_vec)

Read per-vertex ref delta encoding written by `_write_fixwidth_ref_deltas`.
"""
function _read_fixwidth_ref_deltas(r::BitReader, s::Int, window::Int; vlc::Bool=false)
    use_ref_vec = Vector{Bool}(undef, s)
    ref_delta_vec = zeros(UInt32, s)
    if vlc
        for idx in 1:s
            has_ref = read_bit(r)
            use_ref_vec[idx] = has_ref
            if has_ref
                ref_delta_vec[idx] = read_encoded_value(r, :fibonacci, UInt32)
            end
        end
    else
        nbits = max(1, ceil(Int, log2(window)))
        for idx in 1:s
            has_ref = read_bit(r)
            use_ref_vec[idx] = has_ref
            if has_ref
                d = 0
                for _ in 1:nbits
                    d = (d << 1) | Int(read_bit(r))
                end
                ref_delta_vec[idx] = UInt32(d + 1)  # restore delta ∈ [1, window]
            end
        end
    end
    return use_ref_vec, ref_delta_vec
end

"""
    write_compact_int_vector(w, vals; encoding=:golomb)

Writes a small nonnegative integer vector as: [len][each value].
"""
function write_compact_int_vector(w::BitWriter, vals::Vector{T}; encoding::Symbol=:golomb) where {T<:Unsigned}
    write_varint(w, T(length(vals)); encoding=encoding)
    # Some integer codings don't support 0 (e.g., :elias_delta, :fibonacci),
    # in which case we store nonnegative values shifted by +1.
    @inline zero_unsupported(enc::Symbol) = (enc == :elias_delta || enc == :fibonacci || enc == :zeta)
    for v in vals
        local outv = v
        if zero_unsupported(encoding)
            outv = v + one(T)
        end
        write_encoded_value(w, outv, encoding)
    end
end

"""
    write_sparse_int_vector(w::BitWriter, vals::Vector{T}; index_encoding::Symbol=:fibonacci, value_encoding::Symbol=:elias_delta) where {T<:Unsigned}

Encode nonnegative integer vector sparsely:
- [nnz: encoded with `value_encoding` (shifted by +1 if zero-unsafe)]
- positions (1-based) delta-coded with `index_encoding`
- values (non-zero) encoded with `value_encoding` (shifted by +1 if zero-unsafe)
"""
function write_sparse_int_vector(w::BitWriter, vals::Vector{T}; index_encoding::Symbol=:fibonacci, value_encoding::Symbol=:elias_delta) where {T<:Unsigned}
    pos = Int[]
    nzv = T[]
    for i in 1:length(vals)
        vi = vals[i]
        if vi != zero(T)
            push!(pos, i)
            push!(nzv, vi)
        end
    end
    nnz = length(pos)
    # write nnz (shift if encoding doesn't support 0)
    local cntu = UInt(nnz)
    if value_encoding == :elias_delta || value_encoding == :fibonacci || value_encoding == :zeta
        cntu = UInt(nnz + 1)
    end
    write_encoded_value(w, cntu, value_encoding)
    if nnz == 0
        return
    end
    # write positions (>=1) as delta-coded
    write_delta(w, UInt32.(pos), index_encoding)
    # write values (shift if zero-unsafe)
    for v in nzv
        local outv = v
        if value_encoding == :elias_delta || value_encoding == :fibonacci || value_encoding == :zeta
            outv = v + one(T)
        end
        write_encoded_value(w, outv, value_encoding)
    end
end

"""
    write_upper_tri_bitset(w, adjmat)

Write upper-triangular bitset (excluding diagonal) for an undirected cluster.
`adjmat` is an sxs Bool matrix or callable adj accessor.
Bits are emitted row-major for pairs (i<j): i=1..s-1, j=i+1..s.
"""
function write_upper_tri_bitset(w::BitWriter, adjmat::AbstractMatrix{Bool})
    s = size(adjmat, 1)
    @assert size(adjmat, 2) == s
    for i in 1:s-1
        for j in (i+1):s
            write_bit(w, adjmat[i,j])
        end
    end
end

"""
    write_bitset_cluster(w, g, cluster)

For undirected graphs and small clusters, write upper-tri bitset of C's induced subgraph.
Fallback to adjacency-list mode for directed graphs.
"""
function write_bitset_cluster(w::BitWriter, g::AbstractGraph{T}, cluster::Vector{T}) where {T<:Unsigned}
    s = length(cluster)
    # Map global -> local index (1..s)
    local_index = Dict{T,Int}()
    for (i,u) in enumerate(cluster)
        local_index[u] = i
    end
    # Build adjacency matrix upper triangle for undirected graphs
    A = falses(s, s)
    for (i,u) in enumerate(cluster)
        for v in outneighbors(g, Int(u))
            if haskey(local_index, v)
                j = local_index[v]
                A[i,j] = true
                A[j,i] = true
            end
        end
    end
    write_upper_tri_bitset(w, A)
end

"""
    order_from_edges(EAB, A_local, B_local, neighbors_in_B)

Build permutation pi for stub matching based on canonical ordering.
- `A_local`: vertices of A in local order (Vector{T})
- `B_local`: vertices of B in local order (Vector{T})
- `neighbors_in_B[u]`: sorted list of B-local neighbors for each u∈A

Returns pi as Vector{Int} of size m with values in 1..m.
"""
function order_from_edges(A_local::Vector{T}, B_local::Vector{T}, neighbors_in_B::Dict{T,Vector{Int}}) where {T<:Unsigned}
    # Degree vectors and stub counts
    d_out = Int[]
    for u in A_local
        push!(d_out, length(get(neighbors_in_B, u, Int[])))
    end
    d_in_map = Dict{Int,Int}()  # B-local idx => count
    for (_, Ns) in neighbors_in_B
        for vloc in Ns
            d_in_map[vloc] = get(d_in_map, vloc, 0) + 1
        end
    end
    # Total stubs
    m = sum(d_out)
    # Precompute base offsets for IN_STUBS: B vertices in local order 1..|B|
    bases = Dict{Int,Int}()
    off = 1
    for (j,_) in enumerate(B_local)
        bases[j] = off
        off += get(d_in_map, j, 0)
    end
    next_pos = Dict{Int,Int}()
    for (j,_) in enumerate(B_local)
        next_pos[j] = get(bases, j, 0)
    end

    # Walk OUT_STUBS order: for each u in A_local, list Ns in ascending local v
    pi = Vector{Int}(undef, m)
    i = 1
    for u in A_local
        Ns = get(neighbors_in_B, u, Int[])
        for vloc in Ns
            pos = next_pos[vloc]
            if pos == 0
                # No stubs available (should not happen in consistent input)
                error("Stub construction error: no available in-stub for vloc=$vloc")
            end
            pi[i] = pos
            next_pos[vloc] = pos + 1
            i += 1
        end
    end
    return pi
end

"""
    write_permutation(w, pi; strategy=:lehmer, pos_encoding=:fibonacci, count_encoding=:elias_delta, block_sizes=nothing)

Write permutation pi (1..m) using selected strategy.
- :lehmer encodes digits with truncated binary coding in ranges [0..m-i].
- :raw writes m then pi[i] as varints (debug-friendly, not compact).
- :blockpos requires `block_sizes` (d_in), and encodes per-block OUT positions delta-coded.
"""
function write_permutation(w::BitWriter, pi::Vector{Int}; strategy::Symbol=:lehmer, pos_encoding::Symbol=:fibonacci, count_encoding::Symbol=:elias_delta, block_sizes=nothing)
    m = length(pi)
    if strategy == :raw
        write_varint(w, UInt(m))
        for x in pi
            write_varint(w, UInt(x))
        end
        return
    end
    if strategy == :blockpos
        block_sizes === nothing && error("blockpos requires block_sizes")
        bs = Int.(block_sizes)
        nb = length(bs)
        block_of_pos = Vector{Int}(undef, m)
        off = 1
        for j in 1:nb
            s = bs[j]
            for k in 0:(s-1)
                block_of_pos[off + k] = j
            end
            off += s
        end
        pos_lists = [Int[] for _ in 1:nb]
        for p in 1:m
            b = block_of_pos[pi[p]]
            push!(pos_lists[b], p)
        end
        for j in 1:nb
            c = length(pos_lists[j])
            # write count (shift if zero-unsafe)
            local cc = UInt(c)
            if count_encoding == :elias_delta || count_encoding == :fibonacci || count_encoding == :zeta
                cc = UInt(c + 1)
            end
            write_encoded_value(w, cc, count_encoding)
            if c > 0
                # delta-code strictly increasing positions
                write_delta(w, UInt32.(pos_lists[j]), pos_encoding)
            end
        end
        return
    end
    # Lehmer code
    # Maintain an ordered list of remaining items 1..m. Use a Fenwick tree alternative with a simple vector for now.
    remaining = collect(1:m)
    for i in 1:m
        x = pi[i]
        # Find rank of x in remaining (0-based)
        idx = findfirst(==(x), remaining)
        idx === nothing && error("Invalid permutation element $x at position $i")
        digit = idx - 1  # 0-based digit in [0 .. (m-i)]
        # Range size for digit is (m - i + 1)
        nrange = m - i + 1
        write_truncated_binary_coding(w, UInt(digit), nrange)
        # Remove x from remaining
        deleteat!(remaining, idx)
    end
end

# --------------------------
# Core encoder
# --------------------------

"""
    encode_level(w, g, P; params=RCGEParams())

Encode one coarsening level for graph `g` with partition `P` using RCGE.

Inputs:
- w::BitWriter: output bitstream
- g::LightGraphs.AbstractGraph{T<:Unsigned}
- P::Vector{Vector{T}}: list of clusters with global vertex ids
- params::RCGEParams: encoding parameters
"""
function encode_level(w::BitWriter, g::AbstractGraph{T}, P::Vector{Vector{T}}; params::RCGEParams=RCGEParams(), stats::Union{Nothing,RCGEStats}=nothing) where {T<:Unsigned}
    n = nv(g)
    directed = is_directed(g)

    # Sort clusters so local indices match sorted membership order (required for decodability)
    clusters = [sort(copy(C)) for C in P]

    # Precompute membership map: global vertex (Int) -> cluster index
    cluster_of = zeros(Int, n)
    for (ci, C) in enumerate(clusters)
        for u in C
            cluster_of[Int(u)] = ci
        end
    end

    # 1) Cluster membership
    local bits_before = _total_bits(w)
    # write number of clusters
    write_varint(w, T(length(clusters)); encoding=params.varint)
    for C in clusters
        if params.membership == :implicit_ranges
            # Clusters are contiguous ID ranges; store only the size.
            # Decoder reconstructs cluster i as offset+1..offset+|C| from cumulative sizes.
            write_varint(w, T(length(C)); encoding=params.varint)
        else
            Cm = sort(copy(C))
            unique!(Cm)
            if params.membership == :elias_fano
                # Elias-Fano encodes its own length; omit explicit |C| when using EF
                write_elias_fano(w, Cm)
            elseif params.membership == :delta
                # Legacy: explicit length + delta list
                write_varint(w, T(length(Cm)); encoding=params.varint)
                write_gap_coded(w, Cm; encoding=params.gap)
            else
                # STOP-terminated delta list (no explicit length)
                write_stop_delta_list(w, Cm; encoding=params.gap)
            end
        end
    end
    if stats !== nothing
        stats.bits_membership += _total_bits(w) - bits_before
    end

    # 2) Intra-cluster induced graphs
    bits_before = _total_bits(w)
    for C in clusters
        s = length(C)
        # Build local index map for this cluster
        local_index = Dict{T,Int}()
        for (i,u) in enumerate(C)
            local_index[u] = i
        end

        if s <= params.L && !directed
            # Use upper-tri bitset for undirected graphs
            write_bitset_cluster(w, g, C)
        elseif params.intra_mgs
            # Full MGS encoding: inverted-index reference lookup, interval detection,
            # recursive reference — the most powerful encoder available.
            neighbor_lists_mgs = Dict{T, Vector{T}}()
            for (i,u) in enumerate(C)
                nl = T[]
                for v in outneighbors(g, Int(u))
                    if haskey(local_index, v)
                        push!(nl, T(local_index[v]))
                    end
                end
                sort!(nl)
                neighbor_lists_mgs[T(i)] = nl
            end
            Compression.write_compressed_graph_data(w, neighbor_lists_mgs, :children,
                :fibonacci, true, true, true, params.intra_ref_window)
        else
            # Adjacency lists in local order, neighbors restricted to C
            # Optional: try block encoding for the whole cluster and fall back if not beneficial
            if params.intra_block_try
                local use_block = false
                local block_bytes = UInt8[]
                # Build neighbor_lists for this cluster
                neighbor_lists = Dict{T, Vector{T}}()
                for (i,u) in enumerate(C)
                    nl = T[]
                    for v in outneighbors(g, Int(u))
                        if haskey(local_index, v)
                            push!(nl, T(local_index[v]))
                        end
                    end
                    sort!(nl)
                    neighbor_lists[T(i)] = nl
                end
                # Encode block to temp buffer
                local io_blk = IOBuffer(); local w_blk = BitWriter(io_blk)
                Compression.write_compressed_graph_data(w_blk, neighbor_lists, :children, :fibonacci, true, true, true, params.intra_ref_window)
                flush_bitwriter(w_blk; flush_last_bits=true)
                block_bytes = take!(io_blk)
                block_bits = length(block_bytes) * 8
                # Cheap baseline: raw count + delta per vertex
                local io_base = IOBuffer(); local w_base = BitWriter(io_base)
                for u in C
                    nl_int = Int[]
                    for v in outneighbors(g, Int(u))
                        if haskey(local_index, v)
                            push!(nl_int, local_index[v])
                        end
                    end
                    sort!(nl_int)
                    # small-count + delta
                    Compression.write_small_count(w_base, T(length(nl_int)), params.count_varint)
                    if !isempty(nl_int)
                        write_delta(w_base, T.(nl_int), :fibonacci)
                    end
                end
                flush_bitwriter(w_base; flush_last_bits=true)
                baseline_bits = length(take!(io_base)) * 8
                # Decision with small header buffer H=8 bits
                use_block = (block_bits + 8 < baseline_bits)
                # Write cluster mode flag and either block or fall back
                local mode_bits_before = _total_bits(w)
                write_bit(w, use_block)
                if stats !== nothing
                    stats.bits_intra_headers += _total_bits(w) - mode_bits_before
                end
                if use_block
                    # Write block size and bytes
                    local hb0 = _total_bits(w)
                    write_encoded_value(w, T(length(block_bytes)), params.varint)
                    write_bytes(w, block_bytes)
                    if stats !== nothing
                        stats.bits_intra_headers += (_total_bits(w) - hb0) - (length(block_bytes)*8)
                        stats.bits_intra_copy += length(block_bytes) * 8
                    end
                    continue
                end
                # else fall through to per-vertex path
            end
            # Optional intra-cluster referencing to exploit overlap in successive lists
            prev_lists = Vector{Vector{Int}}()
            # First pass: decide refs and compute ref data
            use_ref_vec = Vector{Bool}(undef, s)
            ref_delta_vec = Vector{UInt32}(undef, s)
            ref_positions_list = Vector{Vector{Int}}(undef, s)
            additions_list = Vector{Vector{Int}}(undef, s)
            raw_lists = Vector{Vector{Int}}(undef, s)
            ref_len_list = Vector{Int}(undef, s)
            # Decisions for encoding modes
            pos_use_bitmap_vec = Vector{Bool}(undef, s)
            add_use_intervals_vec = Vector{Bool}(undef, s)
            # Per-vertex mil (greedy search populates this; otherwise fixed)
            mil_vec = fill(params.intra_mil, s)
            _mil_options = [2, 3, 4, 5]
            for (idx_local, u) in enumerate(C)
                # gather local neighbors of u within C, as local indices (1-based)
                nl = Int[]
                for v in outneighbors(g, Int(u))
                    if haskey(local_index, v)
                        push!(nl, local_index[v])
                    end
                end
                sort!(nl)
                # decide reference and mil (greedy by measured bits)
                use_ref = false; ref_positions = Int[]; additions = Int[]; ref_delta_val = UInt32(0)
                if params.intra_greedy_mil
                    # Greedy per-vertex mil search: try all mil values for raw and ref
                    let
                        best_bits = typemax(Int)
                        best_mil_val = params.intra_mil
                        best_is_ref = false
                        best_ref_idx = 0
                        best_ref_pos = Int[]
                        best_ref_add = Int[]

                        # Try raw encoding with each mil
                        _zz_vid = params.intra_zigzag ? T(idx_local) : nothing
                        for mil in _mil_options
                            io_raw = IOBuffer(); w_raw = BitWriter(io_raw)
                            Compression.write_intervals_and_residuals(w_raw, T.(nl), :fibonacci, mil; vertex_id=_zz_vid)
                            flush_bitwriter(w_raw; flush_last_bits=true)
                            raw_bits = length(take!(io_raw)) * 8
                            if raw_bits < best_bits
                                best_bits = raw_bits
                                best_mil_val = mil
                                best_is_ref = false
                            end
                        end

                        # Try ref encoding with each candidate × each mil
                        if params.intra_ref_enabled && !isempty(prev_lists)
                            wstart = max(1, length(prev_lists) - params.intra_ref_window + 1)
                            for rix in wstart:length(prev_lists)
                                ref = prev_lists[rix]
                                # Two-pointer for positions and additions
                                i = 1; j = 1
                                positions = Int[]; adds = Int[]
                                while i <= length(nl) && j <= length(ref)
                                    if nl[i] == ref[j]
                                        push!(positions, j); i += 1; j += 1
                                    elseif nl[i] < ref[j]
                                        push!(adds, nl[i]); i += 1
                                    else
                                        j += 1
                                    end
                                end
                                while i <= length(nl); push!(adds, nl[i]); i += 1; end

                                # Positions cost (doesn't depend on mil)
                                io_pos = IOBuffer(); w_pos = BitWriter(io_pos)
                                if params.intra_copy_blocks
                                    _write_copy_blocks(w_pos, positions, params.varint)
                                else
                                    Compression.write_small_count(w_pos, T(length(positions)), params.count_varint)
                                    if !isempty(positions)
                                        write_delta(w_pos, UInt32.(positions), :fibonacci)
                                    end
                                end
                                flush_bitwriter(w_pos; flush_last_bits=true)
                                pos_bits = length(take!(io_pos)) * 8

                                # Try each mil for additions
                                for mil in _mil_options
                                    io_add = IOBuffer(); w_add = BitWriter(io_add)
                                    Compression.write_intervals_and_residuals(w_add, T.(adds), :fibonacci, mil; vertex_id=_zz_vid)
                                    flush_bitwriter(w_add; flush_last_bits=true)
                                    add_bits = length(take!(io_add)) * 8

                                    ref_bits = pos_bits + add_bits
                                    if ref_bits < best_bits
                                        best_bits = ref_bits
                                        best_mil_val = mil
                                        best_is_ref = true
                                        best_ref_idx = rix
                                        best_ref_pos = positions
                                        best_ref_add = adds
                                    end
                                end
                            end
                        end

                        if best_is_ref
                            use_ref = true
                            ref_positions = best_ref_pos
                            additions = best_ref_add
                            ref_delta_val = UInt32(idx_local - best_ref_idx)
                        end
                        mil_vec[idx_local] = best_mil_val
                    end
                elseif params.intra_ref_enabled && !isempty(prev_lists)
                    # Original reference decision (with or without intervals)
                    let
                        _zz_vid = params.intra_zigzag ? T(idx_local) : nothing
                        io_raw = IOBuffer(); w_raw = BitWriter(io_raw)
                        if params.intra_intervals
                            Compression.write_intervals_and_residuals(w_raw, T.(nl), :fibonacci, params.intra_mil; vertex_id=_zz_vid)
                        elseif params.intra_stop_deltas
                            _write_stop_delta_zigzag(w_raw, T.(nl), :fibonacci, _zz_vid)
                        else
                            Compression.write_small_count(w_raw, T(length(nl)), params.count_varint)
                            if !isempty(nl)
                                write_delta(w_raw, T.(nl), :fibonacci; vertex_id=_zz_vid)
                            end
                        end
                        flush_bitwriter(w_raw; flush_last_bits=true)
                        raw_bits = length(take!(io_raw)) * 8
                        best_bits = raw_bits
                        best_idx = 0
                        best_pos = Int[]
                        best_add = Int[]
                        wstart = max(1, length(prev_lists) - params.intra_ref_window + 1)
                        for rix in wstart:length(prev_lists)
                            ref = prev_lists[rix]
                            # build positions and additions by two-pointer
                            i = 1; j = 1
                            positions = Int[]; adds = Int[]
                            while i <= length(nl) && j <= length(ref)
                                if nl[i] == ref[j]
                                    push!(positions, j); i += 1; j += 1
                                elseif nl[i] < ref[j]
                                    push!(adds, nl[i]); i += 1
                                else
                                    j += 1
                                end
                            end
                            while i <= length(nl); push!(adds, nl[i]); i += 1; end
                            # estimate ref bits
                            io_ref = IOBuffer(); w_ref = BitWriter(io_ref)
                            # positions: copy-blocks or small-count + delta (NO zigzag — positions are indices into ref list)
                            if params.intra_copy_blocks
                                _write_copy_blocks(w_ref, positions, params.varint)
                            else
                                Compression.write_small_count(w_ref, T(length(positions)), params.count_varint)
                                if !isempty(positions)
                                    write_delta(w_ref, UInt32.(positions), :fibonacci)
                                end
                            end
                            # additions: interval-aware or plain delta
                            if params.intra_intervals
                                Compression.write_intervals_and_residuals(w_ref, T.(adds), :fibonacci, params.intra_mil; vertex_id=_zz_vid)
                            elseif params.intra_stop_deltas
                                _write_stop_delta_zigzag(w_ref, T.(adds), :fibonacci, _zz_vid)
                            else
                                Compression.write_small_count(w_ref, T(length(adds)), params.count_varint)
                                if !isempty(adds)
                                    write_delta(w_ref, T.(adds), :fibonacci; vertex_id=_zz_vid)
                                end
                            end
                            flush_bitwriter(w_ref; flush_last_bits=true)
                            ref_bits = length(take!(io_ref)) * 8
                            if ref_bits < best_bits
                                best_bits = ref_bits
                                best_idx = rix
                                best_pos = positions
                                best_add = adds
                            end
                        end
                        if best_idx > 0 && best_bits < raw_bits
                            use_ref = true
                            ref_positions = best_pos
                            additions = best_add
                            ref_delta_val = UInt32(idx_local - best_idx)
                        end
                    end
                end
                use_ref_vec[idx_local] = use_ref
                ref_delta_vec[idx_local] = ref_delta_val
                ref_positions_list[idx_local] = ref_positions
                additions_list[idx_local] = additions
                raw_lists[idx_local] = nl
                if use_ref
                    # store reference list length for bitmap construction
                    # best_ref_idx chosen above
                    # recompute best_ref_idx for storing ref_len
                    # Note: best_ref_idx is the last matching index in loop; recompute quickly
                    # Use nearest previous list matching ref_delta
                    ref_index = idx_local - Int(ref_delta_val)
                    ref_len_list[idx_local] = ref_index >= 1 ? length(prev_lists[ref_index]) : 0
                else
                    ref_len_list[idx_local] = 0
                end
                push!(prev_lists, nl)
            end

            # Write per-cluster ref bitmap and ref deltas
            if params.intra_ref_enabled
                local hb0 = _total_bits(w)
                if params.intra_ref_fixwidth
                    # Fixed-width (or VLC): 1-bit flag + delta encoding per ref vertex
                    _write_fixwidth_ref_deltas(w, use_ref_vec, ref_delta_vec, params.intra_ref_window; vlc=params.intra_ref_vlc)
                else
                    # Legacy: byte-padded bitmap + varint delta list
                    Compression.write_bitpacked_bitmap(w, use_ref_vec)
                    if any(use_ref_vec)
                        deltas = UInt32[]
                        for idx in 1:s
                            if use_ref_vec[idx]
                                push!(deltas, ref_delta_vec[idx])
                            end
                        end
                        if params.intra_ref_rle
                            Compression.write_rle_ones_deltas(w, deltas, params.varint)
                        else
                            write_encoded_value(w, T(length(deltas)), params.varint)
                            for d in deltas
                                write_varint(w, d; encoding=params.varint)
                            end
                        end
                    end
                end
                if stats !== nothing
                    stats.bits_intra_headers += _total_bits(w) - hb0
                end
            end

            # Write per-vertex mil tags if greedy mil search was used
            if params.intra_greedy_mil
                local mil_hdr0 = _total_bits(w)
                # Find default mil (most common)
                mil_counts = Dict{Int,Int}()
                for idx in 1:s
                    mil_counts[mil_vec[idx]] = get(mil_counts, mil_vec[idx], 0) + 1
                end
                default_mil = argmax(mil_counts)
                # Write default mil as 2-bit tag (index into [2,3,4,5])
                default_mil_idx = findfirst(==(default_mil), _mil_options)
                default_mil_idx === nothing && (default_mil_idx = 3)  # fallback to mil=4
                write_bit(w, ((default_mil_idx - 1) >> 1) & 1 == 1)
                write_bit(w, (default_mil_idx - 1) & 1 == 1)
                # Per-vertex: 1-bit flag (0=default, 1=custom+2-bit tag)
                for idx in 1:s
                    if mil_vec[idx] == default_mil
                        write_bit(w, false)
                    else
                        write_bit(w, true)
                        local mi = findfirst(==(mil_vec[idx]), _mil_options)
                        mi === nothing && (mi = 3)
                        write_bit(w, ((mi - 1) >> 1) & 1 == 1)
                        write_bit(w, (mi - 1) & 1 == 1)
                    end
                end
                if stats !== nothing
                    stats.bits_intra_headers += _total_bits(w) - mil_hdr0
                end
            end

            # Now write per-vertex payloads without flags or ref_deltas
            for idx_local in 1:s
                _zz_vid = params.intra_zigzag ? T(idx_local) : nothing

                if use_ref_vec[idx_local]
                    # write copied positions into reference list
                    ref_positions = ref_positions_list[idx_local]
                    additions = additions_list[idx_local]
                    ref_len = ref_len_list[idx_local]
                    local bpos0 = _total_bits(w)
                    if params.intra_copy_blocks
                        _write_copy_blocks(w, ref_positions, params.varint)
                    elseif ref_len > 0
                        copied = fill(false, ref_len)
                        for p in ref_positions
                            if 1 <= p <= ref_len; copied[p] = true; end
                        end
                        Compression.write_bitmap_rle_ones_deltas(w, copied, params.varint)
                    else
                        # empty token list via small count
                        Compression.write_small_count(w, T(0), params.count_varint)
                    end
                    local bpos1 = _total_bits(w)
                    # additions: MGS intervals, custom intervals, or plain delta
                    local ah0 = _total_bits(w)
                    if params.intra_intervals || params.intra_greedy_mil
                        Compression.write_intervals_and_residuals(w, T.(additions), :fibonacci, mil_vec[idx_local]; vertex_id=_zz_vid)
                    elseif params.additions_mode == :intervals
                        # detect runs and write intervals + singles
                        runs = Vector{Tuple{Int,Int}}(); singles = Int[]
                        if !isempty(additions)
                            i3 = 1
                            while i3 <= length(additions)
                                j3 = i3
                                while j3+1 <= length(additions) && additions[j3+1] == additions[j3] + 1
                                    j3 += 1
                                end
                                lenr = j3 - i3 + 1
                                if lenr >= 2
                                    push!(runs, (additions[i3], lenr))
                                else
                                    push!(singles, additions[i3])
                                end
                                i3 = j3 + 1
                            end
                        end
                        Compression.write_small_count(w, T(length(runs)), params.count_varint)
                        for (st, ln) in runs
                            write_encoded_value(w, T(st), params.varint)
                            write_encoded_value(w, T(ln), params.varint)
                        end
                        Compression.write_small_count(w, T(length(singles)), params.count_varint)
                        if !isempty(singles)
                            write_delta(w, T.(singles), :fibonacci; vertex_id=_zz_vid)
                        end
                    elseif params.intra_add_adaptive && params.intra_stop_deltas
                        # Adaptive: per-vertex pick cheaper of STOP-delta vs intervals
                        # Write mode bit inline (0=stop-delta, 1=intervals), then payload
                        local adds_T = T.(additions)
                        local _io_sd = IOBuffer(); local _w_sd = BitWriter(_io_sd)
                        _write_stop_delta_zigzag(_w_sd, adds_T, :fibonacci, _zz_vid)
                        flush_bitwriter(_w_sd; flush_last_bits=true)
                        local _sd_bits = length(take!(_io_sd)) * 8
                        local _io_iv = IOBuffer(); local _w_iv = BitWriter(_io_iv)
                        Compression.write_intervals_and_residuals(_w_iv, adds_T, :fibonacci, params.intra_adapt_mil; vertex_id=_zz_vid)
                        flush_bitwriter(_w_iv; flush_last_bits=true)
                        local _iv_bits = length(take!(_io_iv)) * 8
                        if _iv_bits + 1 < _sd_bits   # +1 for the mode flag itself
                            write_bit(w, true)   # intervals
                            Compression.write_intervals_and_residuals(w, adds_T, :fibonacci, params.intra_adapt_mil; vertex_id=_zz_vid)
                        else
                            write_bit(w, false)  # stop-delta
                            _write_stop_delta_zigzag(w, adds_T, :fibonacci, _zz_vid)
                        end
                    else
                        if params.intra_stop_deltas
                            _write_stop_delta_zigzag(w, T.(additions), :fibonacci, _zz_vid)
                        else
                            Compression.write_small_count(w, T(length(additions)), params.count_varint)
                            if !isempty(additions)
                                write_delta(w, T.(additions), :fibonacci; vertex_id=_zz_vid)
                            end
                        end
                    end
                    local b3 = _total_bits(w)
                    if stats !== nothing
                        # positions bitmap is copy payload; additions include intervals/singles payload
                        stats.bits_intra_copy += (bpos1 - bpos0)
                        stats.bits_intra_add += (b3 - bpos1)
                        stats.intra_ref_used += 1
                    end
                else
                    # raw
                    nl = raw_lists[idx_local]
                    local rb0 = _total_bits(w)
                    if params.intra_intervals || params.intra_greedy_mil
                        Compression.write_intervals_and_residuals(w, T.(nl), :fibonacci, mil_vec[idx_local]; vertex_id=_zz_vid)
                    elseif params.intra_raw_adaptive && params.intra_stop_deltas
                        # Adaptive: per-vertex pick cheaper of STOP-delta vs intervals
                        local nl_T = T.(nl)
                        local _io_sd2 = IOBuffer(); local _w_sd2 = BitWriter(_io_sd2)
                        _write_stop_delta_zigzag(_w_sd2, nl_T, :fibonacci, _zz_vid)
                        flush_bitwriter(_w_sd2; flush_last_bits=true)
                        local _sd_bits2 = length(take!(_io_sd2)) * 8
                        local _io_iv2 = IOBuffer(); local _w_iv2 = BitWriter(_io_iv2)
                        Compression.write_intervals_and_residuals(_w_iv2, nl_T, :fibonacci, params.intra_adapt_mil; vertex_id=_zz_vid)
                        flush_bitwriter(_w_iv2; flush_last_bits=true)
                        local _iv_bits2 = length(take!(_io_iv2)) * 8
                        if _iv_bits2 + 1 < _sd_bits2
                            write_bit(w, true)   # intervals
                            Compression.write_intervals_and_residuals(w, nl_T, :fibonacci, params.intra_adapt_mil; vertex_id=_zz_vid)
                        else
                            write_bit(w, false)  # stop-delta
                            _write_stop_delta_zigzag(w, nl_T, :fibonacci, _zz_vid)
                        end
                    elseif params.intra_stop_deltas
                        _write_stop_delta_zigzag(w, T.(nl), :fibonacci, _zz_vid)
                    else
                        # raw count: use small-count which escapes to varint for large values and handles zero
                        Compression.write_small_count(w, T(length(nl)), params.count_varint)
                        if !isempty(nl)
                            write_delta(w, T.(nl), :fibonacci; vertex_id=_zz_vid)
                        end
                    end
                    if stats !== nothing
                        stats.bits_intra_raw += _total_bits(w) - rb0
                        stats.intra_no_ref += 1
                    end
                end
            end
        end
    end
    if stats !== nothing
        stats.bits_intra += _total_bits(w) - bits_before
    end

    # 3) Inter-cluster bundles via stubs (optimized: only present pairs)
    mclusters = length(clusters)
    # Precompute local index maps for all clusters
    local_index_by_cluster = Vector{Dict{T,Int}}(undef, mclusters)
    for ci in 1:mclusters
        C = clusters[ci]
        li = Dict{T,Int}()
        for (i,u) in enumerate(C)
            li[u] = i
        end
        local_index_by_cluster[ci] = li
    end

    # Build mapping from (ia,ib) -> neighbors_in_B dict: u_in_A => [B-local neighbors]
    pair_map = Dict{Tuple{Int,Int}, Dict{T,Vector{Int}}}()
    for u_int in 1:n
        ia = cluster_of[u_int]
        # Skip vertices not in any cluster (shouldn't happen)
        ia == 0 && continue
        for v in outneighbors(g, u_int)
            ib = cluster_of[Int(v)]
            if ib == 0 || ia == ib
                continue
            end
            key = (ia, ib)
            dict_u = get!(pair_map, key, Dict{T,Vector{Int}}())
            # B-local index for v
            vloc = get(local_index_by_cluster[ib], v, 0)
            if vloc == 0
                continue
            end
            # Append to list for u
            uT = T(u_int)
            lst = get!(dict_u, uT, Int[])
            push!(lst, vloc)
        end
    end

    # Group by source cluster A and emit per-A neighbor lists
    # Build mapping A -> sorted list of neighbor Bs
    neighbors_by_A = Dict{Int, Vector{Int}}()
    for (iaib, _) in pair_map
        ia, ib = iaib
        if directed || !params.undirected_pairs || ia < ib
            push!(get!(neighbors_by_A, ia, Int[]), ib)
        end
    end
    for ia in keys(neighbors_by_A)
        unique!(neighbors_by_A[ia])
        sort!(neighbors_by_A[ia])
    end

    # Emit per-A targets and AB groups in fixed A order (1..K)
    for ia in 1:mclusters
        A = clusters[ia]
        Bs = get(neighbors_by_A, ia, Int[])
        sort!(Bs)
        # A's targets: STOP-terminated delta list of B indices
        bits_before = _total_bits(w)
        write_stop_delta_list(w, T.(Bs); encoding=params.gap)
        if stats !== nothing
            stats.bits_inter_headers += _total_bits(w) - bits_before
        end

        # For each B in targets, write AB group as STOP-terminated list of active u records
        for ib in Bs
            neighbors_in_B = pair_map[(ia,ib)]
            # Ensure neighbor lists are sorted and strictly increasing (B-local ids)
            for lst in values(neighbors_in_B)
                sort!(lst)
                unique!(lst)
            end
            # AB group payload (counted as lists): per active u in ascending local order
            pl0 = _total_bits(w)
            for (idxA, u) in enumerate(A)
                Ns = get(neighbors_in_B, u, Int[])
                if isempty(Ns)
                    continue
                end
                # Record prelude '1', then u_local as varint (non-zero), then STOP-terminated vlist (delta)
                write_bit(w, true)
                write_encoded_value(w, T(idxA), params.varint)
                write_stop_delta_list(w, T.(Ns); encoding=params.gap)
            end
            # End of AB group: STOP bit '0'
            write_bit(w, false)
            if stats !== nothing
                stats.bits_inter_lists += _total_bits(w) - pl0
            end
        end
    end

    return nothing
end

# --------------------------
# Decoder
# --------------------------

"""
    read_stop_delta_list(r::BitReader; encoding::Symbol=:fibonacci, T::Type{<:Unsigned}=UInt32)

Read a STOP-terminated delta-coded list (inverse of `write_stop_delta_list`).
Each element is prefixed by a '1' bit; a '0' bit terminates the list.
"""
function read_stop_delta_list(r::BitReader; encoding::Symbol=:fibonacci, T::Type{<:Unsigned}=UInt32)
    vals = T[]
    prev = zero(T)
    while read_bit(r)  # 1 = more values, 0 = STOP
        delta = read_encoded_value(r, encoding, T)
        prev += delta
        push!(vals, prev)
    end
    return vals
end

"""
    decode_level(r::BitReader, params::RCGEParams; T::Type{<:Unsigned}=UInt32, directed::Bool=true)

Decode one RCGE coarsening level from the bitstream.
Returns a `Dict{T, Vector{T}}` mapping global vertex ID → sorted outneighbors.
"""
function decode_level(r::BitReader, params::RCGEParams; T::Type{<:Unsigned}=UInt32, directed::Bool=true)
    # ----------------------------------------------------------------
    # Section 1: Read cluster membership
    # ----------------------------------------------------------------
    K = Int(read_encoded_value(r, params.varint, T))
    clusters = Vector{Vector{T}}(undef, K)

    if params.membership == :implicit_ranges
        # Clusters are contiguous ID ranges: cluster i = offset+1..offset+size_i
        offset = T(0)
        for ci in 1:K
            sz = Int(read_encoded_value(r, params.varint, T))
            clusters[ci] = collect(offset + T(1) : offset + T(sz))
            offset += T(sz)
        end
    else
        for ci in 1:K
            if params.membership == :elias_fano
                clusters[ci] = read_elias_fano(r, T)
            elseif params.membership == :delta
                len = Int(read_encoded_value(r, params.varint, T))
                clusters[ci] = read_delta(r, params.gap, T; max_elements=len)
            else  # :stop
                clusters[ci] = read_stop_delta_list(r; encoding=params.gap, T=T)
            end
        end
    end

    # Build global neighbor dict
    neighbor_lists = Dict{T, Vector{T}}()

    # ----------------------------------------------------------------
    # Section 2: Read intra-cluster edges
    # ----------------------------------------------------------------
    _mil_options = [2, 3, 4, 5]

    for ci in 1:K
        C = clusters[ci]
        s = length(C)

        if s <= params.L && !directed
            # Upper-tri bitset: s*(s-1)/2 bits
            for i in 1:s-1
                for j in (i+1):s
                    if read_bit(r)
                        u_global = C[i]
                        v_global = C[j]
                        push!(get!(neighbor_lists, u_global, T[]), v_global)
                        push!(get!(neighbor_lists, v_global, T[]), u_global)
                    end
                end
            end
        elseif params.intra_mgs
            # Full MGS block: read_compressed_graph_data needs vertex count
            local_neighbors = read_compressed_graph_data(r, T(s), :children, :fibonacci, T)
            for (local_v, nbs) in local_neighbors
                u_global = C[Int(local_v)]
                for nb in nbs
                    push!(get!(neighbor_lists, u_global, T[]), C[Int(nb)])
                end
            end
        else
            # Per-vertex adjacency list path
            # Handle optional block encoding
            if params.intra_block_try
                use_block = read_bit(r)
                if use_block
                    block_len = Int(read_encoded_value(r, params.varint, T))
                    # Read block_len bytes and decode as MGS compressed graph
                    block_bytes = Vector{UInt8}(undef, block_len)
                    for bi in 1:block_len
                        block_bytes[bi] = UInt8(read_value(r, 8, UInt8))
                    end
                    block_io = IOBuffer(block_bytes)
                    block_r = BitReader(block_io)
                    local_neighbors = read_compressed_graph_data(block_r, T(s), :children, :fibonacci, T)
                    for (local_v, nbs) in local_neighbors
                        u_global = C[Int(local_v)]
                        for nb in nbs
                            push!(get!(neighbor_lists, u_global, T[]), C[Int(nb)])
                        end
                    end
                    continue  # skip per-vertex path for this cluster
                end
                # else fall through to per-vertex path
            end

            # Read reference bitmap and deltas if ref is enabled
            use_ref_vec = fill(false, s)
            ref_delta_vec = zeros(UInt32, s)
            has_any_ref = false

            if params.intra_ref_enabled
                if params.intra_ref_fixwidth
                    # Fixed-width (or VLC): 1-bit flag + delta encoding per ref vertex
                    use_ref_vec, ref_delta_vec = _read_fixwidth_ref_deltas(r, s, params.intra_ref_window; vlc=params.intra_ref_vlc)
                    has_any_ref = any(use_ref_vec)
                else
                    # Legacy: byte-padded bitmap + varint delta list
                    use_ref_vec = Vector{Bool}(read_bitpacked_bitmap(r, s))
                    has_any_ref = any(use_ref_vec)
                    if has_any_ref
                        if params.intra_ref_rle
                            ref_deltas = read_rle_ones_deltas(r, params.varint, UInt32)
                        else
                            ndelt = Int(read_encoded_value(r, params.varint, T))
                            ref_deltas = UInt32[]
                            for _ in 1:ndelt
                                push!(ref_deltas, read_encoded_value(r, params.varint, UInt32))
                            end
                        end
                        di = 1
                        for idx in 1:s
                            if use_ref_vec[idx]
                                ref_delta_vec[idx] = ref_deltas[di]
                                di += 1
                            end
                        end
                    end
                end
            end

            # Read per-vertex mil tags if greedy mil search was used
            mil_vec = fill(params.intra_mil, s)
            if params.intra_greedy_mil
                # Read default mil as 2-bit tag
                b1 = read_bit(r)
                b2 = read_bit(r)
                default_mil_idx = ((Int(b1) << 1) | Int(b2)) + 1
                default_mil = _mil_options[default_mil_idx]
                fill!(mil_vec, default_mil)
                # Per-vertex: 1-bit flag (0=default, 1=custom+2-bit tag)
                for idx in 1:s
                    if read_bit(r)
                        mb1 = read_bit(r)
                        mb2 = read_bit(r)
                        mi = ((Int(mb1) << 1) | Int(mb2)) + 1
                        mil_vec[idx] = _mil_options[mi]
                    end
                end
            end

            # Read per-vertex payloads
            prev_lists = Vector{Vector{Int}}()  # local neighbor lists for reference lookups
            for idx_local in 1:s
                local nl_local::Vector{Int}
                _zz_vid = params.intra_zigzag ? T(idx_local) : nothing

                if use_ref_vec[idx_local]
                    # Reference mode: read positions bitmap then additions
                    ref_index = idx_local - Int(ref_delta_vec[idx_local])
                    ref_list = ref_index >= 1 ? prev_lists[ref_index] : Int[]
                    ref_len = length(ref_list)

                    # Read copied positions from reference
                    if params.intra_copy_blocks
                        copied_positions = _read_copy_blocks(r, params.varint, T)
                        copied_vals = Int[]
                        for p in copied_positions
                            if 1 <= p <= ref_len
                                push!(copied_vals, ref_list[p])
                            end
                        end
                    elseif ref_len > 0
                        copied_bitmap = read_bitmap_rle_ones_deltas(r, params.varint, UInt32)
                        # Extract copied positions from ref
                        copied_vals = Int[]
                        for (pi, flag) in enumerate(copied_bitmap)
                            if flag && pi <= ref_len
                                push!(copied_vals, ref_list[pi])
                            end
                        end
                    else
                        # Empty ref — read small count of 0
                        read_small_count(r, params.count_varint, T)
                        copied_vals = Int[]
                    end

                    # Read additions
                    local additions::Vector{Int}
                    if params.intra_intervals || params.intra_greedy_mil
                        additions = Int.(read_intervals_and_residuals(r, :fibonacci, mil_vec[idx_local], T; vertex_id=_zz_vid))
                    elseif params.additions_mode == :intervals
                        # intervals: runs + singles
                        n_runs = Int(read_small_count(r, params.count_varint, T))
                        add_vals = Int[]
                        for _ in 1:n_runs
                            st = Int(read_encoded_value(r, params.varint, T))
                            ln = Int(read_encoded_value(r, params.varint, T))
                            for k in 0:(ln-1)
                                push!(add_vals, st + k)
                            end
                        end
                        n_singles = Int(read_small_count(r, params.count_varint, T))
                        if n_singles > 0
                            singles = Int.(read_delta(r, :fibonacci, T; max_elements=n_singles, vertex_id=_zz_vid))
                            append!(add_vals, singles)
                        end
                        additions = add_vals
                    elseif params.intra_add_adaptive && params.intra_stop_deltas
                        # Adaptive: read mode bit then decode accordingly
                        if read_bit(r)  # true = intervals
                            additions = Int.(read_intervals_and_residuals(r, :fibonacci, params.intra_adapt_mil, T; vertex_id=_zz_vid))
                        else            # false = stop-delta
                            additions = Int.(_read_stop_delta_zigzag(r, :fibonacci, T, _zz_vid))
                        end
                    elseif params.intra_stop_deltas
                        additions = Int.(_read_stop_delta_zigzag(r, :fibonacci, T, _zz_vid))
                    else
                        n_add = Int(read_small_count(r, params.count_varint, T))
                        if n_add > 0
                            additions = Int.(read_delta(r, :fibonacci, T; max_elements=n_add, vertex_id=_zz_vid))
                        else
                            additions = Int[]
                        end
                    end

                    # Combine copied + additions
                    nl_local = sort(vcat(copied_vals, additions))
                else
                    # Raw mode
                    if params.intra_intervals || params.intra_greedy_mil
                        nl_local = Int.(read_intervals_and_residuals(r, :fibonacci, mil_vec[idx_local], T; vertex_id=_zz_vid))
                    elseif params.intra_raw_adaptive && params.intra_stop_deltas
                        # Adaptive: read mode bit then decode accordingly
                        if read_bit(r)  # true = intervals
                            nl_local = Int.(read_intervals_and_residuals(r, :fibonacci, params.intra_adapt_mil, T; vertex_id=_zz_vid))
                        else            # false = stop-delta
                            nl_local = Int.(_read_stop_delta_zigzag(r, :fibonacci, T, _zz_vid))
                        end
                    elseif params.intra_stop_deltas
                        nl_local = Int.(_read_stop_delta_zigzag(r, :fibonacci, T, _zz_vid))
                    else
                        cnt = Int(read_small_count(r, params.count_varint, T))
                        if cnt > 0
                            nl_local = Int.(read_delta(r, :fibonacci, T; max_elements=cnt, vertex_id=_zz_vid))
                        else
                            nl_local = Int[]
                        end
                    end
                end

                push!(prev_lists, nl_local)

                # Map local indices to global vertex IDs
                u_global = C[idx_local]
                for nb_local in nl_local
                    if 1 <= nb_local <= s
                        push!(get!(neighbor_lists, u_global, T[]), C[nb_local])
                    end
                end
            end
        end
    end

    # ----------------------------------------------------------------
    # Section 3: Read inter-cluster edges
    # ----------------------------------------------------------------
    for ia in 1:K
        A = clusters[ia]
        # Read STOP-terminated delta list of target cluster indices
        Bs = read_stop_delta_list(r; encoding=params.gap, T=T)

        for ib_T in Bs
            ib = Int(ib_T)
            B = clusters[ib]
            sB = length(B)

            # Read AB group: STOP-terminated list of (u_local, neighbor_list) records
            while read_bit(r)  # 1 = more records, 0 = end of group
                idxA = Int(read_encoded_value(r, params.varint, T))
                u_global = A[idxA]
                # Read STOP-terminated delta list of B-local neighbor indices
                Ns_local = read_stop_delta_list(r; encoding=params.gap, T=T)
                for nb_local_T in Ns_local
                    nb_local = Int(nb_local_T)
                    if 1 <= nb_local <= sB
                        push!(get!(neighbor_lists, u_global, T[]), B[nb_local])
                    end
                end
            end
        end
    end

    # Sort and deduplicate all neighbor lists
    for (v, nbs) in neighbor_lists
        sort!(nbs)
        unique!(nbs)
    end

    return neighbor_lists
end

"""
    load_rcge_graph(filepath::AbstractString; params::RCGEParams=RCGEParams(),
                    T::Type{<:Unsigned}=UInt32, directed::Bool=true)

Load an RCGE-compressed graph from a file, decoding all levels.
Returns a `Dict{T, Vector{T}}` mapping vertex → sorted outneighbors.
"""
function load_rcge_graph(filepath::AbstractString; params::RCGEParams=RCGEParams(),
                         T::Type{<:Unsigned}=UInt32, directed::Bool=true)
    data = read(filepath)
    io = IOBuffer(data)
    r = BitReader(io)
    neighbor_lists = Dict{T, Vector{T}}()
    while r.index <= r.length
        level_neighbors = decode_level(r, params; T=T, directed=directed)
        # Merge level edges
        for (v, nbs) in level_neighbors
            existing = get!(neighbor_lists, v, T[])
            append!(existing, nbs)
            sort!(existing)
            unique!(existing)
        end
    end
    return neighbor_lists
end

end # module RCGE
