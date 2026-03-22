#
# Adjacently: Julia Complex Directed Networks Library
# Copyright (C) 2016-2026 Jimmy Dubuisson <jimmy@dubuisson.ch>
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
# CG: Clustered Greedy

module CG

using LightGraphs
using LightGraphs: AbstractGraph, outneighbors, nv, is_directed

using ...IO: BitWriter, BitReader, write_bit, write_bits, write_value, flush_bitwriter,
    read_bit, read_bits, read_value, write_bytes, reset_bitwriter!, bytes_written, get_bytes
import ..Compression
using ..Compression: write_encoded_value, write_delta, write_truncated_binary_coding,
    write_hybrid_mix_encoded_list, delta_encode_vector, write_elias_fano,
    read_encoded_value, read_delta, read_elias_fano,
    read_small_count, read_rle_ones_deltas, read_bitmap_rle_ones_deltas,
    read_bitpacked_bitmap, read_intervals_and_residuals,
    read_compressed_graph_data, compress_intervals,
    estimate_encoded_value_cost, estimate_interval_runlength_encoding_cost

export CGParams, CGStats, encode_level, decode_level, load_cg_graph

"""
    CGStats

Lightweight counters (in bits) for encode_level sections.
"""
mutable struct CGStats
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
    # Profiling: copy-mode selection counts (when intra_copy_adaptive=true)
    intra_copy_bitmap_count::Int
    intra_copy_blocks_count::Int
    intra_copy_complement_count::Int
    # Profiling: overlap fraction histogram (10 buckets: [0%,10%), …, [90%,100%])
    intra_overlap_hist::Vector{Int}
    # Profiling: total additions across all ref vertices
    intra_add_count_total::Int
    CGStats() = new(0,0,0,0,0,0,0,0,0,0,0,0,0,
                      Tuple{Int,Int,Int,Int,Int,Int}[],
                      0, 0, 0, zeros(Int, 10), 0)
end

@inline function _total_bits(w::BitWriter)
    return w.bit_count
end

"""
    _reset!(w::BitWriter)

Reset a BitWriter for reuse without allocation.
"""
@inline function _reset!(w::BitWriter)
    reset_bitwriter!(w)
end

"""
    CostBuffer{T}

Pre-allocated buffers for cost estimation in the reference search loop.
Uses buffer-only BitWriters (no IO target) for fast reset.
Parameterized on T to hold neighbor lists directly as Vector{T}.
"""
mutable struct CostBuffer{T<:Unsigned}
    w1::BitWriter   # positions: copy-blocks / raw encoding
    w2::BitWriter   # complement / stop-delta / additions
    w3::BitWriter   # intervals / second trial
    positions::Vector{T}
    adds::Vector{T}
    skipped::Vector{T}
end

function CostBuffer{T}() where {T<:Unsigned}
    CostBuffer{T}(BitWriter(capacity=1024), BitWriter(capacity=1024), BitWriter(capacity=1024),
               T[], T[], T[])
end

"""
    CGParams(; L=32, varint=:fibonacci, count_varint=:fibonacci, gap=:fibonacci, degree=:golomb, perm_strategy=:lehmer, undirected_pairs=true, membership=:delta, inter_strategy=:perm,
                 intra_ref_enabled::Bool=true, intra_ref_window::Int=16, intra_ref_rle::Bool=true,
                 intra_block_try::Bool=false)

Parameters for CG encoding.
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
Base.@kwdef struct CGParams
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
    intra_copy_adaptive::Bool = false # Per-ref-vertex: pick cheaper of copy-blocks vs dense bitmap over reference degree (requires intra_copy_blocks=true)
    intra_ref_fixwidth::Bool = false  # Fixed-width ref delta encoding (1-bit flag + ceil(log2(window)) bits per ref)
    intra_ref_vlc::Bool = false       # VLC (Fibonacci) for ref delta instead of fixed-width (requires intra_ref_fixwidth=true)
    intra_add_adaptive::Bool = false  # Per-vertex adaptive additions: pick cheaper of STOP-delta vs intervals (requires intra_stop_deltas=true)
    intra_raw_adaptive::Bool = false  # Per-vertex adaptive raw: pick cheaper of STOP-delta vs intervals (requires intra_stop_deltas=true)
    intra_adapt_mil::Int = 2          # MIL for adaptive interval encoding (2=most aggressive)
    intra_lr_split::Bool = false       # Left/right residual split: split residuals at vertex_id after interval extraction
    intra_tight_deltas::Bool = false  # Skip +1 shift on delta gaps in LR residuals (gaps are always ≥ 1 for sorted unique lists)
    cost_model::Int = Compression.DEFAULT_COST_MODEL  # 0=full analytical, 1=fast (skip interval/LR, pure delta cost)
    index_sample_k::Int = 0  # Two-level sampled offsets: 0=full, >0=sample every k-th vertex (must be multiple of 4)
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
        write_encoded_value(w, UInt64(Compression._zigzag_encode(offset) + 1), encoding)
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
        if first && vertex_id !== nothing
            raw64 = read_encoded_value(r, encoding, UInt64)
            val = T(Int64(vertex_id) + Compression._zigzag_decode(raw64 - 1))
            first = false
        else
            raw = read_encoded_value(r, encoding, T)
            val = prev + raw
            first = false
        end
        push!(result, val)
        prev = val
    end
    return result
end

"""
    _write_ir_lr(w, neighbors, encoding, mil, vertex_id)

Write neighbor list using interval+residual encoding with left/right residual split.
Intervals are extracted from the full list, then residuals are split at vertex_id
into left (<vertex_id) and right (>vertex_id) halves. Each half's residuals are
transformed to ascending distances from vertex_id and delta-encoded, yielding
smaller first values than zigzag encoding of the full list.

Format: intervals (zigzag first start) | num_residuals | left_count | left_dists | right_dists
"""
function _write_ir_lr(w::BitWriter, neighbors::Vector{T}, encoding::Symbol,
                       mil::Int, vertex_id; tight_deltas::Bool=false) where {T<:Unsigned}
    if vertex_id === nothing || isempty(neighbors)
        Compression.write_intervals_and_residuals(w, neighbors, encoding, mil; vertex_id=vertex_id)
        return
    end
    vid = T(vertex_id)

    # Extract intervals from full sorted neighbor list
    intervals, residuals = Compression.compress_intervals(neighbors, mil)

    # Write intervals: same format as standard encoding with zigzag first start
    write_encoded_value(w, T(length(intervals)) + T(1), encoding)
    if !isempty(intervals)
        prev_ref = T(0)
        for (idx, (start, len)) in enumerate(intervals)
            if idx == 1
                offset = Int64(start) - Int64(vid)
                write_encoded_value(w, UInt64(Compression._zigzag_encode(offset) + 1), encoding)
            else
                write_encoded_value(w, start - prev_ref, encoding)
            end
            write_encoded_value(w, len - T(mil) + T(1), encoding)
            prev_ref = start
        end
    end

    # Write total residual count (same as standard)
    write_encoded_value(w, T(length(residuals)) + T(1), encoding)
    if !isempty(residuals)
        # Split residuals at vertex_id: left (< vid), right (> vid)
        # Note: val == vid (self-loop) goes to right with distance 0+1 shift
        split_idx = 0
        @inbounds for i in 1:length(residuals)
            if residuals[i] < vid
                split_idx = i
            else
                break
            end
        end
        n_left = split_idx
        n_right = length(residuals) - n_left

        # Write left count (decoder derives right = total - left)
        write_encoded_value(w, T(n_left) + T(1), encoding)

        # Left residuals: reverse to ascending distances from vertex_id
        if n_left > 0
            left_dists = Vector{T}(undef, n_left)
            @inbounds for i in 1:n_left
                left_dists[i] = vid - residuals[n_left - i + 1]
            end
            write_delta(w, left_dists, encoding; positive_gaps=tight_deltas)
        end

        # Right residuals: distances from vertex_id (+1 to handle val==vid edge case)
        if n_right > 0
            right_dists = Vector{T}(undef, n_right)
            @inbounds for i in 1:n_right
                right_dists[i] = residuals[n_left + i] - vid + T(1)
            end
            write_delta(w, right_dists, encoding; positive_gaps=tight_deltas)
        end
    end
end

"""
    _read_ir_lr(r, encoding, mil, T, vertex_id)

Read neighbor list written by `_write_ir_lr`.
"""
function _read_ir_lr(r::BitReader, encoding::Symbol, mil::Int,
                      ::Type{T}, vertex_id; tight_deltas::Bool=false) where {T<:Unsigned}
    if vertex_id === nothing
        return read_intervals_and_residuals(r, encoding, mil, T; vertex_id=vertex_id)
    end
    vid = T(vertex_id)

    # Read intervals (same format as standard)
    num_intervals = Int(read_encoded_value(r, encoding, T)) - 1
    neighbors = T[]

    if num_intervals > 0
        prev_ref = T(0)
        for idx in 1:num_intervals
            if idx == 1
                raw_start = read_encoded_value(r, encoding, UInt64)
                start = T(Int64(vid) + Compression._zigzag_decode(raw_start - 1))
            else
                start_delta = read_encoded_value(r, encoding, T)
                start = prev_ref + start_delta
            end
            len = Int(read_encoded_value(r, encoding, T)) - 1 + mil
            for j in 0:(len-1)
                push!(neighbors, start + T(j))
            end
            prev_ref = start
        end
    end

    # Read total residual count
    num_residuals = Int(read_encoded_value(r, encoding, T)) - 1
    if num_residuals > 0
        # Read left count
        n_left = Int(read_encoded_value(r, encoding, T)) - 1
        n_right = num_residuals - n_left

        # Read left distances → reconstruct left values
        if n_left > 0
            left_dists = read_delta(r, encoding, T; max_elements=n_left, positive_gaps=tight_deltas)
            # Reverse distances back to ascending original values
            for i in n_left:-1:1
                push!(neighbors, vid - left_dists[i])
            end
        end

        # Read right distances → reconstruct right values
        if n_right > 0
            right_dists = read_delta(r, encoding, T; max_elements=n_right, positive_gaps=tight_deltas)
            for d in right_dists
                push!(neighbors, vid + d - T(1))
            end
        end
    end

    sort!(neighbors)
    return neighbors
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
function _write_copy_blocks(w::BitWriter, positions::AbstractVector{<:Integer}, encoding::Symbol)
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

# Minimum number of reference candidates to use threaded search
const MIN_CANDIDATES_FOR_THREADING = 8

# Maximum candidates to fully evaluate after overlap screening (2-phase pruning)
# Full model: evaluate all candidates (low cap was causing BPE regression)
const MAX_REF_CANDIDATES_PHASE2 = 1024
# Fast model: limit evaluation for speed
const MAX_REF_CANDIDATES_PHASE2_FAST = 16

# ---------- Analytical cost estimators (no IOBuffer allocation) ----------

"""
    _estimate_small_count_cost(v, varint) → Int

Analytical cost of `write_small_count`: 2 bits for 0-2, 2+varint for 3+.
"""
@inline function _estimate_small_count_cost(v::Integer, varint::Symbol)::Int
    v <= 2 && return 2
    return 2 + estimate_encoded_value_cost(UInt64(v), varint)
end

"""
    _estimate_delta_list_cost(sorted_vals, encoding; vertex_id=nothing, positive_gaps=false) → Int

Analytical cost of `write_delta` for a sorted list.
"""
function _estimate_delta_list_cost(sorted_vals::Vector{T}, encoding::Symbol;
                                    vertex_id=nothing, positive_gaps::Bool=false)::Int where {T<:Unsigned}
    isempty(sorted_vals) && return 0
    cost = 0
    # First value
    if vertex_id !== nothing
        offset = Int64(sorted_vals[1]) - Int64(vertex_id)
        cost += estimate_encoded_value_cost(UInt64(Compression._zigzag_encode(offset) + 1), encoding)
    else
        cost += estimate_encoded_value_cost(sorted_vals[1], encoding)
    end
    # Remaining gaps
    for i in 2:length(sorted_vals)
        gap = sorted_vals[i] - sorted_vals[i-1]
        if positive_gaps
            cost += estimate_encoded_value_cost(gap, encoding)
        else
            cost += estimate_encoded_value_cost(gap + T(1), encoding)
        end
    end
    return cost
end

"""
    _estimate_stop_delta_zigzag_cost(sorted_vals, encoding, vertex_id) → Int

Analytical cost of `_write_stop_delta_zigzag`.
"""
function _estimate_stop_delta_zigzag_cost(sorted_vals::Vector{T}, encoding::Symbol,
                                           vertex_id)::Int where {T<:Unsigned}
    isempty(sorted_vals) && return 1  # just STOP bit
    cost = 1  # STOP bit
    # First value
    cost += 1  # continuation bit
    if vertex_id !== nothing
        offset = Int64(sorted_vals[1]) - Int64(vertex_id)
        cost += estimate_encoded_value_cost(UInt64(Compression._zigzag_encode(offset) + 1), encoding)
    else
        cost += estimate_encoded_value_cost(sorted_vals[1], encoding)
    end
    # Remaining gaps (>= 1 for sorted unique)
    for i in 2:length(sorted_vals)
        cost += 1 + estimate_encoded_value_cost(sorted_vals[i] - sorted_vals[i-1], encoding)
    end
    return cost
end

"""
    _estimate_copy_blocks_cost(positions, encoding) → Int

Analytical cost of `_write_copy_blocks`: small_count(nblocks) + per-block costs.
"""
function _estimate_copy_blocks_cost(positions::AbstractVector{<:Integer}, encoding::Symbol)::Int
    isempty(positions) && return _estimate_small_count_cost(0, encoding)
    # Compute copy blocks
    nblocks = 0
    i = 1
    block_cost = 0
    prev_end = 0
    while i <= length(positions)
        start = positions[i]
        len = 1
        while i + len <= length(positions) && positions[i + len] == start + len
            len += 1
        end
        nblocks += 1
        if nblocks == 1
            block_cost += estimate_encoded_value_cost(UInt32(start), encoding)
            block_cost += estimate_encoded_value_cost(UInt32(len), encoding)
        else
            gap = start - prev_end
            block_cost += estimate_encoded_value_cost(UInt32(gap), encoding)
            block_cost += estimate_encoded_value_cost(UInt32(len), encoding)
        end
        prev_end = start + len
        i += len
    end
    return _estimate_small_count_cost(nblocks, encoding) + block_cost
end

"""
    _estimate_complement_blocks_cost(positions, ref_len, encoding) → Int

Analytical cost of copy-blocks for the complement (skipped) positions.
Derives complement blocks from positions without allocating a complement vector.
The complement of positions in [1..ref_len] is the set of indices NOT in positions.
"""
function _estimate_complement_blocks_cost(positions::AbstractVector{<:Integer}, ref_len::Int, encoding::Symbol)::Int
    n_skipped = ref_len - length(positions)
    n_skipped <= 0 && return _estimate_small_count_cost(0, encoding)

    # Walk through positions to find complement blocks (gaps)
    nblocks = 0
    block_cost = 0
    prev_complement_end = 0

    # Build position blocks first, then derive complement blocks from gaps
    i = 1
    pos_block_start = 0
    pos_block_end = 0
    complement_cursor = 1  # next expected skipped index

    while i <= length(positions)
        pos_block_start = positions[i]
        blen = 1
        while i + blen <= length(positions) && positions[i + blen] == pos_block_start + blen
            blen += 1
        end
        pos_block_end = pos_block_start + blen - 1

        # Gap before this position block: [complement_cursor, pos_block_start-1]
        if complement_cursor < pos_block_start
            gap_len = pos_block_start - complement_cursor
            nblocks += 1
            if nblocks == 1
                block_cost += estimate_encoded_value_cost(UInt32(complement_cursor), encoding)
            else
                gap_from_prev = complement_cursor - prev_complement_end
                block_cost += estimate_encoded_value_cost(UInt32(gap_from_prev), encoding)
            end
            block_cost += estimate_encoded_value_cost(UInt32(gap_len), encoding)
            prev_complement_end = complement_cursor + gap_len
        end

        complement_cursor = pos_block_end + 1
        i += blen
    end

    # Trailing gap after last position block: [complement_cursor, ref_len]
    if complement_cursor <= ref_len
        gap_len = ref_len - complement_cursor + 1
        nblocks += 1
        if nblocks == 1
            block_cost += estimate_encoded_value_cost(UInt32(complement_cursor), encoding)
        else
            gap_from_prev = complement_cursor - prev_complement_end
            block_cost += estimate_encoded_value_cost(UInt32(gap_from_prev), encoding)
        end
        block_cost += estimate_encoded_value_cost(UInt32(gap_len), encoding)
    end

    return _estimate_small_count_cost(nblocks, encoding) + block_cost
end

"""
    _estimate_ir_lr_cost(neighbors, encoding, mil, vertex_id; tight_deltas) → Int

CG-specific analytical cost of `_write_ir_lr` (intervals + LR-split residuals).
"""
function _estimate_ir_lr_cost(neighbors::Vector{T}, encoding::Symbol, mil::Int,
                               vertex_id; tight_deltas::Bool=false)::Int where {T<:Unsigned}
    if vertex_id === nothing || isempty(neighbors)
        return estimate_interval_runlength_encoding_cost(neighbors, encoding, mil, 3; vertex_id=vertex_id)
    end
    vid = T(vertex_id)
    intervals, residuals = compress_intervals(neighbors, mil)
    cost = 0

    # Interval count
    cost += estimate_encoded_value_cost(T(length(intervals)) + T(1), encoding)

    # Intervals
    prev_ref = T(0)
    for (idx, (start, len)) in enumerate(intervals)
        if idx == 1
            encoded_start = UInt64(Compression._zigzag_encode(Int64(start) - Int64(vid)) + 1)
            cost += estimate_encoded_value_cost(encoded_start, encoding)
        else
            cost += estimate_encoded_value_cost(start - prev_ref, encoding)
        end
        cost += estimate_encoded_value_cost(T(len - mil + 1), encoding)
        prev_ref = start  # CG doesn't use tight_intervals for inter-interval gaps
    end

    # Residual count
    cost += estimate_encoded_value_cost(T(length(residuals)) + T(1), encoding)

    if !isempty(residuals)
        # Split at vid
        split_idx = 0
        @inbounds for i in 1:length(residuals)
            residuals[i] < vid ? (split_idx = i) : break
        end
        n_left = split_idx
        n_right = length(residuals) - n_left

        # Left count
        cost += estimate_encoded_value_cost(T(n_left) + T(1), encoding)

        # Left distances (ascending from vid)
        if n_left > 0
            prev = T(0)
            for i in 1:n_left
                d = vid - residuals[n_left - i + 1]
                if i == 1
                    cost += estimate_encoded_value_cost(d, encoding)
                elseif tight_deltas
                    cost += estimate_encoded_value_cost(d - prev, encoding)
                else
                    cost += estimate_encoded_value_cost(d - prev + T(1), encoding)
                end
                prev = d
            end
        end

        # Right distances (ascending from vid, +1 shift)
        if n_right > 0
            prev = T(0)
            for i in 1:n_right
                d = residuals[n_left + i] - vid + T(1)
                if i == 1
                    cost += estimate_encoded_value_cost(d, encoding)
                elseif tight_deltas
                    cost += estimate_encoded_value_cost(d - prev, encoding)
                else
                    cost += estimate_encoded_value_cost(d - prev + T(1), encoding)
                end
                prev = d
            end
        end
    end

    return cost
end

"""
    _estimate_raw_cost_analytical(nl, params, T, _zz_vid) → (Int, Bool)

Analytical bit cost of encoding `nl` without a reference (raw mode).
No IOBuffer allocations — pure arithmetic cost estimation.
Returns (raw_bits, raw_use_intervals).
"""
function _estimate_raw_cost_analytical(nl::Vector{T}, params::CGParams,
                                        ::Type{T}, _zz_vid)::Tuple{Int,Bool} where {T<:Unsigned}
    raw_use_intervals = false

    # Fast cost model: skip interval/LR estimation, use stop-delta + adaptive interval check
    if params.cost_model == Compression.COST_MODEL_FAST
        if params.intra_stop_deltas
            raw_bits = _estimate_stop_delta_zigzag_cost(nl, :fibonacci, _zz_vid)
        else
            raw_bits = _estimate_small_count_cost(length(nl), params.count_varint)
            if !isempty(nl)
                raw_bits += _estimate_delta_list_cost(nl, :fibonacci; vertex_id=_zz_vid)
            end
        end
        # Keep adaptive decision (cheap interval estimation vs stop-delta)
        if params.intra_raw_adaptive && params.intra_stop_deltas
            raw_iv_bits = estimate_interval_runlength_encoding_cost(nl, :fibonacci, params.intra_adapt_mil, 3; vertex_id=_zz_vid)
            if 1 + raw_iv_bits < 1 + raw_bits
                raw_use_intervals = true
            end
            raw_bits = min(1 + raw_bits, 1 + raw_iv_bits)
        end
        return raw_bits, raw_use_intervals
    end

    if params.intra_intervals && params.intra_lr_split
        raw_bits = _estimate_ir_lr_cost(nl, :fibonacci, params.intra_mil, _zz_vid; tight_deltas=params.intra_tight_deltas)
    elseif params.intra_intervals
        raw_bits = estimate_interval_runlength_encoding_cost(nl, :fibonacci, params.intra_mil, 3; vertex_id=_zz_vid)
    elseif params.intra_stop_deltas
        raw_bits = _estimate_stop_delta_zigzag_cost(nl, :fibonacci, _zz_vid)
    else
        raw_bits = _estimate_small_count_cost(length(nl), params.count_varint)
        if !isempty(nl)
            raw_bits += _estimate_delta_list_cost(nl, :fibonacci; vertex_id=_zz_vid)
        end
    end

    if params.intra_raw_adaptive && params.intra_stop_deltas
        raw_iv_bits = estimate_interval_runlength_encoding_cost(nl, :fibonacci, params.intra_adapt_mil, 3; vertex_id=_zz_vid)
        if 1 + raw_iv_bits < 1 + raw_bits
            raw_use_intervals = true
        end
        raw_bits = min(1 + raw_bits, 1 + raw_iv_bits)
    end

    return raw_bits, raw_use_intervals
end

"""
    _evaluate_candidate_analytical(positions, adds, ref_len, params, T, _zz_vid; best_so_far) → (Int, UInt8, Bool)

Analytical bit cost of encoding with a reference, given pre-computed positions/adds.
Returns (total_bits, copy_mode, add_use_intervals).
Supports early termination: if positions cost alone exceeds `best_so_far`, returns early.
"""
function _evaluate_candidate_analytical(positions::Vector{T}, adds::Vector{T},
                                         ref_len::Int, params::CGParams, ::Type{T},
                                         _zz_vid; best_so_far::Int=typemax(Int))::Tuple{Int,UInt8,Bool} where {T<:Unsigned}
    # --- Positions estimation ---
    local pos_bits::Int
    local copy_mode::UInt8 = 0x00

    if params.cost_model == Compression.COST_MODEL_FAST
        # Fast model: 3-way adaptive without allocating complement vector
        if params.intra_copy_adaptive && params.intra_copy_blocks
            cb_bits = _estimate_copy_blocks_cost(positions, params.varint)
            bm_bits = ref_len
            cc_bits = _estimate_complement_blocks_cost(positions, ref_len, params.varint)
            t_bm = 1 + bm_bits; t_cb = 2 + cb_bits; t_cc = 2 + cc_bits
            if t_bm <= t_cb && t_bm <= t_cc
                pos_bits = t_bm; copy_mode = 0x00
            elseif t_cb <= t_cc
                pos_bits = t_cb; copy_mode = 0x01
            else
                pos_bits = t_cc; copy_mode = 0x02
            end
        elseif params.intra_copy_blocks
            pos_bits = _estimate_copy_blocks_cost(positions, params.varint)
            copy_mode = 0x01
        else
            pos_bits = _estimate_small_count_cost(length(positions), params.count_varint)
            if !isempty(positions)
                pos_bits += _estimate_delta_list_cost(UInt32.(positions), :fibonacci)
            end
        end

        # Early termination
        if pos_bits >= best_so_far
            return pos_bits, copy_mode, false
        end

        # Fast additions: stop-delta + adaptive interval check (cheap)
        local fast_add_use_intervals::Bool = false
        fast_add_bits = if params.intra_stop_deltas
            _estimate_stop_delta_zigzag_cost(adds, :fibonacci, _zz_vid)
        else
            ab = _estimate_small_count_cost(length(adds), params.count_varint)
            if !isempty(adds)
                ab += _estimate_delta_list_cost(adds, :fibonacci; vertex_id=_zz_vid)
            end
            ab
        end
        if params.intra_add_adaptive && params.intra_stop_deltas
            iv_bits = estimate_interval_runlength_encoding_cost(adds, :fibonacci, params.intra_adapt_mil, 3; vertex_id=_zz_vid)
            if 1 + iv_bits < 1 + fast_add_bits
                fast_add_use_intervals = true
            end
            fast_add_bits = min(1 + fast_add_bits, 1 + iv_bits)
        end

        return pos_bits + fast_add_bits, copy_mode, fast_add_use_intervals
    end

    if params.intra_copy_adaptive && params.intra_copy_blocks
        # Copy-blocks cost
        cb_bits = _estimate_copy_blocks_cost(positions, params.varint)
        # Bitmap cost
        bm_bits = ref_len
        # Complement cost
        n_skipped = ref_len - length(positions)
        if n_skipped > 0
            # Build complement inline for cost estimation
            skipped = Vector{T}(undef, n_skipped)
            si = 0; pi = 1
            @inbounds for p in T(1):T(ref_len)
                if pi <= length(positions) && positions[pi] == p
                    pi += 1
                else
                    si += 1
                    skipped[si] = p
                end
            end
            cc_bits = _estimate_copy_blocks_cost(skipped, params.varint)
        else
            cc_bits = _estimate_small_count_cost(0, params.varint)
        end
        # 3-way min with mode flag costs
        t_bm = 1 + bm_bits; t_cb = 2 + cb_bits; t_cc = 2 + cc_bits
        if t_bm <= t_cb && t_bm <= t_cc
            pos_bits = t_bm; copy_mode = 0x00
        elseif t_cb <= t_cc
            pos_bits = t_cb; copy_mode = 0x01
        else
            pos_bits = t_cc; copy_mode = 0x02
        end
    elseif params.intra_copy_blocks
        pos_bits = _estimate_copy_blocks_cost(positions, params.varint)
        copy_mode = 0x01
    else
        pos_bits = _estimate_small_count_cost(length(positions), params.count_varint)
        if !isempty(positions)
            pos_bits += _estimate_delta_list_cost(UInt32.(positions), :fibonacci)
        end
    end

    # Early termination: if positions alone exceed best, skip additions
    if pos_bits >= best_so_far
        return pos_bits, copy_mode, false
    end

    # --- Additions estimation ---
    local add_bits::Int
    local add_use_intervals::Bool = false
    if params.intra_intervals && params.intra_lr_split
        add_bits = _estimate_ir_lr_cost(adds, :fibonacci, params.intra_mil, _zz_vid; tight_deltas=params.intra_tight_deltas)
    elseif params.intra_intervals
        add_bits = estimate_interval_runlength_encoding_cost(adds, :fibonacci, params.intra_mil, 3; vertex_id=_zz_vid)
    elseif params.intra_add_adaptive && params.intra_stop_deltas
        sd_bits = _estimate_stop_delta_zigzag_cost(adds, :fibonacci, _zz_vid)
        iv_bits = estimate_interval_runlength_encoding_cost(adds, :fibonacci, params.intra_adapt_mil, 3; vertex_id=_zz_vid)
        if 1 + iv_bits < 1 + sd_bits
            add_use_intervals = true
        end
        add_bits = min(1 + sd_bits, 1 + iv_bits)
    elseif params.intra_stop_deltas
        add_bits = _estimate_stop_delta_zigzag_cost(adds, :fibonacci, _zz_vid)
    else
        add_bits = _estimate_small_count_cost(length(adds), params.count_varint)
        if !isempty(adds)
            add_bits += _estimate_delta_list_cost(adds, :fibonacci; vertex_id=_zz_vid)
        end
    end

    return pos_bits + add_bits, copy_mode, add_use_intervals
end

"""
    _evaluate_candidate_greedy_analytical(positions, adds, params, T, _zz_vid, mil_options) → (Int, Int)

Analytical greedy MIL variant: try each mil value for additions cost.
Returns (best_total_bits, best_mil).
"""
function _evaluate_candidate_greedy_analytical(positions::Vector{T}, adds::Vector{T},
                                                params::CGParams, ::Type{T},
                                                _zz_vid,
                                                mil_options::Vector{Int})::Tuple{Int,Int} where {T<:Unsigned}
    # Fast model: keep copy-blocks for positions, single MIL with interval vs delta adaptive
    if params.cost_model == Compression.COST_MODEL_FAST
        if params.intra_copy_blocks
            pos_bits = _estimate_copy_blocks_cost(positions, params.varint)
        else
            pos_bits = _estimate_small_count_cost(length(positions), params.count_varint)
            if !isempty(positions)
                pos_bits += _estimate_delta_list_cost(UInt32.(positions), :fibonacci)
            end
        end
        # Single fixed MIL: compare interval vs stop-delta for additions
        fast_mil = params.intra_adapt_mil > 0 ? params.intra_adapt_mil : mil_options[1]
        iv_add = estimate_interval_runlength_encoding_cost(adds, :fibonacci, fast_mil, 3; vertex_id=_zz_vid)
        sd_add = if params.intra_stop_deltas
            _estimate_stop_delta_zigzag_cost(adds, :fibonacci, _zz_vid)
        else
            ab = _estimate_small_count_cost(length(adds), params.count_varint)
            if !isempty(adds)
                ab += _estimate_delta_list_cost(adds, :fibonacci; vertex_id=_zz_vid)
            end
            ab
        end
        add_bits = min(iv_add, sd_add)
        return pos_bits + add_bits, fast_mil
    end

    # Positions cost (doesn't depend on mil) — analytical
    if params.intra_copy_blocks
        pos_bits = _estimate_copy_blocks_cost(positions, params.varint)
    else
        pos_bits = _estimate_small_count_cost(length(positions), params.count_varint)
        if !isempty(positions)
            pos_bits += _estimate_delta_list_cost(UInt32.(positions), :fibonacci)
        end
    end

    # Try each mil for additions — analytical
    best_bits = typemax(Int)
    best_mil = mil_options[1]
    for mil in mil_options
        if params.intra_lr_split
            add_bits = _estimate_ir_lr_cost(adds, :fibonacci, mil, _zz_vid; tight_deltas=params.intra_tight_deltas)
        else
            add_bits = estimate_interval_runlength_encoding_cost(adds, :fibonacci, mil, 3; vertex_id=_zz_vid)
        end
        total = pos_bits + add_bits
        if total < best_bits
            best_bits = total
            best_mil = mil
        end
    end

    return best_bits, best_mil
end

"""
    _merge_positions_adds!(nl, ref, positions, adds)

Two-pointer merge producing positions (1-based indices into ref) and adds (unmatched from nl).
Reuses the provided vectors.
"""
@inline function _merge_positions_adds!(nl::Vector{T}, ref::Vector{T},
                                         positions::Vector{T}, adds::Vector{T}) where {T<:Unsigned}
    empty!(positions); empty!(adds)
    i = 1; j = 1
    @inbounds while i <= length(nl) && j <= length(ref)
        if nl[i] == ref[j]
            push!(positions, T(j)); i += 1; j += 1
        elseif nl[i] < ref[j]
            push!(adds, nl[i]); i += 1
        else
            j += 1
        end
    end
    @inbounds while i <= length(nl); push!(adds, nl[i]); i += 1; end
end

"""
    _sorted_overlap_count(a, b) → Int

Count common elements between two sorted vectors using two-pointer scan.
"""
@inline function _sorted_overlap_count(a::Vector{T}, b::Vector{T})::Int where {T}
    count = 0
    i = 1; j = 1
    @inbounds while i <= length(a) && j <= length(b)
        if a[i] == b[j]
            count += 1; i += 1; j += 1
        elseif a[i] < b[j]
            i += 1
        else
            j += 1
        end
    end
    return count
end

# --------------------------
# Core encoder
# --------------------------

"""
    encode_level(w, g, P; params=CGParams())

Encode one coarsening level for graph `g` with partition `P` using CG.

Inputs:
- w::BitWriter: output bitstream
- g::LightGraphs.AbstractGraph{T<:Unsigned}
- P::Vector{Vector{T}}: list of clusters with global vertex ids
- params::CGParams: encoding parameters
"""
function encode_level(w::BitWriter, g::AbstractGraph{T}, P::Vector{Vector{T}}; params::CGParams=CGParams(), stats::Union{Nothing,CGStats}=nothing, progress::Union{Nothing,Function}=nothing, cluster_offsets::Union{Nothing,Vector{Int}}=nothing) where {T<:Unsigned}
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
    n_clusters = length(clusters)
    for (ci, C) in enumerate(clusters)
        # Record cluster start bit offset for index mode
        if cluster_offsets !== nothing
            cluster_offsets[ci] = _total_bits(w)
        end
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
                local w_blk = BitWriter()
                Compression.write_compressed_graph_data(w_blk, neighbor_lists, :children, :fibonacci, true, true, true, params.intra_ref_window)
                flush_bitwriter(w_blk; flush_last_bits=true)
                block_bytes = copy(get_bytes(w_blk))
                block_bits = length(block_bytes) * 8
                # Cheap baseline: raw count + delta per vertex
                local w_base = BitWriter()
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
                baseline_bits = bytes_written(w_base) * 8
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
            # Pre-build ALL neighbor lists for the cluster upfront (optimization #4)
            all_nl = Vector{Vector{T}}(undef, s)
            for (i, u) in enumerate(C)
                nl = T[]
                for v in outneighbors(g, Int(u))
                    lv = get(local_index, v, 0)
                    if lv != 0
                        push!(nl, T(lv))
                    end
                end
                sort!(nl)
                all_nl[i] = nl
            end

            # First pass: decide refs and compute ref data
            use_ref_vec = Vector{Bool}(undef, s)
            ref_delta_vec = Vector{UInt32}(undef, s)
            ref_positions_list = Vector{Vector{T}}(undef, s)
            additions_list = Vector{Vector{T}}(undef, s)
            ref_len_list = Vector{Int}(undef, s)
            # Decisions for encoding modes (stored during search, replayed in write phase)
            copy_mode_vec = Vector{UInt8}(undef, s)    # 0=bitmap, 1=copy-blocks, 2=complement
            add_mode_vec = Vector{Bool}(undef, s)      # true=intervals, false=stop-delta (ref adds)
            raw_mode_vec = Vector{Bool}(undef, s)      # true=intervals, false=stop-delta (raw)
            # Per-vertex mil (greedy search populates this; otherwise fixed)
            mil_vec = fill(params.intra_mil, s)
            _mil_options = [2, 3, 4, 5]

            # Double-buffer swap: two position/adds vector pairs (optimization #5)
            _pos_a = T[]; _adds_a = T[]
            _pos_b = T[]; _adds_b = T[]

            for idx_local in 1:s
                progress !== nothing && progress(idx_local, s, ci, n_clusters)
                nl = all_nl[idx_local]
                # decide reference and mil
                use_ref = false; ref_delta_val = UInt32(0)
                best_copy_mode = 0x00; best_add_mode = false; best_raw_mode = false
                if params.intra_greedy_mil
                    # Greedy per-vertex mil search: try all mil values for raw and ref (analytical)
                    best_bits = typemax(Int)
                    best_mil_val = params.intra_mil
                    best_is_ref = false
                    best_ref_idx = 0

                    _zz_vid = params.intra_zigzag ? T(idx_local) : nothing
                    if params.cost_model == Compression.COST_MODEL_FAST
                        # Fast model: single MIL, compare interval vs stop-delta
                        fast_mil = params.intra_adapt_mil > 0 ? params.intra_adapt_mil : params.intra_mil
                        iv_raw = estimate_interval_runlength_encoding_cost(nl, :fibonacci, fast_mil, 3; vertex_id=_zz_vid)
                        sd_raw = if params.intra_stop_deltas
                            _estimate_stop_delta_zigzag_cost(nl, :fibonacci, _zz_vid)
                        else
                            ab = _estimate_small_count_cost(length(nl), params.count_varint)
                            if !isempty(nl)
                                ab += _estimate_delta_list_cost(nl, :fibonacci; vertex_id=_zz_vid)
                            end
                            ab
                        end
                        raw_bits = min(iv_raw, sd_raw)
                        if raw_bits < best_bits
                            best_bits = raw_bits
                            best_mil_val = fast_mil
                            best_is_ref = false
                        end
                    else
                    for mil in _mil_options
                        if params.intra_lr_split
                            raw_bits = _estimate_ir_lr_cost(nl, :fibonacci, mil, _zz_vid; tight_deltas=params.intra_tight_deltas)
                        else
                            raw_bits = estimate_interval_runlength_encoding_cost(nl, :fibonacci, mil, 3; vertex_id=_zz_vid)
                        end
                        if raw_bits < best_bits
                            best_bits = raw_bits
                            best_mil_val = mil
                            best_is_ref = false
                        end
                    end
                    end

                    # Try ref encoding with 2-phase pruning + analytical cost
                    if params.intra_ref_enabled && idx_local > 1
                        wstart = max(1, idx_local - params.intra_ref_window)
                        wend = idx_local - 1
                        n_candidates = wend - wstart + 1

                        # Phase 1: overlap screening (cheap O(|nl|+|ref|) per candidate)
                        _max_k_greedy = params.cost_model == Compression.COST_MODEL_FAST ? MAX_REF_CANDIDATES_PHASE2_FAST : MAX_REF_CANDIDATES_PHASE2
                        if n_candidates > _max_k_greedy
                            overlap_scores = Vector{Tuple{Int,Int}}(undef, n_candidates)
                            for (ci2, rix) in enumerate(wstart:wend)
                                ov = _sorted_overlap_count(nl, all_nl[rix])
                                overlap_scores[ci2] = (ov, rix)
                            end
                            sort!(overlap_scores; by = x -> -x[1])
                            phase2_indices = [overlap_scores[k][2] for k in 1:min(_max_k_greedy, n_candidates)]
                        else
                            phase2_indices = collect(wstart:wend)
                        end

                        # Phase 2: full analytical evaluation on top candidates
                        for rix in phase2_indices
                            _merge_positions_adds!(nl, all_nl[rix], _pos_a, _adds_a)
                            bits, mil_val = _evaluate_candidate_greedy_analytical(_pos_a, _adds_a, params, T, _zz_vid, _mil_options)
                            if bits < best_bits
                                best_bits = bits
                                best_mil_val = mil_val
                                best_is_ref = true
                                best_ref_idx = rix
                            end
                        end
                    end

                    if best_is_ref
                        use_ref = true
                        ref_delta_val = UInt32(idx_local - best_ref_idx)
                        # Final merge for the winner — swap into storage
                        _merge_positions_adds!(nl, all_nl[best_ref_idx], _pos_a, _adds_a)
                        ref_positions_list[idx_local] = copy(_pos_a)
                        additions_list[idx_local] = copy(_adds_a)
                    else
                        ref_positions_list[idx_local] = T[]
                        additions_list[idx_local] = T[]
                    end
                    mil_vec[idx_local] = best_mil_val
                elseif params.intra_ref_enabled && idx_local > 1
                    # Analytical reference decision with 2-phase pruning + early termination
                    _zz_vid = params.intra_zigzag ? T(idx_local) : nothing

                    # Raw estimation (analytical)
                    raw_bits, raw_use_iv = _estimate_raw_cost_analytical(nl, params, T, _zz_vid)
                    best_raw_mode = raw_use_iv

                    # Ref delta header overhead
                    ref_overhead = 0
                    if params.intra_ref_fixwidth
                        ref_overhead = max(1, ceil(Int, log2(params.intra_ref_window)))
                    end

                    wstart = max(1, idx_local - params.intra_ref_window)
                    wend = idx_local - 1
                    n_candidates = wend - wstart + 1

                    # Phase 1: overlap screening
                    _max_k_ref = params.cost_model == Compression.COST_MODEL_FAST ? MAX_REF_CANDIDATES_PHASE2_FAST : MAX_REF_CANDIDATES_PHASE2
                    if n_candidates > _max_k_ref
                        overlap_scores = Vector{Tuple{Int,Int}}(undef, n_candidates)
                        for (ci2, rix) in enumerate(wstart:wend)
                            ov = _sorted_overlap_count(nl, all_nl[rix])
                            overlap_scores[ci2] = (ov, rix)
                        end
                        sort!(overlap_scores; by = x -> -x[1])
                        phase2_indices = [overlap_scores[k][2] for k in 1:min(_max_k_ref, n_candidates)]
                    else
                        phase2_indices = collect(wstart:wend)
                    end

                    # Phase 2: analytical evaluation with early termination
                    best_bits = raw_bits
                    best_idx = 0
                    # Use double-buffer swap: _pos_a/_adds_a for current, _pos_b/_adds_b for best
                    for rix in phase2_indices
                        _merge_positions_adds!(nl, all_nl[rix], _pos_a, _adds_a)
                        ref_len = length(all_nl[rix])
                        bits, cm, aim = _evaluate_candidate_analytical(_pos_a, _adds_a, ref_len, params, T, _zz_vid; best_so_far=best_bits - ref_overhead)
                        total = bits + ref_overhead
                        if total < best_bits
                            best_bits = total
                            best_idx = rix
                            best_copy_mode = cm
                            best_add_mode = aim
                            # Swap buffers: _pos_b/_adds_b now hold the best result
                            _pos_a, _pos_b = _pos_b, _pos_a
                            _adds_a, _adds_b = _adds_b, _adds_a
                        end
                    end
                    if best_idx > 0
                        use_ref = true
                        ref_delta_val = UInt32(idx_local - best_idx)
                        # Best result is in _pos_b/_adds_b (after last swap)
                        ref_positions_list[idx_local] = copy(_pos_b)
                        additions_list[idx_local] = copy(_adds_b)
                    end
                    if !use_ref
                        ref_positions_list[idx_local] = T[]
                        additions_list[idx_local] = T[]
                    end
                else
                    ref_positions_list[idx_local] = T[]
                    additions_list[idx_local] = T[]
                end
                use_ref_vec[idx_local] = use_ref
                ref_delta_vec[idx_local] = ref_delta_val
                copy_mode_vec[idx_local] = best_copy_mode
                add_mode_vec[idx_local] = best_add_mode
                raw_mode_vec[idx_local] = best_raw_mode
                if use_ref
                    ref_index = idx_local - Int(ref_delta_val)
                    ref_len_list[idx_local] = ref_index >= 1 ? length(all_nl[ref_index]) : 0
                else
                    ref_len_list[idx_local] = 0
                end
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
            # In index mode (cluster_offsets !== nothing), use two-pass encoding:
            # encode payloads to temp buffer, then write intra-vertex offset table + data
            local _idx_mode = cluster_offsets !== nothing
            local _vtx_payload_bw = _idx_mode ? BitWriter() : nothing
            local _vtx_offsets = _idx_mode ? Vector{Int}(undef, s + 1) : nothing
            local pw = _idx_mode ? _vtx_payload_bw : w  # payload writer

            # Pre-allocate reusable complement vector for write phase
            _wp_skipped = T[]

            for idx_local in 1:s
                if _idx_mode
                    _vtx_offsets[idx_local] = _total_bits(pw)
                end
                _zz_vid = params.intra_zigzag ? T(idx_local) : nothing

                if use_ref_vec[idx_local]
                    # write copied positions into reference list
                    ref_positions = ref_positions_list[idx_local]
                    additions = additions_list[idx_local]
                    ref_len = ref_len_list[idx_local]
                    local bpos0 = _total_bits(pw)
                    if params.intra_copy_adaptive && params.intra_copy_blocks
                        # 3-way adaptive: use stored decision from search phase
                        local _cm = copy_mode_vec[idx_local]
                        if _cm == 0x00
                            # Bitmap mode
                            write_bit(pw, true)   # outer=1: bitmap
                            local _bm_arr = fill(false, ref_len)
                            for p in ref_positions; if 1 <= Int(p) <= ref_len; _bm_arr[Int(p)] = true; end; end
                            for b in _bm_arr; write_bit(pw, b); end
                            if stats !== nothing; stats.intra_copy_bitmap_count += 1; end
                        elseif _cm == 0x01
                            # Copy-blocks mode
                            write_bit(pw, false); write_bit(pw, false)
                            _write_copy_blocks(pw, ref_positions, params.varint)
                            if stats !== nothing; stats.intra_copy_blocks_count += 1; end
                        else
                            # Complement mode — build complement using pre-allocated vector
                            write_bit(pw, false); write_bit(pw, true)
                            empty!(_wp_skipped)
                            local _pi = 1
                            @inbounds for p in T(1):T(ref_len)
                                if _pi <= length(ref_positions) && ref_positions[_pi] == p
                                    _pi += 1
                                else
                                    push!(_wp_skipped, p)
                                end
                            end
                            _write_copy_blocks(pw, _wp_skipped, params.varint)
                            if stats !== nothing; stats.intra_copy_complement_count += 1; end
                        end

                        # Profiling: overlap fraction and addition count
                        if stats !== nothing
                            local _ov = ref_len > 0 ? length(ref_positions) / ref_len : 0.0
                            local _bucket = min(9, floor(Int, _ov * 10))
                            stats.intra_overlap_hist[_bucket + 1] += 1
                        end
                    elseif params.intra_copy_blocks
                        _write_copy_blocks(pw, ref_positions, params.varint)
                    elseif ref_len > 0
                        copied = fill(false, ref_len)
                        for p in ref_positions
                            if 1 <= Int(p) <= ref_len; copied[Int(p)] = true; end
                        end
                        Compression.write_bitmap_rle_ones_deltas(pw, copied, params.varint)
                    else
                        # empty token list via small count
                        Compression.write_small_count(pw, T(0), params.count_varint)
                    end
                    local bpos1 = _total_bits(pw)
                    # additions: MGS intervals, custom intervals, or plain delta
                    local ah0 = _total_bits(pw)
                    if (params.intra_intervals || params.intra_greedy_mil) && params.intra_lr_split
                        _write_ir_lr(pw, additions, :fibonacci, mil_vec[idx_local], _zz_vid; tight_deltas=params.intra_tight_deltas)
                    elseif params.intra_intervals || params.intra_greedy_mil
                        Compression.write_intervals_and_residuals(pw, additions, :fibonacci, mil_vec[idx_local]; vertex_id=_zz_vid)
                    elseif params.additions_mode == :intervals
                        # detect runs and write intervals + singles
                        runs = Vector{Tuple{T,Int}}(); singles = T[]
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
                        Compression.write_small_count(pw, T(length(runs)), params.count_varint)
                        for (st, ln) in runs
                            write_encoded_value(pw, T(st), params.varint)
                            write_encoded_value(pw, T(ln), params.varint)
                        end
                        Compression.write_small_count(pw, T(length(singles)), params.count_varint)
                        if !isempty(singles)
                            write_delta(pw, singles, :fibonacci; vertex_id=_zz_vid)
                        end
                    elseif params.intra_add_adaptive && params.intra_stop_deltas
                        # Adaptive: use stored decision from search phase
                        if add_mode_vec[idx_local]
                            write_bit(pw, true)   # intervals
                            Compression.write_intervals_and_residuals(pw, additions, :fibonacci, params.intra_adapt_mil; vertex_id=_zz_vid)
                        else
                            write_bit(pw, false)  # stop-delta
                            _write_stop_delta_zigzag(pw, additions, :fibonacci, _zz_vid)
                        end
                    else
                        if params.intra_stop_deltas
                            _write_stop_delta_zigzag(pw, additions, :fibonacci, _zz_vid)
                        else
                            Compression.write_small_count(pw, T(length(additions)), params.count_varint)
                            if !isempty(additions)
                                write_delta(pw, additions, :fibonacci; vertex_id=_zz_vid)
                            end
                        end
                    end
                    local b3 = _total_bits(pw)
                    if stats !== nothing
                        # positions bitmap is copy payload; additions include intervals/singles payload
                        stats.bits_intra_copy += (bpos1 - bpos0)
                        stats.bits_intra_add += (b3 - bpos1)
                        stats.intra_ref_used += 1
                        stats.intra_add_count_total += length(additions)
                    end
                else
                    # raw
                    nl = all_nl[idx_local]
                    local rb0 = _total_bits(pw)
                    if (params.intra_intervals || params.intra_greedy_mil) && params.intra_lr_split
                        _write_ir_lr(pw, nl, :fibonacci, mil_vec[idx_local], _zz_vid; tight_deltas=params.intra_tight_deltas)
                    elseif params.intra_intervals || params.intra_greedy_mil
                        Compression.write_intervals_and_residuals(pw, nl, :fibonacci, mil_vec[idx_local]; vertex_id=_zz_vid)
                    elseif params.intra_raw_adaptive && params.intra_stop_deltas
                        # Adaptive: use stored decision from search phase
                        if raw_mode_vec[idx_local]
                            write_bit(pw, true)   # intervals
                            Compression.write_intervals_and_residuals(pw, nl, :fibonacci, params.intra_adapt_mil; vertex_id=_zz_vid)
                        else
                            write_bit(pw, false)  # stop-delta
                            _write_stop_delta_zigzag(pw, nl, :fibonacci, _zz_vid)
                        end
                    elseif params.intra_stop_deltas
                        _write_stop_delta_zigzag(pw, nl, :fibonacci, _zz_vid)
                    else
                        # raw count: use small-count which escapes to varint for large values and handles zero
                        Compression.write_small_count(pw, T(length(nl)), params.count_varint)
                        if !isempty(nl)
                            write_delta(pw, nl, :fibonacci; vertex_id=_zz_vid)
                        end
                    end
                    if stats !== nothing
                        stats.bits_intra_raw += _total_bits(pw) - rb0
                        stats.intra_no_ref += 1
                    end
                end
            end

            # In index mode, write intra-vertex offset table then buffered payload data
            if _idx_mode
                # Record exact payload bit count BEFORE flush (flush adds padding)
                local _payload_total_bits = _total_bits(_vtx_payload_bw)
                _vtx_offsets[s + 1] = _payload_total_bits
                flush_bitwriter(_vtx_payload_bw; flush_last_bits=true)

                local _max_vtx_offset = _vtx_offsets[s + 1]
                local _vtx_ew = _max_vtx_offset > 0 ? max(Int(ceil(log2(_max_vtx_offset + 1))), 1) : 1

                local _sk = params.index_sample_k
                if _sk > 0 && s > _sk
                    # Two-level sampled offsets within cluster
                    write_bit(w, true)  # sampled flag
                    @assert _sk >= 4 && _sk % 4 == 0
                    write_value(w, UInt64(_sk ÷ 4 - 1), 8)
                    write_value(w, UInt64(_vtx_ew), 6)
                    for _vi in 1:_sk:s
                        write_value(w, UInt64(_vtx_offsets[_vi]), _vtx_ew)
                    end
                    write_value(w, UInt64(_vtx_offsets[s + 1]), _vtx_ew)
                else
                    # Full offset table
                    write_bit(w, false)
                    write_value(w, UInt64(_vtx_ew), 6)
                    for _vi in 1:(s + 1)
                        write_value(w, UInt64(_vtx_offsets[_vi]), _vtx_ew)
                    end
                end

                # Write buffered payload data precisely (avoid byte-padding corruption)
                local _vtx_data = get_bytes(_vtx_payload_bw)
                local _full_bytes = _payload_total_bits ÷ 8
                local _remaining_bits = _payload_total_bits % 8
                if _full_bytes > 0
                    write_bytes(w, _vtx_data[1:_full_bytes])
                end
                if _remaining_bits > 0 && _full_bytes < length(_vtx_data)
                    local _last_byte = _vtx_data[_full_bytes + 1]
                    for _bi in 0:(_remaining_bits - 1)
                        write_bit(w, ((_last_byte >> (7 - _bi)) & 0x01) == 0x01)
                    end
                end
            end
        end
    end
    if stats !== nothing
        stats.bits_intra += _total_bits(w) - bits_before
    end

    # Record inter-cluster start bit offset for index mode
    if cluster_offsets !== nothing
        cluster_offsets[n_clusters + 1] = _total_bits(w)
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

    # Record inter-section start for computing per-source-cluster relative offsets
    local _inter_start = _total_bits(w)

    # Emit per-A targets and AB groups in fixed A order (1..K)
    for ia in 1:mclusters
        # Record per-source-cluster inter offset (relative to inter start)
        if cluster_offsets !== nothing
            cluster_offsets[n_clusters + 1 + ia] = _total_bits(w) - _inter_start
        end

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
    decode_level(r::BitReader, params::CGParams; T::Type{<:Unsigned}=UInt32, directed::Bool=true, coding_scheme::Symbol=:children)

Decode one CG coarsening level from the bitstream.
Returns a `Dict{T, Vector{T}}` mapping global vertex ID → sorted outneighbors.
"""
function decode_level(r::BitReader, params::CGParams; T::Type{<:Unsigned}=UInt32, directed::Bool=true, coding_scheme::Symbol=:children)
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

            # In index mode, read intra-vertex offset table (skip for sequential decode)
            if coding_scheme == :index
                local _is_sampled_vtx = read_bit(r)
                if _is_sampled_vtx
                    local _sk_vtx = (Int(read_value(r, 8, UInt64)) + 1) * 4
                    local _vtx_ew = Int(read_value(r, 6, UInt64))
                    local _n_sampled_vtx = length(1:_sk_vtx:s)
                    for _ in 1:(_n_sampled_vtx + 1)
                        read_value(r, _vtx_ew, UInt64)  # skip sampled offsets
                    end
                else
                    local _vtx_ew = Int(read_value(r, 6, UInt64))
                    for _ in 1:(s + 1)
                        read_value(r, _vtx_ew, UInt64)  # skip full offsets
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
                    if params.intra_copy_adaptive && params.intra_copy_blocks
                        # Nested mode bits: outer=1 → bitmap; outer=0,inner=0 → copy-blocks;
                        #                   outer=0,inner=1 → complement (skipped positions)
                        copied_vals = Int[]
                        if read_bit(r)  # outer=1: bitmap
                            for p in 1:ref_len
                                if read_bit(r); push!(copied_vals, ref_list[p]); end
                            end
                        else            # outer=0: copy-blocks or complement
                            if read_bit(r)  # inner=1: complement — read skipped positions
                                local _skip_set = Set{Int}(Int.(collect(_read_copy_blocks(r, params.varint, T))))
                                for p in 1:ref_len
                                    if p ∉ _skip_set; push!(copied_vals, ref_list[p]); end
                                end
                            else            # inner=0: standard copy-blocks
                                for p in _read_copy_blocks(r, params.varint, T)
                                    if 1 <= p <= ref_len; push!(copied_vals, ref_list[p]); end
                                end
                            end
                        end
                    elseif params.intra_copy_blocks
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
                    if (params.intra_intervals || params.intra_greedy_mil) && params.intra_lr_split
                        additions = Int.(_read_ir_lr(r, :fibonacci, mil_vec[idx_local], T, _zz_vid; tight_deltas=params.intra_tight_deltas))
                    elseif params.intra_intervals || params.intra_greedy_mil
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
                    if (params.intra_intervals || params.intra_greedy_mil) && params.intra_lr_split
                        nl_local = Int.(_read_ir_lr(r, :fibonacci, mil_vec[idx_local], T, _zz_vid; tight_deltas=params.intra_tight_deltas))
                    elseif params.intra_intervals || params.intra_greedy_mil
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
    load_cg_graph(filepath::AbstractString; params::CGParams=CGParams(),
                   T::Type{<:Unsigned}=UInt32, directed::Bool=true)

Load a CG-compressed graph from a file, decoding all levels.
Returns a `Dict{T, Vector{T}}` mapping vertex → sorted outneighbors.
"""
function load_cg_graph(filepath::AbstractString; params::CGParams=CGParams(),
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

end # module CG
