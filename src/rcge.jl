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

using ..IO: BitWriter, write_bit, write_bits, write_value
import ..Compression
using ..Compression: write_encoded_value, write_delta, write_truncated_binary_coding, write_hybrid_mix_encoded_list, delta_encode_vector, write_elias_fano

export RCGEParams, RCGEStats, encode_level

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
    intra_ref_used::Int
    intra_no_ref::Int
    bits_inter_headers::Int
    bits_inter_degrees::Int
    bits_inter_perms::Int
    RCGEStats() = new(0,0,0,0,0,0,0,0,0,0,0)
end

@inline function _total_bits(w::BitWriter)
    # Bytes written to underlying IO + buffered bits
    return Int(position(w.io)) * 8 + (w.index - 1)
end

"""
    RCGEParams(; L=32, varint=:fibonacci, count_varint=:fibonacci, gap=:fibonacci, degree=:golomb, perm_strategy=:lehmer, undirected_pairs=true, membership=:delta, inter_strategy=:perm,
                 intra_ref_enabled::Bool=true, intra_ref_window::Int=16, intra_ref_min_overlap::Float64=0.3, intra_ref_rle::Bool=true)

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
    membership::Symbol = :delta
    inter_strategy::Symbol = :perm
    intra_ref_enabled::Bool = true
    intra_ref_window::Int = 16
    intra_ref_min_overlap::Float64 = 0.35
    intra_ref_rle::Bool = true
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

    # Use provided cluster order for intra/inter encoding; keep a sorted copy for membership encoding
    clusters = [copy(C) for C in P]

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
        Cm = sort(copy(C))
        if params.membership == :elias_fano
            # Elias-Fano encodes its own length; omit explicit |C| when using EF
            write_elias_fano(w, Cm)
        else
            write_varint(w, T(length(Cm)); encoding=params.varint)
            write_gap_coded(w, Cm; encoding=params.gap)
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
        else
            # Adjacency lists in local order, neighbors restricted to C
            # Optional intra-cluster referencing to exploit overlap in successive lists
            prev_lists = Vector{Vector{Int}}()
            # First pass: decide refs and compute ref data
            use_ref_vec = Vector{Bool}(undef, s)
            ref_delta_vec = Vector{UInt32}(undef, s)
            ref_positions_list = Vector{Vector{Int}}(undef, s)
            additions_list = Vector{Vector{Int}}(undef, s)
            raw_lists = Vector{Vector{Int}}(undef, s)
            for (idx_local, u) in enumerate(C)
                # gather local neighbors of u within C, as local indices (1-based)
                nl = Int[]
                for v in outneighbors(g, Int(u))
                    if haskey(local_index, v)
                        push!(nl, local_index[v])
                    end
                end
                sort!(nl)
                # decide reference
                use_ref = false; ref_positions = Int[]; additions = Int[]; ref_delta_val = UInt32(0)
                if params.intra_ref_enabled && !isempty(prev_lists) && !isempty(nl)
                    wstart = max(1, length(prev_lists) - params.intra_ref_window + 1)
                    best_overlap = 0; best_ref_idx = 0; best_ref_positions = Int[]; best_additions = Int[]
                    for rix in wstart:length(prev_lists)
                        ref = prev_lists[rix]
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
                        overlap = length(positions)
                        if overlap > best_overlap
                            best_overlap = overlap; best_ref_idx = rix; best_ref_positions = positions; best_additions = adds
                        end
                    end
                    if best_overlap > 0
                        overlap_ratio = best_overlap / max(1, length(nl))
                        # simple cost check: ensure fewer tokens when using reference
                        if overlap_ratio >= params.intra_ref_min_overlap
                            ref_tokens = length(best_ref_positions) + length(best_additions) + 1  # +1 for positions mode/count
                            raw_tokens = length(nl)
                            if ref_tokens < raw_tokens
                                use_ref = true
                                ref_positions = best_ref_positions
                                additions = best_additions
                                ref_delta_val = UInt32(idx_local - best_ref_idx)
                            end
                        end
                    end
                end
                use_ref_vec[idx_local] = use_ref
                ref_delta_vec[idx_local] = ref_delta_val
                ref_positions_list[idx_local] = ref_positions
                additions_list[idx_local] = additions
                raw_lists[idx_local] = nl
                push!(prev_lists, nl)
            end

            # Write per-cluster ref bitmap and ref_delta RLE
            if any(use_ref_vec)
                local hb0 = _total_bits(w)
                Compression.write_bitpacked_bitmap(w, use_ref_vec)
                # collect deltas for ref-only entries in vertex order
                deltas = UInt32[]
                for idx in 1:s
                    if use_ref_vec[idx]
                        push!(deltas, ref_delta_vec[idx])
                    end
                end
                if params.intra_ref_rle
                    Compression.write_rle_ones_deltas(w, deltas, params.varint)
                else
                    # Write plain count and each delta using varint
                    write_encoded_value(w, T(length(deltas)), params.varint)
                    for d in deltas
                        write_varint(w, d; encoding=params.varint)
                    end
                end
                if stats !== nothing
                    stats.bits_intra_headers += _total_bits(w) - hb0
                end
            end

            # Now write per-vertex payloads without flags or ref_deltas
            for idx_local in 1:s
                if use_ref_vec[idx_local]
                    # write copied positions and additions
                    ref_positions = ref_positions_list[idx_local]
                    additions = additions_list[idx_local]
                    # copied positions count (small-count encoding)
                    Compression.write_small_count(w, T(length(ref_positions)), params.count_varint)
                    # choose encoding for positions: runs vs delta
                    local bpos0 = _total_bits(w)
                    if !isempty(ref_positions)
                        # detect runs of consecutive positions
                        runs = Vector{Tuple{Int,Int}}()
                        i2 = 1
                        total_run_len = 0
                        while i2 <= length(ref_positions)
                            j2 = i2
                            while j2+1 <= length(ref_positions) && ref_positions[j2+1] == ref_positions[j2] + 1
                                j2 += 1
                            end
                            lenrun = j2 - i2 + 1
                            if lenrun >= 3
                                push!(runs, (ref_positions[i2], lenrun))
                                total_run_len += lenrun
                            end
                            i2 = j2 + 1
                        end
                        use_runs = total_run_len * 2 >= length(ref_positions)  # heuristic
                        write_bit(w, use_runs)  # mode bit
                        if use_runs
                            # write number of runs (small-count), then pairs
                            Compression.write_small_count(w, T(length(runs)), params.count_varint)
                            for (st, ln) in runs
                                write_encoded_value(w, T(st), params.varint)
                                write_encoded_value(w, T(ln), params.varint)
                            end
                        else
                            # delta-code positions (Fibonacci)
                            write_delta(w, UInt32.(ref_positions), :fibonacci)
                        end
                    else
                        # no positions, write mode bit 0
                        write_bit(w, false)
                    end
                    local bpos1 = _total_bits(w)
                    Compression.write_small_count(w, T(length(additions)), params.count_varint)
                    if !isempty(additions)
                        write_delta(w, T.(additions), :fibonacci)
                    end
                    local b3 = _total_bits(w)
                    if stats !== nothing
                        stats.bits_intra_copy += (bpos1 - bpos0)
                        stats.bits_intra_add += (b3 - bpos1)
                        stats.intra_ref_used += 1
                    end
                else
                    # raw
                    nl = raw_lists[idx_local]
                    local rb0 = _total_bits(w)
                    # raw count can be larger: keep varint
                    write_encoded_value(w, T(length(nl)), params.count_varint)
                    local rb1 = _total_bits(w)
                    if !isempty(nl)
                        write_delta(w, T.(nl), :fibonacci)
                    end
                    if stats !== nothing
                        stats.bits_intra_headers += (rb1 - rb0)
                        stats.bits_intra_raw += _total_bits(w) - rb1
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

    # Write number of A groups
    bits_before = _total_bits(w)
    write_varint(w, T(length(neighbors_by_A)); encoding=params.varint)
    if stats !== nothing
        stats.bits_inter_headers += _total_bits(w) - bits_before
    end

    for ia in sort(collect(keys(neighbors_by_A)))
        A = clusters[ia]
        Bs = neighbors_by_A[ia]
        # Header: A id, and B list (gap-coded with count)
        bits_before = _total_bits(w)
        write_varint(w, T(ia); encoding=params.varint)
        # write k (shift if needed)
        local kcnt = T(length(Bs))
        if params.count_varint == :elias_delta || params.count_varint == :fibonacci || params.count_varint == :zeta
            kcnt = kcnt + one(T)
        end
        write_encoded_value(w, kcnt, params.count_varint)
        # write sorted B ids as gap-coded
        write_gap_coded(w, T.(Bs); encoding=params.gap)
        if stats !== nothing
            stats.bits_inter_headers += _total_bits(w) - bits_before
        end

        # For each B in list, emit m, degrees, permutation
        for ib in Bs
            B = clusters[ib]
            neighbors_in_B = pair_map[(ia,ib)]
            # Sort neighbor lists for determinism
            for lst in values(neighbors_in_B)
                sort!(lst)
            end
            # Degree vectors and totals
            total_edges = 0
            d_out = T[ zero(T) for _ in 1:length(A) ]
            d_in_counts = zeros(Int, length(B))
            for (idxA, u) in enumerate(A)
                Ns = get(neighbors_in_B, u, Int[])
                d_out[idxA] = T(length(Ns))
                total_edges += length(Ns)
                for vloc in Ns
                    d_in_counts[vloc] += 1
                end
            end
            # Header: m
            bits_before = _total_bits(w)
            write_varint(w, T(total_edges); encoding=params.varint)
            if stats !== nothing
                stats.bits_inter_headers += _total_bits(w) - bits_before
            end

            if params.inter_strategy == :perm
                # Degrees sparse
                bits_before = _total_bits(w)
                write_sparse_int_vector(w, d_out; index_encoding=params.gap, value_encoding=params.degree)
                write_sparse_int_vector(w, T.(d_in_counts); index_encoding=params.gap, value_encoding=params.degree)
                if stats !== nothing
                    stats.bits_inter_degrees += _total_bits(w) - bits_before
                end
                # Permutation (blockpos uses d_in_counts)
                pi = order_from_edges(A, B, neighbors_in_B)
                bits_before = _total_bits(w)
                write_permutation(w, pi; strategy=params.perm_strategy, pos_encoding=params.gap, count_encoding=params.count_varint, block_sizes=Int.(d_in_counts))
                if stats !== nothing
                    stats.bits_inter_perms += _total_bits(w) - bits_before
                end
            else
                # Explicit per-u neighbor lists in B (CSR-like); drop degrees/permutation
                # For each u in A local order, write count and gap-coded list
                for u in A
                    Ns = get(neighbors_in_B, u, Int[])
                    # count (shift if zero-unsafe)
                    local cntu::T = T(length(Ns))
                    if params.count_varint == :elias_delta || params.count_varint == :fibonacci || params.count_varint == :zeta
                        cntu = cntu + one(T)
                    end
                    write_encoded_value(w, cntu, params.count_varint)
                    if !isempty(Ns)
                        write_delta(w, T.(Ns), params.gap)
                    end
                end
            end
        end
    end

    return nothing
end

end # module RCGE
