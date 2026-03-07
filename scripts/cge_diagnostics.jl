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

#!/usr/bin/env julia

#
# CGE Diagnostics for CNR-2000
#
# Runs CGE encoding with the best config (zigzag + stop_deltas + copy_blocks)
# and collects detailed per-cluster diagnostics:
#   - Ref delta histogram and statistics
#   - Copy-block count per ref vertex
#   - Average raw / additions list sizes
#   - Per-component bit costs
#

using Pkg
Pkg.activate(normpath(joinpath(@__DIR__, "..")))

using Logging
global_logger(ConsoleLogger(stderr, Logging.Info))

using LightGraphs
using LightGraphs: nv, outneighbors, is_directed
using Adjacently
using Adjacently.CGE: encode_level, CGEParams, CGEStats, decode_level
using Adjacently.IO: load_adjacency_list_from_csv, BitWriter, BitReader, flush_bitwriter
using Adjacently.Clustering: leiden_partition, aggregate_graph
using Adjacently.Relabeling: relabel_vertices_llp
using Adjacently.Graph: subgraph

# ---------------------------------------------------------------------------
# Minimal TestGraph for coarse weighted graph (reused from test_cge_cnr2000.jl)
# ---------------------------------------------------------------------------
struct TestGraph{T<:Unsigned} <: LightGraphs.AbstractGraph{T}
    n::Int
    adj::Dict{Int, Vector{T}}
    radj::Dict{Int, Vector{T}}
    directed::Bool
end
struct TestEdge{T<:Unsigned}
    u::T
    v::T
end
LightGraphs.nv(h::TestGraph) = h.n
LightGraphs.is_directed(h::TestGraph) = h.directed
LightGraphs.outneighbors(h::TestGraph{T}, v::Integer) where {T<:Unsigned} = get(h.adj, Int(v), T[])
LightGraphs.src(e::TestEdge) = e.u
LightGraphs.dst(e::TestEdge) = e.v
function LightGraphs.edges(h::TestGraph{T}) where {T<:Unsigned}
    ed = TestEdge{T}[]
    for (u, lst) in h.adj
        for v in lst
            push!(ed, TestEdge{T}(T(u), v))
        end
    end
    return ed
end
LightGraphs.inneighbors(h::TestGraph{T}, v::Integer) where {T<:Unsigned} = get(h.radj, Int(v), T[])

function coarse_to_testgraph(Gc)::TestGraph{UInt32}
    nC = Gc.n
    adj = Dict{Int, Vector{UInt32}}()
    radj = Dict{Int, Vector{UInt32}}()
    for u in 1:nC
        lst = UInt32[]
        for (v, w) in Gc.out_w[u]
            for _ in 1:Int(round(w))
                push!(lst, UInt32(v))
                push!(get!(radj, v, UInt32[]), UInt32(u))
            end
        end
        !isempty(lst) && (adj[u] = lst)
    end
    return TestGraph{UInt32}(nC, adj, radj, true)
end

count_edges(h) = sum(length(outneighbors(h, v)) for v in 1:nv(h))

# ---------------------------------------------------------------------------
# Reorder within clusters using LLP on induced subgraph
# ---------------------------------------------------------------------------
function reorder_clusters!(clusters, base_g; llp_passes::Int=5)
    for idx in 1:length(clusters)
        C = clusters[idx]
        length(C) <= 2 && continue
        sg, oni, _ = subgraph(base_g, C)
        mapping = relabel_vertices_llp(sg, :sym; passes=llp_passes)
        sort!(C, by = v -> Int(mapping[oni[v]]))
        clusters[idx] = C
    end
    return clusters
end

# ---------------------------------------------------------------------------
# Diagnostics collection structure
# ---------------------------------------------------------------------------
mutable struct DiagnosticCollector
    # Ref delta values (one per ref vertex)
    ref_deltas::Vector{Int}
    # Copy block counts per ref vertex
    copy_block_counts::Vector{Int}
    # Raw list lengths (vertices that did NOT use a reference)
    raw_list_lengths::Vector{Int}
    # Additions list lengths (vertices that DID use a reference)
    additions_list_lengths::Vector{Int}
    # Total list lengths (neighbor count for every vertex, before ref decision)
    total_list_lengths::Vector{Int}
    # Cluster sizes
    cluster_sizes::Vector{Int}

    DiagnosticCollector() = new(Int[], Int[], Int[], Int[], Int[], Int[])
end

# ---------------------------------------------------------------------------
# Diagnostic CGE run: encode + collect per-vertex reference stats
# ---------------------------------------------------------------------------
function run_cge_diagnostics(g, m_original, params; K::Int=8, llp_passes::Int=5)
    max_levels = 5
    min_clusters = 32
    cur_g = g
    total_bytes = 0
    prev_coarse_n = -1

    # Accumulate stats and diagnostics across levels
    all_stats = CGEStats[]
    all_diag = DiagnosticCollector[]
    level_info = []

    println("\nStarting CGE diagnostic encoding (K=$K, max_levels=$max_levels)...")

    for level in 1:max_levels
        ncur = nv(cur_g)
        mcur = count_edges(cur_g)
        println("\n  Level $level: n=$ncur, m=$mcur")

        # Partition
        t1 = time()
        part = if level == 1
            leiden_partition(cur_g; max_passes=8, max_levels=5)
        else
            leiden_partition(cur_g; max_passes=5, max_levels=5)
        end
        nclusters = maximum(part)
        println("    Partitioned into $nclusters clusters ($(round(time()-t1, digits=3))s)")

        # Cap to top-K by size
        counts = Dict{Int,Int}()
        for c in part
            counts[c] = get(counts, c, 0) + 1
        end
        labels_sorted = sort(collect(keys(counts)), by = c -> -counts[c])
        topK = min(K, length(labels_sorted))
        top = labels_sorted[1:topK]
        top_index = Dict{Int,Int}(c => i for (i,c) in enumerate(top))

        Vcur = (typeof(cur_g)).parameters[1]
        clusters = [Vcur[] for _ in 1:(topK + 1)]
        capped_part = similar(part)
        for i in 1:ncur
            c = part[i]
            bucket = get(top_index, c, 0)
            if bucket == 0
                bucket = topK + 1
            end
            capped_part[i] = bucket
            push!(clusters[bucket], Vcur(i))
        end
        clusters = filter(!isempty, clusters)
        println("    Effective clusters: $(length(clusters)) (top=$topK)")

        # Reorder within clusters using LLP
        reorder_clusters!(clusters, cur_g; llp_passes=llp_passes)

        # Prepare next level's K
        K = max(16, min(K, ceil(Int, nclusters / 2)))

        # ---- Encode with stats ----
        io = IOBuffer()
        w = BitWriter(io)
        t2 = time()
        stats = CGEStats()
        encode_level(w, cur_g, clusters; params=params, stats=stats)
        flush_bitwriter(w; flush_last_bits=true)
        bytes = take!(io)
        t3 = time()
        push!(all_stats, stats)

        level_bytes = length(bytes)
        total_bytes += level_bytes
        bpe = 8.0 * level_bytes / max(mcur, 1)
        cum_bpe = 8.0 * total_bytes / max(m_original, 1)

        push!(level_info, (level=level, ncur=ncur, mcur=mcur, level_bytes=level_bytes, bpe=bpe, cum_bpe=cum_bpe))
        println("    Encoded: $(level_bytes) bytes, BPE=$(round(bpe, digits=4)), cumulative=$(round(cum_bpe, digits=4)) ($(round(t3-t2, digits=3))s)")

        # ---- Diagnostic pass: re-run reference matching to collect stats ----
        diag = DiagnosticCollector()
        sorted_clusters = [sort(copy(C)) for C in clusters]

        for C in sorted_clusters
            s = length(C)
            push!(diag.cluster_sizes, s)

            # Build local index map
            local_index = Dict{eltype(C),Int}()
            for (i, u) in enumerate(C)
                local_index[u] = i
            end

            # Skip bitset clusters (small undirected); they don't have ref decisions
            if s <= params.L && !is_directed(cur_g)
                continue
            end

            # Replay per-vertex reference matching (mirrors encode_level logic)
            prev_lists = Vector{Vector{Int}}()

            for (idx_local, u) in enumerate(C)
                # Gather local neighbors
                nl = Int[]
                for v in outneighbors(cur_g, Int(u))
                    if haskey(local_index, v)
                        push!(nl, local_index[v])
                    end
                end
                sort!(nl)

                push!(diag.total_list_lengths, length(nl))

                use_ref = false
                best_positions = Int[]
                best_additions = Int[]
                best_delta = 0

                if params.intra_ref_enabled && !isempty(prev_lists)
                    # Estimate raw cost
                    _zz_vid = params.intra_zigzag ? UInt32(idx_local) : nothing
                    io_raw = IOBuffer(); w_raw = BitWriter(io_raw)
                    if params.intra_stop_deltas
                        Adjacently.CGE._write_stop_delta_zigzag(w_raw, UInt32.(nl), :fibonacci, _zz_vid)
                    else
                        Adjacently.Compression.write_small_count(w_raw, UInt32(length(nl)), params.count_varint)
                        if !isempty(nl)
                            Adjacently.Compression.write_delta(w_raw, UInt32.(nl), :fibonacci)
                        end
                    end
                    flush_bitwriter(w_raw; flush_last_bits=true)
                    raw_bits = length(take!(io_raw)) * 8
                    best_bits = raw_bits
                    best_idx = 0

                    wstart = max(1, length(prev_lists) - params.intra_ref_window + 1)
                    for rix in wstart:length(prev_lists)
                        ref = prev_lists[rix]
                        # Two-pointer intersection
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

                        # Estimate ref cost
                        io_ref = IOBuffer(); w_ref = BitWriter(io_ref)
                        if params.intra_copy_blocks
                            Adjacently.CGE._write_copy_blocks(w_ref, positions, params.varint)
                        else
                            Adjacently.Compression.write_small_count(w_ref, UInt32(length(positions)), params.count_varint)
                            if !isempty(positions)
                                Adjacently.Compression.write_delta(w_ref, UInt32.(positions), :fibonacci)
                            end
                        end
                        if params.intra_stop_deltas
                            Adjacently.CGE._write_stop_delta_zigzag(w_ref, UInt32.(adds), :fibonacci, _zz_vid)
                        else
                            Adjacently.Compression.write_small_count(w_ref, UInt32(length(adds)), params.count_varint)
                            if !isempty(adds)
                                Adjacently.Compression.write_delta(w_ref, UInt32.(adds), :fibonacci)
                            end
                        end
                        flush_bitwriter(w_ref; flush_last_bits=true)
                        ref_bits = length(take!(io_ref)) * 8

                        if ref_bits < best_bits
                            best_bits = ref_bits
                            best_idx = rix
                            best_positions = positions
                            best_additions = adds
                        end
                    end

                    if best_idx > 0 && best_bits < raw_bits
                        use_ref = true
                        best_delta = idx_local - best_idx
                    end
                end

                if use_ref
                    push!(diag.ref_deltas, best_delta)
                    push!(diag.additions_list_lengths, length(best_additions))

                    # Count copy blocks in positions
                    nblocks = 0
                    if !isempty(best_positions)
                        nblocks = 1
                        for k in 2:length(best_positions)
                            if best_positions[k] != best_positions[k-1] + 1
                                nblocks += 1
                            end
                        end
                    end
                    push!(diag.copy_block_counts, nblocks)
                else
                    push!(diag.raw_list_lengths, length(nl))
                end

                push!(prev_lists, nl)
            end
        end
        push!(all_diag, diag)

        # Build coarse graph and check stopping
        t4 = time()
        Gc = aggregate_graph(cur_g, capped_part)
        println("    Coarsened to n=$(Gc.n) ($(round(time()-t4, digits=3))s)")

        if prev_coarse_n == Gc.n
            println("    Stopping: coarse size unchanged")
            break
        end
        prev_coarse_n = Gc.n
        if Gc.n <= min_clusters
            println("    Stopping: coarse size $(Gc.n) <= min_clusters $min_clusters")
            break
        end

        cur_g = coarse_to_testgraph(Gc)
    end

    return all_stats, all_diag, level_info, total_bytes
end

# ---------------------------------------------------------------------------
# Pretty-print diagnostics
# ---------------------------------------------------------------------------
function print_histogram(vals::Vector{Int}; max_bucket::Int=32, title::String="Histogram")
    isempty(vals) && (println("  (no data)"); return)
    println("\n  $title (n=$(length(vals))):")
    bucket_counts = Dict{Int,Int}()
    for v in vals
        bucket = min(v, max_bucket)
        bucket_counts[bucket] = get(bucket_counts, bucket, 0) + 1
    end
    total = length(vals)
    max_count = maximum(Base.values(bucket_counts))
    bar_width = 40

    for k in sort(collect(keys(bucket_counts)))
        c = bucket_counts[k]
        pct = 100.0 * c / total
        bar_len = round(Int, bar_width * c / max_count)
        label = k == max_bucket ? "$k+" : string(k)
        bar = repeat("#", bar_len)
        println("    $(lpad(label, 4)): $(lpad(c, 7))  ($(lpad(round(pct, digits=1), 5))%)  $bar")
    end
end

function print_stats_summary(stats::CGEStats, info)
    mcur = info.mcur
    println("\n  --- Bit Cost Breakdown (Level $(info.level)) ---")
    total_bits = info.level_bytes * 8

    membership_bits = stats.bits_membership
    intra_bits = stats.bits_intra
    intra_header_bits = stats.bits_intra_headers
    intra_copy_bits = stats.bits_intra_copy
    intra_add_bits = stats.bits_intra_add
    intra_raw_bits = stats.bits_intra_raw
    inter_header_bits = stats.bits_inter_headers
    inter_list_bits = stats.bits_inter_lists
    inter_total = inter_header_bits + inter_list_bits

    # The ref_small_headers field tracks ref bitmap + ref delta bits (part of intra_headers)
    ref_bitmap_delta_bits = intra_header_bits

    println("    Total:       $(lpad(total_bits, 10)) bits  ($(round(total_bits / 8, digits=0)) bytes)")
    println("    Membership:  $(lpad(membership_bits, 10)) bits  ($(round(100.0 * membership_bits / total_bits, digits=1))%)")
    println("    Intra:       $(lpad(intra_bits, 10)) bits  ($(round(100.0 * intra_bits / total_bits, digits=1))%)")
    println("      Headers (ref_bitmap + ref_deltas): $(lpad(ref_bitmap_delta_bits, 8)) bits")
    println("      Copy (positions):                  $(lpad(intra_copy_bits, 8)) bits")
    println("      Additions:                         $(lpad(intra_add_bits, 8)) bits")
    println("      Raw:                               $(lpad(intra_raw_bits, 8)) bits")
    println("    Inter:       $(lpad(inter_total, 10)) bits  ($(round(100.0 * inter_total / total_bits, digits=1))%)")
    println("      Headers:   $(lpad(inter_header_bits, 8)) bits")
    println("      Lists:     $(lpad(inter_list_bits, 8)) bits")
    println("    Ref usage:   $(stats.intra_ref_used) ref / $(stats.intra_ref_used + stats.intra_no_ref) total vertices")
    if stats.intra_ref_used + stats.intra_no_ref > 0
        ref_pct = 100.0 * stats.intra_ref_used / (stats.intra_ref_used + stats.intra_no_ref)
        println("    Ref rate:    $(round(ref_pct, digits=1))%")
    end

    # Per-edge cost
    if mcur > 0
        println("    Membership BPE:  $(round(membership_bits / mcur, digits=4))")
        println("    Intra BPE:       $(round(intra_bits / mcur, digits=4))")
        println("    Inter BPE:       $(round(inter_total / mcur, digits=4))")
    end
end

function print_diagnostics(diag::DiagnosticCollector, level::Int)
    println("\n  === Diagnostic Details (Level $level) ===")

    # Cluster sizes
    if !isempty(diag.cluster_sizes)
        sizes = diag.cluster_sizes
        println("\n  Cluster sizes: $(length(sizes)) clusters")
        println("    Min: $(minimum(sizes)), Max: $(maximum(sizes)), Mean: $(round(sum(sizes)/length(sizes), digits=1)), Median: $(sort(sizes)[div(length(sizes),2)+1])")
    end

    # Ref delta histogram
    if !isempty(diag.ref_deltas)
        deltas = diag.ref_deltas
        println("\n  Ref delta statistics ($(length(deltas)) ref vertices):")
        println("    Min: $(minimum(deltas)), Max: $(maximum(deltas)), Mean: $(round(sum(deltas)/length(deltas), digits=2))")

        # Detailed counts for small deltas
        println("\n  Ref delta value counts:")
        for d in 1:min(10, maximum(deltas))
            c = count(==(d), deltas)
            if c > 0
                println("    delta=$d: $c ($(round(100.0*c/length(deltas), digits=1))%)")
            end
        end

        print_histogram(deltas; max_bucket=32, title="Ref delta distribution")
    else
        println("\n  No ref vertices found (all vertices encoded raw)")
    end

    # Copy block counts
    if !isempty(diag.copy_block_counts)
        cbc = diag.copy_block_counts
        println("\n  Copy-block counts per ref vertex ($(length(cbc)) vertices):")
        println("    Min: $(minimum(cbc)), Max: $(maximum(cbc)), Mean: $(round(sum(cbc)/length(cbc), digits=2))")
        print_histogram(cbc; max_bucket=16, title="Copy-block count distribution")
    end

    # Raw list lengths
    if !isempty(diag.raw_list_lengths)
        rll = diag.raw_list_lengths
        println("\n  Raw list lengths (non-ref vertices): $(length(rll)) vertices")
        println("    Min: $(minimum(rll)), Max: $(maximum(rll)), Mean: $(round(sum(rll)/length(rll), digits=2))")
    end

    # Additions list lengths
    if !isempty(diag.additions_list_lengths)
        all_len = diag.additions_list_lengths
        println("\n  Additions list lengths (ref vertices): $(length(all_len)) vertices")
        println("    Min: $(minimum(all_len)), Max: $(maximum(all_len)), Mean: $(round(sum(all_len)/length(all_len), digits=2))")
        # Count how many additions lists are empty (pure copy)
        n_empty = count(==(0), all_len)
        println("    Pure copies (additions=0): $n_empty ($(round(100.0*n_empty/length(all_len), digits=1))%)")
    end

    # Total list lengths
    if !isempty(diag.total_list_lengths)
        tll = diag.total_list_lengths
        println("\n  Total neighbor list lengths (all vertices): $(length(tll)) vertices")
        println("    Min: $(minimum(tll)), Max: $(maximum(tll)), Mean: $(round(sum(tll)/length(tll), digits=2))")
        n_empty_total = count(==(0), tll)
        println("    Empty lists (degree-0 in cluster): $n_empty_total")
    end
end

# ---------------------------------------------------------------------------
# Main
# ---------------------------------------------------------------------------
println("=" ^ 70)
println("CGE CNR-2000 Diagnostics")
println("=" ^ 70)

cnr_csv = joinpath(@__DIR__, "..", "datasets", "webgraph", "cnr-2000", "cnr-2000.csv")
cnr_csv = normpath(cnr_csv)
if !isfile(cnr_csv)
    error("CNR-2000 CSV not found at $cnr_csv")
end

println("\nLoading CNR-2000...")
t0 = time()
g = load_adjacency_list_from_csv(cnr_csv, ',', true)
n = nv(g)
m_original = count_edges(g)
println("  Loaded: $n vertices, $m_original edges ($(round(time()-t0, digits=2))s)")

# Best config: zigzag + stop_deltas + copy_blocks
function make_params(; kwargs...)
    d = Dict{Symbol,Any}(
        :L=>128, :varint=>:fibonacci, :count_varint=>:fibonacci, :gap=>:fibonacci,
        :degree=>:elias_delta, :undirected_pairs=>false, :perm_strategy=>:blockpos,
        :membership=>:elias_fano, :inter_strategy=>:perm, :intra_ref_enabled=>true,
        :intra_ref_window=>32, :intra_ref_rle=>false,
        :intra_block_try=>false, :positions_mode=>:delta, :additions_mode=>:delta,
        :min_cluster_density=>0.0, :intra_intervals=>false, :intra_mil=>4,
        :intra_greedy_mil=>false, :intra_zigzag=>true, :intra_stop_deltas=>true,
        :intra_copy_blocks=>true,
    )
    for (k, v) in kwargs
        d[k] = v
    end
    return CGEParams(; d...)
end

params = make_params()

println("\nConfig: zigzag=true, stop_deltas=true, copy_blocks=true")
println("  ref_window=$(params.intra_ref_window), L=$(params.L), varint=$(params.varint)")
println("  membership=$(params.membership), ref_rle=$(params.intra_ref_rle)")

all_stats, all_diag, level_info, total_bytes = run_cge_diagnostics(g, m_original, params; K=8, llp_passes=5)

# Print per-level results
for (i, info) in enumerate(level_info)
    println("\n" * "=" ^ 70)
    println("LEVEL $(info.level) DIAGNOSTICS")
    println("=" ^ 70)

    print_stats_summary(all_stats[i], info)
    print_diagnostics(all_diag[i], info.level)
end

# Final summary
println("\n" * "=" ^ 70)
println("OVERALL SUMMARY")
println("=" ^ 70)
println("  Total bytes:  $total_bytes ($(round(total_bytes / 1024.0, digits=1)) KB)")
println("  Total BPE:    $(round(8.0 * total_bytes / m_original, digits=4))")
println("  Levels:       $(length(level_info))")

# Aggregate diagnostics across all levels
all_ref_deltas = vcat([d.ref_deltas for d in all_diag]...)
all_copy_blocks = vcat([d.copy_block_counts for d in all_diag]...)
all_raw_lens = vcat([d.raw_list_lengths for d in all_diag]...)
all_add_lens = vcat([d.additions_list_lengths for d in all_diag]...)

if !isempty(all_ref_deltas)
    println("\n  Aggregate ref deltas (all levels): $(length(all_ref_deltas)) ref vertices")
    println("    Mean delta: $(round(sum(all_ref_deltas)/length(all_ref_deltas), digits=2))")
    println("    Median delta: $(sort(all_ref_deltas)[div(length(all_ref_deltas),2)+1])")
end
if !isempty(all_copy_blocks)
    println("  Aggregate copy-blocks/vertex: mean=$(round(sum(all_copy_blocks)/length(all_copy_blocks), digits=2))")
end
if !isempty(all_raw_lens)
    println("  Aggregate raw list length: mean=$(round(sum(all_raw_lens)/length(all_raw_lens), digits=2))")
end
if !isempty(all_add_lens)
    println("  Aggregate additions list length: mean=$(round(sum(all_add_lens)/length(all_add_lens), digits=2))")
end
println("=" ^ 70)
