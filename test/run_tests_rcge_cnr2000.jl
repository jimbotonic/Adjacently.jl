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

#!/usr/bin/env julia

# RCGE multi-level coarsening + encoding on CNR-2000 graph

include("run_tests_main.jl")
using Logging
using Adjacently.RCGE: encode_level, RCGEParams
using Adjacently.IO: load_adjacency_list_from_csv, BitWriter, flush_bitwriter
using Adjacently.Clustering: louvain_partition, leiden_partition, aggregate_graph
using Adjacently.Relabeling: relabel_vertices_rcm, relabel_vertices_llp, relabel_vertices_minhash
using Adjacently.Graph: subgraph
using LightGraphs
using LightGraphs: nv, outneighbors

@testset "RCGE RCGE multi-level encode on CNR-2000" begin
    prev_logger = current_logger()
    global_logger(ConsoleLogger(stderr, Logging.Info))
    try
        cnr_csv = "datasets/webgraph/cnr-2000/cnr-2000.csv"
        if !isfile(cnr_csv)
            @warn "CNR-2000 CSV not found at $cnr_csv; skipping RCGE test"
            @test_skip "CNR-2000 dataset unavailable"
            return
        end

        @info "Loading CNR-2000 graph..."
        g = load_adjacency_list_from_csv(cnr_csv, ',', true)
        n = nv(g)
        @info "Graph loaded: n=$n"

        # Relabeling flags (tune cluster-internal ordering)
        # Available methods: :none, :rcm, :llp, :minhash
        RELABEL_METHOD = Symbol(get(ENV, "RCGE_RELABEL_METHOD", "none"))
        RELABEL_NEIGHBOR = Symbol(get(ENV, "RCGE_RELABEL_NEIGHBOR", "sym"))  # :out or :sym
        LLP_PASSES = (try parse(Int, get(ENV, "RCGE_LLP_PASSES", "5")) catch; 5 end)
        MINHASH_K = (try parse(Int, get(ENV, "RCGE_MINHASH_K", "32")) catch; 32 end)
        @info "Relabeling config: method=$(RELABEL_METHOD), neighbor=$(RELABEL_NEIGHBOR), llp_passes=$(LLP_PASSES), minhash_k=$(MINHASH_K)"

        # Utility: count directed edges for a LightGraphs graph-like API
        count_edges(h) = sum(length(outneighbors(h, v)) for v in 1:nv(h))

        # Convert WeightedCoarseGraph to a minimal LightGraphs-like TestGraph
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
                    # Repeat neighbor v 'w' times (w is float but integer-valued)
                    for _ in 1:Int(round(w))
                        push!(lst, UInt32(v))
                        push!(get!(radj, v, UInt32[]), UInt32(u))
                    end
                end
                !isempty(lst) && (adj[u] = lst)
            end
            return TestGraph{UInt32}(nC, adj, radj, true)
        end

        # Parameters
        max_levels = 5
        min_clusters = 32
        # Use Fibonacci for positive-only fields and Elias-delta(+1) for zero-allowing fields
        INTER_STRATEGY = Symbol(get(ENV, "RCGE_INTER", "perm"))
        params = RCGEParams(L=128, varint=:fibonacci, count_varint=:elias_delta, gap=:fibonacci, degree=:elias_delta, undirected_pairs=false, perm_strategy=:blockpos, membership=:elias_fano, inter_strategy=INTER_STRATEGY)

        # Helper: reorder vertices inside each cluster using RCM on the induced subgraph
        function reorder_clusters!(clusters, base_g)
            for idx in 1:length(clusters)
                C = clusters[idx]
                if length(C) <= 2
                    continue
                end
                # build subgraph induced by C
                sg, oni, _ = subgraph(base_g, C)
                # relabel subgraph vertices according to selected method
                local mapping
                if RELABEL_METHOD == :rcm
                    mapping = relabel_vertices_rcm(sg, RELABEL_NEIGHBOR)
                elseif RELABEL_METHOD == :llp
                    mapping = relabel_vertices_llp(sg, RELABEL_NEIGHBOR; passes=LLP_PASSES)
                elseif RELABEL_METHOD == :minhash
                    mapping = relabel_vertices_minhash(sg, RELABEL_NEIGHBOR; k=MINHASH_K)
                else
                    # no relabeling: identity mapping on subgraph local ids
                    mapping = Dict{eltype(C), eltype(C)}()
                    for lid in values(oni)
                        mapping[lid] = lid
                    end
                end
                # sort original vertex ids in C by their new ids in the subgraph mapping
                sort!(C, by = v -> Int(mapping[oni[v]]))
                clusters[idx] = C
            end
            return clusters
        end

        # Quick sweeps: RLE enabled/disabled for intra ref_deltas (level-1 only)
        INTER = Symbol(get(ENV, "RCGE_SWEEP_INTER", "perm"))
        METHOD = :none
        part1 = leiden_partition(g; max_passes=8, max_levels=5)
        counts = Dict{Int,Int}(); for c in part1; counts[c] = get(counts,c,0)+1; end
        labels_sorted = sort(collect(keys(counts)), by = c -> -counts[c])
        K1 = 8
        topK = min(K1, length(labels_sorted))
        top = labels_sorted[1:topK]
        top_index = Dict{Int,Int}(c => i for (i,c) in enumerate(top))
        Vcur = (typeof(g)).parameters[1]
        clusters1 = [Vcur[] for _ in 1:(topK + 1)]
        other_label = topK + 1
        capped_part = similar(part1)
        for i in 1:n
            c = part1[i]
            bucket = get(top_index, c, other_label)
            capped_part[i] = bucket
            push!(clusters1[bucket], Vcur(i))
        end
        clusters1 = filter(!isempty, clusters1)
        @info "RLE sweep base: K1=$(K1) inter=$(INTER) clusters after capping=$(length(clusters1))"
        for rle in (false, true)
            io = IOBuffer(); w = BitWriter(io)
            params1 = RCGEParams(L=128, varint=:fibonacci, count_varint=:elias_delta, gap=:fibonacci, degree=:elias_delta, undirected_pairs=false, perm_strategy=:blockpos, membership=:elias_fano, inter_strategy=INTER, intra_ref_enabled=true, intra_ref_window=32, intra_ref_min_overlap=0.3, intra_ref_rle=rle)
            stats1 = Adjacently.RCGE.RCGEStats()
            tenc0 = time(); encode_level(w, g, clusters1; params=params1, stats=stats1); flush_bitwriter(w; flush_last_bits=true); bytes1 = take!(io); tenc1 = time()
            bpe1 = 8.0 * length(bytes1) / max(count_edges(g), 1)
            mb = ceil(Int, stats1.bits_membership / 8)
            ib = ceil(Int, stats1.bits_intra / 8)
            ihe = ceil(Int, stats1.bits_intra_headers / 8)
            icp = ceil(Int, stats1.bits_intra_copy / 8)
            iad = ceil(Int, stats1.bits_intra_add / 8)
            irw = ceil(Int, stats1.bits_intra_raw / 8)
            hb = ceil(Int, stats1.bits_inter_headers / 8)
            db = ceil(Int, stats1.bits_inter_degrees / 8)
            pb = ceil(Int, stats1.bits_inter_perms / 8)
            @info "RLE=$(rle) inter=$(INTER): encode_time=$(round(tenc1-tenc0,digits=3))s size=$(length(bytes1)) bytes, bpe=$(round(bpe1,digits=4))"
            @info "RLE=$(rle) inter=$(INTER): sections (bytes): membership=$(mb), intra=$(ib) [headers=$(ihe), copy=$(icp), add=$(iad), raw=$(irw)], inter_headers=$(hb), inter_degrees=$(db), inter_perms=$(pb)"
            @info "RLE=$(rle) inter=$(INTER): intra refs used=$(stats1.intra_ref_used), no_ref=$(stats1.intra_no_ref)"
        end

        # Iterate levels until threshold
        cur_g = g
        m_original = count_edges(g)
        total_bytes = 0
        prev_coarse_n = -1
        level = 1
        # Initial K1 can be overridden via env var RCGE_K1; defaults to 64
        K = try parse(Int, get(ENV, "RCGE_K1", "8")) catch; 8 end  # initial cap on clusters for level 1 (then adaptive)
        @info "Initial K1 for multi-level: $(K)"
        while level <= max_levels
            ncur = nv(cur_g)
            mcur = count_edges(cur_g)
            @info "Level $level: n=$ncur m=$mcur"

            # Partition current graph via Louvain
            t0 = time()
            local part
            if level == 1
                # Stronger refinement at level 1
                part = leiden_partition(cur_g; max_passes=8, max_levels=5)
            else
                # Try Leiden on deeper levels to improve coarsening
                part = leiden_partition(cur_g; max_passes=5, max_levels=5)
            end
            t1 = time()
            nclusters = maximum(part)
            @info "Partitioned level $level into $(nclusters) clusters in $(round(t1-t0,digits=3))s"
            # Cap communities to top-K by size to control inter-cluster pairs
            @info "Capping clusters to K=$(K) for encoding"
            counts = Dict{Int,Int}()
            for c in part
                counts[c] = get(counts, c, 0) + 1
            end
            labels_sorted = sort(collect(keys(counts)), by = c -> -counts[c])
            topK = min(K, length(labels_sorted))
            top = labels_sorted[1:topK]
            top_index = Dict{Int,Int}(c => i for (i,c) in enumerate(top))
            # For remaining small communities, avoid a single giant 'other': spread into M bins
            M_OTHERS = 1  # single 'other' cluster appears best empirically
            other_offset = topK
            # Build capped partition (1..topK+M_OTHERS) and clusters
            Vcur = (typeof(cur_g)).parameters[1]
            clusters = [Vcur[] for _ in 1:(topK + M_OTHERS)]
            capped_part = similar(part)
            rest_order = 0
            for i in 1:ncur
                c = part[i]
                bucket = get(top_index, c, 0)
                if bucket == 0
                    # assign to one of M_OTHERS bins in round-robin to reduce mixing
                    rest_order += 1
                    bucket = other_offset + ((rest_order - 1) % M_OTHERS) + 1
                end
                capped_part[i] = bucket
                push!(clusters[bucket], Vcur(i))
            end
            clusters = filter(!isempty, clusters)
            @info "Effective clusters after binning: $(length(clusters)) (top=$(topK), bins=$(M_OTHERS))"
            # Reorder within clusters
            reorder_clusters!(clusters, cur_g)

            # Prepare next level's K adaptively (halve, but not below 16, and not above current K)
            K = max(16, min(K, ceil(Int, nclusters / 2)))

            # Encode RCGE level and compute stats
            io = IOBuffer(); w = BitWriter(io)
            t2 = time()
            stats = Adjacently.RCGE.RCGEStats()
            encode_level(w, cur_g, clusters; params=params, stats=stats)
            flush_bitwriter(w; flush_last_bits=true)
            bytes = take!(io)
            t3 = time()
            @test length(bytes) > 0
            level_bytes = length(bytes)
            total_bytes += level_bytes
            bpe = 8.0 * level_bytes / max(mcur, 1)
            cum_bpe = 8.0 * total_bytes / max(m_original, 1)
            # Sectional bytes
            memb_b = ceil(Int, stats.bits_membership / 8)
            intra_b = ceil(Int, stats.bits_intra / 8)
            head_b = ceil(Int, stats.bits_inter_headers / 8)
            deg_b  = ceil(Int, stats.bits_inter_degrees / 8)
            perm_b = ceil(Int, stats.bits_inter_perms / 8)
            @info "RCGE Level $(level): encode_time=$(round(t3-t2,digits=3))s size=$(level_bytes) bytes, bits/edge=$(round(bpe, digits=4)), cumulative_bits/edge=$(round(cum_bpe, digits=4))"
            @info "  Sections (bytes): membership=$(memb_b), intra=$(intra_b), inter_headers=$(head_b), inter_degrees=$(deg_b), inter_perms=$(perm_b)"

            # Build coarse graph and check threshold
            t4 = time(); Gc = aggregate_graph(cur_g, capped_part); t5 = time()
            @info "Aggregated to coarse graph: n=$(Gc.n) in $(round(t5-t4,digits=3))s"
            # Early stopping: no change in coarse size
            if prev_coarse_n == Gc.n
                @info "Stopping: coarse size unchanged ($(Gc.n)) from previous level"
                break
            end
            prev_coarse_n = Gc.n
            if Gc.n <= min_clusters
                @info "Stopping: coarse communities $(Gc.n) <= min_clusters $(min_clusters)"
                break
            end
            # Convert coarse weighted graph to TestGraph for next iteration
            t6 = time(); cur_g = coarse_to_testgraph(Gc); t7 = time()
            @info "Converted coarse to TestGraph in $(round(t7-t6,digits=3))s"
            level += 1
        end
    finally
        global_logger(prev_logger)
    end
end
