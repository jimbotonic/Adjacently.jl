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

#
# Compare vertex reordering strategies on CNR-2000:
#  - Lexicographic (built-in)
#  - Label Propagation (LLP) – simple implementation
#  - MinHash shingling on adjacency sets
#  - BFS/RCM (approximate RCM on directed graph treated as undirected)
#

include("run_tests_main.jl")

using Logging, Random
using Adjacently.MGS: write_compressed_mgs3_graph
using Adjacently.IO: load_adjacency_list_from_csv
using Adjacently.Graph: get_neighbor_lists
using Adjacently.Relabeling: relabel_vertices, relabel_graph, relabel_vertices_rcm

@testset "CNR-2000 Ordering Comparison" begin
    # Verbose logging for clarity
    prev_logger = current_logger()
    global_logger(ConsoleLogger(stderr, Logging.Info))
    try
        @info "=== CNR-2000 Ordering Comparison ==="

        # Paths
        cnr_csv_path = "datasets/webgraph/cnr-2000/cnr-2000.csv"
        output_dir = "test_data"

        if !isfile(cnr_csv_path)
            @warn "CNR-2000 CSV not found at: $cnr_csv_path"
            @test_skip "Dataset unavailable"
            return
        end

        # Load graph
        @info "Loading CNR-2000 graph..."
        @time g = load_adjacency_list_from_csv(cnr_csv_path, ',', true)
        n = nv(g)
        m = ne(g)
        @info "Graph: n=$n, m=$m"

        # Helpers
        # Helpers to build robust, total mappings from vertex orderings
        function order_to_mapping(order::Vector{S}, ::Type{V}) where {S<:Integer, V<:Unsigned}
            mapping = Dict{V,V}()
            @assert length(order) == length(unique(order)) "Order must contain each vertex exactly once"
            for (i, v) in enumerate(order)
                mapping[convert(V, v)] = convert(V, i)
            end
            return mapping
        end

        # Strategy 1: Lexicographic (built-in)
        @info "Computing lexicographic mapping..."
        map_lex = relabel_vertices(g, :lexicographic)

        # Strategy 2: Label Propagation (LLP, simple)
        @info "Computing LLP mapping..."
        function lp_order(g; iters::Int=10, seed::Int=42, V::Type{<:Unsigned}=eltype(vertices(g)))
            Random.seed!(seed)
            T = eltype(vertices(g))
            labels = Dict{T,T}()
            vs = collect(vertices(g))
            for v in vs
                labels[v] = v
            end
            for _ in 1:iters
                shuffle!(vs)
                for v in vs
                    # Use neighbors as undirected approximation: outneighbors only
                    neigh = outneighbors(g,v)
                    isempty(neigh) && continue
                    counts = Dict{T,Int}()
                    for u in neigh
                        l = labels[u]
                        counts[l] = get(counts,l,0) + 1
                    end
                    # pick most frequent label, tie-break by smallest label id
                    best = first(sort(collect(counts), by = x -> (-x[2], x[1])))
                    labels[v] = best[1]
                end
            end
            # Build a total order: by (label, vertex id)
            pairs = sort(collect(labels), by = x -> (x[2], x[1]))
            order = Vector{V}(undef, length(pairs))
            for (i, x) in enumerate(pairs)
                order[i] = convert(V, x[1])
            end
            return order
        end
        # V will be defined below before use

        # Strategy 3: MinHash on adjacency sets
        @info "Computing MinHash mapping..."
        function minhash_order(g; k::Int=4, seeds = (0x9e3779b97f4a7c15, 0xc2b2ae3d27d4eb4f, 0x165667b19e3779f9, 0xd6e8feb86659fd93), V::Type{<:Unsigned}=eltype(vertices(g)))
            T = eltype(vertices(g))
            vs = collect(vertices(g))
            sigs = Vector{Tuple{Vararg{UInt64}}}(undef, length(vs))
            for (i,v) in enumerate(vs)
                neigh = outneighbors(g, v)
                mins = fill(typemax(UInt64), length(seeds))
                for u in neigh
                    for (j,sd) in enumerate(seeds)
                        h = hash(u, sd)
                        if h < mins[j]
                            mins[j] = h
                        end
                    end
                end
                sigs[i] = tuple(mins...)
            end
            ordered = sort(collect(zip(vs, sigs)), by = x -> x[2])
            order = Vector{V}(undef, length(ordered))
            for (i, v_sig) in enumerate(ordered)
                order[i] = convert(V, v_sig[1])
            end
            return order
        end
        # V will be defined below before use

        # Strategy 4: BFS/RCM (approximate RCM with undirected view)
        @info "Computing RCM mapping (library :out variant)..."
        # V will be defined below before use

        mkpath(output_dir)

        # Compare compressions
        function compress_bits(path_prefix::AbstractString, g2; use_mix_mode=true, use_reference=true)
            write_compressed_mgs3_graph(g2, path_prefix, :children, :zeta, use_mix_mode, use_reference)
            file = path_prefix * ".mgz"
            sz = filesize(file)
            bpv = (sz * 8) / nv(g2)
            bpe = (sz * 8) / ne(g2)
            return (file=file, bits_per_vertex=bpv, bits_per_edge=bpe, size_bytes=sz)
        end

        # Build mapping lazily per strategy so failures are caught individually
        V = typeof(g).parameters[1]
        order_llp = lp_order(g; V=V)
        order_minhash = minhash_order(g; V=V)
        map_rcm = relabel_vertices_rcm(g, :out)
        configs = [
            ("lexicographic", () -> begin
                map = relabel_vertices(g, :lexicographic) # already Dict{V,V}
                # ensure type matches exactly Dict{V,V}
                Dict{V,V}(map)
            end),
            ("llp",           () -> order_to_mapping(order_llp, V)),
            ("minhash",       () -> order_to_mapping(order_minhash, V)),
            ("rcm",           () -> Dict{V,V}(map_rcm))
        ]

        results = []
        for (name, mkmap) in configs
            try
                @info "\n=== Order: $name ==="
                mapping = mkmap()
                # sanity check
                if length(mapping) != n
                    @warn "  Mapping size $(length(mapping)) != n=$n; attempting to continue"
                end
                @info "  Mapping type: $(typeof(mapping)) key=$(eltype(keys(mapping))) val=$(eltype(values(mapping)))"
                g_ord = relabel_graph(g, mapping)
                out = joinpath(output_dir, "cnr2000_order_" * name)
                @time stats = compress_bits(out, g_ord)
                @info "  File: $(basename(stats.file))"
                @info "  Bits/vertex: $(round(stats.bits_per_vertex, digits=3))"
                @info "  Bits/edge:   $(round(stats.bits_per_edge, digits=3))"
                push!(results, (name=name, bpv=stats.bits_per_vertex, bpe=stats.bits_per_edge, sz=stats.size_bytes))
            catch e
                @error "Ordering '$name' failed: $e"
                bt = stacktrace(catch_backtrace())
                for (i,frm) in enumerate(bt[1:min(end,10)])
                    @error "  [$i] $frm"
                end
            end
        end

        # Simple summary assertion: results exist and are finite
        @test length(results) >= 1
        @test all(isfinite(r.bpe) && isfinite(r.bpv) for r in results)

        # Log best when available
        if !isempty(results)
            best = findmin([(r.bpe, i) for (i,r) in enumerate(results)])[2]
            @info "\nBest ordering by bits/edge: $(results[best].name) with $(round(results[best].bpe,digits=3)) bits/edge"
        else
            @warn "No ordering produced results"
        end
    finally
        global_logger(prev_logger)
    end
end
