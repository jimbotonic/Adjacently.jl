#!/usr/bin/env julia

#
# CNR-2000: Improved RCM ordering + compression (zeta + mix + reference)
# - Symmetrized neighbors (in ∪ out)
# - Pseudo-peripheral start (Gibbs–Poole–Stockmeyer style, few sweeps)
# - Degree-aware BFS with increasing-degree neighbor expansion
# - Reverse order (RCM)
#

include("run_tests_main.jl")

using Logging, Random
using LightGraphs: inneighbors, outneighbors, nv, ne
using Adjacently.IO: load_adjacency_list_from_csv
using Adjacently.Graph: get_in_out_degrees
using Adjacently.MGS: write_compressed_mgs3_graph

@testset "CNR-2000 Improved RCM Ordering" begin
    prev_logger = current_logger()
    global_logger(ConsoleLogger(stderr, Logging.Info))
    try
        @info "=== CNR-2000 Improved RCM Ordering ==="

        # Paths
        cnr_csv_path = "datasets/webgraph/cnr-2000/cnr-2000.csv"
        output_dir = "test_data"
        mkpath(output_dir)

        if !isfile(cnr_csv_path)
            @warn "CNR-2000 CSV not found at: $cnr_csv_path"
            @test_skip "Dataset unavailable"
            return
        end

        # Load graph (typed to custom UInt), use header
        @info "Loading graph..."
        @time g = load_adjacency_list_from_csv(cnr_csv_path, ',', true)
        n = nv(g); m = ne(g)
        V = typeof(g).parameters[1]
        @info "Graph: n=$n, m=$m, V=$(V)"

        # Degrees (symmetrized): indeg+outdeg
        @info "Computing out-degrees (directed) for RCM..."
        # Use out-degree only, as in the better-performing RCM from the comparison test
        _, outdeg = get_in_out_degrees(g)
        deg = Dict{V,Int}()
        for (v, d) in outdeg
            deg[v] = Int(d)
        end

        # Neighbor accessors
        outnbrs(v::V) = outneighbors(g, v)
        symnbrs(v::V) = union(outneighbors(g, v), inneighbors(g, v))

        # Compute BFS levels from a start using symmetrized neighbors
        function bfs_levels(start::V, nbrs::Function)
            visited = falses(Int(n) + 1)
            visited[Int(start)] = true
            levels = Vector{Vector{V}}()
            frontier = V[start]             # 1-element Vector{V}
            push!(levels, copy(frontier))
            while true
                next_frontier = Vector{V}()
                for v in frontier
                    for u in nbrs(v)
                        iu = Int(u)
                        if !visited[iu]
                            visited[iu] = true
                            push!(next_frontier, u)
                        end
                    end
                end
                if isempty(next_frontier)
                    break
                end
                push!(levels, next_frontier)
                frontier = next_frontier
            end
            return levels
        end

        # Variant A: RCM using out-degree and out-neighbors only
        function rcm_order_outdeg()
            vs = collect(vertices(g))
            visited = falses(Int(n) + 1)
            order = V[]
            # global min-outdegree start
            start = vs[argmin(x -> get(deg, x, 0), vs)]
            queue = V[start]
            visited[Int(start)] = true
            push!(order, start)
            qidx = 1
            while qidx <= length(queue)
                v = queue[qidx]
                qidx += 1
                neigh = V[ u for u in outnbrs(v) if !visited[Int(u)] ]
                sort!(neigh, by = x -> get(deg, x, 0))
                for u in neigh
                    iu = Int(u)
                    if !visited[iu]
                        visited[iu] = true
                        push!(queue, u)
                        push!(order, u)
                    end
                end
            end
            # add unvisited in original order
            for v in vs
                if !visited[Int(v)]
                    push!(order, v)
                end
            end
            reverse!(order)
            return order
        end

        # Variant B: RCM using symmetrized neighbors and (in+out) degree, with 2 BFS sweeps
        function rcm_order_sym()
            indeg2, outdeg2 = get_in_out_degrees(g)
            deg2 = Dict{V,Int}()
            for (v,d) in outdeg2; deg2[v] = Int(d) end
            for (v,d) in indeg2; deg2[v] = get(deg2,v,0) + Int(d) end

            visited = falses(Int(n) + 1)
            order = V[]
            remaining = n
            while remaining > 0
                # pick unvisited min-degree as seed
                seed = zero(V); mind = typemax(Int)
                for v in keys(deg2)
                    if !visited[Int(v)] && deg2[v] < mind
                        mind = deg2[v]; seed = v
                    end
                end
                if seed == zero(V); break; end
                # Two sweeps to push towards periphery
                lvls1 = bfs_levels(seed, symnbrs)
                last1 = lvls1[end]
                cand = last1[argmin([(get(deg2,v,0),v) for v in last1])]
                lvls2 = bfs_levels(cand, symnbrs)
                last2 = lvls2[end]
                start = last2[argmin([(get(deg2,v,0),v) for v in last2])]

                # BFS with degree buckets
                queue = [start]
                visited[Int(start)] = true
                push!(order, start)
                remaining -= 1
                qidx = 1
                while qidx <= length(queue)
                    v = queue[qidx]
                    qidx += 1
                    buckets = Dict{Int,Vector{V}}()
                    for u in symnbrs(v)
                        iu = Int(u)
                        if !visited[iu]
                            du = get(deg2, u, 0)
                            push!(get!(buckets, du, V[]), u)
                        end
                    end
                    if !isempty(buckets)
                        ds = sort!(collect(keys(buckets)))
                        for d in ds
                            for u in buckets[d]
                                iu = Int(u)
                                if !visited[iu]
                                    visited[iu] = true
                                    push!(queue, u)
                                    push!(order, u)
                                    remaining -= 1
                                end
                            end
                        end
                    end
                end
            end
            reverse!(order)
            return order
        end

        # Helper to run a variant and report
        function run_variant(order::Vector{V}, tag::AbstractString)
            @info "Order($tag) length: $(length(order)) (expected $n)"
            @test length(order) == n
            mapping = Dict{V,V}()
            for (i,v) in enumerate(order)
                mapping[v] = convert(V, i)
            end
            g_ord = Adjacently.Graph.relabel_graph(g, mapping)
            out = joinpath(output_dir, "cnr2000_rcm_" * tag)
            @info "Compressing $tag (zeta + mix + reference + children)..."
            @time write_compressed_mgs3_graph(g_ord, out, :children, :zeta, true, true)
            file = out * ".mgz"
            sz = filesize(file)
            bpv = (sz * 8) / n
            bpe = (sz * 8) / m
            @info "Compression($tag): file=$(basename(file)), size=$(round(sz/1024/1024; digits=3)) MB, bpv=$(round(bpv; digits=3)), bpe=$(round(bpe; digits=3))"
            @test isfinite(bpv)
            @test isfinite(bpe)
            return (file=file, sz=sz, bpv=bpv, bpe=bpe)
        end

        # Run both approaches
        @info "Computing RCM outdeg-only..."
        orderA = rcm_order_outdeg()
        resA = run_variant(orderA, "outdeg")

        @info "Computing RCM symmetrized..."
        orderB = rcm_order_sym()
        resB = run_variant(orderB, "sym")

        # Summary
        @info "Summary:"
        @info "  Outdeg-only: bpv=$(round(resA.bpv; digits=3)) bpe=$(round(resA.bpe; digits=3)) file=$(basename(resA.file))"
        @info "  Symmetrized: bpv=$(round(resB.bpv; digits=3)) bpe=$(round(resB.bpe; digits=3)) file=$(basename(resB.file))"
    finally
        global_logger(prev_logger)
    end
end
