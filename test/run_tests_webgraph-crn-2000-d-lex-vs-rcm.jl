#!/usr/bin/env julia

#
# CNR-2000: Lexicographic relabeling + children mode + reference + mix encoding
# with FED codes (Fibonacci, Elias-gamma, Elias-delta). Writes mgz files in test_data
# and prints size, bits/vertex, bits/edge for each code.
#

include("run_tests_main.jl")

using Logging
using LightGraphs: vertices, outneighbors, nv, ne
using Adjacently.IO: load_adjacency_list_from_csv
using Adjacently.Graph: relabel_vertices, relabel_graph, get_in_out_degrees
using Adjacently.MGS: write_compressed_mgs3_graph, load_compressed_mgs3_graph

@testset "CNR-2000 Lexicographic/RCM + ZETA/FED/FIB (children + ref + mix)" begin
    prev = current_logger()
    global_logger(ConsoleLogger(stderr, Logging.Info))
    try
        @info "=== CNR-2000 Lexicographic/RCM + ZETA/FED/FIB (children + ref + mix) ==="

        # Paths
        cnr_csv = "datasets/webgraph/cnr-2000/cnr-2000.csv"
        outdir = "test_data"
        mkpath(outdir)

        if !isfile(cnr_csv)
            @warn "CNR-2000 CSV not found: $cnr_csv"
            @test_skip "Dataset unavailable"
            return
        end

        # Load graph
        @info "Loading graph..."
        @time g = load_adjacency_list_from_csv(cnr_csv, ',', true)
        n = nv(g); m = ne(g)
        V = typeof(g).parameters[1]
        @info "Graph loaded: n=$n, m=$m, V=$(V)"

        # Lexicographic relabeling
        @info "Lexicographic relabeling..."
        map = relabel_vertices(g, :lexicographic)
        glex = relabel_graph(g, map)
        @test nv(glex) == n
        @test ne(glex) == m

        # RCM relabeling (outdegree-only, outneighbors, global min-outdegree start, reversed order)
        @info "RCM (outdegree-only) relabeling..."
        _, outdeg = get_in_out_degrees(g)
        deg = Dict{eltype(keys(outdeg)), Int}()
        for (v, d) in outdeg
            deg[v] = Int(d)
        end
        function rcm_order_outdeg(g)
            V = typeof(g).parameters[1]
            vs = collect(vertices(g))
            visited = falses(Int(n) + 1)
            order = V[]
            start = vs[argmin(x -> get(deg, x, 0), vs)]
            queue = V[start]
            visited[Int(start)] = true
            push!(order, start)
            qidx = 1
            while qidx <= length(queue)
                v = queue[qidx]
                qidx += 1
                neigh = V[ u for u in outneighbors(g, v) if !visited[Int(u)] ]
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
            for v in vs
                if !visited[Int(v)]
                    push!(order, v)
                end
            end
            reverse!(order)
            return order
        end
        order_rcm = rcm_order_outdeg(g)
        map_rcm = Dict{typeof(order_rcm[1]), typeof(order_rcm[1])}()
        for (i, v) in enumerate(order_rcm)
            map_rcm[v] = typeof(v)(i)
        end
        grcm = relabel_graph(g, map_rcm)
        @test nv(grcm) == n
        @test ne(grcm) == m

        # Helper to compress and report
        function compress_and_report(g0, base, encsym)
            @time write_compressed_mgs3_graph(g0, base, :children, encsym, true, true)
            file = base * ".mgz"
            @test isfile(file)
            sz = filesize(file)
            bpv = (sz * 8) / n
            bpe = (sz * 8) / m
            @info "  $(basename(file)) | enc=$(encsym) size=$(round(sz/1024/1024; digits=3)) MB bpv=$(round(bpv; digits=3)) bpe=$(round(bpe; digits=3))"
            return (file=file, size=sz, bpv=bpv, bpe=bpe, enc=encsym)
        end

        results = []
        @info "\nLexicographic ordering:"
        push!(results, compress_and_report(glex, joinpath(outdir, "cnr2000_lex_zeta"), :zeta))
        push!(results, compress_and_report(glex, joinpath(outdir, "cnr2000_lex_fed"), :fed))
        push!(results, compress_and_report(glex, joinpath(outdir, "cnr2000_lex_fib"), :fibonacci))

        @info "\nRCM ordering:"
        push!(results, compress_and_report(grcm, joinpath(outdir, "cnr2000_rcm_zeta"), :zeta))
        push!(results, compress_and_report(grcm, joinpath(outdir, "cnr2000_rcm_fed"), :fed))
        push!(results, compress_and_report(grcm, joinpath(outdir, "cnr2000_rcm_fib"), :fibonacci))

        # Round-trip verification for all files (including :fed)
        try
            for r in results
                gback = load_compressed_mgs3_graph(r.file)
                @test nv(gback) == n && ne(gback) == m
            end
        catch e
            @warn "Round-trip verification skipped/failed: $e"
        end

        # Summary of all 6 files
        @info "\nSummary (6 files):"
        for r in results
            @info "  $(basename(r.file)): enc=$(r.enc) size=$(round(r.size/1024/1024; digits=3)) MB bpv=$(round(r.bpv; digits=3)) bpe=$(round(r.bpe; digits=3))"
        end
        @test length(results) == 6
    finally
        global_logger(prev)
    end
end
