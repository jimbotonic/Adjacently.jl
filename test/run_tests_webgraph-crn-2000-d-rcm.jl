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
using Adjacently.Relabeling: relabel_graph, relabel_vertices_rcm
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

        # Use library functions for RCM mappings

        # Helper to run a variant and report
        function run_variant(order::Vector{V}, tag::AbstractString)
            @info "Order($tag) length: $(length(order)) (expected $n)"
            @test length(order) == n
            mapping = Dict{V,V}()
            for (i,v) in enumerate(order)
                mapping[v] = convert(V, i)
            end
            g_ord = relabel_graph(g, mapping)
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
        map_out = relabel_vertices_rcm(g, :out)
        orderA = [v for (v,_) in sort(collect(map_out), by=x->x[2])]
        resA = run_variant(orderA, "outdeg")

        @info "Computing RCM symmetrized..."
        map_sym = relabel_vertices_rcm(g, :sym)
        orderB = [v for (v,_) in sort(collect(map_sym), by=x->x[2])]
        resB = run_variant(orderB, "sym")

        # Summary
        @info "Summary:"
        @info "  Outdeg-only: bpv=$(round(resA.bpv; digits=3)) bpe=$(round(resA.bpe; digits=3)) file=$(basename(resA.file))"
        @info "  Symmetrized: bpv=$(round(resB.bpv; digits=3)) bpe=$(round(resB.bpe; digits=3)) file=$(basename(resB.file))"
    finally
        global_logger(prev_logger)
    end
end
