#!/usr/bin/env julia
#
# Test auto_select_K on CNR-2000 and IN-2004
#

using Pkg
Pkg.activate(normpath(joinpath(@__DIR__, "..")))

using Logging
global_logger(ConsoleLogger(stderr, Logging.Info))

using LightGraphs: nv, ne
using Adjacently
using Adjacently.IO: load_adjacency_list_from_csv
using Adjacently.Clustering: auto_select_K, leiden_partition

datasets = [
    ("cnr-2000", 2),
    ("in-2004", 4),
]

for (ds, expected_K) in datasets
    csv = normpath(joinpath(@__DIR__, "..", "datasets", "webgraph", ds, "$ds.csv"))
    if !isfile(csv)
        @warn "Skipping $ds — CSV not found at $csv"
        continue
    end

    println("\n" * "=" ^ 60)
    println("Dataset: $ds (expected K=$expected_K)")
    println("=" ^ 60)

    @info "Loading $ds..."
    t0 = time()
    g = load_adjacency_list_from_csv(csv, ',', true)
    @info "  Loaded: $(nv(g)) vertices, $(ne(g)) edges ($(round(time()-t0, digits=1))s)"

    @info "Running auto_select_K..."
    t1 = time()
    K, part = auto_select_K(g)
    dt = round(time() - t1, digits=1)

    println("  → K = $K  (expected $expected_K)  [$(dt)s]")
    if K == expected_K
        println("  ✓ PASS")
    else
        println("  ✗ FAIL (got $K, expected $expected_K)")
    end
end

println("\nDone.")
