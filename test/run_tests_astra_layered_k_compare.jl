#!/usr/bin/env julia

# Compare ASTRA-L (new spec) compression for k=2 vs k=3

include("run_tests_main.jl")

using Logging
using LightGraphs: nv, ne
using Adjacently
using Adjacently.IO: load_adjacency_list_from_csv
using Adjacently.Graph: get_reverse_graph
using Adjacently.Compression.ASTRALayered: create_level_decomposition, write_astra_layered_graph

global_logger(ConsoleLogger(stderr, Logging.Info))

println("=" ^ 80)
println("ASTRA-L (New Spec) k-comparison Test")
println("Dataset: CNR-2000 web graph")
println("=" ^ 80)

# Load graph
dataset_csv = joinpath(PROJECT_ROOT, "datasets/webgraph/cnr-2000/cnr-2000.csv")
@info "Loading dataset" dataset=dataset_csv
g = load_adjacency_list_from_csv(dataset_csv, ',', true)
rg = get_reverse_graph(g)
V = nv(g)
E = ne(g)
println("Graph: $V vertices, $E edges")

mkpath(TEST_DIR)

function compress_k(radius::Int)
    @info "Building decomposition" k=radius
    decomp = create_level_decomposition(g, rg; radius=radius, log_every=2000)
    outfile = joinpath(TEST_DIR, "cnr2000_astral_k$(radius).astral")
    @info "Writing ASTRA-L file" file=outfile
    write_astra_layered_graph(decomp, outfile, g; integer_encoding=:fibonacci, radius=radius, log_every=0)
    sz = filesize(outfile)
    bpe = (8 * sz) / max(1, E)
    return outfile, sz, bpe
end

file2, size2, bpe2 = compress_k(2)
file3, size3, bpe3 = compress_k(3)

println("\nResults:")
println("  k=2: $(round(size2/1024/1024, digits=3)) MB, $(round(bpe2, digits=3)) bits/edge")
println("  k=3: $(round(size3/1024/1024, digits=3)) MB, $(round(bpe3, digits=3)) bits/edge")
println("  Δ (k=3 - k=2): $(round(bpe3 - bpe2, digits=3)) bits/edge")

@test size2 > 0 && size3 > 0
@test bpe2 > 0 && bpe3 > 0

println("\nDone.")
