#!/usr/bin/env julia

#
# Greedy encoder: 3-way adaptive copy benchmark
# Tests complement copy-blocks (bitmap/copy-blocks/complement) on CNR-2000
# using the best relabeling (Leiden K=1 + per-group LLP) + w=64 + Fibonacci
#
# Baseline (config 8 from greedy_benchmark_cnr2000.jl):
#   Leiden+LLP + local-window + Fib w=64 (copy-blocks)
# Technique from RCGE:
#   3-way adaptive copy (bitmap/copy-blocks/complement) — same nested-bit scheme
#

using Pkg
Pkg.activate(normpath(joinpath(@__DIR__, "..")))

using Logging
global_logger(ConsoleLogger(stderr, Logging.Info))

using LightGraphs: nv, ne, outneighbors, vertices, edges, src, dst, is_directed, eltype
using Adjacently
using Adjacently.IO: load_adjacency_list_from_csv, BitWriter, BitReader, flush_bitwriter,
                      write_bytes
using Adjacently.Relabeling: relabel_graph, relabel_vertices_llp
using Adjacently.Clustering: leiden_partition
using Adjacently.Graph: subgraph
using Adjacently.MGS: load_compressed_mgs3_graph, create_header_flags,
                       OPTION_GREEDY_BASE, MGS_MAX_SIZE
using Adjacently.Compression: write_greedy_graph_data, read_greedy_graph_data
using Adjacently.Util: infer_uint_custom_type
using Adjacently.IO: BitWriter, BitReader, flush_bitwriter

const PROJECT_ROOT = normpath(joinpath(@__DIR__, ".."))
const CNR_CSV = joinpath(PROJECT_ROOT, "datasets", "webgraph", "cnr-2000", "cnr-2000.csv")
const OUTPUT_DIR = joinpath(PROJECT_ROOT, "datasets", "webgraph", "cnr-2000")

# Two-step relabeling: Leiden K=1 + per-group LLP (mirrors RCGE approach)
function relabel_leiden_k1_llp(g; leiden_max_passes::Int=8, leiden_max_levels::Int=5,
                                   llp_passes::Int=5)
    T = eltype(g)
    n = nv(g)

    print("    Leiden partitioning...")
    t1 = time()
    part = leiden_partition(g; max_passes=leiden_max_passes, max_levels=leiden_max_levels)
    nclusters = maximum(part)
    println(" → $nclusters clusters ($(round(time()-t1, digits=2))s)")

    counts = Dict{Int,Int}()
    for c in part; counts[c] = get(counts, c, 0) + 1; end
    top_label = argmax(counts)

    groups = [T[], T[]]
    for i in 1:n
        push!(groups[part[i] == top_label ? 1 : 2], T(i))
    end
    println("    Group sizes: $(length(groups[1])), $(length(groups[2]))")

    print("    Per-group LLP on induced subgraphs...")
    t2 = time()
    for gi in 1:2
        C = groups[gi]
        length(C) <= 2 && continue
        sg, oni, _ = subgraph(g, C)
        mapping = relabel_vertices_llp(sg, :sym; passes=llp_passes)
        sort!(C, by = v -> Int(mapping[oni[v]]))
        groups[gi] = C
    end
    println(" ($(round(time()-t2, digits=2))s)")

    vertex_mapping = Dict{T,T}()
    new_id = T(1)
    for C in groups
        for v in C
            vertex_mapping[v] = new_id
            new_id += T(1)
        end
    end
    group_sizes = [length(C) for C in groups]
    return vertex_mapping, group_sizes
end

function run_config(name, neighbor_lists::Dict{T,Vector{T}}, m::Int;
        ref_window_size::Int=64, copy_blocks::Bool=true, adaptive_copy::Bool=false,
        cluster_sizes::Vector{Int}=Int[]) where {T<:Unsigned}
    print("  $name ... ")
    io = IOBuffer()
    bw = BitWriter(io)
    t = time()
    write_greedy_graph_data(bw, neighbor_lists, :children, ref_window_size;
        integer_encoding=:fibonacci, stats=nothing,
        copy_blocks=copy_blocks, adaptive_copy=adaptive_copy,
        cluster_sizes=cluster_sizes)
    flush_bitwriter(bw; flush_last_bits=true)
    bytes = take!(io)
    dt = round(time() - t, digits=2)

    nbytes = length(bytes)
    bpe = round(8.0 * nbytes / m, digits=4)
    println("$bpe BPE  ($nbytes bytes, $(dt)s)")

    # Roundtrip
    r = BitReader(IOBuffer(bytes))
    vs = T(length(neighbor_lists))
    decoded = read_greedy_graph_data(r, vs, :children, T;
        integer_encoding=:fibonacci, copy_blocks=copy_blocks, adaptive_copy=adaptive_copy,
        ref_window_size=ref_window_size, cluster_sizes=cluster_sizes)
    dec_edges = sum(length(v) for v in values(decoded))
    ok = dec_edges == m
    println("    Roundtrip: $(ok ? "OK" : "MISMATCH (expected $m, got $dec_edges)!")")

    return nbytes, bpe, ok
end

println("=" ^ 70)
println("Greedy Encoder: 3-Way Adaptive Copy Benchmark")
println("Baseline: best from greedy_benchmark_cnr2000.jl (Leiden+LLP+local+Fib+w64+CB)")
println("=" ^ 70)

isfile(CNR_CSV) || error("CNR-2000 CSV not found at $CNR_CSV")

println("\nLoading CNR-2000...")
t0 = time()
g = load_adjacency_list_from_csv(CNR_CSV, ',', true)
n, m = nv(g), ne(g)
println("  Loaded: $n vertices, $m edges ($(round(time()-t0, digits=2))s)")

println("\nBuilding Leiden+LLP relabeling...")
t_rel = time()
leiden_llp_map, leiden_group_sizes = relabel_leiden_k1_llp(g)
g_rel = relabel_graph(g, leiden_llp_map)
println("  Total: $(round(time()-t_rel, digits=2))s, groups: $(leiden_group_sizes)")

T = eltype(g_rel)
neighbor_lists = Dict{T,Vector{T}}()
for v in vertices(g_rel)
    neighbor_lists[T(v)] = sort([T(o) for o in outneighbors(g_rel, v)])
end

println("\n" * "─" ^ 70)
println("Running configs...")
println("─" ^ 70)

results = []

# Baseline: copy-blocks only (no adaptive), global window w=64
b, bpe, ok = run_config("Baseline: CB w=64 (global)", neighbor_lists, m;
    ref_window_size=64, copy_blocks=true, adaptive_copy=false, cluster_sizes=Int[])
push!(results, ("CB w=64 (global)", b, bpe, ok))

# Baseline: copy-blocks only (no adaptive), local window w=64
b, bpe, ok = run_config("Baseline: CB w=64 (local window)", neighbor_lists, m;
    ref_window_size=64, copy_blocks=true, adaptive_copy=false, cluster_sizes=leiden_group_sizes)
push!(results, ("CB w=64 (local)", b, bpe, ok))

# 3-way adaptive: global window w=64
b, bpe, ok = run_config("3-way adaptive (bitmap/cb/complement) w=64 (global)", neighbor_lists, m;
    ref_window_size=64, copy_blocks=true, adaptive_copy=true, cluster_sizes=Int[])
push!(results, ("Adaptive w=64 (global)", b, bpe, ok))

# 3-way adaptive: local window w=64
b, bpe, ok = run_config("3-way adaptive (bitmap/cb/complement) w=64 (local)", neighbor_lists, m;
    ref_window_size=64, copy_blocks=true, adaptive_copy=true, cluster_sizes=leiden_group_sizes)
push!(results, ("Adaptive w=64 (local)", b, bpe, ok))

println("\n" * "=" ^ 70)
println("SUMMARY")
println("=" ^ 70)
ref_bpe = results[2][3]  # local CB baseline
for (name, nbytes, bpe, ok) in results
    delta = round(bpe - ref_bpe, digits=4)
    sign = delta <= 0 ? "" : "+"
    rt = ok ? "OK" : "FAIL"
    println("  $(rpad(name, 42)) $(round(bpe, digits=4)) BPE  $(sign)$(delta) vs local-CB  RT=$rt")
end
println()
println("  Reference: WebGraph BV                         2.897 BPE")
println("  RCGE best (3-way adaptive, Leiden K=1):        2.4341 BPE")
best_bpe = minimum(r[3] for r in results)
println("  Best this run:                                 $(round(best_bpe, digits=4)) BPE")
println("=" ^ 70)
