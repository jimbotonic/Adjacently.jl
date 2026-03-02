#!/usr/bin/env julia

#
# Test adaptive copy position encoding (bitmap vs copy-blocks per ref vertex).
# The mode bit (1 bit per ref vertex) selects the cheaper encoding:
#   - copy-blocks: good for contiguous run patterns
#   - bitmap: exactly ref_len bits, good for short reference lists
#

using Pkg
Pkg.activate(normpath(joinpath(@__DIR__, "..")))

using Logging
global_logger(ConsoleLogger(stderr, Logging.Info))

using LightGraphs
using LightGraphs: nv, outneighbors
using Adjacently
using Adjacently.RCGE: encode_level, RCGEParams, RCGEStats, decode_level
using Adjacently.IO: load_adjacency_list_from_csv, BitWriter, BitReader, flush_bitwriter
using Adjacently.Clustering: leiden_partition
using Adjacently.Relabeling: relabel_graph

count_edges(h) = sum(length(outneighbors(h, v)) for v in 1:nv(h))

function make_params(; kwargs...)
    d = Dict{Symbol,Any}(
        :L=>128, :varint=>:fibonacci, :count_varint=>:fibonacci, :gap=>:fibonacci,
        :degree=>:elias_delta, :undirected_pairs=>false, :perm_strategy=>:blockpos,
        :membership=>:implicit_ranges, :inter_strategy=>:perm, :intra_ref_enabled=>true,
        :intra_ref_window=>32, :intra_ref_rle=>false,
        :intra_block_try=>false, :positions_mode=>:delta, :additions_mode=>:delta,
        :min_cluster_density=>0.0, :intra_intervals=>false, :intra_mil=>4,
        :intra_greedy_mil=>false, :intra_zigzag=>false, :intra_stop_deltas=>false,
        :intra_copy_blocks=>false, :intra_copy_adaptive=>false, :intra_ref_fixwidth=>false,
        :intra_ref_vlc=>false, :intra_add_adaptive=>false, :intra_raw_adaptive=>false,
        :intra_adapt_mil=>2,
    )
    for (k, v) in kwargs; d[k] = v; end
    return RCGEParams(; d...)
end

# Known best config (2.5501 BPE): adaptive adds+raw, w=64
best64(; kwargs...) = make_params(
    intra_zigzag=true, intra_stop_deltas=true, intra_copy_blocks=true,
    intra_ref_fixwidth=true, intra_add_adaptive=true, intra_raw_adaptive=true,
    intra_adapt_mil=2, intra_ref_window=64; kwargs...)

println("=" ^ 70)
println("RCGE Adaptive Copy Positions Benchmark")
println("Baseline: 2.5501 BPE (adaptive adds+raw, w=64)")
println("=" ^ 70)

cnr_csv = normpath(joinpath(@__DIR__, "..", "datasets", "webgraph", "cnr-2000", "cnr-2000.csv"))
isfile(cnr_csv) || error("CNR-2000 CSV not found at $cnr_csv")

println("\nLoading CNR-2000...")
t0 = time()
g  = load_adjacency_list_from_csv(cnr_csv, ',', true)
n  = nv(g)
m  = count_edges(g)
TV = eltype(g)
println("  Loaded: $n vertices, $m edges ($(round(time()-t0, digits=2))s)")

println("\nBuilding pre-relabeled graph (Leiden K=1, vertex-ID rank relabeling)...")
t_rel = time()
part = leiden_partition(g; max_passes=8, max_levels=5)
counts = Dict{Int,Int}()
for c in part; counts[c] = get(counts, c, 0) + 1; end
top_label = argmax(counts)

clusters_raw = [TV[], TV[]]
for i in 1:n
    push!(clusters_raw[part[i] == top_label ? 1 : 2], TV(i))
end
clusters_by_id = [sort(C) for C in clusters_raw]
S1 = length(clusters_by_id[1])

vertex_map = let new_id = TV(1)
    d = Dict{TV,TV}()
    for C in clusters_by_id
        for v in C
            d[v] = new_id
            new_id += TV(1)
        end
    end
    d
end
g_rel = relabel_graph(g, vertex_map)
clusters_impl = [TV.(1:S1), TV.(S1+1:n)]
println("  Clusters: $S1 + $(n-S1) ($(round(time()-t_rel, digits=2))s)")

function run_config(name, params)
    print("  $name ... ")
    io = IOBuffer(); w = BitWriter(io)
    stats = RCGEStats()
    t = time()
    encode_level(w, g_rel, clusters_impl; params=params, stats=stats)
    flush_bitwriter(w; flush_last_bits=true)
    bytes = take!(io)
    dt = round(time() - t, digits=2)

    nbytes = length(bytes)
    bpe = round(8.0 * nbytes / m, digits=4)
    println("$bpe BPE  ($(nbytes) bytes, $(dt)s)")
    println("    hdrs=$(ceil(Int,stats.bits_intra_headers/8))B  copy=$(ceil(Int,stats.bits_intra_copy/8))B  add=$(ceil(Int,stats.bits_intra_add/8))B  raw=$(ceil(Int,stats.bits_intra_raw/8))B")

    r_v = BitReader(IOBuffer(bytes))
    decoded = decode_level(r_v, params; T=TV, directed=true)
    dec_edges = sum(length(v) for v in values(decoded))
    ok = dec_edges == m
    println("    Roundtrip: $(ok ? "OK" : "MISMATCH (expected $m, got $dec_edges) !")")
    return nbytes, bpe, stats
end

println("\n" * "─" ^ 70)
println("Running configs...")
println("─" ^ 70)

results = []

# Baseline (known)
b, bpe, s = run_config("Baseline (w=64, adaptive adds+raw)", best64())
push!(results, ("Baseline w=64", b, bpe))

# + adaptive copy (bitmap vs copy-blocks per ref vertex)
b, bpe, s = run_config("+ Adaptive copy (bitmap vs copy-blocks)", best64(intra_copy_adaptive=true))
push!(results, ("+ Adaptive copy", b, bpe))

println("\n" * "=" ^ 70)
println("SUMMARY")
println("=" ^ 70)
baseline_bpe = results[1][3]
for (name, nbytes, bpe) in results
    delta = round(bpe - baseline_bpe, digits=4)
    sign = delta <= 0 ? "" : "+"
    println("  $(rpad(name, 32)) $(round(bpe, digits=4)) BPE  $(sign)$(delta) vs baseline")
end
println()
println("  Reference: WebGraph BV         2.897 BPE")
println("  Known best (w=64 adaptive)     2.5501 BPE")
best_bpe = minimum(r[3] for r in results)
println("  Best this run:                 $(round(best_bpe, digits=4)) BPE")
println("=" ^ 70)
