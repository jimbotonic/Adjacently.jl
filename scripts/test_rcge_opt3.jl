#!/usr/bin/env julia

#
# Benchmark three new RCGE optimizations on top of implicit-ranges baseline (2.666 BPE):
#   1. VLC reference distance (Fibonacci instead of fixed-width)
#   2. Interval encoding for additions (per-vertex adaptive)
#   3. Interval encoding for raw lists (per-vertex adaptive)
#
# All three are independent and can be combined.
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
        :intra_ref_window=>32, :intra_ref_min_overlap=>0.3, :intra_ref_rle=>false,
        :intra_block_try=>false, :positions_mode=>:delta, :additions_mode=>:delta,
        :min_cluster_density=>0.0, :intra_intervals=>false, :intra_mil=>4,
        :intra_greedy_mil=>false, :intra_zigzag=>false, :intra_stop_deltas=>false,
        :intra_copy_blocks=>false, :intra_ref_fixwidth=>false,
        :intra_ref_vlc=>false, :intra_add_adaptive=>false, :intra_raw_adaptive=>false,
        :intra_adapt_mil=>2,
    )
    for (k, v) in kwargs; d[k] = v; end
    return RCGEParams(; d...)
end

# Baseline: best known config with implicit ranges (2.666 BPE)
fw64_base(; kwargs...) = make_params(
    intra_zigzag=true, intra_stop_deltas=true, intra_copy_blocks=true,
    intra_ref_window=64, intra_ref_fixwidth=true; kwargs...)

println("=" ^ 70)
println("RCGE Optimization Benchmark (baseline: 2.666 BPE implicit ranges)")
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

# ── Build pre-relabeled graph (shared across all configs) ─────────────────────
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

# ── Helper: encode + verify ───────────────────────────────────────────────────
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
    println("    headers=$(ceil(Int,stats.bits_intra_headers/8))B  copy=$(ceil(Int,stats.bits_intra_copy/8))B  add=$(ceil(Int,stats.bits_intra_add/8))B  raw=$(ceil(Int,stats.bits_intra_raw/8))B  mem=$(ceil(Int,stats.bits_membership/8))B")

    # Roundtrip verify
    r_v = BitReader(IOBuffer(bytes))
    decoded = decode_level(r_v, params; T=TV, directed=true)
    dec_edges = sum(length(v) for v in values(decoded))
    if dec_edges == m
        println("    Roundtrip: OK")
    else
        println("    Roundtrip: MISMATCH (expected $m, got $dec_edges) !")
    end
    return nbytes, bpe, stats
end

println("\n" * "─" ^ 70)
println("Running configs...")
println("─" ^ 70)

results = []

# Baseline
b, bpe, s = run_config("Baseline (implicit, no new opts)", fw64_base())
push!(results, ("Baseline", b, bpe))

# 1. VLC ref distance
b, bpe, s = run_config("+ VLC ref distance", fw64_base(intra_ref_vlc=true))
push!(results, ("+ VLC ref delta", b, bpe))

# 2. Adaptive additions (MIL=2)
b, bpe, s = run_config("+ Adaptive additions (MIL=2)", fw64_base(intra_add_adaptive=true, intra_adapt_mil=2))
push!(results, ("+ Adaptive adds MIL=2", b, bpe))

# 3. Adaptive raw (MIL=2)
b, bpe, s = run_config("+ Adaptive raw (MIL=2)", fw64_base(intra_raw_adaptive=true, intra_adapt_mil=2))
push!(results, ("+ Adaptive raw MIL=2", b, bpe))

# Best: adaptive adds + adaptive raw WITHOUT VLC
b, bpe, s = run_config("+ Adaptive adds+raw (no VLC)", fw64_base(intra_add_adaptive=true, intra_raw_adaptive=true, intra_adapt_mil=2))
push!(results, ("Adaptive adds+raw MIL=2", b, bpe))

# Combined: VLC + adaptive adds + adaptive raw
b, bpe, s = run_config("+ All three combined", fw64_base(intra_ref_vlc=true, intra_add_adaptive=true, intra_raw_adaptive=true, intra_adapt_mil=2))
push!(results, ("All three combined", b, bpe))

println("\n" * "=" ^ 70)
println("SUMMARY")
println("=" ^ 70)
baseline_bpe = results[1][3]
for (name, nbytes, bpe) in results
    delta = round(bpe - baseline_bpe, digits=4)
    sign = delta <= 0 ? "" : "+"
    println("  $(rpad(name, 28)) $(round(bpe, digits=4)) BPE  $(sign)$(delta) vs baseline")
end
println("")
println("  Reference: WebGraph BV                2.897 BPE")
println("  Baseline: implicit ranges             $(baseline_bpe) BPE")
best_bpe = minimum(r[3] for r in results)
println("  New best:                             $(round(best_bpe, digits=4)) BPE")
println("=" ^ 70)
