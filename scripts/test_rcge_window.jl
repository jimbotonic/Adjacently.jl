#!/usr/bin/env julia

#
# Compare reference window sizes 64, 128, 256 on the best RCGE config
# (implicit ranges + adaptive adds+raw, MIL=2).
#
# Per-section cost model:
#   headers:  flag(1 bit) + delta(ceil(log2(window)) bits) per ref vertex
#   copy:     copy-blocks (unchanged by window)
#   add:      adaptive MIL=2 (unchanged by window)
#   raw:      adaptive MIL=2 (raw vertices may convert to ref with larger window)
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
        :intra_copy_blocks=>false, :intra_ref_fixwidth=>false,
        :intra_ref_vlc=>false, :intra_add_adaptive=>false, :intra_raw_adaptive=>false,
        :intra_adapt_mil=>2,
    )
    for (k, v) in kwargs; d[k] = v; end
    return RCGEParams(; d...)
end

# Best config: adaptive adds+raw, fixwidth, zigzag, stop, copy-blocks
best_cfg(w; kwargs...) = make_params(
    intra_zigzag=true, intra_stop_deltas=true, intra_copy_blocks=true,
    intra_ref_fixwidth=true, intra_add_adaptive=true, intra_raw_adaptive=true,
    intra_adapt_mil=2, intra_ref_window=w; kwargs...)

println("=" ^ 70)
println("RCGE Reference Window Comparison (w=64/128/256)")
println("Best config: implicit ranges + adaptive adds+raw MIL=2")
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
    ref_pct = round(100.0 * stats.intra_ref_used / (stats.intra_ref_used + stats.intra_no_ref), digits=1)
    println("$bpe BPE  ($(nbytes) bytes, $(dt)s)")
    println("    hdrs=$(ceil(Int,stats.bits_intra_headers/8))B  copy=$(ceil(Int,stats.bits_intra_copy/8))B  add=$(ceil(Int,stats.bits_intra_add/8))B  raw=$(ceil(Int,stats.bits_intra_raw/8))B  ref%=$(ref_pct)%")

    r_v = BitReader(IOBuffer(bytes))
    decoded = decode_level(r_v, params; T=TV, directed=true)
    dec_edges = sum(length(v) for v in values(decoded))
    ok = dec_edges == m
    println("    Roundtrip: $(ok ? "OK" : "MISMATCH (expected $m, got $dec_edges) !")")
    return nbytes, bpe, stats
end

println("\n" * "─" ^ 70)
println("Running window configs...")
println("─" ^ 70)

results = []

# No-adaptive baseline for each window (to isolate window effect on refs)
for w in [64, 128, 256]
    nbits_delta = max(1, ceil(Int, log2(w)))
    b, bpe, s = run_config("w=$w (no adaptive, $(nbits_delta)b delta)", make_params(
        intra_zigzag=true, intra_stop_deltas=true, intra_copy_blocks=true,
        intra_ref_fixwidth=true, intra_ref_window=w))
    push!(results, ("w=$w no-adaptive", b, bpe))
end

println()

# Best config (adaptive) for each window
for w in [64, 128, 256]
    nbits_delta = max(1, ceil(Int, log2(w)))
    b, bpe, s = run_config("w=$w (adaptive, $(nbits_delta)b delta)", best_cfg(w))
    push!(results, ("w=$w adaptive", b, bpe))
end

println("\n" * "=" ^ 70)
println("SUMMARY")
println("=" ^ 70)
ref_bpe = results[1][3]  # w=64 no-adaptive as reference
for (name, nbytes, bpe) in results
    delta = round(bpe - ref_bpe, digits=4)
    sign = delta <= 0 ? "" : "+"
    println("  $(rpad(name, 24)) $(round(bpe, digits=4)) BPE  $(sign)$(delta) vs w=64-no-adaptive")
end
println()
println("  Reference: WebGraph BV       2.897 BPE")
println("  Known best (w=64 adaptive)   2.5501 BPE")
best_bpe = minimum(r[3] for r in results)
println("  New best:                    $(round(best_bpe, digits=4)) BPE")
println("=" ^ 70)
