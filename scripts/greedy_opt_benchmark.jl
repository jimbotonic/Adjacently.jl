#!/usr/bin/env julia

# Focused benchmark: test new optimizations (compact_copy, tight_intervals, vlc2, window size)
# Previous best baseline: adaptive + stop_d + empty = ~2.93 BPE

using Pkg
Pkg.activate(normpath(joinpath(@__DIR__, "..")))

using Logging
global_logger(ConsoleLogger(stderr, Logging.Info))

using LightGraphs: nv, ne, outneighbors, vertices, eltype
using Adjacently
using Adjacently.IO: load_adjacency_list_from_csv, BitWriter, BitReader, flush_bitwriter
using Adjacently.Relabeling: relabel_graph, relabel_vertices_llp
using Adjacently.Clustering: leiden_partition
using Adjacently.Graph: subgraph
using Adjacently.Compression: write_greedy_graph_data, read_greedy_graph_data

const CNR_CSV = normpath(joinpath(@__DIR__, "..", "datasets", "webgraph", "cnr-2000", "cnr-2000.csv"))

function relabel_leiden_k1_llp(g)
    T = eltype(g); n = nv(g)
    part = leiden_partition(g; max_passes=8, max_levels=5)
    counts = Dict{Int,Int}()
    for c in part; counts[c] = get(counts, c, 0) + 1; end
    top_label = argmax(counts)
    groups = [T[], T[]]
    for i in 1:n; push!(groups[part[i] == top_label ? 1 : 2], T(i)); end
    for gi in 1:2
        C = groups[gi]; length(C) <= 2 && continue
        sg, oni, _ = subgraph(g, C)
        mapping = relabel_vertices_llp(sg, :sym; passes=5)
        sort!(C, by = v -> Int(mapping[oni[v]]))
        groups[gi] = C
    end
    vertex_mapping = Dict{T,T}(); new_id = T(1)
    for C in groups; for v in C; vertex_mapping[v] = new_id; new_id += T(1); end; end
    return vertex_mapping
end

function run_test(name, nls, m, T; window_size=64, compact_copy=false, tight_intervals=false, vlc2=false)
    print("  $name ... ")
    io = IOBuffer(); bw = BitWriter(io)
    t = time()
    write_greedy_graph_data(bw, nls, :children, window_size;
        integer_encoding=:fibonacci, copy_blocks=true, adaptive_copy=true,
        stop_deltas=true, empty_prefix=true,
        compact_copy=compact_copy, tight_intervals=tight_intervals, vlc2=vlc2)
    flush_bitwriter(bw; flush_last_bits=true)
    bytes = take!(io)
    dt = round(time() - t, digits=1)
    bpe = round(8.0 * length(bytes) / m, digits=4)
    println("$bpe BPE  ($(length(bytes)) B, $(dt)s)")
    # Roundtrip
    r = BitReader(IOBuffer(bytes))
    decoded = read_greedy_graph_data(r, T(length(nls)), :children, T;
        integer_encoding=:fibonacci, copy_blocks=true, adaptive_copy=true,
        ref_window_size=window_size, stop_deltas=true, empty_prefix=true,
        compact_copy=compact_copy, tight_intervals=tight_intervals, vlc2=vlc2)
    dec_edges = sum(length(v) for v in values(decoded))
    ok = dec_edges == m
    if !ok; println("    ROUNDTRIP MISMATCH: $dec_edges vs $m"); end
    return bpe, ok
end

println("=" ^ 70)
println("New Optimizations Benchmark — CNR-2000")
println("=" ^ 70)

isfile(CNR_CSV) || error("CNR-2000 CSV not found at $CNR_CSV")

println("\nLoading CNR-2000...")
t0 = time()
g = load_adjacency_list_from_csv(CNR_CSV, ',', true)
n, m = nv(g), ne(g)
println("  $n vertices, $m edges ($(round(time()-t0, digits=1))s)")

println("Leiden+LLP relabeling...")
t_rel = time()
lmap = relabel_leiden_k1_llp(g)
g_rel = relabel_graph(g, lmap)
println("  Done ($(round(time()-t_rel, digits=1))s)")

T = eltype(g_rel)
nls = Dict{T,Vector{T}}()
for v in vertices(g_rel); nls[T(v)] = sort([T(o) for o in outneighbors(g_rel, v)]); end

results = []
println("\n" * "─" ^ 70)

# Baseline: previous best (adaptive + stop_d + empty)
bpe, ok = run_test("Baseline (adapt+stop+empty)", nls, m, T)
push!(results, ("Baseline", bpe, ok))

println("\n" * "─" ^ 70)
println("INDIVIDUAL OPTIMIZATIONS:")
println("─" ^ 70)

# Opt 1: compact_copy
bpe, ok = run_test("+ compact_copy", nls, m, T; compact_copy=true)
push!(results, ("+ compact_copy", bpe, ok))

# Opt 2: tight_intervals
bpe, ok = run_test("+ tight_intervals", nls, m, T; tight_intervals=true)
push!(results, ("+ tight_intervals", bpe, ok))

# Opt 3: vlc2
bpe, ok = run_test("+ vlc2", nls, m, T; vlc2=true)
push!(results, ("+ vlc2", bpe, ok))

println("\n" * "─" ^ 70)
println("COMBINATIONS:")
println("─" ^ 70)

# compact + tight
bpe, ok = run_test("compact + tight", nls, m, T; compact_copy=true, tight_intervals=true)
push!(results, ("compact+tight", bpe, ok))

# compact + vlc2
bpe, ok = run_test("compact + vlc2", nls, m, T; compact_copy=true, vlc2=true)
push!(results, ("compact+vlc2", bpe, ok))

# ALL three
bpe, ok = run_test("ALL: compact + tight + vlc2", nls, m, T; compact_copy=true, tight_intervals=true, vlc2=true)
push!(results, ("ALL 3 opts", bpe, ok))

println("\n" * "─" ^ 70)
println("WINDOW SIZE TUNING (with ALL 3 opts):")
println("─" ^ 70)

for ws in [32, 64, 128, 256]
    bpe, ok = run_test("ALL w=$ws", nls, m, T; compact_copy=true, tight_intervals=true, vlc2=true, window_size=ws)
    push!(results, ("ALL w=$ws", bpe, ok))
end

println("\n" * "=" ^ 70)
println("SUMMARY")
println("=" ^ 70)
baseline_bpe = results[1][2]
for (name, bpe, ok) in results
    d = round(bpe - baseline_bpe, digits=4)
    sign = d <= 0 ? "" : "+"
    rt = ok ? "OK" : "FAIL"
    println("  $(rpad(name, 25)) $(rpad(string(bpe), 8)) BPE  $(sign)$(d)  RT=$rt")
end
println()
println("  WebGraph BV reference:   2.897  BPE")
println("  RCGE best:               2.4341 BPE")
println("=" ^ 70)
