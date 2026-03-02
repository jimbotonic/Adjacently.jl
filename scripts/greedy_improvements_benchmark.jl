#!/usr/bin/env julia

# Benchmark: greedy encoder optimizations on CNR-2000
# Base config: adaptive_copy + stop_deltas (3.0684 BPE known)
# Tests each new optimization individually, then best combos

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

function run_config(name, nls::Dict{T,Vector{T}}, m::Int;
        copy_blocks::Bool=true, adaptive_copy::Bool=true,
        stop_deltas::Bool=true,
        fixwidth_ref::Bool=false,
        adaptive_deltas::Bool=false,
        empty_prefix::Bool=false,
        split_residual::Bool=false,
        bv_blocks::Bool=false,
        compact_copy::Bool=false,
        tight_intervals::Bool=false,
        vlc2::Bool=false) where {T<:Unsigned}

    vs = length(nls)
    println("\n  [$name]")
    println("    Encoding $vs vertices...")
    io = IOBuffer(); bw = BitWriter(io)
    t = time()
    log_interval = max(1, vs ÷ 10)  # log every 10%

    write_greedy_graph_data(bw, nls, :children, 64;
        integer_encoding=:fibonacci, copy_blocks=copy_blocks,
        adaptive_copy=adaptive_copy, fixwidth_ref=fixwidth_ref,
        stop_deltas=stop_deltas, adaptive_deltas=adaptive_deltas,
        empty_prefix=empty_prefix, split_residual=split_residual,
        bv_blocks=bv_blocks, compact_copy=compact_copy,
        tight_intervals=tight_intervals, vlc2=vlc2)
    flush_bitwriter(bw; flush_last_bits=true)
    bytes = take!(io)
    dt = round(time() - t, digits=1)
    bpe = round(8.0 * length(bytes) / m, digits=4)
    println("    Encode done: $bpe BPE  ($(length(bytes)) B, $(dt)s)")

    # Roundtrip verification
    println("    Verifying roundtrip...")
    t2 = time()
    r = BitReader(IOBuffer(bytes))
    decoded = read_greedy_graph_data(r, T(vs), :children, T;
        integer_encoding=:fibonacci, copy_blocks=copy_blocks,
        adaptive_copy=adaptive_copy, ref_window_size=64,
        fixwidth_ref=fixwidth_ref, stop_deltas=stop_deltas,
        adaptive_deltas=adaptive_deltas, empty_prefix=empty_prefix,
        split_residual=split_residual, bv_blocks=bv_blocks,
        compact_copy=compact_copy, tight_intervals=tight_intervals,
        vlc2=vlc2)
    dec_edges = sum(length(v) for v in values(decoded))
    ok = dec_edges == m
    dt2 = round(time() - t2, digits=1)
    println("    Roundtrip: $(ok ? "OK" : "FAIL ($dec_edges vs $m)")  ($(dt2)s)")
    return length(bytes), bpe, ok
end

println("=" ^ 70)
println("Greedy Encoder Optimizations — CNR-2000")
println("Base: adaptive_copy + stop_deltas, w=64, Fibonacci")
println("=" ^ 70)

isfile(CNR_CSV) || error("CNR-2000 CSV not found at $CNR_CSV")

println("\n[1/3] Loading CNR-2000...")
t0 = time()
g = load_adjacency_list_from_csv(CNR_CSV, ',', true)
n, m = nv(g), ne(g)
println("  $n vertices, $m edges ($(round(time()-t0, digits=1))s)")

println("\n[2/3] Leiden+LLP relabeling...")
t_rel = time()
lmap = relabel_leiden_k1_llp(g)
g_rel = relabel_graph(g, lmap)
println("  Done ($(round(time()-t_rel, digits=1))s)")

T = eltype(g_rel)
nls = Dict{T,Vector{T}}()
for v in vertices(g_rel); nls[T(v)] = sort([T(o) for o in outneighbors(g_rel, v)]); end

println("\n[3/3] Running configurations...")
println("─" ^ 70)
results = []

# Winners-only combo: empty + compact + vlc2 + tight
b, bpe, ok = run_config("BEST: empty+compact+vlc2+tight", nls, m;
    empty_prefix=true, compact_copy=true, vlc2=true, tight_intervals=true)
push!(results, ("BEST: empty+compact+vlc2+tight", bpe, ok))

# empty + compact + vlc2 (without tight, in case tight doesn't combine well)
b, bpe, ok = run_config("empty+compact+vlc2", nls, m;
    empty_prefix=true, compact_copy=true, vlc2=true)
push!(results, ("empty+compact+vlc2", bpe, ok))

# empty + compact
b, bpe, ok = run_config("empty+compact", nls, m;
    empty_prefix=true, compact_copy=true)
push!(results, ("empty+compact", bpe, ok))

# empty + vlc2
b, bpe, ok = run_config("empty+vlc2", nls, m;
    empty_prefix=true, vlc2=true)
push!(results, ("empty+vlc2", bpe, ok))

# Summary
println("\n" * "=" ^ 70)
println("SUMMARY")
println("=" ^ 70)
baseline = 2.9559  # known base: adaptive + stop_deltas
for (name, bpe, ok) in results
    d = round(bpe - baseline, digits=4)
    sign = d <= 0 ? "" : "+"
    rt = ok ? "OK" : "FAIL"
    println("  $(rpad(name, 35)) $(rpad(string(bpe), 8)) BPE  $(sign)$(d)  RT=$rt")
end
println()
println("  WebGraph BV reference:        2.897  BPE")
println("  RCGE best:                    2.4341 BPE")
println("=" ^ 70)
