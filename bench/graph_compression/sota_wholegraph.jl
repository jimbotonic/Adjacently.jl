#
# Adjacently: Julia Complex Directed Networks Library
# Copyright (C) 2016-2026 Jimmy Dubuisson <jimmy@dubuisson.ch>
#
# This program is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.
#
# This program is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.
#
# ---------------------------------------------------------------------------
# Reproduces the whole-graph state-of-the-art table: BG / CS / CG at K=1 on the
# :context_range backend, in each dataset's native ordering and under
# Leiden+LLP with the small-cluster merge refinement, against the WebGraph
# BV-HC and Zuckerli baselines.
#
# Two things distinguish this driver from ord_ablation.jl / transfer.jl:
#   1. the backend is :context_range (every structural bit entropy-coded), not
#      Fibonacci;
#   2. the ordered column applies the small-cluster merge
#      (relabel_graph_leiden_llp(g; merge_clusters=tau)) with a per-dataset tau.
# Absolute values are therefore NOT comparable to the ordering-ablation table.
#
# The native column is deterministic and encoded once. The Leiden column is
# stochastic through LLP's visit order, so it is reported as mean +/- std over
# the ordering seeds (leiden_partition itself is deterministic).
#
#   ~/.juliaup/bin/julia --project=. bench/graph_compression/sota_wholegraph.jl [dataset|all]
#
#   dataset in {web-google, amazon-0601, arxiv-hep-ph, eat, cnr-2000,
#               in-2004, enwiki-2013}                        (default: all)
#   ENV: SEEDS=0,1,2   COST_MODEL=0   VERIFY=true   LLP_PASSES=5
#        MERGE=<per-dataset default>   (override the merge threshold; 0 = none)
#        OUT_TSV=<path>                (append one row per cell)
#        WEBGRAPH_CP=<classpath>       (enables the BV-w7 and BV-HC columns)
#        ZUCKERLI_ENCODER=<path>       (enables the Zuckerli max-mode column)
#        EXPORT_DIR=<dir>              (write the ordered graph as .graph-txt/.csr
#                                       for running those tools separately)
#
# BASELINE COLUMNS. BV and Zuckerli are external tools this repository does not
# vendor. Without WEBGRAPH_CP / ZUCKERLI_ENCODER those columns are skipped and
# the driver prints the exact command lines instead, writing the ordered graph
# out (.graph-txt / .csr) when EXPORT_DIR is set. WebGraph .offsets are never
# counted, conservatively against us.
#
# DATASET AVAILABILITY. web-google, amazon-0601, arxiv-hep-ph and EAT run from
# committed data. cnr-2000, in-2004 and enwiki-2013 need the LAW graphs — run
# bench/graph_compression/fetch_datasets.sh first. The committed cnr-2000.mgz
# is pre-reordered, so it CANNOT produce the paper's cnr native row; this
# driver uses the fetched CSV when present and otherwise says so and skips.
#
# Expect differences of up to ~0.08 bpe against the published cells: those used
# slightly different per-cell configurations. See README.md for the per-cell
# comparison on the dataset that was swept end-to-end.
# ---------------------------------------------------------------------------

using Pkg
const REPO_ROOT = normpath(joinpath(@__DIR__, "..", ".."))
Pkg.activate(REPO_ROOT)

using Logging
global_logger(ConsoleLogger(stderr, Logging.Info))
using Printf, Random, Statistics, Dates

using LightGraphs: nv, ne, eltype, outneighbors, vertices
using Adjacently.IO: load_graph_from_pajek, load_adjacency_list_from_csv
using Adjacently.MGS: write_bg_mgs3_graph, write_cs_mgs3_graph,
                      write_cg_mgs3_graph, load_compressed_mgs3_graph
using Adjacently.Graph: get_core
using Adjacently.Relabeling: relabel_graph_leiden_llp
using Adjacently.Compression.CG: CGParams

const COST_MODEL = parse(Int,  get(ENV, "COST_MODEL", "0"))
const LLP_PASSES = parse(Int,  get(ENV, "LLP_PASSES", "5"))
const VERIFY     = parse(Bool, get(ENV, "VERIFY", "true"))
const SEEDS      = parse.(Int, split(get(ENV, "SEEDS", "0,1,2"), ","))
const OUT_TSV    = get(ENV, "OUT_TSV", "")
const WEBGRAPH_CP     = get(ENV, "WEBGRAPH_CP", "")
const ZUCKERLI_ENCODER = get(ENV, "ZUCKERLI_ENCODER", "")
const EXPORT_DIR       = get(ENV, "EXPORT_DIR", "")
const ENCODERS   = ("bg", "cs", "cg")

# Per-dataset spec.
#   core    — reduce to the largest strongly connected component. Set to match
#             the graph each paper row was measured on (its dataset table):
#             web-google 434,818v, amazon-0601 395,234v and EAT 7,754v are
#             cores; arxiv-hep-ph 34,546v is the WHOLE graph (its largest SCC is
#             only 12,711v, which reads ~0.6 bpe lower); the LAW crawls are
#             whole. Spot-checked against the published cells: arxiv native BG
#             8.9537 (paper 8.954), EAT native BG 9.1019 (paper 9.102).
#   merge   — small-cluster merge threshold for the Leiden column: an Int
#             min-size, or nothing for no merge. The LAW crawls are already
#             locality-optimized, so merging is not applied there.
#   bg/cs/cg — encoder parameters, held fixed across both orderings.
const DATASETS = [
    (name="web-google",   loader=:mgz,   file="Web_Google/web-Google.mgz",         core=true,  merge=100,
        bg=(window=64,  lr=false, mr=true),
        cs=(window=256, lr=false),
        cg=(window=8,   intervals=true, lr=true, mil=4, tdd=true)),
    (name="amazon-0601",  loader=:mgz,   file="Amazon_0601/Amazon0601.mgz",        core=true,  merge=100,
        bg=(window=64,  lr=true, mr=true),
        cs=(window=64,  lr=true),
        cg=(window=8,   intervals=true, lr=true, mil=4)),
    (name="arxiv-hep-ph", loader=:mgz,   file="Arxiv_HEP-PH/Cit-HepPh.mgz",        core=false, merge=20,
        bg=(window=64,  lr=true, mr=true),
        cs=(window=256, lr=true),
        cg=(window=256, intervals=true, lr=true, mil=4)),
    (name="eat",          loader=:pajek, file="EAT/EATnew.net",                    core=true,  merge=nothing,
        bg=(window=64,  lr=true, mr=true),
        cs=(window=256, lr=true),
        cg=(window=256, intervals=true, lr=true, mil=4)),
    (name="cnr-2000",     loader=:csv,   file="webgraph/cnr-2000/cnr-2000.csv",    core=false, merge=nothing,
        bg=(window=64,  lr=false, mr=true),
        cs=(window=256, lr=false),
        cg=(window=64,  intervals=true, lr=true, mil=4)),
    (name="in-2004",      loader=:csv,   file="webgraph/in-2004/in-2004.csv",      core=false, merge=nothing,
        bg=(window=64,  lr=true, mr=true),
        cs=(window=256, lr=true),
        cg=(window=8,   intervals=true, lr=true, mil=4)),
    (name="enwiki-2013",  loader=:csv,   file="webgraph/enwiki-2013/enwiki-2013.csv", core=false, merge=nothing,
        bg=(window=64,  lr=true, mr=true),
        cs=(window=128, lr=true),
        cg=(window=64,  intervals=true, lr=true, mil=4)),
]

merge_threshold(spec) =
    haskey(ENV, "MERGE") ? (let m = parse(Int, ENV["MERGE"]); m <= 1 ? nothing : m end) : spec.merge

function build_cg_params(spec)
    return CGParams(
        L=128,
        varint=:fibonacci, count_varint=:fibonacci,
        gap=:fibonacci, degree=:elias_delta,
        undirected_pairs=false,
        perm_strategy=:blockpos, membership=:implicit_ranges,
        inter_strategy=:perm,
        intra_ref_enabled=true, intra_ref_window=spec.window,
        intra_ref_rle=false,
        intra_block_try=false,
        positions_mode=:delta, additions_mode=:delta,
        min_cluster_density=0.0,
        intra_intervals=spec.intervals, intra_mil=spec.mil,
        intra_lr_split=spec.lr,
        intra_tight_deltas = hasproperty(spec, :tdd) ? spec.tdd : false,
        intra_greedy_mil=false,
        intra_zigzag=true, intra_stop_deltas=true,
        intra_copy_blocks=true, intra_copy_adaptive=true,
        intra_ref_fixwidth=true, intra_ref_vlc=false,
        intra_add_adaptive=true, intra_raw_adaptive=true,
        intra_adapt_mil=spec.mil,
    )
end

# Returns the graph, or nothing when the dataset has to be fetched first.
function load_graph(spec)
    path = joinpath(REPO_ROOT, "datasets", spec.file)
    if !isfile(path)
        if spec.loader == :csv
            @warn """$(spec.name): $(spec.file) not found — skipping.
                     Fetch the LAW graphs first:
                       bash bench/graph_compression/fetch_datasets.sh $(spec.name)"""
            spec.name == "cnr-2000" && @warn """The committed cnr-2000.mgz is NOT usable here: it is a
                     pre-reordered graph, so it cannot produce the native row."""
            return nothing
        end
        error("dataset not found (committed?): $path")
    end
    full = spec.loader == :pajek ? load_graph_from_pajek(path) :
           spec.loader == :csv   ? load_adjacency_list_from_csv(path, ',', true) :
                                   load_compressed_mgs3_graph(path)
    g = spec.core ? get_core(full)[1] : full
    @info @sprintf("  %s: %dv/%de%s", spec.name, nv(g), ne(g),
                   spec.core ? " (largest SCC)" : "")
    return g
end

function verify_roundtrip(g_orig, g_dec, label)
    nv(g_orig) == nv(g_dec) || error("$label roundtrip: nv mismatch")
    ne(g_orig) == ne(g_dec) || error("$label roundtrip: ne mismatch")
    for v in vertices(g_orig)
        sort(collect(outneighbors(g_orig, v))) == sort(collect(outneighbors(g_dec, v))) ||
            error("$label roundtrip: neighbor mismatch at v=$v")
    end
end

# Encode with one of our encoders on the context-range backend; BPE = 8*bytes/m.
function encode_bpe(g, enc::AbstractString, spec, tmpdir, m)
    base = joinpath(tmpdir, enc)
    if enc == "bg"
        write_bg_mgs3_graph(g, base; integer_encoding=:context_range,
            ref_window_size=spec.bg.window, copy_blocks=true, stop_deltas=true,
            lr_split=spec.bg.lr, multi_ref=spec.bg.mr, cost_model=COST_MODEL)
    elseif enc == "cs"
        write_cs_mgs3_graph(g, base; integer_encoding=:context_range,
            ref_window_size=spec.cs.window, lr_split=spec.cs.lr, cost_model=COST_MODEL)
    else # cg, K=1
        T = eltype(g)
        write_cg_mgs3_graph(g, base, [collect(T(1):T(nv(g)))];
            params=build_cg_params(spec.cg), integer_encoding=:context_range)
    end
    bpe = 8.0 * filesize(base * ".mgz") / m
    if VERIFY
        verify_roundtrip(g, load_compressed_mgs3_graph(base * ".mgz"), uppercase(enc))
        GC.gc()
    end
    return bpe
end

# ── external baselines ──────────────────────────────────────────────────────

# 0-based ASCIIGraph, the input format of WebGraph's BVGraph converter.
function write_graph_txt(g, path)
    open(path, "w") do io
        println(io, nv(g))
        for v in vertices(g)
            println(io, join((string(Int(u) - 1) for u in sort(collect(outneighbors(g, v)))), " "))
        end
    end
    return path
end

# Zuckerli CSR: uint64 fingerprint, uint32 N, (N+1) uint64 offsets, M uint32 dsts.
function write_zuckerli_csr(g, path)
    n = Int(nv(g))
    open(path, "w") do io
        write(io, UInt64((sizeof(UInt64) << 4) | sizeof(UInt32)))
        write(io, UInt32(n))
        off = UInt64(0)
        write(io, off)
        degs = Vector{Int}(undef, n)
        for v in vertices(g)
            degs[Int(v)] = length(outneighbors(g, v))
        end
        for i in 1:n
            off += UInt64(degs[i])
            write(io, off)
        end
        for v in vertices(g)
            for u in sort(collect(outneighbors(g, v)))
                write(io, UInt32(Int(u) - 1))
            end
        end
    end
    return path
end

# BV variants: "bv-w7" = default, random-access capable; "bv-hc" = high
# compression (w=16, unbounded reference chains, no practical random access).
function bv_bpe(g, variant::AbstractString, tmpdir, m)
    isempty(WEBGRAPH_CP) && return NaN
    txt = joinpath(tmpdir, "bv.graph-txt")
    isfile(txt) || write_graph_txt(g, txt)
    bp = joinpath(tmpdir, "bv_" * variant)
    w, maxref = variant == "bv-hc" ? ("16", "-1") : ("7", "3")
    cmd = `java -Xmx8G -cp $WEBGRAPH_CP it.unimi.dsi.webgraph.BVGraph
           -g it.unimi.dsi.webgraph.ASCIIGraph -w $w -m $maxref $bp $bp`
    try
        run(pipeline(cmd, stdout=devnull, stderr=devnull))
    catch e
        @warn "$variant failed (WebGraph): $e"
        return NaN
    end
    # .offsets deliberately not counted
    return isfile(bp * ".graph") ? 8.0 * filesize(bp * ".graph") / m : NaN
end

function zuckerli_bpe(g, tmpdir, m)
    isempty(ZUCKERLI_ENCODER) && return NaN
    csr = joinpath(tmpdir, "graph.csr")
    isfile(csr) || write_zuckerli_csr(g, csr)
    out = joinpath(tmpdir, "graph.zk")
    try
        run(pipeline(`$ZUCKERLI_ENCODER --input_path $csr --output_path $out`,
                     stdout=devnull, stderr=devnull))
    catch e
        @warn "Zuckerli failed: $e"
        return NaN
    end
    return isfile(out) ? 8.0 * filesize(out) / m : NaN
end

# When a baseline tool is not configured, tell the user how to produce that
# column by hand. Exports the ordered graph only if EXPORT_DIR is set — these
# files are large (enwiki-2013 runs to several GB), so writing them is opt-in;
# otherwise the commands are printed with a placeholder path.
function baseline_hint(spec, ordering, seed, g, m)
    (isempty(WEBGRAPH_CP) || isempty(ZUCKERLI_ENCODER)) || return
    tag = "$(spec.name)_$(ordering)_seed$(seed)"
    exported(ext, writer) = if isempty(EXPORT_DIR)
        "<$ordering-ordered graph>$ext   # set EXPORT_DIR to have this driver write it"
    else
        mkpath(EXPORT_DIR)
        writer(g, joinpath(EXPORT_DIR, tag * ext))
    end
    println("\n  [$(spec.name)/$ordering] baseline columns skipped — to produce them:")
    if isempty(WEBGRAPH_CP)
        txt = exported(".graph-txt", write_graph_txt)
        println("    # BV-w7 / BV-HC  (bpe = 8 * filesize(.graph) / $m, offsets excluded)")
        println("    java -cp '<webgraph-3.6.12 + deps>' it.unimi.dsi.webgraph.BVGraph \\")
        println("         -g it.unimi.dsi.webgraph.ASCIIGraph -w 7  -m 3  $txt <out>")
        println("    java -cp '<webgraph-3.6.12 + deps>' it.unimi.dsi.webgraph.BVGraph \\")
        println("         -g it.unimi.dsi.webgraph.ASCIIGraph -w 16 -m -1 $txt <out>")
    end
    if isempty(ZUCKERLI_ENCODER)
        csr = exported(".csr", write_zuckerli_csr)
        println("    # Zuckerli max-compression (bpe = 8 * filesize(.zk) / $m)")
        println("    <zuckerli>/build/encoder --input_path $csr --output_path <out>.zk")
    end
end

# ── driver ──────────────────────────────────────────────────────────────────

function emit_tsv(rows)
    isempty(OUT_TSV) && return
    new = !isfile(OUT_TSV)
    open(OUT_TSV, "a") do io
        new && println(io, "timestamp\tdataset\tordering\tseed\tmethod\tbpe")
        for r in rows
            @printf(io, "%s\t%s\t%s\t%d\t%s\t%.4f\n",
                    Dates.now(), r.dataset, r.ordering, r.seed, r.method, r.bpe)
        end
    end
end

function run_dataset(spec)
    tau = merge_threshold(spec)
    println("=" ^ 72)
    println("  sota_wholegraph — $(spec.name)   seeds=$(SEEDS) merge=$(tau === nothing ? "none" : tau)")
    println("=" ^ 72)
    g0 = load_graph(spec)
    g0 === nothing && return nothing
    m = ne(g0)
    tmpdir = mktempdir()

    # results[ordering][method] = Vector over seeds
    methods = [ENCODERS...; "bv-w7"; "bv-hc"; "zuckerli"]
    results = Dict(o => Dict(mt => Float64[] for mt in methods) for o in ("native", "leiden"))

    for ordering in ("native", "leiden")
        seeds = ordering == "native" ? [SEEDS[1]] : SEEDS
        for seed in seeds
            g = g0
            if ordering == "leiden"
                Random.seed!(seed)
                g, _ = relabel_graph_leiden_llp(g0; llp_mode=:sym, llp_passes=LLP_PASSES,
                                                merge_clusters=tau)
            end
            seed_dir = mktempdir(tmpdir)
            rows = NamedTuple[]   # flushed per (ordering, seed): a long run keeps its results

            for enc in ENCODERS
                b = encode_bpe(g, enc, spec, seed_dir, m)
                push!(results[ordering][enc], b)
                push!(rows, (dataset=spec.name, ordering=ordering, seed=seed,
                             method=uppercase(enc) * "-ctx", bpe=b))
                @info @sprintf("  [%-6s seed=%d] %s-ctx = %.4f BPE", ordering, seed, uppercase(enc), b)
            end

            for (mt, val) in (("bv-w7", bv_bpe(g, "bv-w7", seed_dir, m)),
                              ("bv-hc", bv_bpe(g, "bv-hc", seed_dir, m)),
                              ("zuckerli", zuckerli_bpe(g, seed_dir, m)))
                isnan(val) && continue
                push!(results[ordering][mt], val)
                push!(rows, (dataset=spec.name, ordering=ordering, seed=seed,
                             method=mt, bpe=val))
                @info @sprintf("  [%-6s seed=%d] %-8s = %.4f BPE", ordering, seed, mt, val)
            end

            emit_tsv(rows)
            seed == seeds[1] && baseline_hint(spec, ordering, seed, g, m)
            ordering == "leiden" && (g = nothing; GC.gc())
        end
    end

    fmt(v) = isempty(v) ? "n/a" :
             length(v) > 1 ? @sprintf("%.4f±%.4f", mean(v), std(v)) : @sprintf("%.4f", v[1])
    println("\n  $(spec.name) — whole-graph BPE (context-range backend), m=$m edges")
    @printf("  %-8s %14s %14s %14s %14s %14s %9s\n",
            "Ord.", "BV-HC", "Zuck.", "BG", "CS", "CG", "vs Zuck.")
    for (o, label) in (("native", "native"), ("leiden", "leiden"))
        ours = [results[o][e] for e in ENCODERS]
        best = minimum(mean(v) for v in ours if !isempty(v))
        z = results[o]["zuckerli"]
        margin = isempty(z) ? "n/a" : @sprintf("%+.1f%%", 100 * (mean(z) - best) / mean(z))
        @printf("  %-8s %14s %14s %14s %14s %14s %9s\n", label,
                fmt(results[o]["bv-hc"]), fmt(z),
                fmt(results[o]["bg"]), fmt(results[o]["cs"]), fmt(results[o]["cg"]), margin)
    end
    println("=" ^ 72, "\n")

    rm(tmpdir; recursive=true)
    return results
end

function main()
    which = isempty(ARGS) ? "all" : lowercase(ARGS[1])
    specs = which == "all" ? DATASETS : filter(s -> s.name == which, DATASETS)
    isempty(specs) && error("unknown dataset '$which'; choose from " *
                            join([s.name for s in DATASETS], ", ") * ", or 'all'")
    isempty(WEBGRAPH_CP) && @info "WEBGRAPH_CP unset — BV columns skipped"
    isempty(ZUCKERLI_ENCODER) && @info "ZUCKERLI_ENCODER unset — Zuckerli column skipped"
    for spec in specs
        run_dataset(spec)
    end
end

main()
