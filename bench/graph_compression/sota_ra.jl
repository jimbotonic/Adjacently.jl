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
# Reproduces the random-access state-of-the-art table: each system measured in
# its own random-access mode, so the comparison is like-for-like on what a
# seekable file costs.
#
#   BG-RA / CS-RA : :context_range + :index, k-vertex chunks (sampled index)
#   CG-RA         : :context_range + :index, cluster = seek chunk
#   BV-w7         : BVGraph w=7, m=3 — the random-access-capable BV variant
#                   (BV-HC has unbounded reference chains and no practical
#                   random access, which is why it is absent here)
#   Zuckerli-RA   : encoder --allow_random_access
#
# SEEK GRANULARITY. All three of our encoders use the same chunk size
# (SAMPLE_K, default 2048 vertices) so the columns are comparable: BG/CS pass
# it as index_sample_k, and CG coalesces its cluster ranges until each chunk
# holds at least that many vertices. That single granularity is what reproduces
# the published cells — on EAT it gives BG 9.1295 / CS 9.0713 / CG native
# 12.5918 / CG leiden 10.3951 against published 9.130 / 9.081 / 12.608 / 10.411.
#
# Verification here is a real seek test, not a full decode: the encoded file is
# reopened with load_indexed_mgs3_graph and a sample of vertices is fetched
# through neighbors(v), exercising the chunk seek and the cross-chunk reference
# resolver.
#
#   ~/.juliaup/bin/julia --project=. bench/graph_compression/sota_ra.jl [dataset|all]
#
#   dataset in {web-google, amazon-0601, arxiv-hep-ph, eat, cnr-2000,
#               in-2004, enwiki-2013}                        (default: all)
#   ENV: SEEDS=0,1,2   COST_MODEL=0   VERIFY=true   LLP_PASSES=5
#        SAMPLE_K=2048                 (seek granularity, multiple of 4)
#        MERGE=<per-dataset default>   (merge threshold for the ordering; 0 = none)
#        OUT_TSV=<path>                (append one row per cell)
#        WEBGRAPH_CP=<classpath>       (enables the BV-w7 column)
#        ZUCKERLI_ENCODER=<path>       (enables the Zuckerli-RA column)
#        EXPORT_DIR=<dir>              (write the ordered graph for those tools)
#
# See sota_wholegraph.jl for the whole-graph half of the comparison, and
# README.md for which datasets need fetch_datasets.sh first.
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
                      write_cg_mgs3_graph, load_compressed_mgs3_graph,
                      load_indexed_mgs3_graph
using Adjacently.Graph: get_core
using Adjacently.Relabeling: relabel_graph_leiden_llp
using Adjacently.Compression.CG: CGParams

const COST_MODEL = parse(Int,  get(ENV, "COST_MODEL", "0"))
const LLP_PASSES = parse(Int,  get(ENV, "LLP_PASSES", "5"))
const VERIFY     = parse(Bool, get(ENV, "VERIFY", "true"))
const SEEDS      = parse.(Int, split(get(ENV, "SEEDS", "0,1,2"), ","))
const SAMPLE_K   = parse(Int,  get(ENV, "SAMPLE_K", "2048"))
const OUT_TSV    = get(ENV, "OUT_TSV", "")
const WEBGRAPH_CP      = get(ENV, "WEBGRAPH_CP", "")
const ZUCKERLI_ENCODER = get(ENV, "ZUCKERLI_ENCODER", "")
const EXPORT_DIR       = get(ENV, "EXPORT_DIR", "")
const ENCODERS   = ("bg", "cs", "cg")
const VERIFY_SAMPLES = parse(Int, get(ENV, "VERIFY_SAMPLES", "200"))

SAMPLE_K % 4 == 0 || error("SAMPLE_K must be a multiple of 4 (got $SAMPLE_K)")

# Same per-dataset specs as sota_wholegraph.jl — see that file for why `core`
# differs per dataset (the paper's rows are not on a uniform reduction).
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

function load_graph(spec)
    path = joinpath(REPO_ROOT, "datasets", spec.file)
    if !isfile(path)
        if spec.loader == :csv
            @warn """$(spec.name): $(spec.file) not found — skipping.
                     Fetch the LAW graphs first:
                       bash bench/graph_compression/fetch_datasets.sh $(spec.name)"""
            return nothing
        end
        error("dataset not found (committed?): $path")
    end
    full = spec.loader == :pajek ? load_graph_from_pajek(path) :
           spec.loader == :csv   ? load_adjacency_list_from_csv(path, ',', true) :
                                   load_compressed_mgs3_graph(path)
    g = spec.core ? get_core(full)[1] : full
    @info @sprintf("  %s: %dv/%de%s", spec.name, nv(g), ne(g), spec.core ? " (largest SCC)" : "")
    return g
end

# CG seek chunks as contiguous vertex ranges of at least SAMPLE_K vertices.
# `segs` are the Leiden cluster sizes in concatenation order (from
# relabel_graph_leiden_llp(...; return_clusters=true)); consecutive clusters are
# coalesced until a chunk reaches the granularity. In the native ordering there
# are no cluster ranges, so fixed SAMPLE_K-vertex chunks are used instead.
function ra_chunks(g, segs::Union{Nothing,Vector{Int}})
    n = Int(nv(g)); T = eltype(g)
    sizes = if segs === nothing
        [min(SAMPLE_K, n - i + 1) for i in 1:SAMPLE_K:n]
    else
        acc = Int[]; cur = 0
        for s in segs
            cur += s
            cur >= SAMPLE_K && (push!(acc, cur); cur = 0)
        end
        cur > 0 && (isempty(acc) ? push!(acc, cur) : (acc[end] += cur))
        acc
    end
    chunks = Vector{Vector{T}}(); off = 0
    for s in sizes
        push!(chunks, collect(T(off + 1):T(off + s)))
        off += s
    end
    return chunks
end

# Random access is the point of this table, so verify by seeking rather than by
# decoding the whole file: reopen through the index and pull a spread of
# vertices, which exercises the chunk seek and the cross-chunk ref resolver.
function verify_random_access(g, path, label)
    idx = load_indexed_mgs3_graph(path)
    n = Int(nv(g))
    Int(idx.n) == n || error("$label RA: nv mismatch ($(idx.n) vs $n)")
    step = max(1, n ÷ max(1, VERIFY_SAMPLES))
    for v in 1:step:n
        got = sort(collect(idx.neighbors(v)))
        want = sort(collect(outneighbors(g, eltype(g)(v))))
        length(got) == length(want) && all(Int(a) == Int(b) for (a, b) in zip(got, want)) ||
            error("$label RA: neighbor mismatch at v=$v")
    end
end

function encode_bpe(g, enc::AbstractString, spec, chunks, tmpdir, m)
    base = joinpath(tmpdir, enc)
    if enc == "bg"
        write_bg_mgs3_graph(g, base; integer_encoding=:context_range, coding_scheme=:index,
            ref_window_size=spec.bg.window, copy_blocks=true, stop_deltas=true,
            lr_split=spec.bg.lr, multi_ref=spec.bg.mr, cost_model=COST_MODEL,
            index_sample_k=SAMPLE_K)
    elseif enc == "cs"
        write_cs_mgs3_graph(g, base; integer_encoding=:context_range, coding_scheme=:index,
            ref_window_size=spec.cs.window, lr_split=spec.cs.lr, cost_model=COST_MODEL,
            index_sample_k=SAMPLE_K)
    else # cg: cluster = seek chunk
        write_cg_mgs3_graph(g, base, chunks; params=build_cg_params(spec.cg),
            integer_encoding=:context_range, coding_scheme=:index)
    end
    bpe = 8.0 * filesize(base * ".mgz") / m
    if VERIFY
        verify_random_access(g, base * ".mgz", uppercase(enc))
        GC.gc()
    end
    return bpe
end

# ── external baselines ──────────────────────────────────────────────────────

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

# BV at w=7, m=3: the random-access-capable variant. .offsets not counted.
function bv_ra_bpe(g, tmpdir, m)
    isempty(WEBGRAPH_CP) && return NaN
    txt = joinpath(tmpdir, "bv.graph-txt")
    isfile(txt) || write_graph_txt(g, txt)
    bp = joinpath(tmpdir, "bv_w7")
    cmd = `java -Xmx8G -cp $WEBGRAPH_CP it.unimi.dsi.webgraph.BVGraph
           -g it.unimi.dsi.webgraph.ASCIIGraph -w 7 -m 3 $bp $bp`
    try
        run(pipeline(cmd, stdout=devnull, stderr=devnull))
    catch e
        @warn "BV-w7 failed (WebGraph): $e"
        return NaN
    end
    return isfile(bp * ".graph") ? 8.0 * filesize(bp * ".graph") / m : NaN
end

function zuckerli_ra_bpe(g, tmpdir, m)
    isempty(ZUCKERLI_ENCODER) && return NaN
    csr = joinpath(tmpdir, "graph.csr")
    isfile(csr) || write_zuckerli_csr(g, csr)
    out = joinpath(tmpdir, "graph_ra.zk")
    try
        run(pipeline(`$ZUCKERLI_ENCODER --input_path $csr --output_path $out --allow_random_access`,
                     stdout=devnull, stderr=devnull))
    catch e
        @warn "Zuckerli-RA failed: $e"
        return NaN
    end
    return isfile(out) ? 8.0 * filesize(out) / m : NaN
end

function baseline_hint(spec, ordering, seed, g, m)
    (isempty(WEBGRAPH_CP) || isempty(ZUCKERLI_ENCODER)) || return
    tag = "$(spec.name)_$(ordering)_seed$(seed)"
    exported(ext, writer) = if isempty(EXPORT_DIR)
        "<$tag$ext>"
    else
        mkpath(EXPORT_DIR)
        writer(g, joinpath(EXPORT_DIR, tag * ext))
    end
    println("\n  [$(spec.name)/$ordering] baseline columns skipped — to produce them:")
    isempty(EXPORT_DIR) &&
        println("    (set EXPORT_DIR to have this driver write the <...> input files)")
    if isempty(WEBGRAPH_CP)
        txt = exported(".graph-txt", write_graph_txt)
        println("    # BV-w7, the RA-capable variant (bpe = 8 * filesize(.graph) / $m)")
        println("    java -cp '<webgraph-3.6.12 + deps>' it.unimi.dsi.webgraph.BVGraph \\")
        println("         -g it.unimi.dsi.webgraph.ASCIIGraph -w 7 -m 3 $txt <out>")
    end
    if isempty(ZUCKERLI_ENCODER)
        csr = exported(".csr", write_zuckerli_csr)
        println("    # Zuckerli random-access mode (bpe = 8 * filesize(.zk) / $m)")
        println("    <zuckerli>/build/encoder --input_path $csr --output_path <out>.zk \\")
        println("         --allow_random_access")
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
    println("  sota_ra — $(spec.name)   seeds=$(SEEDS) merge=$(tau === nothing ? "none" : tau) k=$SAMPLE_K")
    println("=" ^ 72)
    g0 = load_graph(spec)
    g0 === nothing && return nothing
    m = ne(g0)
    tmpdir = mktempdir()

    methods = [ENCODERS...; "bv-w7"; "zuckerli-ra"]
    results = Dict(o => Dict(mt => Float64[] for mt in methods) for o in ("native", "leiden"))

    for ordering in ("native", "leiden")
        seeds = ordering == "native" ? [SEEDS[1]] : SEEDS
        for seed in seeds
            g, segs = g0, nothing
            if ordering == "leiden"
                Random.seed!(seed)
                g, _, segs = relabel_graph_leiden_llp(g0; llp_mode=:sym, llp_passes=LLP_PASSES,
                                                      merge_clusters=tau, return_clusters=true)
            end
            chunks = ra_chunks(g, segs)
            seed_dir = mktempdir(tmpdir)
            rows = NamedTuple[]   # flushed per (ordering, seed)

            for enc in ENCODERS
                b = encode_bpe(g, enc, spec, chunks, seed_dir, m)
                push!(results[ordering][enc], b)
                push!(rows, (dataset=spec.name, ordering=ordering, seed=seed,
                             method=uppercase(enc) * "-RA", bpe=b))
                @info @sprintf("  [%-6s seed=%d] %s-RA = %.4f BPE%s", ordering, seed,
                               uppercase(enc), b, enc == "cg" ? " ($(length(chunks)) chunks)" : "")
            end

            for (mt, val) in (("bv-w7", bv_ra_bpe(g, seed_dir, m)),
                              ("zuckerli-ra", zuckerli_ra_bpe(g, seed_dir, m)))
                isnan(val) && continue
                push!(results[ordering][mt], val)
                push!(rows, (dataset=spec.name, ordering=ordering, seed=seed, method=mt, bpe=val))
                @info @sprintf("  [%-6s seed=%d] %-11s = %.4f BPE", ordering, seed, mt, val)
            end

            emit_tsv(rows)
            seed == seeds[1] && baseline_hint(spec, ordering, seed, g, m)
            ordering == "leiden" && (g = nothing; GC.gc())
        end
    end

    fmt(v) = isempty(v) ? "n/a" :
             length(v) > 1 ? @sprintf("%.4f±%.4f", mean(v), std(v)) : @sprintf("%.4f", v[1])
    println("\n  $(spec.name) — random-access BPE (context-range, k=$SAMPLE_K), m=$m edges")
    @printf("  %-8s %14s %14s %14s %14s %14s %9s\n",
            "Ord.", "BV-w7", "Zuck.-RA", "BG", "CS", "CG", "vs Zuck.")
    for o in ("native", "leiden")
        ours = [results[o][e] for e in ENCODERS]
        best = minimum(mean(v) for v in ours if !isempty(v))
        z = results[o]["zuckerli-ra"]
        margin = isempty(z) ? "n/a" : @sprintf("%+.1f%%", 100 * (mean(z) - best) / mean(z))
        @printf("  %-8s %14s %14s %14s %14s %14s %9s\n", o,
                fmt(results[o]["bv-w7"]), fmt(z),
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
    isempty(WEBGRAPH_CP) && @info "WEBGRAPH_CP unset — BV-w7 column skipped"
    isempty(ZUCKERLI_ENCODER) && @info "ZUCKERLI_ENCODER unset — Zuckerli-RA column skipped"
    for spec in specs
        run_dataset(spec)
    end
end

main()
