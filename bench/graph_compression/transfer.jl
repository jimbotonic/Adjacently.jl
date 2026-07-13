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
# Reproduces the cross-encoder transfer table: the Leiden+LLP-vs-plain-LLP gain
# (bits per edge saved) per encoder, as mean +/- std over three ordering seeds,
# under both the Fibonacci and the context-range backend, for the four
# weakly-ordered SNAP datasets that ship with the repository. Also reports the
# cross-encoder spread of the mean gain (the paper's transfer-invariance
# number). Gains are paired per seed: gain = BPE(LLP) - BPE(Leiden+LLP).
#
#   ~/.juliaup/bin/julia --project=. bench/graph_compression/transfer.jl [dataset|all]
#
#   dataset in {web-google, amazon-0601, arxiv-hep-ph, eat}   (default: all)
#   ENV: SEEDS=0,1,2   BACKENDS=fib,ctx   COST_MODEL=0   VERIFY=true
#
# enwiki-2013 (the fifth weakly-ordered dataset in the paper table) needs the
# LAW-fetch driver; see README.md.
# ---------------------------------------------------------------------------

using Pkg
const REPO_ROOT = normpath(joinpath(@__DIR__, "..", ".."))
Pkg.activate(REPO_ROOT)

using Logging
global_logger(ConsoleLogger(stderr, Logging.Info))
using Printf, Random, Statistics

using LightGraphs: nv, ne, eltype, outneighbors, vertices
using Adjacently.IO: load_graph_from_pajek
using Adjacently.MGS: write_bg_mgs3_graph, write_cs_mgs3_graph,
                      write_cg_mgs3_graph, load_compressed_mgs3_graph
using Adjacently.Graph: get_core
using Adjacently.Relabeling: relabel_graph_llp, relabel_graph_leiden_llp
using Adjacently.Compression.CG: CGParams

const COST_MODEL = parse(Int,  get(ENV, "COST_MODEL", "0"))
const LLP_PASSES = parse(Int,  get(ENV, "LLP_PASSES", "5"))
const VERIFY     = parse(Bool, get(ENV, "VERIFY", "true"))
const SEEDS      = parse.(Int, split(get(ENV, "SEEDS", "0,1,2"), ","))
const BACKENDS   = split(lowercase(get(ENV, "BACKENDS", "fib,ctx")), ",")
const ORDERINGS  = ("llp", "leiden_llp")

# Same per-dataset spec and encoder params as ord_ablation.jl (held fixed across
# orderings and backends). See ord_ablation.jl for the get_core rationale.
const DATASETS = [
    (name="web-google",   loader=:mgz,   file="Web_Google/web-Google.mgz",
        bg=(window=64,  lr=false, mr=true),
        cs=(window=256, lr=false),
        cg=(window=8,   intervals=true, lr=true, mil=4, tdd=true)),
    (name="amazon-0601",  loader=:mgz,   file="Amazon_0601/Amazon0601.mgz",
        bg=(window=64,  lr=true, mr=true),
        cs=(window=64,  lr=true),
        cg=(window=8,   intervals=true, lr=true, mil=4)),
    (name="arxiv-hep-ph", loader=:mgz,   file="Arxiv_HEP-PH/Cit-HepPh.mgz",
        bg=(window=64,  lr=true, mr=true),
        cs=(window=256, lr=true),
        cg=(window=256, intervals=true, lr=true, mil=4)),
    (name="eat",          loader=:pajek, file="EAT/EATnew.net",
        bg=(window=64,  lr=true, mr=true),
        cs=(window=256, lr=true),
        cg=(window=256, intervals=true, lr=true, mil=4)),
]

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

function load_core(spec)
    path = joinpath(REPO_ROOT, "datasets", spec.file)
    isfile(path) || error("dataset not found (committed?): $path")
    full = spec.loader == :pajek ? load_graph_from_pajek(path) :
                                    load_compressed_mgs3_graph(path)
    core = get_core(full)[1]
    @info @sprintf("  %s: core %dv/%de", spec.name, nv(core), ne(core))
    return core
end

function reorder(g, ordering::AbstractString, seed::Int)
    Random.seed!(seed)
    if ordering == "llp"
        g_rel, _ = relabel_graph_llp(g; llp_mode=:sym, passes=LLP_PASSES)
    else
        g_rel, _ = relabel_graph_leiden_llp(g; llp_mode=:sym, llp_passes=LLP_PASSES,
                                            merge_clusters=nothing)
    end
    return g_rel
end

function verify_roundtrip(g_orig, g_dec, label)
    nv(g_orig) == nv(g_dec) || error("$label roundtrip: nv mismatch")
    ne(g_orig) == ne(g_dec) || error("$label roundtrip: ne mismatch")
    for v in vertices(g_orig)
        sort(collect(outneighbors(g_orig, v))) == sort(collect(outneighbors(g_dec, v))) ||
            error("$label roundtrip: neighbor mismatch at v=$v")
    end
end

function encode_bpe(g, enc::AbstractString, spec, ie::Symbol, tmpdir, m)
    base = joinpath(tmpdir, enc)
    if enc == "bg"
        write_bg_mgs3_graph(g, base; integer_encoding=ie,
            ref_window_size=spec.bg.window, copy_blocks=true, stop_deltas=true,
            lr_split=spec.bg.lr, multi_ref=spec.bg.mr, cost_model=COST_MODEL)
    elseif enc == "cs"
        write_cs_mgs3_graph(g, base; integer_encoding=ie,
            ref_window_size=spec.cs.window, lr_split=spec.cs.lr, cost_model=COST_MODEL)
    else
        T = eltype(g)
        write_cg_mgs3_graph(g, base, [collect(T(1):T(nv(g)))];
            params=build_cg_params(spec.cg), integer_encoding=ie)
    end
    bpe = 8.0 * filesize(base * ".mgz") / m
    if VERIFY
        verify_roundtrip(g, load_compressed_mgs3_graph(base * ".mgz"), uppercase(enc))
        GC.gc()
    end
    return bpe
end

function run_dataset(spec)
    println("=" ^ 64)
    println("  transfer — $(spec.name)   seeds=$(SEEDS) backends=$(BACKENDS)")
    println("=" ^ 64)
    g0 = load_core(spec)
    m  = ne(g0)
    tmpdir = mktempdir()

    # bpe[backend][ordering][encoder] = Vector over seeds
    bpe = Dict(bk => Dict(o => Dict(e => Float64[] for e in ("bg","cs","cg"))
                          for o in ORDERINGS) for bk in BACKENDS)
    for seed in SEEDS
        for ordering in ORDERINGS
            g = reorder(g0, ordering, seed)
            for bk in BACKENDS
                ie = bk == "ctx" ? :context_range : :fibonacci
                for enc in ("bg","cs","cg")
                    b = encode_bpe(g, enc, spec, ie, tmpdir, m)
                    push!(bpe[bk][ordering][enc], b)
                    @info @sprintf("  [%s %-10s seed=%d] %s = %.4f", bk, ordering, seed, uppercase(enc), b)
                end
            end
            g = nothing; GC.gc()
        end
    end
    rm(tmpdir; recursive=true)

    println("\n  $(spec.name) — Leiden+LLP-vs-plain-LLP gain (BPE saved), m=$m")
    @printf("  %-5s %13s %13s %13s %9s\n", "Back.", "BG", "CS", "CG", "Spread")
    for bk in BACKENDS
        means = Float64[]
        cells = String[]
        for enc in ("bg","cs","cg")
            a = bpe[bk]["llp"][enc]; b = bpe[bk]["leiden_llp"][enc]
            gains = a .- b                          # paired per seed
            push!(means, mean(gains))
            push!(cells, @sprintf("%.3f±%.3f", mean(gains), std(gains)))
        end
        spread = maximum(means) - minimum(means)
        @printf("  %-5s %13s %13s %13s %9.3f\n", bk, cells[1], cells[2], cells[3], spread)
    end
    println("=" ^ 64, "\n")
end

function main()
    which = isempty(ARGS) ? "all" : lowercase(ARGS[1])
    specs = which == "all" ? DATASETS : filter(s -> s.name == which, DATASETS)
    isempty(specs) && error("unknown dataset '$which'; choose from " *
                            join([s.name for s in DATASETS], ", ") * ", or 'all'")
    for spec in specs
        run_dataset(spec)
    end
end

main()
