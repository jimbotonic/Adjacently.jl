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
# Reproduces the ordering-ablation table (Fibonacci backend) for the four
# weakly-ordered SNAP datasets that ship with the repository. For each dataset
# it extracts the largest strongly-connected core, builds three orderings
# (original / plain LLP / Leiden+LLP), and encodes each with BG, CS and CG
# (K=1), reporting bits-per-edge. The LLP and Leiden+LLP orderings are
# regenerated in-memory under fixed seeds (0,1,2); "original" is the core's
# native ascending-ID order and is deterministic (encoded once).
#
# BPE is a property of the encoded bitstream, hence machine-independent. The
# encoder configuration is held fixed across the three orderings (a true
# ordering ablation); only the vertex order varies.
#
# Verified against the paper's EAT row: every LLP and Leiden+LLP cell and the
# Orig. BG/CS cells reproduce to the printed 3 decimals. The one exception is
# the Orig. CG cell, which in the paper table was populated from the per-dataset
# best-config CG run (tight-deltas + gamma gap coding, which suits the large gaps
# of the unordered graph) rather than this fixed fibonacci config; here it reads
# ~0.08 bpe higher. The ordering-transfer result this driver exists to support
# (LLP / Leiden+LLP rows and the gains between them) reproduces exactly.
#
#   ~/.juliaup/bin/julia --project=. bench/graph_compression/ord_ablation.jl [dataset|all]
#
#   dataset in {web-google, amazon-0601, arxiv-hep-ph, eat}   (default: all)
#   ENV: SEEDS=0,1,2   COST_MODEL=0   VERIFY=true
#
# The two largest LAW graphs (in-2004, enwiki-2013) and cnr-2000 are handled by
# separate drivers; see README.md.
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
const ORDERINGS  = ("original", "llp", "leiden_llp")

# Per-dataset spec. `file` is committed under datasets/; the graph is decoded in
# full and reduced to its largest SCC (get_core) — reproducing the non-committed
# *_core.mgs the research drivers consumed. Encoder params match the ablation
# table and are identical across orderings.
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
    @info @sprintf("  %s: full %dv/%de -> core %dv/%de",
                   spec.name, nv(full), ne(full), nv(core), ne(core))
    return core
end

# Native order for "original"; seeded LLP / Leiden+LLP otherwise. Seeding the
# global RNG is what makes the LLP layers reproducible per seed (leiden_partition
# is deterministic; the variance lives in LLP's shuffled visit order).
function reorder(g, ordering::AbstractString, seed::Int)
    ordering == "original" && return g
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

# Encode g with one encoder; return BPE = 8 * bytes / m.
function encode_bpe(g, enc::AbstractString, spec, tmpdir, m)
    base = joinpath(tmpdir, enc)
    if enc == "bg"
        write_bg_mgs3_graph(g, base; integer_encoding=:fibonacci,
            ref_window_size=spec.bg.window, copy_blocks=true, stop_deltas=true,
            lr_split=spec.bg.lr, multi_ref=spec.bg.mr, cost_model=COST_MODEL)
    elseif enc == "cs"
        write_cs_mgs3_graph(g, base; integer_encoding=:fibonacci,
            ref_window_size=spec.cs.window, lr_split=spec.cs.lr, cost_model=COST_MODEL)
    else # cg, K=1
        T = eltype(g)
        write_cg_mgs3_graph(g, base, [collect(T(1):T(nv(g)))];
            params=build_cg_params(spec.cg), integer_encoding=:fibonacci)
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
    println("  ord_ablation — $(spec.name)   seeds=$(SEEDS) cost_model=$COST_MODEL")
    println("=" ^ 64)
    g0 = load_core(spec)
    m  = ne(g0)
    tmpdir = mktempdir()

    # results[ordering][encoder] = Vector over seeds (length 1 for "original")
    results = Dict(o => Dict(e => Float64[] for e in ("bg","cs","cg")) for o in ORDERINGS)
    for ordering in ORDERINGS
        seeds = ordering == "original" ? [SEEDS[1]] : SEEDS
        for seed in seeds
            g = reorder(g0, ordering, seed)
            for enc in ("bg","cs","cg")
                b = encode_bpe(g, enc, spec, tmpdir, m)
                push!(results[ordering][enc], b)
                @info @sprintf("  [%-10s seed=%d] %s = %.4f BPE", ordering, seed, uppercase(enc), b)
            end
            ordering == "original" || (g = nothing; GC.gc())
        end
    end
    rm(tmpdir; recursive=true)

    fmt(v) = length(v) > 1 ? @sprintf("%.4f±%.4f", mean(v), std(v)) : @sprintf("%.4f", v[1])
    println("\n  $(spec.name) — BPE (Fibonacci backend), m=$m edges")
    @printf("  %-14s %14s %14s %14s\n", "Ordering", "BG", "CS", "CG")
    for (o, label) in (("original","Orig."), ("llp","LLP"), ("leiden_llp","Leiden+LLP"))
        @printf("  %-14s %14s %14s %14s\n", label,
                fmt(results[o]["bg"]), fmt(results[o]["cs"]), fmt(results[o]["cg"]))
    end
    println("=" ^ 64, "\n")
    return results
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
