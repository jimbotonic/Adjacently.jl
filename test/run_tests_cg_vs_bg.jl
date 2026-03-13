#!/usr/bin/env julia
#
# Compare CG K=1 vs BG best on CNR-2000
# Tests whether the BPE gap is from clustering (K=2) or encoder differences.
#
# Usage:
#   julia --project test/run_tests_cge_vs_std.jl [DATASET]

include("run_tests_main.jl")

using Logging
global_logger(ConsoleLogger(stderr, Logging.Info))

using LightGraphs: nv, ne, outneighbors, vertices, eltype
using Adjacently
using Adjacently.IO: load_adjacency_list_from_csv
using Adjacently.Relabeling: relabel_graph, relabel_vertices_llp, load_llp_ordering
using Adjacently.MGS: write_bg_mgs3_graph, write_cg_mgs3_graph, load_compressed_mgs3_graph
using Adjacently.Compression.CG: CGParams

const DATASET = length(ARGS) >= 1 ? ARGS[1] : "cnr-2000"
const PREFIX  = replace(DATASET, "-" => "")
const DS_DIR  = normpath(joinpath(@__DIR__, "..", "datasets", "webgraph", DATASET))
const DS_CSV  = joinpath(DS_DIR, DATASET * ".csv")
const OUT_DIR = DS_DIR

@info "Dataset: $DATASET"
@info "  CSV: $DS_CSV"

# ── Load graph ───────────────────────────────────────────────────────────────
@info "Loading dataset..."
g = load_adjacency_list_from_csv(DS_CSV, ',', true)
@test nv(g) > 0
@test ne(g) > 0
@info "  $(nv(g)) vertices, $(ne(g)) edges"

# Apply LLP reordering if available
llp_file = joinpath(DS_DIR, "$(PREFIX)_llp_p10.bin")
if isfile(llp_file)
    @info "Applying LLP ordering from $llp_file"
    T_v = eltype(g)
    llp_order = load_llp_ordering(llp_file, T_v)
    g = relabel_graph(g, llp_order)
end

total_edges = ne(g)
n = nv(g)
T_v = eltype(g)

function compute_bpe(filename)
    fsize = filesize(filename)
    header_bytes = 12
    data_bytes = fsize - header_bytes
    return data_bytes * 8.0 / total_edges
end

function verify_roundtrip(original_g, loaded_g, label)
    @test nv(loaded_g) == nv(original_g)
    @test ne(loaded_g) == ne(original_g)
    mismatch = 0
    for v in vertices(original_g)
        orig_nbs = sort(collect(outneighbors(original_g, v)))
        load_nbs = sort(collect(outneighbors(loaded_g, v)))
        if orig_nbs != load_nbs
            mismatch += 1
            if mismatch <= 3
                @warn "  Edge mismatch at vertex $v ($label)"
            end
        end
    end
    @test mismatch == 0
    if mismatch > 0
        @error "  $mismatch vertices with mismatched edges ($label)"
    else
        @info "  Roundtrip verified: all edges match ($label)"
    end
end

# ── CG K=1 params (best known settings, but K=1 = single cluster) ──────────
cg_k1_params = CGParams(
    L=128,
    varint=:fibonacci, count_varint=:fibonacci,
    gap=:fibonacci, degree=:elias_delta,
    undirected_pairs=false,
    perm_strategy=:blockpos, membership=:implicit_ranges,
    inter_strategy=:lists,
    intra_ref_enabled=true, intra_ref_window=64,
    intra_ref_rle=false,
    intra_block_try=false,
    positions_mode=:delta, additions_mode=:delta,
    min_cluster_density=0.0,
    intra_intervals=false, intra_mil=4, intra_greedy_mil=false,
    intra_zigzag=true, intra_stop_deltas=true,
    intra_copy_blocks=true, intra_copy_adaptive=true,
    intra_ref_fixwidth=true, intra_ref_vlc=false,
    intra_add_adaptive=true, intra_raw_adaptive=true,
    intra_adapt_mil=4,
)

# K=1: single cluster containing all vertices
clusters_k1 = [collect(T_v(1):T_v(n))]

@testset "CG K=1 vs BG — $DATASET" begin
    results = []

    # ── BG best ─────────────────────────────────────────────────────────────
    @info "───────────────────────────────────────────────"
    @info "Testing: BG best (lr_split + multi_ref)"
    out_base = joinpath(OUT_DIR, "$(PREFIX)_cmp_bg")
    t_write = @elapsed write_bg_mgs3_graph(g, out_base;
        ref_window_size=64, copy_blocks=true, stop_deltas=true,
        lr_split=true, multi_ref=true)
    mgz_file = out_base * ".mgz"
    bpe = compute_bpe(mgz_file)
    @info "  Written in $(round(t_write, digits=2))s — BPE = $(round(bpe, digits=4))"
    t_load = @elapsed loaded_g = load_compressed_mgs3_graph(mgz_file)
    @info "  Loaded in $(round(t_load, digits=2))s"
    verify_roundtrip(g, loaded_g, "BG best")
    push!(results, (label="BG best (lr+mr)", bpe=bpe, write_time=t_write, load_time=t_load))
    rm(mgz_file; force=true)

    # ── CG K=1 ──────────────────────────────────────────────────────────────
    @info "───────────────────────────────────────────────"
    @info "Testing: CG K=1 (single cluster, w=64)"
    out_base = joinpath(OUT_DIR, "$(PREFIX)_cmp_cg_k1")
    t_write = @elapsed write_cg_mgs3_graph(g, out_base, clusters_k1;
        coding_scheme=:children, integer_encoding=:fibonacci,
        params=cg_k1_params)
    mgz_file = out_base * ".mgz"
    bpe = compute_bpe(mgz_file)
    @info "  Written in $(round(t_write, digits=2))s — BPE = $(round(bpe, digits=4))"
    t_load = @elapsed loaded_g = load_compressed_mgs3_graph(mgz_file)
    @info "  Loaded in $(round(t_load, digits=2))s"
    verify_roundtrip(g, loaded_g, "CG K=1")
    push!(results, (label="CG K=1 (w=64)", bpe=bpe, write_time=t_write, load_time=t_load))
    rm(mgz_file; force=true)

    # ── CG K=1 with intervals + lr_split ────────────────────────────────────
    @info "───────────────────────────────────────────────"
    @info "Testing: CG K=1 (intervals + lr_split)"
    cg_k1_iv = CGParams(
        L=128,
        varint=:fibonacci, count_varint=:fibonacci,
        gap=:fibonacci, degree=:elias_delta,
        undirected_pairs=false,
        perm_strategy=:blockpos, membership=:implicit_ranges,
        inter_strategy=:lists,
        intra_ref_enabled=true, intra_ref_window=64,
        intra_ref_rle=false,
        intra_block_try=false,
        positions_mode=:delta, additions_mode=:delta,
        min_cluster_density=0.0,
        intra_intervals=true, intra_mil=4, intra_greedy_mil=false,
        intra_zigzag=true, intra_stop_deltas=true,
        intra_copy_blocks=true, intra_copy_adaptive=true,
        intra_ref_fixwidth=true, intra_ref_vlc=false,
        intra_add_adaptive=true, intra_raw_adaptive=true,
        intra_adapt_mil=4,
        intra_lr_split=true,
    )
    out_base = joinpath(OUT_DIR, "$(PREFIX)_cmp_cg_k1_iv")
    t_write = @elapsed write_cg_mgs3_graph(g, out_base, clusters_k1;
        coding_scheme=:children, integer_encoding=:fibonacci,
        params=cg_k1_iv)
    mgz_file = out_base * ".mgz"
    bpe = compute_bpe(mgz_file)
    @info "  Written in $(round(t_write, digits=2))s — BPE = $(round(bpe, digits=4))"
    t_load = @elapsed loaded_g = load_compressed_mgs3_graph(mgz_file)
    @info "  Loaded in $(round(t_load, digits=2))s"
    verify_roundtrip(g, loaded_g, "CG K=1 iv+lr")
    push!(results, (label="CG K=1 (iv+lr)", bpe=bpe, write_time=t_write, load_time=t_load))
    rm(mgz_file; force=true)

    # ── Summary ──────────────────────────────────────────────────────────────
    @info "\n═══════════════════════════════════════════════"
    @info "Results Summary ($DATASET, $(nv(g)) vertices, $total_edges edges):"
    @info "───────────────────────────────────────────────"
    for r in results
        @info "  $(rpad(r.label, 25)) BPE=$(round(r.bpe, digits=4))  write=$(round(r.write_time, digits=1))s  load=$(round(r.load_time, digits=1))s"
    end
    bg_bpe = results[1].bpe
    for r in results[2:end]
        delta = r.bpe - bg_bpe
        @info "    → $(r.label): $(delta > 0 ? "+" : "")$(round(delta, digits=4)) BPE vs BG"
    end
end
