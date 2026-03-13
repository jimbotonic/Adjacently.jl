#!/usr/bin/env julia
#
# Test: BG compression with LR-split residuals
# Tests roundtrip correctness and measures BPE improvement.
#
# Usage:
#   julia --project test/run_tests_std_adaptive.jl [DATASET]
#   Default DATASET: cnr-2000

include("run_tests_main.jl")

using Logging
global_logger(ConsoleLogger(stderr, Logging.Info))

using LightGraphs: nv, ne, outneighbors, vertices, eltype
using Adjacently
using Adjacently.IO: load_adjacency_list_from_csv
using Adjacently.Relabeling: relabel_graph, relabel_vertices_llp, load_llp_ordering
using Adjacently.MGS: write_bg_mgs3_graph, load_compressed_mgs3_graph

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

# ── Configuration matrix ─────────────────────────────────────────────────────
configs = [
    (label="BG adaptive + lr_split",
     kwargs=Dict(:ref_window_size=>64, :copy_blocks=>true, :stop_deltas=>true,
                 :lr_split=>true, :adaptive_header=>true)),
    (label="BG adaptive + lr + exact",
     kwargs=Dict(:ref_window_size=>64, :copy_blocks=>true, :stop_deltas=>true,
                 :lr_split=>true, :adaptive_header=>true, :exact_costing=>true)),
]

# ── Run tests ────────────────────────────────────────────────────────────────
@testset "BG Compression — $DATASET" begin
    results = []
    for cfg in configs
        @info "───────────────────────────────────────────────"
        @info "Testing: $(cfg.label)"
        out_base = joinpath(OUT_DIR, "$(PREFIX)_bg_test")

        t_write = @elapsed write_bg_mgs3_graph(g, out_base; cfg.kwargs...)
        mgz_file = out_base * ".mgz"
        bpe = compute_bpe(mgz_file)

        @info "  Written in $(round(t_write, digits=2))s — BPE = $(round(bpe, digits=4))"

        t_load = @elapsed loaded_g = load_compressed_mgs3_graph(mgz_file)
        @info "  Loaded in $(round(t_load, digits=2))s"

        verify_roundtrip(g, loaded_g, cfg.label)
        push!(results, (label=cfg.label, bpe=bpe, write_time=t_write, load_time=t_load))

        rm(mgz_file; force=true)
    end

    @info "\n═══════════════════════════════════════════════"
    @info "Results Summary ($DATASET, $(nv(g)) vertices, $total_edges edges):"
    @info "───────────────────────────────────────────────"
    for r in results
        @info "  $(rpad(r.label, 30)) BPE=$(round(r.bpe, digits=4))  write=$(round(r.write_time, digits=1))s  load=$(round(r.load_time, digits=1))s"
    end
    if length(results) >= 2
        baseline = results[1].bpe
        for r in results[2:end]
            delta = r.bpe - baseline
            @info "    → $(r.label): $(delta > 0 ? "+" : "")$(round(delta, digits=4)) BPE vs baseline"
        end
    end
end
