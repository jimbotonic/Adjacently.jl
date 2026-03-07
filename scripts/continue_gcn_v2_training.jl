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

#!/usr/bin/env julia

# Load GCN v2 model (10-feature, hidden=64) and continue training.
#
# Phase 1: 200 more proxy epochs with reduced LR (0.0005) to refine
# Phase 2: 50 REINFORCE epochs using actual RL-greedy BPE as reward
#
# Usage: julia scripts/continue_gcn_v2_training.jl

using Pkg
Pkg.activate(normpath(joinpath(@__DIR__, "..")))

using LightGraphs: nv, ne, outneighbors, vertices, eltype
using Adjacently
using Adjacently.IO: load_adjacency_list_from_csv
using Adjacently.Relabeling: relabel_graph
using Adjacently.GCN: GCNModel, gcn_forward, gcn_ordering, gcn_backward!,
                       train_gcn_proxy!, train_gcn_reinforce!,
                       TrainConfig, save_gcn_model, load_gcn_weights!
using Adjacently.MGS: write_rl_compressed_mgs3_graph, load_compressed_mgs3_graph

const PROJECT_ROOT = normpath(joinpath(@__DIR__, ".."))
const CNR_CSV = joinpath(PROJECT_ROOT, "datasets", "webgraph", "cnr-2000", "cnr-2000.csv")
const V2_WEIGHTS = joinpath(PROJECT_ROOT, "policies", "cnr2000_gcn_v2_ordering.bin")
const OUTPUT_DIR = joinpath(PROJECT_ROOT, "datasets", "webgraph", "cnr-2000")
const MODEL_OUTPUT = joinpath(PROJECT_ROOT, "policies", "cnr2000_gcn_v2_ordering.bin")

function make_compress_fn(g, m)
    return function(ordering::Dict{T,T}) where {T<:Unsigned}
        g_reordered = relabel_graph(g, ordering)
        output_base = joinpath(PROJECT_ROOT, "tmp", "gcn_v2_reinforce_eval")
        mkpath(dirname(output_base))
        write_rl_compressed_mgs3_graph(g_reordered, output_base, nothing, 1;
            coding_scheme=:children, ref_window_size=7, integer_encoding=:fibonacci)
        mgz_file = output_base * ".mgz"
        file_size = filesize(mgz_file)
        bpe = 8.0 * file_size / m
        rm(mgz_file; force=true)
        return bpe
    end
end

function evaluate_bpe(model, g, m)
    V = eltype(g)
    ordering = gcn_ordering(model, V)
    g_reordered = relabel_graph(g, ordering)
    output_base = joinpath(PROJECT_ROOT, "tmp", "gcn_v2_eval_tmp")
    mkpath(dirname(output_base))
    write_rl_compressed_mgs3_graph(g_reordered, output_base, nothing, 1;
        coding_scheme=:children, ref_window_size=7, integer_encoding=:fibonacci)
    mgz_file = output_base * ".mgz"
    file_size = filesize(mgz_file)
    bpe = 8.0 * file_size / m

    # Roundtrip check
    g_loaded = load_compressed_mgs3_graph(mgz_file)
    m2 = ne(g_loaded)
    rm(mgz_file; force=true)

    return bpe, m2 == ne(g_reordered)
end

function main()
    println("=" ^ 70)
    println("Continue Training GCN v2 — CNR-2000")
    println("=" ^ 70)

    # =====================================================================
    # 1. Load graph
    # =====================================================================
    println("\n[1/6] Loading CNR-2000...")
    if !isfile(CNR_CSV)
        error("CNR-2000 CSV not found at $CNR_CSV")
    end
    t0 = time()
    g = load_adjacency_list_from_csv(CNR_CSV, ',', true)
    n, m = nv(g), ne(g)
    println("  Vertices: $n, Edges: $m")
    println("  Load time: $(round(time() - t0, digits=2))s")

    # =====================================================================
    # 2. Build model with same architecture and load saved weights
    # =====================================================================
    println("\n[2/6] Building GCN v2 model and loading saved weights...")
    if !isfile(V2_WEIGHTS)
        error("GCN v2 weights not found at $V2_WEIGHTS — run retrain_gcn_test_cnr2000.jl or compare_cnr2000_bpe.jl first")
    end

    t1 = time()
    model = GCNModel(g; hidden_dim=64, seed=42)
    println("  Model built: $(model.input_dim) -> $(model.hidden_dim) -> 1")
    println("  Feature extraction: $(round(time() - t1, digits=2))s")

    load_gcn_weights!(model, V2_WEIGHTS)
    println("  Weights loaded from: $V2_WEIGHTS")

    # Evaluate starting point
    println("\n  Evaluating loaded model...")
    bpe_start, rt_ok = evaluate_bpe(model, g, m)
    println("  Starting BPE: $(round(bpe_start, digits=4))  Roundtrip: $(rt_ok ? "OK" : "FAIL")")

    # =====================================================================
    # 3. Phase 1: Continue proxy training (200 epochs, reduced LR)
    # =====================================================================
    println("\n[3/6] Phase 1: Proxy loss training (200 epochs, lr=0.0005)...")
    config_proxy = TrainConfig(proxy_epochs=200, proxy_lr=0.0005,
                                reinforce_epochs=0, reinforce_lr=0.0,
                                sigma=0.1, baseline_ema=0.9)

    t2 = time()
    proxy_losses = train_gcn_proxy!(model, g, config_proxy)
    proxy_time = time() - t2
    println("  Initial loss: $(round(proxy_losses[1], digits=2))")
    println("  Final loss:   $(round(proxy_losses[end], digits=2))")
    println("  Training time: $(round(proxy_time, digits=2))s")

    # Evaluate after proxy phase
    println("\n  Evaluating after proxy phase...")
    bpe_proxy, rt_ok = evaluate_bpe(model, g, m)
    println("  BPE after proxy: $(round(bpe_proxy, digits=4))  Roundtrip: $(rt_ok ? "OK" : "FAIL")")

    # Save checkpoint
    save_gcn_model(model, MODEL_OUTPUT)
    println("  Checkpoint saved to: $MODEL_OUTPUT")

    # =====================================================================
    # 4. Phase 2: REINFORCE training (50 epochs, actual BPE as reward)
    # =====================================================================
    println("\n[4/6] Phase 2: REINFORCE training (50 epochs, lr=0.0001)...")
    config_reinforce = TrainConfig(proxy_epochs=0, proxy_lr=0.0,
                                    reinforce_epochs=50, reinforce_lr=0.0001,
                                    sigma=0.1, baseline_ema=0.9)

    compress_fn = make_compress_fn(g, m)
    t3 = time()
    reinforce_results = train_gcn_reinforce!(model, g, compress_fn, config_reinforce)
    reinforce_time = time() - t3

    if !isempty(reinforce_results)
        bpes = [r[2] for r in reinforce_results]
        best_idx = argmin(bpes)
        best_epoch, best_bpe = reinforce_results[best_idx]
        println("  Best BPE: $(round(best_bpe, digits=4)) at epoch $best_epoch")
        println("  Last BPE: $(round(bpes[end], digits=4))")
    end
    println("  REINFORCE time: $(round(reinforce_time, digits=2))s")

    # =====================================================================
    # 5. Save final model
    # =====================================================================
    println("\n[5/6] Saving final model...")
    save_gcn_model(model, MODEL_OUTPUT)
    println("  Saved to: $MODEL_OUTPUT")

    # =====================================================================
    # 6. Final evaluation with roundtrip
    # =====================================================================
    println("\n[6/6] Final evaluation...")
    V = eltype(g)
    ordering = gcn_ordering(model, V)
    g_reordered = relabel_graph(g, ordering)

    output_final = joinpath(OUTPUT_DIR, "cnr2000_rl_greedy_gcn_v2")
    mkpath(OUTPUT_DIR)
    t4 = time()
    write_rl_compressed_mgs3_graph(g_reordered, output_final, nothing, 1;
        coding_scheme=:children, ref_window_size=7, integer_encoding=:fibonacci)
    compress_time = time() - t4

    mgz_final = output_final * ".mgz"
    file_size = filesize(mgz_final)
    bpe_final = 8.0 * file_size / m

    g_loaded = load_compressed_mgs3_graph(mgz_final)
    m2 = ne(g_loaded)
    rt_ok = m2 == ne(g_reordered)

    println("  File size:     $file_size bytes ($(round(file_size / 1024, digits=1)) KB)")
    println("  Final BPE:     $(round(bpe_final, digits=4))")
    println("  Compress time: $(round(compress_time, digits=2))s")
    println("  Roundtrip:     $(rt_ok ? "OK ($m2 edges)" : "MISMATCH (expected $(ne(g_reordered)), got $m2)")")

    # =====================================================================
    # Summary
    # =====================================================================
    println("\n" * "=" ^ 70)
    println("Summary")
    println("=" ^ 70)
    println("  Graph:            CNR-2000 ($n vertices, $m edges)")
    println("  Features:         $(model.input_dim)-dim enhanced")
    println("  Hidden dim:       $(model.hidden_dim)")
    println("  Starting BPE:     $(round(bpe_start, digits=4)) (loaded weights)")
    println("  After proxy:      $(round(bpe_proxy, digits=4)) (+200 epochs, lr=0.0005)")
    println("  Final BPE:        $(round(bpe_final, digits=4)) (+50 REINFORCE epochs)")
    if !isempty(reinforce_results)
        bpes = [r[2] for r in reinforce_results]
        println("  Best REINFORCE:   $(round(minimum(bpes), digits=4))")
    end
    println("  Model saved:      $MODEL_OUTPUT")
    println("  Output .mgz:      $mgz_final")
    println("=" ^ 70)
end

main()
