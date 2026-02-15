#
# Adjacently: Julia Complex Directed Networks Library
# Copyright (C) 2016-2025 Jimmy Dubuisson <jimmy.dubuisson@gmail.com>
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

# Train a GCN vertex ordering model on CNR-2000 and evaluate compression.
#
# Phase 1: Proxy loss (bandwidth minimization) — fast, ~100 epochs
# Phase 2: REINFORCE with actual bpe — slow, ~50 epochs (full compression per epoch)
#
# Usage: julia scripts/train_gcn_ordering.jl [policy_filepath]

using Pkg
Pkg.activate(normpath(joinpath(@__DIR__, "..")))

using LightGraphs: nv, ne, outneighbors, vertices
using Adjacently
using Adjacently.IO: load_adjacency_list_from_csv
using Adjacently.Relabeling: relabel_graph
using Adjacently.GCN: GCNModel, gcn_ordering, train_gcn_proxy!, train_gcn_reinforce!,
                       TrainConfig, save_gcn_model, load_gcn_weights!
using Adjacently.MGS: write_rl_compressed_mgs3_graph, load_compressed_mgs3_graph

const PROJECT_ROOT = normpath(joinpath(@__DIR__, ".."))
const CNR_CSV = joinpath(PROJECT_ROOT, "datasets", "webgraph", "cnr-2000", "cnr-2000.csv")
const DEFAULT_POLICY = joinpath(PROJECT_ROOT, "policies", "cnr2000_policy.qpolicy")
const MODEL_OUTPUT = joinpath(PROJECT_ROOT, "tmp", "gcn_ordering_cnr2000.bin")

function make_compress_fn(g, policy_path, policy_id)
    return function(ordering::Dict{T,T}) where {T<:Unsigned}
        g_reordered = relabel_graph(g, ordering)
        output_base = joinpath(PROJECT_ROOT, "tmp", "gcn_eval_tmp")
        mkpath(dirname(output_base))
        write_rl_compressed_mgs3_graph(g_reordered, output_base, policy_path, policy_id;
            coding_scheme=:children, ref_window_size=7, integer_encoding=:fibonacci)
        mgz_file = output_base * ".mgz"
        file_size = filesize(mgz_file)
        bpe = 8.0 * file_size / ne(g_reordered)
        rm(mgz_file; force=true)
        return bpe
    end
end

function main()
    policy_path = length(ARGS) >= 1 ? ARGS[1] : nothing
    if policy_path == "greedy"
        policy_path = nothing
    elseif policy_path === nothing && isfile(DEFAULT_POLICY)
        policy_path = DEFAULT_POLICY
    end

    println("=" ^ 70)
    println("GCN Vertex Ordering Training — CNR-2000")
    println("=" ^ 70)

    # 1. Load graph
    println("\nLoading CNR-2000...")
    if !isfile(CNR_CSV)
        error("CNR-2000 CSV not found at $CNR_CSV")
    end
    g = load_adjacency_list_from_csv(CNR_CSV, ',', true)
    n, m = nv(g), ne(g)
    println("  Vertices: $n, Edges: $m")

    # 2. Create GCN model
    println("\nCreating GCN model (hidden_dim=32)...")
    model = GCNModel(g; hidden_dim=32, seed=42)
    println("  Model: $(model.input_dim) → $(model.hidden_dim) → 1")
    println("  Parameters: $(model.input_dim * model.hidden_dim + model.hidden_dim) weights")

    # 3. Phase 1: Proxy loss training
    config = TrainConfig(proxy_epochs=100, proxy_lr=0.001,
                         reinforce_epochs=50, reinforce_lr=0.0001,
                         sigma=0.1, baseline_ema=0.9)

    println("\n" * "-" ^ 70)
    println("Phase 1: Proxy Loss Training (bandwidth minimization)")
    println("-" ^ 70)
    proxy_losses = train_gcn_proxy!(model, g, config)
    println("  Initial loss: $(round(proxy_losses[1], digits=2))")
    println("  Final loss:   $(round(proxy_losses[end], digits=2))")

    # Evaluate ordering after Phase 1
    V = typeof(first(vertices(g)))
    ordering = gcn_ordering(model, V)
    compress_fn = make_compress_fn(g, policy_path, 1)
    bpe_phase1 = compress_fn(ordering)
    println("  BPE after Phase 1: $(round(bpe_phase1, digits=4))")

    # 4. Phase 2: REINFORCE training
    println("\n" * "-" ^ 70)
    println("Phase 2: REINFORCE Training (actual compression bpe)")
    println("-" ^ 70)
    reinforce_results = train_gcn_reinforce!(model, g, compress_fn, config)
    if !isempty(reinforce_results)
        best_epoch, best_bpe = reinforce_results[argmin([r[2] for r in reinforce_results])]
        println("  Best BPE: $(round(best_bpe, digits=4)) at epoch $best_epoch")
    end

    # 5. Save model
    mkpath(dirname(MODEL_OUTPUT))
    save_gcn_model(model, MODEL_OUTPUT)
    println("\nModel saved to: $MODEL_OUTPUT")

    # 6. Final evaluation
    println("\n" * "-" ^ 70)
    println("Final Evaluation")
    println("-" ^ 70)
    ordering_final = gcn_ordering(model, V)
    g_reordered = relabel_graph(g, ordering_final)
    output_base = joinpath(PROJECT_ROOT, "tmp", "cnr2000_gcn_final")
    write_rl_compressed_mgs3_graph(g_reordered, output_base, policy_path, 1;
        coding_scheme=:children, ref_window_size=7, integer_encoding=:fibonacci)
    mgz_file = output_base * ".mgz"
    file_size = filesize(mgz_file)
    final_bpe = 8.0 * file_size / ne(g_reordered)
    println("  File size: $(file_size) bytes")
    println("  Final BPE: $(round(final_bpe, digits=4))")

    # Roundtrip verification
    println("\nVerifying roundtrip...")
    g_loaded = load_compressed_mgs3_graph(mgz_file)
    n2, m2 = nv(g_loaded), ne(g_loaded)
    if m2 == ne(g_reordered)
        println("  Roundtrip: OK ($m2 edges match)")
    else
        println("  Roundtrip: MISMATCH (expected $(ne(g_reordered)), got $m2)")
    end
    rm(mgz_file; force=true)

    println("\n" * "=" ^ 70)
end

main()
