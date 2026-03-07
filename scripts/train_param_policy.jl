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

# Train supervised GNN-based CGE parameter predictor.
#
# Supports both fresh training and continuation from saved weights.
# Uses synthetic graphs + subgraphs sampled from real web graphs
# (CNR-2000, IN-2004) when available.
#
# Usage: julia scripts/train_param_policy.jl [epochs] [batch_size] [hidden_dim]

using Pkg
Pkg.activate(normpath(joinpath(@__DIR__, "..")))

using LightGraphs: AbstractGraph, nv, ne, vertices, outneighbors, add_edge!, eltype
using Adjacently
using Adjacently.CustomLightGraphs: SimpleDiGraph
using Adjacently.GraphGenerators: generate_training_batch
using Adjacently.GNN: GNNModel, gnn_forward
using Adjacently.ParamPrediction: ParamPredictor, MasterGNNWeights,
                                   train_param_predictor!, evaluate_param_predictor,
                                   save_param_predictor, load_param_predictor,
                                   predict_params, graph_level_readout,
                                   compute_bpe, grid_search_best_params,
                                   generate_training_labels
using Adjacently.IO: load_adjacency_list_from_csv
using Adjacently.Graph: get_core
using Adjacently.Compression: CGE
using Statistics: mean
using Random: MersenneTwister, shuffle!

const PROJECT_ROOT = normpath(joinpath(@__DIR__, ".."))
const CNR_CSV = joinpath(PROJECT_ROOT, "datasets", "webgraph", "cnr-2000", "cnr-2000.csv")
const IN2004_CSV = joinpath(PROJECT_ROOT, "datasets", "webgraph", "in-2004", "in-2004.csv")
const MODEL_PATH = joinpath(PROJECT_ROOT, "models", "param_predictor_trained.bson")

"""
    sample_bfs_subgraph(g, target_n; seed=42)

Sample a connected subgraph of approximately `target_n` vertices from `g`
using BFS from a random start vertex. Returns a `SimpleDiGraph{UInt32}` with
vertices renumbered 1:n_sub.
"""
function sample_bfs_subgraph(g::AbstractGraph{T}, target_n::Int;
                              seed::Int=42) where {T<:Unsigned}
    rng = MersenneTwister(seed)
    n = Int(nv(g))
    target_n = min(target_n, n)

    # Pick random start vertex
    start = T(rand(rng, 1:n))

    # BFS to collect target_n vertices
    visited = Set{T}()
    queue = T[start]
    push!(visited, start)
    order = T[]

    while !isempty(queue) && length(order) < target_n
        v = popfirst!(queue)
        push!(order, v)
        for u in outneighbors(g, v)
            if !(u in visited) && length(visited) < target_n
                push!(visited, u)
                push!(queue, u)
            end
        end
    end

    # Build renumbered subgraph
    n_sub = length(order)
    old_to_new = Dict{T, UInt32}()
    for (i, v) in enumerate(order)
        old_to_new[v] = UInt32(i)
    end

    g_sub = SimpleDiGraph{UInt32}(UInt32(n_sub))
    for v in order
        new_v = old_to_new[v]
        for u in outneighbors(g, v)
            if haskey(old_to_new, u)
                add_edge!(g_sub, new_v, old_to_new[u])
            end
        end
    end

    return g_sub
end

function main()
    epochs = length(ARGS) >= 1 ? parse(Int, ARGS[1]) : 50
    batch_size = length(ARGS) >= 2 ? parse(Int, ARGS[2]) : 16
    hidden_dim = length(ARGS) >= 3 ? parse(Int, ARGS[3]) : 32
    feature_set = :enhanced

    # Check if we can resume from saved model
    resume = isfile(MODEL_PATH)

    println("=" ^ 70)
    println("Supervised GNN-based CGE Parameter Predictor Training")
    println("  epochs=$epochs, batch_size=$batch_size, hidden_dim=$hidden_dim")
    println("  resume=$(resume ? "yes ($(MODEL_PATH))" : "no (fresh)")")
    println("=" ^ 70)

    # =====================================================================
    # 1. Generate synthetic training graphs
    # =====================================================================
    println("\n[1/7] Generating $batch_size synthetic training graphs...")
    t0 = time()
    synthetic_graphs = generate_training_batch(; n_range=(500, 2000),
                                                 batch_size=batch_size, seed=42)
    gen_time = time() - t0
    for (i, g) in enumerate(synthetic_graphs)
        println("  Synthetic $i: $(nv(g))v, $(ne(g))e")
    end
    println("  Generation time: $(round(gen_time, digits=2))s")

    # =====================================================================
    # 2. Sample subgraphs from real web graphs
    # =====================================================================
    println("\n[2/7] Sampling subgraphs from real web graphs...")
    real_subgraphs = SimpleDiGraph{UInt32}[]
    n_samples_per_graph = 4
    subgraph_sizes = [1000, 2000, 3000, 5000]

    for (name, csv_path, seed_base) in [("CNR-2000", CNR_CSV, 100), ("IN-2004", IN2004_CSV, 200)]
        if !isfile(csv_path)
            println("  $name not found, skipping")
            continue
        end
        println("  Loading $name...")
        t1 = time()
        g_full = load_adjacency_list_from_csv(csv_path, ',', true)
        println("    Full graph: $(nv(g_full))v, $(ne(g_full))e ($(round(time()-t1, digits=1))s)")

        # Extract main SCC core so BFS starts from well-connected vertices
        println("    Extracting main SCC core...")
        t_core = time()
        g_core, _, _ = get_core(g_full)
        println("    Core: $(nv(g_core))v, $(ne(g_core))e ($(round(time()-t_core, digits=1))s)")

        for (i, sz) in enumerate(subgraph_sizes)
            sub = sample_bfs_subgraph(g_core, sz; seed=seed_base + i)
            push!(real_subgraphs, sub)
            println("    $name subgraph $i: $(nv(sub))v, $(ne(sub))e (target=$sz)")
        end
    end

    println("  Total real subgraphs: $(length(real_subgraphs))")

    # =====================================================================
    # 3. Combine all training graphs
    # =====================================================================
    all_graphs = SimpleDiGraph{UInt32}[synthetic_graphs; real_subgraphs]
    println("\n[3/7] Total training graphs: $(length(all_graphs)) ($(length(synthetic_graphs)) synthetic + $(length(real_subgraphs)) real)")

    # =====================================================================
    # 4. Grid-search labels for all graphs
    # =====================================================================
    println("\n[4/7] Grid-searching optimal params (512 combos each)...")
    t2 = time()
    labels = generate_training_labels(all_graphs; verbose=true)
    label_time = time() - t2
    println("  Label generation time: $(round(label_time, digits=1))s")

    # =====================================================================
    # 5. Initialize or load model weights
    # =====================================================================
    if resume
        println("\n[5/7] Loading existing model weights from $(MODEL_PATH)...")
        predictor, master = load_param_predictor(MODEL_PATH)
        println("  Loaded: hidden_dim=$(master.hidden_dim), feat_dim=$(predictor.feat_dim)")
        println("  Heads: $(predictor.head_names)")
    else
        println("\n[5/7] Initializing fresh model weights...")
        ref_gnn = GNNModel(all_graphs[1]; hidden_dim=hidden_dim, seed=42,
                            feature_set=feature_set, use_gat=false)
        master = MasterGNNWeights(ref_gnn)
        input_dim = ref_gnn.input_dim
        println("  Architecture: GCN $(input_dim) → $(hidden_dim) → 1")

        feat_dim = hidden_dim
        predictor = ParamPredictor(feat_dim; lr=0.01, seed=42)
    end
    println("  Feature set: $feature_set")
    for (i, name) in enumerate(predictor.head_names)
        println("    $name: $(predictor.head_options[i])")
    end

    # =====================================================================
    # 6. Train supervised
    # =====================================================================
    println("\n[6/7] Training ($epochs epochs, cross-entropy, $(length(all_graphs)) graphs)...")
    t3 = time()
    loss_history = train_param_predictor!(master, predictor, all_graphs, labels;
                                           epochs=epochs,
                                           verbose=true,
                                           feature_set=feature_set,
                                           gnn_lr=0.0001,
                                           seed=resume ? 1000 : 42)
    train_time = time() - t3
    println("  Training time: $(round(train_time, digits=1))s")

    first_quarter = loss_history[1:max(1, length(loss_history) ÷ 4)]
    last_quarter = loss_history[max(1, end - length(loss_history) ÷ 4 + 1):end]
    println("  Loss first quarter avg: $(round(mean(first_quarter), digits=4))")
    println("  Loss last quarter avg:  $(round(mean(last_quarter), digits=4))")

    # =====================================================================
    # 7. Save model
    # =====================================================================
    mkpath(dirname(MODEL_PATH))
    save_param_predictor(predictor, master, MODEL_PATH)
    println("\n  Model saved to: $MODEL_PATH")

    # =====================================================================
    # Evaluate on full CNR-2000 and IN-2004
    # =====================================================================
    for (name, csv_path) in [("CNR-2000", CNR_CSV), ("IN-2004", IN2004_CSV)]
        if !isfile(csv_path)
            println("\n[Eval] $name CSV not found — skipping")
            continue
        end
        println("\n[Eval] Evaluating on $name...")
        t_eval = time()
        g_eval = load_adjacency_list_from_csv(csv_path, ',', true)
        println("  Loaded: $(nv(g_eval))v, $(ne(g_eval))e")

        bpe_pred, params_pred = evaluate_param_predictor(master, predictor, g_eval;
                                                          feature_set=feature_set, seed=42)
        println("  Predicted BPE (K=1): $(round(bpe_pred, digits=4))")
        println("  Predicted params:")
        println("    window=$(params_pred.intra_ref_window)")
        println("    intervals=$(params_pred.intra_intervals)")
        println("    lr_split=$(params_pred.intra_lr_split)")
        println("    mil=$(params_pred.intra_mil)")
        println("    copy_adaptive=$(params_pred.intra_copy_adaptive)")
        println("    stop_deltas=$(params_pred.intra_stop_deltas)")
        println("    zigzag=$(params_pred.intra_zigzag)")

        # Hand-tuned baseline
        baseline_params = CGE.CGEParams(;
            intra_ref_window=64, intra_intervals=false, intra_lr_split=false,
            intra_mil=4, intra_copy_adaptive=true, intra_stop_deltas=true,
            intra_zigzag=true, intra_copy_blocks=true, intra_ref_fixwidth=true,
            intra_add_adaptive=true, intra_raw_adaptive=true,
            varint=:fibonacci, count_varint=:fibonacci, gap=:fibonacci,
            degree=:elias_delta, membership=:implicit_ranges,
        )
        T_g = eltype(g_eval)
        n_g = Int(nv(g_eval))
        clusters_g = [collect(T_g(1):T_g(n_g))]
        bpe_baseline = compute_bpe(g_eval, baseline_params, clusters_g)
        println("  Hand-tuned BPE (K=1): $(round(bpe_baseline, digits=4))")
        diff = bpe_pred - bpe_baseline
        println("  Difference: $(round(diff, digits=4)) ($(diff > 0 ? "worse" : "better"))")
        println("  Eval time: $(round(time()-t_eval, digits=1))s")
    end

    println("\n" * "=" ^ 70)
    println("Training Summary")
    println("=" ^ 70)
    println("  Resumed:        $resume")
    println("  Epochs:         $epochs")
    println("  Graphs:         $(length(all_graphs)) ($(length(synthetic_graphs)) synthetic + $(length(real_subgraphs)) real)")
    println("  Hidden dim:     $(master.hidden_dim)")
    println("  Feature set:    $feature_set")
    println("  Label time:     $(round(label_time, digits=1))s")
    println("  Training time:  $(round(train_time, digits=1))s")
    println("  Model saved:    $MODEL_PATH")
    println("=" ^ 70)
end

main()
