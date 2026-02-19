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

# Train RL compression policy with GNN-augmented policy on CNR-2000.
#
# Pipeline:
#   1. Load CNR-2000
#   2. Build GNN features (2-layer GNN proxy pass)
#   3. Relabel graph with LLP (community-aware ordering)
#   4. Train RL policy with GNN features (actor-critic)
#   5. Compress with RL header tag (fibonacci encoding)
#   6. Attempt roundtrip
#
# Usage: julia scripts/train_rl_gcn_cnr2000.jl [rl_episodes] [gnn_proxy_epochs] [hidden_dim] [feature_set]

using Pkg
Pkg.activate(normpath(joinpath(@__DIR__, "..")))

using LightGraphs: nv, ne, outneighbors, vertices, eltype
using Adjacently
using Adjacently.IO: load_adjacency_list_from_csv
using Adjacently.Relabeling: relabel_graph, relabel_vertices_llp
using Adjacently.GNN: GNNModel, gnn_forward, TrainConfig as GNNTrainConfig, train_gnn_proxy!, gnn_backward!
using Adjacently.RL: CompressionEnv, QPolicy, save_policy,
                      action_from_index, best_action, update!, decay_epsilon!,
                      get_bits_per_edge, select_action,
                      current_vertex, reset!, step!, Action, NUM_ACTIONS
using StatsBase: fit, Histogram
using Statistics: mean

const PROJECT_ROOT = normpath(joinpath(@__DIR__, ".."))
const CNR_CSV = joinpath(PROJECT_ROOT, "datasets", "webgraph", "cnr-2000", "cnr-2000.csv")

const DEGREE_BINS = 6
const GNN_SCORE_BINS = 10
const GNN_NEIGHBOR_SCORE_BINS = 10
const NUM_Q_STATES = DEGREE_BINS * GNN_SCORE_BINS * GNN_NEIGHBOR_SCORE_BINS

function main()
    rl_episodes = length(ARGS) >= 1 ? parse(Int, ARGS[1]) : 100
    gnn_epochs = length(ARGS) >= 2 ? parse(Int, ARGS[2]) : 50
    hidden_dim = length(ARGS) >= 3 ? parse(Int, ARGS[3]) : 32
    feature_set = length(ARGS) >= 4 ? Symbol(ARGS[4]) : :enhanced

    println("=" ^ 70)
    println("RL (Q-Learning) + GCN Compression Pipeline — CNR-2000")
    println("  hidden_dim=$hidden_dim, feature_set=$feature_set")
    println("=" ^ 70)

    # =====================================================================
    # 1. Load graph
    # =====================================================================
    println("\n[1/6] Loading CNR-2000 dataset...")
    if !isfile(CNR_CSV)
        error("CNR-2000 CSV not found at $CNR_CSV")
    end
    t0 = time()
    g = load_adjacency_list_from_csv(CNR_CSV, ',', true)
    load_time = time() - t0
    n, m = nv(g), ne(g)
    println("  Vertices: $n")
    println("  Edges:    $m")
    println("  Load time: $(round(load_time, digits=2))s")

    # =====================================================================
    # 2. Train GNN proxy (to populate features)
    # =====================================================================
    println("\n[2/6] Training GNN proxy (feature learning, $gnn_epochs epochs)...")
    t1 = time()
    gnn = GNNModel(g; hidden_dim=hidden_dim, seed=42, feature_set=feature_set)
    println("  GNN: $(gnn.input_dim) → $(gnn.hidden_dim) → 1 (features=$feature_set)")
    resume_path = joinpath(PROJECT_ROOT, "policies", "cnr2000_gnn_model.bin")
    if isfile(resume_path)
        try
            Adjacently.GNN.load_gnn_weights!(gnn, resume_path)
            println("  Loaded existing GNN weights from: $(resume_path)")
        catch err
            println("  Warning: could not load GNN weights (", err, ") — starting fresh")
        end
    end
    gnn_cfg = GNNTrainConfig(proxy_epochs=gnn_epochs, proxy_lr=0.001,
                              reinforce_epochs=0, reinforce_lr=0.0,
                              sigma=0.0, baseline_ema=0.9, grad_clip_norm=1.0)
    proxy_losses = train_gnn_proxy!(gnn, g, gnn_cfg)
    gnn_time = time() - t1
    gnn_forward(gnn)  # populate caches H1/Z2
    println("  Proxy loss: start=$(round(proxy_losses[1], digits=2)), end=$(round(proxy_losses[end], digits=2))")
    println("  GNN proxy time: $(round(gnn_time, digits=2))s")

    out_dir = joinpath(PROJECT_ROOT, "policies")
    mkpath(out_dir)
    gnn_path = joinpath(out_dir, "cnr2000_gnn_model.bin")
    Adjacently.GNN.save_gnn_model(gnn, gnn_path)
    println("  GNN model: $gnn_path")

    # =====================================================================
    # 3. Apply LLP relabeling (community-aware ordering)
    # =====================================================================
    println("\n[3/6] Applying LLP relabeling...")
    t2 = time()
    llp_map = relabel_vertices_llp(g, :sym; passes=10)
    g_ord = relabel_graph(g, llp_map)
    relabel_time = time() - t2
    println("  Relabel time: $(round(relabel_time, digits=2))s")
    println("  Reordered graph: $(nv(g_ord)) vertices, $(ne(g_ord)) edges")

    neighbor_lists = Dict{UInt32, Vector{UInt32}}()
    for v in vertices(g_ord)
        nbs = outneighbors(g_ord, v)
        neighbor_lists[UInt32(v)] = sort(UInt32.(collect(nbs)))
    end

    # =====================================================================
    # 4. Define GNN-based Q-learning framework
    # =====================================================================
    println("\n[4/6] Defining GNN-based Q-learning framework...")

    degree_to_bin(d) = (d == 0 && return 1; d <= 3 && return 2; d <= 10 && return 3; d <= 50 && return 4; d <= 200 && return 5; 6)

    gnn_scores = vec(gnn.Z2)
    neighbor_scores = zeros(n)
    for v in vertices(g)
        nbs = outneighbors(g, v)
        if !isempty(nbs)
            neighbor_scores[v] = mean(gnn_scores[u] for u in nbs)
        end
    end

    gnn_score_hist = fit(Histogram, gnn_scores, nbins=GNN_SCORE_BINS)
    neighbor_score_hist = fit(Histogram, neighbor_scores[neighbor_scores .!= 0], nbins=GNN_NEIGHBOR_SCORE_BINS)
    gnn_score_edges = gnn_score_hist.edges[1]
    neighbor_score_edges = neighbor_score_hist.edges[1]

    score_to_bin(score, edges) = max(1, min(length(edges) - 1, searchsortedlast(edges, score)))

    # Create mapping from ordered ID back to original ID for GNN feature lookup
    g_map_inv = Dict(v => k for (k, v) in llp_map)

    function extract_gnn_feature_index(v_id, g_map, gnn_model, neighbor_lists_map)
        v = UInt32(v_id)
        degree = length(get(neighbor_lists_map, v, UInt32[]))
        d_bin = degree_to_bin(degree)

        # Map ordered vertex ID back to original ID to get GNN features
        original_v = g_map[v]
        score = gnn_model.Z2[original_v, 1]
        s_bin = score_to_bin(score, gnn_score_edges)

        orig_neighbors = outneighbors(g, original_v)
        avg_neighbor_score = isempty(orig_neighbors) ? 0.0 : mean(gnn_model.Z2[u, 1] for u in orig_neighbors)
        ns_bin = score_to_bin(avg_neighbor_score, neighbor_score_edges)

        return (d_bin - 1) * (GNN_SCORE_BINS * GNN_NEIGHBOR_SCORE_BINS) +
               (s_bin - 1) * GNN_NEIGHBOR_SCORE_BINS +
               ns_bin
    end

    # =====================================================================
    # 5. Train RL policy (Tabular Q-learning with GNN features)
    # =====================================================================
    println("\n[5/6] Training RL policy with GNN features ($rl_episodes episodes)...")
    env = CompressionEnv(neighbor_lists; ref_window_size=1024)
    policy = QPolicy(; num_states=NUM_Q_STATES, alpha=0.1, gamma=0.9,
                     epsilon=0.5, epsilon_decay=0.97, epsilon_min=0.02)

    t3 = time()
    episode_bpe = Float64[]
    best_bpe = Inf
    for ep in 1:rl_episodes
        reset!(env)
        while !env.done
            v = current_vertex(env)
            s_idx = extract_gnn_feature_index(v, g_map_inv, gnn, neighbor_lists)
            a_idx = select_action(policy, s_idx)
            action = action_from_index(a_idx)

            _, reward, done = step!(env, action)
            
            next_v = current_vertex(env)
            ns_idx = done ? s_idx : extract_gnn_feature_index(next_v, g_map_inv, gnn, neighbor_lists)
            update!(policy, s_idx, a_idx, reward, ns_idx)
        end
        bpe = get_bits_per_edge(env)
        push!(episode_bpe, bpe)
        best_bpe = min(best_bpe, bpe)
        decay_epsilon!(policy)
        if ep == 1 || ep % 5 == 0
            println("  Episode $ep: train_bpe=$(round(bpe, digits=4)), best_bpe=$(round(best_bpe, digits=4)), eps=$(round(policy.epsilon, digits=3))")
        end
    end
    rl_time = time() - t3
    println("\n  RL training time: $(round(rl_time, digits=1))s")

    # =====================================================================
    # 6. Evaluate and save final policy
    # =====================================================================
    println("\n[6/6] Evaluating and saving policy...")
    final_bpe = 0.0
    reset!(env)
    while !env.done
        v = current_vertex(env)
        s_idx = extract_gnn_feature_index(v, g_map_inv, gnn, neighbor_lists)
        a_idx = best_action(policy, s_idx)
        action = action_from_index(a_idx)
        step!(env, action)
    end
    final_bpe = get_bits_per_edge(env)
    println("  Final greedy BPE: $(round(final_bpe, digits=4))")
    
    policy_path = joinpath(out_dir, "cnr2000_gcn_rl_policy.qpolicy")
    save_policy(policy, policy_path)
    println("  Policy saved to: $policy_path")

    println("\n" * "=" ^ 70)
    println("Summary")
    println("=" ^ 70)
    println("  Graph:           CNR-2000 ($n vertices, $m edges)")
    println("  GNN proxy loss:  $(round(proxy_losses[1], digits=2)) → $(round(proxy_losses[end], digits=2))")
    println("  RL episodes:     $rl_episodes")
    println("  Final BPE:       $(round(final_bpe, digits=4))")
    println("=" ^ 70)
end

main()
