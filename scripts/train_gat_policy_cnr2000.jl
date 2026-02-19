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

# Train GAT-based GNN and PolicyNetwork (GnnACPolicy) for compression on CNR-2000.
#
# Pipeline:
#   1. Load CNR-2000
#   2. Build GAT-based GNN and train proxy
#   3. Relabel graph with LLP (community-aware ordering)
#   4. Create and train GnnACPolicy end-to-end with GNN
#   5. Save GNN model and policy
#   6. Test compression with trained policy
#
# Usage: julia scripts/train_gat_policy_cnr2000.jl [gnn_epochs] [rl_episodes] [hidden_dim]

using Pkg
Pkg.activate(normpath(joinpath(@__DIR__, "..")))

using LightGraphs: nv, ne, outneighbors, vertices, eltype
using Adjacently
using Adjacently.IO: load_adjacency_list_from_csv
using Adjacently.Relabeling: relabel_graph, relabel_vertices_llp
using Adjacently.GNN: GNNModel, gnn_forward, TrainConfig as GNNTrainConfig, train_gnn_proxy!,
                      save_gnn_model, load_gnn_weights!, gnn_backward_rl!
using Adjacently.RL: CompressionEnv, GnnACPolicy, PolicyNetwork, save_gnn_ac_policy, load_gnn_ac_policy,
                      action_from_index, best_action_gnn, train_gnn_ac_e2e!, evaluate_gnn_actor_critic,
                      get_bits_per_edge, reset!, step!, Action, NUM_ACTIONS,
                      current_vertex
using Adjacently.MGS: write_rl_compressed_mgs3_graph
using Statistics: mean

const PROJECT_ROOT = normpath(joinpath(@__DIR__, ".."))
const CNR_CSV = joinpath(PROJECT_ROOT, "datasets", "webgraph", "cnr-2000", "cnr-2000.csv")

function main()
    gnn_epochs = length(ARGS) >= 1 ? parse(Int, ARGS[1]) : 100
    rl_episodes = length(ARGS) >= 2 ? parse(Int, ARGS[2]) : 200
    hidden_dim = length(ARGS) >= 3 ? parse(Int, ARGS[3]) : 32

    println("=" ^ 70)
    println("GAT + PolicyNetwork Compression Training — CNR-2000")
    println("  gnn_epochs=$gnn_epochs, rl_episodes=$rl_episodes, hidden_dim=$hidden_dim")
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
    # 2. Train GAT-based GNN proxy
    # =====================================================================
    println("\n[2/6] Training GAT-based GNN proxy ($gnn_epochs epochs)...")
    t1 = time()
    gnn = GNNModel(g; hidden_dim=hidden_dim, seed=42, feature_set=:enhanced, use_gat=true)
    println("  GAT GNN: $(gnn.input_dim) → $(gnn.hidden_dim) → 1")
    println("  Using GAT with $(length(gnn.gat_layers)) layers")

    # Check for existing model
    gnn_path = joinpath(PROJECT_ROOT, "policies", "cnr2000_gat_gnn_model.bin")
    if isfile(gnn_path)
        try
            load_gnn_weights!(gnn, gnn_path)
            println("  Loaded existing GAT GNN weights from: $gnn_path")
        catch err
            println("  Warning: could not load GAT GNN weights (", err, ") — starting fresh")
        end
    end

    gnn_cfg = GNNTrainConfig(proxy_epochs=gnn_epochs, proxy_lr=0.001,
                              reinforce_epochs=0, reinforce_lr=0.0,
                              sigma=0.0, baseline_ema=0.9, grad_clip_norm=1.0)
    proxy_losses = train_gnn_proxy!(gnn, g, gnn_cfg)
    gnn_time = time() - t1
    gnn_forward(gnn)  # populate caches H1/Z2
    println("  Proxy loss: start=$(round(proxy_losses[1], digits=4)), end=$(round(proxy_losses[end], digits=4))")
    println("  GAT GNN training time: $(round(gnn_time, digits=2))s")

    # Save GNN model
    out_dir = joinpath(PROJECT_ROOT, "policies")
    mkpath(out_dir)
    save_gnn_model(gnn, gnn_path)
    println("  GAT GNN model saved to: $gnn_path")

    # =====================================================================
    # 3. Apply LLP relabeling
    # =====================================================================
    println("\n[3/6] Applying LLP relabeling...")
    t2 = time()
    llp_map = relabel_vertices_llp(g, :sym; passes=10)
    g_ord = relabel_graph(g, llp_map)
    relabel_time = time() - t2
    println("  Relabel time: $(round(relabel_time, digits=2))s")

    neighbor_lists = Dict{UInt32, Vector{UInt32}}()
    for v in vertices(g_ord)
        nbs = outneighbors(g_ord, v)
        neighbor_lists[UInt32(v)] = sort(UInt32.(collect(nbs)))
    end

    # Create inverse mapping: ordered ID -> original ID
    g_map_inv = Dict(v => k for (k, v) in llp_map)

    # =====================================================================
    # 4. Create PolicyNetwork (GnnACPolicy)
    # =====================================================================
    println("\n[4/6] Creating PolicyNetwork (GnnACPolicy)...")
    feat_dim = hidden_dim + 1  # H1[v,:] (hidden_dim) + Z2[v] (1)
    policy = GnnACPolicy(feat_dim; num_actions=NUM_ACTIONS,
                         actor_lr=0.05, critic_lr=0.1,
                         gamma=0.0, temperature=1.0)
    println("  Feature dimension: $feat_dim (hidden_dim=$hidden_dim + 1)")
    println("  Num actions: $NUM_ACTIONS")

    # Feature function: extract GNN features for vertex v (in ordered graph)
    function feature_fn(v::UInt32)
        vi = Int(v)
        # Map ordered vertex ID back to original ID for GNN feature lookup
        original_v = g_map_inv[vi]
        h1 = gnn.H1[original_v, :]  # hidden layer features
        z2 = gnn.Z2[original_v, 1]  # output score
        return vcat(h1, z2)  # feature vector of length feat_dim
    end

    # =====================================================================
    # 5. Train PolicyNetwork end-to-end with GNN
    # =====================================================================
    println("\n[5/6] Training PolicyNetwork end-to-end ($rl_episodes episodes)...")
    env = CompressionEnv(neighbor_lists; ref_window_size=7)

    t3 = time()
    episode_bpe = train_gnn_ac_e2e!(env, policy, gnn, feature_fn;
                                     episodes=rl_episodes, gnn_lr=0.0001,
                                     gnn_update_every=5, eval_every=10, verbose=true)
    rl_time = time() - t3
    println("\n  RL training time: $(round(rl_time, digits=1))s")

    # Final evaluation
    reset!(env)
    while !env.done
        v = current_vertex(env)
        x = feature_fn(v)
        a_idx = best_action_gnn(policy, x)
        action = action_from_index(a_idx)
        step!(env, action)
    end
    final_bpe = get_bits_per_edge(env)
    println("  Final greedy BPE: $(round(final_bpe, digits=4))")

    # Save policy
    policy_path = joinpath(out_dir, "cnr2000_gat_policy.bin")
    save_gnn_ac_policy(policy, policy_path)
    println("  Policy saved to: $policy_path")

    # Save updated GNN model (fine-tuned during E2E training)
    save_gnn_model(gnn, gnn_path)
    println("  Updated GAT GNN model saved to: $gnn_path")

    # =====================================================================
    # 6. Test compression with trained policy
    # =====================================================================
    println("\n[6/6] Testing compression with trained policy...")

    # Compute vertex actions using the trained policy
    vertex_actions = Dict{UInt32, Int}()
    for v in vertices(g_ord)
        v_u32 = UInt32(v)
        x = feature_fn(v_u32)
        a_idx = best_action_gnn(policy, x)
        vertex_actions[v_u32] = a_idx
    end

    # Test write
    test_path = joinpath(PROJECT_ROOT, "tmp", "cnr2000_gat_test")
    mkpath(dirname(test_path))

    write_rl_compressed_mgs3_graph(g_ord, test_path, vertex_actions;
                                    coding_scheme=:children, ref_window_size=7,
                                    integer_encoding=:fibonacci)

    # Get file size
    file_size = stat(test_path * ".mgz").size
    achieved_bpe = file_size * 8.0 / m
    println("  Compressed file: $(test_path).mgz")
    println("  File size: $file_size bytes")
    println("  Achieved BPE: $(round(achieved_bpe, digits=4))")

    # Cleanup
    rm(test_path * ".mgz", force=true)

    println("\n" * "=" ^ 70)
    println("Training Summary")
    println("=" ^ 70)
    println("  Graph:           CNR-2000 ($n vertices, $m edges)")
    println("  GNN proxy loss:  $(round(proxy_losses[1], digits=4)) → $(round(proxy_losses[end], digits=4))")
    println("  RL episodes:     $rl_episodes")
    println("  Final RL BPE:    $(round(final_bpe, digits=4))")
    println("  Achieved BPE:    $(round(achieved_bpe, digits=4))")
    println("  GAT GNN model:   $gnn_path")
    println("  Policy:          $policy_path")
    println("=" ^ 70)
end

main()