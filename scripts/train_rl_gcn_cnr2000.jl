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
using Statistics: mean, std
using Adjacently
using Adjacently.IO: load_adjacency_list_from_csv
using Adjacently.Relabeling: relabel_graph, relabel_vertices_llp
using Adjacently.GNN: GNNModel, gnn_forward, TrainConfig as GNNTrainConfig, train_gnn_proxy!
using Adjacently.RL: CompressionEnv,
                      GnnACPolicy, train_gnn_actor_critic!, evaluate_gnn_actor_critic,
                      train_gnn_ac_e2e!, save_gnn_ac_policy, load_gnn_ac_policy,
                      best_action_gnn,
                      current_vertex, action_from_index
using Adjacently.MGS: write_rl_compressed_mgs3_graph, load_compressed_mgs3_graph

const PROJECT_ROOT = normpath(joinpath(@__DIR__, ".."))
const CNR_CSV = joinpath(PROJECT_ROOT, "datasets", "webgraph", "cnr-2000", "cnr-2000.csv")

function main()
    rl_episodes = length(ARGS) >= 1 ? parse(Int, ARGS[1]) : 100
    gnn_epochs = length(ARGS) >= 2 ? parse(Int, ARGS[2]) : 50
    hidden_dim = length(ARGS) >= 3 ? parse(Int, ARGS[3]) : 32
    feature_set = length(ARGS) >= 4 ? Symbol(ARGS[4]) : :enhanced

    println("=" ^ 70)
    println("RL + GCN Compression Pipeline — CNR-2000")
    println("  hidden_dim=$hidden_dim, feature_set=$feature_set")
    println("=" ^ 70)

    # =====================================================================
    # 1. Load graph
    # =====================================================================
    println("\n[1/7] Loading CNR-2000 dataset...")
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
    println("\n[2/7] Training GNN proxy (feature learning, $gnn_epochs epochs)...")
    t1 = time()
    gnn = GNNModel(g; hidden_dim=hidden_dim, seed=42, feature_set=feature_set)
    println("  GNN: $(gnn.input_dim) → $(gnn.hidden_dim) → 1 (features=$feature_set)")
    # Resume from previous weights if available
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
                              sigma=0.0, baseline_ema=0.9)
    proxy_losses = train_gnn_proxy!(gnn, g, gnn_cfg)
    gnn_time = time() - t1
    gnn_forward(gnn)  # populate caches H1/Z2
    println("  Proxy loss: start=$(round(proxy_losses[1], digits=2)), end=$(round(proxy_losses[end], digits=2))")
    println("  GNN proxy time: $(round(gnn_time, digits=2))s")

    # Save GNN model
    println("  Saving GNN model...")
    out_dir = joinpath(PROJECT_ROOT, "policies")
    mkpath(out_dir)
    gnn_path = joinpath(out_dir, "cnr2000_gnn_model.bin")
    Adjacently.GNN.save_gnn_model(gnn, gnn_path)
    println("  GNN model: $gnn_path")

    # Dump proxy losses for convergence analysis
    begin
        losses_dir = joinpath(PROJECT_ROOT, "tmp")
        mkpath(losses_dir)
        out_loss = joinpath(losses_dir, "gnn_proxy_losses.csv")
        open(out_loss, "w") do f
            write(f, "epoch,loss\n")
            for (i, L) in enumerate(proxy_losses)
                write(f, string(i, ",", L, "\n"))
            end
        end
        println("  Wrote losses: $out_loss")
    end

    # =====================================================================
    # 3. Apply LLP relabeling (community-aware ordering)
    # =====================================================================
    println("\n[3/7] Applying LLP relabeling...")
    t2 = time()
    llp_map = relabel_vertices_llp(g, :sym; passes=10)
    g_ord = relabel_graph(g, llp_map)
    relabel_time = time() - t2
    println("  Relabel time: $(round(relabel_time, digits=2))s")
    println("  Reordered graph: $(nv(g_ord)) vertices, $(ne(g_ord)) edges")

    # =====================================================================
    # 4. Train RL policy (GNN-augmented actor-critic) on LLP-ordered graph
    # =====================================================================
    println("\n[4/7] Training RL policy with GNN features ($rl_episodes episodes)...")
    neighbor_lists = Dict{UInt32, Vector{UInt32}}()
    for v in vertices(g_ord)
        nbs = outneighbors(g_ord, v)
        neighbor_lists[UInt32(v)] = sort(UInt32.(collect(nbs)))
    end

    env = CompressionEnv(neighbor_lists; ref_window_size=1024)
    # Feature function from GNN caches (use H1 row concatenated with scalar Z2)
    # Normalize features to zero-mean unit-variance to prevent NaN in actor-critic
    feat_dim = size(gnn.H1, 2) + 1
    raw_features = hcat(gnn.H1, gnn.Z2)  # (n × feat_dim)
    feat_mean = vec(mean(raw_features, dims=1))
    feat_std = vec(std(raw_features, dims=1))
    feat_std[feat_std .< 1e-8] .= 1.0  # avoid division by zero
    println("  Feature stats: mean range [$(round(minimum(feat_mean),digits=2)), $(round(maximum(feat_mean),digits=2))], std range [$(round(minimum(feat_std),digits=2)), $(round(maximum(feat_std),digits=2))]")
    feature_fn = v -> (vcat(gnn.H1[Int(v), :], gnn.Z2[Int(v)]) .- feat_mean) ./ feat_std
    # Learning rates scaled for online per-step updates over ~325K vertices/episode
    policy = GnnACPolicy(feat_dim; num_actions=12, actor_lr=0.0001, critic_lr=0.001,
                         gamma=0.0, temperature=1.0)

    t3 = time()
    _ = train_gnn_actor_critic!(env, policy, feature_fn; episodes=rl_episodes, eval_every=5, verbose=true)
    rl_time = time() - t3
    println("\n  RL training time: $(round(rl_time, digits=1))s")

    # =====================================================================
    # 5. End-to-end GNN+RL fine-tuning
    # =====================================================================
    e2e_episodes = max(rl_episodes ÷ 2, 10)
    println("\n[5/7] End-to-end GNN+RL fine-tuning ($e2e_episodes episodes)...")
    t_e2e = time()
    _ = train_gnn_ac_e2e!(env, policy, gnn, feature_fn;
                          episodes=e2e_episodes, gnn_lr=0.0001,
                          gnn_update_every=5, eval_every=5, verbose=true)
    e2e_time = time() - t_e2e
    println("\n  E2E fine-tuning time: $(round(e2e_time, digits=1))s")

    # Save AC policy
    println("  Saving AC policy...")
    ac_path = joinpath(out_dir, "cnr2000_gnn_ac_policy.bin")
    save_gnn_ac_policy(policy, ac_path)
    println("  AC policy: $ac_path")

    # Save updated GNN model (after e2e fine-tuning)
    gnn_path_e2e = joinpath(out_dir, "cnr2000_gnn_model_e2e.bin")
    Adjacently.GNN.save_gnn_model(gnn, gnn_path_e2e)
    println("  GNN model (e2e): $gnn_path_e2e")

    # Verify save/load roundtrip for policy
    policy_loaded = load_gnn_ac_policy(ac_path; actor_lr=0.05, critic_lr=0.1, gamma=0.0, temperature=1.0)
    println("  Policy load check: feat_dim=$(policy_loaded.feat_dim), num_actions=$(policy_loaded.num_actions)")

    # =====================================================================
    # 6. Compress with RL policy (fibonacci encoding)
    # =====================================================================
    println("\n[6/7] Compressing CNR-2000 (LLP + RL policy + fibonacci)...")
    output_base = joinpath(PROJECT_ROOT, "tmp", "cnr2000_gcn_rl")
    mkpath(dirname(output_base))

    # Pre-compute per-vertex actions from trained GNN policy
    println("  Pre-computing vertex actions from GNN policy...")
    vertex_actions = Dict{UInt32,Int}()
    for v in vertices(g_ord)
        x = feature_fn(UInt32(v))
        a = best_action_gnn(policy, x)
        vertex_actions[UInt32(v)] = a
    end

    # Count action distribution
    action_counts = zeros(Int, 12)
    for a in values(vertex_actions)
        action_counts[a] += 1
    end
    println("  Action distribution:")
    for (i, c) in enumerate(action_counts)
        act = action_from_index(i)
        pct = round(100.0 * c / nv(g_ord), digits=1)
        println("    Action $i (ref=$(act.reference_mode), mil=$(act.min_interval_length)): $c ($pct%)")
    end

    t4 = time()
    write_rl_compressed_mgs3_graph(g_ord, output_base, nothing, 1;
        coding_scheme=:children, ref_window_size=7, integer_encoding=:fibonacci,
        vertex_actions=vertex_actions)
    compress_time = time() - t4

    mgz_file = output_base * ".mgz"
    file_size = filesize(mgz_file)
    bpe = 8.0 * file_size / m
    println("  File size:     $(file_size) bytes")
    println("  Bits/edge:     $(round(bpe, digits=4))")
    println("  Compress time: $(round(compress_time, digits=2))s")

    # =====================================================================
    # 7. Verify roundtrip
    # =====================================================================
    println("\n[7/7] Verifying roundtrip...")
    t5 = time()
    g_loaded = load_compressed_mgs3_graph(mgz_file)
    decompress_time = time() - t5
    n2, m2 = nv(g_loaded), ne(g_loaded)
    println("  Loaded: $n2 vertices, $m2 edges")
    println("  Decompress time: $(round(decompress_time, digits=2))s")

    if m2 == ne(g_ord)
        println("  Roundtrip: OK (edge count matches)")
    else
        println("  Roundtrip: MISMATCH (expected $(ne(g_ord)) edges, got $m2)")
    end

    # Cleanup temp file
    rm(mgz_file; force=true)

    # =====================================================================
    # Summary
    # =====================================================================
    println("\n" * "=" ^ 70)
    println("Summary")
    println("=" ^ 70)
    println("  Graph:           CNR-2000 ($n vertices, $m edges)")
    println("  GNN proxy:       $(round(proxy_losses[1], digits=2)) → $(round(proxy_losses[end], digits=2))")
    println("  Final bpe:       $(round(bpe, digits=4))")
    println("  Encoding:        RL policy + fibonacci")
    println("=" ^ 70)
end

main()
