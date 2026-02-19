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
# Verifies the saved graph by reloading and comparing it to the original.
#
# Usage: julia scripts/train_and_verify_gat_cnr2000.jl [gnn_epochs] [rl_episodes] [hidden_dim] [relabeling_mode]
#   relabeling_mode: bisection (default), llp, or combined (LLP then bisection)

using Pkg
Pkg.activate(normpath(joinpath(@__DIR__, "..")))

using LightGraphs: nv, ne, outneighbors, vertices, eltype, is_directed
using Adjacently
using Adjacently.IO: load_adjacency_list_from_csv
using Adjacently.Relabeling: relabel_graph, relabel_vertices_llp, relabel_vertices_bisection
using Adjacently.GNN: GNNModel, gnn_forward, TrainConfig as GNNTrainConfig, train_gnn_proxy!,
                      save_gnn_model, load_gnn_weights!, gnn_backward_rl!, relabel_vertices_gnn
using Adjacently.RL: CompressionEnv, GnnACPolicy, PolicyNetwork,
                      action_from_index, best_action_gnn, train_gnn_ac_e2e!, evaluate_gnn_actor_critic,
                      get_bits_per_edge, reset!, step!, Action, NUM_ACTIONS,
                      current_vertex, save_gnn_ac_policy
using Adjacently.MGS: write_rl_compressed_mgs3_graph, load_compressed_mgs3_graph
using Statistics: mean, std

const PROJECT_ROOT = normpath(joinpath(@__DIR__, ".."))
const CNR_CSV = joinpath(PROJECT_ROOT, "datasets", "webgraph", "cnr-2000", "cnr-2000.csv")
const REF_WINDOW_SIZE = 7

function main()
    gnn_epochs = length(ARGS) >= 1 ? parse(Int, ARGS[1]) : 100
    rl_episodes = length(ARGS) >= 2 ? parse(Int, ARGS[2]) : 500
    hidden_dim = length(ARGS) >= 3 ? parse(Int, ARGS[3]) : 16
    relabeling_mode = length(ARGS) >= 4 ? Symbol(ARGS[4]) : :bisection

    println("=" ^ 70)
    println("GAT + PolicyNetwork Compression Training & Verification — CNR-2000")
    println("  gnn_epochs=$gnn_epochs, rl_episodes=$rl_episodes, hidden_dim=$hidden_dim")
    println("  relabeling_mode=$relabeling_mode, ref_window=$REF_WINDOW_SIZE")
    println("=" ^ 70)

    # 1. Load graph
    println("\n[1/6] Loading CNR-2000 dataset...")
    t_load = @elapsed begin
        g = load_adjacency_list_from_csv(CNR_CSV, ',', true)
    end
    n, m = nv(g), ne(g)
    println("  Vertices: $n, Edges: $m (loaded in $(round(t_load, digits=2))s)")

    # 2. Train GAT-based GNN proxy with compression features
    println("\n[2/6] Preparing GAT-based GNN proxy...")
    t_gnn = @elapsed begin
        gnn = GNNModel(g; hidden_dim=hidden_dim, seed=42, feature_set=:compression, use_gat=true)

        gnn_path = joinpath(PROJECT_ROOT, "policies", "cnr2000_gat_gnn_model.bin")
        loaded = false
        if isfile(gnn_path)
            println("  Loading existing GNN weights from $gnn_path...")
            try
                load_gnn_weights!(gnn, gnn_path)
                gnn_forward(gnn)  # verify dimensions are compatible
                loaded = true
            catch e
                println("  Warning: Incompatible saved weights ($(e)). Retraining from scratch...")
                # Re-create model to reset weights
                gnn = GNNModel(g; hidden_dim=hidden_dim, seed=42, feature_set=:compression, use_gat=true)
            end
        end
        if !loaded
            println("  Training GNN proxy for $gnn_epochs epochs...")
            gnn_cfg = GNNTrainConfig(proxy_epochs=gnn_epochs, proxy_lr=0.001, reinforce_epochs=0, reinforce_lr=0.0)
            train_gnn_proxy!(gnn, g, gnn_cfg)
            save_gnn_model(gnn, gnn_path)
            gnn_forward(gnn)
        end
    end
    println("  GNN prepared in $(round(t_gnn, digits=2))s")

    # 3. Apply relabeling
    println("\n[3/6] Applying $relabeling_mode relabeling...")
    t_relabel = @elapsed begin
        if relabeling_mode == :bisection
            relabel_map = relabel_vertices_bisection(g)
        elseif relabeling_mode == :llp
            relabel_map = relabel_vertices_llp(g, :sym; passes=10, K=10)
        elseif relabeling_mode == :combined
            # LLP first, then bisection on the reordered graph
            llp_map = relabel_vertices_llp(g, :sym; passes=10, K=10)
            g_llp = relabel_graph(g, llp_map)
            bisect_map = relabel_vertices_bisection(g_llp)
            # Compose mappings: original -> llp -> bisection
            relabel_map = Dict{eltype(g),eltype(g)}()
            for (old_v, llp_v) in llp_map
                relabel_map[old_v] = bisect_map[llp_v]
            end
        else
            error("Unknown relabeling mode: $relabeling_mode. Use :bisection, :llp, or :combined")
        end
        g_ord = relabel_graph(g, relabel_map)
    end
    println("  Relabeling completed in $(round(t_relabel, digits=2))s")

    n_bits_v = convert(UInt8, ceil(log(2, n)))
    V_type = Adjacently.Util.infer_uint_custom_type(n_bits_v)

    neighbor_lists = Dict{V_type, Vector{V_type}}()
    for v in vertices(g_ord)
        neighbor_lists[V_type(v)] = sort(V_type.(collect(outneighbors(g_ord, v))))
    end
    g_map_inv = Dict(V_type(v) => V_type(k) for (k, v) in relabel_map)

    # 4. Create PolicyNetwork
    println("\n[4/6] Creating PolicyNetwork...")
    feat_dim = hidden_dim + 1
    policy = GnnACPolicy(feat_dim; num_actions=NUM_ACTIONS,
                         actor_lr=0.001, critic_lr=0.001,
                         gamma=0.99, temperature=2.0, entropy_coeff=1.0)

    function feature_fn(v)
        v_conv = V_type(Int(v))
        original_v = get(g_map_inv, v_conv, V_type(1))
        h1 = gnn.H1[Int(original_v), :]
        z2 = gnn.Z2[Int(original_v), 1]
        res = vcat(h1, z2)
        return any(isnan, res) ? fill(0.0, length(res)) : res
    end

    # 5. Train end-to-end
    println("\n[5/6] Training PolicyNetwork end-to-end...")
    env = CompressionEnv(neighbor_lists; ref_window_size=REF_WINDOW_SIZE)
    t_train = @elapsed begin
        train_gnn_ac_e2e!(env, policy, gnn, feature_fn;
                          episodes=rl_episodes, gnn_lr=0.00001,
                          gnn_update_every=10, eval_every=10, verbose=true)
    end
    println("  E2E training completed in $(round(t_train, digits=2))s")

    # 6. Test compression and Verify
    println("\n[6/6] Testing compression and verifying...")
    vertex_actions = Dict{V_type, Int}()
    for v in vertices(g_ord)
        v_v = V_type(v)
        vertex_actions[v_v] = best_action_gnn(policy, feature_fn(v_v))
    end

    out_dir = joinpath(PROJECT_ROOT, "policies")
    mkpath(out_dir)
    final_mgz_path = joinpath(out_dir, "cnr2000_gat_trained")

    write_rl_compressed_mgs3_graph(g_ord, final_mgz_path, vertex_actions;
                                    coding_scheme=:children, ref_window_size=REF_WINDOW_SIZE,
                                    integer_encoding=:fibonacci)

    file_size = stat(final_mgz_path * ".mgz").size
    achieved_bpe = file_size * 8.0 / m
    println("  Achieved BPE: $(round(achieved_bpe, digits=4))")

    # Action Stats
    counts = Dict{Int,Int}()
    for a in values(vertex_actions)
        counts[a] = get(counts, a, 0) + 1
    end
    sorted_actions = sort(collect(counts), by=x->x[2], rev=true)
    println("\nAction Statistics:")
    for (a_idx, count) in sorted_actions[1:min(10, length(sorted_actions))]
        action = action_from_index(a_idx)
        pct = round(100.0 * count / n, digits=2)
        println("  Action $a_idx: $(action.reference_mode), $(action.encoding_type), mil=$(action.min_interval_length) - $count vertices ($pct%)")
    end

    println("\n[Verification] Loading saved graph...")
    GC.gc()
    g_decoded = load_compressed_mgs3_graph(final_mgz_path * ".mgz")

    mismatch_count = 0
    for v in vertices(g_ord)
        orig_nbs = sort(V_type.(outneighbors(g_ord, v)))
        deco_nbs = sort(V_type.(outneighbors(g_decoded, v)))
        if orig_nbs != deco_nbs
            mismatch_count += 1
            if mismatch_count <= 1
                println("    FIRST Mismatch at vertex $v")
                println("      Action: $(action_from_index(vertex_actions[V_type(v)]))")
                println("      Original: $(orig_nbs)")
                println("      Decoded:  $(deco_nbs)")
            end
        end
    end

    if mismatch_count == 0
        println("  SUCCESS: Decoded graph is identical to original!")
    else
        println("  FAILED: $mismatch_count vertices have mismatched neighbor lists.")
    end

    # Save artifacts (including post-RL GNN weights)
    println("\n[Saving Artifacts]")
    Adjacently.GNN.save_gnn_model(gnn, joinpath(out_dir, "cnr2000_gat_gnn_model.bin"))
    save_gnn_ac_policy(policy, gnn, joinpath(out_dir, "cnr2000_gat_policy.bin"))
    println("  Artifacts saved to policies/ directory.")
end

main()
