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

# Train an RL policy to learn optimal per-vertex compression encoding
# decisions on the CNR-2000 web graph.
#
# Usage: julia scripts/train_rl_cnr2000.jl [num_episodes]

using Pkg
Pkg.activate(normpath(joinpath(@__DIR__, "..")))

using LightGraphs: nv, ne, outneighbors, vertices
using Adjacently
using Adjacently.IO: load_adjacency_list_from_csv
using Adjacently.Relabeling: relabel_vertices_rcm, relabel_graph
using Adjacently.Graph: get_neighbor_lists
using Adjacently.RL: CompressionEnv, QPolicy, TrainingConfig, train!, evaluate, run_baseline, save_policy

const PROJECT_ROOT = normpath(joinpath(@__DIR__, ".."))
const CNR_CSV = joinpath(PROJECT_ROOT, "datasets", "webgraph", "cnr-2000", "cnr-2000.csv")

function main()
    num_episodes = length(ARGS) >= 1 ? parse(Int, ARGS[1]) : 100

    println("=" ^ 70)
    println("RL Compression Policy Training on CNR-2000")
    println("=" ^ 70)

    # 1. Load graph
    println("\nLoading CNR-2000 dataset...")
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

    # 2. RCM relabeling
    println("\nApplying RCM relabeling...")
    t1 = time()
    rcm_mapping = relabel_vertices_rcm(g, :out)
    g_rcm = relabel_graph(g, rcm_mapping)
    relabel_time = time() - t1
    println("  Relabel time: $(round(relabel_time, digits=2))s")

    # 3. Build neighbor lists (must use Unsigned type for CompressionEnv)
    println("\nBuilding neighbor lists...")
    neighbor_lists = Dict{UInt32, Vector{UInt32}}()
    for v in vertices(g_rcm)
        nbs = outneighbors(g_rcm, v)
        neighbor_lists[UInt32(v)] = sort(UInt32.(collect(nbs)))
    end
    println("  Done: $(length(neighbor_lists)) vertices")

    # 4. Create environment and policy
    println("\nCreating compression environment...")
    env = CompressionEnv(neighbor_lists; ref_window_size=1024)
    policy = QPolicy(; alpha=0.1, gamma=0.99, epsilon=0.5,
                       epsilon_decay=0.97, epsilon_min=0.02)
    println("  State space:  $(policy.num_states) states")
    println("  Action space: $(policy.num_actions) actions")

    # 5. Train
    println("\n" * "=" ^ 70)
    config = TrainingConfig(; num_episodes=num_episodes, eval_every=5, verbose=true)
    t2 = time()
    result = train!(env, policy, config)
    train_time = time() - t2

    println("\nTotal training time: $(round(train_time, digits=1))s")
    println("Average per episode: $(round(train_time / num_episodes, digits=1))s")

    # 6. Save policy
    policy_dir = joinpath(PROJECT_ROOT, "policies")
    mkpath(policy_dir)
    policy_path = joinpath(policy_dir, "cnr2000_policy.qpolicy")
    save_policy(result.policy, policy_path)
    println("\nPolicy saved to: $policy_path")
    println("=" ^ 70)
end

main()
