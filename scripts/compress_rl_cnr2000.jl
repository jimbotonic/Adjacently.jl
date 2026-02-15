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

# Compress CNR-2000 using a trained RL policy and verify roundtrip.
#
# Usage: julia scripts/compress_rl_cnr2000.jl [policy_filepath] [policy_id]

using Pkg
Pkg.activate(normpath(joinpath(@__DIR__, "..")))

using LightGraphs: nv, ne, outneighbors, vertices
using Adjacently
using Adjacently.IO: load_adjacency_list_from_csv
using Adjacently.Relabeling: relabel_vertices_rcm, relabel_vertices_llp, relabel_graph
using Adjacently.MGS: write_rl_compressed_mgs3_graph, load_compressed_mgs3_graph

const PROJECT_ROOT = normpath(joinpath(@__DIR__, ".."))
const CNR_CSV = joinpath(PROJECT_ROOT, "datasets", "webgraph", "cnr-2000", "cnr-2000.csv")
const DEFAULT_POLICY = joinpath(PROJECT_ROOT, "policies", "cnr2000_policy.qpolicy")

function main()
    # Use "greedy" as first arg for greedy mode, otherwise path to policy file
    mode = length(ARGS) >= 1 ? ARGS[1] : "greedy"
    policy_id = length(ARGS) >= 2 ? parse(Int, ARGS[2]) : 1

    println("=" ^ 70)
    println("RL-Compressed CNR-2000")
    println("=" ^ 70)

    # 1. Determine policy
    policy_path = nothing
    if mode == "greedy"
        println("\nMode: greedy per-vertex optimization")
    else
        policy_path = mode
        if !isfile(policy_path)
            error("Policy file not found: $policy_path\nRun scripts/train_rl_cnr2000.jl first to train a policy, or use 'greedy'.")
        end
        println("\nPolicy: $policy_path (ID=$policy_id)")
    end

    # 2. Load graph
    println("\nLoading CNR-2000...")
    if !isfile(CNR_CSV)
        error("CNR-2000 CSV not found at $CNR_CSV")
    end
    g = load_adjacency_list_from_csv(CNR_CSV, ',', true)
    n, m = nv(g), ne(g)
    println("  Vertices: $n, Edges: $m")

    # 3. LLP relabeling (community-aware ordering for better reference compression)
    println("Applying LLP relabeling...")
    llp_mapping = relabel_vertices_llp(g, :sym; passes=10)
    g_rcm = relabel_graph(g, llp_mapping)

    # 4. Compress with RL policy
    output_base = joinpath(PROJECT_ROOT, "tmp", "cnr2000_rl_policy$(policy_id)")
    mkpath(dirname(output_base))
    println("\nCompressing with RL policy...")
    t0 = time()
    write_rl_compressed_mgs3_graph(g_rcm, output_base, policy_path, policy_id;
                                    coding_scheme=:children, ref_window_size=7,
                                    integer_encoding=:fibonacci)
    compress_time = time() - t0

    mgz_file = output_base * ".mgz"
    file_size = filesize(mgz_file)
    bpe = 8.0 * file_size / m
    println("  File size: $(file_size) bytes")
    println("  Bits/edge: $(round(bpe, digits=4))")
    println("  Compress time: $(round(compress_time, digits=2))s")

    # 5. Decompress and verify roundtrip
    println("\nDecompressing and verifying...")
    t1 = time()
    g_loaded = load_compressed_mgs3_graph(mgz_file)
    decompress_time = time() - t1
    n2, m2 = nv(g_loaded), ne(g_loaded)
    println("  Loaded: $n2 vertices, $m2 edges")
    println("  Decompress time: $(round(decompress_time, digits=2))s")

    if m2 == ne(g_rcm)
        println("  Roundtrip: OK (edge count matches)")
    else
        println("  Roundtrip: MISMATCH (expected $(ne(g_rcm)) edges, got $m2)")
    end

    # Cleanup
    rm(mgz_file; force=true)

    println("\n" * "=" ^ 70)
end

main()
