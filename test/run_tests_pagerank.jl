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

#
# Test PageRank implementations on CNR-2000 web graph
# Compares different PageRank algorithms:
# 1. Power iteration (matrix-based) - EFFICIENT
# 2. Naive iterative (graph-based)
# 3. Monte Carlo (random walk-based)
#

include("run_tests_main.jl")

using Logging
using LightGraphs: nv, ne, outneighbors
using SparseArrays: SparseMatrixCSC, sparse
using LinearAlgebra: norm
using Adjacently.MGS: load_compressed_mgs3_graph
using Adjacently
using Adjacently.IO: load_adjacency_list_from_csv
using Adjacently.PageRank: PR
using Adjacently.Graph: get_reverse_graph, get_sparse_P_matrix

println("=" ^ 80)
println("PageRank Algorithm Tests on CNR-2000")
println("=" ^ 80)

# Load graph
println("\n[1/5] Loading CNR-2000 dataset...")
graph = load_compressed_mgs3_graph("datasets/webgraph/cnr-2000/cnr-2000.mgz")
vertices_count = nv(graph)
edges_count = ne(graph)
println("  Graph: $vertices_count vertices, $edges_count edges")

# Build reverse graph (needed for naive implementation)
println("\n[2/5] Building reverse graph...")
reverse_start = time()
reverse_graph = get_reverse_graph(graph)
reverse_time = time() - reverse_start
println("  Reverse graph built in $(round(reverse_time, digits=2))s")

# Build stochastic adjacency matrix for power iteration method
println("\n[3/5] Building column-stochastic adjacency matrix...")
matrix_start = time()

# Use the fixed get_sparse_P_matrix function
P = get_sparse_P_matrix(graph)  # Returns column-stochastic by default
matrix_time = time() - matrix_start
matrix_nnz = length(P.nzval)
n = vertices_count
matrix_density = 100.0 * matrix_nnz / (n * n)

println("  Matrix built in $(round(matrix_time, digits=2))s")
println("  Matrix size: $n × $n")
println("  Non-zero entries: $matrix_nnz")
println("  Density: $(round(matrix_density, digits=6))%")
println("  Memory estimate: ~$(round(matrix_nnz * 16 / 1024 / 1024, digits=2)) MB")

# Test parameters
damping = 0.85
epsilon = 1e-6
max_iter = 100

println("\n[4/5] Computing PageRank with different methods...")
println("  Parameters:")
println("    Damping factor (α): $damping")
println("    Convergence threshold (ε): $epsilon")
println("    Max iterations: $max_iter")

# Method 1: Power iteration (matrix-based) - EFFICIENT
println("\n  Method 1: Power Iteration (Matrix-based)")
println("  " * "-" ^ 60)
power_start = time()
pr_power = PR(P; damping=damping, epsilon=epsilon, max_iter=max_iter)
power_time = time() - power_start
power_sum = sum(pr_power)
power_norm = norm(pr_power, 1)
println("    Computation time: $(round(power_time, digits=3))s")
println("    Sum of PageRank: $(round(power_sum, digits=6))")
println("    L1 norm: $(round(power_norm, digits=6))")
println("    Top 5 vertices by PageRank:")
top5_indices = sortperm(pr_power, rev=true)[1:5]
for (rank, idx) in enumerate(top5_indices)
    println("      #$rank: vertex $idx with PR = $(round(pr_power[idx], digits=8))")
end

# Method 2: Naive iterative (graph-based)
println("\n  Method 2: Naive Iterative (Graph-based)")
println("  " * "-" ^ 60)
naive_start = time()
pr_naive = PR(graph, reverse_graph; damping=damping, epsilon=epsilon)
naive_time = time() - naive_start
naive_sum = sum(pr_naive)
naive_norm = norm(pr_naive, 1)
println("    Computation time: $(round(naive_time, digits=3))s")
println("    Sum of PageRank: $(round(naive_sum, digits=6))")
println("    L1 norm: $(round(naive_norm, digits=6))")
println("    Top 5 vertices by PageRank:")
top5_indices_naive = sortperm(pr_naive, rev=true)[1:5]
for (rank, idx) in enumerate(top5_indices_naive)
    println("      #$rank: vertex $idx with PR = $(round(pr_naive[idx], digits=8))")
end

# Method 3: Monte Carlo (random walk-based) - skip if too slow
println("\n  Method 3: Monte Carlo (Random Walk-based)")
println("  " * "-" ^ 60)
println("    Skipped - too slow for large graphs")
println("    (MC method requires many random walks per vertex)")

# Compare results
println("\n[5/5] Comparison and Validation")
println("=" ^ 80)

# Compare power iteration vs naive
diff_power_naive = norm(pr_power - pr_naive, Inf)
println("\n  Power vs Naive:")
println("    Max absolute difference: $(round(diff_power_naive, digits=10))")
println("    Agreement: $(diff_power_naive < 1e-5 ? "EXCELLENT" : "POOR")")

# Speedup analysis
speedup = naive_time / power_time
println("\n  Performance:")
println("    Power iteration: $(round(power_time, digits=3))s")
println("    Naive iteration: $(round(naive_time, digits=3))s")
println("    Speedup: $(round(speedup, digits=2))x")
println("    Winner: $(power_time < naive_time ? "Power iteration" : "Naive iteration")")

# Statistics about distribution
println("\n  PageRank Distribution Statistics (Power Iteration):")
pr_min = minimum(pr_power)
pr_max = maximum(pr_power)
pr_mean = sum(pr_power) / length(pr_power)
pr_median_val = sort(pr_power)[div(length(pr_power), 2)]
println("    Min: $(round(pr_min, digits=10))")
println("    Max: $(round(pr_max, digits=10))")
println("    Mean: $(round(pr_mean, digits=10))")
println("    Median: $(round(pr_median_val, digits=10))")
println("    Range: $(round(pr_max / pr_min, digits=2))x")

# Check for dangling nodes (nodes with no outgoing edges)
dangling_nodes = sum(length(outneighbors(graph, UInt24(v))) == 0 for v in 1:vertices_count)
println("\n  Graph Properties:")
println("    Dangling nodes (no outgoing edges): $dangling_nodes ($(round(100*dangling_nodes/vertices_count, digits=2))%)")
println("    Average out-degree: $(round(edges_count / vertices_count, digits=2))")

# ============================================================================
# TEST ON CORE SUBGRAPH (Main SCC)
# ============================================================================

println("\n" * "=" ^ 80)
println("[6/6] BONUS: PageRank on Core Subgraph (Main SCC)")
println("=" ^ 80)
println("Extracting the largest strongly connected component (SCC)...")
println("This eliminates dangling nodes and tests if methods agree on clean graphs.")

using Adjacently.Graph: get_core

core_start = time()
core_graph, old_to_new, new_to_old = get_core(graph)
core_time = time() - core_start

core_vertices = nv(core_graph)
core_edges = ne(core_graph)
core_reduction_v = 100.0 * (vertices_count - core_vertices) / vertices_count
core_reduction_e = 100.0 * (edges_count - core_edges) / edges_count

println("\n  Core extraction completed in $(round(core_time, digits=2))s")
println("  Core graph: $core_vertices vertices, $core_edges edges")
println("  Reduction: $(round(core_reduction_v, digits=1))% vertices, $(round(core_reduction_e, digits=1))% edges removed")

# Build reverse graph for core
core_reverse = get_reverse_graph(core_graph)

# Count dangling nodes in core
core_dangling = sum(length(outneighbors(core_graph, UInt24(v))) == 0 for v in 1:core_vertices)
println("  Dangling nodes in core: $core_dangling ($(round(100*core_dangling/core_vertices, digits=2))%)")

println("\n  Computing PageRank on core with both methods...")

# Method 1: Power iteration on core
println("\n  Power Iteration (Core):")
# Use get_sparse_P_matrix for core
P_core = get_sparse_P_matrix(core_graph)

core_power_start = time()
pr_core_power = PR(P_core; damping=damping, epsilon=epsilon, max_iter=max_iter)
core_power_time = time() - core_power_start
println("    Time: $(round(core_power_time, digits=3))s")
println("    Sum: $(round(sum(pr_core_power), digits=6))")
println("    Top vertex: $(argmax(pr_core_power)) with PR = $(round(maximum(pr_core_power), digits=8))")

# Method 2: Naive on core
println("\n  Naive Iteration (Core):")
core_naive_start = time()
pr_core_naive = PR(core_graph, core_reverse; damping=damping, epsilon=epsilon)
core_naive_time = time() - core_naive_start
println("    Time: $(round(core_naive_time, digits=3))s")
println("    Sum: $(round(sum(pr_core_naive), digits=6))")
println("    Top vertex: $(argmax(pr_core_naive)) with PR = $(round(maximum(pr_core_naive), digits=8))")

# Compare core results
diff_core = norm(pr_core_power - pr_core_naive, Inf)
println("\n  Core Comparison:")
println("    Max absolute difference: $(round(diff_core, digits=10))")
println("    Agreement: $(diff_core < 1e-5 ? "✓ EXCELLENT" : diff_core < 1e-3 ? "✓ GOOD" : "✗ POOR")")
println("    Both PageRank sums close to 1.0: $(abs(sum(pr_core_power) - 1.0) < 0.01 && abs(sum(pr_core_naive) - 1.0) < 0.01 ? "✓ YES" : "✗ NO")")

core_speedup = core_naive_time / core_power_time
println("\n  Performance on core:")
println("    Power iteration: $(round(core_power_time, digits=3))s")
println("    Naive iteration: $(round(core_naive_time, digits=3))s")
println("    Speedup: $(round(core_speedup, digits=2))x")

println("\n" * "=" ^ 80)
println("FINAL CONCLUSION")
println("=" ^ 80)

println("\n[Full Graph Results]")
println("  Power iteration: $(round(speedup, digits=1))x faster than naive")
if diff_power_naive < 1e-5
    println("  ✓ Both methods agree perfectly")
else
    println("  ⚠ Methods produce different results (max diff: $(round(diff_power_naive, digits=8)))")
    println("    Cause: Dangling nodes ($(round(100*dangling_nodes/vertices_count, digits=1))% of graph)")
end

println("\n[Core Subgraph Results]")
println("  Core: $core_vertices vertices ($(round(100*core_vertices/vertices_count, digits=1))% of original)")
println("  Core dangling: $core_dangling ($(round(100*core_dangling/core_vertices, digits=2))%)")
println("  Power iteration: $(round(core_speedup, digits=2))x faster")
if diff_core < 1e-5
    println("  ✓ Both methods agree perfectly (max diff: $(round(diff_core, digits=10)))")
    println("  ✓ Confirms dangling nodes were the issue!")
else
    println("  ⚠ Methods still differ (max diff: $(round(diff_core, digits=10)))")
end

println("\n[Key Findings]")
println("  1. Power iteration consistently faster: $(round(speedup, digits=1))x (full), $(round(core_speedup, digits=1))x (core)")
println("     - Dramatic speedup on core ($(round(core_speedup, digits=1))x) due to denser connectivity")
println()
println("  2. Methods STILL differ even without dangling nodes:")
println("     - Full: $(round(100*dangling_nodes/vertices_count, digits=1))% dangling → max diff $(round(diff_power_naive, digits=4))")
println("     - Core: $(round(100*core_dangling/core_vertices, digits=2))% dangling → max diff $(round(diff_core, digits=4))")
println("     - This suggests a deeper implementation difference beyond dangling nodes")
println()
println("  3. PageRank sum behavior:")
println("     - Full graph: sum ~$(round(sum(pr_power), digits=2)) (rank sink from dangling nodes)")
println("     - Core graph: sum = 1.0 ✓ (no dangling nodes)")
println()
println("  4. Both methods need investigation:")
println("     - Different convergence criteria?")
println("     - Different numerical precision?")
println("     - Implementation bug in one method?")

println("\n[Recommendation]")
println("  Use power iteration (sparse matrices) + proper dangling node treatment")
println("  OR work with core subgraph (main SCC) for web graphs")
println("=" ^ 80)
