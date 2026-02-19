#!/usr/bin/env julia

#
# Debug PageRank discrepancy between naive and matrix methods
#

include("run_tests_main.jl")

using LightGraphs: nv, ne, outneighbors, add_edge!, add_vertices!
using SparseArrays: SparseMatrixCSC, sparse
using Adjacently
using Adjacently.CustomLightGraphs: SimpleDiGraph
using Adjacently.IO: load_adjacency_list_from_csv
using Adjacently.PageRank: PR
using Adjacently.Graph: get_reverse_graph

println("=" ^ 80)
println("Debugging PageRank Discrepancy")
println("=" ^ 80)

# Test on a small graph first to understand the issue
println("\n[Test 1] Small graph test (5 nodes)")
println("-" ^ 80)

# Create a simple test graph: 1->2, 2->3, 3->1, 4->1, 5 (dangling)
g_test = SimpleDiGraph{UInt24}()
add_vertices!(g_test, 5)
add_edge!(g_test, UInt24(1), UInt24(2))
add_edge!(g_test, UInt24(2), UInt24(3))
add_edge!(g_test, UInt24(3), UInt24(1))
add_edge!(g_test, UInt24(4), UInt24(1))
# Node 5 is dangling (no outgoing edges)

rg_test = get_reverse_graph(g_test)
n_test = nv(g_test)

# Build matrix
I_test = Int[]
J_test = Int[]
V_test = Float64[]
for j in 1:n_test
    jv = UInt24(j)
    outs = outneighbors(g_test, jv)
    d = length(outs)
    if d == 0
        continue  # Skip dangling nodes - THIS IS THE KEY ISSUE!
    end
    w = 1.0 / d
    for i in outs
        push!(I_test, Int(i))
        push!(J_test, j)
        push!(V_test, w)
    end
end
P_test_int = sparse(I_test, J_test, V_test, n_test, n_test)
P_test = SparseMatrixCSC{Float64, UInt24}(P_test_int)

damping = 0.85

# Compute with both methods
pr_naive_test = PR(g_test, rg_test; damping=damping, epsilon=1e-10)
pr_matrix_test = PR(P_test; damping=damping, epsilon=1e-10)

println("Graph structure:")
for v in 1:n_test
    outs = outneighbors(g_test, UInt24(v))
    println("  Node $v -> $(collect(outs)) (outdegree = $(length(outs)))")
end

println("\nNaive method:")
println("  PR = $pr_naive_test")
println("  Sum = $(sum(pr_naive_test))")

println("\nMatrix method:")
println("  PR = $pr_matrix_test")
println("  Sum = $(sum(pr_matrix_test))")

println("\nDifference:")
println("  Max diff = $(maximum(abs.(pr_naive_test - pr_matrix_test)))")

# Now let's analyze the issue
println("\n" ^ 2 * "=" ^ 80)
println("ANALYSIS OF THE DISCREPANCY")
println("=" ^ 80)

println("\n[Root Cause]")
println("-" ^ 80)
println("The discrepancy comes from how dangling nodes are handled:")
println()
println("1. NAIVE METHOD (line 58 in pr.jl):")
println("   nv += pr[p]/length(outneighbors(g,p))")
println()
println("   When a node p has NO outgoing edges (dangling):")
println("   - length(outneighbors(g,p)) = 0")
println("   - This causes pr[p]/0 = Inf or NaN")
println("   - Julia may skip this in the sum, effectively treating it as 0")
println()
println("2. MATRIX METHOD (line 215 in pr.jl):")
println("   pr2 = damping*P'*pr + (1-damping)*ppr")
println()
println("   When building matrix P, dangling nodes are skipped:")
println("   - if d == 0: continue")
println("   - Column j of P is ZERO for dangling node j")
println("   - This means dangling nodes contribute 0 to all other nodes")
println()
println("3. THE PROBLEM:")
println("   - Dangling nodes accumulate PageRank but don't distribute it")
println("   - This creates a \"rank sink\" - PageRank mass disappears")
println("   - Sum of PageRank < 1.0")
println()

println("\n[Testing with CNR-2000]")
println("-" ^ 80)

# Load CNR-2000
graph = load_adjacency_list_from_csv("datasets/webgraph/cnr-2000/cnr-2000.csv", ',', true)
n = nv(graph)

# Count dangling nodes
dangling = sum(length(outneighbors(graph, UInt24(v))) == 0 for v in 1:n)
dangling_pct = 100.0 * dangling / n

println("CNR-2000: $n nodes, $dangling dangling ($(round(dangling_pct, digits=2))%)")
println()
println("With $(round(dangling_pct, digits=1))% dangling nodes, we expect:")
println("  - PageRank sum ≈ $(round(1.0 - damping * dangling_pct/100, digits=3))")
println("  - Observed sum ≈ 0.71 (close!)")

println("\n[Solutions]")
println("-" ^ 80)
println("To fix this, we need to handle dangling nodes properly:")
println()
println("Option 1: Redistribute dangling mass")
println("  After each iteration, collect mass from dangling nodes")
println("  and redistribute uniformly (or according to ppr)")
println()
println("Option 2: Add self-loops to dangling nodes")
println("  Make dangling nodes point to themselves")
println()
println("Option 3: Teleport from dangling nodes")
println("  Treat dangling nodes as if they teleport uniformly")
println()
println("The standard PageRank algorithm uses Option 1 or 3.")

println("\n" * "=" ^ 80)
