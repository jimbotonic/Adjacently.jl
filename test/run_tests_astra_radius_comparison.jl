#!/usr/bin/env julia

#
# Test ASTRA-L with different ball radii on CNR-2000
# Compare radius=2 vs radius=3
#

include("run_tests_main.jl")

using Logging
using LightGraphs: nv, ne
using Adjacently

# Enable debug logging for ASTRA layered compression
global_logger(ConsoleLogger(stderr, Logging.Info))
using Adjacently.IO: load_adjacency_list_from_csv
using Adjacently.Compression.ASTRALayered: create_level_decomposition, estimate_astra_layered_cost_breakdown, write_astra_layered_graph
using Adjacently.PageRank: PR
using Adjacently.Graph: get_reverse_graph, get_sparse_P_matrix

println("=" ^ 80)
println("ASTRA-L Ball Radius Comparison Test")
println("Testing on CNR-2000 web graph")
println("Comparing radius=2 vs radius=3")
println("=" ^ 80)

# Load graph
println("\n[1/4] Loading CNR-2000 dataset...")
graph = load_adjacency_list_from_csv("datasets/webgraph/cnr-2000/cnr-2000.csv", ',', true)
vertices_count = nv(graph)
edges_count = ne(graph)
println("  Graph: $vertices_count vertices, $edges_count edges")

# Find best starting vertex (highest PageRank)
println("\n[2/4] Finding optimal BFS start vertex using PageRank...")

# Build reverse graph for PageRank computation
reverse_graph = get_reverse_graph(graph)

# Compute PageRank using the new get_sparse_P_matrix
pr_start = time()
P = get_sparse_P_matrix(graph)
pagerank_scores = PR(P; damping=0.85, epsilon=1e-6)
pr_time = time() - pr_start

# Find vertex with highest PageRank among high-degree vertices
min_degree_threshold = 100
max_pr = 0.0
start_vertex = UInt24(1)
candidates = 0
for v in 1:vertices_count
    deg = length(outneighbors(graph, UInt24(v)))
    if deg >= min_degree_threshold && pagerank_scores[v] > max_pr
        global max_pr = pagerank_scores[v]
        global start_vertex = UInt24(v)
        global candidates += 1
    end
end

println("  Start vertex: $start_vertex (PageRank: $(round(max_pr, digits=6)), degree: $(length(outneighbors(graph, start_vertex))))")

# Test both radii
radii = [2, 3]
results = Dict()

for radius in radii
    println("\n" * "=" ^ 80)
    println(" TESTING WITH BALL RADIUS = $radius")
    println("=" ^ 80)

    # Create level decomposition
    println("\n[3/4] Creating level decomposition (radius=$radius)...")
    decomp_start = time()
    decomp = create_level_decomposition(
        graph, reverse_graph;
        radius=radius,
        damping=0.85,
        epsilon=1e-6,
        pr=pagerank_scores
    )
    decomp_time = time() - decomp_start

    println("  Number of levels: $(decomp.num_levels)")
    println("  Decomposition time: $(round(decomp_time, digits=2))s")

    # Analyze edge distribution
    total_intra = sum(length(level.intra_edges) for level in decomp.levels)
    total_forward = sum(length(level.forward_edges) for level in decomp.levels)
    total_backward = sum(length(level.backward_edges) for level in decomp.levels)
    total_skip = sum(length(level.skip_edges) for level in decomp.levels)

    println("\n  Edge Distribution:")
    println("    Intra-level:  $total_intra ($(round(100*total_intra/edges_count, digits=1))%)")
    println("    Forward:      $total_forward ($(round(100*total_forward/edges_count, digits=1))%)")
    println("    Backward:     $total_backward ($(round(100*total_backward/edges_count, digits=1))%)")
    println("    Skip:         $total_skip ($(round(100*total_skip/edges_count, digits=1))%)")

    # Estimate compression cost
    println("\n[4/4] Estimating compression cost...")
    estimate_start = time()
    breakdown = estimate_astra_layered_cost_breakdown(decomp, :fibonacci)
    estimate_time = time() - estimate_start

    estimated_bits = breakdown["total_bits"]
    estimated_bits_per_edge = estimated_bits / edges_count

    println("  Estimated bits per edge: $(round(estimated_bits_per_edge, digits=3))")
    println("  Estimation time: $(round(estimate_time, digits=2))s")

    # Write compressed file
    output_dir = "test_data"
    mkpath(output_dir)
    output_file = joinpath(output_dir, "cnr2000_astra_radius_$(radius).astral")

    println("\n  Writing compressed file...")
    compress_start = time()
    filesize_bytes, actual_bpe = write_astra_layered_graph(decomp, output_file, reverse_graph, graph, :fibonacci)
    compress_time = time() - compress_start

    println("  File: $output_file")
    println("  Size: $(round(filesize_bytes / 1024 / 1024, digits=3)) MB")
    println("  Actual bits per edge: $(round(actual_bpe, digits=3))")
    println("  Compression time: $(round(compress_time, digits=2))s")

    # Store results
    results[radius] = Dict(
        "num_levels" => decomp.num_levels,
        "intra_edges" => total_intra,
        "forward_edges" => total_forward,
        "backward_edges" => total_backward,
        "skip_edges" => total_skip,
        "estimated_bpe" => estimated_bits_per_edge,
        "actual_bpe" => actual_bpe,
        "filesize_mb" => filesize_bytes / 1024 / 1024,
        "decomp_time" => decomp_time,
        "compress_time" => compress_time,
        "breakdown" => breakdown
    )
end

# Final comparison
println("\n" * "=" ^ 80)
println(" COMPARISON SUMMARY")
println("=" ^ 80)

println("\n┌─────────────────────────────┬──────────────┬──────────────┬────────────┐")
println("│ Metric                      │   Radius=2   │   Radius=3   │   Diff     │")
println("├─────────────────────────────┼──────────────┼──────────────┼────────────┤")

# Number of levels
r2_levels = results[2]["num_levels"]
r3_levels = results[3]["num_levels"]
diff_levels = r3_levels - r2_levels
diff_pct_levels = 100 * diff_levels / r2_levels
println("│ Number of levels            │ $(lpad(r2_levels, 12)) │ $(lpad(r3_levels, 12)) │ $(lpad(round(diff_pct_levels, digits=1), 6))%  │")

# Intra-level edges
r2_intra = results[2]["intra_edges"]
r3_intra = results[3]["intra_edges"]
diff_intra = r3_intra - r2_intra
diff_pct_intra = 100 * diff_intra / r2_intra
println("│ Intra-level edges           │ $(lpad(r2_intra, 12)) │ $(lpad(r3_intra, 12)) │ $(lpad(round(diff_pct_intra, digits=1), 6))%  │")

# Forward edges
r2_forward = results[2]["forward_edges"]
r3_forward = results[3]["forward_edges"]
diff_forward = r3_forward - r2_forward
diff_pct_forward = 100 * diff_forward / r2_forward
println("│ Forward edges               │ $(lpad(r2_forward, 12)) │ $(lpad(r3_forward, 12)) │ $(lpad(round(diff_pct_forward, digits=1), 6))%  │")

# Backward edges
r2_backward = results[2]["backward_edges"]
r3_backward = results[3]["backward_edges"]
diff_backward = r3_backward - r2_backward
diff_pct_backward = 100 * diff_backward / r2_backward
println("│ Backward edges              │ $(lpad(r2_backward, 12)) │ $(lpad(r3_backward, 12)) │ $(lpad(round(diff_pct_backward, digits=1), 6))%  │")

# Skip edges
r2_skip = results[2]["skip_edges"]
r3_skip = results[3]["skip_edges"]
diff_skip = r3_skip - r2_skip
diff_pct_skip = 100 * diff_skip / r2_skip
println("│ Skip edges                  │ $(lpad(r2_skip, 12)) │ $(lpad(r3_skip, 12)) │ $(lpad(round(diff_pct_skip, digits=1), 6))%  │")

println("├─────────────────────────────┼──────────────┼──────────────┼────────────┤")

# Estimated bits per edge
r2_est_bpe = results[2]["estimated_bpe"]
r3_est_bpe = results[3]["estimated_bpe"]
diff_est_bpe = r3_est_bpe - r2_est_bpe
diff_pct_est_bpe = 100 * diff_est_bpe / r2_est_bpe
println("│ Estimated bits/edge         │ $(lpad(round(r2_est_bpe, digits=3), 12)) │ $(lpad(round(r3_est_bpe, digits=3), 12)) │ $(lpad(round(diff_pct_est_bpe, digits=1), 6))%  │")

# Actual bits per edge
r2_bpe = results[2]["actual_bpe"]
r3_bpe = results[3]["actual_bpe"]
diff_bpe = r3_bpe - r2_bpe
diff_pct_bpe = 100 * diff_bpe / r2_bpe
println("│ Actual bits/edge            │ $(lpad(round(r2_bpe, digits=3), 12)) │ $(lpad(round(r3_bpe, digits=3), 12)) │ $(lpad(round(diff_pct_bpe, digits=1), 6))%  │")

# File size
r2_size = results[2]["filesize_mb"]
r3_size = results[3]["filesize_mb"]
diff_size = r3_size - r2_size
diff_pct_size = 100 * diff_size / r2_size
println("│ File size (MB)              │ $(lpad(round(r2_size, digits=3), 12)) │ $(lpad(round(r3_size, digits=3), 12)) │ $(lpad(round(diff_pct_size, digits=1), 6))%  │")

println("├─────────────────────────────┼──────────────┼──────────────┼────────────┤")

# Decomposition time
r2_decomp = results[2]["decomp_time"]
r3_decomp = results[3]["decomp_time"]
diff_decomp = r3_decomp - r2_decomp
diff_pct_decomp = 100 * diff_decomp / r2_decomp
println("│ Decomposition time (s)      │ $(lpad(round(r2_decomp, digits=2), 12)) │ $(lpad(round(r3_decomp, digits=2), 12)) │ $(lpad(round(diff_pct_decomp, digits=1), 6))%  │")

# Compression time
r2_comp = results[2]["compress_time"]
r3_comp = results[3]["compress_time"]
diff_comp = r3_comp - r2_comp
diff_pct_comp = 100 * diff_comp / r2_comp
println("│ Compression time (s)        │ $(lpad(round(r2_comp, digits=2), 12)) │ $(lpad(round(r3_comp, digits=2), 12)) │ $(lpad(round(diff_pct_comp, digits=1), 6))%  │")

println("└─────────────────────────────┴──────────────┴──────────────┴────────────┘")

# Determine winner
println("\n CONCLUSION:")
if r2_bpe < r3_bpe
    improvement = 100 * (r3_bpe - r2_bpe) / r3_bpe
    println("  Radius=2 is BETTER by $(round(improvement, digits=1))%")
    println("  $(round(r2_bpe, digits=3)) bits/edge vs $(round(r3_bpe, digits=3)) bits/edge")
elseif r3_bpe < r2_bpe
    improvement = 100 * (r2_bpe - r3_bpe) / r2_bpe
    println("  Radius=3 is BETTER by $(round(improvement, digits=1))%")
    println("  $(round(r3_bpe, digits=3)) bits/edge vs $(round(r2_bpe, digits=3)) bits/edge")
else
    println("  Both radii achieve the SAME compression!")
    println("  $(round(r2_bpe, digits=3)) bits/edge")
end

println("\n  Standard ASTRA baseline: 5.108 bits/edge")
println("  Best ASTRA-L result: $(round(min(r2_bpe, r3_bpe), digits=3)) bits/edge")
best_improvement = 100 * (5.108 - min(r2_bpe, r3_bpe)) / 5.108
println("  Improvement over baseline: $(round(best_improvement, digits=1))%")

println("\n" * "=" ^ 80)
