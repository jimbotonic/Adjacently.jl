#!/usr/bin/env julia

#
# Test ASTRA-L (Layered) compression on CNR-2000
#

include("run_tests_main.jl")

using Logging
using LightGraphs: nv, ne
using Adjacently

# Enable debug logging for ASTRA layered compression
global_logger(ConsoleLogger(stderr, Logging.Debug))
using Adjacently.IO: load_adjacency_list_from_csv
using Adjacently.Compression.ASTRALayered: create_level_decomposition, create_level_decomposition_ppr, create_level_decomposition_ppr_bfs, estimate_astra_layered_cost, estimate_astra_layered_cost_breakdown, write_astra_layered_graph
using Adjacently.PageRank: PR
using Adjacently.Graph: get_reverse_graph

println("=" ^ 80)
println("ASTRA-L (Layered Graph Compression) Test")
println("Testing on CNR-2000 web graph")
println("=" ^ 80)

# Load graph
println("\n Loading CNR-2000 dataset...")
graph = load_adjacency_list_from_csv("datasets/webgraph/cnr-2000/cnr-2000.csv", ',', true)
vertices_count = nv(graph)
edges_count = ne(graph)
println("  Graph: $vertices_count vertices, $edges_count edges")

# Find best starting vertex (highest PageRank)
println("\n Finding optimal BFS start vertex using PageRank...")
println("  Computing PageRank scores...")

# Build reverse graph for PageRank computation
reverse_graph = get_reverse_graph(graph)

# Compute PageRank
pr_start = time()
pagerank_scores = PR(graph, reverse_graph)
pr_time = time() - pr_start

# Find vertex with highest PageRank among high-degree vertices
# Strategy: Only consider vertices with degree >= 100
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

println("  PageRank computation time: $(round(pr_time, digits=2))s")
println("  Candidate vertices (degree >= $min_degree_threshold): $candidates")
println("  Start vertex: $start_vertex (PageRank: $(round(max_pr, digits=6)))")
println("  Degree of start vertex: $(length(outneighbors(graph, start_vertex)))")

# Create level decomposition (BFS radius-2 + PPR filtering)
println("\n" * "=" ^ 80)
println(" PHASE 1: LEVEL DECOMPOSITION (BFS radius-2 + PPR filtering)")
println("=" ^ 80)

decomp_start = time()
# BFS radius-2 + PPR filtering approach with adaptive radius:
# - Explores 2-hop neighborhood from seed
# - Adaptively expands to radius-3, 4, 5 if ball too small
# - Filters using PPR to maintain community cohesion
# - Immediately absorbs sink vertices (out-degree 0)
# - Selects next seed from high-PPR frontier vertices
#
# PARAMETERS FOR FINE-GRAINED LEVELS (many small 2-radius balls):
# - min_ball_size: 50 (small communities for compact local IDs)
# - min_ppr_score: 1e-6 (standard)
# - radius: 2 (2-radius BFS balls centered on seeds)
# Goal: Test baseline with smaller communities
decomp = create_level_decomposition_ppr_bfs(graph, reverse_graph, start_vertex,
                                             alpha=0.15, epsilon=1e-6,
                                             min_ppr_score=1e-6, radius=2,
                                             min_ball_size=50, max_radius=5)
decomp_time = time() - decomp_start

println("\n Decomposition Statistics:")
println("  Number of levels: $(decomp.num_levels)")
println("  Total vertices: $(decomp.total_vertices)")
println("  Total edges: $(decomp.total_edges)")
println("  Decomposition time: $(round(decomp_time, digits=2))s")

# Analyze edge distribution
total_intra = sum(length(level.intra_edges) for level in decomp.levels)
total_forward = sum(length(level.forward_edges) for level in decomp.levels)
total_backward = sum(length(level.backward_edges) for level in decomp.levels)
total_skip = sum(length(level.skip_edges) for level in decomp.levels)

println("\n Edge Distribution:")
println("  Intra-level:  $total_intra ($(round(100*total_intra/edges_count, digits=1))%)")
println("  Forward:      $total_forward ($(round(100*total_forward/edges_count, digits=1))%)")
println("  Backward:     $total_backward ($(round(100*total_backward/edges_count, digits=1))%)")
println("  Skip:         $total_skip ($(round(100*total_skip/edges_count, digits=1))%)")

# Analyze local ID ranges per level
println("\n Local ID Ranges per Level:")
for (idx, level) in enumerate(decomp.levels[1:min(10, length(decomp.levels))])
    max_bits_needed = ceil(Int, log2(max(1, level.num_vertices)))
    println("  Level $(level.level_id): $(level.num_vertices) vertices (needs $max_bits_needed bits)")
end
if length(decomp.levels) > 10
    println("  ... ($(length(decomp.levels) - 10) more levels)")
end

# Estimate compression cost with detailed breakdown
println("\n" * "=" ^ 80)
println(" PHASE 2: COST ESTIMATION & BREAKDOWN")
println("=" ^ 80)

estimate_start = time()
breakdown = estimate_astra_layered_cost_breakdown(decomp, :fibonacci)
estimate_time = time() - estimate_start

estimated_bits = breakdown["total_bits"]
estimated_bits_per_edge = estimated_bits / edges_count
estimated_bytes = estimated_bits / 8

println("\n Estimated Compression (Fibonacci encoding):")
println("  Total bits: $estimated_bits")
println("  Bits per edge: $(round(estimated_bits_per_edge, digits=3))")
println("  Estimated file size: $(round(estimated_bytes / 1024 / 1024, digits=3)) MB")
println("  Estimation time: $(round(estimate_time, digits=2))s")

println("\n Detailed Cost Breakdown by Edge Type:")
println("  ─" ^ 80)
println("  Edge Type    │  Count   │  Bits      │  Bits/Edge  │  % of Total Bits")
println("  ─" ^ 80)

# Header
header_pct = 100 * breakdown["header_bits"] / estimated_bits
println("  Header       │     -    │  $(lpad(breakdown["header_bits"], 9)) │      -      │  $(lpad(round(header_pct, digits=1), 5))%")

# Intra-level
if breakdown["intra_count"] > 0
    intra_bpe = breakdown["intra_bits"] / breakdown["intra_count"]
    intra_pct = 100 * breakdown["intra_bits"] / estimated_bits
    println("  Intra-level  │  $(lpad(breakdown["intra_count"], 6)) │  $(lpad(breakdown["intra_bits"], 9)) │    $(lpad(round(intra_bpe, digits=2), 5))    │  $(lpad(round(intra_pct, digits=1), 5))%")
end

# Forward
if breakdown["forward_count"] > 0
    forward_bpe = breakdown["forward_bits"] / breakdown["forward_count"]
    forward_pct = 100 * breakdown["forward_bits"] / estimated_bits
    println("  Forward      │  $(lpad(breakdown["forward_count"], 6)) │  $(lpad(breakdown["forward_bits"], 9)) │    $(lpad(round(forward_bpe, digits=2), 5))    │  $(lpad(round(forward_pct, digits=1), 5))%")
end

# Backward
if breakdown["backward_count"] > 0
    backward_bpe = breakdown["backward_bits"] / breakdown["backward_count"]
    backward_pct = 100 * breakdown["backward_bits"] / estimated_bits
    println("  Backward     │  $(lpad(breakdown["backward_count"], 6)) │  $(lpad(breakdown["backward_bits"], 9)) │    $(lpad(round(backward_bpe, digits=2), 5))    │  $(lpad(round(backward_pct, digits=1), 5))%")
end

# Skip
if breakdown["skip_count"] > 0
    skip_bpe = breakdown["skip_bits"] / breakdown["skip_count"]
    skip_pct = 100 * breakdown["skip_bits"] / estimated_bits
    println("  Skip         │  $(lpad(breakdown["skip_count"], 6)) │  $(lpad(breakdown["skip_bits"], 9)) │    $(lpad(round(skip_bpe, digits=2), 5))    │  $(lpad(round(skip_pct, digits=1), 5))%")
end

println("  ─" ^ 80)
println("  TOTAL        │  $(lpad(breakdown["total_edges"], 6)) │  $(lpad(estimated_bits, 9)) │    $(lpad(round(estimated_bits_per_edge, digits=2), 5))    │  100.0%")
println("  ─" ^ 80)

# Compare to standard ASTRA
standard_astra_bpe = 5.108  # From previous tests
improvement = 100 * (standard_astra_bpe - estimated_bits_per_edge) / standard_astra_bpe
println("\n Comparison to Standard ASTRA:")
println("  Standard ASTRA: $standard_astra_bpe bits/edge")
println("  ASTRA-L (estimated): $(round(estimated_bits_per_edge, digits=3)) bits/edge")
println("  Improvement: $(round(improvement, digits=1))%")

# Write compressed file
println("\n" * "=" ^ 80)
println(" PHASE 3: COMPRESSION")
println("=" ^ 80)

output_dir = "test_data"
mkpath(output_dir)
output_file = joinpath(output_dir, "cnr2000_astra_layered_ppr_bfs_adaptive.astral")

compress_start = time()
filesize_bytes, actual_bpe = write_astra_layered_graph(decomp, output_file, reverse_graph, graph, :fibonacci)
compress_time = time() - compress_start

println("\n Compression Results:")
println("  File: $output_file")
println("  Size: $(round(filesize_bytes / 1024 / 1024, digits=3)) MB")
println("  Bits per edge: $(round(actual_bpe, digits=3))")
println("  Compression time: $(round(compress_time, digits=2))s")
println("  Compression ratio: $(round(edges_count / filesize_bytes, digits=2)) edges/byte")

# Final summary
println("\n" * "=" ^ 80)
println(" SUMMARY")
println("=" ^ 80)

println("\n ASTRA-L Compression Complete!")
println("  Graph: CNR-2000 ($vertices_count vertices, $edges_count edges)")
println("  Levels: $(decomp.num_levels)")
println("  File size: $(round(filesize_bytes / 1024 / 1024, digits=3)) MB")
println("  Bits per edge: $(round(actual_bpe, digits=3))")
println("")
println("  Comparison:")
println("    Standard ASTRA: 5.108 bits/edge")
println("    ASTRA-L:        $(round(actual_bpe, digits=3)) bits/edge")
if actual_bpe < standard_astra_bpe
    println("    Improvement:    $(round(improvement, digits=1))% better! 🎉")
else
    println("    Note: ASTRA-L needs optimization")
end

println("\n" * "=" ^ 80)
