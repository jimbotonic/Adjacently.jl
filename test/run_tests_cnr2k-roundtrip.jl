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

include("run_tests_main.jl")

#!/usr/bin/env julia

# Complete round-trip test: compress and decompress CNR-2000
# Verifies that encoding/decoding preserves the graph structure

#using Pkg
#Pkg.activate(@__DIR__)

using LightGraphs: nv, ne, outneighbors
using Adjacently
using Adjacently.IO: load_adjacency_list_from_csv
using Adjacently.Relabeling: relabel_vertices_rcm, relabel_vertices_llp, relabel_vertices_bfs, relabel_graph
using Adjacently.MGS: write_compressed_mgs3_graph, load_compressed_mgs3_graph

println("=" ^ 80)
println("CNR-2000 Complete Round-Trip Test with RCM Ordering")
println("=" ^ 80)

# Load and relabel
println("\n Loading CNR-2000 dataset...")
graph_original = load_adjacency_list_from_csv("datasets/webgraph/cnr-2000/cnr-2000.csv", ',', true)
vertices_count = nv(graph_original)
edges_count = ne(graph_original)
println("  Graph: $vertices_count vertices, $edges_count edges")

println("\n🔄 Applying RCM (Reverse Cuthill-McKee) relabeling...")
rcm_mapping = relabel_vertices_rcm(graph_original, :out)
graph = relabel_graph(graph_original, rcm_mapping)
println("  RCM applied successfully")

# Compress
output_dir = "test_data"
mkpath(output_dir)
output_file = joinpath(output_dir, "cnr2000_roundtrip_test")

println("\n" * "=" ^ 80)
println(" STEP 1: COMPRESSION")
println("=" ^ 80)

println("\n Compressing CNR-2000...")
println("  Using RCM ordering + Fibonacci encoding + adaptive block/RLE/raw encoding")
compress_start = time()

write_compressed_mgs3_graph(
    graph,
    output_file,
    :children,                # coding_scheme
    :complex,                 # compression
    :fibonacci,               # integer_encoding (Fibonacci performs best)
    true,                     # use_mix_mode
    true,                     # reference_enabled
    true,                     # recursive_reference
    1024                      # ref_window_size
)

compress_time = time() - compress_start
# The writer adds .mgz extension automatically
output_file_actual = output_file * ".mgz"
file_size = filesize(output_file_actual)
bits_per_edge = 8.0 * file_size / edges_count

println("\n Compression complete!")
println("  Time: $(round(compress_time, digits=2))s")
println("  File size: $(round(file_size / (1024*1024), digits=3)) MB")
println("  Bits per edge: $(round(bits_per_edge, digits=3))")

# Decompress
println("\n" * "=" ^ 80)
println("📂 STEP 2: DECOMPRESSION")
println("=" ^ 80)

println("\n  Decompressing CNR-2000...")
decompress_start = time()

graph_restored = load_compressed_mgs3_graph(output_file_actual)

decompress_time = time() - decompress_start

println("\n Decompression complete!")
println("  Time: $(round(decompress_time, digits=2))s")

# Verify
println("\n" * "=" ^ 80)
println(" STEP 3: VERIFICATION")
println("=" ^ 80)

println("\n Graph statistics:")
println("  Original vertices: $vertices_count")
println("  Restored vertices: $(nv(graph_restored))")
println("  Original edges:    $edges_count")
println("  Restored edges:    $(ne(graph_restored))")

# Check basic stats
vertices_match = nv(graph) == nv(graph_restored)
edges_match = ne(graph) == ne(graph_restored)

println("\n🔬 Basic checks:")
println("  Vertices match: $(vertices_match ? "OK" : "Error")")
println("  Edges match:    $(edges_match ? "OK" : "Error")")

if !vertices_match || !edges_match
    println("\n FAILED: Graph structure doesn't match!")
    exit(1)
end

# Detailed verification: compare all adjacency lists
println("\n Detailed verification (comparing all adjacency lists)...")
verify_start = time()

mismatches = 0
mismatch_details = []

for v in 1:nv(graph)
    global mismatches, mismatch_details
    original_neighbors = sort(collect(outneighbors(graph, v)))
    restored_neighbors = sort(collect(outneighbors(graph_restored, v)))

    if original_neighbors != restored_neighbors
        mismatches += 1
        if mismatches <= 10  # Only store first 10 mismatches for reporting
            push!(mismatch_details, (
                vertex = v,
                original_degree = length(original_neighbors),
                restored_degree = length(restored_neighbors),
                original_sample = original_neighbors[1:min(5, length(original_neighbors))],
                restored_sample = restored_neighbors[1:min(5, length(restored_neighbors))]
            ))
        end
    end

    # Progress indicator
    if v % 50000 == 0
        progress = 100.0 * v / nv(graph)
        println("    Progress: $(round(progress, digits=1))% ($v / $(nv(graph)) vertices)")
    end
end

verify_time = time() - verify_start

println("\n Verification time: $(round(verify_time, digits=2))s")

# Report results
println("\n" * "=" ^ 80)
println(" FINAL RESULTS")
println("=" ^ 80)

if mismatches == 0
    println("\n SUCCESS! Perfect round-trip verification!")
    println("   All $vertices_count vertex adjacency lists match exactly.")
    println("\n Performance Summary:")
    println("   Compression:   $(round(compress_time, digits=2))s")
    println("   Decompression: $(round(decompress_time, digits=2))s")
    println("   Verification:  $(round(verify_time, digits=2))s")
    println("   Total time:    $(round(compress_time + decompress_time + verify_time, digits=2))s")
    println("\n Compression Stats:")
    println("   File size:     $(round(file_size / (1024*1024), digits=3)) MB")
    println("   Bits per edge: $(round(bits_per_edge, digits=3))")
    println("   Compression ratio: $(round(edges_count / file_size, digits=2)) edges/byte")
else
    println("\n FAILED! Found $mismatches vertex mismatches (out of $vertices_count vertices)")
    println("   This means $(round(100.0 * mismatches / vertices_count, digits=2))% of vertices have incorrect adjacency lists")

    if !isempty(mismatch_details)
        println("\n🔍 Sample mismatches (first 10):")
        for (i, detail) in enumerate(mismatch_details)
            println("\n   Mismatch #$i - Vertex $(detail.vertex):")
            println("     Original degree: $(detail.original_degree)")
            println("     Restored degree: $(detail.restored_degree)")
            println("     Original sample: $(detail.original_sample)")
            println("     Restored sample: $(detail.restored_sample)")
        end
    end

    exit(1)
end

println("\n" * "=" ^ 80)
println(" Round-trip test complete!")
println("=" ^ 80)
