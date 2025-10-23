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

# Test CNR-2000 with optimal compression settings:
# - children mode
# - fibonacci encoding  
# - RCM ordering (out-degree variant)
# - reference encoding enabled
# - mix encoding with run-length enabled
# - interval compression enabled

using Pkg
Pkg.activate(".")

using Test
using LightGraphs: nv, ne, outneighbors
using Adjacently.IO: load_adjacency_list_from_csv
using Adjacently.Relabeling: relabel_vertices_rcm, relabel_graph
using Adjacently.MGS: write_compressed_mgs3_graph, load_compressed_mgs3_graph

@testset "CNR-2000 Optimal Compression Test" begin
    # Load CNR-2000 dataset
    cnr_csv = "datasets/webgraph/cnr-2000/cnr-2000.csv"
    
    if !isfile(cnr_csv)
        @warn "CNR-2000 dataset not found at $cnr_csv, skipping test"
        return
    end
    
    println("=== CNR-2000 Optimal Compression Test ===")
    
    # Load graph
    println("Loading CNR-2000...")
    g = load_adjacency_list_from_csv(cnr_csv, ',', true)
    n, m = nv(g), ne(g)
    println("Graph loaded: n=$n, m=$m")
    
    # Apply RCM ordering (out-degree variant for best compression)
    println("Applying RCM ordering (out-degree)...")
    rcm_map = relabel_vertices_rcm(g, :out)
    g_rcm = relabel_graph(g, rcm_map)
    
    # Create test_data directory if it doesn't exist
    mkpath("test_data")
    
    # Test optimal compression settings
    output_file = "test_data/cnr2k_optimal_hybrid"
    
    println("Testing optimal compression: children+fibonacci+rcm+ref+mix+intervals...")
    
    # Compress with optimal settings
    @time write_compressed_mgs3_graph(
        g_rcm, 
        output_file, 
        :children,    # children mode
        :fibonacci,   # fibonacci encoding
        true,         # mix encoding enabled
        true          # reference encoding enabled
    )
    
    # Check file was created
    mgz_file = "$(output_file).mgz"
    @test isfile(mgz_file)
    
    # Get compression statistics
    filesize = stat(mgz_file).size
    bpv = (filesize * 8) / n  # bits per vertex
    bpe = (filesize * 8) / m  # bits per edge
    
    println("\nCompression Results:")
    println("File size: $(round(filesize/1024/1024; digits=3)) MB")
    println("Bits per vertex: $(round(bpv; digits=3))")
    println("Bits per edge: $(round(bpe; digits=3))")
    
    # WebGraph reference for comparison
    webgraph_size = 1164843  # bytes (from WebGraph CNR-2000)
    webgraph_bpe = 2.897     # bits per edge
    gap_to_webgraph = (filesize / webgraph_size - 1) * 100
    
    println("\nComparison to WebGraph:")
    println("WebGraph size: $(round(webgraph_size/1024/1024; digits=3)) MB")
    println("WebGraph BPE: $webgraph_bpe")
    println("Gap: $(round(gap_to_webgraph; digits=1))% $(gap_to_webgraph > 0 ? "larger" : "smaller")")
    
    # Decompress and verify
    println("\nTesting decompression...")
    @time g_decoded = load_compressed_mgs3_graph(mgz_file)
    
    # Verify graph properties
    @test nv(g_decoded) == n
    @test ne(g_decoded) == m
    
    println("Decompressed graph: n=$(nv(g_decoded)), m=$(ne(g_decoded))")
    
    # Verify neighbor lists match (full graph verification)
    println("Verifying neighbor lists for full graph...")
    mismatch_count = 0
    
    for v in 1:n
        original_neighbors = sort(collect(outneighbors(g_rcm, v)))
        decoded_neighbors = sort(collect(outneighbors(g_decoded, v)))
        
        if original_neighbors != decoded_neighbors
            mismatch_count += 1
            if mismatch_count <= 5  # Show first 5 mismatches
                println("Mismatch at vertex $v:")
                println("  Original: $(length(original_neighbors)) neighbors")
                println("  Decoded:  $(length(decoded_neighbors)) neighbors")
            end
        end
    end
    
    @test mismatch_count == 0
    
    if mismatch_count == 0
        println(" All neighbor lists match perfectly!")
    else
        println(" Found $mismatch_count mismatched vertices out of $n tested")
    end
    
    # Performance assessment
    println("\nPerformance Assessment:")
    if gap_to_webgraph < 10
        println(" Excellent! Very close to WebGraph performance")
    elseif gap_to_webgraph < 50
        println(" Good performance, within 50% of WebGraph")
    elseif gap_to_webgraph < 100
        println(" Decent compression, within 100% of WebGraph")
    else
        println(" Need further optimization to reach WebGraph performance")
        println("   Consider: copy-list references, advanced gap encoding")
    end
    
    # Save compressed file in test_data folder
    println("\nCompressed file saved as: $mgz_file")
    println("File size: $(round(filesize/1024/1024; digits=3)) MB")
    
    println("\n✓ Test completed successfully!")
end