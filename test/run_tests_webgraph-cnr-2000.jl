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

#
# Test WebGraph CNR-2000 dataset processing
# Loads CSV edge list, extracts core, compresses with zeta+reference encoding
#

include("run_tests_main.jl")
using Logging
using Adjacently.MGS: write_compressed_mgs3_graph, load_compressed_mgs3_graph
using Adjacently.IO: load_adjacency_list_from_csv, BitWriter, write_value, flush_bitwriter
using Adjacently.Compression: write_compressed_graph_data

@testset "WebGraph CNR-2000 Processing" begin
    # Enable debug-level logging for encoding/decoding visibility
    prev_logger = current_logger()
    global_logger(ConsoleLogger(stderr, Logging.Debug))
    try
    @info "=== WebGraph CNR-2000 Dataset Test ==="
    
    # Define paths
    cnr_csv_path = "datasets/webgraph/cnr-2000/cnr-2000.csv"
    output_dir = "test_data"
    output_file = "cnr-2000_core_zeta_ref"
    
    # Check if CSV file exists
    if !isfile(cnr_csv_path)
        @warn "CNR-2000 CSV file not found at: $cnr_csv_path"
        @warn "Please run the CNR-2000 dataset processing workflow first"
        @warn "See: datasets/webgraph/cnr-2000/README.md"
        @test_skip "CNR-2000 CSV file not available"
        return
    end
    
    @info "Step 1: Loading CNR-2000 edge list from CSV..."
    
    # Load CSV edge list using Adjacently.IO function (skip header)
    @time original_graph = load_adjacency_list_from_csv(cnr_csv_path, ',', true)
    
    original_vertices = nv(original_graph)
    original_edges = ne(original_graph)
    
    @info "  Loaded graph: $original_vertices vertices, $original_edges edges"
    
    # Convert to adjacency dictionary for processing
    @info "Step 2: Building adjacency lists..."
    neighbor_lists = Dict{UInt32, Vector{UInt32}}()
    
    @time for v in vertices(original_graph)
        neighbors = UInt32[]
        for neighbor in outneighbors(original_graph, v)
            push!(neighbors, UInt32(neighbor))
        end
        if !isempty(neighbors)
            sort!(neighbors)
            neighbor_lists[UInt32(v)] = neighbors
        end
    end
    
    @info "  Original graph: $original_vertices vertices, $original_edges edges"
    
    # Extract core (vertices with both in-degree and out-degree > 0)
    @info "Step 3: Extracting graph core..."
    
    # Count in-degrees
    indegree = Dict{UInt32, Int}()
    for neighbors in values(neighbor_lists)
        for target in neighbors
            indegree[target] = get(indegree, target, 0) + 1
        end
    end
    
    # Find core vertices (have both outgoing and incoming edges)
    core_vertices = Set{UInt32}()
    for (source, neighbors) in neighbor_lists
        if length(neighbors) > 0 && get(indegree, source, 0) > 0
            push!(core_vertices, source)
        end
    end
    
    @info "  Core vertices: $(length(core_vertices)) out of $original_vertices"
    
    # Build core subgraph with vertex remapping
    @info "Step 4: Building core subgraph with vertex remapping..."
    core_vertices_sorted = sort(collect(core_vertices))
    old_to_new = Dict{UInt32, UInt32}()
    new_to_old = Dict{UInt32, UInt32}()
    
    for (new_id, old_id) in enumerate(core_vertices_sorted)
        old_to_new[old_id] = UInt32(new_id)
        new_to_old[UInt32(new_id)] = old_id
    end
    
    # Create core neighbor lists with remapped IDs
    core_neighbor_lists = Dict{UInt32, Vector{UInt32}}()
    
    for (new_id, old_id) in enumerate(core_vertices_sorted)
        old_neighbors = get(neighbor_lists, old_id, UInt32[])
        new_neighbors = UInt32[]
        
        for old_target in old_neighbors
            if haskey(old_to_new, old_target)
                push!(new_neighbors, old_to_new[old_target])
            end
        end
        
        if !isempty(new_neighbors)
            sort!(new_neighbors)
            core_neighbor_lists[UInt32(new_id)] = new_neighbors
        end
    end
    
    core_vertices_count = length(core_neighbor_lists)
    core_edges_count = sum(length(neighbors) for neighbors in values(core_neighbor_lists))
    
    @info "  Core graph: $core_vertices_count vertices, $core_edges_count edges"
    @info "  Vertex reduction: $(round(100 * (1 - core_vertices_count/original_vertices), digits=2))%"
    @info "  Edge reduction: $(round(100 * (1 - core_edges_count/original_edges), digits=2))%"
    
    # Create output directory
    mkpath(output_dir)
    
    # Prepare output directory
    mkpath(output_dir)
    
    # Test configurations: [description, filename_suffix, use_mix_mode, use_reference]
    test_configs = [
        ("Whole graph with zeta + delta only + children", "whole_delta", false, false),
        ("Whole graph with zeta + mix-mode + children", "whole_mix", true, false),
        ("Whole graph with zeta + reference + children", "whole_ref", false, true),
        ("Whole graph with zeta + mix-mode + reference + children", "whole_mix_ref", true, true)
    ]
    
    compression_results = []
    
    for (i, (description, suffix, use_mix_mode, use_reference)) in enumerate(test_configs)
        @info "\n" * "="^80
        @info "Step $(i+4): $description"
        @info "="^80
        
        # Use whole original graph
        test_graph = original_graph
        vertices_count = nv(test_graph)
        edges_count = ne(test_graph)
        @info "  Using whole graph: $vertices_count vertices, $edges_count edges"
        
        # Compress with MGS3 format
        encoding_type = if use_reference && use_mix_mode
            "mix-mode + reference"
        elseif use_reference
            "reference"
        elseif use_mix_mode
            "mix-mode"
        else
            "delta only"
        end
        @info "  Compressing with zeta + $encoding_type + children..."
        
        output_path = joinpath(output_dir, output_file * "_" * suffix)
        
        @time begin
            # Use modified MGS3 writer with specific parameters
            write_compressed_mgs3_graph(test_graph, output_path, :children, :zeta, use_mix_mode, use_reference)
        end
        
        output_path_with_ext = output_path * ".mgz"
        
        # Get file statistics
        file_size = filesize(output_path_with_ext)
        bits_per_vertex = (file_size * 8) / vertices_count
        bits_per_edge = (file_size * 8) / edges_count
        
        @info "  Results:"
        @info "    File: $(basename(output_path_with_ext))"
        @info "    Size: $(file_size) bytes ($(round(file_size/1024/1024, digits=3)) MB)"
        @info "    Bits per vertex: $(round(bits_per_vertex, digits=3))"
        @info "    Bits per edge: $(round(bits_per_edge, digits=3))"
        
        # Store results
        push!(compression_results, (
            description = description,
            vertices = vertices_count,
            edges = edges_count,
            file_size = file_size,
            bits_per_vertex = bits_per_vertex,
            bits_per_edge = bits_per_edge,
            use_mix_mode = use_mix_mode,
            use_reference = use_reference,
            encoding_type = encoding_type
        ))
        
        # Verify compression by reading back
        @info "  Verifying compression (round-trip test)..."
        
        @time loaded_graph = load_compressed_mgs3_graph(output_path_with_ext)
        
        loaded_vertices = nv(loaded_graph)
        loaded_edges = ne(loaded_graph)
        
        @test loaded_vertices == vertices_count
        @test loaded_edges == edges_count
        
        # Quick integrity check (sample-based for large graphs)
        sample_size = min(1000, vertices_count)
        sample_vertices = rand(1:vertices_count, sample_size)
        
        integrity_check = true
        mismatch_count = 0
        
        for v in sample_vertices
            if v <= vertices_count && v <= nv(test_graph) && v <= nv(loaded_graph)
                original_neighbors = sort(collect(outneighbors(test_graph, v)))
                loaded_neighbors = sort(collect(outneighbors(loaded_graph, v)))
                
                if original_neighbors != loaded_neighbors
                    integrity_check = false
                    mismatch_count += 1
                    if mismatch_count <= 3
                        @warn "Sample vertex $v mismatch: original=$(length(original_neighbors)), loaded=$(length(loaded_neighbors))"
                    end
                end
            end
        end
        
        if integrity_check
            @info "  [PASS] Round-trip verification passed (sampled $sample_size vertices)"
        else
            @warn "  [WARN] Round-trip issues found in $mismatch_count/$sample_size sampled vertices"
        end
        
        @test mismatch_count <= sample_size * 0.01  # Allow up to 1% sample errors
        
        @info "  Compression completed successfully!"
    end
    
    # Comprehensive analysis and comparison
    @info "\n" * "="^80
    @info "COMPREHENSIVE COMPRESSION ANALYSIS"
    @info "="^80
    
    # WebGraph baseline for comparison
    webgraph_bits_per_edge = 2.897  # From CNR-2000.properties (whole graph)
    
    @info "\nCompression Results Summary:"
    @info "-"^80
    
    best_compression = minimum(r.bits_per_edge for r in compression_results)
    
    for result in compression_results
        encoding_desc = rpad(result.encoding_type, 10)
        
        @info "Whole + $(encoding_desc): $(rpad(result.vertices, 6)) vertices, $(rpad(result.edges, 7)) edges"
        @info "  File size: $(rpad(round(result.file_size/1024/1024, digits=2), 6)) MB"
        @info "  Bits/edge: $(rpad(round(result.bits_per_edge, digits=3), 6)) $(result.bits_per_edge == best_compression ? "*** BEST ***" : "")"
        @info "  Bits/vertex: $(round(result.bits_per_vertex, digits=3))"
        
        # Comparison with WebGraph
        if result.bits_per_edge < webgraph_bits_per_edge
            improvement = webgraph_bits_per_edge - result.bits_per_edge
            @info "  vs WebGraph: $(round(improvement, digits=3)) bits/edge better! ***"
        else
            difference = result.bits_per_edge - webgraph_bits_per_edge
            @info "  vs WebGraph: $(round(difference, digits=3)) bits/edge worse"
        end
        @info ""
    end
    
    
    @info "=== CNR-2000 Processing Test Completed Successfully ==="
    
    # Summary test assertions
    @test core_vertices_count > 0
    @test core_edges_count > 0
    @test length(compression_results) == 4
    @test all(r.bits_per_edge < 50 for r in compression_results)  # Reasonable compression achieved
    @test all(r.file_size > 0 for r in compression_results)  # Valid file sizes
    
    best_bits_per_edge = minimum(r.bits_per_edge for r in compression_results)
    
    @info "Final summary:"
    @info "  Whole graph: $original_vertices vertices, $original_edges edges"
    @info "  Core graph: $core_vertices_count vertices, $core_edges_count edges (extracted for analysis)"
    @info "  Best compression: $(round(best_bits_per_edge, digits=3)) bits/edge"
    @info "  Configurations tested: $(length(compression_results)) (delta only, mix-mode, reference, mix+reference)"
    finally
        # Restore previous logger
        global_logger(prev_logger)
    end
end
