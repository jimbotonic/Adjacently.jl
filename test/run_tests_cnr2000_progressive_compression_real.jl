#!/usr/bin/env julia

#
# Progressive CNR-2000 Compression Test using Real Dataset
# Tests compression with increasing options: delta only, delta + mix, delta + mix + reference, delta + mix + recursive reference
# Uses the actual CNR-2000 dataset from CSV and the write_compressed_mgs3_graph function
#

using Pkg
Pkg.activate(joinpath(@__DIR__, ".."))

using Test
using LightGraphs: nv, ne, vertices, outneighbors
using Adjacently
using Adjacently.IO: load_adjacency_list_from_csv
using Adjacently.MGS: write_compressed_mgs3_graph, load_compressed_mgs3_graph
using Adjacently.Relabeling: relabel_vertices_rcm, relabel_graph
using Logging

@testset "CNR-2000 Real Dataset Progressive Compression Test" begin

    # Enable info-level logging - debug creates too much output
    prev_logger = current_logger()
    global_logger(ConsoleLogger(stderr, Logging.Info))

    try
        @info "=== CNR-2000 Real Dataset Progressive Compression Test ==="

        # Define paths
        cnr_csv_path = "datasets/webgraph/cnr-2000/cnr-2000.csv"
        output_dir = "test_data"
        base_filename = "cnr2000_real_progressive"

        # Check if CSV file exists
        if !isfile(cnr_csv_path)
            @warn "CNR-2000 CSV file not found at: $cnr_csv_path"
            @warn "Please ensure the CNR-2000 dataset is available"
            @test_skip "CNR-2000 CSV file not available"
            return
        end

        @info "Step 1: Loading CNR-2000 graph from CSV..."
        @info "  Loading from: $cnr_csv_path"
        @time graph = load_adjacency_list_from_csv(cnr_csv_path, ',', true)

        vertices_count = nv(graph)
        edges_count = ne(graph)

        @info "  Successfully loaded CNR-2000 graph:"
        @info "    Vertices: $vertices_count"
        @info "    Edges: $edges_count"
        @info "    Average degree: $(round(2*edges_count/vertices_count, digits=2))"

        @info "Step 2: Applying RCM relabeling for better compression..."
        @time begin
            # Apply RCM relabeling using out-degree criterion (best for directed graphs)
            rcm_mapping = relabel_vertices_rcm(graph, :out)
            @info "    RCM mapping computed: $(length(rcm_mapping)) vertices"

            # Apply the relabeling to get the reordered graph
            graph = relabel_graph(graph, rcm_mapping)
            @info "    Graph relabeled with RCM ordering"
        end

        # Verify graph is still the same size after relabeling
        rcm_vertices_count = nv(graph)
        rcm_edges_count = ne(graph)
        @info "    RCM graph verification: $rcm_vertices_count vertices, $rcm_edges_count edges"
        @test rcm_vertices_count == vertices_count
        @test rcm_edges_count == edges_count

        # Create output directory
        mkpath(output_dir)

        # Progressive compression configurations
        # Each test progressively adds more compression features
        compression_configs = [
            (
                name = "Delta Only",
                description = "Delta encoding only",
                filename_suffix = "delta_only",
                compression = :complex,
                use_mix_mode = false,
                reference_enabled = false,
                recursive_reference = false
            ),
            (
                name = "Delta + Mix",
                description = "Delta + mix-mode encoding (run-length + interval)",
                filename_suffix = "delta_mix",
                compression = :complex,
                use_mix_mode = true,
                reference_enabled = false,
                recursive_reference = false
            ),
            (
                name = "Delta + Mix + Reference",
                description = "Delta + mix-mode + reference encoding",
                filename_suffix = "delta_mix_ref",
                compression = :complex,
                use_mix_mode = true,
                reference_enabled = true,
                recursive_reference = false
            ),
            (
                name = "Delta + Mix + Recursive Reference",
                description = "Delta + mix-mode + recursive reference encoding (hybrid+)",
                filename_suffix = "delta_mix_recursive_ref",
                compression = :complex,
                use_mix_mode = true,
                reference_enabled = true,
                recursive_reference = true
            )
        ]

        compression_results = []

        # Test each configuration progressively
        for (i, config) in enumerate(compression_configs)
            @info "\n" * "="^80
            @info "Configuration $(i): Testing $(config.name)"
            @info "Description: $(config.description)"
            @info "="^80

            output_path = joinpath(output_dir, base_filename * "_" * config.filename_suffix)

            @info "  Compression Parameters:"
            @info "    Coding scheme: :children"
            @info "    Compression type: $(config.compression)"
            @info "    Mix mode: $(config.use_mix_mode)"
            @info "    Reference enabled: $(config.reference_enabled)"
            @info "    Recursive reference: $(config.recursive_reference)"

            @info "  Starting compression of CNR-2000 graph..."

            compression_start_time = time()

            @time begin
                write_compressed_mgs3_graph(
                    graph,
                    output_path,
                    :children,                      # coding_scheme
                    config.compression,             # compression (:complex for advanced encoding)
                    :fibonacci,                     # integer_encoding
                    config.use_mix_mode,            # use_mix_mode
                    config.reference_enabled,       # reference_enabled
                    config.recursive_reference      # recursive_reference
                )
            end

            compression_time = time() - compression_start_time
            output_file = output_path * ".mgz"

            # Verify output file was created
            if !isfile(output_file)
                @error "  Compression failed - output file not created: $output_file"
                continue
            end

            # Get compression statistics
            file_size = filesize(output_file)
            original_size_estimate = edges_count * 4  # Rough estimate: 4 bytes per edge
            compression_ratio = original_size_estimate / file_size
            bits_per_vertex = (file_size * 8) / vertices_count
            bits_per_edge = (file_size * 8) / edges_count
            space_saving_pct = (1.0 - 1.0/compression_ratio) * 100

            @info "  Compression Results:"
            @info "    Output file: $(basename(output_file))"
            @info "    File size: $(file_size) bytes ($(round(file_size/1024/1024, digits=3)) MB)"
            @info "    Compression ratio: $(round(compression_ratio, digits=3))"
            @info "    Bits per vertex: $(round(bits_per_vertex, digits=3))"
            @info "    Bits per edge: $(round(bits_per_edge, digits=3))"
            @info "    Space saving: $(round(space_saving_pct, digits=1))%"
            @info "    Compression time: $(round(compression_time, digits=2))s"

            # Store results for comparison
            push!(compression_results, (
                name = config.name,
                description = config.description,
                file_size = file_size,
                compression_ratio = compression_ratio,
                bits_per_vertex = bits_per_vertex,
                bits_per_edge = bits_per_edge,
                space_saving_pct = space_saving_pct,
                compression_time = compression_time,
                use_mix_mode = config.use_mix_mode,
                reference_enabled = config.reference_enabled,
                recursive_reference = config.recursive_reference,
                output_file = output_file
            ))

            # Verify compression integrity by loading back
            @info "  Performing round-trip verification..."

            decompression_start_time = time()

            @time loaded_graph = load_compressed_mgs3_graph(output_file)

            decompression_time = time() - decompression_start_time

            loaded_vertices = nv(loaded_graph)
            loaded_edges = ne(loaded_graph)

            # Basic structure verification
            @test loaded_vertices == vertices_count
            @test loaded_edges == edges_count

            @info "    Decompressed graph: $loaded_vertices vertices, $loaded_edges edges"
            @info "    Decompression time: $(round(decompression_time, digits=2))s"

            # Quick integrity check - just verify basic structure
            if loaded_vertices == vertices_count && loaded_edges == edges_count
                @info "  [✓ PASS] Basic round-trip verification successful (same vertex/edge counts)"
                integrity_pass = true
            else
                @warn "  [⚠ WARN] Vertex/edge count mismatch: loaded vs original"
                integrity_pass = false
            end

            # Skip detailed vertex-by-vertex check to avoid memory issues
            @test integrity_pass

            @info "  Configuration $(config.name) completed successfully!"
            @info "    Total time: $(round(compression_time + decompression_time, digits=2))s"
        end

        # Comprehensive comparison and analysis
        @info "\n" * "="^80
        @info "PROGRESSIVE COMPRESSION ANALYSIS - CNR-2000 REAL DATASET"
        @info "="^80

        @info "\nDataset Summary:"
        @info "  Source: CNR-2000 Web Graph Dataset"
        @info "  Vertices: $(vertices_count)"
        @info "  Edges: $(edges_count)"
        @info "  Average degree: $(round(2*edges_count/vertices_count, digits=2))"
        @info "  Configurations tested: $(length(compression_results))"

        @info "\nCompression Results Comparison:"
        @info "-"^80

        # Find best compression
        best_compression = minimum(r.bits_per_edge for r in compression_results)
        best_config = compression_results[findfirst(r -> r.bits_per_edge == best_compression, compression_results)]

        for (i, result) in enumerate(compression_results)
            @info "$(i). $(rpad(result.name, 40))"
            @info "    Description: $(result.description)"
            @info "    File size: $(rpad(round(result.file_size/1024/1024, digits=3), 8)) MB"
            @info "    Compression ratio: $(rpad(round(result.compression_ratio, digits=2), 6))x"
            @info "    Bits per edge: $(rpad(round(result.bits_per_edge, digits=3), 7)) $(result.bits_per_edge == best_compression ? "*** BEST ***" : "")"
            @info "    Bits per vertex: $(round(result.bits_per_vertex, digits=3))"
            @info "    Space saving: $(round(result.space_saving_pct, digits=1))%"
            @info "    Total time: $(round(result.compression_time, digits=2))s"

            # Show improvement over delta-only baseline
            if i > 1
                baseline = compression_results[1]
                improvement = baseline.bits_per_edge - result.bits_per_edge
                improvement_pct = (improvement / baseline.bits_per_edge) * 100
                if improvement > 0
                    @info "    vs Delta-only: $(round(improvement, digits=3)) bits/edge better ($(round(improvement_pct, digits=1))% improvement)"
                else
                    @info "    vs Delta-only: $(round(abs(improvement), digits=3)) bits/edge worse ($(round(abs(improvement_pct), digits=1))% regression)"
                end
            end
            @info "    Output: $(basename(result.output_file))"
            @info ""
        end

        # Progressive improvement analysis
        @info "Progressive Enhancement Analysis:"
        @info "-"^40
        for i in 2:length(compression_results)
            prev_result = compression_results[i-1]
            curr_result = compression_results[i]
            improvement = prev_result.bits_per_edge - curr_result.bits_per_edge
            improvement_pct = (improvement / prev_result.bits_per_edge) * 100

            enhancement = if curr_result.use_mix_mode && !prev_result.use_mix_mode
                "Added Mix Mode"
            elseif curr_result.reference_enabled && !prev_result.reference_enabled
                "Added Reference Encoding"
            elseif curr_result.recursive_reference && !prev_result.recursive_reference
                "Added Recursive Reference"
            else
                "Unknown Enhancement"
            end

            if improvement > 0
                @info "$enhancement: $(round(improvement, digits=3)) bits/edge improvement ($(round(improvement_pct, digits=1))%)"
            else
                @info "$enhancement: $(round(abs(improvement), digits=3)) bits/edge worse ($(round(abs(improvement_pct), digits=1))%)"
            end
        end

        # Overall summary
        total_improvement = compression_results[1].bits_per_edge - compression_results[end].bits_per_edge
        total_improvement_pct = (total_improvement / compression_results[1].bits_per_edge) * 100

        @info "\nOverall Summary:"
        @info "  Dataset processed with RCM relabeling: ✅ Applied"
        @info "  Best configuration: $(best_config.name)"
        @info "  Best compression: $(round(best_compression, digits=3)) bits/edge"
        @info "  Total improvement: $(round(total_improvement, digits=3)) bits/edge ($(round(total_improvement_pct, digits=1))%)"
        @info "  Compression range: $(round(minimum(r.compression_ratio for r in compression_results), digits=2))x to $(round(maximum(r.compression_ratio for r in compression_results), digits=2))x"

        # Count recursive reference configurations
        recursive_configs = length([r for r in compression_results if r.recursive_reference])
        @info "  Recursive reference configurations tested: $recursive_configs/$(length(compression_results))"

        # WebGraph comparison (if we know WebGraph's performance on CNR-2000)
        webgraph_bpe = 2.897  # bits per edge from WebGraph BV compression
        @info "\nComparison with WebGraph BV:"
        for result in compression_results
            if result.bits_per_edge < webgraph_bpe
                improvement = webgraph_bpe - result.bits_per_edge
                @info "  $(result.name): $(round(improvement, digits=3)) bits/edge better than WebGraph!"
            else
                difference = result.bits_per_edge - webgraph_bpe
                @info "  $(result.name): $(round(difference, digits=3)) bits/edge worse than WebGraph"
            end
        end

        # Test assertions
        @test length(compression_results) == 4
        @test all(r.bits_per_edge > 0 for r in compression_results)
        @test all(r.file_size > 0 for r in compression_results)
        @test all(r.compression_ratio > 0 for r in compression_results)

        # Expect reasonable compression performance
        @test all(r.bits_per_edge < 50 for r in compression_results)  # Should be reasonably compressed
        @test best_compression < 10  # Best should be quite good

        @info "\n=== CNR-2000 Real Dataset Progressive Compression Test Completed Successfully ==="

    finally
        # Restore original logger
        global_logger(prev_logger)
    end
end