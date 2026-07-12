#!/usr/bin/env julia

#
# Adjacently: Julia Complex Directed Networks Library
# Copyright (C) 2016-2026 Anonymous (double-blind review)
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
using Printf

@testset "Graph Compression Rate Analysis" begin
    @info "=== Amazon Graph Compression Rate Analysis ==="
    
    # Load Amazon dataset
    @info "Loading Amazon dataset..."
    amz_g = load_dataset(AMZ_DATASET_IN; separator='\t')
    @test 403394 == convert(Int,nv(amz_g))
    @test 3387388 == ne(amz_g)
    
    # Get core
    @info "Getting core..."
    amz_core,oni,noi = get_core(amz_g)
    @test 395234 == convert(Int,nv(amz_core))
    @test 3301092 == ne(amz_core)
    
    vs = vertices(amz_core)
    gs = convert(UInt64, length(vs))
    n_bits_v = convert(UInt8, ceil(log(2, gs)))
    T = infer_uint_custom_type(n_bits_v) 
    es = ne(amz_core)
    
    @info "Graph statistics:"
    @info "  Vertices: $gs"
    @info "  Edges: $es"
    @info "  Average degree: $(round(2*es/gs, digits=2))"
    @info "  Integer type: $T ($(8*sizeof(T)) bits)"
    
    # Test configurations
    encodings = [:elias_delta, :fibonacci, :zeta]
    relabeling_schemes = [:none, :lexicographic]
    
    # Compression configurations to test
    compression_configs = [
        (name="Delta Only", use_mix=false, reference=false),
        (name="Mix Mode", use_mix=true, reference=false),  
        (name="Mix Mode + Reference", use_mix=true, reference=true)
    ]
    
    # Results storage
    results = []
    
    @info "\nStarting comprehensive compression analysis...\n"
    
    for relabeling_scheme in relabeling_schemes
        @info "="^80
        @info "RELABELING SCHEME: $relabeling_scheme"
        @info "="^80
        
        # Apply relabeling
        current_graph = amz_core
        if relabeling_scheme != :none
            @info "  Applying $relabeling_scheme relabeling..."
            relabel_start = time()
            old_to_new_vertex_ids = relabel_vertices(current_graph, relabeling_scheme)
            current_graph = relabel_graph(current_graph, old_to_new_vertex_ids)
            relabel_time = time() - relabel_start
            @info "  Relabeling completed in $(round(relabel_time, digits=3))s"
        end
        
        # Get neighbor lists
        @info "  Creating neighbor lists..."
        neighbor_start = time()
        neighbor_lists = get_neighbor_lists(current_graph)
        neighbor_time = time() - neighbor_start
        @test length(neighbor_lists) == gs
        @info "  Neighbor lists created in $(round(neighbor_time, digits=3))s"
        
        for config in compression_configs
            @info "\n" * "-"^60
            @info "COMPRESSION CONFIG: $(config.name)"
            @info "-"^60
            
            for encoding in encodings
                @info "\n  Testing encoding: $encoding"
                
                compression_start = time()
                
                # Compress
                io = IOBuffer()
                writer = BitWriter(io)
                
                @info "    Compressing..."
                compress_start = time()
                write_compressed_graph_data(writer, neighbor_lists, encoding, :children, 
                                           config.reference, config.use_mix)
                flush_bitwriter(writer; flush_last_bits=true)
                compress_time = time() - compress_start
                
                compressed_bits = position(io) * 8
                total_time = time() - compression_start
                
                # Calculate compression metrics
                bits_per_vertex = compressed_bits / gs
                bits_per_edge = compressed_bits / es
                
                @info "    Results:"
                @info "      Compressed size: $compressed_bits bits ($(round(compressed_bits/1024/1024, digits=2)) MB)"
                @info "      Compression time: $(round(compress_time, digits=3))s"
                @info "      Bits per vertex: $(round(bits_per_vertex, digits=3))"
                @info "      Bits per edge: $(round(bits_per_edge, digits=3))"
                
                # Store results
                push!(results, Dict(
                    :relabeling => string(relabeling_scheme),
                    :config => config.name,
                    :encoding => string(encoding),
                    :use_mix => config.use_mix,
                    :reference => config.reference,
                    :compressed_bits => compressed_bits,
                    :compress_time => compress_time,
                    :bits_per_vertex => bits_per_vertex,
                    :bits_per_edge => bits_per_edge,
                    :vertices => gs,
                    :edges => es
                ))
                
                # Basic validation - try to decompress first few bytes to ensure format is correct
                @info "    Validating compression format..."
                seekstart(io)
                reader = BitReader(io)
                try
                    # Just read the reference flag to verify format
                    ref_flag = read_bit(reader)
                    @info "      Format validation: ✓ (reference_flag=$ref_flag)"
                catch e
                    @warn "      Format validation failed: $e"
                end
            end
        end
    end
    
    @info "\n" * "="^80
    @info "COMPRESSION ANALYSIS SUMMARY"
    @info "="^80
    
    # Create summary tables
    @info "\nDETAILED RESULTS TABLE:"
    @info "Relabeling | Config        | Encoding     | Bits/Vertex | Bits/Edge | Time(s) | Size(MB)"
    @info "-"^90
    
    # Sort results for better presentation
    sorted_results = sort(results, by=r -> (r[:relabeling], r[:config], r[:encoding]))
    
    for result in sorted_results
        relabel_str = rpad(result[:relabeling], 10)
        config_str = rpad(result[:config], 13) 
        encoding_str = rpad(result[:encoding], 12)
        bits_vertex = @sprintf("%.3f", result[:bits_per_vertex])
        bits_edge = @sprintf("%.3f", result[:bits_per_edge])
        time_str = @sprintf("%.3f", result[:compress_time])
        size_mb = @sprintf("%.2f", result[:compressed_bits]/1024/1024)
        
        @info "$relabel_str | $config_str | $encoding_str | $bits_vertex    | $bits_edge   | $time_str   | $size_mb"
    end
    
    # Analysis by encoding scheme
    @info "\n" * "="^60
    @info "ANALYSIS BY ENCODING SCHEME"
    @info "="^60
    
    for encoding in encodings
        @info "\n$encoding:"
        encoding_results = filter(r -> r[:encoding] == string(encoding), results)
        
        best_result = encoding_results[argmin([r[:bits_per_vertex] for r in encoding_results])]
        worst_result = encoding_results[argmax([r[:bits_per_vertex] for r in encoding_results])]
        
        @info "  Best:  $(best_result[:bits_per_vertex]) bits/vertex ($(best_result[:relabeling]), $(best_result[:config]))"
        @info "  Worst: $(worst_result[:bits_per_vertex]) bits/vertex ($(worst_result[:relabeling]), $(worst_result[:config]))"
        @info "  Range: $(round(worst_result[:bits_per_vertex] - best_result[:bits_per_vertex], digits=3)) bits/vertex"
    end
    
    # Analysis by configuration
    @info "\n" * "="^60  
    @info "ANALYSIS BY COMPRESSION CONFIGURATION"
    @info "="^60
    
    for config in compression_configs
        @info "\n$(config.name):"
        config_results = filter(r -> r[:config] == config.name, results)
        
        avg_bits_vertex = mean(r[:bits_per_vertex] for r in config_results)
        avg_bits_edge = mean(r[:bits_per_edge] for r in config_results)
        avg_time = mean(r[:compress_time] for r in config_results)
        
        @info "  Average bits/vertex: $(round(avg_bits_vertex, digits=3))"
        @info "  Average bits/edge: $(round(avg_bits_edge, digits=3))" 
        @info "  Average time: $(round(avg_time, digits=3))s"
        
        best = config_results[argmin([r[:bits_per_vertex] for r in config_results])]
        @info "  Best combination: $(best[:encoding]) + $(best[:relabeling]) → $(round(best[:bits_per_vertex], digits=3)) bits/vertex"
    end
    
    # Analysis by relabeling
    @info "\n" * "="^60
    @info "ANALYSIS BY RELABELING SCHEME" 
    @info "="^60
    
    for relabeling in relabeling_schemes
        @info "\n$relabeling:"
        relabel_results = filter(r -> r[:relabeling] == string(relabeling), results)
        
        avg_bits_vertex = mean(r[:bits_per_vertex] for r in relabel_results)
        improvement_vs_none = relabeling == :none ? 0.0 : 
            avg_bits_vertex - mean(r[:bits_per_vertex] for r in results if r[:relabeling] == "none")
        
        @info "  Average bits/vertex: $(round(avg_bits_vertex, digits=3))"
        if relabeling != :none
            @info "  Improvement vs none: $(round(-improvement_vs_none, digits=3)) bits/vertex ($(round(-100*improvement_vs_none/avg_bits_vertex, digits=1))%)"
        end
    end
    
    # Overall best and worst results
    @info "\n" * "="^60
    @info "OVERALL BEST AND WORST RESULTS"
    @info "="^60
    
    overall_best = results[argmin([r[:bits_per_vertex] for r in results])]
    overall_worst = results[argmax([r[:bits_per_vertex] for r in results])]
    
    @info "\nBest Overall:"
    @info "  Configuration: $(overall_best[:relabeling]) + $(overall_best[:config]) + $(overall_best[:encoding])"
    @info "  Compression: $(round(overall_best[:bits_per_vertex], digits=3)) bits/vertex, $(round(overall_best[:bits_per_edge], digits=3)) bits/edge"
    @info "  Size: $(round(overall_best[:compressed_bits]/1024/1024, digits=2)) MB"
    @info "  Time: $(round(overall_best[:compress_time], digits=3))s"
    
    @info "\nWorst Overall:"
    @info "  Configuration: $(overall_worst[:relabeling]) + $(overall_worst[:config]) + $(overall_worst[:encoding])"  
    @info "  Compression: $(round(overall_worst[:bits_per_vertex], digits=3)) bits/vertex, $(round(overall_worst[:bits_per_edge], digits=3)) bits/edge"
    @info "  Size: $(round(overall_worst[:compressed_bits]/1024/1024, digits=2)) MB"
    @info "  Time: $(round(overall_worst[:compress_time], digits=3))s"
    
    compression_ratio = overall_worst[:bits_per_vertex] / overall_best[:bits_per_vertex]
    @info "\nOverall compression range: $(round(compression_ratio, digits=2))x difference between best and worst"
    
    # Reference encoding impact analysis  
    @info "\n" * "="^60
    @info "REFERENCE ENCODING IMPACT ANALYSIS"
    @info "="^60
    
    ref_results = filter(r -> r[:reference] == true, results)
    no_ref_results = filter(r -> r[:reference] == false && r[:config] == "Mix Encoding", results)
    
    if !isempty(ref_results) && !isempty(no_ref_results)
        @info "\nReference encoding vs Mix encoding (same conditions):"
        
        for encoding in encodings
            for relabeling in relabeling_schemes
                ref_result = filter(r -> r[:encoding] == string(encoding) && 
                                        r[:relabeling] == string(relabeling) && 
                                        r[:reference] == true, results)
                mix_result = filter(r -> r[:encoding] == string(encoding) && 
                                        r[:relabeling] == string(relabeling) && 
                                        r[:config] == "Mix Encoding", results)
                
                if !isempty(ref_result) && !isempty(mix_result)
                    ref_bits = ref_result[1][:bits_per_vertex]
                    mix_bits = mix_result[1][:bits_per_vertex]
                    improvement = mix_bits - ref_bits
                    improvement_pct = 100 * improvement / mix_bits
                    
                    @info "  $encoding + $relabeling: $(round(improvement, digits=3)) bits/vertex improvement ($(round(improvement_pct, digits=1))%)"
                end
            end
        end
        
        avg_ref = mean(r[:bits_per_vertex] for r in ref_results)
        avg_mix = mean(r[:bits_per_vertex] for r in no_ref_results)
        overall_improvement = avg_mix - avg_ref
        overall_improvement_pct = 100 * overall_improvement / avg_mix
        
        @info "\nOverall reference encoding impact:"
        @info "  Average improvement: $(round(overall_improvement, digits=3)) bits/vertex ($(round(overall_improvement_pct, digits=1))%)"
    end
    
    @info "\n" * "="^80
    @info "COMPRESSION ANALYSIS COMPLETED"
    @info "="^80
    @info "Total configurations tested: $(length(results))"
    @info "Analysis completed successfully!"
    
    # Basic sanity checks
    @test length(results) == length(encodings) * length(relabeling_schemes) * length(compression_configs)
    @test all(r -> r[:compressed_bits] > 0, results)
    @test all(r -> r[:bits_per_vertex] > 0, results)
    @test all(r -> r[:bits_per_edge] > 0, results)
    
end