#!/usr/bin/env julia

#
# Verify the current state of the Amazon graph compression algorithm
# Testing full dataset with proper vertex remapping
#

include("test/run_tests_main.jl")

@testset "Full Amazon Graph Verification Test" begin
    @info "=== Full Amazon Graph Verification Test ==="
    
    # Load Amazon dataset
    @info "Loading full Amazon dataset..."
    amz_g = load_dataset(AMZ_DATASET_IN; separator='\t')
    amz_core, oni, noi = get_core(amz_g)
    full_neighbor_lists = get_neighbor_lists(amz_core)
    
    T = UInt24
    encoding = :zeta
    
    total_vertices = length(full_neighbor_lists)
    total_edges = sum(length(neighbors) for neighbors in values(full_neighbor_lists))
    
    @info "Full Amazon core graph:"
    @info "  Vertices: $total_vertices"
    @info "  Total edges (directed): $total_edges"
    @info "  Average degree: $(round(total_edges / total_vertices, digits=2))"
    
    @info "\n" * "="^80
    @info "CREATING PROPER VERTEX REMAPPING"
    @info "="^80
    
    # Step 1: Get sorted list of all vertex IDs for consistent ordering
    vertex_ids = sort(collect(keys(full_neighbor_lists)))
    @info "Vertex ID range: $(vertex_ids[1]) to $(vertex_ids[end])"
    
    # Step 2: Create old_id -> new_id mapping (1-based sequential)
    old_to_new = Dict{T,T}()
    new_to_old = Dict{T,T}()
    
    for (new_id, old_id) in enumerate(vertex_ids)
        old_to_new[old_id] = T(new_id)
        new_to_old[T(new_id)] = old_id
    end
    
    @info "Created mapping for $total_vertices vertices"
    @info "  Example: vertex $(vertex_ids[1]) -> 1, vertex $(vertex_ids[end]) -> $total_vertices"
    
    # Step 3: Create remapped neighbor lists
    @info "Remapping neighbor lists..."
    remap_start = time()
    
    remapped_neighbor_lists = Dict{T, Vector{T}}()
    vertices_processed = 0
    progress_interval = max(1, total_vertices ÷ 100)  # Progress every 1%
    
    for (new_id, old_id) in enumerate(vertex_ids)
        old_neighbors = full_neighbor_lists[old_id]
        
        # Remap neighbor IDs (all neighbors should be in our mapping since we're using the full graph)
        new_neighbors = T[]
        for old_neighbor in old_neighbors
            if haskey(old_to_new, old_neighbor)
                push!(new_neighbors, old_to_new[old_neighbor])
            else
                @warn "Neighbor $old_neighbor not found in mapping - this shouldn't happen for full graph"
            end
        end
        
        sort!(new_neighbors)  # Keep neighbors sorted
        remapped_neighbor_lists[T(new_id)] = new_neighbors
        
        vertices_processed += 1
        if vertices_processed % progress_interval == 0
            elapsed = time() - remap_start
            @info "  Progress: $vertices_processed/$total_vertices vertices ($(round(100*vertices_processed/total_vertices, digits=1))%), $(round(elapsed, digits=1))s elapsed"
        end
    end
    
    remap_time = time() - remap_start
    @info "Remapping completed in $(round(remap_time, digits=2))s"
    
    # Verify remapping
    @info "Verifying remapped data:"
    @info "  Remapped vertices: $(length(remapped_neighbor_lists))"
    @info "  Vertex IDs: 1 to $(maximum(keys(remapped_neighbor_lists)))"
    @info "  Max neighbor ID: $(maximum(maximum(neighbors) for neighbors in values(remapped_neighbor_lists) if !isempty(neighbors)))"
    
    @info "\n" * "="^80
    @info "COMPRESSION TEST"
    @info "="^80
    
    # Test compression with the properly remapped full graph
    @info "Starting compression of full Amazon graph..."
    compress_start = time()
    
    try
        # Compress
        io = IOBuffer()
        writer = BitWriter(io)
        Adjacently.Compression.write_compressed_graph_data(writer, remapped_neighbor_lists, encoding, :children, true, true)
        flush_bitwriter(writer; flush_last_bits=true)
        
        compress_time = time() - compress_start
        compressed_size = position(io)
        
        @info "Compression successful!"
        @info "  Compressed size: $compressed_size bytes ($(round(compressed_size/1024/1024, digits=2)) MB)"
        @info "  Compression time: $(round(compress_time, digits=2))s"
        @info "  Compression rate: $(round(compressed_size * 8 / total_edges, digits=2)) bits per edge"
        @info "  Throughput: $(round(total_vertices / compress_time, digits=0)) vertices/second"
        
        @info "\n" * "="^80
        @info "DECOMPRESSION TEST"
        @info "="^80
        
        # Test decompression
        @info "Starting decompression..."
        decompress_start = time()
        
        seekstart(io)
        reader = BitReader(io)
        decoded_neighbor_lists = Adjacently.Compression.read_compressed_graph_data(reader, T(total_vertices), encoding, :children, T)
        
        decompress_time = time() - decompress_start
        
        @info "Decompression successful!"
        @info "  Decompression time: $(round(decompress_time, digits=2))s"
        @info "  Throughput: $(round(total_vertices / decompress_time, digits=0)) vertices/second"
        
        @info "\n" * "="^80
        @info "VALIDATION"
        @info "="^80
        
        # Validate the results
        @info "Validating decompressed data..."
        validation_start = time()
        
        all_correct = true
        mismatch_count = 0
        first_error_vertex = nothing
        validation_interval = max(1, total_vertices ÷ 1000)  # Progress every 0.1%
        
        for i in 1:total_vertices
            vertex_id = T(i)
            original = remapped_neighbor_lists[vertex_id]
            decoded = get(decoded_neighbor_lists, vertex_id, T[])
            
            if Set(original) != Set(decoded)
                all_correct = false
                mismatch_count += 1
                
                if first_error_vertex === nothing
                    first_error_vertex = i
                    @error "First error at vertex $i:"
                    if length(original) <= 10
                        @error "  Original ($(length(original))): $original"
                    else
                        @error "  Original ($(length(original))): [$(original[1:5])..., $(original[end-2:end])]"
                    end
                    
                    if length(decoded) <= 10
                        @error "  Decoded  ($(length(decoded))): $decoded"
                    else
                        @error "  Decoded  ($(length(decoded))): [$(decoded[1:5])..., $(decoded[end-2:end])]"
                    end
                    
                    if length(decoded) > length(original) * 10
                        @error "  Massive over-decompression detected!"
                    end
                end
                
                # Stop after finding first few errors to avoid log spam
                if mismatch_count >= 10
                    @warn "Stopping detailed error reporting after 10 mismatches to avoid log spam..."
                    break
                end
            end
            
            # Progress logging
            if i % validation_interval == 0
                elapsed = time() - validation_start
                @info "  Validation progress: $i/$total_vertices vertices ($(round(100*i/total_vertices, digits=1))%), $mismatch_count errors so far, $(round(elapsed, digits=1))s elapsed"
            end
        end
        
        # If we stopped early due to too many errors, count the remaining quickly
        if mismatch_count >= 10
            @info "Counting remaining mismatches quickly..."
            for i in (first_error_vertex+10):total_vertices
                vertex_id = T(i)
                original = remapped_neighbor_lists[vertex_id]
                decoded = get(decoded_neighbor_lists, vertex_id, T[])
                
                if Set(original) != Set(decoded)
                    mismatch_count += 1
                end
            end
        end
        
        validation_time = time() - validation_start
        
        if all_correct
            @info "ALL VERTICES CORRECT!"
            @info "The reference encoding works perfectly on the full Amazon graph."
        else
            @error "Found $mismatch_count mismatches out of $total_vertices vertices"
            @error "First error at vertex: $first_error_vertex"
            @error "Error rate: $(round(100*mismatch_count/total_vertices, digits=3))%"
            @error "Success rate: $(round(100*(total_vertices - mismatch_count)/total_vertices, digits=3))%"
        end
        
        @info "Validation completed in $(round(validation_time, digits=2))s"
        
        @info "\n" * "="^80
        @info "SUMMARY"
        @info "="^80
        
        total_time = compress_time + decompress_time + validation_time
        
        @info "Full Amazon graph test summary:"
        @info "  Vertices: $total_vertices"
        @info "  Edges: $total_edges" 
        @info "  Compressed size: $(round(compressed_size/1024/1024, digits=2)) MB"
        @info "  Compression rate: $(round(compressed_size * 8 / total_edges, digits=2)) bits/edge"
        @info "  Total time: $(round(total_time, digits=2))s"
        if all_correct
            @info "  Result: SUCCESS"
        else
            @info "  Result: FAILED"
        end
        
        if !all_correct && first_error_vertex !== nothing
            @info "  First error at vertex: $first_error_vertex ($(round(100*first_error_vertex/total_vertices, digits=2))% through)"
        end
        
        @test all_correct  # Test should pass only if everything is correct
        
    catch e
        @error "Exception during compression/decompression: $e"
        @error "Stack trace: $(stacktrace())"
        @test false
    end
    
end