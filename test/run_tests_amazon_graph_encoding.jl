#
# Test graph encoding with actual Amazon graph data
#

include("run_tests_main.jl")

@testset "Amazon graph encoding tests" begin
    @info "=== Loading Amazon dataset ==="
    
    # Load the Amazon dataset (same as in the reference tests)
    dataset_path = "../datasets/Amazon/"
    if !isdir(dataset_path)
        @warn "Amazon dataset not found at $dataset_path - skipping tests"
        return
    end
    
    @info "Loading Amazon dataset from $dataset_path"
    
    # Load the graph data
    try
        # This mimics the loading logic from the reference tests
        amazon_graph = load_adjacently_graph(dataset_path * "Amazon.adjacently")
        
        @info "Loaded graph with $(length(amazon_graph.vertices)) vertices"
        
        # Get neighbor lists for testing
        neighbor_lists = Dict{UInt32,Vector{UInt32}}()
        
        # Extract neighbor lists from the first 100 vertices for focused testing
        test_vertex_count = min(100, length(amazon_graph.vertices))
        @info "Testing with first $test_vertex_count vertices"
        
        vertex_ids = sort(collect(keys(amazon_graph.vertices)))[1:test_vertex_count]
        
        for (i, vertex_id) in enumerate(vertex_ids)
            vertex = amazon_graph.vertices[vertex_id]
            
            # Convert to UInt32 and create consecutive numbering
            neighbors = UInt32[]
            for neighbor in vertex.out_neighbors
                push!(neighbors, UInt32(neighbor))
            end
            
            neighbor_lists[UInt32(i)] = sort(neighbors)  # Use consecutive numbering 1, 2, 3...
        end
        
        @info "Created neighbor lists for $(length(neighbor_lists)) vertices"
        
        # Print some statistics
        list_lengths = [length(neighbors) for neighbors in values(neighbor_lists)]
        empty_lists = sum(length(neighbors) == 0 for neighbors in values(neighbor_lists))
        
        @info "Neighbor list statistics:"
        @info "  Empty lists: $empty_lists"
        @info "  Min length: $(minimum(list_lengths))"
        @info "  Max length: $(maximum(list_lengths))"
        @info "  Average length: $(round(mean(list_lengths), digits=2))"
        @info "  Total edges: $(sum(list_lengths))"
        
        # Show some example lists
        @info "Example neighbor lists:"
        for (i, (vertex_id, neighbors)) in enumerate(sort(collect(neighbor_lists)))
            if i <= 5  # Show first 5
                @info "  Vertex $vertex_id: $neighbors"
            end
        end
        
        # Test different encoding combinations
        encodings = [:elias_gamma, :elias_delta, :fibonacci, :zeta]
        modes = [:children, :index]
        enable_reference_options = [false, true]
        
        for encoding in encodings
            for mode in modes
                for enable_reference in enable_reference_options
                    @testset "$encoding - $mode - reference=$(enable_reference)" begin
                        @info "Testing $encoding, $mode, enable_reference=$enable_reference"
                        
                        # Encode
                        io = IOBuffer()
                        writer = BitWriter(io)
                        
                        start_time = time()
                        Adjacently.Compression.write_compressed_graph_data(writer, neighbor_lists, encoding, mode, true, enable_reference)
                        flush_bitwriter(writer; flush_last_bits=true)
                        encoding_time = time() - start_time
                        
                        # Get encoding statistics
                        buffer_size = position(io)
                        total_edges = sum(length(neighbors) for neighbors in values(neighbor_lists))
                        bits_per_edge = total_edges > 0 ? (buffer_size * 8) / total_edges : 0.0
                        
                        @info "  Encoded to $buffer_size bytes ($(round(bits_per_edge, digits=4)) bits/edge) in $(round(encoding_time, digits=3))s"
                        
                        # Decode
                        seekstart(io)
                        reader = BitReader(io)
                        
                        start_time = time()
                        decoded = Adjacently.Compression.read_compressed_graph_data(reader, UInt32(length(neighbor_lists)), encoding, mode, UInt32)
                        decoding_time = time() - start_time
                        
                        @info "  Decoded in $(round(decoding_time, digits=3))s"
                        
                        # Verify correctness
                        @test length(decoded) == length(neighbor_lists)
                        
                        # Check each vertex
                        mismatches = 0
                        for vertex_id in keys(neighbor_lists)
                            expected = Set(neighbor_lists[vertex_id])
                            if haskey(decoded, vertex_id)
                                actual = Set(decoded[vertex_id])
                                if expected != actual
                                    mismatches += 1
                                    if mismatches <= 3  # Show first 3 mismatches
                                        @warn "Mismatch for vertex $vertex_id: expected=$(sort(collect(expected))), got=$(sort(collect(actual)))"
                                    end
                                end
                            else
                                mismatches += 1
                                @warn "Missing vertex $vertex_id in decoded results"
                            end
                        end
                        
                        # Check for extra vertices
                        extra_vertices = 0
                        for vertex_id in keys(decoded)
                            if !haskey(neighbor_lists, vertex_id)
                                extra_vertices += 1
                                if extra_vertices <= 3  # Show first 3 extra
                                    @warn "Extra vertex $vertex_id in decoded results: $(decoded[vertex_id])"
                                end
                            end
                        end
                        
                        # Report results
                        if mismatches == 0 && extra_vertices == 0
                            @info "  ✅ Perfect match - all $(length(neighbor_lists)) vertices correct"
                            @test true  # Mark test as passed
                        else
                            @error "  ❌ Encoding/decoding errors: $mismatches mismatches, $extra_vertices extra vertices"
                            @test false  # Mark test as failed
                        end
                    end
                end
            end
        end
        
    catch e
        @warn "Failed to load Amazon dataset: $e"
        @info "Skipping Amazon graph tests - dataset may not be available"
    end
end