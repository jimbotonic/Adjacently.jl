#!/usr/bin/env julia

#
# Adjacently: Julia Complex Directed Networks Library
# Copyright (C) 2016-2026 Jimmy Dubuisson <jimmy@dubuisson.ch>
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

using Test

# Try to include the full module, but fall back to just index.jl if dependencies are missing
try
    include("../src/Adjacently.jl")
    using .Adjacently
    using .Adjacently.Index
    @info "Using full Adjacently module"
catch e
    @warn "Could not load full Adjacently module due to dependencies, loading Index module directly"
    include("../src/index.jl")
    using .Index
end

@testset "StreamingInvertedIndex Tests" begin
    
    @testset "Basic Index Creation and Operations" begin
        T = UInt32
        
        @testset "Index initialization" begin
            # Test basic index creation
            ncols = 10
            idx = Index.StreamInvertedIndex{T}(ncols)
            
            @test idx isa Index.StreamInvertedIndex{T}
            @test length(idx.postings) == ncols
            @test length(idx.rowdeg) == 0  # No rows added yet
            @test idx.k == 0
            @test idx.ncols == ncols
            
            @info "✓ Basic index creation works"
        end
        
        @testset "Manual candidate addition" begin
            idx = Index.StreamInvertedIndex{T}(10)
            
            # Add arbitrary candidate
            candidate_row = T[2, 5, 6, 9]
            id1 = Index.add_candidate!(idx, candidate_row)
            
            @test id1 == 1
            @test length(idx.rowdeg) == 1
            @test idx.rowdeg[id1] == length(candidate_row)
            @test idx.k == 1
            
            # Add another candidate
            candidate_row2 = T[1, 3, 7]
            id2 = Index.add_candidate!(idx, candidate_row2)
            
            @test id2 == 2
            @test length(idx.rowdeg) == 2
            @test idx.rowdeg[id2] == length(candidate_row2)
            @test idx.k == 2
            
            @info "✓ Manual candidate addition works"
            @info "  Added candidates with degrees: $(idx.rowdeg)"
        end
        
        @testset "Simulated graph-based candidate addition" begin
            # Simulate the graph structure from the original example
            # Graph edges: 1->[2,5,9], 2->[5,7], 3->[4,5,9,10], 4->[2,5,8] 
            graph_neighbors = Dict(
                1 => T[2, 5, 9],
                2 => T[5, 7], 
                3 => T[4, 5, 9, 10],
                4 => T[2, 5, 8]
            )
            
            idx = Index.StreamInvertedIndex{T}(10)
            
            # Add candidates manually using the neighbor lists
            id1 = Index.add_candidate!(idx, graph_neighbors[1])
            id2 = Index.add_candidate!(idx, graph_neighbors[2])
            id3 = Index.add_candidate!(idx, graph_neighbors[3])
            id4 = Index.add_candidate!(idx, graph_neighbors[4])
            
            @test id1 == 1
            @test id2 == 2
            @test id3 == 3
            @test id4 == 4
            @test length(idx.rowdeg) == 4
            
            # Check that neighbors are correctly added
            @test idx.rowdeg[id1] == length(graph_neighbors[1])
            @test idx.rowdeg[id3] == length(graph_neighbors[3])
            
            @info "✓ Simulated graph-based candidate addition works"
            @info "  Added $(length(idx.rowdeg)) candidates with degrees: $(idx.rowdeg)"
        end
    end
    
    @testset "Overlap Computation and Querying" begin
        T = UInt32
        
        @testset "Basic overlap queries" begin
            # Setup using simulated graph structure
            graph_neighbors = Dict(
                1 => T[2, 5, 9],
                2 => T[5, 7], 
                3 => T[4, 5, 9, 10],
                4 => T[2, 5, 8]
            )
            
            idx = Index.StreamInvertedIndex{T}(10)
            
            id1 = Index.add_candidate!(idx, graph_neighbors[1])
            id2 = Index.add_candidate!(idx, graph_neighbors[2])
            id3 = Index.add_candidate!(idx, graph_neighbors[3])
            id4 = Index.add_candidate!(idx, graph_neighbors[4])
            
            # Add manual candidate as in original example
            id5 = Index.add_candidate!(idx, T[2, 5, 6, 9])
            
            # Create workspace
            work = Index.OverlapWorkspace(idx)
            @test work isa Index.OverlapWorkspace{T}
            
            # Query target as in example
            target_children = T[2, 5, 9]
            counts, touched = Index.overlap!(idx, target_children, work)
            
            @test counts isa Vector
            @test touched isa Vector
            @test length(counts) >= length(idx.rowdeg)
            @test length(touched) <= length(idx.rowdeg)
            
            @info "✓ Basic overlap computation works"
            @info "  Query target: $target_children"
            @info "  Touched candidates: $(length(touched))"
            if !isempty(touched)
                @info "  Overlap counts for touched: $([counts[i] for i in touched])"
            end
        end
        
        @testset "Best match and top-K selection" begin
            # Setup same as above
            graph_neighbors = Dict(
                1 => T[2, 5, 9],
                2 => T[5, 7], 
                3 => T[4, 5, 9, 10],
                4 => T[2, 5, 8]
            )
            
            idx = Index.StreamInvertedIndex{T}(10)
            
            Index.add_candidate!(idx, graph_neighbors[1])
            Index.add_candidate!(idx, graph_neighbors[2])
            Index.add_candidate!(idx, graph_neighbors[3])
            Index.add_candidate!(idx, graph_neighbors[4])
            Index.add_candidate!(idx, T[2, 5, 6, 9])
            
            work = Index.OverlapWorkspace(idx)
            target_children = T[2, 5, 9]
            counts, touched = Index.overlap!(idx, target_children, work)
            
            # Test best match
            if !isempty(touched)
                best_i, best_v = Index.argmax_on_touched(counts, touched)
                @test best_i isa Integer
                @test best_v isa Integer
                @test best_v >= 0
                @test best_i in touched
                
                @info "✓ Best match selection works"
                @info "  Best candidate: $best_i with overlap: $best_v"
            end
            
            # Test top-K (even if we have fewer than K candidates)
            k = min(3, length(touched))
            if k > 0
                topk = Index.topk_on_touched(counts, touched, k)
                @test length(topk) <= k
                @test length(topk) <= length(touched)
                
                # Check that top-k is sorted in descending order of overlap
                if length(topk) > 1
                    for i in 1:(length(topk)-1)
                        @test topk[i][2] >= topk[i+1][2]  # (idx, count) tuples
                    end
                end
                
                @info "✓ Top-K selection works"
                @info "  Top-$k candidates: $topk"
            end
        end
        
        @testset "Jaccard similarity computation" begin
            # Setup same as above
            graph_neighbors = Dict(
                1 => T[2, 5, 9],
                2 => T[5, 7], 
                3 => T[4, 5, 9, 10],
                4 => T[2, 5, 8]
            )
            
            idx = Index.StreamInvertedIndex{T}(10)
            
            Index.add_candidate!(idx, graph_neighbors[1])
            Index.add_candidate!(idx, graph_neighbors[2])
            Index.add_candidate!(idx, graph_neighbors[3])
            Index.add_candidate!(idx, graph_neighbors[4])
            Index.add_candidate!(idx, T[2, 5, 6, 9])
            
            work = Index.OverlapWorkspace(idx)
            target_children = T[2, 5, 9]
            counts, touched = Index.overlap!(idx, target_children, work)
            
            # Compute Jaccard similarities
            d_t = T(length(target_children))
            scores = Index.jaccard_on_touched(counts, touched, idx.rowdeg, d_t)
            
            @test scores isa Dict
            
            # Check that Jaccard scores are in [0, 1] for touched candidates
            for i in touched
                @test haskey(scores, i)
                @test 0.0 <= scores[i] <= 1.0
                
                # Verify Jaccard formula: |A ∩ B| / |A ∪ B|
                overlap = counts[i]
                union_size = idx.rowdeg[i] + d_t - overlap
                expected_jaccard = union_size > 0 ? Float32(overlap) / Float32(union_size) : 0.0f0
                @test abs(scores[i] - expected_jaccard) < 1e-6
            end
            
            @info "✓ Jaccard similarity computation works"
            @info "  Target size: $d_t"
            if !isempty(touched)
                @info "  Jaccard scores for touched: $([scores[i] for i in touched])"
            end
        end
    end
    
    @testset "Edge Cases and Error Handling" begin
        T = UInt32
        
        @testset "Empty operations" begin
            idx = Index.StreamInvertedIndex{T}(10)
            work = Index.OverlapWorkspace(idx)
            
            # Query on empty index
            target_children = T[1, 2, 3]
            counts, touched = Index.overlap!(idx, target_children, work)
            
            @test isempty(touched)
            @test length(counts) >= 0
            
            @info "✓ Empty index queries handled correctly"
        end
        
        @testset "Single candidate operations" begin
            idx = Index.StreamInvertedIndex{T}(10)
            
            # Add single candidate
            id1 = Index.add_candidate!(idx, T[2, 5, 7])
            @test id1 == 1
            
            work = Index.OverlapWorkspace(idx)
            
            # Query with exact match
            counts, touched = Index.overlap!(idx, T[2, 5, 7], work)
            @test length(touched) == 1
            @test touched[1] == 1
            @test counts[1] == 3  # Full overlap
            
            # Query with partial match
            counts, touched = Index.overlap!(idx, T[2, 5], work)
            @test length(touched) == 1
            @test counts[1] == 2  # Partial overlap
            
            # Query with no match
            counts, touched = Index.overlap!(idx, T[1, 3, 4], work)
            @test isempty(touched)
            
            @info "✓ Single candidate operations work correctly"
        end
        
        @testset "Large overlap values" begin
            idx = Index.StreamInvertedIndex{T}(100)
            
            # Add candidate with many elements
            large_candidate = T[i for i in 1:50]
            id1 = Index.add_candidate!(idx, large_candidate)
            
            work = Index.OverlapWorkspace(idx)
            
            # Query with large overlap
            large_query = T[i for i in 25:75]  # 26 overlapping elements (25-50)
            counts, touched = Index.overlap!(idx, large_query, work)
            
            @test length(touched) == 1
            @test counts[1] == 26  # Elements 25-50 overlap
            
            # Test Jaccard with large sets
            d_t = T(length(large_query))
            scores = Index.jaccard_on_touched(counts, touched, idx.rowdeg, d_t)
            
            expected_jaccard = Float32(26) / Float32(50 + 51 - 26)  # |A ∩ B| / |A ∪ B|
            @test abs(scores[1] - expected_jaccard) < 1e-6
            
            @info "✓ Large overlap computations work correctly"
            @info "  Large candidate size: $(length(large_candidate))"
            @info "  Large query size: $(length(large_query))"
            @info "  Overlap: $(counts[1])"
            @info "  Jaccard: $(scores[1])"
        end
    end
    
    @testset "Integration with Reference Encoding" begin
        T = UInt32
        
        @testset "StreamingInvertedIndex as reference selection backend" begin
            # Test if the streaming index could be used as a backend for reference encoding
            
            # Create test data similar to reference encoding scenarios
            test_neighbor_lists = Dict{T, Vector{T}}(
                T(1) => T[10, 20, 30, 40],
                T(2) => T[50, 60, 70, 80, 90],
                T(3) => T[10, 20, 30, 40, 100],  # Overlaps with vertex 1
                T(4) => T[50, 60, 70, 80, 90, 110],  # Overlaps with vertex 2
                T(5) => T[15, 25, 35],  # No significant overlap
            )
            
            # Test reference selection for vertex 3 (query a modified version to avoid self-match)
            # Query with vertex 3's neighbors but exclude vertex 3 from candidates
            target_neighbors = test_neighbor_lists[T(3)]
            
            # Create a new index without the target vertex itself
            idx_filtered = Index.StreamInvertedIndex{T}(200)
            candidate_map_filtered = Dict{T, T}()
            ref_v_min_degree = 4
            
            for (vertex_id, neighbors) in test_neighbor_lists
                if length(neighbors) >= ref_v_min_degree && vertex_id != T(3)  # Exclude target vertex
                    stream_id = Index.add_candidate!(idx_filtered, neighbors)
                    candidate_map_filtered[stream_id] = vertex_id
                end
            end
            
            work_filtered = Index.OverlapWorkspace(idx_filtered)
            counts, touched = Index.overlap!(idx_filtered, target_neighbors, work_filtered)
            
            if !isempty(touched)
                best_i, best_overlap = Index.argmax_on_touched(counts, touched)
                best_vertex_id = candidate_map_filtered[best_i]
                
                @info "✓ StreamingInvertedIndex reference selection works"
                @info "  Target vertex: 3 with neighbors: $target_neighbors"
                @info "  Best reference: vertex $best_vertex_id (stream_id $best_i) with overlap: $best_overlap"
                @info "  Reference neighbors: $(test_neighbor_lists[best_vertex_id])"
                
                # Should find vertex 1 with overlap of 4 (elements [10, 20, 30, 40])
                @test best_vertex_id == T(1)
                @test best_overlap == 4
                
                # Test Jaccard similarity
                d_t = T(length(target_neighbors))
                scores = Index.jaccard_on_touched(counts, touched, idx_filtered.rowdeg, d_t)
                jaccard_score = scores[best_i]
                
                expected_jaccard = Float32(4) / Float32(4 + 5 - 4)  # |intersection| / |union|
                @test abs(jaccard_score - expected_jaccard) < 1e-6
                
                @info "  Jaccard similarity: $jaccard_score"
            end
            
            # Test for vertex 4 (should find vertex 2, excluding vertex 4 from candidates)
            target_neighbors = test_neighbor_lists[T(4)]
            
            idx_filtered2 = Index.StreamInvertedIndex{T}(200)
            candidate_map_filtered2 = Dict{T, T}()
            
            for (vertex_id, neighbors) in test_neighbor_lists
                if length(neighbors) >= ref_v_min_degree && vertex_id != T(4)  # Exclude target vertex
                    stream_id = Index.add_candidate!(idx_filtered2, neighbors)
                    candidate_map_filtered2[stream_id] = vertex_id
                end
            end
            
            work_filtered2 = Index.OverlapWorkspace(idx_filtered2)
            counts, touched = Index.overlap!(idx_filtered2, target_neighbors, work_filtered2)
            
            if !isempty(touched)
                best_i, best_overlap = Index.argmax_on_touched(counts, touched)
                best_vertex_id = candidate_map_filtered2[best_i]
                
                @test best_vertex_id == T(2)
                @test best_overlap == 5
                
                @info "  Target vertex: 4 found reference: vertex $best_vertex_id with overlap: $best_overlap"
            end
        end
    end
end

@info "All StreamingInvertedIndex tests completed successfully!"