#
# Adjacently: Julia Complex Directed Networks Library
# Copyright (C) 2016-2024 Jimmy Dubuisson <jimmy.dubuisson@gmail.com>
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

using Test, LightGraphs, Adjacently.io, Adjacently.graph, Adjacently.util

const GRAY_TO_RSA_TABLE_PATH = joinpath(@__DIR__, "data", "gray_to_rsa.txt")

const AMZ_DATASET_IN = joinpath(@__DIR__, "..", "datasets", "Amazon_0601", "Amazon0601.txt")
const AMZ_DATASET_OUT = "Amazon_0601_core"

const GG_DATASET_IN = joinpath(@__DIR__, "..", "datasets", "Web_Google", "web-Google.txt")
const GG_DATASET_OUT = "Web_Google_core"

const ARX_DATASET_IN = joinpath(@__DIR__, "..", "datasets", "Arxiv_HEP-PH", "Cit-HepPh.txt")
const ARX_DATASET_OUT = "Arxiv_HEP-PH_core"

const EAT_DATASET_IN = joinpath(@__DIR__, "..", "datasets", "EAT", "EATnew.net")
const EAT_DATASET_OUT = "EAT_rcore"

"""
	load_dataset(input_path::AbstractString; separator::AbstractChar=',', is_pajek::Bool=false)

Load graph from CSV adjacency list or Pajek file
"""
function load_dataset(input_path::AbstractString; separator::AbstractChar=',', is_pajek::Bool=false)
	if !is_pajek
		g = load_adjacency_list_from_csv(input_path, separator)
	else
	    g = load_graph_from_pajek(input_path)
	end
	return g
end

@testset "Amazon_0601 Graph Tests" begin
	@info("Loading Amazon dataset")
	amz_g = load_dataset(AMZ_DATASET_IN; separator='\t')
	# number of vertices and edges
	@test 403394 == convert(Int,nv(amz_g))
	@test 3387388 == ne(amz_g)
	
	@info("Getting core")
	amz_core,oni,noi = get_core(amz_g)
	@test 395234 == convert(Int,nv(amz_core))
	@test 3301092 == ne(amz_core)

	@info("getting reverse graph")
	amz_rcore = get_reverse_graph(amz_core) 
	@test 395234 == convert(Int,nv(amz_rcore))
	@test 3301092 == ne(amz_rcore)

	@info("Saving Amazon dataset (core) in MGS format")
	write_mgs3_graph(amz_core, AMZ_DATASET_OUT)

	@info("Saving Amazon dataset (core) in MGZ format")
	write_mgs3_huffman_graph(amz_core, amz_rcore, AMZ_DATASET_OUT)
	# serialize_to_jld(amz_core, "core", AMZ_DATASET_OUT)
	
	@info("Loading Amazon dataset (core) from MGS format")
	amz_core_mgs = load_mgs3_graph(AMZ_DATASET_OUT * ".mgs")
	@test 395234 == convert(Int,nv(amz_core_mgs))
	@test 3301092 == ne(amz_core_mgs)

	@info("Loading Amazon dataset (core) from MGZ format")
	amz_core_mgz = load_mgs3_huffman_graph(AMZ_DATASET_OUT * ".mgz")
	@test 395234 == convert(Int,nv(amz_core_mgz))
	@test 3301092 == ne(amz_core_mgz)
end

@testset "Arxiv_HEP-PH Graph Tests" begin
    # Load the test graph
    g = load_adjacency_list_from_csv("../datasets/Arxiv_HEP-PH/Cit-HepPh.txt", '\t')
    
    # Test basic graph properties
    @test nv(g) == 34546
    @test ne(g) == 421578
    @test isapprox(density(g), 0.00035326041393102855, rtol=1e-10)
    
    # Test degree statistics
    @test maximum(outdegree(g)) == 411
    @test minimum(outdegree(g)) == 0
    @test isapprox(mean(outdegree(g)), 12.203380999247381, rtol=1e-10)
    @test isapprox(median(outdegree(g)), 8.0, rtol=1e-10)
    @test isapprox(std(outdegree(g)), 15.224473754546942, rtol=1e-10)
    
    # Test specific vertex neighborhood
    expected_neighbors = UInt16[2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12]
    @test outneighbors(g, vertices(g)[1]) == expected_neighbors
    
    # Test MGS3 format preservation
    write_mgs3_graph(g, "test_graph")
    g2 = load_mgs3_graph("test_graph.mgs")
    
    # Verify graph properties are preserved
    @test nv(g) == nv(g2)
    @test ne(g) == ne(g2)
    @test density(g) ≈ density(g2)
    @test maximum(outdegree(g)) == maximum(outdegree(g2))
    @test minimum(outdegree(g)) == minimum(outdegree(g2))
    @test mean(outdegree(g)) ≈ mean(outdegree(g2))
    @test median(outdegree(g)) ≈ median(outdegree(g2))
    @test std(outdegree(g)) ≈ std(outdegree(g2))
    @test outneighbors(g, vertices(g)[1]) == outneighbors(g2, vertices(g2)[1])
    
    # Clean up test file
    rm("test_graph.mgs", force=true)
end

@testset "Sorting Algorithms" begin
    A = UInt8[3,6,4,7,1,9,2,12,10]
    PA = UInt8[5, 7, 1, 3, 2, 4, 6, 9, 8]
    SA1 = UInt8[1, 2, 3, 4, 6, 7, 9, 10, 12] 
    SA2 = UInt8[12, 10, 9, 7, 6, 4, 3, 2, 1]
    
    # Test mergesort
    R = bottom_up_sort(A)
    @test R == PA
    
    # Test quicksort
    A2 = copy(A)
    R2 = quicksort_iterative_permutation!(A2)
    @test R2 == PA
    
    # Test sorted arrays
    S = get_sorted_array(A, R)
    S2 = get_sorted_array(A, R, false)
    @test S == SA1
    @test S2 == SA2
end

@testset "Binary Search" begin
    v = UInt8[1, 2, 4, 6, 8]
    
    # Test binary search
    @test binary_search(v, UInt8(3)) == 0  # not found
    
    # Test searchsortedfirst
    @test searchsortedfirst(v, UInt8(3)) == 3  # insert position
    @test searchsortedfirst(v, UInt8(4)) == 3  # existing element
    @test searchsortedfirst(v, UInt8(9)) == 6  # beyond end
    @test searchsortedfirst(v, UInt8(0)) == 1  # before start
end

@testset "Huffman Encoding" begin
    v = UInt8[1, 6, 3, 7, 2, 8, 5, 18, 12, 17, 13, 24, 12, 1, 4]
    
    # Test Huffman tree creation
    t = huffman_encoding(v)
    @test t !== nothing
    
    # Test tree encoding/decoding
    S = BitArray{1}()
    D = Array{UInt8,1}()
    encode_tree!(t, S, D)
    
    t2 = decode_tree!(S, D)
    @test t2 !== nothing
    
    # Test Huffman codes
    B = BitArray{1}()
    C = Dict{BitArray{1},UInt8}()
    get_huffman_codes!(t2, C, B)
    
    # Verify unique codes for each value
    unique_values = unique(v)
    @test length(keys(C)) == length(unique_values)
    @test all(val -> any(code_val -> code_val == val, values(C)), unique_values)
end

@testset "Arxiv_HEP-PH Graph Tests 2" begin
    # Test MGS3 graph loading and properties
    g = load_mgs3_graph("../datasets/Arxiv_HEP-PH/Arxiv_HEP-PH_core.mgs")
    initial_vertices = nv(g)
    initial_edges = ne(g)
    
    # Test reverse graph
    rg = get_reverse_graph(g)
    @test ne(rg) == ne(g)
    
    # Test MGS4 format (Huffman compressed) writing and loading
    write_mgs3_huffman_graph(g, rg, "test")
    gb = load_mgs3_huffman_graph("test.mgz")
    
    # Verify graph properties are preserved
    @test nv(gb) == initial_vertices
    @test ne(gb) == initial_edges
    
    # Clean up
    rm("test.mgz", force=true)
end

@testset "Pajek Graph Format" begin
    # Test loading from Pajek format
    g = load_graph_from_pajek("../datasets/EAT/EATnew.net")
    
    # Test basic statistics
    nvs, nes, dens = get_basic_stats(g)
    @test nvs == 23219
    @test nes == 325593
    @test isapprox(dens, 0.0006039580036712458, atol=1e-15)
    
    # Test specific vertex neighborhoods
    @test isempty(outneighbors(g, 1))  # nv1 is empty
    
    # Test vertex 45 neighbors
    expected_nv45 = UInt16[0x002d, 0x0037, 0x00ac, 0x00bf, 0x00c7, 
        0x03ed, 0x042a, 0x085a, 0x08be, 0x0ccc, 0x12a3, 0x14bc, 
        0x2278, 0x25a4, 0x2c86, 0x2ffa, 0x3499, 0x35c4, 0x35e1, 
        0x367c, 0x3737, 0x37bb, 0x3fed, 0x3fef, 0x436c, 0x4760, 
        0x4987, 0x54f5, 0x59f0, 0x5aa2]
    @test sort(outneighbors(g, 45)) == sort(expected_nv45)
    
    # Test vertex 55 neighbors
    expected_nv55 = UInt16[0x0003, 0x0037, 0x003b, 0x0044, 0x0046, 
        0x0068, 0x0079, 0x0091, 0x0097, 0x00a7, 0x00b2, 0x00b7, 
        0x09cf, 0x0add, 0x0e50, 0x0f59, 0x1004, 0x1071, 0x1495, 
        0x14d7, 0x14f2, 0x1833, 0x183a, 0x1a05, 0x1ab4, 0x1e06, 
        0x1e4d, 0x235e, 0x29f9, 0x2fd2, 0x3209, 0x32c8, 0x35c9, 
        0x364d, 0x36ad, 0x36af, 0x36cf, 0x3a2a, 0x3ac6, 0x3c1b, 
        0x3c1e, 0x45f9, 0x4748, 0x4f1f, 0x4ff4, 0x516a, 0x53ea, 
        0x5a51, 0x5a7a, 0x5a7e]
    @test sort(outneighbors(g, 55)) == sort(expected_nv55)
    
    # Test vertex 88 neighbors
    expected_nv88 = UInt16[0x0059, 0x005b, 0x005c, 0x00a9, 0x00b7, 
        0x0368, 0x0ea9, 0x0eb3, 0x103e, 0x18e9, 0x1b7e, 0x1ebf, 
        0x2078, 0x22a8, 0x22ae, 0x2653, 0x26fb, 0x2d31, 0x2d36, 
        0x2d95, 0x2fb2, 0x307a, 0x36ad, 0x3735, 0x37cb, 0x37ef, 
        0x3a2a, 0x3c91, 0x454e, 0x46cf, 0x4a79, 0x50ca, 0x540c, 
        0x54f6, 0x5998, 0x5a4d, 0x5a51]
    @test sort(outneighbors(g, 88)) == sort(expected_nv88)
    
    # Test MGS3 format preservation
    write_mgs3_graph(g, "test_eat")
    g2 = load_mgs3_graph("test_eat.mgs")
    
    # Verify graph properties are preserved
    nvs2, nes2, dens2 = get_basic_stats(g2)
    @test nvs2 == nvs
    @test nes2 == nes
    @test isapprox(dens2, dens, atol=1e-15)
    
    # Verify specific vertex neighborhoods are preserved
    @test isempty(outneighbors(g2, 1))
    @test sort(outneighbors(g2, 45)) == sort(expected_nv45)
    @test sort(outneighbors(g2, 55)) == sort(expected_nv55)
    @test sort(outneighbors(g2, 88)) == sort(expected_nv88)
    
    # Clean up test file
    rm("test_eat.mgs", force=true)
end


