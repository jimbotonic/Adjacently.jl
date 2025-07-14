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

using Test
using Pkg
using Statistics  # Add this import for mean, median, std
Pkg.activate(normpath(joinpath(@__DIR__, "..")))

using Adjacently
using Adjacently.IO: load_adjacency_list_from_csv, load_graph_from_pajek, BitWriter, BitReader, read_bits, flush_bitwriter
using Adjacently.Graph: get_core, get_reverse_graph, get_basic_stats, remap_vertices, relabel_graph
using Adjacently.MGS: write_mgs3_graph, write_compressed_mgs3_graph, load_mgs3_graph, load_compressed_mgs3_graph
using Adjacently.Util: bottom_up_sort, quicksort_iterative_permutation!, get_sorted_array, binary_search
using Adjacently.Compression: huffman_encoding, encode_huffman_tree!, decode_huffman_tree!, get_huffman_codes!, 
write_elias_gamma, write_elias_delta, write_golomb, read_elias_gamma, read_elias_delta, read_golomb, 
write_fibonacci_code, read_fibonacci_code
using LightGraphs: nv, ne, outneighbors, vertices, outdegree, density, add_vertices!, add_edge!

# Get the absolute path to the project root directory
const PROJECT_ROOT = normpath(joinpath(@__DIR__, ".."))

# Helper function to get absolute path for dataset files
function get_dataset_path(filename)
    return normpath(joinpath(PROJECT_ROOT, "datasets", filename))
end

const AMZ_DATASET_IN = get_dataset_path("Amazon_0601/Amazon0601.txt")
const AMZ_DATASET_OUT = "Amazon_0601_core"

const GG_DATASET_IN = get_dataset_path("Web_Google/web-Google.txt")
const GG_DATASET_OUT = "Web_Google_core"

const ARX_DATASET_IN = get_dataset_path("Arxiv_HEP-PH/Cit-HepPh.txt")
const ARX_DATASET_OUT = "Arxiv_HEP-PH_core"

const EAT_DATASET_IN = get_dataset_path("EAT/EATnew.net")
const EAT_DATASET_OUT = "EAT_rcore"

const TEST_DIR = joinpath(PROJECT_ROOT, "test_data")

"""
	load_dataset(input_path::AbstractString; separator::AbstractChar=',', is_pajek::Bool=false)

Load graph from CSV adjacency list or Pajek file
"""
function load_dataset(input_path::String; separator::Char=',', is_pajek::Bool=false)
	if !is_pajek
		g = load_adjacency_list_from_csv(input_path, separator)
	else
	    g = load_graph_from_pajek(input_path)
	end
	return g
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
    @test binary_search(v, UInt8(3)) == -1  # not found
    
    # Test searchsortedfirst
    @test searchsortedfirst(v, UInt8(3)) == 3  # insert position
    @test searchsortedfirst(v, UInt8(4)) == 3  # existing element
    @test searchsortedfirst(v, UInt8(9)) == 6  # beyond end
    @test searchsortedfirst(v, UInt8(0)) == 1  # before start
end

"""
    elias_bits(v::T) where {T<:Unsigned}

    Parameters:
    - v: Value to encode
    - encoding: Encoding scheme (:elias_gamma or :elias_delta)

Return Elias code as a BitVector for a given value.
"""
function elias_bits(v::T, encoding::Symbol) where {T<:Unsigned}
    io = IOBuffer()
    writer = BitWriter(io)
    if encoding == :elias_gamma
        write_elias_gamma(writer, v)
    elseif encoding == :elias_delta
        write_elias_delta(writer, v)
    end
    nbits = writer.index - 1 # capture number of bits before flush resets it
    flush_bitwriter(writer; flush_last_bits=true)
    # read the bits from the buffer
    seekstart(io)
    reader = BitReader(io)
    return read_bits(reader, nbits)
end

"""
    golomb_bits(v::T, k::Int) where {T<:Unsigned}

    Parameters:
    - v: Value to encode
    - k: Number of bits to use for the Golomb coding

Return Golomb code as a BitVector for a given value.
"""
function golomb_bits(v::T, k::Int) where {T<:Unsigned}
    io = IOBuffer()
    writer = BitWriter(io)
    write_golomb(writer, v, k)
    nbits = writer.index - 1 # capture number of bits before flush resets it
    flush_bitwriter(writer; flush_last_bits=true)
    seekstart(io)
    reader = BitReader(io)
    return read_bits(reader, nbits)
end

@testset "Elias Encoding" begin
    @test_throws ArgumentError write_elias_gamma(BitWriter(IOBuffer()), UInt8(0))

    @test elias_bits(UInt8(1), :elias_gamma) == BitVector([1])
    @test elias_bits(UInt8(2), :elias_gamma) == BitVector([0, 1, 0])
    @test elias_bits(UInt8(3), :elias_gamma) == BitVector([0, 1, 1])
    @test elias_bits(UInt8(4), :elias_gamma) == BitVector([0, 0, 1, 0, 0])
    @test elias_bits(UInt8(5), :elias_gamma) == BitVector([0, 0, 1, 0, 1])
    @test elias_bits(UInt8(6), :elias_gamma) == BitVector([0, 0, 1, 1, 0])
    @test elias_bits(UInt8(7), :elias_gamma) == BitVector([0, 0, 1, 1, 1])

    @test elias_bits(UInt8(1), :elias_delta) == BitVector([1])
    @test elias_bits(UInt8(2), :elias_delta) == BitVector([0, 1, 0, 0])
    @test elias_bits(UInt8(3), :elias_delta) == BitVector([0, 1, 0, 1])
    @test elias_bits(UInt8(4), :elias_delta) == BitVector([0, 1, 1, 0, 0])
    @test elias_bits(UInt8(5), :elias_delta) == BitVector([0, 1, 1, 0, 1])
    @test elias_bits(UInt8(6), :elias_delta) == BitVector([0, 1, 1, 1, 0])
    @test elias_bits(UInt8(7), :elias_delta) == BitVector([0, 1, 1, 1, 1])
end

@testset "Golomb Encoding" begin
    @test golomb_bits(UInt8(0), 4) == BitVector([1, 0, 0])
    @test golomb_bits(UInt8(1), 4) == BitVector([1, 0, 1])
    @test golomb_bits(UInt8(2), 4) == BitVector([1, 1, 0])
    @test golomb_bits(UInt8(3), 4) == BitVector([1, 1, 1])
    @test golomb_bits(UInt8(4), 4) == BitVector([0, 1, 0, 0])
    @test golomb_bits(UInt8(5), 4) == BitVector([0, 1, 0, 1])
    @test golomb_bits(UInt8(6), 4) == BitVector([0, 1, 1, 0])
    @test golomb_bits(UInt8(7), 4) == BitVector([0, 1, 1, 1])
end

@testset "Pajek simple graph" begin
    g = load_graph_from_pajek(EAT_DATASET_IN)

    @test nv(g) == 23219
    @test ne(g) == 325593
    @test density(g) == 0.0006039580036712458

    nv1 = UInt16[]

    nv45 = UInt16[0x002d, 0x0037, 0x00ac, 0x00bf, 0x00c7, 
    0x03ed, 0x042a, 0x085a, 0x08be, 0x0ccc, 0x12a3, 0x14bc, 
    0x2278, 0x25a4, 0x2c86, 0x2ffa, 0x3499, 0x35c4, 0x35e1, 
    0x367c, 0x3737, 0x37bb, 0x3fed, 0x3fef, 0x436c, 0x4760, 
    0x4987, 0x54f5, 0x59f0, 0x5aa2]

    nv55 = UInt16[0x0003, 0x0037, 0x003b, 0x0044, 0x0046, 
    0x0068, 0x0079, 0x0091, 0x0097, 0x00a7, 0x00b2, 0x00b7, 
    0x09cf, 0x0add, 0x0e50, 0x0f59, 0x1004, 0x1071, 0x1495, 
    0x14d7, 0x14f2, 0x1833, 0x183a, 0x1a05, 0x1ab4, 0x1e06, 
    0x1e4d, 0x235e, 0x29f9, 0x2fd2, 0x3209, 0x32c8, 0x35c9, 
    0x364d, 0x36ad, 0x36af, 0x36cf, 0x3a2a, 0x3ac6, 0x3c1b, 
    0x3c1e, 0x45f9, 0x4748, 0x4f1f, 0x4ff4, 0x516a, 0x53ea, 
    0x5a51, 0x5a7a, 0x5a7e]

    nv88 = UInt16[0x0059, 0x005b, 0x005c, 0x00a9, 0x00b7, 
    0x0368, 0x0ea9, 0x0eb3, 0x103e, 0x18e9, 0x1b7e, 0x1ebf, 
    0x2078, 0x22a8, 0x22ae, 0x2653, 0x26fb, 0x2d31, 0x2d36, 
    0x2d95, 0x2fb2, 0x307a, 0x36ad, 0x3735, 0x37cb, 0x37ef, 
    0x3a2a, 0x3c91, 0x454e, 0x46cf, 0x4a79, 0x50ca, 0x540c, 
    0x54f6, 0x5998, 0x5a4d, 0x5a51]

    @test outneighbors(g, 1) == nv1
    @test outneighbors(g, 45) == nv45
    @test outneighbors(g, 55) == nv55
    @test outneighbors(g, 88) == nv88

    # save graph in MGSv3 format
    write_mgs3_graph(g, "EAT")

    # load graph in MGSv3 format
    g2 = load_mgs3_graph("EAT.mgs")

    @test nv(g2) == nv(g)
    @test ne(g2) == ne(g)
    @test density(g2) == density(g)

    @test outneighbors(g2, 1) == nv1
    @test outneighbors(g2, 45) == nv45
    @test outneighbors(g2, 55) == nv55
    @test outneighbors(g2, 88) == nv88

    # clean up the test directory
    #rm("EAT.mgs", force=true)
end


@testset "Golomb Encoding (small graph)" begin
    # Create a simple 10-vertex strongly connected directed graph
    # Each vertex connects to the next 3 vertices (wrapping around)
    g = SimpleDiGraph{UInt8}(10)
    
    # create a strongly connected graph where each vertex connects to the next 3 vertices
    for v in 1:10
        for offset in 1:3
            target = mod(v + offset - 1, 10) + 1  # wrap around
            add_edge!(g, v, target)
        end
    end
    
    # verify the graph properties
    @test nv(g) == 10
    @test ne(g) == 30  # each vertex has 3 outgoing edges
    
    # Create test directory if it doesn't exist
    mkpath(TEST_DIR)
    
    # Test both encoding schemes
    encoding_types = [:children, :index]
    
    for encoding_type in encoding_types
        @info("Testing Golomb compression with $encoding_type encoding on small graph")
        
        # create output filename
        output_file = joinpath(TEST_DIR, "small_graph_golomb_$(encoding_type)")
        
        # write the graph using Golomb compression
        write_compressed_mgs3_graph(g, output_file, encoding_type, :golomb)
        
        # load the graph back
        loaded_graph = load_compressed_mgs3_graph(output_file * ".mgz")
        
        # Verify graph properties are preserved
        @test nv(loaded_graph) == nv(g)
        @test ne(loaded_graph) == ne(g)
        @test density(loaded_graph) ≈ density(g)
        
        # Verify specific edges are preserved
        for v in 1:10
            original_neighbors = sort(collect(outneighbors(g, v)))
            loaded_neighbors = sort(collect(outneighbors(loaded_graph, v)))
            @test original_neighbors == loaded_neighbors
        end
        
        # Verify degree statistics
        @test maximum(outdegree(loaded_graph)) == maximum(outdegree(g))
        @test minimum(outdegree(loaded_graph)) == minimum(outdegree(g))
        @test mean(outdegree(loaded_graph)) ≈ mean(outdegree(g))
        @test median(outdegree(loaded_graph)) ≈ median(outdegree(g))
        @test std(outdegree(loaded_graph)) ≈ std(outdegree(g))
        
        @info("Golomb compression with $encoding_type encoding passed all tests")
    end
    
    # Clean up test files
    # rm(joinpath(TEST_DIR, "small_graph_golomb_*"), force=true)
end

@testset "Huffman Encoding" begin
    v = UInt8[1, 6, 3, 7, 2, 8, 5, 18, 12, 17, 13, 24, 12, 1, 4]
    
    # get the Huffman tree for the given values
    t = huffman_encoding(v)
    @test t !== nothing
    
    # encode the Huffman tree
    S = BitArray{1}()
    D = Array{UInt8,1}()
    encode_huffman_tree!(t, S, D)

    # decode the Huffman tree
    t2 = decode_huffman_tree!(S, D)
    @test t2 !== nothing
    
    # get the Huffman codes for the given values
    B = BitArray{1}()
    # dictionary of bitarrays to uint8
    # example: [1, 1, 0, 0, 0, 0] -> 0x05
    C = Dict{BitArray{1},UInt8}()
    get_huffman_codes!(t2, C, B)
    
    # Verify we have one code per unique value (no hidden vertices)
    @test length(keys(C)) == length(unique(v))

    # Verify each unique value has a code
    @test all(val -> any(code_val -> code_val == val, values(C)), unique(v))
    
    # Verify codes are prefix-free (no code is a prefix of another)
    codes = collect(keys(C))
    for i in eachindex(codes)
        for j in eachindex(codes)
            if i != j
                # Convert BitVectors to strings for prefix comparison
                code_i = join(Int.(codes[i]))
                code_j = join(Int.(codes[j]))
                @test !startswith(code_i, code_j)
            end
        end
    end
end

@testset "Remapping Vertices" begin
    g = load_dataset(AMZ_DATASET_IN; separator='\t')
    @test nv(g) == 403394
    @test ne(g) == 3387388
    @test isapprox(density(g), 2.0816473245108078e-5, rtol=1e-10)

    # remap the vertices according to the in-degree
    remapped_vertices = remap_vertices(g, :in_degree)
    remapped_g = relabel_graph(g, remapped_vertices)

    @test nv(remapped_g) == nv(g)
    @test ne(remapped_g) == ne(g)
    @test isapprox(density(remapped_g), density(g), rtol=1e-10)
end

@testset "Amazon Graph Tests" begin
	@info("Loading Amazon dataset")
	amz_g = load_dataset(AMZ_DATASET_IN; separator='\t')
	# number of vertices and edges
	@test 403394 == convert(Int,nv(amz_g))
	@test 3387388 == ne(amz_g)
	
	@info("Getting core")
	amz_core,oni,noi = get_core(amz_g)
	@test 395234 == convert(Int,nv(amz_core))
	@test 3301092 == ne(amz_core)

	@info("Getting reverse graph")
	amz_rcore = get_reverse_graph(amz_core) 
	@test 395234 == convert(Int,nv(amz_rcore))
	@test 3301092 == ne(amz_rcore)
    
    # Create test directory if it doesn't exist
    mkpath(TEST_DIR)

	@info("Saving Amazon dataset (core) in MGS format")
	# NB: the output file is created with extension .mgs
	write_mgs3_graph(amz_core, joinpath(TEST_DIR, AMZ_DATASET_OUT))

	@info("Saving Amazon dataset (core) in MGZ format")
	# NB: the output file is created with extension .mgz
	write_compressed_mgs3_graph(amz_core, joinpath(TEST_DIR, AMZ_DATASET_OUT), :children, :huffman)
	
	@info("Loading Amazon dataset (core) from MGS format")
	# NB: the input file is created with extension .mgs
	amz_core_mgs = load_mgs3_graph(joinpath(TEST_DIR, AMZ_DATASET_OUT * ".mgs"))
	@test 395234 == convert(Int, nv(amz_core_mgs))
	@test 3301092 == ne(amz_core_mgs)

	@info("Loading Amazon dataset (core) from MGZ format")
	# NB: the input file is created with extension .mgz
	amz_core_mgz = load_compressed_mgs3_graph(joinpath(TEST_DIR, AMZ_DATASET_OUT * ".mgz"))
	@test 395234 == convert(Int, nv(amz_core_mgz))
	@test 3301092 == ne(amz_core_mgz)

    # Clean up test files
    #rm(joinpath(TEST_DIR, AMZ_DATASET_OUT * ".mgs"), force=true)
    #rm(joinpath(TEST_DIR, AMZ_DATASET_OUT * ".mgz"), force=true)
    # clean up the test directory
    #rm(TEST_DIR, force=true, recursive=true)  
end

@testset "Amazon Graph Tests (relabeling vertices)" begin
	@info("Loading Amazon dataset")
	amz_g = load_dataset(AMZ_DATASET_IN; separator='\t')
	
	@info("Getting core")
	amz_core,oni,noi = get_core(amz_g)

	@info("Getting reverse graph")
	amz_rcore = get_reverse_graph(amz_core) 

    # relabel core according to in-degree
    remapped_vertices = remap_vertices(amz_core, :in_degree)
    amz_core_relabeled = relabel_graph(amz_core, remapped_vertices)

    # check that the relabeled core has the same number of vertices and edges
    @info("Checking that the relabeled core has the same number of vertices and edges")
    @test nv(amz_core_relabeled) == nv(amz_core)
    @test ne(amz_core_relabeled) == ne(amz_core)

    # relabel reverse graph according to in-degree
    remapped_vertices = remap_vertices(amz_rcore, :in_degree)
    amz_rcore_relabeled = relabel_graph(amz_rcore, remapped_vertices)

    # check that the relabeled reverse core has the same number of vertices and edges
    @info("Checking that the relabeled reverse core has the same number of vertices and edges")
    @test nv(amz_rcore_relabeled) == nv(amz_rcore)
    @test ne(amz_rcore_relabeled) == ne(amz_rcore)

    # save the relabeled graphs
    @info("Saving relabeled core and reverse core (:children, :elias_delta)")
    write_compressed_mgs3_graph(amz_core_relabeled, joinpath(TEST_DIR, AMZ_DATASET_OUT * "_core_elias_delta_relabeled"), :children, :elias_delta)
    write_compressed_mgs3_graph(amz_rcore_relabeled, joinpath(TEST_DIR, AMZ_DATASET_OUT * "_rcore_elias_delta_relabeled"), :children, :elias_delta)
    
    # load the relabeled graphs
    @info("Loading relabeled core and reverse core (:children, :elias_delta)")
    amz_core_relabeled_mgz = load_compressed_mgs3_graph(joinpath(TEST_DIR, AMZ_DATASET_OUT * "_core_elias_delta_relabeled" * ".mgz"))
    amz_rcore_relabeled_mgz = load_compressed_mgs3_graph(joinpath(TEST_DIR, AMZ_DATASET_OUT * "_rcore_elias_delta_relabeled" * ".mgz"))

    # check that the loaded relabeled core and reverse core have the same number of vertices and edges
    @info("Checking that the loaded relabeled core and reverse core have the same number of vertices and edges")
    @test nv(amz_core_relabeled_mgz) == nv(amz_core_relabeled)
    @test ne(amz_core_relabeled_mgz) == ne(amz_core_relabeled)
    @test nv(amz_rcore_relabeled_mgz) == nv(amz_rcore_relabeled)
    @test ne(amz_rcore_relabeled_mgz) == ne(amz_rcore_relabeled)

    # clean up test files
    #rm(joinpath(TEST_DIR, AMZ_DATASET_OUT * "_core_relabeled" * ".mgs"), force=true)
    #rm(joinpath(TEST_DIR, AMZ_DATASET_OUT * "_rcore_relabeled" * ".mgs"), force=true)
end

@testset "Arxiv Graph Tests" begin
	@info("Loading Arxiv_HEP-PH dataset")
	g = load_dataset(ARX_DATASET_IN; separator='\t')
    
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
    
    # Create test directory if it doesn't exist
    mkpath(TEST_DIR)
    
    # Use test directory for output files
    mgs_output_file = joinpath(TEST_DIR, "Arxiv_HEP-PH_test")

    # Test MGS3 format preservation
    # NB: the output file is created with extension .mgs
    write_mgs3_graph(g, mgs_output_file)
    g2 = load_mgs3_graph(mgs_output_file * ".mgs")
    
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
    
    # clean up test file
    #rm(mgs_output_file * ".mgs", force=true)
    # clean up the test directory
    #rm(TEST_DIR, force=true, recursive=true)  
end

@testset "Arxiv_HEP-PH Graph Tests 2" begin
    # First load and save the graph
    @info("Loading Arxiv_HEP-PH dataset")
    g = load_dataset(ARX_DATASET_IN; separator='\t')
    
    # Create test directory if it doesn't exist
    mkpath(TEST_DIR)
    
    # Use test directory for output files
    mgs_output_file = joinpath(TEST_DIR, "Arxiv_HEP-PH_core")
    mgz_output_file = joinpath(TEST_DIR, "Arxiv_HEP-PH_core_compressed")
    
    # Save the graph in MGS format
    # NB: the output file is created with extension .mgs
    write_mgs3_graph(g, mgs_output_file)
    
    # Now test loading the saved file
    @info("Testing MGS3 graph loading and properties")
    # NB: the input file is created with extension .mgs
    g_loaded = load_mgs3_graph(mgs_output_file * ".mgs")
    initial_vertices = nv(g_loaded)
    initial_edges = ne(g_loaded)
    
    # Test reverse graph
    rg = get_reverse_graph(g_loaded)
    @test ne(rg) == ne(g_loaded)
    
    # Test MGS4 format (Huffman compressed) writing and loading
    # NB: the output file is created with extension .mgz
    write_compressed_mgs3_graph(g_loaded, mgz_output_file, :children, :huffman)
    gb = load_compressed_mgs3_graph(mgz_output_file * ".mgz")
    
    # Verify graph properties are preserved
    @test nv(gb) == initial_vertices
    @test ne(gb) == initial_edges
    
    # clean up test files
    #rm(mgs_output_file * ".mgs", force=true)
    #rm(mgz_output_file * ".mgz", force=true)
    # clean up the test directory
    #rm(TEST_DIR, force=true, recursive=true)  
end

@testset "Pajek Graph Format" begin
    # Use proper path construction
    pajek_file = joinpath(PROJECT_ROOT, "datasets", "EAT", "EATnew.net")
    
    # Test loading from Pajek format
    g = load_graph_from_pajek(pajek_file)
    
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

    # create test directory if it doesn't exist
    mkpath(TEST_DIR)
    
    # Use test directory for output files
    mgs_output_file = joinpath(TEST_DIR, "test_eat")
    
    # Test MGS3 format preservation
    # NB: the output file is created with extension .mgs
    write_mgs3_graph(g, mgs_output_file)
    g2 = load_mgs3_graph(mgs_output_file * ".mgs")
    
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
    #rm(mgs_output_file * ".mgs", force=true)
    # clean up the test directory
    #rm(TEST_DIR, force=true, recursive=true)  
end

@testset "Elias Compression" begin
    @info("Loading Amazon dataset")
	amz_g = load_dataset(AMZ_DATASET_IN; separator='\t')
	# number of vertices and edges
	@test 403394 == convert(Int,nv(amz_g))
	@test 3387388 == ne(amz_g)
	
	@info("Getting core")
	amz_core, oni, noi = get_core(amz_g)
	@test 395234 == convert(Int,nv(amz_core))
	@test 3301092 == ne(amz_core)

    # create test directory if it doesn't exist
    mkpath(TEST_DIR)
    
    # Use test directory for output files
    mgs_output_file = joinpath(TEST_DIR, "amz_core")

    encoding_types = [:children, :index]
    compression_types = [:elias_gamma, :elias_delta]

    for encoding_type in encoding_types
        for compression_type in compression_types
            @info("Writing compressed MGS3 graph with $encoding_type encoding and $compression_type compression")
            write_compressed_mgs3_graph(amz_core, mgs_output_file * "_" * string(encoding_type) * "_" * string(compression_type), encoding_type, compression_type)
        end
    end

    # Test MGS3 format preservation
    # NB: the output file is created with extension .mgs
    for encoding_type in encoding_types
        for compression_type in compression_types
            @info("Loading compressed MGS3 graph with $encoding_type encoding and $compression_type compression")
            loaded_graph = load_compressed_mgs3_graph(mgs_output_file * "_" * string(encoding_type) * "_" * string(compression_type) * ".mgz")
            
            @info("Verifying graph properties")
            @test nv(loaded_graph) == nv(amz_core)
            @test ne(loaded_graph) == ne(amz_core)
            @test density(loaded_graph) ≈ density(amz_core)
            @test maximum(outdegree(loaded_graph)) == maximum(outdegree(amz_core))
            @test minimum(outdegree(loaded_graph)) == minimum(outdegree(amz_core))
            @test mean(outdegree(loaded_graph)) ≈ mean(outdegree(amz_core))
            @test median(outdegree(loaded_graph)) ≈ median(outdegree(amz_core))
            @test std(outdegree(loaded_graph)) ≈ std(outdegree(amz_core))
        end
    end

    # clean up the test directory
    # rm(TEST_DIR, force=true, recursive=true)  
end

@testset "Golomb Compression" begin
    @info("Loading Amazon dataset")
	amz_g = load_dataset(AMZ_DATASET_IN; separator='\t')
	# number of vertices and edges
	@test 403394 == convert(Int,nv(amz_g))
	@test 3387388 == ne(amz_g)
	
	@info("Getting core")
	amz_core, oni, noi = get_core(amz_g)
	@test 395234 == convert(Int,nv(amz_core))
	@test 3301092 == ne(amz_core)

    # create test directory if it doesn't exist
    mkpath(TEST_DIR)
    
    # Use test directory for output files
    mgs_output_file = joinpath(TEST_DIR, "amz_core")

    encoding_types = [:children, :index]
    compression_types = [:golomb]

    for encoding_type in encoding_types
        for compression_type in compression_types
            @info("Writing compressed MGS3 graph with $encoding_type encoding and $compression_type compression")
            write_compressed_mgs3_graph(amz_core, mgs_output_file * "_" * string(encoding_type) * "_" * string(compression_type), encoding_type, compression_type)
        end
    end

    # Test MGS3 format preservation
    # NB: the output file is created with extension .mgs
    for encoding_type in encoding_types
        for compression_type in compression_types
            @info("Loading compressed MGS3 graph with $encoding_type encoding and $compression_type compression")
            loaded_graph = load_compressed_mgs3_graph(mgs_output_file * "_" * string(encoding_type) * "_" * string(compression_type) * ".mgz")
            
            @info("Verifying graph properties")
            @test nv(loaded_graph) == nv(amz_core)
            @test ne(loaded_graph) == ne(amz_core)
            @test density(loaded_graph) ≈ density(amz_core)
            @test maximum(outdegree(loaded_graph)) == maximum(outdegree(amz_core))
            @test minimum(outdegree(loaded_graph)) == minimum(outdegree(amz_core))
            @test mean(outdegree(loaded_graph)) ≈ mean(outdegree(amz_core))
            @test median(outdegree(loaded_graph)) ≈ median(outdegree(amz_core))
            @test std(outdegree(loaded_graph)) ≈ std(outdegree(amz_core))
        end
    end

    # clean up the test directory
    # rm(TEST_DIR, force=true, recursive=true)  
end

@testset "Fibonacci Compression" begin
    @info("Loading Amazon dataset")
	amz_g = load_dataset(AMZ_DATASET_IN; separator='\t')
	# number of vertices and edges
	@test 403394 == convert(Int,nv(amz_g))
	@test 3387388 == ne(amz_g)
	
	@info("Getting core")
	amz_core, oni, noi = get_core(amz_g)
	@test 395234 == convert(Int,nv(amz_core))
	@test 3301092 == ne(amz_core)

    # create test directory if it doesn't exist
    mkpath(TEST_DIR)
    
    # Use test directory for output files
    mgs_output_file = joinpath(TEST_DIR, "amz_core")

    encoding_types = [:children, :index]
    compression_types = [:fibonacci]

    for encoding_type in encoding_types
        for compression_type in compression_types
            @info("Writing compressed MGS3 graph with $encoding_type encoding and $compression_type compression")
            write_compressed_mgs3_graph(amz_core, mgs_output_file * "_" * string(encoding_type) * "_" * string(compression_type), encoding_type, compression_type)
        end
    end

    # Test MGS3 format preservation
    # NB: the output file is created with extension .mgs
    for encoding_type in encoding_types
        for compression_type in compression_types
            @info("Loading compressed MGS3 graph with $encoding_type encoding and $compression_type compression")
            loaded_graph = load_compressed_mgs3_graph(mgs_output_file * "_" * string(encoding_type) * "_" * string(compression_type) * ".mgz")
            
            @info("Verifying graph properties")
            @test nv(loaded_graph) == nv(amz_core)
            @test ne(loaded_graph) == ne(amz_core)
            @test density(loaded_graph) ≈ density(amz_core)
            @test maximum(outdegree(loaded_graph)) == maximum(outdegree(amz_core))
            @test minimum(outdegree(loaded_graph)) == minimum(outdegree(amz_core))
            @test mean(outdegree(loaded_graph)) ≈ mean(outdegree(amz_core))
            @test median(outdegree(loaded_graph)) ≈ median(outdegree(amz_core))
            @test std(outdegree(loaded_graph)) ≈ std(outdegree(amz_core))
        end
    end

    # clean up the test directory
    # rm(TEST_DIR, force=true, recursive=true)  
end

@testset "Remapping vertices" begin
    @info("Loading Amazon dataset")
	amz_g = load_dataset(AMZ_DATASET_IN; separator='\t')
	
	@info("Getting core")
	amz_core,oni,noi = get_core(amz_g)

    # create test directory if it doesn't exist
    mkpath(TEST_DIR)
    
    # Use test directory for output files
    mgs_output_file = joinpath(TEST_DIR, "amz_core_rl")

    criterion = [:in_degree, :out_degree, :degree, :pagerank]
    encoding = [:children, :index]
    compression = [:elias_delta, :fibonacci]

    for compression in compression
        for encoding in encoding
            for criterion in criterion
                remapped_vertices = remap_vertices(amz_core, criterion)
                amz_core_rl = relabel_graph(amz_core, remapped_vertices)

                @info("Checking that the relabeled core has the same number of vertices and edges")
                @test nv(amz_core_rl) == nv(amz_core)
                @test ne(amz_core_rl) == ne(amz_core)

                @info("Writing compressed MGS3 graph with $criterion criterion, $encoding encoding and $compression compression")
                write_compressed_mgs3_graph(amz_core_rl, mgs_output_file * "_" * string(criterion) * "_" * string(encoding) * "_" * string(compression), encoding, compression)
            end
        end
    end

    # clean up the test directory
    # rm(TEST_DIR, force=true, recursive=true)  
end

