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
using Adjacently.Graph: get_core, get_reverse_graph, get_basic_stats, relabel_vertices, relabel_graph
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
