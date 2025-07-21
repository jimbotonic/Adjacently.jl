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
