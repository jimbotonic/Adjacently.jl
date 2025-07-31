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
using Statistics 

Pkg.activate(normpath(joinpath(@__DIR__, "..")))

using LightGraphs: nv, ne, outneighbors, vertices, outdegree, density, add_vertices!, add_edge!
using StatsBase
using Plots

using Adjacently

using Adjacently.IO: load_adjacency_list_from_csv, load_graph_from_pajek, BitWriter, BitReader, 
read_bits, flush_bitwriter, write_bit, read_bit

using Adjacently.Graph: get_core, get_reverse_graph, get_basic_stats, relabel_graph, relabel_vertices
using Adjacently.MGS: write_mgs3_graph, write_compressed_mgs3_graph, load_mgs3_graph, load_compressed_mgs3_graph
using Adjacently.Util: bottom_up_sort, quicksort_iterative_permutation!, get_sorted_array, binary_search

using Adjacently.Compression: write_unary_coding, write_truncated_binary_coding, huffman_encoding, 
encode_huffman_tree!, decode_huffman_tree!, get_huffman_codes!, 
write_elias_gamma, write_elias_delta, write_golomb, read_elias_gamma, read_elias_delta, read_golomb, 
write_fibonacci, read_fibonacci, write_zeta, read_zeta

using Adjacently.Distribution: get_degree_entropy, get_entropy, powerlaw_sample

using Adjacently.Constants: GOLOMB_BASE, ZETA_BASE

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
    bitvector_to_bitreader(bits::BitVector)

Create a BitReader from a BitVector for testing purposes.
"""
function bitvector_to_bitreader(bits::BitVector)
    io = IOBuffer()
    writer = BitWriter(io)
    
    # write each bit from the BitVector
    for bit in bits
        write_bit(writer, bit)
    end
    
    # flush and prepare for reading
    flush_bitwriter(writer; flush_last_bits=true)
    seekstart(io)
    
    return BitReader(io)
end