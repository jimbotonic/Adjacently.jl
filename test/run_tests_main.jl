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

using Adjacently.Graph: get_core, get_reverse_graph, get_basic_stats, relabel_graph, relabel_vertices,
get_out_degrees, get_in_degrees, get_in_out_degrees

using Adjacently.MGS: write_mgs3_graph, write_compressed_mgs3_graph, load_mgs3_graph, load_compressed_mgs3_graph

using Adjacently.Util: bottom_up_sort, quicksort_iterative_permutation!, get_sorted_array, binary_search,
infer_uint_custom_type

using Adjacently.Compression: write_unary_coding, write_truncated_binary_coding, huffman_encoding, 
encode_huffman_tree!, decode_huffman_tree!, get_huffman_codes!, 
write_elias_gamma, write_elias_delta, write_golomb, read_elias_gamma, read_elias_delta, read_golomb, 
write_fibonacci, read_fibonacci, write_zeta, read_zeta, delta_encode_vector,
write_encoded_value, read_encoded_value, write_run_length_delta, read_run_length_delta

using Adjacently.Distribution: get_graph_entropy, get_degree_entropy, get_entropy, powerlaw_sample

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

"""
    samples_bits(samples::Vector{T}, encoding::Symbol) where {T<:Unsigned}

    Parameters:
    - samples: Vector of samples to encode
    - encoding: Encoding scheme (:elias_gamma, :elias_delta, :golomb, :fibonacci, :zeta)

Return encoded samples as a BitVector for a given samples.
"""
function samples_bits(samples::Vector{T}, encoding::Symbol) where {T<:Unsigned}
    # Use a simple approach: encode each sample individually and collect the bits
    all_bits = Bool[]
    
    for sample in samples
        io = IOBuffer()
        writer = BitWriter(io)
        
        if encoding == :elias_gamma
            write_elias_gamma(writer, sample)
        elseif encoding == :elias_delta
            write_elias_delta(writer, sample)
        elseif encoding == :golomb
            write_golomb(writer, sample, GOLOMB_BASE)
        elseif encoding == :fibonacci
            write_fibonacci(writer, sample)
        elseif encoding == :zeta
            write_zeta(writer, sample, ZETA_BASE)
        else
            error("Invalid encoding scheme: $encoding")
        end
        
        # capture number of bits for this sample
        sample_nbits = writer.index - 1
        flush_bitwriter(writer; flush_last_bits=true)
        
        # read the bits for this sample
        seekstart(io)
        reader = BitReader(io)
        sample_bits = read_bits(reader, sample_nbits)
        
        # append to all_bits
        append!(all_bits, sample_bits)
    end
    
    return all_bits
end