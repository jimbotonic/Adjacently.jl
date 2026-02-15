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

module MGS

using LightGraphs, DataStructures, Logging
using ..CustomTypes: UInt24, UInt40
using ..NodeTypes: Node, EmptyNode
using ..CustomLightGraphs: SimpleDiGraph, SimpleGraph, SimpleEdge
using ..Util: infer_uint_custom_type, to_bytes
using ..IO: BitWriter, write_bytes, flush_bitwriter, BitReader

using ..Compression: huffman_encoding, get_huffman_codes!, decode_huffman_values,
	delta_encode_vector, write_elias_coding, read_elias_coding,
	write_golomb, read_golomb, write_fibonacci, read_fibonacci,
	write_zeta, read_zeta, write_compressed_graph_data, read_compressed_graph_data,
	write_rl_compressed_graph_data, read_rl_compressed_graph_data

using ..Graph: get_basic_stats, get_in_out_degrees, get_out_degrees
using ..Relabeling: relabel_vertices, relabel_graph
using ..Constants: GOLOMB_BASE, ZETA_BASE

using ..RL: QPolicy, load_policy, best_action, action_from_index, feature_index,
	extract_features, VertexFeatures, Action, NUM_ACTIONS

# constants for new header format
# Header structure: 'MGS' (3 bytes) + major/minor version (2 bytes) + flags (2 bytes) + vertices (5 bytes)
# 
# New flag byte structure:
# Byte 1: graph_type (2 bits) | coding_scheme (2 bits) | integer_encoding (4 bits)  
# Byte 2: option_flags (8 bits)
#
# Graph type constants (2 bits)
const GRAPH_TYPE_DIRECTED = 0b00
const GRAPH_TYPE_UNDIRECTED = 0b01
#
# Coding scheme constants   (2 bits)
const CODING_SCHEME_CHILDREN = 0b00
const CODING_SCHEME_INDEX = 0b01
#
# Varint encoding constants (4 bits, values 1-6)
const INT_ENCODING_ELIAS_GAMMA = 0x1
const INT_ENCODING_ELIAS_DELTA = 0x2
const INT_ENCODING_GOLOMB = 0x3
const INT_ENCODING_FED = 0x4
const INT_ENCODING_ZETA = 0x5
const INT_ENCODING_FIBONACCI = 0x6
#
# Option flag constants (8 bits)
const OPTION_NONE = 0b00000000          # no options
const OPTION_ASTRA = 0b00001111         # ASTRA: Adaptive Streaming Adjacency (greedy cost-based, adaptive bitmaps, recursive references)
const OPTION_HUFFMAN = 0xFF             # Huffman compression (deprecated)
# RL policy range: 0x10-0x8F (128 policy slots)
const OPTION_RL_POLICY_BASE = 0x10
const OPTION_RL_POLICY_MAX = 0x8F

# maximum number of vertices
MGS_MAX_SIZE = 0xffffffffff

# Helper functions for new header format
function create_header_flags(graph_type::Symbol, coding_scheme::Symbol, integer_encoding::Symbol, option_flags::UInt8)
    """Create header flags bytes according to new format."""
    
    # Map symbols to constants
    gt = graph_type == :directed ? GRAPH_TYPE_DIRECTED : GRAPH_TYPE_UNDIRECTED
    cs = coding_scheme == :children ? CODING_SCHEME_CHILDREN : CODING_SCHEME_INDEX
    
    ie = if integer_encoding == :elias_gamma
        INT_ENCODING_ELIAS_GAMMA
    elseif integer_encoding == :elias_delta
        INT_ENCODING_ELIAS_DELTA
    elseif integer_encoding == :golomb
        INT_ENCODING_GOLOMB
    elseif integer_encoding == :fed
        INT_ENCODING_FED
    elseif integer_encoding == :zeta
        INT_ENCODING_ZETA
    elseif integer_encoding == :fibonacci
        INT_ENCODING_FIBONACCI
    else
        error("Unsupported integer encoding: $integer_encoding")
    end
    
    # Construct byte 1: graph_type (2 bits) | coding_scheme (2 bits) | integer_encoding (4 bits)
    byte1 = UInt8((gt << 6) | (cs << 4) | ie)
    
    # Byte 2 is the option flags
    byte2 = option_flags
    
    return (byte1, byte2)
end

function decode_header_flags(byte1::UInt8, byte2::UInt8)
    """Decode header flags bytes according to new format."""
    
    # Extract fields from byte 1
    graph_type_bits = (byte1 >> 6) & 0b11
    coding_scheme_bits = (byte1 >> 4) & 0b11
    integer_encoding_bits = byte1 & 0b1111
    
    # Map bits back to symbols
    graph_type = graph_type_bits == GRAPH_TYPE_DIRECTED ? :directed : :undirected
    coding_scheme = coding_scheme_bits == CODING_SCHEME_CHILDREN ? :children : :index
    
    integer_encoding = if integer_encoding_bits == INT_ENCODING_ELIAS_GAMMA
        :elias_gamma
    elseif integer_encoding_bits == INT_ENCODING_ELIAS_DELTA
        :elias_delta
    elseif integer_encoding_bits == INT_ENCODING_GOLOMB
        :golomb
    elseif integer_encoding_bits == INT_ENCODING_FED
        :fed
    elseif integer_encoding_bits == INT_ENCODING_ZETA
        :zeta
    elseif integer_encoding_bits == INT_ENCODING_FIBONACCI
        :fibonacci
    else
        error("Unknown integer encoding: $integer_encoding_bits")
    end
    
    # Decode option flags
    option_flags = byte2
    
    return (graph_type, coding_scheme, integer_encoding, option_flags)
end

function get_option_flags(use_mix_mode::Bool, reference_enabled::Bool, recursive_reference::Bool, huffman::Bool=false)
    """Get option flags based on compression features."""

    if huffman
        return OPTION_HUFFMAN  # deprecated
    elseif recursive_reference && reference_enabled && use_mix_mode
        return OPTION_ASTRA  # Full ASTRA encoding: adaptive + streaming + greedy references
    else
        return OPTION_NONE  # Basic delta encoding
    end
end

# Export the functions we want to make available
export write_mgs3_graph,
       write_compressed_mgs3_graph,
       load_mgs3_graph,
       load_compressed_mgs3_graph,
       write_huffman_compressed_mgs3_graph,
       write_elias_compressed_mgs3_graph,
       write_golomb_compressed_mgs3_graph,
	   write_fibonacci_compressed_mgs3_graph,
	   write_zeta_compressed_mgs3_graph,
	   write_complex_encoded_compressed_mgs3_graph,
	   load_huffman_compressed_mgs3_graph,
	   load_complex_encoded_compressed_mgs3_graph,
	   write_rl_compressed_mgs3_graph,
	   load_rl_compressed_mgs3_graph

################################################################################
# Write uncompressed MGS v3 graph
################################################################################

"""
    write_mgs3_graph(g::AbstractGraph{T},filename::AbstractString) where {T<:Unsigned}

write graph in format MGS v3

Supported encoding schemes:
- :children
- :index
"""
function write_mgs3_graph(g::AbstractGraph{T}, filename::AbstractString, encoding::Symbol=:children) where {T<:Unsigned}
  	# Header 12 bytes: 
	# -> version: 'MGS' 3 bytes string
	# -> major + minor version: 2 bytes
	# -> flags: 2 bytes
	#	 * Byte 1 (2 bits + 6 bits):
	#    	- graph type (2 bits): 			0x0: directed graph | 0x1: undirected graph
	#	 	- compression scheme (6 bits): 	0x0: no compression
	#	 * Byte 2:
	#		- coding scheme: 		0x0: data section only | 0x1: index+data section with implicit numbering of vertices
	#	 	- reserved flags: 		0x0: reserved
	# -> # vertices: 5 bytes position 
	#
	# <'MGS' string 3 bytes> + <16 bits major|minor version> + <flags 2 bytes> + <# vertices 5 bytes>
	vs = vertices(g)
	# number of vertices
	gs = convert(UInt64, length(vs))

	# if the graph has more than 2^40-1 vertices, `T` should be `UInt64`
	if gs > MGS_MAX_SIZE
		error("Input graph cannot have more than 2^40-1 vertices")
	end
	
	# `n_size_t` is the number of bytes needed to represent the graph vertices
	# NB: `n_size_t` (the computed size) may be lower than `sizeof(T)` of the graph type in parameter
	n_size_t = convert(UInt8, ceil(log(2, gs)/8))
	# NB: `T` should be an unsigned integer of size 1,2,4,8 bytes
	p_size_t = sizeof(T)
	
	#  The difference of size should >= 0 
	diff_size = p_size_t - n_size_t

	# Create header using new format (uncompressed graph = no options)
	option_flags = OPTION_NONE  # No compression options for uncompressed graph
	
	# We don't have integer encoding for uncompressed graphs, so use elias_delta as default
	flag_byte1, flag_byte2 = create_header_flags(:directed, encoding, :elias_delta, option_flags)
	
	# Construct 12-byte header
	header_bytes = UInt8[
		# 'MGS' signature (3 bytes)
		0x4d, 0x47, 0x53,
		# Major version = 3, Minor version = 0 (2 bytes)
		0x03, 0x00,
		# Flag bytes (2 bytes)
		flag_byte1, flag_byte2,
		# Number of vertices (5 bytes, little-endian UInt40)
		(gs & 0xff), ((gs >> 8) & 0xff), ((gs >> 16) & 0xff), ((gs >> 24) & 0xff), ((gs >> 32) & 0xff)
	]

	# create the output file (with extension .mgs)
	f = open(filename * ".mgs", "w")
	
	### write header (12 bytes total with new format)
	write(f, header_bytes)
	
	# coding scheme: data section only with implicit numbering of vertices
	if encoding == :children
		# array of 0x00 of size `n_size_t`
		stop_seq = [0x00 for i in 1:n_size_t]
		for v in vs
			ovs = outneighbors(g, v)
			for c in ovs
				# `c` is of type `T` (e.g. UInt8, UInt16, UInt32, UInt64)
				# `diff_size` is the difference in size between `T` and `n_size_t`
				# `n_size_t` is the number of bytes needed to represent the graph vertices
				# `p_size_t` is the size of the graph type in parameter
				#bytes = reverse(reinterpret(UInt8, [c]))[(diff_size+1):end]
				bytes = to_bytes(c)[(diff_size+1):end]
				write(f, bytes)
			end
			write(f, stop_seq)
		end
	# coding scheme: index+data sections with implicit numbering of vertices
	elseif encoding == :index
		# number of children for each vertex
		# NB: `T` should have a sufficient size to store the number of children
		ods = T[]
		# flattened list of children for all the vertices
		# NB: `T` should have a sufficient size to store the children indices
		children = T[]

		for v in vs
			ovs = outneighbors(g, v)
			push!(ods, convert(T, length(ovs)))
			for o in ovs
				push!(children, o)	
			end
		end

		### write index section
		for o in ods
			# `o` is of type `T` (e.g. UInt8, UInt16, UInt32, UInt64)
			# `diff_size` is the difference in size between `T` and `n_size_t`
			# `n_size_t` is the number of bytes needed to represent the graph vertices
			# `p_size_t` is the size of the graph type in parameter
			#bytes = reverse(reinterpret(UInt8, [o]))[(diff_size+1):end]
			bytes = to_bytes(o)[(diff_size+1):end]
			write(f, bytes)
		end
		### write data section
		for c in children
			# `c` is of type `T` (e.g. UInt8, UInt16, UInt32, UInt64)
			# `diff_size` is the difference in size between `T` and `n_size_t`
			# `n_size_t` is the number of bytes needed to represent the graph vertices
			# `p_size_t` is the size of the graph type in parameter
			#bytes = reverse(reinterpret(UInt8, [c]))[(diff_size+1):end]
			bytes = to_bytes(c)[(diff_size+1):end]
			write(f, bytes)
		end
	end
	close(f)
end

################################################################################
# Load uncompressed MGS v3 graph
################################################################################

""" 
    load_mgs3_graph(filename::AbstractString)

load graph in format MGS v3
"""
function load_mgs3_graph(filename::AbstractString)
	f = open(filename, "r")
	### read header
  	# 7 bytes: <3 bytes string 'MGS'> + <2 bytes major/minor> + <2 bytes flags>
	# 5 bytes (40 bits): number of vertices
	###
	# 7-bytes version
	# e.g. UInt8[0x4d, 0x47, 0x53, 0x03, 0x00, 0x00, 0x00]
	version = read(f,7)
	# number of vertices
	# NB: reinterpret generates an array of length 1 with type reinterpret(UInt64, ::Vector{UInt8})
	# NB: the 40 bits are stored in the last 5 bytes of the header in big endian format
	# NB: reinterpret function reads a byte array in little endian format
	gs = reinterpret(UInt64, vcat(reverse(read(f,5)),[0x00,0x00,0x00]))[1]
	# `n_size_v` is the number of bits needed to represent the graph vertices
	n_bits_v = convert(UInt8, ceil(log(2, gs)))
	# Get appropriate custom UInt type based on number of bits needed
	V = infer_uint_custom_type(n_bits_v)

	# graph type is 2 first bits of 6th byte of header	
	# 0x0: directed graph | 0x1: undirected graph
	graph_type = version[6] >> 6
	# coding scheme is 4 first bits of 7th byte of header	
	encoding = version[7] >> 4 == 0x0 ? :children : :index

	# intialize graph g according to graph type
	g = if graph_type == 0x0
		SimpleDiGraph{V}()
	else
		SimpleGraph{V}()
	end
	
	# vertex set
	vs = range(1, stop=gs)

	# add vertices to graph
    add_vertices!(g, gs)

	# coding scheme: data section only with implicit numbering of vertices
	# NB: each list of children is terminated with 0
	if encoding == :children
		# read data
		children = V[]
		while !eof(f)
			c = read(f, sizeof(V))
			# NB: reinterpret function reads a byte array in little endian format
			# NB: reverse function reverses the byte array in little endian format
			push!(children, reinterpret(V, reverse(c))[1])
		end
		
		# add edges
		pos = 1
		n_children = length(children)

		for i in 1:length(vs)
			if pos <= n_children
				source = convert(V, i)
				while children[pos] != 0x00
					target = children[pos]
					add_edge!(g, source, target)
					pos += 1
				end
				# skip 0x00
				pos += 1
			else
				# if we reached the last child, we are done
				break
			end
		end
	# coding scheme: index+data sections with implicit numbering of vertices
	elseif encoding == :index
		# read index
		# NB: 
		# - position indices are 1-based
		# - each position indicates the index of the first child of a vertex
		ods = V[]
		for i in 1:gs
			p = read(f, sizeof(V))
			push!(ods, reinterpret(V, reverse(p))[1])
		end
		# read data
		children = V[]
		while !eof(f)
			c = read(f, sizeof(V))
			push!(children, reinterpret(V, reverse(c))[1])
		end

		# add edges
		current_pos = 1
		for v in vs
			source = convert(V, v)
			# if we reached the last parent vertex
			if v == length(vs)
				pos1 = current_pos
				pos2 = length(children)
			else
				pos1 = current_pos
				# position of the last child of vertex v
				pos2 = current_pos + ods[v] - 1
				current_pos = pos2 + 1
			end
			# add edges for each child
			for p in pos1:pos2
				target = children[p]
				add_edge!(g, source, target)
			end
		end
	end
	close(f)
	return g::AbstractGraph{V}
end

################################################################################
# Write compressed MGS v3 graph
################################################################################

"""
    write_compressed_mgs3_graph(g::AbstractGraph{T}, filename::AbstractString, compression::Symbol=:huffman, coding_scheme::UInt8=0x00) where {T<:Unsigned}

Write graph in MGS v3 format with specified encoding and compression scheme

Parameters:
- g: Input graph
- filename: Output filename
- compression: Compression scheme to use (default: :huffman)
- encoding: Coding scheme (:children for children section only, :index for index+children sections)

Supported encoding schemes:
- :children
- :index

Supported compression types:
- :huffman - Huffman coding (default)
- :complex - complex encoding (delta + mix + reference + recursive reference)

Supported integer encodings:
- :elias_gamma - Elias gamma coding
- :elias_delta - Elias delta coding
- :golomb - Golomb coding
- :fibonacci - Fibonacci coding
- :zeta - Zeta coding
- :fed - Fibonacci+Elias Delta hybrid coding

Supported complex encoding options:
- :delta - delta encoding only
- :mix - mix encoding (run-length + interval)
- :hybrid - hybrid encoding (delta + mix + reference only)
- :hybrid+ - hybrid+ encoding (delta + mix + recursive reference)

@returns nothing
"""
function write_compressed_mgs3_graph(g::AbstractGraph{T}, filename::AbstractString, coding_scheme::Symbol=:children, compression::Symbol=:huffman, integer_encoding::Symbol=:fibonacci, use_mix_mode::Bool=true, reference_enabled::Bool=true, recursive_reference::Bool=true, ref_window_size::Int=1024) where {T<:Unsigned}
	# supported compression
	supported_compressions = [:huffman, :complex]
	supported_complex_options = [:delta, :mix, :hybrid, :hybrid_plus]
	supported_integer_encodings = [:elias_gamma, :elias_delta, :golomb, :fibonacci, :zeta, :fed]

	if compression == :huffman
        write_huffman_compressed_mgs3_graph(g, filename, coding_scheme)
    elseif compression == :complex
        write_complex_encoded_compressed_mgs3_graph(g, filename, coding_scheme, integer_encoding, use_mix_mode, reference_enabled, recursive_reference, ref_window_size)
    else
        error("Unsupported compression scheme: $compression. Supported schemes are :huffman, :complex")
    end
end

################################################################################
# Huffman compressed MGS v3 graph
################################################################################

"""
    write_huffman_compressed_mgs3_graph(g::AbstractGraph{T}, filename::AbstractString, coding_scheme::UInt8=0x00) where {T<:Unsigned}

write graph in a compressed MGS v3 (Huffman compression scheme)

Parameters:
- g: Input graph
- filename: Output filename
- encoding: Coding scheme (:children for children section only, :index for index+children sections)

@returns nothing
"""
function write_huffman_compressed_mgs3_graph(g::AbstractGraph{T}, filename::AbstractString, coding_scheme::Symbol=:children) where {T<:Unsigned}
	# Header format: see MGS_HEADER.md
	vs = vertices(g)
	# number of vertices
	gs = convert(UInt64, length(vs))

	# if the graph has more than 2^40-1 vertices, `T` should be `UInt64`
	if gs > MGS_MAX_SIZE
		error("Input graph cannot have more than 2^40-1 vertices")
	end
	
	# `n_bits_v` is the number of bits needed to represent the graph vertices
	n_bits_v = convert(UInt8, ceil(log(2, gs)))
	# Get appropriate custom UInt type based on number of bits needed
	V = infer_uint_custom_type(n_bits_v)
	# `n_size_t` is the number of bytes needed to represent the graph vertices
	# NB: `n_size_t` (the computed size) may be lower than `sizeof(T)` of the graph type in parameter
	n_size_t = sizeof(V)
	# NB: `T` should be an unsigned integer of size 1,2,4,8 bytes
	p_size_t = sizeof(T)
	
	#  The difference of size should >= 0 
	diff_size = p_size_t - n_size_t

	# Create header using new format (Huffman uses special option flag)
	option_flags = OPTION_HUFFMAN  # Huffman compression flag
	
	# For Huffman, we still need an integer encoding (doesn't matter since Huffman overrides)
	flag_byte1, flag_byte2 = create_header_flags(:directed, coding_scheme, :elias_delta, option_flags)
	
	# Construct 12-byte header
	header_bytes = UInt8[
		# 'MGS' signature (3 bytes)
		0x4d, 0x47, 0x53,
		# Major version = 3, Minor version = 0 (2 bytes)
		0x03, 0x00,
		# Flag bytes (2 bytes)
		flag_byte1, flag_byte2,
		# Number of vertices (5 bytes, little-endian UInt40)
		(gs & 0xff), ((gs >> 8) & 0xff), ((gs >> 16) & 0xff), ((gs >> 24) & 0xff), ((gs >> 32) & 0xff)
	]

	# create the output file (with extension .mgz)
	f = open(filename * ".mgz", "w")

	# flattened list of children for the whole graph
	children = V[]

	# frequencies of each vertex (in- and out- degrees)
	in_degrees, out_degrees = get_in_out_degrees(g)

	# convert the out-degrees to the custom type `V`
	out_degrees = Dict{V,V}(k => convert(V, v) for (k, v) in out_degrees)
	# convert in-degrees to the custom type `V`
	in_degrees = Dict{V,V}(k => convert(V, v) for (k, v) in in_degrees)
	
	# collect all children and in-degrees
	for v in vs
		ovs = outneighbors(g, v)
		for o in ovs
			push!(children, convert(V, o))
		end
	end

	# if coding scheme is :children, we add a special stop sequence whose frequency is equal to `gs` the number of vertices
	# NB: vertex frequencies are in the range [0, gs-1]
	if coding_scheme == :children
		# add the stop sequence 0 to the frequencies
		# NB: 
		# - vertex IDs are in the range [1, gs]
		# - 0 is the stop sequence and its frequency is set to the total number of vertices
		in_degrees[zero(V)] = convert(V, gs)
	end

	@info("generating Huffman tree")
	# get Huffman encoding tree
	tree = huffman_encoding(in_degrees)

	@info("getting Huffman codes")
	# get Huffman codes in C (C:code -> value)
	C = Dict{BitArray{1}, V}()
	get_huffman_codes!(tree, C, BitArray{1}())

	# reverse dictionary (R: value -> code)
	R = Dict{V, BitArray{1}}()
	[R[value] = key for (key, value) in C]

	@info("writing header section (new format)")
	### write header (12 bytes total with new format)
	write(f, header_bytes)
	
	@info("writing frequency section")
	### write frequency section
	for v in vs
		bytes = to_bytes(in_degrees[v])[(diff_size+1):end]
		write(f, bytes)
	end

	if coding_scheme == :children
		# write stop sequence in frequency section
		bytes = to_bytes(in_degrees[zero(V)])[(diff_size+1):end]
		write(f, bytes)
		
		@info("writing data section with stop sequence")
		### write data section
		cdata = BitArray{1}()
		# stop sequence is the code associated to 0
		# NB: R has a length of `gs+1` because of the stop sequence
		stop_seq_code = R[zero(V)] 
		
		for v in vs
			ovs = outneighbors(g, v)
			for c in ovs
				code = R[convert(V, c)]
				append!(cdata, code)
			end
			append!(cdata, stop_seq_code)
		end
	elseif coding_scheme == :index
		@info("writing index section")
		### write index section
		for v in vs
			bytes = to_bytes(out_degrees[v])[(diff_size+1):end]
			write(f, bytes)
		end
		
		@info("writing data section")
		### write data section
		cdata = BitArray{1}()
		for c in children
			code = R[convert(V, c)]
			append!(cdata, code)
		end
	end
	
	# write the final padding
	# add padding if necessary to make the data section a multiple of 8 bytes
	scd = convert(UInt32, length(cdata))
	sp = 8 - scd%8
	for i in 1:sp
		push!(cdata,0)
	end
	
	# number of bytes to write
	nb = round(Int, length(cdata)/8)

	for i in 1:nb
		# BitArray chunks type is UInt64
		# we only need to keep the last byte of the chunks
		byte = reinterpret(UInt8, cdata[(i-1)*8+1:(i-1)*8+8].chunks)[1]
		write(f, byte)
	end

	# last byte indicates by how many bytes cdata was padded
	b = 0xff
	write(f, b >> sp)

	close(f)
end

################################################################################
# Complex-encoding (delta + run-length + interval + recursive reference) + variable length encoding of MGS v3 graph
################################################################################

"""
    write_complex_encoded_compressed_mgs3_graph(g::AbstractGraph{T}, filename::AbstractString, coding_scheme::Symbol=:children, compression::Symbol=:elias_delta) where {T<:Unsigned}

Write graph in a compressed MGS v3 format (Complex-encoding (delta + run-length + interval + recursive reference) + variable length encoding)

Parameters:
- g: Input graph
- filename: Output filename
- encoding: Coding scheme (:children for children section only, :index for index+children sections)
- compression: Compression scheme to use (default: :elias_delta)
- use_mix_mode: Whether to use mix mode (default: true)
- reference_enabled: Whether to enable reference encoding (default: true)
- recursive_reference: Whether to enable recursive reference encoding (default: true)

@returns nothing
"""
function write_complex_encoded_compressed_mgs3_graph(g::AbstractGraph{T}, filename::AbstractString, coding_scheme::Symbol=:children, integer_encoding::Symbol=:elias_delta, use_mix_mode::Bool=true, reference_enabled::Bool=true, recursive_reference::Bool=true, ref_window_size::Int=1024) where {T<:Unsigned}
	# Header format: see MGS_HEADER.md
	vs = vertices(g)
	# number of vertices
	gs = convert(UInt64, length(vs))

	# if the graph has more than 2^40-1 vertices, `T` should be `UInt64`
	if gs > MGS_MAX_SIZE
		error("Input graph cannot have more than 2^40-1 vertices")
	end

	# `n_bits_v` is the number of bits needed to represent the graph vertices
	n_bits_v = convert(UInt8, ceil(log(2, gs)))
	# Get appropriate custom UInt type based on number of bits needed
	V = infer_uint_custom_type(n_bits_v)

	# Create header using new format
	# Header: 'MGS' (3 bytes) + major/minor version (2 bytes) + flags (2 bytes) + vertices (5 bytes)
	
	# Get option flags based on compression features  
	option_flags = get_option_flags(use_mix_mode, reference_enabled, recursive_reference)
	
	# Create flag bytes using new format
	flag_byte1, flag_byte2 = create_header_flags(:directed, coding_scheme, integer_encoding, option_flags)
	
	# Construct 12-byte header  
	header_bytes = UInt8[
		# 'MGS' signature (3 bytes)
		0x4d, 0x47, 0x53,
		# Major version = 3, Minor version = 0 (2 bytes)
		0x03, 0x00,
		# Flag bytes (2 bytes) 
		flag_byte1, flag_byte2,
		# Number of vertices (5 bytes, little-endian UInt40)
		(gs & 0xff), ((gs >> 8) & 0xff), ((gs >> 16) & 0xff), ((gs >> 24) & 0xff), ((gs >> 32) & 0xff)
	]

	# create the output file (with extension .mgz)
	f = open(filename * ".mgz", "w")

	# create a bitwriter
	bw = BitWriter(f)

	@info("writing header section (new format)")
	### write header (12 bytes total)
	# Write all header bytes at once
	write_bytes(bw, header_bytes)

	# Build neighbor lists from graph for write_compressed_graph_data
	neighbor_lists = Dict{V,Vector{V}}()
	for v in vs
		ovs = outneighbors(g, v)
		neighbor_lists[convert(V, v)] = [convert(V, o) for o in ovs]
	end

	# Use the comprehensive write_compressed_graph_data function
	# This handles mix encoding (delta + run-length) with reference encoding
	@info("writing compressed graph data using write_compressed_graph_data")
	write_compressed_graph_data(bw, neighbor_lists, coding_scheme, integer_encoding, use_mix_mode, reference_enabled, recursive_reference, ref_window_size)

	# flush the bitwriter and close the file
	flush_bitwriter(bw; flush_last_bits=true)
	close(f)
end

################################################################################
# Load compressed MGS v3 graph
################################################################################

"""
    load_compressed_mgs3_graph(filename::AbstractString)

Load graph in MGS v3 format.

Parameters:
- filename: Input filename

Supported compression schemes:
- :huffman
- :elias_gamma
- :elias_delta
- :golomb
- :fibonacci
- :zeta
- :fed

Returns a graph loaded with the compression scheme specified in the header.
"""
function load_compressed_mgs3_graph(filename::AbstractString)
	# Header format: see MGS_HEADER.md
	f = open(filename, "r")
	
	# read header (12 bytes total with new format)
	header_bytes = read(f, 12)
	
	# Verify MGS signature 
	if header_bytes[1:3] != [0x4d, 0x47, 0x53]  # 'MGS'
		error("Invalid MGS file signature")
	end
	
	# Extract version
	major_version = header_bytes[4]
	minor_version = header_bytes[5]
	
	# Extract and decode flags (new format)
	flag_byte1 = header_bytes[6] 
	flag_byte2 = header_bytes[7]
	
	# Decode header flags using new format
	graph_type, encoding, compression, option_flags = decode_header_flags(flag_byte1, flag_byte2)
	
	# Extract number of vertices (5 bytes, little-endian)
	gs_bytes = header_bytes[8:12]
	gs = UInt64(gs_bytes[1]) | (UInt64(gs_bytes[2]) << 8) | (UInt64(gs_bytes[3]) << 16) | 
		 (UInt64(gs_bytes[4]) << 24) | (UInt64(gs_bytes[5]) << 32)

	# supported compression schemes
	supported_schemes = [:elias_gamma, :elias_delta, :golomb, :fed, :zeta, :fibonacci]
	
	# Check if Huffman is enabled via option flags
	if option_flags == OPTION_HUFFMAN
		g = load_huffman_compressed_mgs3_graph(f, graph_type, encoding, gs)
	elseif OPTION_RL_POLICY_BASE <= option_flags <= OPTION_RL_POLICY_MAX
		policy_id = Int(option_flags - OPTION_RL_POLICY_BASE) + 1
		g = load_rl_compressed_mgs3_graph(f, graph_type, encoding, gs, policy_id, compression)
	elseif compression in supported_schemes
		g = load_complex_encoded_compressed_mgs3_graph(f, graph_type, encoding, gs, compression)
	else
		error("Unsupported compression scheme: $compression. Supported integer encodings are: $(join(supported_schemes, ", "))")
    end

	close(f)
	return g
end

""" 
    load_huffman_compressed_mgs3_graph(f::IO, graph_type::Symbol, encoding::Symbol, gs::UInt64)

load graph in compressed MGS v3 format (Huffman compression scheme)

Parameters:
- f: Input stream
- graph_type: graph type (:directed or :undirected)
- encoding: coding scheme (:children or :index)
- gs: number of vertices
"""
function load_huffman_compressed_mgs3_graph(io::IO, graph_type::Symbol, encoding::Symbol, gs::UInt64)
	# `n_size_u` is the number of bits needed to represent the graph vertices
	n_bits_v = convert(UInt8, ceil(log(2, gs)))
	# Get appropriate unsigned int type based on number of bits needed
	V = infer_uint_custom_type(n_bits_v)
	
	# intialize graph g according to graph type
	g = graph_type == :directed ? SimpleDiGraph{V}() : SimpleGraph{V}()
	
	# vertex set
	vs = range(1, stop=gs)
	# NB: stop sequence is the code associated to 0
	stop_seq = zero(V)

	@info("generating graph")
	@info("adding vertices")
	# add vertices to graph
    add_vertices!(g, gs)
	
	# frequencies of each vertex (in-degree)
	in_degrees = Dict{V,V}()

	# read frequency section for Huffman decoding
	for v in vs
		p = read(io, sizeof(V))
		in_degrees[v] = reinterpret(V, reverse(p))[1]
	end

	if encoding == :children
		# read stop sequence from frequency section
		p = read(io, sizeof(V))
		stop_seq_value = reinterpret(V, reverse(p))[1]
		in_degrees[stop_seq] = stop_seq_value
	end

	out_degrees = Dict{V,V}()
	if encoding == :index
		@info("reading index section")
		# read index
		for v in vs
			p = read(io, sizeof(V))
			out_degrees[v] = reinterpret(V, reverse(p))[1]
		end
	end

	@info("reading data section")
	# read data section
	cdata = BitVector()
	while !eof(io)
		# read a byte
		b = read(io, 1)[1]
		# read each bit of the byte
		for j in 0:7
			if ((b >> j) & 0x01) == 1
				push!(cdata, true)
			else
				push!(cdata, false)
			end
		end
	end

	# get last byte number of 0s
	sp = 8 - sum(cdata[end-7:end])
	cdata = cdata[1:end-(7+sp)]
	
	@info("generating Huffman tree")
	tree = huffman_encoding(in_degrees)
	
	@info("decoding values")
	children = decode_huffman_values(tree, cdata)

	# number of children
	n_children = length(children)

	@info("adding edges")
	if encoding == :children
		pos = 1
		for v in vs
			source = convert(V, v)
			while pos <= n_children && children[pos] != stop_seq
				target = children[pos]
				add_edge!(g, source, target)
				pos += 1
			end
			# skip stop sequence
			pos += 1
		end
	elseif encoding == :index
		current_pos = 1
		for v in vs
			source = convert(V,v)
			# if we reached the last parent vertex
			if v == length(vs)
				pos1 = current_pos
				pos2 = n_children
			else
				pos1 = current_pos
				# position of the last child of vertex v
				pos2 = current_pos + out_degrees[v] - 1
				current_pos = pos2 + 1
			end
			# add edges for each child
			for p in pos1:pos2
				target = children[p]
				add_edge!(g, source, target)
			end
		end
	end
	close(io)
	return g
end

"""
    load_complex_encoded_compressed_mgs3_graph(io::IO, graph_type::Symbol, encoding::Symbol, gs::UInt64, compression::Symbol)

Load graph in compressed MGS v3 format (Complex-encoding (delta + run-length + interval + recursive reference) + variable length encoding)

Parameters:
- io: Input stream
- graph_type: Graph type (:directed or :undirected)
- coding_scheme: Coding scheme (:children or :index)
- gs: Number of vertices
- integer_encoding: Integer encoding (:elias_delta, :fibonacci, :zeta, :fed)
- use_mix_mode: Whether to use mix mode (default: true)
- reference_enabled: Whether to enable reference encoding (default: true)

@returns a graph loaded with the compression scheme specified in the header.
"""
function load_complex_encoded_compressed_mgs3_graph(io::IO, graph_type::Symbol, coding_scheme::Symbol, gs::UInt64, integer_encoding::Symbol)
	# `n_size_u` is the number of bits needed to represent the graph vertices
	n_bits_v = convert(UInt8, ceil(log(2, gs)))
	# Get appropriate unsigned int type based on number of bits needed
	V = infer_uint_custom_type(n_bits_v)
	
	# Initialize graph according to graph type
	g = graph_type == :directed ? SimpleDiGraph{V}() : SimpleGraph{V}()
	
	@info("generating graph")
	@info("adding vertices")
	# Add vertices to graph
	add_vertices!(g, gs)
	
	# Create BitReader from the IO stream
	reader = BitReader(io)
	
	@info("reading compressed graph data using read_compressed_graph_data")
	# Use the comprehensive read_compressed_graph_data function
	# This handles mix encoding (delta + run-length) with reference encoding
	neighbor_lists = read_compressed_graph_data(reader, V(gs), coding_scheme, integer_encoding, V)
	
	@info("building graph from neighbor lists")
	# Build graph from the decoded neighbor lists
	for (source_vertex, neighbors) in neighbor_lists
		for target_vertex in neighbors
			add_edge!(g, source_vertex, target_vertex)
		end
	end
	
	@info("graph construction completed: $(nv(g)) vertices, $(ne(g)) edges")
	return g
end

################################################################################
# RL policy-based compressed MGS v3 graph
################################################################################

"""
    write_rl_compressed_mgs3_graph(g, filename, policy_filepath, policy_id; coding_scheme, ref_window_size)

Write graph in compressed MGS v3 format using an RL-learned policy for per-vertex encoding decisions.

The policy determines per-vertex choices of integer encoding, reference mode, and min interval length.
The compressed stream is self-describing: each vertex's data includes a compact encoding tag so the
decoder doesn't need the policy.

Parameters:
- g: Input graph
- filename: Output filename (without .mgz extension)
- policy_filepath: Path to .qpolicy file
- policy_id: Policy ID (1-16) stored in header byte 2
- coding_scheme: Coding scheme (:children or :index)
- ref_window_size: Reference window size
"""
function write_rl_compressed_mgs3_graph(g::AbstractGraph{T}, filename::AbstractString,
		policy_filepath::Union{AbstractString,Nothing}=nothing, policy_id::Int=1;
		coding_scheme::Symbol=:children, ref_window_size::Int=7,
		integer_encoding::Symbol=:fibonacci,
		vertex_actions::Union{Dict,Nothing}=nothing) where {T<:Unsigned}

	if !(1 <= policy_id <= 16)
		error("Policy ID must be between 1 and 16, got: $policy_id")
	end

	# Load the trained policy (or use greedy mode if no policy / vertex_actions)
	policy = nothing
	if vertex_actions !== nothing
		@info("Using pre-computed vertex actions ($(length(vertex_actions)) vertices)")
	elseif policy_filepath !== nothing && isfile(policy_filepath)
		@info("Loading RL policy from $policy_filepath")
		policy = load_policy(policy_filepath)
		@info("Policy loaded: $(policy.num_states) states, $(policy.num_actions) actions")
	else
		@info("Using greedy per-vertex optimization (no policy)")
	end

	vs = vertices(g)
	gs = convert(UInt64, length(vs))

	if gs > MGS_MAX_SIZE
		error("Input graph cannot have more than 2^40-1 vertices")
	end

	n_bits_v = convert(UInt8, ceil(log(2, gs)))
	V = infer_uint_custom_type(n_bits_v)

	# Write integer encoding in header (matches stream encoding tag)
	option_flags = UInt8(OPTION_RL_POLICY_BASE + policy_id - 1)
	flag_byte1, flag_byte2 = create_header_flags(:directed, coding_scheme, integer_encoding, option_flags)

	# Construct 12-byte header
	header_bytes = UInt8[
		0x4d, 0x47, 0x53,  # 'MGS'
		0x03, 0x00,         # Version 3.0
		flag_byte1, flag_byte2,
		(gs & 0xff), ((gs >> 8) & 0xff), ((gs >> 16) & 0xff), ((gs >> 24) & 0xff), ((gs >> 32) & 0xff)
	]

	f = open(filename * ".mgz", "w")
	bw = BitWriter(f)

	@info("writing header section (RL policy mode, policy_id=$policy_id)")
	write_bytes(bw, header_bytes)

	# Build neighbor lists (convert to compact type V)
	neighbor_lists = Dict{V,Vector{V}}()
	for v in vs
		ovs = outneighbors(g, v)
		neighbor_lists[convert(V, v)] = [convert(V, o) for o in ovs]
	end

	# Convert vertex_actions keys to compact type V if provided
	va_compact = nothing
	if vertex_actions !== nothing
		va_compact = Dict{V,Int}()
		for (k, v) in vertex_actions
			va_compact[convert(V, k)] = v
		end
	end

	@info("writing RL-compressed graph data")
	write_rl_compressed_graph_data(bw, neighbor_lists, policy, coding_scheme, ref_window_size;
		integer_encoding=integer_encoding, vertex_actions=va_compact)

	flush_bitwriter(bw; flush_last_bits=true)
	close(f)
end

"""
    load_rl_compressed_mgs3_graph(io, graph_type, coding_scheme, gs, policy_id)

Load graph from RL-compressed MGS v3 format. The data stream is self-describing,
so the policy is not needed for decompression (policy_id is informational).
"""
function load_rl_compressed_mgs3_graph(io::IO, graph_type::Symbol, coding_scheme::Symbol, gs::UInt64, policy_id::Int,
		integer_encoding::Symbol=:fibonacci)
	n_bits_v = convert(UInt8, ceil(log(2, gs)))
	V = infer_uint_custom_type(n_bits_v)

	g = graph_type == :directed ? SimpleDiGraph{V}() : SimpleGraph{V}()

	@info("generating graph (RL policy mode, policy_id=$policy_id)")
	@info("adding vertices")
	add_vertices!(g, gs)

	reader = BitReader(io)

	@info("reading RL-compressed graph data")
	neighbor_lists = read_rl_compressed_graph_data(reader, V(gs), coding_scheme, V;
		integer_encoding=integer_encoding)

	@info("building graph from neighbor lists")
	for (source_vertex, neighbors) in neighbor_lists
		for target_vertex in neighbors
			add_edge!(g, source_vertex, target_vertex)
		end
	end

	@info("graph construction completed: $(nv(g)) vertices, $(ne(g)) edges")
	return g
end

end
