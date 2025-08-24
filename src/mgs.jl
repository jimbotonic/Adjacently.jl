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
	write_zeta, read_zeta, write_compressed_graph_data, read_compressed_graph_data

using ..Graph: get_basic_stats, get_in_out_degrees, get_out_degrees, relabel_vertices, relabel_graph
using ..Constants: GOLOMB_BASE, ZETA_BASE

# constants
# 'MGS' + 0x0300 (major=3, minor=0) 
# - D0 (directed graph + no compression)
# - CS0 (coding scheme 0) = 0x00 (data section only + reserved)
HEADER_MGS3_D0_CS0 = 0x4d475303000000
# 'MGS' + 0x0300 (major=3, minor=0) 
# - D0 (directed graph + no compression)
# - CS1 (coding scheme 1) = 0x10 (index and data sections + reserved)
HEADER_MGS3_D0_CS1 = 0x4d475303000010
# 'MGS' + 0x0300 (major=3, minor=0) 
# - D1 (directed graph 00 + Huffman compression 000001)
# - CS0 (coding scheme 0) = 0x00 (data section only + reserved)
HEADER_MGS3_DH_CS0 = 0x4d475303000100
# 'MGS' + 0x0300 (major=3, minor=0) 
# - D0 (directed graph 00 + Huffman compression 000001)
# - CS1 (coding scheme 1) = 0x10 (index and data sections + reserved)
HEADER_MGS3_DH_CS1 = 0x4d475303000110
# 'MGS' + 0x0300 (major=3, minor=0) 
# - D0 (directed graph 00 + Elias gamma compression 000002)
# - CS0 (coding scheme 0) = 0x00 (data section only + reserved)
HEADER_MGS3_DEG_CS0 = 0x4d475303000200
# 'MGS' + 0x0300 (major=3, minor=0) 
# - D0 (directed graph 00 + Elias gamma compression 000002)
# - CS1 (coding scheme 1) = 0x10 (index and data sections + reserved)
HEADER_MGS3_DEG_CS1 = 0x4d475303000210
# 'MGS' + 0x0300 (major=3, minor=0) 
# - D0 (directed graph 00 + Elias delta compression 000003)
# - CS0 (coding scheme 0) = 0x00 (data section only + reserved)
HEADER_MGS3_DED_CS0 = 0x4d475303000300
# 'MGS' + 0x0300 (major=3, minor=0) 
# - D0 (directed graph 00 + Elias delta compression 000003)
# - CS1 (coding scheme 1) = 0x10 (index and data sections + reserved)
HEADER_MGS3_DED_CS1 = 0x4d475303000310
# 'MGS' + 0x0300 (major=3, minor=0) 
# - D0 (directed graph 00 + Golomb compression 000004)
# - CS0 (coding scheme 0) = 0x00 (data section only + reserved)
HEADER_MGS3_DG_CS0 = 0x4d475303000400
# 'MGS' + 0x0300 (major=3, minor=0) 
# - D0 (directed graph 00 + Golomb compression 000004)
# - CS1 (coding scheme 1) = 0x10 (index and data sections + reserved)
HEADER_MGS3_DG_CS1 = 0x4d475303000410
# 'MGS' + 0x0300 (major=3, minor=0) 
# - D0 (directed graph 00 + Fibonacci compression 000005)
# - CS0 (coding scheme 0) = 0x00 (data section only + reserved)
HEADER_MGS3_DF_CS0 = 0x4d475303000500
# 'MGS' + 0x0300 (major=3, minor=0) 
# - D0 (directed graph 00 + Fibonacci compression 000005)
# - CS1 (coding scheme 1) = 0x10 (index and data sections + reserved)
HEADER_MGS3_DF_CS1 = 0x4d475303000510
# 'MGS' + 0x0300 (major=3, minor=0) 
# - D0 (directed graph 00 + Zeta compression 000006)
# - CS0 (coding scheme 0) = 0x00 (data section only + reserved)
HEADER_MGS3_DZ_CS0 = 0x4d475303000600
# 'MGS' + 0x0300 (major=3, minor=0) 
# - D0 (directed graph 00 + Zeta compression 000006)
# - CS1 (coding scheme 1) = 0x10 (index and data sections + reserved)
HEADER_MGS3_DZ_CS1 = 0x4d475303000610
# 'MGS' + 0x0300 (major=3, minor=0) 
# - D0 (directed graph 00 + Fibonacci+Elias Delta hybrid compression 000007)
# - CS0 (coding scheme 0) = 0x00 (data section only + reserved)
HEADER_MGS3_DFED_CS0 = 0x4d475303000700
# 'MGS' + 0x0300 (major=3, minor=0) 
# - D0 (directed graph 00 + Fibonacci+Elias Delta hybrid compression 000007)
# - CS1 (coding scheme 1) = 0x10 (index and data sections + reserved)
HEADER_MGS3_DFED_CS1 = 0x4d475303000710

# maximum number of vertices
MGS_MAX_SIZE = 0xffffffffff

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
	   write_mix_encoded_compressed_mgs3_graph,
	   load_huffman_compressed_mgs3_graph,
	   load_mix_encoded_compressed_mgs3_graph

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

	if encoding == :children
		# 'MGS' + 0x0300 + 0x00 (directed graph + no compression) + 0x00 (data section only + reserved)
		# encoding: data section only with implicit numbering of vertices
  		version = HEADER_MGS3_D0_CS0
	elseif encoding == :index
		# 'MGS' + 0x0300 + 0x00 (directed graph + no compression) + 0x10 (index and data sections + reserved)
		# encoding: index+data sections with implicit numbering of vertices
  		version = HEADER_MGS3_D0_CS1
	end

	# create the output file (with extension .mgs)
	f = open(filename * ".mgs", "w")
	
	### write header
	# MGS version + parameters (7 bytes)
	# NB: reinterpret generates an array of length 8 even if version has a length of 7 bytes
	bytes = reverse(reinterpret(UInt8, [version]))[2:8]
	write(f, bytes)

	# write the number of vertices (5 bytes)
	bytes = reverse(reinterpret(UInt8, [gs]))[4:8]
	write(f, bytes)
	
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

Supported compression schemes:
- :huffman - Huffman coding (default)
- :elias_gamma - Elias gamma coding
- :elias_delta - Elias delta coding
- :golomb - Golomb coding
- :fibonacci - Fibonacci coding
- :zeta - Zeta coding
- :fed - Fibonacci+Elias Delta hybrid coding

@returns nothing
"""
function write_compressed_mgs3_graph(g::AbstractGraph{T}, filename::AbstractString, encoding::Symbol=:children, compression::Symbol=:huffman, use_mix_mode::Bool=true, reference_enabled::Bool=true) where {T<:Unsigned}
	# supported compression schemes
	supported_schemes = [:huffman, :elias_gamma, :elias_delta, :golomb, :fibonacci, :zeta, :fed]
	
	if compression == :huffman
        write_huffman_compressed_mgs3_graph(g, filename, encoding)
    elseif compression in supported_schemes
        write_mix_encoded_compressed_mgs3_graph(g, filename, encoding, compression, use_mix_mode, reference_enabled)
    else
        error("Unsupported compression scheme: $compression. Supported schemes are :huffman, :elias_gamma, :elias_delta, :golomb, :fibonacci, :zeta, :fed")
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
	# Header 12 bytes: 
	# -> version: 'MGS' 3 bytes string
	# -> major + minor version: 2 bytes
	# -> flags: 2 bytes
	#	 * Byte 1 (2 bits + 6 bits):
	#    	- graph type (2 bits): 			0x0: directed graph | 0x1: undirected graph
	#	 	- compression scheme (6 bits): 	0x1: Huffman | 0x2: Elias gamma | 0x3: Elias delta | 0x4: Golomb | 0x5: Fibonacci | 0x6: Zeta
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

	if coding_scheme == :children
		# 'MGS' + 0x0300 + 0x01 (directed graph + Huffman compression) + 0x00 (data section only + reserved)
		# encoding: data section only with implicit numbering of vertices
  		version = HEADER_MGS3_DH_CS0
	elseif coding_scheme == :index
		# 'MGS' + 0x0300 + 0x01 (directed graph + Huffman compression) + 0x10 (index and data sections + reserved)
		# encoding: index+data sections with implicit numbering of vertices
  		version = HEADER_MGS3_DH_CS1
	end

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

	@info("writing header section")
	### write header
	# MGS version + parameters (7 bytes)
	# NB: reinterpret generates an array of length 8 even if version has a length of 7 bytes
	bytes = reverse(reinterpret(UInt8, [version]))[2:8]
	write(f, bytes)

	# write the number of vertices (5 bytes)
	bytes = reverse(reinterpret(UInt8, [gs]))[4:8]
	write(f, bytes)
	
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
# Mix-encoding (delta + run-length) + reference + variable length encoding of MGS v3 graph
################################################################################

"""
    write_mix_encoded_compressed_mgs3_graph(g::AbstractGraph{T}, filename::AbstractString, coding_scheme::Symbol=:children, compression::Symbol=:elias_delta) where {T<:Unsigned}

Write graph in a compressed MGS v3 format (Mix-encoding (delta + run-length) + reference + variable length encoding)

Parameters:
- g: Input graph
- filename: Output filename
- encoding: Coding scheme (:children for children section only, :index for index+children sections)

@returns nothing
"""
function write_mix_encoded_compressed_mgs3_graph(g::AbstractGraph{T}, filename::AbstractString, coding_scheme::Symbol=:children, compression::Symbol=:elias_delta, use_mix_mode::Bool=true, reference_enabled::Bool=true) where {T<:Unsigned}
	# Header 12 bytes: 
	# -> version: 'MGS' 3 bytes string
	# -> major + minor version: 2 bytes
	# -> flags: 2 bytes
	#	 * Byte 1 (2 bits + 6 bits):
	#    	- graph type (2 bits): 			0x0: directed graph | 0x1: undirected graph
	#	 	- compression scheme (6 bits): 	0x1: Huffman | 0x2: Elias gamma | 0x3: Elias delta | 0x4: Golomb | 0x5: Fibonacci | 0x6: Zeta
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

	# `n_bits_v` is the number of bits needed to represent the graph vertices
	n_bits_v = convert(UInt8, ceil(log(2, gs)))
	# Get appropriate custom UInt type based on number of bits needed
	V = infer_uint_custom_type(n_bits_v)

	if compression == :elias_gamma
		if coding_scheme == :children
			# 'MGS' + 0x0300 + 0x02 (directed graph + compression) + 0x00 (data section only + reserved)
			# encoding: data section only with implicit numbering of vertices
			version = HEADER_MGS3_DEG_CS0
		elseif coding_scheme == :index
			# 'MGS' + 0x0300 + 0x02 (directed graph + compression) + 0x10 (index and data sections + reserved)
			# encoding: index+data sections with implicit numbering of vertices
			version = HEADER_MGS3_DEG_CS1
		end
	elseif compression == :elias_delta
		if coding_scheme == :children
			# 'MGS' + 0x0300 + 0x03 (directed graph + compression) + 0x00 (data section only + reserved)
			# encoding: data section only with implicit numbering of vertices
			version = HEADER_MGS3_DED_CS0
		elseif coding_scheme == :index
			# 'MGS' + 0x0300 + 0x03 (directed graph + compression) + 0x10 (index and data sections + reserved)
			# encoding: index+data sections with implicit numbering of vertices
			version = HEADER_MGS3_DED_CS1
		end
	elseif compression == :golomb
		if coding_scheme == :children
			# 'MGS' + 0x0300 + 0x04 (directed graph + compression) + 0x00 (data section only + reserved)
			# encoding: data section only with implicit numbering of vertices
			version = HEADER_MGS3_DG_CS0
		elseif coding_scheme == :index
			# 'MGS' + 0x0300 + 0x04 (directed graph + compression) + 0x10 (index and data sections + reserved)
			# encoding: index+data sections with implicit numbering of vertices
			version = HEADER_MGS3_DG_CS1
		end
	elseif compression == :fibonacci
		if coding_scheme == :children
			# 'MGS' + 0x0300 + 0x05 (directed graph + compression) + 0x00 (data section only + reserved)
			# encoding: data section only with implicit numbering of vertices
			version = HEADER_MGS3_DF_CS0
		elseif coding_scheme == :index
			# 'MGS' + 0x0300 + 0x05 (directed graph + compression) + 0x10 (index and data sections + reserved)
			# encoding: index+data sections with implicit numbering of vertices
			version = HEADER_MGS3_DF_CS1
		end
	elseif compression == :zeta
		if coding_scheme == :children
			# 'MGS' + 0x0300 + 0x06 (directed graph + compression) + 0x00 (data section only + reserved)
			# encoding: data section only with implicit numbering of vertices
			version = HEADER_MGS3_DZ_CS0
		elseif coding_scheme == :index
			# 'MGS' + 0x0300 + 0x06 (directed graph + compression) + 0x10 (index and data sections + reserved)
			# encoding: index+data sections with implicit numbering of vertices
			version = HEADER_MGS3_DZ_CS1
		end
	elseif compression == :fed
		if coding_scheme == :children
			# 'MGS' + 0x0300 + 0x07 (directed graph + compression) + 0x00 (data section only + reserved)
			# encoding: data section only with implicit numbering of vertices
			version = HEADER_MGS3_DFED_CS0
		elseif coding_scheme == :index
			# 'MGS' + 0x0300 + 0x07 (directed graph + compression) + 0x10 (index and data sections + reserved)
			# encoding: index+data sections with implicit numbering of vertices
			version = HEADER_MGS3_DFED_CS1
		end
	else
		error("Compression scheme not supported")
	end

	# create the output file (with extension .mgz)
	f = open(filename * ".mgz", "w")

	# create a bitwriter
	bw = BitWriter(f)

	@info("writing header section")
	### write header
	# MGS version + parameters (7 bytes)
	# NB: reinterpret generates an array of length 8 even if version has a length of 7 bytes
	bytes = reverse(reinterpret(UInt8, [version]))[2:8]
	write_bytes(bw, bytes)

	# write the number of vertices (5 bytes)
	bytes = reverse(reinterpret(UInt8, [gs]))[4:8]
	write_bytes(bw, bytes)

	# Build neighbor lists from graph for write_compressed_graph_data
	neighbor_lists = Dict{V,Vector{V}}()
	for v in vs
		ovs = outneighbors(g, v)
		neighbor_lists[convert(V, v)] = [convert(V, o) for o in ovs]
	end

	# Use the comprehensive write_compressed_graph_data function
	# This handles mix encoding (delta + run-length) with reference encoding
	@info("writing compressed graph data using write_compressed_graph_data")
	write_compressed_graph_data(bw, neighbor_lists, compression, coding_scheme, use_mix_mode, reference_enabled)

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
	# Header 12 bytes: 
	# -> version: 'MGS' 3 bytes string
	# -> major + minor version: 2 bytes
	# -> flags: 2 bytes
	#	 * Byte 1 (2 bits + 6 bits):
	#    	- graph type (2 bits): 			0x0: directed graph | 0x1: undirected graph
	#	 	- compression scheme (6 bits): 	0x1: Huffman | 0x2: Elias gamma | 0x3: Elias delta | 0x4: Golomb | 0x5: Fibonacci | 0x6: Zeta
	#	 * Byte 2 (4 bits + 4 bits):
	#		- coding scheme: 		0x0: data section only | 0x1: index+data section with implicit numbering of vertices
	#	 	- reserved flags: 		0x0: reserved
	# -> # vertices: 5 bytes position 
	#
	# <'MGS' string 3 bytes> + <16 bits major|minor version> + <flags 2 bytes> + <# vertices 5 bytes>
	f = open(filename, "r")
	
	# read header
	version = read(f,5)
	# major version
	major_version = version[4]
	# minor version
	minor_version = version[5]

	# flags
	flags = read(f,2)
	# graph type
	graph_type = flags[1] >> 6 == 0x0 ? :directed : :undirected
	# compression scheme
	compression_scheme = flags[1] & 0x3F
	# coding scheme
	encoding = flags[2] >> 4 == 0x0 ? :children : :index

	if compression_scheme == 0x1
		compression = :huffman
	elseif compression_scheme == 0x2
		compression = :elias_gamma
	elseif compression_scheme == 0x3
		compression = :elias_delta
	elseif compression_scheme == 0x4
		compression = :golomb
	elseif compression_scheme == 0x5
		compression = :fibonacci
	elseif compression_scheme == 0x6
		compression = :zeta
	elseif compression_scheme == 0x7
		compression = :fed
	else
		error("Unsupported compression scheme: $compression_scheme. Supported schemes are :huffman, :elias_gamma, :elias_delta, :fibonacci, :zeta, :fed")
	end
	# number of vertices
	gs = reinterpret(UInt64, vcat(reverse(read(f,5)),[0x00,0x00,0x00]))[1]

	# supported compression schemes
	supported_schemes = [:huffman, :elias_gamma, :elias_delta, :fibonacci, :zeta, :fed]
	
	if compression == :huffman
		g = load_huffman_compressed_mgs3_graph(f, graph_type, encoding, gs)
	elseif compression in supported_schemes
		g = load_mix_encoded_compressed_mgs3_graph(f, graph_type, encoding, gs, compression)
	else
		error("Unsupported compression scheme: $compression. Supported schemes are :huffman, :elias_gamma, :elias_delta, :fibonacci, :zeta, :fed")
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
    load_mix_encoded_compressed_mgs3_graph(io::IO, graph_type::Symbol, encoding::Symbol, gs::UInt64, compression::Symbol)

Load graph in compressed MGS v3 format (Mix-encoding (delta + run-length) + reference + variable length encoding)

Parameters:
- io: Input stream
- graph_type: Graph type (:directed or :undirected)
- encoding: Coding scheme (:children or :index)
- gs: Number of vertices
- compression: Compression scheme (:elias_delta, :fibonacci, :zeta)

@returns a graph loaded with the compression scheme specified in the header.
"""
function load_mix_encoded_compressed_mgs3_graph(io::IO, graph_type::Symbol, encoding::Symbol, gs::UInt64, compression::Symbol)
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
	neighbor_lists = read_compressed_graph_data(reader, V(gs), compression, encoding, V)
	
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

end