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
	write_golomb, read_golomb, write_fibonacci_code, read_fibonacci_code,
	write_zeta_coding, read_zeta_coding
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
	   load_huffman_compressed_mgs3_graph,
	   load_elias_compressed_mgs3_graph,
	   load_golomb_compressed_mgs3_graph,
	   load_fibonacci_compressed_mgs3_graph,
	   load_zeta_compressed_mgs3_graph

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

@returns nothing
"""
function write_compressed_mgs3_graph(g::AbstractGraph{T}, filename::AbstractString, encoding::Symbol=:children, compression::Symbol=:huffman) where {T<:Unsigned}
    if compression == :huffman
        write_huffman_compressed_mgs3_graph(g, filename, encoding)
    elseif compression == :elias_gamma || compression == :elias_delta
        write_elias_compressed_mgs3_graph(g, filename, encoding, compression)
    elseif compression == :golomb
        write_golomb_compressed_mgs3_graph(g, filename, encoding)
    elseif compression == :fibonacci
        write_fibonacci_compressed_mgs3_graph(g, filename, encoding)
	elseif compression == :zeta
		write_zeta_compressed_mgs3_graph(g, filename, encoding)
    else
        error("Unsupported compression scheme: $compression. Supported schemes are :huffman, :elias, :golomb")
    end
end

"""
    write_huffman_compressed_mgs3_graph(g::AbstractGraph{T}, filename::AbstractString, coding_scheme::UInt8=0x00) where {T<:Unsigned}

write graph in a compressed MGS v3 (Huffman compression scheme)

Parameters:
- g: Input graph
- filename: Output filename
- encoding: Coding scheme (:children for children section only, :index for index+children sections)

@returns nothing
"""
function write_huffman_compressed_mgs3_graph(g::AbstractGraph{T}, filename::AbstractString, encoding::Symbol=:children) where {T<:Unsigned}
	# Header 12 bytes: 
	# -> version: 'MGS' 3 bytes string
	# -> major + minor version: 2 bytes
	# -> flags: 2 bytes
	#	 * Byte 1 (2 bits + 6 bits):
	#    	- graph type (2 bits): 			0x0: directed graph | 0x1: undirected graph
	#	 	- compression scheme (6 bits): 	0x1: Huffman | 0x2: Elias gamma | 0x3: Elias delta | 0x4: Golomb
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

	if encoding == :children
		# 'MGS' + 0x0300 + 0x01 (directed graph + Huffman compression) + 0x00 (data section only + reserved)
		# encoding: data section only with implicit numbering of vertices
  		version = HEADER_MGS3_DH_CS0
	elseif encoding == :index
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
	# convert the degrees to the custom type `V`
	in_degrees = Dict{V,V}(k => convert(V, v) for (k, v) in in_degrees)
	out_degrees = Dict{V,V}(k => convert(V, v) for (k, v) in out_degrees)
	
	# collect all children and in-degrees
	for v in vs
		ovs = outneighbors(g, v)
		for o in ovs
			push!(children, convert(V, o))
		end
	end

	# if coding scheme is :children, we need a special stop sequence equal to the number of vertices
	# NB: frequencies are in the range [0, gs-1] and the stop sequence is the code associated to `gs`
	if encoding == :children
		# add the stop sequence typemax(V) to the frequencies
		# typemax(V) => gs
		in_degrees[typemax(V)] = convert(V, gs)
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

	if encoding == :children
		# write stop sequence in frequency section
		bytes = to_bytes(in_degrees[typemax(V)])[(diff_size+1):end]
		write(f, bytes)
		
		@info("writing data section with stop sequence")
		### write data section
		cdata = BitArray{1}()
		# stop sequence is the code associated to typemax(V)
		# NB: R has a length of `gs+1` because of the stop sequence
		stop_seq_code = R[typemax(V)] 
		
		for v in vs
			ovs = outneighbors(g, v)
			for c in ovs
				code = R[convert(V, c)]
				append!(cdata, code)
			end
			append!(cdata, stop_seq_code)
		end
	elseif encoding == :index
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

"""
    write_elias_compressed_mgs3_graph(g::AbstractGraph{T}, filename::AbstractString, encoding::Symbol=:children) where {T<:Unsigned}

Write graph in a compressed MGS v3 format (Elias gamma or delta compression scheme)

Parameters:
- g: Input graph
- filename: Output filename
- encoding: Coding scheme (:children for children section only, :index for index+children sections)
- compression: Compression scheme (:elias_gamma or :elias_delta)
"""
function write_elias_compressed_mgs3_graph(g::AbstractGraph{T}, filename::AbstractString, encoding::Symbol=:children, compression::Symbol=:elias_gamma) where {T<:Unsigned}
	# Header 12 bytes: 
	# -> version: 'MGS' 3 bytes string
	# -> major + minor version: 2 bytes
	# -> flags: 2 bytes
	#	 * Byte 1 (2 bits + 6 bits):
	#    	- graph type (2 bits): 			0x0: directed graph | 0x1: undirected graph
	#	 	- compression scheme (6 bits): 	0x1: Huffman | 0x2: Elias gamma | 0x3: Elias delta | 0x4: Golomb
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
		if encoding == :children
			# 'MGS' + 0x0300 + 0x02 (directed graph + compression) + 0x00 (data section only + reserved)
			# encoding: data section only with implicit numbering of vertices
			version = HEADER_MGS3_DEG_CS0
		elseif encoding == :index
			# 'MGS' + 0x0300 + 0x02 (directed graph + compression) + 0x10 (index and data sections + reserved)
			# encoding: index+data sections with implicit numbering of vertices
			version = HEADER_MGS3_DEG_CS1
		end
	elseif compression == :elias_delta
		if encoding == :children
			# 'MGS' + 0x0300 + 0x03 (directed graph + compression) + 0x00 (data section only + reserved)
			# encoding: data section only with implicit numbering of vertices
			version = HEADER_MGS3_DED_CS0
		elseif encoding == :index
			# 'MGS' + 0x0300 + 0x03 (directed graph + compression) + 0x10 (index and data sections + reserved)
			# encoding: index+data sections with implicit numbering of vertices
			version = HEADER_MGS3_DED_CS1
		end
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

	if encoding == :children
		@info("writing data section with stop sequence")
		# reserve 1 for the stop sequence
		stop_seq = one(V)

		### write data section
		for v in vs
			# get the outneighbors of the vertex
			# sort the outneighbors in ascending order
			# NB: as we do not deal with muti-graphs, all neighbors are unique
			# so we can sort them in ascending order
			ovs = sort(collect(V, outneighbors(g, v)))
			# delta encode the neighbors
			# NB: the starting value is the first element of the `diffs` vector
			# NB: as we do not deal with muti-graphs, all neighbors are unique
			# and no delta is equal to 0
			diffs = delta_encode_vector(ovs)
			# shift the diffs as 1 is reserved for the stop sequence
			# NB: all values in the original vector are greater than 0
			diffs .+= 1
			# write the diffs using Elias coding
			# NB: the starting value is the first outneighbor of the vertex
			for d in diffs
				write_elias_coding(bw, d, compression)
			end
			# if we did not reach the last parent vertex, write the stop sequence
			if v < gs
				# write the stop sequence
				write_elias_coding(bw, stop_seq, compression)
			end
		end
	elseif encoding == :index
		# frequencies of each vertex (out- degrees)
		out_degrees = get_out_degrees(g)
		# convert the out-degrees to the custom type `V` and shift as 0 is a forbidden value
		out_degrees = Dict{V,V}(k => convert(V, v) + 1 for (k, v) in out_degrees)
		
		@info("writing index section")
		### write index section
		for v in vs
			write_elias_coding(bw, out_degrees[v], compression)
		end
		@info("writing data section")
		### write data section
		for v in vs
			ovs = outneighbors(g, v)
			if !isempty(ovs)
				# sort the outneighbors in ascending order
				ovs = sort(collect(V, ovs))
				# get the delta encoding of the outneighbors
				# NB: all values in the original vector are greater than 0
				diffs = delta_encode_vector(ovs)
				# write the diffs using Elias coding
				# NB: the starting value is the first outneighbor of the vertex
				for d in diffs
					write_elias_coding(bw, d, compression)
				end
			end
		end
	end

	# flush the bitwriter and close the file
	flush_bitwriter(bw; flush_last_bits=true)
	close(f)
end

"""
    write_golomb_compressed_mgs3_graph(g::AbstractGraph{T}, filename::AbstractString, encoding::Symbol=:children) where {T<:Unsigned}

Write graph in a compressed MGS v3 format (Golomb compression scheme)

Parameters:
- g: Input graph
- filename: Output filename
- encoding: Coding scheme (:children for children section only, :index for index+children sections)
"""
function write_golomb_compressed_mgs3_graph(g::AbstractGraph{T}, filename::AbstractString, encoding::Symbol=:children) where {T<:Unsigned}
	# Header 12 bytes: 
	# -> version: 'MGS' 3 bytes string
	# -> major + minor version: 2 bytes
	# -> flags: 2 bytes
	#	 * Byte 1 (2 bits + 6 bits):
	#    	- graph type (2 bits): 			0x0: directed graph | 0x1: undirected graph
	#	 	- compression scheme (6 bits): 	0x1: Huffman | 0x2: Elias gamma | 0x3: Elias delta | 0x4: Golomb
	#	 * Byte 2:
	#		- coding scheme: 		0x0: data section only | 0x1: index+data section with implicit numbering of vertices
	#	 	- reserved flags: 		0x0: reserved
	# -> # vertices: 5 bytes position 
	#
	# <'MGS' string 3 bytes> + <16 bits major|minor version> + <flags 2 bytes> + <# vertices 5 bytes>
	# number of vertices
	gs = convert(UInt64, nv(g))
	vs = vertices(g)

	# if the graph has more than 2^40-1 vertices, `T` should be `UInt64`
	if gs > MGS_MAX_SIZE
		error("Input graph cannot have more than 2^40-1 vertices")
	end
	
	# `n_bits_v` is the number of bits needed to represent the graph vertices
	n_bits_v = convert(UInt8, ceil(log(2, gs)))
	# Get appropriate custom UInt type based on number of bits needed
	V = infer_uint_custom_type(n_bits_v)
	
	if encoding == :children
		version = HEADER_MGS3_DG_CS0
	elseif encoding == :index
		version = HEADER_MGS3_DG_CS1
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

	if encoding == :children
		@info("writing data section with stop sequence")
		# stop sequence is equal to 0
		stop_seq = zero(V)

		### write data section
		for v in vs
			ovs = outneighbors(g, v)
			if !isempty(ovs)
				# sort the outneighbors in ascending order
				ovs = sort(collect(V, ovs))
				# delta encode the neighbors
				# NB: the starting value is the first element of the `diffs` vector
				# NB: as we do not deal with muti-graphs, all neighbors are unique
				# and no delta is equal to 0
				diffs = delta_encode_vector(ovs)
				# write the diffs using Golomb coding
				for d in diffs
					if d == 0
						error("Delta is equal to 0. This should not happen as all neighbors are unique.")
					end
					write_golomb(bw, d, GOLOMB_BASE)
				end
			end
			# if we did not reach the last parent vertex, write the stop sequence
			if v < gs
				# write the stop sequence
				write_golomb(bw, stop_seq, GOLOMB_BASE)
			end
		end
	elseif encoding == :index
		# frequencies of each vertex (out- degrees)
		out_degrees = get_out_degrees(g)
		out_degrees = Dict{V,V}(k => convert(V, v) for (k, v) in out_degrees)

		@info("writing index section")
		### write index section
		for v in vs
			write_golomb(bw, out_degrees[v], GOLOMB_BASE)
		end
		@info("writing data section")
		### write data section
		for v in vs
			ovs = outneighbors(g, v)
			if !isempty(ovs)
				ovs = sort(collect(V, ovs))
				# get the delta encoding of the outneighbors
				# NB: all values in the original vector are greater than 0
				diffs = delta_encode_vector(ovs)
				# write the diffs using Golomb coding
				# NB: the starting value is the first outneighbor of the vertex
				for d in diffs
					if d == 0
						error("Delta is equal to 0. This should not happen as all neighbors are unique.")
					end
					write_golomb(bw, d, GOLOMB_BASE)
				end
			end
		end
	end

	# flush the bitwriter and close the file
	flush_bitwriter(bw; flush_last_bits=true)
	close(f)
end

"""
    write_fibonacci_compressed_mgs3_graph(g::AbstractGraph{T}, filename::AbstractString, encoding::Symbol=:children) where {T<:Unsigned}

Write graph in a compressed MGS v3 format (Fibonacci compression scheme)

Parameters:
- g: Input graph
- filename: Output filename
- encoding: Coding scheme (:children for children section only, :index for index+children sections)
"""
function write_fibonacci_compressed_mgs3_graph(g::AbstractGraph{T}, filename::AbstractString, encoding::Symbol=:children) where {T<:Unsigned}
	# Header 12 bytes: 
	# -> version: 'MGS' 3 bytes string
	# -> major + minor version: 2 bytes
	# -> flags: 2 bytes
	#	 * Byte 1 (2 bits + 6 bits):
	#    	- graph type (2 bits): 			0x0: directed graph | 0x1: undirected graph
	#	 	- compression scheme (6 bits): 	0x1: Huffman | 0x2: Elias gamma | 0x3: Elias delta | 0x4: Golomb
	#	 * Byte 2:
	#		- coding scheme: 		0x0: data section only | 0x1: index+data section with implicit numbering of vertices
	#	 	- reserved flags: 		0x0: reserved
	# -> # vertices: 5 bytes position 
	#
	# <'MGS' string 3 bytes> + <16 bits major|minor version> + <flags 2 bytes> + <# vertices 5 bytes>
	# number of vertices
	gs = convert(UInt64, nv(g))
	vs = vertices(g)

	# if the graph has more than 2^40-1 vertices, `T` should be `UInt64`
	if gs > MGS_MAX_SIZE
		error("Input graph cannot have more than 2^40-1 vertices")
	end
	
	# `n_bits_v` is the number of bits needed to represent the graph vertices
	n_bits_v = convert(UInt8, ceil(log(2, gs)))
	# Get appropriate custom UInt type based on number of bits needed
	V = infer_uint_custom_type(n_bits_v)
	
	if encoding == :children
		version = HEADER_MGS3_DF_CS0
	elseif encoding == :index
		version = HEADER_MGS3_DF_CS1
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

	if encoding == :children
		@info("writing data section with stop sequence")
		# stop sequence is equal to 1
		stop_seq = one(V)

		### write data section
		for v in vs
			ovs = outneighbors(g, v)
			if !isempty(ovs)
				# sort the outneighbors in ascending order
				ovs = sort(collect(V, ovs))
				# delta encode the neighbors
				# NB: the starting value is the first element of the `diffs` vector
				# NB: as we do not deal with muti-graphs, all neighbors are unique
				# and no delta is equal to 0
				diffs = delta_encode_vector(ovs)
				# shift the diffs by 1
				diffs .+= 1
				# write the diffs using Fibonacci coding
				for d in diffs
					if d == 0
						error("Delta is equal to 0. This should not happen as all neighbors are unique.")
					end
					write_fibonacci_code(bw, d)
				end
			end
			# if we did not reach the last parent vertex, write the stop sequence
			if v < gs
				# write the stop sequence
				write_fibonacci_code(bw, stop_seq)
			end
		end
	elseif encoding == :index
		# frequencies of each vertex (out- degrees)
		out_degrees = get_out_degrees(g)
		out_degrees = Dict{V,V}(k => convert(V, v) for (k, v) in out_degrees)

		@info("writing index section")
		### write index section
		for v in vs
			write_fibonacci_code(bw, out_degrees[v])
		end
		@info("writing data section")
		### write data section
		for v in vs
			ovs = outneighbors(g, v)
			if !isempty(ovs)
				ovs = sort(collect(V, ovs))
				# get the delta encoding of the outneighbors
				# NB: all values in the original vector are greater than 0
				diffs = delta_encode_vector(ovs)
				# shift the diffs by 1
				diffs .+= 1
				# write the diffs using Fibonacci coding
				# NB: the starting value is the first outneighbor of the vertex
				for d in diffs
					if d == 0
						error("Delta is equal to 0. This should not happen as all neighbors are unique.")
					end
					write_fibonacci_code(bw, d)
				end
			end
		end
	end

	# flush the bitwriter and close the file
	flush_bitwriter(bw; flush_last_bits=true)
	close(f)
end

"""
    write_zeta_compressed_mgs3_graph(g::AbstractGraph{T}, filename::AbstractString, encoding::Symbol=:children) where {T<:Unsigned}

Write graph in a compressed MGS v3 format (Zeta compression scheme)

Parameters:
- g: Input graph
- filename: Output filename
- encoding: Coding scheme (:children for children section only, :index for index+children sections)
"""
function write_zeta_compressed_mgs3_graph(g::AbstractGraph{T}, filename::AbstractString, encoding::Symbol=:children, k::Int=ZETA_BASE) where {T<:Unsigned}
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
	# number of vertices
	gs = convert(UInt64, nv(g))
	vs = vertices(g)

	# if the graph has more than 2^40-1 vertices, `T` should be `UInt64`
	if gs > MGS_MAX_SIZE
		error("Input graph cannot have more than 2^40-1 vertices")
	end
	
	# `n_bits_v` is the number of bits needed to represent the graph vertices
	n_bits_v = convert(UInt8, ceil(log(2, gs)))
	# Get appropriate custom UInt type based on number of bits needed
	V = infer_uint_custom_type(n_bits_v)
	
	if encoding == :children
		version = HEADER_MGS3_DZ_CS0
	elseif encoding == :index
		version = HEADER_MGS3_DZ_CS1
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

	if encoding == :children
		@info("writing data section with stop sequence")
		# stop sequence is equal to 1
		stop_seq = one(V)

		### write data section
		for v in vs
			ovs = outneighbors(g, v)
			if !isempty(ovs)
				# sort the outneighbors in ascending order
				ovs = sort(collect(V, ovs))
				# delta encode the neighbors
				# NB: the starting value is the first element of the `diffs` vector
				# NB: as we do not deal with muti-graphs, all neighbors are unique
				# and no delta is equal to 0
				diffs = delta_encode_vector(ovs)
				# shift the diffs by 1
				diffs .+= 1
				# write the diffs using Zeta coding
				for d in diffs
					if d == 0
						error("Delta is equal to 0. This should not happen as all neighbors are unique.")
					end
					write_zeta_coding(bw, d, k)
				end
			end
			# if we did not reach the last parent vertex, write the stop sequence
			if v < gs
				# write the stop sequence
				write_zeta_coding(bw, stop_seq, k)
			end
		end
	elseif encoding == :index
		# frequencies of each vertex (out- degrees)
		out_degrees = get_out_degrees(g)
		out_degrees = Dict{V,V}(k => convert(V, v) for (k, v) in out_degrees)

		@info("writing index section")
		### write index section
		for v in vs
			write_zeta_coding(bw, out_degrees[v], k)
		end
		@info("writing data section")
		### write data section
		for v in vs
			ovs = outneighbors(g, v)
			if !isempty(ovs)
				ovs = sort(collect(V, ovs))
				# get the delta encoding of the outneighbors
				# NB: all values in the original vector are greater than 0
				diffs = delta_encode_vector(ovs)
				# shift the diffs by 1
				diffs .+= 1
				# write the diffs using Fibonacci coding
				# NB: the starting value is the first outneighbor of the vertex
				for d in diffs
					if d == 0
						error("Delta is equal to 0. This should not happen as all neighbors are unique.")
					end
					write_zeta_coding(bw, d, k)
				end
			end
		end
	end

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

Returns a graph loaded with the compression scheme specified in the header.
"""
function load_compressed_mgs3_graph(filename::AbstractString)
	# Header 12 bytes: 
	# -> version: 'MGS' 3 bytes string
	# -> major + minor version: 2 bytes
	# -> flags: 2 bytes
	#	 * Byte 1 (2 bits + 6 bits):
	#    	- graph type (2 bits): 			0x0: directed graph | 0x1: undirected graph
	#	 	- compression scheme (6 bits): 	0x1: Huffman | 0x2: Elias gamma | 0x3: Elias delta | 0x4: Golomb
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

	# number of vertices
	gs = reinterpret(UInt64, vcat(reverse(read(f,5)),[0x00,0x00,0x00]))[1]
	
	# `n_size_u` is the number of bits needed to represent the graph vertices
	#n_bits_v = convert(UInt8, ceil(log(2, gs)))
	# Get appropriate unsigned int type based on number of bits needed
	#V = infer_uint_custom_type(n_bits_v)
	
	if compression_scheme == 0x1
		g = load_huffman_compressed_mgs3_graph(f, graph_type, encoding, gs)
	elseif compression_scheme == 0x2
		g = load_elias_compressed_mgs3_graph(f, graph_type, encoding, gs, :elias_gamma)
	elseif compression_scheme == 0x3
		g = load_elias_compressed_mgs3_graph(f, graph_type, encoding, gs, :elias_delta)
	elseif compression_scheme == 0x4
		g =  load_golomb_compressed_mgs3_graph(f, graph_type, encoding, gs)
	elseif compression_scheme == 0x5
		g = load_fibonacci_compressed_mgs3_graph(f, graph_type, encoding, gs)
	elseif compression_scheme == 0x6
		g = load_zeta_compressed_mgs3_graph(f, graph_type, encoding, gs)
	else
		error("Unsupported compression scheme: $compression_scheme. Supported schemes are :huffman, :elias_gamma, :elias_delta, :golomb")
    end

	close(f)
	return g
end

""" 
    load_huffman_compressed_mgs3_graph(f::IO, graph_type::Symbol, encoding::Symbol, gs::UInt64)

load graph in compressed MGS v3 format (Huffman compression scheme)

Parameters:
- f: Input stream
- graph_type: Graph type (:directed or :undirected)
- encoding: Coding scheme (:children or :index)
- gs: Number of vertices
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
	# NB: stop sequence is the code associated to typemax(V)
	stop_seq = typemax(V)

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
    load_elias_compressed_mgs3_graph(io::IO, graph_type::Symbol, encoding::Symbol, gs::UInt64, compression::Symbol)

Load graph in compressed MGS v3 format with Elias coding scheme.

Parameters:
- io: Input stream
- graph_type: Graph type (:directed or :undirected)
- encoding: Coding scheme (:children or :index)
- gs: Number of vertices
- compression: Compression scheme (:elias_gamma or :elias_delta)
"""
function load_elias_compressed_mgs3_graph(io::IO, graph_type::Symbol, encoding::Symbol, gs::UInt64, compression::Symbol)
	# `n_size_u` is the number of bits needed to represent the graph vertices
	n_bits_v = convert(UInt8, ceil(log(2, gs)))
	# Get appropriate unsigned int type based on number of bits needed
	V = infer_uint_custom_type(n_bits_v)
	
	# intialize graph g according to graph type
	g = graph_type == :directed ? SimpleDiGraph{V}() : SimpleGraph{V}()
	
	# vertex set
	vs = range(1, stop=gs)

	@info("generating graph")
	@info("adding vertices")
	# add vertices to graph
    add_vertices!(g, gs)

	reader = BitReader(io)

	if encoding == :children
		# NB: stop sequence is the code associated to one(V)
		stop_seq = one(V)

		@info("reading data section")
		for v in vs
			source = convert(V, v)
			try
				# read the first neighbor value
				first_value = read_elias_coding(reader, compression, V)
				if first_value == stop_seq
					# go to next vertex
					continue
				end
				# NB: the first value is the first outneighbor of the vertex
				# NB: the deltas are shifted as 1 is reserved for the stop sequence
				# Subtract 1 from the delta to undo the shift, then use as vertex
				neighbor = first_value - 1
				add_edge!(g, source, neighbor)
				prev_value = neighbor
				
				# read subsequent neighbors as differences
				while true
					delta = read_elias_coding(reader, compression, V)
					if delta == stop_seq
						# go to next vertex
						break
					end
					# NB: the deltas are shifted as 1 is reserved for the stop sequence
					# Subtract 1 from the delta to undo the shift, then add to previous value
					prev_value += (delta - 1)
					add_edge!(g, source, prev_value)
				end
			catch e
				# do nothing
				#if !(isa(e, EOFError) || isa(e, ArgumentError))
				#	rethrow(e)
				#end
			end
		end
	elseif encoding == :index
		out_degrees = Dict{V,V}()
		@info("reading index section")
		for v in vs
			# NB: the out-degree is shifted as 0 is a forbidden value for Elias coding
			out_degrees[v] = read_elias_coding(reader, compression, V) - 1
		end
		@info("reading data section")
		for v in vs
			# NB: the out-neighbors are delta encoded
			degree = out_degrees[v]
			if degree == 0
				continue
			end
			
			# Read the neighbors for this vertex
			neighbors = V[]
			try
				# Read first value
				first_value = read_elias_coding(reader, compression, V)
				push!(neighbors, first_value)
				
				# Read remaining deltas
				for _ in 2:degree
					delta = read_elias_coding(reader, compression, V)
					push!(neighbors, neighbors[end] + delta)
				end
				
				# Add all edges for this vertex
				for neighbor in neighbors
					if neighbor > 0 && neighbor <= gs
						add_edge!(g, v, neighbor)
					end
				end
			catch
				# If we can't read data for this vertex, skip it and continue
				continue
			end
		end
	end

	close(io)
	return g
end

"""
    load_golomb_compressed_mgs3_graph(io::IO, graph_type::Symbol, encoding::Symbol, gs::UInt64)

Load graph in compressed MGS v3 format with Golomb coding scheme.

Parameters:
- io: Input stream
- graph_type: Graph type (:directed or :undirected)
- encoding: Coding scheme (:children or :index)
- gs: Number of vertices
"""
function load_golomb_compressed_mgs3_graph(io::IO, graph_type::Symbol, encoding::Symbol, gs::UInt64)
	# `n_size_u` is the number of bits needed to represent the graph vertices
	n_bits_v = convert(UInt8, ceil(log(2, gs)))
	# Get appropriate unsigned int type based on number of bits needed
	V = infer_uint_custom_type(n_bits_v)

	# intialize graph g according to graph type
	g = graph_type == :directed ? SimpleDiGraph{V}() : SimpleGraph{V}()

	# vertex set
	vs = range(1, stop=gs)
	
	@info("generating graph")
	@info("adding vertices")
	# add vertices to graph
    add_vertices!(g, gs)
	
	reader = BitReader(io)
	
	if encoding == :children
		# NB: stop sequence is 0
		stop_seq = zero(V)
		
		@info("reading data section")
		for v in vs
			source = convert(V, v)
			try
				# read the first neighbor value
				first_value = read_golomb(reader, GOLOMB_BASE, V)
				if first_value == stop_seq
					# go to next vertex
					continue
				end
				# ensure vertex is within bounds
				if first_value > 0 && first_value <= gs
					add_edge!(g, source, first_value)
					prev_value = first_value
					
					# read subsequent neighbors as differences
					while true
						delta = read_golomb(reader, GOLOMB_BASE, V)
						if delta == stop_seq
							# go to next vertex
							break
						end
						prev_value += delta
						# ensure vertex is within bounds
						if prev_value > 0 && prev_value <= gs
							add_edge!(g, source, prev_value)
						else
							error("Target vertex is out of bounds.")
						end
					end
				else
					error("Initial vertex is out of bounds.")
				end
			catch e
				# do nothing
				#if !(isa(e, EOFError) || isa(e, ArgumentError))
				#	rethrow(e)
				#end
			end
		end
	elseif encoding == :index
		@info("reading index section")
		out_degrees = Dict{V,V}()
		# read the out-degrees
		for v in vs
			out_degrees[v] = read_golomb(reader, GOLOMB_BASE, V)
		end
		@info("reading data section")
		for v in vs
			# NB: the out-neighbors are delta encoded
			degree = out_degrees[v]
			if degree == 0
				continue
			end
			
			# read the neighbors for this vertex
			neighbors = V[]
			try
				# read first value
				first_value = read_golomb(reader, GOLOMB_BASE, V)
				push!(neighbors, first_value)
				prev_value = first_value
				
				# read remaining deltas
				for _ in 2:degree
					delta = read_golomb(reader, GOLOMB_BASE, V)
					prev_value += delta
					push!(neighbors, prev_value)
				end
				
				# add all edges for this vertex
				for neighbor in neighbors
					if neighbor > 0 && neighbor <= gs
						add_edge!(g, v, neighbor)
					else
						error("Target vertex is out of bounds.")
					end
				end
			catch
				# if we can't read data for this vertex, throw an error
				error("Error reading data for vertex $v.")
			end
		end
	end

	close(io)
	return g
end

"""
    load_fibonacci_compressed_mgs3_graph(io::IO, graph_type::Symbol, encoding::Symbol, gs::UInt64)

Load graph in compressed MGS v3 format with Fibonacci coding scheme.

Parameters:
- io: Input stream
- graph_type: Graph type (:directed or :undirected)
- encoding: Coding scheme (:children or :index)
- gs: Number of vertices
"""
function load_fibonacci_compressed_mgs3_graph(io::IO, graph_type::Symbol, encoding::Symbol, gs::UInt64)
	# `n_size_u` is the number of bits needed to represent the graph vertices
	n_bits_v = convert(UInt8, ceil(log(2, gs)))
	# Get appropriate unsigned int type based on number of bits needed
	V = infer_uint_custom_type(n_bits_v)

	# intialize graph g according to graph type
	g = graph_type == :directed ? SimpleDiGraph{V}() : SimpleGraph{V}()

	# vertex set
	vs = range(1, stop=gs)

	@info("generating graph")
	@info("adding vertices")
	# add vertices to graph
    add_vertices!(g, gs)

	reader = BitReader(io)

	if encoding == :children
		# NB: stop sequence is 1
		stop_seq = one(V)
		
		@info("reading data section")
		for v in vs
			source = convert(V, v)
			try
				# read the first neighbor value
				first_value = read_fibonacci_code(reader, V)
				if first_value == stop_seq
					# go to next vertex
					continue
				end
				# NB: the first value is the first outneighbor of the vertex
				neighbor = first_value - 1
				add_edge!(g, source, neighbor)
				prev_value = neighbor

				# read subsequent neighbors as differences
				while true
					delta = read_fibonacci_code(reader, V)
					if delta == stop_seq
						# go to next vertex
						break
					end
					prev_value += (delta - 1)	
					add_edge!(g, source, prev_value)
				end
			catch e
				# do nothing
				#if !(isa(e, EOFError) || isa(e, ArgumentError))
				#	rethrow(e)
				#end
			end
		end
	elseif encoding == :index
		@info("reading index section")
		out_degrees = Dict{V,V}()
		# read the out-degrees
		for v in vs
			out_degrees[v] = read_fibonacci_code(reader, V)
		end
		@info("reading data section")
		for v in vs
			# NB: the out-neighbors are delta encoded
			degree = out_degrees[v]
			if degree == 0
				continue
			end

			# read the neighbors for this vertex
			neighbors = V[]
			try
				# read first value
				first_value = read_fibonacci_code(reader, V)
				push!(neighbors, first_value - 1)
				prev_value = first_value - 1
				
				# read remaining deltas
				for _ in 2:degree
					delta = read_fibonacci_code(reader, V)
					prev_value += (delta - 1)
					push!(neighbors, prev_value)
				end
				
				# add all edges for this vertex
				for neighbor in neighbors
					if neighbor > 0 && neighbor <= gs
						add_edge!(g, v, neighbor)
					else
						error("Target vertex is out of bounds.")
					end
				end
			catch e
				# do nothing
				#if !(isa(e, EOFError) || isa(e, ArgumentError))
				#	rethrow(e)
				#end
			end
		end
	end

	close(io)
	return g
end

"""
    load_zeta_compressed_mgs3_graph(io::IO, graph_type::Symbol, encoding::Symbol, gs::UInt64)

Load graph in compressed MGS v3 format with Zeta coding scheme.

Parameters:
- io: Input stream
- graph_type: Graph type (:directed or :undirected)
- encoding: Coding scheme (:children or :index)
- gs: Number of vertices
- k: Zeta parameter (default: 3)
"""
function load_zeta_compressed_mgs3_graph(io::IO, graph_type::Symbol, encoding::Symbol, gs::UInt64, k::Int=ZETA_BASE)
	# `n_size_u` is the number of bits needed to represent the graph vertices
	n_bits_v = convert(UInt8, ceil(log(2, gs)))
	# Get appropriate unsigned int type based on number of bits needed
	V = infer_uint_custom_type(n_bits_v)

	# intialize graph g according to graph type
	g = graph_type == :directed ? SimpleDiGraph{V}() : SimpleGraph{V}()

	# vertex set
	vs = range(1, stop=gs)

	@info("generating graph")
	@info("adding vertices")
	# add vertices to graph
    add_vertices!(g, gs)

	reader = BitReader(io)

	if encoding == :children
		# NB: stop sequence is 1
		stop_seq = one(V)
		
		@info("reading data section")
		for v in vs
			source = convert(V, v)
			try
				# read the first neighbor value
				first_value = read_zeta_coding(reader, k, V)
				if first_value == stop_seq
					# go to next vertex
					continue
				end
				# NB: the first value is the first outneighbor of the vertex
				neighbor = first_value - 1
				add_edge!(g, source, neighbor)
				prev_value = neighbor

				# read subsequent neighbors as differences
				while true
					delta = read_zeta_coding(reader, k, V)
					if delta == stop_seq
						# go to next vertex
						break
					end
					prev_value += (delta - 1)	
					add_edge!(g, source, prev_value)
				end
			catch e
				# do nothing
				#if !(isa(e, EOFError) || isa(e, ArgumentError))
				#	rethrow(e)
				#end
			end
		end
	elseif encoding == :index
		@info("reading index section")
		out_degrees = Dict{V,V}()
		# read the out-degrees
		for v in vs
			out_degrees[v] = read_zeta_coding(reader, k, V)
		end
		@info("reading data section")
		for v in vs
			# NB: the out-neighbors are delta encoded
			degree = out_degrees[v]
			if degree == 0
				continue
			end

			# read the neighbors for this vertex
			neighbors = V[]
			try
				# read first value
				first_value = read_zeta_coding(reader, k, V)
				push!(neighbors, first_value - 1)
				prev_value = first_value - 1
				
				# read remaining deltas
				for _ in 2:degree
					delta = read_zeta_coding(reader, k, V)
					prev_value += (delta - 1)
					push!(neighbors, prev_value)
				end
				
				# add all edges for this vertex
				for neighbor in neighbors
					if neighbor > 0 && neighbor <= gs
						add_edge!(g, v, neighbor)
					else
						error("Target vertex is out of bounds.")
					end
				end
			catch e
				# do nothing
				#if !(isa(e, EOFError) || isa(e, ArgumentError))
				#	rethrow(e)
				#end
			end
		end
	end

	close(io)
	return g
end

end