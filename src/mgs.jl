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
using ..Util: infer_uint_custom_type, to_bytes, huffman_encoding, encode_tree!, decode_tree!, get_huffman_codes!, decode_values  # Add these imports
using ..Graph: get_basic_stats  # Add specific imports from Graph if needed

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
# - D1 (directed graph + Huffman compression)
# - CS0 (coding scheme 0) = 0x00 (data section only + reserved)
HEADER_MGS3_DH_CS0 = 0x4d475303000100
# 'MGS' + 0x0300 (major=3, minor=0) 
# - D0 (directed graph + Huffman compression)
# - CS1 (coding scheme 1) = 0x10 (index and data sections + reserved)
HEADER_MGS3_DH_CS1 = 0x4d475303000110

MGS_MAX_SIZE = 0xffffffffff

# Export the functions we want to make available
export write_mgs3_graph, 
       write_compressed_mgs3_graph,
       load_mgs3_graph,
       load_compressed_mgs3_graph,
       write_huffman_compressed_mgs3_graph

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
	#	 * Byte 1:
	#    	- graph type: 			0x0: directed graph | 0x1: undirected graph
	#	 	- compression scheme: 	0x0: no compression
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

	f = open("$filename.mgs", "w")
	
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

	# graph type is 4 first bits of 6th byte of header	
	# 0x0: directed graph | 0x1: undirected graph
	graph_type = version[6] >> 4
	# coding scheme is 4 first bits of 7th byte of header	
	coding_scheme = version[7] >> 4

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
	if coding_scheme == 0x00
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
	elseif coding_scheme == 0x01
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
		for i in 1:length(vs)
			source = convert(V, i)
			# if we reached the last parent vertex
			if i == length(vs)
				pos1 = current_pos
				pos2 = length(children)
			else
				pos1 = current_pos
				# position of the last child of vertex i
				pos2 = current_pos + ods[i+1] - 1
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

"""
    write_compressed_mgs3_graph(g::AbstractGraph{T}, rg::AbstractGraph{T}, filename::AbstractString, compression::Symbol=:huffman, coding_scheme::UInt8=0x00) where {T<:Unsigned}

Write graph in MGS v3 format with specified encoding and compression scheme

Supported encoding schemes:
- :children
- :index

Supported compression schemes:
- :huffman - Huffman coding (default)
- :elias - Elias gamma coding
- :golomb - Golomb coding

Parameters:
- g: Input graph
- rg: Reverse graph (required for some compression schemes)
- filename: Output filename
- encoding: Coding scheme (:children for children section only, :index for index+children sections)
- compression: Compression scheme to use (default: :huffman)
"""
function write_compressed_mgs3_graph(g::AbstractGraph{T}, rg::AbstractGraph{T}, filename::AbstractString, encoding::Symbol=:children, compression::Symbol=:huffman) where {T<:Unsigned}
    if compression == :huffman
        write_huffman_compressed_mgs3_graph(g, rg, filename, encoding)
    elseif compression == :elias
        error("Elias gamma coding not yet implemented")
    elseif compression == :golomb
        error("Golomb coding not yet implemented")
    else
        error("Unsupported compression scheme: $compression. Supported schemes are :huffman, :elias, :golomb")
    end
end


"""
    write_huffman_compressed_mgs3_graph(g::AbstractGraph{T}, rg::AbstractGraph{T}, filename::AbstractString, coding_scheme::UInt8=0x00) where {T<:Unsigned}

write graph in a compressed MGS v3 (Huffman compression scheme)

Supported encoding schemes:
- :children
- :index
"""
function write_huffman_compressed_mgs3_graph(g::AbstractGraph{T}, rg::AbstractGraph{T}, filename::AbstractString, encoding::Symbol=:children) where {T<:Unsigned}
	# Header 12 bytes: 
	# -> version: 'MGS' 3 bytes string
	# -> major + minor version: 2 bytes
	# -> flags: 2 bytes
	#	 * Byte 1:
	#    	- graph type: 			0x0: directed graph | 0x1: undirected graph
	#	 	- compression scheme: 	0x1: Huffman compression
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
		# 'MGS' + 0x0300 + 0x00 (directed graph + compression) + 0x00 (data section only + reserved)
		# encoding: data section only with implicit numbering of vertices
  		version = HEADER_MGS3_DH_CS0
	elseif enccoding == :index
		# 'MGS' + 0x0300 + 0x00 (directed graph + compression) + 0x00 (index and data sections + reserved)
		# encoding: index+data sections with implicit numbering of vertices
  		version = HEADER_MGS3_DH_CS1
	end

	f = open("$filename.mgz", "w")

	# vertices in-degree
	in_degrees = T[]
	# number of children for each vertex
	# NB: `T` should have a sufficient size to store the number of children
	ods = T[]
	# flattened list of children for all the vertices
	# NB: `T` should have a sufficient size to store the children indices
	children = T[]

	for v in vs
		ovs = outneighbors(g, v)
		push!(ods, convert(T, length(ovs)))
		push!(in_degrees, length(outneighbors(rg, v)))
		for o in ovs
			push!(children, o)	
		end
	end

	# if coding scheme is :children, we need a special stop sequence for each vertex in the data section
	if encoding == :children
		# this stop sequence will appear `gs` times in the data section
		# let's add to the in-degrees list the number of vertices
		push!(in_degrees, convert(T, gs))
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
	for p in in_degrees
		#bytes = reverse(reinterpret(UInt8, [p]))[diff_size+1:end]
		bytes = to_bytes(p)[(diff_size+1):end]
		write(f, bytes)
	end

	if encoding == :children
		@info("writing data section with stop sequence")
		### write data section
		cdata = BitArray{1}()
		# stop sequence is the code associated to the number of vertices `gs`
		# NB: frequencies are in the range [0, gs-1] and the stop sequence is the code associated to `gs`
		# NB: R has a length of `gs+1` because of the stop sequence
		stop_seq_code = R[gs+1]
		
		for v in vs
			ovs = outneighbors(g, v)
			for c in ovs
				# get code associated to child vertex id
				code = R[c]
				append!(cdata, code)
			end
			append!(cdata, stop_seq_code)
		end
	elseif encoding == :index
		# if coding scheme is 0x01, we need to write the index section
		@info("writing index section")
		### write index section
		for p in ods
			#bytes = reverse(reinterpret(UInt8, [p]))[diff_size+1:end]
			bytes = to_bytes(p)[(diff_size+1):end]
			write(f, bytes)
		end
		
		@info("writing data section")
		### write data section
		cdata = BitArray{1}()
		for c in children
			# get code associated to child id
			code = R[c]
			append!(cdata, code)
		end
	end
	
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
    load_compressed_mgs3_graph(filename::AbstractString, compression::Symbol=:huffman)

Load graph in MGS v3 format with specified compression scheme.
Supported compression schemes:
- :huffman - Huffman coding (default)
- :elias - Elias gamma coding
- :golomb - Golomb coding

Returns a graph loaded with the specified compression scheme.
"""
function load_compressed_mgs3_graph(filename::AbstractString, compression::Symbol=:huffman)
    if compression == :huffman
        return load_huffman_compressed_mgs3_graph(filename)
    elseif compression == :elias
        error("Elias gamma coding not yet implemented")
    elseif compression == :golomb
        error("Golomb coding not yet implemented") 
    else
        error("Unsupported compression scheme: $compression. Supported schemes are :huffman, :elias, :golomb")
    end
end

""" 
    load_huffman_compressed_mgs3_graph(filename::AbstractString)

load graph in compressed MGS v3 format (Huffman compression scheme)
"""
function load_huffman_compressed_mgs3_graph(filename::AbstractString)
	f = open(filename, "r")
	### read header
  	# 7 bytes: <3 bytes string 'MGS'> + <2 bytes major/minor> + <2 bytes flags>
	# 5 bytes (40 bits): number of vertices
	###
	# 7-bytes version
	version = read(f,7)
	# number of vertices
	gs = reinterpret(UInt64, vcat(reverse(read(f,5)),[0x00,0x00,0x00]))[1]
	# `n_size_u` is the number of bits needed to represent the graph vertices
	n_bits_v = convert(UInt8, ceil(log(2, gs)))
	# Get appropriate unsigned int type based on number of bits needed
	V = infer_uint_custom_type(n_bits_v)
	
	# graph type is 4 first bits of 6th byte of header	
	# 0x0: directed graph | 0x1: undirected graph
	graph_type = version[6] >> 4

	# coding scheme is 4 first bits of 7th byte of header	
	coding_scheme = version[7] >> 4

	# intialize graph g according to graph type
	g = if graph_type == 0x0
		SimpleDiGraph{V}()
	else
		SimpleGraph{V}()
	end
	
	# vertex set
	vs = range(1, stop=gs)

	@info("generating graph")
	@info("adding vertices")
	# add vertices to graph
    add_vertices!(g, gs)

	# read frequency section
	F = V[]
	for i in 1:(gs+1)
		p = read(f, sizeof(V))
		push!(F, reinterpret(V, reverse(p))[1])
	end
	
	ods = V[]
	if coding_scheme == 0x01
		@info("reading index section")
		# read index
		# NB: 
		# - position indices are 1-based
		for i in 1:gs
			p = read(f, sizeof(V))
			push!(ods, reinterpret(V, reverse(p))[1])
		end
	end

	@info("reading data section")
	# read data section
	CDATA = BitArray{1}()
	while !eof(f)
		b = read(f, 1)[1]
		for j in 0:7
			if ((b >> j) & 0x01) == 1
				push!(CDATA, true)
			else
				push!(CDATA, false)
			end
		end
	end
	close(f)

	# get last byte number of 0s
	sp = 8 - sum(CDATA[end-7:end])
	CDATA = CDATA[1:end-(7+sp)]
	
	@info("generating Huffman tree")
	# get Huffman encoding tree
	tree = huffman_encoding(F)
	
	@info("decoding values")
	# decode values
	children = decode_values(tree, CDATA)

	@info("adding edges")
	# add edges
	if coding_scheme == 0x00
		pos = 1
		n_children = length(children)
		# NB: stop sequence is the code associated to `gs+1`
		stop_seq = convert(V, gs+1)

		for i in 1:length(vs)
			if pos <= n_children
				source = convert(V,i)
				while children[pos] != stop_seq
					target = children[pos]
					add_edge!(g, source, target)
					pos += 1
				end
				# skip stop sequence
				pos += 1
			else
				# if we reached the last child, we are done
				break
			end
		end
	elseif coding_scheme == 0x01
		current_pos = 1
		for i in 1:length(vs)
			source = convert(V,i)
			# if we reached the last parent vertex
			if i == length(vs)
				pos1 = current_pos
				pos2 = length(children)
			else
				pos1 = current_pos
				# position of the last child of vertex i
				pos2 = current_pos + ods[i+1] - 1
				current_pos = pos2 + 1
			end
			# add edges for each child
			for p in pos1:pos2
				target = children[p]
				add_edge!(g, source, target)
			end
		end
	end
	return g::AbstractGraph{V}
end

end