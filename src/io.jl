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

using LightGraphs, DataStructures, HDF5, JLD

include("util.jl")

# constants
HEADER_MGS3_CS0_C0 = 0x4d475303000000
HEADER_MGS3_CS1_C0 = 0x4d475303000010

MGS3_MAX_SIZE = 0xffffffffff

""" 
    load_jls_serialized(filename::AbstractString)

load serialized JLS data
"""
function load_jls_serialized(filename::AbstractString)
	x = open(filename, "r") do file
		deserialize(file)
	end
	return x
end

""" 
    serialize_to_jls(x::Any, filename::AbstractString)

serialize data to JLS format
"""
function serialize_to_jls(x::Any, filename::AbstractString)
	open("$filename.jls", "w") do file
		serialize(file, x)
	end
end

""" 
    load_jld_serialized(name::AbstractString, filename::AbstractString)

load serialized JLD data

NB: to be favored for long term storage
"""
function load_jld_serialized(name::AbstractString, filename::AbstractString)
	x = jldopen(filename, "r") do file
    		read(file, name)
  	end
	return x
end

""" 
    serialize_to_jld(x::Any, name::AbstractString, filename::AbstractString)

serialize data to JLD format

NB: to be favored for long term storage
"""
function serialize_to_jld(x::Any, name::AbstractString, filename::AbstractString)
	jldopen("$filename.jld", "w") do file
		write(file, name, x)
	end
end

"""
    write_mgs3_graph(g::AbstractGraph{T},filename::AbstractString) where {T<:Unsigned}

write graph in format MGS v3
"""
function write_mgs3_graph(g::AbstractGraph{T}, filename::AbstractString, coding_scheme::UInt8=0x00) where {T<:Unsigned}
  	# Header 12 bytes: 
	# -> version: 'MGS' 3 bytes string
	# -> major + minor version: 2 bytes
	# -> flag: 2 bytes
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
	if gs > MGS3_MAX_SIZE
		error("Input graph cannot have more than 2^40-1 vertices")
	end
	
	# `n_size_t` is the number of bytes needed to represent the graph vertices
	# NB: `n_size_t` (the computed size) may be lower than `sizeof(T)` of the graph type in parameter
	n_size_t = convert(UInt8, ceil(log(2, gs)))
	# NB: `T` should be an unsigned integer of size 1,2,4,8 bytes
	p_size_t = sizeof(T)
	
	#  The difference of size should >= 0 
	diff_size = p_size_st - n_size_t

	if coding_scheme == 0x00
		# 'MGS' + 0x0300 + 0x00 (directed graph + no compression) + 0x00 (data section only + reserved) (7 bytes)
		# encoding: data section only with implicit numbering of vertices
  		version = HEADER_MGS3_CS0_C0
	elseif coding_scheme == 0x01
		# 'MGS' + 0x0300 + 0x00 (directed graph + no compression) + 0x00 (index and data sections + reserved) (7 bytes)
		# encoding: index+data sections with implicit numbering of vertices
  		version = HEADER_MGS3_CS0_C1
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
	if coding_scheme == 0x00
		stop_seq = [0x00 for i in 1:n_size_t]
		for v in vs
			ovs = outneighbors(g, v)
			for c in ovs
				bytes = reverse(reinterpret(UInt8, [c]))[(diff_size+1):p_size_t]
				write(f, bytes)
			end
			write(f, stop_seq)
		end

	# coding scheme: index+data sections with implicit numbering of vertices
	elseif coding_scheme == 0x01
		# number of children for each vertex
		# NB: `T` should have a sufficient size to store the number of children
		ods = T[]
		# flattened list of children for all the vertices
		# NB: `T` should have a sufficient size to store the children indices
		children = T[]

		for v in vs
			ovs = outneighbors(g, v)
			push!(ods, length(ovs))
			for o in ovs
				push!(children, o)	
			end
		end

		### write index section
		for o in ods
			bytes = reverse(reinterpret(UInt8, [o]))[(diff_size+1):p_size_t]
			write(f, bytes)
		end
		### write data section
		for c in children
			bytes = reverse(reinterpret(UInt8, [c]))[(diff_size+1):p_size_t]
			write(f, bytes)
		end
	end

	close(f)
end

""" 
    load_mgs3_graph(filename::AbstractString)::AbstractGraph{U} where {U<:Unsigned}

load graph in format MGS v3
"""
function load_mgs3_graph(filename::AbstractString)::AbstractGraph{U} where {U<:Unsigned}
	f = open(filename, "r")
	### read header
  	# 7 bytes: <3 bytes string 'MGS'> + <2 bytes major/minor> + <2 bytes flags>
	# 5 bytes (40 bits): number of vertices
	###
	# 7-bytes version
	version = read(f,7)
	# number of vertices
	gs = reinterpret(UInt64, vcat(reverse(read(f,5)),[0x00,0x00,0x00]))
	
	# `n_size_u` is the number of bits needed to represent the graph vertices
	n_bits_u = convert(UInt8, ceil(log(2, length(gs))))
	# `U` should be an unsigned integer of size 1,2,4,8 bytes
	U = infer_uint_type(n_bits_u)

	# graph type is 4 first bits of 6th byte of header	
	# 0x0: directed graph | 0x1: undirected graph
	graph_type = version[6] >> 4

	# coding scheme is 4 first bits of 7th byte of header	
	coding_scheme = version[7] >> 4

	# intialize graph g according to graph type
	if graph_type == 0x0
		g = SimpleDiGraph{U}()
	else
		g = SimpleGraph{U}()
	end
	
	# vertex set
	vs = range(1, stop=gs)

	# add vertices to graph
    add_vertices!(g, length(vs))

	# coding scheme: data section only with implicit numbering of vertices
	# NB: each list of children is terminated with 0
	if coding_scheme == 0x00
		# read data
		children = U[]
		while !eof(f)
			c = read(f, sizeof(U))
			push!(children, reinterpret(U,c)[1])
		end
		
		# add edges
		pos = 1
		n_children = length(children)

		for i in 1:length(vs)
			if pos <= n_children
				source = convert(U,i)
				while children[pos] != 0x00
					target = children[pos]
					add_edge!(g,source,target)
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
		pos = U[]
		for i in 1:gs
			p = read(f, sizeof(U))
			push!(pos, reinterpret(U,p)[1])
		end
		# read data
		children = U[]
		while !eof(f)
			c = read(f, sizeof(U))
			push!(children, reinterpret(U,c)[1])
		end

		# add edges
		for i in 1:length(vs)
			source = convert(U,i)
			# if we reached the last parent vertex
			if i == length(vs)
				pos1 = pos[i]
				pos2 = length(children)
			else
				pos1 = pos[i]
				# position of the last child of vertex i
				pos2 = pos[i+1]-1
			end
			# add edges for each child
			for p in pos1:pos2
				target = children[p]
				add_edge!(g,source,target)
			end
		end
	end
		
	close(f)

	return g
end

"""
    write_mgs4_graph(g::AbstractGraph{T},rg::AbstractGraph{T},filename::AbstractString) where {T<:Unsigned}

write graph in format MGS v4
"""
function write_mgs4_graph(g::AbstractGraph{T},rg::AbstractGraph{T},filename::AbstractString) where {T<:Unsigned}
  	# Header 12 bytes: 
	# -> version: 'MGS' 3 bytes string
	# -> major + minor version: 2 bytes
	# -> reserved flag: 2 bytes
	# -> size_t in bytes: 1 byte
	# -> # vertices: 4 bytes position 
	#
	#'MGS' + <16 bits major|minor version> + <8 bits size_t> 
	# 	+ <32 bits offset in size_t of data section>

	# MGS + 0x0400 + 0x0000
  	version = 0x4d475304000000
	# size of type T in bytes
	size_t = convert(UInt8, sizeof(T))
	
	vs = vertices(g)
	# offset of children 
	pos = T[]
	# array of children
	children = T[]
	# vertices in-degree
	in_degrees = T[]
	cpos = convert(T,1)
	for v in vs
		ovs = outneighbors(g,v)
		push!(pos,cpos)
		push!(in_degrees,length(outneighbors(rg,v)))
		for o in ovs
			push!(children,o)	
			cpos += convert(T,1)
		end
	end

	@info("generating Huffman tree")
	# get Huffman encoding tree
	tree = huffman_encoding(in_degrees)

	@info("getting Huffman codes")
	# get Huffman codes in C (C:code -> value)
	C = Dict{BitArray{1},T}()
	get_huffman_codes!(tree, C, BitArray{1}())

	# reverse dictionary (R: value -> code)
	R = Dict{T,BitArray{1}}()
	[R[value] = key for (key, value) in C]
	
	# number of vertices
	gs = convert(UInt32, length(vs))

	f = open("$filename.mgz", "w")
	
	@info("writing header section")
	### write header
	# reinterpret generates an array of length 8 even if version has a length of 7 bytes
	# write version + size of type T (8 bytes)
	bytes = reinterpret(UInt8, [version])[1:7]
	write(f, bytes)
	# size of type T (1 byte)
	bytes = reinterpret(UInt8, [size_t])
	write(f, bytes)
	# write number of vertices (4 bytes)
	bytes = reinterpret(UInt8, [gs])
	write(f, bytes)
	
	@info("writing frequency section")
	### write frequency section
	for p in in_degrees
		bytes = reinterpret(UInt8, [p])
		write(f, bytes)
	end
	
	@info("writing index section")
	### write index section
	for p in pos
		bytes = reinterpret(UInt8, [p])
		write(f, bytes)
	end
	
	@info("writing data section")
	### write data section
	cdata = BitArray{1}()
	for c in children
		# get code associated to child id
		code = R[c]
		append!(cdata,code)
	end
	
	# add padding if necessary
	scd = convert(UInt32, length(cdata))
	sp = 8-scd%8
	for i in 1:sp
		push!(cdata,0)
	end
	
	# number of bytes to write
	nb = round(Int, length(cdata)/8)
	@debug("# bytes to write: $nb")
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
    load_mgs4_graph(g::AbstractGraph{T},filename::AbstractString) where {T<:Unsigned}

load graph in format MGS v4
"""
function load_mgs4_graph(g::AbstractGraph{T},filename::AbstractString) where {T<:Unsigned}
	f = open(filename, "r")
	### read header
  	# 8 bytes: <3 bytes string 'MGS'> + <2 bytes major/minor> + <2 bytes flags> + <1 byte size_t>
	# 4 bytes (32 bits): offset in size_t of graph data section (i.e. # of vertices)
	###
	# 7-bytes version
	version = read(f,7)
	# size in bytes of T
	size_t = convert(UInt8, read(f,1)[1])
	# number of vertices
	gs = reinterpret(UInt32, read(f,4))[1]
	
	# read frequency section
	F = T[]
	for i in 1:gs
		p = read(f,sizeof(T))
		push!(F,reinterpret(T,p)[1])
	end
	
	# read index section
	pos = T[]
	for i in 1:gs
		p = read(f,sizeof(T))
		push!(pos,reinterpret(T,p)[1])
	end

	# read data section
	CDATA = BitArray{1}()
	while !eof(f)
		b = read(f,1)[1]
		for j in 0:7
			if ((b >> j) & 0x01) == 1
				push!(CDATA,true)
			else
				push!(CDATA,false)
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
	
	@info("generating graph")
	# vertex set
	vs = range(1,stop=length(pos))

	@info("adding vertices")
	# add vertices
    add_vertices!(g,length(vs))

	@info("adding edges")
	# add edges
	for i in 1:length(vs)
		source = convert(T,i)
		# if we reached the last parent vertex
		if i == length(vs)
			pos1 = pos[i]
			pos2 = length(children)
		else
			pos1 = pos[i]
			pos2 = pos[i+1]-1
		end
		for p in pos1:pos2
			target = children[p]
			add_edge!(g,source,target)
		end
	end
end

""" 
    load_adjacency_list_from_csv(::Type{T},g::AbstractGraph{T},filename::AbstractString,separator::Char=',') where {T<:Unsigned}

load graph from CSV adjacency list
"""
function load_adjacency_list_from_csv(::Type{T},g::AbstractGraph{T},filename::AbstractString,separator::AbstractChar=',') where {T<:Unsigned}
	f = open(filename,"r")
	oni = Dict{T,T}()
	edges = Array{Tuple{T,T},1}()
	counter = convert(T,1)
	while !eof(f)
		line = strip(readline(f))
		if !startswith(line, "#")
			edge = split(line, separator)
			v1 = parse(T,edge[1])
			v2 = parse(T,edge[2])

			if !haskey(oni, v1)
				oni[v1] = counter
				counter += convert(T,1)
			end
			if !haskey(oni, v2)
				oni[v2] = counter
				counter += convert(T,1)
			end
			push!(edges, (oni[v1], oni[v2]))
		end
	end
	close(f)
	
	# add vertices
	add_vertices!(g,length(keys(oni)))
	
	# add edges
	for edge in edges
		add_edge!(g, edge[1], edge[2])	
	end

	return oni
end

""" 
    load_adjacency_list(::Type{T},g::AbstractGraph{T},adj_list::Array{T,2}) where {T<:Unsigned}

load graph from adjacency list

NB: adjcency list should be represented as a 2-dimensional array with 2 rows and 1 column per edge
"""
function load_adjacency_list(::Type{T},g::AbstractGraph{T},adj_list::Array{T,2}) where {T<:Unsigned}
	oni = Dict{T,T}()
	edges = Array{Tuple{T,T},1}()
	counter = convert(T,1)
	nes = size(adj_list)[2]
	for i in 1:nes
		edge = adj_list[i]
		v1 = parse(T,edge[1])
		v2 = parse(T,edge[2])

		if !haskey(oni, v1)
			oni[v1] = counter
			counter += convert(T,1)
		end
		if !haskey(oni, v2)
			oni[v2] = counter
			counter += convert(T,1)
		end
		push!(edges, (oni[v1], oni[v2]))
	end
	
	# add vertices
	add_vertices!(g,length(keys(oni)))
	
	# add edges
	for edge in edges
		add_edge!(g, edge[1], edge[2])	
	end

	return oni
end

""" 
    load_graph_from_pajek(::Type{T},g::AbstractGraph{T},filename::AbstractString) where {T<:Unsigned}

Load net Pajek file
"""
function load_graph_from_pajek(::Type{T},g::AbstractGraph{T},filename::AbstractString) where {T<:Unsigned}
	f = open(filename,"r")
	inside_vertices_section = false
	inside_edges_section = false
	vdict = Dict{UInt32, UInt32}()
	vcounter = convert(T,1)
	while !eof(f)
		line = lowercase(strip(readline(f)))
		if !startswith(line, "%")
			if startswith(line,"*vertices")
				inside_vertices_section = true
				continue
			elseif startswith(line,"*arcs")
				inside_edges_section = true
				inside_vertices_section = false
				continue
			end
			if inside_vertices_section
				sa = split(line, ' ')
				vdict[parse(T,sa[1])] = vcounter
				vcounter += convert(T,1)	
        		add_vertex!(g)
			elseif inside_edges_section
				sa = split(line, ' ')
				vs = vdict[parse(T,sa[1])]
				vt = vdict[parse(T,sa[2])]
				add_edge!(g,vs,vt)
			end
		end
	end
	close(f)
end

""" 
    load_triangles(::Type{T},filename::AbstractString) where {T<:Unsigned}

load triangles list from CSV text-formatted file 
"""
function load_triangles(::Type{T},filename::AbstractString) where {T<:Unsigned}
	f = open(filename,"r")
	a = (T,T,T)[]
	while !eof(f)
		te = split(readline(f)[2:(end-2)],',')
		v1 = parse(T,te[1])
		v2 = parse(T,te[2])
		v3 = parse(T,te[3])
		push!(a,(v1,v2,v3))
	end
	close(f)
	return a
end
