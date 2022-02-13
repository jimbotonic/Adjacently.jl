#
# Adjacently: Julia Complex Networks Library
# Copyright (C) 2016-2022 Jimmy Dubuisson <jimmy.dubuisson@gmail.com>
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

# MGS adjacency graph format - version 1
#
# unsigned integers stored in big endian format (byte and bit levels)
# positions: 0-based
#
# graph.index: vid1 (1 word), pos1 (1 word) | ...
# graph.data: vid1, vid2, vid3, ...
#
# NB: vertices do not need to be contiguous
#

# MGS adjacency graph format - version 2
#
# unsigned integers stored in big endian format (byte and bit levels)
# positions: 1-based
#
# graph.index: pos1 (1 word) | pos2 | ...
# graph.data: vid1 (1 word), vid2, vid3, ...
#
# NB: vertex ids need to be contiguous. In case, a vertex has no child, pos[k] = pos[k+1]

using LightGraphs, DataStructures, HDF5, JLD

include("util.jl")

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
    load_mgs1_graph_index(ipos::OrderedDict{T,T},filename::AbstractString) where {T<:Unsigned}

get the ordered dictionary vid -> startpos (mgs v1)
"""
function load_mgs1_graph_index(ipos::OrderedDict{T,T},filename::AbstractString) where {T<:Unsigned}
	f = open(filename, "r")
	while !eof(f1)
		id = read(f,UInt8,sizeof(T))
		pos = read(f,UInt8,sizeof(T))
		# Julia is 1-based -> +1
		ipos[reinterpret(T,reverse(id))[1]] = reinterpret(T,reverse(pos))[1] + 1
	end
	close(f)
end

""" 
    load_mgs2_graph_index(pos::Array{T,1},filename::AbstractString) where {T<:Unsigned}

get the set of positions (mgs v2)
"""
function load_mgs2_graph_index(pos::Array{T,1},filename::AbstractString) where {T<:Unsigned}
	f = open(filename, "r")
	while !eof(f)
		p = read(f,UInt8,sizeof(T))
		push!(pos, reinterpret(T,reverse(p))[1])
	end
	close(f)
end

""" 
    write_mgs1_graph_index(ipos::OrderedDict{T,T}, filename::AbstractString) where {T<:Unsigned}

write index file (mgs v1)
"""
function write_mgs1_graph_index(ipos::OrderedDict{T,T}, filename::AbstractString) where {T<:Unsigned}
	f = open(filename, "w")
	for p in ipos
    # reinterpret pair (vid,pos) in an array of bytes
		bytes = reinterpret(UInt8, [p])
		write(f, reverse(bytes))
	end
	close(f)
end

""" 
    write_mgs2_graph_index(pos::Array{T,1}, filename::AbstractString) where {T<:Unsigned}

write index file (mgs v2)
"""
function write_mgs2_graph_index(pos::Array{T,1}, filename::AbstractString) where {T<:Unsigned}
	f = open(filename, "w")
	for p in pos
    # reinterpret pos in an array of bytes
		bytes = reinterpret(UInt8, [p])
		write(f, reverse(bytes))
	end
	close(f)
end

""" 
    load_graph_data(children::Array{T,1},filename::AbstractString) where {T<:Unsigned}

get the array of graph children (MGSv1 & MGSv2)
"""
function load_graph_data(children::Array{T,1},filename::AbstractString) where {T<:Unsigned}
	f = open(filename, "r")
	while !eof(f)
		child = read(f,UInt8,sizeof(T))
		push!(children,reinterpret(T,reverse(child))[1])
	end
	close(f)
end

""" 
    write_graph_data(children::Array{T,1}, filename::AbstractString) where {T<:Unsigned}

write the array of graph children (MGSv1 & MGSv2)
"""
function write_graph_data(children::Array{T,1}, filename::AbstractString) where {T<:Unsigned}
	f = open(filename, "w")
	for c in children
		bytes = reinterpret(UInt8, [c])
		write(f, reverse(bytes))
	end
	close(f)
end

""" 
    load_mgs1_graph(g::AbstractGraph{T},name::AbstractString) where {T<:Unsigned}

load graph in format MGS v1
"""
function load_mgs1_graph(g::AbstractGraph{T},name::AbstractString) where {T<:Unsigned}
	ipos = OrderedDict{T,T}()
	children = T[]

	load_mgs1_graph_index(ipos,"$name.index")
	load_graph_data(children,"$name.data")

	# ipos is an ordered dictionary
	ks = collect(keys(ipos))
	# vertex set
	vs = sort(union(children,ks))

	# add vertices
    add_vertices!(g,length(vs))

	# dictionary old -> new vertex indices
	oni = Dict{T,T}()

	counter = convert(T,1)
	for v in vs
		oni[v] = counter
		counter += convert(T,1)
	end

	# add edges
	for i in 1:length(ks)
		source = oni[ks[i]]
		# if we reached the last parent vertex
		if i == length(ks)
			pos1 = ipos[ks[i]]
			pos2 = length(children)
		else
			pos1 = ipos[ks[i]]
			pos2 = ipos[ks[i+1]]-1
		end
		for p in pos1:pos2
			target = oni[children[p]]
			add_edge!(g,source,target)
		end
	end

	return oni
end

""" 
    load_mgs2_graph(g::AbstractGraph{T},name::AbstractString) where {T<:Unsigned}

load graph in format MGS v2
"""
function load_mgs2_graph(g::AbstractGraph{T},name::AbstractString) where {T<:Unsigned}
	pos = T[]
	children = T[]

	load_mgs2_graph_index(pos,"$name.index")
	load_graph_data(children,"$name.data")

	#@debug("#pos:",length(pos))
	#@debug("#children:",length(children))

	# vertex set
	vs = range(1,length(pos))

	# add vertices
	add_vertices!(g,length(vs))

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
    write_mgs3_graph(g::AbstractGraph{T},filename::AbstractString) where {T<:Unsigned}

load graph in format MGS v3
"""
function write_mgs3_graph(g::AbstractGraph{T},filename::AbstractString) where {T<:Unsigned}
  	# Header 12 bytes: 
	# -> version: 7 bytes string
	# -> size T: 1 byte
	# -> # vertices: 4 bytes position 
	#
	# for example: 'MGSv3  ' + <8bits T size> +  <32bits offset in size_t of data section>
  	version = 0x4d475376332020
	# size of type T in bytes
	size_t = convert(UInt8, sizeof(T))

	pos = T[]
	children = T[]
	vs = vertices(g)
	cpos = convert(T,1)
	for v in vs
		ovs = outneighbors(g,v)
		push!(pos,cpos)
		for o in ovs
			push!(children,o)	
			cpos += convert(T,1)
		end
	end

	# number of vertices
	gs = convert(UInt32, length(vs))

	f = open("$filename.mgs", "w")
	### write header
	# NB: reinterpret generates an array of length 8 even if version has a length of 7 bytes
	bytes = reinterpret(UInt8, [version])[1:7]
	write(f, bytes)
	bytes = reinterpret(UInt8, [size_t])
	write(f, bytes)
	bytes = reinterpret(UInt8, [gs])
	write(f, bytes)
	### write index section
	for p in pos
		bytes = reinterpret(UInt8, [p])
		write(f, bytes)
	end
	### write data section
	for c in children
		bytes = reinterpret(UInt8, [c])
		write(f, bytes)
	end
	close(f)
end

""" 
    load_mgs3_graph(g::AbstractGraph{T},filename::AbstractString) where {T<:Unsigned}

load graph in format MGS v3
"""
function load_mgs3_graph(g::AbstractGraph{T},filename::AbstractString) where {T<:Unsigned}
	f = open(filename, "r")
	### read header
	# 7-bytes version
	version = read(f,7)
	# size in bytes of T
	size_t = convert(UInt8, read(f,1)[1])
	# number of vertices
	gs = reinterpret(UInt32, read(f,4))[1]
	# read index
	pos = T[]
	for i in 1:gs
		p = read(f, sizeof(T))
		push!(pos, reinterpret(T,p)[1])
	end
	# read data
	children = T[]
	while !eof(f)
		c = read(f, sizeof(T))
		push!(children, reinterpret(T,c)[1])
	end
	close(f)
	
	# vertex set
	vs = range(1, stop=length(pos))

	# add vertices
    add_vertices!(g, length(vs))

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
    write_mgs4_graph(g::AbstractGraph{T},rg::AbstractGraph{T},filename::AbstractString) where {T<:Unsigned}

load graph in format MGS v4
"""
function write_mgs4_graph(g::AbstractGraph{T},rg::AbstractGraph{T},filename::AbstractString) where {T<:Unsigned}
  	# 8 bytes: 7 bytes string + 1 byte: 'MGSv4  ' + <8bits T size>
	# 4 bytes (32bits): offset in size_t of graph data section (i.e. # of vertices)
  	version = 0x4d475376342020
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
  	# 8 bytes: 7 bytes string + 1 byte: 'MGSv4  ' + <8bits T size>
	# 4 bytes (32bits): offset in size_t of graph data section (i.e. # of vertices)
	#
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
