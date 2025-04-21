# MGS File Format

MGS stands for Massive Graph Storage. It is a binary file format designed to store the toplogy of massive directed graphs in a compact form. 

In the following sections, we describe the MGS file format in detail. $n$ represents the number of vertices in the graph.

## MGS version 3

MGS version 3 stores graphs in the form of adjacency lists in a single file {graph-name}.mgs.

The graph vertices are encoded as unsigned integers over a variable size of bytes determined by the number $n$ of vertices to be represented (see below). 

The file header is 12 bytes long. It is composed of the following fields:

* File format identifier: 'MGS' string (3 bytes)
* Version number: major 0x03 + minor 0x00 (2 bytes)
* Flags (2 bytes)
  * Graph type  (4 bits) + compression scheme (4 bits)
  * Coding scheme (4 bits) + reserved flags (4 bits)
* Number of vertices (5 bytes)

The possible graph types are:
* 0x0: directed graph
* 0x1: undirected graph

The possible compression algorithms are:
* 0x0: no compression
* 0x1: Huffman encoding of vertex ids

The possible coding schemes are:
* 0x0: data section only with implicit numbering of vertices
* 0x1: index section + data section with explicit numbering of vertices

The minimum number of bits allocated for each vertex id is computed as follows:

$$ nbits = \lceil log_2 (n+1) \rceil$$

The number of bytes allocated for each vertex id is then given by the closest unsigned integer type that can represent the vertex ids among the following:

* `UInt8` for $n \leq 255$
* `UInt16` for $n \leq 65535$
* `UInt24` for $n \leq 16777215$
* `UInt32` for $n \leq 4294967295$
* `UInt40` for $n \leq 1099511627775$

### 0x0 Coding Scheme (data section only)

There is no index section. The data section contains contiguous lists of vertex ids representing the children of each vertex numbered implicitly from $1$ to $n. Each adjacency list is terminated by a special value $0$ of $size\_t$ bytes. 

NB: by reordering the vertices by decreasing out-degree allows, it is possible to spare adding $0$-termination values corresponding to $0$-out-degree vertices. 

### 0x1 Coding Scheme (index + data sections)

The index section starts directly after the header and has a size of $size\_t * n$ bytes. Each entry indicates the number of children of the corresponding vertex numbered implicitly from $1$ to $n$.

The following data section contains contiguous lists of vertex ids representing the children of each vertex.

## Deprecated MGS formats

### MGS version 1 (deprecated)

MGS version 1 stores graphs in the form of adjacency lists in two separate files: 
* The graph index is a set of contiguous (vertex id, start of adjacency list position) pairs.
* The graph data is a set of contiguous vertex ids representing the children of each vertex specified in the index.

NB:
* Vertex ids do not need to be contiguous. 
* Vertex ids and adjacency list positions are stored in big endian format at both bit and byte levels.
* Vertex ids and adjacency list positions are stored as unsigned integers of the specified size `size_t`.
* Adjacency list positions are 1-based.

The two MGS files have the same name with the following extensions: {graph-name}.index and {graph-name}.data.

```julia
""" 
    load_mgs1_graph_index(ipos::OrderedDict{T,T},filename::AbstractString) where {T<:Unsigned}

get the ordered dictionary vid -> startpos (MGS v1) (deprecated)
"""
function load_mgs1_graph_index(ipos::OrderedDict{T,T},filename::AbstractString) where {T<:Unsigned}
	f = open(filename, "r")
	while !eof(f1)
		# read `size_t` bytes
		id = read(f, UInt8, sizeof(T))
		pos = read(f, UInt8, sizeof(T))
		# Add value to the dictionary: ipos[vid] <- pos 
		#
		# NB: `ipos` positions should be 1-based encoded
		# NB: `reinterpet` takes as input an array of bytes in little endian
		# NB: `reinterpret` returns an array of T value -> select the first and only one
		ipos[reinterpret(T, reverse(id))[1]] = reinterpret(T, reverse(pos))[1]
	end
	close(f)
end

""" 
    load_graph_data(children::Array{T,1},filename::AbstractString) where {T<:Unsigned}

get the array of graph children (MGS v1 & MGS v2) (deprecated)
"""
function load_graph_data(children::Array{T,1},filename::AbstractString) where {T<:Unsigned}
	f = open(filename, "r")
	while !eof(f)
		child = read(f, UInt8, sizeof(T))
		# NB: 
		# - `reinterpet` takes as input an array of bytes in little endian
		# - `reinterpret` returns an array of T value -> select the first and only one
		push!(children, reinterpret(T, reverse(child))[1])
	end
	close(f)
end

""" 
    write_mgs1_graph_index(ipos::OrderedDict{T,T}, filename::AbstractString) where {T<:Unsigned}

write index file (MGS v1) (deprecated)
"""
function write_mgs1_graph_index(ipos::OrderedDict{T,T}, filename::AbstractString) where {T<:Unsigned}
	f = open(filename, "w")
	# NB: `ipos` positions should be 1-based encoded
	for p in ipos
    	# reinterpret pair (vid,pos) in an array of bytes
		#
		# See for example: 
		#
		# d = OrderedDic{UInt16,UInt16}()
		# d[1] = 2
		# d[3] = 4
		# for p in d
       	# bytes = reinterpret(UInt8, [p])
       	# println(bytes)
       	# end
		#
		# UInt8[0x01, 0x00, 0x02, 0x00]
		# UInt8[0x03, 0x00, 0x04, 0x00]
		bytes = reinterpret(UInt8, [p])
		# write bytes in big endian
		#
		# See for example:
		# c = UInt16(1)
		# bytes = reinterpret(UInt8, [c])
		# 2-element reinterpret(UInt8, ::Vector{UInt16}):
 		# 0x01
 		# 0x00
		mid_pos = convert(T, length(bytes)/2)
		end_pos = convert(T, lastindex(bytes))
		write(f, reverse(bytes[1:mid_pos]))
		write(f, reverse(bytes[mid_pos+convert(T,1):end_pos]))
	end
	close(f)
end

""" 
    write_graph_data(children::Array{T,1}, filename::AbstractString) where {T<:Unsigned}

write the array of graph children (MGS v1 & MGS v2) (deprecated)
"""
function write_graph_data(children::Array{T,1}, filename::AbstractString) where {T<:Unsigned}
	f = open(filename, "w")
	for c in children
		bytes = reinterpret(UInt8, [c])
		# write bytes in big endian
		#
		# See for example:
		# c = UInt16(1)
		# bytes = reinterpret(UInt8, [c])
		# 2-element reinterpret(UInt8, ::Vector{UInt16}):
 		# 0x01
 		# 0x00
		write(f, reverse(bytes))
	end
	close(f)
end
```

## MGS version 2 (deprecated)

MGS version 2 stores graphs in the form of adjacency lists in two separate files: 
* The graph index is a set of contiguous start of adjacency list positions.
* The graph data is a set of contiguous vertex ids representing the children of each vertex specified in the index.

NB:
* Vertex ids are numbered implicitly from 1 to $n$ in increasing order. 
* Adjacency list positions are stored in big endian format at both bit and byte levels.
* Adjacency list positions are stored as unsigned integers of the specified size `size_t`.
* Adjacency list positions are 1-based.
* In case a vertex has no child, $pos[k] = pos[k+1]$ in the index file. 

The two MGS files have the same name with the following extensions: {graph-name}.index and {graph-name}.data.

```julia
""" 
    load_mgs2_graph_index(pos::Array{T,1},filename::AbstractString) where {T<:Unsigned}

get the set of positions (MGS v2) (deprecated)
"""
function load_mgs2_graph_index(pos::Array{T,1},filename::AbstractString) where {T<:Unsigned}
	f = open(filename, "r")
	while !eof(f)
		# read `size_t` bytes
		p = read(f, UInt8, sizeof(T))
		# Add value to the list: p 
		#
		# NB: `reinterpet` takes as input an array of bytes in little endian
		# NB: `reinterpret` returns an array of T value -> select the first and only one
		push!(pos, reinterpret(T, reverse(p))[1])
	end
	close(f)
end

""" 
    load_mgs2_graph(g::AbstractGraph{T},name::AbstractString) where {T<:Unsigned}

load graph in format MGS v2 (deprecated)
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
    write_mgs2_graph_index(pos::Array{T,1}, filename::AbstractString) where {T<:Unsigned}

write index file (MGS v2) (deprecated)
"""
function write_mgs2_graph_index(pos::Array{T,1}, filename::AbstractString) where {T<:Unsigned}
	f = open(filename, "w")
	for p in pos
    	# reinterpret pos in an array of bytes
		bytes = reinterpret(UInt8, [p])
		# write bytes in big endian
		#
		# See for example:
		# c = UInt16(1)
		# bytes = reinterpret(UInt8, [c])
		# 2-element reinterpret(UInt8, ::Vector{UInt16}):
 		# 0x01
 		# 0x00
		write(f, reverse(bytes))
	end
	close(f)
end
```
