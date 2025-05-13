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

module IO

using LightGraphs, DataStructures, HDF5, JLD
using ..CustomTypes: UInt24, UInt40
using ..NodeTypes: Node, EmptyNode
using ..CustomLightGraphs: SimpleDiGraph, SimpleGraph, SimpleEdge
using ..Util: infer_uint_custom_type

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
    load_adjacency_list_from_csv(filename::AbstractString,separator::Char=',')

load graph from CSV adjacency list
"""
function load_adjacency_list_from_csv(filename::AbstractString, separator::AbstractChar=',')
	f = open(filename,"r")
	oni = Dict{UInt64,UInt64}()
	edges = Array{Tuple{UInt64,UInt64},1}()
	counter = convert(UInt64,1)
	while !eof(f)
		line = strip(readline(f))
		if !startswith(line, "#")
			edge = split(line, separator)
			v1 = parse(UInt64,edge[1])
			v2 = parse(UInt64,edge[2])

			if !haskey(oni, v1)
				oni[v1] = counter
				counter += convert(UInt64,1)
			end
			if !haskey(oni, v2)
				oni[v2] = counter
				counter += convert(UInt64,1)
			end
			push!(edges, (oni[v1], oni[v2]))
		end
	end
	close(f)

	gs = length(keys(oni))
	nbits = convert(UInt8, ceil(log(2, gs)))
	V = infer_uint_custom_type(nbits)

	g = SimpleDiGraph{V}()

	# add vertices
	add_vertices!(g, gs)
	
	# add edges
	for edge in edges
		add_edge!(g, convert(V, edge[1]), convert(V, edge[2]))	
	end

	return g::AbstractGraph{V}
end

""" 
    load_adjacency_list(adj_list::Array{Unsigned,2})

load graph from adjacency list

NB: adjcency list should be represented as a 2-dimensional array with 2 rows and 1 column per edge
"""
function load_adjacency_list(adj_list::Array{Unsigned,2})
	oni = Dict{UInt64,UInt64}()
	edges = Array{Tuple{UInt64,UInt64},1}()
	counter = convert(UInt64,1)
	nes = size(adj_list)[2]

	for i in 1:nes
		edge = adj_list[i]
		v1 = edge[1]
		v2 = edge[2]

		if !haskey(oni, v1)
			oni[v1] = counter
			counter += convert(UInt64,1)
		end
		if !haskey(oni, v2)
			oni[v2] = counter
			counter += convert(UInt64,1)
		end
		push!(edges, (oni[v1], oni[v2]))
	end

	gs = length(keys(oni))
	nbits = convert(UInt8, ceil(log(2, gs)))
	V = infer_uint_custom_type(nbits)

	g = SimpleDiGraph{V}()
	# add vertices
	add_vertices!(g,gs)
	
	# add edges
	for edge in edges
		add_edge!(g, convert(V, edge[1]), convert(V, edge[2]))	
	end

	return g::AbstractGraph{V}
end

""" 
    load_graph_from_pajek(filename::AbstractString)

Load net Pajek file
"""
function load_graph_from_pajek(filename::AbstractString)
	f = open(filename,"r")
	
	inside_vertices_section = false
	inside_edges_section = false
	
	vdict = Dict{UInt64, UInt64}()
	vcounter = convert(UInt64,1)
	edges = Array{Tuple{UInt64,UInt64},1}()

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
			# Example of vertices section:
			#
			# *Vertices 82670
			# 1 "entity"
			# 2 "thing"
			# 3 "anything"
			# 4 "something"
			# 5 "nothing"
			# 6 "whole"
			if inside_vertices_section
				sa = split(line, ' ')
				# Dictionary of vertices ids and their associated counter
				vdict[parse(UInt64, sa[1])] = vcounter
				vcounter += convert(UInt64, 1)	
			# Example of edges section:
			#
			# *Arcs
			# 1 2
			# 1 3
			# 1 4
			# 1 5
			elseif inside_edges_section
				sa = split(line, ' ')
				vs = vdict[parse(UInt64, sa[1])]
				vt = vdict[parse(UInt64, sa[2])]
				push!(edges, (vs, vt))
			end
		end
	end
	close(f)

	gs = length(keys(vdict))
	nbits = convert(UInt8, ceil(log(2, gs)))
	V = infer_uint_custom_type(nbits)

	g = SimpleDiGraph{V}()
	add_vertices!(g, gs)
	for edge in edges
		add_edge!(g, convert(V, edge[1]), convert(V, edge[2]))
	end

	return g::AbstractGraph{V}
end

""" 
    load_triangles(::Type{T},filename::AbstractString) where {T<:Unsigned}

load triangles list from CSV text-formatted file 

TODO: infer type T 
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

end
