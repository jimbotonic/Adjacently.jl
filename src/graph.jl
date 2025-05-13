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

module Graph

using LightGraphs, DataStructures, SparseArrays, LinearAlgebra, Logging
using ..CustomTypes: UInt24, UInt40
using ..CustomLightGraphs: SimpleDiGraph, SimpleGraph, SimpleEdge
using ..Util: infer_uint_custom_type
using ..Algo: pearce_iterative
using ..RandomWalks: RW, RW_aggregated

# Export the functions we want to make available
export get_basic_stats, 
       display_basic_stats, 
       get_out_degree_stats, 
       get_sinks, 
       get_sources, 
       subgraph, 
       subgraph_streamed, 
       get_core, 
       get_core_streamed, 
       get_reverse_graph, 
       get_vertex_in_degrees, 
       get_in_out_degrees, 
       get_avg_out_degree, 
       get_forward_ball, 
       get_clustering_coefficients, 
       get_colink_coefficients, 
       get_inclist_from_adjlist, 
       get_sparse_adj_matrix, 
       get_sparse_P_matrix

###
# Basic stats for directed graphs
###

""" 
    display_basic_stats(g,rg)

display basic stats
"""
function display_basic_stats(g::AbstractGraph,rg::AbstractGraph)
	nvs,nes,dens = get_basic_stats(g)
	nvs2,nes2,dens2 = get_basic_stats(rg)

	avg_od,max_od,sinks = get_out_degree_stats(g)
	avg_id,max_id,sources = get_out_degree_stats(rg)

	info("# vertices: $nvs")
	info("# edges: $nes")
	info("density: $dens")
	info()
	info("avg out-degree: $avg_od")
	info("max out-degree: $max_od")
	nsinks = length(sinks)
	info("# sinks: $nsinks")
	info()
	info("avg in-degree: $avg_id")
	info("max in-degree: $max_id")
	nsources = length(sources)
	info("# sources: $nsources")
end

"""
    get_basic_stats(g::AbstractGraph{T}) where {T<:Unsigned}

get # vertices, # of edges and density
NB: for computing density, we assume a directed graph
"""
function get_basic_stats(g::AbstractGraph{T}) where {T<:Unsigned}
	nvs = nv(g)
	nes = ne(g)
	density = nes/(nvs*(nvs-1))
	return nvs,nes,density
end

"""
    get_out_degree_stats(g::AbstractGraph{T}) where {T<:Unsigned}

get out degree stats

retuns avg out-degree, max out-degree, array of sinks vertices
"""
function get_out_degree_stats(g::AbstractGraph{T}) where {T<:Unsigned}
	sum = 0
	sinks = T[]
	max_degree = 0
	vs = vertices(g)
	for v in vs
		children = outneighbors(g,v)
		od = length(children)
		sum += od
		if od == 0
			push!(sinks,v)
		end
		if od > max_degree
			max_degree = od
		end
	end
	return sum/length(vs),max_degree,sinks
end

""" 
    get_sinks(g::AbstractGraph{T}) where {T<:Unsigned}

get array of sink vertices in the specified graph
"""
function get_sinks(g::AbstractGraph{T}) where {T<:Unsigned}
	sinks = T[]
	vs = vertices(g)
	for v in vs
		children = outneighbors(g,v)
		od = length(children)
		if od == 0
			push!(sinks,v)
		end
	end
	return sinks
end

""" 
    get_sources(g::AbstractGraph{T}) where {T<:Unsigned}

get array of source vertices in the specified graph
"""
function get_sources(g::AbstractGraph{T}) where {T<:Unsigned}
	achildren = T[]
	vs = vertices(g)
	for v in vs
		children = outneighbors(g,v)
		achildren = union(achildren,children)
	end
	return setdiff(vs,achildren)
end

###
# Subgraphs && SCCs
###

""" 
    subgraph(g::AbstractGraph{T},sids::Array{T,1}) where {T<:Unsigned}

get the subgraph of g induced by the set of vertices sids
"""
function subgraph(g::AbstractGraph{T},sids::Array{T,1}) where {T<:Unsigned}
	if typeof(g) == SimpleDiGraph{T}
		ng = SimpleDiGraph{T}()
	else
		ng = SimpleGraph{T}()
	end
	# nvs should be sorted in ascending order
	nvs = sort(sids)

	# add vertices
	add_vertices!(ng,length(nvs))

	# old id -> new id
	oni = Dict{T,T}()
	noi = Dict{T,T}()

	counter = convert(T,1)
	for v in nvs
		oni[v] = counter
		noi[counter] = v
		counter += convert(T,1)
	end

	# add edges
	for v in nvs
		children = outneighbors(g,v)
		source = oni[v]
		for c in children
			#if index_sorted(nvs,c) != 0
			if haskey(oni,c)
				target = oni[c]
				add_edge!(ng,source,target)
			end
		end
	end
	return ng,oni,noi
end

""" 
    subgraph_streamed(g::AbstractGraph{T},sids::Array{T,1},name::String) where {T<:Unsigned}

get the subgraph of g induced by the set of vertices sids
write the subgraph in MGS format v2
"""
function subgraph_streamed(g::AbstractGraph{T},sids::Array{T,1},name::String) where {T<:Unsigned}
	ng = SimpleDiGraph{T}()
	# nvs should be sorted in acending order
	nvs = sort(sids)

	f1 = open("$name.index", "w")
	f2 = open("$name.data", "w")

	# old id -> new id
	oni = Dict{T,T}()

	counter = convert(T,1)
	for v in nvs
		oni[v] = counter
		counter += convert(T,1)
	end

	pos = convert(T,1)
	for v in nvs
		children = outneighbors(g,v)
		bytes = reinterpret(Uint8, [pos])
                write(f1, reverse(bytes))
		for c in children
			if haskey(oni,c)
				target = oni[c]
				bytes = reinterpret(Uint8, [target])
                		write(f2, reverse(bytes))
				# pos += convert(Uint32, 1) -> 8 bytes!?
				pos = convert(T,pos+1)
			end
		end
	end

	@debug("# vertices: ", length(nvs))
	@debug("# edges: ", pos)

	close(f1)
	close(f2)

	return oni
end

""" 
    get_core(g::AbstractGraph{T}) where {T<:Unsigned}

get the main SCC

@returns ng (core subgraph), oni (old->new vertex indices), noi (new->old vertex indices)
"""
function get_core(g::AbstractGraph{T}) where {T<:Unsigned}
	sccs = pearce_iterative(g)
	scc_ids = union(sccs,[])
	id_size = Dict{T,T}()

	# compute the size of each SCC
	for id in scc_ids
		id_size[id] = 0
	end

	for id in sccs
		id_size[id] += 1
	end

	sizes = collect(values(id_size))
	sort!(sizes)
	msize = sizes[end]

	# find the max SCC id
	mid = 0
	for id in keys(id_size)
		id_size[id] == msize && begin mid = id; break end
	end

	sids = T[]
	counter = convert(T,1)
	for id in sccs
		id == mid && push!(sids,counter)
		counter += convert(T,1)
	end

	@debug("# core vids: ", length(sids))
	return subgraph(g,sids)
end

""" 
    get_core_streamed(g::AbstractGraph{T},sccs::Array{T,1},name::String) where {T<:Unsigned}

get the main SCC and write it to the specified file (MGS format)
write the subgraph in MGS format v2
"""
function get_core_streamed(g::AbstractGraph{T},sccs::Array{T,1},name::String) where {T<:Unsigned}
	scc_ids = union(sccs,[])
	id_size = Dict{T,T}()

	for id in scc_ids
		id_size[id] = 0
	end

	for id in sccs
		id_size[id] += 1
	end

	sizes = collect(values(id_size))
	sort!(sizes)
	msize = sizes[end]

	# find the max SCC id
	mid = 0
	for id in keys(id_size)
		id_size[id] == msize && begin mid = id; break end
	end

	sids = T[]
	counter = convert(T,1)
	for id in sccs
		id == mid && push!(sids,counter)
		counter = convert(T,counter+1)
	end

	@debug("# core vids: ", length(sids))
	subgraph_streamed(g,sids,name)
end

""" 
    get_reverse_graph(g::AbstractGraph{T}) where {T<:Unsigned}

get the reverse graph (same graph with all edge directions reversed)

@returns reversed graph
"""
function get_reverse_graph(g::AbstractGraph{T}) where {T<:Unsigned}
	ng = SimpleDiGraph{T}()

	# same set of vertices
	add_vertices!(ng,nv(g))

	# inverse the edge directions
	for v in vertices(g)
		children = outneighbors(g,v)
		for c in children
			add_edge!(ng,c,v)
		end
	end

	return ng
end

""" 
    get_vertex_in_degrees(g::AbstractGraph{T}) where {T<:Unsigned}

get in-degree of g vertices

@returns dictionary (vertex_id -> in-degree)
"""
function get_vertex_in_degrees(g::AbstractGraph{T}) where {T<:Unsigned}
	vin = Dict{T,T}()
	for v in vertices(g)
		ovs = outneighbors(g,v)
		for o in ovs
			if !haskey(vin,o)
				vin[o] = convert(T,1)
			else
				vin[o] = convert(T,vin[o]+1)
			end
		end
	end
	return vin
end

""" 
    get_in_out_degrees(g::AbstractGraph{T}) where {T<:Unsigned}

get in- and out- degree arrays of specified graph
"""
function get_in_out_degrees(g::AbstractGraph{T}) where {T<:Unsigned}
	vin = get_vertex_in_degrees(g)
	in_degrees = T[]
	out_degrees = T[]
	for v in vertices(g)
		push!(in_degrees,vin[v])
		push!(out_degrees,convert(T,length(outneighbors(g,v))))
	end
	return in_degrees,out_degrees
end

""" 
    get_avg_out_degree(g::AbstractGraph{T},visited::Array{T,1},p_avg::Float64=float64(-1),np_steps::UInt64=uint64(0)) where {T<:Unsigned}

get the avg out-degree of the specified set of visited nodes

p_avg: current average
nb_steps: number of points used to compute p_avg

NB: to get the avg in-degree of visited nodes, one can use the reverse graph of g
"""
function get_avg_out_degree(g::AbstractGraph{T},visited::Array{T,1},p_avg::Float64=float64(-1),np_steps::UInt64=uint64(0)) where {T<:Unsigned}
	sum = 0.
	for v in visited
		sum += length(outneighbors(g,v))
	end
	if p_avg != -1
		return (p_avg*np_steps+sum)/(np_steps+length(visited))
	else
		return sum/length(visited)
	end
end

""" 
    get_forward_ball(v::T,g::AbstractGraph{T},radius::Int=2,p::Float64=1) where {T<:Unsigned}

get a forward ball centered at the specified vertex
"""
function get_forward_ball(v::T,g::AbstractGraph{T},radius::Int=2,p::Float64=1) where {T<:Unsigned}
	# vertex ids of the ball
	subids = T[]
	push!(subids,v)
	# set of explored vertices
	explored = T[]
	for t in 1:radius
		so = T[]
		for u in subids
			if !(u in explored)
				append!(so,outneighbors(g,u))
				push!(explored,u)
			end
		end
		if length(so) > 0
			#  if the probability is 1, add the child
			if p == 1
				append!(subids,unique(so))
			else
				# we are at the first level of the BFS
				if t == 1
					append!(subids,unique(so))
				else
					so = unique(so)
					#added = false
					# select children at random with an exponentially decreasing probability
					for s in so
						if rand() <= p^(t-1)
							push!(subids,s)
							#added = true
						end
					end
					# if no child was added, add at least one at random
					#!added && push!(subids,so[rand(1:length(so))])
				end
			end
		end
	end
	return unique(subids)
end

""" 
    get_clustering_coefficients(g::AbstractGraph{T},rg::AbstractGraph{T},ntriangles::Array{T,1},density::Float64=-1.) where {T<:Unsigned}

get the array of clustering coefficients

ncolinks: array of the number of colinks for each vertex
"""
function get_clustering_coefficients(g::AbstractGraph{T},rg::AbstractGraph{T},ntriangles::Array{T,1},density::Float64=-1.) where {T<:Unsigned}
	vs = vertices(g)
	n = nv(g)
	ccs = zeros(Float64,n)
	for v in vs
		parents = outneighbors(rg,v)
		children = outneighbors(g,v)
		p = length(parents)
		c = length(children)
		i = length(intersect(parents,children))
		# maximum number of triangles:
		# = (p_minus_c * c_minus_p) + i * (p_minus_c + c_minus_p + (i - 1))
		# = p*c - i
		npt = p*c-i
		# if there is only 1 colink, no triangle can be formed
		if npt == 0
			ccs[v] = 0
		else
			# correct the result according to the specified density
			if density > 0
				ccs[v] = (ntriangles[v] - density*npt)/(npt*(1-density))
			else
				ccs[v] = ntriangles[v]/npt
			end
		end
	end
	return ccs
end

""" 
    get_colink_coefficients(g::AbstractGraph{T},rg::AbstractGraph{T},ncolinks::Array{T,1},density::Float64=-1.) where {T<:Unsigned}

get the array of colink coefficients

ncolinks: array of the number of colinks for each vertex
"""
function get_colink_coefficients(g::AbstractGraph{T},rg::AbstractGraph{T},ncolinks::Array{T,1},density::Float64=-1.) where {T<:Unsigned}
	vs = vertices(g)
	n = nv(g)
	ccs = zeros(Float64,n)
	for v in vs
		parents = outneighbors(rg,v)
		children = outneighbors(g,v)
		p = length(parents)
		c = length(children)
		i = length(intersect(parents,children))
		# maximum number of colinks:
		# = p+c-i
		npc = p+c-i
		# correct the result according to the specified density
		if density > 0
			ccs[v] = (ncolinks[v] - density*npc)/(npc*(1-density))
		else
			ccs[v] = ncolinks[v]/npc
		end
	end
	return ccs
end

""" 
    get_sparse_adj_matrix(g::AbstractGraph{T}) where {T<:Unsigned}

get sparse adjacency matrix A
"""
function get_sparse_adj_matrix(g::AbstractGraph{T}) where {T<:Unsigned}
	I = Array{T,1}()
	J = Array{T,1}()
	V = Array{Float64,1}()
	for u in vertices(g)
		for v in outneighbors(g,u)
			push!(I,u)
			push!(J,v)
			push!(V,1.)
		end
	end
	return sparse(I,J,V)
end

""" 
    get_adj_matrix(g::AbstractGraph{T}) where {T<:Unsigned}

get adjacency matrix A
"""
function get_adj_matrix(g::AbstractGraph{T}) where {T<:Unsigned}
	n = nv(g) 
	A = zeros(Float64,n,n) 
	for u in vertices(g)
		for v in outneighbors(g,u)
			A[u,v] = 1.
		end
	end
	return A
end

""" 
    get_sparse_P_matrix(g::AbstractGraph{T}) where {T<:Unsigned}

get P = D^-1 * A matrix
"""
function get_sparse_P_matrix(g::AbstractGraph{T}) where {T<:Unsigned}
	n = nv(g)
	A = get_sparse_adj_matrix(g)
	S = dropdims(sum(A, dims=2), dims=2)
	S = 1 ./ S
	range = convert(Array{T,1}, 1:n)
	D = sparse(range, range, S)
	return D*A
end

""" 
    get_sparse_symmetric_P_matrix(g::AbstractGraph{T}) where {T<:Unsigned}

get P = D*^(-1/2) A* D*^(-1/2) matrix with A* = (A+I)

NB: we assume there is no sink in the graph
"""
function get_sparse_symmetric_P_matrix(g::AbstractGraph{T}) where {T<:Unsigned}
	n = length(G.vertices(g))
	A = get_sparse_adj_matrix(g)
	A = A + sparse(I,n,n)
	S = dropdims(sum(A, dims=2), dims=2)
	S = 1 ./ sqrt.(S)
	range = convert(Array{T,1}, 1:n)
	D = sparse(range, range, S)
	return D*A*D
end

end # module Graph

