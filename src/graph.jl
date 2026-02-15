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
using ..PageRank: PR

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
       get_in_degrees,
       get_in_out_degrees,
       get_out_degrees,
	   get_avg_out_degree,
       get_forward_ball,
       get_clustering_coefficients,
       get_colink_coefficients,
       get_inclist_from_adjlist,
       get_sparse_adj_matrix,
       get_sparse_P_matrix,
       get_sparse_symmetric_P_matrix,
	   get_neighbor_lists,
	   # Local clustering functions
	   compute_volume,
	   compute_cut,
	   compute_conductance,
	   compute_modularity_gain,
	   compute_support_score,
	   select_children_by_support,
	   select_children_by_conductance,
	   select_children_by_modularity,
	   # Local PPR
	   push_ppr,
	   select_children_by_ppr

########################################################
# Basic stats for directed graphs
########################################################

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

@param g: the graph
@returns # vertices, # of edges, density
"""
function get_basic_stats(g::AbstractGraph{T}) where {T<:Unsigned}
	nvs = nv(g)
	nes = ne(g)
	if is_directed(g)
		density = nes/(nvs*(nvs-1))
	else
		density = 2*nes/(nvs*(nvs-1))
	end
	return nvs,nes,density
end

"""
    get_out_degree_stats(g::AbstractGraph{T}) where {T<:Unsigned}

get out-degree stats

@param g: the graph
@returns avg out-degree, max out-degree, array of sink vertices
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

@param g: the graph
@returns array of sink vertices
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

@param g: the graph
@returns array of source vertices
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

########################################################
# Subgraphs && SCCs
########################################################

""" 
    subgraph(g::AbstractGraph{T},sids::Array{T,1}) where {T<:Unsigned}

get the subgraph of g induced by the set of vertices sids

@param g: the graph
@param sids: the set of vertices
@returns the subgraph of g induced by the set of vertices sids
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

@param g: the graph
@param sids: the set of vertices
@param name: the name of the file
@returns the subgraph of g induced by the set of vertices sids
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

@param g: the graph
@param sccs: the set of SCCs
@param name: the name of the file
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

@param g: the graph
@param sccs: the set of SCCs
@param name: the name of the file
@returns the subgraph of g induced by the set of vertices sids
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

@param g: the graph
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

########################################################
# Degree stats
########################################################

""" 
    get_in_degrees(g::AbstractGraph{T}, V::Type{<:Unsigned}=T) where {T<:Unsigned}

get in-degree of g vertices

@param g: the graph
@returns dictionary (vertex_id -> in-degree)
"""
function get_in_degrees(g::AbstractGraph{T}) where {T<:Unsigned}
	# if g is a directed graph, compute the in-degrees
	if is_directed(g)
		in_degrees = Dict{T,T}()

		# initialize frequencies for all vertices to ensure every vertex has an entry
		for v in vertices(g)
			if !haskey(in_degrees, v)
				in_degrees[v] = zero(T)
			end
		end

		for v in vertices(g)
			for o in outneighbors(g, v)
				if !haskey(in_degrees, o)
					in_degrees[o] = one(T)
				else
					in_degrees[o] += one(T)
				end
			end
		end
		return in_degrees
	else
		return get_out_degrees(g)
	end
end

"""
    get_out_degrees(g::AbstractGraph{T}) where {T<:Unsigned}

get out-degree of g vertices

@param g: the graph
@returns dictionary (vertex_id -> out-degree)
"""
function get_out_degrees(g::AbstractGraph{T}) where {T<:Unsigned}
	out_degrees = Dict{T,T}()
	for v in vertices(g)
		out_degrees[v] = convert(T, length(outneighbors(g, v)))
	end
	return out_degrees
end

""" 
    get_in_out_degrees(g::AbstractGraph{T}) where {T<:Unsigned}

get in- and out- degree dictionaries of specified graph

@param g: the graph
@returns in-degree dictionary (vertex_id -> in-degree), out-degree dictionary (vertex_id -> out-degree)
"""
function get_in_out_degrees(g::AbstractGraph{T}) where {T<:Unsigned}
	if is_directed(g)
		in_degrees = get_in_degrees(g)
		out_degrees = get_out_degrees(g)
		return in_degrees, out_degrees
	else
		return get_out_degrees(g), get_out_degrees(g)
	end
end

""" 
    get_avg_out_degree(g::AbstractGraph{T},visited::Array{T,1},p_avg::Float64=float64(-1),np_steps::UInt64=uint64(0)) where {T<:Unsigned}

get the avg out-degree of the specified set of visited nodes

NB: to get the avg in-degree of visited nodes, one can use the reverse graph of g

@param g: the graph
@param visited: the set of visited nodes
@param p_avg: current average
@param np_steps: number of points used to compute p_avg

@returns the avg out-degree of the specified set of visited nodes
"""
function get_avg_out_degree(g::AbstractGraph{T}, visited::Array{T,1}, p_avg::Float64=-1., np_steps::UInt64=0) where {T<:Unsigned}
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
    get_neighbor_lists(g::AbstractGraph{T}) where {T<:Unsigned}

get the neighbor lists of the graph

@param g: the graph
@returns the neighbor lists as a dictionary (vertex_id -> list of neighbors)
"""
function get_neighbor_lists(g::AbstractGraph{T}) where {T<:Unsigned}
	neighbor_lists = Dict{T,Vector{T}}()
	for v in vertices(g)
		neighbor_lists[v] = outneighbors(g,v)
	end
	return neighbor_lists
end

""" 
    get_forward_ball(v::T,g::AbstractGraph{T},radius::Int=2,p::Float64=1) where {T<:Unsigned}

get a forward ball centered at the specified vertex

@param v: the vertex
@param g: the graph
@param radius: the radius of the ball
@param p: the probability of adding a child

@returns the set of vertex ids in the ball
"""
function get_forward_ball(v::T, g::AbstractGraph{T}, radius::Int=2, p::Float64=1.) where {T<:Unsigned}
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

@param g: the graph
@param rg: the reverse graph
@param ntriangles: the number of triangles for each vertex
@param density: the density of the graph

@returns the array of clustering coefficients
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

@param g: the graph
@param rg: the reverse graph
@param ncolinks: array of the number of colinks for each vertex
@param density: the density of the graph

@returns the array of colink coefficients
"""
function get_colink_coefficients(g::AbstractGraph{T}, rg::AbstractGraph{T}, ncolinks::Array{T,1}, density::Float64=-1.) where {T<:Unsigned}
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

@param g: the graph
@returns the sparse adjacency matrix
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

@param g: the graph
@returns the adjacency matrix
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
    get_sparse_P_matrix(g::AbstractGraph{T}; column_stochastic::Bool=false) where {T<:Unsigned}

Get column-stochastic transition matrix P where P[i,j] = 1/outdeg(j) for edge j→i.
This is the standard form used for PageRank computation.

For row-stochastic matrix (P = D^-1 * A), use column_stochastic=false.

@param g: the graph
@param column_stochastic: if true, return column-stochastic matrix (default: false for backward compatibility)
@returns the sparse P matrix (column-stochastic by default)

Note: Dangling nodes (outdegree=0) result in zero columns/rows, which is correct for PageRank.
"""
function get_sparse_P_matrix(g::AbstractGraph{T}; column_stochastic::Bool=true) where {T<:Unsigned}
	n = nv(g)
	I_idx = Int[]
	J_idx = Int[]
	V_vals = Float64[]

	if column_stochastic
		# Build column-stochastic matrix: P[i,j] = 1/outdeg(j) for edge j→i
		for j in 1:n
			jv = T(j)
			outs = outneighbors(g, jv)
			d = length(outs)
			if d == 0
				continue  # Skip dangling nodes - zero column is correct
			end
			w = 1.0 / d
			for i in outs
				push!(I_idx, Int(i))
				push!(J_idx, j)
				push!(V_vals, w)
			end
		end
	else
		# Build row-stochastic matrix: P[i,j] = 1/outdeg(i) for edge i→j
		for i in 1:n
			iv = T(i)
			outs = outneighbors(g, iv)
			d = length(outs)
			if d == 0
				continue  # Skip dangling nodes - zero row is correct
			end
			w = 1.0 / d
			for j in outs
				push!(I_idx, i)
				push!(J_idx, Int(j))
				push!(V_vals, w)
			end
		end
	end

	Pint = sparse(I_idx, J_idx, V_vals, n, n)
	return SparseMatrixCSC{Float64, T}(Pint)
end

""" 
    get_sparse_symmetric_P_matrix(g::AbstractGraph{T}) where {T<:Unsigned}

get P = D*^(-1/2) A* D*^(-1/2) matrix with A* = (A+I)

NB: we assume there is no sink in the graph

@param g: the graph
@returns the sparse P matrix
"""
function get_sparse_symmetric_P_matrix(g::AbstractGraph{T}) where {T<:Unsigned}
	n = length(vertices(g))
	A = get_sparse_adj_matrix(g)
	A = A + sparse(I,n,n)
	S = dropdims(sum(A, dims=2), dims=2)
	S = 1 ./ sqrt.(S)
	range = convert(Array{T,1}, 1:n)
	D = sparse(range, range, S)
	return D*A*D
end

########################################################
# LOCAL CLUSTERING / COMMUNITY DETECTION FUNCTIONS
########################################################

"""
    compute_volume(S::Set{T}, graph) where T<:Unsigned

Compute vol(S) = sum of degrees of vertices in S.
For directed graphs, this is the sum of out-degrees.
"""
function compute_volume(S::Set{T}, graph) where T<:Unsigned
    vol = 0
    for v in S
        vol += length(outneighbors(graph, v))
    end
    return vol
end

"""
    compute_cut(S::Set{T}, graph) where T<:Unsigned

Compute cut(S) = number of edges from S to V\\S.
"""
function compute_cut(S::Set{T}, graph) where T<:Unsigned
    cut = 0
    for v in S
        for neighbor in outneighbors(graph, v)
            if !(neighbor in S)
                cut += 1
            end
        end
    end
    return cut
end

"""
    compute_conductance(S::Set{T}, graph, total_volume::Int) where T<:Unsigned

Compute conductance φ(S) = cut(S) / min(vol(S), vol(V\\S)).
"""
function compute_conductance(S::Set{T}, graph, total_volume::Int) where T<:Unsigned
    if isempty(S)
        return Inf
    end

    vol_S = compute_volume(S, graph)
    cut_S = compute_cut(S, graph)

    if vol_S == 0
        return Inf
    end

    vol_complement = total_volume - vol_S
    min_vol = min(vol_S, vol_complement)

    if min_vol == 0
        return Inf
    end

    return cut_S / min_vol
end

"""
    compute_modularity_gain(u::T, S::Set{T}, graph, m::Int, vol_S::Int, k_in::Int) where T<:Unsigned

Compute modularity gain ΔQ from adding vertex u to community S.
ΔQ ≈ (1/2m) * (2*k_in(u) - deg(u)*vol(S)/m)

Parameters:
- u: candidate vertex
- S: current community
- graph: the graph
- m: total number of edges
- vol_S: current volume of S
- k_in: number of edges from u into S
"""
function compute_modularity_gain(u::T, S::Set{T}, graph, m::Int, vol_S::Int, k_in::Int) where T<:Unsigned
    deg_u = length(outneighbors(graph, u))
    return (2 * k_in - (deg_u * vol_S) / m) / (2 * m)
end

"""
    compute_support_score(k_in::Int, deg::Int, lambda::Float64=0.5)

Compute support score: k_in - λ*(deg - k_in)
Higher score means stronger connection to community.

Parameters:
- k_in: edges into community
- deg: total degree
- lambda: trade-off parameter (0 = only k_in matters, 1 = balanced)
"""
function compute_support_score(k_in::Int, deg::Int, lambda::Float64=0.5)
    return k_in - lambda * (deg - k_in)
end

"""
    select_children_by_support(level_vertices::Set{T}, graph,
                                rho::Float64=0.4, tau::Int=2) where T<:Unsigned

Select children using support/tie-strength threshold (Method #1 - fastest).
Returns Set of selected children.

Parameters:
- level_vertices: vertices in current level
- graph: the graph
- rho: relative support threshold (fraction of degree)
- tau: absolute support threshold (minimum edges into S)
"""
function select_children_by_support(level_vertices::Set{T}, graph,
                                     rho::Float64=0.4, tau::Int=2) where T<:Unsigned
    # Get all potential children
    all_children = Set{T}()
    for v in level_vertices
        for neighbor in outneighbors(graph, v)
            if !(neighbor in level_vertices)
                push!(all_children, neighbor)
            end
        end
    end

    # Compute support for each child
    selected = Set{T}()
    for child in all_children
        k_in = 0  # edges from child into level_vertices
        for neighbor in outneighbors(graph, child)
            if neighbor in level_vertices
                k_in += 1
            end
        end

        # Also count incoming edges (for directed graphs)
        for v in level_vertices
            for neighbor in outneighbors(graph, v)
                if neighbor == child
                    k_in += 1
                end
            end
        end

        deg = length(outneighbors(graph, child))

        # Accept if passes thresholds
        if k_in >= tau || (deg > 0 && k_in / deg >= rho)
            push!(selected, child)
        end
    end

    return selected
end

"""
    select_children_by_conductance(level_vertices::Set{T}, graph,
                                    max_select::Int=typemax(Int)) where T<:Unsigned

Select children that improve conductance (Method #2).
Returns Set of selected children.

Parameters:
- level_vertices: vertices in current level
- graph: the graph
- max_select: maximum number of children to select
"""
function select_children_by_conductance(level_vertices::Set{T}, graph,
                                        max_select::Int=typemax(Int)) where T<:Unsigned
    S = copy(level_vertices)
    vol_S = compute_volume(S, graph)
    cut_S = compute_cut(S, graph)

    # Get all potential children with their k_in values
    frontier = Dict{T, Int}()  # child => k_in
    for v in level_vertices
        for neighbor in outneighbors(graph, v)
            if !(neighbor in S)
                frontier[neighbor] = get(frontier, neighbor, 0) + 1
            end
        end
    end

    selected = Set{T}()
    num_selected = 0

    # Greedily add children that improve conductance
    while !isempty(frontier) && num_selected < max_select
        best_child = nothing
        best_improvement = 0.0
        current_conductance = vol_S > 0 ? cut_S / vol_S : Inf

        for (child, k_in) in frontier
            deg_child = length(outneighbors(graph, child))
            k_out = deg_child - k_in

            # Compute conductance after adding child
            new_vol = vol_S + deg_child
            new_cut = cut_S - k_in + k_out
            new_conductance = new_vol > 0 ? new_cut / new_vol : Inf

            improvement = current_conductance - new_conductance
            if improvement > best_improvement
                best_improvement = improvement
                best_child = child
            end
        end

        # Stop if no improvement
        if best_improvement <= 0
            break
        end

        # Add best child
        push!(selected, best_child)
        push!(S, best_child)
        num_selected += 1

        # Update vol and cut
        deg_best = length(outneighbors(graph, best_child))
        k_in_best = frontier[best_child]
        vol_S += deg_best
        cut_S = cut_S - k_in_best + (deg_best - k_in_best)

        # Remove from frontier and update neighbors
        delete!(frontier, best_child)
        for neighbor in outneighbors(graph, best_child)
            if !(neighbor in S)
                frontier[neighbor] = get(frontier, neighbor, 0) + 1
            end
        end
    end

    return selected
end

"""
    select_children_by_modularity(level_vertices::Set{T}, graph, m::Int,
                                   max_select::Int=typemax(Int)) where T<:Unsigned

Select children using modularity gain (Method #3).
Returns Set of selected children.

Parameters:
- level_vertices: vertices in current level
- graph: the graph
- m: total number of edges in graph
- max_select: maximum number of children to select
"""
function select_children_by_modularity(level_vertices::Set{T}, graph, m::Int,
                                       max_select::Int=typemax(Int)) where T<:Unsigned
    S = copy(level_vertices)
    vol_S = compute_volume(S, graph)

    # Get all potential children with their k_in values
    frontier = Dict{T, Int}()
    for v in level_vertices
        for neighbor in outneighbors(graph, v)
            if !(neighbor in S)
                frontier[neighbor] = get(frontier, neighbor, 0) + 1
            end
        end
    end

    selected = Set{T}()
    num_selected = 0

    # Greedily add children with positive modularity gain
    candidates = collect(frontier)
    sort!(candidates, by=x->begin
        deg = length(outneighbors(graph, x[1]))
        compute_modularity_gain(x[1], S, graph, m, vol_S, x[2])
    end, rev=true)

    for (child, k_in) in candidates
        if num_selected >= max_select
            break
        end

        deg_child = length(outneighbors(graph, child))
        delta_Q = compute_modularity_gain(child, S, graph, m, vol_S, k_in)

        if delta_Q > 0
            push!(selected, child)
            push!(S, child)
            vol_S += deg_child
            num_selected += 1
        else
            break  # Sorted by gain, so no more positive gains
        end
    end

    return selected
end

########################################################
# PUSH-BASED LOCAL PERSONALIZED PAGERANK (PPR)
########################################################

"""
    push_ppr(seeds::Union{T,Vector{T},Set{T}}, graph, reverse_graph;
             alpha::Float64=0.15, epsilon::Float64=1e-6,
             max_pushes::Int=typemax(Int)) where T<:Unsigned

Compute local Personalized PageRank using push-based algorithm.

This is the state-of-the-art for local community detection. It only explores
a small subgraph around the seed vertices, making it much faster than global
PageRank computation.

Algorithm:
1. Start with residual mass at seed vertices
2. Push mass from high-residual vertices to neighbors
3. Stop when all residuals are below epsilon

The push operation:
- For vertex u with residual r[u]:
  - Add α*r[u] to PPR estimate p[u]
  - Distribute (1-α)*r[u] equally to out-neighbors
  - Set r[u] = 0

Returns:
- p: PPR vector (Dict{T, Float64}) - only non-zero entries
- support: Set of vertices explored (subgraph touched by algorithm)

Parameters:
- seeds: seed vertex/vertices (can be single vertex, Vector, or Set)
- graph: the forward graph
- reverse_graph: the reverse graph (for directed graphs)
- alpha: teleport probability (default 0.15, standard PageRank damping = 1-alpha)
- epsilon: residual threshold (default 1e-6)
- max_pushes: maximum number of push operations (prevents infinite loops)

Reference:
- Andersen, Chung, Lang. "Local Graph Partitioning using PageRank Vectors" (FOCS 2006)
- http://www.cs.cmu.edu/afs/cs/user/glmiller/public/Scientific-Computing/F-11/RelatedWork/local-partitioning-by-pagerank.pdf
"""
function push_ppr(seeds::Union{T,Vector{T},Set{T}}, graph, reverse_graph;
                  alpha::Float64=0.15, epsilon::Float64=1e-6,
                  max_pushes::Int=typemax(Int)) where T<:Unsigned

    # Normalize seeds to Set for uniform handling
    seed_set = if seeds isa T
        Set{T}([seeds])
    elseif seeds isa Vector{T}
        Set{T}(seeds)
    else
        seeds
    end

    # Initialize PPR estimate and residual vectors
    p = Dict{T, Float64}()  # PPR estimate
    r = Dict{T, Float64}()  # Residual mass to be pushed

    # Initialize: put all mass at seeds
    seed_mass = 1.0 / length(seed_set)
    for seed in seed_set
        r[seed] = seed_mass
        p[seed] = 0.0
    end

    # Track vertices we've seen
    support = Set{T}(seed_set)

    # Priority queue: (residual/degree, vertex)
    # We want to push from vertices with highest residual/degree ratio
    pq = PriorityQueue{T, Float64}()
    for seed in seed_set
        deg = length(outneighbors(graph, seed))
        if deg > 0
            pq[seed] = -(get(r, seed, 0.0) / deg)  # Negative for max-heap
        else
            pq[seed] = -get(r, seed, 0.0)
        end
    end

    push_count = 0

    # Push loop
    while !isempty(pq) && push_count < max_pushes
        # Get vertex with highest residual/degree
        u = dequeue!(pq)

        # Skip if residual below threshold
        deg_u = length(outneighbors(graph, u))
        if deg_u == 0
            continue
        end

        if get(r, u, 0.0) / deg_u < epsilon
            continue
        end

        push_count += 1

        # Push operation
        residual_u = get(r, u, 0.0)

        # Add teleport mass to PPR estimate
        p[u] = get(p, u, 0.0) + alpha * residual_u

        # Distribute random walk mass to neighbors
        mass_to_distribute = (1.0 - alpha) * residual_u / deg_u

        for v in outneighbors(graph, u)
            # Add to residual of neighbor
            old_residual = get(r, v, 0.0)
            new_residual = old_residual + mass_to_distribute
            r[v] = new_residual

            # Track that we've seen this vertex
            push!(support, v)

            # Update priority queue
            deg_v = length(outneighbors(graph, v))
            if deg_v > 0
                pq[v] = -(new_residual / deg_v)
            else
                pq[v] = -new_residual
            end
        end

        # Clear residual at u
        r[u] = 0.0
    end

    return p, support
end

"""
    select_children_by_ppr(level_vertices::Set{T}, graph, reverse_graph;
                           alpha::Float64=0.15, epsilon::Float64=1e-6,
                           min_score::Float64=1e-4, max_select::Int=typemax(Int)) where T<:Unsigned

Select children using local Personalized PageRank (PPR) ranking.

This is the state-of-the-art method for local community finding. It:
1. Computes local PPR seeded at current level vertices
2. Ranks frontier nodes by π(u)/deg(u) (PPR per degree)
3. Accepts nodes with score above threshold

The PPR/degree score identifies vertices that are:
- Strongly connected to the current community (high PPR)
- Not too high-degree (avoids hub vertices that connect to everything)

Returns Set of selected children.

Parameters:
- level_vertices: vertices in current level (seed set)
- graph: the forward graph
- reverse_graph: the reverse graph
- alpha: PPR teleport probability (default 0.15)
- epsilon: PPR convergence threshold (default 1e-6)
- min_score: minimum PPR/degree score to accept (default 1e-4)
- max_select: maximum number of children to select

Reference:
- Andersen, Chung, Lang. "Local Graph Partitioning using PageRank Vectors" (FOCS 2006)
"""
function select_children_by_ppr(level_vertices::Set{T}, graph, reverse_graph;
                                alpha::Float64=0.15, epsilon::Float64=1e-6,
                                min_score::Float64=1e-4, max_select::Int=typemax(Int)) where T<:Unsigned

    # Compute local PPR seeded at level_vertices
    ppr_scores, support = push_ppr(level_vertices, graph, reverse_graph,
                                   alpha=alpha, epsilon=epsilon)

    # Get frontier: vertices in support but not in current level
    frontier = Dict{T, Float64}()  # vertex => ppr/degree score

    for v in support
        if !(v in level_vertices)
            deg = length(outneighbors(graph, v))
            if deg > 0
                ppr_v = get(ppr_scores, v, 0.0)
                score = ppr_v / deg
                if score >= min_score
                    frontier[v] = score
                end
            end
        end
    end

    # Sort by score (descending) and select top candidates
    candidates = sort(collect(frontier), by=x->x[2], rev=true)

    selected = Set{T}()
    for (i, (vertex, score)) in enumerate(candidates)
        if i > max_select
            break
        end
        push!(selected, vertex)
    end

    return selected
end

end # module Graph

