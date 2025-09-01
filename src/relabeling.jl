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

module Relabeling

using LightGraphs, DataStructures, Logging, SparseArrays
using ..CustomTypes: UInt24, UInt40
using ..CustomLightGraphs: SimpleDiGraph, SimpleGraph, SimpleEdge
using ..Graph: get_in_degrees, get_out_degrees, get_in_out_degrees, get_reverse_graph
using ..PageRank: PR

export relabel_graph, relabel_vertices, relabel_vertices_score, relabel_vertices_lexicographic, relabel_vertices_rcm, relabel_vertices_webgraph_lex

"""
    relabel_graph(g::AbstractGraph{T}, vertex_mapping::Vector{T}) where {T<:Unsigned}

relabel the graph vertices according to the specified mapping

@param g: the graph
@param vertex_mapping: the mapping of the vertices (vertex_id[old_id] -> new_id)

@returns the relabeled graph
"""
function relabel_graph(g::AbstractGraph{T}, vertex_mapping::Dict{T,T}) where {T<:Unsigned}
	if is_directed(g)
		ng = SimpleDiGraph{T}()
	else
		ng = SimpleGraph{T}()
	end
	add_vertices!(ng, length(vertex_mapping))
	for e in edges(g)
		src_id = vertex_mapping[src(e)]
		dst_id = vertex_mapping[dst(e)]
		add_edge!(ng, src_id, dst_id)
	end
	return ng
end 

"""
    relabel_vertices(g::AbstractGraph{T}, criterion::Symbol=:in_degree) where {T<:Unsigned}

relabel the vertices of the graph according to the specified criterion

@param g: the graph
@param criterion: the criterion to relabel the vertices

Vertices are first assigned a score according to the specified criterion.
Then, they are sorted by score in descending order and relabeled accordingly.

Possible modes:
- :score
- :lexicographic (default)
- :rcm (out-degree or sym)

Possible score-based criterion:
- :in_degree
- :out_degree
- :degree
- :pagerank

Possible RCM criterion:
- :out
- :sym (symmetric)

@returns the mapping of the vertices (vertex_id[old_id] -> new_id)
"""
function relabel_vertices(g::AbstractGraph{T}, mode::Symbol=:lexicographic, criterion::Symbol=:in_degree) where {T<:Unsigned}
	if mode == :lexicographic
		return relabel_vertices_lexicographic(g)
	elseif mode == :webgraph_lex
		return relabel_vertices_webgraph_lex(g)
	elseif mode == :score
		return relabel_vertices_score(g, criterion)
	elseif mode == :rcm
		return relabel_vertices_rcm(g, criterion)
	else
		@error("Invalid mode: $mode")
	end
end

"""
    relabel_vertices_score(g::AbstractGraph{T}, criterion::Symbol=:in_degree) where {T<:Unsigned}

relabel the vertices of the graph according to the specified score-based criterion

@param g: the graph
@param criterion: the criterion to relabel the vertices

Possible score-based criterion:
- :in_degree
- :out_degree
- :degree
- :pagerank

@returns the mapping of the vertices (vertex_id[old_id] -> new_id)
"""
function relabel_vertices_score(g::AbstractGraph{T}, criterion::Symbol=:in_degree) where {T<:Unsigned}
	if criterion == :in_degree
		vertex_scores = get_in_degrees(g)
	elseif criterion == :out_degree
		vertex_scores = get_out_degrees(g)
	elseif criterion == :degree
		if is_directed(g)
			vertex_scores = Dict{T,T}()
			in_degrees, out_degrees = get_in_out_degrees(g)
			for v in vertices(g)
				vertex_scores[v] = in_degrees[v] + out_degrees[v]
			end
		else
			vertex_scores = get_out_degrees(g)
		end
	elseif criterion == :pagerank
		rg = get_reverse_graph(g)
		vertex_scores = PR(g, rg)
	else
		@error("Invalid criterion: $criterion")
	end

	# create a vector of (vertex, degree) pairs and sort by degree
	vertex_scores_pairs = [(v, vertex_scores[v]) for v in vertices(g)]
	# sort by score in descending order
	sort!(vertex_scores_pairs, by=x->x[2], rev=true)
	
	# create the mapping dictionary: old_id -> new_id
	vertex_mapping = Dict{T,T}()
	for (new_id, (old_id, _)) in enumerate(vertex_scores_pairs)
		vertex_mapping[old_id] = convert(T, new_id)
	end
	
	return vertex_mapping
end

"""
    relabel_vertices_lexicographic(g::AbstractGraph{T}) where {T<:Unsigned}

relabel the vertices of the graph lexicographically

@param g: the graph
@returns the mapping of the vertices (vertex_id[old_id] -> new_id)
"""
function relabel_vertices_lexicographic(g::AbstractGraph{T}) where {T<:Unsigned}
	vs = vertices(g)
	n = nv(g)
	vertex_mapping = Dict{T,T}()
	
	# create sparse bit vectors for each vertex's out-neighbors
	vertex_bitvectors = Vector{Tuple{T, SparseVector{Bool, T}}}()
	
	for v in vs
		out_neighbors = outneighbors(g, v)
		# create sparse bit vector: 1 for each out-neighbor
		# indices where bit is set to 1
		I = T[]  
		# values (all true)
		V = Bool[]
		# for each out-neighbor, set the bit to 1
		for neighbor in out_neighbors
			push!(I, neighbor)
			push!(V, true)
		end
		
		# create sparse bit vector of size n
		bitvector = sparsevec(I, V, n)
		push!(vertex_bitvectors, (v, bitvector))
	end
	
	# sort vertices by their sparse bit vectors in lexicographic order
	# optimized comparison function for sparse bit vectors
	function lex_compare(bv1::SparseVector{Bool, T}, bv2::SparseVector{Bool, T})
		# get sorted unique indices from both vectors
		all_indices = sort(union(bv1.nzind, bv2.nzind))
		
		# compare only at positions where at least one vector has a 1
		for i in all_indices
			val1 = i in bv1.nzind
			val2 = i in bv2.nzind
			if val1 != val2
				# false < true
				return val1 < val2
			end
		end
		# equal vectors
		return false
	end
	
	sort!(vertex_bitvectors, by=x->x[2], lt=lex_compare)
	
	# create the mapping dictionary: old_id -> new_id
	for (new_id, (old_id, _)) in enumerate(vertex_bitvectors)
		vertex_mapping[old_id] = convert(T, new_id)
	end
	
	return vertex_mapping
end

"""
    relabel_vertices_rcm(g::AbstractGraph{T}, criterion::Symbol=:sym) where {T<:Unsigned}

relabel the vertices of the graph using the RCM algorithm

@param g: the graph
@param criterion: the criterion of the RCM algorithm (:sym for symmetric, :out for out-degree)

@returns the mapping of the vertices (vertex_id[old_id] -> new_id)
"""
function relabel_vertices_rcm(g::AbstractGraph{T}, criterion::Symbol=:sym) where {T<:Unsigned}
    n = nv(g)
    vs = collect(vertices(g))

    # variant A: RCM using out-degree only, with 1 BFS sweep
    if criterion == :out
        # out-degree only
        deg = get_out_degrees(g)

        # build RCM order (global min-out-degree start, directed out-neighbors)
        order = T[]
        # visited vertices
        visited = falses(n + 1)
        # start at global min-out-degree
        start = vs[argmin(x -> get(deg, x, zero(T)), vs)]
        # queue
        queue = T[start]
        # visited
        visited[Int(start)] = true
        # add start to order
        push!(order, start)
        qidx = 1
        # BFS sweep
        while qidx <= length(queue)
            v = queue[qidx]; qidx += 1
            # neighbors not yet visited, sorted by increasing out-degree
            neigh = T[ u for u in outneighbors(g, v) if !visited[Int(u)] ]
            # sort by increasing out-degree
            sort!(neigh, by = x -> get(deg, x, zero(T)))
            # add neighbors to queue and order
            for u in neigh
                iu = Int(u)
                if !visited[iu]
                    visited[iu] = true
                    push!(queue, u)
                    push!(order, u)
                end
            end
        end
        # append any remaining vertices in original order (isolates or other components)
        for v in vs
            if !visited[Int(v)]
                push!(order, v)
            end
        end
        # reverse order
        reverse!(order)

        # build mapping: old->new
        mapping = Dict{T,T}()
        for (i, v) in enumerate(order)
            mapping[v] = T(i)
        end
        return mapping
    # variant B: RCM using symmetrized neighbors and (in+out) degree, with 2 BFS sweeps
    elseif criterion == :sym
        # symmetric neighbors (in ∪ out), degree = in+out
        indeg, outdeg = get_in_out_degrees(g)
        deg = Dict{T,T}()
        for (v, d) in outdeg
            deg[v] = d
        end
        for (v, d) in indeg
            deg[v] = get(deg, v, zero(T)) + d
        end

        # RCM using symmetrized neighbors and increasing degree expansion
        order = T[]
        # visited vertices
        visited = falses(n + 1)
        # start at global min-degree
        start = vs[argmin(x -> get(deg, x, zero(T)), vs)]
        # queue
        queue = T[start]
        # visited
        visited[Int(start)] = true
        # add start to order
        push!(order, start)
        qidx = 1
        # BFS sweep
        while qidx <= length(queue)
            v = queue[qidx]; qidx += 1
            # neighbors: union(in,out)
            neigh_set = union(outneighbors(g, v), inneighbors(g, v))
            # neighbors not yet visited, sorted by increasing degree
            neigh = T[ u for u in neigh_set if !visited[Int(u)] ]
            # sort by increasing degree
            sort!(neigh, by = x -> get(deg, x, zero(T)))
            # add neighbors to queue and order
            for u in neigh
                iu = Int(u)
                if !visited[iu]
                    visited[iu] = true
                    push!(queue, u)
                    push!(order, u)
                end
            end
        end
        # append unvisited in original order
        for v in vs
            if !visited[Int(v)]
                push!(order, v)
            end
        end
        reverse!(order)

        mapping = Dict{T,T}()
        for (i, v) in enumerate(order)
            mapping[v] = T(i)
        end
        return mapping
    else
        error("Unsupported RCM mode: $criterion. Use :out or :sym")
    end
end

"""
    relabel_vertices_webgraph_lex(g::AbstractGraph{T}) where {T<:Unsigned}

WebGraph-style lexicographic relabeling optimized for better compression.
Uses a more efficient bit-vector comparison that prioritizes compression.

@param g: the graph
@returns the mapping of the vertices (vertex_id[old_id] -> new_id)
"""
function relabel_vertices_webgraph_lex(g::AbstractGraph{T}) where {T<:Unsigned}
	vs = vertices(g)
	n = nv(g)
	vertex_mapping = Dict{T,T}()
	
	# Create compact representations optimized for web graph structure
	vertex_sigs = Vector{Tuple{T, Vector{T}}}()
	
	for v in vs
		# Get out-neighbors, already sorted by vertex ID for web graphs
		out_neighbors = sort!(collect(outneighbors(g, v)))
		
		# For web graphs, shorter neighbor lists often indicate different page types
		# Priority: empty lists first, then by lexicographic comparison
		push!(vertex_sigs, (v, out_neighbors))
	end
	
	# Sort with WebGraph-inspired comparison:
	# 1. Empty neighbor lists first (sinks/leaves)
	# 2. Length-based grouping for better reference encoding
	# 3. Lexicographic within groups
	function webgraph_compare(sig1::Tuple{T, Vector{T}}, sig2::Tuple{T, Vector{T}})
		_, neighbors1 = sig1
		_, neighbors2 = sig2
		
		# Empty lists first
		if isempty(neighbors1) && !isempty(neighbors2)
			return true
		elseif !isempty(neighbors1) && isempty(neighbors2)
			return false
		elseif isempty(neighbors1) && isempty(neighbors2)
			return false  # equal
		end
		
		# Group by length for better reference window efficiency
		len1, len2 = length(neighbors1), length(neighbors2)
		if len1 != len2 && (len1 <= 3 || len2 <= 3)
			return len1 < len2  # Small lists first for better references
		end
		
		# Lexicographic comparison
		for i in 1:min(len1, len2)
			if neighbors1[i] != neighbors2[i]
				return neighbors1[i] < neighbors2[i]
			end
		end
		
		return len1 < len2
	end
	
	sort!(vertex_sigs, lt=webgraph_compare)
	
	# Create the mapping dictionary: old_id -> new_id
	for (new_id, (old_id, _)) in enumerate(vertex_sigs)
		vertex_mapping[old_id] = convert(T, new_id)
	end
	
	return vertex_mapping
end

end # module
