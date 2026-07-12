#
# Adjacently: Julia Complex Directed Networks Library
# Copyright (C) 2016-2026 Anonymous (double-blind review)
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

using LightGraphs, DataStructures, Logging, SparseArrays, Random
using ..CustomTypes: UInt24, UInt40
using ..CustomLightGraphs: SimpleDiGraph, SimpleGraph, SimpleEdge
using ..Graph: get_in_degrees, get_out_degrees, get_in_out_degrees, get_reverse_graph, subgraph
using ..Clustering: leiden_partition
using ..Distribution: get_graph_entropy
using ..PageRank: PR

export relabel_graph, relabel_vertices, relabel_vertices_score, relabel_vertices_lexicographic, relabel_vertices_rcm, relabel_vertices_webgraph_lex,
       relabel_vertices_llp, relabel_vertices_minhash, relabel_vertices_bisection,
       relabel_graph_llp, relabel_graph_leiden_llp, merge_small_clusters,
       relabel_graph_leiden_greedy, relabel_graph_rgb, relabel_graph_leiden_rgb,
       save_llp_ordering, load_llp_ordering,
       ordering_quality_metrics, print_ordering_metrics, compare_ordering_metrics

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
	elseif mode == :llp
		# criterion selects neighbor type: :out or :sym
		return relabel_vertices_llp(g, criterion)
	elseif mode == :bisection
		return relabel_vertices_bisection(g)
	elseif mode == :minhash
		# criterion selects neighbor type: :out or :sym
		return relabel_vertices_minhash(g, criterion)
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
    _apm_label_propagation(g, neighbor_mode, gamma, passes)

Absolute Potts Model label propagation with resolution parameter gamma.
Update rule: score(l) = k_l - gamma * (size[l] - k_l) where k_l is the number
of neighbors with label l, and size[l] is the total number of vertices with label l.

Returns a label vector (Int array indexed by vertex Int id).
"""
function _apm_label_propagation(g::AbstractGraph{T}, neighbor_mode::Symbol,
                                  gamma::Float64, passes::Int) where {T<:Unsigned}
    n = Int(nv(g))
    vs = collect(vertices(g))

    # Initialize: each vertex gets its own label
    labels = collect(1:n)

    # Track label sizes
    label_size = Dict{Int,Int}()
    for i in 1:n
        label_size[i] = 1
    end

    neigh_iter(v) = neighbor_mode == :out ? outneighbors(g, v) :
                    union(outneighbors(g, v), inneighbors(g, v))

    perm = collect(1:n)
    for pass in 1:max(passes, 1)
        # Process vertices in random order
        Random.shuffle!(perm)
        changed = false
        for idx in perm
            v = vs[idx]
            vi = Int(v)
            nbs = neigh_iter(v)
            isempty(nbs) && continue

            # Count neighbor labels
            neigh_label_count = Dict{Int,Int}()
            for u in nbs
                lu = labels[Int(u)]
                neigh_label_count[lu] = get(neigh_label_count, lu, 0) + 1
            end

            # Evaluate score for each candidate label
            old_label = labels[vi]
            best_label = old_label
            # Score for current label
            k_old = get(neigh_label_count, old_label, 0)
            best_score = k_old - gamma * (get(label_size, old_label, 0) - k_old)

            for (lab, k_l) in neigh_label_count
                sz = get(label_size, lab, 0)
                # If we move to this label, its size grows by 1 (but we left old, so old shrinks)
                # For scoring: score(l) = k_l - gamma * (size[l] - k_l)
                score = k_l - gamma * (sz - k_l)
                if score > best_score || (score == best_score && lab < best_label)
                    best_score = score
                    best_label = lab
                end
            end

            if best_label != old_label
                # Update label sizes
                label_size[old_label] = get(label_size, old_label, 1) - 1
                if label_size[old_label] <= 0
                    delete!(label_size, old_label)
                end
                label_size[best_label] = get(label_size, best_label, 0) + 1
                labels[vi] = best_label
                changed = true
            end
        end
        !changed && break
    end

    return labels
end

"""
    relabel_vertices_llp(g::AbstractGraph{T}, neighbor_mode::Symbol=:sym; passes::Int=5, K::Int=10) where {T<:Unsigned}

Multi-resolution Layered Label Propagation (Boldi et al.). Runs Absolute Potts Model
label propagation at multiple resolution parameters gamma in {0} ∪ {2^(-i) : i=0,...,K},
then composes the orderings to build a final vertex mapping.

@param g: the graph
@param neighbor_mode: :out (directed out-neighbors) or :sym (union of in/out)
@param passes: number of label-propagation passes per resolution
@param K: number of resolution levels (gamma = 2^(-i) for i=0,...,K)
@returns Dict{T,T}: mapping old_id -> new_id
"""
function relabel_vertices_llp(g::AbstractGraph{T}, neighbor_mode::Symbol=:sym;
                               passes::Int=5, K::Int=10) where {T<:Unsigned}
    n = Int(nv(g))
    vs = collect(vertices(g))

    # Generate gamma values: {0} ∪ {2^(-i) : i=0,...,K}
    gammas = Float64[0.0]
    for i in 0:K
        push!(gammas, 2.0^(-i))
    end
    @info "LLP: running APM at $(length(gammas)) resolution levels (K=$K, passes=$passes)"

    # Start with identity ordering
    current_order = collect(1:n)  # current_order[position] = vertex_int_id

    for (gi, gamma) in enumerate(gammas)
        # Run APM label propagation at this gamma
        labels = _apm_label_propagation(g, neighbor_mode, gamma, passes)

        # Group vertices by label, preserving current order
        groups = Dict{Int, Vector{Int}}()
        for pos in 1:n
            vid = current_order[pos]
            lab = labels[vid]
            push!(get!(groups, lab, Int[]), vid)
        end

        # Sort group keys by size descending, tie-break by smallest member position
        group_keys = collect(keys(groups))
        sort!(group_keys, by = k -> (-length(groups[k]), k))

        # Compose: within each group, vertices keep their prior relative order
        new_order = Int[]
        sizehint!(new_order, n)
        for k in group_keys
            append!(new_order, groups[k])
        end
        current_order = new_order
    end

    # Build final mapping: old_id -> new_id
    mapping = Dict{T,T}()
    for (new_id, old_vid) in enumerate(current_order)
        mapping[T(old_vid)] = T(new_id)
    end
    return mapping
end

"""
    relabel_vertices_minhash(g::AbstractGraph{T}, neighbor_mode::Symbol=:out; k::Int=32, seed::UInt64=0x9e3779b97f4a7c15) where {T<:Unsigned}

MinHash-based relabeling: compute a k-dimensional MinHash signature of each
vertex's neighbor set (out-neighbors by default), then sort vertices by their
signature to cluster similar adjacency lists before relabeling.

@param g: the graph
@param neighbor_mode: :out (default) or :sym (union of in/out)
@param k: number of hash functions (signature length)
@param seed: RNG seed to derive hash functions
@returns Dict{T,T}: mapping old_id -> new_id
"""
function relabel_vertices_minhash(g::AbstractGraph{T}, neighbor_mode::Symbol=:out; k::Int=32, seed::UInt64=0x9e3779b97f4a7c15) where {T<:Unsigned}
    n = nv(g)
    k = max(1, k)
    vs = collect(vertices(g))

    # derive k pairs (a,b) for universal hashing h(x) = a*x + b in UInt64 space
    rng = Random.MersenneTwister(seed)
    a = [rand(rng, UInt64) | 1 for _ in 1:k]  # force odd to reduce trivial cycles
    b = [rand(rng, UInt64) for _ in 1:k]

    # neighbor iterator
    neigh_iter(v) = neighbor_mode == :out ? outneighbors(g, v) : union(outneighbors(g, v), inneighbors(g, v))

    # compute signatures
    UMAX = typemax(UInt64)
    signatures = Vector{Tuple{T, Vector{UInt64}}}(undef, length(vs))
    idx = 0
    for v in vs
        idx += 1
        sig = fill(UMAX, k)
        for u in neigh_iter(v)
            xu = UInt64(u)
            @inbounds for i in 1:k
                # mix then take min
                hv = a[i] * (xu ⊻ 0x9e3779b97f4a7c15) + b[i]
                if hv < sig[i]
                    sig[i] = hv
                end
            end
        end
        signatures[idx] = (v, sig)
    end

    # sort by signature lexicographically; tie-break by out-degree and id
    outdeg = get_out_degrees(g)
    sort!(signatures, lt = (x,y)->(x[2] < y[2] || (x[2] == y[2] && (get(outdeg, x[1], zero(T)) < get(outdeg, y[1], zero(T)) || (get(outdeg, x[1], zero(T)) == get(outdeg, y[1], zero(T)) && x[1] < y[1])))))

    mapping = Dict{T,T}()
    for (i, (v, _)) in enumerate(signatures)
        mapping[v] = T(i)
    end
    return mapping
end

"""
    relabel_vertices_bisection(g::AbstractGraph{T}; max_iters::Int=20, leaf_size::Int=32) where {T<:Unsigned}

Recursive bisection relabeling (Dhulipala et al., KDD 2016). Recursively splits the
vertex set into two halves and refines each split via move-gain optimization using
a BiMLogA cost function. Produces a vertex ordering that groups structurally similar
vertices together for better compression.

@param g: the graph
@param max_iters: max refinement iterations per bisection level
@param leaf_size: stop recursion when partition size <= leaf_size
@returns Dict{T,T}: mapping old_id -> new_id
"""
function relabel_vertices_bisection(g::AbstractGraph{T}; max_iters::Int=20, leaf_size::Int=32) where {T<:Unsigned}
    n = Int(nv(g))
    vs = collect(vertices(g))

    # Pre-build adjacency as Vector{Vector{Int}} for O(1) lookups
    adj = Vector{Vector{Int}}(undef, n)
    for v in vs
        adj[Int(v)] = Int.(collect(outneighbors(g, v)))
    end
    # Also build reverse adjacency for symmetric neighbors
    radj = [Int[] for _ in 1:n]
    for v in vs
        vi = Int(v)
        for u in adj[vi]
            push!(radj[u], vi)
        end
    end
    # Symmetric adjacency: union of out and in neighbors
    sym_adj = Vector{Vector{Int}}(undef, n)
    for i in 1:n
        sym_adj[i] = sort(unique(vcat(adj[i], radj[i])))
    end

    @info "Bisection: starting recursive bisection on $n vertices (max_iters=$max_iters, leaf_size=$leaf_size)"

    # Run recursive bisection to get ordering
    verts = collect(1:n)
    order = _recursive_bisect(verts, sym_adj, max_iters, leaf_size)

    # Build mapping: old_id -> new_id
    mapping = Dict{T,T}()
    for (new_id, old_vid) in enumerate(order)
        mapping[T(old_vid)] = T(new_id)
    end
    return mapping
end

"""
    _move_gain(v, partition, deg1, deg2, size1, size2, adj)

Compute the BiMLogA gain of moving vertex v from its current partition to the other.
Cost for v in partition p: deg_p(v) * log(n_p / (deg_p(v)+1)) + deg_q(v) * log(n_q / (deg_q(v)+1))
Gain = cost_before - cost_after (positive = beneficial move).
"""
function _move_gain(v::Int, partition::Vector{Int}, deg1::Vector{Int}, deg2::Vector{Int},
                    size1::Int, size2::Int)
    from = partition[v]
    d1 = deg1[v]
    d2 = deg2[v]

    if from == 1
        n_from, n_to = size1, size2
        d_from, d_to = d1, d2
    else
        n_from, n_to = size2, size1
        d_from, d_to = d2, d1
    end

    # Cost before move
    cost_before = d_from * log(max(n_from, 1) / (d_from + 1)) +
                  d_to * log(max(n_to, 1) / (d_to + 1))
    # Cost after move (v moves from -> to, so sizes change by 1)
    cost_after = d_from * log(max(n_to + 1, 1) / (d_from + 1)) +
                 d_to * log(max(n_from - 1, 1) / (d_to + 1))

    return cost_before - cost_after
end

"""
    _apply_move!(v, partition, deg1, deg2, adj)

Move vertex v to the other partition and update neighbor degree counts.
"""
function _apply_move!(v::Int, partition::Vector{Int}, deg1::Vector{Int}, deg2::Vector{Int}, adj::Vector{Vector{Int}})
    old_part = partition[v]
    new_part = old_part == 1 ? 2 : 1
    partition[v] = new_part

    # Update neighbor degree counts
    for u in adj[v]
        if old_part == 1
            # v left partition 1: neighbors lose a partition-1 neighbor, gain a partition-2 neighbor
            deg1[u] -= 1
            deg2[u] += 1
        else
            deg2[u] -= 1
            deg1[u] += 1
        end
    end
end

"""
    _recursive_bisect(verts, adj, max_iters, leaf_size)

Recursively bisect the vertex set and return the ordering.
"""
function _recursive_bisect(verts::Vector{Int}, adj::Vector{Vector{Int}}, max_iters::Int, leaf_size::Int)
    n = length(verts)
    if n <= leaf_size
        return verts
    end

    # Initial partition: first half = 1, second half = 2
    half = div(n, 2)
    vert_set = Set(verts)

    # Map global vertex IDs to local partition assignment
    partition = zeros(Int, length(adj))  # 0 = not in this subproblem
    for (i, v) in enumerate(verts)
        partition[v] = i <= half ? 1 : 2
    end

    # Compute initial degree counts within each partition
    deg1 = zeros(Int, length(adj))  # number of neighbors in partition 1
    deg2 = zeros(Int, length(adj))  # number of neighbors in partition 2
    for v in verts
        for u in adj[v]
            if partition[u] == 1
                deg1[v] += 1
            elseif partition[u] == 2
                deg2[v] += 1
            end
        end
    end

    size1 = half
    size2 = n - half

    # Refinement iterations
    for iter in 1:max_iters
        # Compute gains for all vertices in the subproblem
        gains = [(v, _move_gain(v, partition, deg1, deg2, size1, size2)) for v in verts]
        sort!(gains, by = x -> -x[2])  # descending by gain

        moved = false
        # Greedily swap pairs while net gain > 0
        # Process vertices one by one
        for (v, g) in gains
            if g <= 0.0
                break
            end
            # Check balance constraint: don't let either partition become too small
            from = partition[v]
            if from == 1 && size1 <= div(n, 4)
                continue
            end
            if from == 2 && size2 <= div(n, 4)
                continue
            end

            _apply_move!(v, partition, deg1, deg2, adj)
            if from == 1
                size1 -= 1; size2 += 1
            else
                size2 -= 1; size1 += 1
            end
            moved = true
        end

        !moved && break
    end

    # Collect vertices in each partition preserving input order
    part1 = Int[]
    part2 = Int[]
    for v in verts
        if partition[v] == 1
            push!(part1, v)
        else
            push!(part2, v)
        end
    end

    # Recurse on each half
    left_order = _recursive_bisect(part1, adj, max_iters, leaf_size)
    right_order = _recursive_bisect(part2, adj, max_iters, leaf_size)

    return vcat(left_order, right_order)
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

"""
    relabel_vertices_bfs(g::AbstractGraph{T}, start_vertex::Union{Nothing,T}=nothing) where {T<:Unsigned}

Relabel vertices using Breadth-First Search (BFS) traversal order.
This is the ordering strategy used by WebGraph for optimal compression.

The BFS ordering tends to group vertices that are close in the graph structure,
which improves reference encoding effectiveness since nearby vertices in the
traversal order are likely to have similar adjacency lists.

@param g::AbstractGraph{T}: the graph to relabel
@param start_vertex::Union{Nothing,T}: starting vertex for BFS (if nothing, uses highest out-degree vertex)
@return vertex_mapping::Dict{T,T}: mapping from old vertex IDs to new vertex IDs
"""
function relabel_vertices_bfs(g::AbstractGraph{T}, start_vertex::Union{Nothing,T}=nothing) where {T<:Unsigned}
    n = nv(g)
    vertex_mapping = Dict{T,T}()

    # If no start vertex specified, choose vertex with highest out-degree
    if start_vertex === nothing
        max_degree = -1
        start_vertex = T(1)
        for v in vertices(g)
            deg = length(outneighbors(g, v))
            if deg > max_degree
                max_degree = deg
                start_vertex = v
            end
        end
    end

    # BFS traversal
    visited = Set{T}()
    queue = T[]
    bfs_order = T[]

    # Start BFS from the chosen vertex
    push!(queue, start_vertex)
    push!(visited, start_vertex)

    while !isempty(queue)
        v = popfirst!(queue)
        push!(bfs_order, v)

        # Visit neighbors in sorted order for deterministic results
        neighbors = sort(collect(outneighbors(g, v)))
        for neighbor in neighbors
            if !(neighbor in visited)
                push!(visited, neighbor)
                push!(queue, neighbor)
            end
        end
    end

    # Handle disconnected components: add unvisited vertices in degree order
    if length(visited) < n
        unvisited = T[]
        for v in vertices(g)
            if !(v in visited)
                push!(unvisited, v)
            end
        end

        # Sort unvisited by out-degree (descending)
        sort!(unvisited, by=v -> -length(outneighbors(g, v)))

        # Continue BFS from each unvisited component
        for start_v in unvisited
            if start_v in visited
                continue
            end

            push!(queue, start_v)
            push!(visited, start_v)

            while !isempty(queue)
                v = popfirst!(queue)
                push!(bfs_order, v)

                neighbors = sort(collect(outneighbors(g, v)))
                for neighbor in neighbors
                    if !(neighbor in visited)
                        push!(visited, neighbor)
                        push!(queue, neighbor)
                    end
                end
            end
        end
    end

    # Create mapping: old_id -> new_id
    for (new_id, old_id) in enumerate(bfs_order)
        vertex_mapping[old_id] = T(new_id)
    end

    return vertex_mapping
end

"""
    relabel_graph_llp(g::AbstractGraph{T}; llp_mode::Symbol=:sym, passes::Int=5) where {T<:Unsigned}

Apply global LLP reordering to graph `g` and return the relabeled graph and vertex mapping.

@param g: the graph
@param llp_mode: LLP neighbor mode (:sym or :out)
@param passes: number of LLP passes
@returns (relabeled_graph, vertex_mapping::Dict{T,T})
"""
function relabel_graph_llp(g::AbstractGraph{T}; llp_mode::Symbol=:sym, passes::Int=5) where {T<:Unsigned}
    mapping = relabel_vertices_llp(g, llp_mode; passes=passes)
    return relabel_graph(g, mapping), mapping
end

"""
    merge_small_clusters(g, part, min_size) -> Vector{Int}

Merge every cluster smaller than `min_size` into the neighbouring cluster it
shares the most edges with (one pass, smallest-first). `leiden_partition` here is
a lightweight Louvain approximation that over-fragments into many tiny clusters;
grouping the small ones improves within-cluster locality and hence compression
across all reference-based encoders. Returns a new label vector; `min_size <= 1`
is a no-op.
"""
function merge_small_clusters(g::AbstractGraph{T}, part::Vector{Int}, min_size::Int) where {T<:Unsigned}
    min_size <= 1 && return part
    n = Int(nv(g))
    lab = copy(part)
    verts = Dict{Int,Vector{Int}}()
    for v in 1:n; push!(get!(verts, lab[v], Int[]), v); end
    smalls = sort([l for (l, vs) in verts if length(vs) < min_size]; by = l -> length(verts[l]))
    for l in smalls
        haskey(verts, l) || continue
        vs = verts[l]
        length(vs) < min_size || continue
        cnt = Dict{Int,Int}()
        for v in vs
            for u in outneighbors(g, T(v)); cl = lab[Int(u)]; cl == l || (cnt[cl] = get(cnt, cl, 0) + 1); end
            for u in inneighbors(g, T(v));  cl = lab[Int(u)]; cl == l || (cnt[cl] = get(cnt, cl, 0) + 1); end
        end
        isempty(cnt) && continue
        tgt = argmax(cnt)
        for v in vs; lab[v] = tgt; end
        append!(verts[tgt], vs); delete!(verts, l)
    end
    return lab
end

# Build the concatenated per-cluster-LLP vertex order for a partition.
function _leiden_llp_order(g::AbstractGraph{T}, part; llp_mode, llp_passes, sort_clusters) where {T<:Unsigned}
    label_to_idx = Dict{Int,Int}()
    fine_clusters = Vector{Vector{T}}()
    for v in vertices(g)
        l = part[Int(v)]
        if !haskey(label_to_idx, l)
            label_to_idx[l] = length(label_to_idx) + 1
            push!(fine_clusters, T[])
        end
        push!(fine_clusters[label_to_idx[l]], T(v))
    end
    if sort_clusters == :size_desc
        sort!(fine_clusters, by=length, rev=true)
    end
    new_order = T[]
    sizehint!(new_order, nv(g))
    for C in fine_clusters
        if length(C) <= 2
            append!(new_order, sort(C))
        else
            sg, oni, _ = subgraph(g, C)
            mapping = relabel_vertices_llp(sg, llp_mode; passes=llp_passes)
            sort!(C, by = v -> Int(mapping[oni[v]]))
            append!(new_order, C)
        end
    end
    return new_order, length(fine_clusters)
end

# Pick the small-cluster-merge threshold that minimises the cheap delta-encoding
# entropy proxy (no encoder run). The proxy's minimum coincides with the true
# BPE minimum (validated on Web-Google/Amazon), and the optimum is
# dataset-dependent, so a single fixed threshold will not do.
function _auto_merge_threshold(g::AbstractGraph{T}, part; llp_mode, llp_passes, sort_clusters,
                               grid=[0, 10, 20, 50, 100, 200, 400, 800]) where {T<:Unsigned}
    best_ms, best_de = 0, Inf
    for ms in grid
        p = ms <= 1 ? part : merge_small_clusters(g, part, ms)
        order, _ = _leiden_llp_order(g, p; llp_mode, llp_passes, sort_clusters)
        vmap = Dict{T,T}(old => T(i) for (i, old) in enumerate(order))
        de = get_graph_entropy(relabel_graph(g, vmap), :bits_per_edge, :delta)
        de < best_de && (best_de = de; best_ms = ms)
    end
    return best_ms
end

"""
    relabel_graph_leiden_llp(g; llp_mode=:sym, llp_passes=5, sort_clusters=:size_desc,
                             merge_clusters=nothing)

`merge_clusters` controls the small-cluster merge (see [`merge_small_clusters`]):
`nothing` (default) keeps the legacy behaviour with no merge; an `Int` uses that
fixed min-size threshold; `:auto` sweeps a threshold grid and picks the one that
minimises the delta-entropy proxy (dataset-adaptive, but orders the graph once
per grid point — expensive on very large graphs, so prefer a fixed `Int` there).
"""
function relabel_graph_leiden_llp(g::AbstractGraph{T}; llp_mode::Symbol=:sym, llp_passes::Int=5,
                                  sort_clusters::Symbol=:size_desc,
                                  merge_clusters::Union{Nothing,Integer,Symbol}=nothing,
                                  return_clusters::Bool=false) where {T<:Unsigned}
    # Step 1: Leiden partition → fine clusters (label vector)
    part = leiden_partition(g)

    # Step 2 (optional): merge small clusters
    if merge_clusters === :auto
        ms = _auto_merge_threshold(g, part; llp_mode, llp_passes, sort_clusters)
        @info "Leiden+LLP: auto-selected merge min_size=$ms"
        part = merge_small_clusters(g, part, ms)
    elseif merge_clusters isa Integer && merge_clusters > 1
        part = merge_small_clusters(g, part, Int(merge_clusters))
    end

    # Step 3: per-cluster LLP ordering
    new_order, n_clusters = _leiden_llp_order(g, part; llp_mode, llp_passes, sort_clusters)
    szc = Dict{Int,Int}(); for l in part; szc[l] = get(szc, l, 0) + 1; end
    top_sizes = sort(collect(values(szc)), rev=true)[1:min(5, n_clusters)]
    @info "Leiden+LLP: $n_clusters clusters, largest: $top_sizes"

    # Step 4: Build vertex mapping and relabel
    vertex_map = Dict{T,T}()
    for (new_id, old_id) in enumerate(new_order)
        vertex_map[old_id] = T(new_id)
    end
    if return_clusters
        # Cluster sizes in concatenation order: in the NEW labeling each Leiden
        # cluster occupies one contiguous ID range (new_order concatenates whole
        # clusters), so the sizes fully describe the partition (implicit_ranges).
        seg_sizes = Int[]
        prev_label = -1
        for old_id in new_order
            l = part[Int(old_id)]
            if l != prev_label
                push!(seg_sizes, 0)
                prev_label = l
            end
            seg_sizes[end] += 1
        end
        return relabel_graph(g, vertex_map), vertex_map, seg_sizes
    end
    return relabel_graph(g, vertex_map), vertex_map
end

# ─────────────────────────────────────────────────────────────────────────────
# Leiden + Greedy Gap Minimization
# ─────────────────────────────────────────────────────────────────────────────

"""
    relabel_graph_leiden_greedy(g; sort_clusters=:size_desc)

Reorder vertices using Leiden community detection followed by a greedy
within-cluster ordering that maximizes consecutive neighbor overlap.

Unlike Leiden+LLP which uses label propagation (optimizing modularity),
this method directly optimizes for reference-copy quality by greedily
placing the vertex with maximum neighbor overlap to the last-placed vertex.
This preserves both first-order locality (high consecutive Jaccard) and
second-order locality (compact residuals after reference copy).

Algorithm:
1. Leiden partition → fine clusters
2. Sort clusters by size descending
3. Within each cluster (>2 vertices): build induced subgraph, then
   greedily order starting from the highest-degree vertex, always
   picking the unplaced vertex with maximum overlap to the last placed.
4. Concatenate ordered clusters → global vertex mapping

Returns `(relabeled_graph, vertex_mapping)`.
"""
function relabel_graph_leiden_greedy(g::AbstractGraph{T};
        sort_clusters::Symbol=:size_desc) where {T<:Unsigned}
    n = nv(g)

    # Step 1: Leiden partition
    part = leiden_partition(g)
    label_to_idx = Dict{Int,Int}()
    fine_clusters = Vector{Vector{T}}()
    for v in vertices(g)
        l = part[Int(v)]
        if !haskey(label_to_idx, l)
            label_to_idx[l] = length(label_to_idx) + 1
            push!(fine_clusters, T[])
        end
        push!(fine_clusters[label_to_idx[l]], T(v))
    end

    if sort_clusters == :size_desc
        sort!(fine_clusters, by=length, rev=true)
    end

    n_clusters = length(fine_clusters)
    top_sizes = sort([length(C) for C in fine_clusters], rev=true)[1:min(5, n_clusters)]
    @info "Leiden+Greedy: $n_clusters fine clusters, largest: $top_sizes"

    # Step 2: Greedy ordering within each cluster
    new_order = T[]
    sizehint!(new_order, n)

    for C in fine_clusters
        if length(C) <= 2
            append!(new_order, sort(C))
        else
            ordered = _greedy_overlap_order(g, C)
            append!(new_order, ordered)
        end
    end

    # Step 3: Build mapping
    vertex_map = Dict{T,T}()
    for (new_id, old_id) in enumerate(new_order)
        vertex_map[old_id] = T(new_id)
    end
    return relabel_graph(g, vertex_map), vertex_map
end

"""Greedy ordering of vertices in `cluster`: start from highest-degree vertex,
then always pick the unplaced vertex with maximum neighbor overlap to last placed."""
function _greedy_overlap_order(g::AbstractGraph{T}, cluster::Vector{T}) where {T<:Unsigned}
    # Build neighbor sets for cluster members
    member_set = Set{T}(cluster)
    adj = Dict{T, Set{T}}()
    for v in cluster
        adj[v] = Set{T}(u for u in outneighbors(g, v))
    end

    # Start from highest-degree vertex in cluster
    remaining = Set{T}(cluster)
    start_v = cluster[1]
    best_deg = 0
    for v in cluster
        d = length(adj[v])
        if d > best_deg
            best_deg = d
            start_v = v
        end
    end

    order = T[start_v]
    delete!(remaining, start_v)
    last_nbrs = adj[start_v]

    while !isempty(remaining)
        best_v = first(remaining)
        best_ov = -1

        for v in remaining
            ov = length(intersect(adj[v], last_nbrs))
            if ov > best_ov
                best_ov = ov
                best_v = v
            end
        end

        push!(order, best_v)
        delete!(remaining, best_v)
        last_nbrs = adj[best_v]
    end

    return order
end

# ─────────────────────────────────────────────────────────────────────────────
# Recursive Graph Bisection (RGB)
# ─────────────────────────────────────────────────────────────────────────────

"""
    relabel_graph_rgb(g; max_iters=20, min_partition=64)

Reorder vertices using Recursive Graph Bisection (Dhulipala et al., KDD 2016).

Directly optimizes the log-gap compression cost: Σ log₂(gap) across all
adjacency lists, where gap = difference between consecutive sorted neighbors.
At each recursion level, vertices are swapped between two halves if the swap
reduces the total log-gap cost.

Algorithm:
1. Split vertices into two halves (initialized by degree sorting)
2. Iteratively swap vertices between halves to minimize log-gap cost
3. Recurse on each half until partition size < min_partition
4. Concatenate final orderings

Parameters:
- `max_iters`: maximum swap iterations per recursion level (default: 20)
- `min_partition`: stop recursing below this size (default: 64)

Returns `(relabeled_graph, vertex_mapping)`.

Reference: Dhulipala, Kabiljo, Karrer, Ottaviano, Pupyrev, Shalita.
"Compressing Graphs and Indexes with Recursive Graph Bisection." KDD 2016.
"""
function relabel_graph_rgb(g::AbstractGraph{T};
        max_iters::Int=20, min_partition::Int=64) where {T<:Unsigned}
    n = Int(nv(g))

    # Build sorted neighbor lists (using Int for speed)
    adj = Vector{Vector{Int}}(undef, n)
    for v in vertices(g)
        adj[Int(v)] = sort(Int.(outneighbors(g, v)))
    end

    # Initial order: sort by degree descending (heuristic initialization)
    order = sortperm([length(adj[v]) for v in 1:n], rev=true)

    # Recursive bisection
    _rgb_recurse!(adj, order, 1, n, max_iters, min_partition)

    # Build mapping: order[position] = old_vertex_id
    vertex_map = Dict{T,T}()
    for (new_id, old_id) in enumerate(order)
        vertex_map[T(old_id)] = T(new_id)
    end
    return relabel_graph(g, vertex_map), vertex_map
end

"""Recursive bisection: optimize log-gap cost for order[lo:hi]."""
function _rgb_recurse!(adj::Vector{Vector{Int}}, order::Vector{Int},
                       lo::Int, hi::Int, max_iters::Int, min_partition::Int)
    size = hi - lo + 1
    size <= min_partition && return

    mid = lo + size ÷ 2

    # Build position lookup: vertex_id → current position
    pos = Dict{Int,Int}()
    for i in lo:hi
        pos[order[i]] = i
    end

    # Iterative swapping to minimize log-gap cost
    for iter in 1:max_iters
        n_swaps = 0

        for i in lo:mid
            v_left = order[i]

            # Find best swap partner in right half
            best_j = 0
            best_gain = 0.0

            for j in (mid+1):hi
                v_right = order[j]
                gain = _rgb_swap_gain(adj, pos, v_left, v_right, mid)
                if gain > best_gain
                    best_gain = gain
                    best_j = j
                end
            end

            if best_j > 0
                # Perform swap
                v_right = order[best_j]
                order[i], order[best_j] = order[best_j], order[i]
                pos[v_right] = i
                pos[v_left] = best_j
                n_swaps += 1
            end
        end

        n_swaps == 0 && break
    end

    # Recurse on each half
    _rgb_recurse!(adj, order, lo, mid, max_iters, min_partition)
    _rgb_recurse!(adj, order, mid + 1, hi, max_iters, min_partition)
end

"""Compute the log-gap cost reduction from swapping v_left (in left half)
with v_right (in right half). Positive = beneficial swap."""
function _rgb_swap_gain(adj::Vector{Vector{Int}}, pos::Dict{Int,Int},
                        v_left::Int, v_right::Int, mid::Int)
    gain = 0.0

    # For each neighbor of v_left: moving v_left to right half
    # means its neighbors that are in the left half get a larger gap cost
    for u in adj[v_left]
        haskey(pos, u) || continue
        p_u = pos[u]
        if p_u <= mid
            gain += 0.5  # neighbor in same half as v_right (after swap)
        else
            gain -= 0.5  # neighbor in different half
        end
    end

    # For each neighbor of v_right: moving v_right to left half
    for u in adj[v_right]
        haskey(pos, u) || continue
        p_u = pos[u]
        if p_u <= mid
            gain -= 0.5  # was in same half, now different
        else
            gain += 0.5  # was in different half, now same
        end
    end

    return gain
end

# ─────────────────────────────────────────────────────────────────────────────
# Leiden + Recursive Graph Bisection
# ─────────────────────────────────────────────────────────────────────────────

"""
    relabel_graph_leiden_rgb(g; max_iters=20, min_partition=32, sort_clusters=:size_desc)

Reorder vertices using Leiden community detection followed by Recursive Graph
Bisection within each cluster. Combines community-aware grouping with
log-gap-optimal local ordering.

Algorithm:
1. Leiden partition → fine clusters
2. Sort clusters by size descending
3. Within each cluster (>min_partition vertices): apply RGB to optimize
   the log-gap cost locally on the induced subgraph
4. Small clusters: sort by vertex ID
5. Concatenate ordered clusters → global vertex mapping

This hybrid approach preserves the community structure (good for inter-cluster
encoding) while using RGB's direct cost optimization within each community
(better residual locality than LLP).

Returns `(relabeled_graph, vertex_mapping)`.
"""
function relabel_graph_leiden_rgb(g::AbstractGraph{T};
        max_iters::Int=20, min_partition::Int=32,
        sort_clusters::Symbol=:size_desc) where {T<:Unsigned}
    n = nv(g)

    # Step 1: Leiden partition
    part = leiden_partition(g)
    label_to_idx = Dict{Int,Int}()
    fine_clusters = Vector{Vector{T}}()
    for v in vertices(g)
        l = part[Int(v)]
        if !haskey(label_to_idx, l)
            label_to_idx[l] = length(label_to_idx) + 1
            push!(fine_clusters, T[])
        end
        push!(fine_clusters[label_to_idx[l]], T(v))
    end

    if sort_clusters == :size_desc
        sort!(fine_clusters, by=length, rev=true)
    end

    n_clusters = length(fine_clusters)
    top_sizes = sort([length(C) for C in fine_clusters], rev=true)[1:min(5, n_clusters)]
    @info "Leiden+RGB: $n_clusters fine clusters, largest: $top_sizes"

    # Step 2: RGB within each cluster
    new_order = T[]
    sizehint!(new_order, n)

    for C in fine_clusters
        if length(C) <= 2
            append!(new_order, sort(C))
        else
            # Build local adjacency for this cluster
            member_set = Set{Int}(Int.(C))
            local_n = length(C)
            old_to_local = Dict{Int,Int}()
            local_to_old = Vector{Int}(undef, local_n)
            for (i, v) in enumerate(C)
                old_to_local[Int(v)] = i
                local_to_old[i] = Int(v)
            end

            # Build local neighbor lists (only intra-cluster edges)
            local_adj = Vector{Vector{Int}}(undef, local_n)
            for i in 1:local_n
                v = T(local_to_old[i])
                nbrs = Int[]
                for u in outneighbors(g, v)
                    u_int = Int(u)
                    if haskey(old_to_local, u_int)
                        push!(nbrs, old_to_local[u_int])
                    end
                end
                local_adj[i] = sort(nbrs)
            end

            # Apply RGB on local ordering
            local_order = sortperm([length(local_adj[i]) for i in 1:local_n], rev=true)
            _rgb_recurse!(local_adj, local_order, 1, local_n, max_iters, min_partition)

            # Map back to global vertex IDs
            for pos in 1:local_n
                push!(new_order, T(local_to_old[local_order[pos]]))
            end
        end
    end

    # Step 3: Build mapping
    vertex_map = Dict{T,T}()
    for (new_id, old_id) in enumerate(new_order)
        vertex_map[old_id] = T(new_id)
    end
    return relabel_graph(g, vertex_map), vertex_map
end

"""
    save_llp_ordering(path, mapping, n)

Save an LLP vertex ordering to a binary file.
Format: Int32 header (number of vertices), followed by n UInt32 values
where index = old vertex ID and value = new vertex ID.
"""
function save_llp_ordering(path::String, mapping::Dict{T,T}, n::Int) where {T}
    ordering = Vector{UInt32}(undef, n)
    for (old_id, new_id) in mapping
        ordering[Int(old_id)] = UInt32(new_id)
    end
    open(path, "w") do f
        write(f, Int32(n))
        write(f, ordering)
    end
    @info "Saved LLP ordering to $path ($(filesize(path)) bytes)"
end

"""
    load_llp_ordering(path, T)

Load an LLP vertex ordering from a binary file.
Returns a Dict{T,T} mapping old_id → new_id.
"""
function load_llp_ordering(path::String, ::Type{T}) where {T}
    data = open(path, "r") do f
        n = read(f, Int32)
        ordering = Vector{UInt32}(undef, n)
        read!(f, ordering)
        ordering
    end
    mapping = Dict{T,T}()
    for (old_id, new_id) in enumerate(data)
        mapping[T(old_id)] = T(new_id)
    end
    @info "Loaded LLP ordering from $path ($(length(data)) vertices)"
    return mapping
end

"""
    _sorted_overlap_count(a::Vector{Int}, b::Vector{Int})

Count elements common to both sorted integer vectors using merge scan. O(|a| + |b|).
"""
function _sorted_overlap_count(a::Vector{Int}, b::Vector{Int})
    count = 0
    i, j = 1, 1
    na, nb = length(a), length(b)
    @inbounds while i <= na && j <= nb
        if a[i] == b[j]
            count += 1; i += 1; j += 1
        elseif a[i] < b[j]
            i += 1
        else
            j += 1
        end
    end
    return count
end

"""
    _sorted_setdiff(a::Vector{Int}, b::Vector{Int})

Return elements in sorted vector `a` that are not in sorted vector `b`. O(|a| + |b|).
"""
function _sorted_setdiff(a::Vector{Int}, b::Vector{Int})
    result = Int[]
    i, j = 1, 1
    na, nb = length(a), length(b)
    @inbounds while i <= na
        if j > nb || a[i] < b[j]
            push!(result, a[i])
            i += 1
        elseif a[i] == b[j]
            i += 1; j += 1
        else
            j += 1
        end
    end
    return result
end

"""
    ordering_quality_metrics(g::AbstractGraph{T}; window::Int=7, mil::Int=4, ref_sample::Int=0) where {T<:Unsigned}

Compute ordering quality metrics that characterize how well a vertex ordering
supports different compression strategies.

Gap/interval metrics are computed over all vertices. Reference overlap metrics
can be sampled for large graphs by setting `ref_sample > 0`.

Returns a Dict{Symbol, Any} with:
- Gap metrics (lower = better for CS/BG delta/interval encoding):
  - :avg_log_gap — mean log₂(gap+1) between consecutive sorted neighbors
  - :avg_successor_dist — mean |v - u| for all edges (v,u)
  - :avg_successor_dist_norm — avg_successor_dist / n
- Interval metrics (higher = better for CS/BG):
  - :interval_coverage — fraction of neighbors in consecutive runs ≥ mil
  - :avg_intervals_per_vertex — mean interval count per vertex
- Reference metrics (higher = better for BV/CG copy-list encoding):
  - :avg_best_overlap — mean best overlap count within window
  - :avg_best_jaccard — mean best Jaccard similarity within window
  - :ref_hit_rate — fraction of vertices with any ref overlap > 0
  - :avg_copy_fraction — mean fraction of neighbors copyable from best ref
"""
function ordering_quality_metrics(g::AbstractGraph{T}; window::Int=7, mil::Int=4, ref_sample::Int=0) where {T<:Unsigned}
    n = Int(nv(g))

    # Pre-build sorted neighbor lists
    adj = Vector{Vector{Int}}(undef, n)
    total_edges = 0
    for v in vertices(g)
        vi = Int(v)
        adj[vi] = sort(Int.(collect(outneighbors(g, v))))
        total_edges += length(adj[vi])
    end

    # === Gap statistics (between consecutive sorted neighbors) ===
    total_log_gap = 0.0
    n_gaps = 0
    for vi in 1:n
        nbrs = adj[vi]
        @inbounds for j in 2:length(nbrs)
            gap = nbrs[j] - nbrs[j-1]
            total_log_gap += log2(Float64(gap) + 1.0)
            n_gaps += 1
        end
    end
    avg_log_gap = n_gaps > 0 ? total_log_gap / n_gaps : 0.0

    # === Successor distance (|v - neighbor|) ===
    total_succ_dist = 0.0
    for vi in 1:n
        @inbounds for u in adj[vi]
            total_succ_dist += abs(vi - u)
        end
    end
    avg_succ_dist = total_edges > 0 ? total_succ_dist / total_edges : 0.0

    # === Interval statistics ===
    total_in_intervals = 0
    total_interval_count = 0
    for vi in 1:n
        nbrs = adj[vi]
        deg = length(nbrs)
        deg < mil && continue
        run_start = 1
        @inbounds for j in 2:deg
            if nbrs[j] != nbrs[j-1] + 1
                run_len = j - run_start
                if run_len >= mil
                    total_in_intervals += run_len
                    total_interval_count += 1
                end
                run_start = j
            end
        end
        run_len = deg - run_start + 1
        if run_len >= mil
            total_in_intervals += run_len
            total_interval_count += 1
        end
    end
    interval_coverage = total_edges > 0 ? Float64(total_in_intervals) / total_edges : 0.0
    avg_intervals = n > 0 ? Float64(total_interval_count) / n : 0.0

    # === Reference overlap statistics ===
    if ref_sample > 0 && ref_sample < n
        rng = Random.MersenneTwister(42)
        sample_indices = sort(Random.randperm(rng, n)[1:ref_sample])
    else
        sample_indices = collect(1:n)
    end

    total_best_overlap = 0
    total_best_jaccard = 0.0
    total_copy_fraction = 0.0
    ref_hits = 0
    sampled_with_neighbors = 0
    # Residual stats: after removing copied elements, measure residual encoding cost
    total_residual_log_gap = 0.0
    n_residual_gaps = 0
    total_residuals_in_intervals = 0
    total_residual_count = 0

    for vi in sample_indices
        nbrs = adj[vi]
        deg = length(nbrs)
        deg == 0 && continue
        sampled_with_neighbors += 1

        best_overlap = 0
        best_jaccard = 0.0
        best_ri = 0

        for ri in max(1, vi - window):(vi - 1)
            ref_nbrs = adj[ri]
            isempty(ref_nbrs) && continue
            overlap = _sorted_overlap_count(nbrs, ref_nbrs)
            if overlap > best_overlap
                best_overlap = overlap
                union_size = deg + length(ref_nbrs) - overlap
                best_jaccard = Float64(overlap) / union_size
                best_ri = ri
            end
        end

        total_best_overlap += best_overlap
        total_best_jaccard += best_jaccard
        total_copy_fraction += Float64(best_overlap) / deg
        if best_overlap > 0
            ref_hits += 1
        end

        # Compute residual stats
        if best_ri > 0 && best_overlap > 0
            residuals = _sorted_setdiff(nbrs, adj[best_ri])
        else
            residuals = nbrs
        end
        nr = length(residuals)
        total_residual_count += nr
        # Gap stats on residuals
        @inbounds for j in 2:nr
            gap = residuals[j] - residuals[j-1]
            total_residual_log_gap += log2(Float64(gap) + 1.0)
            n_residual_gaps += 1
        end
        # Interval coverage of residuals
        if nr >= mil
            run_start = 1
            @inbounds for j in 2:nr
                if residuals[j] != residuals[j-1] + 1
                    if j - run_start >= mil
                        total_residuals_in_intervals += j - run_start
                    end
                    run_start = j
                end
            end
            if nr - run_start + 1 >= mil
                total_residuals_in_intervals += nr - run_start + 1
            end
        end
    end

    sn = max(sampled_with_neighbors, 1)
    avg_residual_log_gap = n_residual_gaps > 0 ? total_residual_log_gap / n_residual_gaps : 0.0
    residual_interval_coverage = total_residual_count > 0 ? Float64(total_residuals_in_intervals) / total_residual_count : 0.0

    return Dict{Symbol, Any}(
        :n_vertices => n,
        :n_edges => total_edges,
        :window => window,
        :mil => mil,
        :ref_sample => ref_sample > 0 ? ref_sample : n,
        :avg_log_gap => avg_log_gap,
        :avg_successor_dist => avg_succ_dist,
        :avg_successor_dist_norm => avg_succ_dist / max(n, 1),
        :interval_coverage => interval_coverage,
        :avg_intervals_per_vertex => avg_intervals,
        :avg_best_overlap => Float64(total_best_overlap) / sn,
        :avg_best_jaccard => total_best_jaccard / sn,
        :ref_hit_rate => Float64(ref_hits) / sn,
        :avg_copy_fraction => total_copy_fraction / sn,
        :avg_residual_log_gap => avg_residual_log_gap,
        :residual_interval_coverage => residual_interval_coverage,
    )
end

"""
    print_ordering_metrics(metrics::Dict{Symbol, Any}; label::String="")

Pretty-print ordering quality metrics.
"""
function print_ordering_metrics(metrics::Dict{Symbol, Any}; label::String="")
    _pad(s, w) = s * " "^max(0, w - length(s))
    _rpad(s, w) = " "^max(0, w - length(s)) * s
    _f(v, d) = d == 1 ? string(round(v; digits=1)) :
               d == 3 ? string(round(v; digits=3)) :
               d == 4 ? string(round(v; digits=4)) :
               string(round(v; sigdigits=6))

    hdr = isempty(label) ? "Ordering Quality Metrics" : "Ordering Quality Metrics [$label]"
    println(hdr)
    println("-"^62)
    println("  $(_pad("vertices", 32)) $(metrics[:n_vertices])")
    println("  $(_pad("edges", 32)) $(metrics[:n_edges])")
    println("  $(_pad("window", 32)) $(metrics[:window])")
    println("  $(_pad("mil", 32)) $(metrics[:mil])")
    if haskey(metrics, :ref_sample)
        println("  $(_pad("ref_sample", 32)) $(metrics[:ref_sample])")
    end
    println("-"^62)
    println("  Gap metrics (lower = better locality, helps CS/BG):")
    println("    $(_pad("avg_log_gap", 30)) $(_rpad(_f(metrics[:avg_log_gap], 4), 10))")
    println("    $(_pad("avg_successor_dist", 30)) $(_rpad(_f(metrics[:avg_successor_dist], 1), 10))")
    println("    $(_pad("avg_succ_dist_norm", 30)) $(_rpad(_f(metrics[:avg_successor_dist_norm], 6), 10))")
    println("  Interval metrics (higher = more intervals, helps CS/BG):")
    ic = metrics[:interval_coverage]
    println("    $(_pad("interval_coverage", 30)) $(_rpad(_f(ic, 4), 10))  ($(round(100*ic; digits=1))%)")
    println("    $(_pad("avg_intervals/vertex", 30)) $(_rpad(_f(metrics[:avg_intervals_per_vertex], 4), 10))")
    println("  Reference metrics (higher = better ref copy, helps BV/CG):")
    println("    $(_pad("avg_best_overlap", 30)) $(_rpad(_f(metrics[:avg_best_overlap], 3), 10))")
    println("    $(_pad("avg_best_jaccard", 30)) $(_rpad(_f(metrics[:avg_best_jaccard], 4), 10))")
    rh = metrics[:ref_hit_rate]
    println("    $(_pad("ref_hit_rate", 30)) $(_rpad(_f(rh, 4), 10))  ($(round(100*rh; digits=1))%)")
    cf = metrics[:avg_copy_fraction]
    println("    $(_pad("avg_copy_fraction", 30)) $(_rpad(_f(cf, 4), 10))  ($(round(100*cf; digits=1))%)")
    println("  Residual metrics (after best-ref copy):")
    println("    $(_pad("avg_residual_log_gap", 30)) $(_rpad(_f(metrics[:avg_residual_log_gap], 4), 10))")
    ric = metrics[:residual_interval_coverage]
    println("    $(_pad("residual_interval_cov", 30)) $(_rpad(_f(ric, 4), 10))  ($(round(100*ric; digits=1))%)")
end

"""
    compare_ordering_metrics(m1::Dict{Symbol,Any}, m2::Dict{Symbol,Any}; label1::String="Order 1", label2::String="Order 2")

Print two sets of ordering metrics side by side for comparison.
Delta column shows m2 - m1.
"""
function compare_ordering_metrics(m1::Dict{Symbol,Any}, m2::Dict{Symbol,Any}; label1::String="Order 1", label2::String="Order 2")
    _rpad(s, w) = " "^max(0, w - length(s)) * s
    _f(v, d) = d == 1 ? string(round(v; digits=1)) :
               d == 3 ? string(round(v; digits=3)) :
               d == 4 ? string(round(v; digits=4)) :
               string(round(v; sigdigits=6))
    _pad(s, w) = s * " "^max(0, w - length(s))

    function _row(name, key, digits)
        v1 = Float64(m1[key])
        v2 = Float64(m2[key])
        d = v2 - v1
        s1 = _rpad(_f(v1, digits), 12)
        s2 = _rpad(_f(v2, digits), 12)
        sd = _rpad((d >= 0 ? "+" : "") * _f(d, digits), 12)
        println("  $(_pad(name, 28)) $s1 $s2 $sd")
    end

    println("\n  $(_pad("Metric", 28)) $(_rpad(label1, 12)) $(_rpad(label2, 12)) $(_rpad("Delta(2-1)", 12))")
    println("  ", "-"^66)
    println("  Gap metrics (lower = better locality, helps CS/BG):")
    _row("avg_log_gap", :avg_log_gap, 4)
    _row("avg_successor_dist", :avg_successor_dist, 1)
    _row("avg_succ_dist_norm", :avg_successor_dist_norm, 6)
    println("  Interval metrics (higher = more intervals, helps CS/BG):")
    _row("interval_coverage", :interval_coverage, 4)
    _row("avg_intervals/vertex", :avg_intervals_per_vertex, 4)
    println("  Reference metrics (higher = better ref copy, helps BV/CG):")
    _row("avg_best_overlap", :avg_best_overlap, 3)
    _row("avg_best_jaccard", :avg_best_jaccard, 4)
    _row("ref_hit_rate", :ref_hit_rate, 4)
    _row("avg_copy_fraction", :avg_copy_fraction, 4)
    println("  Residual metrics (after best-ref copy):")
    _row("avg_residual_log_gap", :avg_residual_log_gap, 4)
    _row("residual_interval_cov", :residual_interval_coverage, 4)
end

end # module
