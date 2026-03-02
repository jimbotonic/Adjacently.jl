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

using LightGraphs, DataStructures, Logging, SparseArrays, Random
using ..CustomTypes: UInt24, UInt40
using ..CustomLightGraphs: SimpleDiGraph, SimpleGraph, SimpleEdge
using ..Graph: get_in_degrees, get_out_degrees, get_in_out_degrees, get_reverse_graph
using ..PageRank: PR

export relabel_graph, relabel_vertices, relabel_vertices_score, relabel_vertices_lexicographic, relabel_vertices_rcm, relabel_vertices_webgraph_lex,
       relabel_vertices_llp, relabel_vertices_minhash, relabel_vertices_bisection,
       save_llp_ordering, load_llp_ordering

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

end # module
