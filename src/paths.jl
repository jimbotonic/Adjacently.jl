#
# Adjacently: Julia Complex Directed Networks Library
# Copyright (C) 2016-2026 Jimmy Dubuisson <jimmy@dubuisson.ch>
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

module Paths

using LightGraphs, DataStructures, Logging
using ..CustomTypes: UInt24, UInt40
using ..CustomLightGraphs: SimpleDiGraph, SimpleGraph, SimpleEdge
using ..PageRank: PPR
using ..Graph: get_reverse_graph

export ppr_path_scores, extract_path, extract_subgraph, bfs_shortest_path,
       discover_paths, discover_paths_bottleneck, discover_paths_fp,
       wppr, wppr_path_scores

"""
    ppr_path_scores(g::AbstractGraph{T}, s::T, t::T; damping=0.85, epsilon=1e-6) where {T<:Unsigned}

Compute bidirectional Personalized PageRank vectors that "illuminate" paths from `s` to `t`.

Returns `(π_s, π_t)` where:
- `π_s`: PPR from source `s` on `g` (high values = reachable from s)
- `π_t`: PPR from target `t` on the reverse of `g` (high values = can reach t)

The element-wise product `π_s .* π_t` scores vertices by how strongly they lie on paths from `s` to `t`.
"""
function ppr_path_scores(g::AbstractGraph{T}, s::T, t::T;
                         damping::Float64=0.85, epsilon::Float64=1e-6) where {T<:Unsigned}
    rg = get_reverse_graph(g)
    return ppr_path_scores(g, rg, s, t; damping=damping, epsilon=epsilon)
end

"""
    ppr_path_scores(g::AbstractGraph{T}, rg::AbstractGraph{T}, s::T, t::T; damping=0.85, epsilon=1e-6)

Same as above but accepts a precomputed reverse graph `rg` to avoid recomputing it.
"""
function ppr_path_scores(g::AbstractGraph{T}, rg::AbstractGraph{T}, s::T, t::T;
                         damping::Float64=0.85, epsilon::Float64=1e-6) where {T<:Unsigned}
    @info "Computing forward PPR from source $s"
    π_s = PPR(s, g, rg; damping=damping, epsilon=epsilon)

    @info "Computing backward PPR from target $t"
    π_t = PPR(t, rg, g; damping=damping, epsilon=epsilon)

    return (π_s, π_t)
end

"""
    extract_path(g::AbstractGraph{T}, s::T, t::T, scores::Vector{Float64}; max_hops::Int=1000) where {T<:Unsigned}

Greedily extract a path from `s` to `t` by following the highest-scoring outneighbor at each step.

Returns a vector of vertex IDs forming the path (including `s` and `t`), or an empty vector if
no path is found within `max_hops` steps.
"""
function extract_path(g::AbstractGraph{T}, s::T, t::T, scores::Vector{Float64};
                      max_hops::Int=1000) where {T<:Unsigned}
    path = T[s]
    visited = Set{T}(s)
    current = s

    for _ in 1:max_hops
        current == t && return path

        neighbors = outneighbors(g, current)
        isempty(neighbors) && return T[]  # dead end

        # pick unvisited neighbor with highest score
        best_score = -1.0
        best_next = zero(T)
        for n in neighbors
            if n ∉ visited && scores[n] > best_score
                best_score = scores[n]
                best_next = n
            end
        end

        best_next == zero(T) && return T[]  # all neighbors visited, stuck
        push!(path, best_next)
        push!(visited, best_next)
        current = best_next
    end

    return T[]  # max_hops exceeded
end

"""
    extract_subgraph(g::AbstractGraph{T}, s::T, t::T, scores::Vector{Float64}; threshold::Float64=0.0) where {T<:Unsigned}

Extract the induced subgraph on vertices whose score exceeds `threshold`, always including `s` and `t`.

The `threshold` controls which vertices are kept: only those with `score > threshold` are included.
If `threshold` is 0.0 (default), uses the mean of nonzero scores as an automatic threshold.

Returns `(subg, vmap)` where:
- `subg` is a new `SimpleDiGraph` on the high-scoring vertices
- `vmap` maps subgraph vertex IDs back to original graph vertex IDs
"""
function extract_subgraph(g::AbstractGraph{T}, s::T, t::T, scores::Vector{Float64};
                          threshold::Float64=0.0) where {T<:Unsigned}
    # auto-threshold: mean of nonzero scores
    if threshold == 0.0
        nonzero = filter(x -> x > 0.0, scores)
        threshold = isempty(nonzero) ? 0.0 : sum(nonzero) / length(nonzero)
    end

    # select vertices above threshold, always including s and t
    vset = Set{T}(T(v) for v in 1:length(scores) if scores[v] > threshold)
    push!(vset, s)
    push!(vset, t)
    vmap = sort!(collect(vset))

    # build reverse mapping: original ID → subgraph ID
    id_map = Dict{T,T}()
    for (i, v) in enumerate(vmap)
        id_map[v] = T(i)
    end

    # build induced subgraph
    subg = SimpleDiGraph{T}()
    add_vertices!(subg, length(vmap))
    for (i, v) in enumerate(vmap)
        for n in outneighbors(g, v)
            if haskey(id_map, n)
                add_edge!(subg, T(i), id_map[n])
            end
        end
    end

    return (subg, vmap)
end

"""
    bfs_shortest_path(g::AbstractGraph{T}, s::T, t::T) where {T<:Unsigned}

Find the shortest path from `s` to `t` in an unweighted directed graph using BFS.

Returns a vector of vertex IDs forming the path (including `s` and `t`), or an empty vector
if no path exists.
"""
function bfs_shortest_path(g::AbstractGraph{T}, s::T, t::T) where {T<:Unsigned}
    s == t && return T[s]

    visited = Set{T}(s)
    parent = Dict{T,T}()
    queue = Deque{T}()
    push!(queue, s)

    while !isempty(queue)
        current = popfirst!(queue)
        for n in outneighbors(g, current)
            if n ∉ visited
                parent[n] = current
                if n == t
                    # reconstruct path
                    path = T[t]
                    v = t
                    while v != s
                        v = parent[v]
                        pushfirst!(path, v)
                    end
                    return path
                end
                push!(visited, n)
                push!(queue, n)
            end
        end
    end

    return T[]  # no path found
end

"""
    discover_paths(g::AbstractGraph{T}, s::T, t::T, scores::Vector{Float64};
                   max_paths::Int=10, decay::Float64=0.1, max_hops::Int=1000) where {T<:Unsigned}

Discover multiple alternative paths from `s` to `t` by iteratively penalizing
previously discovered vertices.

Unlike naive zeroing (which deletes all vertices of a found path), this function
multiplies the scores of intermediate vertices by `decay` (default 0.1) after each
extraction. This means alternative paths that share most of their structure with
an earlier path — differing only in a short detour — can still be found.

Returns a vector of paths (each a `Vector{T}`), ordered by discovery. Duplicate
paths are skipped.
"""
function discover_paths(g::AbstractGraph{T}, s::T, t::T, scores::Vector{Float64};
                        max_paths::Int=10, decay::Float64=0.1,
                        max_hops::Int=1000) where {T<:Unsigned}
    current_scores = copy(scores)
    paths = Vector{T}[]
    seen = Set{Vector{T}}()

    for _ in 1:max_paths * 3  # extra iterations to account for duplicates
        length(paths) >= max_paths && break

        path = extract_path(g, s, t, current_scores; max_hops=max_hops)
        isempty(path) && break

        if path ∉ seen
            push!(seen, path)
            push!(paths, path)
        end

        # penalize intermediate vertices (keep s and t at full strength)
        for v in path
            if v != s && v != t
                current_scores[v] *= decay
            end
        end
    end

    return paths
end

"""
    discover_paths_bottleneck(g::AbstractGraph{T}, s::T, t::T, scores::Vector{Float64};
                              max_paths::Int=10, max_hops::Int=1000) where {T<:Unsigned}

Discover multiple alternative paths from `s` to `t` by iteratively removing
the bottleneck edge of each discovered path.

For each found path, the bottleneck edge is the edge `(u, v)` with the lowest
edge score `min(scores[u], scores[v])`. Removing only this single edge forces
the next extraction to find a detour around the weakest point while preserving
the rest of the shared structure.

This is more precise than vertex penalization: if two paths share 8 of 10
vertices and differ only in a 2-vertex detour, this method will find both.

Returns a vector of `(path, removed_edge)` tuples, where `removed_edge` is
the `(u, v)` pair that was cut after finding that path (`(0, 0)` for the last path).
"""
function discover_paths_bottleneck(g::AbstractGraph{T}, s::T, t::T, scores::Vector{Float64};
                                    max_paths::Int=10, max_hops::Int=1000) where {T<:Unsigned}
    # work on a copy so we don't modify the original graph
    gc = SimpleDiGraph{T}()
    add_vertices!(gc, nv(g))
    for e in edges(g)
        add_edge!(gc, src(e), dst(e))
    end

    results = Tuple{Vector{T}, Tuple{T,T}}[]
    seen = Set{Vector{T}}()

    for _ in 1:max_paths * 5  # extra budget for skipped duplicates
        length(results) >= max_paths && break

        path = extract_path(gc, s, t, scores; max_hops=max_hops)
        isempty(path) && break
        length(path) < 2 && break  # s == t, nothing to cut

        # find the bottleneck edge: lowest min(scores[u], scores[v])
        worst_score = Inf
        worst_idx = 1
        for i in 1:length(path)-1
            edge_score = min(scores[path[i]], scores[path[i+1]])
            if edge_score < worst_score
                worst_score = edge_score
                worst_idx = i
            end
        end
        bottleneck = (path[worst_idx], path[worst_idx+1])

        # record the path if new
        if path ∉ seen
            push!(seen, path)
            push!(results, (path, bottleneck))
        end

        # always remove the bottleneck edge (even for duplicates, to make progress)
        rem_edge!(gc, bottleneck[1], bottleneck[2])
    end

    return results
end

"""
    discover_paths_fp(g::AbstractGraph{T}, rg::AbstractGraph{T}, s::T, t::T;
                      damping=0.85, epsilon=1e-6, d::Int=0, max_paths::Int=10,
                      max_hops::Int=1000,
                      edge_weights::Union{Nothing,Dict{Tuple{T,T},Float64}}=nothing,
                      hub_penalty_alpha::Float64=0.0,
                      initial_path::Union{Nothing,Vector{T}}=nothing) where {T<:Unsigned}

Discover multiple paths from `s` to `t` using per-branch PPR fingerprints.

Algorithm:
1. Compute backward fingerprint fp_t = PPR(t, rg, g) once.
2. Compute forward fingerprint fp_s = PPR(s, g, rg), form fp_st = fp_s .* fp_t.
3. Find initial greedy path of length n. Mark its vertices as known.
4. Walk along the initial path. At each vertex v_i (position i from s):
   - Enumerate children of v_i, sorted by fp_st score descending, known children last.
   - For each non-known child c (stop when hitting the known child on the initial path):
     * Compute a fresh fp_c = PPR(c, g, rg) and fp_ct = fp_c .* fp_t
     * Greedy-extract a path from c to t using fp_ct
     * If path length ≤ n - i + d, record the full path s→...→v_i→c→...→t
     * Mark new vertices as known
5. Parameter `d` controls slack: d=0 means alternatives must be ≤ n steps total,
   d=1 allows 1 extra step, etc.

If `edge_weights` is provided, weighted PPR (`wppr`) is used in place of unweighted
PPR for all forward / backward / per-branch fingerprints, biasing the walk toward
high-weight edges (e.g. FBA-flux-weighted reactions). Reverse weights are derived
automatically.

`hub_penalty_alpha` (default 0.0 = no change) divides the score vector used for
greedy extraction by `deg(v)^alpha`, where `deg(v) = outdegree(g, v) + outdegree(rg, v)`.
This combats a pathology of weighted PPR on heavy-tailed graphs: high-flux edges
concentrate mass on hub metabolites, and the greedy walk gets stuck in the dense
amino-acid / central-carbon hub region instead of reaching the target. `alpha=1.0`
gives a Croes-style (2005) hub penalty; 0.5 is a gentler mix.

`initial_path` (default nothing = greedy-extract from `fp_st`) lets the caller supply
a pre-computed source→target skeleton. Use this to seed FP-branch with a high-quality
backbone (e.g. Yen's 1-shortest path under a chemistry-aware cost) and let the
per-branch fingerprints enumerate alternatives around it. The supplied path is used
verbatim — no validity checking — so the caller is responsible for making sure it
starts at `s`, ends at `t`, and only traverses existing edges.

Returns a vector of paths (each a `Vector{T}`).
"""
function discover_paths_fp(g::AbstractGraph{T}, rg::AbstractGraph{T}, s::T, t::T;
                           damping::Float64=0.85, epsilon::Float64=1e-6,
                           d::Int=0, max_paths::Int=10,
                           max_hops::Int=1000,
                           edge_weights::Union{Nothing,Dict{Tuple{T,T},Float64}}=nothing,
                           hub_penalty_alpha::Float64=0.0,
                           initial_path::Union{Nothing,Vector{T}}=nothing) where {T<:Unsigned}
    # Pick (weighted or unweighted) forward / backward PPR closures.
    # Anonymous lambdas (not `function` defs) so they close over the right locals.
    local fwd_ppr, bwd_ppr
    if edge_weights === nothing
        fwd_ppr = src -> PPR(src, g, rg; damping=damping, epsilon=epsilon)
        bwd_ppr = src -> PPR(src, rg, g; damping=damping, epsilon=epsilon)
    else
        rev_weights = Dict{Tuple{T,T},Float64}()
        for ((u, v), w) in edge_weights
            rev_weights[(v, u)] = w
        end
        fwd_ppr = src -> wppr(src, g, rg, edge_weights; damping=damping, epsilon=epsilon)
        bwd_ppr = src -> wppr(src, rg, g, rev_weights; damping=damping, epsilon=epsilon)
    end

    # Optional hub-penalty multiplier: 1 / deg(v)^alpha (or ones if alpha == 0).
    # deg(v) = outdegree(g, v) + outdegree(rg, v) (i.e. total connectivity).
    # Precomputed once; applied to every score vector used for greedy extraction.
    hub_mult = if hub_penalty_alpha > 0.0
        dv = [Float64(length(outneighbors(g, v)) + length(outneighbors(rg, v)))
              for v in T(1):T(nv(g))]
        # guard against deg=0 (isolated vertex), yields a zero multiplier which
        # also makes that vertex unreachable via greedy — fine since it has no edges.
        [x > 0 ? 1.0 / x^hub_penalty_alpha : 0.0 for x in dv]
    else
        nothing
    end

    apply_penalty(scores) = hub_mult === nothing ? scores : scores .* hub_mult

    # Step 1: compute backward fingerprint (fixed for all branches)
    fp_t = bwd_ppr(t)

    # Forward fingerprint from s
    fp_s = fwd_ppr(s)
    fp_st = apply_penalty(fp_s .* fp_t)

    # Initial path: caller-supplied seed (e.g. Yen's 1-shortest with a chemistry cost),
    # or the default greedy walk on fp_st when no seed is provided. Seeding with an
    # externally-computed skeleton lets us decouple the quality of the backbone from
    # the quality of the per-branch fingerprint scoring.
    local initial_path_used::Vector{T}
    if initial_path !== nothing
        initial_path_used = initial_path
    else
        initial_path_used = extract_path(g, s, t, fp_st; max_hops=max_hops)
    end
    if isempty(initial_path_used)
        @warn "No initial path found from $s to $t"
        return Vector{T}[]
    end
    initial_path = initial_path_used

    n = length(initial_path) - 1  # number of edges
    known = Set{T}(initial_path)
    paths = [initial_path]
    seen = Set{Vector{T}}([initial_path])

    @info "Initial path: $n steps, exploring branches (d=$d, budget=$n+$d=$(n+d))"

    # Steps 2-4: walk along the initial path, explore branches at each vertex
    for (i, v) in enumerate(initial_path[1:end-1])  # positions 1..n (don't explore from t)
        budget = n - i + d  # max allowed length for child→t subpath
        budget <= 0 && continue

        # Children of v, sorted: non-known first by fp_st desc, then known by fp_st desc
        children = collect(outneighbors(g, v))
        isempty(children) && continue

        sort!(children, by=c -> (c ∈ known ? 1 : 0, -fp_st[c]))

        # The known child on the initial path at this position
        next_on_path = initial_path[i + 1]

        for c in children
            length(paths) >= max_paths && break

            # Stop when we reach the known child from the initial path
            c == next_on_path && break

            # Skip if c is already known (they're sorted last, so this also stops)
            c ∈ known && break

            # Compute fresh fingerprint from c
            fp_c = fwd_ppr(c)
            fp_ct = apply_penalty(fp_c .* fp_t)

            # Greedy path from c to t using the branch-specific scores
            c_path = extract_path(g, c, t, fp_ct; max_hops=max_hops)
            isempty(c_path) && continue

            c_path_len = length(c_path) - 1
            c_path_len > budget && continue

            # Build full path: initial_path[1:i] → c_path
            # initial_path[1:i] = s, v1, ..., v_{i-1} = the prefix up to v
            full_path = vcat(initial_path[1:i], c_path)

            # Reject paths with revisited vertices (spur path re-entered the prefix)
            length(Set(full_path)) == length(full_path) || continue

            if full_path ∉ seen
                push!(seen, full_path)
                push!(paths, full_path)
                union!(known, Set(full_path))
            end
        end

        length(paths) >= max_paths && break
    end

    return paths
end

"""
    wppr(src::T, g::AbstractGraph{T}, rg::AbstractGraph{T},
         edge_weights::Dict{Tuple{T,T}, Float64};
         damping=0.85, epsilon=1e-6) where {T<:Unsigned}

Weighted Personalized PageRank. Like PPR, but the random walk follows edges
with probability proportional to their weight instead of uniformly.

`edge_weights` maps `(u, v) => weight` for edges in `g`.
Edges not in the dict get weight 1.0 (unweighted fallback).
"""
function wppr(src::T, g::AbstractGraph{T}, rg::AbstractGraph{T},
              edge_weights::Dict{Tuple{T,T}, Float64};
              damping::Float64=0.85, epsilon::Float64=1e-6) where {T<:Unsigned}
    n = nv(g)
    pr = zeros(Float64, n)
    pr[src] = 1.0
    pr2 = zeros(Float64, n)

    # precompute weighted out-degrees
    w_outdeg = zeros(Float64, n)
    for v in vertices(g)
        for c in outneighbors(g, v)
            w_outdeg[v] += get(edge_weights, (v, c), 1.0)
        end
    end

    while true
        for v in vertices(g)
            nv_val = 0.0
            # in-neighbors of v = out-neighbors of v in rg
            for p in outneighbors(rg, v)
                w_outdeg[p] > 0 || continue
                w_pv = get(edge_weights, (p, v), 1.0)
                nv_val += pr[p] * w_pv / w_outdeg[p]
            end
            pr2[v] = v == src ? (1 - damping) + damping * nv_val : damping * nv_val
        end

        maxdiff = 0.0
        for i in 1:n
            maxdiff = max(maxdiff, abs(pr2[i] - pr[i]))
        end
        if maxdiff <= epsilon
            return pr2
        end
        pr, pr2 = pr2, pr  # swap buffers
    end
end

"""
    wppr_path_scores(g::AbstractGraph{T}, rg::AbstractGraph{T}, s::T, t::T,
                     edge_weights::Dict{Tuple{T,T}, Float64};
                     damping=0.85, epsilon=1e-6) where {T<:Unsigned}

Compute bidirectional weighted PPR vectors for path discovery.

Returns `(π_s, π_t)` where:
- `π_s`: weighted PPR from `s` on `g` (forward, using `edge_weights`)
- `π_t`: weighted PPR from `t` on `rg` (backward, using reverse weights)

The reverse weights are computed automatically: `rw[(v,u)] = edge_weights[(u,v)]`.
"""
function wppr_path_scores(g::AbstractGraph{T}, rg::AbstractGraph{T}, s::T, t::T,
                          edge_weights::Dict{Tuple{T,T}, Float64};
                          damping::Float64=0.85, epsilon::Float64=1e-6) where {T<:Unsigned}
    # Build reverse weights for the backward PPR
    rev_weights = Dict{Tuple{T,T}, Float64}()
    for ((u, v), w) in edge_weights
        rev_weights[(v, u)] = w
    end

    @info "Computing forward weighted PPR from source $s"
    π_s = wppr(s, g, rg, edge_weights; damping=damping, epsilon=epsilon)

    @info "Computing backward weighted PPR from target $t"
    π_t = wppr(t, rg, g, rev_weights; damping=damping, epsilon=epsilon)

    return (π_s, π_t)
end

end # module Paths
