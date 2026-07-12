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

module Clustering

using Random
using LightGraphs
using SparseArrays
using LinearAlgebra
using ..CustomLightGraphs: SimpleGraph, SimpleDiGraph
using LightGraphs: src, dst, AbstractGraph, outneighbors, edges, degree, nv, is_directed, add_vertices!, inneighbors

export heavy_edge_matching,
       contract_graph,
       initial_kway_partition,
       project_partition_up,
       refine_fm!,
       metis_partition,
       # Community detection
       louvain_partition,
       leiden_partition,
       auto_select_K,
       # AMG exports
       Smoother,
       Jacobi,
       AMGHierarchy,
       jacobi_setup,
       smooth!,
       strength_of_connection,
       greedy_aggregate,
       piecewise_constant_P,
       smoothed_P,
       amg_setup,
       vcycle!

# -----------------------------------------------------------------------------
# METIS
# -----------------------------------------------------------------------------

# ---------- Types ----------
"""
Lightweight wrapper to carry edge weights without a dedicated weighted graph type.
Provide weight(u,v) to all routines. Default weight = 1.0 for all edges.
"""
const WeightFun = Function

# ---------- Heavy-Edge Matching ----------
"""
    heavy_edge_matching(g; weight=(u,v)->1.0, rng=Random.default_rng())

Return:
- matches::Vector{Tuple{T,T}} : pairs (u,v), singletons as (u,u)
"""
function heavy_edge_matching(g::AbstractGraph{T}; weight::WeightFun = (u,v)->1.0, rng::AbstractRNG=Random.default_rng()) where {T<:Unsigned}
    n = nv(g)
    order = collect(1:n); shuffle!(rng, order)
    marked = falses(n)
    matches = Vector{Tuple{T,T}}()

    for u in order
        if marked[u]; continue; end
        best_found = false
        bestv = zero(T)
        bestw = -Inf
        for vT in outneighbors(g, u)
            if !marked[vT]
                w = weight(convert(T, u), vT)
                if w > bestw
                    bestw = w; bestv = vT; best_found = true
                end
            end
        end
        if !best_found
            push!(matches, (convert(T,u), convert(T,u)))
            marked[u] = true
        else
            push!(matches, (convert(T,u), bestv))
            marked[u] = true
            marked[bestv] = true
        end
    end
    return matches
end

# ---------- Contraction (quotient graph) ----------
"""
    contract_graph(g, matches; weight)

- Each (u,v) becomes a supernode (singleton allowed).
- Edges between supernodes are summed by weight.
Return:
- gc::SimpleGraph or SimpleDiGraph matching directedness of g
- coarse_w::Dict{Tuple{T,T},Float64} edge weights on gc
- coarse_of::Vector{Int} fine->coarse map
"""
function contract_graph(g::AbstractGraph{T}, matches::Vector{Tuple{T,T}}; weight::WeightFun=(u,v)->1.0) where {T<:Unsigned}
    n = nv(g)
    coarse_of = zeros(Int, n)

    # assign coarse ids
    cid = 0
    for (uT,vT) in matches
        cid += 1
        coarse_of[uT] = cid
        if uT != vT
            coarse_of[vT] = cid
        end
    end

    # Build weighted edges on coarse graph
    # Use Dict to accumulate weights
    wdict = Dict{Tuple{T,T}, Float64}()
    directed = is_directed(g)
    for e in edges(g)
        uT = src(e); vT = dst(e)
        cu = coarse_of[uT]; cv = coarse_of[vT]
        if cu == cv; continue; end
        if directed
            aT = convert(T, cu); bT = convert(T, cv)
            wdict[(aT,bT)] = get(wdict, (aT,bT), 0.0) + weight(uT, vT)
        else
            if cu <= cv
                aT = convert(T, cu); bT = convert(T, cv)
            else
                aT = convert(T, cv); bT = convert(T, cu)
            end
            wdict[(aT,bT)] = get(wdict, (aT,bT), 0.0) + weight(uT, vT)
        end
    end

    m = maximum(coarse_of)
    gc = directed ? SimpleDiGraph{T}() : SimpleGraph{T}()
    add_vertices!(gc, m)
    for ((aT,bT), _) in wdict
        add_edge!(gc, aT, bT)
    end
    return gc, wdict, coarse_of
end

# ---------- Initial k-way (greedy grow) ----------
"""
Greedy graph growing into k parts with capacity constraint.

Returns part::Vector{Int} of length nv(g) with labels 1..k
"""
function initial_kway_partition(g::AbstractGraph{T}, k::Int; imbalance_tol=0.03) where {T<:Unsigned}
    n = nv(g)
    target = ceil(Int, (1.0 + imbalance_tol) * n / k)
    part = zeros(Int, n)
    # Seeds: pick k highest degrees (or random if small)
    degs = degree(g)
    seeds = sortperm(degs, rev=true)[1:min(k,n)]
    for (i,s) in enumerate(seeds)
        part[s] = i
    end
    # BFS-like growth around seeds
    frontier = [T[] for _ in 1:k]
    for (i,s) in enumerate(seeds)
        push!(frontier[i], convert(T, s))
    end
    assigned = count(!=(0), part)
    while assigned < n
        for p in 1:k
            if sum(part .== p) >= target; continue; end
            if isempty(frontier[p])
                # find an unassigned vertex touching current part
                expand = T[]
                for u in 1:n
                    if part[u] == p
                        for vT in outneighbors(g, u)
                            if part[vT] == 0; push!(expand, vT); end
                        end
                    end
                end
                if isempty(expand)
                    # fallback: assign any unassigned vertex
                    v = findfirst(==(0), part)
                    v === nothing && break
                    part[v] = p; assigned += 1
                else
                    for v in expand
                        if part[v] == 0 && sum(part .== p) < target
                            part[v] = p; assigned += 1; push!(frontier[p], v)
                        end
                    end
                end
            else
                # expand from frontier
                newfront = T[]
                for u in frontier[p]
                    for v in outneighbors(g,u)
                        if part[v] == 0 && sum(part .== p) < target
                            part[v] = p; assigned += 1; push!(newfront, v)
                        end
                    end
                end
                frontier[p] = newfront
            end
            assigned >= n && break
        end
    end
    return part
end

# ---------- Project partition up ----------
"""
    project_partition_up(Pc, coarse_of): map coarse labels to fine vertices.

Pc is Vector{Int} of length nv(G_coarse), coarse_of maps fine->coarse
"""
project_partition_up(Pc::Vector{Int}, coarse_of::Vector{Int}) = Pc[coarse_of]

# ---------- Simple k-way FM refinement (outline + working heuristic) ----------
"""
    refine_fm!(g, part; imbalance_tol, max_passes)

A light, working heuristic inspired by FM:
- For each vertex, consider move to best neighboring part that reduces cut
- Enforce balance hard
"""
function refine_fm!(g::AbstractGraph{T}, part::Vector{Int}; imbalance_tol=0.03, max_passes=5) where {T<:Unsigned}
    n = nv(g); k = maximum(part)
    target = ceil(Int, (1.0 + imbalance_tol) * n / k)
    sizes = [count(==(p), part) for p in 1:k]

    function boundary_cost(u, pnew)
        # cut difference if u moves to pnew
        pold = part[u]
        gain = 0
        for v in outneighbors(g,u)
            if part[v] == pold; gain += 1
            elseif part[v] == pnew; gain -= 1
            end
        end
        return gain # positive → reduces cut
    end

    for _ in 1:max_passes
        improved = false
        for u in 1:n
            pold = part[u]
            bestp = pold
            bestgain = 0
            # try moving to any neighboring part
            for v in outneighbors(g,u)
                pnew = part[v]
                if pnew == pold; continue; end
                if sizes[pnew] + 1 > target; continue; end
                gain = boundary_cost(u, pnew)
                if gain > bestgain
                    bestgain = gain; bestp = pnew
                end
            end
            if bestp != pold && bestgain > 0
                part[u] = bestp
                sizes[pold] -= 1; sizes[bestp] += 1
                improved = true
            end
        end
        !improved && break
    end
    return part
end

# ---------- Top-level multilevel partition ----------
"""
    metis_partition(g, k; ...)

Returns final part::Vector{Int}
"""
function metis_partition(g::AbstractGraph{T}, k::Int;
    min_coarse_size::Int=1000, coarsen_ratio::Float64=0.7, max_levels::Int=50,
    imbalance_tol::Float64=0.03, rng::AbstractRNG=Random.default_rng(),
    weight::WeightFun=(u,v)->1.0) where {T<:Unsigned}

    # Build multilevel hierarchy
    graphs = typeof(g)[g]
    maps   = Vector{Int}[]           # fine->coarse per level
    weights = Dict{Tuple{T,T},Float64}[]

    level = 1
    while nv(graphs[level]) > min_coarse_size && level < max_levels
        matches = heavy_edge_matching(graphs[level]; weight, rng)
        gc, wdict, coarse_of = contract_graph(graphs[level], matches; weight)
        push!(graphs, gc)
        push!(maps, coarse_of)
        push!(weights, wdict)
        level += 1
        if nv(gc) > coarsen_ratio * nv(graphs[level-1])
            # stop if coarsening stalls
            break
        end
    end

    # Initial partition on coarsest
    P = initial_kway_partition(graphs[end], k; imbalance_tol)

    # Uncoarsen + refine
    for ℓ in reverse(1:length(maps))
        P = project_partition_up(P, maps[ℓ])
        refine_fm!(graphs[ℓ], P; imbalance_tol)
    end

    return P
end

# -----------------------------------------------------------------------------
# AMG (smoothed aggregation)
# -----------------------------------------------------------------------------

# ---------- Smoothers ----------
abstract type Smoother end

struct Jacobi <: Smoother
    omega::Float64
    Dinv::Vector{Float64}  # cached inverse diagonal per level (filled later)
end

function jacobi_setup(A::SparseMatrixCSC{Float64,T}, omega::Float64=2/3) where {T<:Unsigned}
    n = size(A,1)
    D = Vector{Float64}(undef, n)
    @inbounds for i in 1:n
        D[i] = A[i,i]
    end
    Dinv = map(x -> x == 0 ? 0.0 : 1.0/x, D)
    return Jacobi(omega, Dinv)
end

function smooth!(A::SparseMatrixCSC{Float64,T}, b::Vector{Float64}, x::Vector{Float64}, S::Jacobi) where {T<:Unsigned}
    # x ← x + ω D^{-1}(b - A x)
    r = b - A * x
    @inbounds @simd for i in eachindex(x)
        x[i] += S.omega * S.Dinv[i] * r[i]
    end
    return x
end

# ---------- Strength of connection ----------
"""
    strength_of_connection(A, θ)

i strongly depends on j if -A[i,j] ≥ θ * max_{k≠i}(-A[i,k]).
Returns a Sparse Bool pattern S (n×n) marking strong edges.
"""
function strength_of_connection(A::SparseMatrixCSC{Float64,T}, theta::Float64) where {T<:Unsigned}
    n = size(A,1)
    rows = Int[]; cols = Int[]
    @inbounds for i in 1:n
        # work on row i: get column indices and values for A[i,:]
        col_idx, vals = findnz(A[i, :])  # returns (indices, values) for SparseVector
        # compute alpha = max_{k≠i} (-A[i,k])
        alpha = 0.0
        for aij in vals
            alpha = max(alpha, -aij)
        end
        alpha == 0.0 && continue
        for (jT, aij) in zip(col_idx, vals)
            j = Int(jT)
            j == i && continue
            if -aij >= theta * alpha
                push!(rows, i); push!(cols, j)
            end
        end
    end
    S_int = sparse(rows, cols, ones(Bool,length(rows)), n, n)
    return SparseMatrixCSC{Bool,T}(S_int)
end

# ---------- Aggregation ----------
"""
    greedy_aggregate(S::SparseMatrixCSC{Bool})

Greedy seeds, then attach strong neighbors. Returns agg_id::Vector{Int} (fine -> coarse id).
"""
function greedy_aggregate(S::SparseMatrixCSC{Bool,T}) where {T<:Unsigned}
    n = size(S,1)
    unassigned = trues(n)
    agg_id = zeros(Int, n)
    cid = 0
    strong_deg = [nnz(S[i,:]) + nnz(S[:,i]) for i in 1:n]

    while any(unassigned)
        # pick seed with max strong degree among unassigned
        scores = map(i -> unassigned[i] ? strong_deg[i] : -1, 1:n)
        s = argmax(scores)
        if !unassigned[s]
            s = findfirst(unassigned)
        end
        cid += 1
        agg_id[s] = cid
        unassigned[s] = false
        # attach neighbors v if S[s,v] || S[v,s]
        Si = findnz(S[s,:])[1]  # column indices in row s (SparseVector)
        Sj = findnz(S[:,s])[1]  # row indices in column s (SparseVector)
        for v in union(Si, Sj)
            if unassigned[v]
                agg_id[v] = cid
                unassigned[v] = false
            end
        end
    end
    return agg_id, cid
end

# ---------- Interpolation ----------
"""
    piecewise_constant_P(agg_id, nAgg)
"""
function piecewise_constant_P(::Type{T}, agg_id::Vector{Int}, nAgg::Int) where {T<:Unsigned}
    n = length(agg_id)
    I = Int[]; J = Int[]; V = Float64[]
    @inbounds for i in 1:n
        c = agg_id[i]
        push!(I, i); push!(J, c); push!(V, 1.0)
    end
    Pint = sparse(I, J, V, n, nAgg)
    return SparseMatrixCSC{Float64,T}(Pint)
end

"""
    smoothed_P(A, P, ω=4/3) : P ← (I - ω D^{-1} A) P
"""
function smoothed_P(A::SparseMatrixCSC{Float64,T}, P::SparseMatrixCSC{Float64,T}, omega::Float64=1.0) where {T<:Unsigned}
    n = size(A,1)
    D = Vector{Float64}(undef, n)
    @inbounds for i in 1:n
        D[i] = A[i,i]
    end
    Dinv = map(x-> x==0 ? 0.0 : 1.0/x, D)
    AP = A * P
    rows, cols, vals = findnz(AP)
    @inbounds for k in eachindex(vals)
        vals[k] = Dinv[rows[k]] * vals[k]
    end
    AP_scaled_int = sparse(rows, cols, vals, n, size(P,2))
    AP_scaled = SparseMatrixCSC{Float64,T}(AP_scaled_int)
    # P_smoothed = P - omega * (D^{-1} A) P
    return P - omega * AP_scaled
end

# ---------- Hierarchy ----------
struct AMGHierarchy{T<:Unsigned}
    A::Vector{SparseMatrixCSC{Float64,T}}
    P::Vector{SparseMatrixCSC{Float64,T}}
    R::Vector{SparseMatrixCSC{Float64,T}}
    smoother::Vector{Jacobi}
end

"""
    amg_setup(A0; θ=0.08, max_levels=25, min_coarse=50, smoothed=true)

Builds SA-AMG hierarchy.
"""
function amg_setup(A0::SparseMatrixCSC{Float64,T};
                   theta::Float64=0.08, max_levels::Int=25, min_coarse::Int=50, smoothed::Bool=true) where {T<:Unsigned}
    A = SparseMatrixCSC{Float64,T}[A0]
    P = SparseMatrixCSC{Float64,T}[]
    R = SparseMatrixCSC{Float64,T}[]
    smoother = Jacobi[]
    push!(smoother, jacobi_setup(A0))

    level = 1
    while level <= max_levels && size(A[level],1) > min_coarse
        S = strength_of_connection(A[level], theta)
        agg_id, nAgg = greedy_aggregate(S)
        P_l = piecewise_constant_P(T, agg_id, nAgg)
        if smoothed
            P_l = smoothed_P(A[level], P_l)
        end
        R_l = transpose(P_l)
        Acoarse = R_l * A[level] * P_l

        push!(P, P_l); push!(R, R_l)
        push!(A, Acoarse)
        push!(smoother, jacobi_setup(Acoarse))
        level += 1
        if size(Acoarse,1) >= size(A[level-1],1) # stall guard
            break
        end
    end

    return AMGHierarchy{T}(A, P, R, smoother)
end

# ---------- V-cycle ----------
function vcycle!(H::AMGHierarchy{T}, b::Vector{Float64}, x::Vector{Float64};
                 nu_pre::Int=1, nu_post::Int=1, level::Int=1) where {T<:Unsigned}
    if level == length(H.A)
        # coarse solve (small) — use a direct factorization
        # Convert to Int-indexed CSC for UMFPACK
        Ac = SparseMatrixCSC{Float64,Int}(H.A[level])
        x[:] = Ac \ b
        return x
    end

    # pre-smooth
    for _ in 1:nu_pre
        smooth!(H.A[level], b, x, H.smoother[level])
    end

    # residual → coarse RHS
    r = b - H.A[level] * x
    b_c = H.R[level] * r

    # recursive coarse solve
    e_c = zeros(size(H.A[level+1],1))
    vcycle!(H, b_c, e_c; nu_pre=nu_pre, nu_post=nu_post, level=level+1)

    # correction and post-smooth
    x .+= H.P[level] * e_c
    for _ in 1:nu_post
        smooth!(H.A[level], b, x, H.smoother[level])
    end
    return x
end

# -----------------------------------------------------------------------------
# Louvain + Leiden (community detection)
# -----------------------------------------------------------------------------

"""
    modularity(g, part)

Compute modularity of a partition `part` (Int labels) for directed/undirected graphs.
Directed uses Leicht–Newman formulation.
"""
function modularity(g::AbstractGraph{T}, part::Vector{Int}) where {T<:Unsigned}
    n = nv(g)
    if is_directed(g)
        # directed modularity: Q = sum_c (e_c - a_out_c * a_in_c)
        m = 0.0
        kout = zeros(Float64, n)
        kin  = zeros(Float64, n)
        # count edges and degrees
        for u in 1:n
            du = 0
            for v in outneighbors(g, u)
                du += 1
                kin[v] += 1
            end
            kout[u] += du
            m += du
        end
        m == 0 && return 0.0
        # per-community sums
        maxc = maximum(part)
        e_c = zeros(Float64, maxc)
        aout = zeros(Float64, maxc)
        ain  = zeros(Float64, maxc)
        for c in 1:maxc
            aout[c] = sum(kout[i] for i in 1:n if part[i]==c) / m
            ain[c]  = sum(kin[i]  for i in 1:n if part[i]==c) / m
        end
        # count internal directed edges
        for u in 1:n
            cu = part[u]
            for v in outneighbors(g, u)
                if part[v] == cu
                    e_c[cu] += 1
                end
            end
        end
        for c in 1:maxc
            e_c[c] /= m
        end
        return sum(e_c[c] - aout[c]*ain[c] for c in 1:maxc)
    else
        # undirected modularity: Q = sum_c (e_c - a_c^2)
        m = 0.0
        deg = zeros(Float64, n)
        for u in 1:n
            du = length(outneighbors(g,u))
            deg[u] = du
            m += du
        end
        m /= 2
        m == 0 && return 0.0
        maxc = maximum(part)
        e_c = zeros(Float64, maxc)
        a_c = zeros(Float64, maxc)
        for c in 1:maxc
            a_c[c] = sum(deg[i] for i in 1:n if part[i]==c) / (2m)
        end
        # count internal edges once
        for u in 1:n
            cu = part[u]
            for v in outneighbors(g,u)
                if part[v] == cu && v > u
                    e_c[cu] += 1
                end
            end
        end
        for c in 1:maxc
            e_c[c] /= m
        end
        return sum(e_c[c] - a_c[c]^2 for c in 1:maxc)
    end
end

"""
    louvain_partition(g; max_passes=10)

Basic Louvain-style greedy modularity optimization (single-level).
Returns `part::Vector{Int}` with labels 1..k.
"""
function louvain_partition(g::AbstractGraph{T}; max_passes::Int=10, max_levels::Int=10) where {T<:Unsigned}
    # Multilevel Louvain with local moving and weighted coarsening
    parts = Vector{Vector{Int}}()
    # Level 0 local moving
    p0 = nothing
    if is_directed(g)
        p0 = _louvain_directed_optimized(g; max_passes=max_passes)
    else
        p0 = _louvain_undirected_optimized(g; max_passes=max_passes)
    end
    push!(parts, p0)
    # Build first coarse level
    Gc = aggregate_graph(g, p0)
    prev_n = nv(g)
    level = 1
    while level <= max_levels && Gc.n < prev_n
        # Local moving on coarse graph
        pc = louvain_local_move_coarse(Gc; max_passes=max_passes)
        push!(parts, pc)
        # Aggregate further
        prev_n = Gc.n
        Gc = aggregate_graph(Gc, pc)
        level += 1
    end
    # Compose labels back to original
    combined = parts[end]
    for ℓ in reverse(1:length(parts)-1)
        combined = project_partition_up(combined, parts[ℓ])
    end
    # Relabel compact 1..k
    labels = unique(combined)
    remap = Dict{Int,Int}(l => i for (i,l) in enumerate(sort(labels)))
    @inbounds for i in 1:length(combined)
        combined[i] = remap[combined[i]]
    end
    return combined
end

function _louvain_undirected_optimized(g::AbstractGraph{T}; max_passes::Int=10) where {T<:Unsigned}
    n = nv(g)
    part = collect(1:n)
    # degree and total weight
    deg = [length(outneighbors(g,u)) for u in 1:n]
    two_m = sum(deg)
    two_m == 0 && return ones(Int, n) # empty graph
    # community total degree Σ_tot (by community id)
    ctot = copy(deg)  # initially each node is its own community
    # main loop
    for _ in 1:max_passes
        improved = false
        for u in 1:n
            cu = part[u]
            k_i = deg[u]
            # remove u from its community temporarily
            ctot[cu] -= k_i
            # accumulate k_i_in per neighboring community
            neigh = Dict{Int,Int}()
            for v in outneighbors(g,u)
                cv = part[v]
                neigh[cv] = get(neigh, cv, 0) + 1
            end
            best_c = cu
            best_gain = 0.0
            for (c, k_i_in_c) in neigh
                # gain proportional (constant 1/(2m) omitted for comparison)
                gain = k_i_in_c - (k_i * ctot[c]) / two_m
                if gain > best_gain + 1e-12
                    best_gain = gain
                    best_c = c
                end
            end
            # move if positive gain and different
            if best_c != cu && best_gain > 0
                part[u] = best_c
                ctot[best_c] += k_i
                improved = true
            else
                # put back
                ctot[cu] += k_i
            end
        end
        !improved && break
    end
    # relabel compact 1..k
    labels = unique(part)
    remap = Dict{Int,Int}(l => i for (i,l) in enumerate(sort(labels)))
    @inbounds for i in 1:n
        part[i] = remap[part[i]]
    end
    return part
end

function _louvain_directed_optimized(g::AbstractGraph{T}; max_passes::Int=10) where {T<:Unsigned}
    n = nv(g)
    part = collect(1:n)
    # degrees
    k_out = [length(outneighbors(g,u)) for u in 1:n]
    k_in  = [length(inneighbors(g,u)) for u in 1:n]
    m = sum(k_out)
    m == 0 && return ones(Int, n)
    # community totals
    K_out = copy(k_out)
    K_in  = copy(k_in)
    for _ in 1:max_passes
        improved = false
        for u in 1:n
            cu = part[u]
            ku_out = k_out[u]
            ku_in  = k_in[u]
            # Temporarily remove u from its community totals
            K_out[cu] -= ku_out
            K_in[cu]  -= ku_in
            # accumulate weights to/from neighbor communities
            neigh_out = Dict{Int,Int}()
            for v in outneighbors(g,u)
                c = part[v]
                neigh_out[c] = get(neigh_out, c, 0) + 1
            end
            neigh_in = Dict{Int,Int}()
            for v in inneighbors(g,u)
                c = part[v]
                neigh_in[c] = get(neigh_in, c, 0) + 1
            end
            # current internal edges counts
            s_current = get(neigh_out, cu, 0) + get(neigh_in, cu, 0)
            best_c = cu
            best_gain = 0.0
            cand = union(keys(neigh_out), keys(neigh_in))
            for c in cand
                s_c = get(neigh_out, c, 0) + get(neigh_in, c, 0)
                gain = (s_c - s_current) / m - (ku_out * (K_in[c] - K_in[cu]) + ku_in * (K_out[c] - K_out[cu])) / (m*m)
                if gain > best_gain + 1e-12
                    best_gain = gain
                    best_c = c
                end
            end
            if best_c != cu && best_gain > 0
                part[u] = best_c
                K_out[best_c] += ku_out
                K_in[best_c]  += ku_in
                improved = true
            else
                # put back
                K_out[cu] += ku_out
                K_in[cu]  += ku_in
            end
        end
        !improved && break
    end
    # compact labels
    labs = unique(part)
    rem = Dict{Int,Int}(l => i for (i,l) in enumerate(sort(labs)))
    @inbounds for i in 1:n
        part[i] = rem[part[i]]
    end
    return part
end

# -------------------------------
# Weighted coarse graph structure
# -------------------------------
struct WeightedCoarseGraph
    directed::Bool
    n::Int
    out_w::Vector{Dict{Int,Float64}}
    in_w::Vector{Dict{Int,Float64}}
    kout::Vector{Float64}
    kin::Vector{Float64}
    m::Float64
end

function aggregate_graph(g::AbstractGraph{T}, part::Vector{Int}) where {T<:Unsigned}
    C = maximum(part)
    directed = is_directed(g)
    out_w = [Dict{Int,Float64}() for _ in 1:C]
    in_w  = [Dict{Int,Float64}() for _ in 1:C]
    m = 0.0
    if directed
        for e in edges(g)
            cu = part[src(e)]; cv = part[dst(e)]
            out_w[cu][cv] = get(out_w[cu], cv, 0.0) + 1.0
            in_w[cv][cu]  = get(in_w[cv],  cu, 0.0) + 1.0
            m += 1.0
        end
    else
        # add both directions
        for e in edges(g)
            u = Int(src(e)); v = Int(dst(e))
            cu = part[u]; cv = part[v]
            out_w[cu][cv] = get(out_w[cu], cv, 0.0) + 1.0
            in_w[cv][cu]  = get(in_w[cv],  cu, 0.0) + 1.0
            out_w[cv][cu] = get(out_w[cv], cu, 0.0) + 1.0
            in_w[cu][cv]  = get(in_w[cu],  cv, 0.0) + 1.0
            m += 2.0
        end
    end
    kout = [sum(values(out_w[c])) for c in 1:C]
    kin  = [sum(values(in_w[c]))  for c in 1:C]
    return WeightedCoarseGraph(directed, C, out_w, in_w, kout, kin, m)
end

function aggregate_graph(G::WeightedCoarseGraph, part::Vector{Int})
    C2 = maximum(part)
    out_w = [Dict{Int,Float64}() for _ in 1:C2]
    in_w  = [Dict{Int,Float64}() for _ in 1:C2]
    m = 0.0
    for u in 1:G.n
        su = part[u]
        for (v, w) in G.out_w[u]
            sv = part[v]
            out_w[su][sv] = get(out_w[su], sv, 0.0) + w
            in_w[sv][su]  = get(in_w[sv],  su, 0.0) + w
            m += w
        end
    end
    kout = [sum(values(out_w[c])) for c in 1:C2]
    kin  = [sum(values(in_w[c]))  for c in 1:C2]
    return WeightedCoarseGraph(G.directed, C2, out_w, in_w, kout, kin, m)
end

function louvain_local_move_coarse(G::WeightedCoarseGraph; max_passes::Int=10)
    n = G.n
    part = collect(1:n)
    # community totals start as node degrees
    K_out = copy(G.kout)
    K_in  = copy(G.kin)
    m = G.m
    for _ in 1:max_passes
        improved = false
        for u in 1:n
            cu = part[u]
            ku_out = G.kout[u]
            ku_in  = G.kin[u]
            K_out[cu] -= ku_out
            K_in[cu]  -= ku_in
            # candidate communities from neighbors
            neigh_out = Dict{Int,Float64}()
            for (v, w) in G.out_w[u]
                c = part[v]
                neigh_out[c] = get(neigh_out, c, 0.0) + w
            end
            neigh_in = Dict{Int,Float64}()
            for (v, w) in G.in_w[u]
                c = part[v]
                neigh_in[c] = get(neigh_in, c, 0.0) + w
            end
            s_current = get(neigh_out, cu, 0.0) + get(neigh_in, cu, 0.0)
            best_c = cu
            best_gain = 0.0
            cand = union(keys(neigh_out), keys(neigh_in))
            for c in cand
                s_c = get(neigh_out, c, 0.0) + get(neigh_in, c, 0.0)
                gain = (s_c - s_current) / m - (ku_out * (K_in[c] - K_in[cu]) + ku_in * (K_out[c] - K_out[cu])) / (m*m)
                if gain > best_gain + 1e-12
                    best_gain = gain
                    best_c = c
                end
            end
            if best_c != cu && best_gain > 0
                part[u] = best_c
                K_out[best_c] += ku_out
                K_in[best_c]  += ku_in
                improved = true
            else
                K_out[cu] += ku_out
                K_in[cu]  += ku_in
            end
        end
        !improved && break
    end
    # compact labels
    labs = unique(part)
    rem = Dict{Int,Int}(l => i for (i,l) in enumerate(sort(labs)))
    @inbounds for i in 1:n
        part[i] = rem[part[i]]
    end
    return part
end

"""
    leiden_partition(g; max_passes=10)

Lightweight Leiden-like refinement: run Louvain, split disconnected communities,
then a second Louvain pass on refined labels.
"""
function leiden_partition(g::AbstractGraph{T}; max_passes::Int=10, max_levels::Int=10) where {T<:Unsigned}
    parts = Vector{Vector{Int}}()
    # local moving on original
    p0 = louvain_local_move(g; max_passes=max_passes)
    p0r = refine_partition(g, p0)
    push!(parts, p0r)
    Gc = aggregate_graph(g, p0r)
    prev_n = nv(g)
    level = 1
    while level <= max_levels && Gc.n < prev_n
        pc = louvain_local_move(Gc; max_passes=max_passes)
        pcr = refine_partition(Gc, pc)
        push!(parts, pcr)
        prev_n = Gc.n
        Gc = aggregate_graph(Gc, pcr)
        level += 1
    end
    combined = parts[end]
    for ℓ in reverse(1:length(parts)-1)
        combined = project_partition_up(combined, parts[ℓ])
    end
    labs = unique(combined)
    rem = Dict{Int,Int}(l => i for (i,l) in enumerate(sort(labs)))
    @inbounds for i in 1:length(combined)
        combined[i] = rem[combined[i]]
    end
    return combined
end

function louvain_local_move(g::AbstractGraph{T}; max_passes::Int=10) where {T<:Unsigned}
    if is_directed(g)
        return _louvain_directed_optimized(g; max_passes=max_passes)
    else
        return _louvain_undirected_optimized(g; max_passes=max_passes)
    end
end

function louvain_local_move(G::WeightedCoarseGraph; max_passes::Int=10)
    return louvain_local_move_coarse(G; max_passes=max_passes)
end

function refine_partition(g::AbstractGraph{T}, part::Vector{Int}) where {T<:Unsigned}
    n = nv(g)
    current = copy(part)
    lab = maximum(current)
    for c in 1:lab
        verts = [i for i in 1:n if current[i]==c]
        isempty(verts) && continue
        cid = 0
        visited = falses(n)
        comp_id = Dict{Int,Int}()
        for s in verts
            if visited[s]; continue; end
            cid += 1
            q = [s]; visited[s] = true; comp_id[s] = cid
            while !isempty(q)
                u = pop!(q)
                if is_directed(g)
                    for v in outneighbors(g,u)
                        if current[v]==c && !visited[v]
                            visited[v] = true; comp_id[v] = cid; push!(q,v)
                        end
                    end
                    for v in inneighbors(g,u)
                        if current[v]==c && !visited[v]
                            visited[v] = true; comp_id[v] = cid; push!(q,v)
                        end
                    end
                else
                    for v in outneighbors(g,u)
                        if current[v]==c && !visited[v]
                            visited[v] = true; comp_id[v] = cid; push!(q,v)
                        end
                    end
                end
            end
        end
        if cid > 1
            for (v,cc) in comp_id
                if cc > 1
                    current[v] = lab + cc - 1
                end
            end
            lab = maximum(current)
        end
    end
    return current
end

function refine_partition(G::WeightedCoarseGraph, part::Vector{Int})
    n = G.n
    current = copy(part)
    lab = maximum(current)
    for c in 1:lab
        verts = [i for i in 1:n if current[i]==c]
        isempty(verts) && continue
        cid = 0
        visited = falses(n)
        comp_id = Dict{Int,Int}()
        for s in verts
            if visited[s]; continue; end
            cid += 1
            q = [s]; visited[s] = true; comp_id[s] = cid
            while !isempty(q)
                u = pop!(q)
                for (v, _) in G.out_w[u]
                    if current[v]==c && !visited[v]
                        visited[v] = true; comp_id[v] = cid; push!(q,v)
                    end
                end
                for (v, _) in G.in_w[u]
                    if current[v]==c && !visited[v]
                        visited[v] = true; comp_id[v] = cid; push!(q,v)
                    end
                end
            end
        end
        if cid > 1
            for (v,cc) in comp_id
                if cc > 1
                    current[v] = lab + cc - 1
                end
            end
            lab = maximum(current)
        end
    end
    return current
end

"""
    auto_select_K(g; max_K=8, min_community_frac=0.005, min_granularity=0.01, partition=nothing)

Automatically select the optimal number of CG clusters K.

First runs a cheap single-pass Louvain to estimate community granularity (L/n).
If the graph has coarse community structure (L/n < `min_granularity`), returns K=1
immediately without running full Leiden — K=1 with global LLP is preferred for such
graphs. Otherwise, runs full Leiden and selects K based on graph size and community
structure:
- K grows logarithmically with vertex count: `ceil(log2(n / 100_000))`
  (graphs ≤ 100K vertices get K=1; larger graphs benefit from more clusters)
- K is capped by the number of significant Leiden communities (those with
  ≥ `min_community_frac` fraction of total vertices)

Returns `(K, partition)`:
- `K`: recommended number of clusters for `leiden_partition_k`
- `partition`: the Leiden partition vector (for reuse, avoids recomputing);
  `nothing` when K=1 (partition not needed for K=1 path)
"""
function auto_select_K(g::AbstractGraph{T};
        max_K::Int=8,
        min_community_frac::Float64=0.005,
        min_granularity::Float64=0.01,
        partition::Union{Nothing,Vector{Int}}=nothing) where {T<:Unsigned}
    n = nv(g)

    # Quick granularity check: 2-pass Louvain to estimate community count.
    # 2 passes converge enough to distinguish fine-grained graphs (many small communities)
    # from coarse ones (few large communities), while being much cheaper than full Leiden
    # (2 passes vs 8 passes + refinement + multi-level aggregation).
    quick_part = louvain_local_move(g; max_passes=2)
    L_quick = length(Set(quick_part))
    granularity = L_quick / n

    # Community granularity: L/n = communities per vertex.
    # Fine-grained (e.g. 0.08–0.10): many small tight communities;
    #   K>1 splitting groups them effectively and crawl-order locality within clusters is good.
    # Coarse (e.g. < 0.01): few large communities;
    #   K>1 splitting is ineffective, K=1 with global LLP is preferred.
    if granularity < min_granularity
        @info "auto_select_K: n=$n, ~$L_quick communities (quick), granularity=$(round(granularity, digits=4)) < $(min_granularity) (coarse) → K=1"
        return 1, nothing
    end

    # Fine-grained: run full Leiden for accurate partition (needed for K>1 path)
    part = partition !== nothing ? partition : leiden_partition(g; max_passes=8, max_levels=5)

    # Count community sizes and sort descending
    counts = Dict{Int,Int}()
    for c in part; counts[c] = get(counts, c, 0) + 1; end
    sorted_labels = sort(collect(keys(counts)); by=l->counts[l], rev=true)
    L = length(sorted_labels)

    # Graph-size baseline: K grows logarithmically with n.
    # Graphs ≤ 100K vertices → K=1 (one cluster suffices for reference window).
    # Larger graphs benefit from more clusters as community diversity increases.
    K_base = n > 100_000 ? max(1, ceil(Int, log2(n / 100_000))) : 1

    # Cap by number of significant communities (≥ min_community_frac of vertices)
    min_size = max(1, ceil(Int, n * min_community_frac))
    n_significant = count(l -> counts[l] >= min_size, sorted_labels)

    K = clamp(K_base, 1, min(max_K, n_significant, L))

    # Log top community sizes for diagnostics
    n_show = min(K + 2, length(sorted_labels))
    top_sizes = [counts[sorted_labels[i]] for i in 1:n_show]
    @info "auto_select_K: n=$n, $L communities, granularity=$(round(L / n, digits=4)), $n_significant significant (≥$(min_size)v), K_base=$K_base → K=$K  top_sizes=$top_sizes"

    return K, part
end

end # module Clustering
