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

"""
    Generators

Random directed graph generators for training GNN-based compression parameter
prediction. All generators produce `SimpleDiGraph{UInt32}`.
"""
module Generators

using LightGraphs: nv, ne, add_edge!, add_vertices!, vertices, outneighbors,
                   watts_strogatz, edges, src, dst
using Random: MersenneTwister, rand, randperm, shuffle!

using ..CustomLightGraphs: SimpleDiGraph

export random_erdos_renyi_digraph, random_barabasi_albert_digraph,
       random_sbm_digraph, random_modular_hub_digraph,
       random_lfr_digraph, random_watts_strogatz_digraph,
       random_web_digraph, generate_training_batch

"""
    random_erdos_renyi_digraph(n, p; seed=42)

Erdős–Rényi directed graph: each directed edge (i→j), i≠j, exists with
probability `p`. Expected edges: n*(n-1)*p.
"""
function random_erdos_renyi_digraph(n::Int, p::Float64; seed::Int=42)
    rng = MersenneTwister(seed)
    g = SimpleDiGraph{UInt32}(UInt32(n))
    for i in UInt32(1):UInt32(n)
        for j in UInt32(1):UInt32(n)
            i == j && continue
            if rand(rng) < p
                add_edge!(g, i, j)
            end
        end
    end
    return g
end

"""
    random_barabasi_albert_digraph(n, m; seed=42)

Barabási–Albert preferential attachment directed graph. Starts with a clique of
`m+1` vertices, then each new vertex attaches to `m` existing vertices chosen
proportional to in-degree+1. Edges are directed from new → existing.
"""
function random_barabasi_albert_digraph(n::Int, m::Int; seed::Int=42)
    @assert m >= 1 && n > m
    rng = MersenneTwister(seed)
    g = SimpleDiGraph{UInt32}(UInt32(n))

    # Initial clique of m+1 vertices
    init = m + 1
    for i in UInt32(1):UInt32(init)
        for j in UInt32(1):UInt32(init)
            i == j && continue
            add_edge!(g, i, j)
        end
    end

    # In-degree tracker (for preferential attachment)
    in_deg = zeros(Int, n)
    for v in UInt32(1):UInt32(init)
        in_deg[Int(v)] = init - 1  # each vertex in clique has in-degree (init-1)
    end

    # Repeated edge list for fast weighted sampling
    targets = Int[]
    for v in 1:init
        for _ in 1:(in_deg[v] + 1)  # weight = in_deg + 1
            push!(targets, v)
        end
    end

    for new_v in UInt32(init + 1):UInt32(n)
        chosen = Set{Int}()
        attempts = 0
        while length(chosen) < m && attempts < m * 20
            t = targets[rand(rng, 1:length(targets))]
            if t != Int(new_v) && !(t in chosen)
                push!(chosen, t)
            end
            attempts += 1
        end
        for t in chosen
            add_edge!(g, new_v, UInt32(t))
            in_deg[t] += 1
            push!(targets, t)  # update weights
        end
        # New vertex gets weight 1 (in_deg=0, +1 baseline)
        push!(targets, Int(new_v))
    end

    return g
end

"""
    random_sbm_digraph(sizes, p_matrix; seed=42)

Stochastic Block Model directed graph. `sizes` is a vector of block sizes,
`p_matrix` is a k×k matrix where `p_matrix[a,b]` is the probability of an
edge from a vertex in block `a` to a vertex in block `b`.
"""
function random_sbm_digraph(sizes::Vector{Int}, p_matrix::Matrix{Float64}; seed::Int=42)
    rng = MersenneTwister(seed)
    k = length(sizes)
    @assert size(p_matrix) == (k, k)
    n = sum(sizes)
    g = SimpleDiGraph{UInt32}(UInt32(n))

    # Assign blocks
    block_start = cumsum([0; sizes[1:end-1]]) .+ 1
    block_end = cumsum(sizes)

    for a in 1:k
        for b in 1:k
            p = p_matrix[a, b]
            p <= 0.0 && continue
            for i in block_start[a]:block_end[a]
                for j in block_start[b]:block_end[b]
                    i == j && continue
                    if rand(rng) < p
                        add_edge!(g, UInt32(i), UInt32(j))
                    end
                end
            end
        end
    end

    return g
end

"""
    random_modular_hub_digraph(n, k; p_intra=0.1, p_inter=0.01, hub_fraction=0.05, seed=42)

Modular directed graph with clear community structure connected by hub vertices.

Creates `k` communities of roughly equal size. Within each community, edges exist
with probability `p_intra`. A fraction `hub_fraction` of vertices in each community
are designated as hubs; hubs connect to vertices in other communities with probability
`p_inter`. Non-hub vertices have no inter-community edges.

This produces graphs with high modularity and identifiable bridge vertices,
useful for testing community-aware compression (e.g. Leiden+LLP ordering).
"""
function random_modular_hub_digraph(n::Int, k::Int;
                                    p_intra::Float64=0.1,
                                    p_inter::Float64=0.01,
                                    hub_fraction::Float64=0.05,
                                    seed::Int=42)
    @assert k >= 2 && n >= k
    @assert 0.0 < hub_fraction <= 1.0
    rng = MersenneTwister(seed)
    g = SimpleDiGraph{UInt32}(UInt32(n))

    # Partition vertices into k communities
    base_size = n ÷ k
    sizes = fill(base_size, k)
    sizes[end] += n - sum(sizes)

    comm_start = cumsum([0; sizes[1:end-1]]) .+ 1
    comm_end = cumsum(sizes)

    # Select hub vertices per community
    hub_set = Set{Int}()
    for c in 1:k
        cs, ce = comm_start[c], comm_end[c]
        comm_size = ce - cs + 1
        n_hubs = max(1, round(Int, hub_fraction * comm_size))
        perm = randperm(rng, comm_size)
        for h in 1:n_hubs
            push!(hub_set, cs + perm[h] - 1)
        end
    end

    # Intra-community edges
    for c in 1:k
        cs, ce = comm_start[c], comm_end[c]
        for i in cs:ce
            for j in cs:ce
                i == j && continue
                if rand(rng) < p_intra
                    add_edge!(g, UInt32(i), UInt32(j))
                end
            end
        end
    end

    # Inter-community edges (hubs only)
    for c1 in 1:k
        for c2 in 1:k
            c1 == c2 && continue
            cs1, ce1 = comm_start[c1], comm_end[c1]
            cs2, ce2 = comm_start[c2], comm_end[c2]
            for i in cs1:ce1
                i in hub_set || continue
                for j in cs2:ce2
                    if rand(rng) < p_inter
                        add_edge!(g, UInt32(i), UInt32(j))
                    end
                end
            end
        end
    end

    return g
end

"""
    random_lfr_digraph(n; tau1=2.5, tau2=1.5, mu=0.1, avg_degree=10,
                        max_degree=0, min_community=0, max_community=0, seed=42)

Lancichinetti-Fortunato-Radicchi (LFR) benchmark directed graph with planted
community structure. Produces realistic modular graphs with:
- Power-law degree distribution (exponent `tau1`)
- Power-law community size distribution (exponent `tau2`)
- Tunable mixing parameter `mu` (fraction of inter-community edges per node)

Parameters:
- `n`: number of vertices
- `tau1`: degree distribution power-law exponent (>1, typically 2-3)
- `tau2`: community size distribution power-law exponent (>1, typically 1-2)
- `mu`: mixing parameter in [0,1] — fraction of each node's edges going outside its community
- `avg_degree`: target average out-degree
- `max_degree`: maximum degree (default: n÷10)
- `min_community`: minimum community size (default: avg_degree)
- `max_community`: maximum community size (default: n÷5)
- `seed`: random seed

Returns `(graph, communities)` where `communities` is a `Vector{Vector{UInt32}}`
of ground-truth community membership (sorted vertex IDs per community).

Reference: Lancichinetti, Fortunato, Radicchi (2008) "Benchmark graphs for
testing community detection algorithms" Phys Rev E 78:046110.
"""
function random_lfr_digraph(n::Int;
        tau1::Float64=2.5, tau2::Float64=1.5, mu::Float64=0.1,
        avg_degree::Int=10, max_degree::Int=0,
        min_community::Int=0, max_community::Int=0,
        seed::Int=42)
    @assert n >= 10
    @assert tau1 > 1.0 && tau2 > 1.0
    @assert 0.0 <= mu <= 1.0
    @assert avg_degree >= 1

    rng = MersenneTwister(seed)

    max_degree = max_degree > 0 ? max_degree : max(avg_degree * 3, n ÷ 10)
    min_community = min_community > 0 ? min_community : max(avg_degree, 3)
    max_community = max_community > 0 ? max_community : max(n ÷ 5, min_community + 1)

    # ── Step 1: Generate power-law degree sequence ──────────────────────
    degrees = _sample_powerlaw_sequence(rng, n, tau1, avg_degree, max_degree)

    # ── Step 2: Generate community sizes ────────────────────────────────
    comm_sizes = _sample_community_sizes(rng, n, tau2, min_community, max_community)
    k = length(comm_sizes)

    # ── Step 3: Assign vertices to communities ──────────────────────────
    # Sort vertices by degree descending, assign to communities that need
    # high-degree nodes (largest communities first)
    vertex_order = sortperm(degrees; rev=true)
    membership = zeros(Int, n)  # vertex → community index
    communities = [UInt32[] for _ in 1:k]

    # Fill communities round-robin by sorted degree
    comm_fill = zeros(Int, k)
    ci = 1
    for v in vertex_order
        # Find next community that isn't full
        attempts = 0
        while comm_fill[ci] >= comm_sizes[ci]
            ci = mod1(ci + 1, k)
            attempts += 1
            attempts > k && break
        end
        if comm_fill[ci] < comm_sizes[ci]
            membership[v] = ci
            push!(communities[ci], UInt32(v))
            comm_fill[ci] += 1
        else
            # All full — shouldn't happen, but assign to largest
            ci_largest = argmax(comm_sizes)
            membership[v] = ci_largest
            push!(communities[ci_largest], UInt32(v))
        end
        ci = mod1(ci + 1, k)
    end

    # Sort each community by vertex ID
    for C in communities
        sort!(C)
    end

    # Build community lookup: for each community, set of member vertex IDs
    comm_members = [Set{Int}(Int.(C)) for C in communities]

    # ── Step 4: Wire edges using configuration model ────────────────────
    g = SimpleDiGraph{UInt32}(UInt32(n))

    for v in 1:n
        deg = degrees[v]
        deg == 0 && continue
        c = membership[v]

        n_external = round(Int, mu * deg)
        n_internal = deg - n_external

        # Internal edges: random targets within same community
        my_comm = communities[c]
        if length(my_comm) > 1 && n_internal > 0
            for _ in 1:n_internal
                for attempt in 1:20
                    t = Int(my_comm[rand(rng, 1:length(my_comm))])
                    if t != v
                        add_edge!(g, UInt32(v), UInt32(t))
                        break
                    end
                end
            end
        end

        # External edges: random targets outside community
        if n_external > 0 && k > 1
            for _ in 1:n_external
                for attempt in 1:20
                    t = rand(rng, 1:n)
                    if t != v && !(t in comm_members[c])
                        add_edge!(g, UInt32(v), UInt32(t))
                        break
                    end
                end
            end
        end
    end

    return g, communities
end

"""Sample n values from a discrete power-law distribution P(x) ∝ x^(-tau),
   x ∈ [x_min, x_max], rescaled to achieve target average."""
function _sample_powerlaw_sequence(rng, n::Int, tau::Float64,
                                    avg_target::Int, x_max::Int)
    x_min = max(1, avg_target ÷ 3)
    degrees = zeros(Int, n)

    # Inverse CDF sampling from truncated power-law
    # P(x) ∝ x^(-tau) for x in [x_min, x_max]
    for i in 1:n
        u = rand(rng)
        lo = Float64(x_min)^(1.0 - tau)
        hi = Float64(x_max)^(1.0 - tau)
        x = (lo + u * (hi - lo))^(1.0 / (1.0 - tau))
        degrees[i] = clamp(round(Int, x), x_min, x_max)
    end

    # Rescale to achieve target average degree
    current_avg = sum(degrees) / n
    if current_avg > 0
        scale = avg_target / current_avg
        for i in 1:n
            degrees[i] = clamp(round(Int, degrees[i] * scale), 1, x_max)
        end
    end

    return degrees
end

"""Sample community sizes from a truncated power-law distribution that sum to n."""
function _sample_community_sizes(rng, n::Int, tau::Float64,
                                  min_size::Int, max_size::Int)
    sizes = Int[]
    remaining = n

    while remaining > 0
        if remaining <= max_size
            # Last community gets the remainder
            if remaining >= min_size
                push!(sizes, remaining)
            elseif !isempty(sizes)
                # Too small — merge with last community
                sizes[end] += remaining
            else
                push!(sizes, remaining)
            end
            remaining = 0
        else
            # Sample from truncated power-law
            u = rand(rng)
            lo = Float64(min_size)^(1.0 - tau)
            hi = Float64(max_size)^(1.0 - tau)
            x = (lo + u * (hi - lo))^(1.0 / (1.0 - tau))
            s = clamp(round(Int, x), min_size, min(max_size, remaining))
            push!(sizes, s)
            remaining -= s
        end
    end

    # Sort descending (largest communities first)
    sort!(sizes; rev=true)
    return sizes
end

"""
    random_web_digraph(n; avg_degree=8, seed=42)

Simulates web-like graph structure: sequential pages with locality-heavy
linking. Each vertex links to nearby predecessors (geometric decay) and
random long-range targets. Produces graphs with web-like characteristics:
high locality, skewed degree distribution, and some hub structure.
"""
function random_web_digraph(n::Int; avg_degree::Int=8, seed::Int=42)
    rng = MersenneTwister(seed)
    g = SimpleDiGraph{UInt32}(UInt32(n))

    for v in 2:n
        v_u32 = UInt32(v)
        n_links = max(1, avg_degree + rand(rng, -2:2))

        for _ in 1:n_links
            if rand(rng) < 0.7 && v > 1
                # Local back-link: geometric decay distance
                max_dist = min(v - 1, 100)
                dist = min(max_dist, 1 + floor(Int, abs(randn(rng)) * 5))
                target = v - dist
            else
                # Random long-range link
                target = rand(rng, 1:n)
            end
            if target != v && target >= 1 && target <= n
                add_edge!(g, v_u32, UInt32(target))
            end
        end
    end

    return g
end

"""
    random_watts_strogatz_digraph(n, k=8, β=0.05; seed=42)

Watts–Strogatz small-world *directed* graph. Wraps
`LightGraphs.watts_strogatz(n, k, β; is_directed=true, seed=seed)` and
re-exports the result as a `SimpleDiGraph{UInt32}` for consistency with
the rest of the `Adjacently.Generators` API.

`k` is the (mean) total degree on the underlying ring lattice (must be
even); each vertex is connected to its `k/2` nearest neighbors on each
side before rewiring. `β ∈ [0,1]` is the rewire probability — `β=0`
gives the pure ring lattice (high clustering, long average path),
`β=1` gives a near-random graph. Small β (≈ 0.01–0.1) gives the
small-world regime: high clustering plus short paths.
"""
function random_watts_strogatz_digraph(n::Int, k::Int=8, β::Float64=0.05;
                                       seed::Int=42)
    iseven(k) || throw(ArgumentError("Watts–Strogatz expects an even k, got $k"))
    src_g = watts_strogatz(n, k, β; is_directed=true, seed=seed)
    g = SimpleDiGraph{UInt32}(UInt32(n))
    for e in edges(src_g)
        add_edge!(g, UInt32(src(e)), UInt32(dst(e)))
    end
    return g
end

"""
    generate_training_batch(; n_range=(500, 2000), batch_size=16, seed=42)

Generate a mixed batch of directed graphs with diverse topologies for training.
Returns `Vector{SimpleDiGraph{UInt32}}`.
"""
function generate_training_batch(; n_range::Tuple{Int,Int}=(500, 2000),
                                   batch_size::Int=16, seed::Int=42)
    rng = MersenneTwister(seed)
    graphs = SimpleDiGraph{UInt32}[]
    n_lo, n_hi = n_range

    for i in 1:batch_size
        n = n_lo + rand(rng, 0:(n_hi - n_lo))
        gseed = seed + i * 1000

        type_idx = mod(i - 1, 4)
        if type_idx == 0
            # Erdős–Rényi with varying density
            p = 2.0 / n + rand(rng) * 10.0 / n
            g = random_erdos_renyi_digraph(n, p; seed=gseed)
        elseif type_idx == 1
            # Barabási–Albert with varying attachment
            m = rand(rng, 2:6)
            g = random_barabasi_albert_digraph(n, m; seed=gseed)
        elseif type_idx == 2
            # SBM with 2-4 blocks
            k = rand(rng, 2:4)
            base_size = n ÷ k
            sizes = fill(base_size, k)
            sizes[end] += n - sum(sizes)  # fix rounding
            p_intra = 0.01 + rand(rng) * 0.04
            p_inter = p_intra * (0.05 + rand(rng) * 0.15)
            p_mat = fill(p_inter, k, k)
            for j in 1:k; p_mat[j, j] = p_intra; end
            g = random_sbm_digraph(sizes, p_mat; seed=gseed)
        else
            # Web-like
            avg_deg = rand(rng, 4:12)
            g = random_web_digraph(n; avg_degree=avg_deg, seed=gseed)
        end
        push!(graphs, g)
    end

    return graphs
end

end # module Generators
