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

"""
    GraphGenerators

Random directed graph generators for training GNN-based compression parameter
prediction. All generators produce `SimpleDiGraph{UInt32}`.
"""
module GraphGenerators

using LightGraphs: nv, ne, add_edge!, add_vertices!, vertices, outneighbors
using Random: MersenneTwister, rand, randperm, shuffle!

using ..CustomLightGraphs: SimpleDiGraph

export random_erdos_renyi_digraph, random_barabasi_albert_digraph,
       random_sbm_digraph, random_web_digraph, generate_training_batch

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

end # module GraphGenerators
