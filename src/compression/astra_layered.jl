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

module ASTRALayered

using LightGraphs: AbstractGraph, nv, ne, vertices, outneighbors
using DataStructures: PriorityQueue, dequeue!
using ..Compression
using ...IO: BitWriter, BitReader, flush_bitwriter
using ...CustomTypes: UInt24, UInt40
using ...Graph: get_reverse_graph, get_forward_ball, subgraph, get_sparse_P_matrix
using SparseArrays: SparseMatrixCSC, sparse
using ...PageRank: PR
using ...Algo: pearce_iterative

export create_level_decomposition,
       write_astra_layered_graph,
       read_astra_layered_graph,
       estimate_astra_layered_cost,
       estimate_astra_layered_cost_breakdown

# =====================================================================================
# Data structures
# =====================================================================================

struct Level{T<:Unsigned}
    level_id::Int
    seed::T
    local_to_global::Vector{T}
    global_to_local::Dict{T,Int}
    intra_edges::Vector{Tuple{Int,Int}}
    parents_next_seed::Vector{Int}
end

struct Decomposition{T<:Unsigned}
    levels::Vector{Level{T}}
    reconciliation::Vector{Tuple{T,Tuple{Int,Int},Tuple{Int,Int}}}
    remaining_cross_edges::Vector{Tuple{T,T}}
    encoded_edges::Set{Tuple{T,T}}
end

# =====================================================================================
# Layering builder — full-ball mode with PageRank-guided seed selection
# =====================================================================================

function create_level_decomposition(g::AbstractGraph{T}, rg::AbstractGraph{T};
        radius::Int=2, damping::Float64=0.85, epsilon::Float64=1e-6,
        pr::Union{Vector{Float64},Nothing}=nothing, log_every::Int=1000,
        use_scc::Bool=false, min_level_size::Int=1) where {T<:Unsigned}
    n = nv(g)
    t0 = time()
    local_pr = pr === nothing ? PR(get_sparse_P_matrix(g); damping=damping, epsilon=epsilon) : pr
    @info "PR computed" seconds=round(time()-t0, digits=2)

    # Precompute PR order (highest first)
    order = collect(1:n)
    sort!(order, by=i -> local_pr[i], rev=true)
    order_ptr = 1

    assigned = Set{T}()
    levels = Level{T}[]
    remaining_cross_edges = Tuple{T,T}[]
    encoded_edges = Set{Tuple{T,T}}()
    pq = PriorityQueue{T,Float64}()
    in_pq = Set{T}()

    function pick_seed()::Union{T,Nothing}
        while !isempty(pq)
            v = dequeue!(pq)
            if !(v in assigned)
                return v
            end
        end
        while order_ptr <= n
            v = T(order[order_ptr])
            order_ptr += 1
            if !(v in assigned)
                return v
            end
        end
        return nothing
    end

    level_id = 0
    while true
        seed = pick_seed()
        seed === nothing && break

        # Explore radius-k ball
        ball = get_forward_ball(seed::T, g, radius)

        # Add ball vertices to PQ
        for v in ball
            if !(v in assigned) && !(v in in_pq)
                pq[v] = -local_pr[Int(v)]
                push!(in_pq, v)
            end
        end

        # Also add frontier vertices (1-hop beyond ball) to PQ for better seed coverage
        ball_set = Set(ball)
        for v in ball
            for w in outneighbors(g, Int(v))
                wT = T(w)
                if !(wT in assigned) && !(wT in in_pq) && !(wT in ball_set)
                    pq[wT] = -local_pr[Int(w)]
                    push!(in_pq, wT)
                end
            end
        end

        if use_scc
            # Legacy SCC-based approach
            sg, old2new, new2old = subgraph(g, ball)
            scc_ids = pearce_iterative(sg)
            seed_new = old2new[seed]
            seed_scc_id = scc_ids[seed_new]
            scc_new_vertices = T[]
            for (idx, cid) in enumerate(scc_ids)
                cid == seed_scc_id && push!(scc_new_vertices, T(idx))
            end
            level_vertices = [v for v in [new2old[v] for v in scc_new_vertices] if !(v in assigned)]
        else
            # Full-ball approach: use all unassigned vertices in the ball
            level_vertices = [v for v in ball if !(v in assigned)]
        end

        if isempty(level_vertices)
            continue
        end

        # Merge small levels into previous level
        if length(level_vertices) < min_level_size && !isempty(levels)
            prev = levels[end]
            for gid in level_vertices
                if gid in assigned; continue; end
                newlid = length(prev.local_to_global) + 1
                push!(prev.local_to_global, gid)
                prev.global_to_local[gid] = newlid
                push!(assigned, gid)
                # Outgoing edges from new vertex
                for dst in outneighbors(g, gid)
                    if haskey(prev.global_to_local, dst)
                        push!(prev.intra_edges, (newlid, prev.global_to_local[dst]))
                        push!(encoded_edges, (gid, dst))
                    else
                        push!(remaining_cross_edges, (gid, dst))
                    end
                end
                # Incoming edges from existing level members to new vertex
                for src in outneighbors(rg, gid)
                    srcT = T(src)
                    if haskey(prev.global_to_local, srcT) && !((srcT, gid) in encoded_edges)
                        push!(prev.intra_edges, (prev.global_to_local[srcT], newlid))
                        push!(encoded_edges, (srcT, gid))
                    end
                end
            end
            continue
        end

        # Sort by PR (descending) for better delta encoding
        sort!(level_vertices, by=v -> (local_pr[v], v), rev=true)
        local_to_global = T[level_vertices...]
        global_to_local = Dict{T,Int}()
        for (lid, gid) in enumerate(local_to_global)
            global_to_local[gid] = lid
        end

        # Build intra-level edges
        intra_edges = Tuple{Int,Int}[]
        for gid in local_to_global
            src_local = global_to_local[gid]
            for dst in outneighbors(g, gid)
                if haskey(global_to_local, dst)
                    push!(intra_edges, (src_local, global_to_local[dst]))
                    push!(encoded_edges, (gid, dst))
                else
                    push!(remaining_cross_edges, (gid, dst))
                end
            end
        end

        if log_every > 0 && (level_id % log_every == 0)
            @info "Level built" level=level_id size=length(local_to_global) intra=length(intra_edges) cross_added=count(e -> e[1] in Set(local_to_global), remaining_cross_edges)
        end

        lvl = Level{T}(level_id, seed, local_to_global, global_to_local, intra_edges, Int[])
        push!(levels, lvl)
        for v in local_to_global
            push!(assigned, v)
        end
        level_id += 1

        if length(assigned) == n
            break
        end
    end

    return Decomposition{T}(levels, Tuple{T,Tuple{Int,Int},Tuple{Int,Int}}[], remaining_cross_edges, encoded_edges)
end

# =====================================================================================
# Writer — uses full MGS compression for intra-level, delta-coded cross-level
# =====================================================================================

function write_astra_layered_graph(decomp::Decomposition{T}, filename::String, g::AbstractGraph{T};
        integer_encoding::Symbol=:fibonacci, ref_window_size::Int=7,
        log_every::Int=1000, radius::Int=2, damping::Float64=0.85,
        epsilon::Float64=1e-6) where {T<:Unsigned}
    io = open(filename, "w")
    w = BitWriter(io)
    try
        # Header: number of levels
        Compression.write_encoded_value(w, T(length(decomp.levels)), integer_encoding)

        # Forward pass levels with full MGS compression
        for (li, lvl) in enumerate(decomp.levels)
            # Seed ID
            Compression.write_encoded_value(w, T(lvl.seed), integer_encoding)
            # Level size
            nloc = length(lvl.local_to_global)
            Compression.write_encoded_value(w, T(nloc), integer_encoding)

            # Build neighbor lists in local IDs (1..nloc)
            neighbor_lists = Dict{T,Vector{T}}()
            for v in T(1):T(nloc)
                neighbor_lists[v] = T[]
            end
            for (s,d) in lvl.intra_edges
                push!(neighbor_lists[T(s)], T(d))
            end
            for v in T(1):T(nloc)
                sort!(neighbor_lists[v])
            end

            # Write compressed adjacency using full reference + interval encoding
            Compression.write_compressed_graph_data(w, neighbor_lists, :children,
                integer_encoding, true, true, true, ref_window_size)

            if log_every > 0 && (li % log_every == 0)
                @info "Level written" level=li size=nloc
            end
        end

        # Cross-level edges: grouped by source, delta-encoded
        # Build adjacency map: src_global -> sorted [dst_global...]
        cross_adj = Dict{T, Vector{T}}()
        for (u, v) in decomp.remaining_cross_edges
            lst = get!(cross_adj, u, T[])
            push!(lst, v)
        end
        # Sort and deduplicate
        for (u, lst) in cross_adj
            sort!(lst)
            unique!(lst)
        end
        sources = sort(collect(keys(cross_adj)))

        # Write number of active cross-level sources
        Compression.write_encoded_value(w, T(length(sources)), integer_encoding)

        if !isempty(sources)
            # Delta-encode source IDs
            prev_src = T(0)
            for src in sources
                delta = src - prev_src
                Compression.write_encoded_value(w, T(delta), integer_encoding)
                prev_src = src
                # Write target count + delta-encoded targets
                targets = cross_adj[src]
                Compression.write_encoded_value(w, T(length(targets)), integer_encoding)
                if !isempty(targets)
                    deltas = Compression.delta_encode_vector(targets)
                    Compression.write_hybrid_mix_encoded_list(w, deltas, integer_encoding, true, 4, false)
                end
            end
        end

        # End marker
        Compression.write_encoded_value(w, T(1), integer_encoding)

    finally
        flush_bitwriter(w; flush_last_bits=false)
        close(io)
    end
end

# =====================================================================================
# Reader (skeleton)
# =====================================================================================

function read_astra_layered_graph(filename::String)
    return nothing
end

# =====================================================================================
# Cost estimators
# =====================================================================================

function estimate_astra_layered_cost(decomp::Decomposition{T}; integer_encoding::Symbol=:fibonacci) where {T<:Unsigned}
    ints = 1 # level count
    for lvl in decomp.levels
        ints += 2 # seed + size
        nloc = length(lvl.local_to_global)
        for s in 1:nloc
            local_deg = count(e -> e[1] == s, lvl.intra_edges)
            ints += max(1, local_deg + 1)
        end
    end
    # Cross-level edges
    cross_adj = Dict{T, Int}()
    for (u, v) in decomp.remaining_cross_edges
        cross_adj[u] = get(cross_adj, u, 0) + 1
    end
    ints += 1 # source count
    ints += length(cross_adj) * 2 # delta-src + target-count per source
    ints += length(decomp.remaining_cross_edges) # target values
    ints += 1 # end marker
    return ints
end

function estimate_astra_layered_cost_breakdown(decomp::Decomposition{T}; integer_encoding::Symbol=:fibonacci) where {T<:Unsigned}
    intra_count = sum(length(l.intra_edges) for l in decomp.levels)
    cross_count = length(decomp.remaining_cross_edges)
    total_edges = intra_count + cross_count
    return Dict(
        "intra_count" => intra_count,
        "cross_count" => cross_count,
        "total_edges" => total_edges,
        "num_levels" => length(decomp.levels),
    )
end

end # module ASTRALayered
