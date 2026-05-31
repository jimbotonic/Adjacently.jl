"""
    topologies

Per-topology builders for the Mycelial Polis multiplex (roadmap §7).

`build_topology(topology, n; seed, params) → Multiplex` constructs all six
layers `{G_S, G_C, G_E, G_T, G_G, G_O}` for one of three baseline
topologies:

- `:modular_cells`   — cell-of-cells with **no hubs** (uniform inter-cell
  edges). Cells of size ~n/30.
- `:federated_hubs`  — Barabási–Albert backbone of `n_hubs ≈ √n` hubs,
  each hub the centre of a small ER cell of leaves.
- `:p2p_mesh`        — Watts–Strogatz small-world digraph (`k_ring=8,
  β=0.05`).

Layer projection rules (item 3 spec):

- `G_C` (communication): edge-thinned copy of `G_S` (`keep_C` of edges).
- `G_E` (economic / mutual-aid): edge-thinned copy of `G_S` (`keep_E`).
- `G_T` (technical infra): edge-thinned copy of `G_C` (`keep_T`).
- `G_G` (governance): top-k highest-out-degree neighbours of each node on
  `G_S` (proxy for "highest-trust"; revisited in item 6/7).
- `G_O` (host observation): noisy view of `G_S ∪ G_C` with leakage
  probability `p_leak`.
"""

"""
    build_topology(topology, n; seed=42, p_leak=0.10,
                   keep_C=0.7, keep_E=0.5, keep_T=0.5, governance_k=3)

Build the full multiplex for the given topology. Returns a `Multiplex`.
All randomness goes through `MersenneTwister(seed)` so two calls with the
same `(topology, n, seed, …)` produce bit-identical multiplexes.
"""
function build_topology(topology::Symbol, n::Int;
                        seed::Int=42,
                        p_leak::Float64=0.10,
                        keep_C::Float64=0.7,
                        keep_E::Float64=0.5,
                        keep_T::Float64=0.5,
                        governance_k::Int=3)
    rng = MersenneTwister(seed)
    g_s = _build_GS(topology, n, seed)
    g_c = _derive_GC(g_s, rng, keep_C)
    g_e = _derive_GE(g_s, rng, keep_E)
    g_t = _derive_GT(g_c, rng, keep_T)
    g_g = _derive_GG(g_s, governance_k)
    g_o = _derive_GO(g_s, g_c, rng, p_leak)
    return Multiplex(Dict(:S => g_s, :C => g_c, :E => g_e,
                          :T => g_t, :G => g_g, :O => g_o))
end

# --- per-topology G_S builders -----------------------------------------------

function _build_GS(topology::Symbol, n::Int, seed::Int)
    if topology === :modular_cells
        k = max(2, ceil(Int, n / 30))
        # p_inter is the per-edge-pair probability — but #cell-pairs grows as
        # k² while #cells (and thus total intra) grows as k. Holding p_inter
        # at the §28 spec value 0.01 makes inter dominate intra for k ≳ 10
        # and modularity collapses. Scaling by 1/(k-1) keeps the *expected
        # number of inter-edges per cell* roughly constant across n.
        p_inter = 0.01 / max(k - 1, 1)
        return _GS_modular_cells(n, k; p_intra=0.20, p_inter=p_inter, seed=seed)
    elseif topology === :federated_hubs
        n_hubs = max(2, ceil(Int, sqrt(n)))
        return _GS_federated_hubs(n, n_hubs; hub_m=2, leaf_intra_p=0.05, seed=seed)
    elseif topology === :p2p_mesh
        return Adjacently.Generators.random_watts_strogatz_digraph(
            n, 8, 0.05; seed=seed)
    else
        throw(ArgumentError("unknown topology: $topology " *
                            "(expected :modular_cells, :federated_hubs, :p2p_mesh)"))
    end
end

# :modular_cells — no hubs, uniform inter-cell edges. The roadmap (§19.1
# Non-Domination, §21 Centralization failure mode) treats hubs as an
# anti-pattern for this topology; we therefore can't use the existing
# `random_modular_hub_digraph` which asserts hub_fraction > 0.
function _GS_modular_cells(n::Int, k::Int;
                           p_intra::Float64, p_inter::Float64, seed::Int)
    rng = MersenneTwister(seed)
    g = SimpleDiGraph{UInt32}(UInt32(n))
    sizes = fill(div(n, k), k)
    sizes[end] += n - sum(sizes)
    starts = cumsum([0; sizes[1:end-1]]) .+ 1
    ends   = cumsum(sizes)
    # intra-cell ER
    @inbounds for c in 1:k
        cs, ce = starts[c], ends[c]
        for i in cs:ce, j in cs:ce
            i == j && continue
            if rand(rng) < p_intra
                add_edge!(g, UInt32(i), UInt32(j))
            end
        end
    end
    # inter-cell ER (any node, no hub restriction)
    @inbounds for c1 in 1:k, c2 in 1:k
        c1 == c2 && continue
        cs1, ce1 = starts[c1], ends[c1]
        cs2, ce2 = starts[c2], ends[c2]
        for i in cs1:ce1, j in cs2:ce2
            if rand(rng) < p_inter
                add_edge!(g, UInt32(i), UInt32(j))
            end
        end
    end
    return g
end

# :federated_hubs — Barabási–Albert backbone of `n_hubs` hubs; each hub
# gets a small ER cell of leaves. Vertex IDs 1..n_hubs are hubs; the
# remainder are leaves assigned round-robin to hubs.
function _GS_federated_hubs(n::Int, n_hubs::Int;
                            hub_m::Int, leaf_intra_p::Float64, seed::Int)
    rng = MersenneTwister(seed)
    g = SimpleDiGraph{UInt32}(UInt32(n))

    # BA backbone over the hub vertices 1..n_hubs.
    if n_hubs >= 2
        ba = Adjacently.Generators.random_barabasi_albert_digraph(n_hubs, hub_m; seed=seed)
        for e in edges(ba)
            add_edge!(g, UInt32(src(e)), UInt32(dst(e)))
        end
    end

    # Assign each leaf to a hub (round-robin), connect leaf↔hub both directions,
    # and add intra-cell ER among co-assigned leaves.
    cells = [Int[] for _ in 1:n_hubs]
    @inbounds for v in (n_hubs+1):n
        h = 1 + mod(v - n_hubs - 1, n_hubs)
        push!(cells[h], v)
        add_edge!(g, UInt32(h), UInt32(v))
        add_edge!(g, UInt32(v), UInt32(h))
    end
    @inbounds for h in 1:n_hubs
        leaves = cells[h]
        for i in leaves, j in leaves
            i == j && continue
            if rand(rng) < leaf_intra_p
                add_edge!(g, UInt32(i), UInt32(j))
            end
        end
    end
    return g
end

# --- layer derivations -------------------------------------------------------

# Drop each edge of `g` with probability `1-keep`. Same vertex set, same eltype.
function _edge_thin(g::AbstractGraph, keep::Float64, rng::AbstractRNG)
    T = eltype(g)
    out = SimpleDiGraph{T}(T(nv(g)))
    @inbounds for e in edges(g)
        rand(rng) < keep && add_edge!(out, src(e), dst(e))
    end
    return out
end

_derive_GC(g_s, rng, keep) = _edge_thin(g_s, keep, rng)
_derive_GE(g_s, rng, keep) = _edge_thin(g_s, keep, rng)
_derive_GT(g_c, rng, keep) = _edge_thin(g_c, keep, rng)

# Governance: keep the `top_k` out-neighbours of each vertex by *their* total
# degree on `g_s` (proxy for trust standing). Item 6 will replace this with
# the real `q_i` trust scoring.
function _derive_GG(g_s::AbstractGraph, top_k::Int)
    T = eltype(g_s)
    out = SimpleDiGraph{T}(T(nv(g_s)))
    @inbounds for v in 1:nv(g_s)
        nbrs = collect(outneighbors(g_s, T(v)))
        isempty(nbrs) && continue
        if length(nbrs) > top_k
            scores = [outdegree(g_s, u) + indegree(g_s, u) for u in nbrs]
            order  = sortperm(scores; rev=true)
            nbrs   = nbrs[order[1:top_k]]
        end
        for u in nbrs
            add_edge!(out, T(v), u)
        end
    end
    return out
end

# Host observation: each edge of (G_S ∪ G_C) is leaked with probability p_leak.
function _derive_GO(g_s::AbstractGraph, g_c::AbstractGraph,
                    rng::AbstractRNG, p_leak::Float64)
    T = eltype(g_s)
    out = SimpleDiGraph{T}(T(nv(g_s)))
    seen = Set{Tuple{Int,Int}}()
    @inbounds for e in edges(g_s)
        push!(seen, (Int(src(e)), Int(dst(e))))
        rand(rng) < p_leak && add_edge!(out, src(e), dst(e))
    end
    @inbounds for e in edges(g_c)
        key = (Int(src(e)), Int(dst(e)))
        key ∈ seen && continue
        push!(seen, key)
        rand(rng) < p_leak && add_edge!(out, src(e), dst(e))
    end
    return out
end

# --- topology summary --------------------------------------------------------

"""
    topology_summary(mp; cell_partition=nothing, diameter_cap=1500)
        -> Vector{NamedTuple}

One row per layer with `(layer, n, m, mean_degree, modularity, clustering,
diameter)`. `modularity` is computed only if `cell_partition` is supplied
(natural choice for `:modular_cells`); otherwise `NaN`. `diameter` is
computed only for `n ≤ diameter_cap` (otherwise `NaN`); it is the largest
shortest-path length on the underlying undirected projection.
"""
function topology_summary(mp::Multiplex;
                          cell_partition::Union{Nothing,Vector{Int}}=nothing,
                          diameter_cap::Int=1500)
    rows = NamedTuple[]
    for (k, g) in mp.layers
        n_v = nv(g); n_e = ne(g)
        mean_deg = n_v == 0 ? 0.0 : Float64(2 * n_e) / Float64(n_v)
        mod_val  = cell_partition === nothing ? NaN :
                   _directed_modularity(g, cell_partition)
        clust    = _avg_clustering(g)
        diam     = n_v <= diameter_cap ? _approx_diameter(g) : NaN
        push!(rows, (layer=k, n=Int(n_v), m=Int(n_e),
                     mean_degree=mean_deg, modularity=mod_val,
                     clustering=clust, diameter=diam))
    end
    return sort!(rows, by=r -> string(r.layer))
end

# Newman directed modularity:
#   Q = (1/m) Σ_c [ L_c − (D^out_c · D^in_c) / m ]
# where L_c = number of directed edges (u→v) with both endpoints in c,
# D^out_c / D^in_c = total out-/in-degree summed over c. Sketch in
# Leicht & Newman (2008), "Community structure in directed networks".
function _directed_modularity(g::AbstractGraph, partition::Vector{Int})
    n = nv(g)
    length(partition) == n || throw(ArgumentError("partition size mismatch"))
    m = ne(g)
    m == 0 && return 0.0
    comms = unique(partition)
    L  = Dict{Int,Float64}(c => 0.0 for c in comms)
    Do = Dict{Int,Float64}(c => 0.0 for c in comms)
    Di = Dict{Int,Float64}(c => 0.0 for c in comms)
    T = eltype(g)
    @inbounds for v in 1:n
        cv = partition[v]
        Do[cv] += Float64(outdegree(g, T(v)))
        Di[cv] += Float64(indegree(g,  T(v)))
        for u in outneighbors(g, T(v))
            partition[Int(u)] == cv && (L[cv] += 1.0)
        end
    end
    q = 0.0
    @inbounds for c in comms
        q += L[c] - Do[c] * Di[c] / m
    end
    return q / m
end

# Average undirected local clustering coefficient: per vertex, ratio of
# actual edges among its (undirected) neighbours to the maximum possible.
function _avg_clustering(g::AbstractGraph)
    n = nv(g); n == 0 && return 0.0
    total = 0.0
    counted = 0
    T = eltype(g)
    @inbounds for v in 1:n
        nbrs = Set(Int.(outneighbors(g, T(v))))
        union!(nbrs, Int.(inneighbors(g, T(v))))
        delete!(nbrs, v)
        nb = collect(nbrs)
        deg = length(nb)
        deg < 2 && continue
        nbset = Set(nb)
        triangles = 0
        for u in nb
            for w in outneighbors(g, T(u))
                w == v && continue
                Int(w) ∈ nbset && (triangles += 1)
            end
        end
        # `triangles` counts each (u,w) twice (out from u to w and out from w to u
        # if both exist); normalise by deg*(deg-1).
        total += triangles / (deg * (deg - 1))
        counted += 1
    end
    return counted == 0 ? 0.0 : total / counted
end

# Approximate diameter via two-pivot BFS on the undirected projection.
# Exact diameter on n=1000+ is too slow for a per-run summary; this returns
# the longer of two double-sweep estimates and is within 0–1 of exact on
# most graphs.
function _approx_diameter(g::AbstractGraph)
    n = nv(g); n == 0 && return 0
    T = eltype(g)
    function bfs_far(s)
        dist = fill(-1, n); dist[s] = 0
        q = Int[s]; head = 1
        far = s
        while head <= length(q)
            v = q[head]; head += 1
            for u in outneighbors(g, T(v))
                dist[Int(u)] == -1 && (dist[Int(u)] = dist[v] + 1; push!(q, Int(u));
                                      dist[Int(u)] > dist[far] && (far = Int(u)))
            end
            for u in inneighbors(g, T(v))
                dist[Int(u)] == -1 && (dist[Int(u)] = dist[v] + 1; push!(q, Int(u));
                                      dist[Int(u)] > dist[far] && (far = Int(u)))
            end
        end
        return far, dist[far]
    end
    pivot1, _ = bfs_far(1)
    pivot2, d2 = bfs_far(pivot1)
    pivot3, d3 = bfs_far(pivot2)
    return max(d2, d3)
end

"""
    natural_partition(topology, n) -> Vector{Int}

Returns the "natural" cell partition for a topology (vertex → cell id), or
`nothing` if the topology has no obvious partition. Used by
`topology_summary` to compute modularity meaningfully.
"""
function natural_partition(topology::Symbol, n::Int)
    if topology === :modular_cells
        k = max(2, ceil(Int, n / 30))
        sizes = fill(div(n, k), k); sizes[end] += n - sum(sizes)
        out = Vector{Int}(undef, n)
        idx = 1
        for (c, s) in enumerate(sizes), _ in 1:s
            out[idx] = c; idx += 1
        end
        return out
    elseif topology === :federated_hubs
        n_hubs = max(2, ceil(Int, sqrt(n)))
        out = Vector{Int}(undef, n)
        @inbounds for v in 1:n
            out[v] = v <= n_hubs ? v : 1 + mod(v - n_hubs - 1, n_hubs)
        end
        return out
    else
        return nothing
    end
end
