"""
    metrics

Viability metrics for the Mycelial Polis (roadmap §14). Item 7 replaces
the item-2 skeleton with the full §14 form:

| Metric        | Range    | Source       | Φ-weight (default) |
|---------------|----------|--------------|--------------------|
| `core_frac`   | [0,1]    | committed / N | 0.25              |
| `giant_frac`  | [0,1]    | weakly-connected component on G_S, active vertices only / N | 0.25 |
| `kcore_norm`  | [0,1]    | k_core_depth / `params.kcore_scale` | 0.10  |
| `Ψ_T`         | [0,1]    | §12 weighted infra indicator           | 0.20  |
| `Λ_sat`       | [0,1)    | Λ / (1+Λ), Λ = mean(L_int) / mean(L_ext) | 0.20 |

`params.phi_weights` is an `NamedTuple` keyed by `(:core, :giant, :kcore,
:psi_T, :lambda_sat)`. Override at construction to retune.

All individual metric functions are exported so they can be inspected and
tested independently.
"""

# --- counts ------------------------------------------------------------------

active_count(world::World)    = count(is_active,    world.agents)
committed_count(world::World) = count(is_committed, world.agents)

mean_field(world::World, f::Function) =
    isempty(world.agents) ? 0f0 :
    Float32(sum(f, world.agents) / length(world.agents))

# --- core-fraction-style metrics --------------------------------------------

"""
    committed_frac(world) -> Float32 in [0,1]

Fraction of agents with role ≥ `:user`. Numerator excludes `:removed`,
`:defector`, `:infiltrator`, `:outsider`.
"""
committed_frac(world::World) = length(world.agents) == 0 ? 0f0 :
    Float32(committed_count(world) / length(world.agents))

# --- structural metrics on G_S | active -------------------------------------

# Boolean mask: which agent ids are currently "active" in the polis-internal
# sense (not removed, not outsider, not infiltrator, not defector). This is
# stricter than `is_active` because infiltrators and defectors are technically
# in the network but not on the polis's side.
function _polis_mask(world::World)
    n = nv(layer(world.multiplex, :S))
    mask = falses(Int(n))
    for a in world.agents
        if a.role !== :removed && a.role !== :outsider &&
           a.role !== :infiltrator && a.role !== :defector && a.id <= Int(n)
            mask[a.id] = true
        end
    end
    return mask
end

"""
    giant_active_frac(world) -> Float32 in [0,1]

Size of the largest weakly-connected component on `G_S` restricted to
polis-active agents, divided by the total number of agents.
"""
function giant_active_frac(world::World)
    g = layer(world.multiplex, :S)
    n = Int(nv(g))
    n == 0 && return 0f0
    mask = _polis_mask(world)
    sum(mask) == 0 && return 0f0
    visited = falses(n)
    best = 0; q = Int[]; T = eltype(g)
    @inbounds for v in 1:n
        (mask[v] && !visited[v]) || continue
        empty!(q); push!(q, v); visited[v] = true
        size = 0
        while !isempty(q)
            u = pop!(q); size += 1
            for w in outneighbors(g, T(u))
                iw = Int(w); iw <= n && mask[iw] && !visited[iw] &&
                    (visited[iw] = true; push!(q, iw))
            end
            for w in inneighbors(g, T(u))
                iw = Int(w); iw <= n && mask[iw] && !visited[iw] &&
                    (visited[iw] = true; push!(q, iw))
            end
        end
        size > best && (best = size)
    end
    return Float32(best / n)
end

"""
    kcore_depth(world) -> Int

Largest `k` for which a non-empty `k`-core exists on the *active-vertex*
projection of `G_S`. Uses LightGraphs' `core_number` on the full G_S
then zeros the entries of non-polis-active vertices; the resulting max
core number is the depth.

**Convention**: LightGraphs.core_number on a directed graph counts
total degree (out + in). A k-clique therefore has core number 2(k−1),
not k−1. This matches what we want for the polis since both directions
of trust count toward "membership in the dense core".
"""
function kcore_depth(world::World)
    g = layer(world.multiplex, :S)
    nv(g) == 0 && return 0
    mask = _polis_mask(world)
    sum(mask) == 0 && return 0
    cn = LightGraphs.core_number(g)
    best = 0
    @inbounds for i in eachindex(cn)
        if i <= length(mask) && mask[i]
            cn[i] > best && (best = cn[i])
        end
    end
    return Int(best)
end

"""
    modularity_active(world; partition=nothing) -> Float64

Newman directed modularity on the active-vertex projection of `G_S`.
If `partition === nothing`, returns `NaN` (no natural cell grouping).
Topologies built via `build_topology` should pass
`natural_partition(topology, n)` to get a meaningful number.
"""
function modularity_active(world::World; partition::Union{Nothing,Vector{Int}}=nothing)
    partition === nothing && return NaN
    g = layer(world.multiplex, :S)
    n = Int(nv(g))
    n == 0 && return 0.0
    mask = _polis_mask(world)
    sum(mask) == 0 && return 0.0
    # Restrict to active vertices in-place by skipping edges with inactive endpoints.
    m = 0
    comms = unique(partition[mask])
    L  = Dict(c => 0.0 for c in comms)
    Do = Dict(c => 0.0 for c in comms)
    Di = Dict(c => 0.0 for c in comms)
    T = eltype(g)
    @inbounds for v in 1:n
        mask[v] || continue
        cv = partition[v]
        ods = 0; ids = 0
        for u in outneighbors(g, T(v))
            iu = Int(u); iu <= n && mask[iu] || continue
            ods += 1
            partition[iu] == cv && (L[cv] += 1.0)
        end
        for u in inneighbors(g, T(v))
            iu = Int(u); iu <= n && mask[iu] || continue
            ids += 1
        end
        Do[cv] += Float64(ods)
        Di[cv] += Float64(ids)
        m += ods   # each directed edge counted once at its tail
    end
    m == 0 && return 0.0
    q = 0.0
    @inbounds for c in comms
        q += L[c] - Do[c] * Di[c] / m
    end
    return q / m
end

# --- legibility -------------------------------------------------------------

"""
    lambda(world; eps=1f-3) -> Float32

§11 legibility ratio Λ(t) = mean(L_int) / max(mean(L_ext), ε), restricted
to polis-active agents. Can exceed 1; see `lambda_sat`.
"""
function lambda(world::World; eps::Float32=1f-3)
    active = filter(a -> a.role !== :removed && a.role !== :outsider &&
                         a.role !== :infiltrator && a.role !== :defector,
                    world.agents)
    isempty(active) && return 0f0
    l_int = sum(a -> a.L_int, active) / length(active)
    l_ext = sum(a -> a.L_ext, active) / length(active)
    return Float32(l_int / max(l_ext, eps))
end

"""
    lambda_sat(world) -> Float32 in [0,1)

Squashed legibility ratio Λ / (1 + Λ) — bounded so legibility cannot
dominate `Φ`. (Spec requirement from the *Metrics normalization spec*
section of RESEARCH_TASKS.md.)
"""
lambda_sat(world::World) = (Λ = lambda(world); Float32(Λ / (1f0 + Λ)))

# --- infrastructure (§12) ---------------------------------------------------

# Replica latency lookup — returns `true` if the agent is still in the
# warm-up window for this function (and therefore should not count toward
# Ψ_T yet). Currently stored only when latency > 0; absence means ready.
@inline _is_warming(world::World, k::Symbol, id::Int) =
    haskey(world.infra_latency, k) && haskey(world.infra_latency[k], id)

# An agent is an "active replica" if it is a polis-active committed
# agent (not removed, not defector, not infiltrator) AND has finished
# its replica-latency warm-up.
function _is_effective_replica(world::World, k::Symbol, id::Int)
    a = world.agents[id]
    is_committed(a) || return false
    return !_is_warming(world, k, id)
end

"""
    psi_T(world) -> Float32 in [0,1]

§12 technical-resilience score, weighted sum over critical functions of
the indicator `|R_k_active| ≥ m_k`. Skeleton uses uniform weights
`ω_k = 1/|K|`; future versions tied to `params.phi_weights` could
distinguish "communication" from "archives".
"""
function psi_T(world::World)
    K = collect(keys(world.infra))
    isempty(K) && return 0f0
    score = 0.0
    for k in K
        active = count(id -> _is_effective_replica(world, k, id), world.infra[k])
        score += (active >= world.infra_min[k]) ? 1.0 : 0.0
    end
    return Float32(score / length(K))
end

"""
    infra_outages(world) -> Int

Number of critical functions currently below their `m_k` threshold of
effective replicas.
"""
function infra_outages(world::World)
    n_out = 0
    for k in keys(world.infra)
        active = count(id -> _is_effective_replica(world, k, id), world.infra[k])
        active < world.infra_min[k] && (n_out += 1)
    end
    return n_out
end

# --- governance ------------------------------------------------------------

"""
    gamma(world) -> Float32 in [0,1]

§13 governance score placeholder. Skeleton: fraction of agents who hold
the `:steward` role. Item 9 (principles) will broaden this with the full
`F(participation, legitimacy, speed, capture_risk, ...)` form.
"""
gamma(world::World) =
    isempty(world.agents) ? 0f0 :
    Float32(count(is_steward, world.agents) / length(world.agents))

# --- Φ composite -----------------------------------------------------------

"""
    default_phi_weights() -> NamedTuple

Default Φ-component weights. Sum to 1.0 so `Φ ∈ [0,1]` even before any
saturation. Override by passing
`params = merge(default_params(), (phi_weights = (...),))`.
"""
default_phi_weights() = (
    core        = 0.25f0,
    giant       = 0.25f0,
    kcore       = 0.10f0,
    psi_T       = 0.20f0,
    lambda_sat  = 0.20f0,
)

"""
    phi_components(world) -> NamedTuple

Per-step raw inputs to the Φ composite, each already normalised to `[0, 1]`
(`kcore` is clamped by `kcore_scale`). Returned separately so the composite is
transparent and Φ can be recomputed offline under arbitrary `phi_weights`
(Path 1 R_NH sensitivity). `phi` is the `phi_weights`-weighted sum of these.
"""
function phi_components(world::World)
    kcore_scale = Float32(get(world.params, :kcore_scale, 10))
    return (core       = committed_frac(world),
            giant      = giant_active_frac(world),
            kcore      = min(Float32(kcore_depth(world) / kcore_scale), 1f0),
            psi_T      = psi_T(world),
            lambda_sat = lambda_sat(world))
end

"""
    phi(world; partition=nothing) -> Float32 in [0,1]

§14 viability composite. Uses `world.params.phi_weights` if present,
else `default_phi_weights()`. K-core depth is normalised by
`params.kcore_scale` (default 10) so the normalised value lives in
[0,1].
"""
function phi(world::World; partition::Union{Nothing,Vector{Int}}=nothing)
    n = length(world.agents)
    n == 0 && return 0f0
    w = get(world.params, :phi_weights, default_phi_weights())
    c = phi_components(world)
    return clamp(
        w.core       * c.core +
        w.giant      * c.giant +
        w.kcore      * c.kcore +
        w.psi_T      * c.psi_T +
        w.lambda_sat * c.lambda_sat,
        0f0, 1f0)
end

# --- hierarchy score H(t) — paper 2 ---------------------------------------
#
# H(t) ∈ [0, 1] is the polis's *hierarchy / domination / minority-control*
# score, deliberately kept ORTHOGONAL to Φ so that the politically central
# failure mode "resilient domination" (high Φ + high H) is detectable.
#
# Components (mean of available signals):
#
#   steward_concentration  = top-3 stewards' G_S outdegree / total committed G_S edges
#                            ∈ [0, 1], high = a few stewards dominate the trust graph.
#   infra_concentration    = max replicas per agent / max(total replicas, 1)
#                            ∈ [0, 1], high = one agent holds all infrastructure keys.
#   faction_dominance      = max faction size / max(n_committed, 1)
#                            ∈ [0, 1], high = one faction controls the committed core.
#   steward_fraction       = n_steward / max(n_committed, 1) capped at the "few-
#                            stewards" tail: small steward population governing
#                            many is hierarchy; large rotating set is not.
#                            Maps n_steward/n_committed = 0.05 → H=1, ≥ 0.5 → H=0.
#
# Each component is reported separately for paper-2 transparency; the
# scalar H(t) is the equally-weighted mean of the components actually
# defined for the world's current state.

"""
    hierarchy_components(world) -> NamedTuple

Per-step components of the hierarchy score, in `[0, 1]` (higher = more
hierarchical). Returns each signal separately so the composite is
transparent.
"""
function hierarchy_components(world::World)
    g_S = layer(world.multiplex, :S)
    committed_ids = [a.id for a in world.agents if is_committed(a)]
    n_committed = length(committed_ids)
    n_steward = count(is_steward, world.agents)

    # 1. steward_concentration
    sc = 0f0
    if n_committed > 0 && n_steward > 0
        steward_ids = [a.id for a in world.agents if a.role === :steward]
        # Total committed-to-committed edges on G_S
        total_committed_edges = 0
        for u in committed_ids
            for v in outneighbors(g_S, eltype(g_S)(u))
                is_committed(world.agents[Int(v)]) && (total_committed_edges += 1)
            end
        end
        if total_committed_edges > 0
            steward_out = [outdegree(g_S, eltype(g_S)(s)) for s in steward_ids]
            sort!(steward_out; rev = true)
            top = sum(steward_out[1:min(3, length(steward_out))])
            sc = clamp(Float32(top / total_committed_edges), 0f0, 1f0)
        end
    end

    # 2. infra_concentration
    ic = 0f0
    # `init = 0` because a world with no infrastructure functions at all
    # (e.g. an adoption-only fixture) leaves `world.infra` empty, and an
    # empty generator has no zero for `sum` to fall back on.
    total_replicas = sum(length(s) for (_, s) in world.infra; init = 0)
    if total_replicas > 0
        per_agent = Dict{Int, Int}()
        for (_, s) in world.infra
            for id in s
                per_agent[id] = get(per_agent, id, 0) + 1
            end
        end
        max_per = isempty(per_agent) ? 0 : maximum(values(per_agent))
        ic = clamp(Float32(max_per / total_replicas), 0f0, 1f0)
    end

    # 3. faction_dominance
    fd = 0f0
    if n_committed > 0
        fcounts = Dict{Symbol, Int}()
        for a in world.agents
            is_committed(a) || continue
            a.faction === :none && continue
            fcounts[a.faction] = get(fcounts, a.faction, 0) + 1
        end
        if !isempty(fcounts)
            fd = clamp(Float32(maximum(values(fcounts)) / n_committed), 0f0, 1f0)
        end
    end

    # 4. steward_fraction — "few stewards governing many" maps to high H.
    #    n_steward/n_committed of 0.05 or below → H=1; ≥ 0.5 → H=0.
    sf = 0f0
    if n_committed > 0
        ratio = n_steward / n_committed
        sf = clamp(Float32((0.5 - ratio) / 0.45), 0f0, 1f0)
    end

    return (steward_concentration = sc,
            infra_concentration   = ic,
            faction_dominance     = fd,
            steward_fraction      = sf)
end

"""
    hierarchy_score(world) -> Float32

Scalar H(t) ∈ [0, 1]: equally-weighted mean of the components in
`hierarchy_components(world)`. Higher = more hierarchical / captured.

A "captured" polis is one where H exceeds a threshold (e.g., 0.5)
for many consecutive steps; this single-step scalar is the input to
that downstream detector.
"""
function hierarchy_score(world::World)
    c = hierarchy_components(world)
    return Float32(mean((c.steward_concentration, c.infra_concentration,
                          c.faction_dominance, c.steward_fraction)))
end

# --- informal hierarchy H_informal (paper 2 — Path 3) -----------------------
#
# The structural H above is role-based: it can only see hierarchy that is
# *declared* (steward roles, infra keys, faction labels). It is blind to the
# "tyranny of structurelessness" — an agent that dominates the realised
# interaction network while holding no title. H_informal measures domination
# from the realised trust graph G_S itself, independent of any role:
#
#   H_informal = Gini of the eigenvector-centrality distribution over the
#                polis-active induced subgraph of G_S.
#
# High Gini = influence concentrated in a few agents (informal elite);
# low Gini = flat / genuinely non-hierarchical. Because it reads the network,
# not the role labels, a captured-but-untitled elite scores high here even when
# structural H reads low.

"""
    _gini(x) -> Float32 in [0,1]

Gini coefficient of a non-negative vector (concentration; 0 = perfectly flat).
"""
function _gini(x::AbstractVector{<:Real})
    n = length(x)
    n == 0 && return 0f0
    s = sum(x)
    s <= 0 && return 0f0
    xs = sort(x)
    cum = 0.0
    @inbounds for i in 1:n
        cum += i * xs[i]
    end
    g = (2cum) / (n * s) - (n + 1) / n
    return Float32(clamp(g, 0f0, 1f0))
end

"""
    influence_centrality(world; mode=:eigenvector) -> Vector{Float32}

Centrality of each polis-active agent on the realised trust graph `G_S`,
restricted to the active-vertex induced subgraph. `:eigenvector` (default) uses
eigenvector centrality on the undirected projection (Perron–Frobenius stable);
`:pagerank` uses PageRank on the directed subgraph (robust to disconnection).
Returns `Float32[]` when fewer than 3 active agents.
"""
function influence_centrality(world::World; mode::Symbol=:eigenvector)
    g = layer(world.multiplex, :S)
    mask = _polis_mask(world)
    ids = findall(mask)
    length(ids) < 3 && return Float32[]
    sub, _ = LightGraphs.induced_subgraph(g, ids)
    if mode === :pagerank
        return Float32.(LightGraphs.pagerank(sub))
    end
    ug = LightGraphs.Graph(sub)                     # undirected projection
    LightGraphs.ne(ug) == 0 && return zeros(Float32, length(ids))
    try
        return Float32.(LightGraphs.eigenvector_centrality(ug))
    catch                                            # non-convergence fallback
        d = Float32.(LightGraphs.degree(ug)); sd = sum(d)
        return sd > 0 ? d ./ sd : zeros(Float32, length(ids))
    end
end

"""
    informal_hierarchy(world; mode=:eigenvector) -> NamedTuple

Realised-influence concentration on `G_S`, role-independent.
`H_informal` = Gini of `influence_centrality`; `top_decile_share` = fraction of
total centrality held by the top 10% of active agents; `n` = active count.
"""
function informal_hierarchy(world::World; mode::Symbol=:eigenvector)
    c = influence_centrality(world; mode=mode)
    isempty(c) && return (H_informal = 0f0, top_decile_share = 0f0, n = 0)
    xs = sort(c; rev=true)
    k = max(1, ceil(Int, 0.1 * length(xs)))
    s = sum(c)
    tds = s > 0 ? Float32(sum(xs[1:k]) / s) : 0f0
    return (H_informal = _gini(c), top_decile_share = tds, n = length(c))
end

"""
    H_informal(world) -> Float32 in [0,1]

Informal-hierarchy order parameter: Gini of realised-influence centrality on
`G_S`. See [`informal_hierarchy`](@ref).
"""
H_informal(world::World) = informal_hierarchy(world).H_informal

"""
    influence_share(world, ids; mode=:eigenvector) -> Float32 in [0,1]

Fraction of total realised-influence centrality (see [`influence_centrality`])
held by the agents in `ids`. Compared against `length(ids)/n_active`, this is a
monotonic, group-level measure of informal capture: a small group whose
influence share far exceeds its population share is an informal elite,
independent of any formal role. Robust to the non-monotonicity of the global
Gini `H_informal` w.r.t. single-tier vs nested hub structure.
"""
function influence_share(world::World, ids; mode::Symbol=:eigenvector)
    g = layer(world.multiplex, :S)
    mask = _polis_mask(world)
    active = findall(mask)
    length(active) < 3 && return 0f0
    c = influence_centrality(world; mode=mode)      # aligned to `active` order
    s = sum(c)
    s <= 0 && return 0f0
    pos = Dict(v => i for (i, v) in enumerate(active))
    tot = 0f0
    for id in ids
        p = get(pos, Int(id), 0)
        p > 0 && (tot += c[p])
    end
    return Float32(tot / s)
end

# --- null-corrected informal hierarchy (Path 3 deferred refinement) ---------
#
# Raw H_informal has a non-zero baseline set by the trust topology's intrinsic
# hub structure: a graph with power-law degrees scores high Gini even with no
# *designed* elite. H_informal_excess subtracts the mean H_informal of degree-
# preserving rewirings of the same active subgraph, isolating influence
# concentration BEYOND what the degree sequence alone forces. ~0 = influence is
# exactly as concentrated as the degrees imply (topological baseline only);
# > 0 = concentration in excess of degree (nested / rich-club hierarchy). This
# gives the clean zero point the raw Gini lacks and makes the order parameter
# comparable across topologies with different intrinsic hub structure.

function _degree_preserving_null_gini(ug, rng; n_null::Int, swaps_mult::Int)
    m = LightGraphs.ne(ug)
    m < 2 && return 0f0
    n = LightGraphs.nv(ug)
    ginis = Float32[]
    for _ in 1:n_null
        h = LightGraphs.Graph(n)
        E = Vector{Tuple{Int,Int}}(undef, m)
        i = 0
        for e in LightGraphs.edges(ug)
            i += 1
            E[i] = (LightGraphs.src(e), LightGraphs.dst(e))
            LightGraphs.add_edge!(h, E[i][1], E[i][2])
        end
        for _ in 1:(swaps_mult * m)          # double-edge swaps preserve degrees
            p = rand(rng, 1:m); q = rand(rng, 1:m)
            p == q && continue
            a, b = E[p]
            c, d = E[q]
            rand(rng, Bool) && ((c, d) = (d, c))
            length(Set((a, b, c, d))) < 4 && continue
            (LightGraphs.has_edge(h, a, d) || LightGraphs.has_edge(h, c, b)) && continue
            LightGraphs.rem_edge!(h, a, b); LightGraphs.rem_edge!(h, c, d)
            LightGraphs.add_edge!(h, a, d); LightGraphs.add_edge!(h, c, b)
            E[p] = (a, d); E[q] = (c, b)
        end
        cent = try
            Float32.(LightGraphs.eigenvector_centrality(h))
        catch
            deg = Float32.(LightGraphs.degree(h)); sd = sum(deg)
            sd > 0 ? deg ./ sd : zeros(Float32, n)
        end
        push!(ginis, _gini(cent))
    end
    return mean(ginis)
end

"""
    H_informal_excess(world; mode=:eigenvector, n_null=6, swaps_mult=8,
                      rng=MersenneTwister(0)) -> Float32

Null-corrected informal-hierarchy order parameter: `H_informal(world)` minus the
mean `H_informal` of `n_null` degree-preserving rewirings of the same active
induced subgraph of `G_S`. Isolates realised-influence concentration in excess
of the degree sequence, giving a clean ~0 zero point on a run whose only
concentration is topological (see [`informal_hierarchy`](@ref)). Deterministic
given `rng`; does not touch `world.rng`. May be slightly negative (sampling noise
around a true zero).
"""
function H_informal_excess(world::World; mode::Symbol=:eigenvector,
                           n_null::Int=6, swaps_mult::Int=8,
                           rng::AbstractRNG=MersenneTwister(0))
    g = layer(world.multiplex, :S)
    mask = _polis_mask(world)
    ids = findall(mask)
    length(ids) < 3 && return 0f0
    sub, _ = LightGraphs.induced_subgraph(g, ids)
    ug = LightGraphs.Graph(sub)
    LightGraphs.ne(ug) == 0 && return 0f0
    observed = _gini(influence_centrality(world; mode=mode))
    null = _degree_preserving_null_gini(ug, rng; n_null=n_null, swaps_mult=swaps_mult)
    return Float32(observed - null)
end

# --- per-step snapshot ------------------------------------------------------

"""
    snapshot(world; t, attacked, defected, backfired, partition=nothing)

Per-step metrics record (one row of `metrics.tsv`). Returns a NamedTuple
with the full column set from RESEARCH_TASKS.md item 7 spec, plus the
extras introduced in items 5–6 (`n_infiltrators`, `n_defectors`).

`partition` is optional but required for a meaningful `modularity_S`
value; pass `natural_partition(topology, n)` from `topologies.jl`.
"""
function snapshot(world::World; t::Int, attacked::Int=0,
                  defected::Int=0, backfired::Int=0,
                  partition::Union{Nothing,Vector{Int}}=nothing)
    Λ   = lambda(world)
    Λs  = lambda_sat(world)
    pc  = phi_components(world)        # 5 Φ inputs, each ∈ [0,1]
    hc  = hierarchy_components(world)  # 4 H inputs, each ∈ [0,1]
    ih  = informal_hierarchy(world)    # realised-influence concentration (Path 3)
    return (
        t              = t,
        n_active       = active_count(world),
        n_committed    = committed_count(world),
        n_steward      = count(is_steward, world.agents),
        n_infiltrators = count(is_infiltrator, world.agents),
        n_defectors    = count(is_defector,    world.agents),
        k_core_depth   = kcore_depth(world),
        giant_comp_S   = giant_active_frac(world),
        modularity_S   = Float32(modularity_active(world; partition=partition)),
        mean_fear      = mean_field(world, a -> a.fear),
        mean_backfire  = mean_field(world, a -> a.backfire),
        mean_L_int     = mean_field(world, a -> a.L_int),
        mean_L_ext     = mean_field(world, a -> a.L_ext),
        Lambda         = Λ,
        Lambda_sat     = Λs,
        Psi_T          = psi_T(world),
        infra_outages  = infra_outages(world),
        Gamma          = gamma(world),
        Phi            = phi(world; partition=partition),
        H              = hierarchy_score(world),
        # --- Φ components (Path 1: offline R_NH sensitivity to phi_weights) ---
        phi_core       = pc.core,
        phi_giant      = pc.giant,
        phi_kcore      = pc.kcore,
        phi_psi_T      = pc.psi_T,
        phi_lambda_sat = pc.lambda_sat,
        # --- H components (offline sensitivity to H-weighting / thresholds) ---
        h_steward_conc = hc.steward_concentration,
        h_infra_conc   = hc.infra_concentration,
        h_faction_dom  = hc.faction_dominance,
        h_steward_frac = hc.steward_fraction,
        # --- informal hierarchy (Path 3): realised-influence concentration ---
        H_informal          = ih.H_informal,
        informal_top_decile = ih.top_decile_share,
        n_attacked     = attacked,
        n_defected     = defected,
        n_backfired    = backfired,
    )
end

"""
    write_metrics_tsv(path, snaps; prefix=nothing)

Serialise a vector of `snapshot` NamedTuples to a TSV at `path`, one row per
snapshot, header derived dynamically from the snapshot keys (so the Φ/H
component columns added for Path 1 are always included — no hand-maintained
column list to fall out of sync).

`prefix`, if given, is a vector of NamedTuples of the same length as `snaps`
providing leading run-identifier columns (e.g. `(scenario=..., seed=...)`);
its keys are prepended to the header and its values to each row. All `prefix`
entries must share the same keys.
"""
function write_metrics_tsv(path::AbstractString, snaps::AbstractVector;
                           prefix::Union{Nothing,AbstractVector}=nothing)
    isempty(snaps) && (touch(path); return path)
    prefix !== nothing && length(prefix) != length(snaps) &&
        throw(ArgumentError("prefix length $(length(prefix)) != snaps length $(length(snaps))"))
    pcols = prefix === nothing ? () : keys(prefix[1])
    scols = keys(snaps[1])
    open(path, "w") do io
        println(io, join(string.((pcols..., scols...)), '\t'))
        for (i, s) in enumerate(snaps)
            pvals = prefix === nothing ? () : values(prefix[i])
            println(io, join(string.((pvals..., values(s)...)), '\t'))
        end
    end
    return path
end
