"""
    HostStrategy

Abstract type for host adversary strategies (roadmap §10).

A host implements `act!(host, world, t) -> Vector{Int}`. The return value
is the list of agent ids the host *acted on* this step — used by
`metrics.snapshot` to populate `n_attacked`. Strategies select targets in
their own way (random, by-degree, by-betweenness, by-legibility, localized,
infiltration-first, attrition) and then spend their per-step budget on
one or more attack types from `ATTACK_COST`.

Item 5 coverage (vs roadmap §10):

| §10 strategy                    | Concrete type           | Status     |
|---------------------------------|-------------------------|------------|
| random                          | `RandomHost`            | item 2     |
| targeted by degree centrality   | `DegreeHost`            | item 5     |
| targeted by betweenness         | `BetweennessHost`       | item 5     |
| targeted by external legibility | `LegibilityHost`        | item 5     |
| localized attack                | `LocalizedHost`         | item 5     |
| infiltration-first              | `InfiltrationFirstHost` | item 5     |
| attrition                       | `AttritionHost`         | item 5     |
| co-optation                     | —                       | item 6+    |
| legal-bureaucratic strangulation| —                       | Phase 8    |
| shock repression                | —                       | (deferred) |

Attack types (vs §10 attack types):

| §10 attack type      | Status                                            |
|----------------------|---------------------------------------------------|
| node removal         | `_node_remove!`        (item 2)                   |
| edge disruption      | `_edge_disrupt!`       (item 5)                   |
| resource attack      | `_resource_attack!`    (item 5)                   |
| legibility attack    | `_legibility_attack!`  (item 5)                   |
| infiltration         | `_infiltrate!`         (item 5)                   |
| narrative attack     | item 6+                                           |
| capture attack       | item 9                                            |
| fragmentation attack | (deferred)                                        |
"""
abstract type HostStrategy end

# Per-step cost for each attack type. Set so that arrests are the
# baseline (1), infiltration is more expensive (operationally complex),
# edge disruption is cheap-per-edge but the strategy decides how many
# edges it cuts per step.
const ATTACK_COST = Dict{Symbol,Int}(
    :node_removal      => 1,
    :legibility_attack => 1,
    :resource_attack   => 1,
    :edge_disruption   => 1,   # one disruption ≈ removing a small bundle of edges
    :infiltration      => 2,
)

# Which graph layer the host "sees" by default. The roadmap §10 model
# `H_t : (G_O, X_O) → A_t` makes the host see only the leaked observation
# layer; strategies fall back to G_S if the host is given `view = :S`
# (useful for testing topology fragility without confounds from leakage).
const DEFAULT_HOST_VIEW = :O

# --- attack-action helpers ---------------------------------------------------
# All helpers are idempotent on already-removed/infiltrated/disrupted nodes
# so a host can blindly apply them without pre-filtering.

function _node_remove!(world::World, id::Int)
    a = world.agents[id]
    a.role === :removed && return false
    a.role       = :removed
    a.awareness  = 0f0
    a.commitment = 0f0
    for (_, replicas) in world.infra
        delete!(replicas, id)
    end
    return true
end

function _legibility_attack!(world::World, id::Int; bump::Float32=0.4f0)
    a = world.agents[id]
    a.role === :removed && return false
    a.L_ext = clamp(a.L_ext + bump, 0f0, 1f0)
    return true
end

function _resource_attack!(world::World, id::Int)
    a = world.agents[id]
    a.role === :removed && return false
    a.capacity = clamp(a.capacity - 0.25f0, 0f0, 1f0)
    # Drop the target from at most one infra function per step.
    for (_, replicas) in world.infra
        if id ∈ replicas
            delete!(replicas, id)
            break
        end
    end
    return true
end

function _infiltrate!(world::World, id::Int)
    a = world.agents[id]
    a.role === :removed     && return false
    a.role === :infiltrator && return false
    a.role     = :infiltrator
    # Infiltrators report on their neighbours — bump their L_ext.
    g_s = layer(world.multiplex, :S)
    for v in outneighbors(g_s, a.id)
        nb = world.agents[Int(v)]
        nb.role === :removed && continue
        nb.L_ext = clamp(nb.L_ext + 0.15f0, 0f0, 1f0)
    end
    return true
end

function _edge_disrupt!(world::World, ids::Vector{Int}; layer_sym::Symbol=:C, per_node::Int=2)
    g = layer(world.multiplex, layer_sym)
    disrupted = 0
    @inbounds for id in ids
        edges_out = collect(outneighbors(g, eltype(g)(id)))
        isempty(edges_out) && continue
        # Drop up to `per_node` outgoing edges (idempotent if already removed).
        for v in @view edges_out[1:min(per_node, length(edges_out))]
            rem_edge!(g, eltype(g)(id), v) && (disrupted += 1)
        end
    end
    return disrupted
end

# --- 1. RandomHost (already part of item 2) ---------------------------------

"""
    RandomHost(budget=5; budget_decay=0.0)

Uniform random arrests up to `budget` units per step. `budget_decay` is
the per-step exponential decay rate applied to `budget` after each
step — `0.0` = constant budget (item 5 behaviour); higher values let
the host's repression capacity erode over time, enabling item-8-v2
classifications like `:host_accommodation`. Applied by `step!` via
`_decay_budget!`.
"""
Base.@kwdef mutable struct RandomHost <: HostStrategy
    budget::Int          = 5
    rho_H::Float64       = 0.0
    budget_decay::Float64       = 0.0
    accommodation_rate::Float64 = 0.0
    _rho_acc::Float64           = 0.0
    _decay_buffer::Float64      = 0.0
    _committed_ewma::Float32    = -1f0
end

function act!(host::RandomHost, world::World, t::Int)
    pool = [a.id for a in world.agents if is_active(a)]
    isempty(pool) && return Int[]
    n = min(div(host.budget, ATTACK_COST[:node_removal]), length(pool))
    n <= 0 && return Int[]
    order = randperm(world.rng, length(pool))
    targets = pool[order[1:n]]
    for id in targets; _node_remove!(world, id); end
    return targets
end

# --- 2. DegreeHost ----------------------------------------------------------

"""
    DegreeHost(budget=5, view=:O)

Arrest the highest-degree active agents on the chosen view layer (`:O`
for the realistic "host sees only leaks" mode; `:S` for the worst-case
"oracle" mode used in stress tests). Implements §10's
"targeted attack by degree centrality".
"""
Base.@kwdef mutable struct DegreeHost <: HostStrategy
    budget::Int           = 5
    view::Symbol          = DEFAULT_HOST_VIEW
    rho_H::Float64        = 0.0
    budget_decay::Float64       = 0.0
    accommodation_rate::Float64 = 0.0
    _rho_acc::Float64           = 0.0
    _decay_buffer::Float64      = 0.0
    _committed_ewma::Float32    = -1f0
end

function act!(host::DegreeHost, world::World, t::Int)
    g = layer(world.multiplex, host.view)
    pool = [(a.id, outdegree(g, eltype(g)(a.id)) + indegree(g, eltype(g)(a.id)))
            for a in world.agents if is_active(a)]
    isempty(pool) && return Int[]
    sort!(pool; by = x -> -x[2])
    n = min(div(host.budget, ATTACK_COST[:node_removal]), length(pool))
    targets = [p[1] for p in pool[1:n]]
    for id in targets; _node_remove!(world, id); end
    return targets
end

# --- 3. BetweennessHost -----------------------------------------------------

"""
    BetweennessHost(budget=5, view=:O, recompute_every=10)

Like `DegreeHost` but ranks targets by betweenness centrality on the
chosen view layer. Betweenness is `O(V·E)` per recompute, so the
cache is refreshed only every `recompute_every` steps; in between, the
last-known top-`k` ranking is used (with already-removed nodes skipped).
"""
Base.@kwdef mutable struct BetweennessHost <: HostStrategy
    budget::Int           = 5
    view::Symbol          = DEFAULT_HOST_VIEW
    recompute_every::Int  = 10
    rho_H::Float64        = 0.0
    budget_decay::Float64 = 0.0
    accommodation_rate::Float64 = 0.0
    _cache_t::Int         = -1
    _cache_order::Vector{Int} = Int[]
    _rho_acc::Float64        = 0.0
    _decay_buffer::Float64   = 0.0
    _committed_ewma::Float32 = -1f0
end

function act!(host::BetweennessHost, world::World, t::Int)
    if isempty(host._cache_order) || (t - host._cache_t) >= host.recompute_every
        g = layer(world.multiplex, host.view)
        # Convert to a normal undirected projection for betweenness — the
        # signal we want is "agent sits on many shortest paths", direction-
        # agnostic. LightGraphs.betweenness_centrality works on digraphs
        # directly though, so we just call it on g.
        bc = betweenness_centrality(g)
        order = sortperm(bc; rev=true)
        host._cache_order = Int.(order)
        host._cache_t = t
    end
    pool = [id for id in host._cache_order if is_active(world.agents[id])]
    n = min(div(host.budget, ATTACK_COST[:node_removal]), length(pool))
    targets = pool[1:n]
    for id in targets; _node_remove!(world, id); end
    return targets
end

# --- 4. LegibilityHost ------------------------------------------------------

"""
    LegibilityHost(budget=5, escalate_first=true, threshold=0.7)

Targets the most-legible agents (highest `L_ext`). When `escalate_first`
is true, the strategy spends one budget unit on a `:legibility_attack`
(raising the target's `L_ext` toward 1.0) before spending another on a
`:node_removal`. With it false, the strategy arrests directly.
"""
Base.@kwdef mutable struct LegibilityHost <: HostStrategy
    budget::Int          = 5
    escalate_first::Bool = true
    threshold::Float32   = 0.7f0
    rho_H::Float64       = 0.0
    budget_decay::Float64       = 0.0
    accommodation_rate::Float64 = 0.0
    _rho_acc::Float64           = 0.0
    _decay_buffer::Float64      = 0.0
    _committed_ewma::Float32    = -1f0
end

function act!(host::LegibilityHost, world::World, t::Int)
    pool = [(a.id, a.L_ext) for a in world.agents if is_active(a)]
    isempty(pool) && return Int[]
    sort!(pool; by = x -> -x[2])

    budget = host.budget
    targets = Int[]
    for (id, l) in pool
        budget <= 0 && break
        a = world.agents[id]
        if host.escalate_first && l < host.threshold && budget >= ATTACK_COST[:legibility_attack]
            _legibility_attack!(world, id)
            budget -= ATTACK_COST[:legibility_attack]
            push!(targets, id)
            continue
        end
        if budget >= ATTACK_COST[:node_removal]
            _node_remove!(world, id)
            budget -= ATTACK_COST[:node_removal]
            push!(targets, id)
        end
    end
    return targets
end

# --- 5. LocalizedHost -------------------------------------------------------

"""
    LocalizedHost(budget=5, region_size=8, view=:S, refresh_every=20)

Picks a connected region of `~region_size` agents on the view layer
(by BFS from a random seed) and concentrates this step's budget there.
Mixes arrests with edge disruption inside the region to model a
"crackdown on one neighbourhood / cell".
"""
Base.@kwdef mutable struct LocalizedHost <: HostStrategy
    budget::Int           = 5
    region_size::Int      = 8
    view::Symbol          = :S
    refresh_every::Int    = 20
    rho_H::Float64        = 0.0
    budget_decay::Float64 = 0.0
    accommodation_rate::Float64 = 0.0
    _region::Vector{Int}  = Int[]
    _region_t::Int        = -1
    _rho_acc::Float64        = 0.0
    _decay_buffer::Float64   = 0.0
    _committed_ewma::Float32 = -1f0
end

function _bfs_region(g::AbstractGraph, seed::Int, want::Int)
    T = eltype(g); n = nv(g)
    seen = falses(n); seen[seed] = true
    q = Int[seed]; head = 1
    out = Int[seed]
    while head <= length(q) && length(out) < want
        v = q[head]; head += 1
        for u in outneighbors(g, T(v))
            !seen[Int(u)] && (seen[Int(u)] = true; push!(q, Int(u)); push!(out, Int(u)))
            length(out) >= want && break
        end
        for u in inneighbors(g, T(v))
            !seen[Int(u)] && (seen[Int(u)] = true; push!(q, Int(u)); push!(out, Int(u)))
            length(out) >= want && break
        end
    end
    return out
end

function act!(host::LocalizedHost, world::World, t::Int)
    g = layer(world.multiplex, host.view)
    if isempty(host._region) || (t - host._region_t) >= host.refresh_every ||
       all(id -> !is_active(world.agents[id]), host._region)
        actives = [a.id for a in world.agents if is_active(a)]
        isempty(actives) && return Int[]
        seed = actives[rand(world.rng, 1:length(actives))]
        host._region = _bfs_region(g, seed, host.region_size)
        host._region_t = t
    end
    targets = filter(id -> is_active(world.agents[id]), host._region)
    isempty(targets) && return Int[]

    budget = host.budget
    arrested = Int[]
    # Split budget: half arrests, half edge disruption on the region.
    arrests_budget = div(budget, 2)
    n_arr = min(arrests_budget, length(targets))
    for id in targets[1:n_arr]
        _node_remove!(world, id)
        push!(arrested, id)
    end
    budget -= n_arr * ATTACK_COST[:node_removal]
    if budget >= ATTACK_COST[:edge_disruption]
        survivors = filter(id -> is_active(world.agents[id]), host._region)
        _edge_disrupt!(world, survivors; layer_sym=:C, per_node=2)
    end
    return arrested
end

# --- 6. InfiltrationFirstHost ----------------------------------------------

"""
    InfiltrationFirstHost(budget=4, plant_phase=5, escalate_phase=10)

For the first `plant_phase` steps, spends the whole budget infiltrating
random committed agents (sets their role to `:infiltrator` and raises
their `:S`-neighbours' `L_ext`). Then for the next `escalate_phase`
steps it arrests the *most-legible* agents (i.e. the ones the
infiltrators outed). After that, repeats the cycle.
"""
Base.@kwdef mutable struct InfiltrationFirstHost <: HostStrategy
    budget::Int          = 4
    plant_phase::Int     = 5
    escalate_phase::Int  = 10
    rho_H::Float64       = 0.0
    budget_decay::Float64       = 0.0
    accommodation_rate::Float64 = 0.0
    _rho_acc::Float64           = 0.0
    _decay_buffer::Float64      = 0.0
    _committed_ewma::Float32    = -1f0
end

function act!(host::InfiltrationFirstHost, world::World, t::Int)
    cycle = host.plant_phase + host.escalate_phase
    phase = mod(t - 1, cycle)
    if phase < host.plant_phase
        # Plant phase — infiltrate committed agents.
        pool = [a.id for a in world.agents
                if is_active(a) && is_committed(a) && a.role !== :infiltrator]
        isempty(pool) && return Int[]
        n = min(div(host.budget, ATTACK_COST[:infiltration]), length(pool))
        n <= 0 && return Int[]
        order = randperm(world.rng, length(pool))
        targets = pool[order[1:n]]
        for id in targets; _infiltrate!(world, id); end
        return targets
    else
        # Escalate phase — arrest the most legible (i.e. infiltrator-outed).
        pool = [(a.id, a.L_ext) for a in world.agents if is_active(a)]
        isempty(pool) && return Int[]
        sort!(pool; by = x -> -x[2])
        n = min(div(host.budget, ATTACK_COST[:node_removal]), length(pool))
        targets = [p[1] for p in pool[1:n]]
        for id in targets; _node_remove!(world, id); end
        return targets
    end
end

# --- 7. AttritionHost -------------------------------------------------------

"""
    AttritionHost(budget=2, fear_bias=true)

Slow, continuous pressure. Per-step budget is small. When `fear_bias`
is true, prioritises arresting already-fearful agents — the model is
that attrition compounds, picking off the demoralised first. Also
applies one resource attack per step (drains capacity, removes one
infra-replica slot).
"""
Base.@kwdef mutable struct AttritionHost <: HostStrategy
    budget::Int          = 2
    fear_bias::Bool      = true
    rho_H::Float64       = 0.0
    budget_decay::Float64       = 0.0
    accommodation_rate::Float64 = 0.0
    _rho_acc::Float64           = 0.0
    _decay_buffer::Float64      = 0.0
    _committed_ewma::Float32    = -1f0
end

# Resolve per-capita pressure into an integer attack budget for this
# step. When `host.rho_H > 0`, the host runs in "pressure mode": the
# desired budget per step is `rho_H × N`, accumulated as a fractional
# debt in `host._rho_acc` so non-integer per-capita rates are
# represented faithfully without per-step rounding noise. After
# resolution `host.budget` holds the integer count of attacks for the
# current step (overwritten next step). When `rho_H == 0` this is a
# no-op and the legacy fixed-`budget::Int` semantics apply.
function _resolve_rho_pressure!(host::HostStrategy, world::World)
    host.rho_H > 0 || return host
    N = length(world.agents)         # standing population (includes :removed)
    host._rho_acc += host.rho_H * N
    n = floor(Int, host._rho_acc)
    host._rho_acc -= n
    host.budget = max(n, 0)
    return host
end

# Generic decay helper: called by `step!` after `act!`. In legacy mode
# (`rho_H == 0`) decays the integer budget via a fractional accumulator
# (so a budget of 1 with decay 0.05 still erodes — loses 1 unit every
# ~20 steps). In pressure mode (`rho_H > 0`) decays `rho_H` itself,
# preserving the per-capita interpretation across `N`.
function _decay_budget!(host::HostStrategy)
    host.budget_decay > 0 || return host
    if host.rho_H > 0
        host.rho_H *= (1 - host.budget_decay)
    elseif host.budget > 0
        host._decay_buffer += host.budget_decay * host.budget
        if host._decay_buffer >= 1.0
            n_lost = floor(Int, host._decay_buffer)
            host.budget = max(0, host.budget - n_lost)
            host._decay_buffer -= n_lost
        end
    end
    return host
end

# EWMA smoothing constant for the polis-size signal used by accommodation.
const ACCOM_EWMA_ALPHA = 0.3f0
# Slack: polis must be at most this fraction below its EWMA before the
# host stops accommodating. 0.95 = "5% shrinkage counts as success".
const ACCOM_TOLERANCE  = 0.95f0

"""
    _accommodate!(host, n_committed)

Endogenous budget decay (the "AccommodatingHost" mechanism): when the
host observes that the polis's committed core isn't actually shrinking,
it lowers its own repression budget. Modelled as: maintain an EWMA of
`n_committed`; if the current count is ≥ `ACCOM_TOLERANCE × ewma`, the
host "sees its attacks failing" and pays a budget cost proportional to
`accommodation_rate × budget` (drained via the same fractional
accumulator as `_decay_budget!`).

Unlike `budget_decay` — which is a fixed schedule — accommodation is
*revealed-preference withdrawal*: the host stops spending when the
target isn't getting smaller. Set `accommodation_rate = 0` to disable;
both mechanisms can be active simultaneously.
"""
function _accommodate!(host::HostStrategy, n_committed::Int)
    rate = host.accommodation_rate
    rate <= 0 && return host
    pressure_mode = host.rho_H > 0
    if !pressure_mode && host.budget <= 0
        return host
    end
    if host._committed_ewma < 0
        # First observation — no decision yet, just initialise.
        host._committed_ewma = Float32(n_committed)
        return host
    end
    # Polis still stable or growing → host accommodates.
    if Float32(n_committed) >= host._committed_ewma * ACCOM_TOLERANCE
        if pressure_mode
            host.rho_H *= (1 - rate)
        else
            host._decay_buffer += rate * host.budget
            if host._decay_buffer >= 1.0
                n_lost = floor(Int, host._decay_buffer)
                host.budget = max(0, host.budget - n_lost)
                host._decay_buffer -= n_lost
            end
        end
    end
    # Update EWMA for next step.
    host._committed_ewma = ACCOM_EWMA_ALPHA * Float32(n_committed) +
                           (1 - ACCOM_EWMA_ALPHA) * host._committed_ewma
    return host
end

function act!(host::AttritionHost, world::World, t::Int)
    pool = [(a.id, a.fear) for a in world.agents if is_active(a)]
    isempty(pool) && return Int[]
    if host.fear_bias
        sort!(pool; by = x -> -x[2])
    else
        shuffle!(world.rng, pool)
    end
    budget = host.budget
    arrested = Int[]
    n = min(div(budget, ATTACK_COST[:node_removal]), length(pool))
    for (id, _) in pool[1:n]
        _node_remove!(world, id)
        push!(arrested, id)
    end
    budget -= n * ATTACK_COST[:node_removal]
    # One resource attack per step if there's still budget for it.
    if budget >= ATTACK_COST[:resource_attack]
        remaining = filter(p -> is_active(world.agents[p[1]]), pool)
        !isempty(remaining) && _resource_attack!(world, remaining[1][1])
    end
    return arrested
end

# --- 8. AdaptiveHost (Path 5 — co-evolving adversary) ------------------------

# Self-contained factory for a named static arm (host.jl must not depend on
# basin.jl's `_make_host`).
function _arm_host(name::Symbol, budget::Int)
    name === :random             ? RandomHost(budget=budget)             :
    name === :degree             ? DegreeHost(budget=budget)             :
    name === :betweenness        ? BetweennessHost(budget=budget)        :
    name === :legibility         ? LegibilityHost(budget=budget)         :
    name === :localized          ? LocalizedHost(budget=budget)          :
    name === :infiltration_first ? InfiltrationFirstHost(budget=budget)  :
    name === :attrition          ? AttritionHost(budget=budget)          :
    throw(ArgumentError("unknown adaptive-host arm: $name"))
end

"""
    AdaptiveHost(; budget=5,
                 arms=[:random, :degree, :infiltration_first, :attrition],
                 epsilon=0.15, lr=0.25)

Non-stationary ε-greedy **bandit** host (Path 5 — the co-evolving adversary).
Each step it selects one of its sub-strategy `arms` and delegates the attack to
it; the reward is the polis-size loss (drop in committed count) over that step,
credited to the arm played, and its per-arm value estimate `Q` is updated by an
EWMA with rate `lr`. With probability `epsilon` it explores a random arm,
otherwise it plays the current best. Unlike any single static host, it learns
online which fixed strategy most damages *this* polis on *this* topology — the
adversary the "constitution as boundary-shifter" claim must survive.

`plays` logs the chosen arm index per step; `Q`/`counts` expose the learned
values. Carries the standard host budget/pressure fields so `step!`'s generic
`_resolve_rho_pressure!` / `_decay_budget!` / `_accommodate!` apply unchanged.
"""
mutable struct AdaptiveHost <: HostStrategy
    arms::Vector{HostStrategy}
    arm_names::Vector{Symbol}
    Q::Vector{Float64}
    counts::Vector{Int}
    plays::Vector{Int}
    epsilon::Float64
    lr::Float64
    last_arm::Int
    prev_committed::Int
    budget::Int
    rho_H::Float64
    budget_decay::Float64
    accommodation_rate::Float64
    _rho_acc::Float64
    _decay_buffer::Float64
    _committed_ewma::Float32
end

function AdaptiveHost(; budget::Int=5,
                       arms::Vector{Symbol}=[:random, :degree,
                                             :infiltration_first, :attrition],
                       epsilon::Float64=0.15, lr::Float64=0.25,
                       rho_H::Float64=0.0, budget_decay::Float64=0.0,
                       accommodation_rate::Float64=0.0)
    inst = HostStrategy[_arm_host(a, budget) for a in arms]
    return AdaptiveHost(inst, collect(arms), zeros(Float64, length(arms)),
                        zeros(Int, length(arms)), Int[], epsilon, lr, 0, -1,
                        budget, rho_H, budget_decay, accommodation_rate,
                        0.0, 0.0, -1f0)
end

function act!(host::AdaptiveHost, world::World, t::Int)
    cur = count(is_committed, world.agents)
    # credit the previous step's arm with the polis loss it produced
    if host.last_arm > 0 && host.prev_committed >= 0
        r = Float64(host.prev_committed - cur)
        a = host.last_arm
        host.counts[a] += 1
        host.Q[a] += host.lr * (r - host.Q[a])
    end
    k = length(host.arms)
    a = rand(world.rng) < host.epsilon ? rand(world.rng, 1:k) : argmax(host.Q)
    host.last_arm = a
    host.prev_committed = cur
    push!(host.plays, a)
    arm = host.arms[a]
    arm.budget = host.budget            # spend the host's current budget via the chosen arm
    return act!(arm, world, t)
end
