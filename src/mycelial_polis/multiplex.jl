"""
    Multiplex

Container for the six multiplex layers of the Mycelial Polis (roadmap §7):

- `:S` — social trust network `G_S`
- `:C` — communication network `G_C`
- `:E` — economic / mutual-aid network `G_E`
- `:T` — technical infrastructure network `G_T`
- `:G` — governance network `G_G`
- `:O` — observation layer available to the host `G_O`

Each layer is a `SimpleDiGraph` over the same vertex set `1:N`. The host
sees only `:O` (a partial, noisy projection of the other layers); the
asymmetry between what is seen internally vs externally is the central
design variable.
"""
struct Multiplex
    layers::Dict{Symbol,SimpleDiGraph}
end

Multiplex() = Multiplex(Dict{Symbol,SimpleDiGraph}())

layer(m::Multiplex, k::Symbol) = m.layers[k]
nv_mp(m::Multiplex) = nv(first(values(m.layers)))

"""
    World

The full simulation state. Holds agents, the multiplex, the host adversary,
infrastructure replicas, parameter blocks, and an RNG.

Infrastructure functions `K` (§12) and per-function replica sets `R_k` are
stored here so `step!` can update them as host attacks remove or compromise
stewards.
"""
mutable struct World
    agents::Vector{Agent}
    multiplex::Multiplex
    infra::Dict{Symbol,Set{Int}}                # function-name → replica agent ids
    infra_min::Dict{Symbol,Int}                 # function-name → m_k threshold
    rng::AbstractRNG
    params::NamedTuple                          # dynamics hyperparameters
    t::Int                                      # current step
    recent_backfire_count::Vector{Int}          # rolling Δt_narrative window (item 6)
    infra_latency::Dict{Symbol,Dict{Int,Int}}   # function → (agent_id → remaining warm-up) (item 7)
    dcs::Any                                    # ::Union{Nothing, DCSState}; left untyped
                                                # so multiplex.jl doesn't need to import dcs.jl
    # H-13 — scratch counters for scaling-coalition modes. Set when the
    # coalition is seeded; read each step to compute pool growth deltas.
    # Default 0 (no-op for fixed coalition or non-coalition runs).
    _coalition_seed_committed::Int
    _coalition_prev_committed::Int
    # E3 — treasury layer. Per-cell treasury balances (Float32). Empty
    # when E3 is off; populated by `_treasury_step!` once
    # `params.treasury_enabled == true`.
    treasury::Dict{Int,Float32}
    # N-1 — counters for the faction_diversity_floor mechanism audit.
    # `faction_change_attempts` counts every call to
    # `_faction_change_rejected`; `faction_change_blocks` counts those
    # that returned `true`. Both reset by `reset_faction_counters!`.
    faction_change_attempts::Int
    faction_change_blocks::Int

    # Inner constructor with backward-compatible defaults for both
    # narrative buffer and latency map — the 7-arg positional form is
    # still in use from test helpers and `build_world`.
    World(agents, multiplex, infra, infra_min, rng, params, t,
          recent_backfire_count::Vector{Int}=Int[],
          infra_latency::Dict{Symbol,Dict{Int,Int}}=Dict{Symbol,Dict{Int,Int}}(),
          dcs::Any = nothing) =
        new(agents, multiplex, infra, infra_min, rng, params, t,
            recent_backfire_count, infra_latency, dcs, 0, 0,
            Dict{Int,Float32}(), 0, 0)
end

# Convenience constructors -----------------------------------------------------

"""
    build_world(; topology, n, seed=42, params=default_params())

Construct a `World` with `n` agents on the requested topology. Skeleton:
uses `:modular_cells` (the only one implemented at item 2) via
`Adjacently.Generators.random_modular_hub_digraph`. Items 3 will add
`:federated_hubs` and `:p2p_mesh` and richer per-layer construction.

The same trust graph `G_S` is currently reused as a starting point for
`G_C` and `G_E` (item 3 will derive each layer with edge-thinning /
projection rules). `G_O` is a noisy projection — by default the host sees
a random `p_leak`-fraction of `G_S` edges.
"""
function build_world(; topology::Symbol=:modular_cells, n::Int=200,
                       seed::Int=42, params::NamedTuple=default_params(),
                       p_leak::Float64=0.10)
    rng = MersenneTwister(seed)
    # Item 3 owns the per-topology builders + layer projections.
    mp = build_topology(topology, n; seed=seed, p_leak=p_leak)

    agents = [Agent(id=i) for i in 1:n]
    _seed_initial_committed!(agents, params.initial_committed_frac, rng)

    # Infrastructure: §12 critical functions. Skeleton assigns the initial
    # stewards/contributors round-robin to each function. Item 7 refines.
    K = (:identity, :communication, :storage, :governance,
         :payments, :reputation, :archives, :discovery)
    infra = Dict{Symbol,Set{Int}}(k => Set{Int}() for k in K)
    infra_min = Dict{Symbol,Int}(k => 2 for k in K)
    committed = [a.id for a in agents if is_committed(a)]
    for (idx, aid) in enumerate(committed)
        push!(infra[K[1 + mod(idx-1, length(K))]], aid)
    end

    return World(agents, mp, infra, infra_min, rng, params, 0)
end

function _seed_initial_committed!(agents::Vector{Agent}, frac::Float64, rng::AbstractRNG)
    n = length(agents)
    n_seed = max(1, round(Int, frac * n))
    seeds = randperm(rng, n)[1:n_seed]
    for s in seeds
        a = agents[s]
        a.role       = :contributor
        a.awareness  = 1f0
        a.commitment = 0.8f0
        a.identity   = 0.8f0
        a.L_int      = 0.7f0
        a.L_ext      = 0.3f0
    end
    return agents
end

"""
    default_params()

Baseline dynamics parameters; item 4 (adoption), item 5 (host), item 6
(repression) will own the precise tuning of each block.
"""
function default_params()
    return (
        # adoption (§8) — rule selection + stage thresholds + weight coefficients
        adoption_rule = :threshold,  # one of :threshold | :complex_contagion
        stage_thresholds = (observer=0.05f0, sympathizer=0.15f0,
                            user=0.30f0, contributor=0.50f0, steward=0.75f0),
        alpha    = 1.0f0,   # trusted-neighbor weight
        beta     = 0.3f0,   # communication-exposure weight
        gamma    = 0.2f0,   # utility weight
        delta    = 0.8f0,   # fear weight (negative)
        eta      = 0.6f0,   # backfire weight (positive)
        chi      = 0.3f0,   # identity weight
        # complex contagion (§25.3 Centola–Macy) — only used when
        # adoption_rule === :complex_contagion
        k_cc     = 2,       # required distinct adopted in-neighbours
        theta_cc = 1.0f0,   # required cumulative trust weight (Σ commitment)
        # repression / backfire (§9, full form in item 6)
        rho      = 0.3f0,   # fear gain per repression event
        kappa    = 0.1f0,   # protection-driven fear decay
        # Four-component backfire: b_i = clip(w_R·R_i + w_J·J_i + w_V·V_i + w_A·A_i)
        backfire_R_weight = 0.3f0,   # experienced/observed repression (flow)
        backfire_J_weight = 0.3f0,   # cumulative perceived injustice (stock)
        backfire_V_weight = 0.2f0,   # global visibility × broadcast amplifier
        backfire_A_weight = 0.2f0,   # narrative availability (rolling)
        # visibility on/off knob for Experiment 4. Set to 0.0 to model
        # "complete media blackout": V_i ≡ 0, only fear matters under repression.
        visibility_scaling = 1.0f0,
        narrative_window   = 5,      # Δt_narrative — rolling buffer length
        narrative_scale    = 20.0f0, # n_backfired per step that saturates A_i to 1
        # defection (§9 → :defector role transition, item 6)
        defect_fear_threshold       = 0.7f0,  # fear ≥ this AND
        defect_protection_threshold = 0.3f0,  # protection < this → committed→:defector
        # legibility (§11)
        leak_rate = 0.02f0, # per-step L_ext drift up for unprotected agents
        hide_rate = 0.05f0, # per-step L_ext drift down for active agents
        # infrastructure (§12, item 7) — replica warm-up + Φ weights
        replica_latency  = 1,        # steps before a new steward counts toward Ψ_T
        kcore_scale      = 10,       # k_core_depth / kcore_scale → [0,1] for Φ
        phi_weights = (core=0.25f0, giant=0.25f0, kcore=0.10f0,
                       psi_T=0.20f0, lambda_sat=0.20f0),
        # principles (§19, item 9) — ten monotone hooks in [0,1]
        principles = default_principles(),
        # AI capability asymmetry (§17, item 11) — both in [0, 1]
        # A_L: latent-society AI capability (boosts adoption coordination,
        #      compresses narrative_scale so backfire saturates faster).
        # A_H: host AI capability (boosts attack budget effectiveness,
        #      boosts G_O leakage rate).
        A_L = 0.0f0,
        A_H = 0.0f0,
        # seeding
        initial_committed_frac = 0.05,
        # open-population recruitment: each step, every :removed agent is
        # revived as a fresh :outsider with this per-step probability. Models
        # the social position of an arrested member being filled by a new
        # individual. At 0.0 the population is closed (default; reproduces
        # all pre-recruitment behaviour). Mass-balance interpretation:
        # steady-state attrition fraction = ρ_H / (recruitment_rate + ρ_H).
        recruitment_rate = 0.0,
        # Legacy boolean: when true, sets `recruitment_mode = :sponsored`.
        # Prefer `recruitment_mode` directly.
        recruitment_sponsored = false,
        # `:ambient | :sponsored | :bridge | :cellular`. See dynamics.jl
        # for the per-mode admission filter. The recruitment ecology
        # experiment in paper2 sweeps across these.
        recruitment_mode = :ambient,
        # Probability a `:cellular` recruitment crosses a cell boundary
        # (per :removed agent, per step). Only used when mode = :cellular.
        recruitment_cross_cell_p = 0.10,
        # E2 endogenous conflict (paper 2). Off by default
        # (disagreement_rate = 0 → no factions, no schisms).
        disagreement_rate    = 0.0,      # per-step prob of an event
        disagreement_decay   = 0.95,     # per-step decay of per-agent disagreement
        schism_threshold     = 0.6f0,    # flip faction + edge decay above this
        trust_decay_on_fight = 0.10,     # per-edge prob of removal at schism
        n_factions_cap       = 4,
        # Reversible governance: steward outdegree-fraction above this
        # triggers a potential capture-demotion (prob = principle strength).
        capture_centrality_thresh = 0.20,
        # Forkability: faction must have ≥ this many committed members to
        # trigger a clean fork.
        fork_min_cell_size = 5,
        # Per-schism probability that the schism-ing agent defects
        # (loses committed status). Scaled by (1 − forkability):
        # forkability suppresses this defection. Active only when
        # disagreement_rate > 0.
        factional_defection_prob = 0.15,
    )
end
