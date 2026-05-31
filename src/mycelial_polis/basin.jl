"""
    basin

Basin-of-attraction estimator for the Mycelial Polis (item 8 v1 of
`research/mycelial_polis/RESEARCH_TASKS.md`). Splits the (L_0, fear_0,
host_budget_0, p_leak_0, θ_0) initial-condition space into one
Latin-hypercube sample per "cell", runs each forward, then classifies the
trajectory into one of three §15.1 equilibria (plus an `:unclassified`
catch-all). v2 deferred until items 6 / 9 / AI variables.

The testable pieces (`lhs`, `classify`, `time_to_attractor`, `map_lhs`,
`run_sample`) live here; `scripts/mp_basin.jl` is a thin CLI wrapper
that calls `estimate_basins`.
"""

# --- classifier thresholds --------------------------------------------------

const EPS_EXT      = 0.05f0     # Φ below this counts as extinction
const TAIL_FRAC    = 0.25       # use last 25% of steps for stability check
const STABLE_STD   = 0.05f0     # tail-std below this = "stabilised"
const PARALLEL_MIN = 0.35f0     # tail-mean Φ above this = parallel coexistence

# v2 thresholds (item 8 v2)
const TRANSFORM_FRAC      = 0.50    # committed/N above this = transformative replacement
const ATTACK_DECAY_THRESH = 0.50    # attack-density drop (last vs first quarter) for accommodation/hybridization
const POLIS_GROWTH_THRESH = 0.20    # relative growth in committed = accommodation

# --- tier profiles (item 0) --------------------------------------------------

const BASIN_TIERS = Dict{Symbol,NamedTuple}(
    :smoke    => (n_samples=200,  n_world=100, t_max=40),
    :standard => (n_samples=500,  n_world=300, t_max=80),
    :nightly  => (n_samples=2000, n_world=500, t_max=200),
)

# --- Latin-hypercube sampler -------------------------------------------------

"""
    lhs(n_samples, n_dims; rng) -> Matrix{Float64}

Standard Latin-hypercube: split `[0,1]` into `n_samples` strata per
dimension; draw one uniform point from each stratum; shuffle strata
independently per dimension. Returns `n_samples × n_dims` matrix in
`[0,1]`.
"""
function lhs(n_samples::Int, n_dims::Int; rng::AbstractRNG)
    M = Matrix{Float64}(undef, n_samples, n_dims)
    @inbounds for d in 1:n_dims
        perm = randperm(rng, n_samples)
        for i in 1:n_samples
            M[i, d] = (perm[i] - 1 + rand(rng)) / n_samples
        end
    end
    return M
end

# --- sample-parameter mapping ------------------------------------------------

struct SampleParams
    L_0::Float64         # initial_committed_frac  ∈ [0.02, 0.30]
    fear_0::Float32      # injected initial fear   ∈ [0.0,  0.5]
    host_budget::Int     # RandomHost budget       ∈ {0..8}
    p_leak::Float64      # G_O leak probability    ∈ [0.0,  0.5]
    theta_scale::Float32 # stage-threshold scaling ∈ [0.5,  1.5]
end

map_lhs(u::AbstractVector{<:Real}) = SampleParams(
    0.02 + 0.28 * u[1],
    Float32(0.0  + 0.5  * u[2]),
    round(Int,   0  +   8  * u[3]),
    0.0  + 0.5  * u[4],
    Float32(0.5  + 1.0  * u[5]),
)

function build_sample_params(sp::SampleParams;
                             visibility::Float32=1.0f0,
                             principles::Principles=default_principles(),
                             A_L::Float32=0f0, A_H::Float32=0f0)
    base = default_params()
    thr  = base.stage_thresholds
    scaled = (observer    = thr.observer    * sp.theta_scale,
              sympathizer = thr.sympathizer * sp.theta_scale,
              user        = thr.user        * sp.theta_scale,
              contributor = thr.contributor * sp.theta_scale,
              steward     = thr.steward     * sp.theta_scale)
    return merge(base, (initial_committed_frac = sp.L_0,
                        stage_thresholds       = scaled,
                        visibility_scaling     = visibility,
                        principles             = principles,
                        A_L                    = A_L,
                        A_H                    = A_H))
end

# --- classifier --------------------------------------------------------------

"""
    classify(history::Vector{Float32}, n_committed_final::Int) -> Symbol

Map a per-step Φ trajectory + final committed count to one of:
`:extinction`, `:dormant_persistence`, `:parallel_coexistence`,
`:unclassified`. Discrimination rules in the module docstring.
"""
function classify(history::Vector{Float32}, n_committed_final::Int)
    T = length(history)
    T == 0 && return :unclassified
    if history[end] < EPS_EXT || n_committed_final == 0
        return :extinction
    end
    tail_start = max(1, ceil(Int, (1 - TAIL_FRAC) * T))
    tail = @view history[tail_start:T]
    tail_std  = length(tail) > 1 ? std(tail) : 0f0
    tail_mean = mean(tail)
    is_stable = tail_std < STABLE_STD
    is_stable || return :unclassified
    return tail_mean < PARALLEL_MIN ? :dormant_persistence : :parallel_coexistence
end

"""
    classify_v2(phi_history, n_committed_history, n_attacked_history,
                n_total) -> Symbol

Extended classification adding three v2 equilibria from §15.1:
`:transformative_replacement`, `:host_accommodation`,
`:institutional_hybridization`. Falls back to v1 labels when the v2
predicates don't fire. `:ai_irrelevance` is not implemented — needs
the §17 AI variables (deferred).

Discrimination rules (in order):
1. **Extinction** — `Φ(end) < EPS_EXT` OR `n_committed(end) == 0`.
2. **Transformative replacement** — `n_committed(end) ≥ TRANSFORM_FRAC × N`.
3. If the tail of `Φ` is not stable: `:unclassified`.
4. **Host accommodation** — attack-density in last quarter dropped by
   `≥ ATTACK_DECAY_THRESH` vs first quarter AND polis grew (relative
   change in mean committed ≥ POLIS_GROWTH_THRESH). The host gave up;
   polis won.
5. **Institutional hybridization** — same attack decay, but polis
   *stable* (no growth). Both sides settled into a lower-intensity
   equilibrium.
6. **Dormant persistence / parallel coexistence** — v1 categories on
   tail-mean Φ.

Requires `n_attacked_history` per step; the v1 `classify` is preserved
unchanged for callers that don't carry that channel.
"""
function classify_v2(phi_history::Vector{Float32},
                     n_committed_history::Vector{Int},
                     n_attacked_history::Vector{Int},
                     n_total::Int;
                     A_L::Float32=0f0, A_H::Float32=0f0)
    T = length(phi_history)
    T == 0 && return :unclassified
    final_committed = n_committed_history[end]
    if phi_history[end] < EPS_EXT || final_committed == 0
        return :extinction
    end
    # §17 / §15.1 Equilibrium 6 — `:ai_irrelevance`. When both sides are
    # heavily AI-mediated (A_L ≥ 0.5 AND A_H ≥ 0.5) and the polis hasn't
    # collapsed, the dynamics are dominated by AI capabilities rather
    # than human agency on either side. Coarse rule for v1; future work
    # can refine with explicit human-agency metrics.
    if A_L >= 0.5f0 && A_H >= 0.5f0
        return :ai_irrelevance
    end
    # v2: polis dominates the network → transformative replacement.
    if final_committed >= TRANSFORM_FRAC * n_total
        return :transformative_replacement
    end
    # Stability gate
    tail_start = max(1, ceil(Int, (1 - TAIL_FRAC) * T))
    head_end   = max(1, floor(Int, TAIL_FRAC * T))
    tail_phi   = @view phi_history[tail_start:T]
    tail_std   = length(tail_phi) > 1 ? std(tail_phi) : 0f0
    tail_phi_mean = mean(tail_phi)
    tail_std < STABLE_STD || return :unclassified

    # v2: attack-density decay
    head_attacks = mean(@view n_attacked_history[1:head_end])
    tail_attacks = mean(@view n_attacked_history[tail_start:T])
    attack_decay = head_attacks > 0 ? 1.0 - tail_attacks / head_attacks : 0.0

    if attack_decay >= ATTACK_DECAY_THRESH
        # Host gave up. Did polis grow or stay flat?
        head_comm = mean(@view n_committed_history[1:head_end])
        tail_comm = mean(@view n_committed_history[tail_start:T])
        polis_growth = head_comm > 0 ? (tail_comm - head_comm) / head_comm : 0.0
        if polis_growth >= POLIS_GROWTH_THRESH
            return :host_accommodation
        else
            return :institutional_hybridization
        end
    end

    # Fall through to v1 dormant / parallel.
    return tail_phi_mean < PARALLEL_MIN ? :dormant_persistence : :parallel_coexistence
end

"""
    time_to_attractor(history::Vector{Float32}) -> Int

First step `t` for which the rolling std over a `max(5, T÷10)`-step
window ending at `t` drops below `STABLE_STD`. Returns `length(history)`
if it never stabilises.
"""
function time_to_attractor(history::Vector{Float32})
    T = length(history)
    window = max(5, div(T, 10))
    @inbounds for t in window:T
        s = std(@view history[(t - window + 1):t])
        s < STABLE_STD && return t
    end
    return T
end

# --- per-sample runner -------------------------------------------------------

"""
    _make_host(strategy, budget)

Construct a `HostStrategy` instance from a symbol name. Adds the
basin-estimator coverage of item 5's seven host strategies without
having to pass live host objects through the LHS sampler.
"""
function _make_host(strategy::Symbol, budget::Int)
    if strategy === :random
        return RandomHost(budget=budget)
    elseif strategy === :degree
        return DegreeHost(budget=budget)             # default view=:O
    elseif strategy === :betweenness
        return BetweennessHost(budget=budget)
    elseif strategy === :legibility
        return LegibilityHost(budget=budget)
    elseif strategy === :localized
        return LocalizedHost(budget=budget)
    elseif strategy === :infiltration_first
        return InfiltrationFirstHost(budget=budget)
    elseif strategy === :attrition
        return AttritionHost(budget=budget)
    else
        throw(ArgumentError("unknown host_strategy: $strategy"))
    end
end

function run_sample(sp::SampleParams; topology::Symbol, n_world::Int,
                    t_max::Int, cell_seed::Int, visibility::Float32=1.0f0,
                    principles::Principles=default_principles(),
                    host_strategy::Symbol=:random,
                    budget_decay::Float64=0.0,
                    accommodation_rate::Float64=0.0,
                    A_L::Float32=0f0, A_H::Float32=0f0)
    params = build_sample_params(sp; visibility=visibility, principles=principles,
                                  A_L=A_L, A_H=A_H)
    world  = build_world(; topology=topology, n=n_world, seed=cell_seed,
                           params=params, p_leak=sp.p_leak)
    if sp.fear_0 > 0
        for a in world.agents
            a.fear = sp.fear_0
        end
    end
    host = _make_host(host_strategy, sp.host_budget)
    if budget_decay > 0
        host.budget_decay = budget_decay
    end
    if accommodation_rate > 0
        host.accommodation_rate = accommodation_rate
    end
    part = natural_partition(topology, n_world)
    phi_hist  = Vector{Float32}(undef, t_max)
    comm_hist = Vector{Int}(undef, t_max)
    atk_hist  = Vector{Int}(undef, t_max)
    for t in 1:t_max
        snap = step!(world, host; partition=part)
        phi_hist[t]  = snap.Phi
        comm_hist[t] = snap.n_committed
        atk_hist[t]  = snap.n_attacked
    end
    return (
        equilibrium    = classify(phi_hist, committed_count(world)),
        t_to_attractor = time_to_attractor(phi_hist),
        final_phi      = phi_hist[end],
        history        = phi_hist,
        comm_history   = comm_hist,
        atk_history    = atk_hist,
        A_L            = A_L,
        A_H            = A_H,
    )
end

# --- main entry point --------------------------------------------------------

"""
    estimate_basins(; tier=:smoke, topology=:modular_cells, seed=0,
                      visibility=1.0f0, out_dir=nothing)

Runs `BASIN_TIERS[tier].n_samples` Latin-hypercube samples, classifies
each, accumulates counts. Returns a NamedTuple with:

- `counts::Dict{Symbol,Int}` — per-equilibrium count
- `t_attractors::Dict{Symbol,Vector{Int}}` — list of T-to-attractor
- `final_phis::Dict{Symbol,Vector{Float32}}`
- `n_samples`, `tier`, `topology`, `visibility`, `seed`, `elapsed_s`

If `out_dir` is provided, also writes `<out_dir>/basins.tsv`.
"""
function estimate_basins(; tier::Symbol=:smoke,
                           topology::Symbol=:modular_cells,
                           seed::Int=0,
                           visibility::Float32=1.0f0,
                           principles::Principles=default_principles(),
                           host_strategy::Symbol=:random,
                           budget_decay::Float64=0.0,
                           accommodation_rate::Float64=0.0,
                           A_L::Float32=0f0, A_H::Float32=0f0,
                           classifier::Symbol=:v1,
                           out_dir::Union{Nothing,String}=nothing)
    classifier ∈ (:v1, :v2) ||
        error("classifier must be :v1 or :v2 (got $classifier)")
    haskey(BASIN_TIERS, tier) ||
        error("unknown basin tier: $tier (expected :smoke, :standard, :nightly)")
    cfg = BASIN_TIERS[tier]
    rng = MersenneTwister(seed)
    samples = lhs(cfg.n_samples, 5; rng=rng)

    EQ = classifier === :v1 ?
         (:extinction, :dormant_persistence, :parallel_coexistence, :unclassified) :
         (:extinction, :dormant_persistence, :parallel_coexistence,
          :host_accommodation, :institutional_hybridization,
          :transformative_replacement, :ai_irrelevance, :unclassified)
    counts       = Dict{Symbol,Int}(eq => 0 for eq in EQ)
    t_attractors = Dict{Symbol,Vector{Int}}(eq => Int[] for eq in EQ)
    final_phis   = Dict{Symbol,Vector{Float32}}(eq => Float32[] for eq in EQ)

    t0 = time()
    for i in 1:cfg.n_samples
        sp = map_lhs(@view samples[i, :])
        r  = run_sample(sp; topology=topology, n_world=cfg.n_world,
                        t_max=cfg.t_max, cell_seed=seed + 1000 + i,
                        visibility=visibility, principles=principles,
                        host_strategy=host_strategy,
                        budget_decay=budget_decay,
                        accommodation_rate=accommodation_rate,
                        A_L=A_L, A_H=A_H)
        eq = classifier === :v1 ? r.equilibrium :
             classify_v2(r.history, r.comm_history, r.atk_history, cfg.n_world;
                          A_L=A_L, A_H=A_H)
        counts[eq] += 1
        push!(t_attractors[eq], r.t_to_attractor)
        push!(final_phis[eq],   r.final_phi)
    end
    elapsed = time() - t0

    if out_dir !== nothing
        mkpath(out_dir)
        out_path = joinpath(out_dir, "basins.tsv")
        open(out_path, "w") do io
            println(io, "equilibrium\tcount\tfraction\tmean_T_to_attractor\tmean_final_phi")
            for eq in EQ
                c    = counts[eq]
                frac = c / cfg.n_samples
                mt   = isempty(t_attractors[eq]) ? NaN : mean(t_attractors[eq])
                mphi = isempty(final_phis[eq])   ? NaN : mean(final_phis[eq])
                print(io, eq, "\t", c, "\t")
                print(io, round(frac; digits=4), "\t")
                print(io, round(mt;   digits=2), "\t")
                print(io, round(mphi; digits=4), "\n")
            end
        end
    end

    return (
        counts       = counts,
        t_attractors = t_attractors,
        final_phis   = final_phis,
        n_samples    = cfg.n_samples,
        tier         = tier,
        topology     = topology,
        visibility   = visibility,
        seed         = seed,
        elapsed_s    = elapsed,
    )
end
