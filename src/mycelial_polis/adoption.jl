"""
    adoption

Adoption rules for the Mycelial Polis (roadmap §8). Two rules are exposed,
selectable at runtime via `params.adoption_rule`:

- `:threshold`         — Granovetter / Watts. Continuous weighted exposure
  accumulates into a per-step "drive" signal; role advances when the drive
  clears the stage-specific threshold.
- `:complex_contagion` — Centola & Macy (2007). Risky / costly behaviours
  spread only via *redundant* social reinforcement. Role advancement is
  gated on the agent having ≥ `params.k_cc` distinct in-neighbours
  simultaneously at role `≥ :user` with cumulative trust weight `≥
  params.theta_cc`. When the gate is open, the same drive computation as
  `:threshold` is used (so on a saturated neighbourhood the two rules
  converge); when the gate is closed, no role advancement happens, only
  `awareness` accumulates from continuous exposure.

`apply_adoption!(world, eS, eC)` is the dispatch entry point called by
`step!`. It mutates `world.agents` in place.

Roadmap refs: §8 (adoption dynamics), §25.1 (Granovetter), §25.2 (Watts),
§25.3 (Centola–Macy), §25.4 (decade review).
"""

# Dispatch entry point — picks the rule from `world.params.adoption_rule`.
function apply_adoption!(world::World, exposures_S::Vector{Float32},
                         exposures_C::Vector{Float32})
    rule = world.params.adoption_rule
    if rule === :threshold
        threshold_adopt!(world, exposures_S, exposures_C)
    elseif rule === :complex_contagion
        complex_contagion_adopt!(world, exposures_S, exposures_C)
    else
        throw(ArgumentError("unknown adoption_rule: $rule " *
                            "(expected :threshold or :complex_contagion)"))
    end
    return world
end

# --- threshold rule ----------------------------------------------------------

"""
    threshold_adopt!(world, exposures_S, exposures_C)

Granovetter-style staged threshold adoption (§8). The "drive" signal is

    drive_i = (α · e_S + β · e_C) / deg_norm − δ·r_i + η·b_i + χ·I_i

with `deg_norm = max(outdeg(g_S, i) + indeg(g_S, i), 1)` so per-stage
thresholds are scale-invariant. Each stage transition gates on its own
threshold via `_maybe_advance_role!`.
"""
function threshold_adopt!(world::World, exposures_S::Vector{Float32},
                          exposures_C::Vector{Float32})
    p   = world.params
    pr_global = get_principles(p)
    g_s = layer(world.multiplex, :S)
    dcs_active = world.dcs !== nothing
    # Saturation cap: if `growth_ceiling_frac < 1.0` and the committed
    # pool has already hit the ceiling, gate further promotions into
    # the committed bracket (computed once per step from start-of-step
    # state; the one-step lag is benign).
    gc = Float64(get(p, :growth_ceiling_frac, 1.0))
    cap_reached = false
    if gc < 1.0
        n_committed_start = count(is_committed, world.agents)
        cap_reached = n_committed_start >= round(Int, gc * length(world.agents))
    end
    disable_nm = Bool(get(p, :disable_nm_effect, false))
    @inbounds for a in world.agents
        a.role === :removed && continue
        # Per-cell principle resolution when DCS is active.
        pr = dcs_active ? cell_principles_for(world, world.dcs, a.id) : pr_global
        drive = _drive(a, g_s, exposures_S, exposures_C, p, pr)
        a.awareness = clamp(0.7f0 * a.awareness + 0.3f0 * drive, 0f0, 1f0)
        _maybe_advance_role!(a, drive, p.stage_thresholds; principles=pr,
                              disable_nm_effect = disable_nm,
                              block_committed_promotion = cap_reached)
    end
    return world
end

# --- complex-contagion rule --------------------------------------------------

"""
    complex_contagion_adopt!(world, exposures_S, exposures_C)

Centola–Macy complex contagion (§8). Role advancement requires the agent
to be in the *socially-reinforcing regime*:

    m_i      = number of distinct in-neighbours with role ≥ :user
    trust_i  = Σ commitment_j over those in-neighbours
    gate     = (m_i ≥ k_cc) ∧ (trust_i ≥ θ_cc)

When `gate` is open, the drive is computed exactly as `:threshold` (so on
a saturated neighbourhood the two rules behave identically). When closed,
role does not advance, but `awareness` still accumulates from continuous
exposure — agents on the wrong side of the gate slowly become aware of the
polis without committing to it.

The qualitative Centola–Macy effect (complex stalls on tree-like
topologies where threshold spreads) is reproduced by
`test/run_tests_mp_adoption.jl`.
"""
function complex_contagion_adopt!(world::World, exposures_S::Vector{Float32},
                                  exposures_C::Vector{Float32})
    p   = world.params
    pr_global = get_principles(p)
    g_s = layer(world.multiplex, :S)
    k_cc     = p.k_cc
    theta_cc = p.theta_cc
    user_rank = ROLE_RANK[:user]
    dcs_active = world.dcs !== nothing
    gc = Float64(get(p, :growth_ceiling_frac, 1.0))
    cap_reached = false
    if gc < 1.0
        n_committed_start = count(is_committed, world.agents)
        cap_reached = n_committed_start >= round(Int, gc * length(world.agents))
    end
    disable_nm = Bool(get(p, :disable_nm_effect, false))
    @inbounds for a in world.agents
        a.role === :removed && continue
        pr = dcs_active ? cell_principles_for(world, world.dcs, a.id) : pr_global
        drive = _drive(a, g_s, exposures_S, exposures_C, p, pr)
        a.awareness = clamp(0.7f0 * a.awareness + 0.3f0 * drive, 0f0, 1f0)
        # Centola–Macy gate.
        m_i = 0
        trust = 0f0
        for v in inneighbors(g_s, a.id)
            nb = world.agents[Int(v)]
            ROLE_RANK[nb.role] >= user_rank || continue
            m_i  += 1
            trust += nb.commitment
        end
        (m_i >= k_cc && trust >= theta_cc) || continue
        _maybe_advance_role!(a, drive, p.stage_thresholds; principles=pr,
                              disable_nm_effect = disable_nm,
                              block_committed_promotion = cap_reached)
    end
    return world
end

# --- shared drive computation -----------------------------------------------

# Per-agent "drive" used by both rules. Kept inline-able for the hot loop.
# §19.9 Controlled permeability — scales α so a more permeable polis
# absorbs more exposure per unit degree. Hook is neutral at s=0 (α
# unchanged) and at s=1 amplifies α by 1.5×.
@inline function _drive(a::Agent, g_s::AbstractGraph,
                        exposures_S::Vector{Float32},
                        exposures_C::Vector{Float32},
                        p::NamedTuple, pr::Principles)
    # §17 A_L compounds with §19.9 controlled_permeability — AI on the
    # latent side amplifies whatever neighbour exposure is already there.
    A_L = get(p, :A_L, 0f0)
    α_eff = p.alpha * (1f0 + 0.5f0 * pr.controlled_permeability) *
                       (1f0 + 0.5f0 * A_L)
    deg_s = max(outdegree(g_s, a.id) + indegree(g_s, a.id), 1)
    base  = (α_eff * exposures_S[a.id] + p.beta * exposures_C[a.id]) /
            Float32(deg_s)
    return base - p.delta * a.fear + p.eta * a.backfire + p.chi * a.identity
end
