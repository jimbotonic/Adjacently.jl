"""
    Principles

The ten minimal internal attractor principles of roadmap §19, as
hyperparameters in `[0,1]`. Convention everywhere: **0 = absent, 1 =
strong**. Each principle is **monotone** in the metric it most directly
affects — see the table below.

| Field | Hook (where it's applied) | Direction (↑ strength →) |
|---|---|---|
| `non_domination`                | `step!` — rotates stewards back to contributors every `max(3, round(30/(1+9·s)))` steps | `Γ` ↓ slightly, `Φ` stable (anti-capture) |
| `functional_autonomy`           | `_maybe_advance_role!` — boosts commitment of new entrants by `1 + s · stage_index` | `Φ` ↑ (more durable commitment) |
| `bounded_legibility`            | `drift_legibility!` — multiplies `hide_rate` by `1 + 4·s` | `Λ_sat` ↑ (privacy tools) |
| `reversible_governance`         | `step!` — defectors revert to `:contributor` with prob `s` per step | `Φ` ↑ (less permanent attrition) |
| `forkability`                   | `update_repression!` — defection is *prevented* with prob `s` (agent stays committed) | `Φ` ↑ (no schism) |
| `redundant_stewardship`         | `refill_infrastructure!` — refill target = `infra_min + round(3·s)`, while `Ψ_T` threshold stays at `infra_min` | `Ψ_T` ↑ (spare capacity) |
| `normative_minimalism`          | `_maybe_advance_role!` — stage thresholds scaled by `1 − 0.5·s` (thinner constitution = easier advancement) | `Φ` ↑ (more recruits) |
| `adaptive_memory`               | `update_repression!` — `narrative_scale` scaled by `1 − 0.5·s` (saturates faster) | `Φ` ↑ (better backfire response) |
| `controlled_permeability`       | `apply_adoption!` — `α` scaled by `0.7 + 0.6·s`; `drift_legibility!` — `leak_rate` scaled by `1.0 + 0.5·s` | interior optimum; mostly `Φ` ↑ |
| `transformative_non_absorption` | `update_repression!` — backfire `b_i` scaled by `0.5 + 1.5·s` | `Φ` ↑ (transformation amplifies backfire) |

All hooks gracefully degrade to the item-2/-4/-6/-7 behaviour when
`get(world.params, :principles, default_principles())` returns the
defaults (all 0.5).

NAMING HAZARD (flagged Tier-1 E4, 2026-07-30): two distinct mechanisms share
the name "adaptive_memory". (A) the `adaptive_memory` field below — the
static AM *knob*: a global principle strength that compresses narrative
saturation (hook in `dynamics.jl`) and, additionally, gates mechanism (B).
(B) the DCS per-cell AM *memory layer* (`DCSState.memory` in `dcs.jl`): a
decaying per-cell, per-principle strength floor for reactivation.
Ablations/matrix sweeps act on (A) via `ALL_PRINCIPLE_NAMES`; (B) is
internal to the DCS.
"""
Base.@kwdef struct Principles
    non_domination::Float32                = 0.5f0
    functional_autonomy::Float32           = 0.5f0
    bounded_legibility::Float32            = 0.5f0
    reversible_governance::Float32         = 0.5f0
    forkability::Float32                   = 0.5f0
    redundant_stewardship::Float32         = 0.5f0
    normative_minimalism::Float32          = 0.5f0
    adaptive_memory::Float32               = 0.5f0
    controlled_permeability::Float32       = 0.5f0
    transformative_non_absorption::Float32 = 0.5f0
    # N-1 — defence against narrative coalition. When > 0, faction-change
    # attempts (coalition narrative conversion + schism reassignment) are
    # rejected with probability proportional to the *target faction's*
    # share of committed agents in the target's cell. Logic: the bigger
    # a faction already is, the harder it is for it to grow further —
    # which forces narrative attackers to spread their conversions across
    # cells instead of saturating one.
    faction_diversity_floor::Float32       = 0f0
end

"""
    default_principles()

Defaults all at 0.5 so removing one (setting to 0) is a meaningful
baseline ablation in item 10's OAT sweep.
"""
default_principles() = Principles()

"""
    zero_principles()

All principles set to 0 (absent). Useful when a test wants the bare
§9 / §12 / §8 dynamics without any principle modulation — e.g. for
isolating Experiment 4 (§23) from the baseline `0.5` modulation that
item 9 introduces in `default_principles()`.
"""
zero_principles() = Principles(
    non_domination                = 0f0,
    functional_autonomy           = 0f0,
    bounded_legibility            = 0f0,
    reversible_governance         = 0f0,
    forkability                   = 0f0,
    redundant_stewardship         = 0f0,
    normative_minimalism          = 0f0,
    adaptive_memory               = 0f0,
    controlled_permeability       = 0f0,
    transformative_non_absorption = 0f0,
    faction_diversity_floor       = 0f0,
)

"""
    get_principles(params) -> Principles

Read the `Principles` block out of `params`, returning the defaults if
not present. All hooks use this so they continue to work on worlds built
before item 9 was introduced (any existing `default_params()` merge that
doesn't explicitly carry a `principles` field).
"""
get_principles(params) = get(params, :principles, default_principles())

"""
    principle_strength(params, name::Symbol) -> Float32

Shortcut for `getfield(get_principles(params), name)` clamped to `[0,1]`.
"""
principle_strength(params, name::Symbol) =
    clamp(getfield(get_principles(params), name), 0f0, 1f0)
