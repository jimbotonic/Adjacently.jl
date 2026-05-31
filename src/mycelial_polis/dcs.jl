# Distributed constitutional sensing (DCS) — paper 2.
#
# Adds a per-cell sensor + activation layer that lets each cell of the
# polis infer its threat regime from local observations and activate
# the appropriate constitutional principle as a reflex, without
# requiring a global oracle.
#
# Phase 1 (this file): infrastructure only. CellSensors + CellPrinciples
# structs, per-step update step that computes the 6 local signals, and
# a per-cell activation step that maps signals to principle strengths.
# No simulator behaviour change yet — the per-cell principle strengths
# are computed but NOT read by the dynamics hooks. Phase 2 wires NM,
# Phase 3 wires the remaining hooks.

# --- Per-cell sensor struct ------------------------------------------------

"""
    CellSensors

Per-cell rolling sensor signals. Each field is an EMA of an observable
that the cell can detect from its own members' state, without needing
a global view.

- `recruitment_slow`: 1 - (fraction of new outsiders in this cell that
  reached :user this step) — high when newcomers fail to advance.
- `attack_intensity`: per-cell removal rate per step / N_cell.
- `faction_diversity`: number of distinct factions in this cell /
  n_factions_cap.
- `trust_decay`: edges lost from this cell's G_S subgraph / initial.
- `steward_dominance`: max steward outdegree / total cell edges.
- `resource_pressure`: low treasury / threshold (E3 hook; 0 for now).

EMA decay is `sensor_ema_decay` (default 0.85) per step.
"""
Base.@kwdef mutable struct CellSensors
    cell_id::Int
    members::Vector{Int}
    initial_edges::Int             = 0
    recruitment_slow::Float32      = 0f0
    attack_intensity::Float32      = 0f0
    faction_diversity::Float32     = 0f0
    trust_decay::Float32           = 0f0
    steward_dominance::Float32     = 0f0
    resource_pressure::Float32     = 0f0
    # H-5: infrastructure-capture sensor — local Herfindahl of replicas
    # held by stewards in this cell. 1.0 = one steward in this cell
    # holds all replicas held by this cell's stewards; 0 = uniform.
    infra_concentration::Float32   = 0f0
    # Counters reset each step
    _new_outsiders_this_step::Int  = 0
    _advanced_this_step::Int       = 0
    _attacks_this_step::Int        = 0
end

"""
    CellPrinciples

Per-cell principle strengths, in `[0, 1]`. Activated by `cell_activate!`
from the cell's sensors. The dynamics hooks read these instead of the
global `params.principles` when `params.dcs_enabled == true`.

Phase 1: computed but not yet read by hooks (NM hook reads global).
"""
Base.@kwdef mutable struct CellPrinciples
    cell_id::Int
    normative_minimalism::Float32          = 0f0
    transformative_non_absorption::Float32 = 0f0
    forkability::Float32                   = 0f0
    reversible_governance::Float32         = 0f0
    non_domination::Float32                = 0f0
    controlled_permeability::Float32       = 0f0
    bounded_legibility::Float32            = 0f0
    functional_autonomy::Float32           = 0f0
    redundant_stewardship::Float32         = 0f0
    adaptive_memory::Float32               = 0f0
end

"""
    DCSState

World-level container for the DCS layer.

- `cell_of[agent_id] = cell_id`
- `sensors[cell_id]`, `principles[cell_id]`
- `communication_reliability ∈ [0, 1]`: fraction of cell observations
  the activation step "sees" (degraded communication). 1.0 = perfect.
"""
Base.@kwdef mutable struct DCSState
    cell_of::Vector{Int}
    sensors::Dict{Int, CellSensors}
    principles::Dict{Int, CellPrinciples}
    communication_reliability::Float32 = 1f0
    activation_mode::Symbol = :adaptive
    # `:adaptive`  — sensors → activation logic per cell
    # `:global`    — every cell uses params.principles
    # `:universal` — every cell sets NM=TNA=1 always (pre-positioned default)
    # `:hybrid`    — universal NM+TNA floor of 0.5 + adaptive top-up (DCS-4
    #                tests this as a sense-attack defence: even if `q` collapses
    #                and the adaptive layer goes silent, the NM/TNA floor
    #                remains pre-positioned and active.)
    # `:off`       — sensors run but principles stay 0 (debug)
    sensor_ema_decay::Float32 = 0.85f0
    # H-5: when true, the adaptive activation logic uses
    # `steward_dominance` + `infra_concentration` to fire the anti-capture
    # reflex strongly (redundant_stewardship + non_domination +
    # reversible_governance). When false (default), centralization sensing
    # is weak and only `steward_dominance` modulates non_domination via
    # the pre-H-5 mapping. Decoupling lets DCS-3-style runs measure the
    # pre-H-5 baseline cleanly.
    h_sensing::Bool = false
    # DCS-2: fraction of cells that activate the WRONG principle each
    # step. Misclassified cells set `controlled_permeability = 1` and
    # zero all other principles — this is the "tighten openness when
    # conflict hits" misclassification, the cleanest case to test
    # because CP is known to HURT under conflict (matrix Δ = −0.03).
    misclassification_rate::Float32 = 0f0
    # When true, misclassified cells use the WORST principle for the
    # cell's dominant signal (per the principle-threat matrix). This
    # models an adversary that can corrupt the local diagnosis. When
    # false, all misclassified cells simply set CP = 1.
    misclassification_adversarial::Bool = false
end

"""
    build_dcs_state(world, partition; activation_mode=:adaptive,
                     communication_reliability=1.0)

Initialise DCS state for an existing world given a partition vector
(agent_id → cell_id). Computes initial per-cell edge counts.
"""
function build_dcs_state(world::World, partition::Vector{Int};
                          activation_mode::Symbol = :adaptive,
                          communication_reliability::Real = 1f0,
                          sensor_ema_decay::Real = 0.85f0,
                          misclassification_rate::Real = 0f0,
                          misclassification_adversarial::Bool = false,
                          h_sensing::Bool = false)
    g_S = world.multiplex.layers[:S]
    cells = unique(partition)
    sensors    = Dict{Int, CellSensors}()
    principles = Dict{Int, CellPrinciples}()
    members_by_cell = Dict{Int, Vector{Int}}()
    for c in cells
        members_by_cell[c] = Int[]
    end
    for (id, c) in enumerate(partition)
        push!(members_by_cell[c], id)
    end
    for c in cells
        # Count edges fully internal to the cell
        n_int = 0
        for u in members_by_cell[c]
            for v in outneighbors(g_S, eltype(g_S)(u))
                partition[Int(v)] == c && (n_int += 1)
            end
        end
        sensors[c]    = CellSensors(cell_id = c,
                                     members = copy(members_by_cell[c]),
                                     initial_edges = n_int)
        principles[c] = CellPrinciples(cell_id = c)
    end
    return DCSState(cell_of = copy(partition),
                     sensors = sensors,
                     principles = principles,
                     communication_reliability = Float32(communication_reliability),
                     activation_mode = activation_mode,
                     sensor_ema_decay = Float32(sensor_ema_decay),
                     misclassification_rate = Float32(misclassification_rate),
                     misclassification_adversarial = misclassification_adversarial,
                     h_sensing = h_sensing)
end

# --- Sensor update ---------------------------------------------------------

"""
    _update_dcs_sensors!(world, dcs, snapshot)

Called each step after `step!` finishes. Computes per-cell signals from
the world's current state and the just-finished snapshot, then EMAs them
into the sensors.

Communication reliability degrades observations: with probability
`1 - q` each per-cell signal is replaced by its previous-step value
(stale data instead of fresh).
"""
function _update_dcs_sensors!(world::World, dcs::DCSState, snap::NamedTuple)
    g_S = world.multiplex.layers[:S]
    q   = dcs.communication_reliability
    α   = 1f0 - dcs.sensor_ema_decay   # EMA mixing rate
    n_factions_cap = Float32(get(world.params, :n_factions_cap, 4))

    for (c, sens) in dcs.sensors
        n_mem = length(sens.members)
        n_mem == 0 && continue

        # Active/committed counts in this cell
        n_active = 0
        n_committed = 0
        n_advanced = 0
        n_attacks = 0
        max_steward_out = 0
        factions_in_cell = Set{Symbol}()
        n_edges_int = 0
        for id in sens.members
            a = world.agents[id]
            if ROLE_RANK[a.role] >= 1     # observer+
                n_active += 1
            end
            if ROLE_RANK[a.role] >= ROLE_RANK[:user]
                n_committed += 1
                if a.faction !== :none
                    push!(factions_in_cell, a.faction)
                end
            end
            if a.role === :user || a.role === :contributor || a.role === :steward
                # Approx: count members who reached :user-or-above
                n_advanced += 1
            end
            if a.role === :steward
                od = outdegree(g_S, eltype(g_S)(id))
                max_steward_out = max(max_steward_out, od)
            end
            if a.role === :removed
                n_attacks += 1
            end
            # Internal edges count (out only — symmetric undercount but
            # consistent)
            for v in outneighbors(g_S, eltype(g_S)(id))
                dcs.cell_of[Int(v)] == c && (n_edges_int += 1)
            end
        end

        # --- compute fresh signals ---
        # recruitment_slow: 1 - n_advanced/n_mem; high when few have advanced
        sig_recruit_slow = Float32(1 - n_advanced / max(n_mem, 1))
        # attack_intensity: cumulative :removed in this cell / n_mem
        sig_attack_int   = Float32(n_attacks / max(n_mem, 1))
        # faction_diversity: factions in cell / cap
        sig_faction_div  = Float32(length(factions_in_cell) / n_factions_cap)
        # trust_decay: 1 - current_edges / initial_edges
        sig_trust_decay  = sens.initial_edges == 0 ? 0f0 :
                           clamp(Float32(1 - n_edges_int / sens.initial_edges), 0f0, 1f0)
        # steward_dominance: max_steward_out / max(total_cell_edges, 1)
        sig_steward_dom  = sens.initial_edges == 0 ? 0f0 :
                           Float32(max_steward_out / max(sens.initial_edges, 1))
        # resource_pressure: hook for E3 (treasury layer); 0 for now
        sig_resource     = 0f0
        # infra_concentration: local Herfindahl-style — among stewards in
        # this cell, what fraction of replicas held by them is held by the
        # single most-loaded steward. Fires when one steward in this cell
        # hoards replicas across functions (the H-2 shock signature).
        stewards_in_cell = [id for id in sens.members
                             if world.agents[id].role === :steward]
        sig_infra_conc = 0f0
        if !isempty(stewards_in_cell)
            rcount = Dict{Int, Int}(id => 0 for id in stewards_in_cell)
            for (_, replicas) in world.infra
                for sid in replicas
                    if haskey(rcount, sid)
                        rcount[sid] += 1
                    end
                end
            end
            total_in_cell = sum(values(rcount))
            if total_in_cell > 0
                sig_infra_conc = Float32(maximum(values(rcount)) / total_in_cell)
            end
        end

        # Apply communication noise: with probability (1-q), skip this
        # cell's fresh update (use last-step value).
        if rand(world.rng) < q
            sens.recruitment_slow   = (1-α) * sens.recruitment_slow + α * sig_recruit_slow
            sens.attack_intensity   = (1-α) * sens.attack_intensity + α * sig_attack_int
            sens.faction_diversity  = (1-α) * sens.faction_diversity + α * sig_faction_div
            sens.trust_decay        = (1-α) * sens.trust_decay + α * sig_trust_decay
            sens.steward_dominance  = (1-α) * sens.steward_dominance + α * sig_steward_dom
            sens.resource_pressure  = (1-α) * sens.resource_pressure + α * sig_resource
            sens.infra_concentration = (1-α) * sens.infra_concentration + α * sig_infra_conc
        end
    end
end

# --- Activation logic ------------------------------------------------------

"""
    _cell_activate!(dcs)

Maps the EMA'd sensor signals to per-cell principle strengths. Each
principle's strength is the clamp of one or two normalised signals.
Conservative thresholds so principles activate only when a signal
is meaningfully present.

In `:universal` mode: always set NM=TNA=1 (the pre-positioned default
class). In `:global` mode: copy from `params.principles`. In `:off`
mode: leave all at zero.
"""
function _cell_activate!(world::World, dcs::DCSState)
    if dcs.activation_mode === :off
        return
    end
    if dcs.activation_mode === :universal
        for (_, p) in dcs.principles
            p.normative_minimalism = 1f0
            p.transformative_non_absorption = 1f0
        end
        return
    end
    if dcs.activation_mode === :global
        pr = get_principles(world.params)
        for (_, p) in dcs.principles
            for f in fieldnames(CellPrinciples)
                f === :cell_id && continue
                setfield!(p, f, getfield(pr, f))
            end
        end
        return
    end
    # :adaptive and :hybrid both run the sensor-mapping logic below;
    # hybrid additionally enforces an NM+TNA floor at the end of each
    # cell's update.
    # :hybrid mode — universal NM+TNA floor + adaptive top-up. Compute
    # the adaptive activations first by setting mode = :adaptive
    # transiently, then enforce the floor.
    hybrid = dcs.activation_mode === :hybrid
    # :adaptive mode — signals → principles, with optional misclassification
    mc_rate = dcs.misclassification_rate
    for (c, sens) in dcs.sensors
        p = dcs.principles[c]
        if mc_rate > 0 && rand(world.rng) < mc_rate
            # Zero all principles, then activate the wrong one(s).
            p.normative_minimalism = 0f0
            p.transformative_non_absorption = 0f0
            p.forkability = 0f0
            p.reversible_governance = 0f0
            p.non_domination = 0f0
            p.functional_autonomy = 0f0
            p.redundant_stewardship = 0f0
            p.bounded_legibility = 0f0
            p.controlled_permeability = 0f0
            p.adaptive_memory = 0f0
            if dcs.misclassification_adversarial
                # Adversarial: pick the WORST principle for the
                # dominant local signal, per the matrix.
                dom_signal = argmax((sens.faction_diversity,
                                      sens.attack_intensity,
                                      sens.recruitment_slow,
                                      sens.steward_dominance))
                if dom_signal == 1            # conflict → CP (Δ -0.03)
                    p.controlled_permeability = 1f0
                elseif dom_signal == 2        # ext pressure → adaptive_memory (Δ -0.10)
                    p.adaptive_memory = 1f0
                elseif dom_signal == 3        # recruitment_slow → bounded_legibility (Δ ~0)
                    p.bounded_legibility = 1f0
                else                           # steward dominance → CP
                    p.controlled_permeability = 1f0
                end
            else
                # Random-wrong: always CP=1 (default DCS-2 model).
                p.controlled_permeability = 1f0
            end
            continue
        end
        # Renewal bottleneck → NM (the universal default)
        p.normative_minimalism = clamp(sens.recruitment_slow, 0f0, 1f0)
        # External pressure → TNA
        p.transformative_non_absorption = clamp(sens.attack_intensity * 3f0, 0f0, 1f0)
        # Internal conflict → forkability + reversible_governance
        p.forkability = clamp(sens.faction_diversity * 1.5f0, 0f0, 1f0)
        p.reversible_governance = clamp((sens.faction_diversity * 1.5f0 +
                                          sens.steward_dominance * 2f0) / 2, 0f0, 1f0)
        # Capture → non_domination
        p.non_domination = clamp(sens.steward_dominance * 2f0, 0f0, 1f0)
        # Resources → functional_autonomy + redundant_stewardship (E3)
        p.functional_autonomy = clamp(sens.resource_pressure, 0f0, 1f0)
        p.redundant_stewardship = clamp(sens.resource_pressure * 0.7f0, 0f0, 1f0)
        # H-5: hierarchy-sensing reflex. When `h_sensing`, the
        # `steward_dominance` (G_S edges) and `infra_concentration`
        # (replica hoarding) sensors fire the anti-capture reflex per
        # H-2: redundant_stewardship is strongest, non_domination and
        # reversible_governance secondary. infra_concentration only fires
        # in the cell that contains the captured steward, so the global
        # aggregation (max across cells) lifts the principle for refill.
        if dcs.h_sensing
            h_signal = max(sens.steward_dominance, sens.infra_concentration)
            p.redundant_stewardship = max(p.redundant_stewardship,
                                          clamp(h_signal * 2f0, 0f0, 1f0))
            p.non_domination = max(p.non_domination,
                                   clamp(h_signal * 1.5f0, 0f0, 1f0))
            p.reversible_governance = max(p.reversible_governance,
                                          clamp(h_signal * 1.5f0, 0f0, 1f0))
        end
        # Infiltration proxy (cell trust decay) → bounded_legibility / CP
        p.bounded_legibility = clamp(sens.trust_decay, 0f0, 1f0)
        p.controlled_permeability = clamp(0.5f0 - sens.trust_decay, 0f0, 1f0)
        # Adaptive memory: stub for now (needs cycle-detection logic)
        p.adaptive_memory = 0f0
        # Hybrid mode: enforce NM and TNA floor of 0.5 even when sensors
        # are silent. Adaptive activation may push them higher.
        if hybrid
            p.normative_minimalism = max(p.normative_minimalism, 0.5f0)
            p.transformative_non_absorption = max(p.transformative_non_absorption, 0.5f0)
        end
    end
end

"""
    global_principle_max(dcs, field) -> Float32

Aggregator for principles that have a global mechanical effect (e.g.,
`redundant_stewardship` modulates the *global* infrastructure refill
target; `non_domination` triggers a single global rotation event each
step). When the polis is governed by per-cell principles, we take the
**max** across cells: if any cell senses capture, the global anti-
capture reflex fires. Returns 0 when DCS is off.
"""
function global_principle_max(dcs::Union{Nothing, DCSState}, field::Symbol)
    dcs === nothing && return 0f0
    dcs.activation_mode === :off && return 0f0
    m = 0f0
    for (_, p) in dcs.principles
        v = getfield(p, field)
        v > m && (m = v)
    end
    return m
end

"""
    cell_principles_for(world, dcs, agent_id) -> Principles

Returns the effective Principles struct for an agent. Used by dynamics
hooks that want per-agent (= per-cell) principle strengths. Falls back
to global `params.principles` when `dcs === nothing`.
"""
function cell_principles_for(world::World, dcs::Union{Nothing, DCSState},
                              agent_id::Int)
    dcs === nothing && return get_principles(world.params)
    c = dcs.cell_of[agent_id]
    cp = dcs.principles[c]
    return Principles(;
        normative_minimalism = cp.normative_minimalism,
        transformative_non_absorption = cp.transformative_non_absorption,
        forkability = cp.forkability,
        reversible_governance = cp.reversible_governance,
        non_domination = cp.non_domination,
        controlled_permeability = cp.controlled_permeability,
        bounded_legibility = cp.bounded_legibility,
        functional_autonomy = cp.functional_autonomy,
        redundant_stewardship = cp.redundant_stewardship,
        adaptive_memory = cp.adaptive_memory,
    )
end
