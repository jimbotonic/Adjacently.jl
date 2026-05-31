"""
    step!(world, host) -> snapshot

Advance the world by one time step and return a metrics snapshot. The six
sub-steps follow item 2's spec / roadmap §5–§12 ordering:

    1. compute exposures on G_S and G_C
    2. apply host action H_t        (§10)
    3. update r_i (fear) and b_i (backfire) per §9
    4. apply adoption rule and recompute roles    (§8)
    5. update infrastructure replicas R_k(t)      (§12)
    6. snapshot metrics                            (§14)

Skeleton implementations are deliberately minimal. The eventual home of
the full §8 / §9 / §12 dynamics is items 4, 6, 7.
"""
function step!(world::World, host::HostStrategy;
               partition::Union{Nothing,Vector{Int}}=nothing)
    world.t += 1
    t = world.t
    p = world.params
    pr = get_principles(p)

    # Non-domination (§19.1) — rotate one steward back to :contributor every
    # `max(3, round(30/(1+9·s)))` steps. At s=0 effectively off (every ~30
    # steps with negligible effect); at s=1 a constant churn of stewards.
    # DCS aggregation: lift `s` to the per-cell max so a single cell that
    # senses capture triggers the rotation, matching the H-2 anti-capture
    # reflex semantics.
    nd_strength = max(pr.non_domination,
                       global_principle_max(world.dcs, :non_domination))
    nd_period = max(3, round(Int, 30 / (1 + 9 * nd_strength)))
    if nd_strength > 0 && t % nd_period == 0
        stewards = [a for a in world.agents if a.role === :steward]
        if !isempty(stewards)
            victim = stewards[1 + (t % length(stewards))]
            victim.role = :contributor
        end
    end

    # Open-population recruitment (default off): each :removed agent is
    # revived as a fresh :outsider with per-step probability
    # `params.recruitment_rate`. Without this the population is closed and
    # the host's cumulative attack budget eventually depletes everyone —
    # producing a 1/T artifactual phase boundary rather than a true
    # dynamical attractor.
    recruitment_rate = get(p, :recruitment_rate, 0.0)
    # `recruitment_mode` ∈ (:ambient, :sponsored, :bridge, :cellular).
    # Defaults to :ambient (random outsider revival), preserving
    # CN26 baseline behaviour. Legacy `recruitment_sponsored = true`
    # is auto-translated to :sponsored.
    recruitment_mode = get(p, :recruitment_mode, :ambient)
    if get(p, :recruitment_sponsored, false) && recruitment_mode === :ambient
        recruitment_mode = :sponsored
    end
    recruitment_cross_cell_p = get(p, :recruitment_cross_cell_p, 0.10)
    if recruitment_rate > 0
        g_S = world.multiplex.layers[:S]
        # For bridge mode we want to revive :removed agents with FEW
        # committed in-neighbours; for sponsored we want MANY; cellular
        # combines local concentration with a small cross-cell probability.
        # Count committed in-neighbours per :removed agent — a single pass.
        for a in world.agents
            a.role === :removed || continue
            n_committed_nbrs = 0
            if recruitment_mode !== :ambient
                for v in inneighbors(g_S, eltype(g_S)(a.id))
                    nb = world.agents[Int(v)]
                    if ROLE_RANK[nb.role] >= ROLE_RANK[:user]
                        n_committed_nbrs += 1
                    end
                end
            end
            # Per-mode admission filter.
            admit = if recruitment_mode === :ambient
                true
            elseif recruitment_mode === :sponsored
                n_committed_nbrs >= 1
            elseif recruitment_mode === :bridge
                # Targets depleted zones: only revive :removed with ≤ 1
                # committed in-neighbour. Repopulates cleared regions
                # without raising local L_ext (the new agent enters as
                # :outsider, not as a known sympathiser).
                n_committed_nbrs <= 1
            elseif recruitment_mode === :cellular
                # Local concentration (sponsored-like) plus a small cross-
                # cell seeding probability. Without partition info we
                # approximate "cross-cell" by reviving isolated :removed
                # agents with probability `recruitment_cross_cell_p`.
                n_committed_nbrs >= 1 ||
                    rand(world.rng) < recruitment_cross_cell_p
            else
                error("unknown recruitment_mode: $recruitment_mode")
            end
            admit || continue
            if rand(world.rng) < recruitment_rate
                a.role       = :outsider
                a.awareness  = 0f0
                a.commitment = 0f0
                a.fear       = 0f0
                a.backfire   = 0f0
                a.capacity   = 1f0
                a.L_ext      = 0f0
                a.L_int      = 0f0
            end
        end
    end

    # E2 endogenous conflict (paper 2). Off when disagreement_rate = 0.
    # Three sub-steps: (a) decay all per-agent disagreement; (b) maybe
    # fire one disagreement event (faction birth or propagation, raising
    # disagreement on the involved agents); (c) schism for any agent whose
    # disagreement crosses `schism_threshold` (flip faction, trim cross-
    # faction trust edges with probability `trust_decay_on_fight`).
    _disagreement_step!(world, p)

    # 1. exposures — count adopted neighbors on G_S and G_C.
    exposures_S, exposures_C = compute_exposures(world)

    # 2. host action + exogenous budget decay + endogenous accommodation
    # (item 8 v2: enables `:host_accommodation` and
    # `:institutional_hybridization` either schedule-driven or in response
    # to the host's revealed-preference withdrawal).
    # §17 A_H boost: host AI temporarily amplifies the budget that this
    # step's `act!` sees, without permanently raising the budget (so
    # decay still operates on the underlying value). All seven host
    # strategies read `host.budget` so this is sufficient.
    _resolve_rho_pressure!(host, world)
    A_H = get(p, :A_H, 0f0)
    saved_budget = host.budget
    if A_H > 0
        host.budget = round(Int, host.budget * (1f0 + 0.5f0 * A_H))
    end
    attacked_ids = act!(host, world, t)
    host.budget = saved_budget
    repressed_set = Set(attacked_ids)
    _decay_budget!(host)
    _accommodate!(host, committed_count(world))

    # 3. fear + backfire update (§9 full form, item 6) + defection trigger.
    n_defected, n_backfired = update_repression!(world, repressed_set)

    # Reversible governance (§19.4) — defectors revert to :contributor with
    # probability `s` per step. Implementation: scan defectors created so
    # far and randomly restore. (Cheap: at most `n_defectors` agents.)
    if pr.reversible_governance > 0
        for a in world.agents
            a.role === :defector || continue
            if rand(world.rng) < pr.reversible_governance
                a.role       = :contributor
                a.commitment = 0.6f0
                a.fear       = 0.3f0       # post-recovery still cautious
            end
        end
    end

    # 4. adoption + role advancement (§8 threshold skeleton).
    apply_adoption!(world, exposures_S, exposures_C)

    # 5. infrastructure replicas (§12 skeleton: promote committed agents
    #    to fill missing infra slots if any function is below m_k).
    refill_infrastructure!(world)

    # 6. legibility drift (§11) and snapshot.
    drift_legibility!(world)
    snap = snapshot(world; t=t, attacked=length(attacked_ids),
                    defected=n_defected, backfired=n_backfired,
                    partition=partition)

    # 7. Distributed constitutional sensing (paper 2): update per-cell
    # sensors from this step's events, then run the activation logic.
    # No-op when `world.dcs === nothing`.
    if world.dcs !== nothing
        _update_dcs_sensors!(world, world.dcs, snap)
        _cell_activate!(world, world.dcs)
    end
    return snap
end

# 1 — exposures ---------------------------------------------------------------

"""
    compute_exposures(world) -> (exposures_S, exposures_C)

Per-agent count of in-neighbors on `G_S` and `G_C` that are themselves at
or above the `:user` role (i.e. "actively adopted"). Item 4's complex
contagion will weight by trust; the skeleton uses a uniform count.
"""
function compute_exposures(world::World)
    g_s = layer(world.multiplex, :S)
    g_c = layer(world.multiplex, :C)
    n   = length(world.agents)
    eS  = zeros(Float32, n)
    eC  = zeros(Float32, n)
    @inbounds for a in world.agents
        is_committed(a) || continue
        for v in outneighbors(g_s, a.id)
            v <= n && (eS[v] += 1f0)
        end
        for v in outneighbors(g_c, a.id)
            v <= n && (eC[v] += 1f0)
        end
    end
    return (eS, eC)
end

# 3 — repression / backfire ---------------------------------------------------

"""
    update_repression!(world, repressed_set) -> (n_defected, n_backfired)

Full §9 dynamics (item 6). Per non-removed agent, compute four
independent backfire components and the fear update:

    R_i = max(𝟙[i ∈ repressed], fraction of :S in-neighbours repressed this step)
    J_i = fraction of :S in-neighbours currently :removed (cumulative injustice)
    V_i = visibility_scaling · V_global · (1 + broadcast_amp)
    A_i = narrative availability — rolling mean of n_backfired / narrative_scale

    r_i(t+1) = clip( r_i(t) + ρ · R_i − κ · P_i )
    b_i(t+1) = clip( w_R·R_i + w_J·J_i + w_V·V_i + w_A·A_i )

After the update, a committed agent (`role ≥ :user`) whose
`fear ≥ defect_fear_threshold` AND whose protection (§S-neighbour
committed density) is below `defect_protection_threshold` transitions to
`:defector` (rank −3). Defectors leave their infra replica seats.

The step's `n_backfired` (count of agents with `b_i > 0.5`) is pushed onto
`world.recent_backfire_count`, truncated to `params.narrative_window` —
that buffer is the only state behind `A_i` for the next step.
"""
function update_repression!(world::World, repressed_set::Set{Int})
    p   = world.params
    g_s = layer(world.multiplex, :S)
    g_o = layer(world.multiplex, :O)
    n   = length(world.agents)

    # Global visibility — fraction of N attacked this step, amplified by how
    # much the host's observation layer covers G_S.
    V_global = n == 0 ? 0f0 : Float32(length(repressed_set) / n)
    edges_S  = max(ne(g_s), 1)
    broadcast_amp = Float32(ne(g_o) / edges_S)

    # Narrative availability — rolling mean of recent n_backfired counts.
    # Adaptive memory (§19.8) compresses `narrative_scale`: stronger memory
    # = faster A_i saturation = better backfire response. §17 A_L
    # compounds the compression — latent AI accelerates narrative spread.
    pr = get_principles(p)
    A_L = get(p, :A_L, 0f0)
    nm_scale = p.narrative_scale * (1f0 - 0.5f0 * pr.adaptive_memory) *
                                    (1f0 - 0.5f0 * A_L)
    history = world.recent_backfire_count
    window  = p.narrative_window
    A_global = if isempty(history) || nm_scale <= 0
        0f0
    else
        recent = @view history[max(1, end - window + 1):end]
        clamp(Float32(sum(recent) / length(recent)) / nm_scale, 0f0, 1f0)
    end

    n_defected  = 0
    n_backfired = 0
    @inbounds for a in world.agents
        a.role === :removed && continue
        nbrs = inneighbors(g_s, a.id)
        if !isempty(nbrs)
            committed_nbrs = 0
            removed_nbrs   = 0
            unjust_nbrs    = 0
            for v in nbrs
                nb = world.agents[Int(v)]
                is_committed(nb) && (committed_nbrs += 1)
                nb.role === :removed && (removed_nbrs += 1)
                Int(v) ∈ repressed_set && (unjust_nbrs += 1)
            end
            protection   = Float32(committed_nbrs / length(nbrs))
            local_unjust = Float32(unjust_nbrs    / length(nbrs))
            J_i          = Float32(removed_nbrs   / length(nbrs))
        else
            protection = 0f0; local_unjust = 0f0; J_i = 0f0
        end

        was_repressed = a.id ∈ repressed_set
        R_i = max(was_repressed ? 1f0 : 0f0, local_unjust)
        V_i = p.visibility_scaling * V_global * (1f0 + broadcast_amp)

        a.fear = clamp(a.fear + p.rho * R_i - p.kappa * protection, 0f0, 1f0)
        # Transformative non-absorption (§19.10) — backfire amplification.
        # Neutral at s=0 (×1.0), doubled at s=1 (×2.0).
        b_scale = 1f0 + pr.transformative_non_absorption
        a.backfire = clamp(b_scale * (p.backfire_R_weight * R_i +
                                      p.backfire_J_weight * J_i +
                                      p.backfire_V_weight * V_i +
                                      p.backfire_A_weight * A_global),
                           0f0, 1f0)
        a.backfire > 0.5f0 && (n_backfired += 1)

        # Defection: committed-and-isolated-and-scared → :defector.
        # Forkability (§19.5) — with probability `s`, the agent splits into
        # a fork instead of defecting (stays :contributor with reduced
        # commitment). Models the §19.5 "smooth divergence without schism".
        if is_committed(a) &&
           a.fear >= p.defect_fear_threshold &&
           protection < p.defect_protection_threshold
            if rand(world.rng) < pr.forkability
                a.commitment = clamp(a.commitment * 0.5f0, 0f0, 1f0)
                # stay :contributor — no role transition
            else
                a.role       = :defector
                a.commitment = 0f0
                n_defected  += 1
                for (_, replicas) in world.infra
                    delete!(replicas, a.id)
                end
            end
        end
    end

    push!(world.recent_backfire_count, n_backfired)
    if length(world.recent_backfire_count) > p.narrative_window
        popfirst!(world.recent_backfire_count)
    end
    return (n_defected, n_backfired)
end

# 4 — adoption / role advancement --------------------------------------------
# Implementation lives in `adoption.jl` (two rules: `:threshold` /
# `:complex_contagion`, selected by `params.adoption_rule`). The role state
# machine (`STAGE_SEQUENCE`, `_maybe_advance_role!`) lives here because it
# is used by both rules and by item 6's defection / infiltration logic.

const STAGE_SEQUENCE = (:outsider, :observer, :sympathizer, :user, :contributor, :steward)

function _maybe_advance_role!(a::Agent, drive::Float32, thresholds::NamedTuple;
                              principles::Principles=default_principles())
    rank = ROLE_RANK[a.role]
    # Negative-rank roles (:removed, :infiltrator, :defector) are not on the
    # advancement ladder. They get filtered upstream too, but the guard makes
    # the function safe to call on any agent.
    rank < 0 && return a
    # §19.7 Normative minimalism — thin constitution = lower stage thresholds.
    thr_scale = 1f0 - 0.5f0 * principles.normative_minimalism
    # §19.2 Functional autonomy — commitment boost per stage reached.
    fa_boost  = 1f0 + principles.functional_autonomy
    while rank < ROLE_RANK[:steward]
        next = STAGE_SEQUENCE[rank + 2]   # +2 because outsider=0, sympathizer=1, ...
        thr = getfield(thresholds, next) * thr_scale
        drive >= thr || break
        a.role = next
        stage_idx = Float32(ROLE_RANK[next] / ROLE_RANK[:steward])
        a.commitment = clamp(max(a.commitment, stage_idx) * fa_boost, 0f0, 1f0)
        rank += 1
    end
    return a
end

# 5 — infrastructure refill ---------------------------------------------------

"""
    refill_infrastructure!(world)

§12 dynamics (item 7). Each step:
  1. Decrement replica warm-up timers; remove zero entries.
  2. Drop any replica whose role is no longer committed (removed,
     defector, infiltrator) — they no longer hold the function.
  3. For each function whose *effective* replica count is below `m_k`,
     promote the most-committed unassigned agents and start their
     warm-up timer at `params.replica_latency`.
"""
function refill_infrastructure!(world::World)
    p  = world.params
    pr = get_principles(p)
    latency_default = get(p, :replica_latency, 1)
    # §19.6 Redundant stewardship — refill *target* is bumped above the
    # `infra_min` threshold used by Ψ_T. Stronger principle = more spare
    # replicas standing by, so the function's threshold stays met for
    # longer when the host arrests one.
    # DCS aggregation: if any cell senses infrastructure capture, lift the
    # redundancy bonus to that cell's principle strength. Mirrors the
    # global aggregation used for non_domination.
    rs_strength = max(pr.redundant_stewardship,
                      global_principle_max(world.dcs, :redundant_stewardship))
    redundancy_bonus = round(Int, 3 * rs_strength)

    # 1. tick down warm-ups.
    for (k, latmap) in world.infra_latency
        for (id, lat) in latmap
            lat - 1 <= 0 ? delete!(latmap, id) : (latmap[id] = lat - 1)
        end
    end

    # 2. drop replicas who are no longer committed (defected / infiltrated / removed).
    for (k, replicas) in world.infra
        for id in collect(replicas)
            is_committed(world.agents[id]) || begin
                delete!(replicas, id)
                haskey(world.infra_latency, k) && delete!(world.infra_latency[k], id)
            end
        end
    end

    # 3. refill where we're short.
    candidates = sort([a for a in world.agents if is_committed(a)];
                      by = a -> -a.commitment)
    for (k, replicas) in world.infra
        target = world.infra_min[k] + redundancy_bonus
        # Count only *effective* replicas (off warm-up) for the threshold check;
        # but the replicas Set still gets new entries up to `target` even when
        # some existing ones are still warming.
        effective = count(id -> id ∈ replicas &&
                                !haskey(get(world.infra_latency, k, Dict{Int,Int}()), id),
                          replicas)
        needed = target - effective
        needed <= 0 && continue
        for a in candidates
            length(replicas) >= target && break
            a.id ∈ replicas && continue
            push!(replicas, a.id)
            # Start the warm-up timer (no-op if latency_default == 0).
            if latency_default > 0
                latmap = get!(world.infra_latency, k, Dict{Int,Int}())
                latmap[a.id] = latency_default
            end
        end
    end
    return world
end

# 6 — legibility drift --------------------------------------------------------

"""
    drift_legibility!(world)

Skeleton §11 dynamics. Active agents that use privacy tools (proxied by
`commitment`) reduce their `L_ext`; everyone else accumulates legibility at
the global `leak_rate`. Stewards and contributors who curate the polis
(proxied by `commitment`) raise their `L_int` so internal trust can form.
Item 5 (host) will couple `L_ext` to surveillance pressure properly.
"""
function drift_legibility!(world::World)
    p  = world.params
    pr = get_principles(p)
    # §19.3 Bounded legibility — privacy tools scale the hide rate.
    hide_scale = 1f0 + 4f0 * pr.bounded_legibility
    # §19.9 Controlled permeability — opens the polis up, but more
    # permeable = more leakage to G_O. §17 A_H amplifies leakage
    # multiplicatively (host AI sees more).
    A_H = get(p, :A_H, 0f0)
    leak_scale = (1f0 + 0.5f0 * pr.controlled_permeability) *
                 (1f0 + 0.5f0 * A_H)
    hide_eff = p.hide_rate * hide_scale
    leak_eff = p.leak_rate * leak_scale
    @inbounds for a in world.agents
        a.role === :removed && continue
        if is_active(a)
            a.L_ext = clamp(a.L_ext + leak_eff - hide_eff * a.commitment, 0f0, 1f0)
            a.L_int = clamp(a.L_int + 0.05f0 * a.commitment, 0f0, 1f0)
        else
            a.L_ext = clamp(a.L_ext + leak_eff, 0f0, 1f0)
        end
    end
    return world
end

# --- E2 endogenous conflict (paper 2) --------------------------------------

const FACTION_LABELS = (:alpha, :beta, :gamma, :delta)

function _next_faction_label(rng, cap::Int, exclude::Symbol = :none)
    @assert cap <= length(FACTION_LABELS)
    pool = collect(FACTION_LABELS[1:cap])
    exclude !== :none && filter!(!=(exclude), pool)
    isempty(pool) && return FACTION_LABELS[1]
    return pool[1 + rand(rng, 0:length(pool)-1)]
end

"""
    _disagreement_step!(world, params)

E2 conflict mechanism. No-op when `params.disagreement_rate == 0`.

(1) Decay every agent's `disagreement` accumulator by
`params.disagreement_decay`.

(2) With probability `disagreement_rate` (per step), fire one
disagreement event: pick a random committed agent. If they have no
faction yet, assign them a fresh label. Otherwise, propagate the
faction to a random committed in-neighbour and raise both agents'
`disagreement` by 0.05.

(3) For any committed agent whose `disagreement` exceeds
`params.schism_threshold`, flip their faction to a different label
and trim a fraction of their incoming trust edges from
differently-factioned committed neighbours (probability
`params.trust_decay_on_fight` per edge).
"""
function _disagreement_step!(world::World, p::NamedTuple)
    rate = get(p, :disagreement_rate, 0.0)
    rate <= 0 && return world
    decay   = Float32(get(p, :disagreement_decay, 0.95))
    thresh  = Float32(get(p, :schism_threshold, 0.6f0))
    edge_p  = get(p, :trust_decay_on_fight, 0.10)
    cap     = get(p, :n_factions_cap, 4)
    g_S     = world.multiplex.layers[:S]

    # (1) Decay
    for a in world.agents
        a.disagreement *= decay
    end

    # (2) Fire k events this step, where k = floor(rate) plus a
    # Bernoulli draw on the fractional part (so the per-step mean
    # event count equals `rate`). Each event is a "small earthquake"
    # that bumps disagreement on the seed agent AND all committed
    # in-neighbours, propagating the seed's faction to unaligned ones.
    base_k = floor(Int, rate)
    frac = rate - base_k
    k_events = base_k + (rand(world.rng) < frac ? 1 : 0)
    if k_events > 0
        committed_ids = [a.id for a in world.agents
                         if ROLE_RANK[a.role] >= ROLE_RANK[:user]]
        if !isempty(committed_ids)
            for _ in 1:k_events
                i = world.agents[committed_ids[1 + rand(world.rng,
                                                        0:length(committed_ids)-1)]]
                i.disagreement += 0.20f0
                if i.faction === :none
                    i.faction = _next_faction_label(world.rng, cap)
                end
                for v in inneighbors(g_S, eltype(g_S)(i.id))
                    nb = world.agents[Int(v)]
                    ROLE_RANK[nb.role] >= ROLE_RANK[:user] || continue
                    nb.disagreement += 0.10f0
                    if nb.faction === :none
                        nb.faction = i.faction
                    end
                end
            end
        end
    end

    # (3) Forkability runs BEFORE schism so the schism step's
    # disagreement-reset doesn't pre-empt the fork trigger. If a
    # faction with at least `fork_min` committed members has any
    # member at-or-above the schism threshold, with probability
    # `forkability` trigger a clean fork: bulk-remove the cross-
    # faction edges for all faction members at once and reset their
    # disagreement. Substitutes for the messier per-agent schism
    # below, which would otherwise also fire.
    pr = get_principles(p)
    cap_thresh = get(p, :capture_centrality_thresh, 0.20)
    fork_min = get(p, :fork_min_cell_size, 5)
    forked_faction = :none
    if pr.forkability > 0
        faction_counts = Dict{Symbol, Int}()
        faction_max_d = Dict{Symbol, Float32}()
        for a in world.agents
            (a.faction === :none) && continue
            ROLE_RANK[a.role] >= ROLE_RANK[:user] || continue
            faction_counts[a.faction] = get(faction_counts, a.faction, 0) + 1
            faction_max_d[a.faction] = max(get(faction_max_d, a.faction, 0f0),
                                           a.disagreement)
        end
        for (f, count) in faction_counts
            count >= fork_min || continue
            get(faction_max_d, f, 0f0) >= thresh || continue
            if rand(world.rng) < pr.forkability
                forked_faction = f
                members = [a for a in world.agents
                           if a.faction == f &&
                              ROLE_RANK[a.role] >= ROLE_RANK[:user]]
                for a in members
                    vid = eltype(g_S)(a.id)
                    to_rem = eltype(g_S)[]
                    for v in inneighbors(g_S, vid)
                        nb = world.agents[Int(v)]
                        ROLE_RANK[nb.role] >= ROLE_RANK[:user] || continue
                        nb.faction === f && continue
                        push!(to_rem, v)
                    end
                    for v in to_rem; rem_edge!(g_S, v, vid); end
                    a.disagreement = 0f0
                end
                break    # at most one fork per step
            end
        end
    end

    # (4) Schisms: flip faction + trim trust edges for any over-threshold
    # agent NOT already in the forked faction. A small fraction of
    # schism-ing agents defect outright ("fed up, leave"); forkability
    # reduces this defection probability since the agent had a clean
    # path out via fork.
    defect_prob = Float32(get(p, :factional_defection_prob, 0.15)) *
                  (1f0 - pr.forkability)
    for a in world.agents
        ROLE_RANK[a.role] >= ROLE_RANK[:user] || continue
        a.faction === forked_faction && continue
        a.disagreement >= thresh || continue
        old_faction = a.faction
        a.faction = _next_faction_label(world.rng, cap, old_faction)
        a.disagreement = 0f0
        vid = eltype(g_S)(a.id)
        to_remove = eltype(g_S)[]
        for v in inneighbors(g_S, vid)
            nb = world.agents[Int(v)]
            ROLE_RANK[nb.role] >= ROLE_RANK[:user] || continue
            nb.faction === a.faction && continue
            if rand(world.rng) < edge_p
                push!(to_remove, v)
            end
        end
        for v in to_remove
            rem_edge!(g_S, v, vid)
        end
        # Endogenous-conflict defection: schism-ing agents may exit
        # the polis. Forkability suppresses this.
        if defect_prob > 0 && rand(world.rng) < defect_prob
            a.role = :defector
            # Drop them from any infrastructure they held.
            for (_, replicas) in world.infra
                delete!(replicas, a.id)
            end
        end
    end

    # (5) Reversible governance: a steward is "captured" when their
    # active-committed-neighbour count exceeds `cap_thresh × n_committed`
    # (default 20% of the committed population trusts this single
    # steward). With probability `reversible_governance` demote them
    # to contributor (recovery from panic-centralisation).
    if pr.reversible_governance > 0
        n_committed = count(a -> ROLE_RANK[a.role] >= ROLE_RANK[:user],
                            world.agents)
        if n_committed > 0
            capture_n = ceil(Int, cap_thresh * n_committed)
            stewards = [a for a in world.agents if a.role === :steward]
            for s in stewards
                # Count committed in-neighbours (trust ties pointing AT s).
                n_followers = 0
                for v in inneighbors(g_S, eltype(g_S)(s.id))
                    nb = world.agents[Int(v)]
                    if ROLE_RANK[nb.role] >= ROLE_RANK[:user]
                        n_followers += 1
                    end
                end
                if n_followers >= capture_n &&
                   rand(world.rng) < pr.reversible_governance
                    s.role = :contributor
                    s.disagreement = 0f0
                end
            end
        end
    end

    return world
end

# Expose for test/diagnostic scripts.
export _disagreement_step!
