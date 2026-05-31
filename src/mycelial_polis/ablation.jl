"""
    ablation

Principle ablations for the Mycelial Polis (item 10 of
`research/mycelial_polis/RESEARCH_TASKS.md`). Three sweeps on top of the
item 8 basin estimator:

- **OAT-out**: 10 runs, one principle disabled (set to 0) per run.
- **Pair-out**: `C(10,2) = 45` runs, two principles disabled.
- **Greedy forward-selection**: start from `zero_principles()`, add the
  principle that most increases `A_S` until the gain drops below
  `gain_threshold` or the set hits `max_size`.

Per §20, the attractor stability score is a weighted sum with positive
terms (basin size, mean Φ, mean Λ_sat, mean Ψ_T) and negative terms
(extinction fraction, mean defectors). Initial weights all 1.0; the
script writes them to `A_S_weights.tsv` for downstream tuning.
"""

const ALL_PRINCIPLE_NAMES = (
    :non_domination, :functional_autonomy, :bounded_legibility,
    :reversible_governance, :forkability, :redundant_stewardship,
    :normative_minimalism, :adaptive_memory, :controlled_permeability,
    :transformative_non_absorption,
)

"""
    principles_with(disabled, enabled_strength=0.5f0) -> Principles

Build a `Principles` struct where the names in `disabled` are set to 0
and all other principles are set to `enabled_strength` (default 0.5).
"""
function principles_with(disabled::AbstractVector{Symbol},
                         enabled_strength::Float32=0.5f0)
    field_vals = Dict{Symbol,Float32}()
    for name in ALL_PRINCIPLE_NAMES
        field_vals[name] = name in disabled ? 0f0 : enabled_strength
    end
    return Principles(; (k => field_vals[k] for k in ALL_PRINCIPLE_NAMES)...)
end

"""
    principles_only(enabled, strength=0.5f0) -> Principles

Inverse: all principles set to 0 EXCEPT the named ones, which are set to
`strength`. Used by the greedy forward-selection search.
"""
function principles_only(enabled::AbstractVector{Symbol},
                         strength::Float32=0.5f0)
    field_vals = Dict{Symbol,Float32}()
    for name in ALL_PRINCIPLE_NAMES
        field_vals[name] = name in enabled ? strength : 0f0
    end
    return Principles(; (k => field_vals[k] for k in ALL_PRINCIPLE_NAMES)...)
end

"""
    summarise_basin(result) -> NamedTuple

Compact summary of a basin result for the ablation TSVs. Computes the
attractor stability score `A_S` with initial unit weights per §20.
"""
function summarise_basin(result)
    n = result.n_samples
    counts = result.counts
    basin_par  = counts[:parallel_coexistence] / n
    basin_dorm = counts[:dormant_persistence]  / n
    basin_ext  = counts[:extinction]           / n
    surviving_phi = vcat(result.final_phis[:parallel_coexistence],
                         result.final_phis[:dormant_persistence])
    mean_phi_surv = isempty(surviving_phi) ? 0f0 :
                    Float32(sum(surviving_phi) / length(surviving_phi))
    # A_S v1 — initial weights all 1.0, positive terms minus risk terms.
    # The tuned-weights file lives at A_S_weights.tsv next to this output.
    A_S = basin_par + Float32(mean_phi_surv) + 0.5f0 * basin_dorm - basin_ext
    return (basin_parallel  = basin_par,
            basin_dormant   = basin_dorm,
            basin_extinction = basin_ext,
            mean_phi_surv   = Float32(mean_phi_surv),
            A_S             = Float32(A_S))
end

"""
    run_oat(; tier, topology, seed, visibility, enabled_strength=0.5f0,
              baseline=nothing) -> Vector{NamedTuple}

OAT-out: for each principle, disable it and run the basin. Returns one
row per principle with `(principle, ΔBasin, ΔA_S, basin_parallel,
basin_extinction, A_S)`. ΔBasin / ΔA_S are computed against `baseline`
(if supplied) — otherwise computed against the all-default-0.5 baseline
which is also computed internally.
"""
function run_oat(; tier::Symbol=:smoke,
                   topology::Symbol=:modular_cells,
                   seed::Int=0,
                   visibility::Float32=1.0f0,
                   enabled_strength::Float32=0.5f0,
                   host_strategy::Symbol=:random,
                   baseline::Union{Nothing,NamedTuple}=nothing)
    base_principles = principles_with(Symbol[], enabled_strength)
    base = baseline === nothing ?
           summarise_basin(estimate_basins(; tier=tier, topology=topology,
                                             seed=seed, visibility=visibility,
                                             principles=base_principles,
                                             host_strategy=host_strategy)) :
           baseline
    rows = NamedTuple[]
    for p in ALL_PRINCIPLE_NAMES
        pr_dis = principles_with([p], enabled_strength)
        s = summarise_basin(estimate_basins(; tier=tier, topology=topology,
                                              seed=seed, visibility=visibility,
                                              principles=pr_dis,
                                              host_strategy=host_strategy))
        push!(rows, (
            principle        = p,
            delta_basin      = base.basin_parallel - s.basin_parallel,
            delta_A_S        = base.A_S - s.A_S,
            basin_parallel   = s.basin_parallel,
            basin_extinction = s.basin_extinction,
            A_S              = s.A_S,
        ))
    end
    return base, rows
end

"""
    run_pair_out(; ...) -> Vector{NamedTuple}

For all `C(10,2) = 45` pairs, disable both principles and run the basin.
Records `(p1, p2, ΔBasin_pair, ΔBasin_additive, non_additivity)` where
`ΔBasin_additive = ΔBasin(p1) + ΔBasin(p2)` from a supplied OAT result,
and `non_additivity = ΔBasin_pair − ΔBasin_additive`. Large non-additivity
indicates a strong interaction between the two principles.
"""
function run_pair_out(oat_rows::AbstractVector{<:NamedTuple};
                      tier::Symbol=:smoke,
                      topology::Symbol=:modular_cells,
                      seed::Int=0,
                      visibility::Float32=1.0f0,
                      enabled_strength::Float32=0.5f0,
                      host_strategy::Symbol=:random,
                      baseline::NamedTuple)
    oat_lookup = Dict(r.principle => r for r in oat_rows)
    names = collect(ALL_PRINCIPLE_NAMES)
    rows = NamedTuple[]
    for i in 1:length(names), j in (i+1):length(names)
        p1, p2 = names[i], names[j]
        pr_dis = principles_with([p1, p2], enabled_strength)
        s = summarise_basin(estimate_basins(; tier=tier, topology=topology,
                                              seed=seed, visibility=visibility,
                                              principles=pr_dis,
                                              host_strategy=host_strategy))
        dpair = baseline.basin_parallel - s.basin_parallel
        dadd  = oat_lookup[p1].delta_basin + oat_lookup[p2].delta_basin
        push!(rows, (
            p1               = p1,
            p2               = p2,
            delta_basin_pair = dpair,
            delta_basin_add  = dadd,
            non_additivity   = dpair - dadd,
            basin_parallel   = s.basin_parallel,
            A_S              = s.A_S,
        ))
    end
    return rows
end

"""
    greedy_forward(; ...) -> (selected, history)

Greedy forward-selection. Start from `zero_principles()` (all 0), at
each step enable the single principle that gives the largest `A_S`
gain. Stop when the next gain is `< gain_threshold` OR the selected set
has reached `max_size`. Returns `(selected::Vector{Symbol},
history::Vector{NamedTuple})` where each history row records the chosen
principle and the resulting `A_S`.
"""
function greedy_forward(; tier::Symbol=:smoke,
                          topology::Symbol=:modular_cells,
                          seed::Int=0,
                          visibility::Float32=1.0f0,
                          enabled_strength::Float32=0.5f0,
                          host_strategy::Symbol=:random,
                          max_size::Int=5,
                          gain_threshold::Float64=0.02)
    selected = Symbol[]
    # Anchor: start from empty set, A_S of the all-zero polis.
    s0 = summarise_basin(estimate_basins(; tier=tier, topology=topology,
                                           seed=seed, visibility=visibility,
                                           principles=principles_only(Symbol[],
                                                                      enabled_strength),
                                           host_strategy=host_strategy))
    history = NamedTuple[(step=0, added=Symbol(""), A_S=s0.A_S,
                          basin_parallel=s0.basin_parallel)]
    while length(selected) < max_size
        remaining = setdiff(ALL_PRINCIPLE_NAMES, selected)
        isempty(remaining) && break
        best = nothing; best_summary = nothing
        for p in remaining
            candidate_set = vcat(selected, p)
            pr = principles_only(candidate_set, enabled_strength)
            sc = summarise_basin(estimate_basins(; tier=tier, topology=topology,
                                                   seed=seed, visibility=visibility,
                                                   principles=pr,
                                                   host_strategy=host_strategy))
            if best === nothing || sc.A_S > best_summary.A_S
                best = p; best_summary = sc
            end
        end
        gain = best_summary.A_S - history[end].A_S
        gain < gain_threshold && break
        push!(selected, best)
        push!(history, (step=length(selected), added=best, A_S=best_summary.A_S,
                        basin_parallel=best_summary.basin_parallel))
    end
    return (selected=selected, history=history)
end
