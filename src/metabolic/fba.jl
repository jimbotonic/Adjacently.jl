#
# Flux Balance Analysis (FBA) utilities.
#
# Solves: max v_biomass  s.t. S·v = 0, lb ≤ v ≤ ub.
# Used for producing edge weights that bias path-finding toward flux-active
# regions of the metabolic graph.
#

# ── FBA condition presets: exchange-reaction bound overrides ──
const FBA_CONDITIONS = Dict{Symbol, Vector{Tuple{String,Float64,Float64}}}(
    :aerobic_glucose   => Tuple{String,Float64,Float64}[],  # default M9 + glucose + O2

    :anaerobic_glucose => [
        ("EX_o2_e", 0.0, 1000.0),
    ],

    :aerobic_acetate   => [
        ("EX_glc__D_e", 0.0, 1000.0),
        ("EX_ac_e",    -10.0, 1000.0),
    ],

    :aerobic_glycerol  => [
        ("EX_glc__D_e", 0.0, 1000.0),
        ("EX_glyc_e",  -10.0, 1000.0),
    ],

    :nitrogen_limited  => [
        ("EX_nh4_e",    -1.0, 1000.0),
    ],
)

"""
    load_fba_model(model_path)

Parse a BiGG JSON into (`S`, `lbs`, `ubs`, `rxn_ids`, `met_ids`, `rxn_idx`,
`met_idx`). `S` is the sparse stoichiometric matrix.
"""
function load_fba_model(model_path::String)
    d = JSON.parsefile(model_path)

    all_mets = [m["id"] for m in d["metabolites"]]
    met_idx = Dict(m => i for (i, m) in enumerate(all_mets))
    n_mets = length(all_mets)

    all_rxns = [r["id"] for r in d["reactions"]]
    rxn_idx = Dict(r => i for (i, r) in enumerate(all_rxns))
    n_rxns = length(all_rxns)

    I_s, J_s, V_s = Int[], Int[], Float64[]
    lbs = Float64[]
    ubs = Float64[]
    for (j, r) in enumerate(d["reactions"])
        push!(lbs, r["lower_bound"])
        push!(ubs, r["upper_bound"])
        for (met_id, coeff) in r["metabolites"]
            i = met_idx[met_id]
            push!(I_s, i); push!(J_s, j); push!(V_s, coeff)
        end
    end
    S = sparse(I_s, J_s, V_s, n_mets, n_rxns)

    return (S=S, lbs=lbs, ubs=ubs,
            rxn_ids=all_rxns, met_ids=all_mets,
            rxn_idx=rxn_idx, met_idx=met_idx)
end

"""
    apply_condition(model, condition::Symbol) -> (lbs, ubs)

Return a copy of the model's bounds with `FBA_CONDITIONS[condition]`
overrides applied. The base model is not mutated.
"""
function apply_condition(model, condition::Symbol)
    haskey(FBA_CONDITIONS, condition) || error("Unknown FBA condition: $condition")
    lbs = copy(model.lbs)
    ubs = copy(model.ubs)
    for (rxn_id, lb, ub) in FBA_CONDITIONS[condition]
        haskey(model.rxn_idx, rxn_id) || (@warn "Exchange $rxn_id not found — skipping"; continue)
        j = model.rxn_idx[rxn_id]
        lbs[j] = lb
        ubs[j] = ub
    end
    return (lbs=lbs, ubs=ubs)
end

"""
    solve_fba(model, lbs, ubs; biomass_id)

Maximize biomass flux subject to `S·v = 0` and `lbs ≤ v ≤ ubs`.
Returns `flux::Vector{Float64}` (all zeros if infeasible).
"""
function solve_fba(model, lbs, ubs;
                   biomass_id::String="BIOMASS_Ec_iML1515_core_75p37M")
    n_rxns = length(lbs)
    bio_idx = model.rxn_idx[biomass_id]

    jm = Model(HiGHS.Optimizer); set_silent(jm)
    @variable(jm, lbs[j] <= v[j=1:n_rxns] <= ubs[j])
    @constraint(jm, model.S * v .== 0)
    @objective(jm, Max, v[bio_idx])
    optimize!(jm)

    status = termination_status(jm)
    if status != OPTIMAL
        @warn "FBA infeasible: $status"
        return zeros(n_rxns)
    end
    return value.(v)
end

"""
    run_fba_condition(model_path, condition; biomass_id) -> (flux, rxn_ids, met_ids)

Load the model, apply condition-specific exchange bounds, solve FBA.
"""
function run_fba_condition(model_path::String, condition::Symbol;
                           biomass_id::String="BIOMASS_Ec_iML1515_core_75p37M")
    m = load_fba_model(model_path)
    @info "FBA[$condition]: $(length(m.met_ids)) metabolites, $(length(m.rxn_ids)) reactions"
    c = apply_condition(m, condition)
    flux = solve_fba(m, c.lbs, c.ubs; biomass_id=biomass_id)
    active = count(x -> abs(x) > 1e-9, flux)
    bio_flux = flux[m.rxn_idx[biomass_id]]
    @info "FBA[$condition]: biomass=$(round(bio_flux, digits=4)), active=$active/$(length(flux))"
    return (flux, m.rxn_ids, m.met_ids)
end

"""
    run_fba(model_path; biomass_id)

Aerobic-glucose FBA. Thin wrapper over `run_fba_condition(..., :aerobic_glucose; ...)`.
"""
function run_fba(model_path::String;
                 biomass_id::String="BIOMASS_Ec_iML1515_core_75p37M")
    return run_fba_condition(model_path, :aerobic_glucose; biomass_id=biomass_id)
end

"""
    run_fba_sampled(model_path, n_samples; condition, biomass_id, perturb, seed)

Sample multiple near-optimal FBA solutions by random objective perturbations.
Returns `(mean_flux, rxn_ids, met_ids)` where `mean_flux[j] = mean(|flux_j|)`
across the N feasible samples. Only samples with biomass ≥ 0.95·optimum are
kept. See scripts/test_sampled_flux.jl for the (null) experimental result —
sampled flux under-performs single-FBA on our benchmark.
"""
function run_fba_sampled(model_path::String, n_samples::Int;
                         condition::Symbol=:aerobic_glucose,
                         biomass_id::String="BIOMASS_Ec_iML1515_core_75p37M",
                         perturb::Float64=1e-2,
                         seed::Int=42)
    rng = Random.MersenneTwister(seed)

    m = load_fba_model(model_path)
    c = apply_condition(m, condition)
    n_rxns = length(m.rxn_ids)
    bio_idx = m.rxn_idx[biomass_id]

    bio_opt_flux = solve_fba(m, c.lbs, c.ubs; biomass_id=biomass_id)
    bio_opt = bio_opt_flux[bio_idx]
    if bio_opt <= 0
        @warn "Biomass infeasible or zero — returning zero flux"
        return (zeros(n_rxns), m.rxn_ids, m.met_ids)
    end

    lbs = copy(c.lbs); ubs = copy(c.ubs)
    lbs[bio_idx] = 0.95 * bio_opt

    flux_sum = zeros(n_rxns)
    n_kept = 0
    for i in 1:n_samples
        coeffs = perturb * (2 * rand(rng, n_rxns) .- 1)
        coeffs[bio_idx] = 1.0
        jm = Model(HiGHS.Optimizer); set_silent(jm)
        @variable(jm, lbs[j] <= v[j=1:n_rxns] <= ubs[j])
        @constraint(jm, m.S * v .== 0)
        @objective(jm, Max, sum(coeffs[j] * v[j] for j in 1:n_rxns))
        optimize!(jm)
        status = termination_status(jm)
        if status == OPTIMAL
            fluxes = value.(v)
            flux_sum .+= abs.(fluxes)
            n_kept += 1
        end
    end
    n_kept == 0 && (@warn "No feasible samples"; return (zeros(n_rxns), m.rxn_ids, m.met_ids))
    mean_flux = flux_sum ./ n_kept
    active = count(x -> abs(x) > 1e-9, mean_flux)
    @info "Sampled FBA ($condition, N=$n_kept/$n_samples): mean-|flux| active = $active / $n_rxns"
    return (mean_flux, m.rxn_ids, m.met_ids)
end

"""
    flux_to_edge_weights(flux, model_path, id_to_idx;
                         compartment="c", currency=CURRENCY_METABOLITES, transform=:log)

Convert an FBA flux vector into per-edge weights. `transform` ∈ {:raw, :log,
:binary}. Returns `Dict{Tuple{UInt32,UInt32}, Float64}` indexed by the graph's
vertex IDs (from `build_metabolic_graph`'s `id_to_idx`).
"""
function flux_to_edge_weights(flux::Vector{Float64},
                               model_path::String,
                               id_to_idx::Dict{String,Int};
                               compartment::String="c",
                               currency::Set{String}=CURRENCY_METABOLITES,
                               transform::Symbol=:log)
    d = JSON.parsefile(model_path)
    weights = Dict{Tuple{UInt32,UInt32}, Float64}()

    for (j, r) in enumerate(d["reactions"])
        f_raw = abs(flux[j])
        f_raw < 1e-9 && continue
        f = if transform == :log
            log(1 + f_raw)
        elseif transform == :binary
            1.0
        else
            f_raw
        end

        mets = r["metabolites"]
        lb = get(r, "lower_bound", 0.0)

        substrates = UInt32[]
        products = UInt32[]
        for (met_id, coeff) in mets
            endswith(met_id, "_$compartment") || continue
            base_id = replace(met_id, r"_[a-z]$" => "")
            base_id in currency && continue
            haskey(id_to_idx, met_id) || continue

            idx = UInt32(id_to_idx[met_id])
            if coeff < 0
                push!(substrates, idx)
            elseif coeff > 0
                push!(products, idx)
            end
        end

        isempty(substrates) && continue
        isempty(products) && continue

        for s in substrates, p in products
            s == p && continue
            key = (s, p)
            weights[key] = get(weights, key, 0.0) + f
        end

        if lb < 0 && flux[j] < -1e-9
            for p in products, s in substrates
                s == p && continue
                key = (p, s)
                weights[key] = get(weights, key, 0.0) + abs(flux[j])
            end
        end
    end

    @info "Edge weights: $(length(weights)) weighted edges"
    return weights
end
