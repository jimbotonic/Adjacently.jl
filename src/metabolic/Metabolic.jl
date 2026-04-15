#
# Adjacently: Julia Complex Directed Networks Library
# Copyright (C) 2016-2026 Jimmy Dubuisson <jimmy@dubuisson.ch>
#
# This program is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.
#

#
# Adjacently.Metabolic — metabolic pathway discovery infrastructure.
#
# - BiGG JSON model loading → directed metabolite graph
# - Flux Balance Analysis (JuMP + HiGHS) with 5 condition presets
# - Chemistry-aware edge weighting (atom similarity from BiGG formulas)
# - MetaRoute-style reaction-context ranking
# - RDT atom-tracing edge-weight loader
#
# Primary entry points:
#   build_metabolic_graph(model_path; compartment, keep_currency)
#   run_fba_condition(model_path, condition; biomass_id)
#   flux_to_edge_weights(flux, model_path, id_to_idx; transform)
#   build_atom_weights(model_path, g, id_to_idx; transform)
#   build_metaroute_weights(model_path, g, id_to_idx)
#   load_rdt_weights(tsv_path, id_to_idx)
#

module Metabolic

using JSON
using JuMP
using HiGHS
using SparseArrays
using Random
using LightGraphs

include("atom_weights.jl")    # parse_formula, atom_similarity, build_atom_weights, atom_weight_fn
include("graph.jl")           # build_metabolic_graph, load_bigg_metabolic, CURRENCY_METABOLITES
include("fba.jl")             # FBA_CONDITIONS, run_fba, run_fba_condition, run_fba_sampled, flux_to_edge_weights
include("metaroute.jl")       # build_metaroute_weights (depends on atom_weights functions)
include("rdt.jl")             # load_rdt_weights

export
    # Graph loading
    build_metabolic_graph, load_bigg_metabolic, save_metabolic_graph,
    CURRENCY_METABOLITES,

    # FBA
    FBA_CONDITIONS, load_fba_model, apply_condition, solve_fba,
    run_fba, run_fba_condition, run_fba_sampled, flux_to_edge_weights,

    # Atom-similarity weights
    ATOM_ELEMENTS, parse_formula, atom_similarity, heavy_atom_count,
    build_atom_weights, atom_weight_fn,

    # MetaRoute-style
    build_metaroute_weights,

    # RDT atom-tracing
    load_rdt_weights

end # module Metabolic
