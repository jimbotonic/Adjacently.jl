#
# Adjacently: Julia Complex Directed Networks Library
# Copyright (C) 2016-2026 Jimmy Dubuisson <jimmy@dubuisson.ch>
#

"""
    MycelialPolis

Skeleton agent-based simulation of non-hierarchical latent virtual societies
in reactive host environments. Implements roadmap §5–§14 at a level
sufficient for the smoke test of item 2 in `research/mycelial_polis/RESEARCH_TASKS.md`.

Filled-in components vs roadmap (deliberately minimal):

| Subsystem              | Status        | Owning TODO item |
|------------------------|---------------|------------------|
| Agent state (§6)       | full          | item 2           |
| Multiplex (§7)         | scaffold      | item 3 (per-layer builders) |
| Topologies (§7)        | one real, two stubs | item 3      |
| Adoption (§8)          | threshold only       | item 4      |
| Host adversary (§10)   | random only          | item 5      |
| Repression (§9)        | skeleton equations   | item 6      |
| Infrastructure (§12)   | simple refill        | item 7      |
| Metrics (§14)          | Φ/Ψ_T/Λ/Γ skeletons  | item 7      |
| Basin estimation (§15) | n/a                  | item 8      |
| Principles (§19)       | n/a                  | item 9      |
"""
module MycelialPolis

using Random
using Random: MersenneTwister, AbstractRNG, randperm
using Statistics: mean, std
using LightGraphs: AbstractGraph, SimpleDiGraph, nv, ne, edges, src, dst,
                   add_edge!, rem_edge!, outneighbors, inneighbors,
                   outdegree, indegree, betweenness_centrality, core_number
import LightGraphs

import ..Adjacently   # used by build_world to reach Adjacently.Generators

include("agent.jl")
include("principles.jl")
include("multiplex.jl")
include("topologies.jl")
include("host.jl")
include("metrics.jl")
include("dynamics.jl")
include("adoption.jl")
include("basin.jl")
include("ablation.jl")
include("dcs.jl")              # Distributed constitutional sensing (paper 2)

export Agent, ROLE_RANK, is_active, is_committed, is_steward,
       is_infiltrator, is_defector,
       Multiplex, World, build_world, default_params, layer,
       build_topology, topology_summary, natural_partition,
       HostStrategy, RandomHost, DegreeHost, BetweennessHost,
       LegibilityHost, LocalizedHost, InfiltrationFirstHost, AttritionHost,
       ATTACK_COST, act!,
       step!, snapshot,
       phi, psi_T, lambda, lambda_sat, gamma,
       hierarchy_score, hierarchy_components,
       committed_frac, giant_active_frac, kcore_depth, modularity_active,
       infra_outages, default_phi_weights,
       active_count, committed_count,
       Principles, default_principles, zero_principles, get_principles, principle_strength,
       apply_adoption!, threshold_adopt!, complex_contagion_adopt!,
       lhs, classify, classify_v2, time_to_attractor, SampleParams, map_lhs,
       build_sample_params, run_sample, estimate_basins, BASIN_TIERS,
       ALL_PRINCIPLE_NAMES, principles_with, principles_only,
       summarise_basin, run_oat, run_pair_out, greedy_forward,
       # DCS (paper 2)
       CellSensors, CellPrinciples, DCSState, build_dcs_state,
       cell_principles_for, global_principle_max

end # module MycelialPolis
