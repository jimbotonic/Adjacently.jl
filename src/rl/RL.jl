#
# Adjacently: Julia Complex Directed Networks Library
# Copyright (C) 2016-2025 Jimmy Dubuisson <jimmy.dubuisson@gmail.com>
#
# This program is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.
#
# This program is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.
#

"""
    RL

Reinforcement learning framework for learning compression encoding policies.
Uses tabular Q-learning to find per-vertex encoding decisions that minimize
bits per edge on a given graph.

The action space covers:
- Integer encoding: Fibonacci, Zeta, Elias gamma, Elias delta
- Reference mode: none, reference, recursive reference
- Minimum interval length: 2, 3, 4, 5

State features are extracted per-vertex:
- Degree bucket
- Interval density (fraction of neighbors in consecutive intervals)
- Maximum delta gap
- Best reference overlap ratio
- Reference window density
"""
module RL

using Random
using ..Compression: compress_intervals, delta_encode_vector, estimate_encoded_value_cost, estimate_block_encoding_cost
using ..GNN: gnn_backward_rl!, gnn_forward

include("features.jl")
include("environment.jl")
include("policy.jl")
include("agent.jl")
include("actor_critic.jl")
include("gnn_head.jl")

export VertexFeatures, extract_features, feature_index,
       Action, action_to_index, action_from_index,
       CompressionEnv, reset!, step!, get_bits_per_edge, encode_vertex,
       QPolicy, select_action, best_action, update!, decay_epsilon!,
       get_policy_summary, save_policy, load_policy,
       TrainingConfig, TrainingResult, train!, evaluate, run_baseline,
       NUM_STATES, NUM_ACTIONS,
       ACPolicy, select_action_ac, best_action_ac, train_actor_critic!, evaluate_actor_critic,
       GnnACPolicy, select_action_gnn, best_action_gnn, train_gnn_actor_critic!, evaluate_gnn_actor_critic,
       train_gnn_ac_e2e!, save_gnn_ac_policy, load_gnn_ac_policy

end # module RL
