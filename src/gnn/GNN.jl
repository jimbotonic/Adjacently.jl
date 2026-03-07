#
# Adjacently: Julia Complex Directed Networks Library
# Copyright (C) 2016-2026 Jimmy Dubuisson <jimmy@dubuisson.ch>
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
    GNN

Graph Neural Network utilities for learning compression-aware signals.
Currently provides a lightweight 2-layer GNN that produces per-vertex scores
and supports proxy + REINFORCE training. This module replaces the former GCN
module and removes vertex renumbering/ordering APIs.

No external ML dependencies — manual gradients through a 2-layer architecture.
"""
module GNN

using LightGraphs: AbstractGraph, nv, ne, vertices, edges, outneighbors, inneighbors, is_directed, src, dst, add_edge!, add_vertices!
using SparseArrays: sparse, SparseMatrixCSC
using LinearAlgebra: I, norm, dot
using Random: MersenneTwister, randn
using Statistics: var, mean
using Logging

using ..CustomLightGraphs: SimpleDiGraph, SimpleGraph
using ..Graph: get_sparse_P_matrix
using ..PageRank: PR
using ..CompressionUtils: compress_intervals

include("gat.jl")
include("model.jl")
include("training.jl")
include("relabeling.jl")

export GNNModel, gnn_forward, gnn_backward!, gnn_backward_rl!,
       save_gnn_model, load_gnn_weights!,
       TrainConfig, train_gnn_proxy!, train_gnn_reinforce!,
       relabel_vertices_gnn

end # module GNN

