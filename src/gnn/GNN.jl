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

