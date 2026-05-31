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
    Fingerprints

Neural diffusion fingerprints (NDF) — APPNP-style generalization of the
linear "Diffusion Fingerprints" construction (Dubuisson, Eckmann, Agazzi,
arXiv:1408.4966v2). The original DF is recovered as the linear, no-features,
K→∞ limit of NDF.
"""
module Fingerprints

using BSON
using JLD2
using ChainRulesCore
using Flux
using LinearAlgebra
using NPZ
using SparseArrays: SparseMatrixCSC, sparse, spdiagm
using Statistics: mean
using Random
using Random: MersenneTwister
using LightGraphs: AbstractGraph, nv, edges, src, dst, is_directed,
                   vertices, outdegree, indegree

using SparseArrays: findnz

include("ndf.jl")
include("training.jl")
include("text.jl")
include("text_gcn.jl")

export NDF, NDFEncoder, propagate, prepare_adjacency,
       prepare_adjacency_directed_ppr, default_node_features,
       top_central_nodes,
       Document, TrainConfig, train!, fingerprints, load_ndf_state,
       BlogDoc, read_blog_corpus, build_vocab,
       compute_ppmi, compute_v1_collocation,
       build_domain_graph, blogdocs_to_documents,
       porter_stem_lite,
       TextGcnDoc, read_text_gcn_corpus,
       textgcn_to_documents, build_text_gcn_vocab,
       load_bert_word_features, read_blog_posts

end # module
