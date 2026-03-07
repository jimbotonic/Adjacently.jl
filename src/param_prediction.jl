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
    ParamPrediction

Supervised GNN-based CGE parameter predictor. A GNN produces per-vertex
embeddings which are mean-pooled into a graph-level feature vector. Seven
independent linear heads then predict parameter configurations via softmax.
Training uses cross-entropy loss against grid-search-optimal labels.

Master GNN weights are shared across varying-size training graphs via a
copy-in / copy-out pattern.
"""
module ParamPrediction

using LightGraphs: AbstractGraph, nv, ne, vertices, outneighbors, is_directed
using SparseArrays: sparse, SparseMatrixCSC
using LinearAlgebra: norm, dot
using Random: MersenneTwister, rand, randn, randperm
using Statistics: mean
using Logging
using BSON: @save, @load

using ..CustomLightGraphs: SimpleDiGraph
using ..GNN: GNNModel, gnn_forward, gnn_backward_rl!,
             save_gnn_model, load_gnn_weights!
using ..IO: BitWriter, flush_bitwriter
using ..Compression: CGE

export ParamPredictor, MasterGNNWeights,
       graph_level_readout, predict_params, actions_to_cge_params,
       compute_bpe, grid_search_best_params, generate_training_labels,
       supervised_loss_backward!, train_param_predictor!,
       evaluate_param_predictor, save_param_predictor, load_param_predictor

# ─────────────────────────────────────────────────────────────────────────────
# Structs
# ─────────────────────────────────────────────────────────────────────────────

"""
    ParamPredictor

Multi-head supervised classifier over 7 CGE parameter groups.
Each head is a linear layer from graph features to a softmax over options,
trained with cross-entropy loss.
"""
mutable struct ParamPredictor
    heads::Vector{Matrix{Float64}}    # 7 heads: each (feat_dim × n_options_i)
    biases::Vector{Vector{Float64}}   # 7 bias vectors
    head_names::Vector{Symbol}
    head_options::Vector{Vector}      # concrete option values per head
    lr::Float64
    feat_dim::Int
end

const PARAM_HEADS = [
    (:window,        [8, 16, 32, 64]),
    (:intervals,     [true, false]),
    (:lr_split,      [true, false]),
    (:mil,           [2, 3, 4, 5]),
    (:copy_adaptive, [true, false]),
    (:stop_deltas,   [true, false]),
    (:zigzag,        [true, false]),
]

function ParamPredictor(feat_dim::Int; lr::Float64=0.01, seed::Int=42)
    rng = MersenneTwister(seed)
    heads = Matrix{Float64}[]
    biases = Vector{Float64}[]
    names = Symbol[]
    options = Vector[]

    for (name, opts) in PARAM_HEADS
        n_opts = length(opts)
        W = randn(rng, feat_dim, n_opts) .* 0.01
        b = zeros(Float64, n_opts)
        push!(heads, W)
        push!(biases, b)
        push!(names, name)
        push!(options, opts)
    end

    return ParamPredictor(heads, biases, names, options, lr, feat_dim)
end

"""
    MasterGNNWeights

Holds master copies of GNN weights that are shared across varying-size graphs.
Before each forward pass, master weights are copied into a per-graph GNNModel.
After backward, gradients are accumulated back to master weights.
"""
mutable struct MasterGNNWeights
    W1::Matrix{Float64}
    W2::Matrix{Float64}
    input_dim::Int
    hidden_dim::Int
    use_gat::Bool
    gat_W1::Matrix{Float64}
    gat_a1::Matrix{Float64}
    gat_W2::Matrix{Float64}
    gat_a2::Matrix{Float64}
end

function MasterGNNWeights(gnn::GNNModel)
    if gnn.use_gat
        return MasterGNNWeights(
            Matrix{Float64}(undef, 0, 0), Matrix{Float64}(undef, 0, 0),
            gnn.input_dim, gnn.hidden_dim, true,
            copy(gnn.gat_layers[1].W), copy(gnn.gat_layers[1].a),
            copy(gnn.gat_layers[2].W), copy(gnn.gat_layers[2].a)
        )
    else
        return MasterGNNWeights(
            copy(gnn.W1), copy(gnn.W2),
            gnn.input_dim, gnn.hidden_dim, false,
            Matrix{Float64}(undef, 0, 0), Matrix{Float64}(undef, 0, 0),
            Matrix{Float64}(undef, 0, 0), Matrix{Float64}(undef, 0, 0)
        )
    end
end

function copy_master_to_gnn!(gnn::GNNModel, master::MasterGNNWeights)
    if master.use_gat
        gnn.gat_layers[1].W .= master.gat_W1
        gnn.gat_layers[1].a .= master.gat_a1
        gnn.gat_layers[2].W .= master.gat_W2
        gnn.gat_layers[2].a .= master.gat_a2
    else
        gnn.W1 .= master.W1
        gnn.W2 .= master.W2
    end
end

function copy_gnn_to_master!(master::MasterGNNWeights, gnn::GNNModel)
    if master.use_gat
        master.gat_W1 .= gnn.gat_layers[1].W
        master.gat_a1 .= gnn.gat_layers[1].a
        master.gat_W2 .= gnn.gat_layers[2].W
        master.gat_a2 .= gnn.gat_layers[2].a
    else
        master.W1 .= gnn.W1
        master.W2 .= gnn.W2
    end
end

# ─────────────────────────────────────────────────────────────────────────────
# Softmax utilities
# ─────────────────────────────────────────────────────────────────────────────

function _softmax(z::AbstractVector{Float64})
    m = maximum(z)
    ex = exp.(z .- m)
    s = sum(ex)
    if s <= 0.0 || isnan(s)
        return fill(1.0 / length(z), length(z))
    end
    return ex ./ s
end

# ─────────────────────────────────────────────────────────────────────────────
# Core functions
# ─────────────────────────────────────────────────────────────────────────────

"""
    graph_level_readout(gnn::GNNModel)

Mean-pool hidden layer activations H1 into a single graph-level feature vector.
Returns a vector of length `hidden_dim`.
"""
function graph_level_readout(gnn::GNNModel)::Vector{Float64}
    return vec(mean(gnn.H1, dims=1))
end

"""
    predict_params(predictor, h_graph)

Predict parameter configuration from all 7 heads via argmax.
Returns `(actions, cge_params)`.
"""
function predict_params(predictor::ParamPredictor, h_graph::Vector{Float64})
    n_heads = length(predictor.heads)
    actions = zeros(Int, n_heads)

    for i in 1:n_heads
        logits = predictor.heads[i]' * h_graph .+ predictor.biases[i]
        clamp!(logits, -20.0, 20.0)
        probs = _softmax(logits)
        actions[i] = argmax(probs)
    end

    params = actions_to_cge_params(predictor, actions)
    return actions, params
end

"""
    actions_to_cge_params(predictor, actions)

Convert action indices to a concrete `CGEParams` struct.
Non-predicted parameters use best-known defaults.
"""
function actions_to_cge_params(predictor::ParamPredictor, actions::Vector{Int})
    window = predictor.head_options[1][actions[1]]::Int
    intervals = predictor.head_options[2][actions[2]]::Bool
    lr_split_raw = predictor.head_options[3][actions[3]]::Bool
    mil = predictor.head_options[4][actions[4]]::Int
    copy_adaptive = predictor.head_options[5][actions[5]]::Bool
    stop_deltas = predictor.head_options[6][actions[6]]::Bool
    zigzag = predictor.head_options[7][actions[7]]::Bool

    # Enforce dependency: lr_split requires intervals
    lr_split = intervals ? lr_split_raw : false

    return CGE.CGEParams(;
        intra_ref_window = window,
        intra_intervals = intervals,
        intra_lr_split = lr_split,
        intra_mil = mil,
        intra_adapt_mil = mil,
        intra_copy_adaptive = copy_adaptive,
        intra_stop_deltas = stop_deltas,
        intra_zigzag = zigzag,
        varint = :fibonacci,
        count_varint = :fibonacci,
        gap = :fibonacci,
        degree = :elias_delta,
        intra_copy_blocks = true,
        intra_ref_fixwidth = true,
        intra_ref_vlc = false,
        intra_add_adaptive = true,
        intra_raw_adaptive = true,
        membership = :implicit_ranges,
        intra_ref_enabled = true,
    )
end

"""
    compute_bpe(g, cge_params, clusters)

Run CGE `encode_level` on graph `g` with given parameters and cluster
partition. Returns bits per edge (BPE) as a Float64.
"""
function compute_bpe(g::AbstractGraph{T}, cge_params::CGE.CGEParams,
                     clusters::Vector{Vector{T}}) where {T<:Unsigned}
    m = ne(g)
    if m == 0
        return 0.0
    end

    io = IOBuffer()
    w = BitWriter(io)
    CGE.encode_level(w, g, clusters; params=cge_params)
    flush_bitwriter(w)

    total_bytes = position(io)
    return 8.0 * total_bytes / m
end

# ─────────────────────────────────────────────────────────────────────────────
# Grid search
# ─────────────────────────────────────────────────────────────────────────────

"""
    grid_search_best_params(g; clusters=nothing, verbose=false)

Exhaustive search over all 512 parameter combinations. Returns
`(best_actions, best_bpe, all_bpes)`.
"""
function grid_search_best_params(g::AbstractGraph{T};
                                  clusters::Union{Nothing, Vector{Vector{T}}}=nothing,
                                  verbose::Bool=false) where {T<:Unsigned}
    n = Int(nv(g))
    if clusters === nothing
        clusters = [collect(T(1):T(n))]
    end

    # Enumerate all 512 combinations
    n_heads = length(PARAM_HEADS)
    head_sizes = [length(opts) for (_, opts) in PARAM_HEADS]
    n_combos = prod(head_sizes)  # 4*2*2*4*2*2*2 = 512

    best_bpe = Inf
    best_actions = zeros(Int, n_heads)
    all_bpes = Float64[]

    # Build a dummy predictor to reuse actions_to_cge_params
    dummy = ParamPredictor(1; lr=0.0)

    combo_idx = 0
    for w_idx in 1:head_sizes[1]
        for int_idx in 1:head_sizes[2]
            for lr_idx in 1:head_sizes[3]
                for mil_idx in 1:head_sizes[4]
                    for ca_idx in 1:head_sizes[5]
                        for sd_idx in 1:head_sizes[6]
                            for zz_idx in 1:head_sizes[7]
                                combo_idx += 1
                                actions = [w_idx, int_idx, lr_idx, mil_idx,
                                           ca_idx, sd_idx, zz_idx]
                                params = actions_to_cge_params(dummy, actions)
                                bpe = compute_bpe(g, params, clusters)
                                push!(all_bpes, bpe)

                                if bpe < best_bpe
                                    best_bpe = bpe
                                    best_actions .= actions
                                end

                                if verbose && combo_idx % 64 == 0
                                    println("  Grid search: $combo_idx/$n_combos, best=$(round(best_bpe, digits=4))")
                                end
                            end
                        end
                    end
                end
            end
        end
    end

    if verbose
        println("  Grid search complete: best BPE=$(round(best_bpe, digits=4)) ($n_combos combos)")
    end

    return best_actions, best_bpe, all_bpes
end

"""
    generate_training_labels(graphs; verbose=false)

Grid-search each graph to find optimal parameter actions. Returns a vector of
best-action vectors (one per graph).
"""
function generate_training_labels(graphs::Vector{<:AbstractGraph{<:Unsigned}};
                                   verbose::Bool=false)
    labels = Vector{Int}[]
    for (i, g) in enumerate(graphs)
        if verbose
            println("  Labeling graph $i/$(length(graphs)): $(nv(g))v, $(ne(g))e")
        end
        best_actions, best_bpe, _ = grid_search_best_params(g; verbose=false)
        push!(labels, best_actions)
        if verbose
            println("    best BPE=$(round(best_bpe, digits=4)), actions=$best_actions")
        end
    end
    return labels
end

# ─────────────────────────────────────────────────────────────────────────────
# Supervised backward pass
# ─────────────────────────────────────────────────────────────────────────────

"""
    supervised_loss_backward!(predictor, gnn, h_graph, target_actions; gnn_lr)

Cross-entropy gradient for each head + backprop into GNN via mean-pool.

For each head i, the loss is: -log(softmax(logits)[target_i]).
Gradient w.r.t. logits: softmax(logits) - one_hot(target_i).
"""
function supervised_loss_backward!(predictor::ParamPredictor, gnn::GNNModel,
                                    h_graph::Vector{Float64},
                                    target_actions::Vector{Int};
                                    gnn_lr::Float64=0.0001)
    n_heads = length(predictor.heads)
    hidden_dim = length(h_graph)

    dL_dh = zeros(Float64, hidden_dim)
    total_loss = 0.0

    for i in 1:n_heads
        W = predictor.heads[i]    # (feat_dim × n_opts)
        b = predictor.biases[i]   # (n_opts,)
        n_opts = length(b)
        target = target_actions[i]

        # Forward: logits = W' * h + b, probs = softmax(logits)
        logits = W' * h_graph .+ b
        clamp!(logits, -20.0, 20.0)
        probs = _softmax(logits)

        # Cross-entropy loss contribution
        total_loss -= log(probs[target] + 1e-10)

        # Gradient w.r.t. logits: probs - one_hot(target)
        dlogits = copy(probs)
        dlogits[target] -= 1.0

        # Clip gradient per head
        dlogits_norm = norm(dlogits)
        if dlogits_norm > 5.0
            dlogits .*= 5.0 / dlogits_norm
        end

        # Update head weights: W -= lr * h_graph * dlogits'
        for k in 1:n_opts
            predictor.heads[i][:, k] .-= predictor.lr .* h_graph .* dlogits[k]
        end
        predictor.biases[i] .-= predictor.lr .* dlogits

        # Gradient w.r.t. h_graph: dL/dh += W * dlogits
        dL_dh .+= W * dlogits
    end

    # Clip aggregate gradient
    dh_norm = norm(dL_dh)
    if dh_norm > 1.0
        dL_dh .*= 1.0 / dh_norm
    end

    # Backprop into GNN: h_graph = mean(H1, dims=1)
    n = gnn.n
    dL_dH1 = repeat(dL_dh' ./ n, n, 1)  # (n × hidden_dim)
    dL_dZ2 = zeros(Float64, n)

    gnn_backward_rl!(gnn, dL_dZ2, dL_dH1; lr=gnn_lr)

    return total_loss
end

# ─────────────────────────────────────────────────────────────────────────────
# Training loop
# ─────────────────────────────────────────────────────────────────────────────

"""
    train_param_predictor!(master, predictor, graphs, labels; epochs, lr,
                           feature_set, gnn_lr, verbose, seed)

Supervised training over the dataset for multiple epochs.
Each epoch iterates over all graphs, computing cross-entropy loss against
grid-search-optimal labels and updating predictor heads + GNN weights.
"""
function train_param_predictor!(master::MasterGNNWeights,
                                 predictor::ParamPredictor,
                                 graphs::Vector{<:AbstractGraph{<:Unsigned}},
                                 labels::Vector{Vector{Int}};
                                 epochs::Int=50,
                                 verbose::Bool=true,
                                 feature_set::Symbol=:enhanced,
                                 gnn_lr::Float64=0.0001,
                                 seed::Int=42)
    @assert length(graphs) == length(labels)
    n_graphs = length(graphs)

    loss_history = Float64[]
    best_avg_loss = Inf

    for epoch in 1:epochs
        # Shuffle order each epoch
        rng = MersenneTwister(seed + epoch)
        perm = randperm(rng, n_graphs)
        epoch_loss = 0.0

        for idx in perm
            g = graphs[idx]
            target = labels[idx]
            n = Int(nv(g))

            if ne(g) == 0
                continue
            end

            # Create per-graph GNNModel, copy master weights in
            gnn = GNNModel(g; hidden_dim=master.hidden_dim, seed=seed + idx,
                            feature_set=feature_set, use_gat=master.use_gat)
            copy_master_to_gnn!(gnn, master)

            # Forward pass
            gnn_forward(gnn)
            h_graph = graph_level_readout(gnn)

            # Supervised backward
            loss = supervised_loss_backward!(predictor, gnn, h_graph, target;
                                              gnn_lr=gnn_lr)
            epoch_loss += loss

            # Copy updated GNN weights back to master
            copy_gnn_to_master!(master, gnn)
        end

        avg_loss = epoch_loss / n_graphs
        push!(loss_history, avg_loss)

        if avg_loss < best_avg_loss
            best_avg_loss = avg_loss
        end

        if verbose && (epoch % max(1, epochs ÷ 20) == 0 || epoch == 1)
            println("  Epoch $epoch/$epochs: avg_loss=$(round(avg_loss, digits=4)), best=$(round(best_avg_loss, digits=4))")
        end
    end

    if verbose
        println("\nTraining complete. Best avg loss: $(round(best_avg_loss, digits=4))")
    end

    return loss_history
end

# ─────────────────────────────────────────────────────────────────────────────
# Evaluation
# ─────────────────────────────────────────────────────────────────────────────

"""
    evaluate_param_predictor(master, predictor, g; feature_set, seed)

Evaluate the trained predictor on a graph using argmax prediction.
Returns `(bpe, cge_params)`.
"""
function evaluate_param_predictor(master::MasterGNNWeights,
                                   predictor::ParamPredictor,
                                   g::AbstractGraph{T};
                                   feature_set::Symbol=:enhanced,
                                   seed::Int=42) where {T<:Unsigned}
    gnn = GNNModel(g; hidden_dim=master.hidden_dim, seed=seed,
                    feature_set=feature_set, use_gat=master.use_gat)
    copy_master_to_gnn!(gnn, master)
    gnn_forward(gnn)

    h_graph = graph_level_readout(gnn)
    actions, params = predict_params(predictor, h_graph)

    clusters = [collect(T(1):T(nv(g)))]
    bpe = compute_bpe(g, params, clusters)

    return bpe, params
end

# ─────────────────────────────────────────────────────────────────────────────
# Save / Load
# ─────────────────────────────────────────────────────────────────────────────

function save_param_predictor(predictor::ParamPredictor, master::MasterGNNWeights,
                               filepath::AbstractString)
    @save filepath predictor master
end

function load_param_predictor(filepath::AbstractString)
    @load filepath predictor master
    return predictor, master
end

end # module ParamPrediction
