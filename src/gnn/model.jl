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

mutable struct GNNModel
    W1::Matrix{Float64}                    # (input_dim × hidden_dim)
    W2::Matrix{Float64}                    # (hidden_dim × 1)
    A_hat::SparseMatrixCSC{Float64,Int}    # Normalized adjacency (precomputed)
    n::Int                                 # Number of vertices
    input_dim::Int
    hidden_dim::Int
    use_gat::Bool                          # Use GAT instead of GCN
    gat_layers::Vector{GATLayer}           # GAT layers
    # Caches for backprop
    X::Matrix{Float64}                     # Input features (N × input_dim)
    Z1::Matrix{Float64}                    # Pre-activation layer 1
    H1::Matrix{Float64}                    # Post-activation layer 1
    Z2::Matrix{Float64}                    # Output scores (N × 1)
    # Gradients
    dW1::Matrix{Float64}
    dW2::Matrix{Float64}
end

"""
    GNNModel(g; hidden_dim=32, seed=42, X=nothing, feature_set=:enhanced)

Build a 2-layer GNN from graph `g`. Computes normalized adjacency Â = D^{-1/2}(A+I)D^{-1/2},
extracts node features, and Xavier-initializes weights.

Feature sets:
- `:basic` → 4-D features (degrees, clustering)
- `:enhanced` → 10-D features (+ PageRank, neighbor stats, spectral)
- `:extended` → 12-D features (+ reciprocity, avg neighbor clustering)
- `:compression` → 19-D features (+ gap stats, HITS, core number, interval fraction)

If `X` is provided, it is used directly as the feature matrix (N × d) instead of computing
features internally.
"""
function GNNModel(g::AbstractGraph{T}; hidden_dim::Int=32, seed::Int=42,
                  X::Union{Matrix{Float64},Nothing}=nothing,
                  feature_set::Symbol=:enhanced, use_gat::Bool=false) where {T<:Unsigned}
    rng = MersenneTwister(seed)
    n = Int(nv(g))

    # Build normalized adjacency Â = D^{-1/2}(A+I)D^{-1/2}
    A = _build_sparse_adj(g)
    A_plus_I = A + sparse(1:n, 1:n, ones(Float64, n), n, n)
    deg_vec = vec(sum(A_plus_I, dims=2))
    inv_sqrt_deg = 1.0 ./ sqrt.(max.(deg_vec, 1.0))
    D_inv_sqrt = sparse(1:n, 1:n, inv_sqrt_deg, n, n)
    A_hat = D_inv_sqrt * A_plus_I * D_inv_sqrt

    # Extract or use provided node features
    if X === nothing
        if feature_set == :basic
            X = _extract_node_features(g, n)
        elseif feature_set == :extended
            X = _extract_extended_features(g, n, A_hat)
        elseif feature_set == :compression
            X = _extract_compression_features(g, n, A_hat)
        else
            X = _extract_enhanced_features(g, n, A_hat)
        end
    end
    input_dim = size(X, 2)

    W1, W2 = Matrix{Float64}(undef,0,0), Matrix{Float64}(undef,0,0)
    gat_layers = GATLayer[]
    if use_gat
        gat_layer1 = GATLayer(input_dim, hidden_dim; seed=seed)
        gat_layer2 = GATLayer(hidden_dim, 1; concat=false, seed=seed+1)
        push!(gat_layers, gat_layer1)
        push!(gat_layers, gat_layer2)
    else
        # Xavier initialization for GCN
        W1 = randn(rng, input_dim, hidden_dim) .* sqrt(2.0 / (input_dim + hidden_dim))
        W2 = randn(rng, hidden_dim, 1) .* sqrt(2.0 / (hidden_dim + 1))
    end

    # Allocate caches
    Z1 = zeros(n, hidden_dim)
    H1 = zeros(n, hidden_dim)
    Z2 = zeros(n, 1)
    dW1 = zeros(input_dim, hidden_dim)
    dW2 = zeros(hidden_dim, 1)

    return GNNModel(W1, W2, A_hat, n, input_dim, hidden_dim, use_gat, gat_layers, X, Z1, H1, Z2, dW1, dW2)
end

function _build_sparse_adj(g::AbstractGraph{T}) where {T<:Unsigned}
    n = Int(nv(g))
    rows = Int[]
    cols = Int[]
    for v in vertices(g)
        for u in outneighbors(g, v)
            push!(rows, Int(v))
            push!(cols, Int(u))
        end
    end
    return sparse(rows, cols, ones(Float64, length(rows)), n, n)
end

function _extract_node_features(g::AbstractGraph{T}, n::Int) where {T<:Unsigned}
    X = zeros(n, 4)

    # Compute in-degrees in O(E) single pass
    in_degrees = zeros(Int, n)
    for v in vertices(g)
        for u in outneighbors(g, v)
            in_degrees[Int(u)] += 1
        end
    end

    # Compute per-vertex features
    for v in vertices(g)
        vi = Int(v)
        out_deg = length(outneighbors(g, v))
        in_deg = in_degrees[vi]
        total_deg = out_deg + in_deg
        cc = _local_clustering_coeff(g, v)
        X[vi, 1] = log(1.0 + out_deg)
        X[vi, 2] = log(1.0 + in_deg)
        X[vi, 3] = log(1.0 + total_deg)
        X[vi, 4] = cc
    end
    return X
end

"""
    _extract_enhanced_features(g, n, A_hat)

Compute 10-dimensional enhanced node features:
  1. log(1 + out_deg)
  2. log(1 + in_deg)
  3. log(1 + total_deg)
  4. clustering_coeff
  5. log(1 + pagerank × n)
  6. log(1 + mean_neighbor_out_deg)
  7. log(1 + max_neighbor_out_deg)
  8-10. spectral coordinates (top-3 eigenvectors of Â)
"""
function _extract_enhanced_features(g::AbstractGraph{T}, n::Int,
                                     A_hat::SparseMatrixCSC{Float64,Int}) where {T<:Unsigned}
    X = zeros(n, 10)

    # Compute in-degrees and out-degrees in O(E)
    in_degrees = zeros(Int, n)
    out_degrees = zeros(Int, n)
    for v in vertices(g)
        od = length(outneighbors(g, v))
        out_degrees[Int(v)] = od
        for u in outneighbors(g, v)
            in_degrees[Int(u)] += 1
        end
    end

    # Features 1-4: degree + clustering (same as before)
    for v in vertices(g)
        vi = Int(v)
        out_deg = out_degrees[vi]
        in_deg = in_degrees[vi]
        total_deg = out_deg + in_deg
        cc = _local_clustering_coeff(g, v)
        X[vi, 1] = log(1.0 + out_deg)
        X[vi, 2] = log(1.0 + in_deg)
        X[vi, 3] = log(1.0 + total_deg)
        X[vi, 4] = cc
    end

    # Feature 5: PageRank via power iteration on transition matrix
    @info "Computing PageRank features..."
    P = get_sparse_P_matrix(g)
    pr = PR(P)
    for vi in 1:n
        X[vi, 5] = log(1.0 + pr[vi] * n)
    end

    # Features 6-7: neighbor degree statistics
    @info "Computing neighbor degree features..."
    for v in vertices(g)
        vi = Int(v)
        nbs = outneighbors(g, v)
        if length(nbs) > 0
            neighbor_out_degs = [Float64(out_degrees[Int(u)]) for u in nbs]
            X[vi, 6] = log(1.0 + sum(neighbor_out_degs) / length(neighbor_out_degs))
            X[vi, 7] = log(1.0 + maximum(neighbor_out_degs))
        end
    end

    # Features 8-10: spectral coordinates (top-3 eigenvectors of Â)
    @info "Computing spectral features (3 eigenvectors, 30 iterations each)..."
    spectral = _compute_spectral_features(A_hat, 3)
    X[:, 8:10] = spectral

    @info "Enhanced features: $(size(X, 2)) dimensions for $n vertices"
    return X
end

"""
    _extract_extended_features(g, n, A_hat)

Compute 12-dimensional extended node features (enhanced + 2 extra):
  1-10. Same as enhanced features
  11. reciprocity ratio = |{v : (u,v) ∈ E ∧ (v,u) ∈ E}| / out_deg
  12. log(1 + avg_neighbor_clustering_coeff)
"""
function _extract_extended_features(g::AbstractGraph{T}, n::Int,
                                     A_hat::SparseMatrixCSC{Float64,Int}) where {T<:Unsigned}
    # Start with the 10 enhanced features
    X_base = _extract_enhanced_features(g, n, A_hat)
    X = zeros(n, 12)
    X[:, 1:10] = X_base

    # Build reverse edge set for reciprocity
    @info "Computing extended features (reciprocity, neighbor clustering)..."
    out_sets = Dict{Int, Set{Int}}()
    for v in vertices(g)
        out_sets[Int(v)] = Set{Int}(Int(u) for u in outneighbors(g, v))
    end

    for v in vertices(g)
        vi = Int(v)
        nbs = outneighbors(g, v)
        out_deg = length(nbs)

        # Feature 11: reciprocity ratio
        if out_deg > 0
            reciprocal_count = count(u -> vi in get(out_sets, Int(u), Set{Int}()), nbs)
            X[vi, 11] = reciprocal_count / out_deg
        end

        # Feature 12: log(1 + avg neighbor clustering coeff)
        if out_deg > 0
            cc_sum = 0.0
            for u in nbs
                cc_sum += _local_clustering_coeff(g, u)
            end
            X[vi, 12] = log(1.0 + cc_sum / out_deg)
        end
    end

    @info "Extended features: $(size(X, 2)) dimensions for $n vertices"
    return X
end

"""
    _extract_compression_features(g, n, A_hat)

Compute 19-dimensional compression-aware node features (extended + 7 extra):
  1-12. Same as extended features
  13. Gap variance: var(log.(1 .+ gaps)) of sorted neighbor ID gaps
  14. Reference overlap potential: avg Jaccard similarity with k nearest predecessors
  15. Hub score (10-iter HITS)
  16. Authority score (10-iter HITS)
  17. Core number (k-core decomposition)
  18. Log gap mean: mean(log.(1 .+ gaps))
  19. Interval fraction: fraction of neighbors in consecutive intervals >= 2
"""
function _extract_compression_features(g::AbstractGraph{T}, n::Int,
                                        A_hat::SparseMatrixCSC{Float64,Int}) where {T<:Unsigned}
    # Start with the 12 extended features
    X_base = _extract_extended_features(g, n, A_hat)
    X = zeros(n, 19)
    X[:, 1:12] = X_base

    @info "Computing compression-aware features (gap stats, HITS, core, intervals)..."

    # Features 15-16: HITS hub and authority scores
    hub, auth = _compute_hits(g, 10)
    for vi in 1:n
        X[vi, 15] = hub[vi]
        X[vi, 16] = auth[vi]
    end

    # Feature 17: core numbers
    core = _compute_core_numbers(g)
    for vi in 1:n
        X[vi, 17] = Float64(core[vi])
    end

    # Features 13, 14, 18, 19: gap-based and interval features (per-vertex)
    # Pre-collect sorted out-neighbor lists for Jaccard reference overlap
    sorted_nbs = Vector{Vector{Int}}(undef, n)
    for v in vertices(g)
        sorted_nbs[Int(v)] = sort(Int.(collect(outneighbors(g, v))))
    end

    for v in vertices(g)
        vi = Int(v)
        nbs = sorted_nbs[vi]
        deg = length(nbs)

        if deg >= 2
            # Compute gaps between consecutive sorted neighbors
            gaps = Float64[nbs[i] - nbs[i-1] for i in 2:deg]
            log_gaps = log.(1.0 .+ gaps)
            X[vi, 13] = length(log_gaps) >= 2 ? var(log_gaps) : 0.0
            X[vi, 18] = mean(log_gaps)    # log gap mean
        end

        # Feature 14: reference overlap potential (avg Jaccard with k=5 nearest predecessors)
        if vi > 1 && deg > 0
            nbs_set = Set(nbs)
            k_pred = min(5, vi - 1)
            jaccard_sum = 0.0
            for offset in 1:k_pred
                pred_vi = vi - offset
                pred_nbs = sorted_nbs[pred_vi]
                if isempty(pred_nbs)
                    continue
                end
                pred_set = Set(pred_nbs)
                isect = length(intersect(nbs_set, pred_set))
                union_size = length(union(nbs_set, pred_set))
                if union_size > 0
                    jaccard_sum += isect / union_size
                end
            end
            X[vi, 14] = jaccard_sum / k_pred
        end

        # Feature 19: interval fraction (fraction of neighbors in consecutive runs >= 2)
        if deg >= 2
            in_interval = 0
            i = 1
            while i <= deg
                run_len = 1
                while i + run_len <= deg && nbs[i + run_len] == nbs[i] + run_len
                    run_len += 1
                end
                if run_len >= 2
                    in_interval += run_len
                end
                i += run_len
            end
            X[vi, 19] = in_interval / deg
        end
    end

    @info "Compression features: $(size(X, 2)) dimensions for $n vertices"
    return X
end

"""
    _compute_hits(g, num_iters)

HITS algorithm returning (hub_scores, authority_scores) as length-n Float64 vectors.
"""
function _compute_hits(g::AbstractGraph{T}, num_iters::Int) where {T<:Unsigned}
    n = Int(nv(g))
    hub = ones(Float64, n)
    auth = ones(Float64, n)

    for _ in 1:num_iters
        # Authority update: auth[v] = sum of hub[u] for all u -> v
        new_auth = zeros(Float64, n)
        for u in vertices(g)
            ui = Int(u)
            for v in outneighbors(g, u)
                new_auth[Int(v)] += hub[ui]
            end
        end
        # Normalize
        auth_norm = sqrt(sum(new_auth .^ 2))
        if auth_norm > 0
            new_auth ./= auth_norm
        end

        # Hub update: hub[u] = sum of auth[v] for all u -> v
        new_hub = zeros(Float64, n)
        for u in vertices(g)
            ui = Int(u)
            for v in outneighbors(g, u)
                new_hub[ui] += new_auth[Int(v)]
            end
        end
        # Normalize
        hub_norm = sqrt(sum(new_hub .^ 2))
        if hub_norm > 0
            new_hub ./= hub_norm
        end

        hub = new_hub
        auth = new_auth
    end

    return hub, auth
end

"""
    _compute_core_numbers(g)

K-core decomposition based on out-degree. Returns a length-n Int vector of core numbers.
"""
function _compute_core_numbers(g::AbstractGraph{T}) where {T<:Unsigned}
    n = Int(nv(g))
    # Compute out-degrees
    deg = zeros(Int, n)
    for v in vertices(g)
        deg[Int(v)] = length(outneighbors(g, v))
    end

    core = copy(deg)
    # Sort vertices by degree
    order = sortperm(core)
    pos = zeros(Int, n)  # position of each vertex in order
    for (i, v) in enumerate(order)
        pos[v] = i
    end

    # Bin-sort based peeling
    for i in 1:n
        v = order[i]
        for u in outneighbors(g, T(v))
            ui = Int(u)
            if core[ui] > core[v]
                # Decrease core number of neighbor
                core[ui] -= 1
                # Move neighbor forward in the order (swap with first vertex of its old degree)
                # This is the standard k-core peeling algorithm
            end
        end
    end

    return core
end

"""
    _compute_spectral_features(A_hat, k; num_iters=30, seed=42)

Compute top-k approximate eigenvectors of Â via power iteration with deflation.
Returns an N × k matrix of spectral coordinates.
"""
function _compute_spectral_features(A_hat::SparseMatrixCSC{Float64,Int}, k::Int;
                                     num_iters::Int=30, seed::Int=42)
    rng = MersenneTwister(seed)
    n = size(A_hat, 1)
    vectors = Vector{Vector{Float64}}()

    for i in 1:k
        v = randn(rng, n)
        v ./= norm(v)
        for iter in 1:num_iters
            v = A_hat * v
            # Deflate: project out previous eigenvectors
            for prev in vectors
                v .-= dot(v, prev) .* prev
            end
            v ./= norm(v)
        end
        push!(vectors, v)
    end

    return hcat(vectors...)
end

function _local_clustering_coeff(g::AbstractGraph, v)
    neighbors = outneighbors(g, v)
    k = length(neighbors)
    if k < 2
        return 0.0
    end
    if k > 200
        return 0.0
    end
    triangles = 0
    for i in 1:k
        ni_out_set = Set(outneighbors(g, neighbors[i]))
        for j in (i+1):k
            if neighbors[j] in ni_out_set
                triangles += 1
            end
        end
    end
    return 2.0 * triangles / (k * (k - 1))
end

"""
    gnn_forward(model::GNNModel)::Vector{Float64}

Run forward pass through the 2-layer GNN. Returns N-dimensional score vector.
"""
function gnn_forward(model::GNNModel; training::Bool=true)::Vector{Float64}
    if model.use_gat
        model.H1 = gat_forward(model.gat_layers[1], model.X, model.A_hat; training=training)
        model.Z2 = gat_forward(model.gat_layers[2], model.H1, model.A_hat; training=training)
    else
        model.Z1 = model.A_hat * model.X * model.W1
        model.H1 = max.(model.Z1, 0.0)
        model.Z2 = model.A_hat * model.H1 * model.W2
    end
    return vec(model.Z2)
end

"""
    gnn_backward!(model::GNNModel, dL_dZ2::Vector{Float64}; lr::Float64=0.001)

Backpropagate through the 2-layer GNN and update weights.
`dL_dZ2` is the gradient of the loss w.r.t. the output scores (length N).
"""
function gnn_backward!(model::GNNModel, dL_dZ2::Vector{Float64}; lr::Float64=0.001)
    if model.use_gat
        dZ2 = reshape(dL_dZ2, model.n, 1)
        dH1 = gat_backward!(model.gat_layers[2], dZ2, model.H1, model.A_hat)
        gat_backward!(model.gat_layers[1], dH1, model.X, model.A_hat)

        # SGD update for GAT layers
        model.gat_layers[1].W .-= lr .* model.gat_layers[1].dW
        model.gat_layers[1].a .-= lr .* model.gat_layers[1].da
        model.gat_layers[2].W .-= lr .* model.gat_layers[2].dW
        model.gat_layers[2].a .-= lr .* model.gat_layers[2].da
    else
        dZ2 = reshape(dL_dZ2, model.n, 1)

        # dL/dW2 = (Â · H1)ᵀ · dL/dZ2
        AH1 = model.A_hat * model.H1
        model.dW2 = AH1' * dZ2

        # dL/dH1 = Âᵀ · dL/dZ2 · W2ᵀ
        dH1 = (model.A_hat' * dZ2) * model.W2'

        # dL/dZ1 = dL/dH1 .* (Z1 .> 0)  (ReLU gradient)
        dZ1 = dH1 .* (model.Z1 .> 0)

        # dL/dW1 = (Â · X)ᵀ · dL/dZ1
        AX = model.A_hat * model.X
        model.dW1 = AX' * dZ1

        # Gradient clipping (max norm per parameter matrix)
        max_grad_norm = 1.0
        dW1_norm = sqrt(sum(model.dW1 .^ 2))
        if dW1_norm > max_grad_norm
            model.dW1 .*= max_grad_norm / dW1_norm
        end
        dW2_norm = sqrt(sum(model.dW2 .^ 2))
        if dW2_norm > max_grad_norm
            model.dW2 .*= max_grad_norm / dW2_norm
        end

        # SGD update
        model.W1 .-= lr .* model.dW1
        model.W2 .-= lr .* model.dW2
    end
end

"""
    gnn_backward_rl!(model, dL_dZ2_ext, dL_dH1_ext; lr=0.001)

Backpropagate RL gradients through the 2-layer GNN and update weights.
Accepts external gradients w.r.t. both Z2 (output, N×1) and H1 (hidden, N×hidden_dim),
since the RL feature vector is [H1[v,:]; Z2[v]].
"""
function gnn_backward_rl!(model::GNNModel, dL_dZ2_ext::Vector{Float64},
                          dL_dH1_ext::Matrix{Float64}; lr::Float64=0.001)
    if model.use_gat
        dZ2 = reshape(dL_dZ2_ext, model.n, 1)
        dH1 = gat_backward!(model.gat_layers[2], dZ2, model.H1, model.A_hat)
        dH1 .+= dL_dH1_ext
        gat_backward!(model.gat_layers[1], dH1, model.X, model.A_hat)

        # SGD update for GAT layers
        model.gat_layers[1].W .-= lr .* model.gat_layers[1].dW
        model.gat_layers[1].a .-= lr .* model.gat_layers[1].da
        model.gat_layers[2].W .-= lr .* model.gat_layers[2].dW
        model.gat_layers[2].a .-= lr .* model.gat_layers[2].da
    else
        dZ2 = reshape(dL_dZ2_ext, model.n, 1)

        # dL/dW2 = (Â · H1)ᵀ · dL/dZ2
        AH1 = model.A_hat * model.H1
        model.dW2 = AH1' * dZ2

        # Combine Z2-chain gradient with direct H1 gradient from RL
        dH1 = (model.A_hat' * dZ2) * model.W2' .+ dL_dH1_ext

        # ReLU gradient
        dZ1 = dH1 .* (model.Z1 .> 0)

        # dL/dW1 = (Â · X)ᵀ · dL/dZ1
        AX = model.A_hat * model.X
        model.dW1 = AX' * dZ1

        # Gradient clipping (max norm per parameter matrix)
        max_grad_norm = 1.0
        dW1_norm = sqrt(sum(model.dW1 .^ 2))
        if dW1_norm > max_grad_norm
            model.dW1 .*= max_grad_norm / dW1_norm
        end
        dW2_norm = sqrt(sum(model.dW2 .^ 2))
        if dW2_norm > max_grad_norm
            model.dW2 .*= max_grad_norm / dW2_norm
        end

        # SGD update
        model.W1 .-= lr .* model.dW1
        model.W2 .-= lr .* model.dW2
    end
end

"""
    save_gnn_model(model::GNNModel, filepath::AbstractString)

Save GNN weights (W1, W2) in binary format. Â is recomputed from graph on load.
"""
function save_gnn_model(model::GNNModel, filepath::AbstractString)
    open(filepath, "w") do f
        write(f, model.use_gat)
        if model.use_gat
            write(f, Int32(length(model.gat_layers)))
            for layer in model.gat_layers
                write(f, Int32(layer.in_dim))
                write(f, Int32(layer.out_dim))
                write(f, layer.W)
                write(f, layer.a)
            end
        else
            write(f, Int32(model.input_dim))
            write(f, Int32(model.hidden_dim))
            write(f, Int32(model.n))
            write(f, model.W1)
            write(f, model.W2)
        end
    end
end

"""
    load_gnn_weights!(model::GNNModel, filepath::AbstractString)

Load GNN weights from file into an existing model.
"""
function load_gnn_weights!(model::GNNModel, filepath::AbstractString)
    open(filepath, "r") do f
        use_gat = read(f, Bool)
        model.use_gat = use_gat
        if use_gat
            num_layers = read(f, Int32)
            for i in 1:num_layers
                in_dim = read(f, Int32)
                out_dim = read(f, Int32)
                W = reshape(reinterpret(Float64, read(f, sizeof(Float64) * in_dim * out_dim)),
                              in_dim, out_dim)
                a = reshape(reinterpret(Float64, read(f, sizeof(Float64) * 2 * out_dim)),
                              2 * out_dim, 1)
                if i > length(model.gat_layers)
                    # This case happens if the model is not initialized with use_gat=true
                    # We can create a new layer, but the seed will be different
                    push!(model.gat_layers, GATLayer(in_dim, out_dim))
                end
                model.gat_layers[i].W = W
                model.gat_layers[i].a = a
            end
        else
            input_dim = read(f, Int32)
            hidden_dim = read(f, Int32)
            _n = read(f, Int32)
            if input_dim != model.input_dim || hidden_dim != model.hidden_dim
                error("Model dimensions mismatch: file has ($input_dim, $hidden_dim), model has ($(model.input_dim), $(model.hidden_dim))")
            end
            model.W1 = reshape(reinterpret(Float64, read(f, sizeof(Float64) * model.input_dim * model.hidden_dim)),
                              model.input_dim, model.hidden_dim)
            model.W2 = reshape(reinterpret(Float64, read(f, sizeof(Float64) * model.hidden_dim * 1)),
                              model.hidden_dim, 1)
        end
    end
end

