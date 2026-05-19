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
    NDFEncoder(d_in, hidden; dropout=0.1f0)

Per-vertex MLP that maps a `d_in`-dimensional feature vector to a
`hidden`-dimensional representation. Applied row-wise to a feature matrix
of shape (n_vertices × d_in) or batched array (n_vertices × d_in × B).
"""
function NDFEncoder(d_in::Int, hidden::Int; dropout::Float32=0.1f0)
    return Chain(
        Dense(d_in => hidden, relu),
        Dropout(dropout),
        Dense(hidden => hidden),
    )
end

"""
    NDF(d_in, n_classes; hidden=64, K=10, α=0.15f0, dropout=0.1f0, readout=:mean)

Neural Diffusion Fingerprint — APPNP-style generalization of the 2014
"Diffusion Fingerprints" construction.

Pipeline:
    H_0 = encoder(Φ)                 # per-vertex MLP
    H_0 = H_0 ⊙ seed_mask            # optional: zero non-seed rows
    Z_{k+1} = (1-α) Â Z_k + α H_0     # K APPNP propagation steps
    f       = readout(Z_K)            # pool over vertices
    ŷ       = classifier(f)           # optional

The original linear DF is recovered when `d_in=1`, `hidden=1`, encoder is the
identity, `K→∞`, and Φ is the personalization vector v_k.

# Arguments
- `d_in`: per-vertex input feature dim (e.g. 1 for seed-only, 1+d_struct for seed+features).
- `n_classes`: number of output classes; pass 0 to return the fingerprint directly.
- `hidden`: hidden dim shared by encoder and propagation.
- `K`: number of propagation steps. v1 DF used K→∞; APPNP defaults to 10.
- `α`: teleportation constant (v1 used 0.15).
- `readout`: how to pool the propagated state Z before the classifier.
    - `:mean` / `:sum`: pool over all vertices (permutation-invariant; good for
      graph-level tasks where total/average diffusion matters).
    - `:seed_mean`: pool over seed vertices only (requires `seed_mask` at call time).
    - `:flatten`: flatten Z to `(n_vertices * hidden,)`. Use when the SPATIAL
      pattern of diffusion is what discriminates (per-vertex weights; matches
      v1 DF's setup where the full n-dim PPR vector was fed to a classifier).
      Requires `n_vertices` so the classifier head can be sized.
- `n_vertices`: number of vertices in the domain graph. Required when
      `readout=:flatten` so the classifier has the right input dim; ignored otherwise.
"""
struct NDF{E, C}
    encoder::E
    classifier::C
    K::Int
    α::Float32
    readout::Symbol
end
Flux.@layer NDF

function NDF(d_in::Int, n_classes::Int;
             hidden::Int=64, K::Int=10, α::Float32=0.15f0,
             dropout::Float32=0.1f0, readout::Symbol=:mean,
             n_vertices::Int=0)
    readout ∈ (:mean, :sum, :seed_mean, :flatten) ||
        throw(ArgumentError("readout must be :mean, :sum, :seed_mean, or :flatten"))
    if readout === :flatten && n_vertices <= 0
        throw(ArgumentError("readout=:flatten requires n_vertices > 0"))
    end
    encoder = NDFEncoder(d_in, hidden; dropout=dropout)
    classifier_in = readout === :flatten ? n_vertices * hidden : hidden
    classifier = n_classes > 0 ?
        Chain(Dropout(dropout), Dense(classifier_in => n_classes)) :
        identity
    return NDF(encoder, classifier, K, α, readout)
end

# Sparse-aware SpMV wrapper. Zygote's generic `*` adjoint for `A * X` computes
# `∂L/∂A = ΔY * X^T`, which for our 60K-vertex APPNP graph materializes a
# 60K × 60K dense (~14 GB) matrix — even though Â is not a parameter. The
# custom rrule below returns `NoTangent()` for A and only the cheap
# `A^T * ΔY` (sparse * dense) gradient for X. Drops peak memory from 17 GB
# to <1 GB on the modern-pipeline 60K-vertex graph.
_spmm(A::AbstractMatrix, X::AbstractMatrix) = A * X
function ChainRulesCore.rrule(::typeof(_spmm), A::AbstractMatrix, X::AbstractMatrix)
    Y = A * X
    function _spmm_pullback(ΔY)
        return ChainRulesCore.NoTangent(), ChainRulesCore.NoTangent(), A' * ΔY
    end
    return Y, _spmm_pullback
end

"""
    propagate(Â, H_0, K, α)

Run K APPNP propagation steps `Z ← (1-α) Â Z + α H_0` starting from Z_0 = H_0.
Accepts `H_0` of shape `(n × hidden)` for a single document, or
`(n × hidden × B)` for a batch of B documents.
"""
function propagate(Â::AbstractMatrix, H_0::AbstractMatrix, K::Int, α::Float32)
    Z = H_0
    one_minus_α = 1.0f0 - α
    for _ in 1:K
        Z = one_minus_α .* _spmm(Â, Z) .+ α .* H_0
    end
    return Z
end

function propagate(Â::AbstractMatrix, H_0::AbstractArray{T,3}, K::Int, α::Float32) where T
    n, h, B = size(H_0)
    H_0_flat = reshape(H_0, n, h * B)
    Z = H_0_flat
    one_minus_α = 1.0f0 - α
    for _ in 1:K
        Z = one_minus_α .* _spmm(Â, Z) .+ α .* H_0_flat
    end
    return reshape(Z, n, h, B)
end

# Differentiable masked mean. mask is converted to Float32 and the result is
# `sum(Z .* mask) / sum(mask)`, avoiding non-differentiable logical indexing.
function _masked_mean(Z::AbstractMatrix, mask::AbstractVector{Bool})
    mask_f = Float32.(mask)
    norm = max(sum(mask_f), 1.0f0)
    return vec(sum(Z .* reshape(mask_f, :, 1); dims=1)) ./ norm
end

function _masked_mean(Z::AbstractArray{T,3}, mask::AbstractMatrix{Bool}) where T
    n, _, B = size(Z)
    mask_f = Float32.(mask)
    mask_3d = reshape(mask_f, n, 1, B)
    weighted = dropdims(sum(Z .* mask_3d; dims=1); dims=1)          # (h × B)
    norm = reshape(max.(sum(mask_f; dims=1), 1.0f0), 1, B)          # (1 × B)
    return weighted ./ norm
end

"""
    (m::NDF)(Φ, Â; seed_mask=nothing)

Forward pass. `Φ` is the per-vertex feature input — `(n × d_in)` for a single
document or `(n × d_in × B)` for a batch. `seed_mask` is a Bool vector of
length n for single documents or a `(n × B)` Bool matrix for batches; when
provided, non-seed rows of H_0 are zeroed (preserving DF's "only seed
vertices teleport" semantics) and required when `readout == :seed_mean`.

Returns logits `(n_classes,)` or `(n_classes × B)` when a classifier is
configured, otherwise the pooled fingerprint of length `hidden` or shape
`(hidden × B)`.
"""
function (m::NDF)(Φ::AbstractMatrix, Â::AbstractMatrix;
                  seed_mask::Union{Nothing,AbstractVector{Bool}}=nothing)
    H_0 = transpose(m.encoder(transpose(Φ)))                       # (n × hidden)

    if seed_mask !== nothing
        H_0 = H_0 .* reshape(Float32.(seed_mask), :, 1)
    end

    Z = propagate(Â, H_0, m.K, m.α)

    f = if m.readout === :mean
        vec(mean(Z; dims=1))
    elseif m.readout === :sum
        vec(sum(Z; dims=1))
    elseif m.readout === :flatten
        vec(Z)
    else  # :seed_mean
        seed_mask === nothing &&
            throw(ArgumentError("seed_mask is required for readout=:seed_mean"))
        _masked_mean(Z, seed_mask)
    end

    return m.classifier(f)
end

function (m::NDF)(Φ::AbstractArray{T,3}, Â::AbstractMatrix;
                  seed_mask::Union{Nothing,AbstractMatrix{Bool}}=nothing) where T
    n, d_in, B = size(Φ)
    # Flatten (n × d_in × B) → (d_in × n*B) for the Dense encoder, then reshape
    # back to (n × hidden × B). permutedims is required because Dense operates
    # on dim 1; can't avoid the transpose without rewriting the MLP layer.
    Φ_for_dense = reshape(permutedims(Φ, (2, 1, 3)), d_in, n * B)
    H_flat = m.encoder(Φ_for_dense)                                # (hidden × n*B)
    h = size(H_flat, 1)
    H_0 = permutedims(reshape(H_flat, h, n, B), (2, 1, 3))         # (n × hidden × B)

    if seed_mask !== nothing
        H_0 = H_0 .* reshape(Float32.(seed_mask), n, 1, B)
    end

    Z = propagate(Â, H_0, m.K, m.α)

    f = if m.readout === :mean
        dropdims(mean(Z; dims=1); dims=1)                          # (hidden × B)
    elseif m.readout === :sum
        dropdims(sum(Z; dims=1); dims=1)
    elseif m.readout === :flatten
        reshape(Z, n * size(Z, 2), B)                              # (n*hidden × B)
    else  # :seed_mean
        seed_mask === nothing &&
            throw(ArgumentError("seed_mask is required for readout=:seed_mean"))
        _masked_mean(Z, seed_mask)
    end

    return m.classifier(f)
end

"""
    prepare_adjacency(g; self_loops=true) -> SparseMatrixCSC{Float32}

Build the symmetric-normalized adjacency Â = D^{-1/2}(A + I) D^{-1/2} used by
APPNP. Directed input graphs are symmetrized; diffusion needs an undirected
operator.
"""
function prepare_adjacency(g::AbstractGraph; self_loops::Bool=true)
    n = Int(nv(g))
    rows = Int[]
    cols = Int[]
    for e in edges(g)
        u, v = Int(src(e)), Int(dst(e))
        push!(rows, u); push!(cols, v)
        if is_directed(g) && u != v
            push!(rows, v); push!(cols, u)
        end
    end
    A = sparse(rows, cols, ones(Float32, length(rows)), n, n)
    return prepare_adjacency(A; self_loops=self_loops)
end

"""
    prepare_adjacency(A::SparseMatrixCSC; self_loops=true)

Same normalization as the graph version, but takes a precomputed adjacency
matrix directly. `A` must be square and symmetric (will be used as-is). Use
this when the adjacency is built outside of LightGraphs — for instance after
thresholding a PPMI matrix with `build_domain_graph`.
"""
function prepare_adjacency(A::SparseMatrixCSC; self_loops::Bool=true)
    n = size(A, 1)
    size(A, 2) == n || throw(ArgumentError("A must be square"))
    if self_loops
        A = A + sparse(1:n, 1:n, ones(Float32, n), n, n)
    end
    deg = max.(vec(sum(A; dims=2)), 1.0f0)
    D_inv_sqrt = spdiagm(0 => 1.0f0 ./ sqrt.(deg))
    return D_inv_sqrt * A * D_inv_sqrt
end

"""
    prepare_adjacency_directed_ppr(A::SparseMatrixCSC; self_loops=true)

Build the propagation operator for directed personalized PageRank: returns
`Â = A^T D⁻¹` where `D` is the diagonal out-degree matrix of `A`. With this
`Â`, the APPNP iteration `Z = (1-α) Â Z + α H_0` is equivalent to the
directed-PPR power iteration

    π[u] = α v[u] + (1-α) ∑_{w: w→u} π[w] / out_deg(w).

Use this when the domain graph is genuinely directed (e.g., the v1
collocation graph where u→v means v follows u). Dangling nodes (out-degree
0) are handled by adding self-loops when `self_loops=true`.
"""
function prepare_adjacency_directed_ppr(A::SparseMatrixCSC; self_loops::Bool=true)
    n = size(A, 1)
    size(A, 2) == n || throw(ArgumentError("A must be square"))
    if self_loops
        A = A + sparse(1:n, 1:n, ones(Float32, n), n, n)
    end
    out_deg = max.(vec(sum(A; dims=2)), 1.0f0)
    D_inv = spdiagm(0 => 1.0f0 ./ out_deg)
    return transpose(A) * D_inv
end

"""
    default_node_features(g) -> Matrix{Float32}

Minimal per-vertex feature matrix `(n × 4)`: normalized out-degree,
normalized in-degree, log-scaled out-degree, log-scaled in-degree. Undirected
graphs use the same value for in- and out-degree columns. Intended as a
sensible default to concatenate with the seed vector before passing to `NDF`.
"""
function default_node_features(g::AbstractGraph)
    n = Int(nv(g))
    od = Float32[outdegree(g, v) for v in vertices(g)]
    id = is_directed(g) ? Float32[indegree(g, v) for v in vertices(g)] : copy(od)
    return _features_from_degrees(od, id, n)
end

"""
    default_node_features(A::SparseMatrixCSC) -> Matrix{Float32}

Same feature set as the graph version, computed directly from a (symmetric)
adjacency matrix. Out- and in-degree columns are identical for symmetric
inputs.
"""
function default_node_features(A::SparseMatrixCSC)
    n = size(A, 1)
    od = Float32.(vec(sum(A; dims=2)))
    id = Float32.(vec(sum(A; dims=1)))
    return _features_from_degrees(od, id, n)
end

function _features_from_degrees(od::Vector{Float32}, id::Vector{Float32}, n::Int)
    max_od = max(maximum(od), 1.0f0)
    max_id = max(maximum(id), 1.0f0)
    log1p_n = log1p(Float32(n))
    return hcat(
        od ./ max_od,
        id ./ max_id,
        log1p.(od) ./ log1p_n,
        log1p.(id) ./ log1p_n,
    )
end
