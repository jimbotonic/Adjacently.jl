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
    - `:opc_flatten`: v1 paper's OPC (Optimized PageRank Centrality) idea —
      flatten Z restricted to a precomputed top-`n_central` central node set.
      Classifier sees `n_central × hidden` features instead of `n_vertices ×
      hidden`. Massively shrinks the classifier head (e.g. 60K → 1K) and is
      the right move when the diffusion mass concentrates on a small "core"
      of central vertices. Caller passes the central indices as a kwarg
      (`central_nodes=`) at every forward call.
- `n_vertices`: number of vertices in the domain graph. Required when
      `readout=:flatten` so the classifier has the right input dim; ignored
      for `:opc_flatten` (which uses `n_central`) and for the pool readouts.
- `n_central`: size of the central-node set for `:opc_flatten`. Required
      when `readout=:opc_flatten`.
"""
struct NDF{E, C, G, W}
    encoder::E
    classifier::C
    K::Int
    α::Float32
    readout::Symbol
    propagation::Symbol      # :appnp (fixed) | :gpr (learnable γ_k, still LINEAR) | :nonlin (ReLU MP)
    gamma::G                 # length-(K+1) learnable per-hop weights (:gpr); else Float32[]
    Wprop::W                 # hidden×hidden learnable transform for :nonlin ReLU MP; else 0×0
end
Flux.@layer NDF trainable=(encoder, classifier, gamma, Wprop)

function NDF(d_in::Int, n_classes::Int;
             hidden::Int=64, K::Int=10, α::Float32=0.15f0,
             dropout::Float32=0.1f0, readout::Symbol=:mean,
             n_vertices::Int=0, n_central::Int=0,
             head::Symbol=:linear, head_hidden::Int=256,
             head_style::Symbol=:widemlp, propagation::Symbol=:appnp,
             n_gates::Int=4)
    readout ∈ (:mean, :sum, :seed_mean, :flatten, :opc_flatten, :flatten_seed_residual, :seed_multipool) ||
        throw(ArgumentError("readout must be :mean, :sum, :seed_mean, :flatten, :opc_flatten, :flatten_seed_residual, or :seed_multipool"))
    if readout ∈ (:flatten, :flatten_seed_residual) && n_vertices <= 0
        throw(ArgumentError("readout=$readout requires n_vertices > 0"))
    end
    if readout === :opc_flatten && n_central <= 0
        throw(ArgumentError("readout=:opc_flatten requires n_central > 0"))
    end
    head ∈ (:linear, :mlp) || throw(ArgumentError("head must be :linear or :mlp"))
    head_style ∈ (:widemlp, :pre_dropout) ||
        throw(ArgumentError("head_style must be :widemlp or :pre_dropout"))
    propagation ∈ (:appnp, :gpr, :nonlin, :gate) ||
        throw(ArgumentError("propagation must be :appnp, :gpr, :nonlin, or :gate"))
    if propagation === :gate && n_vertices <= 0
        throw(ArgumentError("propagation=:gate requires n_vertices > 0"))
    end
    # idea C: learnable per-hop weights γ_0..γ_K (GPR-GNN). Init to the APPNP
    # geometric weights α(1-α)^k so :gpr CONTAINS the fixed-α diffusion and
    # learns away from it only if that lowers loss.
    #
    # THE KEY (:gate): multiplicative seed-masked gate.
    #   Z[w] = v[w] · (1 + Σ_j g_j·σ(a_j·(Âv)[w] + b_j))
    # Reuse `gamma` to hold the 3r gate params [a₁..aᵣ ‖ b₁..bᵣ ‖ g₁..gᵣ]. Init
    # g=0 ⇒ Z=v EXACTLY (BoW floor is exact); the gate is a masked, absent-word-
    # vanishing, degree-2 graph correction that ESCAPES the linear-BoW ceiling.
    gamma = if propagation === :gpr
        Float32[α * (1f0 - α)^k for k in 0:K]
    elseif propagation === :gate
        vcat(fill(0.1f0, n_gates), zeros(Float32, n_gates), zeros(Float32, n_gates))
    else
        Float32[]
    end
    # NONLINEAR message passing: Z_{k+1}=ReLU((1-α)Â Z_k W + α H0). W init=I so
    # the first pass ≈ linear-APPNP-with-ReLU; ReLU breaks Z=Pv (escapes the
    # BoW-reparametrization ceiling). Needs hidden>1 + a pooling readout.
    Wprop = zeros(Float32, 0, 0)
    if propagation === :nonlin
        Wprop = zeros(Float32, hidden, hidden)
        @inbounds for i in 1:hidden; Wprop[i, i] = 1f0; end
    end
    encoder = NDFEncoder(d_in, hidden; dropout=dropout)
    classifier_in = if propagation === :gate
        n_vertices                     # gated seed Z is V-dim (BoW-shaped)
    elseif readout === :flatten
        n_vertices * hidden
    elseif readout === :flatten_seed_residual
        # flatten(Z_K) ⊕ raw seed vector — the diffused fingerprint plus the
        # verbatim BoW lexical channel WideMLP keeps. (idea B)
        n_vertices * hidden + n_vertices
    elseif readout === :seed_multipool
        # learned document aggregator: [tfidf-weighted-sum ‖ max ‖ mean] over the
        # doc's words — a NONLINEAR pool (max) that replaces near-linear seed_mean.
        3 * hidden
    elseif readout === :opc_flatten
        n_central * hidden
    else
        hidden
    end
    classifier = if n_classes <= 0
        identity
    elseif head === :mlp && head_style === :widemlp
        # WideMLP-EXACT head (Galke & Scherp): Dense→Dropout→Dense, NO input
        # dropout. Makes idea B a strict superset of WideMLP (head can ignore Z).
        Chain(Dense(classifier_in => head_hidden, relu),
              Dropout(dropout), Dense(head_hidden => n_classes))
    elseif head === :mlp
        # :pre_dropout — extra input dropout (heavier regularization variant).
        Chain(Dropout(dropout), Dense(classifier_in => head_hidden, relu),
              Dropout(dropout), Dense(head_hidden => n_classes))
    else
        Chain(Dropout(dropout), Dense(classifier_in => n_classes))
    end
    return NDF(encoder, classifier, K, α, readout, propagation, gamma, Wprop)
end

# Sparse-aware SpMV wrapper. Zygote's generic `*` adjoint for `A * X` computes
# `∂L/∂A = ΔY * X^T`, which for our 60K-vertex APPNP graph materializes a
# 60K × 60K dense (~14 GB) matrix — even though Â is not a parameter. The
# custom rrule below returns `NoTangent()` for A and only the cheap
# `A^T * ΔY` (sparse * dense) gradient for X. Drops peak memory from 17 GB
# to <1 GB on the modern-pipeline 60K-vertex graph.
_spmm(A::AbstractMatrix, X::AbstractMatrix) = A * X

# Zygote may hand this pullback a CPU Matrix cotangent even when the primal
# sparse multiply returned a CuArray. Materialize ΔY on the same array backend
# as the primal output Y before computing A' * ΔY — otherwise
# `CuSparseMatrix' * Matrix{CPU}` falls back to generic CPU matmul and
# scalar-indexes the CuSparse object (disallowed on GPU).
function _spmm_tangent_like(ΔY, Y::AbstractMatrix)
    ΔY = ChainRulesCore.unthunk(ΔY)
    ΔY isa AbstractMatrix || return ΔY
    typeof(ΔY) === typeof(Y) && return ΔY
    out = similar(Y, eltype(Y), size(ΔY))
    copyto!(out, ΔY)
    return out
end

function ChainRulesCore.rrule(::typeof(_spmm), A::AbstractMatrix, X::AbstractMatrix)
    Y = A * X
    function _spmm_pullback(ΔY)
        ΔY_dev = _spmm_tangent_like(ΔY, Y)
        return ChainRulesCore.NoTangent(), ChainRulesCore.NoTangent(), A' * ΔY_dev
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

"""
    propagate_gpr(Â, H_0, K, γ)

Generalized-PageRank propagation (idea C, GPR-GNN): `Z = Σ_{k=0}^{K} γ_k Â^k H_0`
with LEARNABLE per-hop weights `γ` (length K+1). Recovers APPNP when
`γ_k = α(1-α)^k`. `γ` is sliced as 1-element arrays (`γ[i:i]`) so the weighting
broadcasts on GPU without scalar indexing, and stays differentiable.
"""
function propagate_gpr(Â::AbstractMatrix, H_0::AbstractMatrix, K::Int, γ::AbstractVector)
    acc = reshape(γ[1:1], 1, 1) .* H_0
    Hk = H_0
    for k in 1:K
        Hk = _spmm(Â, Hk)
        acc = acc .+ reshape(γ[k+1:k+1], 1, 1) .* Hk
    end
    return acc
end

function propagate_gpr(Â::AbstractMatrix, H_0::AbstractArray{T,3}, K::Int, γ::AbstractVector) where T
    n, h, B = size(H_0)
    H_flat = reshape(H_0, n, h * B)
    acc = reshape(γ[1:1], 1, 1) .* H_flat
    Hk = H_flat
    for k in 1:K
        Hk = _spmm(Â, Hk)
        acc = acc .+ reshape(γ[k+1:k+1], 1, 1) .* Hk
    end
    return reshape(acc, n, h, B)
end

"""
    propagate_nonlin(Â, H_0, K, W, α)

NONLINEAR message passing: `Z_{k+1} = ReLU((1-α) Â Z_k W + α H_0)`, with a
learnable `hidden×hidden` transform `W` and a ReLU between hops. Unlike APPNP/GPR
(both linear ⇒ `Z=Pv`), this is a nonlinear function of the seed that requires
the graph — it escapes the BoW-reparametrization ceiling (AND-like adjacency-
gated conjunctions). `W` init = I ⇒ first pass ≈ linear-APPNP-with-ReLU.
"""
function propagate_nonlin(Â::AbstractMatrix, H_0::AbstractMatrix, K::Int, W::AbstractMatrix, α::Float32)
    oma = 1f0 - α
    Z = H_0
    for _ in 1:K
        Z = relu.(oma .* (_spmm(Â, Z) * W) .+ α .* H_0)
    end
    return Z
end

function propagate_nonlin(Â::AbstractMatrix, H_0::AbstractArray{T,3}, K::Int, W::AbstractMatrix, α::Float32) where T
    V, h, B = size(H_0)
    oma = 1f0 - α
    Z = H_0
    for _ in 1:K
        AZ = reshape(_spmm(Â, reshape(Z, V, h * B)), V, h, B)     # Â-message over vertices
        M = reshape(permutedims(AZ, (2, 1, 3)), h, V * B)         # (h × V*B)
        trans = permutedims(reshape(W' * M, h, V, B), (2, 1, 3))  # apply W to hidden dim
        Z = relu.(oma .* trans .+ α .* H_0)
    end
    return Z
end

"""
    _seed_gate(Â, v, γ)

THE KEY — multiplicative seed-masked gate:
    Z[w] = v[w] · (1 + Σ_{j=1..r} g_j · σ(a_j·(Âv)[w] + b_j))
with `γ = [a₁..aᵣ ‖ b₁..bᵣ ‖ g₁..gᵣ]` (learnable). `v` is the raw seed (V,) or
(V×B). Escapes the linear-BoW ceiling (degree-2 word×neighbor product ≠ P·v),
VANISHES on absent words (v[w]=0 ⇒ Z[w]=0, keeps the per-word BoW signal clean +
no capacity on absent dims), and is BoW-exact at init (g=0 ⇒ Z=v). γ sliced as
1-element arrays for GPU safety. Returns the same shape as `v`.
"""
function _seed_gate(Â::AbstractMatrix, v::AbstractVecOrMat, γ::AbstractVector)
    r = length(γ) ÷ 3
    vmat = v isa AbstractVector ? reshape(v, :, 1) : v          # (V×1) or (V×B)
    Av = _spmm(Â, vmat)                                          # neighbor context
    corr = zero(Av)
    for j in 1:r
        aj = reshape(γ[j:j], 1, 1); bj = reshape(γ[r+j:r+j], 1, 1); gj = reshape(γ[2r+j:2r+j], 1, 1)
        corr = corr .+ gj .* sigmoid.(aj .* Av .+ bj)
    end
    Zmat = vmat .* (1f0 .+ corr)
    return v isa AbstractVector ? vec(Zmat) : Zmat
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

# Learned document aggregator: concat [tfidf-weighted-sum ‖ max ‖ mean] over the
# doc's seed words. `max` is a NONLINEAR pool (cannot collapse to BoW); the
# weighted sum uses the tfidf seed weights. Replaces near-linear seed_mean.
function _seed_multipool(Z::AbstractMatrix, mask::AbstractVector{Bool}, seed_w::AbstractVector)
    mask_f = Float32.(mask)
    wsum = vec(sum(Z .* reshape(seed_w, :, 1); dims=1))
    mean_p = _masked_mean(Z, mask)
    Zneg = Z .+ reshape((1f0 .- mask_f) .* -1f9, :, 1)             # non-seed rows → -inf
    maxp = vec(maximum(Zneg; dims=1))
    return vcat(wsum, maxp, mean_p)
end

function _seed_multipool(Z::AbstractArray{T,3}, mask::AbstractMatrix{Bool}, seed_w::AbstractMatrix) where T
    V, h, B = size(Z)
    mask_f = Float32.(mask)
    wsum = dropdims(sum(Z .* reshape(seed_w, V, 1, B); dims=1); dims=1)     # (h × B)
    mean_p = _masked_mean(Z, mask)                                          # (h × B)
    Zneg = Z .+ reshape((1f0 .- mask_f) .* -1f9, V, 1, B)
    maxp = dropdims(maximum(Zneg; dims=1); dims=1)                          # (h × B)
    return vcat(wsum, maxp, mean_p)
end

"""
    (m::NDF)(Φ, Â; seed_mask=nothing, central_nodes=nothing)

Forward pass. `Φ` is the per-vertex feature input — `(n × d_in)` for a single
document or `(n × d_in × B)` for a batch. `seed_mask` is a Bool vector of
length n for single documents or a `(n × B)` Bool matrix for batches; when
provided, non-seed rows of H_0 are zeroed (preserving DF's "only seed
vertices teleport" semantics) and required when `readout == :seed_mean`.
`central_nodes` is a `Vector{Int}` of vertex indices to keep when
`readout == :opc_flatten` — its length must equal the NDF's `n_central`
constructor arg.

Returns logits `(n_classes,)` or `(n_classes × B)` when a classifier is
configured, otherwise the pooled fingerprint of length `hidden` or shape
`(hidden × B)`.
"""
function (m::NDF)(Φ::AbstractMatrix, Â::AbstractMatrix;
                  seed_mask::Union{Nothing,AbstractVector{Bool}}=nothing,
                  central_nodes::Union{Nothing,AbstractVector{<:Integer}}=nothing)
    if m.propagation === :gate     # THE KEY: multiplicative seed-masked gate (V-dim, BoW-exact at init)
        return m.classifier(_seed_gate(Â, Φ[:, 1], m.gamma))
    end
    H_0 = transpose(m.encoder(transpose(Φ)))                       # (n × hidden)

    if seed_mask !== nothing
        H_0 = H_0 .* reshape(Float32.(seed_mask), :, 1)
    end

    Z = m.propagation === :gpr    ? propagate_gpr(Â, H_0, m.K, m.gamma) :
        m.propagation === :nonlin ? propagate_nonlin(Â, H_0, m.K, m.Wprop, m.α) :
                                    propagate(Â, H_0, m.K, m.α)

    f = if m.readout === :mean
        vec(mean(Z; dims=1))
    elseif m.readout === :sum
        vec(sum(Z; dims=1))
    elseif m.readout === :flatten
        vec(Z)
    elseif m.readout === :flatten_seed_residual
        vcat(vec(Z), Φ[:, 1])                                      # (n*hidden + n,)
    elseif m.readout === :seed_multipool
        seed_mask === nothing &&
            throw(ArgumentError("seed_mask is required for readout=:seed_multipool"))
        _seed_multipool(Z, seed_mask, Φ[:, 1])                     # (3*hidden,)
    elseif m.readout === :opc_flatten
        central_nodes === nothing &&
            throw(ArgumentError("central_nodes required for readout=:opc_flatten"))
        vec(Z[central_nodes, :])                                   # (n_central * hidden,)
    else  # :seed_mean
        seed_mask === nothing &&
            throw(ArgumentError("seed_mask is required for readout=:seed_mean"))
        _masked_mean(Z, seed_mask)
    end

    return m.classifier(f)
end

function (m::NDF)(Φ::AbstractArray{T,3}, Â::AbstractMatrix;
                  seed_mask::Union{Nothing,AbstractMatrix{Bool}}=nothing,
                  central_nodes::Union{Nothing,AbstractVector{<:Integer}}=nothing) where T
    if m.propagation === :gate     # THE KEY: multiplicative seed-masked gate (V-dim, BoW-exact at init)
        return m.classifier(_seed_gate(Â, Φ[:, 1, :], m.gamma))
    end
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

    Z = m.propagation === :gpr    ? propagate_gpr(Â, H_0, m.K, m.gamma) :
        m.propagation === :nonlin ? propagate_nonlin(Â, H_0, m.K, m.Wprop, m.α) :
                                    propagate(Â, H_0, m.K, m.α)

    f = if m.readout === :mean
        dropdims(mean(Z; dims=1); dims=1)                          # (hidden × B)
    elseif m.readout === :sum
        dropdims(sum(Z; dims=1); dims=1)
    elseif m.readout === :flatten
        reshape(Z, n * size(Z, 2), B)                              # (n*hidden × B)
    elseif m.readout === :flatten_seed_residual
        vcat(reshape(Z, n * size(Z, 2), B), Φ[:, 1, :])           # (n*hidden + n × B)
    elseif m.readout === :seed_multipool
        seed_mask === nothing &&
            throw(ArgumentError("seed_mask is required for readout=:seed_multipool"))
        _seed_multipool(Z, seed_mask, Φ[:, 1, :])                 # (3*hidden × B)
    elseif m.readout === :opc_flatten
        central_nodes === nothing &&
            throw(ArgumentError("central_nodes required for readout=:opc_flatten"))
        # Z[central_nodes, :, :] → (n_central × hidden × B); flatten to
        # (n_central*hidden × B). The gather of central_nodes is a non-
        # differentiable index op which Zygote handles fine.
        Z_c = Z[central_nodes, :, :]
        reshape(Z_c, length(central_nodes) * size(Z_c, 2), B)
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
    top_central_nodes(Â, d; K=50, α=0.15f0) -> Vector{Int}

Return the indices of the top-`d` vertices ranked by uniform-personalization
PageRank centrality on the normalized propagation operator `Â`. Used to
select the central-node set for `readout=:opc_flatten` (the v1 paper's OPC
construction).

Implementation: run the APPNP iteration `π = (1-α) Â π + α (1/n)·1` for
`K` steps from a uniform initial vector. With `K=50, α=0.15` this converges
to within float precision on graphs we've tested up to V=60K (l2 error
≤ 1e-6 vs K=200). Falls back gracefully on GPU `Â` — the centrality result
is moved back to CPU before `partialsortperm`.
"""
function top_central_nodes(Â::AbstractMatrix, d::Int; K::Int=50,
                           α::Float32=0.15f0)
    n = size(Â, 1)
    d > 0 || throw(ArgumentError("d must be positive"))
    d <= n || throw(ArgumentError("d ($d) exceeds graph size ($n)"))
    # Called on a CPU Â — caller is responsible for passing a CPU version
    # when on GPU, since centrality only needs the graph topology (no
    # gradients, no batching) and the sparse matmul falls back to a slow
    # CPU/GPU mixed path otherwise.
    v = fill(1.0f0 / n, n, 1)
    π = propagate(Â, v, K, α)
    π_cpu = Array(π)
    return partialsortperm(vec(π_cpu), 1:d; rev=true)
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
