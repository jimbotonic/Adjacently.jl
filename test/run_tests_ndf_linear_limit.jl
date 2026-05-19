#!/usr/bin/env julia
#
# Test: NDF in the linear, no-features, K→∞ limit reproduces closed-form
# personalized PageRank.
#
# This is the anchor experiment for the v2 fingerprints paper claim:
# "Diffusion Fingerprints are the linear, no-features, K→∞ limit of NDF."
#
# Construction tested:
#   d_in    = 1   (seed vector only, no structural features)
#   hidden  = 1   (no expansion)
#   encoder = identity (Dense weights set to 1, biases to 0; ReLU is identity
#             on the non-negative seed vector)
#   K       = 200 (large enough for truncated APPNP to be ε-close to fixed point)
#
# Closed form: π_k = α(I - (1-α)Â)^{-1} v_k
#
# Usage:
#   julia --project test/run_tests_ndf_linear_limit.jl

include("run_tests_main.jl")

using LinearAlgebra: norm, I as LinIdentity
using Random: MersenneTwister
using Flux
using SparseArrays: SparseMatrixCSC, nnz

using Adjacently.Fingerprints: NDF, NDFEncoder, propagate, prepare_adjacency

@info "NDF linear-limit test"

# Tiny but non-trivial directed graph with two loose communities.
function build_test_graph(n_per::Int=20, p_intra::Float64=0.25,
                          p_inter::Float64=0.02; seed::Int=0)
    rng = MersenneTwister(seed)
    n = 2 * n_per
    g = SimpleDiGraph(n)
    for u in 1:n, v in 1:n
        u == v && continue
        same = (u <= n_per) == (v <= n_per)
        p = same ? p_intra : p_inter
        rand(rng) < p && add_edge!(g, u, v)
    end
    return g
end

g = build_test_graph()
n = nv(g)
Â = prepare_adjacency(g)
@info "graph: n=$n  m=$(ne(g))  Â nnz=$(nnz(Â))"

# Closed-form PPR solver. For a Float32 dense Â of modest n, the inverse is
# fine; we only use this as the ground-truth oracle for the test.
function closed_form_ppr(Â::SparseMatrixCSC{Float32}, v::Vector{Float32}, α::Float32)
    n = size(Â, 1)
    M = LinIdentity(n) - (1.0f0 - α) .* Matrix(Â)
    return α .* (M \ v)
end

# Build a linear-limit NDF: weights manually set to identity, no dropout.
function linear_limit_ndf(α::Float32, K::Int)
    encoder = NDFEncoder(1, 1; dropout=0.0f0)
    # encoder = Chain(Dense(1=>1, relu), Dropout(0), Dense(1=>1))
    # Set weights/biases to identity. Both Dense layers are 1×1 matrices.
    encoder[1].weight .= 1.0f0; encoder[1].bias .= 0.0f0
    encoder[3].weight .= 1.0f0; encoder[3].bias .= 0.0f0
    classifier = identity
    return NDF(encoder, classifier, K, α, :flatten)
end

α = 0.15f0
K = 200
model = linear_limit_ndf(α, K)
Flux.testmode!(model)

# Several seed sets to test, including single-seed (the standard PPR setup),
# multi-seed sets within / spanning the planted communities, and an isolated
# seed at the boundary.
seed_sets = [
    [1],
    [3, 7, 11],
    [21, 25, 29, 33],
    [5, 25],
    [1, 2, 3, 19, 20],
]

@info "Comparing NDF(K=$K) vs closed-form PPR for $(length(seed_sets)) seed sets"
global max_diff_over_sets = 0.0f0
@testset "NDF linear limit equals closed-form PPR" begin
    for (i, idx) in enumerate(seed_sets)
        v_k = zeros(Float32, n)
        v_k[idx] .= 1.0f0 / length(idx)
        Φ = reshape(v_k, n, 1)

        # NDF forward (per-vertex output via :flatten readout, no classifier).
        z_ndf = model(Φ, Â)                           # (n,) vector
        # Closed-form PPR.
        z_exact = closed_form_ppr(Â, v_k, α)

        max_abs = maximum(abs.(z_ndf .- z_exact))
        rel = norm(z_ndf .- z_exact) / norm(z_exact)
        global max_diff_over_sets = max(max_diff_over_sets, max_abs)
        @info "seed_set $i (|T'|=$(length(idx))): max_abs=$(round(max_abs; sigdigits=3)) rel_l2=$(round(rel; sigdigits=3))"

        @test max_abs < 1.0f-4
        @test rel < 1.0f-3
    end
end

# Convergence rate: APPNP at K iterations has residual O((1-α)^K). For α=0.15,
# (1-α)^K = 0.85^K. The residual should drop by ~0.85 per iteration on
# average; verify the residual decreases monotonically.
@info "Convergence rate sweep"
v_k = zeros(Float32, n); v_k[[3, 7, 11]] .= 1.0f0 / 3
Φ = reshape(v_k, n, 1)
exact = closed_form_ppr(Â, v_k, α)

@testset "K-sweep monotone convergence" begin
    last_err = Inf32
    for K_try in [1, 5, 10, 25, 50, 100, 200]
        m_try = linear_limit_ndf(α, K_try)
        Flux.testmode!(m_try)
        z = m_try(Φ, Â)
        err = norm(z .- exact)
        @info "  K=$K_try  l2_err=$(round(err; sigdigits=3))"
        @test err <= last_err + 1.0f-7
        last_err = err
    end
    @test last_err < 1.0f-5
end

@info "Linear-limit test complete (max abs diff across all seed sets: $(round(max_diff_over_sets; sigdigits=3)))"
