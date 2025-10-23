#
# Basic tests for AMG (Smoothed Aggregation) routines
#

using Test
using SparseArrays, LinearAlgebra
using Adjacently
using Adjacently.Clustering
using Adjacently.CustomLightGraphs: SimpleDiGraph
using Adjacently.Graph
using LightGraphs: add_vertices!, add_edge!, nv

function path_laplacian(n::Int)
    I = Int[]; J = Int[]; V = Float64[]
    for i in 1:n
        di = 0
        if i > 1
            push!(I, i); push!(J, i-1); push!(V, -1.0)
            di += 1
        end
        if i < n
            push!(I, i); push!(J, i+1); push!(V, -1.0)
            di += 1
        end
        push!(I, i); push!(J, i); push!(V, di)
    end
    return sparse(I, J, V, n, n)
end

@testset "AMG Setup and V-cycle (UInt24 indices)" begin
    n = 64
    A = path_laplacian(n)
    # Convert to custom Unsigned index type (UInt24)
    A = SparseMatrixCSC{Float64,Adjacently.CustomTypes.UInt24}(A)

    # Zero-sum RHS to avoid Laplacian nullspace (L*1 = 0)
    b = zeros(n); b[1] = 1.0; b[end] = -1.0
    x0 = zeros(n)
    H = amg_setup(A; theta=0.08, max_levels=10, min_coarse=8, smoothed=true)

    r0 = norm(b - A * x0)
    x = copy(x0)
    vcycle!(H, b, x; nu_pre=1, nu_post=1)
    r1 = norm(b - A * x)
    @test r1 < r0

    for _ in 1:3
        vcycle!(H, b, x; nu_pre=1, nu_post=1)
    end
    r2 = norm(b - A * x)
    @test r2 < r1
end

@testset "AMG Setup and V-cycle" begin
    n = 64
    A = path_laplacian(n)
    # Convert to Unsigned index type (e.g., UInt32)
    A = SparseMatrixCSC{Float64,UInt32}(A)
    # Zero-sum RHS to avoid Laplacian nullspace (L*1 = 0)
    b = zeros(n); b[1] = 1.0; b[end] = -1.0
    x0 = zeros(n)
    H = amg_setup(A; theta=0.08, max_levels=10, min_coarse=8, smoothed=true)
    # One V-cycle should reduce residual
    r0 = norm(b - A * x0)
    x = copy(x0)
    vcycle!(H, b, x; nu_pre=1, nu_post=1)
    r1 = norm(b - A * x)
    @test r1 < r0
    # A few cycles should further reduce residual
    for _ in 1:3
        vcycle!(H, b, x; nu_pre=1, nu_post=1)
    end
    r2 = norm(b - A * x)
    @test r2 < r1
end

@testset "AMG from directed graph adjacency" begin
    # Build a small directed graph with reciprocal edges (undirected structure)
    n = 64
    T = Adjacently.CustomTypes.UInt24
    g = SimpleDiGraph{T}()
    add_vertices!(g, n)
    for i in 1:n-1
        add_edge!(g, convert(T, i), convert(T, i+1))
        add_edge!(g, convert(T, i+1), convert(T, i))
    end

    # Export adjacency and form Laplacian-like M = D - A (indices type T)
    A = Graph.get_sparse_adj_matrix(g)
    @test isa(A, SparseMatrixCSC{Float64, T})
    S = dropdims(sum(A, dims=2), dims=2)
    range = convert(Vector{T}, collect(1:n))
    D = sparse(range, range, S)
    M = D - A

    # Zero-sum RHS to avoid nullspace; initial guess zero
    b = zeros(n); b[1] = 1.0; b[end] = -1.0
    x0 = zeros(n)

    H = amg_setup(M; theta=0.08, max_levels=10, min_coarse=8, smoothed=true)
    r0 = norm(b - M * x0)
    x = copy(x0)
    vcycle!(H, b, x; nu_pre=1, nu_post=1)
    r1 = norm(b - M * x)
    @test r1 < r0
    for _ in 1:3
        vcycle!(H, b, x; nu_pre=1, nu_post=1)
    end
    r2 = norm(b - M * x)
    @test r2 < r1
end
