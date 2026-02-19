#!/usr/bin/env julia

#
# Test if _stochastic_matrix is equivalent to transpose(get_sparse_P_matrix)
#

include("run_tests_main.jl")

using LightGraphs: nv, add_edge!, add_vertices!, outneighbors, AbstractGraph
using SparseArrays: SparseMatrixCSC, sparse
using LinearAlgebra: norm
using Adjacently
using Adjacently.CustomLightGraphs: SimpleDiGraph
using Adjacently.Graph: get_sparse_P_matrix

println("=" ^ 80)
println("Testing Matrix Equivalence")
println("=" ^ 80)

# Build _stochastic_matrix (from astra_layered.jl)
function _stochastic_matrix(g::AbstractGraph{T}) where {T<:Unsigned}
    n = nv(g)
    I = Int[]
    J = Int[]
    V = Float64[]
    for j in 1:n
        jv = T(j)
        outs = outneighbors(g, jv)
        d = length(outs)
        if d == 0
            continue
        end
        w = 1.0 / d
        for i in outs
            push!(I, Int(i))
            push!(J, j)
            push!(V, w)
        end
    end
    Aint = sparse(I, J, V, n, n)
    return SparseMatrixCSC{Float64,T}(Aint)
end

# Create test graph (no dangling nodes - all vertices have outgoing edges)
g = SimpleDiGraph{UInt24}()
add_vertices!(g, 4)
add_edge!(g, UInt24(1), UInt24(2))
add_edge!(g, UInt24(1), UInt24(3))
add_edge!(g, UInt24(2), UInt24(3))
add_edge!(g, UInt24(3), UInt24(1))
add_edge!(g, UInt24(4), UInt24(1))

println("\nTest graph: 4 vertices (no dangling nodes)")
println("Edges: 1→2, 1→3, 2→3, 3→1, 4→1")

# Method 1: _stochastic_matrix (column-stochastic)
P_stochastic = _stochastic_matrix(g)

# Method 2: transpose(get_sparse_P_matrix) (should be column-stochastic)
P_existing = get_sparse_P_matrix(g)
P_existing_transpose = sparse(P_existing')

println("\n_stochastic_matrix (column-stochastic):")
println(Matrix(P_stochastic))

println("\nget_sparse_P_matrix (row-stochastic):")
println(Matrix(P_existing))

println("\ntranspose(get_sparse_P_matrix) (should be column-stochastic):")
println(Matrix(P_existing_transpose))

# Compare
diff = norm(Matrix(P_stochastic) - Matrix(P_existing_transpose), Inf)
println("\nMax difference: $diff")

if diff < 1e-10
    println("✓ IDENTICAL! Can use transpose(get_sparse_P_matrix(g)) instead of _stochastic_matrix(g)")
else
    println("✗ DIFFERENT!")
end

# Verify stochasticity
println("\nColumn sums of _stochastic_matrix:")
col_sums = vec(sum(P_stochastic, dims=1))
println(col_sums)
println("All close to 1.0? $(all(abs.(col_sums .- 1.0) .< 0.01))")

println("\nColumn sums of transpose(get_sparse_P_matrix):")
col_sums_existing = vec(sum(P_existing_transpose, dims=1))
println(col_sums_existing)
println("All close to 1.0? $(all(abs.(col_sums_existing .- 1.0) .< 0.01))")

println("\n" ^ 2 * "=" ^ 80)
println("RECOMMENDATION")
println("=" ^ 80)

if diff < 1e-10
    println("YES! The _stochastic_matrix function is redundant.")
    println()
    println("Replace:")
    println("  P = _stochastic_matrix(g)")
    println()
    println("With:")
    println("  P = sparse(get_sparse_P_matrix(g)')")
    println()
    println("Or define a helper:")
    println("  get_column_stochastic_matrix(g) = sparse(get_sparse_P_matrix(g)')")
    println()
    println("This:")
    println("  • Reduces code duplication")
    println("  • Uses existing, tested functions")
    println("  • Makes relationship between row/column stochastic clear")
else
    println("NO - The functions produce different results.")
    println("Keep _stochastic_matrix as a separate implementation.")
end

println("=" ^ 80)
