#
# Adjacently: Julia Complex Directed Networks Library
# Copyright (C) 2016-2025 Jimmy Dubuisson <jimmy.dubuisson@gmail.com>
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

include("run_tests_main.jl")

#
# Basic tests for Clustering (METIS-like) routines
#

using Test
using Adjacently
using Adjacently.CustomLightGraphs: SimpleGraph, SimpleDiGraph
using Adjacently.Clustering
using LightGraphs: add_edge!, add_vertex!, add_vertices!, nv

@testset "Clustering-METIS UInt24" begin
    # Build a small graph with two dense halves and a light cut
    n = 12
    g = SimpleDiGraph{Adjacently.CustomTypes.UInt24}()
    add_vertices!(g, n)

    # Bidirectional edges to emulate undirected cliques
    # Clique on 1..6
    for i in 1:6, j in i+1:6
        add_edge!(g, convert(Adjacently.CustomTypes.UInt24, i), convert(Adjacently.CustomTypes.UInt24, j))
        add_edge!(g, convert(Adjacently.CustomTypes.UInt24, j), convert(Adjacently.CustomTypes.UInt24, i))
    end
    # Clique on 7..12
    for i in 7:12, j in i+1:12
        add_edge!(g, convert(Adjacently.CustomTypes.UInt24, i), convert(Adjacently.CustomTypes.UInt24, j))
        add_edge!(g, convert(Adjacently.CustomTypes.UInt24, j), convert(Adjacently.CustomTypes.UInt24, i))
    end
    # Few cross edges (both directions)
    add_edge!(g, convert(Adjacently.CustomTypes.UInt24, 3), convert(Adjacently.CustomTypes.UInt24, 8))
    add_edge!(g, convert(Adjacently.CustomTypes.UInt24, 8), convert(Adjacently.CustomTypes.UInt24, 3))
    add_edge!(g, convert(Adjacently.CustomTypes.UInt24, 5), convert(Adjacently.CustomTypes.UInt24, 10))
    add_edge!(g, convert(Adjacently.CustomTypes.UInt24, 10), convert(Adjacently.CustomTypes.UInt24, 5))

    @test nv(g) == n

    # Heavy-edge matching
    matches = heavy_edge_matching(g)
    covered = Set{Int}()
    for (u,v) in matches
        push!(covered, Int(u))
        push!(covered, Int(v))
    end
    @test length(covered) == n

    # Contraction
    gc, wdict, coarse_of = contract_graph(g, matches)
    @test maximum(coarse_of) == nv(gc)
    @test nv(gc) == length(matches)
    @test all(1 .<= coarse_of .<= nv(gc))

    # Initial partition + refinement
    k = 2
    P0 = initial_kway_partition(g, k)
    @test length(P0) == n
    @test all(1 .<= P0 .<= k)

    P1 = copy(P0)
    refine_fm!(g, P1)
    @test length(P1) == n
    @test all(1 .<= P1 .<= k)

    # End-to-end
    P = metis_partition(g, k)
    @test length(P) == n
    @test all(1 .<= P .<= k)
end

@testset "Clustering-METIS UInt40 small" begin
    n = 6
    g = SimpleDiGraph{Adjacently.CustomTypes.UInt40}()
    add_vertices!(g, n)
    for i in 1:n-1
        add_edge!(g, convert(Adjacently.CustomTypes.UInt40, i), convert(Adjacently.CustomTypes.UInt40, i+1))
        add_edge!(g, convert(Adjacently.CustomTypes.UInt40, i+1), convert(Adjacently.CustomTypes.UInt40, i))
    end
    @test nv(g) == n
    P = metis_partition(g, 2; min_coarse_size=2)  # keep small to exercise pipeline
    @test length(P) == n
    @test all(1 .<= P .<= 2)
end

@testset "Clustering-METIS Directed UInt24" begin
    # Two strongly connected groups with sparse inter-links
    n = 10
    T = Adjacently.CustomTypes.UInt24
    g = SimpleDiGraph{T}()
    add_vertices!(g, n)

    # Strongly connect 1..5 bidirectionally
    for i in 1:5, j in 1:5
        i == j && continue
        add_edge!(g, convert(T, i), convert(T, j))
    end
    # Strongly connect 6..10 bidirectionally
    for i in 6:10, j in 6:10
        i == j && continue
        add_edge!(g, convert(T, i), convert(T, j))
    end
    # A couple of directed cross edges
    add_edge!(g, convert(T, 3), convert(T, 8))
    add_edge!(g, convert(T, 9), convert(T, 2))

    @test nv(g) == n

    # Matching should cover all nodes
    matches = heavy_edge_matching(g)
    covered = Set{Int}()
    for (u,v) in matches
        push!(covered, Int(u))
        push!(covered, Int(v))
    end
    @test length(covered) == n

    # Contraction should preserve directedness
    gc, wdict, coarse_of = contract_graph(g, matches)
    @test maximum(coarse_of) == nv(gc)
    @test nv(gc) == length(matches)
    @test all(1 .<= coarse_of .<= nv(gc))
    @test gc isa SimpleDiGraph{T}

    # End-to-end partitioning
    k = 2
    P = metis_partition(g, k)
    @test length(P) == n
    @test all(1 .<= P .<= k)
end
