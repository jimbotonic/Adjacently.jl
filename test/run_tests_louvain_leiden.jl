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

include("run_tests_main.jl")


#
# Tests for Louvain and Leiden community detection
#

using Test
using Logging
using Adjacently.MGS: load_compressed_mgs3_graph
using Adjacently
using Adjacently.CustomLightGraphs: SimpleDiGraph, SimpleGraph
using Adjacently.Clustering
using Adjacently.IO: load_adjacency_list_from_csv
using LightGraphs: add_vertices!, add_edge!, nv, ne
using Random

@testset "Louvain/Leiden on two-clique graph (UInt24)" begin
    T = Adjacently.CustomTypes.UInt24
    n = 12
    g = SimpleDiGraph{T}()
    add_vertices!(g, n)
    # two cliques 1..6 and 7..12, reciprocal edges
    for i in 1:6, j in i+1:6
        add_edge!(g, convert(T,i), convert(T,j))
        add_edge!(g, convert(T,j), convert(T,i))
    end
    for i in 7:12, j in i+1:12
        add_edge!(g, convert(T,i), convert(T,j))
        add_edge!(g, convert(T,j), convert(T,i))
    end
    # sparse cross edges
    add_edge!(g, convert(T,3), convert(T,8)); add_edge!(g, convert(T,8), convert(T,3))
    add_edge!(g, convert(T,5), convert(T,10)); add_edge!(g, convert(T,10), convert(T,5))

    P_louvain = louvain_partition(g)
    @test length(P_louvain) == n
    # Expect 2 communities
    @test maximum(P_louvain) == 2
    # First 6 likely same community
    c1 = P_louvain[1]
    @test all(P_louvain[i] == c1 for i in 1:6)

    P_leiden = leiden_partition(g)
    @test length(P_leiden) == n
    @test maximum(P_leiden) == 2
    c1l = P_leiden[1]
    @test all(P_leiden[i] == c1l for i in 1:6)
end

@testset "Louvain/Leiden on full CNR-2000" begin
    prev_logger = current_logger()
    global_logger(ConsoleLogger(stderr, Logging.Info))
    try
    cnr_csv_path = "datasets/webgraph/cnr-2000/cnr-2000.mgz"
    if !isfile(cnr_csv_path)
        @warn "CNR-2000 CSV not found at $cnr_csv_path; skipping full-graph clustering test"
        @test_skip "CNR-2000 dataset unavailable"
    else
        # Load directed webgraph
        t0 = time(); g = load_compressed_mgs3_graph(cnr_csv_path); t1 = time()
        n = nv(g)
        @info "CNR-2000 loaded: n=$(n), m=$(ne(g)) in $(round(t1-t0,digits=3))s"

        # Run Louvain/Leiden with moderate limits for runtime
        t2 = time(); P_l = louvain_partition(g; max_passes=3, max_levels=5); t3 = time()
        @info "Louvain produced $(maximum(P_l)) communities in $(round(t3-t2,digits=3))s"
        t4 = time(); P_d = leiden_partition(g; max_passes=3, max_levels=5); t5 = time()
        @info "Leiden produced $(maximum(P_d)) communities in $(round(t5-t4,digits=3))s"

        @test length(P_l) == n && length(P_d) == n

        # Sanity: communities are within 1..k and nontrivial
        @test maximum(P_l) >= 2 && maximum(P_d) >= 2

        # Quality: Leiden modularity should not be worse than Louvain by more than a small tolerance
        t6 = time(); Ql = Adjacently.Clustering.modularity(g, P_l); t7 = time()
        t8 = time(); Qd = Adjacently.Clustering.modularity(g, P_d); t9 = time()
        @info "Modularity: Louvain=$(round(Ql, digits=4)) in $(round(t7-t6,digits=3))s; Leiden=$(round(Qd, digits=4)) in $(round(t9-t8,digits=3))s"
        @test Qd + 1e-3 >= Ql - 1e-3
    end
    finally
        global_logger(prev_logger)
    end
end

@testset "Louvain/Leiden on random directed graph (UInt24)" begin
    T = Adjacently.CustomTypes.UInt24
    n = 1000
    g = SimpleDiGraph{T}()
    add_vertices!(g, n)
    rng = MersenneTwister(24)
    p = 0.0015
    for i in 1:n
        for j in 1:n
            i == j && continue
            if rand(rng) < p
                add_edge!(g, convert(T,i), convert(T,j))
            end
        end
    end
    P_l = louvain_partition(g; max_passes=5, max_levels=5)
    P_d = leiden_partition(g; max_passes=5, max_levels=5)
    @test length(P_l) == n && length(P_d) == n
    Ql = Adjacently.Clustering.modularity(g, P_l)
    Qd = Adjacently.Clustering.modularity(g, P_d)
    # Leiden may split communities conservatively on directed random graphs;
    # allow small tolerance against Louvain's modularity.
    @test Qd >= Ql - 1e-2
end

@testset "Louvain/Leiden on undirected two-clique graph (UInt24)" begin
    T = Adjacently.CustomTypes.UInt24
    n = 12
    g = SimpleDiGraph{T}()
    add_vertices!(g, n)
    # two cliques 1..6 and 7..12 (reciprocal edges to emulate undirected)
    for i in 1:6, j in i+1:6
        add_edge!(g, convert(T,i), convert(T,j))
        add_edge!(g, convert(T,j), convert(T,i))
    end
    for i in 7:12, j in i+1:12
        add_edge!(g, convert(T,i), convert(T,j))
        add_edge!(g, convert(T,j), convert(T,i))
    end
    # sparse cross edges (reciprocal)
    add_edge!(g, convert(T,3), convert(T,8))
    add_edge!(g, convert(T,8), convert(T,3))
    add_edge!(g, convert(T,5), convert(T,10))
    add_edge!(g, convert(T,10), convert(T,5))

    P_louvain = louvain_partition(g)
    @test length(P_louvain) == n
    @test maximum(P_louvain) == 2
    c1 = P_louvain[1]
    @test all(P_louvain[i] == c1 for i in 1:6)

    P_leiden = leiden_partition(g)
    @test length(P_leiden) == n
    @test maximum(P_leiden) == 2
    c1l = P_leiden[1]
    @test all(P_leiden[i] == c1l for i in 1:6)
end

@testset "Louvain/Leiden on random undirected graph (UInt24)" begin
    T = Adjacently.CustomTypes.UInt24
    n = 1000
    g = SimpleDiGraph{T}()
    add_vertices!(g, n)
    rng = MersenneTwister(42)
    p = 0.002
    for i in 1:n-1
        for j in i+1:n
            if rand(rng) < p
                add_edge!(g, convert(T,i), convert(T,j))
                add_edge!(g, convert(T,j), convert(T,i))
            end
        end
    end
    P_l = louvain_partition(g; max_passes=5, max_levels=5)
    P_d = leiden_partition(g; max_passes=5, max_levels=5)
    @test length(P_l) == n && length(P_d) == n
    # Compare modularity quality (Leiden should be no worse)
    Ql = Adjacently.Clustering.modularity(g, P_l)
    Qd = Adjacently.Clustering.modularity(g, P_d)
    @test Qd + 1e-8 >= Ql
end
