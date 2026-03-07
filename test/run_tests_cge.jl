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

using Test
using Pkg
using Logging
Pkg.activate(normpath(joinpath(@__DIR__, "..")))

using Adjacently
using Adjacently.CGE
using Adjacently.IO: BitWriter, flush_bitwriter, load_adjacency_list_from_csv
using Adjacently.Clustering: louvain_partition, leiden_partition
using LightGraphs
import LightGraphs: nv, ne, is_directed, outneighbors

# Minimal graph type that subtypes LightGraphs.AbstractGraph{T}
struct TestGraph{T<:Unsigned} <: LightGraphs.AbstractGraph{T}
    n::Int
    adj::Dict{Int, Vector{T}}
    directed::Bool
end

function LightGraphs.nv(g::TestGraph)
    g.n
end
function LightGraphs.is_directed(g::TestGraph)
    g.directed
end
function LightGraphs.outneighbors(g::TestGraph{T}, v::Integer) where {T<:Unsigned}
    get(g.adj, Int(v), T[])
end

@testset "CGE basic encode_level and helpers" begin
    # Build a tiny undirected graph with 4 vertices and a few edges (UInt32 typed)
    T = UInt32
    undirected_adj = Dict{Int, Vector{T}}()
    function add_undirected!(u::Int, v::Int)
        push!(get!(undirected_adj, u, T[]), T(v))
        push!(get!(undirected_adj, v, T[]), T(u))
    end
    add_undirected!(1, 2)
    add_undirected!(2, 3)
    add_undirected!(3, 4)
    add_undirected!(1, 4)
    g = TestGraph{T}(4, undirected_adj, false)

    # Partition into two clusters
    P = [T[1, 2], T[3, 4]]

    # Encode one level
    io = IOBuffer()
    w = BitWriter(io)
    encode_level(w, g, P; params=CGEParams(L=128, varint=:fibonacci, count_varint=:elias_delta, gap=:fibonacci, degree=:elias_delta, perm_strategy=:blockpos, membership=:elias_fano))
    flush_bitwriter(w; flush_last_bits=true)
    bytes = take!(io)
    @test length(bytes) > 0

    # Test stub ordering -> permutation construction
    A_local = UInt32[1, 2]
    B_local = UInt32[10, 11]
    neighbors_in_B = Dict{UInt32,Vector{Int}}(
        UInt32(1) => [1, 2],
        UInt32(2) => Int[]
    )
    pi = CGE.order_from_edges(A_local, B_local, neighbors_in_B)
    @test pi == [1, 2]

    # Test permutation writers (lehmer and raw) produce some output
    io2 = IOBuffer(); w2 = BitWriter(io2)
    CGE.write_permutation(w2, [1,2,3]; strategy=:lehmer)
    flush_bitwriter(w2; flush_last_bits=true)
    @test length(take!(io2)) > 0

    io3 = IOBuffer(); w3 = BitWriter(io3)
    CGE.write_permutation(w3, [2,1]; strategy=:raw)
    flush_bitwriter(w3; flush_last_bits=true)
    @test length(take!(io3)) > 0
end

@testset "CGE encode_level on CNR-2000 with Louvain" begin
    prev_logger = current_logger()
    global_logger(ConsoleLogger(stderr, Logging.Info))
    try
    cnr_csv_path = "datasets/webgraph/cnr-2000/cnr-2000.csv"
    if !isfile(cnr_csv_path)
        @warn "CNR-2000 CSV not found at $cnr_csv_path; skipping CGE full-graph test"
        @test_skip "CNR-2000 dataset unavailable"
    else
        # Load full directed graph
        t0 = time(); g = load_adjacency_list_from_csv(cnr_csv_path, ',', true); t1 = time()
        n = nv(g)
        @info "Loaded CNR-2000: n=$(n) in $(round(t1-t0,digits=3))s"

        # Compute clustering (use moderate passes to keep test time reasonable)
        t2 = time(); part = louvain_partition(g; max_passes=3, max_levels=5); t3 = time()
        @info "Louvain communities: $(maximum(part)) in $(round(t3-t2,digits=3))s"
        @test length(part) == n

        # Collapse to at most K+1 clusters: top K by size, rest merged into last
        K = 7
        counts = Dict{Int,Int}()
        for c in part
            counts[c] = get(counts, c, 0) + 1
        end
        # sort labels by decreasing size
        labels_sorted = sort(collect(keys(counts)), by = c -> -counts[c])
        top = labels_sorted[1:min(K, length(labels_sorted))]
        top_index = Dict{Int,Int}(c => i for (i,c) in enumerate(top))

        # infer vertex type V from graph type parameter
        V = (typeof(g)).parameters[1]
        clusters = [Vector{V}() for _ in 1:(min(K, length(labels_sorted)) + 1)]
        other_label = length(clusters)  # last bucket for remaining communities
        for i in 1:n
            c = part[i]
            bucket = get(top_index, c, other_label)
            push!(clusters[bucket], convert(V, i))
        end

        # Ensure all clusters are non-empty (merging may leave empties if len<=K)
        clusters = filter(!isempty, clusters)

        # Encode CGE level
        io = IOBuffer(); w = BitWriter(io)
        params = CGEParams(L=128, varint=:fibonacci, count_varint=:elias_delta, gap=:fibonacci, degree=:elias_delta, undirected_pairs=false, perm_strategy=:blockpos, membership=:elias_fano)
        t4 = time(); encode_level(w, g, clusters; params=params); flush_bitwriter(w; flush_last_bits=true); bytes = take!(io); t5 = time()
        @info "CGE encode: size=$(length(bytes)) bytes, time=$(round(t5-t4,digits=3))s"
        @test length(bytes) > 0
    finally
        global_logger(prev_logger)
    end
end
