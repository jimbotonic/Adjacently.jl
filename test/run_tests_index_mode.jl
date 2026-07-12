#
# Adjacently: Julia Complex Directed Networks Library
# Copyright (C) 2016-2026 Anonymous (double-blind review)
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

#
# Test per-vertex index mode roundtrip for BG, CS, and CG.
# Verifies that encoding with coding_scheme=:index and decoding produces
# identical edge lists to the children mode encoding.
#
# Usage:
#   julia --project test/run_tests_index_mode.jl [DATASET]
#   DATASET defaults to "cnr-2000"

include("run_tests_main.jl")

using Logging
global_logger(ConsoleLogger(stderr, Logging.Info))

using LightGraphs: nv, ne, outneighbors, vertices, eltype
using Adjacently
using Adjacently.IO: load_adjacency_list_from_csv
using Adjacently.Clustering: leiden_partition
using Adjacently.MGS: write_bg_mgs3_graph, write_cs_mgs3_graph, write_cg_mgs3_graph, load_compressed_mgs3_graph
using Adjacently.Compression.CG: CGParams

const DATASET = length(ARGS) >= 1 ? ARGS[1] : "cnr-2000"
const PREFIX  = replace(DATASET, "-" => "")
const DS_DIR  = normpath(joinpath(@__DIR__, "..", "datasets", "webgraph", DATASET))
const DS_CSV  = joinpath(DS_DIR, DATASET * ".csv")

function extract_neighbor_lists(g)
    T = eltype(g)
    nls = Dict{T,Vector{T}}()
    for v in vertices(g)
        nls[T(v)] = sort([T(o) for o in outneighbors(g, v)])
    end
    return nls
end

function verify_roundtrip(original_nls, decoded_nls, m)
    dec_edges = sum(length(v) for v in values(decoded_nls))
    dec_edges != m && return false
    for (v, orig) in original_nls
        haskey(decoded_nls, v) || return false
        sort(decoded_nls[v]) == sort(orig) || return false
    end
    return true
end

if !isfile(DS_CSV)
    @warn "$DATASET CSV not found at $DS_CSV; skipping"
    @testset "Index Mode - $DATASET" begin
        @test_skip "$DATASET dataset unavailable"
    end
else

@info "Loading $DATASET graph..."
g = load_adjacency_list_from_csv(DS_CSV, ',', true)
n_v = nv(g)
m = ne(g)
@info "  Loaded: $n_v vertices, $m edges"

nls = extract_neighbor_lists(g)

# ============================================================================
# 1. BG Index Mode
# ============================================================================

@testset "BG Index Mode — $DATASET roundtrip" begin
    @info "--- BG Index Mode: write .mgz → load .mgz → verify ---"

    bg_base = joinpath(DS_DIR, "$(PREFIX)_bg_index")
    bg_mgz = bg_base * ".mgz"

    t_enc = time()
    write_bg_mgs3_graph(g, bg_base;
        coding_scheme=:index, integer_encoding=:fibonacci,
        ref_window_size=64, copy_blocks=true, stop_deltas=true)
    dt_enc = round(time() - t_enc, digits=2)

    @test isfile(bg_mgz)
    bpe = round(8.0 * filesize(bg_mgz) / m, digits=4)
    @info "  Encoded: $bpe BPE ($(filesize(bg_mgz)) bytes, $(dt_enc)s)"

    t_dec = time()
    g_loaded = load_compressed_mgs3_graph(bg_mgz)
    dt_dec = round(time() - t_dec, digits=2)

    @test nv(g_loaded) == n_v
    @test ne(g_loaded) == m
    loaded_nls = extract_neighbor_lists(g_loaded)
    @test verify_roundtrip(nls, loaded_nls, m)
    @info "  BG Index Mode roundtrip: PASSED ($bpe BPE, decode $(dt_dec)s)"
end

# ============================================================================
# 2. CS Index Mode
# ============================================================================

@testset "CS Index Mode — $DATASET roundtrip" begin
    @info "--- CS Index Mode: write .mgz → load .mgz → verify ---"

    cs_base = joinpath(DS_DIR, "$(PREFIX)_cs_index")
    cs_mgz = cs_base * ".mgz"

    t_enc = time()
    write_cs_mgs3_graph(g, cs_base;
        coding_scheme=:index, integer_encoding=:fibonacci,
        ref_window_size=64, compact_copy=true, tight_intervals=true)
    dt_enc = round(time() - t_enc, digits=2)

    @test isfile(cs_mgz)
    bpe = round(8.0 * filesize(cs_mgz) / m, digits=4)
    @info "  Encoded: $bpe BPE ($(filesize(cs_mgz)) bytes, $(dt_enc)s)"

    t_dec = time()
    g_loaded = load_compressed_mgs3_graph(cs_mgz)
    dt_dec = round(time() - t_dec, digits=2)

    @test nv(g_loaded) == n_v
    @test ne(g_loaded) == m
    loaded_nls = extract_neighbor_lists(g_loaded)
    @test verify_roundtrip(nls, loaded_nls, m)
    @info "  CS Index Mode roundtrip: PASSED ($bpe BPE, decode $(dt_dec)s)"
end

# ============================================================================
# 3. CG Index Mode (K=2)
# ============================================================================

@testset "CG Index Mode K=2 — $DATASET roundtrip" begin
    @info "--- CG Index Mode K=2: write .mgz → load .mgz → verify ---"

    TV = eltype(g)
    n_orig = Int(nv(g))

    # Leiden K=2 split
    part = leiden_partition(g; max_passes=8, max_levels=5)
    counts = Dict{Int,Int}()
    for c in part; counts[c] = get(counts, c, 0) + 1; end
    sorted_labels = sort(collect(keys(counts)); by = l -> counts[l], rev=true)
    top_label = sorted_labels[1]

    clusters = [TV[], TV[]]
    for i in 1:n_orig
        if part[i] == top_label
            push!(clusters[1], TV(i))
        else
            push!(clusters[2], TV(i))
        end
    end
    for C in clusters; sort!(C); end
    @info "  Clusters: $(length(clusters[1])) + $(length(clusters[2])) vertices"

    cg_params = CGParams(
        L=128,
        varint=:fibonacci, count_varint=:fibonacci,
        gap=:fibonacci, degree=:elias_delta,
        undirected_pairs=false,
        perm_strategy=:blockpos, membership=:implicit_ranges,
        inter_strategy=:perm,
        intra_ref_enabled=true, intra_ref_window=64,
        intra_ref_rle=false, intra_block_try=false,
        positions_mode=:delta, additions_mode=:delta,
        min_cluster_density=0.0,
        intra_intervals=false, intra_mil=4, intra_greedy_mil=false,
        intra_zigzag=true, intra_stop_deltas=true,
        intra_copy_blocks=true, intra_copy_adaptive=true,
        intra_ref_fixwidth=true, intra_ref_vlc=false,
        intra_add_adaptive=true, intra_raw_adaptive=true,
        intra_adapt_mil=4,   # must equal intra_mil; header encodes mil ∈ {3,4,5}
    )

    # Children mode reference (same clusters)
    ref_base = joinpath(DS_DIR, "$(PREFIX)_cg_k2_children")
    ref_mgz = ref_base * ".mgz"
    write_cg_mgs3_graph(g, ref_base, clusters;
        coding_scheme=:children, integer_encoding=:fibonacci,
        params=cg_params)
    g_children = load_compressed_mgs3_graph(ref_mgz)
    children_nls = extract_neighbor_lists(g_children)

    # Index mode
    cg_base = joinpath(DS_DIR, "$(PREFIX)_cg_k2_index")
    cg_mgz = cg_base * ".mgz"

    t_enc = time()
    write_cg_mgs3_graph(g, cg_base, clusters;
        coding_scheme=:index, integer_encoding=:fibonacci,
        params=cg_params)
    dt_enc = round(time() - t_enc, digits=2)

    @test isfile(cg_mgz)
    bpe = round(8.0 * filesize(cg_mgz) / m, digits=4)
    @info "  Encoded: $bpe BPE ($(filesize(cg_mgz)) bytes, $(dt_enc)s)"

    t_dec = time()
    g_loaded = load_compressed_mgs3_graph(cg_mgz)
    dt_dec = round(time() - t_dec, digits=2)

    @test nv(g_loaded) == n_v
    @test ne(g_loaded) == ne(g_children)
    loaded_nls = extract_neighbor_lists(g_loaded)
    m_children = ne(g_children)
    @test verify_roundtrip(children_nls, loaded_nls, m_children)
    @info "  CG Index Mode K=2 roundtrip: PASSED ($bpe BPE, decode $(dt_dec)s)"
end

end # isfile check
