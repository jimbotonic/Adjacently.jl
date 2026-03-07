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

#
# Roundtrip tests for legacy uncompressed MGS and Huffman compression
# on CNR-2000 (both children and index modes).
# Compares BPE with best known results.
#
# Usage:
#   julia --project test/run_tests_legacy_huffman.jl
#

using Pkg
Pkg.activate(normpath(joinpath(@__DIR__, "..")))

using Test
using Logging
global_logger(ConsoleLogger(stderr, Logging.Info))
using Printf

using LightGraphs: nv, ne, outneighbors, vertices, eltype
using Adjacently
using Adjacently.IO: load_adjacency_list_from_csv
using Adjacently.MGS: write_mgs3_graph, load_mgs3_graph,
    write_huffman_compressed_mgs3_graph, load_compressed_mgs3_graph

const DS_DIR = normpath(joinpath(@__DIR__, "..", "datasets", "webgraph", "cnr-2000"))
const DS_CSV = joinpath(DS_DIR, "cnr-2000.csv")

function extract_neighbor_lists(g)
    T = eltype(g)
    nls = Dict{T,Vector{T}}()
    for v in vertices(g)
        nls[T(v)] = sort([T(o) for o in outneighbors(g, v)])
    end
    return nls
end

function verify_roundtrip(orig_nls, loaded_nls, expected_edges)
    total = sum(length(v) for v in values(loaded_nls))
    total != expected_edges && return false
    for (v, orig) in orig_nls
        haskey(loaded_nls, v) || return false
        sort(loaded_nls[v]) == sort(orig) || return false
    end
    return true
end

if !isfile(DS_CSV)
    @warn "CNR-2000 CSV not found at $DS_CSV; skipping"
    @testset "Legacy & Huffman" begin
        @test_skip "CNR-2000 dataset unavailable"
    end
else

@info "Loading CNR-2000 graph..."
t0 = time()
g = load_adjacency_list_from_csv(DS_CSV, ',', true)
n_v = nv(g)
m = ne(g)
@info "  Loaded in $(round(time()-t0, digits=1))s: $n_v vertices, $m edges"

orig_nls = extract_neighbor_lists(g)

# Known best BPE for reference
const BEST_CGE_BPE = 2.3191

results = Dict{String, Float64}()

# ==========================================================================
# 1. Legacy Uncompressed MGS — children mode
# ==========================================================================
@testset "Legacy Uncompressed — children mode" begin
    @info "--- Legacy children: write .mgs → load .mgs → verify ---"
    base = joinpath(DS_DIR, "cnr2000_legacy_children")
    mgs_file = base * ".mgs"

    t = time()
    write_mgs3_graph(g, base, :children)
    dt = round(time() - t, digits=2)
    @test isfile(mgs_file)
    sz = filesize(mgs_file)
    bpe = round(8.0 * sz / m, digits=4)
    @info "  Encoded: $bpe BPE ($sz bytes, $(dt)s)"
    results["Legacy children"] = bpe

    t = time()
    g_loaded = load_mgs3_graph(mgs_file)
    dt = round(time() - t, digits=2)
    @test nv(g_loaded) == n_v
    @test ne(g_loaded) == m
    loaded_nls = extract_neighbor_lists(g_loaded)
    @test verify_roundtrip(orig_nls, loaded_nls, m)
    @info "  Roundtrip: PASSED ($bpe BPE, decode $(dt)s)"

    rm(mgs_file)
end

# ==========================================================================
# 2. Legacy Uncompressed MGS — index mode
# ==========================================================================
@testset "Legacy Uncompressed — index mode" begin
    @info "--- Legacy index: write .mgs → load .mgs → verify ---"
    base = joinpath(DS_DIR, "cnr2000_legacy_index")
    mgs_file = base * ".mgs"

    t = time()
    write_mgs3_graph(g, base, :index)
    dt = round(time() - t, digits=2)
    @test isfile(mgs_file)
    sz = filesize(mgs_file)
    bpe = round(8.0 * sz / m, digits=4)
    @info "  Encoded: $bpe BPE ($sz bytes, $(dt)s)"
    results["Legacy index"] = bpe

    t = time()
    g_loaded = load_mgs3_graph(mgs_file)
    dt = round(time() - t, digits=2)
    @test nv(g_loaded) == n_v
    @test ne(g_loaded) == m
    loaded_nls = extract_neighbor_lists(g_loaded)
    @test verify_roundtrip(orig_nls, loaded_nls, m)
    @info "  Roundtrip: PASSED ($bpe BPE, decode $(dt)s)"

    rm(mgs_file)
end

# ==========================================================================
# 3. Huffman — children mode
# ==========================================================================
@testset "Huffman — children mode" begin
    @info "--- Huffman children: write .mgz → load .mgz → verify ---"
    base = joinpath(DS_DIR, "cnr2000_huffman_children")
    mgz_file = base * ".mgz"

    t = time()
    write_huffman_compressed_mgs3_graph(g, base, :children)
    dt = round(time() - t, digits=2)
    @test isfile(mgz_file)
    sz = filesize(mgz_file)
    bpe = round(8.0 * sz / m, digits=4)
    @info "  Encoded: $bpe BPE ($sz bytes, $(dt)s)"
    results["Huffman children"] = bpe

    t = time()
    g_loaded = load_compressed_mgs3_graph(mgz_file)
    dt = round(time() - t, digits=2)
    @test nv(g_loaded) == n_v
    @test ne(g_loaded) == m
    loaded_nls = extract_neighbor_lists(g_loaded)
    @test verify_roundtrip(orig_nls, loaded_nls, m)
    @info "  Roundtrip: PASSED ($bpe BPE, decode $(dt)s)"

    rm(mgz_file)
end

# ==========================================================================
# 4. Huffman — index mode
# ==========================================================================
@testset "Huffman — index mode" begin
    @info "--- Huffman index: write .mgz → load .mgz → verify ---"
    base = joinpath(DS_DIR, "cnr2000_huffman_index")
    mgz_file = base * ".mgz"

    t = time()
    write_huffman_compressed_mgs3_graph(g, base, :index)
    dt = round(time() - t, digits=2)
    @test isfile(mgz_file)
    sz = filesize(mgz_file)
    bpe = round(8.0 * sz / m, digits=4)
    @info "  Encoded: $bpe BPE ($sz bytes, $(dt)s)"
    results["Huffman index"] = bpe

    t = time()
    g_loaded = load_compressed_mgs3_graph(mgz_file)
    dt = round(time() - t, digits=2)
    @test nv(g_loaded) == n_v
    @test ne(g_loaded) == m
    loaded_nls = extract_neighbor_lists(g_loaded)
    @test verify_roundtrip(orig_nls, loaded_nls, m)
    @info "  Roundtrip: PASSED ($bpe BPE, decode $(dt)s)"

    rm(mgz_file)
end

# ==========================================================================
# Summary table
# ==========================================================================
println()
println("=" ^ 72)
println("  CNR-2000 BPE Summary: Legacy & Huffman vs Best Known")
println("  ($n_v vertices, $m edges)")
println("=" ^ 72)
println()
@printf("  %-22s %10s %12s\n", "Mode", "BPE", "vs CGE best")
println("  " * "-" ^ 46)

for label in ["Legacy children", "Legacy index", "Huffman children", "Huffman index"]
    bpe = results[label]
    ratio = round(bpe / BEST_CGE_BPE, digits=2)
    @printf("  %-22s %10.4f %10.1fx\n", label, bpe, ratio)
end

@printf("\n  %-22s %10.4f %10s\n", "CGE best (K=2)", BEST_CGE_BPE, "1.0x")
println()

end # isfile check
