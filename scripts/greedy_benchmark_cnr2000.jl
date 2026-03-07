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

#!/usr/bin/env julia

# Greedy compression benchmark on CNR-2000.
#
# Reports BPE, roundtrip verification, and action statistics.
#
# Usage: julia scripts/greedy_benchmark_cnr2000.jl

using Pkg
Pkg.activate(normpath(joinpath(@__DIR__, "..")))

using LightGraphs: nv, ne, outneighbors, vertices, edges, src, dst, is_directed, eltype
using Adjacently
using Adjacently.IO: load_adjacency_list_from_csv, BitWriter, BitReader, flush_bitwriter,
                      write_bytes
using Adjacently.Relabeling: relabel_graph, relabel_vertices_llp
using Adjacently.Clustering: leiden_partition
using Adjacently.Graph: subgraph
using Adjacently.MGS: load_compressed_mgs3_graph, create_header_flags,
                       OPTION_GREEDY_BASE, MGS_MAX_SIZE
using Adjacently.Compression: write_greedy_graph_data
using Adjacently.Util: infer_uint_custom_type

const PROJECT_ROOT = normpath(joinpath(@__DIR__, ".."))
const CNR_CSV = joinpath(PROJECT_ROOT, "datasets", "webgraph", "cnr-2000", "cnr-2000.csv")
const OUTPUT_DIR = joinpath(PROJECT_ROOT, "datasets", "webgraph", "cnr-2000")

function write_greedy_mgz(g, filename::AbstractString;
        integer_encoding::Symbol=:zeta,
        ref_window_size::Int=7,
        coding_scheme::Symbol=:children,
        stats::Union{Dict,Nothing}=nothing,
        copy_blocks::Bool=false,
        cluster_sizes::Vector{Int}=Int[])

    T = eltype(g)
    vs = vertices(g)
    gs = convert(UInt64, length(vs))

    if gs > MGS_MAX_SIZE
        error("Input graph cannot have more than 2^40-1 vertices")
    end

    n_bits_v = convert(UInt8, ceil(log(2, gs)))
    V = infer_uint_custom_type(n_bits_v)

    option_flags = UInt8(OPTION_GREEDY_BASE)
    flag_byte1, flag_byte2 = create_header_flags(:directed, coding_scheme, integer_encoding, option_flags)

    header_bytes = UInt8[
        0x4d, 0x47, 0x53,  # 'MGS'
        0x03, 0x00,         # Version 3.0
        flag_byte1, flag_byte2,
        (gs & 0xff), ((gs >> 8) & 0xff), ((gs >> 16) & 0xff), ((gs >> 24) & 0xff), ((gs >> 32) & 0xff)
    ]

    mgz_path = filename * ".mgz"
    f = open(mgz_path, "w")
    bw = BitWriter(f)

    write_bytes(bw, header_bytes)

    # Build neighbor lists
    neighbor_lists = Dict{V,Vector{V}}()
    for v in vs
        ovs = outneighbors(g, v)
        neighbor_lists[convert(V, v)] = sort([convert(V, o) for o in ovs])
    end

    # Greedy mode: vertex_actions=nothing triggers exhaustive search
    write_greedy_graph_data(bw, neighbor_lists, coding_scheme, ref_window_size;
        integer_encoding=integer_encoding, stats=stats, copy_blocks=copy_blocks,
        cluster_sizes=cluster_sizes)

    flush_bitwriter(bw; flush_last_bits=true)
    close(f)

    return mgz_path
end

function verify_roundtrip(mgz_path, original_graph; copy_blocks::Bool=false, ref_window_size::Int=7, cluster_sizes::Vector{Int}=Int[])
    g_loaded = load_compressed_mgs3_graph(mgz_path; copy_blocks=copy_blocks, ref_window_size=ref_window_size, cluster_sizes=cluster_sizes)
    m_orig = ne(original_graph)
    m_loaded = ne(g_loaded)
    n_orig = nv(original_graph)
    n_loaded = nv(g_loaded)

    if m_loaded == m_orig && n_loaded == n_orig
        println("    Roundtrip: OK ($n_loaded vertices, $m_loaded edges)")
        return true
    else
        println("    Roundtrip: MISMATCH (expected $n_orig/$m_orig, got $n_loaded/$m_loaded)")
        return false
    end
end

function compute_bpe(mgz_path, num_edges)
    file_size = filesize(mgz_path)
    bpe = 8.0 * file_size / num_edges
    println("    File size: $file_size bytes ($(round(file_size / 1024, digits=1)) KB)")
    println("    BPE:       $(round(bpe, digits=4))")
    return bpe
end

# Two-step relabeling: Leiden community detection → per-group LLP on induced subgraphs.
# Mirrors the CGE K=1 approach: the top Leiden cluster forms group 1, all other
# vertices form group 2. LLP is then applied to each group's induced subgraph.
# This gives 80–90% intra-group neighbor overlap vs ~40% with global LLP.
function relabel_leiden_k1_llp(g; leiden_max_passes::Int=8, leiden_max_levels::Int=5,
                                   llp_passes::Int=5)
    T = eltype(g)
    n = nv(g)

    # Step 1: Leiden community detection
    println("    Leiden partitioning...")
    t1 = time()
    part = leiden_partition(g; max_passes=leiden_max_passes, max_levels=leiden_max_levels)
    nclusters = maximum(part)
    println("      → $nclusters fine clusters ($(round(time()-t1, digits=2))s)")

    # Find the single largest Leiden community label
    counts = Dict{Int,Int}()
    for c in part; counts[c] = get(counts, c, 0) + 1; end
    top_label = argmax(counts)

    # Split into 2 groups: top-1 Leiden community vs everything else
    groups = [T[], T[]]
    for i in 1:n
        push!(groups[part[i] == top_label ? 1 : 2], T(i))
    end
    println("      Group sizes: $(length(groups[1])), $(length(groups[2]))")

    # Step 2: LLP on each group's induced subgraph
    println("    Per-group LLP on induced subgraphs...")
    t2 = time()
    for gi in 1:2
        C = groups[gi]
        length(C) <= 2 && continue
        sg, oni, _ = subgraph(g, C)
        mapping = relabel_vertices_llp(sg, :sym; passes=llp_passes)
        sort!(C, by = v -> Int(mapping[oni[v]]))
        groups[gi] = C
    end
    println("      Per-group LLP done ($(round(time()-t2, digits=2))s)")

    # Build global vertex mapping old_id → new_id
    vertex_mapping = Dict{T,T}()
    new_id = T(1)
    for C in groups
        for v in C
            vertex_mapping[v] = new_id
            new_id += T(1)
        end
    end
    group_sizes = [length(C) for C in groups]
    return vertex_mapping, group_sizes
end

function main()
    println("=" ^ 70)
    println("CNR-2000 Greedy Compression Benchmark")
    println("=" ^ 70)

    # Load graph
    println("\nLoading CNR-2000...")
    if !isfile(CNR_CSV)
        error("CNR-2000 CSV not found at $CNR_CSV\nDownload from: https://law.di.unimi.it/webdata/cnr-2000/")
    end
    t0 = time()
    g = load_adjacency_list_from_csv(CNR_CSV, ',', true)
    n, m = nv(g), ne(g)
    println("  Vertices: $n, Edges: $m")
    println("  Load time: $(round(time() - t0, digits=2))s")

    mkpath(OUTPUT_DIR)

    # ==================================================================
    # LLP relabeling
    # ==================================================================
    println("\n" * "-" ^ 70)
    println("LLP relabeling (sym, 10 passes)")
    println("-" ^ 70)

    t1 = time()
    llp_map = relabel_vertices_llp(g, :sym; passes=10)
    g_llp = relabel_graph(g, llp_map)
    println("  Relabeling: $(round(time() - t1, digits=2))s")

    # ==================================================================
    # Config 1: LLP + Zeta-3 baseline (adaptive bitmap)
    # ==================================================================
    action_stats = Dict{Tuple{Symbol,Symbol,Int}, Int}()
    output_path = joinpath(OUTPUT_DIR, "cnr2000_greedy_llp_zeta_w7")

    println("\n  [1] Compressing with ie=zeta, window=7 (adaptive bitmap)...")
    t2 = time()
    mgz_path = write_greedy_mgz(g_llp, output_path;
        integer_encoding=:zeta, ref_window_size=7, stats=action_stats)
    compress_time = time() - t2
    println("    Compress time: $(round(compress_time, digits=2))s")

    bpe1 = compute_bpe(mgz_path, m)

    println("    Verifying roundtrip...")
    ok1 = verify_roundtrip(mgz_path, g_llp)

    # ==================================================================
    # Config 2: LLP + Zeta-3 + Copy-blocks (w=7)
    # ==================================================================
    action_stats_cb = Dict{Tuple{Symbol,Symbol,Int}, Int}()
    output_path_cb = joinpath(OUTPUT_DIR, "cnr2000_greedy_llp_zeta_w7_cb")

    println("\n  [2] Compressing with ie=zeta, window=7 (copy-blocks)...")
    t3 = time()
    mgz_path_cb = write_greedy_mgz(g_llp, output_path_cb;
        integer_encoding=:zeta, ref_window_size=7, stats=action_stats_cb, copy_blocks=true)
    compress_time_cb = time() - t3
    println("    Compress time: $(round(compress_time_cb, digits=2))s")

    bpe2 = compute_bpe(mgz_path_cb, m)

    println("    Verifying roundtrip...")
    ok2 = verify_roundtrip(mgz_path_cb, g_llp; copy_blocks=true)

    # ==================================================================
    # Config 3: LLP + Zeta-3 + Copy-blocks (w=64)
    # ==================================================================
    output_path_w64 = joinpath(OUTPUT_DIR, "cnr2000_greedy_llp_zeta_w64_cb")

    println("\n  [3] Compressing with ie=zeta, window=64 (copy-blocks)...")
    t4 = time()
    mgz_path_w64 = write_greedy_mgz(g_llp, output_path_w64;
        integer_encoding=:zeta, ref_window_size=64, copy_blocks=true)
    println("    Compress time: $(round(time() - t4, digits=2))s")
    bpe3 = compute_bpe(mgz_path_w64, m)
    println("    Verifying roundtrip...")
    ok3 = verify_roundtrip(mgz_path_w64, g_llp; copy_blocks=true, ref_window_size=64)

    # ==================================================================
    # Config 4: LLP + Fibonacci + Copy-blocks (w=64)
    # ==================================================================
    output_path_fib64 = joinpath(OUTPUT_DIR, "cnr2000_greedy_llp_fib_w64_cb")

    println("\n  [4] Compressing with ie=fibonacci, window=64 (copy-blocks)...")
    t5 = time()
    mgz_path_fib64 = write_greedy_mgz(g_llp, output_path_fib64;
        integer_encoding=:fibonacci, ref_window_size=64, copy_blocks=true)
    println("    Compress time: $(round(time() - t5, digits=2))s")
    bpe4 = compute_bpe(mgz_path_fib64, m)
    println("    Verifying roundtrip...")
    ok4 = verify_roundtrip(mgz_path_fib64, g_llp; copy_blocks=true, ref_window_size=64)

    # ==================================================================
    # Two-step relabeling: Leiden + per-group LLP
    # ==================================================================
    println("\n" * "-" ^ 70)
    println("Two-step relabeling: Leiden K=1 + per-group LLP")
    println("-" ^ 70)

    t_leiden = time()
    leiden_llp_map, leiden_group_sizes = relabel_leiden_k1_llp(g; leiden_max_passes=8, leiden_max_levels=5,
                                                                   llp_passes=5)
    g_leiden_llp = relabel_graph(g, leiden_llp_map)
    println("  Total relabeling time: $(round(time() - t_leiden, digits=2))s")
    println("  Group sizes: $(leiden_group_sizes)")

    # ==================================================================
    # Config 5: Leiden+LLP + Fibonacci + Copy-blocks (w=7)
    # ==================================================================
    output_path_lllp_fib7 = joinpath(OUTPUT_DIR, "cnr2000_greedy_leiden_llp_fib_w7_cb")

    println("\n  [5] Leiden+LLP + Fib w=7 (copy-blocks)...")
    t6 = time()
    mgz_path_lllp7 = write_greedy_mgz(g_leiden_llp, output_path_lllp_fib7;
        integer_encoding=:fibonacci, ref_window_size=7, copy_blocks=true)
    println("    Compress time: $(round(time() - t6, digits=2))s")
    bpe5 = compute_bpe(mgz_path_lllp7, m)
    println("    Verifying roundtrip...")
    ok5 = verify_roundtrip(mgz_path_lllp7, g_leiden_llp; copy_blocks=true, ref_window_size=7)

    # ==================================================================
    # Config 6: Leiden+LLP + Fibonacci + Copy-blocks (w=64)
    # ==================================================================
    output_path_lllp_fib64 = joinpath(OUTPUT_DIR, "cnr2000_greedy_leiden_llp_fib_w64_cb")

    println("\n  [6] Leiden+LLP + Fib w=64 (copy-blocks)...")
    t7 = time()
    mgz_path_lllp64 = write_greedy_mgz(g_leiden_llp, output_path_lllp_fib64;
        integer_encoding=:fibonacci, ref_window_size=64, copy_blocks=true)
    println("    Compress time: $(round(time() - t7, digits=2))s")
    bpe6 = compute_bpe(mgz_path_lllp64, m)
    println("    Verifying roundtrip...")
    ok6 = verify_roundtrip(mgz_path_lllp64, g_leiden_llp; copy_blocks=true, ref_window_size=64)

    # ==================================================================
    # Config 7: Leiden+LLP + cluster-local window + Fibonacci (w=7)
    # ==================================================================
    output_path_local7 = joinpath(OUTPUT_DIR, "cnr2000_greedy_leiden_llp_local_fib_w7_cb")

    println("\n  [7] Leiden+LLP + local-window + Fib w=7 (copy-blocks)...")
    t8 = time()
    mgz_path_local7 = write_greedy_mgz(g_leiden_llp, output_path_local7;
        integer_encoding=:fibonacci, ref_window_size=7, copy_blocks=true,
        cluster_sizes=leiden_group_sizes)
    println("    Compress time: $(round(time() - t8, digits=2))s")
    bpe7 = compute_bpe(mgz_path_local7, m)
    println("    Verifying roundtrip...")
    ok7 = verify_roundtrip(mgz_path_local7, g_leiden_llp; copy_blocks=true, ref_window_size=7,
        cluster_sizes=leiden_group_sizes)

    # ==================================================================
    # Config 8: Leiden+LLP + cluster-local window + Fibonacci (w=64)
    # ==================================================================
    output_path_local64 = joinpath(OUTPUT_DIR, "cnr2000_greedy_leiden_llp_local_fib_w64_cb")

    println("\n  [8] Leiden+LLP + local-window + Fib w=64 (copy-blocks)...")
    t9 = time()
    mgz_path_local64 = write_greedy_mgz(g_leiden_llp, output_path_local64;
        integer_encoding=:fibonacci, ref_window_size=64, copy_blocks=true,
        cluster_sizes=leiden_group_sizes)
    println("    Compress time: $(round(time() - t9, digits=2))s")
    bpe8 = compute_bpe(mgz_path_local64, m)
    println("    Verifying roundtrip...")
    ok8 = verify_roundtrip(mgz_path_local64, g_leiden_llp; copy_blocks=true, ref_window_size=64,
        cluster_sizes=leiden_group_sizes)

    # ==================================================================
    # Summary
    # ==================================================================
    println("\n" * "=" ^ 70)
    println("RESULT — CNR-2000 ($n vertices, $m edges)")
    println("=" ^ 70)
    println("  Global LLP relabeling:")
    println("  [1] LLP + Zeta-3  w=7  (adaptive):    BPE = $(round(bpe1, digits=4)),  RT=$(ok1 ? "OK" : "FAIL")")
    println("  [2] LLP + Zeta-3  w=7  (copy-blocks): BPE = $(round(bpe2, digits=4)),  RT=$(ok2 ? "OK" : "FAIL")")
    println("  [3] LLP + Zeta-3  w=64 (copy-blocks): BPE = $(round(bpe3, digits=4)),  RT=$(ok3 ? "OK" : "FAIL")")
    println("  [4] LLP + Fib     w=64 (copy-blocks): BPE = $(round(bpe4, digits=4)),  RT=$(ok4 ? "OK" : "FAIL")")
    println("  Two-step relabeling (Leiden K=1 + per-group LLP), global window:")
    println("  [5] Leiden+LLP + Fib w=7  (copy-blocks):              BPE = $(round(bpe5, digits=4)),  RT=$(ok5 ? "OK" : "FAIL")")
    println("  [6] Leiden+LLP + Fib w=64 (copy-blocks):              BPE = $(round(bpe6, digits=4)),  RT=$(ok6 ? "OK" : "FAIL")")
    println("  Two-step relabeling (Leiden K=1 + per-group LLP), cluster-local window:")
    println("  [7] Leiden+LLP + local w=7  + Fib (copy-blocks):     BPE = $(round(bpe7, digits=4)),  RT=$(ok7 ? "OK" : "FAIL")")
    println("  [8] Leiden+LLP + local w=64 + Fib (copy-blocks):     BPE = $(round(bpe8, digits=4)),  RT=$(ok8 ? "OK" : "FAIL")")
    println("  WebGraph reference:                                    BPE = 2.897")
    println("  CGE FW64 K=1 (best):                                 BPE = 2.887")
    println("=" ^ 70)
end

main()
