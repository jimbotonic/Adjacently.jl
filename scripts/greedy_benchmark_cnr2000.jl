#!/usr/bin/env julia

# Greedy compression benchmark on CNR-2000 — LLP + Zeta-3.
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
using Adjacently.MGS: load_compressed_mgs3_graph, create_header_flags,
                       OPTION_RL_POLICY_BASE, MGS_MAX_SIZE
using Adjacently.Compression: write_rl_compressed_graph_data
using Adjacently.Util: infer_uint_custom_type

const PROJECT_ROOT = normpath(joinpath(@__DIR__, ".."))
const CNR_CSV = joinpath(PROJECT_ROOT, "datasets", "webgraph", "cnr-2000", "cnr-2000.csv")
const OUTPUT_DIR = joinpath(PROJECT_ROOT, "datasets", "webgraph", "cnr-2000")

function write_greedy_mgz(g, filename::AbstractString;
        integer_encoding::Symbol=:zeta,
        ref_window_size::Int=7,
        coding_scheme::Symbol=:children,
        stats::Union{Dict,Nothing}=nothing,
        copy_blocks::Bool=false)

    T = eltype(g)
    vs = vertices(g)
    gs = convert(UInt64, length(vs))

    if gs > MGS_MAX_SIZE
        error("Input graph cannot have more than 2^40-1 vertices")
    end

    n_bits_v = convert(UInt8, ceil(log(2, gs)))
    V = infer_uint_custom_type(n_bits_v)

    option_flags = UInt8(OPTION_RL_POLICY_BASE)
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
    write_rl_compressed_graph_data(bw, neighbor_lists, coding_scheme, ref_window_size;
        integer_encoding=integer_encoding, vertex_actions=nothing, stats=stats, copy_blocks=copy_blocks)

    flush_bitwriter(bw; flush_last_bits=true)
    close(f)

    return mgz_path
end

function verify_roundtrip(mgz_path, original_graph; copy_blocks::Bool=false)
    g_loaded = load_compressed_mgs3_graph(mgz_path; copy_blocks=copy_blocks)
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

function main()
    println("=" ^ 70)
    println("CNR-2000 Greedy Compression Benchmark — LLP + Zeta-3")
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
    # Config 2: LLP + Zeta-3 + Copy-blocks
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
    # Summary
    # ==================================================================
    println("\n" * "=" ^ 70)
    println("RESULT — CNR-2000 ($n vertices, $m edges)")
    println("=" ^ 70)
    rt1 = ok1 ? "OK" : "FAIL"
    rt2 = ok2 ? "OK" : "FAIL"
    println("  LLP + Zeta-3 (adaptive):    BPE = $(round(bpe1, digits=4)),  Roundtrip = $rt1")
    println("  LLP + Zeta-3 (copy-blocks): BPE = $(round(bpe2, digits=4)),  Roundtrip = $rt2")
    println("  WebGraph reference:         BPE = 2.90")
    println("=" ^ 70)
end

main()
