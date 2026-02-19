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
        stats::Union{Dict,Nothing}=nothing)

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
        integer_encoding=integer_encoding, vertex_actions=nothing, stats=stats)

    flush_bitwriter(bw; flush_last_bits=true)
    close(f)

    return mgz_path
end

function verify_roundtrip(mgz_path, original_graph)
    g_loaded = load_compressed_mgs3_graph(mgz_path)
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
    # LLP + Zeta-3 with action statistics
    # ==================================================================
    action_stats = Dict{Tuple{Symbol,Symbol,Int}, Int}()
    output_path = joinpath(OUTPUT_DIR, "cnr2000_greedy_llp_zeta_w7")

    println("\n  Compressing with ie=zeta, window=7 (VLC headers)...")
    t2 = time()
    mgz_path = write_greedy_mgz(g_llp, output_path;
        integer_encoding=:zeta, ref_window_size=7, stats=action_stats)
    compress_time = time() - t2
    println("    Compress time: $(round(compress_time, digits=2))s")

    bpe = compute_bpe(mgz_path, m)

    println("    Verifying roundtrip...")
    ok = verify_roundtrip(mgz_path, g_llp)

    # ==================================================================
    # Summary
    # ==================================================================
    println("\n" * "=" ^ 70)
    println("RESULT — CNR-2000 ($n vertices, $m edges)")
    println("=" ^ 70)
    rt = ok ? "OK" : "FAIL"
    println("  LLP + Zeta-3 (VLC):  BPE = $(round(bpe, digits=4)),  Roundtrip = $rt")
    println("  WebGraph reference:  BPE = 2.90")
    println("=" ^ 70)

    # ==================================================================
    # Action Statistics
    # ==================================================================
    println("\n" * "=" ^ 70)
    println("ACTION STATISTICS — LLP + Zeta-3")
    println("=" ^ 70)
    total_vertices = sum(values(action_stats))
    sorted_stats = sort(collect(action_stats), by=x -> -x[2])

    entropy = 0.0
    for (_, count) in sorted_stats
        p = count / total_vertices
        if p > 0; entropy -= p * log2(p); end
    end
    println("\n  Total vertices: $total_vertices")
    println("  Distinct action combos: $(length(sorted_stats))")
    println("  Shannon entropy:        $(round(entropy, digits=4)) bits")
    println()

    println("  $(rpad("Ref Mode", 12)) $(rpad("Enc Type", 12)) $(rpad("MIL", 6)) $(rpad("Count", 10)) $(rpad("Pct %", 10)) Cum %")
    println("  " * "-" ^ 62)
    cum_pct = 0.0
    for (key, count) in sorted_stats
        ref_mode, enc_type, mil = key
        pct = 100.0 * count / total_vertices
        cum_pct += pct
        mil_str = enc_type == :delta ? "-" : string(mil)
        println("  $(rpad(ref_mode, 12)) $(rpad(enc_type, 12)) $(rpad(mil_str, 6)) $(rpad(count, 10)) $(rpad(round(pct, digits=2), 10)) $(round(cum_pct, digits=1))")
    end
    println()
end

main()
