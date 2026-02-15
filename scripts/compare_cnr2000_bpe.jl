#!/usr/bin/env julia

# Compare 3 compression approaches on CNR-2000:
#   1. Legacy ASTRA (RCM ordering + complex/fibonacci + ref window 1024)
#   2. RL greedy + LLP ordering
#   3. RL greedy + GCN v2 ordering (enhanced 10-feature)
#
# All .mgz files are saved under datasets/webgraph/cnr-2000/
# and roundtrip-verified.
#
# Usage: julia scripts/compare_cnr2000_bpe.jl

using Pkg
Pkg.activate(normpath(joinpath(@__DIR__, "..")))

using LightGraphs: nv, ne, outneighbors, vertices, eltype
using Adjacently
using Adjacently.IO: load_adjacency_list_from_csv
using Adjacently.Relabeling: relabel_vertices_rcm, relabel_vertices_llp, relabel_graph
using Adjacently.GCN: GCNModel, gcn_ordering, train_gcn_proxy!, TrainConfig, save_gcn_model
using Adjacently.MGS: write_compressed_mgs3_graph, write_rl_compressed_mgs3_graph,
                       load_compressed_mgs3_graph

const PROJECT_ROOT = normpath(joinpath(@__DIR__, ".."))
const CNR_CSV = joinpath(PROJECT_ROOT, "datasets", "webgraph", "cnr-2000", "cnr-2000.csv")
const OUTPUT_DIR = joinpath(PROJECT_ROOT, "datasets", "webgraph", "cnr-2000")

function verify_roundtrip(mgz_path, expected_edges)
    g_loaded = load_compressed_mgs3_graph(mgz_path)
    m2 = ne(g_loaded)
    if m2 == expected_edges
        println("    Roundtrip: OK ($m2 edges)")
        return true
    else
        println("    Roundtrip: MISMATCH (expected $expected_edges, got $m2)")
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
    println("CNR-2000 Compression Comparison (3 approaches)")
    println("=" ^ 70)

    # Load graph
    println("\nLoading CNR-2000...")
    if !isfile(CNR_CSV)
        error("CNR-2000 CSV not found at $CNR_CSV")
    end
    t0 = time()
    g = load_adjacency_list_from_csv(CNR_CSV, ',', true)
    n, m = nv(g), ne(g)
    println("  Vertices: $n, Edges: $m")
    println("  Load time: $(round(time() - t0, digits=2))s")

    mkpath(OUTPUT_DIR)
    results = Dict{String, Float64}()

    # ==================================================================
    # 1. Legacy ASTRA: RCM ordering + complex/fibonacci + ref_window=1024
    # ==================================================================
    println("\n" * "-" ^ 70)
    println("[1/3] Legacy ASTRA (RCM + fibonacci + ref_window=1024)")
    println("-" ^ 70)

    t1 = time()
    rcm_mapping = relabel_vertices_rcm(g, :out)
    g_rcm = relabel_graph(g, rcm_mapping)
    println("  RCM relabeling: $(round(time() - t1, digits=2))s")

    output_legacy = joinpath(OUTPUT_DIR, "cnr2000_legacy_astra")
    t2 = time()
    write_compressed_mgs3_graph(
        g_rcm,
        output_legacy,
        :children,      # coding_scheme
        :complex,       # compression
        :fibonacci,     # integer_encoding
        true,           # use_mix_mode
        true,           # reference_enabled
        true,           # recursive_reference
        1024            # ref_window_size
    )
    println("  Compress time: $(round(time() - t2, digits=2))s")

    mgz_legacy = output_legacy * ".mgz"
    bpe_legacy = compute_bpe(mgz_legacy, m)
    verify_roundtrip(mgz_legacy, ne(g_rcm))
    results["Legacy ASTRA (RCM)"] = bpe_legacy

    # ==================================================================
    # 2. RL greedy + LLP ordering
    # ==================================================================
    println("\n" * "-" ^ 70)
    println("[2/3] RL Greedy + LLP ordering (sym, 5 passes)")
    println("-" ^ 70)

    t3 = time()
    llp_mapping = relabel_vertices_llp(g, :sym; passes=5)
    g_llp = relabel_graph(g, llp_mapping)
    println("  LLP relabeling: $(round(time() - t3, digits=2))s")

    output_rl_llp = joinpath(OUTPUT_DIR, "cnr2000_rl_greedy_llp")
    t4 = time()
    write_rl_compressed_mgs3_graph(
        g_llp,
        output_rl_llp,
        nothing,        # no policy = greedy mode
        1;
        coding_scheme=:children,
        ref_window_size=7,
        integer_encoding=:fibonacci
    )
    println("  Compress time: $(round(time() - t4, digits=2))s")

    mgz_rl_llp = output_rl_llp * ".mgz"
    bpe_rl_llp = compute_bpe(mgz_rl_llp, m)
    verify_roundtrip(mgz_rl_llp, ne(g_llp))
    results["RL Greedy (LLP)"] = bpe_rl_llp

    # ==================================================================
    # 3. RL greedy + GCN v2 ordering (enhanced 10-feature)
    # ==================================================================
    println("\n" * "-" ^ 70)
    println("[3/3] RL Greedy + GCN v2 ordering (10 features, hidden=64, 200 epochs)")
    println("-" ^ 70)

    t5 = time()
    println("  Building enhanced GCN model...")
    model = GCNModel(g; hidden_dim=64, seed=42)
    println("  Model: $(model.input_dim) -> $(model.hidden_dim) -> 1")
    println("  Feature extraction: $(round(time() - t5, digits=2))s")

    t6 = time()
    config = TrainConfig(proxy_epochs=200, proxy_lr=0.001,
                         reinforce_epochs=0, reinforce_lr=0.0001,
                         sigma=0.1, baseline_ema=0.9)
    proxy_losses = train_gcn_proxy!(model, g, config)
    println("  Training time: $(round(time() - t6, digits=2))s")
    println("  Proxy loss: $(round(proxy_losses[1], digits=2)) -> $(round(proxy_losses[end], digits=2))")

    V = eltype(g)
    gcn_map = gcn_ordering(model, V)
    g_gcn = relabel_graph(g, gcn_map)

    output_rl_gcn = joinpath(OUTPUT_DIR, "cnr2000_rl_greedy_gcn_v2")
    t7 = time()
    write_rl_compressed_mgs3_graph(
        g_gcn,
        output_rl_gcn,
        nothing,        # no policy = greedy mode
        1;
        coding_scheme=:children,
        ref_window_size=7,
        integer_encoding=:fibonacci
    )
    println("  Compress time: $(round(time() - t7, digits=2))s")

    mgz_rl_gcn = output_rl_gcn * ".mgz"
    bpe_rl_gcn = compute_bpe(mgz_rl_gcn, m)
    verify_roundtrip(mgz_rl_gcn, ne(g_gcn))
    results["RL Greedy (GCN v2)"] = bpe_rl_gcn

    # ==================================================================
    # Summary
    # ==================================================================
    println("\n" * "=" ^ 70)
    println("COMPARISON SUMMARY — CNR-2000 ($n vertices, $m edges)")
    println("=" ^ 70)
    println()
    println("  Approach                    BPE")
    println("  " * "-" ^ 42)

    # Sort by BPE ascending (best first)
    sorted = sort(collect(results), by=x -> x[2])
    for (name, bpe) in sorted
        marker = name == first(sorted)[1] ? " <-- best" : ""
        println("  $(rpad(name, 28)) $(round(bpe, digits=4))$marker")
    end
    println()
    println("  Output files saved in: $OUTPUT_DIR")
    println("=" ^ 70)
end

main()
