#!/usr/bin/env julia

# Retrain GCN with enhanced features (10-dim) and test CNR-2000 compression.
#
# Uses the existing RL policy (cnr2000_gcn_rl_policy.qpolicy) to isolate
# the effect of improved vertex ordering from richer node features.
#
# Enhanced features: degrees (3), clustering (1), PageRank (1),
#                    neighbor degree stats (2), spectral coordinates (3)
#
# Usage: julia scripts/retrain_gcn_test_cnr2000.jl

using Pkg
Pkg.activate(normpath(joinpath(@__DIR__, "..")))

using LightGraphs: nv, ne, outneighbors, vertices, eltype
using Adjacently
using Adjacently.IO: load_adjacency_list_from_csv
using Adjacently.Relabeling: relabel_graph
using Adjacently.GCN: GCNModel, gcn_ordering, train_gcn_proxy!,
                       TrainConfig, save_gcn_model
using Adjacently.MGS: write_rl_compressed_mgs3_graph, load_compressed_mgs3_graph

const PROJECT_ROOT = normpath(joinpath(@__DIR__, ".."))
const CNR_CSV = joinpath(PROJECT_ROOT, "datasets", "webgraph", "cnr-2000", "cnr-2000.csv")
const RL_POLICY = joinpath(PROJECT_ROOT, "policies", "cnr2000_gcn_rl_policy.qpolicy")
const MODEL_OUTPUT = joinpath(PROJECT_ROOT, "policies", "cnr2000_gcn_v2_ordering.bin")

function main()
    println("=" ^ 70)
    println("Enhanced GCN Ordering (v2) — CNR-2000")
    println("=" ^ 70)

    # =====================================================================
    # 1. Load graph
    # =====================================================================
    println("\n[1/6] Loading CNR-2000 dataset...")
    if !isfile(CNR_CSV)
        error("CNR-2000 CSV not found at $CNR_CSV")
    end
    t0 = time()
    g = load_adjacency_list_from_csv(CNR_CSV, ',', true)
    load_time = time() - t0
    n, m = nv(g), ne(g)
    println("  Vertices: $n")
    println("  Edges:    $m")
    println("  Load time: $(round(load_time, digits=2))s")

    # Verify RL policy exists
    if !isfile(RL_POLICY)
        error("RL policy not found at $RL_POLICY — train it first with train_rl_gcn_cnr2000.jl")
    end
    println("  RL policy: $RL_POLICY")

    # =====================================================================
    # 2. Create enhanced GCN model (10 features, hidden_dim=64)
    # =====================================================================
    println("\n[2/6] Creating enhanced GCN model (10 features, hidden_dim=64)...")
    t1 = time()
    model = GCNModel(g; hidden_dim=64, seed=42)
    build_time = time() - t1
    println("  Model: $(model.input_dim) → $(model.hidden_dim) → 1")
    println("  Parameters: $(model.input_dim * model.hidden_dim + model.hidden_dim) weights")
    println("  Feature extraction time: $(round(build_time, digits=2))s")

    # =====================================================================
    # 3. Train with proxy loss (200 epochs)
    # =====================================================================
    println("\n[3/6] Training proxy loss (200 epochs, lr=0.001)...")
    config = TrainConfig(proxy_epochs=200, proxy_lr=0.001,
                         reinforce_epochs=0, reinforce_lr=0.0001,
                         sigma=0.1, baseline_ema=0.9)

    t2 = time()
    proxy_losses = train_gcn_proxy!(model, g, config)
    train_time = time() - t2
    println("  Initial loss: $(round(proxy_losses[1], digits=2))")
    println("  Final loss:   $(round(proxy_losses[end], digits=2))")
    println("  Loss ratio:   $(round(proxy_losses[end] / proxy_losses[1], digits=4))")
    println("  Training time: $(round(train_time, digits=2))s")

    # =====================================================================
    # 4. Save enhanced GCN model
    # =====================================================================
    println("\n[4/6] Saving enhanced GCN model...")
    mkpath(dirname(MODEL_OUTPUT))
    save_gcn_model(model, MODEL_OUTPUT)
    println("  Saved to: $MODEL_OUTPUT")

    # =====================================================================
    # 5. Apply GCN ordering and compress with existing RL policy
    # =====================================================================
    println("\n[5/6] Compressing with enhanced GCN ordering + existing RL policy...")
    V = eltype(g)
    gcn_map = gcn_ordering(model, V)
    g_reordered = relabel_graph(g, gcn_map)
    println("  Reordered graph: $(nv(g_reordered)) vertices, $(ne(g_reordered)) edges")

    output_base = joinpath(PROJECT_ROOT, "tmp", "cnr2000_gcn_v2")
    mkpath(dirname(output_base))

    t3 = time()
    write_rl_compressed_mgs3_graph(g_reordered, output_base, RL_POLICY, 1;
        coding_scheme=:children, ref_window_size=7, integer_encoding=:fibonacci)
    compress_time = time() - t3

    mgz_file = output_base * ".mgz"
    file_size = filesize(mgz_file)
    bpe = 8.0 * file_size / m
    println("  File size:     $(file_size) bytes")
    println("  Bits/edge:     $(round(bpe, digits=4))")
    println("  Compress time: $(round(compress_time, digits=2))s")

    # =====================================================================
    # 6. Verify roundtrip
    # =====================================================================
    println("\n[6/6] Verifying roundtrip...")
    t4 = time()
    g_loaded = load_compressed_mgs3_graph(mgz_file)
    decompress_time = time() - t4
    n2, m2 = nv(g_loaded), ne(g_loaded)
    println("  Loaded: $n2 vertices, $m2 edges")
    println("  Decompress time: $(round(decompress_time, digits=2))s")

    if m2 == ne(g_reordered)
        println("  Roundtrip: OK ($(m2) edges match)")
    else
        println("  Roundtrip: MISMATCH (expected $(ne(g_reordered)) edges, got $m2)")
    end

    # Cleanup temp file
    rm(mgz_file; force=true)

    # =====================================================================
    # Summary
    # =====================================================================
    println("\n" * "=" ^ 70)
    println("Summary")
    println("=" ^ 70)
    println("  Graph:          CNR-2000 ($n vertices, $m edges)")
    println("  Features:       $(model.input_dim)-dim (degrees, clustering, PageRank, neighbor stats, spectral)")
    println("  Hidden dim:     $(model.hidden_dim)")
    println("  Proxy loss:     $(round(proxy_losses[1], digits=2)) → $(round(proxy_losses[end], digits=2))")
    println("  Final bpe:      $(round(bpe, digits=4))")
    println("  Previous bpe:   7.80 (4-feature GCN)")
    println("  Baseline bpe:   5.19 (LLP ordering)")
    println("  Model saved:    $MODEL_OUTPUT")
    println("=" ^ 70)
end

main()
