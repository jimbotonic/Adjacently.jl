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
# Load Web_Google core & rcore, predict CGE params with trained GNN
# (on subgraphs), grid-search for best params, further train GNN,
# save best .mgz files.
#
# Usage:
#   julia --project scripts/compress_web_google.jl
#

using Pkg
Pkg.activate(normpath(joinpath(@__DIR__, "..")))

using Logging
global_logger(ConsoleLogger(stderr, Logging.Info))
using Printf

using LightGraphs: AbstractGraph, nv, ne, outneighbors, vertices, eltype, add_edge!
using Adjacently
using Adjacently.CustomLightGraphs: SimpleDiGraph
using Adjacently.MGS: load_mgs3_graph, write_cge_mgs3_graph, load_compressed_mgs3_graph
using Adjacently.Relabeling: relabel_graph
using Adjacently.GNN: GNNModel, gnn_forward
using Adjacently.ParamPrediction: ParamPredictor, MasterGNNWeights,
    predict_params, graph_level_readout, compute_bpe,
    grid_search_best_params, generate_training_labels,
    actions_to_cge_params, save_param_predictor, load_param_predictor,
    train_param_predictor!, PARAM_HEADS
using Adjacently.GraphGenerators: generate_training_batch
using Adjacently.Compression.CGE: CGEParams
using Adjacently.Clustering: leiden_partition
using Random: MersenneTwister

const PROJECT_ROOT = normpath(joinpath(@__DIR__, ".."))
const DS_DIR = joinpath(PROJECT_ROOT, "datasets", "Web_Google")
const CORE_MGS = joinpath(DS_DIR, "Web_Google_core.mgs")
const RCORE_MGS = joinpath(DS_DIR, "Web_Google_rcore.mgs")
const MODEL_PATH = joinpath(PROJECT_ROOT, "models", "param_predictor_trained.bson")

fp(args...) = (println(args...); flush(stdout))

function sample_bfs_subgraph(g::AbstractGraph{T}, target_n::Int;
                              seed::Int=42) where {T<:Unsigned}
    rng = MersenneTwister(seed)
    n = Int(nv(g))
    target_n = min(target_n, n)
    start = T(rand(rng, 1:n))
    visited = Set{T}()
    queue = T[start]
    push!(visited, start)
    order = T[]
    while !isempty(queue) && length(order) < target_n
        v = popfirst!(queue)
        push!(order, v)
        for u in outneighbors(g, v)
            if !(u in visited) && length(visited) < target_n
                push!(visited, u)
                push!(queue, u)
            end
        end
    end
    n_sub = length(order)
    old_to_new = Dict{T, UInt32}()
    for (i, v) in enumerate(order)
        old_to_new[v] = UInt32(i)
    end
    g_sub = SimpleDiGraph{UInt32}(UInt32(n_sub))
    for v in order
        new_v = old_to_new[v]
        for u in outneighbors(g, v)
            if haskey(old_to_new, u)
                add_edge!(g_sub, new_v, old_to_new[u])
            end
        end
    end
    return g_sub
end

function extract_neighbor_lists(g)
    T = eltype(g)
    nls = Dict{T,Vector{T}}()
    for v in vertices(g)
        nls[T(v)] = sort([T(o) for o in outneighbors(g, v)])
    end
    return nls
end

"""
    relabel_for_cge(g, clusters) → (g_relabeled, clusters_impl)

Relabel graph so that clusters become contiguous ID ranges.
Required for membership=:implicit_ranges to work correctly.
Returns the relabeled graph and new contiguous cluster vectors.
"""
function relabel_for_cge(g, clusters)
    TV = eltype(g)
    vertex_map = Dict{TV,TV}()
    new_id = TV(1)
    for C in clusters
        for v in C
            vertex_map[v] = new_id
            new_id += TV(1)
        end
    end
    g_rel = relabel_graph(g, vertex_map)
    # Build implicit-range cluster vectors
    clusters_impl = Vector{Vector{TV}}()
    offset = TV(0)
    for C in clusters
        sz = TV(length(C))
        push!(clusters_impl, TV.(offset+1:offset+sz))
        offset += sz
    end
    return g_rel, clusters_impl
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

function main()
    fp("=" ^ 72)
    fp("  Web_Google CGE Compression Pipeline")
    fp("=" ^ 72)

    # =====================================================================
    # 1. Load graphs
    # =====================================================================
    fp("\n[1/7] Loading Web_Google core and rcore from .mgs files...")

    t0 = time()
    g_core = load_mgs3_graph(CORE_MGS)
    fp("  Core:  $(nv(g_core)) vertices, $(ne(g_core)) edges ($(round(time()-t0, digits=1))s)")

    t0 = time()
    g_rcore = load_mgs3_graph(RCORE_MGS)
    fp("  RCore: $(nv(g_rcore)) vertices, $(ne(g_rcore)) edges ($(round(time()-t0, digits=1))s)")

    ne_core = ne(g_core)
    ne_rcore = ne(g_rcore)

    # =====================================================================
    # 2. Load trained GNN predictor & evaluate on subgraphs
    # =====================================================================
    fp("\n[2/7] GNN prediction on subgraphs (full graphs too large for GNN)...")

    predictor, master = load_param_predictor(MODEL_PATH)
    fp("  Loaded model: hidden_dim=$(master.hidden_dim)")

    # Sample 2000-vertex subgraphs for GNN evaluation
    sub_core = sample_bfs_subgraph(g_core, 2000; seed=301)
    sub_rcore = sample_bfs_subgraph(g_rcore, 2000; seed=401)
    fp("  Core subgraph:  $(nv(sub_core))v, $(ne(sub_core))e")
    fp("  RCore subgraph: $(nv(sub_rcore))v, $(ne(sub_rcore))e")

    for (name, sub) in [("Core", sub_core), ("RCore", sub_rcore)]
        gnn = GNNModel(sub; hidden_dim=master.hidden_dim, seed=42,
                        feature_set=:enhanced, use_gat=master.use_gat)
        Adjacently.ParamPrediction.copy_master_to_gnn!(gnn, master)
        gnn_forward(gnn)
        h_graph = graph_level_readout(gnn)
        actions, params = predict_params(predictor, h_graph)

        T_sub = eltype(sub)
        clusters_sub = [collect(T_sub(1):T_sub(nv(sub)))]
        bpe = compute_bpe(sub, params, clusters_sub)
        fp("\n  $name subgraph — GNN prediction: $(round(bpe, digits=4)) BPE")
        fp("    window=$(params.intra_ref_window), intervals=$(params.intra_intervals)")
        fp("    lr_split=$(params.intra_lr_split), mil=$(params.intra_mil)")
        fp("    stop_deltas=$(params.intra_stop_deltas), zigzag=$(params.intra_zigzag)")
        fp("    copy_adaptive=$(params.intra_copy_adaptive)")
    end

    # =====================================================================
    # 3. Grid search on subgraphs
    # =====================================================================
    fp("\n[3/7] Grid-searching best params on subgraphs (512 combos)...")

    t0 = time()
    best_ac_core, best_bpe_core_sub, _ = grid_search_best_params(sub_core; verbose=false)
    fp("  Core subgraph grid search:  $(round(best_bpe_core_sub, digits=4)) BPE ($(round(time()-t0, digits=1))s)")

    t0 = time()
    best_ac_rcore, best_bpe_rcore_sub, _ = grid_search_best_params(sub_rcore; verbose=false)
    fp("  RCore subgraph grid search: $(round(best_bpe_rcore_sub, digits=4)) BPE ($(round(time()-t0, digits=1))s)")

    dummy_pred = ParamPredictor(1; lr=0.0)
    gs_params_core = actions_to_cge_params(dummy_pred, best_ac_core)
    gs_params_rcore = actions_to_cge_params(dummy_pred, best_ac_rcore)

    fp("\n  Best params (Core sub):  w=$(gs_params_core.intra_ref_window) int=$(gs_params_core.intra_intervals) lr=$(gs_params_core.intra_lr_split) mil=$(gs_params_core.intra_mil) sd=$(gs_params_core.intra_stop_deltas) zz=$(gs_params_core.intra_zigzag) ca=$(gs_params_core.intra_copy_adaptive)")
    fp("  Best params (RCore sub): w=$(gs_params_rcore.intra_ref_window) int=$(gs_params_rcore.intra_intervals) lr=$(gs_params_rcore.intra_lr_split) mil=$(gs_params_rcore.intra_mil) sd=$(gs_params_rcore.intra_stop_deltas) zz=$(gs_params_rcore.intra_zigzag) ca=$(gs_params_rcore.intra_copy_adaptive)")

    # =====================================================================
    # 4. Evaluate on full graphs with K=1
    # =====================================================================
    fp("\n[4/7] Evaluating best params on full graphs (K=1)...")

    TV_core = eltype(g_core)
    TV_rcore = eltype(g_rcore)
    clusters_core_k1 = [collect(TV_core(1):TV_core(nv(g_core)))]
    clusters_rcore_k1 = [collect(TV_rcore(1):TV_rcore(nv(g_rcore)))]

    baseline_params = CGEParams(;
        intra_ref_window=64, intra_intervals=false, intra_lr_split=false,
        intra_mil=4, intra_adapt_mil=4, intra_copy_adaptive=true, intra_stop_deltas=true,
        intra_zigzag=true, intra_copy_blocks=true, intra_ref_fixwidth=true,
        intra_add_adaptive=true, intra_raw_adaptive=true,
        varint=:fibonacci, count_varint=:fibonacci, gap=:fibonacci,
        degree=:elias_delta, membership=:implicit_ranges,
    )

    t0 = time()
    bpe_core_k1_gs = compute_bpe(g_core, gs_params_core, clusters_core_k1)
    fp("  Core K=1 (grid-search):  $(round(bpe_core_k1_gs, digits=4)) BPE ($(round(time()-t0, digits=1))s)")

    t0 = time()
    bpe_core_k1_bl = compute_bpe(g_core, baseline_params, clusters_core_k1)
    fp("  Core K=1 (baseline):     $(round(bpe_core_k1_bl, digits=4)) BPE ($(round(time()-t0, digits=1))s)")

    t0 = time()
    bpe_rcore_k1_gs = compute_bpe(g_rcore, gs_params_rcore, clusters_rcore_k1)
    fp("  RCore K=1 (grid-search): $(round(bpe_rcore_k1_gs, digits=4)) BPE ($(round(time()-t0, digits=1))s)")

    t0 = time()
    bpe_rcore_k1_bl = compute_bpe(g_rcore, baseline_params, clusters_rcore_k1)
    fp("  RCore K=1 (baseline):    $(round(bpe_rcore_k1_bl, digits=4)) BPE ($(round(time()-t0, digits=1))s)")

    # =====================================================================
    # 5. Try K=2 with Leiden clustering
    # =====================================================================
    fp("\n[5/7] Testing K=2 Leiden clustering...")

    function make_k2_clusters(g, K=2)
        TV = eltype(g)
        n_orig = Int(nv(g))
        part = leiden_partition(g; max_passes=8, max_levels=5)
        counts = Dict{Int,Int}()
        for c in part; counts[c] = get(counts, c, 0) + 1; end
        sorted_labels = sort(collect(keys(counts)); by = l -> counts[l], rev=true)
        cls = [TV[] for _ in 1:K]
        l2c = Dict{Int,Int}()
        ci = 1
        for l in sorted_labels
            if ci < K && l in Set(sorted_labels[1:K-1])
                l2c[l] = ci; ci += 1
            else
                l2c[l] = K
            end
        end
        for i in 1:n_orig
            push!(cls[l2c[part[i]]], TV(i))
        end
        for C in cls; sort!(C); end
        return cls
    end

    t0 = time()
    clusters_core_k2 = make_k2_clusters(g_core)
    fp("  Core K=2: $(join(length.(clusters_core_k2), " + ")) ($(round(time()-t0, digits=1))s)")

    t0 = time()
    clusters_rcore_k2 = make_k2_clusters(g_rcore)
    fp("  RCore K=2: $(join(length.(clusters_rcore_k2), " + ")) ($(round(time()-t0, digits=1))s)")

    t0 = time()
    bpe_core_k2_gs = compute_bpe(g_core, gs_params_core, clusters_core_k2)
    fp("  Core K=2 (grid-search):  $(round(bpe_core_k2_gs, digits=4)) BPE ($(round(time()-t0, digits=1))s)")

    t0 = time()
    bpe_core_k2_bl = compute_bpe(g_core, baseline_params, clusters_core_k2)
    fp("  Core K=2 (baseline):     $(round(bpe_core_k2_bl, digits=4)) BPE ($(round(time()-t0, digits=1))s)")

    t0 = time()
    bpe_rcore_k2_gs = compute_bpe(g_rcore, gs_params_rcore, clusters_rcore_k2)
    fp("  RCore K=2 (grid-search): $(round(bpe_rcore_k2_gs, digits=4)) BPE ($(round(time()-t0, digits=1))s)")

    t0 = time()
    bpe_rcore_k2_bl = compute_bpe(g_rcore, baseline_params, clusters_rcore_k2)
    fp("  RCore K=2 (baseline):    $(round(bpe_rcore_k2_bl, digits=4)) BPE ($(round(time()-t0, digits=1))s)")

    # Pick overall best for each graph
    core_candidates = [
        ("K=1 grid-search", bpe_core_k1_gs, gs_params_core, clusters_core_k1, 1),
        ("K=1 baseline",    bpe_core_k1_bl, baseline_params, clusters_core_k1, 1),
        ("K=2 grid-search", bpe_core_k2_gs, gs_params_core,  clusters_core_k2, 2),
        ("K=2 baseline",    bpe_core_k2_bl, baseline_params, clusters_core_k2, 2),
    ]
    rcore_candidates = [
        ("K=1 grid-search", bpe_rcore_k1_gs, gs_params_rcore, clusters_rcore_k1, 1),
        ("K=1 baseline",    bpe_rcore_k1_bl, baseline_params,  clusters_rcore_k1, 1),
        ("K=2 grid-search", bpe_rcore_k2_gs, gs_params_rcore,  clusters_rcore_k2, 2),
        ("K=2 baseline",    bpe_rcore_k2_bl, baseline_params,  clusters_rcore_k2, 2),
    ]

    best_core = sort(core_candidates; by=x->x[2])[1]
    best_rcore = sort(rcore_candidates; by=x->x[2])[1]

    fp("\n  === Best Core:  $(best_core[1]) → $(round(best_core[2], digits=4)) BPE ===")
    fp("  === Best RCore: $(best_rcore[1]) → $(round(best_rcore[2], digits=4)) BPE ===")

    # =====================================================================
    # 6. Further train GNN with Web_Google subgraphs
    # =====================================================================
    fp("\n[6/7] Further training GNN with Web_Google subgraphs...")

    all_subs = SimpleDiGraph{UInt32}[sub_core, sub_rcore]

    fp("  Generating grid-search labels...")
    t0 = time()
    new_labels = generate_training_labels(all_subs; verbose=true)
    fp("  Label time: $(round(time()-t0, digits=1))s")

    fp("  Training for 30 additional epochs...")
    t0 = time()
    loss_history = train_param_predictor!(master, predictor, all_subs, new_labels;
                                           epochs=30, verbose=true,
                                           feature_set=:enhanced,
                                           gnn_lr=0.0001, seed=2000)
    fp("  Training time: $(round(time()-t0, digits=1))s")

    save_param_predictor(predictor, master, MODEL_PATH)
    fp("  Updated model saved to: $MODEL_PATH")

    # =====================================================================
    # 7. Save best .mgz files and verify roundtrips
    # =====================================================================
    fp("\n[7/7] Saving best compressed .mgz files...")

    final_results = Dict{String, Any}()

    for (label, g, best_tuple, ne_g) in [
        ("core", g_core, best_core, ne_core),
        ("rcore", g_rcore, best_rcore, ne_rcore)
    ]
        config_name, bpe, params, clusters, K = best_tuple
        base = joinpath(DS_DIR, "Web_Google_$(label)_cge_k$(K)")
        mgz_file = base * ".mgz"

        fp("\n  --- $label ($config_name, K=$K) ---")

        # Relabel graph for implicit_ranges when K>1 (clusters must be contiguous)
        if K > 1
            g_enc, clusters_enc = relabel_for_cge(g, clusters)
            fp("  Relabeled for K=$K implicit_ranges encoding")
        else
            g_enc = g
            clusters_enc = clusters
        end

        t0 = time()
        write_cge_mgs3_graph(g_enc, base, clusters_enc;
            coding_scheme=:children, integer_encoding=:fibonacci, params=params)
        sz = filesize(mgz_file)
        actual_bpe = round(8.0 * sz / ne_g, digits=4)
        fp("  Encoded: $actual_bpe BPE ($sz bytes, $(round(time()-t0, digits=1))s)")

        t0 = time()
        g_loaded = load_compressed_mgs3_graph(mgz_file)
        fp("  Decoded: $(nv(g_loaded))v, $(ne(g_loaded))e ($(round(time()-t0, digits=1))s)")

        orig_nls = extract_neighbor_lists(g_enc)
        loaded_nls = extract_neighbor_lists(g_loaded)
        ok = verify_roundtrip(orig_nls, loaded_nls, ne_g)
        fp("  Roundtrip: $(ok ? "PASSED" : "FAILED")")

        final_results[label] = (actual_bpe, sz, K, config_name, params, ok)
    end

    # Print data for JSON update
    fp("\n  === JSON update data ===")
    for (label, g, best_tuple, ne_g) in [
        ("core", g_core, best_core, ne_core),
        ("rcore", g_rcore, best_rcore, ne_rcore)
    ]
        config_name, bpe, params, clusters, K = best_tuple
        mgz_file = joinpath(DS_DIR, "Web_Google_$(label)_cge_k$(K).mgz")
        actual_bpe = round(8.0 * filesize(mgz_file) / ne_g, digits=4)
        nv_g = Int(nv(g))
        avg_deg = round(ne_g / nv_g, digits=3)

        fp("\n  web-google-$(label):")
        fp("    K=$K, nv=$nv_g, ne=$ne_g, avg_deg=$avg_deg, best_bpe=$actual_bpe")
        fp("    window=$(params.intra_ref_window) intervals=$(params.intra_intervals) lr_split=$(params.intra_lr_split)")
        fp("    mil=$(params.intra_mil) adapt_mil=$(params.intra_adapt_mil)")
        fp("    stop_deltas=$(params.intra_stop_deltas) zigzag=$(params.intra_zigzag)")
        fp("    copy_adaptive=$(params.intra_copy_adaptive) copy_blocks=$(params.intra_copy_blocks)")
        fp("    ref_fixwidth=$(params.intra_ref_fixwidth) ref_vlc=$(params.intra_ref_vlc)")
        fp("    add_adaptive=$(params.intra_add_adaptive) raw_adaptive=$(params.intra_raw_adaptive)")
        fp("    ref_enabled=$(params.intra_ref_enabled)")
        fp("    varint=$(params.varint) degree=$(params.degree) membership=$(params.membership)")
    end

    # Summary table
    fp("\n" * "=" ^ 72)
    fp("  Final Results Summary")
    fp("=" ^ 72)
    fp()
    @printf("  %-14s %8s %8s %10s %10s  %s\n", "Graph", "Vertices", "Edges", "BPE", "Size", "Config")
    flush(stdout)
    fp("  " * "-" ^ 66)
    for label in ["core", "rcore"]
        bpe, sz, K, config, params, ok = final_results[label]
        g_ref = label == "core" ? g_core : g_rcore
        ne_g = label == "core" ? ne_core : ne_rcore
        @printf("  %-14s %8d %8d %10.4f %10d  %s\n",
            label, Int(nv(g_ref)), ne_g, bpe, sz, config)
        flush(stdout)
    end
    fp()
    fp("  Files saved in: $DS_DIR")
    fp("=" ^ 72)
end

main()
