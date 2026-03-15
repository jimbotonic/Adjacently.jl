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

# Best-known compression roundtrip tests for any WebGraph dataset.
# Tests: BG (Standard greedy), CS (Command Stream), and CG
# Each method: write .mgz → load .mgz → verify graph restored
#
# Usage:
#   julia --project test/run_tests_webgraph_best_compression.jl [DATASET] [METHOD] [K] [MAX_CLUSTER_SIZE] [LLP_PASSES] [COST_MODEL]
#
# DATASET           name of a subdirectory under datasets/webgraph/ that contains
#                   a CSV file named DATASET.csv (e.g. "cnr-2000", "in-2004").
#                   Defaults to "cnr-2000".
#
# METHOD            one of "bg", "cs", "cg", or "all" (default: "all").
#
# K                 top-level Leiden split: K-1 largest communities become
#                   separate clusters, the rest merge into cluster K.
#                   Use "auto" to automatically select K via auto_select_K().
#                   Default: 1 (→ 1 cluster, i.e. largest community + rest).
#
# MAX_CLUSTER_SIZE  optional. If provided, any cluster larger than this is
#                   recursively split via Leiden until all clusters ≤ this size.
#                   Default: 0 (disabled — no recursive splitting).
#
# LLP_PASSES        number of label-propagation passes per resolution level
#                   for LLP reordering (BG/CS). Default: 10.
#
# COST_MODEL        cost estimation model: 0=full (all encoding options),
#                   1=fast (skip RLE, single MIL, simpler copy cost).
#                   Default: 0.
#
# Examples:
#   julia --project test/run_tests_webgraph_best_compression.jl
#   julia --project test/run_tests_webgraph_best_compression.jl cnr-2000
#   julia --project test/run_tests_webgraph_best_compression.jl in-2004 cg auto
#   julia --project test/run_tests_webgraph_best_compression.jl in-2004 cg 4
#   julia --project test/run_tests_webgraph_best_compression.jl in-2004 cg 4 500000
#   julia --project test/run_tests_webgraph_best_compression.jl enwiki-2013 cs 1 0 20

include("run_tests_main.jl")

using Logging
global_logger(ConsoleLogger(stderr, Logging.Info))

using LightGraphs: nv, ne, outneighbors, vertices, eltype
using Adjacently
using Adjacently.IO: load_adjacency_list_from_csv
using Adjacently.Relabeling: relabel_graph, relabel_vertices_llp, save_llp_ordering, load_llp_ordering
using Adjacently.Clustering: leiden_partition, metis_partition, auto_select_K
using Adjacently.Graph: subgraph
using Adjacently.MGS: write_bg_mgs3_graph, write_cs_mgs3_graph, write_cg_mgs3_graph, load_compressed_mgs3_graph
using Adjacently.Compression.CG: CGParams

# ---------------------------------------------------------------------------
# CLI argument parsing
# ---------------------------------------------------------------------------
const DATASET          = length(ARGS) >= 1 ? ARGS[1] : "cnr-2000"
const METHOD           = length(ARGS) >= 2 ? lowercase(ARGS[2]) : "all"
const K_ARG            = length(ARGS) >= 3 ? ARGS[3] : "1"
const AUTO_K           = lowercase(K_ARG) == "auto"
const K                = AUTO_K ? 0 : parse(Int, K_ARG)   # 0 = placeholder, resolved after loading graph
const MAX_CLUSTER_SIZE = length(ARGS) >= 4 ? parse(Int, ARGS[4]) : 0
const LLP_PASSES       = length(ARGS) >= 5 ? parse(Int, ARGS[5]) : 10
const COST_MODEL_ARG   = length(ARGS) >= 6 ? parse(Int, ARGS[6]) : 0  # 0=full, 1=fast

# Sanitised short name for file prefixes (e.g. "cnr-2000" → "cnr2000", "in-2004" → "in2004")
const PREFIX  = replace(DATASET, "-" => "")
const DS_DIR  = normpath(joinpath(@__DIR__, "..", "datasets", "webgraph", DATASET))
const DS_CSV  = joinpath(DS_DIR, DATASET * ".csv")

const RUN_BG  = METHOD in ("all", "bg")
const RUN_CS   = METHOD in ("all", "cs")
const RUN_CG = METHOD in ("all", "cg")

@info "Dataset: $DATASET  (prefix=$PREFIX)"
@info "Method:  $METHOD  (BG=$RUN_BG, CS=$RUN_CS, CG=$RUN_CG)"
@info "CG: K=$(AUTO_K ? "auto" : K), max_cluster_size=$(MAX_CLUSTER_SIZE == 0 ? "disabled" : MAX_CLUSTER_SIZE)"
@info "LLP: passes=$LLP_PASSES"
@info "Cost model: $(COST_MODEL_ARG == 0 ? "full" : "fast") ($COST_MODEL_ARG)"
@info "  CSV:   $DS_CSV"
@info "  Out:   $DS_DIR"

# ---------------------------------------------------------------------------
# Per-dataset parameter dictionaries
# ---------------------------------------------------------------------------
const DATASET_PARAMS = Dict(
    "cnr-2000" => (
        bg_write = (coding_scheme=:children, integer_encoding=:fibonacci,
                     ref_window_size=64, copy_blocks=true,
                     stop_deltas=true, lr_split=false, multi_ref=true),
        cs_write  = (coding_scheme=:children, integer_encoding=:fibonacci,
                     ref_window_size=256, compact_copy=true, tight_intervals=true),
        cg_params = CGParams(
            L=128,
            varint=:fibonacci, count_varint=:fibonacci,
            gap=:fibonacci, degree=:elias_delta,
            undirected_pairs=false,
            perm_strategy=:blockpos, membership=:implicit_ranges,
            inter_strategy=:lists,
            intra_ref_enabled=true, intra_ref_window=64,
            intra_ref_rle=false,
            intra_block_try=false,
            positions_mode=:delta, additions_mode=:delta,
            min_cluster_density=0.0,
            intra_intervals=false, intra_mil=4, intra_greedy_mil=false,
            intra_zigzag=true, intra_stop_deltas=true,
            intra_copy_blocks=true, intra_copy_adaptive=true,
            intra_ref_fixwidth=true, intra_ref_vlc=false,
            intra_add_adaptive=true, intra_raw_adaptive=true,
            intra_adapt_mil=4,
        ),
    ),
    "in-2004" => (
        bg_write = (coding_scheme=:children, integer_encoding=:fibonacci,
                     ref_window_size=64, copy_blocks=true,
                     stop_deltas=true, lr_split=true, multi_ref=true),
        cs_write  = (coding_scheme=:children, integer_encoding=:fibonacci,
                     ref_window_size=64, compact_copy=true, tight_intervals=true),
        cg_params = CGParams(
            L=128,
            varint=:fibonacci, count_varint=:fibonacci,
            gap=:fibonacci, degree=:elias_delta,
            undirected_pairs=false,
            perm_strategy=:blockpos, membership=:implicit_ranges,
            inter_strategy=:lists,
            intra_ref_enabled=true, intra_ref_window=8,
            intra_ref_rle=false,
            intra_block_try=false,
            positions_mode=:delta, additions_mode=:delta,
            min_cluster_density=0.0,
            intra_intervals=false, intra_mil=4, intra_greedy_mil=false,
            intra_zigzag=true, intra_stop_deltas=true,
            intra_copy_blocks=true, intra_copy_adaptive=true,
            intra_ref_fixwidth=true, intra_ref_vlc=false,
            intra_add_adaptive=true, intra_raw_adaptive=true,
            intra_adapt_mil=4,
        ),
    ),
    "enwiki-2013" => (
        bg_write = (coding_scheme=:children, integer_encoding=:zeta,
                     ref_window_size=8, copy_blocks=true,
                     stop_deltas=true, lr_split=true, multi_ref=true),
        cs_write  = (coding_scheme=:children, integer_encoding=:zeta,
                     ref_window_size=8, compact_copy=true, tight_intervals=true),
        cg_params = CGParams(
            L=128,
            varint=:fibonacci, count_varint=:fibonacci,
            gap=:fibonacci, degree=:elias_delta,
            undirected_pairs=false,
            perm_strategy=:blockpos, membership=:implicit_ranges,
            inter_strategy=:lists,
            intra_ref_enabled=true, intra_ref_window=64,
            intra_ref_rle=false,
            intra_block_try=false,
            positions_mode=:delta, additions_mode=:delta,
            min_cluster_density=0.0,
            intra_intervals=true, intra_mil=4, intra_greedy_mil=false,
            intra_zigzag=true, intra_stop_deltas=true,
            intra_copy_blocks=true, intra_copy_adaptive=true,
            intra_ref_fixwidth=true, intra_ref_vlc=false,
            intra_add_adaptive=true, intra_raw_adaptive=true,
            intra_adapt_mil=4,
            intra_lr_split=true,
        ),
    ),
)

# Look up params for this dataset (fall back to cnr-2000 defaults)
const PARAMS = get(DATASET_PARAMS, DATASET, DATASET_PARAMS["cnr-2000"])

# ---------------------------------------------------------------------------
# Helpers
# ---------------------------------------------------------------------------

# Compute or load cached global LLP ordering
function get_global_llp(g, ds_dir::String, prefix::String; passes::Int=10)
    T = eltype(g)
    n = Int(nv(g))
    llp_path = joinpath(ds_dir, "$(prefix)_llp_p$(passes).bin")
    if isfile(llp_path)
        @info "Found cached LLP ordering: $llp_path"
        return load_llp_ordering(llp_path, T)
    else
        @info "Computing global LLP (passes=$passes)..."
        mapping = relabel_vertices_llp(g, :sym; passes=passes)
        save_llp_ordering(llp_path, mapping, n)
        return mapping
    end
end

# Generalized Leiden K-split: K-1 largest clusters + rest merged into cluster K.
# Returns Vector{Vector{T}} of K clusters, each ordered by LLP for locality.
function leiden_partition_k(g, K::Int)
    T = eltype(g); n = nv(g)
    part = leiden_partition(g; max_passes=8, max_levels=5)

    # Count vertices per Leiden label
    counts = Dict{Int,Int}()
    for c in part; counts[c] = get(counts, c, 0) + 1; end

    # Sort labels by count descending
    sorted_labels = sort(collect(keys(counts)); by = l -> counts[l], rev=true)

    # The K-1 largest labels each get their own cluster; everything else merges into cluster K
    top_labels = Set(sorted_labels[1:min(K-1, length(sorted_labels))])

    # Build a mapping: Leiden label → cluster index (1..K)
    label_to_cluster = Dict{Int,Int}()
    ci = 1
    for l in sorted_labels
        if l in top_labels
            label_to_cluster[l] = ci
            ci += 1
        else
            label_to_cluster[l] = K
        end
    end

    # Assign vertices to clusters
    clusters = [T[] for _ in 1:K]
    for i in 1:n
        push!(clusters[label_to_cluster[part[i]]], T(i))
    end

    # Sort each cluster by vertex ID
    for C in clusters; sort!(C); end

    return clusters
end

# Balanced K-way partition via METIS. Returns Vector{Vector{T}} of K clusters.
function metis_partition_k(g, K::Int)
    T = eltype(g); n = nv(g)
    part = metis_partition(g, K)

    clusters = [T[] for _ in 1:K]
    for i in 1:n
        push!(clusters[part[i]], T(i))
    end

    # Sort each cluster by vertex ID
    for C in clusters; sort!(C); end

    return clusters
end

# Recursive Leiden: split large clusters until all are below max_cluster_size.
# Returns Vector{Vector{T}} sorted by cluster size descending.
function recursive_leiden_partition(g, max_cluster_size::Int; max_depth::Int=3)
    T = eltype(g); n = nv(g)
    part = leiden_partition(g; max_passes=8, max_levels=5)

    # Group vertices by Leiden label
    label_groups = Dict{Int,Vector{T}}()
    for i in 1:n
        push!(get!(label_groups, part[i], T[]), T(i))
    end

    final_clusters = Vector{Vector{T}}()

    for (_, verts) in label_groups
        sort!(verts)
        if length(verts) > max_cluster_size && max_depth > 0
            # Extract subgraph and recurse
            sg, _, noi = subgraph(g, verts)
            sub_clusters = recursive_leiden_partition(sg, max_cluster_size; max_depth=max_depth-1)
            # Map subgraph vertex IDs back to original
            for sc in sub_clusters
                push!(final_clusters, [noi[v] for v in sc])
            end
        else
            push!(final_clusters, verts)
        end
    end

    # Sort clusters by size descending, then sort vertices within each
    sort!(final_clusters; by=length, rev=true)
    for C in final_clusters; sort!(C); end

    return final_clusters
end

# Extract sorted neighbor lists from a graph
function extract_neighbor_lists(g)
    T = eltype(g)
    nls = Dict{T,Vector{T}}()
    for v in vertices(g)
        nls[T(v)] = sort([T(o) for o in outneighbors(g, v)])
    end
    return nls
end

# Verify two neighbor-list dicts match
function verify_roundtrip(original_nls::Dict{T,Vector{T}}, decoded_nls::Dict{T,Vector{T}}, m::Int) where {T}
    dec_edges = sum(length(v) for v in values(decoded_nls))
    dec_edges != m && return false
    for (v, orig) in original_nls
        haskey(decoded_nls, v) || return false
        sort(decoded_nls[v]) == sort(orig) || return false
    end
    return true
end

# ============================================================================
# Load and prepare graph
# ============================================================================

if !isfile(DS_CSV)
    @warn "$DATASET CSV not found at $DS_CSV; skipping all compression tests"
    @testset "$DATASET Best-Known Compression" begin
        @test_skip "$DATASET dataset unavailable"
    end
else

@info "Loading $DATASET graph..."
t_load = time()
g_original = load_adjacency_list_from_csv(DS_CSV, ',', true)
n_orig = nv(g_original)
m_orig = ne(g_original)
@info "  Loaded in $(round(time()-t_load, digits=1))s: $n_orig vertices, $m_orig edges"

# Compute (or load cached) ordering for BG/CS/CG-K=1.
# Prefer Leiden+LLP ordering if available (produces best BPE for BG/CS).
# Falls back to plain LLP, or original ordering if LLP_PASSES=0.
leiden_llp_path = joinpath(DS_DIR, "$(PREFIX)_leiden_llp_order.bin")
if isfile(leiden_llp_path)
    @info "Using cached Leiden+LLP ordering: $leiden_llp_path"
    t_rel = time()
    TV = eltype(g_original)
    vertex_map = load_llp_ordering(leiden_llp_path, TV)
    g_rel = relabel_graph(g_original, vertex_map)
    @info "  Leiden+LLP ready in $(round(time()-t_rel, digits=1))s"
elseif LLP_PASSES > 0
    @info "Global LLP ordering (passes=$LLP_PASSES)..."
    t_rel = time()
    vertex_map = get_global_llp(g_original, DS_DIR, PREFIX; passes=LLP_PASSES)
    g_rel = relabel_graph(g_original, vertex_map)
    @info "  LLP ready in $(round(time()-t_rel, digits=1))s"
else
    @info "Skipping global LLP (passes=0), using original ordering"
    g_rel = g_original
end

T_rel = eltype(g_rel)
m = ne(g_rel)

# Free original graph if not needed for CG K>1
# (but never free if g_rel === g_original, i.e. no LLP was applied)
needs_g_original = RUN_CG && (AUTO_K || K > 1)
if !needs_g_original && g_rel !== g_original
    @info "Freeing original graph to save memory..."
    g_original = nothing
end
vertex_map = nothing
GC.gc()

# Extract neighbor lists only if needed for BG/CS roundtrip verification
if RUN_BG || RUN_CS
    nls = extract_neighbor_lists(g_rel)
end

# ============================================================================
# 1. BG (Standard greedy)
# ============================================================================

if RUN_BG
@testset "BG — $DATASET roundtrip" begin
    @info "--- BG: write .mgz → load .mgz → verify ---"

    bg_base = joinpath(DS_DIR, "$(PREFIX)_bg")
    bg_mgz = bg_base * ".mgz"

    # Write .mgz
    t_enc = time()
    write_bg_mgs3_graph(g_rel, bg_base; PARAMS.bg_write..., cost_model=COST_MODEL_ARG)
    dt_enc = round(time() - t_enc, digits=2)

    @test isfile(bg_mgz)
    bpe = round(8.0 * filesize(bg_mgz) / m, digits=4)
    @info "  Encoded: $bpe BPE ($(filesize(bg_mgz)) bytes, $(dt_enc)s)"

    # Load .mgz (v3.1: params auto-decoded from header)
    t_dec = time()
    g_loaded = load_compressed_mgs3_graph(bg_mgz)
    dt_dec = round(time() - t_dec, digits=2)

    # Verify
    @test nv(g_loaded) == nv(g_rel)
    @test ne(g_loaded) == ne(g_rel)
    loaded_nls = extract_neighbor_lists(g_loaded)
    @test verify_roundtrip(nls, loaded_nls, m)
    @info "  BG roundtrip: PASSED ($bpe BPE, decode $(dt_dec)s)"
end
end # RUN_BG

# ============================================================================
# 2. CS (Command Stream)
# ============================================================================

if RUN_CS
@testset "CS — $DATASET roundtrip" begin
    @info "--- CS: write .mgz → load .mgz → verify ---"

    cs_base = joinpath(DS_DIR, "$(PREFIX)_cs")
    cs_mgz = cs_base * ".mgz"

    # Write .mgz
    t_enc = time()
    write_cs_mgs3_graph(g_rel, cs_base; PARAMS.cs_write..., cost_model=COST_MODEL_ARG)
    dt_enc = round(time() - t_enc, digits=2)

    @test isfile(cs_mgz)
    bpe = round(8.0 * filesize(cs_mgz) / m, digits=4)
    @info "  Encoded: $bpe BPE ($(filesize(cs_mgz)) bytes, $(dt_enc)s)"

    # Load .mgz (v3.1: params auto-decoded from header)
    t_dec = time()
    g_loaded = load_compressed_mgs3_graph(cs_mgz)
    dt_dec = round(time() - t_dec, digits=2)

    # Verify
    @test nv(g_loaded) == nv(g_rel)
    @test ne(g_loaded) == ne(g_rel)
    loaded_nls = extract_neighbor_lists(g_loaded)
    @test verify_roundtrip(nls, loaded_nls, m)
    @info "  CS roundtrip: PASSED ($bpe BPE, decode $(dt_dec)s)"
end
end # RUN_CS

# ============================================================================
# 3. CG (Clustered Graph Encoding) — loop over K values
# ============================================================================

if RUN_CG
    # Merge cost_model from CLI into CG params
    cg_params = let p = PARAMS.cg_params
        CGParams(; (fn => getfield(p, fn) for fn in fieldnames(CGParams) if fn != :cost_model)..., cost_model=COST_MODEL_ARG)
    end

    # Resolve K: auto-select or use CLI value
    CG_K = if AUTO_K
        K_auto, _ = auto_select_K(g_original)
        K_auto
    else
        K
    end

    @testset "CG K=$CG_K — $DATASET roundtrip" begin
        @info "--- CG K=$CG_K: write .mgz → load .mgz → verify ---"

        if CG_K == 1
            # K=1: use global LLP-reordered graph with single cluster
            @info "  Using global LLP ordering (K=1)"
            TV = eltype(g_rel)
            n_v = Int(nv(g_rel))
            g_cg = g_rel
            m_cg = ne(g_cg)
            clusters_impl = [TV.(1:n_v)]
        else
            # K>1: Leiden K-split on original graph, then LLP within each cluster
            TV = eltype(g_original)
            n_v = Int(nv(g_original))
            clusters_raw = leiden_partition_k(g_original, CG_K)
            @info "  Leiden K=$CG_K split: $(length.(clusters_raw)) vertices per cluster"

            # Optionally recurse on large clusters
            if MAX_CLUSTER_SIZE > 0
                expanded = Vector{Vector{TV}}()
                for C in clusters_raw
                    if length(C) > MAX_CLUSTER_SIZE
                        sg, _, noi = subgraph(g_original, C)
                        sub = recursive_leiden_partition(sg, MAX_CLUSTER_SIZE)
                        for sc in sub
                            push!(expanded, [noi[v] for v in sc])
                        end
                    else
                        push!(expanded, C)
                    end
                end
                clusters_raw = expanded
                sort!(clusters_raw; by=length, rev=true)
                for C in clusters_raw; sort!(C); end
                @info "  After recursive split (max=$MAX_CLUSTER_SIZE): $(length(clusters_raw)) clusters, top sizes=$(length.(clusters_raw[1:min(10,length(clusters_raw))]))"
            end

            # Apply LLP within each cluster for better intra-cluster ordering
            if LLP_PASSES > 0
                @info "  Applying LLP within each cluster (passes=$LLP_PASSES)..."
                t_llp = time()
                for (ci, C) in enumerate(clusters_raw)
                    length(C) <= 2 && continue
                    sg, oni, _ = subgraph(g_original, C)
                    mapping = relabel_vertices_llp(sg, :sym; passes=LLP_PASSES)
                    sort!(C, by = v -> Int(mapping[oni[v]]))
                    clusters_raw[ci] = C
                end
                @info "  Within-cluster LLP done ($(round(time()-t_llp, digits=1))s)"
            else
                @info "  Skipping within-cluster LLP (passes=0)"
                for C in clusters_raw; sort!(C); end
            end

            # Build vertex mapping: concatenate clusters in order
            cg_vertex_map = let new_id = TV(1)
                d = Dict{TV,TV}()
                for C in clusters_raw
                    for v in C
                        d[v] = new_id
                        new_id += TV(1)
                    end
                end
                d
            end

            g_cg = relabel_graph(g_original, cg_vertex_map)
            m_cg = ne(g_cg)

            # Build implicit-range cluster vectors
            clusters_impl = Vector{Vector{TV}}()
            offset = TV(0)
            for C in clusters_raw
                sz = TV(length(C))
                push!(clusters_impl, TV.(offset+1:offset+sz))
                offset += sz
            end
        end

        mcs_suffix = MAX_CLUSTER_SIZE > 0 ? "_mcs$(MAX_CLUSTER_SIZE)" : ""
        k_suffix = CG_K == 1 ? "" : "_k$CG_K"
        cg_base = joinpath(DS_DIR, "$(PREFIX)_cg$(k_suffix)$(mcs_suffix)")
        cg_mgz = cg_base * ".mgz"

        # Write .mgz
        t_enc = time()
        _last_progress_time = Ref(time())
        function _cg_progress(idx_local::Int, cluster_size::Int, ci::Int, n_clusters::Int)
            now = time()
            if now - _last_progress_time[] >= 10.0 || idx_local == cluster_size
                elapsed = round(now - t_enc, digits=1)
                pct = round(100.0 * idx_local / cluster_size, digits=1)
                println(stderr, "  [CG] cluster $ci/$n_clusters: vertex $idx_local/$cluster_size ($pct%) — $(elapsed)s elapsed")
                flush(stderr)
                _last_progress_time[] = now
            end
        end
        write_cg_mgs3_graph(g_cg, cg_base, clusters_impl; params=cg_params, progress=_cg_progress)
        dt_enc = round(time() - t_enc, digits=2)

        @test isfile(cg_mgz)
        bpe = round(8.0 * filesize(cg_mgz) / m_cg, digits=4)
        @info "  Encoded: $bpe BPE ($(filesize(cg_mgz)) bytes, $(dt_enc)s)"

        # Load .mgz (v3.1: params auto-decoded from header)
        t_dec = time()
        g_loaded = load_compressed_mgs3_graph(cg_mgz)
        dt_dec = round(time() - t_dec, digits=2)

        # Verify (lightweight: compare vertex-by-vertex without building full dicts)
        @test nv(g_loaded) == nv(g_cg)
        @test ne(g_loaded) == ne(g_cg)
        match = true
        for v in vertices(g_cg)
            orig = sort(collect(outneighbors(g_cg, v)))
            decoded = sort(collect(outneighbors(g_loaded, v)))
            if orig != decoded
                match = false
                @warn "  Mismatch at vertex $v: expected $(length(orig)) neighbors, got $(length(decoded))"
                break
            end
        end
        @test match
        @info "  CG roundtrip: PASSED ($bpe BPE, decode $(dt_dec)s)"
    end
end # RUN_CG

# ============================================================================
# Summary
# ============================================================================

@info "=== $DATASET Compression Summary ==="
if RUN_BG || RUN_CS
    for (label, suffix) in [("BG", "bg"), ("CS", "cs")]
        (label == "BG" && !RUN_BG) && continue
        (label == "CS" && !RUN_CS) && continue
        path = joinpath(DS_DIR, "$(PREFIX)_$(suffix).mgz")
        if isfile(path)
            sz = filesize(path)
            bpe = round(8.0 * sz / m, digits=4)
            @info "  $label: $bpe BPE  ($sz bytes)"
        end
    end
end
if RUN_CG
    mcs_suffix = MAX_CLUSTER_SIZE > 0 ? "_mcs$(MAX_CLUSTER_SIZE)" : ""
    k_suffix = K == 1 ? "" : "_k$K"
    path = joinpath(DS_DIR, "$(PREFIX)_cg$(k_suffix)$(mcs_suffix).mgz")
    if isfile(path)
        sz = filesize(path)
        bpe = round(8.0 * sz / m_orig, digits=4)
        @info "  CG K=$K: $bpe BPE  ($sz bytes)"
    end
end

end # if isfile(DS_CSV)
