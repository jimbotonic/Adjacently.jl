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
#
# Compress enwiki-2013 with CGE K=1 (no clustering, pure reference coding)
#

using Pkg
Pkg.activate(normpath(joinpath(@__DIR__, "..")))

using Logging
global_logger(ConsoleLogger(stderr, Logging.Info))

using LightGraphs: nv, ne, outneighbors, vertices, eltype
using Adjacently
using Adjacently.IO: load_adjacency_list_from_csv
using Adjacently.MGS: write_cge_mgs3_graph, load_compressed_mgs3_graph
using Adjacently.Compression.CGE: CGEParams

const DS = "enwiki-2013"
const DS_DIR = normpath(joinpath(@__DIR__, "..", "datasets", "webgraph", DS))
const DS_CSV = joinpath(DS_DIR, "$DS.csv")
const PREFIX = replace(DS, "-" => "")

isfile(DS_CSV) || error("CSV not found: $DS_CSV")

println("=" ^ 70)
println("CGE compression — $DS  K=1 (no clustering)")
println("=" ^ 70)

# ── Load graph ───────────────────────────────────────────────────────────────

@info "Loading $DS..."
t0 = time()
g = load_adjacency_list_from_csv(DS_CSV, ',', true)
n = nv(g); m = ne(g)
T = eltype(g)
@info "  Loaded: $n vertices, $m edges ($(round(time()-t0, digits=1))s)"

# ── K=1: single cluster with all vertices ─────────────────────────────────────

clusters = [T.(1:n)]

# ── CGE parameters ──────────────────────────────────────────────────────────

cge_params = CGEParams(
    L=128,
    varint=:fibonacci, count_varint=:fibonacci,
    gap=:fibonacci, degree=:elias_delta,
    undirected_pairs=false,
    perm_strategy=:blockpos, membership=:implicit_ranges,
    inter_strategy=:perm,
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
    intra_adapt_mil=2,
)

# ── Progress bar ─────────────────────────────────────────────────────────────

const PROGRESS_INTERVAL = 2.0
last_progress_time = Ref(time())
vertices_done = Ref(0)
t_encode_start = Ref(0.0)

function progress_callback(idx_local::Int, cluster_size::Int, cluster_idx::Int, num_clusters::Int)
    now = time()
    if idx_local == cluster_size || now - last_progress_time[] >= PROGRESS_INTERVAL
        total_done = vertices_done[] + idx_local
        elapsed = now - t_encode_start[]
        rate = total_done / max(elapsed, 0.001)
        remaining = (n - total_done) / max(rate, 0.001)

        pct = round(100.0 * total_done / n, digits=1)
        bar_width = 40
        filled = round(Int, bar_width * total_done / n)
        bar = "█" ^ filled * "░" ^ (bar_width - filled)

        print("\r  [$bar] $pct%  $(round(Int, total_done))/$n vertices  $(round(rate, digits=0)) v/s  ETA $(round(Int, remaining))s   ")
        last_progress_time[] = now

        if idx_local == cluster_size
            vertices_done[] += cluster_size
        end
    end
end

# ── Encode ───────────────────────────────────────────────────────────────────

cge_base = joinpath(DS_DIR, "$(PREFIX)_cge")
cge_mgz = cge_base * ".mgz"

@info "Encoding CGE K=1 → $cge_mgz"
t_encode_start[] = time()

write_cge_mgs3_graph(g, cge_base, clusters;
    params=cge_params, progress=progress_callback)

println()
dt_enc = round(time() - t_encode_start[], digits=1)
bpe = round(8.0 * filesize(cge_mgz) / m, digits=4)
@info "  Encoded: $bpe BPE ($(filesize(cge_mgz)) bytes, $(dt_enc)s)"

# ── Verify roundtrip ─────────────────────────────────────────────────────────

@info "Verifying roundtrip..."
t_dec = time()
g_loaded = load_compressed_mgs3_graph(cge_mgz; cge_params=cge_params)
dt_dec = round(time() - t_dec, digits=1)

if nv(g_loaded) == nv(g) && ne(g_loaded) == ne(g)
    @info "  Roundtrip: PASSED ($(nv(g_loaded)) vertices, $(ne(g_loaded)) edges, decode $(dt_dec)s)"
else
    @warn "  Roundtrip: MISMATCH! Expected $(nv(g))v/$(ne(g))e, got $(nv(g_loaded))v/$(ne(g_loaded))e"
end

# ── Summary ──────────────────────────────────────────────────────────────────

println("\n" * "=" ^ 70)
println("SUMMARY — $DS CGE K=1")
println("=" ^ 70)
println("  CGE K=1:    $bpe BPE  ($(filesize(cge_mgz)) bytes)")
println("  CGE K=6:    16.9861 BPE (previous)")
println("  WebGraph BV: 12.639 BPE (reference)")
println("  Encode time: $(dt_enc)s")
println("  Decode time: $(dt_dec)s")
println("=" ^ 70)
