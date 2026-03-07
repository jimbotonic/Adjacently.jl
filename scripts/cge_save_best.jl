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
# Save CNR-2000 compressed with best known CGE parameters.
# Best known: 2.4341 BPE (3-way adaptive copy, w=64, zigzag, stop-deltas, adaptive adds/raw)
#
# Output: datasets/webgraph/cnr-2000/cnr2000_cge_best.cge
#

using Pkg
Pkg.activate(normpath(joinpath(@__DIR__, "..")))

using Logging
global_logger(ConsoleLogger(stderr, Logging.Info))

using LightGraphs
using LightGraphs: nv, outneighbors
using Adjacently
using Adjacently.CGE: encode_level, CGEParams, CGEStats, decode_level
using Adjacently.IO: load_adjacency_list_from_csv, BitWriter, BitReader, flush_bitwriter
using Adjacently.Clustering: leiden_partition
using Adjacently.Relabeling: relabel_graph

count_edges(h) = sum(length(outneighbors(h, v)) for v in 1:nv(h))

# Best known parameters
best_params = CGEParams(
    L=128,
    varint=:fibonacci,
    count_varint=:fibonacci,
    gap=:fibonacci,
    degree=:elias_delta,
    undirected_pairs=false,
    perm_strategy=:blockpos,
    membership=:implicit_ranges,
    inter_strategy=:perm,
    intra_ref_enabled=true,
    intra_ref_window=64,
    intra_ref_rle=false,
    intra_block_try=false,
    positions_mode=:delta,
    additions_mode=:delta,
    min_cluster_density=0.0,
    intra_intervals=false,
    intra_mil=4,
    intra_greedy_mil=false,
    intra_zigzag=true,
    intra_stop_deltas=true,
    intra_copy_blocks=true,
    intra_copy_adaptive=true,   # 3-way: bitmap / copy-blocks / complement
    intra_ref_fixwidth=true,
    intra_ref_vlc=false,
    intra_add_adaptive=true,
    intra_raw_adaptive=true,
    intra_adapt_mil=2,
)

println("=" ^ 70)
println("CGE Best-Params Save: CNR-2000")
println("=" ^ 70)
println("Parameters:")
println("  membership:      implicit_ranges")
println("  intra_ref_window: 64")
println("  intra_zigzag:    true")
println("  intra_stop_deltas: true")
println("  intra_copy_adaptive: true  (3-way: bitmap/cb/complement)")
println("  intra_ref_fixwidth: true")
println("  intra_add_adaptive: true")
println("  intra_raw_adaptive: true")
println("  intra_adapt_mil: 2")
println()

cnr_csv = normpath(joinpath(@__DIR__, "..", "datasets", "webgraph", "cnr-2000", "cnr-2000.csv"))
isfile(cnr_csv) || error("CNR-2000 CSV not found at $cnr_csv")

out_path = normpath(joinpath(@__DIR__, "..", "datasets", "webgraph", "cnr-2000", "cnr2000_cge_best.cge"))

println("Loading CNR-2000...")
t0 = time()
g  = load_adjacency_list_from_csv(cnr_csv, ',', true)
n  = nv(g)
m  = count_edges(g)
TV = eltype(g)
println("  Loaded: $n vertices, $m edges ($(round(time()-t0, digits=2))s)")

println("\nBuilding pre-relabeled graph (Leiden K=1, vertex-ID rank relabeling)...")
t_rel = time()
part = leiden_partition(g; max_passes=8, max_levels=5)
counts = Dict{Int,Int}()
for c in part; counts[c] = get(counts, c, 0) + 1; end
top_label = argmax(counts)

clusters_raw = [TV[], TV[]]
for i in 1:n
    push!(clusters_raw[part[i] == top_label ? 1 : 2], TV(i))
end
clusters_by_id = [sort(C) for C in clusters_raw]
S1 = length(clusters_by_id[1])

vertex_map = let new_id = TV(1)
    d = Dict{TV,TV}()
    for C in clusters_by_id
        for v in C
            d[v] = new_id
            new_id += TV(1)
        end
    end
    d
end
g_rel = relabel_graph(g, vertex_map)
clusters_impl = [TV.(1:S1), TV.(S1+1:n)]
println("  Clusters: $S1 + $(n-S1) ($(round(time()-t_rel, digits=2))s)")

println("\nEncoding...")
t_enc = time()
io = IOBuffer()
w = BitWriter(io)
stats = CGEStats()
encode_level(w, g_rel, clusters_impl; params=best_params, stats=stats)
flush_bitwriter(w; flush_last_bits=true)
bytes = take!(io)
dt_enc = round(time() - t_enc, digits=2)

nbytes = length(bytes)
bpe    = round(8.0 * nbytes / m, digits=4)
println("  Encoded: $bpe BPE  ($nbytes bytes, $(dt_enc)s)")
println("  hdrs=$(ceil(Int,stats.bits_intra_headers/8))B  copy=$(ceil(Int,stats.bits_intra_copy/8))B  add=$(ceil(Int,stats.bits_intra_add/8))B  raw=$(ceil(Int,stats.bits_intra_raw/8))B")

bm  = stats.intra_copy_bitmap_count
cb  = stats.intra_copy_blocks_count
cc  = stats.intra_copy_complement_count
tot = bm + cb + cc
if tot > 0
    println("  Copy modes ($tot ref vertices): bitmap=$(bm) ($(round(100bm/tot,digits=1))%), cb=$(cb) ($(round(100cb/tot,digits=1))%), complement=$(cc) ($(round(100cc/tot,digits=1))%)")
end

println("\nVerifying roundtrip...")
t_dec = time()
r_v = BitReader(IOBuffer(bytes))
decoded = decode_level(r_v, best_params; T=TV, directed=true)
dec_edges = sum(length(v) for v in values(decoded))
dt_dec = round(time() - t_dec, digits=2)

if dec_edges == m
    println("  Roundtrip: OK ($dec_edges edges, $(dt_dec)s)")
else
    error("  Roundtrip MISMATCH: expected $m edges, got $dec_edges!")
end

println("\nSaving to: $out_path")
open(out_path, "w") do f
    write(f, bytes)
end
println("  Saved: $nbytes bytes")

println("\n" * "=" ^ 70)
println("RESULT")
println("=" ^ 70)
println("  BPE:           $bpe")
println("  Size:          $nbytes bytes  ($(round(nbytes/1024, digits=1)) KB)")
println("  Vertices:      $n")
println("  Edges:         $m")
println("  Output:        $out_path")
println("  Roundtrip:     OK")
println("=" ^ 70)
