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
# Compare children vs index coding scheme BPE for STD, CS, CGE on CNR-2000.
# Saves compressed files with descriptive names and prints a summary table.
#
# Usage:
#   julia --project scripts/compare_children_vs_index.jl

using Pkg
Pkg.activate(normpath(joinpath(@__DIR__, "..")))

using Logging
global_logger(ConsoleLogger(stderr, Logging.Info))
using Printf

using LightGraphs: nv, ne, outneighbors, vertices, eltype
using Adjacently
using Adjacently.IO: load_adjacency_list_from_csv
using Adjacently.Clustering: leiden_partition
using Adjacently.MGS: write_std_mgs3_graph, write_cs_mgs3_graph, write_cge_mgs3_graph, load_compressed_mgs3_graph
using Adjacently.Compression.CGE: CGEParams

const DS_DIR = normpath(joinpath(@__DIR__, "..", "datasets", "webgraph", "cnr-2000"))
const DS_CSV = joinpath(DS_DIR, "cnr-2000.csv")

# ---------------------------------------------------------------------------
# Best known parameters for cnr-2000
# ---------------------------------------------------------------------------

const STD_PARAMS = (
    integer_encoding=:fibonacci,
    ref_window_size=64, copy_blocks=true, adaptive_copy=true,
    stop_deltas=true, empty_prefix=true, compact_copy=true,
    tight_intervals=true, vlc2=true,
)

const CS_PARAMS = (
    integer_encoding=:fibonacci,
    ref_window_size=64, compact_copy=true, tight_intervals=true,
)

const CGE_K = 2

const CGE_PARAMS = CGEParams(
    L=128,
    varint=:fibonacci, count_varint=:fibonacci,
    gap=:fibonacci, degree=:elias_delta,
    undirected_pairs=false,
    perm_strategy=:blockpos, membership=:implicit_ranges,
    inter_strategy=:perm,
    intra_ref_enabled=true, intra_ref_window=64,
    intra_ref_rle=false, intra_block_try=false,
    positions_mode=:delta, additions_mode=:delta,
    min_cluster_density=0.0,
    intra_intervals=false, intra_mil=4, intra_greedy_mil=false,
    intra_zigzag=true, intra_stop_deltas=true,
    intra_copy_blocks=true, intra_copy_adaptive=true,
    intra_ref_fixwidth=true, intra_ref_vlc=false,
    intra_add_adaptive=true, intra_raw_adaptive=true,
    intra_adapt_mil=2,
)

# ---------------------------------------------------------------------------
# Load graph
# ---------------------------------------------------------------------------

@info "Loading cnr-2000 graph..."
t0 = time()
g = load_adjacency_list_from_csv(DS_CSV, ',', true)
n_v = nv(g)
m = ne(g)
@info "  Loaded in $(round(time()-t0, digits=1))s: $n_v vertices, $m edges"

# ---------------------------------------------------------------------------
# Prepare CGE clusters (Leiden K=2)
# ---------------------------------------------------------------------------

@info "Computing Leiden K=$CGE_K partition..."
TV = eltype(g)
part = leiden_partition(g; max_passes=8, max_levels=5)
counts = Dict{Int,Int}()
for c in part; counts[c] = get(counts, c, 0) + 1; end
sorted_labels = sort(collect(keys(counts)); by = l -> counts[l], rev=true)
top_label = sorted_labels[1]

clusters = let
    cls = [TV[] for _ in 1:CGE_K]
    l2c = Dict{Int,Int}()
    ci = 1
    for l in sorted_labels
        if ci < CGE_K && l in Set(sorted_labels[1:CGE_K-1])
            l2c[l] = ci; ci += 1
        else
            l2c[l] = CGE_K
        end
    end
    for i in 1:Int(n_v)
        push!(cls[l2c[part[i]]], TV(i))
    end
    for C in cls; sort!(C); end
    cls
end
@info "  Clusters: $(join(length.(clusters), " + ")) vertices"

# ---------------------------------------------------------------------------
# Encode & measure
# ---------------------------------------------------------------------------

results = Dict{String, NamedTuple{(:file, :bytes, :bpe), Tuple{String, Int, Float64}}}()

function encode_and_measure(label, mgz_path)
    bytes = filesize(mgz_path)
    bpe = round(8.0 * bytes / m, digits=4)
    results[label] = (file=mgz_path, bytes=bytes, bpe=bpe)
    @info "  $label: $bpe BPE ($bytes bytes)"
end

# --- STD children ---
@info "Encoding STD children..."
t = time()
p = joinpath(DS_DIR, "cnr2000_std_children")
write_std_mgs3_graph(g, p; coding_scheme=:children, STD_PARAMS...)
@info "  $(round(time()-t, digits=1))s"
encode_and_measure("STD children", p * ".mgz")

# --- STD index ---
@info "Encoding STD index..."
t = time()
p = joinpath(DS_DIR, "cnr2000_std_index")
write_std_mgs3_graph(g, p; coding_scheme=:index, STD_PARAMS...)
@info "  $(round(time()-t, digits=1))s"
encode_and_measure("STD index", p * ".mgz")

# --- CS children ---
@info "Encoding CS children..."
t = time()
p = joinpath(DS_DIR, "cnr2000_cs_children")
write_cs_mgs3_graph(g, p; coding_scheme=:children, CS_PARAMS...)
@info "  $(round(time()-t, digits=1))s"
encode_and_measure("CS children", p * ".mgz")

# --- CS index ---
@info "Encoding CS index..."
t = time()
p = joinpath(DS_DIR, "cnr2000_cs_index")
write_cs_mgs3_graph(g, p; coding_scheme=:index, CS_PARAMS...)
@info "  $(round(time()-t, digits=1))s"
encode_and_measure("CS index", p * ".mgz")

# --- CGE children ---
@info "Encoding CGE K=$CGE_K children..."
t = time()
p = joinpath(DS_DIR, "cnr2000_cge_k$(CGE_K)_children")
write_cge_mgs3_graph(g, p, clusters; coding_scheme=:children, integer_encoding=:fibonacci, params=CGE_PARAMS)
@info "  $(round(time()-t, digits=1))s"
encode_and_measure("CGE children", p * ".mgz")

# --- CGE index ---
@info "Encoding CGE K=$CGE_K index..."
t = time()
p = joinpath(DS_DIR, "cnr2000_cge_k$(CGE_K)_index")
write_cge_mgs3_graph(g, p, clusters; coding_scheme=:index, integer_encoding=:fibonacci, params=CGE_PARAMS)
@info "  $(round(time()-t, digits=1))s"
encode_and_measure("CGE index", p * ".mgz")

# ---------------------------------------------------------------------------
# Summary table
# ---------------------------------------------------------------------------

println()
println("=" ^ 70)
println("  CNR-2000 BPE Comparison: children vs index coding scheme")
println("  ($n_v vertices, $m edges)")
println("=" ^ 70)
println()
@printf("  %-16s %10s %10s %10s\n", "Algorithm", "Children", "Index", "Delta")
println("  " * "-" ^ 48)

for algo in ["STD", "CS", "CGE"]
    c = results["$algo children"]
    i = results["$algo index"]
    delta = round(i.bpe - c.bpe, digits=4)
    sign = delta >= 0 ? "+" : ""
    @printf("  %-16s %10.4f %10.4f %10s\n", algo, c.bpe, i.bpe, "$sign$delta")
end

println()
println("  Files saved in: $DS_DIR")
println()
for (label, r) in sort(collect(results); by = x -> x.first)
    println("    $(basename(r.file))")
end
println()
