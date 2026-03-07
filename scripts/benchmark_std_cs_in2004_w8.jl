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
# Benchmark STD and CS on IN-2004 with ref_window_size=8
#

using Pkg
Pkg.activate(normpath(joinpath(@__DIR__, "..")))

using Logging
global_logger(ConsoleLogger(stderr, Logging.Info))

using LightGraphs: nv, ne, outneighbors, vertices, eltype
using Adjacently
using Adjacently.IO: load_adjacency_list_from_csv
using Adjacently.Relabeling: relabel_graph, relabel_vertices_llp
using Adjacently.Clustering: leiden_partition
using Adjacently.Graph: subgraph
using Adjacently.MGS: write_std_mgs3_graph, write_cs_mgs3_graph, load_compressed_mgs3_graph

const DS = "in-2004"
const DS_DIR = normpath(joinpath(@__DIR__, "..", "datasets", "webgraph", DS))
const DS_CSV = joinpath(DS_DIR, "$DS.csv")

isfile(DS_CSV) || error("CSV not found: $DS_CSV")

println("=" ^ 60)
println("STD & CS benchmark — $DS — window=8")
println("=" ^ 60)

# Load graph
@info "Loading $DS..."
t0 = time()
g = load_adjacency_list_from_csv(DS_CSV, ',', true)
n = nv(g); m = ne(g)
@info "  Loaded: $n vertices, $m edges ($(round(time()-t0, digits=1))s)"

# Leiden K=1 + LLP relabeling (same as best-compression test)
@info "Leiden + LLP relabeling..."
t_rel = time()
T = eltype(g)
part = leiden_partition(g; max_passes=8, max_levels=5)
counts = Dict{Int,Int}()
for c in part; counts[c] = get(counts, c, 0) + 1; end
top_label = argmax(counts)
groups = [T[], T[]]
for i in 1:n; push!(groups[part[i] == top_label ? 1 : 2], T(i)); end
for gi in 1:2
    C = groups[gi]; length(C) <= 2 && continue
    sg, oni, _ = subgraph(g, C)
    mapping = relabel_vertices_llp(sg, :sym; passes=5)
    sort!(C, by = v -> Int(mapping[oni[v]]))
    groups[gi] = C
end
vertex_mapping = let new_id = T(1)
    d = Dict{T,T}()
    for C in groups; for v in C; d[v] = new_id; new_id += T(1); end; end
    d
end
g_rel = relabel_graph(g, vertex_mapping)
m_rel = ne(g_rel)
@info "  Relabeled in $(round(time()-t_rel, digits=1))s"

# Helper: extract neighbor lists
function extract_nls(g)
    T = eltype(g)
    nls = Dict{T,Vector{T}}()
    for v in vertices(g)
        nls[T(v)] = sort([T(o) for o in outneighbors(g, v)])
    end
    return nls
end

nls = extract_nls(g_rel)

# ── STD window=8 ─────────────────────────────────────────────────────────────

println("\n" * "─" ^ 60)
println("STD — window=8")
println("─" ^ 60)

std_base = joinpath(DS_DIR, "in2004_std_w8")
std_mgz = std_base * ".mgz"

t_enc = time()
write_std_mgs3_graph(g_rel, std_base;
    coding_scheme=:children, integer_encoding=:fibonacci,
    ref_window_size=8, copy_blocks=true, adaptive_copy=true,
    stop_deltas=true, empty_prefix=true, compact_copy=true,
    tight_intervals=true, vlc2=true)
dt_enc = round(time() - t_enc, digits=2)

bpe_std = round(8.0 * filesize(std_mgz) / m_rel, digits=4)
@info "  STD w=8 encoded: $bpe_std BPE ($(filesize(std_mgz)) bytes, $(dt_enc)s)"

# Verify roundtrip
t_dec = time()
g_std = load_compressed_mgs3_graph(std_mgz;
    copy_blocks=true, adaptive_copy=true, ref_window_size=8,
    stop_deltas=true, empty_prefix=true, compact_copy=true,
    tight_intervals=true, vlc2=true)
dt_dec = round(time() - t_dec, digits=2)

loaded_nls = extract_nls(g_std)
dec_ok = nv(g_std) == nv(g_rel) && ne(g_std) == ne(g_rel) &&
    all(sort(get(loaded_nls, v, T[])) == sort(orig) for (v, orig) in nls)
@info "  STD roundtrip: $(dec_ok ? "PASSED" : "FAILED") (decode $(dt_dec)s)"

# ── CS window=8 ──────────────────────────────────────────────────────────────

println("\n" * "─" ^ 60)
println("CS — window=8")
println("─" ^ 60)

cs_base = joinpath(DS_DIR, "in2004_cs_w8")
cs_mgz = cs_base * ".mgz"

t_enc = time()
write_cs_mgs3_graph(g_rel, cs_base;
    coding_scheme=:children, integer_encoding=:fibonacci,
    ref_window_size=8, compact_copy=true, tight_intervals=true)
dt_enc = round(time() - t_enc, digits=2)

bpe_cs = round(8.0 * filesize(cs_mgz) / m_rel, digits=4)
@info "  CS w=8 encoded: $bpe_cs BPE ($(filesize(cs_mgz)) bytes, $(dt_enc)s)"

# Verify roundtrip
t_dec = time()
g_cs = load_compressed_mgs3_graph(cs_mgz;
    compact_copy=true, tight_intervals=true, ref_window_size=8)
dt_dec = round(time() - t_dec, digits=2)

loaded_nls = extract_nls(g_cs)
dec_ok = nv(g_cs) == nv(g_rel) && ne(g_cs) == ne(g_rel) &&
    all(sort(get(loaded_nls, v, T[])) == sort(orig) for (v, orig) in nls)
@info "  CS roundtrip: $(dec_ok ? "PASSED" : "FAILED") (decode $(dt_dec)s)"

# ── Summary ──────────────────────────────────────────────────────────────────

println("\n" * "=" ^ 60)
println("SUMMARY — $DS window=8")
println("=" ^ 60)
println("  STD w=8:   $bpe_std BPE")
println("  CS  w=8:   $bpe_cs BPE")
println("  CGE K=1:  1.7513 BPE  (best known)")
println("  WebGraph:  1.767 BPE   (reference)")
println("=" ^ 60)
