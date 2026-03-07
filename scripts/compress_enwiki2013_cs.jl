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
# Compress enwiki-2013 with CS (Command Stream) — Leiden K=1 + LLP reordering
#

using Pkg
Pkg.activate(normpath(joinpath(@__DIR__, "..")))

using Logging
global_logger(ConsoleLogger(stderr, Logging.Info))

using LightGraphs: nv, ne, eltype, outneighbors, vertices
using Adjacently
using Adjacently.IO: load_adjacency_list_from_csv
using Adjacently.Relabeling: relabel_graph, relabel_vertices_llp
using Adjacently.Clustering: leiden_partition
using Adjacently.Graph: subgraph
using Adjacently.MGS: write_cs_mgs3_graph, load_compressed_mgs3_graph

const DS = "enwiki-2013"
const DS_DIR = normpath(joinpath(@__DIR__, "..", "datasets", "webgraph", DS))
const DS_CSV = joinpath(DS_DIR, "$DS.csv")
const PREFIX = replace(DS, "-" => "")

isfile(DS_CSV) || error("CSV not found: $DS_CSV")

println("=" ^ 70)
println("CS compression — $DS (Leiden+LLP reordering)")
println("=" ^ 70)

# ── Load graph ───────────────────────────────────────────────────────────────

@info "Loading $DS..."
t0 = time()
g = load_adjacency_list_from_csv(DS_CSV, ',', true)
n = nv(g); m = ne(g)
T = eltype(g)
@info "  Loaded: $n vertices, $m edges ($(round(time()-t0, digits=1))s)"

# ── Leiden K=1 + LLP reordering ──────────────────────────────────────────────

@info "Running Leiden partition..."
t1 = time()
part = leiden_partition(g; max_passes=8, max_levels=5)
counts = Dict{Int,Int}()
for c in part; counts[c] = get(counts, c, 0) + 1; end
top_label = argmax(counts)
groups = [T[], T[]]
for i in 1:n; push!(groups[part[i] == top_label ? 1 : 2], T(i)); end
@info "  Leiden: group sizes $(length.(groups)) ($(round(time()-t1, digits=1))s)"

@info "Running LLP within each group..."
t2 = time()
for gi in 1:2
    C = groups[gi]; length(C) <= 2 && continue
    sg, oni, _ = subgraph(g, C)
    @info "  LLP group $gi ($(nv(sg)) vertices)..."
    mapping = relabel_vertices_llp(sg, :sym; passes=5)
    sort!(C, by = v -> Int(mapping[oni[v]]))
    groups[gi] = C
end
@info "  LLP done ($(round(time()-t2, digits=1))s)"

@info "Building relabeled graph..."
t3 = time()
vertex_mapping = let new_id = T(1)
    d = Dict{T,T}()
    for C in groups
        for v in C
            d[v] = new_id
            new_id += T(1)
        end
    end
    d
end
g_rel = relabel_graph(g, vertex_mapping)
@info "  Relabeled in $(round(time()-t3, digits=1))s"

# ── CS encode ────────────────────────────────────────────────────────────────

cs_base = joinpath(DS_DIR, "$(PREFIX)_cs")
cs_mgz = cs_base * ".mgz"

@info "Encoding CS → $cs_mgz"
t_enc = time()
write_cs_mgs3_graph(g_rel, cs_base;
    coding_scheme=:children, integer_encoding=:fibonacci,
    ref_window_size=8, compact_copy=true, tight_intervals=true)
dt_enc = round(time() - t_enc, digits=1)

bpe = round(8.0 * filesize(cs_mgz) / m, digits=4)
@info "  Encoded: $bpe BPE ($(filesize(cs_mgz)) bytes, $(dt_enc)s)"

# ── Verify roundtrip ─────────────────────────────────────────────────────────

@info "Verifying roundtrip..."
t_dec = time()
g_loaded = load_compressed_mgs3_graph(cs_mgz;
    compact_copy=true, tight_intervals=true, ref_window_size=8)
dt_dec = round(time() - t_dec, digits=1)

if nv(g_loaded) == nv(g_rel) && ne(g_loaded) == ne(g_rel)
    @info "  Roundtrip: PASSED ($(nv(g_loaded)) vertices, $(ne(g_loaded)) edges, decode $(dt_dec)s)"
else
    @warn "  Roundtrip: MISMATCH! Expected $(nv(g_rel))v/$(ne(g_rel))e, got $(nv(g_loaded))v/$(ne(g_loaded))e"
end

# ── Summary ──────────────────────────────────────────────────────────────────

println("\n" * "=" ^ 70)
println("SUMMARY — $DS CS (Leiden+LLP)")
println("=" ^ 70)
println("  CS w=8 LLP:  $bpe BPE  ($(filesize(cs_mgz)) bytes)")
println("  CS w=8 nat:  16.6219 BPE (no reordering)")
println("  CGE K=1:    16.5337 BPE (no reordering)")
println("  WebGraph BV: 12.639 BPE (reference)")
println("  Encode time: $(dt_enc)s")
println("  Decode time: $(dt_dec)s")
println("=" ^ 70)
