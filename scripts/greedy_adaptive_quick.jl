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

# Quick benchmark: just the 3-way adaptive (global + local) on Leiden+LLP + w=64 + Fibonacci
# Baseline already known: 3.3484 BPE (copy-blocks only)

using Pkg
Pkg.activate(normpath(joinpath(@__DIR__, "..")))

using Logging
global_logger(ConsoleLogger(stderr, Logging.Info))

using LightGraphs: nv, ne, outneighbors, vertices, eltype
using Adjacently
using Adjacently.IO: load_adjacency_list_from_csv, BitWriter, BitReader, flush_bitwriter
using Adjacently.Relabeling: relabel_graph, relabel_vertices_llp
using Adjacently.Clustering: leiden_partition
using Adjacently.Graph: subgraph
using Adjacently.Compression: write_greedy_graph_data, read_greedy_graph_data

const CNR_CSV = normpath(joinpath(@__DIR__, "..", "datasets", "webgraph", "cnr-2000", "cnr-2000.csv"))

function relabel_leiden_k1_llp(g)
    T = eltype(g); n = nv(g)
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
    vertex_mapping = Dict{T,T}(); new_id = T(1)
    for C in groups; for v in C; vertex_mapping[v] = new_id; new_id += T(1); end; end
    return vertex_mapping, [length(C) for C in groups]
end

println("Loading CNR-2000...")
g = load_adjacency_list_from_csv(CNR_CSV, ',', true)
n, m = nv(g), ne(g)
println("  $n vertices, $m edges")

println("Leiden+LLP relabeling...")
t_rel = time()
map, gsizes = relabel_leiden_k1_llp(g)
g_rel = relabel_graph(g, map)
println("  Groups: $(gsizes) ($(round(time()-t_rel, digits=1))s)")

T = eltype(g_rel)
nls = Dict{T,Vector{T}}()
for v in vertices(g_rel); nls[T(v)] = sort([T(o) for o in outneighbors(g_rel, v)]); end

# Config: 3-way adaptive, global window w=64
println("\n3-way adaptive w=64 (global window)...")
io = IOBuffer(); bw = BitWriter(io)
t = time()
write_greedy_graph_data(bw, nls, :children, 64;
    integer_encoding=:fibonacci, copy_blocks=true, adaptive_copy=true)
flush_bitwriter(bw; flush_last_bits=true)
bytes = take!(io)
dt = round(time() - t, digits=2)
bpe = round(8.0 * length(bytes) / m, digits=4)
println("  $bpe BPE  ($(length(bytes)) bytes, $(dt)s)")

# Roundtrip
r = BitReader(IOBuffer(bytes))
decoded = read_greedy_graph_data(r, T(length(nls)), :children, T;
    integer_encoding=:fibonacci, copy_blocks=true, adaptive_copy=true, ref_window_size=64)
dec_edges = sum(length(v) for v in values(decoded))
println("  Roundtrip: $(dec_edges == m ? "OK" : "MISMATCH ($dec_edges vs $m)")")

# Config: 3-way adaptive, local window w=64
println("\n3-way adaptive w=64 (local window)...")
io2 = IOBuffer(); bw2 = BitWriter(io2)
t2 = time()
write_greedy_graph_data(bw2, nls, :children, 64;
    integer_encoding=:fibonacci, copy_blocks=true, adaptive_copy=true, cluster_sizes=gsizes)
flush_bitwriter(bw2; flush_last_bits=true)
bytes2 = take!(io2)
dt2 = round(time() - t2, digits=2)
bpe2 = round(8.0 * length(bytes2) / m, digits=4)
println("  $bpe2 BPE  ($(length(bytes2)) bytes, $(dt2)s)")

r2 = BitReader(IOBuffer(bytes2))
decoded2 = read_greedy_graph_data(r2, T(length(nls)), :children, T;
    integer_encoding=:fibonacci, copy_blocks=true, adaptive_copy=true, ref_window_size=64, cluster_sizes=gsizes)
dec_edges2 = sum(length(v) for v in values(decoded2))
println("  Roundtrip: $(dec_edges2 == m ? "OK" : "MISMATCH ($dec_edges2 vs $m)")")

println("\n" * "=" ^ 50)
println("SUMMARY (Leiden+LLP + Fib + w=64)")
println("=" ^ 50)
println("  Baseline (CB only):       3.3484 BPE")
println("  3-way adaptive (global):  $bpe BPE  ($(bpe < 3.3484 ? "-" : "+")$(round(abs(bpe - 3.3484), digits=4)))")
println("  3-way adaptive (local):   $bpe2 BPE  ($(bpe2 < 3.3484 ? "-" : "+")$(round(abs(bpe2 - 3.3484), digits=4)))")
println("=" ^ 50)
