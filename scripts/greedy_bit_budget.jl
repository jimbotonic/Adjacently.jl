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

# Bit budget analysis for the greedy encoder on CNR-2000
# Measures exactly where bits go for the best configuration

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
using Adjacently.Compression: write_greedy_graph_data, read_greedy_graph_data,
    _greedy_vertex_search, _vlc_header_cost, ENCODING_OPTIONS, _estimate_base_cost,
    _estimate_adaptive_copy_cost, _estimate_copy_blocks_cost, _try_find_reference,
    _write_vertex_header_vlc, _write_adaptive_copy, _write_copy_blocks,
    _write_stop_delta, write_intervals_and_residuals, write_delta,
    write_encoded_value, write_bitmap_adaptive, _estimate_stop_delta_cost,
    estimate_interval_runlength_encoding_cost, MIN_INTERVAL_LENGTH

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
    return vertex_mapping
end

println("=" ^ 70)
println("Greedy Encoder Bit Budget Analysis — CNR-2000")
println("=" ^ 70)

isfile(CNR_CSV) || error("CNR-2000 CSV not found at $CNR_CSV")

println("\nLoading CNR-2000...")
t0 = time()
g = load_adjacency_list_from_csv(CNR_CSV, ',', true)
n, m = nv(g), ne(g)
println("  $n vertices, $m edges ($(round(time()-t0, digits=1))s)")

println("Leiden+LLP relabeling...")
t_rel = time()
lmap = relabel_leiden_k1_llp(g)
g_rel = relabel_graph(g, lmap)
println("  Done ($(round(time()-t_rel, digits=1))s)")

T = eltype(g_rel)
nls = Dict{T,Vector{T}}()
for v in vertices(g_rel); nls[T(v)] = sort([T(o) for o in outneighbors(g_rel, v)]); end

# Run the best config first to get total
println("\nEncoding with best config (adaptive + stop_d + empty)...")
io = IOBuffer(); bw = BitWriter(io)
t = time()
stats = Dict{Any,Int}()
write_greedy_graph_data(bw, nls, :children, 64;
    integer_encoding=:fibonacci, copy_blocks=true,
    adaptive_copy=true, stop_deltas=true, empty_prefix=true,
    stats=stats)
flush_bitwriter(bw; flush_last_bits=true)
bytes = take!(io)
dt = round(time() - t, digits=1)
total_bits = 8 * length(bytes)
bpe = round(total_bits / m, digits=4)
println("  Total: $bpe BPE  ($(length(bytes)) B, $(dt)s)")

println("\n" * "─" ^ 70)
println("Action distribution:")
println("─" ^ 70)
total_v = sum(values(stats))
sorted_stats = sort(collect(stats), by=x->-x[2])
for (key, count) in sorted_stats
    ref_mode, enc_type, mil = key
    pct = round(100.0 * count / total_v, digits=1)
    hdr_bits = _vlc_header_cost(ref_mode, enc_type, mil)
    println("  $(rpad(string(key), 40)) $(rpad(count, 8)) ($(pct)%)  hdr=$(hdr_bits) bits")
end

println("\n" * "─" ^ 70)
println("Estimated bit budget by component:")
println("─" ^ 70)

# Now do a manual per-vertex analysis
ie = :fibonacci
ref_window_size = 64
reference_window = T[]

function analyze_bit_budget(nls, m, ie, ref_window_size, T)

bits_empty_prefix = 0
bits_vlc_header = 0
bits_ref_distance = 0
bits_copy_positions = 0
bits_body = 0
count_ref = 0
count_noref = 0
count_empty = 0
count_interval = 0
count_delta = 0
count_rle = 0

for v_idx in 1:Int(length(nls))
    v = T(v_idx)
    current_neighbors = sort(get(nls, v, T[]))

    # Empty prefix
    bits_empty_prefix += 1
    if isempty(current_neighbors)
        count_empty += 1
        push!(reference_window, v)
        if length(reference_window) > ref_window_size
            popfirst!(reference_window)
        end
        continue
    end

    # Greedy search
    _, ref_mode, mil, ref_result, enc_type, use_stop, res_enc_type, res_mil = _greedy_vertex_search(
        current_neighbors, nls, reference_window, ie;
        vertex_id=v, copy_blocks=true, adaptive_copy=true,
        stop_deltas=true)

    # VLC header
    bits_vlc_header += _vlc_header_cost(ref_mode, enc_type, mil)

    if enc_type == :interval; count_interval += 1
    elseif enc_type == :delta; count_delta += 1
    elseif enc_type == :rle; count_rle += 1
    end

    target_list = current_neighbors
    if ref_mode != :none && ref_result !== nothing
        count_ref += 1
        ref_distance, copy_bitmap, residuals = ref_result
        # Reference distance
        bits_ref_distance += Adjacently.Compression.estimate_encoded_value_cost(ref_distance, ie)
        # Copy positions (adaptive)
        bits_copy_positions += _estimate_adaptive_copy_cost(copy_bitmap, ie)
        target_list = residuals
    else
        count_noref += 1
    end

    # Body encoding
    body_cost = _estimate_base_cost(target_list, ie, mil, enc_type; vertex_id=v, stop_deltas=use_stop)
    bits_body += body_cost

    push!(reference_window, v)
    if length(reference_window) > ref_window_size
        popfirst!(reference_window)
    end
end

return (bits_empty_prefix, bits_vlc_header, bits_ref_distance, bits_copy_positions, bits_body,
        count_ref, count_noref, count_empty, count_interval, count_delta, count_rle)
end  # function analyze_bit_budget

bits_empty_prefix, bits_vlc_header, bits_ref_distance, bits_copy_positions, bits_body,
    count_ref, count_noref, count_empty, count_interval, count_delta, count_rle =
    analyze_bit_budget(nls, m, :fibonacci, 64, T)

total_estimated = bits_empty_prefix + bits_vlc_header + bits_ref_distance + bits_copy_positions + bits_body
# Stream header: 4 bits (coding_scheme flag + encoding tag)
bits_stream_header = 4
total_estimated += bits_stream_header

println("  Stream header:      $(rpad(bits_stream_header, 12)) bits  $(round(bits_stream_header / m, digits=4)) BPE")
println("  Empty-prefix flags: $(rpad(bits_empty_prefix, 12)) bits  $(round(bits_empty_prefix / m, digits=4)) BPE")
println("  VLC headers:        $(rpad(bits_vlc_header, 12)) bits  $(round(bits_vlc_header / m, digits=4)) BPE")
println("  Ref distances:      $(rpad(bits_ref_distance, 12)) bits  $(round(bits_ref_distance / m, digits=4)) BPE")
println("  Copy positions:     $(rpad(bits_copy_positions, 12)) bits  $(round(bits_copy_positions / m, digits=4)) BPE")
println("  Body (enc data):    $(rpad(bits_body, 12)) bits  $(round(bits_body / m, digits=4)) BPE")
println("  ─────────────────────────────────────────────────────")
println("  Estimated total:    $(rpad(total_estimated, 12)) bits  $(round(total_estimated / m, digits=4)) BPE")
println("  Actual total:       $(rpad(total_bits, 12)) bits  $(round(total_bits / m, digits=4)) BPE")
println("  Estimation error:   $(total_bits - total_estimated) bits ($(round((total_bits - total_estimated) / m, digits=4)) BPE)")

println("\n" * "─" ^ 70)
println("Vertex counts:")
println("─" ^ 70)
println("  Empty (degree 0):   $count_empty")
println("  Referenced:         $count_ref  ($(round(100.0*count_ref/(count_ref+count_noref-count_empty), digits=1))%)")
println("  Non-referenced:     $(count_noref - count_empty)")
println("  Interval enc:       $count_interval")
println("  Delta enc:          $count_delta")
println("  RLE enc:            $count_rle")

println("\n" * "─" ^ 70)
println("WebGraph BV comparison:")
println("─" ^ 70)
println("  Component          BV BPE    Greedy BPE   Difference")
println("  Outdegrees         0.516     (implicit)   —")
println("  References         0.243     $(round(bits_ref_distance / m, digits=3))        $(round(bits_ref_distance/m - 0.243, digits=3))")
println("  Copy blocks        0.421     $(round(bits_copy_positions / m, digits=3))        $(round(bits_copy_positions/m - 0.421, digits=3))")
println("  Intervals          0.258     (in body)    —")
println("  Residuals          1.460     (in body)    —")
println("  Structure/headers  —         $(round((bits_vlc_header + bits_empty_prefix) / m, digits=3))        —")
println("  TOTAL              2.897     $(round(total_bits / m, digits=3))        $(round(total_bits/m - 2.897, digits=3))")
println("=" ^ 70)
