#!/usr/bin/env julia
#
# Diagnostic: Compare BG vs CG K=1 bit budgets on CNR-2000
# Instruments BG to measure per-section bit costs matching CG's breakdown.

include("run_tests_main.jl")

using Logging
global_logger(ConsoleLogger(stderr, Logging.Info))

using LightGraphs: nv, ne, outneighbors, vertices, eltype
using Adjacently
using Adjacently.IO: load_adjacency_list_from_csv, BitWriter, BitReader, write_bit, flush_bitwriter
using Adjacently.Relabeling: relabel_graph, load_llp_ordering
using Adjacently.MGS: write_bg_mgs3_graph, write_cg_mgs3_graph, load_compressed_mgs3_graph
using Adjacently.Compression: write_greedy_graph_data, _write_vertex_header_vlc,
    _vlc_header_cost, write_encoded_value, estimate_encoded_value_cost,
    _write_copy_blocks, _write_adaptive_copy, _write_stop_delta,
    write_intervals_and_residuals, _write_intervals_lr, write_delta,
    delta_encode_vector, write_hybrid_mix_encoded_list,
    ADAPTIVE_ENCODING_OPTIONS, ADAPTIVE_MIL, _greedy_vertex_search,
    _total_bits, GreedyCostBuffer, _greedy_reset!
using Adjacently.Compression.CG: CGParams, CGStats

const DATASET = length(ARGS) >= 1 ? ARGS[1] : "cnr-2000"
const PREFIX  = replace(DATASET, "-" => "")
const DS_DIR  = normpath(joinpath(@__DIR__, "..", "datasets", "webgraph", DATASET))
const DS_CSV  = joinpath(DS_DIR, DATASET * ".csv")

@info "Dataset: $DATASET"
g = load_adjacency_list_from_csv(DS_CSV, ',', true)
@info "  $(nv(g)) vertices, $(ne(g)) edges"

llp_file = joinpath(DS_DIR, "$(PREFIX)_llp_p10.bin")
if isfile(llp_file)
    @info "Applying LLP ordering"
    T_v = eltype(g)
    llp_order = load_llp_ordering(llp_file, T_v)
    g = relabel_graph(g, llp_order)
end

total_edges = ne(g)
n = nv(g)
T_v = eltype(g)

# --- CG K=1 with stats ---
@info "═══ CG K=1 ═══"
cg_params = CGParams(
    intra_ref_enabled=true, intra_ref_window=64,
    intra_ref_rle=false, intra_zigzag=true, intra_stop_deltas=true,
    intra_copy_blocks=true, intra_copy_adaptive=true,
    intra_ref_fixwidth=true, intra_ref_vlc=false,
    intra_add_adaptive=true, intra_raw_adaptive=true,
    intra_adapt_mil=4, intra_intervals=false, intra_mil=4,
    membership=:implicit_ranges,
)
clusters_k1 = [collect(T_v(1):T_v(n))]
out_base = joinpath(DS_DIR, "$(PREFIX)_diag_cg")
write_cg_mgs3_graph(g, out_base, clusters_k1; params=cg_params)
mgz_file = out_base * ".mgz"
cg_bpe = (filesize(mgz_file) - 12) * 8.0 / total_edges
@info "  CG K=1 BPE = $(round(cg_bpe, digits=4))"
rm(mgz_file; force=true)

# --- BG adaptive + lr_split with bit counting ---
@info "═══ BG adaptive + lr_split (instrumented) ═══"

# Build neighbor lists
nls = Dict{T_v, Vector{T_v}}()
for v in vertices(g)
    nls[T_v(v)] = sort([T_v(o) for o in outneighbors(g, v)])
end
vs = n
ie = :fibonacci
ref_window_size = 64
copy_blocks = true
adaptive_copy = true
fixwidth_ref = true
stop_deltas = true
lr_split = true
adaptive_header = true
ref_dist_bits = ceil(Int, log2(max(2, ref_window_size)))

function run_instrumented_bg(nls, vs, T_v, ie, ref_window_size, ref_dist_bits, total_edges; use_lr_split::Bool=true)
    header_bits_total = 0
    ref_dist_bits_total = 0
    copy_bits_total = 0
    add_bits_total = 0
    raw_bits_total = 0
    ref_count = 0
    noref_count = 0
    empty_count = 0

    reference_window = T_v[]
    function add_to_ref_window!(vertex)
        push!(reference_window, vertex)
        if length(reference_window) > ref_window_size
            popfirst!(reference_window)
        end
    end

    bw_out = BitWriter()
    write_bit(bw_out, false)
    for _ in 1:3; write_bit(bw_out, false); end

    for v_idx in 1:vs
        v = T_v(v_idx)
        current_neighbors = sort(get(nls, v, T_v[]))

        if isempty(current_neighbors)
            header_bits_total += 2
            raw_bits_total += 1
            empty_count += 1
            noref_count += 1
            add_to_ref_window!(v)
            _write_vertex_header_vlc(bw_out, :none, :delta, 0; adaptive=true)
            write_bit(bw_out, false)
            continue
        end

        _, actual_ref_mode, mil, ref_result, enc_type, use_stop, res_enc_type, res_mil =
            _greedy_vertex_search(current_neighbors, nls, reference_window, ie;
                vertex_id=v, copy_blocks=true, adaptive_copy=true,
                fixwidth_ref=true, ref_dist_bits=ref_dist_bits,
                stop_deltas=true, lr_split=use_lr_split, adaptive_header=true)

        header_bits_total += 2
        _write_vertex_header_vlc(bw_out, actual_ref_mode, enc_type, mil; adaptive=true)

        target_list = current_neighbors
        write_enc = enc_type
        write_mil = mil

        if actual_ref_mode == :reference && ref_result !== nothing
            ref_distance, copy_bitmap, residuals = ref_result
            ref_count += 1

            b0 = _total_bits(bw_out)
            d = Int(ref_distance) - 1
            for b in (ref_dist_bits-1):-1:0
                write_bit(bw_out, ((d >> b) & 1) == 1)
            end
            ref_dist_bits_total += _total_bits(bw_out) - b0

            b0 = _total_bits(bw_out)
            _write_adaptive_copy(bw_out, copy_bitmap, ie; compact_copy=true)
            copy_bits_total += _total_bits(bw_out) - b0

            target_list = residuals
        else
            noref_count += 1
        end

        b0 = _total_bits(bw_out)
        if write_enc == :interval
            if use_lr_split
                _write_intervals_lr(bw_out, target_list, ie, write_mil; vertex_id=v, tight_intervals=true)
            else
                write_intervals_and_residuals(bw_out, target_list, ie, write_mil; vertex_id=v, tight_intervals=true)
            end
        elseif write_enc == :delta
            _write_stop_delta(bw_out, target_list, ie; vertex_id=v)
        end
        payload_bits = _total_bits(bw_out) - b0
        if actual_ref_mode == :reference
            add_bits_total += payload_bits
        else
            raw_bits_total += payload_bits
        end

        add_to_ref_window!(v)
    end

    flush_bitwriter(bw_out; flush_last_bits=true)
    total_data_bits = bw_out.bit_count

    return (header_bits_total=header_bits_total, ref_dist_bits_total=ref_dist_bits_total,
            copy_bits_total=copy_bits_total, add_bits_total=add_bits_total,
            raw_bits_total=raw_bits_total, ref_count=ref_count, noref_count=noref_count,
            empty_count=empty_count, total_data_bits=total_data_bits)
end

@info "--- BG adaptive + lr_split ---"
r = run_instrumented_bg(nls, vs, T_v, ie, ref_window_size, ref_dist_bits, total_edges)
header_bits_total = r.header_bits_total
ref_dist_bits_total = r.ref_dist_bits_total
copy_bits_total = r.copy_bits_total
add_bits_total = r.add_bits_total
raw_bits_total = r.raw_bits_total
ref_count = r.ref_count
noref_count = r.noref_count
empty_count = r.empty_count
total_data_bits = r.total_data_bits
stream_header_bits = 4

@info "BG bit budget ($total_edges edges):"
@info "  stream header:  $stream_header_bits bits"
@info "  headers:        $header_bits_total bits ($(round(header_bits_total/total_edges, digits=4)) BPE)"
@info "  ref distances:  $ref_dist_bits_total bits ($(round(ref_dist_bits_total/total_edges, digits=4)) BPE)"
@info "  copy:           $copy_bits_total bits ($(round(copy_bits_total/total_edges, digits=4)) BPE)"
@info "  additions:      $add_bits_total bits ($(round(add_bits_total/total_edges, digits=4)) BPE)"
@info "  raw:            $raw_bits_total bits ($(round(raw_bits_total/total_edges, digits=4)) BPE)"
@info "  total:          $total_data_bits bits ($(round(total_data_bits/total_edges, digits=4)) BPE)"
@info "  ref/no_ref:     $ref_count/$noref_count (empty=$empty_count)"

# Compare with CG
cg_headers = 1484769
cg_copy = 1603551
cg_add = 3540386
cg_raw = 1854855
@info ""
@info "Section comparison (BG - CG):"
bg_hdr = header_bits_total + ref_dist_bits_total
@info "  headers+ref:    $bg_hdr vs $cg_headers = $(bg_hdr - cg_headers) bits ($(round((bg_hdr - cg_headers)/total_edges, digits=4)) BPE)"
@info "  copy:           $copy_bits_total vs $cg_copy = $(copy_bits_total - cg_copy) bits ($(round((copy_bits_total - cg_copy)/total_edges, digits=4)) BPE)"
@info "  add:            $add_bits_total vs $cg_add = $(add_bits_total - cg_add) bits ($(round((add_bits_total - cg_add)/total_edges, digits=4)) BPE)"
@info "  raw:            $raw_bits_total vs $cg_raw = $(raw_bits_total - cg_raw) bits ($(round((raw_bits_total - cg_raw)/total_edges, digits=4)) BPE)"

# --- BG adaptive WITHOUT lr_split (for ref count comparison) ---
@info ""
@info "--- BG adaptive WITHOUT lr_split ---"
r2 = run_instrumented_bg(nls, vs, T_v, ie, ref_window_size, ref_dist_bits, total_edges; use_lr_split=false)
@info "BG no-lr_split bit budget:"
@info "  headers:        $(r2.header_bits_total) bits ($(round(r2.header_bits_total/total_edges, digits=4)) BPE)"
@info "  ref distances:  $(r2.ref_dist_bits_total) bits ($(round(r2.ref_dist_bits_total/total_edges, digits=4)) BPE)"
@info "  copy:           $(r2.copy_bits_total) bits ($(round(r2.copy_bits_total/total_edges, digits=4)) BPE)"
@info "  additions:      $(r2.add_bits_total) bits ($(round(r2.add_bits_total/total_edges, digits=4)) BPE)"
@info "  raw:            $(r2.raw_bits_total) bits ($(round(r2.raw_bits_total/total_edges, digits=4)) BPE)"
@info "  total:          $(r2.total_data_bits) bits ($(round(r2.total_data_bits/total_edges, digits=4)) BPE)"
@info "  ref/no_ref:     $(r2.ref_count)/$(r2.noref_count) (empty=$(r2.empty_count))"
@info ""
@info "Ref count comparison:"
@info "  CG K=1 (no lr_split): 193202 refs"
@info "  BG + lr_split:        $ref_count refs"
@info "  BG no lr_split:       $(r2.ref_count) refs"

# --- Degree distribution analysis ---
deg_counts = zeros(Int, 11)  # 0, 1, 2, ..., 9, 10+
for v in 1:n
    d = length(get(nls, T_v(v), T_v[]))
    if d <= 9; deg_counts[d+1] += 1; else; deg_counts[11] += 1; end
end
@info ""
@info "Degree distribution:"
for i in 0:9
    @info "  degree $i: $(deg_counts[i+1]) vertices"
end
@info "  degree 10+: $(deg_counts[11]) vertices"
low_deg = deg_counts[1] + deg_counts[2] + deg_counts[3]  # degree 0, 1, 2
@info "  Total degree 0-2: $low_deg ($(round(100*low_deg/n, digits=1))%)"

# --- BG matching CG settings: no lr_split, no compact_copy, no tight_intervals ---
# This tests if the ref count difference comes from those flags
@info ""
@info "--- BG matching CG settings (no lr, no compact, no tight) ---"
function run_instrumented_bg_cg_style(nls, vs, T_v, ie, ref_window_size, ref_dist_bits, total_edges)
    header_bits_total = 0
    ref_dist_bits_total = 0
    copy_bits_total = 0
    add_bits_total = 0
    raw_bits_total = 0
    ref_count = 0
    noref_count = 0
    empty_count = 0

    reference_window = T_v[]
    function add_to_ref_window!(vertex)
        push!(reference_window, vertex)
        if length(reference_window) > ref_window_size
            popfirst!(reference_window)
        end
    end

    bw_out = BitWriter()
    write_bit(bw_out, false)
    for _ in 1:3; write_bit(bw_out, false); end

    for v_idx in 1:vs
        v = T_v(v_idx)
        current_neighbors = sort(get(nls, v, T_v[]))

        if isempty(current_neighbors)
            header_bits_total += 2
            raw_bits_total += 1
            empty_count += 1
            noref_count += 1
            add_to_ref_window!(v)
            _write_vertex_header_vlc(bw_out, :none, :delta, 0; adaptive=true)
            write_bit(bw_out, false)
            continue
        end

        # Match CG settings: no lr_split, no compact_copy, no tight_intervals
        _, actual_ref_mode, mil, ref_result, enc_type, use_stop, res_enc_type, res_mil =
            _greedy_vertex_search(current_neighbors, nls, reference_window, ie;
                vertex_id=v, copy_blocks=true, adaptive_copy=true,
                fixwidth_ref=true, ref_dist_bits=ref_dist_bits,
                stop_deltas=true, lr_split=false, adaptive_header=true,
                compact_copy=false, tight_intervals=false)

        header_bits_total += 2
        _write_vertex_header_vlc(bw_out, actual_ref_mode, enc_type, mil; adaptive=true)

        target_list = current_neighbors
        if actual_ref_mode == :reference && ref_result !== nothing
            ref_distance, copy_bitmap, residuals = ref_result
            ref_count += 1
            b0 = _total_bits(bw_out)
            d = Int(ref_distance) - 1
            for b in (ref_dist_bits-1):-1:0
                write_bit(bw_out, ((d >> b) & 1) == 1)
            end
            ref_dist_bits_total += _total_bits(bw_out) - b0
            b0 = _total_bits(bw_out)
            _write_adaptive_copy(bw_out, copy_bitmap, ie)  # non-compact
            copy_bits_total += _total_bits(bw_out) - b0
            target_list = residuals
        else
            noref_count += 1
        end

        b0 = _total_bits(bw_out)
        if enc_type == :interval
            write_intervals_and_residuals(bw_out, target_list, ie, mil; vertex_id=v)
        elseif enc_type == :delta
            _write_stop_delta(bw_out, target_list, ie; vertex_id=v)
        end
        payload_bits = _total_bits(bw_out) - b0
        if actual_ref_mode == :reference
            add_bits_total += payload_bits
        else
            raw_bits_total += payload_bits
        end
        add_to_ref_window!(v)
    end
    flush_bitwriter(bw_out; flush_last_bits=true)
    total_data_bits = bw_out.bit_count
    return (ref_count=ref_count, noref_count=noref_count, empty_count=empty_count,
            total_data_bits=total_data_bits, header_bits_total=header_bits_total,
            ref_dist_bits_total=ref_dist_bits_total, copy_bits_total=copy_bits_total,
            add_bits_total=add_bits_total, raw_bits_total=raw_bits_total)
end
r3 = run_instrumented_bg_cg_style(nls, vs, T_v, ie, ref_window_size, ref_dist_bits, total_edges)
@info "BG CG-style bit budget:"
@info "  total: $(r3.total_data_bits) bits ($(round(r3.total_data_bits/total_edges, digits=4)) BPE)"
@info "  ref/no_ref: $(r3.ref_count)/$(r3.noref_count) (empty=$(r3.empty_count))"
@info "  headers+ref: $(r3.header_bits_total + r3.ref_dist_bits_total) bits"
@info "  copy: $(r3.copy_bits_total) bits"
@info "  add:  $(r3.add_bits_total) bits"
@info "  raw:  $(r3.raw_bits_total) bits"
