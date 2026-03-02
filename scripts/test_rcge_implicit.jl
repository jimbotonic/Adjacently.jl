#!/usr/bin/env julia

#
# Fair comparison: RCGE FW64 K=1 with EF vs implicit_ranges membership
#
# Key insight: encode_level() re-sorts clusters by vertex ID internally (line 565).
# Local indices are always based on vertex ID rank within the cluster, NOT LLP rank.
#
# To match, the relabeling must also use vertex ID rank:
#   new_id(v) = rank of v in ID-sorted cluster
# This makes the relabeled graph's local index space identical to the EF case.
# The only difference becomes the membership encoding (88K bytes vs 7 bytes).
#

using Pkg
Pkg.activate(normpath(joinpath(@__DIR__, "..")))

using Logging
global_logger(ConsoleLogger(stderr, Logging.Info))

using LightGraphs
using LightGraphs: nv, outneighbors
using Adjacently
using Adjacently.RCGE: encode_level, RCGEParams, RCGEStats, decode_level
using Adjacently.IO: load_adjacency_list_from_csv, BitWriter, BitReader, flush_bitwriter
using Adjacently.Clustering: leiden_partition
using Adjacently.Relabeling: relabel_graph

count_edges(h) = sum(length(outneighbors(h, v)) for v in 1:nv(h))

function make_params(; kwargs...)
    d = Dict{Symbol,Any}(
        :L=>128, :varint=>:fibonacci, :count_varint=>:fibonacci, :gap=>:fibonacci,
        :degree=>:elias_delta, :undirected_pairs=>false, :perm_strategy=>:blockpos,
        :membership=>:elias_fano, :inter_strategy=>:perm, :intra_ref_enabled=>true,
        :intra_ref_window=>32, :intra_ref_rle=>false,
        :intra_block_try=>false, :positions_mode=>:delta, :additions_mode=>:delta,
        :min_cluster_density=>0.0, :intra_intervals=>false, :intra_mil=>4,
        :intra_greedy_mil=>false, :intra_zigzag=>false, :intra_stop_deltas=>false,
        :intra_copy_blocks=>false, :intra_ref_fixwidth=>false,
    )
    for (k, v) in kwargs; d[k] = v; end
    return RCGEParams(; d...)
end

fw64_base(; kwargs...) = make_params(
    intra_zigzag=true, intra_stop_deltas=true, intra_copy_blocks=true,
    intra_ref_window=64, intra_ref_fixwidth=true; kwargs...)

println("=" ^ 70)
println("RCGE FW64 K=1: EF vs Implicit Membership (fair comparison)")
println("=" ^ 70)

cnr_csv = normpath(joinpath(@__DIR__, "..", "datasets", "webgraph", "cnr-2000", "cnr-2000.csv"))
isfile(cnr_csv) || error("CNR-2000 CSV not found at $cnr_csv")

println("\nLoading CNR-2000...")
t0 = time()
g  = load_adjacency_list_from_csv(cnr_csv, ',', true)
n  = nv(g)
m  = count_edges(g)
TV = eltype(g)
println("  Loaded: $n vertices, $m edges, eltype=$TV ($(round(time()-t0, digits=2))s)")

# ── Step 1: Leiden K=1 partition ──────────────────────────────────────────────
println("\nStep 1: Leiden K=1 partition...")
t1 = time()
part = leiden_partition(g; max_passes=8, max_levels=5)
nclusters = maximum(part)
println("  Partitioned into $nclusters clusters ($(round(time()-t1, digits=3))s)")

counts = Dict{Int,Int}()
for c in part; counts[c] = get(counts, c, 0) + 1; end
top_label = argmax(counts)

clusters = [TV[], TV[]]
for i in 1:n
    push!(clusters[part[i] == top_label ? 1 : 2], TV(i))
end
println("  Clusters: $(length(clusters[1])) (top) + $(length(clusters[2])) (rest)")

# ── Step 2: Sort clusters by vertex ID (what encode_level does internally) ────
# Note: encode_level re-sorts by vertex ID at line 565. Local indices are always
# based on vertex-ID rank. We do the same sort here so local indices match.
clusters_by_id = [sort(copy(C)) for C in clusters]
println("\nStep 2: Clusters sorted by vertex ID (matches encode_level internal sort)")
println("  Cluster 1 ID range: [$(clusters_by_id[1][1])..$(clusters_by_id[1][end])]")
println("  Cluster 2 ID range: [$(clusters_by_id[2][1])..$(clusters_by_id[2][end])]")

# ── Step 3: Build relabeled graph — relabel by vertex-ID rank ─────────────────
# new_id(v) = rank of v in its ID-sorted cluster (offset by previous cluster size)
# This matches the local index space used by encode_level with EF.
println("\nStep 3: Build relabeled graph (relabel by vertex-ID rank)...")
t3 = time()
vertex_map = let new_id = TV(1)
    d = Dict{TV,TV}()
    for C in clusters_by_id      # each C is sorted by original vertex ID
        for v in C
            d[v] = new_id
            new_id += TV(1)
        end
    end
    d
end
g_rel = relabel_graph(g, vertex_map)
S1 = length(clusters_by_id[1])
S2 = length(clusters_by_id[2])
clusters_contiguous = [TV.(1:S1), TV.(S1+1:n)]
println("  Relabeled: cluster 1 → 1..$S1 ($S1), cluster 2 → $(S1+1)..$n ($S2) ($(round(time()-t3, digits=2))s)")

# ── Encoding A: EF membership (original graph) ────────────────────────────────
# encode_level receives LLP-sorted clusters but re-sorts internally → ID-sorted local indices
println("\nEncoding A: elias_fano membership (original graph, ID-sorted clusters)...")
fw64_ef = fw64_base(membership=:elias_fano)
io_ef = IOBuffer(); w_ef = BitWriter(io_ef)
t_a = time()
stats_ef = RCGEStats()
# Pass ID-sorted clusters (LLP would be overridden anyway)
encode_level(w_ef, g, clusters_by_id; params=fw64_ef, stats=stats_ef)
flush_bitwriter(w_ef; flush_last_bits=true)
bytes_ef = take!(io_ef)
bpe_ef   = 8.0 * length(bytes_ef) / m
println("  Encoded: $(length(bytes_ef)) bytes, BPE=$(round(bpe_ef, digits=4)) ($(round(time()-t_a, digits=2))s)")
println("  Stats: membership=$(ceil(Int, stats_ef.bits_membership/8))B, intra=$(ceil(Int, stats_ef.bits_intra/8))B [headers=$(ceil(Int, stats_ef.bits_intra_headers/8))B, copy=$(ceil(Int, stats_ef.bits_intra_copy/8))B, add=$(ceil(Int, stats_ef.bits_intra_add/8))B, raw=$(ceil(Int, stats_ef.bits_intra_raw/8))B], inter=$(ceil(Int, (stats_ef.bits_inter_headers+stats_ef.bits_inter_lists)/8))B")

# ── Encoding B: Implicit ranges (relabeled graph, contiguous clusters) ─────────
# The relabeled graph has:
#   new_id(v) = vertex-ID rank in its cluster → same local index as EF
#   clusters_contiguous = [1..S1, S1+1..N] → encode_level's internal sort leaves them unchanged
# → intra encoding should be bit-for-bit identical to EF
println("\nEncoding B: implicit_ranges membership (relabeled graph, contiguous clusters)...")
fw64_impl = fw64_base(membership=:implicit_ranges)
io_impl = IOBuffer(); w_impl = BitWriter(io_impl)
t_b = time()
stats_impl = RCGEStats()
encode_level(w_impl, g_rel, clusters_contiguous; params=fw64_impl, stats=stats_impl)
flush_bitwriter(w_impl; flush_last_bits=true)
bytes_impl = take!(io_impl)
bpe_impl   = 8.0 * length(bytes_impl) / m
println("  Encoded: $(length(bytes_impl)) bytes, BPE=$(round(bpe_impl, digits=4)) ($(round(time()-t_b, digits=2))s)")
println("  Stats: membership=$(ceil(Int, stats_impl.bits_membership/8))B, intra=$(ceil(Int, stats_impl.bits_intra/8))B [headers=$(ceil(Int, stats_impl.bits_intra_headers/8))B, copy=$(ceil(Int, stats_impl.bits_intra_copy/8))B, add=$(ceil(Int, stats_impl.bits_intra_add/8))B, raw=$(ceil(Int, stats_impl.bits_intra_raw/8))B], inter=$(ceil(Int, (stats_impl.bits_inter_headers+stats_impl.bits_inter_lists)/8))B")

# Assert intra bits match
if stats_ef.bits_intra == stats_impl.bits_intra
    println("  ✓ Intra encoding is bit-for-bit identical to EF (as expected)")
else
    diff = stats_impl.bits_intra - stats_ef.bits_intra
    println("  ✗ Intra differs by $(diff) bits ($(round(diff/8, digits=0)) bytes)")
end

# ── Verify roundtrip for implicit ─────────────────────────────────────────────
println("\nVerifying implicit roundtrip...")
t_dec = time()
r_impl = BitReader(IOBuffer(bytes_impl))
decoded = decode_level(r_impl, fw64_impl; T=TV, directed=true)
dec_edges = sum(length(v) for v in values(decoded))
dt = round(time() - t_dec, digits=3)
if dec_edges == m
    println("  Roundtrip: OK ($dec_edges edges, $(dt)s)")
    local bad = 0
    for v in TV.(1:min(n, 200))
        orig = sort(TV.(outneighbors(g_rel, v)))
        dec  = sort(get(decoded, v, TV[]))
        if orig != dec; bad += 1; end
    end
    println("  Spot-check (first 200 vertices): $(bad == 0 ? "all match" : "$bad mismatches")")
else
    println("  Roundtrip: MISMATCH (expected $m, got $dec_edges)")
end

# Save implicit file
out_path = normpath(joinpath(@__DIR__, "..", "datasets", "webgraph", "cnr-2000", "cnr2000_rcge_gllp_implicit.rcge"))
open(out_path, "w") do f; write(f, bytes_impl); end
println("  Saved: $out_path")

# ── Summary ───────────────────────────────────────────────────────────────────
println("\n" * "=" ^ 70)
println("SUMMARY (isolates membership overhead — same intra encoding)")
println("=" ^ 70)
ef_bytes   = length(bytes_ef)
impl_bytes = length(bytes_impl)
savings_b   = ef_bytes - impl_bytes
savings_bpe = bpe_ef - bpe_impl
println("  EF membership:       $(round(bpe_ef, digits=4)) BPE  ($(ef_bytes) bytes)")
println("  Implicit membership: $(round(bpe_impl, digits=4)) BPE  ($(impl_bytes) bytes)")
println("  Savings:             $(savings_b) bytes, $(round(savings_bpe, digits=4)) BPE")
println("")
println("  Membership only: EF=$(ceil(Int, stats_ef.bits_membership/8))B vs Implicit=$(ceil(Int, stats_impl.bits_membership/8))B → saves $(ceil(Int, stats_ef.bits_membership/8) - ceil(Int, stats_impl.bits_membership/8)) bytes")
println("  Intra:           EF=$(ceil(Int, stats_ef.bits_intra/8))B vs Implicit=$(ceil(Int, stats_impl.bits_intra/8))B")
println("")
println("  Reference: RCGE FW64 K=1 (saved file):  2.887 BPE  (1160561 bytes)")
println("=" ^ 70)
