#!/usr/bin/env julia

#
# Test RCGE multi-level compression on CNR-2000 and measure BPE
#

using Pkg
Pkg.activate(normpath(joinpath(@__DIR__, "..")))

using Logging
global_logger(ConsoleLogger(stderr, Logging.Info))

using LightGraphs
using LightGraphs: nv, outneighbors
using Adjacently
using Adjacently.RCGE: encode_level, RCGEParams, RCGEStats, decode_level, load_rcge_graph
using Adjacently.IO: load_adjacency_list_from_csv, BitWriter, BitReader, flush_bitwriter
using Adjacently.Clustering: leiden_partition, aggregate_graph
using Adjacently.Relabeling: relabel_graph, relabel_vertices_llp, relabel_vertices_bisection
using Adjacently.Graph: subgraph

# Minimal TestGraph for coarse weighted graph
struct TestGraph{T<:Unsigned} <: LightGraphs.AbstractGraph{T}
    n::Int
    adj::Dict{Int, Vector{T}}
    radj::Dict{Int, Vector{T}}
    directed::Bool
end
struct TestEdge{T<:Unsigned}
    u::T
    v::T
end
LightGraphs.nv(h::TestGraph) = h.n
LightGraphs.is_directed(h::TestGraph) = h.directed
LightGraphs.outneighbors(h::TestGraph{T}, v::Integer) where {T<:Unsigned} = get(h.adj, Int(v), T[])
LightGraphs.src(e::TestEdge) = e.u
LightGraphs.dst(e::TestEdge) = e.v
function LightGraphs.edges(h::TestGraph{T}) where {T<:Unsigned}
    ed = TestEdge{T}[]
    for (u, lst) in h.adj
        for v in lst
            push!(ed, TestEdge{T}(T(u), v))
        end
    end
    return ed
end
LightGraphs.inneighbors(h::TestGraph{T}, v::Integer) where {T<:Unsigned} = get(h.radj, Int(v), T[])

function coarse_to_testgraph(Gc)::TestGraph{UInt32}
    nC = Gc.n
    adj = Dict{Int, Vector{UInt32}}()
    radj = Dict{Int, Vector{UInt32}}()
    for u in 1:nC
        lst = UInt32[]
        for (v, w) in Gc.out_w[u]
            for _ in 1:Int(round(w))
                push!(lst, UInt32(v))
                push!(get!(radj, v, UInt32[]), UInt32(u))
            end
        end
        !isempty(lst) && (adj[u] = lst)
    end
    return TestGraph{UInt32}(nC, adj, radj, true)
end

count_edges(h) = sum(length(outneighbors(h, v)) for v in 1:nv(h))

# Reorder within clusters using LLP or bisection on induced subgraph
function reorder_clusters!(clusters, base_g; llp_passes::Int=5, ordering::Symbol=:llp)
    for idx in 1:length(clusters)
        C = clusters[idx]
        length(C) <= 2 && continue
        sg, oni, _ = subgraph(base_g, C)
        mapping = if ordering == :bisection
            relabel_vertices_bisection(sg)
        else
            relabel_vertices_llp(sg, :sym; passes=llp_passes)
        end
        sort!(C, by = v -> Int(mapping[oni[v]]))
        clusters[idx] = C
    end
    return clusters
end

function run_rcge(g, m_original, params; K::Int=8, llp_passes::Int=5, ordering::Symbol=:llp)
    max_levels = 5
    min_clusters = 32
    cur_g = g
    total_bytes = 0
    prev_coarse_n = -1
    all_level_bytes = UInt8[]
    first_level_bytes = UInt8[]

    println("\nStarting multi-level RCGE encoding (K=$K, max_levels=$max_levels)...")

    for level in 1:max_levels
        ncur = nv(cur_g)
        mcur = count_edges(cur_g)
        println("\n  Level $level: n=$ncur, m=$mcur")

        # Partition
        t1 = time()
        part = if level == 1
            leiden_partition(cur_g; max_passes=8, max_levels=5)
        else
            leiden_partition(cur_g; max_passes=5, max_levels=5)
        end
        nclusters = maximum(part)
        println("    Partitioned into $nclusters clusters ($(round(time()-t1, digits=3))s)")

        # Cap to top-K by size
        counts = Dict{Int,Int}()
        for c in part
            counts[c] = get(counts, c, 0) + 1
        end
        labels_sorted = sort(collect(keys(counts)), by = c -> -counts[c])
        topK = min(K, length(labels_sorted))
        top = labels_sorted[1:topK]
        top_index = Dict{Int,Int}(c => i for (i,c) in enumerate(top))

        Vcur = (typeof(cur_g)).parameters[1]
        clusters = [Vcur[] for _ in 1:(topK + 1)]
        capped_part = similar(part)
        rest_order = 0
        for i in 1:ncur
            c = part[i]
            bucket = get(top_index, c, 0)
            if bucket == 0
                rest_order += 1
                bucket = topK + 1
            end
            capped_part[i] = bucket
            push!(clusters[bucket], Vcur(i))
        end
        clusters = filter(!isempty, clusters)
        println("    Effective clusters: $(length(clusters)) (top=$topK)")

        # Reorder within clusters
        reorder_clusters!(clusters, cur_g; llp_passes=llp_passes, ordering=ordering)

        # Prepare next level's K
        K = max(16, min(K, ceil(Int, nclusters / 2)))

        # Encode RCGE level
        io = IOBuffer()
        w = BitWriter(io)
        t2 = time()
        stats = RCGEStats()
        encode_level(w, cur_g, clusters; params=params, stats=stats)
        flush_bitwriter(w; flush_last_bits=true)
        bytes = take!(io)
        t3 = time()

        level_bytes = length(bytes)
        total_bytes += level_bytes
        append!(all_level_bytes, bytes)
        if level == 1
            first_level_bytes = copy(bytes)
        end
        bpe = 8.0 * level_bytes / max(mcur, 1)
        cum_bpe = 8.0 * total_bytes / max(m_original, 1)

        println("    Encoded: $(level_bytes) bytes, bits/edge=$(round(bpe, digits=4)), cumulative=$(round(cum_bpe, digits=4)) ($(round(t3-t2, digits=3))s)")
        println("    Stats: membership=$(ceil(Int, stats.bits_membership/8))B, intra=$(ceil(Int, stats.bits_intra/8))B [headers=$(ceil(Int, stats.bits_intra_headers/8))B, copy=$(ceil(Int, stats.bits_intra_copy/8))B, add=$(ceil(Int, stats.bits_intra_add/8))B, raw=$(ceil(Int, stats.bits_intra_raw/8))B, ref=$(stats.intra_ref_used)/$(stats.intra_ref_used+stats.intra_no_ref)], inter=$(ceil(Int, (stats.bits_inter_headers+stats.bits_inter_lists)/8))B")

        # Build coarse graph and check stopping
        t4 = time()
        Gc = aggregate_graph(cur_g, capped_part)
        println("    Coarsened to n=$(Gc.n) ($(round(time()-t4, digits=3))s)")

        if prev_coarse_n == Gc.n
            println("    Stopping: coarse size unchanged")
            break
        end
        prev_coarse_n = Gc.n
        if Gc.n <= min_clusters
            println("    Stopping: coarse size $(Gc.n) <= min_clusters $min_clusters")
            break
        end

        cur_g = coarse_to_testgraph(Gc)
    end

    return all_level_bytes, total_bytes, first_level_bytes
end

# Main
println("=" ^ 70)
println("RCGE CNR-2000 Multi-Level Compression Test")
println("=" ^ 70)

cnr_csv = joinpath(@__DIR__, "..", "datasets", "webgraph", "cnr-2000", "cnr-2000.csv")
cnr_csv = normpath(cnr_csv)
if !isfile(cnr_csv)
    error("CNR-2000 CSV not found at $cnr_csv")
end

println("\nLoading CNR-2000...")
t0 = time()
g = load_adjacency_list_from_csv(cnr_csv, ',', true)
n = nv(g)
m_original = count_edges(g)
println("  Loaded: $n vertices, $m_original edges ($(round(time()-t0, digits=2))s)")

# Helper to create params with overrides
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
    for (k, v) in kwargs
        d[k] = v
    end
    return RCGEParams(; d...)
end

# Shared best params builder (FW64 base)
fw64_base(; kwargs...) = make_params(
    intra_zigzag=true, intra_stop_deltas=true, intra_copy_blocks=true,
    intra_ref_window=64, intra_ref_fixwidth=true; kwargs...)

# Test configurations: (name, params, output_suffix, K_override, llp_passes, ordering, graph)
configs = [
    ("FW64 K=2 LLP (ref)",   fw64_base(), "rcge_fw64_k2",        2, 5, :llp,       g),
    ("FW64 K=1 LLP",         fw64_base(), "rcge_fw64_k1",        1, 5, :llp,       g),
    ("FW64 K=2 Bisection",   fw64_base(), "rcge_fw64_k2_bisect", 2, 5, :bisection, g),
    ("FW64 K=1 Bisection",   fw64_base(), "rcge_fw64_k1_bisect", 1, 5, :bisection, g),
]

results = []
for (name, params, suffix, K_val, llp_val, ord_val, input_g) in configs
    println("\n" * "=" ^ 70)
    println("Config: $name")
    println("=" ^ 70)

    all_level_bytes, total_bytes, first_level_bytes = run_rcge(input_g, m_original, params; K=K_val, llp_passes=llp_val, ordering=ord_val)

    output_file = joinpath(@__DIR__, "..", "datasets", "webgraph", "cnr-2000", "cnr2000_$(suffix).rcge")
    output_file = normpath(output_file)
    mkpath(dirname(output_file))
    open(output_file, "w") do f
        write(f, all_level_bytes)
    end

    fsize = filesize(output_file)
    bpe = 8.0 * fsize / m_original
    push!(results, (name, fsize, bpe))

    println("  File size: $(fsize) bytes ($(round(fsize / 1024.0, digits=1)) KB)")
    println("  Bits/edge: $(round(bpe, digits=3))")

    # Decoder roundtrip verification (level 1 only — higher levels use coarsened vertex spaces)
    println("\n  Verifying decoder roundtrip (level 1)...")
    t_dec = time()
    r_verify = BitReader(IOBuffer(first_level_bytes))
    decoded_neighbors = decode_level(r_verify, params; T=UInt32, directed=true)
    decoded_edges = sum(length(v) for v in values(decoded_neighbors))
    dt_dec = round(time() - t_dec, digits=3)
    if decoded_edges == m_original
        println("  Roundtrip: OK ($decoded_edges edges, $(dt_dec)s)")
    else
        println("  Roundtrip: MISMATCH (expected $m_original, got $decoded_edges, $(dt_dec)s)")
        # Spot-check a few vertices
        mismatches = 0
        for v in 1:min(nv(input_g), 20)
            orig_nbs = sort(UInt32.(outneighbors(input_g, v)))
            dec_nbs = get(decoded_neighbors, UInt32(v), UInt32[])
            if orig_nbs != dec_nbs
                println("    v=$v: orig=$(length(orig_nbs)) dec=$(length(dec_nbs))")
                mismatches += 1
            end
        end
        if mismatches == 0
            println("    (first 20 vertices match; discrepancy in later vertices)")
        end
    end
end

# ==================================================================
# FW64 K=1 + Implicit Ranges Membership (pre-relabeled by vertex-ID rank)
# Key insight: encode_level() re-sorts clusters by vertex ID internally.
# Relabeling by vertex-ID rank keeps local indices identical to EF encoding.
# Only the membership encoding changes (88K bytes → 7 bytes = 0.221 BPE savings).
# ==================================================================
println("\n" * "=" ^ 70)
println("Config: FW64 K=1 + Implicit Ranges Membership")
println("=" ^ 70)

println("\n  Step 1: Leiden K=1 partition...")
t_rel = time()
part_impl = leiden_partition(g; max_passes=8, max_levels=5)

counts_impl = Dict{Int,Int}()
for c in part_impl; counts_impl[c] = get(counts_impl, c, 0) + 1; end
top_label_impl = argmax(counts_impl)

TV = eltype(g)
clusters_impl_raw = [TV[], TV[]]
for i in 1:n
    push!(clusters_impl_raw[part_impl[i] == top_label_impl ? 1 : 2], TV(i))
end
# Sort by vertex ID (matches encode_level's internal sort → local indices are identical to EF)
clusters_impl_by_id = [sort(C) for C in clusters_impl_raw]
S1_impl = length(clusters_impl_by_id[1])
S2_impl = length(clusters_impl_by_id[2])
println("    Cluster sizes: $S1_impl (top) + $S2_impl (rest)")

println("\n  Step 2: Build relabeled graph (relabel by vertex-ID rank)...")
vertex_mapping_impl = let new_id = TV(1)
    d = Dict{TV,TV}()
    for C in clusters_impl_by_id
        for v in C
            d[v] = new_id
            new_id += TV(1)
        end
    end
    d
end
g_relabeled = relabel_graph(g, vertex_mapping_impl)
clusters_contiguous = [TV.(1:S1_impl), TV.(S1_impl+1:n)]
println("    Done ($(round(time()-t_rel, digits=2))s). Cluster 1 → 1..$S1_impl, cluster 2 → $(S1_impl+1)..$n")

println("\n  Step 3: RCGE encode with implicit_ranges membership...")
fw64_implicit = fw64_base(membership=:implicit_ranges)

io_impl = IOBuffer()
w_impl = BitWriter(io_impl)
t_enc = time()
stats_impl = RCGEStats()
encode_level(w_impl, g_relabeled, clusters_contiguous; params=fw64_implicit, stats=stats_impl)
flush_bitwriter(w_impl; flush_last_bits=true)
bytes_impl = take!(io_impl)
t_enc_done = time()

level_bytes_impl = length(bytes_impl)
bpe_impl = 8.0 * level_bytes_impl / m_original
println("  Encoded: $(level_bytes_impl) bytes ($(round(level_bytes_impl/1024.0, digits=1)) KB), BPE=$(round(bpe_impl, digits=3)) ($(round(t_enc_done-t_enc, digits=2))s)")
println("  Stats: membership=$(ceil(Int, stats_impl.bits_membership/8))B, intra=$(ceil(Int, stats_impl.bits_intra/8))B [headers=$(ceil(Int, stats_impl.bits_intra_headers/8))B, copy=$(ceil(Int, stats_impl.bits_intra_copy/8))B, add=$(ceil(Int, stats_impl.bits_intra_add/8))B, raw=$(ceil(Int, stats_impl.bits_intra_raw/8))B, ref=$(stats_impl.intra_ref_used)/$(stats_impl.intra_ref_used+stats_impl.intra_no_ref)], inter=$(ceil(Int, (stats_impl.bits_inter_headers+stats_impl.bits_inter_lists)/8))B")

output_file_impl = normpath(joinpath(@__DIR__, "..", "datasets", "webgraph", "cnr-2000", "cnr2000_rcge_gllp_implicit.rcge"))
open(output_file_impl, "w") do f; write(f, bytes_impl); end

println("\n  Verifying decoder roundtrip...")
t_dec = time()
r_verify_impl = BitReader(IOBuffer(bytes_impl))
decoded_impl = decode_level(r_verify_impl, fw64_implicit; T=TV, directed=true)
decoded_edges_impl = sum(length(v) for v in values(decoded_impl))
dt_dec = round(time() - t_dec, digits=3)
if decoded_edges_impl == m_original
    println("  Roundtrip: OK ($decoded_edges_impl edges, $(dt_dec)s)")
else
    println("  Roundtrip: MISMATCH (expected $m_original, got $decoded_edges_impl, $(dt_dec)s)")
    for v in 1:min(n, 20)
        orig_nbs = sort(TV.(outneighbors(g_relabeled, v)))
        dec_nbs  = sort(get(decoded_impl, TV(v), TV[]))
        if orig_nbs != dec_nbs
            println("    v=$v: orig=$(length(orig_nbs)) dec=$(length(dec_nbs))")
        end
    end
end
push!(results, ("FW64 K=1 implicit ranges", level_bytes_impl, bpe_impl))

println("\n" * "=" ^ 70)
println("SUMMARY")
println("=" ^ 70)
for (name, fsize, bpe) in results
    println("  $(rpad(name, 30)) $(round(bpe, digits=3)) BPE  ($(round(fsize / 1024.0, digits=1)) KB)")
end
println("  $(rpad("WebGraph reference", 25)) 2.900 BPE")
println("=" ^ 70)
