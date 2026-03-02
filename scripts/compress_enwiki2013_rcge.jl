#!/usr/bin/env julia
#
# Compress enwiki-2013 with RCGE best known parameters
#

using Pkg
Pkg.activate(normpath(joinpath(@__DIR__, "..")))

using Logging
global_logger(ConsoleLogger(stderr, Logging.Info))

using LightGraphs: nv, ne, outneighbors, vertices, eltype
using Adjacently
using Adjacently.IO: load_adjacency_list_from_csv
using Adjacently.Relabeling: relabel_graph
using Adjacently.Clustering: leiden_partition, auto_select_K
using Adjacently.Graph: subgraph
using Adjacently.MGS: write_rcge_mgs3_graph, load_compressed_mgs3_graph
using Adjacently.Compression.RCGE: RCGEParams

const DS = "enwiki-2013"
const DS_DIR = normpath(joinpath(@__DIR__, "..", "datasets", "webgraph", DS))
const DS_CSV = joinpath(DS_DIR, "$DS.csv")
const PREFIX = replace(DS, "-" => "")

isfile(DS_CSV) || error("CSV not found: $DS_CSV")

println("=" ^ 70)
println("RCGE compression — $DS")
println("=" ^ 70)

# ── Load graph ───────────────────────────────────────────────────────────────

@info "Loading $DS..."
t0 = time()
g = load_adjacency_list_from_csv(DS_CSV, ',', true)
n = nv(g); m = ne(g)
T = eltype(g)
@info "  Loaded: $n vertices, $m edges ($(round(time()-t0, digits=1))s)"

# ── Auto-select K ────────────────────────────────────────────────────────────

@info "Running auto_select_K..."
t1 = time()
K, part = auto_select_K(g)
@info "  auto K=$K ($(round(time()-t1, digits=1))s)"

# ── Build K-partition ────────────────────────────────────────────────────────

@info "Building Leiden K=$K partition..."
t2 = time()

# Count community sizes and sort descending
counts = Dict{Int,Int}()
for c in part; counts[c] = get(counts, c, 0) + 1; end
sorted_labels = sort(collect(keys(counts)); by=l->counts[l], rev=true)

# Top K-1 communities get their own cluster, rest merge into cluster K
top_labels = Set(sorted_labels[1:min(K-1, length(sorted_labels))])
label_to_cluster = let ci = 1
    d = Dict{Int,Int}()
    for l in sorted_labels
        if l in top_labels
            d[l] = ci; ci += 1
        else
            d[l] = K
        end
    end
    d
end

clusters_raw = [T[] for _ in 1:K]
for i in 1:n
    push!(clusters_raw[label_to_cluster[part[i]]], T(i))
end
for C in clusters_raw; sort!(C); end
@info "  Cluster sizes: $(length.(clusters_raw)) ($(round(time()-t2, digits=1))s)"

# ── Build vertex mapping and relabel ─────────────────────────────────────────

@info "Relabeling graph..."
t3 = time()
rcge_vertex_map = let new_id = T(1)
    d = Dict{T,T}()
    for C in clusters_raw
        for v in C
            d[v] = new_id
            new_id += T(1)
        end
    end
    d
end

g_rcge = relabel_graph(g, rcge_vertex_map)
m_rcge = ne(g_rcge)

# Build implicit-range cluster vectors
clusters_impl = let off = T(0)
    result = Vector{Vector{T}}()
    for C in clusters_raw
        sz = T(length(C))
        push!(result, T.(off+1:off+sz))
        off += sz
    end
    result
end
@info "  Relabeled in $(round(time()-t3, digits=1))s"

# ── RCGE parameters ──────────────────────────────────────────────────────────

# Use best known parameters (window=8 like in-2004 since this is a large graph)
rcge_params = RCGEParams(
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

const PROGRESS_INTERVAL = 2.0  # seconds between updates
last_progress_time = Ref(time())
vertices_done = Ref(0)
t_encode_start = Ref(0.0)

function progress_callback(idx_local::Int, cluster_size::Int, cluster_idx::Int, num_clusters::Int)
    now = time()
    if idx_local == 1
        # New cluster starting
        vertices_done[] += 0  # no-op, just mark
    end
    if idx_local == cluster_size || now - last_progress_time[] >= PROGRESS_INTERVAL
        total_done = vertices_done[] + idx_local
        elapsed = now - t_encode_start[]
        rate = total_done / max(elapsed, 0.001)
        remaining = (n - total_done) / max(rate, 0.001)

        pct = round(100.0 * total_done / n, digits=1)
        bar_width = 40
        filled = round(Int, bar_width * total_done / n)
        bar = "█" ^ filled * "░" ^ (bar_width - filled)

        print("\r  [$bar] $pct%  cluster $cluster_idx/$num_clusters  $(round(Int, total_done))/$n vertices  $(round(rate, digits=0)) v/s  ETA $(round(Int, remaining))s   ")
        last_progress_time[] = now

        if idx_local == cluster_size
            vertices_done[] += cluster_size
        end
    end
end

# ── Encode ───────────────────────────────────────────────────────────────────

k_suffix = K == 1 ? "" : "_k$K"
rcge_base = joinpath(DS_DIR, "$(PREFIX)_rcge$(k_suffix)")
rcge_mgz = rcge_base * ".mgz"

@info "Encoding RCGE K=$K → $rcge_mgz"
t_encode_start[] = time()

write_rcge_mgs3_graph(g_rcge, rcge_base, clusters_impl;
    params=rcge_params, progress=progress_callback)

println()  # newline after progress bar
dt_enc = round(time() - t_encode_start[], digits=1)
bpe = round(8.0 * filesize(rcge_mgz) / m_rcge, digits=4)
@info "  Encoded: $bpe BPE ($(filesize(rcge_mgz)) bytes, $(dt_enc)s)"

# ── Verify roundtrip ─────────────────────────────────────────────────────────

@info "Verifying roundtrip..."
t_dec = time()
g_loaded = load_compressed_mgs3_graph(rcge_mgz; rcge_params=rcge_params)
dt_dec = round(time() - t_dec, digits=1)

if nv(g_loaded) == nv(g_rcge) && ne(g_loaded) == ne(g_rcge)
    @info "  Roundtrip: PASSED ($(nv(g_loaded)) vertices, $(ne(g_loaded)) edges, decode $(dt_dec)s)"
else
    @warn "  Roundtrip: MISMATCH! Expected $(nv(g_rcge))v/$(ne(g_rcge))e, got $(nv(g_loaded))v/$(ne(g_loaded))e"
end

# ── Summary ──────────────────────────────────────────────────────────────────

println("\n" * "=" ^ 70)
println("SUMMARY — $DS RCGE K=$K")
println("=" ^ 70)
println("  RCGE K=$K:   $bpe BPE  ($(filesize(rcge_mgz)) bytes)")
println("  WebGraph BV: 12.639 BPE (reference)")
println("  Encode time: $(dt_enc)s")
println("  Decode time: $(dt_dec)s")
println("=" ^ 70)
