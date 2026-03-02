#!/usr/bin/env julia
#
# Verify roundtrip of enwiki2013_rcge.mgz (RCGE K=1 with LLP ordering)
#

using Pkg
Pkg.activate(normpath(joinpath(@__DIR__, "..")))

using Logging
global_logger(ConsoleLogger(stderr, Logging.Info))

using LightGraphs: nv, ne, outneighbors, vertices, eltype
using Adjacently
using Adjacently.IO: load_adjacency_list_from_csv
using Adjacently.Relabeling: relabel_graph, load_llp_ordering
using Adjacently.MGS: load_compressed_mgs3_graph
using Adjacently.Compression.RCGE: RCGEParams

const DS = "enwiki-2013"
const DS_DIR = normpath(joinpath(@__DIR__, "..", "datasets", "webgraph", DS))
const DS_CSV = joinpath(DS_DIR, "$DS.csv")
const PREFIX = replace(DS, "-" => "")

rcge_mgz = joinpath(DS_DIR, "$(PREFIX)_rcge.mgz")
llp_path = joinpath(DS_DIR, "$(PREFIX)_llp_p10.bin")

isfile(rcge_mgz) || error("RCGE file not found: $rcge_mgz")
isfile(llp_path) || error("LLP ordering not found: $llp_path")

# ── Load original graph and apply LLP ────────────────────────────────────────

@info "Loading $DS..."
t0 = time()
g = load_adjacency_list_from_csv(DS_CSV, ',', true)
T = eltype(g)
n = nv(g); m = ne(g)
@info "  Loaded: $n vertices, $m edges ($(round(time()-t0, digits=1))s)"

@info "Loading LLP ordering..."
vertex_map = load_llp_ordering(llp_path, T)
g_rel = relabel_graph(g, vertex_map)
@info "  Relabeled graph: $(nv(g_rel)) vertices, $(ne(g_rel)) edges"

# Free original graph
g = nothing
GC.gc()

# ── Load compressed file ─────────────────────────────────────────────────────

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

@info "Loading compressed file: $rcge_mgz"
t_dec = time()
g_loaded = load_compressed_mgs3_graph(rcge_mgz; rcge_params=rcge_params)
dt_dec = round(time() - t_dec, digits=1)
@info "  Decoded: $(nv(g_loaded)) vertices, $(ne(g_loaded)) edges ($(dt_dec)s)"

# ── Verify vertex-by-vertex ──────────────────────────────────────────────────

@info "Verifying roundtrip (vertex-by-vertex)..."
t_ver = time()
n_checked = 0
mismatches = 0
for v in vertices(g_rel)
    orig = sort(collect(outneighbors(g_rel, v)))
    decoded = sort(collect(outneighbors(g_loaded, v)))
    if orig != decoded
        mismatches += 1
        if mismatches <= 5
            @warn "  Mismatch at vertex $v: expected $(length(orig)) neighbors, got $(length(decoded))"
        end
    end
    n_checked += 1
    if n_checked % 1_000_000 == 0
        @info "  Checked $n_checked / $(nv(g_rel)) vertices..."
    end
end
dt_ver = round(time() - t_ver, digits=1)

bpe = round(8.0 * filesize(rcge_mgz) / m, digits=4)

if mismatches == 0
    @info "  Roundtrip: PASSED ($n_checked vertices, $bpe BPE, verify $(dt_ver)s)"
else
    @warn "  Roundtrip: FAILED ($mismatches / $n_checked vertices mismatched)"
end

println("\n" * "=" ^ 70)
println("SUMMARY — $DS RCGE K=1 (LLP p=10)")
println("=" ^ 70)
println("  RCGE K=1 LLP: $bpe BPE  ($(filesize(rcge_mgz)) bytes)")
println("  CS zeta LLP:  13.072 BPE (previous)")
println("  WebGraph BV:  13.114 BPE (reference)")
println("  Decode time:  $(dt_dec)s")
println("  Verify time:  $(dt_ver)s")
println("  Mismatches:   $mismatches")
println("=" ^ 70)
