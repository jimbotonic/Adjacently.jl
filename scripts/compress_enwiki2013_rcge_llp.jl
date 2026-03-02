#!/usr/bin/env julia
#
# Compress enwiki-2013 with RCGE K=1 + global LLP ordering (memory-lean)
#

using Pkg
Pkg.activate(normpath(joinpath(@__DIR__, "..")))

using Logging
global_logger(ConsoleLogger(stderr, Logging.Info))

using LightGraphs: nv, ne, outneighbors, vertices, eltype
using Adjacently
using Adjacently.IO: load_adjacency_list_from_csv
using Adjacently.Relabeling: relabel_graph, load_llp_ordering
using Adjacently.MGS: write_rcge_mgs3_graph, load_compressed_mgs3_graph
using Adjacently.Compression.RCGE: RCGEParams

const DS = "enwiki-2013"
const DS_DIR = normpath(joinpath(@__DIR__, "..", "datasets", "webgraph", DS))
const DS_CSV = joinpath(DS_DIR, "$DS.csv")
const PREFIX = replace(DS, "-" => "")
const LLP_PASSES = parse(Int, get(ENV, "LLP_PASSES", "10"))
const WINDOW = parse(Int, get(ENV, "WINDOW", "8"))
const INTERVALS = parse(Bool, get(ENV, "INTERVALS", "false"))

llp_path = joinpath(DS_DIR, "$(PREFIX)_llp_p$(LLP_PASSES).bin")
isfile(DS_CSV) || error("CSV not found: $DS_CSV")
isfile(llp_path) || error("LLP ordering not found: $llp_path — run save_llp_enwiki2013.jl first")

# ── Load graph ───────────────────────────────────────────────────────────────

@info "Loading $DS..."
t0 = time()
g_original = load_adjacency_list_from_csv(DS_CSV, ',', true)
T = eltype(g_original)
n = Int(nv(g_original)); m = Int(ne(g_original))
@info "  Loaded: $n vertices, $m edges ($(round(time()-t0, digits=1))s)"

# ── Apply LLP ordering ──────────────────────────────────────────────────────

@info "Loading LLP ordering (passes=$LLP_PASSES)..."
vertex_map = load_llp_ordering(llp_path, T)

@info "Relabeling graph..."
t1 = time()
g = relabel_graph(g_original, vertex_map)
@info "  Relabeled in $(round(time()-t1, digits=1))s"

# Free original graph and mapping
g_original = nothing
vertex_map = nothing
GC.gc()
@info "  Freed original graph ($(round(Sys.free_memory()/1024^3, digits=1)) GB free)"

# ── RCGE parameters ──────────────────────────────────────────────────────────

rcge_params = RCGEParams(
    L=128,
    varint=:fibonacci, count_varint=:fibonacci,
    gap=:fibonacci, degree=:elias_delta,
    undirected_pairs=false,
    perm_strategy=:blockpos, membership=:implicit_ranges,
    inter_strategy=:perm,
    intra_ref_enabled=true, intra_ref_window=WINDOW,
    intra_ref_rle=false,
    intra_block_try=false,
    positions_mode=:delta, additions_mode=:delta,
    min_cluster_density=0.0,
    intra_intervals=INTERVALS, intra_mil=4, intra_greedy_mil=false,
    intra_zigzag=true, intra_stop_deltas=true,
    intra_copy_blocks=true, intra_copy_adaptive=true,
    intra_ref_fixwidth=true, intra_ref_vlc=false,
    intra_add_adaptive=true, intra_raw_adaptive=true,
    intra_adapt_mil=2,
)

# ── Encode ───────────────────────────────────────────────────────────────────

clusters = [T.(1:n)]

rcge_base = joinpath(DS_DIR, "$(PREFIX)_rcge")
rcge_mgz = rcge_base * ".mgz"

@info "Encoding RCGE K=1 → $rcge_mgz"
t_enc = time()
write_rcge_mgs3_graph(g, rcge_base, clusters; params=rcge_params)
dt_enc = round(time() - t_enc, digits=1)

bpe = round(8.0 * filesize(rcge_mgz) / m, digits=4)
@info "  Encoded: $bpe BPE ($(filesize(rcge_mgz)) bytes, $(dt_enc)s)"

# ── Verify roundtrip ─────────────────────────────────────────────────────────

@info "Verifying roundtrip..."
t_dec = time()
g_loaded = load_compressed_mgs3_graph(rcge_mgz; rcge_params=rcge_params)
dt_dec = round(time() - t_dec, digits=1)

@info "  Decoded: $(nv(g_loaded)) vertices, $(ne(g_loaded)) edges ($(dt_dec)s)"

mismatches = 0
for v in vertices(g)
    orig = sort(collect(outneighbors(g, v)))
    decoded = sort(collect(outneighbors(g_loaded, v)))
    if orig != decoded
        mismatches += 1
        if mismatches <= 5
            @warn "  Mismatch at vertex $v"
        end
    end
end

if mismatches == 0
    @info "  Roundtrip: PASSED"
else
    @warn "  Roundtrip: FAILED ($mismatches mismatches)"
end

# ── Summary ──────────────────────────────────────────────────────────────────

println("\n" * "=" ^ 70)
println("SUMMARY — $DS RCGE K=1 (LLP p=$LLP_PASSES)")
println("=" ^ 70)
println("  RCGE K=1 LLP: $bpe BPE  ($(filesize(rcge_mgz)) bytes)")
println("  CS zeta LLP:  13.072 BPE")
println("  WebGraph BV:  13.114 BPE")
println("  WebGraph hc:  12.639 BPE")
println("  Encode time:  $(dt_enc)s")
println("  Decode time:  $(dt_dec)s")
println("=" ^ 70)
