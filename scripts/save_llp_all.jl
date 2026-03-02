#!/usr/bin/env julia
#
# Compute and save global LLP ordering for a given dataset
# Usage: julia --project scripts/save_llp_all.jl DATASET [PASSES]
#

using Pkg
Pkg.activate(normpath(joinpath(@__DIR__, "..")))

using Logging
global_logger(ConsoleLogger(stderr, Logging.Info))

using LightGraphs: nv, ne, eltype
using Adjacently
using Adjacently.IO: load_adjacency_list_from_csv
using Adjacently.Relabeling: relabel_vertices_llp, save_llp_ordering

const DS = length(ARGS) >= 1 ? ARGS[1] : error("Usage: julia save_llp_all.jl DATASET [PASSES]")
const PASSES = length(ARGS) >= 2 ? parse(Int, ARGS[2]) : 10
const DS_DIR = normpath(joinpath(@__DIR__, "..", "datasets", "webgraph", DS))
const DS_CSV = joinpath(DS_DIR, "$DS.csv")
const PREFIX = replace(DS, "-" => "")

isfile(DS_CSV) || error("CSV not found: $DS_CSV")

llp_path = joinpath(DS_DIR, "$(PREFIX)_llp_p$(PASSES).bin")
if isfile(llp_path)
    @info "Already exists: $llp_path ($(filesize(llp_path)) bytes) — skipping"
    exit(0)
end

@info "Loading $DS..."
t0 = time()
g = load_adjacency_list_from_csv(DS_CSV, ',', true)
n = Int(nv(g))
T = eltype(g)
@info "  Loaded: $n vertices, $(ne(g)) edges ($(round(time()-t0, digits=1))s)"

@info "Computing global LLP (passes=$PASSES)..."
t1 = time()
mapping = relabel_vertices_llp(g, :sym; passes=PASSES)
dt = round(time() - t1, digits=1)
@info "  LLP done in $(dt)s"

save_llp_ordering(llp_path, mapping, n)
@info "Done: $llp_path"
