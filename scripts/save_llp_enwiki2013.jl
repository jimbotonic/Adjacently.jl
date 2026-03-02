#!/usr/bin/env julia
#
# Compute and save global LLP ordering for enwiki-2013 (5 passes)
#

using Pkg
Pkg.activate(normpath(joinpath(@__DIR__, "..")))

using Logging
global_logger(ConsoleLogger(stderr, Logging.Info))

using LightGraphs: nv, ne, eltype
using Adjacently
using Adjacently.IO: load_adjacency_list_from_csv
using Adjacently.Relabeling: relabel_vertices_llp, save_llp_ordering

const DS = "enwiki-2013"
const DS_DIR = normpath(joinpath(@__DIR__, "..", "datasets", "webgraph", DS))
const DS_CSV = joinpath(DS_DIR, "$DS.csv")
const PREFIX = replace(DS, "-" => "")
const PASSES = 5

isfile(DS_CSV) || error("CSV not found: $DS_CSV")

@info "Loading $DS..."
t0 = time()
g = load_adjacency_list_from_csv(DS_CSV, ',', true)
n = Int(nv(g))
T = eltype(g)
@info "  Loaded: $n vertices, $(ne(g)) edges ($(round(time()-t0, digits=1))s)"

llp_path = joinpath(DS_DIR, "$(PREFIX)_llp_p$(PASSES).bin")
@info "Computing global LLP (passes=$PASSES)..."
t1 = time()
mapping = relabel_vertices_llp(g, :sym; passes=PASSES)
dt = round(time() - t1, digits=1)
@info "  LLP done in $(dt)s"

save_llp_ordering(llp_path, mapping, n)
@info "Done: $llp_path"
