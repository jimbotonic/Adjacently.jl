# feature_ablation.jl — tab:ablation.
#
# STATUS: the published feature-ablation table does NOT reproduce and its claims
# are withdrawn. Measured on cnr-2000 under the LLP ordering with the shipped
# encoder (2026-09-06), every row differs from the published value and three rows
# have the WRONG SIGN:
#     published  + VLC header   -0.195   measured  +0.074  (costs bits)
#     published  + LR-split     -0.064   measured  +0.026  (costs bits)
#     published  + multi-ref    -0.015   measured   0.000  (no effect at all)
# The published endpoint (2.493, "+ low-degree reference search", -0.30 bpe) is
# not reachable: it required REF_ENCODING_TH=1, and that constant sits in
# find_best_reference_greedy_cost, a code path the BG encoder does not call
# (BG uses _greedy_vertex_search). Setting it to 1 leaves output bit-identical.
#
# This driver measures the ladder that IS reproducible, so the table can be
# rebuilt honestly or dropped on evidence rather than assumption.
using Pkg; Pkg.activate(normpath(joinpath(@__DIR__, "..", "..")); io=devnull)
using Adjacently
using Adjacently.MGS: write_bg_mgs3_graph
include(joinpath(@__DIR__, "config.jl")); using .ReproConfig
using LightGraphs: ne
using Printf

const DATASET = get(ENV, "DATASET", "cnr-2000")
g0 = load_dataset(DATASET)
g  = apply_ordering(g0, "llp", DATASET; seed=0)
m  = ne(g); tmp = mktempdir()

rows = [
 ("baseline (zeta-3, w=7)",   (ie=:zeta,      w=7,  cb=false, sd=false, lr=false, mr=false, ah=false)),
 ("+ copy-blocks",            (ie=:zeta,      w=7,  cb=true,  sd=false, lr=false, mr=false, ah=false)),
 ("+ larger window (w=64)",   (ie=:zeta,      w=64, cb=true,  sd=false, lr=false, mr=false, ah=false)),
 ("+ Fibonacci encoding",     (ie=:fibonacci, w=64, cb=true,  sd=false, lr=false, mr=false, ah=false)),
 ("+ STOP-terminated deltas", (ie=:fibonacci, w=64, cb=true,  sd=true,  lr=false, mr=false, ah=false)),
 ("+ VLC/adaptive header",    (ie=:fibonacci, w=64, cb=true,  sd=true,  lr=false, mr=false, ah=true )),
 ("+ LR-split residuals",     (ie=:fibonacci, w=64, cb=true,  sd=true,  lr=true,  mr=false, ah=true )),
 ("+ multi-reference",        (ie=:fibonacci, w=64, cb=true,  sd=true,  lr=true,  mr=true,  ah=true )),
]
@printf("%s under LLP ordering, %d edges\n\n", DATASET, m)
@printf("%-28s %9s %9s\n", "cumulative feature", "bpe", "delta"); println("-"^50)
prev = NaN
for (i, (name, c)) in enumerate(rows)
    b = joinpath(tmp, "r$i")
    write_bg_mgs3_graph(g, b; integer_encoding=c.ie, ref_window_size=c.w, copy_blocks=c.cb,
        stop_deltas=c.sd, lr_split=c.lr, multi_ref=c.mr, adaptive_header=c.ah, cost_model=0)
    v = 8 * filesize(b * ".mgz") / m
    @printf("%-28s %9.4f %9s\n", name, v, isnan(prev) ? "---" : @sprintf("%+.4f", v - prev))
    global prev = v
end
println("\nNOTE: rows that increase bpe are genuine regressions of the current encoder,")
println("not measurement error. See the header of this file.")
