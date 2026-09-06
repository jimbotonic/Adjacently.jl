# ord_ablation_all.jl — tab:ord-ablation for ALL datasets, including the LAW crawls.
#
# The previous ord_ablation.jl covered only the four committed SNAP-style graphs.
# cnr-2000 and in-2004 had no reproduction driver, and their published cells were
# consequently wrong for two paper versions: the "Orig." row for those two graphs
# was measured on an LLP re-ordering, not on the distributed ordering, because the
# table caption assumed the two coincide. They do not. This driver measures all
# three orderings separately for every dataset.
#
#   DATASETS=cnr-2000,in-2004  SEEDS=0  julia ord_ablation_all.jl
using Pkg; Pkg.activate(normpath(joinpath(@__DIR__, "..", "..")); io=devnull)
using Adjacently
include(joinpath(@__DIR__, "config.jl")); using .ReproConfig
using LightGraphs: nv, ne
using Printf, Dates, Statistics

const SEEDS = parse.(Int, split(get(ENV, "SEEDS", "0"), ","))
const DS    = split(get(ENV, "DATASETS", "eat,arxiv-hep-ph,amazon-0601,web-google,cnr-2000,in-2004"), ",")
const ORDS  = split(get(ENV, "ORDERINGS", "original,llp,leiden_llp"), ",")
const OUT   = joinpath(@__DIR__, get(ENV, "OUT", "ord_ablation_all.tsv"))

isfile(OUT) || open(OUT, "w") do io
    println(io, "timestamp\tdataset\tordering\tencoder\tbpe\tseed")
end
for ds in DS
    g0 = load_dataset(ds); tmp = mktempdir()
    @printf("\n=== %s : %d v / %d e ===\n", ds, nv(g0), ne(g0))
    @printf("  %-14s %9s %9s %9s\n", "ordering", "BG", "CS", "CG(K=2)")
    for ord in ORDS, seed in SEEDS
        g = apply_ordering(g0, ord, ds; seed=seed)
        bg = encode_bpe(g, "bg", ds, joinpath(tmp, "bg_$(ord)_$(seed)"))
        cs = encode_bpe(g, "cs", ds, joinpath(tmp, "cs_$(ord)_$(seed)"))
        cg = cg_k2_bpe(g, joinpath(tmp, "cg_$(ord)_$(seed)"))
        @printf("  %-14s %9.4f %9.4f %9.4f\n", ord, bg, cs, cg); flush(stdout)
        open(OUT, "a") do io
            for (e, v) in (("bg", bg), ("cs", cs), ("cg", cg))
                @printf(io, "%s\t%s\t%s\t%s\t%.6f\t%d\n", now(), ds, ord, e, v, seed)
            end
        end
    end
    rm(tmp; recursive=true, force=true)
end
println("\nWrote $OUT")
