# verify_cg_config.jl — regression guard for the recovered CG configuration.
#
# The paper's cnr-2000 CG cell (2.329) was for months reproducible only from a
# configuration buried in test/run_tests_index_mode.jl. This asserts that the
# configuration in configs/encoders.json still reproduces it, so the recipe can
# never silently drift out of the reproduction path again.
using Pkg; Pkg.activate(normpath(joinpath(@__DIR__, "..", "..")); io=devnull)
using Adjacently
include(joinpath(@__DIR__, "config.jl")); using .ReproConfig
using Printf

const EXPECTED = 2.3286
const TOL      = 0.001

g = load_dataset("cnr-2000")
tmp = mktempdir()
bpe = cg_k2_bpe(g, joinpath(tmp, "cg_k2"))
@printf("cnr-2000 CG K=2 (paper config): %.4f   expected %.4f\n", bpe, EXPECTED)
if abs(bpe - EXPECTED) <= TOL
    println("PASS — configs/encoders.json reproduces the published CG cell.")
else
    println("FAIL — configuration has drifted; the published CG cell is no longer reproducible.")
    exit(1)
end
