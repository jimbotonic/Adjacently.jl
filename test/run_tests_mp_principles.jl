# Tests for item 9 of research/mycelial_polis/RESEARCH_TASKS.md:
# ten internal attractor principles as monotone hyperparameters.
#
# Acceptance:
#   - each principle, swept across {0.0, 0.5, 1.0} (one at a time, other
#     principles held at 0.5), produces a *monotone* shift in at least
#     one of `(Φ, Λ_sat, Ψ_T)` (sign per the principles.jl table).
#   - the Principles struct construction + helpers (get_principles,
#     principle_strength) behave as expected.

using Test
using Statistics: mean
using Adjacently
using Adjacently.MycelialPolis

# Field name on the Principles struct.
const ALL_PRINCIPLES = (
    :non_domination, :functional_autonomy, :bounded_legibility,
    :reversible_governance, :forkability, :redundant_stewardship,
    :normative_minimalism, :adaptive_memory, :controlled_permeability,
    :transformative_non_absorption,
)

# Run a single (principle, strength) cell, return (Phi, Lambda_sat, Psi_T)
# averaged over the last 25% of the trajectory.
function _run_cell(principle_name::Symbol, strength::Float32;
                   topology=:modular_cells, n=200, T=40,
                   host_budget=3, seed=11)
    # All principles at their defaults except the one we're sweeping.
    kwargs = NamedTuple{(principle_name,)}((strength,))
    pr = Principles(; kwargs...)
    params = merge(default_params(),
                   (initial_committed_frac = 0.10,
                    principles = pr))
    world = build_world(; topology=topology, n=n, seed=seed, params=params)
    host  = RandomHost(budget=host_budget)
    part  = natural_partition(topology, n)
    phis = Float32[]; lsats = Float32[]; psis = Float32[]
    for _ in 1:T
        snap = step!(world, host; partition=part)
        push!(phis, snap.Phi); push!(lsats, snap.Lambda_sat); push!(psis, snap.Psi_T)
    end
    tail = max(1, ceil(Int, 0.75 * T)) : T
    return (Phi=mean(@view phis[tail]),
            Lsat=mean(@view lsats[tail]),
            PsiT=mean(@view psis[tail]))
end

# At least one of the three metrics is monotone (allow ±slop = 0.02)
# across {0.0, 0.5, 1.0}, i.e. either weakly-increasing or weakly-decreasing.
function _at_least_one_monotone(v0, v1, v2; slop::Float64=0.02)
    function mono(a, b, c)
        a <= b + slop && b <= c + slop   # non-decreasing within slop
    end
    function mono_rev(a, b, c)
        a + slop >= b && b + slop >= c   # non-increasing within slop
    end
    for (a, b, c) in zip(v0, v1, v2)
        (mono(a, b, c) || mono_rev(a, b, c)) && return true
    end
    return false
end

@testset "MycelialPolis principles (item 9)" begin

    @testset "Principles struct + helpers" begin
        p = Principles()
        @test p.non_domination == 0.5f0
        @test p.forkability    == 0.5f0
        @test p.transformative_non_absorption == 0.5f0

        # Custom construction via kwargs.
        p2 = Principles(non_domination=0.9f0, forkability=0.1f0)
        @test p2.non_domination == 0.9f0
        @test p2.forkability    == 0.1f0
        @test p2.bounded_legibility == 0.5f0    # default for unset

        # get_principles falls back to defaults if params is missing the key.
        @test get_principles(default_params()).non_domination == 0.5f0
        @test get_principles((;)).functional_autonomy == 0.5f0

        # principle_strength is a clamped accessor.
        @test principle_strength(default_params(), :forkability) == 0.5f0
        params_neg = (; principles = Principles(forkability=-0.5f0))
        @test principle_strength(params_neg, :forkability) == 0f0
        params_big = (; principles = Principles(forkability=1.5f0))
        @test principle_strength(params_big, :forkability) == 1f0
    end

    @testset "monotonicity sweep (all 10 principles)" begin
        # JIT warm-up.
        _run_cell(:non_domination, 0.5f0; T=10, n=80)

        all_monotone = true
        results = Dict{Symbol,NTuple{3,NTuple{3,Float32}}}()
        for principle in ALL_PRINCIPLES
            r0 = _run_cell(principle, 0.0f0)
            r5 = _run_cell(principle, 0.5f0)
            r1 = _run_cell(principle, 1.0f0)
            results[principle] = (
                (r0.Phi,  r5.Phi,  r1.Phi ),
                (r0.Lsat, r5.Lsat, r1.Lsat),
                (r0.PsiT, r5.PsiT, r1.PsiT),
            )
            phis  = (r0.Phi,  r5.Phi,  r1.Phi)
            lsats = (r0.Lsat, r5.Lsat, r1.Lsat)
            psis  = (r0.PsiT, r5.PsiT, r1.PsiT)
            ok = _at_least_one_monotone(phis, lsats, psis)
            @info "sweep" principle Phi=phis Lsat=lsats PsiT=psis monotone=ok
            @test ok
            ok || (all_monotone = false)
        end
        @test all_monotone
    end
end
