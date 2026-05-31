# Tests for item 10 of research/mycelial_polis/RESEARCH_TASKS.md:
# principle ablations (OAT + pair-out + greedy forward selection).
#
# Acceptance:
#   - principles_with / principles_only build the expected Principles structs;
#   - summarise_basin returns a NamedTuple with the expected fields and A_S
#     is finite;
#   - one-cell OAT pass (just 2 principles disabled, smoke tier) returns
#     a sensible NamedTuple per principle;
#   - greedy_forward terminates (≤ max_size picks) and the A_S history is
#     monotone non-decreasing.

using Test
using Adjacently
using Adjacently.MycelialPolis

@testset "MycelialPolis ablations (item 10)" begin

    @testset "principles_with / principles_only construction" begin
        # principles_with disables the named ones, leaves the rest at 0.5
        p1 = principles_with([:non_domination, :forkability])
        @test p1.non_domination == 0f0
        @test p1.forkability    == 0f0
        @test p1.functional_autonomy == 0.5f0
        @test p1.bounded_legibility  == 0.5f0

        # principles_only enables only the named ones, rest at 0
        p2 = principles_only([:functional_autonomy, :bounded_legibility])
        @test p2.functional_autonomy == 0.5f0
        @test p2.bounded_legibility  == 0.5f0
        @test p2.non_domination      == 0f0
        @test p2.forkability         == 0f0

        # ALL_PRINCIPLE_NAMES is the 10-name tuple
        @test length(ALL_PRINCIPLE_NAMES) == 10
        @test :non_domination       in ALL_PRINCIPLE_NAMES
        @test :transformative_non_absorption in ALL_PRINCIPLE_NAMES
    end

    @testset "summarise_basin produces sane numbers" begin
        # Toy basin result NamedTuple — same shape as estimate_basins().
        fake = (
            n_samples = 100,
            counts = Dict(:extinction => 40,
                          :dormant_persistence => 20,
                          :parallel_coexistence => 30,
                          :unclassified => 10),
            t_attractors = Dict(:extinction => [10, 12], :dormant_persistence => [20],
                                :parallel_coexistence => [8], :unclassified => [15]),
            final_phis = Dict(:extinction => Float32[0.01, 0.02],
                              :dormant_persistence => Float32[0.20],
                              :parallel_coexistence => Float32[0.70, 0.75],
                              :unclassified => Float32[0.45]),
            tier = :smoke, topology = :modular_cells,
            visibility = 1.0f0, seed = 0, elapsed_s = 1.0,
        )
        s = summarise_basin(fake)
        @test s.basin_parallel   ≈ 0.30 atol=1e-6
        @test s.basin_dormant    ≈ 0.20 atol=1e-6
        @test s.basin_extinction ≈ 0.40 atol=1e-6
        @test isfinite(s.A_S)
        # A_S = basin_par + mean_phi_surv + 0.5*basin_dorm - basin_ext
        # mean_phi_surv = mean([0.20, 0.70, 0.75]) ≈ 0.55
        # A_S = 0.30 + 0.55 + 0.10 - 0.40 = 0.55
        @test s.A_S ≈ 0.55f0 atol=5e-3
    end

    @testset "1-OAT smoke (single principle disabled)" begin
        # Run a single OAT cell at smoke tier with a fixed-small budget.
        # Just verify the result NamedTuple has the expected fields and
        # the basin parallel fraction is sane.
        baseline = summarise_basin(estimate_basins(; tier=:smoke,
                                                     topology=:modular_cells,
                                                     seed=99, visibility=1.0f0,
                                                     principles=principles_with(Symbol[])))
        for p in (:non_domination, :transformative_non_absorption)
            pr_dis = principles_with([p])
            s = summarise_basin(estimate_basins(; tier=:smoke,
                                                  topology=:modular_cells,
                                                  seed=99, visibility=1.0f0,
                                                  principles=pr_dis))
            @info "OAT cell" principle=p basin_par=s.basin_parallel A_S=s.A_S
            @test 0.0 <= s.basin_parallel <= 1.0
            @test isfinite(s.A_S)
        end
        @test isfinite(baseline.A_S)
    end

    @testset "greedy_forward terminates and A_S non-decreasing" begin
        gf = greedy_forward(; tier=:smoke, topology=:modular_cells,
                              seed=3, visibility=1.0f0,
                              max_size=3, gain_threshold=0.02)
        @info "greedy result" selected=gf.selected n_history=length(gf.history)
        @test length(gf.selected) <= 3
        # Anchor is at history[1] (step 0). Each subsequent step has a
        # NamedTuple with the same fields.
        @test gf.history[1].step == 0
        @test all(haskey(r, :A_S) for r in gf.history)
        # A_S is non-decreasing along the history (we only accept gains
        # >= 0.02, so history is monotone within ±1e-6).
        as_vals = [r.A_S for r in gf.history]
        for i in 2:length(as_vals)
            @test as_vals[i] >= as_vals[i-1] - 1f-6
        end
    end
end
