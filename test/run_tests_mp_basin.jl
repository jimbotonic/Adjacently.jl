# Tests for item 8 v1 of research/mycelial_polis/RESEARCH_TASKS.md:
# basin-of-attraction estimator (3 equilibria).
#
# Acceptance:
#   - classify() returns the expected label on hand-built synthetic histories
#     for each of the four labels (:extinction, :dormant_persistence,
#     :parallel_coexistence, :unclassified);
#   - time_to_attractor() detects a stable plateau quickly and reports T
#     when the trajectory never stabilises;
#   - LHS sampler has correct stratification + reproducibility per seed;
#   - estimate_basins() at :smoke tier runs end-to-end, total counts equal
#     n_samples, basin TSV is written with the expected columns;
#   - reproducibility: same `seed` → identical counts.

using Test
using Random: MersenneTwister
using Statistics: mean, std
using Adjacently
using Adjacently.MycelialPolis

@testset "MycelialPolis basin estimator (item 8 v1)" begin

    @testset "classify on synthetic histories" begin
        # All-zero trajectory → :extinction (final Φ < EPS_EXT)
        @test classify(zeros(Float32, 40), 0) === :extinction
        # n_committed_final == 0 forces :extinction even if Φ tail is non-zero
        @test classify(fill(0.6f0, 40), 0)   === :extinction
        # Stable low-Φ trajectory → :dormant_persistence
        @test classify(fill(0.20f0, 40), 5)  === :dormant_persistence
        # Stable high-Φ trajectory → :parallel_coexistence
        @test classify(fill(0.70f0, 40), 50) === :parallel_coexistence
        # Big oscillation in the tail → :unclassified
        h = Float32[0.3 + 0.5 * sin(t * 0.7) for t in 1:40]
        # Make sure tail mean > EPS_EXT (else it'd be :extinction)
        h[end] = 0.4f0
        @test classify(h, 30) === :unclassified
    end

    @testset "time_to_attractor" begin
        # Quickly stable: starts noisy, settles. Window is max(5, T÷10) = 5.
        h = Float32[1.0; 0.5; 0.3; 0.2; 0.1; fill(0.1f0, 35)...]
        t = time_to_attractor(h)
        @info "T_attr on quickly-stable" t=t
        @test t <= 15
        # Never stable — large oscillation; T_attr == T.
        h2 = Float32[0.2 + 0.6 * sin(i * 0.9) for i in 1:40]
        @test time_to_attractor(h2) == 40
    end

    @testset "lhs sampler stratification + reproducibility" begin
        # Stratification: each column should have exactly one sample in each
        # of the n_samples [k/n, (k+1)/n) bins. By construction of LHS.
        rng = MersenneTwister(7)
        n = 50; d = 5
        M = lhs(n, d; rng=rng)
        @test size(M) == (n, d)
        @test all(0.0 .<= M .< 1.0)
        for col in 1:d
            bins = Int.(floor.(M[:, col] .* n)) .+ 1
            @test sort(bins) == 1:n     # one per stratum
        end
        # Reproducibility: same seed → identical matrix.
        rng_a = MersenneTwister(99); rng_b = MersenneTwister(99)
        @test lhs(20, 3; rng=rng_a) == lhs(20, 3; rng=rng_b)
    end

    @testset "smoke-tier estimate_basins end-to-end" begin
        # Smoke tier params per BASIN_TIERS[:smoke] but we'll just run
        # estimate_basins and check the invariants.
        result = estimate_basins(; tier=:smoke, topology=:modular_cells,
                                   seed=42, visibility=1.0f0)
        @info "smoke basin counts" counts=result.counts elapsed=result.elapsed_s

        total = sum(values(result.counts))
        @test total == result.n_samples == BASIN_TIERS[:smoke].n_samples
        # Across the 5-d LHS space, at minimum we should see SOME extinction
        # AND some parallel_coexistence — empirically at smoke tier the
        # distribution is dominated by extinction (high host_budget arm) and
        # parallel coexistence (low-budget arm).
        @test result.counts[:extinction] > 0
        @test result.counts[:parallel_coexistence] > 0
        # Wall budget: smoke must finish in ≤ 30 s.
        @test result.elapsed_s < 30.0
    end

    @testset "reproducibility: same seed → identical counts" begin
        a = estimate_basins(; tier=:smoke, topology=:modular_cells,
                              seed=1234, visibility=1.0f0)
        b = estimate_basins(; tier=:smoke, topology=:modular_cells,
                              seed=1234, visibility=1.0f0)
        @test a.counts == b.counts
    end

    @testset "TSV writer emits expected columns" begin
        tmp_dir = mktempdir()
        result = estimate_basins(; tier=:smoke, topology=:modular_cells,
                                   seed=5, visibility=1.0f0, out_dir=tmp_dir)
        tsv_path = joinpath(tmp_dir, "basins.tsv")
        @test isfile(tsv_path)
        # Header line + 4 equilibrium rows.
        lines = readlines(tsv_path)
        @test length(lines) == 5
        header = split(lines[1], "\t")
        @test header == ["equilibrium", "count", "fraction",
                         "mean_T_to_attractor", "mean_final_phi"]
    end
end
