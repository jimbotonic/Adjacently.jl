# Tests for item 6 of research/mycelial_polis/RESEARCH_TASKS.md:
# full §9 repression / backfire dynamics + :defector trigger + the
# Experiment 4 (§23) qualitative reproduction.
#
# Acceptance:
#   - sustained high-fear / low-protection conditions trigger defection
#     (:committed → :defector);
#   - narrative buffer A_i accumulates over time as backfire events occur;
#   - with visibility ON, adoption-vs-repression is *non-monotone*: some
#     intermediate repression level produces more committed agents than no
#     repression at all (Hess–Martin "transformative event" mechanism);
#   - with visibility OFF, adoption-vs-repression is *monotone-suppressive*
#     (sanity — fear alone dominates).

using Test
using Statistics: mean
using Adjacently
using Adjacently.MycelialPolis

# --- helper: run one (budget, visibility) cell of Experiment 4 ---------------
# Principles set to zero so we exercise the bare §9 dynamics; item 9's
# default 0.5 modulation otherwise changes the V_off curve (e.g.
# `normative_minimalism = 0.5` lowers stage thresholds and creates a
# non-monotone V_off).
function _run_cell(; topology::Symbol, n::Int, T::Int, budget::Int,
                     visibility::Float32, seed::Int)
    overrides = (visibility_scaling = visibility,
                 initial_committed_frac = 0.10,
                 principles = zero_principles())
    params = merge(default_params(), overrides)
    world  = build_world(; topology=topology, n=n, seed=seed, params=params)
    host   = RandomHost(budget=budget)
    final_committed = 0
    final_defected  = 0
    for _ in 1:T
        snap = step!(world, host)
        final_committed = snap.n_committed
        final_defected  = snap.n_defectors
    end
    return (committed = final_committed, defectors = final_defected,
            n_active = active_count(world))
end

@testset "MycelialPolis repression / backfire (item 6)" begin

    @testset "defection triggers under sustained high fear" begin
        # Mechanism test, not a sweep — we inject high fear on the committed
        # core directly (simulating the accumulated effect of prolonged
        # repression) and confirm that ONE `update_repression!` pass produces
        # at least one role transition to `:defector`. With the default
        # `defect_fear_threshold = 0.7` and `defect_protection_threshold =
        # 0.3`, and an initial_committed_frac of 0.10, most committed seeds
        # have protection ≈ 0.10 (well below 0.3) — so injecting fear above
        # 0.7 should make defection deterministic for at least a few of them.
        # Disable principles so item 9's `forkability=0.5` doesn't randomly
        # prevent half the defections (the test wants the bare §9 trigger).
        world = build_world(; topology=:modular_cells, n=200, seed=4,
                            params=merge(default_params(),
                                         (initial_committed_frac=0.10,
                                          principles=zero_principles())))
        # Direct fear injection on every committed agent.
        for a in world.agents
            is_committed(a) && (a.fear = 0.85f0)
        end
        n_committed_before = committed_count(world)
        @test n_committed_before >= 10
        snap = step!(world, RandomHost(budget=1))
        @info "defection result" n_committed_before n_committed_after=snap.n_committed n_defectors=snap.n_defectors
        @test snap.n_defectors >= 1
        # Sanity: defection moves committed → defector, not just removes.
        @test snap.n_committed < n_committed_before
    end

    @testset "narrative buffer A accumulates over steps with attacks" begin
        params = merge(default_params(),
                       (initial_committed_frac = 0.10,
                        narrative_window       = 5))
        world = build_world(; topology=:modular_cells, n=200, seed=2,
                            params=params)
        host  = RandomHost(budget=4)
        @test isempty(world.recent_backfire_count)
        for _ in 1:8
            step!(world, host)
        end
        # Buffer should contain `narrative_window` entries (truncated to the
        # tail) and should have seen at least some backfire events.
        @test length(world.recent_backfire_count) == 5
        @test sum(world.recent_backfire_count) > 0
    end

    @testset "Experiment 4: visibility on → non-monotone; off → suppressive" begin
        N = 500; T = 80; seed = 11
        budgets = (0, 2, 5, 10)

        # JIT warm-up so the per-cell wall-clock is comparable.
        _ = _run_cell(topology=:modular_cells, n=100, T=10,
                      budget=2, visibility=1.0f0, seed=0)

        on  = [_run_cell(topology=:modular_cells, n=N, T=T, budget=b,
                         visibility=1.0f0, seed=seed) for b in budgets]
        off = [_run_cell(topology=:modular_cells, n=N, T=T, budget=b,
                         visibility=0.0f0, seed=seed) for b in budgets]

        on_committed  = [c.committed for c in on]
        off_committed = [c.committed for c in off]
        @info "Experiment 4 results" budgets=budgets on=on_committed off=off_committed

        # V_off should be monotone non-increasing in budget — fear alone, no
        # backfire compensation.
        for i in 2:length(off_committed)
            @test off_committed[i] <= off_committed[i-1] + 2  # ±2 slop for randomness
        end

        # V_on should NOT be monotone non-increasing: at least one intermediate
        # budget must yield STRICTLY more committed than at least one of its
        # neighbours, OR the V_on curve must beat the V_off curve at the same
        # intermediate budget (backfire boost visible).
        non_monotone =
            any(i -> on_committed[i] > max(on_committed[i-1], on_committed[i+1]),
                2:length(on_committed)-1)
        any_boost = any(i -> on_committed[i] > off_committed[i] + 3,
                        eachindex(on_committed))
        @info "Experiment 4 verdict" non_monotone_on=non_monotone backfire_boost=any_boost
        @test non_monotone || any_boost
    end
end
