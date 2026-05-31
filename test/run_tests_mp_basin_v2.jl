# Tests for item 8 v2 of research/mycelial_polis/RESEARCH_TASKS.md:
# extended basin classifier with :transformative_replacement,
# :host_accommodation, :institutional_hybridization.
# :ai_irrelevance is still deferred (needs §17 AI variables).
#
# Acceptance:
#   - classify_v2 returns the expected v2 labels on hand-built synthetic
#     trajectories that exercise each predicate;
#   - v2 falls back to v1 categories cleanly when none of the new
#     predicates fire;
#   - estimate_basins with classifier=:v2 + budget_decay > 0 actually
#     produces some of the new equilibria (smoke regime).

using Test
using Adjacently
using Adjacently.MycelialPolis

@testset "MycelialPolis basin v2 (item 8 v2)" begin

    @testset "classify_v2 on synthetic trajectories" begin
        T = 40; N = 100
        # Extinction — Φ→0 at end.
        phi = vcat(fill(0.5f0, 20), fill(0.01f0, 20))
        comm = vcat(fill(20, 20), fill(0, 20))
        atk  = fill(2, 40)
        @test classify_v2(phi, comm, atk, N) === :extinction

        # Transformative replacement — committed >= 50% of N.
        phi  = fill(0.7f0, T);  comm = vcat(fill(20, 20), fill(60, 20))
        atk  = fill(2, T)
        @test classify_v2(phi, comm, atk, N) === :transformative_replacement

        # Host accommodation — attacks decay AND polis grows.
        phi  = vcat(range(0.4f0, 0.6f0; length=T)...)
        phi  = Float32[0.40 + 0.005*t for t in 1:T]
        comm = [round(Int, 15 + 0.5*t) for t in 1:T]     # grows 15 → 35
        atk  = vcat(fill(5, 20), fill(0, 20))            # first quarter 5/step, last quarter 0/step
        @test classify_v2(phi, comm, atk, N) === :host_accommodation

        # Institutional hybridization — attacks decay BUT polis flat.
        phi  = fill(0.45f0, T)
        comm = fill(20, T)
        atk  = vcat(fill(5, 20), fill(1, 20))
        @test classify_v2(phi, comm, atk, N) === :institutional_hybridization

        # Falls through to v1 :parallel_coexistence when nothing v2 fires.
        phi  = fill(0.60f0, T); comm = fill(25, T); atk = fill(2, T)
        @test classify_v2(phi, comm, atk, N) === :parallel_coexistence

        # Falls through to v1 :dormant_persistence at lower Φ.
        phi  = fill(0.20f0, T); comm = fill(8, T); atk = fill(2, T)
        @test classify_v2(phi, comm, atk, N) === :dormant_persistence
    end

    @testset "estimate_basins :v2 with budget_decay produces new equilibria" begin
        # With strong budget decay, the host's repression erodes over
        # time. At least some of the LHS samples should fall into the
        # v2 equilibria; the v1 estimator would only report the three
        # original categories.
        result_v1 = estimate_basins(; tier=:smoke, topology=:modular_cells,
                                      seed=42, classifier=:v1,
                                      budget_decay=0.05)
        result_v2 = estimate_basins(; tier=:smoke, topology=:modular_cells,
                                      seed=42, classifier=:v2,
                                      budget_decay=0.05)
        # v1 counts should sum to n_samples and contain only the v1 labels.
        @test sum(values(result_v1.counts)) == result_v1.n_samples
        @test haskey(result_v1.counts, :extinction)
        @test !haskey(result_v1.counts, :host_accommodation)

        # v2 counts should sum to n_samples and include the v2 labels.
        @test sum(values(result_v2.counts)) == result_v2.n_samples
        @test haskey(result_v2.counts, :host_accommodation)
        @test haskey(result_v2.counts, :institutional_hybridization)
        @test haskey(result_v2.counts, :transformative_replacement)

        # With budget_decay = 0.05 (significant per-step erosion), at
        # least ONE sample should land in a v2 category.
        v2_only = result_v2.counts[:host_accommodation] +
                  result_v2.counts[:institutional_hybridization] +
                  result_v2.counts[:transformative_replacement]
        @info "v2 estimate" v1=result_v1.counts v2=result_v2.counts v2_only=v2_only
        @test v2_only >= 1
    end

    @testset "endogenous accommodation_rate produces v2 equilibria" begin
        # accommodation_rate > 0 is the "AccommodatingHost" mechanism:
        # the host lowers its own budget when the polis isn't shrinking.
        # Distinct from the schedule-based budget_decay tested above.
        r_baseline = estimate_basins(; tier=:smoke, topology=:modular_cells,
                                       seed=42, classifier=:v2)
        r_accom    = estimate_basins(; tier=:smoke, topology=:modular_cells,
                                       seed=42, classifier=:v2,
                                       accommodation_rate=0.10)
        v2_only_base  = r_baseline.counts[:host_accommodation] +
                        r_baseline.counts[:institutional_hybridization]
        v2_only_accom = r_accom.counts[:host_accommodation] +
                        r_accom.counts[:institutional_hybridization]
        @info "accommodation_rate effect" baseline_v2=v2_only_base accom_v2=v2_only_accom
        @test v2_only_accom > v2_only_base    # accommodation enabled → more v2 cells

        # Sanity: with accommodation_rate = 0 but no decay, the host
        # never backs off → no accommodation/hybridization at all.
        @test r_baseline.counts[:host_accommodation] <= 5
        @test r_baseline.counts[:institutional_hybridization] <= 5
    end

    @testset "_accommodate! direct mechanism check" begin
        # Build a host with high accommodation_rate, manually feed it
        # decreasing then stable n_committed values, and confirm budget
        # only decays once the polis stops shrinking.
        host = RandomHost(budget=10, accommodation_rate=0.5)
        # Polis shrinking → no accommodation, budget intact.
        Adjacently.MycelialPolis._accommodate!(host, 100)
        Adjacently.MycelialPolis._accommodate!(host, 90)
        Adjacently.MycelialPolis._accommodate!(host, 80)
        @test host.budget == 10
        # Polis stops shrinking (3 stable steps with rate=0.5 should erode budget).
        for _ in 1:5
            Adjacently.MycelialPolis._accommodate!(host, 80)
        end
        @info "accommodation budget trace" final_budget=host.budget
        @test host.budget < 10
    end

    @testset "AI capability asymmetry (A_L, A_H) — :ai_irrelevance classifier" begin
        # Synthetic trajectories with A_L=A_H=0.5 should classify as
        # :ai_irrelevance (both sides AI-mediated).
        T = 40; N = 100
        phi  = fill(0.5f0, T); comm = fill(25, T); atk = fill(2, T)
        @test classify_v2(phi, comm, atk, N; A_L=0.6f0, A_H=0.6f0) === :ai_irrelevance
        # A_L low, A_H high — not :ai_irrelevance; fall through to v2/v1 rules
        @test classify_v2(phi, comm, atk, N; A_L=0.0f0, A_H=0.9f0) === :parallel_coexistence
        # Both above 0.5 but polis extinct — still :extinction (checked first)
        phi_dead = fill(0.01f0, T); comm_dead = fill(0, T)
        @test classify_v2(phi_dead, comm_dead, atk, N; A_L=0.7f0, A_H=0.7f0) === :extinction
    end

    @testset "A_L boosts adoption, A_H boosts surveillance" begin
        # Run two cells: one with A_L=0.8, one with A_H=0.8. The A_L cell
        # should produce more committed agents on a moderately-hostile
        # regime; the A_H cell should produce higher mean L_ext.
        params_AL = merge(default_params(), (initial_committed_frac=0.10, A_L=0.8f0))
        params_AH = merge(default_params(), (initial_committed_frac=0.10, A_H=0.8f0))
        for (label, p) in (("A_L", params_AL), ("A_H", params_AH))
            world = build_world(; topology=:modular_cells, n=200, seed=3, params=p)
            host  = RandomHost(budget=2)
            for _ in 1:30
                step!(world, host)
            end
            @info "AI hook check" label committed=committed_count(world) mean_L_ext=Adjacently.MycelialPolis.mean_field(world, a -> a.L_ext)
        end
        # A_L should produce more committed than A_H baseline (both with budget=2)
        w_L = build_world(; topology=:modular_cells, n=200, seed=3, params=params_AL)
        w_H = build_world(; topology=:modular_cells, n=200, seed=3, params=params_AH)
        h = RandomHost(budget=2)
        h2 = RandomHost(budget=2)
        for _ in 1:30
            step!(w_L, h); step!(w_H, h2)
        end
        @test committed_count(w_L) > committed_count(w_H)   # A_L wins recruitment
    end

    @testset "estimate_basins :v2 with A_L=A_H=0.7 yields :ai_irrelevance" begin
        r = estimate_basins(; tier=:smoke, topology=:modular_cells, seed=0,
                              classifier=:v2,
                              A_L=0.7f0, A_H=0.7f0)
        @info "AI basin" counts=r.counts
        # Most non-extinct samples should classify as :ai_irrelevance
        # (since A_L=A_H=0.7 ≥ 0.5 threshold).
        @test haskey(r.counts, :ai_irrelevance)
        @test r.counts[:ai_irrelevance] > 0
    end

    @testset "v2 with budget_decay=0 reduces to v1 categories" begin
        # Without decay, host attacks at constant rate — no accommodation,
        # no hybridization possible. transformative_replacement is still
        # achievable if polis grows past 50%.
        result = estimate_basins(; tier=:smoke, topology=:modular_cells,
                                   seed=0, classifier=:v2, budget_decay=0.0)
        # accommodation + hybridization should both be near-zero.
        n_acc = result.counts[:host_accommodation]
        n_hyb = result.counts[:institutional_hybridization]
        @info "v2 no-decay counts" host_accommodation=n_acc institutional_hybridization=n_hyb
        @test n_acc <= 5    # essentially zero (allow tiny edge effects)
        @test n_hyb <= 5
    end
end
