# Tests for item 5 of research/mycelial_polis/RESEARCH_TASKS.md:
# seven host strategies + the Experiment 2 (§23) qualitative reproduction.
#
# Acceptance:
#   - all seven HostStrategy subtypes run end-to-end through step!(world, host)
#     without throwing on a small world (smoke);
#   - on the federated_hubs topology, by-degree attack fragments the
#     giant component substantially more than uniform-random attack
#     (Albert–Jeong–Barabási 2000);
#   - on the modular_cells topology, by-degree and random produce roughly
#     comparable fragmentation (no scale-free hubs to exploit).
#
# All host views are set to :S (the underlying social-trust graph) so the
# test isolates *topology* robustness from G_O leakage effects. The
# realistic :O view is the default in `host.jl` and is exercised by the
# CLI driver, not by this test.

using Test
using Random
using LightGraphs: nv, outneighbors, inneighbors
using Adjacently
using Adjacently.MycelialPolis

# --- helper: largest weakly-connected component on G_S | active -------------
function _giant_active(world)
    g = layer(world.multiplex, :S)
    n = Int(nv(g))
    active = falses(n)
    for a in world.agents
        if a.role !== :removed && a.role !== :outsider && a.id <= n
            active[a.id] = true
        end
    end
    sum(active) == 0 && return 0
    visited = falses(n)
    best = 0; q = Int[]
    T = eltype(g)
    for v in 1:n
        (active[v] && !visited[v]) || continue
        empty!(q); push!(q, v); visited[v] = true
        size = 0
        while !isempty(q)
            u = pop!(q); size += 1
            for w in outneighbors(g, T(u))
                iw = Int(w); iw <= n && active[iw] && !visited[iw] &&
                    (visited[iw] = true; push!(q, iw))
            end
            for w in inneighbors(g, T(u))
                iw = Int(w); iw <= n && active[iw] && !visited[iw] &&
                    (visited[iw] = true; push!(q, iw))
            end
        end
        size > best && (best = size)
    end
    return best
end

@testset "MycelialPolis host strategies (item 5)" begin

    @testset "all seven strategies run on a small world" begin
        for (name, host) in (
                ("RandomHost",            RandomHost(budget=2)),
                ("DegreeHost",            DegreeHost(budget=2, view=:S)),
                ("BetweennessHost",       BetweennessHost(budget=2, view=:S, recompute_every=3)),
                ("LegibilityHost",        LegibilityHost(budget=2)),
                ("LocalizedHost",         LocalizedHost(budget=2, region_size=6)),
                ("InfiltrationFirstHost", InfiltrationFirstHost(budget=4, plant_phase=3, escalate_phase=3)),
                ("AttritionHost",         AttritionHost(budget=2)),
            )
            world = build_world(; topology=:modular_cells, n=100, seed=11,
                                params=merge(default_params(),
                                             (initial_committed_frac=0.15,)))
            local ok = true
            local n_attacked_total = 0
            try
                for t in 1:8
                    snap = step!(world, host)
                    n_attacked_total += snap.n_attacked
                end
            catch e
                @error "host strategy crashed" name=name exception=(e, catch_backtrace())
                ok = false
            end
            @info "host smoke" name=name attacked=n_attacked_total
            @test ok
        end
    end

    @testset "Experiment 2: scale-free fragile to degree; modular flatter" begin
        N = 400
        T = 25
        seed = 7
        base_params = merge(default_params(), (initial_committed_frac=0.15,))

        function run_cell(topology, host)
            world = build_world(; topology=topology, n=N, seed=seed, params=base_params)
            for _ in 1:T
                step!(world, host)
            end
            return _giant_active(world)
        end

        # JIT warm-up so the test's wall-clock doesn't penalize the first cell.
        _ = run_cell(:modular_cells, RandomHost(budget=1))

        fed_random = run_cell(:federated_hubs, RandomHost(budget=4))
        fed_degree = run_cell(:federated_hubs, DegreeHost(budget=4, view=:S))
        mod_random = run_cell(:modular_cells, RandomHost(budget=4))
        mod_degree = run_cell(:modular_cells, DegreeHost(budget=4, view=:S))

        @info "Experiment 2 final giant (active)" fed_random fed_degree mod_random mod_degree

        # 1. Federated hubs: by-degree fragments substantially more than random.
        #    Albert–Jeong–Barabási qualitative: scale-free networks are fragile
        #    to targeted hub removal.
        @test fed_degree < 0.7 * fed_random

        # 2. Modular cells: no scale-free hubs to exploit, so by-degree should
        #    not dramatically outperform random. Allow ±30 % spread.
        ratio_mod = mod_degree / max(mod_random, 1)
        @test 0.7 < ratio_mod < 1.3
    end

    @testset "infiltration plants :infiltrator role" begin
        # InfiltrationFirstHost should put at least some agents into the
        # :infiltrator role during its plant phase, before arrests start.
        world = build_world(; topology=:modular_cells, n=200, seed=42,
                            params=merge(default_params(),
                                         (initial_committed_frac=0.15,)))
        host = InfiltrationFirstHost(budget=6, plant_phase=4, escalate_phase=4)
        for t in 1:4   # just the plant phase
            step!(world, host)
        end
        n_infil = count(is_infiltrator, world.agents)
        @info "infiltration plant phase result" n_infiltrators=n_infil
        @test n_infil >= 2
    end

    @testset "legibility attack escalates L_ext before removing" begin
        world = build_world(; topology=:modular_cells, n=200, seed=42,
                            params=merge(default_params(),
                                         (initial_committed_frac=0.15,)))
        host = LegibilityHost(budget=3, escalate_first=true, threshold=0.9f0)
        # One step of escalation (threshold=0.9 means most targets get an
        # L_ext bump rather than an arrest in step 1).
        step!(world, host)
        # Some agents should now have L_ext bumped above their default.
        bumped = count(a -> a.L_ext > 0.5f0 && a.role !== :removed, world.agents)
        @info "legibility escalation result" n_bumped=bumped
        @test bumped >= 1
    end
end
