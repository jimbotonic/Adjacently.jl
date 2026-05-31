# Tests for item 7 of research/mycelial_polis/RESEARCH_TASKS.md:
# infrastructure continuity metrics with per-metric regression checks on
# toy worlds with hand-computed expected values.
#
# Acceptance:
#   - every metric has a documented range and a regression test on a small
#     world with a known correct value;
#   - `Lambda_sat ∈ [0,1)` holds on every snapshot of a real run;
#   - per-step cost ≤ 50 ms at N=1000 (already verified by smoke test);
#   - replica latency: a freshly-promoted steward does NOT count toward
#     Ψ_T until its warm-up period elapses.

using Test
using LightGraphs: SimpleDiGraph, add_edge!, nv
using Adjacently
using Adjacently.MycelialPolis

# Build a tiny custom world directly (bypassing build_topology) so we can
# reason about expected metric values by hand.
function _toy_world(; n::Int, edges::Vector{Tuple{Int,Int}},
                     roles::Vector{Symbol}=fill(:outsider, n),
                     params::NamedTuple = default_params())
    g = SimpleDiGraph{UInt32}(UInt32(n))
    for (u, v) in edges
        add_edge!(g, UInt32(u), UInt32(v))
    end
    layers = Dict(
        :S => g,
        :C => SimpleDiGraph{UInt32}(UInt32(n)),
        :E => SimpleDiGraph{UInt32}(UInt32(n)),
        :T => SimpleDiGraph{UInt32}(UInt32(n)),
        :G => SimpleDiGraph{UInt32}(UInt32(n)),
        :O => SimpleDiGraph{UInt32}(UInt32(n)),
    )
    mp = Multiplex(layers)
    agents = [Agent(id=i) for i in 1:n]
    for (i, r) in enumerate(roles)
        agents[i].role = r
        if r === :contributor || r === :steward
            agents[i].commitment = 0.8f0
            agents[i].L_int = 0.6f0
            agents[i].L_ext = 0.2f0
        end
    end
    return World(agents, mp, Dict{Symbol,Set{Int}}(), Dict{Symbol,Int}(),
                 Adjacently.MycelialPolis.MersenneTwister(0), params, 0)
end

@testset "MycelialPolis metrics (item 7)" begin

    @testset "committed_frac" begin
        w = _toy_world(n=10, edges=Tuple{Int,Int}[],
                       roles=[:contributor, :contributor, :user,
                              :outsider, :outsider, :outsider,
                              :outsider, :outsider, :outsider, :sympathizer])
        # 3 committed (rank ≥ :user) out of 10.
        @test committed_frac(w) ≈ 0.30 atol=1e-6
        @test committed_count(w) == 3
    end

    @testset "giant_active_frac on a known path" begin
        # 6-node bidirectional path: 1-2-3-4-5-6
        edges = vcat([(i, i+1) for i in 1:5], [(i+1, i) for i in 1:5])
        roles = [:contributor for _ in 1:6]
        w = _toy_world(n=6, edges=edges, roles=roles)
        # All 6 active and connected → giant fraction = 6/6 = 1.0
        @test giant_active_frac(w) ≈ 1.0 atol=1e-6
        # Remove node 3 → splits into {1,2} and {4,5,6}; giant = 3/6 = 0.5
        w.agents[3].role = :removed
        @test giant_active_frac(w) ≈ 0.5 atol=1e-6
    end

    @testset "kcore_depth on a known clique-plus-tail" begin
        # 5-clique on nodes 1..5 + a tail 5↔6↔7. LightGraphs.core_number on
        # a directed graph counts *total* degree (in + out), so each
        # clique node has 4 in + 4 out = 8 → depth = 8. The tail nodes
        # have lower core numbers and are peeled first.
        edges = Tuple{Int,Int}[]
        for i in 1:5, j in 1:5
            i == j && continue
            push!(edges, (i, j))
        end
        push!(edges, (5, 6)); push!(edges, (6, 5))
        push!(edges, (6, 7)); push!(edges, (7, 6))
        w = _toy_world(n=7, edges=edges, roles=fill(:contributor, 7))
        d = kcore_depth(w)
        @info "k-core depth" depth=d
        @test d == 8
    end

    @testset "lambda + lambda_sat: range invariants" begin
        # Construct a deliberately high-Λ case: L_int = 0.9, L_ext = 0.01 → Λ ≈ 90.
        w = _toy_world(n=5, edges=Tuple{Int,Int}[], roles=fill(:contributor, 5))
        for a in w.agents
            a.L_int = 0.9f0
            a.L_ext = 0.01f0
        end
        Λ  = lambda(w)
        Λs = lambda_sat(w)
        @info "Λ check" Lambda=Λ Lambda_sat=Λs
        @test Λ > 10f0                # large by construction
        @test 0f0 <= Λs < 1f0         # saturated form is bounded
    end

    @testset "psi_T with explicit infra setup + replica latency" begin
        # Two functions: :comm needs 2 replicas, :store needs 1.
        w = _toy_world(n=4, edges=Tuple{Int,Int}[],
                       roles=[:contributor, :contributor, :user, :outsider])
        w.infra     = Dict(:comm => Set([1, 2]), :store => Set([3]))
        w.infra_min = Dict(:comm => 2,           :store => 1)
        # No latency — both functions are fully staffed → Ψ_T = 1.0.
        @test psi_T(w) ≈ 1.0 atol=1e-6
        @test infra_outages(w) == 0
        # Put agent 1 on a 2-step warm-up → :comm has only 1 effective
        # replica, below threshold 2 → Ψ_T = 0.5 (1/2 functions OK).
        w.infra_latency[:comm] = Dict(1 => 2)
        @test psi_T(w) ≈ 0.5 atol=1e-6
        @test infra_outages(w) == 1
    end

    @testset "phi range + reweighting" begin
        # A modest world — fully connected 8-node ring, all committed.
        edges = vcat([(i, mod1(i+1, 8)) for i in 1:8],
                     [(mod1(i+1, 8), i) for i in 1:8])
        w = _toy_world(n=8, edges=edges,
                       roles=fill(:contributor, 8),
                       params=default_params())
        p1 = phi(w)
        @info "phi default" p1
        @test 0f0 <= p1 <= 1f0
        # Reweight to "lambda_sat counts double, giant counts nothing" — the
        # output should change. Just check it stays in [0,1] and differs.
        wts = (core=0.20f0, giant=0.0f0, kcore=0.0f0,
               psi_T=0.0f0, lambda_sat=0.80f0)
        w.params = merge(w.params, (phi_weights=wts,))
        p2 = phi(w)
        @info "phi reweighted" p2
        @test 0f0 <= p2 <= 1f0
        @test p1 != p2
    end

    @testset "Lambda_sat ∈ [0,1) on a real-run snapshot trajectory" begin
        params = merge(default_params(), (initial_committed_frac=0.10,))
        w = build_world(; topology=:modular_cells, n=200, seed=42, params=params)
        host = RandomHost(budget=2)
        for _ in 1:25
            snap = step!(w, host)
            @test 0f0 <= snap.Lambda_sat < 1f0
        end
    end
end
