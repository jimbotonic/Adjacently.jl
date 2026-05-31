# Tests for item 4 of research/mycelial_polis/RESEARCH_TASKS.md:
# threshold + complex-contagion adoption rules.
#
# Acceptance:
#   - on a clique with multiple seeds, threshold and complex converge to the
#     same committed-count fixed point;
#   - on a path graph with one seed, threshold spreads to most nodes while
#     complex stalls; the Centola–Macy ratio bound is
#         committed(complex) < 0.5 × committed(threshold) on the path;
#   - on the :modular_cells topology with the host disabled, both rules
#     produce monotone non-decreasing committed counts (no role back-slide
#     in item 4's scope — defection / infiltration are item 6).

using Test
using Random
using LightGraphs: SimpleDiGraph, nv, ne, add_edge!
using Adjacently
using Adjacently.MycelialPolis

# ---- test-only world builders ----------------------------------------------
# Built directly (not via `build_topology`) so we can construct the exact
# graph shapes the Centola–Macy gap requires.

function _make_world_from_graph(g_s, rule::Symbol;
                                seed::Int=42, overrides::NamedTuple=NamedTuple())
    n = Int(nv(g_s))
    layers = Dict(:S => g_s,
                  :C => SimpleDiGraph{eltype(g_s)}(eltype(g_s)(n)),
                  :E => SimpleDiGraph{eltype(g_s)}(eltype(g_s)(n)),
                  :T => SimpleDiGraph{eltype(g_s)}(eltype(g_s)(n)),
                  :G => SimpleDiGraph{eltype(g_s)}(eltype(g_s)(n)),
                  :O => SimpleDiGraph{eltype(g_s)}(eltype(g_s)(n)))
    mp = Multiplex(layers)
    agents = [Agent(id=i) for i in 1:n]
    infra = Dict{Symbol,Set{Int}}()
    infra_min = Dict{Symbol,Int}()
    params = merge(default_params(), (adoption_rule = rule,), overrides)
    return World(agents, mp, infra, infra_min, MersenneTwister(seed), params, 0)
end

function _make_path_world(n::Int, rule::Symbol; overrides=NamedTuple())
    g = SimpleDiGraph{UInt32}(UInt32(n))
    for i in 1:(n-1)
        add_edge!(g, UInt32(i),   UInt32(i+1))
        add_edge!(g, UInt32(i+1), UInt32(i))   # bidirectional path
    end
    return _make_world_from_graph(g, rule; overrides=overrides)
end

function _make_clique_world(n::Int, rule::Symbol; overrides=NamedTuple())
    g = SimpleDiGraph{UInt32}(UInt32(n))
    for i in 1:n, j in 1:n
        i == j && continue
        add_edge!(g, UInt32(i), UInt32(j))
    end
    return _make_world_from_graph(g, rule; overrides=overrides)
end

function _seed!(world, ids::Vector{Int}; role::Symbol=:contributor, commitment=0.8f0)
    for id in ids
        a = world.agents[id]
        a.role       = role
        a.awareness  = 1f0
        a.commitment = commitment
        a.identity   = 0.8f0
    end
    return world
end

# ---- 1. clique convergence (both rules → same fixed point) ------------------

@testset "clique convergence: threshold ≡ complex" begin
    n = 12
    seeds = [1, 2, 3]   # 3 of 12 — opens the k_cc=2 gate for every non-seed

    # α=4 boosts the drive so the per-stage thresholds can actually be
    # cleared from a single-step exposure pattern on a 12-clique. Both
    # rules see the same boosted drive; only the gate differs.
    overrides = (alpha = 4f0,)

    wt = _make_clique_world(n, :threshold;        overrides=overrides); _seed!(wt, seeds)
    wc = _make_clique_world(n, :complex_contagion; overrides=overrides); _seed!(wc, seeds)
    host = RandomHost(budget=0)

    for _ in 1:20
        step!(wt, host)
        step!(wc, host)
    end

    ct = committed_count(wt)
    cc = committed_count(wc)
    @info "clique result" n=n seeds=length(seeds) threshold_committed=ct complex_committed=cc
    @test ct == cc                          # exact same fixed point
    @test ct == n                           # full convergence on a clique
end

# ---- 2. path graph: threshold spreads, complex stalls ----------------------

@testset "path-graph Centola–Macy gap" begin
    n = 20
    seeds = [1]   # one seed at the end → complex with k_cc=2 has no path to spread

    # α=4 again — without it, a single-committed-neighbour signal can't
    # clear the :user stage threshold (0.30) on a path (drive = α·1/4).
    # With α=4 the threshold rule pushes a one-node-per-step wave from
    # the seed; complex still stalls because every step the wave-front
    # has exactly one committed in-neighbour, not k_cc=2.
    overrides = (alpha = 4f0,)

    wt = _make_path_world(n, :threshold;         overrides=overrides); _seed!(wt, seeds)
    wc = _make_path_world(n, :complex_contagion; overrides=overrides); _seed!(wc, seeds)
    host = RandomHost(budget=0)

    for _ in 1:30
        step!(wt, host)
        step!(wc, host)
    end

    ct = committed_count(wt)
    cc = committed_count(wc)
    @info "path result" n=n seeds=length(seeds) threshold_committed=ct complex_committed=cc
    @test ct >= div(3n, 4)                  # threshold reaches ≥ 75% of the path
    @test cc <= length(seeds)               # complex did not spread past the seed(s)
    @test cc < 0.5 * ct                     # Centola–Macy gap (the spec-named bound)
end

# ---- 3. modular_cells: both monotone non-decreasing when host disabled ------

@testset "modular_cells monotone adoption (host disabled)" begin
    for rule in (:threshold, :complex_contagion)
        params = merge(default_params(),
                       (adoption_rule = rule, initial_committed_frac = 0.05))
        world = build_world(; topology=:modular_cells, n=200, seed=7, params=params)
        # Start the complex rule with a slightly larger seed fraction so its
        # gate has somewhere to open; the test is about MONOTONICITY, not
        # spread rate, so this is fair.
        if rule === :complex_contagion
            params2 = merge(default_params(),
                            (adoption_rule = rule, initial_committed_frac = 0.15))
            world = build_world(; topology=:modular_cells, n=200, seed=7, params=params2)
        end
        host = RandomHost(budget=0)
        history = Int[committed_count(world)]
        for _ in 1:30
            step!(world, host)
            push!(history, committed_count(world))
        end
        # Monotone non-decreasing committed count over the whole trajectory.
        violations = sum(1 for i in 2:length(history) if history[i] < history[i-1]; init=0)
        @info "monotone check" rule=rule first=history[1] last=history[end] violations=violations
        @test violations == 0
        @test history[end] >= history[1]    # at least no net regression
    end
end
