# Tests for item 3 of research/mycelial_polis/RESEARCH_TASKS.md:
# three baseline topologies + layer projections + summary printer.
#
# Acceptance:
#   - each topology builds for n ∈ {200, 1000, 5000} in < 2 s
#   - topology summaries match expected qualitative fingerprints:
#     * modular: high modularity on the natural cell partition
#     * federated: high-skew degree (hubs have much higher degree than leaves)
#     * p2p:      moderate clustering (small-world regime); low diameter
#   - layer projections respect the spec inclusions:
#     G_T ⊆ G_C ⊆ G_S, G_E ⊆ G_S, G_O is non-empty under default p_leak

using Test
using Statistics: mean, std
using Random: MersenneTwister
using LightGraphs: nv, ne, edges, src, dst, outneighbors, inneighbors,
                   outdegree, indegree
using Adjacently
using Adjacently.MycelialPolis

# Defined here (before any @testset) so the testsets can see it.
function median_safe(v::AbstractVector{<:Real})
    isempty(v) && return 0.0
    s = sort(v); n = length(s)
    isodd(n) ? Float64(s[(n+1) ÷ 2]) : Float64(s[n ÷ 2] + s[n ÷ 2 + 1]) / 2
end

@testset "MycelialPolis topologies (item 3)" begin

    @testset "all three build for n ∈ {200, 1000, 5000} in < 2 s each" begin
        # JIT warm-up — without this, the first topology build for the first n
        # absorbs the compile cost of `build_topology` + the layer derivations
        # + the cluster builders. On the bench that's ~1.5–2 s, which can
        # bump the smallest case over the budget on a cold run.
        for topology in (:modular_cells, :federated_hubs, :p2p_mesh)
            build_topology(topology, 64; seed=0)
        end

        for topology in (:modular_cells, :federated_hubs, :p2p_mesh)
            for n in (200, 1000, 5000)
                t0 = time()
                mp = build_topology(topology, n; seed=42)
                elapsed = time() - t0
                @info "topology build" topology=topology n=n elapsed_s=round(elapsed; digits=2)
                @test elapsed < 2.0
                @test length(mp.layers) == 6
                for sym in (:S, :C, :E, :T, :G, :O)
                    @test haskey(mp.layers, sym)
                    @test nv(mp.layers[sym]) == n
                end
            end
        end
    end

    @testset "layer projection inclusions and shapes" begin
        mp = build_topology(:modular_cells, 500; seed=42,
                            keep_C=0.7, keep_E=0.5, keep_T=0.5, p_leak=0.10)
        gS = mp.layers[:S]; gC = mp.layers[:C]; gE = mp.layers[:E]
        gT = mp.layers[:T]; gG = mp.layers[:G]; gO = mp.layers[:O]

        # G_C ⊆ G_S (edge-thinned).
        edges_S = Set((Int(src(e)), Int(dst(e))) for e in edges(gS))
        edges_C = Set((Int(src(e)), Int(dst(e))) for e in edges(gC))
        edges_E = Set((Int(src(e)), Int(dst(e))) for e in edges(gE))
        edges_T = Set((Int(src(e)), Int(dst(e))) for e in edges(gT))
        @test issubset(edges_C, edges_S)
        @test issubset(edges_E, edges_S)
        @test issubset(edges_T, edges_C)

        # Thinning rates are close to nominal (0.7 / 0.5 / 0.5 of upstream).
        @test 0.55 < ne(gC) / max(ne(gS), 1) < 0.85
        @test 0.35 < ne(gE) / max(ne(gS), 1) < 0.65
        @test 0.35 < ne(gT) / max(ne(gC), 1) < 0.65

        # G_G keeps ≤ governance_k out-neighbours per vertex.
        @test all(outdegree(gG, v) <= 3 for v in 1:nv(gG))

        # G_O is non-empty under default p_leak (probabilistic; very unlikely empty).
        @test ne(gO) > 0
    end

    @testset "fingerprints: modular, federated, p2p" begin
        n = 800

        # --- modular_cells: modularity on natural partition should be high.
        mp_mod = build_topology(:modular_cells, n; seed=42)
        part   = natural_partition(:modular_cells, n)
        rows_m = topology_summary(mp_mod; cell_partition=part)
        sm_S   = rows_m[findfirst(r -> r.layer === :S, rows_m)]
        @info "modular_cells fingerprint" modularity=sm_S.modularity
        @test sm_S.modularity > 0.5     # high-mod by construction

        # --- federated_hubs: degree distribution should be heavy-tailed
        # (hubs hold most connections). Use the ratio of max degree to median
        # degree as a rough skew proxy.
        mp_fed = build_topology(:federated_hubs, n; seed=42)
        gS     = mp_fed.layers[:S]
        degs   = Float64[outdegree(gS, v) + indegree(gS, v) for v in 1:nv(gS)]
        skew   = maximum(degs) / max(median_safe(degs), 1.0)
        @info "federated_hubs fingerprint" max_deg=Int(maximum(degs)) median_deg=median_safe(degs) skew=skew
        @test skew > 5.0                # hubs vastly outweigh leaves

        # --- p2p_mesh: small-world regime — moderate clustering and small diameter.
        mp_p2p = build_topology(:p2p_mesh, n; seed=42)
        rows_p = topology_summary(mp_p2p)
        sp_S   = rows_p[findfirst(r -> r.layer === :S, rows_p)]
        @info "p2p_mesh fingerprint" clustering=sp_S.clustering diameter=sp_S.diameter
        @test sp_S.clustering > 0.05    # WS β=0.05 retains some lattice triangles
        @test sp_S.clustering < 0.6     # ... but is not a pure lattice
        @test sp_S.diameter < 50        # small-world: diameter ≪ n
    end
end


# Regression test for thin_to! (TIER1_EXPERIMENT_PLAN.md §E1, E1b matched-3
# low anchor): uniform edge deletion down to a target mean out-degree.
@testset "thin_to! (E1b matched-3 anchor)" begin
    n = 500
    edge_set(g) = Set((Int(src(e)), Int(dst(e))) for e in edges(g))

    # modular_cells native mean out-degree ≈ 5.94 → thin to 3.
    g = build_topology(:modular_cells, n; seed=42).layers[:S]
    es_before = edge_set(g)
    m_before = ne(g) / nv(g)
    @test m_before > 3.0          # thinning actually has work to do
    thin_to!(g, 3, MersenneTwister(42 + 13))
    @test nv(g) == n              # node count preserved
    @test abs(ne(g) / nv(g) - 3.0) <= 0.05   # mean out-degree on target
    @test edge_set(g) ⊆ es_before # thinning only removes edges

    # No-op at or below target (federated_hubs native ≈ 2.97 < 3).
    g2 = build_topology(:federated_hubs, n; seed=42).layers[:S]
    ne2 = ne(g2)
    thin_to!(g2, 3, MersenneTwister(42 + 13))
    @test ne(g2) == ne2

    # Reproducible per seed (same rng discipline as densify_to!).
    g3 = build_topology(:modular_cells, n; seed=42).layers[:S]
    g4 = build_topology(:modular_cells, n; seed=42).layers[:S]
    thin_to!(g3, 3, MersenneTwister(42 + 13))
    thin_to!(g4, 3, MersenneTwister(42 + 13))
    @test edge_set(g3) == edge_set(g4)
end
