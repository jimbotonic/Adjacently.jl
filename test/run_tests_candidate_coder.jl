# Round-trip tests for the candidate-bit coder (src/compression/candidate_coder.jl).
using Test, Random
using Adjacently
const C = Adjacently.Compression
isdefined(C, :cand_encode) || Base.include(C, joinpath(@__DIR__, "..", "src", "compression", "candidate_coder.jl"))

function randrows(rng, n, p; loops=true, emptyfrac=0.0)
    rows = Vector{Vector{UInt32}}(undef, n)
    for v in 1:n
        if rand(rng) < emptyfrac; rows[v] = UInt32[]; continue; end
        rows[v] = UInt32[x for x in 1:n if rand(rng) < p && (loops || x != v)]
    end
    rows
end
# clustered graph with reciprocity and triangles, so that every generator fires
function clustered(rng, n, k, deg)
    rows = [Set{UInt32}() for _ in 1:n]
    for v in 1:n, _ in 1:deg
        x = rand(rng) < 0.8 ? clamp(v + rand(rng, -k:k), 1, n) : rand(rng, 1:n)
        push!(rows[v], UInt32(x)); rand(rng) < 0.3 && push!(rows[x], UInt32(v))
    end
    [sort!(collect(s)) for s in rows]
end
const FROZEN_DEFAULT = (14044, 0xddea2d643410987c)   # (bytes, FNV-1a) of the default encode of the fixed graph below
roundtrip(rows; params=C.CandParams()) = (b = C.cand_encode(rows; params=params); (C.cand_decode(b, length(rows); params=params) == rows, b))

@testset "candidate coder" begin
    rng = MersenneTwister(20261120)
    @testset "parameters" begin
        p = C.CandParams(ref_window=7, gen_innbr=true, remainder=:cell, mix_lr=11, seed=false)
        @test C.cand_params_from_bytes(C.cand_params_bytes(p)) == p
        @test C.cand_params_from_bytes(C.cand_params_bytes(C.CandParams())) == C.CandParams()
        b = C.cand_params_bytes(p); b[1] = 0x7f
        @test_throws ArgumentError C.cand_params_from_bytes(b)
        @test_throws ArgumentError C.cand_params_from_bytes(C.cand_params_bytes(p)[1:end-1])
        @test_throws ArgumentError C.cand_encode([UInt32[]]; params=C.CandParams(remainder=:other))
    end
    variants = [C.CandParams(), C.CandParams(remainder=:cell), C.CandParams(ref_window=0), C.CandParams(seed=false, gen_innbr=true),
                C.CandParams(gen_recip=false, gen_cocite=false, gen_2hop=false), C.CandParams(cap_a=1, cap_in=1, cap_row=2, budget_base=0, budget_per_arc=1),
                C.CandParams(ref_window=3, ref_gain=100, ref_bias=0, min_score=2, remainder=:cell, cell_ds=1, cell_kappa=0)]
    @testset "edge cases" begin
        for p in variants
            @test C.cand_decode(C.cand_encode(Vector{UInt32}[]; params=p), 0; params=p) == Vector{UInt32}[]
            for rows in ([UInt32[]], [UInt32[1]], [UInt32[], UInt32[]], [UInt32[1, 2], UInt32[1, 2]], [UInt32[2], UInt32[1]],
                         [collect(UInt32, 1:40) for _ in 1:40], [UInt32[] for _ in 1:50])
                @test roundtrip(rows; params=p)[1]
            end
        end
        @test_throws ArgumentError C.cand_encode([UInt32[2, 1], UInt32[]])
        @test_throws ArgumentError C.cand_encode([UInt32[1, 1], UInt32[]])
        @test_throws ArgumentError C.cand_encode([UInt32[3], UInt32[]])
    end
    @testset "random graphs" begin
        for p in variants, (n, q, ef) in ((1, 0.5, 0.0), (2, 0.5, 0.0), (17, 0.3, 0.3), (60, 0.9, 0.0), (200, 0.03, 0.5), (300, 0.1, 0.1))
            @test roundtrip(randrows(rng, n, q; emptyfrac=ef); params=p)[1]
        end
        for p in variants, (n, k, dg) in ((400, 5, 4), (1500, 20, 8), (3000, 3, 12))
            @test roundtrip(clustered(rng, n, k, dg); params=p)[1]
        end
    end
    @testset "determinism and mismatch guard" begin
        rows = clustered(rng, 2000, 10, 8)
        ok, b = roundtrip(rows)
        @test ok
        @test C.cand_encode(rows) == b
        # the coder is useful on this graph: well below the trivial log2(n) bits per arc
        @test 8 * length(b) / sum(length, rows) < 9.0
        # a different graph must not decode to `rows`, and wrong parameters must not silently give `rows`
        rows2 = deepcopy(rows); i = findfirst(r -> length(r) > 1, rows2); pop!(rows2[i])
        @test C.cand_decode(C.cand_encode(rows2), length(rows2)) == rows2
        @test C.cand_decode(C.cand_encode(rows2), length(rows2)) != rows
        wrong = try; C.cand_decode(b, length(rows); params=C.CandParams(ref_window=0, mix_lr=9)); catch; nothing; end
        @test wrong != rows
        bad = copy(b); bad[length(bad) ÷ 2] ⊻= 0x55
        dec = try; C.cand_decode(bad, length(rows)); catch; nothing; end
        @test dec != rows
    end
    @testset "popgap remainder" begin
        pv = [C.CandParams(remainder=:popgap), C.CandParams(remainder=:popgap, pop_split=0), C.CandParams(remainder=:popgap, pop_split=3, gap_deg_buckets=4),
              C.CandParams(remainder=:popgap, pop_prior=64, gap_deg_buckets=16, ref_window=0),
              C.CandParams(remainder=:popgap, gen_recip=false, gen_cocite=false, gen_2hop=false, ref_window=0, pop_split=0),
              C.CandParams(gap_deg_buckets=8), C.CandParams(remainder=:cell, gap_deg_buckets=3)]
        for p in pv
            @test C.cand_params_from_bytes(C.cand_params_bytes(p)) == p
            @test C.cand_decode(C.cand_encode(Vector{UInt32}[]; params=p), 0; params=p) == Vector{UInt32}[]
            for rows in ([UInt32[]], [UInt32[1]], [UInt32[2], UInt32[1]], [collect(UInt32, 1:40) for _ in 1:40], [UInt32[] for _ in 1:50],
                         [UInt32[300], UInt32[]], vcat([UInt32[1, 500]], [UInt32[] for _ in 1:499]))
                n = max(length(rows), maximum(r -> isempty(r) ? 0 : Int(maximum(r)), rows))
                full = vcat(rows, [UInt32[] for _ in 1:n-length(rows)])
                @test roundtrip(full; params=p)[1]
            end
            for (n, q, ef) in ((1, 0.5, 0.0), (17, 0.3, 0.3), (60, 0.9, 0.0), (700, 0.01, 0.4), (300, 0.1, 0.1))
                @test roundtrip(randrows(rng, n, q; emptyfrac=ef); params=p)[1]
            end
            for (n, k, dg) in ((400, 5, 4), (5000, 400, 6), (3000, 3, 12))
                @test roundtrip(clustered(rng, n, k, dg); params=p)[1]
            end
        end
        # hub-heavy graph: popularity weighting must not lose to plain gaps by more than noise, and must round-trip
        hub = [sort!(unique!(UInt32[rand(rng) < 0.7 ? rand(rng, 1:40) * 97 : rand(rng, 1:4000) for _ in 1:6])) for v in 1:4000]
        okp, bp = roundtrip(hub; params=C.CandParams(remainder=:popgap, ref_window=0, gen_recip=false, gen_cocite=false, gen_2hop=false))
        okg, bg = roundtrip(hub; params=C.CandParams(ref_window=0, gen_recip=false, gen_cocite=false, gen_2hop=false))
        @test okp && okg
        @test length(bp) < length(bg)
        @test_throws ArgumentError C.cand_encode([UInt32[]]; params=C.CandParams(remainder=:popgap, pop_split=7))
        @test_throws ArgumentError C.cand_encode([UInt32[]]; params=C.CandParams(gap_deg_buckets=0))
    end
    @testset "default bytes are frozen" begin
        # Default parameters keep the first shipped layout: 21 fields, and the same stream as before :popgap existed.
        pb = C.cand_params_bytes(C.CandParams())
        @test length(pb) == 2 + 4 * 21 && pb[1] == 0x01 && pb[2] == 21
        @test C.cand_params_bytes(C.CandParams(remainder=:cell))[2] == 21
        @test C.cand_params_bytes(C.CandParams(remainder=:popgap))[2] == 24   # extended layout: v1 readers must reject it
        @test C.cand_params_bytes(C.CandParams(gap_deg_buckets=2))[2] == 24
        @test C.cand_params_from_bytes(pb) == C.CandParams()
        b3 = C.cand_params_bytes(C.CandParams(remainder=:popgap)); b3[2 + 4 * 18 + 1] = 0x03
        @test_throws ArgumentError C.cand_params_from_bytes(b3)
        # deterministic graph (own LCG, independent of Julia's RNG streams)
        st = UInt64(88172645463325252); nx() = (st = st * 6364136223846793005 + 1442695040888963407; Int(st >> 33))
        n = 3000; sets = [Set{UInt32}() for _ in 1:n]
        for v in 1:n, _ in 1:(nx() % 9)
            x = nx() % 10 < 7 ? clamp(v + nx() % 41 - 20, 1, n) : nx() % n + 1
            push!(sets[v], UInt32(x)); nx() % 3 == 0 && push!(sets[x], UInt32(v))
        end
        rows = [sort!(collect(q)) for q in sets]
        b = C.cand_encode(rows)
        fnv = foldl((h, x) -> (h ⊻ UInt64(x)) * 0x00000100000001b3, b; init=0xcbf29ce484222325)
        @test C.cand_decode(b, n) == rows
        @test (length(b), fnv) == FROZEN_DEFAULT
    end
    @testset "popgap through the MGS container" begin
        rows = clustered(rng, 1500, 30, 6)
        for r in rows; filter!(x -> true, r); end
        g = Adjacently.MGS.LightGraphs.SimpleDiGraph{UInt32}(length(rows))
        for v in 1:length(rows), u in rows[v]; Adjacently.MGS.LightGraphs.add_edge!(g, UInt32(v), u); end
        for (p, kw) in ((C.CandParams(remainder=:popgap), (;)), (C.CandParams(remainder=:popgap, gap_deg_buckets=4, pop_split=3), (;)),
                        (C.CandParams(remainder=:popgap), (; transpose=:on)), (C.CandParams(remainder=:popgap), (; reciprocal=:on)), (C.CandParams(), (;)))
            base = tempname()
            Adjacently.MGS.write_cand_mgs3_graph(g, base; params=p, kw...)
            @test Adjacently.MGS.load_adjacency_mgs3_graph(base * ".mgz") == rows
            @test Adjacently.MGS.load_adjacency_mgs3_graph(base * ".mgz"; optimize_decode=true) == rows
            rm(base * ".mgz"; force=true)
        end
    end
    @testset "tuning graphs" begin
        dir = get(ENV, "ADJACENTLY_CSR_DIR", "")   # optional: directory with <graph>_community.csr inputs
        for g in ("wiki-Vote", "eat")
            f = joinpath(dir, g * "_community.csr")
            if !isfile(f); @test_skip false; continue; end
            rows = open(f) do io
                read(io, UInt64); n = Int(read(io, UInt32)); off = Vector{UInt64}(undef, n + 1); read!(io, off)
                [(a = Vector{UInt32}(undef, off[i+1] - off[i]); read!(io, a); a .+= UInt32(1); a) for i in 1:n]
            end
            @test roundtrip(rows)[1]
        end
    end
@testset "popgap parameters are unreadable by v1 readers" begin
    C2 = Adjacently.Compression
    @test C2.cand_params_bytes(C2.CandParams())[2] == 21
    b = C2.cand_params_bytes(C2.CandParams(remainder=:popgap))
    @test b[2] > 21                       # a v1 reader accepts only 21 fields
    @test C2.cand_params_from_bytes(b) == C2.CandParams(remainder=:popgap)
end
end
