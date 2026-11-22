using Test, Random, Adjacently, LightGraphs, SHA
using Adjacently.MGS: write_bg_mgs3_graph, write_cs_mgs3_graph, write_cg_mgs3_graph,
    load_adjacency_mgs3_graph, load_compressed_mgs3_graph, load_indexed_mgs3_graph
using Adjacently.Compression.CG: CGParams
const M = Adjacently.MGS

function rand_digraph(rng, n, m; symmetric=false)
    g = SimpleDiGraph{UInt32}(n)
    while ne(g) < m
        u = rand(rng, 1:n); v = clamp(u + rand(rng, -40:40), 1, n)
        add_edge!(g, u, v); symmetric && add_edge!(g, v, u)
    end
    g
end
adjof(g) = [sort!(collect(UInt32, outneighbors(g, v))) for v in 1:nv(g)]
cgp(w) = CGParams(optimize_costs=true, L=128, varint=:fibonacci, count_varint=:fibonacci, gap=:fibonacci,
    degree=:elias_delta, undirected_pairs=false, perm_strategy=:blockpos, membership=:implicit_ranges,
    inter_strategy=:perm, intra_ref_enabled=true, intra_ref_window=w, intra_ref_rle=false, intra_block_try=false,
    positions_mode=:delta, additions_mode=:delta, min_cluster_density=0.0, intra_intervals=true, intra_mil=4,
    intra_lr_split=true, intra_tight_deltas=false, intra_greedy_mil=false, intra_zigzag=true, intra_stop_deltas=true,
    intra_copy_blocks=true, intra_copy_adaptive=true, intra_ref_fixwidth=true, intra_ref_vlc=false,
    intra_add_adaptive=true, intra_raw_adaptive=true, intra_adapt_mil=4)
one(n) = [collect(UInt32(1):UInt32(n))]
minorbyte(path) = read(path)[5]

@testset "MGS header extension" begin
    rng = MersenneTwister(7)
    dir = mktempdir()
    g = rand_digraph(rng, 3000, 24000)
    a = adjof(g)

    @testset "no extension => byte-identical to the core writer" begin
        b1 = joinpath(dir, "core"); b2 = joinpath(dir, "pub")
        M._write_cg_mgs3_graph_core(g, b1, one(3000); integer_encoding=:context_range, params=cgp(256))
        write_cg_mgs3_graph(g, b2, one(3000); integer_encoding=:context_range, params=cgp(256))
        @test read(b1 * ".mgz") == read(b2 * ".mgz")
        @test minorbyte(b2 * ".mgz") & 0x80 == 0
        for (w, f) in ((M._write_bg_mgs3_graph_core, write_bg_mgs3_graph), (M._write_cs_mgs3_graph_core, write_cs_mgs3_graph))
            w(g, b1; integer_encoding=:context_range); f(g, b2; integer_encoding=:context_range)
            @test read(b1 * ".mgz") == read(b2 * ".mgz")
        end
    end

    @testset "CG windows beyond the byte2 table" begin
        prev = typemax(Int)
        for w in (512, 1000, 1024, 4096)
            b = joinpath(dir, "cgw$w")
            write_cg_mgs3_graph(g, b, one(3000); integer_encoding=:context_range, params=cgp(w))
            f = b * ".mgz"
            @test minorbyte(f) & 0x80 != 0
            @test load_adjacency_mgs3_graph(f) == a
            @test load_adjacency_mgs3_graph(f; optimize_decode=true) == a
            @test adjof(load_compressed_mgs3_graph(f)) == a
            # extension = varint(len) + tag + varint(len) + varint(window): 4-6 bytes
            b256 = joinpath(dir, "cg256")
            write_cg_mgs3_graph(g, b256, one(3000); integer_encoding=:context_range, params=cgp(256))
            @test filesize(f) - filesize(b256 * ".mgz") in -100000:100000   # sizes comparable
        end
        # the overhead itself is exactly the extension block
        @test length(M._mgs_ext_block([M.MGS_EXT_CG_WINDOW => M._mgs_varint_bytes(1024)])) == 5   # len + tag + len + 2-byte varint
    end

    @testset "symmetric mode" begin
        gs = rand_digraph(rng, 3000, 30000; symmetric=true)
        as = adjof(gs)
        for (name, wr) in (("bg", (h, b; kw...) -> write_bg_mgs3_graph(h, b; integer_encoding=:context_range, kw...)),
                           ("cs", (h, b; kw...) -> write_cs_mgs3_graph(h, b; integer_encoding=:context_range, kw...)),
                           ("cg", (h, b; kw...) -> write_cg_mgs3_graph(h, b, one(3000); integer_encoding=:context_range, params=cgp(256), kw...)),
                           ("cg1024", (h, b; kw...) -> write_cg_mgs3_graph(h, b, one(3000); integer_encoding=:context_range, params=cgp(1024), kw...)))
            bp = joinpath(dir, "sym_$(name)_plain"); bs = joinpath(dir, "sym_$name")
            wr(gs, bp); wr(gs, bs; symmetric=:auto)
            fs = bs * ".mgz"
            @test minorbyte(fs) & 0x80 != 0
            @test filesize(fs) < 0.7 * filesize(bp * ".mgz")
            @test load_adjacency_mgs3_graph(fs) == as
            @test load_adjacency_mgs3_graph(fs; optimize_decode=true) == as
            @test adjof(load_compressed_mgs3_graph(fs)) == as
            # :auto on a non-symmetric graph changes nothing
            bn = joinpath(dir, "nsym_$name"); wr(g, bn; symmetric=:auto); wr(g, bp)
            @test read(bn * ".mgz") == read(bp * ".mgz")
            @test_throws ArgumentError wr(g, bn; symmetric=:on)
        end
        # self-loops survive (stored once, not mirrored twice)
        gl = copy(gs); add_edge!(gl, 5, 5); add_edge!(gl, 77, 77)
        b = joinpath(dir, "loops"); write_cs_mgs3_graph(gl, b; integer_encoding=:context_range, symmetric=:on)
        @test load_adjacency_mgs3_graph(b * ".mgz") == adjof(gl)
        @test_throws ArgumentError write_bg_mgs3_graph(gs, joinpath(dir, "ix"); integer_encoding=:context_range,
            coding_scheme=:index, symmetric=:on)
    end

    @testset "reciprocal mode" begin
        # partially reciprocal graph (~50% of arcs returned) + a few self-loops
        gr = rand_digraph(rng, 3000, 20000)
        for u in 1:3000, v in collect(outneighbors(gr, u)); rand(rng) < 0.5 && add_edge!(gr, v, u); end
        add_edge!(gr, 9, 9)
        ar = adjof(gr)
        for (name, wr) in (("bg", (h, b; kw...) -> write_bg_mgs3_graph(h, b; integer_encoding=:context_range, kw...)),
                           ("cs", (h, b; kw...) -> write_cs_mgs3_graph(h, b; integer_encoding=:context_range, kw...)),
                           ("cg1024", (h, b; kw...) -> write_cg_mgs3_graph(h, b, one(3000); integer_encoding=:context_range, params=cgp(1024), kw...)))
            bp = joinpath(dir, "rec_$(name)_plain"); br = joinpath(dir, "rec_$name"); ba = joinpath(dir, "rec_$(name)_auto")
            wr(gr, bp); wr(gr, br; reciprocal=:on); wr(gr, ba; reciprocal=:auto)
            fr = br * ".mgz"
            @test minorbyte(fr) & 0x80 != 0
            @test load_adjacency_mgs3_graph(fr) == ar
            @test load_adjacency_mgs3_graph(fr; optimize_decode=true) == ar
            @test adjof(load_compressed_mgs3_graph(fr)) == ar
            @test filesize(ba * ".mgz") == min(filesize(bp * ".mgz"), filesize(fr))
            @test load_adjacency_mgs3_graph(ba * ".mgz") == ar
            @test !isfile(ba * ".recip_tmp.mgz")
        end
        # graph with no reciprocity: :on still round-trips; :auto keeps the plain file
        b1 = joinpath(dir, "rec_none"); b2 = joinpath(dir, "rec_none_plain")
        write_cs_mgs3_graph(g, b1; integer_encoding=:context_range, reciprocal=:on)
        @test load_adjacency_mgs3_graph(b1 * ".mgz") == a
        # symmetric graphs: :auto symmetric wins, reciprocal skipped; :on + symmetric refused
        gs = rand_digraph(rng, 2000, 12000; symmetric=true)
        bs = joinpath(dir, "rec_sym"); write_bg_mgs3_graph(gs, bs; integer_encoding=:context_range, symmetric=:auto, reciprocal=:auto)
        @test load_adjacency_mgs3_graph(bs * ".mgz") == adjof(gs)
        @test_throws ArgumentError write_bg_mgs3_graph(gs, bs; integer_encoding=:context_range, symmetric=:on, reciprocal=:on)
        @test_throws ArgumentError write_cs_mgs3_graph(gr, bs; integer_encoding=:context_range, coding_scheme=:index, reciprocal=:on)
    end

    @testset "transposed storage" begin
        # skewed graph: many sources pointing into a few hub blocks, so G and G^T differ in cost
        gt = SimpleDiGraph{UInt32}(3000)
        for u in 1:3000, k in 1:6; add_edge!(gt, u, 1 + (u * 7 + k * 13) % 200); end
        for u in collect(1:3000)[1:40:end], v in collect(outneighbors(gt, u)); add_edge!(gt, v, u); end
        at = adjof(gt)
        for (name, wr) in (("bg", (h, b; kw...) -> write_bg_mgs3_graph(h, b; integer_encoding=:context_range, kw...)),
                           ("cs", (h, b; kw...) -> write_cs_mgs3_graph(h, b; integer_encoding=:context_range, kw...)),
                           ("cg", (h, b; kw...) -> write_cg_mgs3_graph(h, b, one(3000); integer_encoding=:context_range, params=cgp(256), kw...)))
            b0 = joinpath(dir, "tr_$(name)_plain"); b1 = joinpath(dir, "tr_$name"); b2 = joinpath(dir, "tr_$(name)_auto")
            b3 = joinpath(dir, "tr_$(name)_recip")
            wr(gt, b0); wr(gt, b1; transpose=:on); wr(gt, b2; transpose=:auto); wr(gt, b3; transpose=:on, reciprocal=:on)
            for f in (b1, b2, b3) .* ".mgz"
                @test load_adjacency_mgs3_graph(f) == at
                @test load_adjacency_mgs3_graph(f; optimize_decode=true) == at
                @test adjof(load_compressed_mgs3_graph(f)) == at
            end
            @test minorbyte(b1 * ".mgz") & 0x80 != 0
            @test filesize(b2 * ".mgz") == min(filesize(b0 * ".mgz"), filesize(b1 * ".mgz"))
            @test !isfile(b2 * ".transpose_tmp.mgz")
        end
        # symmetric graph: :auto does not transpose (identical), symmetric mode still applies
        gs = rand_digraph(rng, 1500, 9000; symmetric=true)
        bs = joinpath(dir, "tr_sym"); write_cs_mgs3_graph(gs, bs; integer_encoding=:context_range, symmetric=:auto, transpose=:auto)
        @test load_adjacency_mgs3_graph(bs * ".mgz") == adjof(gs)
        @test_throws ArgumentError write_cs_mgs3_graph(gt, bs; integer_encoding=:context_range, coding_scheme=:index, transpose=:on)
    end

    @testset "candidate coder (ALG_CAND)" begin
        C = Adjacently.Compression
        gr = rand_digraph(rng, 2500, 15000)
        for u in 1:2500, v in collect(outneighbors(gr, u)); rand(rng) < 0.4 && add_edge!(gr, v, u); end
        ar = adjof(gr)
        b = joinpath(dir, "cand"); M.write_cand_mgs3_graph(gr, b); f = b * ".mgz"
        @test read(f)[7] == M.ALG_CAND
        @test minorbyte(f) == 0x02                       # no extension
        @test load_adjacency_mgs3_graph(f) == ar
        @test load_adjacency_mgs3_graph(f; optimize_decode=true) == ar
        @test adjof(load_compressed_mgs3_graph(f)) == ar
        # file = header + varint + params + payload, nothing else
        p = C.CandParams(); pb = C.cand_params_bytes(p)
        @test filesize(f) == 12 + length(M._mgs_varint_bytes(length(pb))) + length(pb) + length(C.cand_encode(ar; params=p))
        # non-default params are carried by the file
        p2 = C.CandParams(gen_2hop=false, ref_window=16, remainder=:cell)
        b2 = joinpath(dir, "cand_p2"); M.write_cand_mgs3_graph(gr, b2; params=p2)
        @test load_adjacency_mgs3_graph(b2 * ".mgz") == ar
        # storage modes
        for (name, kw) in (("tr", (; transpose=:on)), ("rec", (; reciprocal=:on)), ("auto", (; transpose=:auto, reciprocal=:auto)))
            bm = joinpath(dir, "cand_$name"); M.write_cand_mgs3_graph(gr, bm; kw...)
            @test minorbyte(bm * ".mgz") & 0x80 != 0 || name == "auto"
            @test load_adjacency_mgs3_graph(bm * ".mgz") == ar
            @test adjof(load_compressed_mgs3_graph(bm * ".mgz")) == ar
        end
        gs = rand_digraph(rng, 2000, 14000; symmetric=true)
        bs = joinpath(dir, "cand_sym"); M.write_cand_mgs3_graph(gs, bs; symmetric=:auto)
        @test minorbyte(bs * ".mgz") & 0x80 != 0
        @test load_adjacency_mgs3_graph(bs * ".mgz") == adjof(gs)
        @test adjof(load_compressed_mgs3_graph(bs * ".mgz")) == adjof(gs)
        @test_throws ArgumentError M.write_cand_mgs3_graph(gr, bs; coding_scheme=:index)
        d = read(f); d[5] = 0x03; bad = joinpath(dir, "cand_bad.mgz"); write(bad, d)
        @test_throws ArgumentError load_adjacency_mgs3_graph(bad)
        @test_throws ArgumentError load_adjacency_mgs3_graph(f; fused_lr=true)
    end

    @testset "malformed extensions are refused" begin
        b = joinpath(dir, "cgw512"); f = b * ".mgz"
        write_cg_mgs3_graph(g, b, one(3000); integer_encoding=:context_range, params=cgp(512))
        d = read(f); d[14] = 0x7e   # unknown tag
        bad = joinpath(dir, "bad.mgz"); write(bad, d)
        @test_throws ArgumentError load_adjacency_mgs3_graph(bad)
        @test_throws ArgumentError load_compressed_mgs3_graph(bad)
    end
end
