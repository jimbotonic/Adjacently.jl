# Focused, dataset-independent regression suite for packed readers and decode workspaces.
# julia --project --startup-file=no test/run_tests_decode_workspace.jl
using Test, Random, Logging, Adjacently, LightGraphs
const DIO = Adjacently.IO
const DC = Adjacently.Compression
const DM = Adjacently.MGS
global_logger(ConsoleLogger(stderr, Logging.Warn))

@testset "Packed bit reader" begin
    rng = MersenneTwister(53)
    bytes = rand(rng, UInt8, 40)
    bits = [((b >> i) & 1) != 0 for b in bytes for i in 7:-1:0]
    r = DIO.BitReader(bytes)
    @test r.buffer === bytes
    @test length(r.buffer) == 40
    @test r.length == 320
    @test DIO.peek_bit(r) == bits[1]
    @test r.bit_count == 0
    @test DIO.read_bits(r, 320) == bits
    @test_throws ErrorException DIO.read_bit(r)
    @test_throws ErrorException DIO.peek_bit(r)
    @test DIO.read_value(r, 0, UInt64) == 0
    @test_throws ArgumentError DIO.read_value(r, -1, UInt64)
    for T in (UInt8, UInt16, UInt32, UInt64, UInt128,
              Adjacently.CustomTypes.UInt24, Adjacently.CustomTypes.UInt40)
        for start in 0:15, n in 0:min(8sizeof(T), 128)
            r = DIO.BitReader(bytes)
            r.index = start + 1; r.bit_count = start
            want = zero(T)
            for b in bits[start+1:start+n]; want = (want << 1) | T(b); end
            @test DIO.read_value(r, n, T) == want
            @test r.index == start + n + 1
            @test r.bit_count == start + n
        end
    end
    a = DIO.BitReader(bytes); b = DIO.BitReader(bytes)
    DIO.read_value(a, 13, UInt16)
    @test b.index == 1
    @test a.buffer === b.buffer
    io = IOBuffer(bytes); seek(io, 3)
    @test DIO.read_bits(DIO.BitReader(io), 296) == bits[25:end]
    @test_throws ErrorException DIO.read_value(DIO.BitReader(UInt8[0]), 9, UInt16)
    @test_throws ErrorException DIO.read_bit(DIO.BitReader(UInt8[]))
end

@testset "Integer code roundtrips" begin
    for T in (UInt16, UInt32, UInt64, Adjacently.CustomTypes.UInt24),
        (writefn, readfn) in ((DC.write_fibonacci, DC.read_fibonacci),
                             (DC.write_elias_gamma, DC.read_elias_gamma),
                             (DC.write_elias_delta, DC.read_elias_delta))
        vals = T.([1, 2, 3, 4, 7, 8, 255, 256, 1023, 16384, 32767])
        io = IOBuffer(); w = DIO.BitWriter(io)
        for v in vals; writefn(w, v); end
        DIO.flush_bitwriter(w; flush_last_bits=true)
        r = DIO.BitReader(take!(io))
        @test [readfn(r, T) for _ in vals] == vals
    end
end

@testset "Range spans and model reuse" begin
    rng = MersenneTwister(54)
    values = [UInt64[], UInt64[0], rand(rng, UInt64(0):UInt64(1_000_000), 1024), zeros(UInt64, 300_000)]
    blobs = map(values) do vs
        e = DC.CtxRangeEncoder()
        for v in vs; DC.rc_encode_value!(e, v); end
        DC.rc_finish!(e)
    end
    backing = vcat(UInt8[0xff, 0xff], blobs..., fill(0xff, 16))
    d = DC.CtxRangeDecoder(UInt8[]); freq = d.freq; tot = d.tot
    offset = 2
    for (vs, blob) in zip(values, blobs)
        DC.reset_range_decoder!(d, backing, offset + 1, offset + length(blob))
        @test d.freq === freq && d.tot === tot && d.inp === backing
        @test [DC.rc_decode_value!(d) for _ in vs] == vs
        fresh = DC.CtxRangeDecoder(backing, offset + 1, offset + length(blob))
        @test [DC.rc_decode_value!(fresh) for _ in vs] == vs
        offset += length(blob)
    end
    for first in 1:5
        @test DC.CtxRangeDecoder(fill(0xff, 5), first, first - 1).code == 0
        @test DC.BinRangeDecoder(fill(0xff, 5), first, first - 1).code == 0
    end
    for D in (DC.CtxRangeDecoder, DC.BinRangeDecoder)
        @test_throws ArgumentError D(UInt8[0], 0, 1)
        @test_throws ArgumentError D(UInt8[0], 1, 2)
    end
    chunks = [Bool[], [true], rand(rng, 100_000) .< 0.85, falses(100_000)]
    bb = map(chunks) do bs
        e = DC.BinRangeEncoder()
        for b in bs; DC.brc_encode_bit!(e, b); end
        DC.brc_finish!(e)
    end
    bytes = vcat(bb..., fill(0xff, 16)); offset = 0
    d = DC.BinRangeDecoder(UInt8[]); c0 = d.c0; c1 = d.c1
    for (bs, blob) in zip(chunks, bb)
        DC.reset_range_decoder!(d, bytes, offset + 1, offset + length(blob))
        @test d.c0 === c0 && d.c1 === c1 && d.inp === bytes
        @test [DC.brc_decode_bit!(d) for _ in bs] == bs
        offset += length(blob)
    end
end

@testset "Indexed continuation validation" begin
    state = DC.IndexedRangeState(UInt16)
    bytes = UInt8[]
    @test_throws ArgumentError DC._validate_indexed_resume(state,nothing,true,16,192,bytes)
    @test_throws ArgumentError DC._validate_indexed_resume(state,1:0,true,16,192,bytes)
    @test_throws ArgumentError DC._validate_indexed_resume(state,1:8,false,16,192,bytes)
    @test_throws ArgumentError DC._validate_indexed_resume(state,2:8,true,16,192,bytes)
    @test_throws ArgumentError DC._validate_indexed_resume(state,1:17,true,16,192,bytes)
    @test DC._validate_indexed_resume(state,1:8,true,16,192,bytes) === nothing
    state.next_vertex=9; state.resid=DC.CtxRangeDecoder(bytes)
    @test DC._validate_indexed_resume(state,9:16,true,16,192,bytes) === nothing
    @test_throws ArgumentError DC._validate_indexed_resume(state,8:16,true,16,192,bytes)
    @test_throws ArgumentError DC._validate_indexed_resume(state,9:16,true,16,192,UInt8[])
    state.next_vertex=17
    @test_throws ArgumentError DC._validate_indexed_resume(state,17:20,true,16,192,bytes)
end

@testset "Graph roundtrips and nested random access" begin
    # Empty lists, intervals, high degree, repeated lists, and cross-chunk references.
    g = Adjacently.SimpleDiGraph{UInt16}()
    add_vertices!(g, 192)
    for u in 1:184, v in 1:192
        if u == 1 || (u % 11 != 0 && (20 <= v <= 45 || 110 <= v <= 130 || v == mod(u + 60, 192) + 1))
            u != v && add_edge!(g, u, v)
        end
    end
    mktempdir() do dir
        for writer in (DM.write_bg_mgs3_graph, DM.write_cs_mgs3_graph)
            @test_throws ArgumentError writer(g, joinpath(dir,"invalid"); chunk_local_refs=true)
        end
        for alg in (:bg, :cs, :cg), backend in (:fibonacci, :context_range), mode in (:children, :index)
            @testset "$alg/$backend/$mode" begin
                base = joinpath(dir, "$(alg)_$(backend)_$(mode)")
                if alg == :bg
                    DM.write_bg_mgs3_graph(g, base; integer_encoding=backend, coding_scheme=mode,
                        ref_window_size=64, lr_split=true, multi_ref=true, cost_model=0, index_sample_k=16)
                elseif alg == :cs
                    DM.write_cs_mgs3_graph(g, base; integer_encoding=backend, coding_scheme=mode,
                        ref_window_size=64, lr_split=true, cost_model=0, index_sample_k=16)
                else
                    clusters = [UInt16.(i:min(i+31,192)) for i in 1:32:192]
                    # Use exactly the configuration represented by the file header.
                    p = DM.decode_cg_params(DM.encode_cg_params(DC.CG.CGParams(
                        intra_ref_window=64, intra_mil=4, intra_adapt_mil=4,
                        intra_intervals=true, intra_lr_split=true)))
                    DM.write_cg_mgs3_graph(g, base, clusters; params=p,
                        integer_encoding=backend, coding_scheme=mode)
                end
                path = base * ".mgz"
                loaded = DM.load_compressed_mgs3_graph(path)
                @test nv(loaded) == nv(g)
                @test ne(loaded) == ne(g)
                @test all(outneighbors(loaded, u) == outneighbors(g, u) for u in 1:192)
                if mode == :index && backend == :context_range
                    idx = DM.load_indexed_mgs3_graph(path)
                    for u in randperm(MersenneTwister(55),192)
                        idx.reset_fn()
                        @test idx.neighbors(u) == outneighbors(g,u)
                        @test idx.neighbors(u) == outneighbors(g,u)
                    end
                    @test DC._RESID_SOURCE[] === nothing
                    @test DC._FLAG_SOURCE[] === nothing
                    if alg != :cg
                        writer = alg == :bg ? DM.write_bg_mgs3_graph : DM.write_cs_mgs3_graph
                        bounded = base * "_bounded"
                        writer(g, bounded; integer_encoding=backend, coding_scheme=mode,
                            ref_window_size=64, lr_split=true, cost_model=0,
                            index_sample_k=16, chunk_local_refs=true)
                        bidx = DM.load_indexed_mgs3_graph(bounded*".mgz")
                        for u in randperm(MersenneTwister(56),192)
                            bidx.reset_fn()
                            @test bidx.neighbors(u) == outneighbors(g,u)
                            @test bidx.stats.chunks_decoded[] == 1
                            @test bidx.stats.vertices_decoded[] == mod(u-1,16)+1
                        end
                        # Ascending requests resume live state without replay;
                        # backward/memo hits don't decode anything or expose owned lists.
                        bidx.reset_fn()
                        for u in 1:192
                            @test bidx.neighbors(u) == outneighbors(g,u)
                            @test bidx.stats.vertices_decoded[] == u
                            @test bidx.stats.active_workspaces() == (u % 16 == 0 ? 0 : 1)
                        end
                        @test bidx.stats.arcs_decoded[] == ne(g)
                        @test bidx.stats.chunks_decoded[] == 12
                        @test bidx.stats.decode_calls[] == 192
                        returned = bidx.neighbors(1); empty!(returned)
                        @test bidx.neighbors(1) == outneighbors(g,1)
                        @test bidx.stats.vertices_decoded[] == 192
                        bidx.reset_fn()
                        @test bidx.stats.active_workspaces() == 0
                        for u in randperm(MersenneTwister(59),192)
                            @test bidx.neighbors(u) == outneighbors(g,u)
                        end
                        @test bidx.stats.vertices_decoded[] == 192
                        @test bidx.stats.arcs_decoded[] == ne(g)
                        @test bidx.stats.active_workspaces() == 0
                        @test_throws ErrorException bidx.neighbors(0)
                        @test_throws ErrorException bidx.neighbors(193)
                        fullidx = DM.load_indexed_mgs3_graph(bounded*".mgz"; decode_prefix=false)
                        @test fullidx.neighbors(1) == outneighbors(g,1)
                        @test fullidx.stats.vertices_decoded[] == 16
                        @test fullidx.stats.active_workspaces() == 0
                        @test_throws ArgumentError DM.load_indexed_mgs3_graph(bounded*".mgz";max_live_chunks=0)
                        tinycache = DM.load_indexed_mgs3_graph(bounded*".mgz";max_live_chunks=1)
                        for u in (1,17,2,18,16,32,33,49,34,50,48,64)
                            @test tinycache.neighbors(u) == outneighbors(g,u)
                            @test tinycache.stats.active_workspaces() <= 1
                        end
                        @test tinycache.stats.chunks_decoded[] == 4
                        @test tinycache.stats.vertices_decoded[] > 64 # evicted prefixes replay
                        # Legacy global-reference files: nested calls must preserve
                        # their caller's live models, including across many resumes.
                        idx.reset_fn()
                        for u in randperm(MersenneTwister(60),192)
                            @test idx.neighbors(u) == outneighbors(g,u)
                        end
                        @test idx.stats.vertices_decoded[] == 192
                        @test idx.stats.arcs_decoded[] == ne(g)
                        @test idx.stats.active_workspaces() == 0
                        @test DC._RESID_SOURCE[] === nothing
                        @test DC._FLAG_SOURCE[] === nothing
                        # Eviction can remove an ancestor still active on the
                        # resolver stack. Such models must not enter the idle pool.
                        nested = DM.load_indexed_mgs3_graph(path; max_live_chunks=1)
                        for u in randperm(MersenneTwister(74),192)
                            nested.reset_fn()
                            @test nested.neighbors(u) == outneighbors(g,u)
                            @test nested.stats.active_workspaces() <= 1
                        end
                        nested.reset_fn()
                        @test nested.stats.scratch_bytes() > 0
                    end
                end
            end
        end
    end
end

@testset "Owned reconstruction and counted residuals" begin
    rng = MersenneTwister(57)
    for T in (UInt16, UInt32, UInt64, Adjacently.CustomTypes.UInt24), _ in 1:100
        a = sort(T.(rand(rng,1:200,rand(rng,0:80))))
        b = sort(T.(rand(rng,1:200,rand(rng,0:80))))
        original_b = copy(b); want = sort(vcat(a,b))
        @test DC._merge_sorted_owned!(a,b) === a
        @test a == want
        @test b == original_b
        copied = sort!(unique!(T.(rand(rng,1:200,rand(rng,0:80)))))
        residual = sort!(setdiff(unique(T.(rand(rng,1:200,rand(rng,0:80)))),copied))
        ranks = DC._rankgap_forward(residual,copied)
        original_ranks = copy(ranks)
        @test DC._rankgap_inverse(ranks,copied) == residual
        @test ranks == original_ranks
        @test DC._rankgap_inverse!(ranks,copied) === ranks
        @test ranks == residual
    end
    x = UInt16[1,2]
    @test_throws ArgumentError DC._merge_sorted_owned!(x,x)
    r = DIO.BitReader(UInt8[])
    @test isempty(DC.read_delta(r,:fibonacci,UInt16;max_elements=0))
    @test r.index == 1
    @test_throws ArgumentError DC.read_delta(r,:fibonacci,UInt16;max_elements=-1)
    @test isempty(DC.read_delta(r,:fibonacci,UInt16;max_elements=10))
    io = IOBuffer(); w = DIO.BitWriter(io)
    DC.write_delta(w,UInt16[10,20,25],:fibonacci)
    DIO.flush_bitwriter(w;flush_last_bits=true)
    bytes = take!(io)
    @test DC.read_delta(DIO.BitReader(bytes),:fibonacci,UInt16;max_elements=3) == UInt16[10,20,25]
    # Truncated counted input retains the historical partial-list behavior.
    @test DC.read_delta(DIO.BitReader(bytes),:fibonacci,UInt16;max_elements=100)[1:3] == UInt16[10,20,25]
    ref = UInt16[10,20,30,40]; flags = Bool[true,false,true,false]
    e = DC.BinRangeEncoder()
    for f in flags; DC.brc_encode_bit!(e,f); end
    @test DC._read_copied_reference(DC.BinRangeDecoder(DC.brc_finish!(e)),ref) == UInt16[10,30]
    @test ref == UInt16[10,20,30,40]
end
