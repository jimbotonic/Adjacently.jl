using Test,Random,Logging,LightGraphs,Adjacently,SHA,JSON
global_logger(ConsoleLogger(stderr,Logging.Warn))
const C=Adjacently.Compression;const M=Adjacently.MGS;const AIO=Adjacently.IO
function fixture_graph(::Type{T}=UInt16) where T
    g=SimpleDiGraph{T}();add_vertices!(g,120);rng=MersenneTwister(83)
    for u in 1:120,v in 1:120
        if u!=v && ((u-1)÷12==(v-1)÷12 ? rand(rng)<.3 : rand(rng)<.004)
            add_edge!(g,T(u),T(v))
        end
    end
    return g
end
function stream(g,alg;parallel=false,workers=4,batch=17,ie=:fibonacci,cost=0,flags=(;))
    T=eltype(g);nls=Dict(T(v)=>sort(collect(outneighbors(g,v))) for v in 1:nv(g))
    w=AIO.BitWriter();stats=Dict()
    result=if alg==:bg
        C.write_greedy_graph_data(w,nls,:children,16;integer_encoding=ie,
            cost_model=cost,stats=stats,parallel_search=parallel,search_workers=workers,
            search_batch_size=batch,flags...)
    else
        C.write_cmdstream_graph_data(w,nls,:children,16;integer_encoding=ie,
            cost_model=cost,parallel_search=parallel,search_workers=workers,
            search_batch_size=batch,flags...)
    end
    AIO.flush_bitwriter(w;flush_last_bits=true)
    return (collect(AIO.get_bytes(w)),result,stats)
end
@testset "Parallel search, byte-identical decisions and streams" begin
    for T in (UInt16,Adjacently.CustomTypes.UInt24)
        g=fixture_graph(T)
        for alg in (:bg,:cs),ie in (:fibonacci,:context_range),cost in (0,1)
            flags=alg==:bg ? (;copy_blocks=true,adaptive_copy=true,compact_copy=true,
                stop_deltas=true,lr_split=true,multi_ref=true,tight_intervals=true) : (;lr_split=true)
            serial=stream(g,alg;ie,cost,flags)
            for (workers,batch) in ((1,1),(2,7),(4,31),(100,4096))
                @test stream(g,alg;parallel=true,workers,batch,ie,cost,flags)==serial
            end
        end
    end
    g=fixture_graph()
    for flags in ((;), (;adaptive_header=true,stop_deltas=true),
                  (;copy_blocks=true,adaptive_copy=true,bv_blocks=true,split_residual=true,adaptive_deltas=true),
                  (;fixwidth_ref=true,adaptive_copy=true,copy_blocks=true,compact_copy=true))
        for ie in (:fibonacci,:context_range)
            @test stream(g,:bg;parallel=true,ie,flags)==stream(g,:bg;ie,flags)
        end
    end
    for n in (0,1,2,33)
        empty=SimpleDiGraph{UInt16}(n)
        for alg in (:bg,:cs),ie in (:fibonacci,:context_range)
            @test stream(empty,alg;parallel=true,ie)==stream(empty,alg;ie)
        end
    end
    # Search results must not alias the next vertex's reused workspace.
    @test stream(g,:bg;parallel=true,batch=120,flags=(;multi_ref=true))==stream(g,:bg;flags=(;multi_ref=true))
end
@testset "Unsupported modes fail before file creation" begin
    g=fixture_graph();tmp=mktempdir()
    for writer in (M.write_bg_mgs3_graph,M.write_cs_mgs3_graph)
        @test_throws ArgumentError writer(g,joinpath(tmp,"bad");parallel_search=true,coding_scheme=:index)
        @test_throws ArgumentError writer(g,joinpath(tmp,"bad");parallel_search=true,search_workers=0)
        @test_throws ArgumentError writer(g,joinpath(tmp,"bad");parallel_search=true,search_batch_size=0)
    end
    @test_throws ArgumentError M.write_bg_mgs3_graph(g,joinpath(tmp,"bad");parallel_search=true,exact_costing=true)
    @test !isfile(joinpath(tmp,"bad.mgz"))
    # Worker errors propagate, with no emission of an incompletely searched batch.
    nls=Dict(UInt16(1)=>UInt16[1],UInt16(2)=>UInt16[1]);seen=Int[]
    @test_throws Exception C._parallel_vertex_batches!((args...)->error("worker failure"),
        (v,args...)->push!(seen,Int(v)),nls,2,1;workers=4,batch_size=2)
    @test isempty(seen)
end
@testset "Frozen pre-change EAT bytes and roundtrip" begin
    legacy=joinpath(@__DIR__,"fixtures/codec_legacy_v1.json")
    @test isfile(legacy)
    if isfile(legacy)
        hashes=JSON.parsefile(legacy)["files"]
        g,_,_=Adjacently.Graph.get_core(AIO.load_graph_from_pajek(joinpath(@__DIR__,"../datasets/EAT/EATnew.net")))
        tmp=mktempdir()
        for alg in (:bg,:cs),ie in (:fibonacci,:context_range),parallel in (false,true)
            base=joinpath(tmp,"eat_$(alg)_$(ie)")
            if alg==:bg
                M.write_bg_mgs3_graph(g,base;integer_encoding=ie,ref_window_size=64,lr_split=true,multi_ref=true,cost_model=0,parallel_search=parallel,search_batch_size=127)
            else
                M.write_cs_mgs3_graph(g,base;integer_encoding=ie,ref_window_size=256,lr_split=true,cost_model=0,parallel_search=parallel,search_batch_size=127)
            end
            @test bytes2hex(sha256(read(base*".mgz")))==hashes[basename(base)*".mgz"]
            decoded=M.load_compressed_mgs3_graph(base*".mgz")
            @test nv(decoded)==nv(g) && ne(decoded)==ne(g)
            @test all(outneighbors(decoded,v)==outneighbors(g,v) for v in 1:nv(g))
        end
    end
end
