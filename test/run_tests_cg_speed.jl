using Test, Adjacently, LightGraphs, Random, Logging
const CG=Adjacently.Compression.CG
const CIO=Adjacently.IO
const CM=Adjacently.MGS
global_logger(ConsoleLogger(stderr,Logging.Warn))
function cgfixture(::Type{T}) where T
    g=SimpleDiGraph{T}();add_vertices!(g,120);rng=MersenneTwister(391)
    for u in 1:120,v in 1:120
        u!=v && rand(rng)<(abs(u-v)<12 ? .3 : .004) && add_edge!(g,T(u),T(v))
    end
    g
end

function cgparams(;kw...)
    CG.CGParams(;membership=:implicit_ranges,intra_ref_window=16,intra_adapt_mil=4,
        intra_intervals=true,intra_lr_split=true,intra_zigzag=true,
        intra_ref_fixwidth=true,intra_copy_blocks=true,intra_copy_adaptive=true,
        intra_stop_deltas=true,cost_model=0,kw...)
end
function cgstream(g,P,p,ctx,index,workers)
    w=CIO.BitWriter();stats=CG.CGStats();events=Any[]
    offsets=index ? zeros(Int,2length(P)+1) : nothing
    streams=CG.encode_level(w,g,P;params=p,ctx_range=ctx,cluster_offsets=offsets,
        stats=stats,progress=(args...)->push!(events,(args...,Threads.threadid())),
        parallel_search=workers>0,search_workers=max(1,workers))
    CIO.flush_bitwriter(w;flush_last_bits=true)
    (collect(CIO.get_bytes(w)),streams,offsets,[getfield(stats,k) for k in fieldnames(typeof(stats))],events)
end
@testset "CG parallel cost search" begin
    for T in (UInt16,Adjacently.CustomTypes.UInt24), k in (1,4),
            ctx in (false,true), index in (false,true), greedy in (false,true)
        g=cgfixture(T);P=[collect(T(a):T(min(a+120÷k-1,120))) for a in 1:120÷k:120]
        p=cgparams(;intra_greedy_mil=greedy)
        serial=cgstream(g,P,p,ctx,index,0)
        for workers in (1,2,4,100)
            @test cgstream(g,P,p,ctx,index,workers)==serial
        end
    end
    for n in (1,2)
        g=SimpleDiGraph{UInt16}();add_vertices!(g,n);P=[collect(UInt16(1):UInt16(n))]
        @test cgstream(g,P,cgparams(),true,false,4)==cgstream(g,P,cgparams(),true,false,0)
    end
    g=cgfixture(UInt16);P=[collect(UInt16(1):UInt16(120))];dir=mktempdir()
    @test_throws ArgumentError CM.write_cg_mgs3_graph(g,joinpath(dir,"bad"),P;parallel_search=true,search_workers=0)
    @test_throws ArgumentError CM.write_cg_mgs3_graph(g,joinpath(dir,"bad"),P;parallel_search=true,params=cgparams(intra_block_try=true))
    @test !isfile(joinpath(dir,"bad.mgz"))
end

@testset "CG explicit-source decoding and cluster parallelism" begin
    for T in (UInt16,Adjacently.CustomTypes.UInt24), k in (1,4),
            ctx in (false,true), index in (false,true), variant in (1,2,3)
        g=cgfixture(T);P=[collect(T(a):T(min(a+120÷k-1,120))) for a in 1:120÷k:120]
        p=variant==1 ? cgparams() : variant==2 ? cgparams(intra_greedy_mil=true,intra_tight_deltas=true) :
            CG.CGParams(membership=:implicit_ranges,intra_ref_fixwidth=true,intra_stop_deltas=true,
                intra_add_adaptive=true,intra_raw_adaptive=true,intra_copy_blocks=true,intra_copy_adaptive=true)
        bytes,streams,offs,_,_=cgstream(g,P,p,ctx,index,0)
        chunks=ctx && index ? (streams[4],streams[5],streams[6]) : nothing
        rc,rd,cp=ctx ? streams[1:3] : (UInt8[],UInt8[],UInt8[])
        for I in (Int,T), fused in (false,true), workers in (1,4)
            kw=ctx && index ? (;parallel_clusters=true,decode_workers=workers) : (;)
            nls=CG.decode_level(CIO.BitReader(bytes),p;T=T,ctx_range=ctx,
                coding_scheme=index ? :index : :children, resid_bytes=rc,refdist_bytes=rd,
                copy_bytes=cp,cg_offsets=offs,chunk_offsets=chunks,local_type=I,fused_lr=fused,kw...)
            @test all(get(nls,T(v),T[])==outneighbors(g,v) for v in 1:nv(g))
        end
        if ctx && index
            sentinel=Adjacently.Compression.CtxRangeDecoder(UInt8[])
            Adjacently.Compression._RESID_SOURCE[]=sentinel
            try
                nls=CG.decode_level(CIO.BitReader(bytes),p;T=T,ctx_range=true,coding_scheme=:index,
                    resid_bytes=rc,refdist_bytes=rd,copy_bytes=cp,cg_offsets=offs,
                    chunk_offsets=chunks,parallel_clusters=true)
                @test Adjacently.Compression._RESID_SOURCE[]===sentinel
                @test all(get(nls,T(v),T[])==outneighbors(g,v) for v in 1:nv(g))
            finally
                Adjacently.Compression._RESID_SOURCE[]=nothing
            end
        end
    end
    g=cgfixture(UInt16);P=[collect(UInt16(a):UInt16(a+29)) for a in 1:30:120]
    dir=mktempdir();base=joinpath(dir,"cg")
    CM.write_cg_mgs3_graph(g,base,P;params=cgparams(),integer_encoding=:context_range,
        coding_scheme=:index,parallel_search=true)
    for workers in (1,2,4,100)
        decoded=CM.load_compressed_mgs3_graph(base*".mgz";parallel_decode=true,decode_workers=workers)
        @test all(outneighbors(decoded,v)==outneighbors(g,v) for v in 1:nv(g))
    end
    @test_throws ArgumentError CM.load_compressed_mgs3_graph(base*".mgz";parallel_decode=true,decode_workers=0)
    CM.write_cg_mgs3_graph(g,base,P;params=cgparams(),integer_encoding=:context_range)
    @test_throws ArgumentError CM.load_compressed_mgs3_graph(base*".mgz";parallel_decode=true)
end
