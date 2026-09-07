include("run_tests_cg_speed.jl")
const C=Adjacently.Compression

@testset "CG opt-in decode switches" begin
    for T in (UInt16,Adjacently.CustomTypes.UInt24), k in (1,4), ctx in (false,true), index in (false,true), variant in 1:5
        g=cgfixture(T);P=[collect(T(a):T(min(a+120÷k-1,120))) for a in 1:120÷k:120]
        p=variant==1 ? cgparams() :
          variant==2 ? cgparams(intra_greedy_mil=true,intra_tight_deltas=true) :
          variant==3 ? CG.CGParams(membership=:implicit_ranges,intra_ref_fixwidth=true,intra_stop_deltas=true,
              intra_add_adaptive=true,intra_raw_adaptive=true,intra_copy_blocks=true,intra_copy_adaptive=true) :
          variant==4 ? cgparams(intra_zigzag=false) :
          CG.CGParams(membership=:implicit_ranges,intra_ref_fixwidth=true,additions_mode=:intervals)
        bytes,streams,offs,_,_=cgstream(g,P,p,ctx,index,0)
        chunks=ctx && index ? (streams[4],streams[5],streams[6]) : nothing
        rc,rd,cp=ctx ? streams[1:3] : (UInt8[],UInt8[],UInt8[])
        function run(a,b,c,I=T)
            r=CIO.BitReader(bytes)
            nls=CG.decode_level(r,p;T=T,ctx_range=ctx,coding_scheme=index ? :index : :children,
                resid_bytes=rc,refdist_bytes=rd,copy_bytes=cp,cg_offsets=offs,chunk_offsets=chunks,
                local_type=I,fused_lr=true,reuse_identity=a,sorted_fastpath=b,reuse_scratch=c)
            (nls,r.index,r.bit_count)
        end
        baseline=run(false,false,false)
        @test all(get(baseline[1],T(v),T[])==outneighbors(g,v) for v in 1:nv(g))
        for a in (false,true),b in (false,true),c in (false,true)
            @test run(a,b,c)==baseline
        end
        @test run(true,true,true,Int)==baseline
        # A returned list is never scratch for a later decode or another vertex.
        owned=run(true,true,true)[1]
        if !isempty(owned)
            firstkey=first(keys(owned));empty!(owned[firstkey])
            @test all(v==firstkey || owned[v]==baseline[1][v] for v in keys(owned))
            @test run(true,true,true)==baseline
        end
    end
    for n in (1,2,20), edges in (false,true)
        g=SimpleDiGraph{UInt16}();add_vertices!(g,n)
        if edges
            for u in 1:n,v in 1:n;u!=v && add_edge!(g,u,v);end
        end
        p=cgparams();P=[UInt16.(1:n)]
        bytes,streams,_,_,_=cgstream(g,P,p,true,false,0)
        nls=CG.decode_level(CIO.BitReader(bytes),p;T=UInt16,ctx_range=true,
            resid_bytes=streams[1],refdist_bytes=streams[2],copy_bytes=streams[3],
            compact_lists=true,fused_lr=true,reuse_identity=true,sorted_fastpath=true,reuse_scratch=true)
        @test all(get(nls,UInt16(v),UInt16[])==outneighbors(g,v) for v in 1:n)
    end
end

@testset "Sorted uniqueness and reusable merge" begin
    rng=MersenneTwister(53)
    for T in (UInt16,Adjacently.CustomTypes.UInt24), n in 0:150
        a=T.(rand(rng,1:20,n));expected=sort!(unique(a))
        @test CG._cg_sorted_unique!(copy(a))==expected
        @test CG._cg_sorted_unique!(sort(a))==expected
        scratch=T[]
        for split in unique([0,n,n÷2,min(1,n),max(0,n-1)])
            input=vcat(sort(a[1:split]),sort(a[split+1:end]))
            @test C._merge_lr_runs!(input,split,scratch)==sort(a)
            @test input !== scratch
        end
    end
end

@testset "Identity mapping fallback and undirected bitset" begin
    g=cgfixture(UInt16);p=cgparams();P=[UInt16.(1:120)]
    bytes,streams,_,_,_=cgstream(g,P,p,true,false,0)
    function mapped(fast)
        r=CIO.BitReader(bytes)
        CG._read_membership(r,p,UInt16)
        nls=CG.decode_level(r,p;T=UInt16,ctx_range=true,
            resid_bytes=streams[1],refdist_bytes=streams[2],copy_bytes=streams[3],
            preparsed_clusters=[UInt16.(2:2:240)],compact_lists=true,fused_lr=true,reuse_identity=fast,
            sorted_fastpath=fast,reuse_scratch=fast)
        (nls,r.index,r.bit_count)
    end
    @test mapped(true)==mapped(false)
    @test all(get(mapped(true)[1],UInt16(2v),UInt16[])==2outneighbors(g,v) for v in 1:nv(g))
    u=SimpleGraph{UInt16}();add_vertices!(u,20)
    for v in 2:20;add_edge!(u,1,v);end
    w=CIO.BitWriter();CG.encode_level(w,u,[UInt16.(1:20)];params=p)
    CIO.flush_bitwriter(w;flush_last_bits=true)
    nls=CG.decode_level(CIO.BitReader(collect(CIO.get_bytes(w))),p;T=UInt16,directed=false,
        reuse_identity=true,sorted_fastpath=true,reuse_scratch=true)
    @test all(get(nls,UInt16(v),UInt16[])==outneighbors(u,v) for v in 1:nv(u))
end
