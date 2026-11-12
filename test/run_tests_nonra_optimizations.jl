using Test, Adjacently, LightGraphs, Random, SHA, JSON, Logging
const M = Adjacently.MGS
const C = Adjacently.Compression
const CG = C.CG
const AIO = Adjacently.IO
global_logger(ConsoleLogger(stderr, Logging.Warn))

function params(; kw...)
    config = JSON.parsefile(joinpath(@__DIR__, "../bench/graph_compression/configs/cg_speed_v1.json"))["encoder"]
    CG.CGParams(; merge(Dict(Symbol(k) => (v isa String ? Symbol(v) : v) for (k,v) in config), Dict(kw))...)
end
function graph_case(n, pattern)
    g = SimpleDiGraph{UInt32}(); add_vertices!(g, n)
    if pattern == :chain
        for v in 1:n; add_edge!(g,v,min(n,v+1)); end
    elseif pattern == :overlap
        for v in 1:n, u in max(1,v-7):min(n,v+7); add_edge!(g,v,u); end
    elseif pattern == :random
        rng = MersenneTwister(n)
        for _ in 1:8n; add_edge!(g,rand(rng,1:n),rand(rng,1:n)); end
    end
    g
end
function write_case(g,base,alg,backend; k=1, optimized=false, workers=1,
        mode=:children, lr=true, tight=false, greedy=false, model=0)
    if alg == :bg
        M.write_bg_mgs3_graph(g,base;integer_encoding=backend,coding_scheme=mode,
            lr_split=lr,multi_ref=true,index_sample_k=mode==:index ? 16 : 0,cost_model=model)
    elseif alg == :cs
        M.write_cs_mgs3_graph(g,base;integer_encoding=backend,coding_scheme=mode,
            lr_split=lr,index_sample_k=mode==:index ? 16 : 0,cost_model=model)
    else
        p=params(membership=k==1 ? :implicit_ranges : :stop,
            optimize_costs=optimized,intra_tight_deltas=tight,intra_lr_split=lr,
            intra_greedy_mil=greedy,cost_model=model)
        P=k==1 ? [UInt32.(1:nv(g))] : [UInt32.(i:k:nv(g)) for i in 1:k]
        M.write_cg_mgs3_graph(g,base,P;params=p,integer_encoding=backend,
            coding_scheme=mode,parallel_search=workers>1,search_workers=workers)
    end
    base*".mgz"
end
function check_owned(a,g)
    @test a == [outneighbors(g,v) for v in 1:nv(g)]
    @test length(Set(objectid(x) for x in a)) == length(a)
    b=deepcopy(a); empty!(a[1])
    @test a[2:end] == b[2:end]
end

@testset "Opt-in full adjacency and staged reconstruction" begin
    @test !CG.CGParams().optimize_costs
    loader_source=read(joinpath(@__DIR__,"../src/mgs_adjacency.jl"),String)
    @test occursin("optimize_decode::Bool=false",loader_source)
    @test occursin("parallel_reconstruct::Bool=false",loader_source)
    decoder_source=read(joinpath(@__DIR__,"../src/compression/cge.jl"),String)
    @test occursin("skip_normalize::Bool=false, direct_output::Bool=false",decoder_source)
    mktempdir() do dir
        for n in (1,2,65,257), pattern in (:empty,:chain,:overlap,:random),
                alg in (:bg,:cs,:cg), backend in (:context_range,:fibonacci)
            g=graph_case(n,pattern); path=write_case(g,joinpath(dir,"case"),alg,backend)
            before=read(path)
            legacy=M.load_adjacency_mgs3_graph(path)
            @test legacy == M.load_adjacency_mgs3_graph(path;optimize_decode=false,parallel_reconstruct=false)
            for compact in (false,true)
                a=M.load_adjacency_mgs3_graph(path;optimize_decode=true,compact_output=compact)
                @test a == legacy
                check_owned(a,g)
            end
            if backend == :context_range && alg != :bg
                for workers in (1,8)
                    a=M.load_adjacency_mgs3_graph(path;parallel_reconstruct=true,decode_workers=workers)
                    @test a == legacy
                    check_owned(a,g)
                end
            else
                @test_throws ArgumentError M.load_adjacency_mgs3_graph(path;parallel_reconstruct=true)
            end
            @test read(path) == before
        end
        g=graph_case(129,:overlap)
        for alg in (:bg,:cs,:cg)
            path=write_case(g,joinpath(dir,"indexed"),alg,:context_range;mode=:index)
            @test_throws ArgumentError M.load_adjacency_mgs3_graph(path;optimize_decode=true)
            @test_throws ArgumentError M.load_adjacency_mgs3_graph(path;parallel_reconstruct=true)
            path=write_case(g,joinpath(dir,"plain"),alg,:context_range;lr=false)
            @test M.load_adjacency_mgs3_graph(path;optimize_decode=true) == M.load_adjacency_mgs3_graph(path)
            @test_throws ArgumentError M.load_adjacency_mgs3_graph(path;parallel_reconstruct=true)
        end
        path=write_case(g,joinpath(dir,"clusters"),:cg,:context_range;k=2)
        check_owned(M.load_adjacency_mgs3_graph(path;optimize_decode=true),g)
        @test_throws ArgumentError M.load_adjacency_mgs3_graph(path;parallel_reconstruct=true)
        @test_throws ArgumentError M.load_adjacency_mgs3_graph(path;parallel_reconstruct=true,decode_workers=0)
        for tight in (true,)
            path=write_case(g,joinpath(dir,"fallback"),:cg,:context_range;tight)
            check_owned(M.load_adjacency_mgs3_graph(path;optimize_decode=true),g)
            @test_throws ArgumentError M.load_adjacency_mgs3_graph(path;parallel_reconstruct=true)
        end
    end
end

@testset "Exact CG cost scans and unchanged encoding decisions" begin
    rng=MersenneTwister(53)
    for T in (UInt8,UInt16,Adjacently.CustomTypes.UInt24,UInt32), trial in 1:500,
            mil in (2,3,4,5,8), tight in (false,true)
        upper=T===UInt8 ? 200 : 2000
        values=sort!(unique(T.(rand(rng,1:upper,rand(rng,0:80)))))
        if trial%3==0; append!(values,T.(20:40)); sort!(unique!(values)); end
        vid=T(rand(rng,1:upper))
        for enc in (:fibonacci,:elias_delta,:elias_gamma)
            @test CG._estimate_ir_lr_cost(values,enc,mil,vid;tight_deltas=tight,optimize=true) ==
                CG._estimate_ir_lr_cost(values,enc,mil,vid;tight_deltas=tight)
        end
        len=rand(rng,0:upper)
        pos=sort!(unique(T.(rand(rng,1:max(1,len),len==0 ? 0 : rand(rng,0:80)))))
        @test CG._estimate_complement_blocks_cost(pos,len,:fibonacci) ==
            CG._estimate_copy_blocks_cost(setdiff(T.(1:len),pos),:fibonacci)
    end
    for T in (UInt16,Adjacently.CustomTypes.UInt24), values in (T[],T[1],T[1,2,3,4],T[1,3,5,7]), vid in (nothing,T(1),T(4),T(100))
        @test CG._estimate_ir_lr_cost(values,:fibonacci,4,vid;optimize=true) ==
            CG._estimate_ir_lr_cost(values,:fibonacci,4,vid)
    end
    mktempdir() do dir
        g=graph_case(129,:overlap)
        for k in (1,2), backend in (:fibonacci,:context_range), mode in (:children,:index),
                model in (0,1), tight in (false,true)
            path=write_case(g,joinpath(dir,"old"),:cg,backend;k,mode,model,tight)
            original=read(path)
            for workers in (1,8)
                path=write_case(g,joinpath(dir,"new"),:cg,backend;k,mode,model,tight,optimized=true,workers)
                @test read(path)==original
                if mode==:children
                    check_owned(M.load_adjacency_mgs3_graph(path;optimize_decode=true),g)
                else
                    decoded=M.load_compressed_mgs3_graph(path)
                    @test all(outneighbors(decoded,v)==outneighbors(g,v) for v in 1:nv(g))
                end
            end
        end
    end
    # Greedy MIL is a raw CG grammar, not a field in the MGS byte2 header.
    # Supply the same explicit parameters to its reader.
    g=graph_case(129,:overlap); P=[UInt32.(1:nv(g))]
    for tight in (false,true), model in (0,1)
        p=params(intra_greedy_mil=true,intra_tight_deltas=tight,cost_model=model)
        w=AIO.BitWriter(); streams=CG.encode_level(w,g,P;params=p,ctx_range=true)
        AIO.flush_bitwriter(w;flush_last_bits=true); original=collect(AIO.get_bytes(w))
        for workers in (1,8)
            q=params(intra_greedy_mil=true,intra_tight_deltas=tight,cost_model=model,optimize_costs=true)
            w=AIO.BitWriter(); out=CG.encode_level(w,g,P;params=q,ctx_range=true,
                parallel_search=workers>1,search_workers=workers)
            AIO.flush_bitwriter(w;flush_last_bits=true); structural=collect(AIO.get_bytes(w))
            @test structural==original
            @test out==streams
            rows=CG.decode_level(AIO.BitReader(structural),p;T=UInt32,ctx_range=true,
                resid_bytes=out[1],refdist_bytes=out[2],copy_bytes=out[3],
                compact_lists=true,fused_lr=true,reuse_identity=true,reuse_scratch=true,
                sorted_fastpath=true,skip_normalize=true,direct_output=true)
            check_owned(M._own_mgs3_adjacency(rows,Int(nv(g)),false,true),g)
        end
    end
end

@testset "Independent entropy sources across calls" begin
    mktempdir() do dir
        g=graph_case(257,:overlap)
        paths=[write_case(g,joinpath(dir,string(alg)),alg,:context_range) for alg in (:bg,:cs,:cg)]
        resid=C.CtxRangeDecoder(UInt8[]); flag=C.BinRangeDecoder(UInt8[])
        oldr=C._RESID_SOURCE[]; oldf=C._FLAG_SOURCE[]
        try
            C._RESID_SOURCE[]=resid; C._FLAG_SOURCE[]=flag
            jobs=[Threads.@spawn M.load_adjacency_mgs3_graph(paths[mod1(i,3)];optimize_decode=true) for i in 1:24]
            for job in jobs; check_owned(fetch(job),g); end
            jobs=[Threads.@spawn M.load_adjacency_mgs3_graph(paths[mod1(i,2)+1];parallel_reconstruct=true,decode_workers=2) for i in 1:8]
            for job in jobs; check_owned(fetch(job),g); end
            @test C._RESID_SOURCE[]===resid
            @test C._FLAG_SOURCE[]===flag
        finally
            C._RESID_SOURCE[]=oldr; C._FLAG_SOURCE[]=oldf
        end
    end
end

@testset "Absolute and relative hybrid residual coordinates" begin
    rng=MersenneTwister(79); absolute=0
    patterns=[UInt32[7],UInt32.(3:60),UInt32.(3:5:303),
        UInt32[2,5,9,10,11,12,13,14,30,35,40,45,50,51,52,53,54,55,80]]
    append!(patterns,[sort!(unique(rand(rng,UInt32(1):UInt32(600),n))) for n in (2,8,64,256) for _ in 1:10])
    for ns in patterns, vid in (UInt32(1),UInt32(200),UInt32(700)), hybrid in (false,true), mil in (3,5)
        rc=C.CtxRangeEncoder(); flag=C.BinRangeEncoder(); w=AIO.BitWriter()
        oldr=C._RESID_SINK[]; oldf=C._FLAG_SINK[]
        try
            C._RESID_SINK[]=rc; C._FLAG_SINK[]=flag
            C.write_hybrid_mix_encoded_list(w,C.delta_encode_vector(ns),:context_range,hybrid,mil,false;vertex_id=vid)
        finally
            C._RESID_SINK[]=oldr; C._FLAG_SINK[]=oldf
        end
        rb=C.rc_finish!(rc); fb=C.brc_finish!(flag)
        coords=M._nonra_cs_hybrid!(UInt64[],length(ns),C.CtxRangeDecoder(rb),C.BinRangeDecoder(fb))
        anchor=M._NONRA_ANCHOR
        @test [x>anchor÷2 ? Int64(x)-Int64(anchor)+Int64(vid) : Int64(x) for x in coords]==ns
        absolute+=any(x->x<anchor÷2,coords)
        for T in (UInt16,Adjacently.CustomTypes.UInt24)
            @test C._read_nonra_hybrid(AIO.BitReader(UInt8[]),:index,:context_range,T,
                C.CtxRangeDecoder(rb),C.BinRangeDecoder(fb);max_elements=length(ns),vertex_id=T(vid))==ns
        end
    end
    @test absolute>0
end
