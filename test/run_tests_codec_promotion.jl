using SHA, JSON

const DECODE_OFF=(;compact_lists=false,fused_lr=false,reuse_identity=false,
    sorted_fastpath=false,reuse_scratch=false)
const DECODE_ON=(;compact_lists=true,fused_lr=true,reuse_identity=true,
    sorted_fastpath=true,reuse_scratch=true)

@testset "Public CG options and threaded decode preserve exact output" begin
    g=cgfixture(UInt16);tmp=mktempdir()
    for k in (1,4), mode in (:children,:index), backend in (:fibonacci,:context_range)
        P=[UInt16.(a:min(a+120÷k-1,120)) for a in 1:120÷k:120]
        base=joinpath(tmp,"public_$(k)_$(mode)_$(backend)")
        CM.write_cg_mgs3_graph(g,base,P;params=cgparams(),coding_scheme=mode,integer_encoding=backend)
        path=base*".mgz";before=read(path)
        baseline=CM.load_compressed_mgs3_graph(path)
        off=CM.load_compressed_mgs3_graph(path;DECODE_OFF...,parallel_decode=false)
        @test all(outneighbors(off,v)==outneighbors(baseline,v)==outneighbors(g,v) for v in 1:nv(g))
        options=[DECODE_ON]
        for name in keys(DECODE_OFF)
            push!(options,merge(DECODE_OFF,NamedTuple{(name,)}((true,))))
        end
        for opts in options
            decoded=CM.load_compressed_mgs3_graph(path;opts...)
            @test all(outneighbors(decoded,v)==outneighbors(g,v) for v in 1:nv(g))
            idx=CM.load_indexed_mgs3_graph(path;opts...)
            for v in (1,17,120)
                if hasproperty(idx,:reset_fn)
                    idx.reset_fn()
                end
                @test idx.neighbors(v)==outneighbors(g,v)
                if mode==:index && backend==:context_range
                    answer=idx.neighbors(v);empty!(answer)
                    @test idx.neighbors(v)==outneighbors(g,v)
                end
            end
            if mode==:index && backend==:context_range
                for workers in (1,4)
                    par=CM.load_compressed_mgs3_graph(path;opts...,parallel_decode=true,decode_workers=workers)
                    @test all(outneighbors(par,v)==outneighbors(g,v) for v in 1:nv(g))
                end
            end
        end
        @test read(path)==before
    end
    CM.write_bg_mgs3_graph(g,joinpath(tmp,"bg"))
    @test_throws ArgumentError CM.load_compressed_mgs3_graph(joinpath(tmp,"bg.mgz");reuse_identity=true)
    @test_throws ArgumentError CM.load_indexed_mgs3_graph(joinpath(tmp,"bg.mgz");reuse_scratch=true)
end

@testset "Frozen CG EAT files, independent of research checkout" begin
    config=JSON.parsefile(joinpath(@__DIR__,"../bench/graph_compression/configs/cg_speed_v1.json"))
    p=CG.CGParams(;[Symbol(k)=>(v isa String ? Symbol(v) : v) for (k,v) in config["encoder"]]...)
    expected=Dict(
        (:fibonacci,:children)=>"f259b29793ad8f94c63fff6a2852aee71c046f2a8eb1f396a7adc7c56c686c18",
        (:context_range,:children)=>"bd4c5fe687f60bff97262527480a73504957795a90c1ad0408611eb44cf35f64",
        (:context_range,:index)=>"0cc0a1d937b4c8606523dc4e9a39c26346873a1cdaaab5f334c156d38740a6ab")
    g=Adjacently.Graph.get_core(CIO.load_graph_from_pajek(joinpath(@__DIR__,"../datasets/EAT/EATnew.net")))[1]
    P=[eltype(g).(1:nv(g))];tmp=mktempdir()
    for ((backend,mode),digest) in expected, parallel in (false,true)
        base=joinpath(tmp,"golden_$(backend)_$(mode)_$(parallel)")
        CM.write_cg_mgs3_graph(g,base,P;params=p,integer_encoding=backend,coding_scheme=mode,
            parallel_search=parallel,search_workers=4)
        path=base*".mgz"
        @test bytes2hex(sha256(read(path)))==digest
        for opts in (DECODE_OFF,DECODE_ON)
            decoded=CM.load_compressed_mgs3_graph(path;opts...)
            @test all(outneighbors(decoded,v)==outneighbors(g,v) for v in 1:nv(g))
        end
    end
end

@testset "Promotion switches are visibly default-off" begin
    root=normpath(joinpath(@__DIR__,".."))
    cfg=JSON.parsefile(joinpath(root,"bench/graph_compression/configs/codec_acceleration_v1.json"))
    @test cfg["cg_decode"]["recommended"]==Dict(string(k)=>v for (k,v) in pairs(DECODE_ON))
    @test cfg["cg_decode"]["baseline"]==Dict(string(k)=>v for (k,v) in pairs(DECODE_OFF))
    @test all(v===false || (k=="cluster_seed" && v===nothing) for (k,v) in cfg["api_defaults"])
    for file in ("src/mgs.jl","src/compression/cge.jl"),
            flag in ("compact_lists","fused_lr","reuse_identity","sorted_fastpath","reuse_scratch")
        source=read(joinpath(root,file),String)
        @test occursin(Regex(flag*"::Bool=false"),source)
        @test !occursin(Regex(flag*"::Bool=true"),source)
    end
    for (file,flag) in (("src/mgs.jl","parallel_search"),("src/mgs.jl","parallel_decode"),
            ("src/compression.jl","parallel_search"),("src/relabeling.jl","parallel_clusters"))
        source=read(joinpath(root,file),String)
        @test occursin(Regex(flag*"::Bool=false"),source)
        @test !occursin(Regex(flag*"::Bool=true"),source)
    end
end

@testset "Default BG/CS indexed encoding remains unchanged" begin
    g=cgfixture(UInt16);tmp=mktempdir()
    for writer in (CM.write_bg_mgs3_graph,CM.write_cs_mgs3_graph), backend in (:fibonacci,:context_range)
        base=joinpath(tmp,string(nameof(writer))*string(backend))
        writer(g,base;coding_scheme=:index,integer_encoding=backend,index_sample_k=16)
        writer(g,base*"_off";coding_scheme=:index,integer_encoding=backend,index_sample_k=16,parallel_search=false)
        @test read(base*".mgz")==read(base*"_off.mgz")
        decoded=CM.load_compressed_mgs3_graph(base*".mgz")
        @test all(outneighbors(decoded,v)==outneighbors(g,v) for v in 1:nv(g))
    end
end
