using Test, Adjacently, LightGraphs, JSON, SHA, Logging
const M = Adjacently.MGS
const CG = Adjacently.Compression.CG
const OFF = (;compact_lists=false,fused_lr=false,reuse_identity=false,
    sorted_fastpath=false,reuse_scratch=false)
const ON = (;compact_lists=true,fused_lr=true,reuse_identity=true,
    sorted_fastpath=true,reuse_scratch=true)
global_logger(ConsoleLogger(stderr, Logging.Warn))

function fixture(n, pattern)
    g=SimpleDiGraph{UInt16}();add_vertices!(g,n)
    if pattern == :overlap
        for v in 1:n, u in max(1,v-8):min(n,v+8)
            add_edge!(g,v,u) # Includes self-loops.
        end
    elseif pattern == :chain
        for v in 1:n-1;add_edge!(g,v,v+1);end
    end
    g
end

function cgparams(;kw...)
    cfg=JSON.parsefile(joinpath(@__DIR__,"../bench/graph_compression/configs/cg_speed_v1.json"))["encoder"]
    CG.CGParams(;merge(Dict(Symbol(k)=>(v isa String ? Symbol(v) : v) for (k,v) in cfg),Dict(kw))...)
end

function write_fixture(g, base, alg, backend; mode=:children)
    if alg == :bg
        M.write_bg_mgs3_graph(g,base;integer_encoding=backend,coding_scheme=mode,lr_split=true,multi_ref=true,
            index_sample_k=mode==:index ? 16 : 0)
    elseif alg == :cs
        M.write_cs_mgs3_graph(g,base;integer_encoding=backend,coding_scheme=mode,lr_split=true,
            index_sample_k=mode==:index ? 16 : 0)
    else
        n=Int(nv(g));T=eltype(g)
        # Non-contiguous K=2 membership exercises permutation/inter-cluster arcs.
        P=alg == :cg1 ? [T.(1:n)] : [T.(1:2:n),T.(2:2:n)]
        M.write_cg_mgs3_graph(g,base,P;params=cgparams(membership=alg==:cg1 ? :implicit_ranges : :stop),integer_encoding=backend,coding_scheme=mode)
    end
    base*".mgz"
end

@testset "Owned non-RA adjacency, unchanged files and legacy API" begin
    mktempdir() do dir
        for n in (2,129,257), pattern in (:empty,:chain,:overlap),
                alg in (:bg,:cs,:cg1,:cg2), backend in (:fibonacci,:context_range)
            g=fixture(n,pattern)
            path=write_fixture(g,joinpath(dir,"case"),alg,backend);before=read(path)
            options=alg in (:cg1,:cg2) ? (OFF,ON) : (OFF,)
            a=M.load_adjacency_mgs3_graph(path)
            legacy=M.load_compressed_mgs3_graph(path)
            @test a isa Vector{<:Vector{<:Unsigned}}
            @test length(a)==n
            @test all(a[v]==outneighbors(g,v)==outneighbors(legacy,v) for v in 1:n)
            for flags in options, compact in (false,true)
                b=M.load_adjacency_mgs3_graph(path;compact_output=compact,flags...)
                @test b==a
                empty!(b[1]);push!(b[1],eltype(b[1])(1))
                @test all(b[v]==a[v] for v in 2:n)
                @test M.load_adjacency_mgs3_graph(path;compact_output=compact,flags...)==a
            end
            @test read(path)==before
            if alg in (:bg,:cs)
                for flag in keys(OFF)
                    @test_throws ArgumentError M.load_adjacency_mgs3_graph(path;NamedTuple{(flag,)}((true,))...)
                end
            end
        end
    end
end

@testset "Default algorithm IDs" begin
    mktempdir() do dir
        g=fixture(129,:overlap)
        for alg in (:bg,:cs,:cg),backend in (:fibonacci,:context_range)
            base=joinpath(dir,"default");p=alg==:bg ? M._bg_default_params() : alg==:cs ? M._cs_default_params() : M._cg_default_params()
            if alg==:bg
                M.write_bg_mgs3_graph(g,base;integer_encoding=backend,
                    ref_window_size=p.ref_window_size,copy_blocks=p.copy_blocks,
                    stop_deltas=p.stop_deltas,lr_split=p.lr_split,adaptive_header=p.adaptive_header)
            elseif alg==:cs
                M.write_cs_mgs3_graph(g,base;integer_encoding=backend,
                    ref_window_size=p.ref_window_size,lr_split=p.lr_split)
            else
                M.write_cg_mgs3_graph(g,base,[UInt16.(1:129)];integer_encoding=backend,params=p)
            end
            bytes=read(base*".mgz");bytes[7]=alg==:bg ? M.ALG_BG : alg==:cs ? M.ALG_CS : M.ALG_CG
            open(io->write(io,bytes),base*".mgz","w")
            a=M.load_adjacency_mgs3_graph(base*".mgz")
            @test all(a[v]==outneighbors(g,v) for v in 1:129)
        end
    end
end

@testset "Unsupported headers and truncated containers" begin
    mktempdir() do dir
        g=fixture(129,:chain);path=write_fixture(g,joinpath(dir,"valid"),:bg,:context_range)
        original=read(path);bad=joinpath(dir,"bad.mgz")
        for len in (0,2,11,12,20,51)
            open(io->write(io,original[1:len]),bad,"w")
            @test_throws EOFError M.load_adjacency_mgs3_graph(bad)
        end
        for (pos,value) in ((1,0),(4,4),(5,2),(6,0x47),(6,0x17),(6,0x27),(6,0x03),(7,0x01),(7,0x05))
            b=copy(original);b[pos]=value;open(io->write(io,b),bad,"w")
            @test_throws ArgumentError M.load_adjacency_mgs3_graph(bad)
        end
        b=copy(original);b[8:12].=0;open(io->write(io,b),bad,"w")
        @test_throws ArgumentError M.load_adjacency_mgs3_graph(bad)
        b=copy(original);b[13:20].=0xff;open(io->write(io,b),bad,"w")
        @test_throws EOFError M.load_adjacency_mgs3_graph(bad)
        for alg in (:bg,:cs,:cg1), backend in (:fibonacci,:context_range)
            path=write_fixture(g,joinpath(dir,"indexed"),alg,backend;mode=:index)
            @test_throws ArgumentError M.load_adjacency_mgs3_graph(path)
        end
    end
end

@testset "Typed materialization" begin
    for T in (UInt8,UInt16,Adjacently.CustomTypes.UInt24)
        rows=Dict(T(1)=>T[3,2,2],T(2)=>T[1])
        for compact in (false,true)
            fresh=deepcopy(rows);a=M._own_mgs3_adjacency(fresh,3,compact)
            @test a==[T[2,3],T[1],T[]]
            @test a[1] !== a[2] !== a[3]
        end
    end
end

@testset "Exact-power vertex universes" begin
    mktempdir() do dir
        for n in (256,65536)
            # CG can emit these universes directly. BG/CS writers have their
            # own legacy ceil(log2(n)) boundary limitation, outside this API.
            g=SimpleDiGraph{UInt32}();add_vertices!(g,n)
            path=write_fixture(g,joinpath(dir,"power"),:cg1,:fibonacci)
            adj=M.load_adjacency_mgs3_graph(path)
            @test length(adj)==n && all(isempty,adj)
            @test eltype(eltype(adj))== (n==256 ? UInt16 : Adjacently.CustomTypes.UInt24)
        end
    end
end

@testset "Frozen EAT CG file and default-off contract" begin
    mktempdir() do dir
        g=Adjacently.Graph.get_core(Adjacently.IO.load_graph_from_pajek(joinpath(@__DIR__,"../datasets/EAT/EATnew.net")))[1]
        base=joinpath(dir,"eat");M.write_cg_mgs3_graph(g,base,[eltype(g).(1:nv(g))];params=cgparams(),integer_encoding=:context_range)
        path=base*".mgz";hash=bytes2hex(sha256(read(path)))
        @test hash=="bd4c5fe687f60bff97262527480a73504957795a90c1ad0408611eb44cf35f64"
        for flags in (OFF,ON), compact in (false,true)
            a=M.load_adjacency_mgs3_graph(path;flags...,compact_output=compact)
            @test all(a[v]==outneighbors(g,v) for v in 1:nv(g))
        end
        @test bytes2hex(sha256(read(path)))==hash
    end
    source=read(joinpath(@__DIR__,"../src/mgs_adjacency.jl"),String)
    for flag in (keys(OFF)...,:compact_output)
        @test occursin(string(flag,"::Bool=false"),source)
        @test !occursin(string(flag,"::Bool=true"),source)
    end
end
