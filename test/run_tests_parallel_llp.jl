using Test,Random,Logging,LightGraphs,Adjacently,JSON
global_logger(ConsoleLogger(stderr,Logging.Warn))
const R=Adjacently.Relabeling
function fixture_graph()
    g=SimpleDiGraph{UInt16}(120);rng=MersenneTwister(83)
    for u in 1:120,v in 1:120
        if u!=v && ((u-1)÷12==(v-1)÷12 ? rand(rng)<.3 : rand(rng)<.004)
            add_edge!(g,UInt16(u),UInt16(v))
        end
    end
    return g
end
@testset "Legacy cluster order and RNG consumption preserved" begin
    g=fixture_graph();part=[(v-1)÷12+1 for v in 1:120]
    Random.seed!(391)
    order,_=R._leiden_llp_order(g,part;llp_mode=:sym,llp_passes=3,sort_clusters=:size_desc)
    tail=string.(rand(UInt64,4))
    path=joinpath(@__DIR__,"fixtures/codec_legacy_v1.json")
    @test isfile(path)
    if isfile(path)
        frozen=JSON.parsefile(path)
        @test Int.(order)==frozen["cluster_order"]
        @test tail==frozen["rng_tail"]
    end
end
@testset "Schedule-independent cluster LLP" begin
    g=fixture_graph();part=[(v-1)÷12+1 for v in 1:120]
    before=[copy(outneighbors(g,v)) for v in 1:nv(g)];saved=copy(part)
    for mode in (:sym,:out),sort_clusters in (:size_desc,:none),seed in (0,53)
        serial=R._leiden_llp_order(g,part;llp_mode=mode,llp_passes=3,sort_clusters,cluster_seed=seed)
        for workers in (1,2,4,100),repeat in 1:2
            parallel=R._leiden_llp_order(g,part;llp_mode=mode,llp_passes=3,sort_clusters,
                cluster_seed=seed,parallel_clusters=true,cluster_workers=workers)
            @test serial==parallel
        end
        @test sort(serial[1])==collect(UInt16(1):UInt16(nv(g)))
    end
    @test part==saved
    @test all(outneighbors(g,v)==before[v] for v in 1:nv(g))
    @test_throws ArgumentError R.relabel_graph_leiden_llp(g;parallel_clusters=true)
    @test_throws ArgumentError R.relabel_graph_leiden_llp(g;parallel_clusters=true,cluster_seed=53,cluster_workers=0)
    for merge in (nothing,10,:auto)
        Random.seed!(912)
        a=R.relabel_graph_leiden_llp(g;llp_passes=2,cluster_seed=53,return_clusters=true,merge_clusters=merge)
        Random.seed!(912)
        b=R.relabel_graph_leiden_llp(g;llp_passes=2,cluster_seed=53,parallel_clusters=true,return_clusters=true,merge_clusters=merge)
        @test a[2:3]==b[2:3]
        @test all(outneighbors(a[1],v)==outneighbors(b[1],v) for v in 1:nv(g))
    end
    for n in (0,1,2,7)
        h=SimpleDiGraph{UInt16}(n);p=collect(1:n)
        a=R._leiden_llp_order(h,p;llp_mode=:sym,llp_passes=2,sort_clusters=:size_desc,cluster_seed=53)
        @test a==R._leiden_llp_order(h,p;llp_mode=:sym,llp_passes=2,sort_clusters=:size_desc,cluster_seed=53,parallel_clusters=true)
    end
end
