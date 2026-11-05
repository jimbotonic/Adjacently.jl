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
    # The frozen fixture predates the first-appearance default.
    prev=R.LLP_GROUP_ORDER[];R.LLP_GROUP_ORDER[]=:size_label
    order,_=try R._leiden_llp_order(g,part;llp_mode=:sym,llp_passes=3,sort_clusters=:size_desc)
    finally R.LLP_GROUP_ORDER[]=prev end
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

@testset "LLP group order" begin
    @test R.LLP_GROUP_ORDER[]===:first_appearance
    g=fixture_graph()
    for mode in (:sym,:out)
        a,la=R.relabel_vertices_llp(g,mode;passes=3,return_labels=true,rng=MersenneTwister(7))
        b,lb=R.relabel_vertices_llp(g,mode;passes=3,return_labels=true,rng=MersenneTwister(7),group_order=:size_label)
        @test la==lb                      # labels do not depend on the composition rule
        @test sort(collect(values(a)))==collect(UInt16(1):UInt16(nv(g)))
        @test sort(collect(values(b)))==collect(UInt16(1):UInt16(nv(g)))
    end
    @test_throws ArgumentError R.relabel_vertices_llp(g;group_order=:size)
    # Two components, the smaller one holding vertex 1: first appearance keeps it
    # first and contiguous; the legacy rule moves the larger component ahead.
    h=SimpleDiGraph{UInt16}(12)
    for (lo,hi) in ((1,3),(4,12)),u in lo:hi,v in lo:hi
        u!=v && add_edge!(h,UInt16(u),UInt16(v))
    end
    fa=R.relabel_vertices_llp(h;passes=3,rng=MersenneTwister(1))
    sl=R.relabel_vertices_llp(h;passes=3,rng=MersenneTwister(1),group_order=:size_label)
    @test sort([Int(fa[UInt16(v)]) for v in 1:3])==[1,2,3]
    @test sort([Int(sl[UInt16(v)]) for v in 1:3])==[10,11,12]
end

@testset "Leiden+LLP cluster order" begin
    @test R.LEIDEN_SORT_CLUSTERS[]===:none
    g=fixture_graph()
    @test_throws ArgumentError R.relabel_graph_leiden_llp(g;sort_clusters=:size)
    Random.seed!(5); a=R.relabel_graph_leiden_llp(g;llp_passes=2)
    Random.seed!(5); b=R.relabel_graph_leiden_llp(g;llp_passes=2,sort_clusters=:none)
    @test a[2]==b[2]
    # Under :none the cluster holding vertex 1 starts the order.
    part=[(v-1)÷12+1 for v in 1:120]; part[1:12].=99   # cluster of vertex 1 has the largest label
    order,_=R._leiden_llp_order(g,part;llp_mode=:sym,llp_passes=2,sort_clusters=:none,cluster_seed=3)
    @test Set(Int.(order[1:12]))==Set(1:12)
end

@testset "Guarded LLP" begin
    # BVGraph R=0 row cost: γ(5)+γ(0)+γ(1 interval)+[γ(0)+γ(nat 1)]+ζ₃(nat 10) = 5+1+3+4+8
    @test R._r0_row_bits(0,[1,2,3,4,10])==21
    @test R._r0_row_bits(3,Int[])==1
    g=fixture_graph(); n=Int(nv(g)); ctx=R._R0Ctx(g)
    b0=R._r0_bytes(ctx,collect(1:n))
    @test_throws ArgumentError R.relabel_vertices_llp(g;guard=:all)
    @test_throws ArgumentError R.relabel_vertices_llp(g;group_order=:size_label,guard=:block)
    Random.seed!(2); a=R.relabel_vertices_llp(g;passes=2)
    Random.seed!(2); b=R.relabel_vertices_llp(g;passes=2,guard=:none)
    @test a==b
    for guard in (:stage,:block), seed in 1:3
        Random.seed!(seed); m=R.relabel_vertices_llp(g;passes=2,guard=guard)
        order=zeros(Int,n); for (o,nw) in m; order[Int(nw)]=Int(o); end
        @test sort(order)==collect(1:n)
        @test R._r0_bytes(ctx,order)<=b0
    end
end

@testset "Guarded Leiden+LLP" begin
    g=fixture_graph(); n=Int(nv(g)); ctx=R._R0Ctx(g)
    part=[(v-1)÷12+1 for v in 1:n]; part[[3,50,77]].=1
    # the unguarded composite chain reproduces the per-cluster order, both RNG policies
    for cs in (nothing,7)
        Random.seed!(4); a,_=R._leiden_llp_order(g,part;llp_mode=:sym,llp_passes=2,sort_clusters=:none,cluster_seed=cs)
        Random.seed!(4); b,_=R._leiden_llp_order(g,part;llp_mode=:sym,llp_passes=2,sort_clusters=:none,cluster_seed=cs,guard=:chain)
        @test a==b
    end
    @test_throws ArgumentError R.relabel_graph_leiden_llp(g;guard=:chain)
    @test_throws ArgumentError R._leiden_llp_order(g,part;llp_mode=:sym,llp_passes=2,sort_clusters=:size_desc,guard=:block)
    b0=R._r0_bytes(ctx,collect(1:n))
    for guard in (:stage,:block), seed in 1:3
        Random.seed!(seed); h,m=R.relabel_graph_leiden_llp(g;llp_passes=2,guard=guard)
        order=zeros(Int,n); for (o,nw) in m; order[Int(nw)]=Int(o); end
        @test sort(order)==collect(1:n)
        @test R._r0_bytes(ctx,order)<=b0
    end
end

@testset "Median group order" begin
    g=fixture_graph(); n=Int(nv(g)); ctx=R._R0Ctx(g)
    # lab by vertex: class 1 = {5,1,4} at positions 1,2,5 (median 2); class 2 = {2,3} at 3,4 (median 3)
    @test R._stable_group([5,1,2,3,4],[1,2,2,1,1],:median)==[5,1,4,2,3]
    # an early outlier leader does not pull its class forward
    @test R._stable_group([1,2,3,4,5,6],[1,2,2,1,1,1],:median)==[2,3,1,4,5,6]
    @test_throws ArgumentError R.relabel_vertices_llp(g;group_order=:mean)
    part=[(v-1)÷12+1 for v in 1:n]; part[[3,50,77]].=1
    old=R.LLP_GROUP_ORDER[]
    try
        R.LLP_GROUP_ORDER[]=:median
        for cs in (nothing,7)
            Random.seed!(4); a,_=R._leiden_llp_order(g,part;llp_mode=:sym,llp_passes=2,sort_clusters=:median,cluster_seed=cs)
            Random.seed!(4); b,_=R._leiden_llp_order(g,part;llp_mode=:sym,llp_passes=2,sort_clusters=:median,cluster_seed=cs,guard=:chain)
            @test a==b
        end
        b0=R._r0_bytes(ctx,collect(1:n))
        for seed in 1:3
            Random.seed!(seed); m=R.relabel_vertices_llp(g;passes=2,guard=:block)
            order=zeros(Int,n); for (o,nw) in m; order[Int(nw)]=Int(o); end
            @test R._r0_bytes(ctx,order)<=b0
            Random.seed!(seed); _,m2=R.relabel_graph_leiden_llp(g;llp_passes=2,sort_clusters=:median,guard=:block)
            order=zeros(Int,n); for (o,nw) in m2; order[Int(nw)]=Int(o); end
            @test sort(order)==collect(1:n) && R._r0_bytes(ctx,order)<=b0
        end
    finally
        R.LLP_GROUP_ORDER[]=old
    end
end
