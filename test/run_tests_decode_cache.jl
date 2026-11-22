using Test, Random, Logging, Adjacently, LightGraphs
const CM=Adjacently.MGS
global_logger(ConsoleLogger(stderr,Logging.Warn))

@testset "Byte-budgeted LRU accounting" begin
    block(i)=CM._CachedAdjacencyChunk(i,Int[1,4],UInt16[2,3,4])
    a,b,c=block(1),block(2),block(3)
    entries=Dict(1=>CM._DecodedCacheEntry(a,UInt64(1)),2=>CM._DecodedCacheEntry(b,UInt64(2)))
    cache=CM._DecodedChunkCache(UInt16,Base.summarysize(entries))
    CM._decoded_cache_admit!(cache,1,a); CM._decoded_cache_admit!(cache,2,b)
    @test Set(keys(cache.entries)) == Set([1,2])
    @test CM._decoded_cache_get(cache,1) === a
    CM._decoded_cache_admit!(cache,3,c)
    @test Set(keys(cache.entries)) == Set([1,3])
    @test cache.evictions == 1
    @test cache.resident == Base.summarysize(cache.entries) <= cache.budget
    huge=CM._CachedAdjacencyChunk(4,Int[1,10001],ones(UInt16,10000))
    CM._decoded_cache_admit!(cache,4,huge)
    @test Set(keys(cache.entries)) == Set([1,3])
    @test cache.bypasses == 1
    for _ in 1:1000; CM._decoded_cache_get(cache,1); end
    @test cache.resident == Base.summarysize(cache.entries) <= cache.budget
    emptyblock=CM._CachedAdjacencyChunk(5,Int[1,1],UInt16[])
    @test isempty(CM._cached_neighbors(emptyblock,5))
end

@testset "Budgeted graph queries" begin
    g=Adjacently.SimpleDiGraph{UInt16}(); add_vertices!(g,192)
    for u in 1:184, v in 1:192
        if u%11 != 0 && (20<=v<=45 || 110<=v<=130 || v==mod(u+60,192)+1)
            u != v && add_edge!(g,u,v)
        end
    end
    mktempdir() do dir
        for alg in (:bg,:cs), localrefs in (false,true)
            writer=alg==:bg ? CM.write_bg_mgs3_graph : CM.write_cs_mgs3_graph
            opts=alg==:bg ? (;multi_ref=true) : NamedTuple()
            base=joinpath(dir,"$(alg)_$(localrefs)")
            writer(g,base;integer_encoding=:context_range,coding_scheme=:index,
                ref_window_size=64,lr_split=true,cost_model=0,index_sample_k=16,
                chunk_local_refs=localrefs,opts...)
            path=base*".mgz"
            @test_throws ArgumentError CM.load_indexed_mgs3_graph(path;cache_budget_bytes=0)
            @test_throws ArgumentError CM.load_indexed_mgs3_graph(path;decode_prefix=false,cache_budget_bytes=-1)
            for budget in (0,1,8192,65536)
                idx=CM.load_indexed_mgs3_graph(path;decode_prefix=false,cache_budget_bytes=budget)
                @test idx.stats.cache_resident_bytes()==0
                for q in randperm(MersenneTwister(62),192)
                    @test idx.neighbors(q)==outneighbors(g,q)
                    @test idx.stats.cache_resident_bytes() <= budget
                end
                @test idx.stats.cache_hits[]+idx.stats.cache_misses[] == 192
                @test idx.stats.cache_peak_bytes() <= budget
                if budget <= 1
                    @test idx.stats.cached_chunks()==0
                    @test idx.stats.cache_hits[]==0
                    @test idx.stats.oversized_bypasses()>0
                end
                idx.reset_fn()
                @test idx.stats.cached_chunks()==0 && idx.stats.cache_resident_bytes()==0
                @test idx.stats.cache_hits[]==0 && idx.stats.cache_misses[]==0
                @test idx.stats.decode_calls[]==0
                @test_throws ErrorException idx.neighbors(0)
                @test_throws ErrorException idx.neighbors(193)
                if budget==65536
                    got=idx.neighbors(1); empty!(got)
                    @test idx.neighbors(1)==outneighbors(g,1)
                    before=idx.stats.vertices_decoded[]
                    @test idx.neighbors(16)==outneighbors(g,16)
                    @test idx.stats.vertices_decoded[]==before
                    @test idx.stats.cache_hits[]==2
                end
            end
        end
        base=joinpath(dir,"not_indexed")
        CM.write_bg_mgs3_graph(g,base;integer_encoding=:fibonacci,coding_scheme=:children)
        @test_throws ArgumentError CM.load_indexed_mgs3_graph(base*".mgz";decode_prefix=false,cache_budget_bytes=1024)
    end
end
