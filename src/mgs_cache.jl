# Optional whole-chunk decoded cache. Included inside MGS.
# Compressed input, live decode work and returned copies are outside the retained
# cache budget. Nonempty cache residency includes its dictionaries and arrays.
struct _CachedAdjacencyChunk{T<:Unsigned}
    first_vertex::Int
    offsets::Vector{Int}
    targets::Vector{T}
end
mutable struct _DecodedCacheEntry{T<:Unsigned}
    chunk::_CachedAdjacencyChunk{T}
    used::UInt64
end
mutable struct _DecodedChunkCache{T<:Unsigned}
    entries::Dict{Int,_DecodedCacheEntry{T}}
    budget::Int
    clock::UInt64
    resident::Int
    peak::Int
    evictions::Int
    bypasses::Int
end
_DecodedChunkCache(::Type{T}, budget) where {T<:Unsigned} =
    _DecodedChunkCache{T}(Dict{Int,_DecodedCacheEntry{T}}(),budget,0,0,0,0,0)

function _pack_decoded_chunk(nls::Dict{T,Vector{T}}, range) where T
    offsets=Vector{Int}(undef,length(range)+1); offsets[1]=1
    for (i,v) in enumerate(range)
        offsets[i+1]=offsets[i]+length(nls[T(v)])
    end
    targets=Vector{T}(undef,offsets[end]-1)
    for (i,v) in enumerate(range)
        ns=nls[T(v)]; sort!(ns)
        copyto!(targets,offsets[i],ns,1,length(ns))
    end
    return _CachedAdjacencyChunk{T}(first(range),offsets,targets)
end
function _cached_neighbors(chunk::_CachedAdjacencyChunk{T}, v::Int) where T
    i=v-chunk.first_vertex+1
    first=chunk.offsets[i]; len=chunk.offsets[i+1]-first
    result=Vector{T}(undef,len)
    copyto!(result,1,chunk.targets,first,len)
    return result
end
function _decoded_cache_get(cache, c)
    entry=get(cache.entries,c,nothing)
    entry === nothing && return nothing
    cache.clock += 1; entry.used=cache.clock
    return entry.chunk
end
function _decoded_cache_admit!(cache::_DecodedChunkCache{T}, c, block) where T
    cache.clock += 1
    entry=_DecodedCacheEntry{T}(block,cache.clock)
    # An oversized chunk must not evict other useful entries merely to fail admission.
    if cache.budget == 0 || Base.summarysize(Dict(c=>entry)) > cache.budget
        cache.bypasses += 1
        return
    end
    cache.entries[c]=entry
    bytes=Base.summarysize(cache.entries)
    while bytes > cache.budget
        victim=first(keys(cache.entries))
        for (id,e) in cache.entries
            e.used < cache.entries[victim].used && (victim=id)
        end
        delete!(cache.entries,victim); cache.evictions += 1
        # Rebuild to release hash-table capacity left behind by evicted entries.
        cache.entries=Dict(id=>e for (id,e) in cache.entries)
        bytes=isempty(cache.entries) ? 0 : Base.summarysize(cache.entries)
    end
    cache.resident=bytes; cache.peak=max(cache.peak,bytes)
    return nothing
end

function _budgeted_chunk_reader(decode!, ::Type{T}, n, k, algorithm, budget) where T
    cache=_DecodedChunkCache(T,budget)
    # Query-local scratch prevents repeated recursive decoding when dependencies
    # exceed the persistent budget. It is released after the outermost query.
    scratch=Ref{Union{Nothing,Dict{Int,_CachedAdjacencyChunk{T}}}}(nothing)
    hits=Ref(0); misses=Ref(0); decodes=Ref(0); verts=Ref(0); arcs=Ref(0)
    function get_chunk(c)
        block=_decoded_cache_get(cache,c)
        block !== nothing && return block
        scratch[] === nothing && (scratch[]=Dict{Int,_CachedAdjacencyChunk{T}}())
        temporary=scratch[]::Dict{Int,_CachedAdjacencyChunk{T}}
        haskey(temporary,c) && return temporary[c]
        range=((c-1)*k+1):min(c*k,n)
        state=IndexedRangeState(T)
        nls=decode!(range,state)
        block=_pack_decoded_chunk(nls,range)
        temporary[c]=block
        decodes[] += 1; verts[] += length(range); arcs[] += length(block.targets)
        _decoded_cache_admit!(cache,c,block)
        return block
    end
    resolve(u)=_cached_neighbors(get_chunk((Int(u)-1)÷k+1),Int(u))
    function neighbors(v::Integer)
        1 <= v <= n || error("Vertex $v out of range [1, $n]")
        c=(Int(v)-1)÷k+1
        haskey(cache.entries,c) ? (hits[] += 1) : (misses[] += 1)
        try
            return resolve(v)
        finally
            scratch[]=nothing
        end
    end
    function reset!()
        cache.entries=Dict{Int,_DecodedCacheEntry{T}}()
        cache.clock=0; cache.resident=0; cache.peak=0; cache.evictions=0; cache.bypasses=0
        scratch[]=nothing
        hits[]=0; misses[]=0; decodes[]=0; verts[]=0; arcs[]=0
        return nothing
    end
    stats=(cache_hits=hits,cache_misses=misses,decode_calls=decodes,
        vertices_decoded=verts,arcs_decoded=arcs,total_chunks=cld(n,k),samp_k=k,
        cache_budget_bytes=budget,cache_resident_bytes=()->cache.resident,
        cache_peak_bytes=()->cache.peak,cache_evictions=()->cache.evictions,
        oversized_bypasses=()->cache.bypasses,cached_chunks=()->length(cache.entries),
        cache_actual_heap_bytes=()->Base.summarysize(cache))
    return (;resolve,api=(n=n,m=-1,neighbors=neighbors,algorithm=algorithm,stats=stats,reset_fn=reset!))
end
