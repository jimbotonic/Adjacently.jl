using Test, Random, Adjacently
const C = Adjacently.Compression

@testset "Contiguous decoder frequencies preserve logical layout" begin
    f=C._DecoderFrequencies(ones(UInt32,C._RC_NSYM,C._RC_NCTX))
    reference=ones(UInt32,C._RC_NCTX,C._RC_NSYM)
    rng=MersenneTwister(73)
    for _ in 1:1000
        i=rand(rng,eachindex(reference)); x=rand(rng,UInt32)
        f[i]=x; reference[i]=x
    end
    @test f == reference
    @test collect(vec(f)) == vec(reference)
    @test collect(eachindex(f)) == collect(eachindex(reference))
    @test stride(f.data,1) == 1
    @test fill!(f,1) === f
    @test all(==(1),f)
    @test_throws BoundsError f[0]
    @test_throws BoundsError f[length(f)+1]
    # Identical updates, including rescaling, through both physical layouts.
    a=ones(UInt32,C._RC_NCTX,C._RC_NSYM); ta=fill(UInt32(C._RC_NSYM),C._RC_NCTX)
    b=C._DecoderFrequencies(copy(permutedims(a))); tb=copy(ta)
    for _ in 1:300_000
        C._rc_update!(a,ta,1,5); C._rc_update!(b,tb,1,5)
    end
    @test a == b && ta == tb
end

@testset "Consecutive prefix lookup" begin
    p=C._PrefixLists{UInt16}()
    @test isempty(p)
    for v in UInt16(200):UInt16(400); p[v]=UInt16[v]; end
    @test length(p) == 201
    @test Dict(p) == Dict(v=>UInt16[v] for v in UInt16(200):UInt16(400))
    @test !haskey(p,199) && !haskey(p,401)
    @test_throws KeyError p[199]
    @test_throws ArgumentError (p[UInt16(402)]=UInt16[])
    ns=UInt16[7,8]; p[UInt16(250)]=ns
    @test p[250] === ns
    empty!(p); p[UInt16(600)]=UInt16[]
    @test length(p) == 1 && haskey(p,600) && isempty(p[600])
end

@testset "Circular reference order and reuse" begin
    for window in (0,1,7,64,256)
        ring=C._ReferenceRing{UInt32}(); ordinary=UInt32[]
        for v in UInt32(1):UInt32(1000)
            push!(ring,v); push!(ordinary,v)
            if length(ring)>window
                @test popfirst!(ring) == popfirst!(ordinary)
            end
            @test collect(ring) == ordinary
        end
        @test_throws BoundsError ring[0]
        @test_throws BoundsError ring[length(ring)+1]
        data=ring.data; empty!(ring)
        @test isempty(ring) && ring.data === data
        @test_throws ArgumentError popfirst!(ring)
        push!(ring,UInt32(42)); @test ring[1] == 42
    end
end

@testset "Recycled scratch detaches answers and resets models" begin
    for dense in (false,true)
        st=C.IndexedRangeState(UInt16;dense)
        st.neighbor_lists[UInt16(20)]=UInt16[1,5]
        saved=st.neighbor_lists
        push!(st.reference_window,UInt16(20)); st.next_vertex=21
        enc=C.CtxRangeEncoder(); values=UInt64[1,8,2,300,7]
        for v in values; C.rc_encode_value!(enc,v); end
        bytes=C.rc_finish!(enc); st.resid=C.CtxRangeDecoder(bytes)
        @test C.rc_decode_value!(st.resid) == values[1]
        freq=st.resid.freq
        C._recycle_indexed_state!(st)
        @test saved[UInt16(20)] == UInt16[1,5]
        @test isempty(st.neighbor_lists) && isempty(st.reference_window) && st.next_vertex == 0
        d=C._indexed_source(false,st.resid,C.CtxRangeDecoder,bytes,[0,length(bytes)],1)
        @test d.freq === freq
        @test [C.rc_decode_value!(d) for _ in values] == values
    end
end
