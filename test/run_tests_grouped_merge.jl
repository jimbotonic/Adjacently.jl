using Test, Random, Adjacently
const C=Adjacently.Compression
@testset "Exact grouped model sums" begin
    rng=MersenneTwister(81)
    f=C._DecoderFrequencies(ones(UInt32,C._RC_NSYM,C._RC_NCTX))
    function checkgroups(f)
        for ctx in 1:C._RC_NCTX, g in 1:C._RC_NGROUP
            @test f.groups[g,ctx] == sum(f[ctx,s] for s in (1+(g-1)*8):min(g*8,C._RC_NSYM))
        end
    end
    checkgroups(f)
    for _ in 1:1000; f[rand(rng,1:length(f))]=rand(rng,UInt32(1):UInt32(10000)); end
    checkgroups(f); fill!(f,1); checkgroups(f)
    totals=fill(UInt32(C._RC_NSYM),C._RC_NCTX)
    for _ in 1:300_000; C._rc_update!(f,totals,1,150); end
    checkgroups(f)
    @test totals == vec(sum(f;dims=2))
end
@testset "Sorted LR run merge" begin
    rng=MersenneTwister(82)
    for T in (UInt8,UInt16,UInt32,Adjacently.CustomTypes.UInt24), n in 0:80, split in 0:n
        a=rand(rng,T.(0:120),n)
        want=sort(a); sort!(view(a,1:split)); sort!(view(a,split+1:n))
        @test C._merge_lr_runs!(a,split) === a
        @test a == want
    end
    @test_throws ArgumentError C._merge_lr_runs!(UInt16[],1)
    @test_throws ArgumentError C._merge_lr_runs!(UInt16[1],-1)
end
