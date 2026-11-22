using Test, Random, Adjacently
const SC=Adjacently.Compression
const SIO=Adjacently.IO
@testset "Borrowed LR/STOP scratch is consumed into owned output" begin
    rng=MersenneTwister(93)
    for T in (UInt16,UInt32,Adjacently.CustomTypes.UInt24), mode in (:lr,:stop)
        scratch=T[999]; retained=Vector{T}[]; expected=Vector{T}[]
        for count in (0,1,2,20,2000,3,0,50)
            ns=sort!(unique!(T.(rand(rng,1:6000,count))))
            enc=SC.CtxRangeEncoder();w=SIO.BitWriter()
            old=SC._RESID_SINK[]
            try
                SC._RESID_SINK[]=enc
                if mode==:lr
                    SC._write_intervals_lr(w,ns,:context_range,4;vertex_id=T(3000),tight_intervals=true)
                else
                    SC._write_stop_delta(w,ns,:context_range;vertex_id=T(3000))
                end
            finally
                SC._RESID_SINK[]=old
            end
            SIO.flush_bitwriter(w;flush_last_bits=true)
            r=SIO.BitReader(collect(SIO.get_bytes(w)))
            d=SC.CtxRangeDecoder(SC.rc_finish!(enc))
            got=mode==:lr ? SC._read_intervals_lr_source(r,:context_range,4,T,T(3000),true,d,scratch) :
                SC._read_stop_delta_source(r,:context_range,T,T(3000),d,scratch)
            @test got===scratch
            @test got==ns
            owned=SC._merge_sorted_owned!(T[],got)
            @test owned!==scratch
            push!(retained,owned);push!(expected,ns)
            @test retained==expected
        end
    end
end
