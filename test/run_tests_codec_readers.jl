# Packed-reader and bounded-model prerequisites for opt-in CG decoding.
# julia --project --startup-file=no test/run_tests_codec_readers.jl
using Test, Random, Logging, Adjacently, LightGraphs
const DIO = Adjacently.IO
const DC = Adjacently.Compression
const DM = Adjacently.MGS
global_logger(ConsoleLogger(stderr, Logging.Warn))

@testset "Packed bit reader" begin
    rng = MersenneTwister(53)
    bytes = rand(rng, UInt8, 40)
    bits = [((b >> i) & 1) != 0 for b in bytes for i in 7:-1:0]
    r = DIO.BitReader(bytes)
    @test r.buffer === bytes
    @test length(r.buffer) == 40
    @test r.length == 320
    @test DIO.peek_bit(r) == bits[1]
    @test r.bit_count == 0
    @test DIO.read_bits(r, 320) == bits
    @test_throws ErrorException DIO.read_bit(r)
    @test_throws ErrorException DIO.peek_bit(r)
    @test DIO.read_value(r, 0, UInt64) == 0
    @test_throws ArgumentError DIO.read_value(r, -1, UInt64)
    for T in (UInt8, UInt16, UInt32, UInt64, UInt128,
              Adjacently.CustomTypes.UInt24, Adjacently.CustomTypes.UInt40)
        for start in 0:15, n in 0:min(8sizeof(T), 128)
            r = DIO.BitReader(bytes)
            r.index = start + 1; r.bit_count = start
            want = zero(T)
            for b in bits[start+1:start+n]; want = (want << 1) | T(b); end
            @test DIO.read_value(r, n, T) == want
            @test r.index == start + n + 1
            @test r.bit_count == start + n
        end
    end
    a = DIO.BitReader(bytes); b = DIO.BitReader(bytes)
    DIO.read_value(a, 13, UInt16)
    @test b.index == 1
    @test a.buffer === b.buffer
    io = IOBuffer(bytes); seek(io, 3)
    @test DIO.read_bits(DIO.BitReader(io), 296) == bits[25:end]
    @test_throws ErrorException DIO.read_value(DIO.BitReader(UInt8[0]), 9, UInt16)
    @test_throws ErrorException DIO.read_bit(DIO.BitReader(UInt8[]))
end

@testset "Integer code roundtrips" begin
    for T in (UInt16, UInt32, UInt64, Adjacently.CustomTypes.UInt24),
        (writefn, readfn) in ((DC.write_fibonacci, DC.read_fibonacci),
                             (DC.write_elias_gamma, DC.read_elias_gamma),
                             (DC.write_elias_delta, DC.read_elias_delta))
        vals = T.([1, 2, 3, 4, 7, 8, 255, 256, 1023, 16384, 32767])
        io = IOBuffer(); w = DIO.BitWriter(io)
        for v in vals; writefn(w, v); end
        DIO.flush_bitwriter(w; flush_last_bits=true)
        r = DIO.BitReader(take!(io))
        @test [readfn(r, T) for _ in vals] == vals
    end
end

@testset "Range spans and model reuse" begin
    rng = MersenneTwister(54)
    values = [UInt64[], UInt64[0], rand(rng, UInt64(0):UInt64(1_000_000), 1024), zeros(UInt64, 300_000)]
    blobs = map(values) do vs
        e = DC.CtxRangeEncoder()
        for v in vs; DC.rc_encode_value!(e, v); end
        DC.rc_finish!(e)
    end
    backing = vcat(UInt8[0xff, 0xff], blobs..., fill(0xff, 16))
    d = DC.CtxRangeDecoder(UInt8[]); freq = d.freq; tot = d.tot
    offset = 2
    for (vs, blob) in zip(values, blobs)
        DC.reset_range_decoder!(d, backing, offset + 1, offset + length(blob))
        @test d.freq === freq && d.tot === tot && d.inp === backing
        @test [DC.rc_decode_value!(d) for _ in vs] == vs
        fresh = DC.CtxRangeDecoder(backing, offset + 1, offset + length(blob))
        @test [DC.rc_decode_value!(fresh) for _ in vs] == vs
        offset += length(blob)
    end
    for first in 1:5
        @test DC.CtxRangeDecoder(fill(0xff, 5), first, first - 1).code == 0
        @test DC.BinRangeDecoder(fill(0xff, 5), first, first - 1).code == 0
    end
    for D in (DC.CtxRangeDecoder, DC.BinRangeDecoder)
        @test_throws ArgumentError D(UInt8[0], 0, 1)
        @test_throws ArgumentError D(UInt8[0], 1, 2)
    end
    chunks = [Bool[], [true], rand(rng, 100_000) .< 0.85, falses(100_000)]
    bb = map(chunks) do bs
        e = DC.BinRangeEncoder()
        for b in bs; DC.brc_encode_bit!(e, b); end
        DC.brc_finish!(e)
    end
    bytes = vcat(bb..., fill(0xff, 16)); offset = 0
    d = DC.BinRangeDecoder(UInt8[]); c0 = d.c0; c1 = d.c1
    for (bs, blob) in zip(chunks, bb)
        DC.reset_range_decoder!(d, bytes, offset + 1, offset + length(blob))
        @test d.c0 === c0 && d.c1 === c1 && d.inp === bytes
        @test [DC.brc_decode_bit!(d) for _ in bs] == bs
        offset += length(blob)
    end
end
