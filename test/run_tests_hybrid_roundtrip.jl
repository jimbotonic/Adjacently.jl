#!/usr/bin/env julia

using Pkg
Pkg.activate(".")

using Test

import Adjacently.Compression: write_hybrid_mix_encoded_list, read_hybrid_mix_encoded_list,
    delta_encode_vector, write_encoded_value
using Adjacently.IO: BitWriter, BitReader, write_bit, flush_bitwriter
using Adjacently.Constants: MIN_INTERVAL_LENGTH

const T = UInt32
const ENC = :fibonacci

function roundtrip_children(neighs::Vector{T}; use_mix::Bool=true)
    # delta encode and children-shift by +1
    deltas = delta_encode_vector(neighs) .+ T(1)

    # write
    buf = IOBuffer()
    w = BitWriter(buf)
    write_hybrid_mix_encoded_list(w, deltas, ENC, use_mix, MIN_INTERVAL_LENGTH, true)
    # emulate trailing stop exactly like children mode writer: vertex-flag 0 then stop value
    if use_mix
        write_bit(w, false)
    end
    write_encoded_value(w, T(1), ENC)
    flush_bitwriter(w; flush_last_bits=true)

    # read
    seekstart(buf)
    r = BitReader(buf)
    decoded = read_hybrid_mix_encoded_list(r, ENC, :children, T; stop_value=T(1))
    return decoded
end

@testset "Hybrid Children Roundtrip" begin
    @testset "delta-only A" begin
        local inp = T[5,7,10]
        local outp = roundtrip_children(inp)
        @info "delta-only A" inp outp
        @test outp == inp
    end
    @testset "delta-only B" begin
        local inp = T[2,3]
        local outp = roundtrip_children(inp)
        @info "delta-only B" inp outp
        @test outp == inp
    end
    @testset "run-length" begin
        local inp = T[5,10,15,20]
        local outp = roundtrip_children(inp)
        @info "run-length" inp outp
        @test outp == inp
    end
    @testset "interval" begin
        local inp = T[100,101,102,103]
        local outp = roundtrip_children(inp)
        @info "interval" inp outp
        @test outp == inp
    end
    @testset "mixed A" begin
        local inp = T[1,2,3,4,10,12,15]
        local outp = roundtrip_children(inp)
        @info "mixed A" inp outp
        @test outp == inp
    end
    @testset "mixed B" begin
        local inp = T[100,101,150,151,152,153,200,205,210]
        local outp = roundtrip_children(inp)
        @info "mixed B" inp outp
        @test outp == inp
    end
    @testset "single" begin
        local inp = T[42]
        # Use no per-value flag before stop to exercise alternate termination path
        local outp = roundtrip_children(inp; use_mix=false)
        @info "single" inp outp
        @test outp == inp
    end
end
