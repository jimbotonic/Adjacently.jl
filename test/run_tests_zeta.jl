#
# Adjacently: Julia Complex Directed Networks Library
# Copyright (C) 2016-2025 Jimmy Dubuisson <jimmy.dubuisson@gmail.com>
#
# This program is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.
#
# This program is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.
#

include("run_tests_main.jl")

"""
    unary_bits(v::T) where {T<:Unsigned}

Return unary code as a BitVector for a given value.
"""
function unary_bits(v::T, invert::Bool = false) where {T<:Unsigned}
    io = IOBuffer()
    writer = BitWriter(io)
    write_unary_coding(writer, v, invert)
    nbits = writer.index - 1 # capture number of bits before flush resets it
    flush_bitwriter(writer; flush_last_bits=true)
    seekstart(io)
    reader = BitReader(io)
    return read_bits(reader, nbits)
end

function truncated_binary_bits(v::T, n::Int) where {T<:Unsigned}
    io = IOBuffer()
    writer = BitWriter(io)
    write_truncated_binary_coding(writer, v, n)
    nbits = writer.index - 1 # capture number of bits before flush resets it
    flush_bitwriter(writer; flush_last_bits=true)
    seekstart(io)
    reader = BitReader(io)
    return read_bits(reader, nbits)
end

function zeta_bits(v::T, k::Int) where {T<:Unsigned}
    io = IOBuffer()
    writer = BitWriter(io)
    write_zeta_coding(writer, v, k)
    nbits = writer.index - 1 # capture number of bits before flush resets it
    flush_bitwriter(writer; flush_last_bits=true)
    seekstart(io)
    reader = BitReader(io)
    return read_bits(reader, nbits)
end

@testset "Unary Coding" begin
    @test unary_bits(UInt8(0), true) == BitVector([0])
    @test unary_bits(UInt8(1), true) == BitVector([1, 0])
    @test unary_bits(UInt8(2), true) == BitVector([1, 1, 0])
    @test unary_bits(UInt8(3), true) == BitVector([1, 1, 1, 0])
    @test unary_bits(UInt8(4), true) == BitVector([1, 1, 1, 1, 0])
end

@testset "Truncated Binary Coding" begin
    # special case: n = 0
    @test_throws ArgumentError truncated_binary_bits(UInt8(0), 0)

    # n = 1
    @test truncated_binary_bits(UInt8(0), 1) == BitVector([])
    @test_throws ArgumentError truncated_binary_bits(UInt8(1), 1)

    # n = 5
    @test truncated_binary_bits(UInt8(0), 5) == BitVector([0, 0])
    @test truncated_binary_bits(UInt8(1), 5) == BitVector([0, 1])
    @test truncated_binary_bits(UInt8(2), 5) == BitVector([1, 0])
    @test truncated_binary_bits(UInt8(3), 5) == BitVector([1, 1, 0])
    @test truncated_binary_bits(UInt8(4), 5) == BitVector([1, 1, 1])

    # n = 7
    @test truncated_binary_bits(UInt8(0), 7) == BitVector([0, 0])
    @test truncated_binary_bits(UInt8(1), 7) == BitVector([0, 1, 0])
    @test truncated_binary_bits(UInt8(2), 7) == BitVector([0, 1, 1])
    @test truncated_binary_bits(UInt8(3), 7) == BitVector([1, 0, 0])
    @test truncated_binary_bits(UInt8(4), 7) == BitVector([1, 0, 1])
    @test truncated_binary_bits(UInt8(5), 7) == BitVector([1, 1, 0])
    @test truncated_binary_bits(UInt8(6), 7) == BitVector([1, 1, 1])
end

@testset "Zeta Coding" begin
    # for k = 1, it is the same as Gamma code
    @test zeta_bits(UInt8(1), 1) == BitVector([1])
    @test zeta_bits(UInt8(2), 1) == BitVector([0, 1, 0])
    @test zeta_bits(UInt8(3), 1) == BitVector([0, 1, 1])
    @test zeta_bits(UInt8(4), 1) == BitVector([0, 0, 1, 0, 0])
    @test zeta_bits(UInt8(5), 1) == BitVector([0, 0, 1, 0, 1])
    @test zeta_bits(UInt8(6), 1) == BitVector([0, 0, 1, 1, 0])
    @test zeta_bits(UInt8(7), 1) == BitVector([0, 0, 1, 1, 1])
    @test zeta_bits(UInt8(8), 1) == BitVector([0, 0, 0, 1, 0, 0, 0])
    
    # for k = 2
    @test zeta_bits(UInt8(1), 2) == BitVector([1, 0])
    @test zeta_bits(UInt8(2), 2) == BitVector([1, 1, 0])
    @test zeta_bits(UInt8(3), 2) == BitVector([1, 1, 1])
    @test zeta_bits(UInt8(4), 2) == BitVector([0, 1, 0, 0, 0])
    @test zeta_bits(UInt8(5), 2) == BitVector([0, 1, 0, 0, 1])
    @test zeta_bits(UInt8(6), 2) == BitVector([0, 1, 0, 1, 0])
    @test zeta_bits(UInt8(7), 2) == BitVector([0, 1, 0, 1, 1])
    @test zeta_bits(UInt8(8), 2) == BitVector([0, 1, 1, 0, 0, 0])
    
    # for k = 3
    @test zeta_bits(UInt8(1), 3) == BitVector([1, 0, 0])
    @test zeta_bits(UInt8(2), 3) == BitVector([1, 0, 1, 0])
    @test zeta_bits(UInt8(3), 3) == BitVector([1, 0, 1, 1])
    @test zeta_bits(UInt8(4), 3) == BitVector([1, 1, 0, 0])
    @test zeta_bits(UInt8(5), 3) == BitVector([1, 1, 0, 1])
    @test zeta_bits(UInt8(6), 3) == BitVector([1, 1, 1, 0])
    @test zeta_bits(UInt8(7), 3) == BitVector([1, 1, 1, 1])
    @test zeta_bits(UInt8(8), 3) == BitVector([0, 1, 0, 0, 0, 0, 0])
    
    # for k = 4
    @test zeta_bits(UInt8(1), 4) == BitVector([1, 0, 0, 0])
    @test zeta_bits(UInt8(2), 4) == BitVector([1, 0, 0, 1, 0])
    @test zeta_bits(UInt8(3), 4) == BitVector([1, 0, 0, 1, 1])
    @test zeta_bits(UInt8(4), 4) == BitVector([1, 0, 1, 0, 0])
    @test zeta_bits(UInt8(5), 4) == BitVector([1, 0, 1, 0, 1])
    @test zeta_bits(UInt8(6), 4) == BitVector([1, 0, 1, 1, 0])
    @test zeta_bits(UInt8(7), 4) == BitVector([1, 0, 1, 1, 1])
    @test zeta_bits(UInt8(8), 4) == BitVector([1, 1, 0, 0, 0])
end