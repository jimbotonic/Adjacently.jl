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

@testset "Sorting Algorithms" begin
    A = UInt8[3,6,4,7,1,9,2,12,10]
    PA = UInt8[5, 7, 1, 3, 2, 4, 6, 9, 8]
    SA1 = UInt8[1, 2, 3, 4, 6, 7, 9, 10, 12] 
    SA2 = UInt8[12, 10, 9, 7, 6, 4, 3, 2, 1]
    
    # Test mergesort
    R = bottom_up_sort(A)
    @test R == PA
    
    # Test quicksort
    A2 = copy(A)
    R2 = quicksort_iterative_permutation!(A2)
    @test R2 == PA
    
    # Test sorted arrays
    S = get_sorted_array(A, R)
    S2 = get_sorted_array(A, R, false)
    @test S == SA1
    @test S2 == SA2
end

@testset "Binary Search" begin
    v = UInt8[1, 2, 4, 6, 8]
    
    # Test binary search
    @test binary_search(v, UInt8(3)) == -1  # not found
    
    # Test searchsortedfirst
    @test searchsortedfirst(v, UInt8(3)) == 3  # insert position
    @test searchsortedfirst(v, UInt8(4)) == 3  # existing element
    @test searchsortedfirst(v, UInt8(9)) == 6  # beyond end
    @test searchsortedfirst(v, UInt8(0)) == 1  # before start
end

"""
    elias_bits(v::T) where {T<:Unsigned}

    Parameters:
    - v: Value to encode
    - encoding: Encoding scheme (:elias_gamma or :elias_delta)

Return Elias code as a BitVector for a given value.
"""
function elias_bits(v::T, encoding::Symbol) where {T<:Unsigned}
    io = IOBuffer()
    writer = BitWriter(io)
    if encoding == :elias_gamma
        write_elias_gamma(writer, v)
    elseif encoding == :elias_delta
        write_elias_delta(writer, v)
    end
    nbits = writer.index - 1 # capture number of bits before flush resets it
    flush_bitwriter(writer; flush_last_bits=true)
    # read the bits from the buffer
    seekstart(io)
    reader = BitReader(io)
    return read_bits(reader, nbits)
end

"""
    golomb_bits(v::T, k::Int) where {T<:Unsigned}

    Parameters:
    - v: Value to encode
    - k: Number of bits to use for the Golomb coding

Return Golomb code as a BitVector for a given value.
"""
function golomb_bits(v::T, k::Int) where {T<:Unsigned}
    io = IOBuffer()
    writer = BitWriter(io)
    write_golomb(writer, v, k)
    nbits = writer.index - 1 # capture number of bits before flush resets it
    flush_bitwriter(writer; flush_last_bits=true)
    seekstart(io)
    reader = BitReader(io)
    return read_bits(reader, nbits)
end

@testset "Elias Encoding" begin
    @test_throws ArgumentError write_elias_gamma(BitWriter(IOBuffer()), UInt8(0))

    @test elias_bits(UInt8(1), :elias_gamma) == BitVector([1])
    @test elias_bits(UInt8(2), :elias_gamma) == BitVector([0, 1, 0])
    @test elias_bits(UInt8(3), :elias_gamma) == BitVector([0, 1, 1])
    @test elias_bits(UInt8(4), :elias_gamma) == BitVector([0, 0, 1, 0, 0])
    @test elias_bits(UInt8(5), :elias_gamma) == BitVector([0, 0, 1, 0, 1])
    @test elias_bits(UInt8(6), :elias_gamma) == BitVector([0, 0, 1, 1, 0])
    @test elias_bits(UInt8(7), :elias_gamma) == BitVector([0, 0, 1, 1, 1])

    @test elias_bits(UInt8(1), :elias_delta) == BitVector([1])
    @test elias_bits(UInt8(2), :elias_delta) == BitVector([0, 1, 0, 0])
    @test elias_bits(UInt8(3), :elias_delta) == BitVector([0, 1, 0, 1])
    @test elias_bits(UInt8(4), :elias_delta) == BitVector([0, 1, 1, 0, 0])
    @test elias_bits(UInt8(5), :elias_delta) == BitVector([0, 1, 1, 0, 1])
    @test elias_bits(UInt8(6), :elias_delta) == BitVector([0, 1, 1, 1, 0])
    @test elias_bits(UInt8(7), :elias_delta) == BitVector([0, 1, 1, 1, 1])
end

@testset "Golomb Encoding" begin
    @test golomb_bits(UInt8(0), 4) == BitVector([1, 0, 0])
    @test golomb_bits(UInt8(1), 4) == BitVector([1, 0, 1])
    @test golomb_bits(UInt8(2), 4) == BitVector([1, 1, 0])
    @test golomb_bits(UInt8(3), 4) == BitVector([1, 1, 1])
    @test golomb_bits(UInt8(4), 4) == BitVector([0, 1, 0, 0])
    @test golomb_bits(UInt8(5), 4) == BitVector([0, 1, 0, 1])
    @test golomb_bits(UInt8(6), 4) == BitVector([0, 1, 1, 0])
    @test golomb_bits(UInt8(7), 4) == BitVector([0, 1, 1, 1])
end

