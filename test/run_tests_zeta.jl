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

"""
    truncated_binary_bits(v::T, n::Int) where {T<:Unsigned}

Return truncated binary code as a BitVector for a given value.
"""
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

"""
    zeta_bits(v::T, k::Int) where {T<:Unsigned}

Return zeta code as a BitVector for a given value.
"""
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

"""
    bitvector_to_bitreader(bits::BitVector)

Create a BitReader from a BitVector for testing purposes.
"""
function bitvector_to_bitreader(bits::BitVector)
    io = IOBuffer()
    writer = BitWriter(io)
    
    # Write each bit from the BitVector
    for bit in bits
        write_bit(writer, bit)
    end
    
    # Flush and prepare for reading
    flush_bitwriter(writer; flush_last_bits=true)
    seekstart(io)
    
    return BitReader(io)
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

@testset "Zeta Decoding" begin
    @test read_zeta_coding(bitvector_to_bitreader(BitVector([1])), 1, UInt8) == UInt8(1)
    @test read_zeta_coding(bitvector_to_bitreader(BitVector([0, 1, 0])), 1, UInt8) == UInt8(2)
    @test read_zeta_coding(bitvector_to_bitreader(BitVector([0, 1, 1])), 1, UInt8) == UInt8(3)
    @test read_zeta_coding(bitvector_to_bitreader(BitVector([0, 0, 1, 0, 0])), 1, UInt8) == UInt8(4)
    @test read_zeta_coding(bitvector_to_bitreader(BitVector([0, 0, 1, 0, 1])), 1, UInt8) == UInt8(5)
    @test read_zeta_coding(bitvector_to_bitreader(BitVector([0, 0, 1, 1, 0])), 1, UInt8) == UInt8(6)
    @test read_zeta_coding(bitvector_to_bitreader(BitVector([0, 0, 1, 1, 1])), 1, UInt8) == UInt8(7)
    @test read_zeta_coding(bitvector_to_bitreader(BitVector([0, 0, 0, 1, 0, 0, 0])), 1, UInt8) == UInt8(8)
end

@testset "Zeta Compression" begin
    @info("Loading Amazon dataset")
	amz_g = load_dataset(AMZ_DATASET_IN; separator='\t')
	# number of vertices and edges
	@test 403394 == convert(Int,nv(amz_g))
	@test 3387388 == ne(amz_g)
	
	@info("Getting core")
	amz_core, oni, noi = get_core(amz_g)
	@test 395234 == convert(Int,nv(amz_core))
	@test 3301092 == ne(amz_core)

    # create test directory if it doesn't exist
    mkpath(TEST_DIR)
    
    # Use test directory for output files
    mgs_output_file = joinpath(TEST_DIR, "amz_core")

    encoding_types = [:children, :index]
    compression_types = [:zeta]

    for encoding_type in encoding_types
        for compression_type in compression_types
            @info("Writing compressed MGS3 graph with $encoding_type encoding and $compression_type compression")
            write_compressed_mgs3_graph(amz_core, mgs_output_file * "_" * string(encoding_type) * "_" * string(compression_type), encoding_type, compression_type)
        end
    end

    # Test MGS3 format preservation
    # NB: the output file is created with extension .mgs
    for encoding_type in encoding_types
        for compression_type in compression_types
            @info("Loading compressed MGS3 graph with $encoding_type encoding and $compression_type compression")
            loaded_graph = load_compressed_mgs3_graph(mgs_output_file * "_" * string(encoding_type) * "_" * string(compression_type) * ".mgz")
            
            @info("Verifying graph properties")
            @test nv(loaded_graph) == nv(amz_core)
            @test ne(loaded_graph) == ne(amz_core)
            @test density(loaded_graph) ≈ density(amz_core)
            @test maximum(outdegree(loaded_graph)) == maximum(outdegree(amz_core))
            @test minimum(outdegree(loaded_graph)) == minimum(outdegree(amz_core))
            @test mean(outdegree(loaded_graph)) ≈ mean(outdegree(amz_core))
            @test median(outdegree(loaded_graph)) ≈ median(outdegree(amz_core))
            @test std(outdegree(loaded_graph)) ≈ std(outdegree(amz_core))
        end
    end

    # clean up the test directory
    # rm(TEST_DIR, force=true, recursive=true)  
end