#
# Adjacently: Julia Complex Directed Networks Library
# Copyright (C) 2016-2026 Jimmy Dubuisson <jimmy@dubuisson.ch>
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

using Adjacently.MGS: write_bg_mgs3_graph

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
    nbits = writer.bits_in # capture number of bits before flush resets it
    flush_bitwriter(writer; flush_last_bits=true)
    # read the bits from the buffer
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

@testset "Elias Compression" begin
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
    compression_types = [:elias_gamma, :elias_delta]

    for encoding_type in encoding_types
        for compression_type in compression_types
            @info("Writing compressed MGS3 graph with $encoding_type encoding and $compression_type compression")
            write_bg_mgs3_graph(amz_core, mgs_output_file * "_" * string(encoding_type) * "_" * string(compression_type); coding_scheme=encoding_type, integer_encoding=compression_type)
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
