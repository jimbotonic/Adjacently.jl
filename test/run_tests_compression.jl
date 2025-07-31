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
    samples_bits(samples::Vector{T}, encoding::Symbol) where {T<:Unsigned}

    Parameters:
    - samples: Vector of samples to encode
    - encoding: Encoding scheme (:elias_gamma, :elias_delta, :golomb, :fibonacci, :zeta)

Return Elias code as a BitVector for a given samples.
"""
function samples_bits(samples::Vector{T}, encoding::Symbol) where {T<:Unsigned}
    # Use a simple approach: encode each sample individually and collect the bits
    all_bits = Bool[]
    
    for sample in samples
        io = IOBuffer()
        writer = BitWriter(io)
        
        if encoding == :elias_gamma
            write_elias_gamma(writer, sample)
        elseif encoding == :elias_delta
            write_elias_delta(writer, sample)
        elseif encoding == :golomb
            write_golomb(writer, sample, GOLOMB_BASE)
        elseif encoding == :fibonacci
            write_fibonacci(writer, sample)
        elseif encoding == :zeta
            write_zeta(writer, sample, ZETA_BASE)
        else
            error("Invalid encoding scheme: $encoding")
        end
        
        # capture number of bits for this sample
        sample_nbits = writer.index - 1
        flush_bitwriter(writer; flush_last_bits=true)
        
        # read the bits for this sample
        seekstart(io)
        reader = BitReader(io)
        sample_bits = read_bits(reader, sample_nbits)
        
        # append to all_bits
        append!(all_bits, sample_bits)
    end
    
    return all_bits
end

"""
    decode_samples(bits::BitVector, n_samples::Int, algorithm::Symbol)

Parameters:
- bits: BitVector to decode
- n_samples: Number of samples to decode
- algorithm: Algorithm to use for decoding

Return decoded samples as a Vector{UInt64}.
"""
function decode_samples(bits::BitVector, n_samples::Int, algorithm::Symbol)
    @info "Decoding $n_samples samples with $algorithm from $(length(bits)) bits"
    reader = bitvector_to_bitreader(bits)
    decoded = Vector{UInt64}(undef, n_samples)
    for i in 1:n_samples
        try
            if algorithm == :elias_gamma
                decoded[i] = read_elias_gamma(reader, UInt64)
            elseif algorithm == :elias_delta
                decoded[i] = read_elias_delta(reader, UInt64)
            elseif algorithm == :golomb
                decoded[i] = read_golomb(reader, GOLOMB_BASE, UInt64)
            elseif algorithm == :fibonacci
                decoded[i] = read_fibonacci(reader, UInt64)
            elseif algorithm == :zeta
                decoded[i] = read_zeta(reader, ZETA_BASE, UInt64)
            else
                error("Invalid algorithm: $algorithm")
            end
        catch e
            @error "Error decoding sample $i with $algorithm: $e"
            rethrow(e)
        end
    end
    return decoded
end

@testset "Compression" begin
    # generate 1M samples from a powerlaw distribution
    k = 1.2
    min_value = 1
    max_value = 1000
    n_samples = 1000000
    values = powerlaw_sample(k, n_samples, min_value, max_value, UInt64)

    @test length(values) == n_samples
    @test all(values .>= 1)
    @test all(values .<= 1000)

    algorithms = [:elias_gamma, :elias_delta, :golomb, :fibonacci, :zeta]

    # encode the samples
    for algorithm in algorithms
        bits = samples_bits(values, algorithm)
        @info "$algorithm: $(length(bits)) bits ($(length(bits) / n_samples) bits/sample)"
    end

    # compute entropy of the samples
    samples_entropy = get_entropy(values)
    @info "Samples entropy: $samples_entropy"

    # decode the samples
    for algorithm in algorithms
        bits = samples_bits(values, algorithm)
        decoded = decode_samples(BitVector(bits), n_samples, algorithm)
        @test decoded == values
    end
    

end

