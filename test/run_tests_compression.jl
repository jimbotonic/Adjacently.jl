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

@testset "Plot Compression Rate" begin
    # Parameters for the compression rate analysis
    k_values = 1.1:0.1:3.0
    min_value = 1
    max_value = 1000
    n_samples = 100000
    # average over 3 runs
    n_repetitions = 3
    algorithms = [:elias_gamma, :elias_delta, :golomb, :fibonacci, :zeta]
    
    # Initialize results storage
    compression_rates = Dict{Symbol, Vector{Float64}}()
    for algorithm in algorithms
        compression_rates[algorithm] = Float64[]
    end
    
    # Storage for entropy values
    entropy_values = Float64[]
    
    @info "Computing compression rates for different power law exponents..."
    
    # test each k value
    for k in k_values
        @info "Testing k = $k"
        
        # storage for this k value
        k_rates = Dict{Symbol, Vector{Float64}}()
        for algorithm in algorithms
            k_rates[algorithm] = Float64[]
        end
        
        # storage for entropy values for this k
        k_entropies = Float64[]
        
        # run multiple repetitions for averaging
        for rep in 1:n_repetitions
            # generate power law samples
            values = powerlaw_sample(k, n_samples, min_value, max_value, UInt64)
            
            # compute entropy for this sample set
            sample_entropy = get_entropy(values)
            push!(k_entropies, sample_entropy)
            
            # test each algorithm
            for algorithm in algorithms
                bits = samples_bits(values, algorithm)
                rate = length(bits) / n_samples
                push!(k_rates[algorithm], rate)
            end
        end
        
        # compute averages for this k value
        for algorithm in algorithms
            avg_rate = sum(k_rates[algorithm]) / n_repetitions
            push!(compression_rates[algorithm], avg_rate)
        end
        
        # compute average entropy for this k value
        avg_entropy = sum(k_entropies) / n_repetitions
        push!(entropy_values, avg_entropy)
    end
    
    # Create the plot
    @info "Creating compression rate plot..."
    
    # Convert k_values to vector for plotting
    k_vec = collect(k_values)
    
    # Create plot with different colors for each algorithm
    plot_obj = plot(
        title="Compression Rate vs Power Law Exponent",
        xlabel="Power Law Exponent (k)",
        ylabel="Compression Rate (bits/sample)",
        legend=:topright,
        size=(800, 600),
        dpi=300
    )
    
    # Color palette for algorithms and entropy
    colors = [:blue, :red, :green, :purple, :orange, :black]
    
    # Add line for each algorithm
    for (i, algorithm) in enumerate(algorithms)
        plot!(plot_obj, k_vec, compression_rates[algorithm],
              label=string(algorithm),
              color=colors[i],
              linewidth=2,
              marker=:circle,
              markersize=4)
    end
    
    # Add entropy line
    plot!(plot_obj, k_vec, entropy_values,
          label="entropy",
          color=colors[6],
          linewidth=2,
          marker=:diamond,
          markersize=4,
          linestyle=:dash)
    
    # save the plot
    mkpath("test_data")
    savefig(plot_obj, "test_data/compression_rates.png")
    @info "Plot saved to test_data/compression_rates.png"
    
    # print summary statistics
    @info "Compression Rate Summary:"
    for algorithm in algorithms
        rates = compression_rates[algorithm]
        min_rate = minimum(rates)
        max_rate = maximum(rates)
        avg_rate = sum(rates) / length(rates)
        @info "$algorithm: min=$min_rate, max=$max_rate, avg=$avg_rate bits/sample"
    end
    
    # Print entropy summary
    min_entropy = minimum(entropy_values)
    max_entropy = maximum(entropy_values)
    avg_entropy = sum(entropy_values) / length(entropy_values)
    @info "entropy: min=$min_entropy, max=$max_entropy, avg=$avg_entropy bits/sample"
    
    # basic validation tests
    @test length(k_vec) == length(k_values)
    @test length(entropy_values) == length(k_values)
    @test all(entropy -> entropy > 0, entropy_values)
    for algorithm in algorithms
        @test length(compression_rates[algorithm]) == length(k_values)
        @test all(rate -> rate > 0, compression_rates[algorithm])
    end
    
    @info "Compression rate analysis completed successfully"
end

