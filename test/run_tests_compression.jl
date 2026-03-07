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

"""
    decode_samples(bits::BitVector, n_samples::Int, algorithm::Symbol)

Parameters:
- bits: BitVector to decode
- n_samples: Number of samples to decode
- encoding: Encoding to use for decoding

Return decoded samples as a Vector{UInt64}.
"""
function decode_samples(bits::BitVector, n_samples::Int, encoding::Symbol)
    @info "Decoding $n_samples samples with $encoding from $(length(bits)) bits"
    reader = bitvector_to_bitreader(bits)
    decoded = Vector{UInt64}(undef, n_samples)
    for i in 1:n_samples
        try
            if encoding == :elias_gamma
                decoded[i] = read_elias_gamma(reader, UInt64)
            elseif encoding == :elias_delta
                decoded[i] = read_elias_delta(reader, UInt64)
            elseif encoding == :golomb
                decoded[i] = read_golomb(reader, GOLOMB_BASE, UInt64)
            elseif encoding == :fibonacci
                decoded[i] = read_fibonacci(reader, UInt64)
            elseif encoding == :zeta
                decoded[i] = read_zeta(reader, ZETA_BASE, UInt64)
            else
                error("Invalid encoding: $encoding")
            end
        catch e
            @error "Error decoding sample $i with $encoding: $e"
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

    encodings = [:elias_gamma, :elias_delta, :golomb, :fibonacci, :zeta]

    # encode the samples
    for encoding in encodings
        bits = samples_bits(values, encoding)
        @info "$encoding: $(length(bits)) bits ($(length(bits) / n_samples) bits/sample)"
    end

    # compute entropy of the samples
    samples_entropy = get_entropy(values)
    @info "Samples entropy: $samples_entropy"

    # decode the samples
    for encoding in encodings
        bits = samples_bits(values, encoding)
        decoded = decode_samples(BitVector(bits), n_samples, encoding)
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
    encodings = [:elias_gamma, :elias_delta, :golomb, :fibonacci, :zeta]
    
    # Initialize results storage
    compression_rates = Dict{Symbol, Vector{Float64}}()
    for encoding in encodings
        compression_rates[encoding] = Float64[]
    end
    
    # Storage for entropy values
    entropy_values = Float64[]
    
    @info "Computing compression rates for different power law exponents..."
    
    # test each k value
    for k in k_values
        @info "Testing k = $k"
        
        # storage for this k value
        k_rates = Dict{Symbol, Vector{Float64}}()
        for encoding in encodings
            k_rates[encoding] = Float64[]
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
            for encoding in encodings
                bits = samples_bits(values, encoding)
                rate = length(bits) / n_samples
                push!(k_rates[encoding], rate)
            end
        end
        
        # compute averages for this k value
        for encoding in encodings
            avg_rate = sum(k_rates[encoding]) / n_repetitions
            push!(compression_rates[encoding], avg_rate)
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
    for (i, encoding) in enumerate(encodings)
        plot!(plot_obj, k_vec, compression_rates[encoding],
              label=string(encoding),
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
    for encoding in encodings
        rates = compression_rates[encoding]
        min_rate = minimum(rates)
        max_rate = maximum(rates)
        avg_rate = sum(rates) / length(rates)
        @info "$encoding: min=$min_rate, max=$max_rate, avg=$avg_rate bits/sample"
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
    for encoding in encodings
        @test length(compression_rates[encoding]) == length(k_values)
        @test all(rate -> rate > 0, compression_rates[encoding])
    end
    
    @info "Compression rate analysis completed successfully"
end

@testset "Plot Compression Rate (delta)" begin
    # Parameters for the delta compression rate analysis
    k_values = 1.1:0.1:3.0
    min_value = 1
    max_value = 1000
    n_samples = 100000
    # average over 3 runs
    n_repetitions = 3
    encodings = [:elias_gamma, :elias_delta, :golomb, :fibonacci, :zeta]
    
    # Initialize results storage for delta distribution
    delta_compression_rates = Dict{Symbol, Vector{Float64}}()
    for encoding in encodings
        delta_compression_rates[encoding] = Float64[]
    end
    
    # Storage for delta entropy values
    delta_entropy_values = Float64[]
    
    @info "Computing compression rates for delta distribution of power law samples..."
    
    # test each k value
    for k in k_values
        @info "Testing k = $k"
        
        # storage for this k value
        k_rates = Dict{Symbol, Vector{Float64}}()
        for encoding in encodings
            k_rates[encoding] = Float64[]
        end
        
        # storage for entropy values for this k
        k_entropies = Float64[]
        
        # run multiple repetitions for averaging
        for rep in 1:n_repetitions
            # generate power law samples
            raw_values = powerlaw_sample(k, n_samples, min_value, max_value, UInt64)
            
            # sort the samples and compute delta encoding
            sorted_values = sort(raw_values)
            delta_values = delta_encode_vector(sorted_values)
            # shift the delta values by 1 to avoid 0
            delta_values = delta_values .+ 1
            
            # compute entropy for the delta distribution
            delta_entropy = get_entropy(delta_values)
            push!(k_entropies, delta_entropy)
            
            # test each algorithm on delta values
            for encoding in encodings
                bits = samples_bits(delta_values, encoding)
                rate = length(bits) / n_samples
                push!(k_rates[encoding], rate)
            end
        end
        
        # compute averages for this k value
        for encoding in encodings
            avg_rate = sum(k_rates[encoding]) / n_repetitions
            push!(delta_compression_rates[encoding], avg_rate)
        end
        
        # compute average entropy for this k value
        avg_entropy = sum(k_entropies) / n_repetitions
        push!(delta_entropy_values, avg_entropy)
    end
    
    # Create the delta plot
    @info "Creating delta compression rate plot..."
    
    # Convert k_values to vector for plotting
    k_vec = collect(k_values)
    
    # Create plot with different colors for each algorithm
    delta_plot_obj = plot(
        title="Delta Compression Rate vs Power Law Exponent",
        xlabel="Power Law Exponent (k)",
        ylabel="Compression Rate (bits/sample)",
        legend=:topright,
        size=(800, 600),
        dpi=300
    )
    
    # Color palette for algorithms and entropy
    colors = [:blue, :red, :green, :purple, :orange, :black]
    
    # Add line for each algorithm
    for (i, encoding) in enumerate(encodings)
        plot!(delta_plot_obj, k_vec, delta_compression_rates[encoding],
              label=string(encoding),
              color=colors[i],
              linewidth=2,
              marker=:circle,
              markersize=4)
    end
    
    # Add delta entropy line
    plot!(delta_plot_obj, k_vec, delta_entropy_values,
          label="entropy",
          color=colors[6],
          linewidth=2,
          marker=:diamond,
          markersize=4,
          linestyle=:dash)
    
    # save the delta plot
    mkpath("test_data")
    savefig(delta_plot_obj, "test_data/delta_compression_rates.png")
    @info "Delta plot saved to test_data/delta_compression_rates.png"
    
    # print delta summary statistics
    @info "Delta Compression Rate Summary:"
    for encoding in encodings
        rates = delta_compression_rates[encoding]
        min_rate = minimum(rates)
        max_rate = maximum(rates)
        avg_rate = sum(rates) / length(rates)
        @info "$encoding: min=$min_rate, max=$max_rate, avg=$avg_rate bits/sample"
    end
    
    # Print delta entropy summary
    min_entropy = minimum(delta_entropy_values)
    max_entropy = maximum(delta_entropy_values)
    avg_entropy = sum(delta_entropy_values) / length(delta_entropy_values)
    @info "delta entropy: min=$min_entropy, max=$max_entropy, avg=$avg_entropy bits/sample"
    
    # basic validation tests for delta
    @test length(k_vec) == length(k_values)
    @test length(delta_entropy_values) == length(k_values)
    @test all(entropy -> entropy > 0, delta_entropy_values)
    for encoding in encodings
        @test length(delta_compression_rates[encoding]) == length(k_values)
        @test all(rate -> rate > 0, delta_compression_rates[encoding])
    end
    
    @info "Delta compression rate analysis completed successfully"
end

