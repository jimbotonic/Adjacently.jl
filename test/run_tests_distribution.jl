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

@testset "Distribution" begin
    @info("Loading Amazon dataset")
	amz_g = load_dataset(AMZ_DATASET_IN; separator='\t')
	# number of vertices and edges
	@test 403394 == convert(Int,nv(amz_g))
	@test 3387388 == ne(amz_g)
	
	@info("Getting core")
	amz_core,oni,noi = get_core(amz_g)
	@test 395234 == convert(Int,nv(amz_core))
	@test 3301092 == ne(amz_core)

	@info("Getting reverse graph")
	amz_rcore = get_reverse_graph(amz_core) 
	@test 395234 == convert(Int,nv(amz_rcore))
	@test 3301092 == ne(amz_rcore)

    @info("Computing in-degree entropy of core")
    in_degree_entropy = get_degree_entropy(amz_core, :in_degree)
    @test isapprox(in_degree_entropy, 0.8850680456612653, atol=1e-10)

    @info("Computing out-degree entropy of core")
    out_degree_entropy = get_degree_entropy(amz_core, :out_degree)
    @test isapprox(out_degree_entropy, 0.6033916390343523, atol=1e-10)

    @info("Computing in-degree entropy of reverse core")
    in_degree_entropy = get_degree_entropy(amz_rcore, :in_degree)
    @test isapprox(in_degree_entropy, 0.6033916390343523, atol=1e-10)

    @info("Computing out-degree entropy of reverse core")
    out_degree_entropy = get_degree_entropy(amz_rcore, :out_degree)
    @test isapprox(out_degree_entropy, 0.8850680456612653, atol=1e-10)
end

@testset "Powerlaw Sampling" begin
    @info("Sampling from powerlaw distribution")
    k = 1.2
    min_value = 1
    max_value = 1000
    n_samples = 1000000
    values = powerlaw_sample(k, n_samples, min_value, max_value, UInt64)
    @test length(values) == n_samples
    @test all(values .>= 1)
    @test all(values .<= 1000)

    # Plot histogram (normalized)
    hist = histogram(values;
        bins=1000,
        normalize=true,
        title="Power-law Histogram (k=$k)",
        xlabel="Value",
        ylabel="Probability",
        legend=false)

    # create test directory if it doesn't exist
    mkpath(TEST_DIR)
    
    # Use test directory for output files
    powerlaw_output_file = joinpath(TEST_DIR, "powerlaw_histogram_k$k.png")

    # Save the plot as PNG
    savefig(powerlaw_output_file)

    # Clean up test file
    #rm(powerlaw_output_file, force=true)
    # clean up the test directory
    #rm(TEST_DIR, force=true, recursive=true)  
end