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

@testset "Fibonacci Compression" begin
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
    compression_types = [:fibonacci]

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
