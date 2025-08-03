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

@testset "Amazon Graph Tests" begin
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

    # Create test directory if it doesn't exist
    mkpath(TEST_DIR)

	@info("Saving Amazon dataset (core) in MGS format")
	# NB: the output file is created with extension .mgs
	write_mgs3_graph(amz_core, joinpath(TEST_DIR, AMZ_DATASET_OUT))

	@info("Saving Amazon dataset (core) in MGZ format")
	# NB: the output file is created with extension .mgz
	write_compressed_mgs3_graph(amz_core, joinpath(TEST_DIR, AMZ_DATASET_OUT), :children, :huffman)
	
	@info("Loading Amazon dataset (core) from MGS format")
	# NB: the input file is created with extension .mgs
	amz_core_mgs = load_mgs3_graph(joinpath(TEST_DIR, AMZ_DATASET_OUT * ".mgs"))
	@test 395234 == convert(Int, nv(amz_core_mgs))
	@test 3301092 == ne(amz_core_mgs)

	@info("Loading Amazon dataset (core) from MGZ format")
	# NB: the input file is created with extension .mgz
	amz_core_mgz = load_compressed_mgs3_graph(joinpath(TEST_DIR, AMZ_DATASET_OUT * ".mgz"))
	@test 395234 == convert(Int, nv(amz_core_mgz))
	@test 3301092 == ne(amz_core_mgz)

    # Clean up test files
    #rm(joinpath(TEST_DIR, AMZ_DATASET_OUT * ".mgs"), force=true)
    #rm(joinpath(TEST_DIR, AMZ_DATASET_OUT * ".mgz"), force=true)
    # clean up the test directory
    #rm(TEST_DIR, force=true, recursive=true)  
end

@testset "Arxiv Graph Tests" begin
	@info("Loading Arxiv_HEP-PH dataset")
	g = load_dataset(ARX_DATASET_IN; separator='\t')
    
    # Test basic graph properties
    @test nv(g) == 34546
    @test ne(g) == 421578
    @test isapprox(density(g), 0.00035326041393102855, rtol=1e-10)
    
    # Test degree statistics
    @test maximum(outdegree(g)) == 411
    @test minimum(outdegree(g)) == 0
    @test isapprox(mean(outdegree(g)), 12.203380999247381, rtol=1e-10)
    @test isapprox(median(outdegree(g)), 8.0, rtol=1e-10)
    @test isapprox(std(outdegree(g)), 15.224473754546942, rtol=1e-10)
    
    # Test specific vertex neighborhood
    expected_neighbors = UInt16[2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12]
    @test outneighbors(g, vertices(g)[1]) == expected_neighbors
    
    # Create test directory if it doesn't exist
    mkpath(TEST_DIR)
    
    # Use test directory for output files
    mgs_output_file = joinpath(TEST_DIR, "Arxiv_HEP-PH_test")

    # Test MGS3 format preservation
    # NB: the output file is created with extension .mgs
    write_mgs3_graph(g, mgs_output_file)
    g2 = load_mgs3_graph(mgs_output_file * ".mgs")
    
    # Verify graph properties are preserved
    @test nv(g) == nv(g2)
    @test ne(g) == ne(g2)
    @test density(g) ≈ density(g2)
    @test maximum(outdegree(g)) == maximum(outdegree(g2))
    @test minimum(outdegree(g)) == minimum(outdegree(g2))
    @test mean(outdegree(g)) ≈ mean(outdegree(g2))
    @test median(outdegree(g)) ≈ median(outdegree(g2))
    @test std(outdegree(g)) ≈ std(outdegree(g2))
    @test outneighbors(g, vertices(g)[1]) == outneighbors(g2, vertices(g2)[1])
    
    # clean up test file
    #rm(mgs_output_file * ".mgs", force=true)
    # clean up the test directory
    #rm(TEST_DIR, force=true, recursive=true)  
end

@testset "Arxiv_HEP-PH Graph Tests 2" begin
    # First load and save the graph
    @info("Loading Arxiv_HEP-PH dataset")
    g = load_dataset(ARX_DATASET_IN; separator='\t')
    
    # Create test directory if it doesn't exist
    mkpath(TEST_DIR)
    
    # Use test directory for output files
    mgs_output_file = joinpath(TEST_DIR, "Arxiv_HEP-PH_core")
    mgz_output_file = joinpath(TEST_DIR, "Arxiv_HEP-PH_core_compressed")
    
    # Save the graph in MGS format
    # NB: the output file is created with extension .mgs
    write_mgs3_graph(g, mgs_output_file)
    
    # Now test loading the saved file
    @info("Testing MGS3 graph loading and properties")
    # NB: the input file is created with extension .mgs
    g_loaded = load_mgs3_graph(mgs_output_file * ".mgs")
    initial_vertices = nv(g_loaded)
    initial_edges = ne(g_loaded)
    
    # Test reverse graph
    rg = get_reverse_graph(g_loaded)
    @test ne(rg) == ne(g_loaded)
    
    # Test MGS4 format (Huffman compressed) writing and loading
    # NB: the output file is created with extension .mgz
    write_compressed_mgs3_graph(g_loaded, mgz_output_file, :children, :huffman)
    gb = load_compressed_mgs3_graph(mgz_output_file * ".mgz")
    
    # Verify graph properties are preserved
    @test nv(gb) == initial_vertices
    @test ne(gb) == initial_edges
    
    # clean up test files
    #rm(mgs_output_file * ".mgs", force=true)
    #rm(mgz_output_file * ".mgz", force=true)
    # clean up the test directory
    #rm(TEST_DIR, force=true, recursive=true)  
end
