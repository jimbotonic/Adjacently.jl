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

@testset "Relabeling Vertices 1" begin
    g = load_dataset(AMZ_DATASET_IN; separator='\t')
    @test nv(g) == 403394
    @test ne(g) == 3387388
    @test isapprox(density(g), 2.0816473245108078e-5, rtol=1e-10)

    # relabel the vertices according to the in-degree
    relabeled_vertices = relabel_vertices(g, :in_degree)
    relabeled_g = relabel_graph(g, relabeled_vertices)

    @test nv(relabeled_g) == nv(g)
    @test ne(relabeled_g) == ne(g)
    @test isapprox(density(relabeled_g), density(g), rtol=1e-10)
end

@testset "Amazon Graph Tests (relabeling vertices)" begin
	@info("Loading Amazon dataset")
	amz_g = load_dataset(AMZ_DATASET_IN; separator='\t')
	
	@info("Getting core")
	amz_core,oni,noi = get_core(amz_g)

	@info("Getting reverse graph")
	amz_rcore = get_reverse_graph(amz_core) 

    # relabel core according to in-degree
    relabeled_vertices = relabel_vertices(amz_core, :in_degree)
    amz_core_relabeled = relabel_graph(amz_core, relabeled_vertices)

    # check that the relabeled core has the same number of vertices and edges
    @info("Checking that the relabeled core has the same number of vertices and edges")
    @test nv(amz_core_relabeled) == nv(amz_core)
    @test ne(amz_core_relabeled) == ne(amz_core)

    # relabel reverse graph according to in-degree
    relabeled_vertices = relabel_vertices(amz_rcore, :in_degree)
    amz_rcore_relabeled = relabel_graph(amz_rcore, relabeled_vertices)

    # check that the relabeled reverse core has the same number of vertices and edges
    @info("Checking that the relabeled reverse core has the same number of vertices and edges")
    @test nv(amz_rcore_relabeled) == nv(amz_rcore)
    @test ne(amz_rcore_relabeled) == ne(amz_rcore)

    # save the relabeled graphs
    @info("Saving relabeled core and reverse core (:children, :elias_delta)")
    write_compressed_mgs3_graph(amz_core_relabeled, joinpath(TEST_DIR, AMZ_DATASET_OUT * "_core_elias_delta_relabeled"), :children, :elias_delta)
    write_compressed_mgs3_graph(amz_rcore_relabeled, joinpath(TEST_DIR, AMZ_DATASET_OUT * "_rcore_elias_delta_relabeled"), :children, :elias_delta)
    
    # load the relabeled graphs
    @info("Loading relabeled core and reverse core (:children, :elias_delta)")
    amz_core_relabeled_mgz = load_compressed_mgs3_graph(joinpath(TEST_DIR, AMZ_DATASET_OUT * "_core_elias_delta_relabeled" * ".mgz"))
    amz_rcore_relabeled_mgz = load_compressed_mgs3_graph(joinpath(TEST_DIR, AMZ_DATASET_OUT * "_rcore_elias_delta_relabeled" * ".mgz"))

    # check that the loaded relabeled core and reverse core have the same number of vertices and edges
    @info("Checking that the loaded relabeled core and reverse core have the same number of vertices and edges")
    @test nv(amz_core_relabeled_mgz) == nv(amz_core_relabeled)
    @test ne(amz_core_relabeled_mgz) == ne(amz_core_relabeled)
    @test nv(amz_rcore_relabeled_mgz) == nv(amz_rcore_relabeled)
    @test ne(amz_rcore_relabeled_mgz) == ne(amz_rcore_relabeled)

    # clean up test files
    #rm(joinpath(TEST_DIR, AMZ_DATASET_OUT * "_core_relabeled" * ".mgs"), force=true)
    #rm(joinpath(TEST_DIR, AMZ_DATASET_OUT * "_rcore_relabeled" * ".mgs"), force=true)
end

@testset "Relabeling vertices 2" begin
    @info("Loading Amazon dataset")
	amz_g = load_dataset(AMZ_DATASET_IN; separator='\t')
	
	@info("Getting core")
	amz_core,oni,noi = get_core(amz_g)

    # create test directory if it doesn't exist
    mkpath(TEST_DIR)
    
    # Use test directory for output files
    mgs_output_file = joinpath(TEST_DIR, "amz_core_rl")

    criterion = [:in_degree, :out_degree, :degree, :pagerank, :lexicographic]
    encoding = [:children, :index]
    compression = [:elias_delta, :fibonacci, :zeta]

    for compression in compression
        for encoding in encoding
            for criterion in criterion
                relabeled_vertices = relabel_vertices(amz_core, criterion)
                amz_core_rl = relabel_graph(amz_core, relabeled_vertices)

                @info("Checking that the relabeled core has the same number of vertices and edges")
                @test nv(amz_core_rl) == nv(amz_core)
                @test ne(amz_core_rl) == ne(amz_core)

                @info("Writing compressed MGS3 graph with $criterion criterion, $encoding encoding and $compression compression")
                write_compressed_mgs3_graph(amz_core_rl, mgs_output_file * "_" * string(criterion) * "_" * string(encoding) * "_" * string(compression), encoding, compression)
            end
        end
    end

    # clean up the test directory
    # rm(TEST_DIR, force=true, recursive=true)  
end

