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