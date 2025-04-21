#
# Adjacently: Julia Complex Directed Networks Library
# Copyright (C) 2016-2024 Jimmy Dubuisson <jimmy.dubuisson@gmail.com>
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

using Test, LightGraphs, Adjacently.io, Adjacently.graph, Adjacently.util

const GRAY_TO_RSA_TABLE_PATH = joinpath(@__DIR__, "data", "gray_to_rsa.txt")

const AMZ_DATASET_IN = joinpath(@__DIR__, "..", "datasets", "Amazon_0601", "Amazon0601.txt")
const AMZ_DATASET_OUT = "Amazon_0601_core"

const GG_DATASET_IN = joinpath(@__DIR__, "..", "datasets", "Web_Google", "web-Google.txt")
const GG_DATASET_OUT = "Web_Google_core"

const ARX_DATASET_IN = joinpath(@__DIR__, "..", "datasets", "Arxiv_HEP-PH", "Cit-HepPh.txt")
const ARX_DATASET_OUT = "Arxiv_HEP-PH_core"

const EAT_DATASET_IN = joinpath(@__DIR__, "..", "datasets", "EAT", "EATnew.net")
const EAT_DATASET_OUT = "EAT_rcore"

"""
	load_dataset(input_path::AbstractString; separator::AbstractChar=',', is_pajek::Bool=false)

Load graph from CSV adjacency list or Pajek file
"""
function load_dataset(input_path::AbstractString; separator::AbstractChar=',', is_pajek::Bool=false)
	if !is_pajek
		g = load_adjacency_list_from_csv(input_path, separator)
	else
	    g = load_graph_from_pajek(input_path)
	end
	return g
end

@testset "io.jl" begin
	@info("Loading Amazon dataset")
	amz_g = load_dataset(AMZ_DATASET_IN; separator='\t')
	# number of vertices and edges
	@test 403394 == convert(Int,nv(amz_g))
	@test 3387388 == ne(amz_g)
	
	@info("Getting core")
	amz_core,oni,noi = get_core(amz_g)
	@test 395234 == convert(Int,nv(amz_core))
	@test 3301092 == ne(amz_core)

	@info("getting reverse graph")
	amz_rcore = get_reverse_graph(amz_core) 
	@test 395234 == convert(Int,nv(amz_rcore))
	@test 3301092 == ne(amz_rcore)

	@info("Saving Amazon dataset (core) in MGS format")
	write_mgs3_graph(amz_core, AMZ_DATASET_OUT)

	@info("Saving Amazon dataset (core) in MGZ format")
	write_mgs3_huffman_graph(amz_core, amz_rcore, AMZ_DATASET_OUT)
	# serialize_to_jld(amz_core, "core", AMZ_DATASET_OUT)
	
	@info("Loading Amazon dataset (core) from MGS format")
	amz_core_mgs = load_mgs3_graph(AMZ_DATASET_OUT * ".mgs")
	@test 395234 == convert(Int,nv(amz_core_mgs))
	@test 3301092 == ne(amz_core_mgs)

	@info("Loading Amazon dataset (core) from MGZ format")
	amz_core_mgz = load_mgs3_huffman_graph(AMZ_DATASET_OUT * ".mgz")
	@test 395234 == convert(Int,nv(amz_core_mgz))
	@test 3301092 == ne(amz_core_mgz)
end


