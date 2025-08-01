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

@testset "Graph Compression" begin
    @info("Loading Amazon dataset")
	amz_g = load_dataset(AMZ_DATASET_IN; separator='\t')
	# number of vertices and edges
	@test 403394 == convert(Int,nv(amz_g))
	@test 3387388 == ne(amz_g)
	
	@info("Getting core")
	amz_core,oni,noi = get_core(amz_g)
	@test 395234 == convert(Int,nv(amz_core))
	@test 3301092 == ne(amz_core)

	@info("Getting reverse core")
	amz_rcore = get_reverse_graph(amz_core) 
	@test 395234 == convert(Int,nv(amz_rcore))
	@test 3301092 == ne(amz_rcore)

    # relabeling schemes
    relabel_schemes = [:in_degree, :out_degree, :degree, :pagerank, :lexicographic]

    # compression schemes
    # encodings = [:elias_gamma, :elias_delta, :fibonacci, :zeta]
    encodings = [:elias_gamma, :elias_delta, :fibonacci, :zeta]

    # first, no relabeling
    # get in- and out-degrees
    in_degrees, out_degrees = get_in_out_degrees(amz_core)

    vs = vertices(amz_core)
	# number of vertices
	gs = convert(UInt64, length(vs))

	# `n_bits_v` is the number of bits needed to represent the graph vertices
	n_bits_v = convert(UInt8, ceil(log(2, gs)))
	# Get appropriate custom UInt type based on number of bits needed
	T = infer_uint_custom_type(n_bits_v) 

	# number of edges
	es = ne(amz_core)

    # get the list of all children and delta children
    children = Vector{T}(undef, es + gs)
	delta_children = Vector{T}(undef, es + gs)

    for v in vs
        onei = outneighbors(amz_core, v)
        for nei in onei
            push!(children, T(nei))
        end
		# sort the children
		sort!(onei)
		# compute the delta
		if !isempty(onei)
			first_child = onei[1]
			for i in 2:length(onei)
				delta_children[i] = onei[i] - first_child
			end
		end
		# stop value
		push!(children, T(0))
		push!(delta_children, T(0))
    end

	# shift the children and delta children by 1
	children .+= 1
	delta_children .+= 1

	# compute theoretical entropy of the graph
	adj_entropy = get_graph_entropy(amz_core, :bits_per_edge)
	@info "Graph adjacency entropy: $adj_entropy bits/edge"
	deg_entropy = get_graph_entropy(amz_core, :bits_per_vertex)
	@info "Graph degree entropy: $deg_entropy bits/vertex"

	@info "--------------------------------"
	@info "No encoding: $(es * n_bits_v) bits ($(es * n_bits_v / es) bits/edge)"
	@info "--------------------------------"

	# encode the list of children + stop values	
	for encoding in encodings
		children_bits = samples_bits(children, encoding)
		children_delta_bits = samples_bits(delta_children, encoding)

		@info "Direct encoding $encoding: $(length(children_bits)) bits ($(length(children_bits) / es) bits/edge)"
		@info "Delta encoding $encoding: $(length(children_delta_bits)) bits ($(length(children_delta_bits) / es) bits/edge)"
		@info "--------------------------------"
	end
	
end