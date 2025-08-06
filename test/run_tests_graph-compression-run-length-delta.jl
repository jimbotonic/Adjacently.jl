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

@testset "Run-length + delta encoding (actual graph)" begin
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
    #relabeling_schemes = [:none, :in_degree, :out_degree, :degree, :pagerank, :lexicographic]
    relabeling_schemes = [:none, :lexicographic]

    # compression schemes
    #encodings = [:elias_gamma, :elias_delta, :fibonacci, :zeta]
    encodings = [:elias_delta, :fibonacci, :zeta]  

    vs = vertices(amz_core)
	# number of vertices
	gs = convert(UInt64, length(vs))

	# `n_bits_v` is the number of bits needed to represent the graph vertices
	n_bits_v = convert(UInt8, ceil(log(2, gs)))
	# Get appropriate custom UInt type based on number of bits needed
	T = infer_uint_custom_type(n_bits_v) 

	# number of edges
	es = ne(amz_core)

	# Test run-length + delta encoding on actual neighbor list
	@info "Testing run-length + delta encoding on neighbor lists..."	

	for relabeling_scheme in relabeling_schemes
		@info("Relabeling scheme: $relabeling_scheme")
		if relabeling_scheme != :none
			old_to_new_vertex_ids = relabel_vertices(amz_core, relabeling_scheme)
			amz_core = relabel_graph(amz_core, old_to_new_vertex_ids)
		end

		@info("Getting neighbor lists")
		amz_neighbor_lists = get_neighbor_lists(amz_core)
		@test length(amz_neighbor_lists) == 395234
		@info "Created $(length(amz_neighbor_lists)) neighbor lists"

		for encoding in encodings
			@info "Testing run-length + delta with $encoding, $relabeling_scheme..."
		
			# Test run-length + delta encoding on individual neighbor lists
			total_bits = 0
			successful_decodes = 0
	
			for v in keys(amz_neighbor_lists)
				original_list = amz_neighbor_lists[v]
				if !isempty(original_list)
					# Convert to appropriate type and sort
					list_T = T.(original_list)
					sort!(list_T)
			
					# Encode using run-length + delta
					io = IOBuffer()
					writer = BitWriter(io)
					write_run_length_delta(writer, encoding, list_T)
					flush_bitwriter(writer; flush_last_bits=true)
					
					total_bits += position(io) * 8
					
					# Decode
					seekstart(io)
					reader = BitReader(io)
					try
						decoded_list = read_run_length_delta(reader, encoding, T)
						if decoded_list == list_T
							successful_decodes += 1
						end
					catch e
						@info "Failed to decode list for vertex $v: $e"
					end
				end
			end
			@info "Run-length + delta $encoding: $total_bits bits ($(total_bits / es) bits/edge)"
			@info "Successfully decoded $successful_decodes out of $(length(amz_neighbor_lists)) lists"
			@test successful_decodes == length(amz_neighbor_lists)  # All lists should decode successfully
			@info "Run-length + delta $encoding, $relabeling_scheme round-trip test: PASSED"
			@info "--------------------------------"
		end # encoding
	end # relabeling_scheme
end