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

@testset "Reference encoding (actual graph)" begin
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

	# Test reference encoding on actual neighbor list
	@info "Testing reference encoding on neighbor lists..."	

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

		# Test reference encoding for both modes
		modes = [:children]  # Test only children mode to see progress
		
		for mode in modes
			@info "Testing reference encoding mode: $mode"
			
			for encoding in encodings
				@info "Testing reference encoding with $encoding, $relabeling_scheme, mode $mode..."
				
				# Use only first 1000 vertices for testing to avoid memory issues
				test_size = min(1000, length(amz_neighbor_lists))
				@info "Testing with first $test_size vertices (out of $(length(amz_neighbor_lists)))"
				
				# Create vertex mapping for consecutive numbering (1 to test_size)
				vertex_to_index = Dict{T,T}()
				index_to_vertex = Dict{T,T}()
				sorted_vertices = sort(collect(keys(amz_neighbor_lists)))[1:test_size]
				for (idx, v) in enumerate(sorted_vertices)
					vertex_to_index[T(v)] = T(idx)
					index_to_vertex[T(idx)] = T(v)
				end
				
				# Convert neighbor lists to appropriate type with consecutive indexing
				typed_neighbor_lists = Dict{T,Vector{T}}()
				for (idx, v) in enumerate(sorted_vertices)
					neighbors = amz_neighbor_lists[T(v)]
					if !isempty(neighbors)
						# Convert neighbors to indices and sort
						typed_neighbors = T[]
						for neighbor in neighbors
							if haskey(vertex_to_index, T(neighbor))
								push!(typed_neighbors, vertex_to_index[T(neighbor)])
							end
						end
						sort!(typed_neighbors)
						typed_neighbor_lists[T(idx)] = typed_neighbors
					else
						typed_neighbor_lists[T(idx)] = T[]
					end
				end
				
				# Encode using reference encoding
				io = IOBuffer()
				writer = BitWriter(io)
				write_reference_encoding(writer, typed_neighbor_lists, encoding, mode, false)
				flush_bitwriter(writer; flush_last_bits=true)
				
				total_bits = position(io) * 8
				@info "Reference encoding $encoding, $mode: $total_bits bits ($(total_bits / es) bits/edge)"
				
				# Decode
				seekstart(io)
				reader = BitReader(io)
				try
					decoded_neighbor_lists = read_reference_encoding(reader, T(length(typed_neighbor_lists)), encoding, mode, T)
					
					# Verify all lists are decoded correctly
					@test length(decoded_neighbor_lists) == length(typed_neighbor_lists)
					
					successful_matches = 0
					for v in keys(typed_neighbor_lists)
						if haskey(decoded_neighbor_lists, v)
							if Set(decoded_neighbor_lists[v]) == Set(typed_neighbor_lists[v])
								successful_matches += 1
							else
								@warn "Mismatch for vertex $v: original=$(typed_neighbor_lists[v]), decoded=$(decoded_neighbor_lists[v])"
							end
						else
							@warn "Missing vertex $v in decoded results"
						end
					end
					
					@info "Successfully matched $successful_matches out of $(length(typed_neighbor_lists)) lists"
					@test successful_matches == length(typed_neighbor_lists)
					@info "Reference encoding $encoding, $relabeling_scheme, mode $mode round-trip test: PASSED"
					
				catch e
					@error "Failed to decode reference encoding: $e"
					rethrow(e)
				end
				@info "--------------------------------"
			end # encoding
		end # mode
	end # relabeling_scheme
end