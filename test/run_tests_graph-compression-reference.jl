#
# Adjacently: Julia Complex Directed Networks Library
# Copyright (C) 2016-2026 Anonymous (double-blind review)
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
using Printf

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

	# Test 1: Specific empty lists patterns to validate our recent fixes
	@info "Testing reference encoding with specific empty list patterns..."
	
	@testset "Empty lists patterns with reference encoding" begin
		test_patterns = [
			("Mixed empty pattern", Dict{UInt16,Vector{UInt16}}(
				UInt16(1) => UInt16[],
				UInt16(2) => UInt16[10, 20],
				UInt16(3) => UInt16[]
			)),
			("Empty at start and middle", Dict{UInt16,Vector{UInt16}}(
				UInt16(1) => UInt16[],
				UInt16(2) => UInt16[],
				UInt16(3) => UInt16[5, 7],
				UInt16(4) => UInt16[1, 2, 3]
			)),
			("Empty at end", Dict{UInt16,Vector{UInt16}}(
				UInt16(1) => UInt16[5, 10],
				UInt16(2) => UInt16[1, 3],
				UInt16(3) => UInt16[],
				UInt16(4) => UInt16[]
			)),
			("All empty", Dict{UInt16,Vector{UInt16}}(
				UInt16(1) => UInt16[],
				UInt16(2) => UInt16[],
				UInt16(3) => UInt16[]
			))
		]
		
		for (pattern_name, neighbor_lists) in test_patterns
			@info "Testing pattern: $pattern_name"
			
			#encodings = [:elias_gamma, :elias_delta, :fibonacci, :zeta]
			encodings = [:fibonacci, :zeta]
			for encoding in encodings
				for reference_enabled in [false, true]
					@info "  Encoding: $encoding, Reference: $reference_enabled"
					
					# Encode
					io = IOBuffer()
					writer = BitWriter(io)
					write_compressed_graph_data(writer, neighbor_lists, encoding, :children, true, reference_enabled)
					flush_bitwriter(writer; flush_last_bits=true)
					
					# Decode
					seekstart(io)
					reader = BitReader(io)
					decoded = read_compressed_graph_data(reader, UInt16(length(neighbor_lists)), encoding, :children, UInt16)
					
					# Verify
					@test length(decoded) == length(neighbor_lists)
					for v in keys(neighbor_lists)
						@test haskey(decoded, v)
						@test Set(decoded[v]) == Set(neighbor_lists[v])
					end
					
					@info "    Pattern $pattern_name with $encoding (ref=$reference_enabled): PASSED"
				end
			end
		end
	end

	# Test 2: Reference encoding on actual neighbor lists
	@info "Testing reference encoding on actual Amazon graph data..."	

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
		modes = [:children, :index]  # Test only children mode for large data
		
		for mode in modes
			@info "Testing encoding with mode: $mode"
			
			# Test both reference enabled and disabled
			for reference_enabled in [false, true]
				@info "Graph size: $(length(amz_neighbor_lists)) vertices"
				
				compression_start_time = time()
			
				for encoding in encodings
					@info "Testing with encoding = $encoding, mode = $mode, relabeling_scheme = $relabeling_scheme, reference_enabled = $reference_enabled"
					encoding_start_time = time()
					
					@info "    Starting compression..."
					# Encode
					io = IOBuffer()
					writer = BitWriter(io)
					
					# Add detailed timing for the compression step
					write_start_time = time()
					write_compressed_graph_data(writer, amz_neighbor_lists, encoding, mode, true, reference_enabled)
					write_time = time() - write_start_time
					
					flush_start_time = time()
					flush_bitwriter(writer; flush_last_bits=true)
					flush_time = time() - flush_start_time
					
					total_encoding_time = time() - encoding_start_time
					compressed_size = position(io) * 8

					@info "    Compression completed!"
					@info "    Write time: $(@sprintf("%.3f", write_time))s"
					@info "    Flush time: $(@sprintf("%.3f", flush_time))s" 
					@info "    Total encoding time: $(@sprintf("%.3f", total_encoding_time))s"
					@info "    Compressed size: $compressed_size bits ($(@sprintf("%.2f", compressed_size/ (1024*1024*8))) MB)"
					@info "    Compression rate: $(@sprintf("%.3f", compressed_size / gs)) bits per vertex"
					@info "    Compression rate: $(@sprintf("%.3f", compressed_size / es)) bits per edge"
					@info "--------------------------------"
				end # encoding
				
				total_compression_time = time() - compression_start_time
				@info "Total time for reference_enabled=$reference_enabled: $(@sprintf("%.3f", total_compression_time))s"
			end # reference_enabled
		end # mode
	end # relabeling_scheme
end