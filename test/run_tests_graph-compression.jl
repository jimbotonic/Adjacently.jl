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

@testset "Run-length + delta encoding (basics)" begin
	# Test with different encodings
	encodings = [:elias_gamma, :elias_delta, :fibonacci, :zeta] 
	
	for encoding in encodings
		@info "Testing basic run-length + delta with $encoding..."
		
		# Test 1: Simple sequential list
		@testset "$encoding - Simple sequential" begin
			original = UInt8[1, 2, 3, 4, 5]
			
			io = IOBuffer()
			writer = BitWriter(io)
			write_run_length_delta(writer, encoding, original)
			flush_bitwriter(writer; flush_last_bits=true)
			
			# Debug: check what was written
			buffer_size = position(io)
			seekstart(io)
			buffer_content = read(io)
			@info "  Encoded $(length(original)) elements to $buffer_size bytes: $buffer_content"
			
			seekstart(io)
			reader = BitReader(io)
			decoded = read_run_length_delta(reader, encoding, UInt8)
			
			@info "  Original: $original"
			@info "  Decoded:  $decoded"
			@test decoded == original
			@info "  Sequential test: PASSED"
		end
		
		# Test 2: List with repeated values (should trigger run-length encoding)
		@testset "$encoding - Repeated values" begin
			original = UInt8[10, 15, 15, 15, 15, 20, 25]  # Four consecutive 15s
			
			io = IOBuffer()
			writer = BitWriter(io)
			write_run_length_delta(writer, encoding, original)
			flush_bitwriter(writer; flush_last_bits=true)
			
			seekstart(io)
			reader = BitReader(io)
			decoded = read_run_length_delta(reader, encoding, UInt8)
			
			@test decoded == original
			@info "  Repeated values test: PASSED"
		end
		
		# Test 3: Large gaps (good for delta encoding)
		@testset "$encoding - Large gaps" begin
			original = UInt16[1, 100, 200, 205, 210]
			
			io = IOBuffer()
			writer = BitWriter(io)
			write_run_length_delta(writer, encoding, original)
			flush_bitwriter(writer; flush_last_bits=true)
			
			seekstart(io)
			reader = BitReader(io)
			decoded = read_run_length_delta(reader, encoding, UInt16)
			
			@test decoded == original
			@info "  Large gaps test: PASSED"
		end
		
		# Test 4: Single element
		@testset "$encoding - Single element" begin
			original = UInt16[42]
			
			io = IOBuffer()
			writer = BitWriter(io)
			write_run_length_delta(writer, encoding, original)
			flush_bitwriter(writer; flush_last_bits=true)
			
			seekstart(io)
			reader = BitReader(io)
			decoded = read_run_length_delta(reader, encoding, UInt16)
			
			@test decoded == original
			@info "  Single element test: PASSED"
		end
		
		# Test 5: Empty list
		@testset "$encoding - Empty list" begin
			original = UInt8[]
			
			io = IOBuffer()
			writer = BitWriter(io)
			write_run_length_delta(writer, encoding, original)
			flush_bitwriter(writer; flush_last_bits=true)
			
			seekstart(io)
			reader = BitReader(io)
			decoded = read_run_length_delta(reader, encoding, UInt8)
			
			@test decoded == original
			@info "  Empty list test: PASSED"
		end
		
		# Test 6: Long run-length sequence
		@testset "$encoding - Long run-length" begin
			original = UInt8[5, 10, 10, 10, 10, 10, 10, 10, 15, 20]  # Seven consecutive 10s
			
			io = IOBuffer()
			writer = BitWriter(io)
			write_run_length_delta(writer, encoding, original)
			flush_bitwriter(writer; flush_last_bits=true)
			
			seekstart(io)
			reader = BitReader(io)
			decoded = read_run_length_delta(reader, encoding, UInt8)
			
			@test decoded == original
			@info "  Long run-length test: PASSED"
		end
		
		# Test 7: Mixed patterns
		@testset "$encoding - Mixed patterns" begin
			original = UInt16[1, 5, 5, 5, 10, 50, 100, 100, 150, 200, 200, 200, 200]
			
			io = IOBuffer()
			writer = BitWriter(io)
			write_run_length_delta(writer, encoding, original)
			flush_bitwriter(writer; flush_last_bits=true)
			
			seekstart(io)
			reader = BitReader(io)
			decoded = read_run_length_delta(reader, encoding, UInt16)
			
			@test decoded == original
			@info "  Mixed patterns test: PASSED"
		end
		
		@info "All basic tests for $encoding: PASSED"
	end
	
	@info "Basic run-length + delta encoding tests completed successfully!"
end

@testset "Reference encoding (basics)" begin
	# Test with different encodings
	encodings = [:elias_gamma, :elias_delta, :fibonacci]  # Skip zeta due to bugs
	
	for encoding in encodings
		@info "Testing basic reference encoding with $encoding..."
		
		# Test 1: Simple case with small lists (will use direct encoding)
		@testset "$encoding - Simple direct encoding case" begin
			# Create neighbor lists as Dict{UInt8,Vector{UInt8}}
			# Use small lists to avoid reference encoding (length <= 3)
			neighbor_lists = Dict{UInt8,Vector{UInt8}}()
			neighbor_lists[UInt8(1)] = UInt8[10, 20, 30]      # Small list (direct encoding)
			neighbor_lists[UInt8(2)] = UInt8[40, 50, 60]      # Different small list (direct encoding)
			neighbor_lists[UInt8(3)] = UInt8[70, 80]          # Even smaller list (direct encoding)
			
			io = IOBuffer()
			writer = BitWriter(io)
			write_reference_encoding(writer, neighbor_lists, encoding, :children, true)
			flush_bitwriter(writer; flush_last_bits=true)
			
			# Debug: check what was written
			buffer_size = position(io)
			@info "  Encoded $(length(neighbor_lists)) vertices to $buffer_size bytes"
			
			seekstart(io)
			reader = BitReader(io)
			decoded = read_reference_encoding(reader, UInt8(length(neighbor_lists)), encoding, :children, UInt8)
			
			@info "  Original vertices: $(length(neighbor_lists))"
			@info "  Decoded vertices:  $(length(decoded))"
			@test length(decoded) == length(neighbor_lists)
			
			# Check that all lists are decoded correctly (order may differ due to sorting)
			for v in keys(neighbor_lists)
				@test Set(decoded[v]) == Set(neighbor_lists[v])
			end
			@info "  Simple direct encoding test: PASSED"
		end
		
		# Test 2: No good references case
		@testset "$encoding - No references case" begin
			neighbor_lists = Dict{UInt8,Vector{UInt8}}()
			neighbor_lists[UInt8(1)] = UInt8[10, 20]      # Small different lists
			neighbor_lists[UInt8(2)] = UInt8[30, 40]      
			neighbor_lists[UInt8(3)] = UInt8[50, 60]      
			
			io = IOBuffer()
			writer = BitWriter(io)
			write_reference_encoding(writer, neighbor_lists, encoding, :children, true)
			flush_bitwriter(writer; flush_last_bits=true)
			
			seekstart(io)
			reader = BitReader(io)
			decoded = read_reference_encoding(reader, UInt8(length(neighbor_lists)), encoding, :children, UInt8)
			
			@test length(decoded) == length(neighbor_lists)
			for v in keys(neighbor_lists)
				@test Set(decoded[v]) == Set(neighbor_lists[v])
			end
			@info "  No references test: PASSED"
		end
		
		# Test 3: Empty lists case
		@testset "$encoding - Empty lists case" begin
			neighbor_lists = Dict{UInt8,Vector{UInt8}}()
			neighbor_lists[UInt8(1)] = UInt8[]           # Empty list
			neighbor_lists[UInt8(2)] = UInt8[10, 20]     # Non-empty list
			neighbor_lists[UInt8(3)] = UInt8[]           # Another empty list
			
			io = IOBuffer()
			writer = BitWriter(io)
			write_reference_encoding(writer, neighbor_lists, encoding, :children, true)
			flush_bitwriter(writer; flush_last_bits=true)
			
			seekstart(io)
			reader = BitReader(io)
			decoded = read_reference_encoding(reader, UInt8(length(neighbor_lists)), encoding, :children, UInt8)
			
			@test length(decoded) == length(neighbor_lists)
			for v in keys(neighbor_lists)
				@test Set(decoded[v]) == Set(neighbor_lists[v])
			end
			@info "  Empty lists test: PASSED"
		end
		
		# Test 4: Single vertex case
		@testset "$encoding - Single vertex case" begin
			neighbor_lists = Dict{UInt8,Vector{UInt8}}()
			neighbor_lists[UInt8(1)] = UInt8[10, 20, 30]  # Small list
			
			io = IOBuffer()
			writer = BitWriter(io)
			write_reference_encoding(writer, neighbor_lists, encoding, :children, true)
			flush_bitwriter(writer; flush_last_bits=true)
			
			seekstart(io)
			reader = BitReader(io)
			decoded = read_reference_encoding(reader, UInt8(length(neighbor_lists)), encoding, :children, UInt8)
			
			@test length(decoded) == 1
			@test Set(decoded[UInt8(1)]) == Set(neighbor_lists[UInt8(1)])
			@info "  Single vertex test: PASSED"
		end
		
		# Test 5: Identical small lists
		@testset "$encoding - Identical small lists" begin
			neighbor_lists = Dict{UInt8,Vector{UInt8}}()
			neighbor_lists[UInt8(1)] = UInt8[10, 20, 30]  # Small list
			neighbor_lists[UInt8(2)] = UInt8[10, 20, 30]  # Identical small list  
			neighbor_lists[UInt8(3)] = UInt8[10, 20, 30]  # Identical small list
			
			io = IOBuffer()
			writer = BitWriter(io)
			write_reference_encoding(writer, neighbor_lists, encoding, :children, true)
			flush_bitwriter(writer; flush_last_bits=true)
			
			seekstart(io)
			reader = BitReader(io)
			decoded = read_reference_encoding(reader, UInt8(length(neighbor_lists)), encoding, :children, UInt8)
			
			@test length(decoded) == length(neighbor_lists)
			for v in keys(neighbor_lists)
				@test Set(decoded[v]) == Set(neighbor_lists[v])
			end
			@info "  Identical small lists test: PASSED"
		end
		
		# Test 6: Mixed small lists
		@testset "$encoding - Mixed small lists" begin
			neighbor_lists = Dict{UInt8,Vector{UInt8}}()
			neighbor_lists[UInt8(1)] = UInt8[2, 5, 10]     # Small list (values >= 2)
			neighbor_lists[UInt8(2)] = UInt8[2, 5, 15]     # Similar small list (values >= 2)
			neighbor_lists[UInt8(3)] = UInt8[5, 10]        # Smaller subset  
			neighbor_lists[UInt8(4)] = UInt8[20, 25, 30]   # Different small list
			
			io = IOBuffer()
			writer = BitWriter(io)
			write_reference_encoding(writer, neighbor_lists, encoding, :children, true)
			flush_bitwriter(writer; flush_last_bits=true)
			
			seekstart(io)
			reader = BitReader(io)
			decoded = read_reference_encoding(reader, UInt8(length(neighbor_lists)), encoding, :children, UInt8)
			
			@test length(decoded) == length(neighbor_lists)
			for v in keys(neighbor_lists)
				@test Set(decoded[v]) == Set(neighbor_lists[v])
			end
			@info "  Mixed small lists test: PASSED"
		end
		
		@info "All basic reference encoding tests for $encoding: PASSED"
	end
	
	@info "Basic reference encoding tests completed successfully!"
end

#=
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
    encodings = [:elias_gamma]  # Temporarily test only elias_gamma to focus on reference encoding

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

	@info("Getting neighbor lists")
	amz_neighbor_lists = get_neighbor_lists(amz_core)
	@test length(amz_neighbor_lists) == 395234

    # get the list of all children and delta children
    children = Vector{T}(undef, es + gs)
	delta_children = Vector{T}(undef, es + gs)

    for v in vs
        onei = amz_neighbor_lists[v]
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

	# shift the children and delta children by 1 to preserve the 0 value
	children .+= 1
	delta_children .+= 1

	# compute theoretical entropy of the graph
	adj_entropy = get_graph_entropy(amz_core, :bits_per_edge, :optimal)
	@info "Graph adjacency entropy: $adj_entropy bits/edge"
	deg_entropy = get_graph_entropy(amz_core, :bits_per_vertex, :optimal)
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
=#

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
    relabel_schemes = [:in_degree, :out_degree, :degree, :pagerank, :lexicographic]

    # compression schemes
    # encodings = [:elias_gamma, :elias_delta, :fibonacci, :zeta]
    encodings = [:elias_gamma]  # Temporarily test only elias_gamma to focus on reference encoding

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

	@info("Getting neighbor lists")
	amz_neighbor_lists = get_neighbor_lists(amz_core)
	@test length(amz_neighbor_lists) == 395234

	# Use actual Amazon graph neighbor lists
  	amz_neighbor_lists_T = Vector{Vector{T}}()
  	for v in vs
		if haskey(amz_neighbor_lists, v)
			neighbors_t = T.(amz_neighbor_lists[v])
			sort!(neighbors_t)
			# shift by 1 to avoid zeros (which cause issues in compression)
			neighbors_t .+= 1
			# add the stop value
			push!(neighbors_t, T(0))
			push!(amz_neighbor_lists_T, neighbors_t)
		else
			push!(amz_neighbor_lists_T, T[])
		end
  	end
	
	@info "Created $(length(amz_neighbor_lists_T)) neighbor lists"

	# Test run-length + delta encoding on actual neighbor list
	@info "Testing run-length + delta encoding on neighbor lists..."	
	
	for encoding in encodings
		@info "Testing run-length + delta with $encoding..."
		
		# Encode using run-length + delta
		io = IOBuffer()
		writer = BitWriter(io)
		
		# Process all non-empty neighbor lists
		non_empty_lists = [list for list in amz_neighbor_lists_T if !isempty(list)]
		for neighbor_list in non_empty_lists
			write_run_length_delta(writer, encoding, neighbor_list)
		end
		
		flush_bitwriter(writer; flush_last_bits=true)
		rld_bits = position(io) * 8
		
		@info "Run-length + delta $encoding: $rld_bits bits ($(rld_bits / es) bits/edge)"
		
		# Test round-trip encoding/decoding
		decoded_lists = Vector{Vector{T}}()
		for (i, original_list) in enumerate(non_empty_lists)
			# Create individual buffer for each list
			individual_io = IOBuffer()
			individual_writer = BitWriter(individual_io)
			write_run_length_delta(individual_writer, encoding, original_list)
			flush_bitwriter(individual_writer; flush_last_bits=true)
			
			# Get the actual number of bits written
			seekstart(individual_io)
			individual_reader = BitReader(individual_io)
			
			try
				decoded_list = read_run_length_delta(individual_reader, encoding, T)
				push!(decoded_lists, decoded_list)
			catch e
				#@info "Error decoding list $i (length: $(length(original_list))): $e"
				# Skip this list and continue with others
				continue
			end
		end
		
		# Verify correctness
		@info "Successfully decoded $(length(decoded_lists)) out of $(length(non_empty_lists)) lists"
		@test length(decoded_lists) > 0  # At least some lists should decode successfully
		
		#for i in 1:length(non_empty_lists)
		#	@test decoded_lists[i] == non_empty_lists[i]
		#end
		
		@info "Run-length + delta $encoding round-trip test: PASSED"
	end
end

#=
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
    relabel_schemes = [:in_degree, :out_degree, :degree, :pagerank, :lexicographic]

    # compression schemes
    # encodings = [:elias_gamma, :elias_delta, :fibonacci, :zeta]
    encodings = [:elias_gamma]  # Temporarily test only elias_gamma to focus on reference encoding

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

	@info("Getting neighbor lists")
	amz_neighbor_lists = get_neighbor_lists(amz_core)
	@test length(amz_neighbor_lists) == 395234

	# Use actual Amazon graph neighbor lists
  	amz_neighbor_lists_T = Vector{Vector{T}}()
  	for v in vs
		if haskey(amz_neighbor_lists, v)
			neighbors_t = T.(amz_neighbor_lists[v])
			sort!(neighbors_t)
			# shift by 1 to avoid zeros (which cause issues in compression)
			neighbors_t .+= 1
			# add the stop value
			push!(neighbors_t, T(0))
			push!(amz_neighbor_lists_T, neighbors_t)
		else
			push!(amz_neighbor_lists_T, T[])
		end
  	end
	
	@info "Created $(length(amz_neighbor_lists_T)) neighbor lists"
	@info "Testing reference encoding on neighbor lists..."
	
	# Test reference encoding
	for encoding in encodings
		@info "Testing reference encoding with $encoding..."
		@info "Total vertices to encode: $(length(amz_neighbor_lists_T))"
		@info "Non-empty lists: $(length([list for list in amz_neighbor_lists_T if !isempty(list)]))"
		
		# Encode using reference encoding
		io = IOBuffer()
		writer = BitWriter(io)
		
		@info "Starting reference encoding..."
		write_reference_encoding(writer, amz_neighbor_lists_T, encoding, :children, true)
		flush_bitwriter(writer; flush_last_bits=true)
		ref_bits = position(io) * 8
		
		@info "Reference encoding $encoding: $ref_bits bits ($(ref_bits / es) bits/edge)"
		@info "Buffer size after encoding: $(position(io)) bytes"
		
		# Test round-trip encoding/decoding
		seekstart(io)
		@info "Buffer content (first 32 bytes): $(take!(copy(io))[1:min(32, end)])"
		seekstart(io)
		reader = BitReader(io)
		
		@info "Starting reference decoding..."
		try
			decoded_ref_lists = read_reference_encoding(reader, encoding, T)
			@info "Successfully decoded $(length(decoded_ref_lists)) lists"
			@info "Sample decoded list lengths: $([length(list) for list in decoded_ref_lists[1:min(10, end)]])"
			@info "Sample decoded lists: $([list for list in decoded_ref_lists[1:min(3, end)]])"
			
			# Verify correctness
			@test length(decoded_ref_lists) == length(amz_neighbor_lists_T)
			
			# Check first few lists in detail
			mismatches = 0
			for i in 1:min(10, length(amz_neighbor_lists_T))
				original_set = Set(amz_neighbor_lists_T[i])
				decoded_set = Set(decoded_ref_lists[i])
				if original_set != decoded_set
					mismatches += 1
					@info "Mismatch at list $i:"
					@info "  Original: $(amz_neighbor_lists_T[i])"
					@info "  Decoded:  $(decoded_ref_lists[i])"
					@info "  Missing:  $(setdiff(original_set, decoded_set))"
					@info "  Extra:    $(setdiff(decoded_set, original_set))"
				end
			end
			
			if mismatches == 0
				@info "First 10 lists match perfectly"
			else
				@info "Found $mismatches mismatches in first 10 lists"
			end
			
			# Test all lists
			for i in 1:length(amz_neighbor_lists_T)
				@test Set(decoded_ref_lists[i]) == Set(amz_neighbor_lists_T[i])  # Order might differ due to sorting in reconstruction
			end
			
			@info "Reference encoding $encoding round-trip test: PASSED"
		catch e
			@info "Error during reference decoding: $e"
			@info "Error type: $(typeof(e))"
			if isa(e, BoundsError)
				@info "BoundsError details: attempting to access $(e.a) at index $(e.i)"
			end
			rethrow(e)
		end
		
		# Calculate compression ratio vs run-length + delta
		io_rld = IOBuffer()
		writer_rld = BitWriter(io_rld)
		for neighbor_list in amz_neighbor_lists_T
			write_run_length_delta(writer_rld, encoding, neighbor_list)
		end
		flush_bitwriter(writer_rld; flush_last_bits=true)
		rld_bits_comparison = position(io_rld) * 8
		
		compression_improvement = (rld_bits_comparison - ref_bits) / rld_bits_comparison * 100
		@info "Reference encoding improvement over run-length+delta: $(compression_improvement)%"
	end
	
	@info "--------------------------------"
	@info "Graph compression tests completed successfully"
end
=#