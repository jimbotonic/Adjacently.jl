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

@testset "Basic delta encoding" begin
	# Test with different encodings
	encodings = [:elias_gamma, :elias_delta, :fibonacci, :zeta] 
	
	for encoding in encodings
		@info "Testing basic delta encoding with $encoding..."
		
		@testset "$encoding - Simple sequential" begin
			original = UInt16[1, 2, 3, 4, 5]
			io = IOBuffer()
			writer = BitWriter(io)
			write_delta(writer, original, encoding)
			flush_bitwriter(writer; flush_last_bits=true)
			
			seekstart(io)
			reader = BitReader(io)
			decoded = read_delta(reader, encoding, UInt16)
			@test decoded == original
		end
		
		@testset "$encoding - Large gaps" begin
			original = UInt16[1, 100, 200, 205, 210]
			io = IOBuffer()
			writer = BitWriter(io)
			write_delta(writer, original, encoding)
			flush_bitwriter(writer; flush_last_bits=true)
			
			seekstart(io)
			reader = BitReader(io)
			decoded = read_delta(reader, encoding, UInt16)
			@test decoded == original
		end
		
		@testset "$encoding - Single element" begin
			original = UInt16[42]
			io = IOBuffer()
			writer = BitWriter(io)
			write_delta(writer, original, encoding)
			flush_bitwriter(writer; flush_last_bits=true)
			
			seekstart(io)
			reader = BitReader(io)
			decoded = read_delta(reader, encoding, UInt16)
			@test decoded == original
		end
		
		@testset "$encoding - Empty list" begin
			original = UInt16[]
			io = IOBuffer()
			writer = BitWriter(io)
			write_delta(writer, original, encoding)
			flush_bitwriter(writer; flush_last_bits=true)
			
			seekstart(io)
			reader = BitReader(io)
			decoded = read_delta(reader, encoding, UInt16)
			@test decoded == original
		end
		
		@info "All basic delta tests for $encoding: PASSED"
	end
	
	@info "Basic delta encoding tests completed successfully!"
end

@testset "Reference encoding (basics)" begin
	# Test with different encodings
	encodings = [:elias_gamma, :elias_delta, :fibonacci, :zeta] 
	
	for encoding in encodings
		@info "Testing basic reference encoding with $encoding..."
		
		# Test 1: Simple case with small lists (will use direct encoding)
		@testset "$encoding - Simple direct encoding case" begin
			# Create neighbor lists as Dict{UInt8,Vector{UInt8}}
			# Use small lists to avoid reference encoding (length <= 3)
			neighbor_lists = Dict{UInt16,Vector{UInt16}}()
			neighbor_lists[UInt16(1)] = UInt16[10, 20, 30]      # Small list (direct encoding)
			neighbor_lists[UInt16(2)] = UInt16[40, 50, 60]      # Different small list (direct encoding)
			neighbor_lists[UInt16(3)] = UInt16[70, 80]          # Even smaller list (direct encoding)
			
			io = IOBuffer()
			writer = BitWriter(io)
			write_compressed_graph_data(writer, neighbor_lists, :children, encoding, true, true)
			flush_bitwriter(writer; flush_last_bits=true)
			
			# Debug: check what was written
			buffer_size = position(io)
			@info "  Encoded $(length(neighbor_lists)) vertices to $buffer_size bytes"
			
			seekstart(io)
			reader = BitReader(io)
			decoded = read_compressed_graph_data(reader, UInt16(length(neighbor_lists)), :children, encoding, UInt16)
			
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
			neighbor_lists = Dict{UInt16,Vector{UInt16}}()
			neighbor_lists[UInt16(1)] = UInt16[10, 20]      # Small different lists
			neighbor_lists[UInt16(2)] = UInt16[30, 40]      
			neighbor_lists[UInt16(3)] = UInt16[50, 60]      
			
			io = IOBuffer()
			writer = BitWriter(io)
			write_compressed_graph_data(writer, neighbor_lists, :children, encoding, true, true)
			flush_bitwriter(writer; flush_last_bits=true)
			
			seekstart(io)
			reader = BitReader(io)
			decoded = read_compressed_graph_data(reader, UInt16(length(neighbor_lists)), :children, encoding, UInt16)
			
			@test length(decoded) == length(neighbor_lists)
			for v in keys(neighbor_lists)
				@test Set(decoded[v]) == Set(neighbor_lists[v])
			end
			@info "  No references test: PASSED"
		end
		
		# Test 3: Empty lists case
		@testset "$encoding - Empty lists case" begin
			neighbor_lists = Dict{UInt16,Vector{UInt16}}()
			neighbor_lists[UInt16(1)] = UInt16[]           # Empty list
			neighbor_lists[UInt16(2)] = UInt16[10, 20]     # Non-empty list
			neighbor_lists[UInt16(3)] = UInt16[]           # Another empty list
			
			io = IOBuffer()
			writer = BitWriter(io)
			write_compressed_graph_data(writer, neighbor_lists, :children, encoding, true, true)
			flush_bitwriter(writer; flush_last_bits=true)
			
			seekstart(io)
			reader = BitReader(io)
			decoded = read_compressed_graph_data(reader, UInt16(length(neighbor_lists)), :children, encoding, UInt16)
			
			@test length(decoded) == length(neighbor_lists)
			for v in keys(neighbor_lists)
				@test Set(decoded[v]) == Set(neighbor_lists[v])
			end
			@info "  Empty lists test: PASSED"
		end
		
		# Test 4: Single vertex case
		@testset "$encoding - Single vertex case" begin
			neighbor_lists = Dict{UInt16,Vector{UInt16}}()
			neighbor_lists[UInt16(1)] = UInt16[10, 20, 30]  # Small list
			
			io = IOBuffer()
			writer = BitWriter(io)
			write_compressed_graph_data(writer, neighbor_lists, :children, encoding, true, true)
			flush_bitwriter(writer; flush_last_bits=true)
			
			seekstart(io)
			reader = BitReader(io)
			decoded = read_compressed_graph_data(reader, UInt16(length(neighbor_lists)), :children, encoding, UInt16)
			
			@test length(decoded) == 1
			@test Set(decoded[UInt16(1)]) == Set(neighbor_lists[UInt16(1)])
			@info "  Single vertex test: PASSED"
		end
		
		# Test 5: Identical small lists
		@testset "$encoding - Identical small lists" begin
			neighbor_lists = Dict{UInt16,Vector{UInt16}}()
			neighbor_lists[UInt16(1)] = UInt16[10, 20, 30]  # Small list
			neighbor_lists[UInt16(2)] = UInt16[10, 20, 30]  # Identical small list  
			neighbor_lists[UInt16(3)] = UInt16[10, 20, 30]  # Identical small list
			
			io = IOBuffer()
			writer = BitWriter(io)
			write_compressed_graph_data(writer, neighbor_lists, :children, encoding, true, true)
			flush_bitwriter(writer; flush_last_bits=true)
			
			seekstart(io)
			reader = BitReader(io)
			decoded = read_compressed_graph_data(reader, UInt16(length(neighbor_lists)), :children, encoding, UInt16)
			
			@test length(decoded) == length(neighbor_lists)
			for v in keys(neighbor_lists)
				@test Set(decoded[v]) == Set(neighbor_lists[v])
			end
			@info "  Identical small lists test: PASSED"
		end
		
		# Test 6: Mixed small lists
		@testset "$encoding - Mixed small lists" begin
			neighbor_lists = Dict{UInt16,Vector{UInt16}}()
			neighbor_lists[UInt16(1)] = UInt16[2, 5, 10]     # Small list (values >= 2)
			neighbor_lists[UInt16(2)] = UInt16[2, 5, 15]     # Similar small list (values >= 2)
			neighbor_lists[UInt16(3)] = UInt16[5, 10]        # Smaller subset  
			neighbor_lists[UInt16(4)] = UInt16[20, 25, 30]   # Different small list
			
			io = IOBuffer()
			writer = BitWriter(io)
			write_compressed_graph_data(writer, neighbor_lists, :children, encoding, true, true)
			flush_bitwriter(writer; flush_last_bits=true)
			
			seekstart(io)
			reader = BitReader(io)
			decoded = read_compressed_graph_data(reader, UInt16(length(neighbor_lists)), :children, encoding, UInt16)
			
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