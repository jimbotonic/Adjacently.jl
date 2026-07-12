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

@testset "Huffman Encoding" begin
    v = UInt8[1, 6, 3, 7, 2, 8, 5, 18, 12, 17, 13, 24, 12, 1, 4]
    
    # get the Huffman tree for the given values
    t = huffman_encoding(v)
    @test t !== nothing
    
    # encode the Huffman tree
    S = BitArray{1}()
    D = Array{UInt8,1}()
    encode_huffman_tree!(t, S, D)

    # decode the Huffman tree
    t2 = decode_huffman_tree!(S, D)
    @test t2 !== nothing
    
    # get the Huffman codes for the given values
    B = BitArray{1}()
    # dictionary of bitarrays to uint8
    # example: [1, 1, 0, 0, 0, 0] -> 0x05
    C = Dict{BitArray{1},UInt8}()
    get_huffman_codes!(t2, C, B)
    
    # Verify we have one code per unique value (no hidden vertices)
    @test length(keys(C)) == length(unique(v))

    # Verify each unique value has a code
    @test all(val -> any(code_val -> code_val == val, values(C)), unique(v))
    
    # Verify codes are prefix-free (no code is a prefix of another)
    codes = collect(keys(C))
    for i in eachindex(codes)
        for j in eachindex(codes)
            if i != j
                # Convert BitVectors to strings for prefix comparison
                code_i = join(Int.(codes[i]))
                code_j = join(Int.(codes[j]))
                @test !startswith(code_i, code_j)
            end
        end
    end
end
