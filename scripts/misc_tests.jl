#
# Adjacently: Julia Complex Networks Library
# Copyright (C) 2016-2020 Jimmy Dubuisson <jimmy.dubuisson@gmail.com>
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

using Pkg
Pkg.activate(normpath(joinpath(@__DIR__, "..")))

include("../util.jl")
include("../io.jl")
include("../graph.jl")

#using Base.Test

###
# sorting algorithms
###

@info("############ Testing sorting functions")

A = Int[3,6,4,7,1,9,2,12,10]

@info("mergesort algorithm")
R = bottom_up_sort(A)	

@info("quicksort algorithm")
A2 = copy(A)
R2 = quicksort_iterative!(A2)	

# permutation array
PA = [5,7,1,3,2,4,6,9,8]

@assert R == R2 == PA

# sorted arrays
SA1 = [1,2,3,4,6,7,9,10,12] 
SA2 = [12,10,9,7,6,4,3,2,1]

S = get_sorted_array(A,R)
S2 = get_sorted_array(A,R,false)

@info("original array: ", A)
@info("permutation array (merge sort): ", R)
@info("permutation array (quicksort): ", R2)
@info("sorted array (increasing): ", S)
@info("sorted array (decreasing): ", S2)

@assert S == SA1
@assert S2 == SA2

@info("##########")

###
# search functions
###

@info("############ Testing search functions")

v = Int[1,2,4,6,8]
i = binary_search(v,3)
i2 = searchsortedfirst(v,3)
i3 = searchsortedfirst(v,4)
i4 = searchsortedfirst(v,9)
i5 = searchsortedfirst(v,0)

@info("original array: ", v)
@info("position of 3: ", i)
@info("insert position of 3: ", i2)
@info("insert position of 4: ", i3)
@info("insert position of 9: ", i4)
@info("insert position of 0: ", i5)

@info("##########")

###
# Huffman encoding
###

@info("############ Testing Huffman encoding functions")

v = UInt8[1,6,3,7,2,8,5,18,12,17,13,24,12,1,4]
#v = UInt8[1,1,1,1]

t = huffman_encoding(v)

@info("original array: ", v)
@info("Huffman tree: ", t)

S = BitArray{1}()
D = Array{UInt8,1}()
encode_tree!(t, S, D)

@info("lenght S: ", length(S))
@info("length D: ", length(D))

t2 = decode_tree!(S, D)

@info("decoded Huffman tree: ", t2)
		
B = BitArray{1}()
C = Dict{BitArray{1},UInt8}()
get_huffman_codes!(t2, C, B)

@info("Huffman codes: ", C)
@info("# codes: ", length(keys(C)))

@info("##########")

###
# loading MGS3 and writing MGS4
###

@info("############ Testing graph serialization functions")

# load graph
@info("loading MGS3 graph")
g = adjlist(UInt32, is_directed=true)
load_mgs3_graph(g,"../datasets/Arxiv_HEP-PH/Arxiv_HEP-PH_core.mgs")
@info("# vertices:", length(vertices(g)))
@info("# edges:", num_edges(g))

@info("getting reverse graph")
rg = get_reverse_graph(g) 
@info("# edges (rg):", num_edges(rg))

@info("writing MGS4 graph")
write_mgs4_graph(g, rg, "test")

@info("loading MGS4 graph")
g2 = adjlist(UInt32, is_directed=true)
load_mgs4_graph(g2,"test.mgz")
@info("# vertices:", length(vertices(g2)))
@info("# edges:", num_edges(g2))

@info("##########")
