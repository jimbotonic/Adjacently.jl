#
# Adjacently: Julia Complex Dorected Networks Library
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

using Pkg
Pkg.activate(normpath(joinpath(@__DIR__, "..")))

include("../src/io.jl")
include("../src/graph.jl")

using Statistics

# load CSV adjacency list
g = load_adjacency_list_from_csv("../datasets/Arxiv_HEP-PH/Cit-HepPh.txt", '\t')

nvs,nes,dens = get_basic_stats(g)

# save graph in MGSv3 format
write_mgs3_graph(g, "Arxiv_HEP-PH")

# load graph in MGSv3 format
gb = load_mgs3_graph("Arxiv_HEP-PH.mgs")

println("Original graph:")
println("# vertices: ", nv(g))
println("# edges: ", ne(g))
println("density: ", density(g))
println("max degree: ", maximum(outdegree(g)))
println("min degree: ", minimum(outdegree(g)))
println("avg degree: ", mean(outdegree(g)))
println("median degree: ", median(outdegree(g)))
println("std degree: ", std(outdegree(g)))
println("outneighbors of vertex 1: ", outneighbors(g,vertices(g)[1]))

println("--------------------------------")

println("MGSv3 graph:")
println("# vertices: ", nv(gb))
println("# edges: ", ne(gb))
println("density: ", density(gb))
println("max degree: ", maximum(outdegree(gb)))
println("min degree: ", minimum(outdegree(gb)))
println("avg degree: ", mean(outdegree(gb)))
println("median degree: ", median(outdegree(gb)))
println("std degree: ", std(outdegree(gb)))
println("outneighbors of vertex 1: ", outneighbors(gb,vertices(gb)[1]))

rm("Arxiv_HEP-PH.mgs", force=true)
