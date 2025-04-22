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

filename = ARGS[1]

# load CSV adjacency list
g = load_adjacency_list_from_csv(filename, '\t')

nvs,nes,dens = get_basic_stats(g)

# save graph in MGSv3 format
write_mgs3_graph(g, "Arxiv_HEP-PH")

# load graph in MGSv3 format
gb = load_mgs3_graph("Arxiv_HEP-PH.mgs")

println(nv(g))
println(outneighbors(g,vertices(g)[1]))
println(nv(gb))
println(outneighbors(gb,vertices(gb)[1]))
#println(ne(gb))
