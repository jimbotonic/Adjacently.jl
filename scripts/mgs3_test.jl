#
# Adjacently: Julia Complex Networks Library
# Copyright (C) 2016-2019  Jimmy Dubuisson <jimmy.dubuisson@gmail.com>
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

include("../io.jl")
include("../graph.jl")

g = adjlist(UInt32, is_directed=true)

filename = ARGS[1]

# load CSV adjacency list
load_adjacency_list_from_csv(UInt32, g, filename, '\t')

nvs,nes,dens = get_basic_stats(g)

# save graph in MGSv3 format
write_mgs3_graph(g, "Arxiv_HEP-PH")

# load graph in MGSv3 format
g2 = adjlist(UInt32, is_directed=true)
load_mgs3_graph(g2, "Arxiv_HEP-PH.mgs")

println(length(vertices(g)))
println(out_neighbors(vertices(g)[1],g))
println(length(vertices(g2)))
println(out_neighbors(vertices(g2)[1],g2))
#println(length(edges(g2)))
