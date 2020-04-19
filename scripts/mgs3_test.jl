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

include("../src/io.jl")
include("../src/graph.jl")

g = adjlist(UInt32, is_directed=true)
g2 = SimpleDiGraph(UInt32)

filename = ARGS[1]

# load CSV adjacency list
load_adjacency_list_from_csv(UInt32, g, filename, '\t')
load_adjacency_list_from_csv(UInt32, g2, filename, '\t')

nvs,nes,dens = get_basic_stats(g)

# save graph in MGSv3 format
write_mgs3_graph(g, "Arxiv_HEP-PH")
write_mgs3_graph(g2, "Arxiv_HEP-PH2")

# load graph in MGSv3 format
ga = adjlist(UInt32, is_directed=true)
load_mgs3_graph(ga, "Arxiv_HEP-PH.mgs")

println(length(G.vertices(g)))
println(out_neighbors(G.vertices(g)[1],g))
println(length(G.vertices(ga)))
println(out_neighbors(G.vertices(ga)[1],g2))
#println(length(edges(ga)))

gb = SimpleDiGraph(UInt32)
load_mgs3_graph(gb, "Arxiv_HEP-PH2.mgs")

println(nv(g2))
println(outneighbors(g2,LG.vertices(g2)[1]))
println(nv(gb))
println(outneighbors(gb,LG.vertices(gb)[1]))
#println(ne(gb))
