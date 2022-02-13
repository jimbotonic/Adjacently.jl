#
# Adjacently: Julia Complex Networks Library
# Copyright (C) 2016-2022 Jimmy Dubuisson <jimmy.dubuisson@gmail.com>
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

g = SimpleDiGraph{UInt32}()

filename = ARGS[1]

# load CSV adjacency list
load_graph_from_pajek(UInt32, g, filename)

nvs,nes,dens = get_basic_stats(g)
println(nvs, " ", nes," ", dens)

# save graph in MGSv3 format
write_mgs3_graph(g, "EAT")

# load graph in MGSv3 format
g2 = SimpleDiGraph{UInt32}()
load_mgs3_graph(g2, "EAT.mgs")

nvs,nes,dens = get_basic_stats(g2)
println(nvs, " ", nes," ", dens)

vs = vertices(g2)
for v in vertices(g2)[1:100]
	println(outneighbors(g2,v))
end
