#
# JCNL: Julia Complex Networks Library
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

include("../io.jl")
include("../graph.jl")

g = adjlist(UInt32, is_directed=true)

filename = ARGS[1]

# load CSV adjacency list
load_graph_from_pajek(UInt32, g, filename)

nvs,nes,dens = get_basic_stats(g)
println(nvs, " ", nes," ", dens)

# save graph in MGSv3 format
write_mgs3_graph(g, "EAT")

# load graph in MGSv3 format
g2 = adjlist(UInt32, is_directed=true)
load_mgs3_graph(g2, "EAT.mgs")

nvs,nes,dens = get_basic_stats(g2)
println(nvs, " ", nes," ", dens)

vs = vertices(g2)
for v in vertices(g2)[1:1000]
	println(out_neighbors(v,g2))
end
