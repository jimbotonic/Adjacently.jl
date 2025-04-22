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

using Pkg, Match
Pkg.activate(normpath(joinpath(@__DIR__, "..")))

include("../src/io.jl")
include("../src/graph.jl")

filename = ARGS[1]

# Load graph
g = load_graph_from_pajek(filename)

# Get stats
nvs, nes, dens = get_basic_stats(g)
println("Vertices: $nvs, Edges: $nes, Density: $dens")

# Process vertices
for v in vertices(g)[1:100]
    println("Neighbors of $v: ", LightGraphs.outneighbors(g, v))
end

# save graph in MGSv3 format
write_mgs3_graph(g, "EAT")

# load graph in MGSv3 format
g2 = load_mgs3_graph("EAT.mgs")

# Get stats for loaded graph
nvs, nes, dens = get_basic_stats(g2)
println("Vertices: $nvs, Edges: $nes, Density: $dens")

# Process vertices
for v in vertices(g2)[1:100]
    println("Neighbors of $v: ", LightGraphs.outneighbors(g2, v))
end
