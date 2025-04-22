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

# Load graph
g = load_graph_from_pajek("../datasets/EAT/EATnew.net")

# Get stats
# expected stats:
# Vertices: 23219, Edges: 325593, Density: 0.0006039580036712458
nvs, nes, dens = get_basic_stats(g)
println("Vertices: $nvs, Edges: $nes, Density: $dens")

nv1 = UInt16[]

nv45 = UInt16[0x002d, 0x0037, 0x00ac, 0x00bf, 0x00c7, 
0x03ed, 0x042a, 0x085a, 0x08be, 0x0ccc, 0x12a3, 0x14bc, 
0x2278, 0x25a4, 0x2c86, 0x2ffa, 0x3499, 0x35c4, 0x35e1, 
0x367c, 0x3737, 0x37bb, 0x3fed, 0x3fef, 0x436c, 0x4760, 
0x4987, 0x54f5, 0x59f0, 0x5aa2]

nv55 = UInt16[0x0003, 0x0037, 0x003b, 0x0044, 0x0046, 
0x0068, 0x0079, 0x0091, 0x0097, 0x00a7, 0x00b2, 0x00b7, 
0x09cf, 0x0add, 0x0e50, 0x0f59, 0x1004, 0x1071, 0x1495, 
0x14d7, 0x14f2, 0x1833, 0x183a, 0x1a05, 0x1ab4, 0x1e06, 
0x1e4d, 0x235e, 0x29f9, 0x2fd2, 0x3209, 0x32c8, 0x35c9, 
0x364d, 0x36ad, 0x36af, 0x36cf, 0x3a2a, 0x3ac6, 0x3c1b, 
0x3c1e, 0x45f9, 0x4748, 0x4f1f, 0x4ff4, 0x516a, 0x53ea, 
0x5a51, 0x5a7a, 0x5a7e]

nv88 = UInt16[0x0059, 0x005b, 0x005c, 0x00a9, 0x00b7, 
0x0368, 0x0ea9, 0x0eb3, 0x103e, 0x18e9, 0x1b7e, 0x1ebf, 
0x2078, 0x22a8, 0x22ae, 0x2653, 0x26fb, 0x2d31, 0x2d36, 
0x2d95, 0x2fb2, 0x307a, 0x36ad, 0x3735, 0x37cb, 0x37ef, 
0x3a2a, 0x3c91, 0x454e, 0x46cf, 0x4a79, 0x50ca, 0x540c, 
0x54f6, 0x5998, 0x5a4d, 0x5a51]

# Process vertices
println("Neighbors of vertex 1: ", LightGraphs.outneighbors(g, 1))
println("Neighbors of vertex 45: ", LightGraphs.outneighbors(g, 45))
println("Neighbors of vertex 55: ", LightGraphs.outneighbors(g, 55))
println("Neighbors of vertex 88: ", LightGraphs.outneighbors(g, 88))

# save graph in MGSv3 format
write_mgs3_graph(g, "EAT")

# load graph in MGSv3 format
g2 = load_mgs3_graph("EAT.mgs")

# Get stats for loaded graph
# expected stats:
# Vertices: 23219, Edges: 325593, Density: 0.0006039580036712458
nvs, nes, dens = get_basic_stats(g2)
println("Vertices: $nvs, Edges: $nes, Density: $dens")

# Process vertices
println("Neighbors of vertex 1: ", LightGraphs.outneighbors(g2, 1))
println("Neighbors of vertex 45: ", LightGraphs.outneighbors(g2, 45))
println("Neighbors of vertex 55: ", LightGraphs.outneighbors(g2, 55))
println("Neighbors of vertex 88: ", LightGraphs.outneighbors(g2, 88))

rm("EAT.mgs", force=true)