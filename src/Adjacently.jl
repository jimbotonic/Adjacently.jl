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

module Adjacently

# First include and export the base types that other modules depend on
include("node_types.jl")
using .NodeTypes: Node, EmptyNode, AbstractNode
export Node, EmptyNode, AbstractNode

include("custom_types.jl")
using .CustomTypes: UInt24, UInt40
export UInt24, UInt40

include("custom_lightgraphs.jl")
using .CustomLightGraphs: SimpleDiGraph, SimpleGraph, SimpleEdge
export SimpleDiGraph, SimpleGraph, SimpleEdge

include("util.jl")
using .Util
export Util

include("io.jl")
using .IO
export IO

include("algo.jl")
using .Algo
export Algo

include("rw.jl")  # Load RandomWalks before PageRank
using .RandomWalks
export RandomWalks

include("pr.jl")  # Now PageRank will have access to both CustomTypes and RandomWalks
using .PageRank
export PageRank

include("graph.jl")
using .Graph
export Graph

include("cycles.jl")
using .Cycles
export Cycles

include("mgs.jl")
using .MGS
export MGS

end # module Adjacently
