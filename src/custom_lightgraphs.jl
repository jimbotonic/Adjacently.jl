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

import LightGraphs
import LightGraphs.SimpleGraphs: SimpleDiGraph, SimpleGraph, SimpleEdge

# Implement nv for our custom types
LightGraphs.nv(g::SimpleDiGraph{UInt24}) = convert(Int, g.ne)
LightGraphs.nv(g::SimpleDiGraph{UInt40}) = convert(Int, g.ne)

# Implement ne for our custom types
LightGraphs.ne(g::SimpleDiGraph{UInt24}) = sum(length.(g.fadjlist))
LightGraphs.ne(g::SimpleDiGraph{UInt40}) = sum(length.(g.fadjlist))

# Implement add_vertex! for our custom types
function LightGraphs.add_vertex!(g::SimpleDiGraph{UInt24})
    nvg = LightGraphs.nv(g)
    if nvg >= typemax(UInt24)
        return false
    end
    g.ne = convert(UInt24, nvg + one(Int))
    resize!(g.fadjlist, LightGraphs.nv(g))
    resize!(g.badjlist, LightGraphs.nv(g))
    # Initialize the new vertex's adjacency lists
    g.fadjlist[LightGraphs.nv(g)] = UInt24[]
    g.badjlist[LightGraphs.nv(g)] = UInt24[]
    return true
end

function LightGraphs.add_vertex!(g::SimpleDiGraph{UInt40})
    nvg = LightGraphs.nv(g)
    if nvg >= typemax(UInt40)
        return false
    end
    g.ne = convert(UInt40, nvg + one(Int))
    resize!(g.fadjlist, LightGraphs.nv(g))
    resize!(g.badjlist, LightGraphs.nv(g))
    # Initialize the new vertex's adjacency lists
    g.fadjlist[LightGraphs.nv(g)] = UInt40[]
    g.badjlist[LightGraphs.nv(g)] = UInt40[]
    return true
end

# Implement add_edge! for our custom types
function LightGraphs.add_edge!(g::SimpleDiGraph{UInt24}, x::UInt24, y::UInt24)
    # Ensure vertices exist
    (x <= LightGraphs.nv(g) && y <= LightGraphs.nv(g)) || return false
    
    # Initialize adjacency lists if they don't exist
    if length(g.fadjlist) < LightGraphs.nv(g)
        resize!(g.fadjlist, LightGraphs.nv(g))
        g.fadjlist = [isnothing(l) ? UInt24[] : l for l in g.fadjlist]  # Changed from loop to list comprehension
    end
    if length(g.badjlist) < LightGraphs.nv(g)
        resize!(g.badjlist, LightGraphs.nv(g))
        g.badjlist = [isnothing(l) ? UInt24[] : l for l in g.badjlist]  # Changed from loop to list comprehension
    end
    
    push!(g.fadjlist[Int(x)], y)
    push!(g.badjlist[Int(y)], x)
    
    return true
end

# Same for UInt40
function LightGraphs.add_edge!(g::SimpleDiGraph{UInt40}, x::UInt40, y::UInt40)
    # Ensure vertices exist
    (x <= LightGraphs.nv(g) && y <= LightGraphs.nv(g)) || return false
    
    # Initialize adjacency lists if they don't exist
    if length(g.fadjlist) < LightGraphs.nv(g)
        resize!(g.fadjlist, LightGraphs.nv(g))
        g.fadjlist = [isnothing(l) ? UInt40[] : l for l in g.fadjlist]  # Changed from loop to list comprehension
    end
    if length(g.badjlist) < LightGraphs.nv(g)
        resize!(g.badjlist, LightGraphs.nv(g))
        g.badjlist = [isnothing(l) ? UInt40[] : l for l in g.badjlist]  # Changed from loop to list comprehension
    end
    
    push!(g.fadjlist[Int(x)], y)
    push!(g.badjlist[Int(y)], x)
    
    return true
end

# Add edge method for SimpleEdge
function LightGraphs.add_edge!(g::SimpleDiGraph{UInt24}, e::SimpleEdge{UInt24})
    if length(g.fadjlist) < LightGraphs.nv(g)
        resize!(g.fadjlist, LightGraphs.nv(g))
        g.fadjlist = [isnothing(l) ? UInt24[] : l for l in g.fadjlist]
    end
    if length(g.badjlist) < LightGraphs.nv(g)
        resize!(g.badjlist, LightGraphs.nv(g))
        g.badjlist = [isnothing(l) ? UInt24[] : l for l in g.badjlist]
    end
    
    push!(g.fadjlist[Int(src(e))], dst(e))
    push!(g.badjlist[Int(dst(e))], src(e))
    return true
end

# Same for UInt40
function LightGraphs.add_edge!(g::SimpleDiGraph{UInt40}, e::SimpleEdge{UInt40})
    if length(g.fadjlist) < LightGraphs.nv(g)
        resize!(g.fadjlist, LightGraphs.nv(g))
        g.fadjlist = [isnothing(l) ? UInt40[] : l for l in g.fadjlist]
    end
    if length(g.badjlist) < LightGraphs.nv(g)
        resize!(g.badjlist, LightGraphs.nv(g))
        g.badjlist = [isnothing(l) ? UInt40[] : l for l in g.badjlist]
    end
    
    push!(g.fadjlist[Int(src(e))], dst(e))
    push!(g.badjlist[Int(dst(e))], src(e))
    return true
end

# nv implementations for SimpleGraph
LightGraphs.nv(g::SimpleGraph{UInt24}) = convert(Int, g.ne)
LightGraphs.nv(g::SimpleGraph{UInt40}) = convert(Int, g.ne)

# ne implementations for SimpleGraph
LightGraphs.ne(g::SimpleGraph{UInt24}) = sum(length.(g.adjlist)) ÷ 2  # divide by 2 because edges are counted twice in undirected graphs
LightGraphs.ne(g::SimpleGraph{UInt40}) = sum(length.(g.adjlist)) ÷ 2

# add_vertex! implementations for SimpleGraph
function LightGraphs.add_vertex!(g::SimpleGraph{UInt24})
    nvg = LightGraphs.nv(g)
    if nvg >= typemax(UInt24)
        return false
    end
    g.ne = convert(UInt24, nvg + one(Int))
    resize!(g.adjlist, LightGraphs.nv(g))
    # Initialize the new vertex's adjacency list
    g.adjlist[LightGraphs.nv(g)] = UInt24[]
    return true
end

function LightGraphs.add_vertex!(g::SimpleGraph{UInt40})
    nvg = LightGraphs.nv(g)
    if nvg >= typemax(UInt40)
        return false
    end
    g.ne = convert(UInt40, nvg + one(Int))
    resize!(g.adjlist, LightGraphs.nv(g))
    # Initialize the new vertex's adjacency list
    g.adjlist[LightGraphs.nv(g)] = UInt40[]
    return true
end

# add_edge! implementations for SimpleGraph
function LightGraphs.add_edge!(g::SimpleGraph{UInt24}, x::UInt24, y::UInt24)
    # Ensure vertices exist and edge is valid
    (x <= LightGraphs.nv(g) && y <= LightGraphs.nv(g) && x != y) || return false
    
    # Initialize adjacency lists if they don't exist
    if length(g.adjlist) < LightGraphs.nv(g)
        resize!(g.adjlist, LightGraphs.nv(g))
        g.adjlist = [isnothing(l) ? UInt24[] : l for l in g.adjlist]
    end
    
    # Add edge in both directions (undirected graph)
    push!(g.adjlist[Int(x)], y)
    push!(g.adjlist[Int(y)], x)
    
    return true
end

function LightGraphs.add_edge!(g::SimpleGraph{UInt40}, x::UInt40, y::UInt40)
    # Ensure vertices exist and edge is valid
    (x <= LightGraphs.nv(g) && y <= LightGraphs.nv(g) && x != y) || return false
    
    # Initialize adjacency lists if they don't exist
    if length(g.adjlist) < LightGraphs.nv(g)
        resize!(g.adjlist, LightGraphs.nv(g))
        g.adjlist = [isnothing(l) ? UInt40[] : l for l in g.adjlist]
    end
    
    # Add edge in both directions (undirected graph)
    push!(g.adjlist[Int(x)], y)
    push!(g.adjlist[Int(y)], x)
    
    return true
end

# Add edge methods for SimpleEdge
function LightGraphs.add_edge!(g::SimpleGraph{UInt24}, e::SimpleEdge{UInt24})
    if length(g.adjlist) < LightGraphs.nv(g)
        resize!(g.adjlist, LightGraphs.nv(g))
        g.adjlist = [isnothing(l) ? UInt24[] : l for l in g.adjlist]
    end
    
    # Add edge in both directions
    push!(g.adjlist[Int(src(e))], dst(e))
    push!(g.adjlist[Int(dst(e))], src(e))
    return true
end

function LightGraphs.add_edge!(g::SimpleGraph{UInt40}, e::SimpleEdge{UInt40})
    if length(g.adjlist) < LightGraphs.nv(g)
        resize!(g.adjlist, LightGraphs.nv(g))
        g.adjlist = [isnothing(l) ? UInt40[] : l for l in g.adjlist]
    end
    
    # Add edge in both directions
    push!(g.adjlist[Int(src(e))], dst(e))
    push!(g.adjlist[Int(dst(e))], src(e))
    return true
end

import LightGraphs: outneighbors

# For SimpleDiGraph - outneighbors returns the forward adjacency list
function LightGraphs.outneighbors(g::SimpleDiGraph{UInt24}, v::Integer)
    return g.fadjlist[v]
end

function LightGraphs.outneighbors(g::SimpleDiGraph{UInt40}, v::Integer)
    return g.fadjlist[v]
end

# For SimpleGraph - outneighbors returns the full adjacency list since edges are undirected
function LightGraphs.outneighbors(g::SimpleGraph{UInt24}, v::Integer)
    return g.adjlist[v]
end

function LightGraphs.outneighbors(g::SimpleGraph{UInt40}, v::Integer)
    return g.adjlist[v]
end

# You might also want to implement inneighbors for SimpleDiGraph
function LightGraphs.inneighbors(g::SimpleDiGraph{UInt24}, v::Integer)
    return g.badjlist[v]
end

function LightGraphs.inneighbors(g::SimpleDiGraph{UInt40}, v::Integer)
    return g.badjlist[v]
end

# LightGraphs integration
import LightGraphs: AbstractGraphFormat
LightGraphs.is_directed(::Type{UInt24}) = true
LightGraphs.is_directed(::Type{UInt40}) = true

#####



