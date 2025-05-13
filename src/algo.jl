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

module Algo

using LightGraphs, DataStructures, Logging
using ..CustomTypes: UInt24, UInt40

# Export the functions we want to make available
export tarjan, pearce, pearce_iterative

"""
    tarjan(g::AbstractGraph{T}) where {T<:Unsigned}

Tarjan algorithm (recursive version)

NB: successfully tested with FA core
NB: the recursive calls may create a stack overflow error
"""
function tarjan(g::AbstractGraph{T}) where {T<:Unsigned}
	sccs = Vector{Vector{T}}()
	n = nv(g)
	indices = zeros(T,n)
	lowlinks = zeros(T,n)
	S = T[]
	index = 1

	# recursive function
	function visit(v)
		indices[v] = index
		lowlinks[v] = index
		index += 1
		push!(S,v)
		children = outneighbors(g,v)
		for w in children
			# w was not visited yet
			if indices[w] == 0
				# recursive call
				visit(w)
				# back propagate the lowlink
				lowlinks[v] = min(lowlinks[v], lowlinks[w])
			# w is in S and hence in the current SCC
			elseif findfirst(S,w) != 0
				lowlinks[v] = min(lowlinks[v], indices[w])
			end
		end

		if lowlinks[v] == indices[v]
			component = T[]
			while true
				w = pop!(S)
				push!(component,w)
				w == v && break
			end
			push!(sccs,component)
		end
	end

	for v in vertices(g)
		indices[v] == 0 && visit(v)
	end

	return sccs
end

"""
    pearce(g::AbstractGraph{T}) where {T<:Unsigned}

Pearce version of Tarjan algorithm

NB: successfully tested with FA core
NB: the recursive calls may create a stack overflow error
"""
function pearce(g::AbstractGraph{T}) where {T<:Unsigned}
	n = nv(g)
	rindex = zeros(T,n)
	S = T[]
	index = 1
	c = n-1

	# recursive function
	function visit(v)
		root = true
		rindex[v] = index
		index += 1
		children = outneighbors(g,v)
		for w in children
			rindex[w] == 0 && visit(w)
			rindex[w] < rindex[v] && begin rindex[v] = rindex[w]; root = false end
		end
		if root
			index -= 1
			while !isempty(S) && rindex[v] <= rindex[S[end]]
				w = pop!(S)
				rindex[w] = c
				index -= 1
			end
			rindex[v] = c
			c -= 1
		else
			push!(S,v)
		end
	end

	for v in vertices(g)
		if rindex[v] == 0
			visit(v)
		end
	end

	return rindex
end

# Pearce algorithm state
mutable struct State{T}
	v::T
	stage::T
	root::Bool
end

# Make State constructor handle type conversion for custom types UInt24 and UInt40
# Direct constructor - no recursion
function State(x::T, y::T, visited::Bool) where T <: Unsigned
    return State{T}(x, y, visited)
end

# Type conversion constructor
function State(x::Integer, y::T, visited::Bool) where T <: Unsigned
    return State{T}(convert(T, x), y, visited)
end

"""
    pearce_iterative(g::AbstractGraph{T) where {T<:Unsigned}

Pearce version of Tarjan algorithm - iterative version
"""
function pearce_iterative(g::AbstractGraph{T}) where {T<:Unsigned}
	n = nv(g)
	rindex = zeros(T,n)
	S = T[]
	index = 1
	c = n-1

	function visit(v)
		states = State[]
		current_state = State(v,convert(T,0),true)
		push!(states, current_state)

		@label start
		while !isempty(states)
			current_state = pop!(states)
			active_loop = false

			if current_state.stage == 0
				rindex[current_state.v] = index
				index += 1
				active_loop = true
			end
			children = outneighbors(g,current_state.v)
			for w in children
				if active_loop
					if rindex[w] == 0
						current_state.stage = w
						push!(states,current_state)

						new_state = State(w,convert(T,0),true)
						push!(states,new_state)
						@goto start
					end
					rindex[w] < rindex[current_state.v] && begin rindex[current_state.v] = rindex[w]; current_state.root = false end
				elseif current_state.stage == w
					active_loop = true
					rindex[w] < rindex[current_state.v] && begin rindex[current_state.v] = rindex[w]; current_state.root = false end
				end
			end
			if current_state.root
				index -= 1
				while !isempty(S) && rindex[current_state.v] <= rindex[S[end]]
					w = pop!(S)
					rindex[w] = c
					index -= 1
				end
				rindex[current_state.v] = c
				c -= 1
			else
				push!(S,current_state.v)
			end
		end
	end

	for v in vertices(g)
		if rindex[v] == 0
			visit(v)
		end
	end

	return rindex
end

end # module Algo
