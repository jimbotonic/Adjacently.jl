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

module RandomWalks

using LightGraphs, DataStructures, Logging
using ..CustomTypes: UInt24, UInt40
using ..CustomLightGraphs: SimpleDiGraph, SimpleGraph, SimpleEdge

# Export the functions we want to make available
export RW, 
       RW_aggregated, 
       US, 
       ARW, 
       ARW_flying, 
       MHRW, 
       MHRW_flying, 
       CC_MHRW_flying

""" 
    RW(g::AbstractGraph{T},n_steps::UInt64,starting_v::T=convert(T,1)) where {T<:Unsigned}

Random Walk

proceeds for n_steps starting from the specified vertex id
@return the list of visited nodes
"""
function RW(g::AbstractGraph{T},n_steps::UInt64,starting_v::T=convert(T,1)) where {T<:Unsigned}
	visited_nodes =  T[]
	v = starting_v
	for i in 1:n_steps
		nei = outneighbors(g,v)
		nv = nei[rand(1:length(nei))]
		push!(visited_nodes,v)
		v = nv
	end
	return visited_nodes
end

"""
    RW_aggregated(g::AbstractGraph{T},jumping_constant::Float64,starting_v::T=convert(T,1)) where {T<:Unsigned}

Random Walk

proceeds until a sink node is reached or if rand()>jumping_constant from the specified vertex id
@return an array vid position -> # visits
"""
function RW_aggregated(g::AbstractGraph{T},jumping_constant::Float64,starting_v::T=convert(T,1)) where {T<:Unsigned}
	vv =  zeros(UInt32,nv(g))
	v = starting_v
	vv[v] = 1
	while rand() > jumping_constant
		nei = outneighbors(g,v)
		length(nei) == 0 && break
		nv = nei[rand(1:length(nei))]
		vv[nv] += 1
		v = nv
	end
	return vv
end

""" 
    RW_aggregated(g::AbstractGraph{T},n_steps::UInt64,starting_v::T=convert(T,1)) where {T<:Unsigned}

Random Walk

proceeds for n_steps starting from the specified vertex id
@return an array vid position -> # visits
"""
function RW_aggregated(g::AbstractGraph{T},n_steps::UInt64,starting_v::T=convert(T,1)) where {T<:Unsigned}
	vv =  zeros(UInt32,nv(g))
	v = starting_v
	vv[v] = 1
	for i in 1:n_steps
		nei = outneighbors(g,v)
		nv = nei[rand(1:length(nei))]
		vv[nv] += 1
		v = nv
	end
	return vv
end

""" 
    US(g::AbstractGraph{T},n_steps::UInt64) where {T<:Unsigned}

Uniform sampling

proceeds for n_steps 
@return the list of sampled vertices
"""
function US(g::AbstractGraph{T},n_steps::UInt64) where {T<:Unsigned}
	visited_nodes =  T[]
	n = nv(g)
	for i in 1:n_steps
		pv = convert(T,rand(1:n))
		push!(visited_nodes,pv)
	end
	return visited_nodes
end

""" 
    get_flying_index(a::Array{Float64,1})

get the flying index (select a child with a probability equals to its score)

@input array of child scores
"""
function get_flying_index(a::Array{Float64,1})
	r = rand()
	sum = 0.
	pos = 1
	for i in a
		p_sum = sum
		sum += i
		if p_sum <= r && r <= sum
			return pos
		end
		pos += 1
	end
	return -1
end

"""
    ARW(g::AbstractGraph{T},n_steps::UInt64,rd::Array{Float64,1},starting_v::T=convert(T,1)) where {T<:Unsigned}

Avoiding Random Walk

rd: stochastic repulsive vector
@return an array vid position -> # visits
"""
function ARW(g::AbstractGraph{T},n_steps::UInt64,rd::Array{Float64,1},starting_v::T=convert(T,1)) where {T<:Unsigned}
	visited_nodes =  T[]
	v = starting_v
	push!(visited_nodes,v)
	rejected = 0
	for i in 1:n_steps
		nei = outneighbors(g,v)
		if length(nei) > 1
			# select a child randomly
			rpos = rand(1:length(nei))
			pv = nei[rpos]
			# compute the repulsive score (higher if selected child is less repulsive than current node)
			score = rd[v]/rd[pv]
			if rand() <= score
				push!(visited_nodes,pv)
				v = pv
			else
				rejected += 1
			end
		else
			# if there is only one child, follow the link
			pv = nei[1]
			push!(visited_nodes,pv)
			v = pv
		end
	end
	return visited_nodes, rejected
end

""" 
    ARW_flying(g::AbstractGraph{T},n_steps::UInt64,rd::Array{Float64,1},starting_v::T=convert(T,1)) where {T<:Unsigned}

Avoiding Random Walk (flying mode)

NB: flying mode avoids to get stuck at a given node to sample nodes more efficiently

@input rd: stochastic repulsive vector
@return an array vid position -> # visits
"""
function ARW_flying(g::AbstractGraph{T},n_steps::UInt64,rd::Array{Float64,1},starting_v::T=convert(T,1)) where {T<:Unsigned}
	visited_nodes =  T[]
	v = starting_v
	push!(visited_nodes,v)
	for i in 1:n_steps
		nei = outneighbors(g,v)
		if length(nei) > 1
			# compute repulsive scores of neighbors
			scores = stochastic([rd[v]/rd[i] for i in nei])
			pos = get_flying_index(scores)
			pv = nei[pos]
		else
			# if there is only one child, follow the link
			pv = nei[1]
		end
		push!(visited_nodes,pv)
		v = pv
	end
	return visited_nodes
end

""" 
    MHRW(g::AbstractGraph{T},n_steps::UInt64,in_degrees::Array{T,1},out_degrees::Array{T,1},starting_v::T=convert(T,1)) where {T<:Unsigned}

Metropolis-Hasting Random Walk

NB: in the case of undirected graphs, MHRW corrects the sampling bias by avoiding high degree nodes
NB: but in the case of directed graphs, the strategy to adopt is less clear

@return an array vid position -> # visits
"""
function MHRW(g::AbstractGraph{T},n_steps::UInt64,in_degrees::Array{T,1},out_degrees::Array{T,1},starting_v::T=convert(T,1)) where {T<:Unsigned}
	visited_nodes =  T[]
	v = starting_v
	push!(visited_nodes,v)
	rejected = 0
	for i in 1:n_steps
		nei = outneighbors(g,v)
		if length(nei) > 1
			# choose a child uniformly at random
			rpos = rand(1:length(nei))
			pv = nei[rpos]
			# get the probability to move
			v_stat = out_degrees[v]
			pv_stat = out_degrees[pv]
			mp = minimum([1.,v_stat/pv_stat])
			if rand() <= mp
				push!(visited_nodes,pv)
				v = pv
			else
				rejected += 1
			end
		else
			pv = nei[1]
			push!(visited_nodes,pv)
			v = pv
		end
	end
	return visited_nodes, rejected
end

"""
    MHRW_flying(g::AbstractGraph{T},n_steps::UInt64,in_degrees::Array{T,1},out_degrees::Array{T,1},starting_v::T=convert(T,1)) where {T<:Unsigned}

Metropolis-Hasting Random Walk (flying mode)
@return an array vid position -> # visits
"""
function MHRW_flying(g::AbstractGraph{T},n_steps::UInt64,in_degrees::Array{T,1},out_degrees::Array{T,1},starting_v::T=convert(T,1)) where {T<:Unsigned}
	visited_nodes =  T[]
	v = starting_v
	push!(visited_nodes,v)
	for i in 1:n_steps
		nei = outneighbors(g,v)
		if length(nei) > 1
			# compute scores of neighbors
			in_v = in_degrees[v]
			out_v = out_degrees[v]
			#scores = stochastic([(in_v+out_v)/(in_degrees[i]+out_degrees[i]) for i in nei])
			scores = stochastic([sqrt(in_v*out_v)/sqrt(in_degrees[i]*out_degrees[i]) for i in nei])
			#scores = stochastic([out_v/out_degrees[i] for i in nei])
			pos = get_flying_index(scores)
			pv = nei[pos]
		else
			pv = nei[1]
		end
		push!(visited_nodes,pv)
		v = pv
	end
	return visited_nodes
end

"""
    CC_MHRW_flying(g::AbstractGraph{T},n_steps::UInt64,ccs::Array{Float64,1},starting_v::T=convert(T,1)) where {T<:Unsigned}

RW in flying mode guided by the colink (C1C) or clustering (C2C) coefficient

Twitter dataset: nchains: 100, burning_time: 100, nsteps: 5000
"""
function CC_MHRW_flying(g::AbstractGraph{T},n_steps::UInt64,ccs::Array{Float64,1},starting_v::T=convert(T,1)) where {T<:Unsigned}
	visited_nodes =  T[]
	v = starting_v
	push!(visited_nodes,v)
	visited_distr = Dict{T,T}()
	visited_distr[v] = 1
	for i in 1:n_steps
		nei = outneighbors(g,v)
		if length(nei) > 1
			# compute scores of neighbors
			scores = Float64[]
			for i in nei
				if ccs[i] != 0
					push!(scores,ccs[v]/ccs[i])
				else
					# null number of colinks: this is an interesting direction
					push!(scores,-1)
				end
			end
			if maximum(scores) > 0
				# get the maximum value
				# NB: mv could be set to the minimum, 0 or 1
				mv = maximum(scores)
				# set -1 entries to mv value
				scores = map(x->x==-1 ? mv : x, scores)
				# normalize vector
				scores = stochastic(scores)
			# all entries are equal to 0
			else
				# set all entries to the same value
				scores = float64(ones(length(nei)))/length(nei)
			end
			pos = get_flying_index(scores)
			pv = nei[pos]
		else
			pv = nei[1]
		end
		push!(visited_nodes,pv)
		if haskey(visited_distr,pv)
			visited_distr[pv] += 1
		else
			visited_distr[pv] = 1
		end
		v = pv
	end
	return visited_nodes,visited_distr
end

""" 
    get_CC_MHRW_flying_ball(g::AbstractGraph{T},n_vertices::T,ccs::Array{Float64,1},starting_v::T=convert(T,1),jumping_constant::Float64=0.) where {T<:Unsigned}

Explore a ball centered around starting vertex with a modified random walk
"""
function get_CC_MHRW_flying_ball(g::AbstractGraph{T},n_vertices::T,ccs::Array{Float64,1},starting_v::T=convert(T,1),jumping_constant::Float64=0.) where {T<:Unsigned}
	visited_nodes =  Set{T}()
	v = starting_v
	push!(visited_nodes, v)
	while length(visited_nodes) < n_vertices
		nei = outneighbors(g,v)
		#@debug("# children $v: ", length(nei))
		diff = collect(setdiff(Set(nei),visited_nodes))
		#@debug("# unvisited children $v: ", length(diff))
		if length(diff) > 1
			# compute scores of neighbors
			scores = Float64[]
			for i in diff
				if ccs[i] != 0
					push!(scores,ccs[v]/ccs[i])
				else
					# to be set to the max value
					push!(scores,-1.)
				end
			end
			if maximum(scores) > 0
				# get the maximum value
				# NB: mv could be set to the minimum, 0 or 1
				mv = maximum(scores)
				# set -1 entries to mv value
				scores = map(x->x==-1 ? mv : x, scores)
			else
				# all entries are equal
				scores = float64(ones(length(diff)))/length(diff)
			end
			# normalize score vector
			scores = stochastic(scores)
			scores = scores/(1-jumping_constant)
			# add starting_v as a virtual child
			push!(scores,jumping_constant)
			pos = get_flying_index(scores)
			# if the last pos was selected
			if pos == length(scores)
				v = starting_v
			else
				v = diff[pos]
			end
		elseif length(diff) == 1
			# check if a restart is needed, or if stuck at the starting node, move forward
			#if rand() > jumping_constant || v == starting_v
			if rand() > jumping_constant
				v = diff[1]
			else
				v = starting_v
			end
		else
			vns = Float64[]
			avn = collect(visited_nodes)
			for u in avn
				if ccs[u] != 0
					push!(vns,1/ccs[u])
				else
					# to be set to the max value
					push!(vns,-1.)
				end
			end
			if maximum(vns) > 0
				# get the maximum value
				# NB: mv could be set to the minimum, 0 or 1
				mv = maximum(vns)
				# set -1 entries to mv value
				vns = map(x->x==-1 ? mv : x, vns)
			else
				# all entries are equal
				vns = float64(ones(length(vns)))/length(vns)
			end
			vns = stochastic(vns)
			pos = get_flying_index(vns)
			v = avn[pos]
		end
		#@debug("Adding node: $v")
		#@debug("Total # nodes: ", length(visited_nodes))
		#@debug("---")
		push!(visited_nodes,v)
	end
	return collect(visited_nodes)
end

end # module RandomWalks
