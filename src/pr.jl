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

include("rw.jl")

using LightGraphs, DataStructures, Logging, Distances, SparseArrays

DAMPING_FACTOR = 0.85
EPSILON = 1e-6
MAX_ITER = 100

###
# Simple graph based Pagerank algorithm
###

""" 
    PR(g::AbstractGraph{T},rg::AbstractGraph{T};init_pr::Array{Float64,1}=Float64[],damping::Float64=DAMPING_FACTOR, epsilon::Float64=EPSILON) where {T<:Unsigned}

Naive inefficient implementation of Pagerank algorithm
"""
function PR(g::AbstractGraph{T},rg::AbstractGraph{T};init_pr::Array{Float64,1}=Float64[],damping::Float64=DAMPING_FACTOR, epsilon::Float64=EPSILON) where {T<:Unsigned}
	vs = vertices(g)
	n = nv(g)
	@info("computing Pagerank (size of graph $n)")
	# initialize pagerank vector
	if length(init_pr) == n
		pr = init_pr
	else
		pr = Float64[1/n for k in vs]
	end
	pr2 = zeros(Float64,n)
	while true
		for v in vs
			nv = 0.
			# get v children in the reverse graph
			in_nei = outneighbors(rg,v)
			if length(in_nei) > 0
				for p in in_nei
					nv +=  pr[p]/length(outneighbors(g,p))
				end
			end
			pr2[v] = (1-damping)/n+damping*nv
		end
		d = chebyshev(pr,pr2)
		if d <= epsilon
			pr = pr2
			break
		end
		pr = copy(pr2)
	end
	return pr
end

""" 
    PR(rg::AbstractGraph{T},out_degrees::Array{T,1};init_pr::Array{Float64,1}=Float64[],damping::Float64=DAMPING_FACTOR,epsilon::Float64=EPSILON,save_pr::Bool=False) where {T<:Unsigned}

Naive inefficient implementation of Pagerank algorithm (designed for large graphs)
"""
function PR(rg::AbstractGraph{T},out_degrees::Array{T,1};init_pr::Array{Float64,1}=Float64[],damping::Float64=DAMPING_FACTOR,epsilon::Float64=EPSILON,save_pr::Bool=False) where {T<:Unsigned}
	vs = vertices(rg)
	n = nv(g)
	@info("computing Pagerank (size of graph $n)")
	# initialize pagerank vector
	if length(init_pr) == n
		pr = init_pr
	else
		pr = Float64[1/n for k in vs]
	end
	pr2 = zeros(Float64, n)
	# iteration number
	ic = 1
	# progression
	thp = ceil(Int,n/100)
	while true
		for v in vs
			nv = 0.
			# get v children in the reverse graph
			in_nei = outneighbors(rg,v)
			if length(in_nei) > 0
				for p in in_nei
					nv +=  pr[p]/out_degrees[p]
				end
			end
			pr2[v] = (1-damping)/n+damping*nv
			if v % thp == 0
				# log progression
				@info("iteration $ic: ", (v/n)*100, " %")
			end
		end
		d = chebyshev(pr,pr2)
		@info("distance(t,t+1) = $d")
		if save_pr
			serialize_to_file(pr, "pr-iter-$ic.jld")
		end
		ic = ic + 1
		if d <= epsilon
			pr = pr2
			break
		end
		pr = copy(pr2)
	end
	return pr
end

""" 
    PPR(src::T,g::AbstractGraph{T},rg::AbstractGraph{T};damping::Float64=DAMPING_FACTOR,epsilon::Float64=EPSILON) where {T<:Unsigned}

Naive inefficient implementation of personalized Pagerank with a single source vertex
"""
function PPR(src::T,g::AbstractGraph{T},rg::AbstractGraph{T};damping::Float64=DAMPING_FACTOR,epsilon::Float64=EPSILON) where {T<:Unsigned}
	vs = vertices(g)
	n = nv(g)
	@info("computing personalized Pagerank (size of graph $n, source $src)")
	# initialize pagerank vector
	pr = [0. for k in vs]
	pr[src] = 1.
	pr2 = zeros(Float64, n)
	while true
		for v in vs
			nv = 0.
			# get v children in the reverse graph
			in_nei = outneighbors(rg,v)
			if length(in_nei) > 0
				for p in in_nei
					nv +=  pr[p]/length(outneighbors(g,p))
				end
			end
			if v == src
				pr2[v] = (1-damping)+damping*nv
			else
				pr2[v] = damping*nv
			end
		end
		d = chebyshev(pr,pr2)
		if d <= epsilon
			pr = pr2
			break
		end
		pr = copy(pr2)
	end
	return pr
end

###
# Monte-Carlo Pagerank algorithm
###

""" 
    PR(g::AbstractGraph{T},n_cycles::Int;damping::Float64=DAMPING_FACTOR,epsilon::Float64=EPSILON) where {T<:Unsigned}

MC Pagerank with cyclic start of complete path stopping at sink nodes
http://www-sop.inria.fr/members/Konstantin.Avratchenkov/pubs/mc.pdf
"""
function PR(g::AbstractGraph{T},n_cycles::Int;damping::Float64=DAMPING_FACTOR,epsilon::Float64=EPSILON) where {T<:Unsigned}
	vs = vertices(g)
	n = nv(g)
	vv = zeros(Float64, n)
	@info("computing Monte-Carlo Pagerank (size of graph $n)")
	for i in 1:n_cycles
		for v in vs
			vv += RW_aggregated(g, (1-damping), v)
		end
	end
	return vv/sum(vv)
end

###
# Power iteration method
###

"""
    PR(P::SparseMatrixCSC{Float64,T};ppr::Array{Float64,1}=Float64[],init_pr::Array{Float64,1}=Float64[], damping::Float64=DAMPING_FACTOR, epsilon::Float64=EPSILON, max_iter::Int=MAX_ITER) where {T<:Unsigned}

Compute Pagerank (Power iteration method)

P: sparse adjacency matrix
ppr: personalized distribution of random jumps
init_pr: initial Pagerank distribution
"""
function PR(P::SparseMatrixCSC{Float64,T};ppr::Array{Float64,1}=Float64[],init_pr::Array{Float64,1}=Float64[], damping::Float64=DAMPING_FACTOR, epsilon::Float64=EPSILON, max_iter::Int=MAX_ITER) where {T<:Unsigned}
	n = size(P)[1]
	# initialize personalized vector
	if length(ppr) != n
		ppr = Float64[1/n for i in 1:n]
	end
	# initialize Pagerank vector
	if length(init_pr) != n
		pr = Float64[1/n for i in 1:n]
	else
		pr = init_pr
	end
	pr2 = zeros(Float64, n)
    i = 0
	while i < max_iter
		#pr2 = (damping*pr'*P + (1-damping)*ppr')'
		pr2 = damping*P'*pr + (1-damping)*ppr
		d = chebyshev(pr,pr2)
		if d <= epsilon
			pr = pr2
			break
		end
		pr = copy(pr2)
        i += 1
	end
	return pr
end

""" 
    PR(P::SparseMatrixCSC{Float64,T},fun::Function;ppr::Array{Float64,1}=Float64[],init_pr::Array{Float64,1}=Float64[], damping::Float64=DAMPING_FACTOR, epsilon::Float64=EPSILON, max_iter::Int=MAX_ITER) where {T<:Unsigned}

Compute non-linear Pagerank (Power iteration method)

P: sparse adjacency matrix
fun: activation function
ppr: personalized distribution of random jumps
init_pr: initial Pagerank distribution
damping: damping factor 
epsilon: 
max_iter: maximum number of iterations
"""
function PR(P::SparseMatrixCSC{Float64,T},fun::Function;ppr::Array{Float64,1}=Float64[],init_pr::Array{Float64,1}=Float64[], damping::Float64=DAMPING_FACTOR, epsilon::Float64=EPSILON, max_iter::Int=MAX_ITER) where {T<:Unsigned}
	n = size(P)[1]
	# initialize personalized vector if needed
	if length(ppr) != n
		ppr = Float64[1/n for i in 1:n]
	end
	# initialize Pagerank vector
	if length(init_pr) != n
		pr = Float64[1/n for i in 1:n]
	else
		pr = init_pr
	end
	pr2 = zeros(Float64, n)
    i = 0
	while i < max_iter
		pr2 = fun.(damping*P'*pr + (1-damping)*ppr)
		d = chebyshev(pr,pr2)
		if d <= epsilon
			pr = pr2
			break
		end
		pr = copy(pr2)
        i += 1
	end
	return pr
end
