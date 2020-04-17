#
# Adjacently: Julia Complex Networks Library
# Copyright (C) 2016-2020 Jimmy Dubuisson <jimmy.dubuisson@gmail.com>
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

include("../util.jl")
include("../io.jl")
include("../graph.jl")
include("../pr.jl")
include("../rw.jl")

function get_sp_greedy(s::T, t::T, g::GenericAdjacencyList{T,Array{T,1},Array{Array{T,1},1}}, df::Array{Float64,1}, max_length::Int64) where {T<:Unsigned}
	@info("searching shortest path (greedy approach) between vertices $s and $t")
	# current vertex
	cv = s
	# current path
	sp = T[]
	while length(sp) < max_length
		push!(sp,cv)
		nei = out_neighbors(cv,g)
		nnei = setdiff(nei,sp)
		if length(nnei) == 0
			@info("--- exploration reached a dead end")
			@info("--- explored path: ", sp)
			break
		end
		cv = nnei[findmax(df[nnei])[2]]

		# a path was found
		if cv == t
			push!(sp,t)
			@info("path found: ", sp)
			break
		end
	end
	return sp
end

function get_sp_proba(s::T, t::T, g::GenericAdjacencyList{T,Array{T,1},Array{Array{T,1},1}}, df::Array{Float64,1}, max_length::Int64, max_iter::Int64) where {T<:Unsigned}
	@info("searching shortest path (probabilistic approach) between vertices $s and $t")
	# array of paths
	sps = Array{Array{T,1},1}()
	# current path
	sp = T[]
	# search in a probabilistic way shortest paths between s and t
	cv = s
	c = 0
	# shortest length found so far
	min_length = Inf
	while c < max_iter
		sp = T[]
		cv = s
		# path length
		pl = 1
		while pl < max_length && pl < (min_length-1)
			push!(sp,cv)
			nei = out_neighbors(cv, g)
			nnei = setdiff(nei,sp)
			if length(nnei) == 0
				@debug("--- exploration reached a dead end")
				@debug("--- explored path: ", sp)
				break
			end
			pos = get_flying_index(df[nnei] / sum(df[nnei]))
			cv = nnei[pos]

			# we found a path
			if cv == t
				push!(sp,t)
				pl += 1
				if !(sp in sps)
					@info("path found between $s and $t (length $pl): ", sp)
					push!(sps,sp)
					# update min_length if necessary
					if pl < min_length
						min_length = pl
					end
				end
				break
			end
			pl += 1
		end
		c += 1
	end
	return sps
end

BOOST = true
PROBA = false

@info("############ Testing shortest paths functions")

core = adjlist(UInt32, is_directed=true)
load_mgs3_graph(core, "../datasets/Arxiv_HEP-PH/Arxiv_HEP-PH_core.mgs")

@info("getting rcore")
rcore = get_reverse_graph(core)

@info("transforming core from adj to inc list")
core2 = get_inclist_from_adjlist(core)

@info("computing Djikstra from vertex 1")
dists = ones(num_edges(core2))
s = convert(UInt32,1)
n = num_vertices(core)
@time r = dijkstra_shortest_paths(core2, dists, s)

@info("getting P for core")
P_core = get_sparse_P_matrix(core)
P_rcore = get_sparse_P_matrix(rcore)

@info("computing Pagerank of core and rcore")
EPSILON = 1e-8
@time pr_core = PR(P_core, epsilon=EPSILON)
@time pr_rcore = PR(P_rcore, epsilon=EPSILON)

pr_boost = Float64[]
if BOOST
	pr_boost = pr_core .* pr_rcore
end

@info("computing personalized Pageranks for nodes s and t")
ppr_s = zeros(Float64,n)
ppr_s[s] = 1.
@time pr_s = PR(P_core, ppr=ppr_s, init_pr=pr_core, epsilon=EPSILON)

# max paths length
#ml = floor(Int,log(length(vertices(core))))
ml = 20

sp = UInt32[]

# counters
tc = 0
sc = 0
sc2 = 0

for t in UInt32(2):UInt32(100)
	# Djikstra
	paths = enumerate_paths(vertices(core2), r.parent_indices, Int64[t])

	ppr_t = zeros(Float64,n)
	ppr_t[t] = 1.
	@time pr_t = PR(P_rcore, ppr=ppr_t, init_pr=pr_rcore, epsilon=EPSILON)

	# DF
	df = pr_s .* pr_t

	if BOOST
		df = df ./ pr_boost
	end

	if PROBA
		sps = get_sp_proba(s, t, core, df, ml, convert(Int64,1e3))
	else
		sp = get_sp_greedy(s, t, core, df, ml)
	end
	
	@info("----- sp (dj): ", convert(Array{Int64,1}, paths[1]))

	if PROBA
		if length(sps) > 0
			mi = findmin(Int64[length(a) for a in sps])[2]
			sp = sps[mi]
		else
			sp = UInt32[]
		end
	end
	@info("----- sp (df): ", convert(Array{Int64,1}, sp))

	if paths[1] == sp
		sc += 1
	end
	if length(paths[1]) == length(sp)
		sc2 += 1
	end
	tc += 1
end

@info("% same SP: ", sc/tc)
@info("% same length: ", sc2/tc)

