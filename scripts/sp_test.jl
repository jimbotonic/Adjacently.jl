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

using Statistics

# override mkindx in a_star_spath.jl
#mkindx(t) = typeof(t) == UInt32 ? t : t.index

###
# shortest distance
###

@info("############ Testing shortest paths functions")

core = adjlist(UInt32, is_directed=true)
load_mgs3_graph(core, "../datasets/Arxiv_HEP-PH/Arxiv_HEP-PH_core.mgs")
#load_mgs4_graph(core, "../datasets/Arxiv_HEP-PH/Arxiv_HEP-PH_core.mgz")

@info("getting rcore")
rcore = get_reverse_graph(core)

@info("transforming core from adj to inc list")
core2 = get_inclist_from_adjlist(core)
#serialize_to_jld(core2, "core", "Arxiv_HEP-PH_core_inclist")

@info("computing Djikstra from vertex 1")
dists = ones(num_edges(core2))
s = convert(UInt32,1)
t = convert(UInt32,100)
n = num_vertices(core)
@time r = dijkstra_shortest_paths(core2, dists, s)

@info("mean distance from vertex 1: ", mean(r.dists))
paths = enumerate_paths(vertices(core2), r.parent_indices, [10,100])
@info("shortest paths from vertex 1 to vertices 10 and 100: ", convert(Array{Array{Int64,1},1}, paths))

@info("computing Bellman-Ford from vertex 1")
@time r = bellman_ford_shortest_paths(core2, dists, UInt32[1])

@info("mean distance from vertex 1: ", mean(r.dists))
paths = enumerate_paths(vertices(core2), r.parent_indices, [10,100])
@info("shortest paths from vertex 1 to vertices 10 and 100: ", convert(Array{Array{Int64,1},1}, paths))

#@info("computing A* from node 1 to vertex 100")
#@time r = shortest_path(core2, dists, s, t)
#@info("shortest path from vertex 1 to vertex 100: ", r)

@info("getting P for core and rcore")
P_core = get_sparse_P_matrix(core)
P_rcore = get_sparse_P_matrix(rcore)

@info("computing Pagerank of core and rcore")
eps = 1e-4
#@time pr_core = PR(core, rcore, epsilon=eps)
#@time pr_rcore = PR(rcore, core, epsilon=eps)
@time pr_core = PR(P_core, epsilon=eps)
@time pr_rcore = PR(P_rcore, epsilon=eps)

@info("computing personalized Pageranks for nodes 1 and 100")
#@time pr_s = PPR(s, core, rcore, epsilon=eps)
#@time pr_t = PPR(t, rcore, core, epsilon=eps)
ppr_s = zeros(Float64,n)
ppr_s[s] = 1.
@time pr_s = PR(P_core, ppr=ppr_s, epsilon=1e-10)
ppr_t = zeros(Float64,n)
ppr_t[t] = 1.
@time pr_t = PR(P_rcore, ppr=ppr_t, epsilon=1e-10)

@info("computing diffusion fingerprints")
df = pr_s .* pr_t
pr_boost = pr_core .* pr_rcore
df_boost = df ./ pr_boost

# max length
ml = log(length(vertices(core)))
# current vertex
cv = s
# current path
sp = UInt32[]

@info("finding shortest path (greedy approach)")

# greedy approach
while length(sp) < ml
	global cv
	@debug("adding vertex: ", cv)
	# add new vertex to the path
	push!(sp,cv)
	# get neighbors not already in the path
	nei = out_neighbors(cv,core)
	nnei = setdiff(nei,sp)
	# we reached a dead end...
	if length(nnei) == 0
		@info("--- exploration reached a dead end")
		@info("--- explored path: ", sp)
		break
	end
	# find the neighbor with highest DF value
	cv = nnei[findmax(df[nnei])[2]]
	#cv = nei[findmax(df_boost[nnei])[2]]

	# we found a path
	if cv == t
		@info("path found: ", sp)
		push!(sp,t)
		break
	end
end

@info("shortest path: ", convert(Array{Int64,1}, sp))
@info("finding shortest path (probabilistic approach)")

# array of paths
sps = Array{Array{UInt32,1},1}()
# current path
sp = UInt32[]
cv = s

# search in a probabilistic way shortest paths between s and t
c = 0
max_iter = 1e3

while c < max_iter
	global c
	while length(sp) < ml 
		global sp, cv
		@debug("adding vertex: ", cv)
		push!(sp,cv)
		nei = out_neighbors(cv,core)
		nnei = setdiff(nei,sp)
		if length(nnei) == 0
			@debug("--- exploration reached a dead end")
			@debug("--- explored path: ", sp)
			sp = UInt32[]
			cv = s
			break
		end
		pos = get_flying_index(df[nnei] / sum(df[nnei]))
		cv = nnei[pos]

		# we found a path
		if cv == t
			@debug("path found: ", sp)
			push!(sp,t)
			if !(sp in sps)
				push!(sps,sp)
			end
			sp = UInt32[]
			cv = s
			break
		end
	end
	c += 1
end

@info("shortest paths: ", convert(Array{Array{Int64,1},1}, sps))


