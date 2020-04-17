#
# Adjacently: Julia Complex Networks Library
# Copyright (C) 2016-2020  Jimmy Dubuisson <jimmy.dubuisson@gmail.com>
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

using Graphs: adjlist, num_vertices
using Distances: chebyshev

using Adjacently.io: load_mgs3_graph
using Adjacently.pr: PR, PPR
using Adjacently.graph: get_reverse_graph, get_sparse_P_matrix

# load graph in MGSv3 format
core = adjlist(UInt32, is_directed=true)
load_mgs3_graph(core, "../datasets/Arxiv_HEP-PH/Arxiv_HEP-PH_core.mgs")
#load_mgs4_graph(core, "../datasets/Arxiv_HEP-PH/Arxiv_HEP-PH_core.mgz")

@info("getting rcore")
rcore = get_reverse_graph(core)

# constants
n = num_vertices(core)
damping = .85
epsilon = 1e-8

@info("computing Pagerank of core and rcore")
@time pr_core = PR(core, rcore, epsilon=epsilon)
#@time pr_rcore = PR(rcore, core, epsilon=epsilon)

@info("computing personalized Pageranks for nodes 1 and 100")
# source and target node
s = convert(UInt32,1)
t = convert(UInt32,100)
@time pr_s = PPR(s, core, rcore, epsilon=epsilon)
#@time pr_t = PPR(t, rcore, core, epsilon=epsilon)

# compute Monte Carlo PR
#niter = 400
#@info("compute Pagerank (Monte Carlo)")
#@time pr_mc = PR(core, niter)
#
#@info("pr_core <-> pr_mc: ", chebyshev(pr_core, pr_mc))

# matrix computation
@info("getting P matrix")
@time P = get_sparse_P_matrix(core)

@info("compute Pagerank (power iteration)")
ppr = zeros(Float64,n)
ppr[s] = 1.
@time pr_pi = PR(P, epsilon=epsilon)

@info("computing personalized Pageranks for node 1 (power iteration)")
@time pr_pi_s = PR(P, ppr=ppr, epsilon=epsilon)

@info("pr_core <-> pr_pi: ", chebyshev(pr_core, pr_pi))
@info("pr_s <-> pr_pi_s: ", chebyshev(pr_s, pr_pi_s))
