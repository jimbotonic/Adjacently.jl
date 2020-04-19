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

include("../src/util.jl")
include("../src/io.jl")
include("../src/graph.jl")

using Graphs
using LightGraphs

const G = Graphs
const LG = LightGraphs

#Logging.configure(level=INFO)

###
# loading and exporting datasets
###

# amazon_0601, web_google, arxiv_hep-ph, eat
dataset = "eat"

if dataset == "amazon"
	@info("loading Amazon_0601 graph")
	g = adjlist(UInt32, is_directed=true)
	load_adjacency_list_from_csv(UInt32, g, "../datasets/Amazon_0601/Amazon0601.txt")
	@info("# vertices:", length(G.vertices(g)))
	g = adjlist(UInt32, is_directed=true)
	@info("# edges:", num_edges(g))

	@info("getting core")
	core,oni,noi = get_core(g)
	@info("# vertices:", length(G.vertices(core)))
	@info("# edges:", num_edges(core))

	@info("getting reverse graph")
	rcore = get_reverse_graph(core) 
	@info("# edges (rcore):", num_edges(rcore))

	write_mgs3_graph(core, "Amazon_0601_core")
	write_mgs4_graph(core, rcore, "Amazon_0601_core")
	#serialize_to_jld(core, "core", "Amazon_0601_core")
	
	write_mgs3_graph(rcore, "Amazon_0601_rcore")
	write_mgs4_graph(rcore, core, "Amazon_0601_rcore")
	#serialize_to_jld(rcore, "rcore", "Amazon_0601_rcore")
elseif dataset == "google"
	@info("loading Web_Google graph")
	g = adjlist(UInt32, is_directed=true)
	load_adjacency_list_from_csv(UInt32, g, "../datasets/Web_Google/web-Google.txt")
	@info("# vertices:", length(G.vertices(g)))
	@info("# edges:", num_edges(g))

	@info("getting core")
	core,oni,noi = get_core(g)
	@info("# vertices:", length(G.vertices(core)))
	@info("# edges:", num_edges(core))

	@info("getting reverse graph")
	rcore = get_reverse_graph(core) 
	@info("# edges (rcore):", num_edges(rcore))

	write_mgs3_graph(core, "Web_Google_core")
	write_mgs4_graph(core, rcore, "Web_Google_core")
	#serialize_to_jld(core, "core", "Web_Google_core")
	
	write_mgs3_graph(rcore, "Web_Google_rcore")
	write_mgs4_graph(rcore, core, "Web_Google_rcore")
	#serialize_to_jld(rcore, "rcore", "Web_Google_rcore")
elseif dataset == "arxiv"
	@info("loading Arxiv_HEP-PH graph")
	g = adjlist(UInt32, is_directed=true)
	load_adjacency_list_from_csv(UInt32, g, "../datasets/Arxiv_HEP-PH/Cit-HepPh.txt")
	@info("# vertices:", length(G.vertices(g)))
	@info("# edges:", num_edges(g))

	@info("getting core")
	core,oni,noi = get_core(g)
	@info("# vertices:", length(G.vertices(core)))
	@info("# edges:", num_edges(core))

	@info("getting reverse graph")
	rcore = get_reverse_graph(core) 
	@info("# edges (rcore):", num_edges(rcore))

	write_mgs3_graph(core, "Arxiv_HEP-PH_core")
	write_mgs4_graph(core, rcore, "Arxiv_HEP-PH_core")
	#serialize_to_jld(core, "core", "Arxiv_HEP-PH_core")
	
	write_mgs3_graph(rcore, "Arxiv_HEP-PH_rcore")
	write_mgs4_graph(rcore, core, "Arxiv_HEP-PH_rcore")
	#serialize_to_jld(rcore, "rcore", "Arxiv_HEP-PH_rcore")
elseif dataset == "eat"
	@info("loading EAT (new) graph")
	g = adjlist(UInt32, is_directed=true)
	load_graph_from_pajek(UInt32, g, "../datasets/EAT/EATnew.net")
	@info("# vertices:", length(G.vertices(g)))
	@info("# edges:", num_edges(g))

	@info("getting core")
	core,oni,noi = get_core(g)
	@info("# vertices:", length(G.vertices(core)))
	@info("# edges:", num_edges(core))

	@info("getting reverse graph")
	rcore = get_reverse_graph(core) 
	@info("# edges (rcore):", num_edges(rcore))

	write_mgs3_graph(core, "EAT_core")
	write_mgs4_graph(core, rcore, "EAT_core")
	#serialize_to_jld(core, "core", "EAT_core")
	
	write_mgs3_graph(rcore, "EAT_rcore")
	write_mgs4_graph(rcore, core, "EAT_rcore")
	#serialize_to_jld(rcore, "rcore", "EAT_rcore")
end
