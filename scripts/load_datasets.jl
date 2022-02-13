#
# Adjacently: Julia Complex Networks Library
# Copyright (C) 2016-2022 Jimmy Dubuisson <jimmy.dubuisson@gmail.com>
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

#Logging.configure(level=Debug)

using GraphPlot

function load_dataset(input_path::AbstractString,output_filename::AbstractString;is_pajek=false) where {T<:Unsigned}
	g = SimpleDiGraph{UInt32}()
	if !is_pajek
		load_adjacency_list_from_csv(UInt32, g, input_path)
	else
	    load_graph_from_pajek(UInt32, g, input_path)
	end
	@info("# vertices:", convert(Int,nv(g)))
	@info("# edges:", ne(g))

	@info("getting core")
	core,oni,noi = get_core(g)
	@info("# vertices:", convert(Int,nv(core)))
	@info("# edges:", ne(core))

	@info("getting reverse graph")
	rcore = get_reverse_graph(core) 
	@info("# edges (rcore):", ne(rcore))

	write_mgs3_graph(core, output_filename)
	write_mgs4_graph(core, rcore, output_filename)
	#serialize_to_jld(core, "core", output_filename)
	
	write_mgs3_graph(rcore, output_filename)
	write_mgs4_graph(rcore, core, output_filename)
	#serialize_to_jld(rcore, "rcore", output_filename)
end

###
# loading and exporting datasets
###

# amazon_0601, web_google, arxiv_hep-ph, eat
dataset = "eat"

if dataset == "amazon"
	@info("loading Amazon_0601 graph")
	load_dataset("../datasets/Amazon_0601/Amazon0601.txt", "Amazon_0601_core")
elseif dataset == "google"
	@info("loading Web_Google graph")
	load_dataset("../datasets/Web_Google/web-Google.txt","Web_Google_core")
elseif dataset == "arxiv"
	@info("loading Arxiv_HEP-PH graph")
	load_dataset("../datasets/Arxiv_HEP-PH/Cit-HepPh.txt", "Arxiv_HEP-PH_core")
elseif dataset == "eat"
	@info("loading EAT (new) graph")
	load_dataset("../datasets/EAT/EATnew.net", "EAT_rcore", is_pajek=true)
end
