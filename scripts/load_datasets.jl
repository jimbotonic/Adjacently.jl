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

using Pkg, Match
Pkg.activate(normpath(joinpath(@__DIR__, "..")))

using Adjacently
using Adjacently.IO: load_adjacency_list_from_csv, load_graph_from_pajek
using Adjacently.Graph: get_core, get_reverse_graph, get_basic_stats
using Adjacently.MGS: write_mgs3_graph, write_mgs3_huffman_graph
using LightGraphs: nv, ne

#Logging.configure(level=Debug)

using GraphPlot

function manage_dataset(input_path::AbstractString, output_filename::AbstractString; is_pajek=false, separator::AbstractChar=',')
	if !is_pajek
		g = load_adjacency_list_from_csv(input_path, separator)
	else
	    g = load_graph_from_pajek(input_path)
	end

	# display basic stats
	@info("Full graph #v:", convert(Int, nv(g)))
	@info("Full graph #e:", ne(g))

	@info("getting core")
	core,oni,noi = get_core(g)
	@info("Core #v:", convert(Int, nv(core)))
	@info("Core #e:", ne(core))

	@info("getting reverse core")
	rcore = get_reverse_graph(core) 
	@info("RCore #e:", ne(rcore))

	# export graph data of core and reverse core
	write_mgs3_graph(core, output_filename)
	write_mgs3_huffman_graph(core, rcore, output_filename)
	#serialize_to_jld(core, "core", output_filename)
	
	write_mgs3_graph(rcore, replace(output_filename, "core" => "rcore"))
	write_mgs3_huffman_graph(rcore, core, replace(output_filename, "core" => "rcore"))
	#serialize_to_jld(rcore, "rcore", output_filename)
end

###
# loading and exporting datasets
###

# amazon_0601, web_google, arxiv_hep-ph, eat
if length(ARGS) == 0 || ARGS[1] == "-h"
    println("Usage: julia load_datasets.jl [dataset]")
    println("Available datasets:")
    println("  amazon - Amazon_0601 graph")
    println("  google - Web_Google graph")
    println("  arxiv  - Arxiv_HEP-PH graph")
    println("  eat    - EAT (new) graph")
    exit(0)
end

dataset = ARGS[1]

@match dataset begin
    "amazon" => begin
        @info("loading Amazon_0601 graph")
        manage_dataset("../datasets/Amazon_0601/Amazon0601.txt", "Amazon_0601_core", separator='\t')
    end
    "google" => begin
        @info("loading Web_Google graph")
        manage_dataset("../datasets/Web_Google/web-Google.txt","Web_Google_core", separator='\t')
    end
    "arxiv" => begin
        @info("loading Arxiv_HEP-PH graph")
        manage_dataset("../datasets/Arxiv_HEP-PH/Cit-HepPh.txt", "Arxiv_HEP-PH_core", separator='\t')
    end
    "eat" => begin
        @info("loading EAT (new) graph")
        manage_dataset("../datasets/EAT/EATnew.net", "EAT_core", is_pajek=true)
    end
end
