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

module Adjacently

module algo

export tarjan, pearce, pearce_iterative

include("algo.jl")
end

module io

export load_jls_serialized, serialize_to_jls, load_jld_serialized, 
    serialize_to_jld, write_mgs3_graph, load_mgs3_graph, 
    write_mgs3_huffman_graph, load_mgs3_huffman_graph, load_adjacency_list_from_csv, 
    load_graph_from_pajek

include("io.jl")
end

module pr

export PR, PPR

include("pr.jl")
end

module util

export quicksort_iterative_permutation!, bottom_up_sort, binary_search, huffman_encoding

include("util.jl")
end

module rw

export RW, RW_aggregated, RW_aggregated, 
    US, ARW, ARW_flying, MHRW, MHRW_flying, CC_MHRW_flying

include("rw.jl")
end

module graph

export get_basic_stats, display_basic_stats, get_out_degree_stats, 
    get_sinks, get_sources, subgraph, 
    subgraph_streamed, get_core, get_core_streamed, 
    get_reverse_graph, get_vertex_in_degrees, get_in_out_degrees, 
    get_avg_out_degree, get_forward_ball, get_clustering_coefficients, 
    get_colink_coefficients, get_inclist_from_adjlist, get_sparse_adj_matrix, 
    get_sparse_P_matrix

include("graph.jl")
end

module cycles

export list_colinks, list_triangles, load_colinks_distribution, 
    load_triangles_distribution

include("cycles.jl")
end

end
