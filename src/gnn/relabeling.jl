#
# Adjacently: Julia Complex Directed Networks Library
# Copyright (C) 2016-2026 Jimmy Dubuisson <jimmy@dubuisson.ch>
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

"""
    relabel_vertices_gnn(g::AbstractGraph{T}, model::GNNModel) where {T<:Unsigned}

Relabel the vertices of the graph based on the GNN output scores.
Sorts vertices by their GNN score (ascending) and returns the mapping.

@param g: the graph
@param model: the trained GNN model

@returns the mapping of the vertices (vertex_id[old_id] -> new_id)
"""
function relabel_vertices_gnn(g::AbstractGraph{T}, model::GNNModel) where {T<:Unsigned}
    scores = gnn_forward(model)
    perm = sortperm(scores)
    
    mapping = Dict{T,T}()
    for (new_id, old_id) in enumerate(perm)
        mapping[T(old_id)] = T(new_id)
    end
    
    return mapping
end
