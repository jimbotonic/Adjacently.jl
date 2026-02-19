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
