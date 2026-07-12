#
# Adjacently: Julia Complex Directed Networks Library
# Copyright (C) 2016-2026 Anonymous (double-blind review)
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

module InterEncoding

import ...Compression
using LightGraphs: outneighbors, is_directed
using ....IO: BitWriter

export encode_ab_pairlists,
       cluster_density,
       should_stop_coarsening_by_density

"""
    encode_ab_pairlists(w, A, neighbors_in_B, count_encoding, integer_encoding)

Encode an (A,B) inter-block by writing, for each u in A (cluster-local order),
the count of neighbors in B followed by a delta(Fibonacci)-encoded list of
cluster-local indices into B.
"""
function encode_ab_pairlists(w::BitWriter, A::Vector{T}, neighbors_in_B::Dict{T,Vector{Int}}, count_encoding::Symbol, integer_encoding::Symbol) where {T<:Unsigned}
    for u in A
        Ns = get(neighbors_in_B, u, Int[])
        Compression.write_small_count(w, T(length(Ns)), count_encoding)
        if !isempty(Ns)
            Compression.write_delta(w, T.(Ns), integer_encoding)
        end
    end
end

"""
    cluster_density(g, C::Vector{T}) where {T<:Unsigned}

Compute the internal density of the cluster induced by C.
- Directed graphs: m_C / (|C| * (|C| - 1)) where m_C counts directed arcs within C.
- Undirected graphs: m_Cu / (|C| * (|C| - 1) / 2) where m_Cu counts undirected pairs within C once.
"""
function cluster_density(g, C::Vector{T}) where {T<:Unsigned}
    n = length(C)
    n <= 1 && return 0.0
    s = Set{Int}(Int.(C))
    directed = is_directed(g)
    m = 0
    for u in Int.(C)
        for v in outneighbors(g, u)
            vi = Int(v)
            if vi in s
                if directed
                    m += 1  # count arcs
                else
                    # count each undirected pair once using u < v gate
                    if u < vi
                        m += 1
                    end
                end
            end
        end
    end
    denom = directed ? (n * (n - 1)) : (n * (n - 1) / 2)
    return denom > 0 ? (m / denom) : 0.0
end

"""
    should_stop_coarsening_by_density(g, clusters, min_density::Float64)

Return true if all cluster densities are >= min_density.
"""
function should_stop_coarsening_by_density(g, clusters::Vector{Vector{T}}, min_density::Float64) where {T<:Unsigned}
    for C in clusters
        if cluster_density(g, C) < min_density
            return false
        end
    end
    return true
end

end # module InterEncoding
