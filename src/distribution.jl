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

module Distribution

using LightGraphs
using ..Graph: get_in_degrees, get_out_degrees, get_in_out_degrees

export get_degree_entropy

"""
    get_degree_entropy(g::AbstractGraph{T}, type::Symbol=:in_degree) where {T<:Unsigned}

Compute the entropy of the degree distribution of the graph in bits.
"""
function get_degree_entropy(g::AbstractGraph{T}, type::Symbol=:in_degree) where {T<:Unsigned}
    if type == :in_degree
        degrees = get_in_degrees(g)
    elseif type == :out_degree
        degrees = get_out_degrees(g)
    else
        @error("Invalid distribution type: $type")
    end

    # compute the entropy of the distribution
    entropy = 0.0
    # degree distribution (degree => count)
    degree_distribution = Dict{T, Float64}()
    total_degree = 0
    for v in keys(degrees)
        degree = degrees[v]
        degree_distribution[degree] = get(degree_distribution, degree, 0.0) + 1.0
        total_degree += degree
    end
    for degree in keys(degree_distribution)
        probability = degree_distribution[degree] / total_degree
        entropy += probability * log2(probability)
    end
    return -entropy
end


end # module Distribution