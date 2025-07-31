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

using Random
using LightGraphs
using ..Graph: get_in_degrees, get_out_degrees, get_in_out_degrees

export get_degree_entropy, get_entropy, powerlaw_sample

"""
    get_degree_entropy(g::AbstractGraph{T}, type::Symbol=:in_degree) where {T<:Unsigned}

Compute the entropy of the degree distribution of the graph in bits.

@param g::AbstractGraph{T}: the graph
@param type::Symbol: the type of degree distribution to compute (:in_degree or :out_degree)
@return::Float64: the entropy of the degree distribution in bits
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

"""
    get_entropy(samples::Vector{T}) where {T<:Unsigned}

Compute the entropy of the samples.

@param samples::Vector{T}: the samples
@return::Float64: the entropy of the samples
"""
function get_entropy(samples::Vector{T}) where {T<:Unsigned}
    # compute the entropy of the distribution using Shannon entropy formula
    entropy = 0.0
    # sample distribution (value => count)
    sample_dist = Dict{T, Float64}()
    total_sample = length(samples)
    for sample in samples
        sample_dist[sample] = get(sample_dist, sample, 0.0) + 1.0
    end
    for sample in keys(sample_dist)
        probability = sample_dist[sample] / total_sample
        entropy += probability * log2(probability)
    end
    return -entropy
end

"""
    powerlaw_sample(k::Float64, n::Int, min_val::Int=1, max_val::Int, T::Type=UInt8) where {T<:Unsigned}

Generate `n` values of type `T` drawn from a power-law distribution P(x) ~ x^(-k),
with support in [min_val, max_val]. Use inverse transform sampling.

@param k::Float64: the exponent of the power-law distribution
@param n::Int: the number of values to generate
@param min_val::Int: the minimum value of the distribution
@param max_val::Int: the maximum value of the distribution
@param T::Type: the type of the values to generate
@return::Vector{T}: the generated values
"""
function powerlaw_sample(k::Float64, n::Int, min_val::Int, max_val::Int, ::Type{T}=UInt8) where {T<:Unsigned}
    @assert min_val >= 1 "min_val must be >= 1"
    @assert max_val >= min_val "max_val must be >= min_val"
    @assert k > 1.0 "k must be > 1 for normalization"

    # precompute the cumulative distribution function (CDF)
    values = min_val:max_val
    probs = cumsum([x^(-k) for x in values])
    # normalize to [0, 1]
    cdf = probs ./ last(probs)

    # inverse transform sampling: draw n uniform samples and binary search in CDF
    rand_vals = rand(n)
    results = Vector{T}(undef, n)
    for i in 1:n
        # binary search in CDF
        idx = searchsortedfirst(cdf, rand_vals[i])
        # convert to type T
        results[i] = T(values[idx])
    end
    return results
end



end # module Distribution