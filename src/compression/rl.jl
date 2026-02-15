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

module RL

"""
    Compression.RL

A lightweight scaffolding for learning-based encoding decisions.
Provides helpers to extract features, estimate bit costs using the
existing Compression routines, and a greedy selector that can be
plugged into encoders for reference and mode decisions.
"""

import ..Compression
using ...IO: BitWriter, flush_bitwriter

export EnvConfig,
       FeatureVector,
       estimate_raw_bits,
       estimate_ref_bits,
       greedy_ref_select,
       choose_positions_mode

"""
    EnvConfig

Holds environment parameters relevant to decisions.
"""
Base.@kwdef struct EnvConfig
    integer_encoding::Symbol = :fibonacci
    count_encoding::Symbol = :fibonacci
    ref_window::Int = 32
end

"""
    FeatureVector

Simple feature container for per-vertex decisions.
"""
struct FeatureVector
    nl_len::Int
    ref_len::Int
    overlap::Int
    adds_len::Int
end

"""
    estimate_raw_bits(neighbors::Vector{T}, cfg::EnvConfig) where {T<:Unsigned}

Estimate bits for the raw path: small-count(len) + delta(Fibonacci) list.
"""
function estimate_raw_bits(neighbors::Vector{T}, cfg::EnvConfig) where {T<:Unsigned}
    io = IOBuffer(); w = BitWriter(io)
    Compression.write_small_count(w, T(length(neighbors)), cfg.count_encoding)
    if !isempty(neighbors)
        Compression.write_delta(w, neighbors, cfg.integer_encoding)
    end
    flush_bitwriter(w; flush_last_bits=true)
    return length(take!(io)) * 8
end

"""
    estimate_ref_bits(positions::Vector{Int}, additions::Vector{T}, cfg::EnvConfig) where {T<:Unsigned}

Estimate bits for the reference path using the current write-path formats:
- positions: small-count + delta(Fibonacci) of positions
- additions: small-count + delta(Fibonacci) of additions
"""
function estimate_ref_bits(positions::Vector{Int}, additions::Vector{T}, cfg::EnvConfig) where {T<:Unsigned}
    io = IOBuffer(); w = BitWriter(io)
    # positions
    Compression.write_small_count(w, T(length(positions)), cfg.count_encoding)
    if !isempty(positions)
        Compression.write_delta(w, T.(positions), cfg.integer_encoding)
    end
    # additions
    Compression.write_small_count(w, T(length(additions)), cfg.count_encoding)
    if !isempty(additions)
        Compression.write_delta(w, additions, cfg.integer_encoding)
    end
    flush_bitwriter(w; flush_last_bits=true)
    return length(take!(io)) * 8
end

"""
    greedy_ref_select(nl::Vector{Int}, ref_lists::Vector{Vector{Int}}, cfg::EnvConfig)

Given current neighbors (nl, sorted local indices) and a list of previous
neighbor lists (ref_lists), pick the ref index that minimizes estimated bits.
Returns (use_ref::Bool, best_ref_idx::Int, positions::Vector{Int}, additions::Vector{Int}).
"""
function greedy_ref_select(nl::Vector{Int}, ref_lists::Vector{Vector{Int}}, cfg::EnvConfig)
    isempty(nl) && return (false, 0, Int[], Int[])
    raw_bits = estimate_raw_bits(Int.(nl), cfg)
    best_bits = raw_bits
    best_idx = 0
    best_pos = Int[]
    best_add = Int[]
    for (rix, ref) in pairs(ref_lists)
        i = 1; j = 1
        positions = Int[]; adds = Int[]
        while i <= length(nl) && j <= length(ref)
            if nl[i] == ref[j]
                push!(positions, j); i += 1; j += 1
            elseif nl[i] < ref[j]
                push!(adds, nl[i]); i += 1
            else
                j += 1
            end
        end
        while i <= length(nl)
            push!(adds, nl[i]); i += 1
        end
        ref_bits = estimate_ref_bits(positions, Int.(adds), cfg)
        if ref_bits < best_bits
            best_bits = ref_bits
            best_idx = rix
            best_pos = positions
            best_add = adds
        end
    end
    return (best_idx > 0, best_idx, best_pos, best_add)
end

"""
    choose_positions_mode(ref_len::Int, positions::Vector{Int}, cfg::EnvConfig)

Return :delta or :bitmap based on total-bit estimates for positions only.
"""
function choose_positions_mode(ref_len::Int, positions::Vector{Int}, cfg::EnvConfig)
    # bitmap
    bm_bits = 0
    if ref_len > 0
        copied = fill(false, ref_len)
        for p in positions
            if 1 <= p <= ref_len; copied[p] = true; end
        end
        io = IOBuffer(); w = BitWriter(io)
        Compression.write_bitmap_rle_ones_deltas(w, copied, cfg.integer_encoding)
        flush_bitwriter(w; flush_last_bits=true)
        bm_bits = length(take!(io)) * 8
    end
    # delta
    io2 = IOBuffer(); w2 = BitWriter(io2)
    Compression.write_small_count(w2, UInt32(length(positions)), cfg.count_encoding)
    if !isempty(positions)
        Compression.write_delta(w2, UInt32.(positions), cfg.integer_encoding)
    end
    flush_bitwriter(w2; flush_last_bits=true)
    dl_bits = length(take!(io2)) * 8
    return (bm_bits <= dl_bits) ? :bitmap : :delta
end

end # module RL
