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

module CompressionUtils

using ..Constants: MIN_INTERVAL_LENGTH

export compress_intervals, delta_encode_vector, estimate_encoded_value_cost, estimate_block_encoding_cost,
       reconstruct_from_delta, analyze_delta_patterns_hybrid

function compress_intervals(neighbors::Vector{T}, min_interval_length::Int=4) where {T<:Unsigned}
    isempty(neighbors) && return (Tuple{T,T}[], T[])
    
    intervals = Tuple{T,T}[]  # (start, length)
    residuals = T[]
    i = 1
    
    while i <= length(neighbors)
        # Check for consecutive sequence starting at i
        consecutive_len = 1
        while i + consecutive_len <= length(neighbors) && 
              neighbors[i + consecutive_len] == neighbors[i] + T(consecutive_len)
            consecutive_len += 1
        end
        
        if consecutive_len >= min_interval_length
            # Create interval: (start, length)
            push!(intervals, (neighbors[i], T(consecutive_len)))
            i += consecutive_len
        else
            # Add to residuals
            for j in 0:(consecutive_len-1)
                push!(residuals, neighbors[i + j])
            end
            i += consecutive_len
        end
    end
    
    return (intervals, residuals)
end

function delta_encode_vector(lst::Vector{T})::Vector{T} where {T<:Unsigned}
    # if the list is empty, return an empty list
    isempty(lst) && return T[]
    # initialize the differences with the first element
    diffs = [T(lst[firstindex(lst)])]
    # for each element in the list, compute the difference with the previous element
    for i in eachindex(lst)[2:end]
        push!(diffs, T(lst[i] - lst[i-1]))
    end
    return diffs
end

"""
Precomputed Fibonacci numbers F(2)..F(86) for exact Fibonacci coding bit-length.
Fibonacci coding of positive integer n uses exactly k bits, where F(k) ≤ n < F(k+1).
"""
const _FIB_TABLE = let
    fibs = UInt64[1, 2]  # F(2)=1, F(3)=2
    while fibs[end] < typemax(UInt64) ÷ 2
        push!(fibs, fibs[end] + fibs[end-1])
    end
    fibs
end

"""
Exact bit length of Fibonacci coding for value n ≥ 1.
"""
@inline function _fibonacci_bit_length(n::Integer)::Int
    n <= 0 && return 1
    # Binary search for k such that _FIB_TABLE[k] ≤ n < _FIB_TABLE[k+1]
    # _FIB_TABLE[i] = F(i+1), so bit length = i + 1
    lo, hi = 1, length(_FIB_TABLE)
    while lo < hi
        mid = (lo + hi + 1) >> 1
        if _FIB_TABLE[mid] <= UInt64(n)
            lo = mid
        else
            hi = mid - 1
        end
    end
    return lo + 1  # _FIB_TABLE index lo corresponds to F(lo+1), bit length = lo + 1
end

function estimate_encoded_value_cost(value::T, encoding::Symbol) where {T<:Unsigned}
    if value == 0
        return 1  # Special case
    end

    if encoding == :fibonacci
        return _fibonacci_bit_length(value)
    elseif encoding == :elias_gamma
        # Elias gamma: 2⌊log2(n)⌋ + 1 bits  
        return 2 * floor(Int, log(2, max(1, value))) + 1
    elseif encoding == :elias_delta
        # Elias delta: roughly log2(n) + 2log2(log2(n)) bits
        log_val = max(1, log(2, max(1, value)))
        return ceil(Int, log_val + 2 * log(2, max(1, log_val)))
    elseif encoding == :zeta
        # Zeta coding with k=3: h in unary (h+1 bits) + remainder in truncated binary (~k bits)
        # For k=3: values 1-7 → h=0, cost ≈ 1 + 3 = 4 bits; values 8-63 → h=1, cost ≈ 2 + 6 = 8 bits
        k = 3 # ZETA_BASE
        log2_v = max(0, floor(Int, log(2, max(1, value))))
        h = div(log2_v, k)
        return (h + 1) + k * (h + 1)  # unary(h) + truncated binary bits
    else
        # Default fallback
        return ceil(Int, log(2, max(1, value))) + 2
    end
end

function estimate_block_encoding_cost(copy_bitmap::Vector{Bool}, varint::Symbol)::Int
    if isempty(copy_bitmap)
        return estimate_encoded_value_cost(UInt32(1), varint)
    end

    # Extract blocks to count them
    blocks = UInt32[]
    i = 1
    n = length(copy_bitmap)
    expecting_copy = true

    while i <= n
        block_len = 0

        if expecting_copy
            while i <= n && copy_bitmap[i]
                block_len += 1
                i += 1
            end
        else
            while i <= n && !copy_bitmap[i]
                block_len += 1
                i += 1
            end
        end

        if block_len > 0
            push!(blocks, UInt32(block_len))
            expecting_copy = !expecting_copy
        else
            if expecting_copy && i <= n && !copy_bitmap[i]
                push!(blocks, UInt32(0))
                expecting_copy = false
            else
                break
            end
        end
    end

    # Cost = block count + sum of block lengths
    cost = estimate_encoded_value_cost(UInt32(length(blocks)) + UInt32(1), varint)

    for block_len in blocks
        cost += estimate_encoded_value_cost(block_len + UInt32(1), varint)
    end

    return cost
end

"""
    reconstruct_from_delta(delta_values::Vector{T}) where {T<:Unsigned}

Reconstruct original values from delta encoded list.
"""
function reconstruct_from_delta(delta_values::Vector{T}) where {T<:Unsigned}
    isempty(delta_values) && return T[]
    result = Vector{T}(undef, length(delta_values))
    result[1] = delta_values[1]
    for i in 2:length(delta_values)
        result[i] = result[i-1] + delta_values[i]
    end
    return result
end

function find_consecutive_length(neighbors::Vector{T}, start::Int) where {T<:Unsigned}
    if start > length(neighbors)
        return 0
    end
    
    len = 1
    while start + len <= length(neighbors) && 
          neighbors[start + len] == neighbors[start] + T(len)
        len += 1
    end
    
    return len
end

function count_consecutive(values::Vector{T}, start::Int) where {T<:Unsigned}
    if start > length(values)
        return 0
    end
    
    count = 1
    val = values[start]
    while start + count <= length(values) && values[start + count] == val
        count += 1
    end
    
    return count
end

function find_run_length_patterns(delta_values::Vector{T}) where {T<:Unsigned}
    sections = []
    i = 1
    
    while i <= length(delta_values)
        # Count consecutive occurrences
        current_val = delta_values[i]
        run_len = 1
        while i + run_len <= length(delta_values) && 
              delta_values[i + run_len] == current_val
            run_len += 1
        end
        
        if run_len >= 3  # Use run-length for 3+ consecutive values
            push!(sections, (type=:run_length, data=[current_val, T(run_len)]))
            i += run_len
        else
            # Create delta section - collect individual values until next run
            delta_start = i
            delta_end = i
            
            # Extend until we find another run-length opportunity
            while delta_end < length(delta_values)
                next_run = count_consecutive(delta_values, delta_end + 1)
                if next_run >= 3
                    break
                end
                delta_end += 1
            end
            
            push!(sections, (type=:delta, data=delta_values[delta_start:delta_end]))
            i = delta_end + 1
        end
    end
    
    return sections
end

function analyze_delta_patterns_hybrid(delta_values::Vector{T}, original_neighbors::Vector{T}, min_interval_length::Int=MIN_INTERVAL_LENGTH) where {T<:Unsigned}
    if isempty(delta_values)
        return []
    end
    
    sections = []
    delta_i = 1  # index into delta_values
    orig_i = 1   # index into original_neighbors
    
    while delta_i <= length(delta_values) && orig_i <= length(original_neighbors)
        # Check for consecutive interval starting at current position in original neighbors
        interval_len = find_consecutive_length(original_neighbors, orig_i)
        
        if interval_len >= min_interval_length
            push!(sections, (type=:interval, data=[original_neighbors[orig_i], T(interval_len)]))
            delta_i += interval_len
            orig_i += interval_len
            continue
        end
        
        # No interval found, create delta section for this region
        delta_start = delta_i
        delta_end = delta_i
        orig_end = orig_i
        
        # Extend delta section until we find a good interval
        while orig_end < length(original_neighbors)
            next_interval_len = find_consecutive_length(original_neighbors, orig_end + 1)
            if next_interval_len >= min_interval_length
                break
            end
            delta_end += 1
            orig_end += 1
        end
        
        # Create delta section from the delta values
        if delta_end >= delta_start && delta_end <= length(delta_values)
            section_deltas = delta_values[delta_start:delta_end]
            
            # Check if delta values have run-length patterns worth encoding
            run_length_sections = find_run_length_patterns(section_deltas)
            append!(sections, run_length_sections)
        end
        
        delta_i = delta_end + 1
        orig_i = orig_end + 1
    end
    
    return sections
end

end
