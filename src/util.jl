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

module Util

using Base.Order
using DataStructures
using Logging

using ..NodeTypes: AbstractNode, Node, EmptyNode
using ..CustomTypes: UInt24, UInt40

# Export the functions we want to make available
export get_next_smallest,
       infer_uint_custom_type,
       infer_uint_std_type,
       to_bytes

################################################################################    
# custom implementation of QuickSort
################################################################################

"""
    swap(A::Vector{T}, i::T, j::T) where {T<:Unsigned}

swap two elements in an array

@param A: array
@param i: index of first element
@param j: index of second element
"""
function swap(A::Vector{T}, i::T, j::T) where {T<:Unsigned}
    t = A[i]
    A[i] = A[j]
    A[j] = t
end

"""
    partition_permutation!(A::Vector{T}, R::Vector{T}, l::T, h::T) where {T<:Unsigned}

Partition function for quicksort that maintains a permutation array.

@param A: array to partition
@param R: permutation array to maintain
@param l: lower index of partition range
@param h: higher index of partition range
@return index of the pivot element after partitioning
"""
function partition_permutation!(A::Vector{T}, R::Vector{T}, l::T, h::T) where {T<:Unsigned}
    pivot = A[h]
    i = convert(T, l)
    
    for j in l:(h-1)
        if A[j] < pivot
            if i != j
                # Convert indices to the appropriate type
                i_t = convert(T, i)
                j_t = convert(T, j)
                swap(A, i_t, j_t)
                swap(R, i_t, j_t)
            end
            i += 1
        end
    end
    
    if i != h
        # Convert indices for final swap
        i_t = convert(T, i)
        h_t = convert(T, h)
        swap(A, i_t, h_t)
        swap(R, i_t, h_t)
    end
    
    return convert(T, i)
end

"""
    quicksort_iterative_permutation!(A::Vector{T}) where T <: Unsigned

iterative quicksort with permutation array

@return permutation array and sort A in ascending order
"""
function quicksort_iterative_permutation!(A::Vector{T}) where T <: Unsigned
    n = convert(T, length(A))
    P = Vector{T}(1:n)  # Initialize permutation array with same type T
    
    # Create a stack for storing the partition range
    stack = Tuple{T, T}[]
    
    # Initialize the range with first and last indices
    push!(stack, (one(T), n))
    
    # Process the stack
    while !isempty(stack)
        low, high = pop!(stack)
        
        # If there are elements to be sorted
        if low < high
            # Get pivot element
            pivot = partition_permutation!(A, P, low, high)
            
            # Push subarrays to stack
            if pivot > low
                push!(stack, (low, pivot - one(T)))
            end
            if pivot < high
                push!(stack, (pivot + one(T), high))
            end
        end
    end
    return P
end

"""
    partition!(A::Vector{T}, low::T, high::T) where T <: Unsigned

partition function

@param A: array
@param low: lower index
@param high: higher index
"""
function partition!(A::Vector{T}, low::T, high::T) where T <: Unsigned
    pivot = A[high]
    i = low - one(T)
    
    for j in low:(high - one(T))
        if A[j] <= pivot
            i += one(T)
            A[i], A[j] = A[j], A[i]
        end
    end
    
    A[i+one(T)], A[high] = A[high], A[i + one(T)]
    return i + one(T)
end

"""
    quicksort_iterative!(A::Vector{T}) where T <: Unsigned

iterative quicksort

@return permutation array and sort A in ascending order
"""
function quicksort_iterative!(A::Vector{T}) where T <: Unsigned
    # Create a stack for storing the partition range
    stack = Tuple{T, T}[]
    
    # Initialize the range with first and last indices
    push!(stack, (one(T), convert(T, length(A))))
    
    # Process the stack
    while !isempty(stack)
        low, high = pop!(stack)
        
        # If there are elements to be sorted
        if low < high
            # Get pivot element
            pivot = partition!(A, low, high)
            
            # Push subarrays to stack
            if pivot > low
                push!(stack, (low, pivot - one(T)))
            end
            if pivot < high
                push!(stack, (pivot + one(T), high))
            end
        end
    end
    return A
end

################################################################################
# custom implementation of MergeSort
################################################################################

"""
    bottom_up_sort(A::Vector{T}) where {T<:Unsigned}

merge sort (ascending order)

@returns permutation array R

NB: to get the sorted array sA: sA[i] = A[R[i]] 
     R[i] is thus the index in the original array of element at new index i in the sorted array
"""
function bottom_up_sort(A::Vector{T}) where {T<:Unsigned}
	n = convert(T, length(A))
	B = zeros(T, n)
	# permutation array
	R = collect(convert(T, 1):n)
	S = zeros(T, n)
	width = 1
	while width < n
		i = 0
		while i < n
			bottom_up_merge(A, convert(T,i), convert(T,min(i+width,n)), convert(T,min(i+2*width,n)), B, R, S)
			i += 2*width
		end
		A = copy(B)
		R = copy(S)
		width = 2*width
	end
	return R
end

"""
    bottom_up_merge(A::Vector{T}, iLeft::T, iRight::T, iEnd::T, B::Vector{T}, R::Vector{T}, S::Vector{T}) where {T<:Unsigned}

merge two sorted arrays

@param A: array
@param iLeft: left index
@param iRight: right index
@param iEnd: end index
@param B: array
@param R: permutation array
@param S: sorted array
"""
function bottom_up_merge(A::Vector{T}, iLeft::T, iRight::T, iEnd::T, B::Vector{T}, R::Vector{T}, S::Vector{T}) where {T<:Unsigned}
	i0 = iLeft
	i1 = iRight
	for j in (iLeft+1):iEnd
		if i0 < iRight && (i1 >= iEnd || A[i0+1] <= A[i1+1])
			B[j] = A[i0+1]
			S[j] = R[i0+1]
			i0 += 1
		else
			B[j] = A[i1+1]
			S[j] = R[i1+1]
			i1 += 1
		end
	end
end

"""
    binary_search(A::Vector{T}, x::T) where {T<:Unsigned}

custom binary search

search x position in array A

NB: array A is assumed to be sorted in ascending order
"""
function binary_search(A::Vector{T}, x::T) where {T<:Unsigned}
	low = 1 
	high = length(A)
	while true
		if high == low
			if A[low] == x
				return low
			else
				return -1
			end
		else
			p = low + floor(Int, (high-low)/2)
			if x == A[p]
				return p
			elseif x > A[p]
				if p < high
					low = p+1
				else
					return -1
				end
			else
				if p > low
					high = p-1
				else
					return -1
				end
			end
		end
	end
end

"""
    get_sorted_array(A::Vector{T}, R::Vector{T}, asc::Bool=true) where {T}

get sorted array

@param A: initial array
@param R: permutation array
@param asc: ascending order
"""
function get_sorted_array(A::Vector{T}, R::Vector{T}, asc::Bool=true) where {T}
	S = zeros(Int,length(A))
	n = length(A)
	# ascending order
	if asc
		for i in 1:n
			S[i] = A[R[i]]
		end
	# decreasing order
	else
		for i in 1:n
			S[n+1-i] = A[R[i]]
		end
	end
	return S
end

################################################################################
# custom unsigned integer types
################################################################################

"""
	infer_uint_custom_type(n::UInt8)

Infer the unsigned integer type needed to store the specified size

NB: standard types are UInt8, UInt16, UInt32, UInt64
NB: custom types are UInt24, UInt40
"""
function infer_uint_custom_type(n_bits::UInt8)
    if n_bits <= 8
        return UInt8
    elseif n_bits <= 16
        return UInt16
    elseif n_bits <= 24
		return UInt24
    elseif n_bits <= 32
        return UInt32
    elseif n_bits <= 40
		return UInt40
    else
        throw(ArgumentError("List size is too large. Maximum supported bits are 40."))
    end
end

"""
	infer_uint_std_type(n::UInt8)

Infer the unsigned integer type needed to store the specified size

NB: standard types are UInt8, UInt16, UInt32, UInt64
"""
function infer_uint_std_type(n_bits::UInt8)
    if n_bits <= 8
        return UInt8
    elseif n_bits <= 16
        return UInt16
    elseif n_bits <= 32
        return UInt32
    elseif n_bits <= 64
		return UInt64
    else
        throw(ArgumentError("List size is too large. Maximum supported bits are 64."))
    end
end

"""
    to_bytes(x::T) where T <: Unsigned

Convert unsigned integer to array of bytes in big-endian order.
"""
function to_bytes(x::T) where T <: Unsigned
    n = sizeof(T)
    bytes = Vector{UInt8}(undef, n)
    for i in 1:n
        bytes[i] = (x >> (8 * (n - i))) & 0xFF
    end
    return bytes
end

"""
    to_bytes(x::UInt24)

Convert UInt24 to array of 3 bytes in big-endian order.
"""
function to_bytes(x::UInt24)
    val = reinterpret(UInt32, x) & 0xFFFFFF
    return UInt8[
        (val >> 16) & 0xFF,
        (val >> 8) & 0xFF,
        val & 0xFF
    ]
end

"""
    to_bytes(x::UInt40)

Convert UInt40 to array of 5 bytes in big-endian order.
"""
function to_bytes(x::UInt40)
    val = reinterpret(UInt64, x) & 0xFFFFFFFFFF
    return UInt8[
        (val >> 32) & 0xFF,
        (val >> 24) & 0xFF,
        (val >> 16) & 0xFF,
        (val >> 8) & 0xFF,
        val & 0xFF
    ]
end

end # module Util

