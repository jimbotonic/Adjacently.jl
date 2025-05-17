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
export swap, 
       partition_permutation!, 
       quicksort_iterative_permutation!,
       partition!,
       quicksort_iterative!,
       bottom_up_sort,
       bottom_up_merge,
       binary_search,
       get_sorted_array,
       encode_tree!,
       get_huffman_codes!,
       decode_tree!,
       decode_values,
       huffman_encoding,
       get_next_smallest,
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
# custom implementation of Huffman encoding 
################################################################################

"""
    _sorted_pairs(frequencies)

Helper – turn the Dict into a sorted Vector of (sym,weight) pairs
"""
function _sorted_pairs(freqs::Dict{T, T}) where {T<:Unsigned}
    pairs = collect(freqs)                         # Vector{Pair{T,T}}
    sort!(pairs, by = x -> (x.second, x.first))     
    return pairs
end

"""
    huffman_encoding(data::Vector{T}) where {T<:Unsigned}

Build a Huffman tree from an array of symbols (any Unsigned type)
"""
function huffman_encoding(data::Vector{T}) where {T<:Unsigned}
    isempty(data) && throw(ArgumentError("Input array must not be empty"))

    # Build frequency dictionary
    freqs = Dict{T, T}()
    for s in data
        freqs[s] = get(freqs, s, zero(T)) + one(T)     
    end

    # Use the dictionary-based implementation
    return huffman_encoding(freqs)
end

"""
    huffman_encoding(frequencies::Dict{T,T}; use_heap::Bool = true) where {T<:Unsigned}

Build and return a Huffman tree from a dictionary that maps each symbol ID
to its frequency.

* **`use_heap = true` (default)** – use an incremental binary min-heap.
  Good for medium-to-huge tables because each push / pop is *O(log n)* and no
  full vector re-sorts are needed.

* **`use_heap = false`** – fall back to the fully sorted-vector algorithm you
  already had.  Usually fastest for very small dictionaries because the heap
  setup costs more than one quick sort.

The function handles all corner cases (0, 1, 2 symbols) just like the other
overloads.
"""
function huffman_encoding(freqs::Dict{T,T}; use_heap::Bool = true) where {T<:Unsigned}
    isempty(freqs) &&
        throw(ArgumentError("Frequency table must contain at least one symbol"))

    # Handle trivial cases first
    if length(freqs) == 1
        (sym, _) = first(freqs)
        return Node{T}(sym, EmptyNode, EmptyNode)
    elseif length(freqs) == 2
        (sym1, w1), (sym2, w2) = collect(freqs)
        left, right = w1 < w2 || (w1 == w2 && sym1 < sym2) ?
                      (Node{T}(sym1, EmptyNode, EmptyNode), Node{T}(sym2, EmptyNode, EmptyNode))  :
                      (Node{T}(sym2, EmptyNode, EmptyNode), Node{T}(sym1, EmptyNode, EmptyNode))
        return Node{T}(zero(T), left, right)
    end

    # fast path: using heap
    if use_heap
        # Each heap element is (weight, seq, Node).
        # • weight : the Huffman weight (frequency sum)
        # • seq    : a strictly increasing counter → guarantees uniqueness,
        #            so Node is never looked at when two weights are equal.
        #            That keeps the default tuple ordering well-defined.
        # • Node   : the actual tree node
        #
        HeapElem{T} = Tuple{T, Int, Node{T}}
        heap        = BinaryMinHeap{HeapElem{T}}()      # default ForwardOrdering

        seq = 0
        for (sym, w) in freqs
            seq += 1
            push!(heap, (w, seq, Node{T}(sym, EmptyNode, EmptyNode)))
        end

        while length(heap) > 1
            (w1, _, n1) = pop!(heap)
            (w2, _, n2) = pop!(heap)

            # deterministic left / right choice
            left, right = w1 < w2 || (w1 == w2 && n1.key < n2.key) ? (n1, n2) : (n2, n1)

            seq += 1
            push!(heap, (w1 + w2, seq, Node{T}(zero(T), left, right)))
        end

        return pop!(heap)[3]            # the single remaining Node
    else
        # slow-but-simple path: fully sort once, then vector queue
        pairs = _sorted_pairs(freqs)                   # sorted vector of (symbol, weight) pairs
        pq = [(p.second, Node{T}(p.first, EmptyNode, EmptyNode)) for p in pairs]   # priority queue

        while length(pq) > 1
            (w1, n1) = popfirst!(pq)
            (w2, n2) = popfirst!(pq)

            left, right = w1 < w2 || (w1 == w2 && n1.key < n2.key) ? (n1, n2) : (n2, n1)
            new = (w1 + w2, Node{T}(zero(T), left, right))

            # Insert in sorted order
            pos = searchsortedfirst(pq, new; by = x -> x[1])
            insert!(pq, pos, new)
        end
        return pq[1][2]
    end
end

"""
    encode_tree!(root::AbstractNode, S::BitArray{1}, D::Vector{T}) where {T}

Encode a binary tree into a bitarray and a vector of leaf node values

@param root: root of the tree
@param S: bitarray to store the encoded tree
@param D: vector to store the leaf node values
"""
function encode_tree!(root::AbstractNode, S::BitArray{1}, D::Vector{T}) where {T}
    root === EmptyNode && return
    if root.left === EmptyNode && root.right === EmptyNode      # leaf
        push!(S, true);   push!(D, root.key)
    else                                                        # internal
        push!(S, false)
        encode_tree!(root.left,  S, D)
        encode_tree!(root.right, S, D)
    end
end

function decode_tree!(S::BitArray{1}, D::Vector{T}) where {T}
    isempty(S) && return EmptyNode
    b = popfirst!(S)
    if b == 1                           # leaf
        return Node{T}(popfirst!(D), EmptyNode, EmptyNode)
    else                                # internal
        left  = decode_tree!(S, D)
        right = decode_tree!(S, D)
        return Node{T}(zero(T), left, right)
    end
end

"""
    get_huffman_codes!(root::AbstractNode, C::Dict{BitArray{1},T}, B::BitArray{1}) where {T}

Produce the canonical "code → symbol" dictionary

@param root: root of the tree
@param C: dictionary (bitarray -> value::T)
@param B: current bitarray
"""
function get_huffman_codes!(root::AbstractNode, C::Dict{BitArray{1},T}, B::BitArray{1}) where {T}
    root === EmptyNode && return
    if root.left === EmptyNode && root.right === EmptyNode      # leaf
        C[copy(B)] = root.key
        return
    end
    push!(B, false);  get_huffman_codes!(root.left,  C, B);  pop!(B)
    push!(B, true);   get_huffman_codes!(root.right, C, B);  pop!(B)
end

"""
    decode_values(tree::Node{T}, bits::BitArray{1}) where {T}

Decode a bit-stream of data with the tree

@param tree: root of the tree
@param bits: bitarray to decode
"""
function decode_values(tree::Node{T}, bits::BitArray{1}) where {T}
    out  = T[];        node = tree
    for bit in bits
        node = (bit == 0) ? node.left : node.right
        if node.left === EmptyNode && node.right === EmptyNode
            push!(out, node.key)
            node = tree                       # restart for next symbol
        end
    end
    # the last symbol is still in 'node' if the stream ended on a leaf
    if node.left === EmptyNode && node.right === EmptyNode && node ≠ tree
        push!(out, node.key)
    end
    return out
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

