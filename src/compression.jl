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

module Compression

using DataStructures
using ..NodeTypes: Node, EmptyNode, AbstractNode
using ..CustomTypes: UInt24, UInt40

# Export the functions we want to make available
export huffman_encoding,
       encode_tree!,
       decode_tree!,
       get_huffman_codes!,
       decode_values

################################################################################
# Huffman encoding
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

end # module Compression
