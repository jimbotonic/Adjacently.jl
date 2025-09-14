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
using SparseArrays
using ..NodeTypes: Node, EmptyNode, AbstractNode
using ..CustomTypes: UInt24, UInt40

using ..IO: BitWriter, BitReader, write_bit, write_bits, read_bit, read_bits, peek_bit, 
flush, write_value, read_value, flush_bitwriter

using ..Constants: FIB_NUMBERS, BUFFER_SIZE, ZETA_H_BOUNDS, ZETA_POWER_BASES, ZETA_BASE, 
GOLOMB_BASE, REF_ENCODING_TH, REF_V_MIN_DEGREE, FED_BLOCK_SIZE, MIN_INTERVAL_LENGTH, REF_WINDOW_SIZE

using ..Index

# Export the functions we want to make available
export write_unary_coding,
       write_truncated_binary_coding,
       write_encoded_value,
       read_encoded_value,
       get_encoded_value,
       huffman_encoding,
       encode_huffman_tree!,
       decode_huffman_tree!,
       get_huffman_codes!,
       decode_huffman_values,
       delta_encode_vector,
       write_delta,
       read_delta,
       write_elias_gamma,
       read_elias_gamma,
       write_elias_coding,
       read_elias_coding,
       write_elias_fano,
       read_elias_fano,
       write_golomb,
       read_golomb,
       write_fibonacci,
       read_fibonacci,
       write_zeta,
       read_zeta,
       write_fed,
       read_fed,
       write_run_length_delta,
       read_run_length_delta,
       compress_intervals,
       write_intervals_and_residuals,
       read_intervals_and_residuals,
       write_compressed_graph_data,
       read_compressed_graph_data

# Lightweight workspace types for reference building
struct RefBuildWorkspace{T<:Unsigned}
    copy_bitmap::Vector{Bool}
    residuals::Vector{T}
end

RefBuildWorkspace{T}() where {T<:Unsigned} = RefBuildWorkspace{T}(Bool[], T[])

################################################################################
# Basic encoding / decoding
################################################################################

"""
    write_unary_coding(w::BitWriter, v::T) where {T<:Unsigned}

Write a unary coding to the bitwriter.

@param w::BitWriter: the bitwriter
@param v::T: the value to write
"""
function write_unary_coding(w::BitWriter, v::T, invert::Bool = false) where {T<:Unsigned}
    for i in 1:v
        write_bit(w, invert)
    end
    write_bit(w, !invert)
end

"""
    read_unary_coding(w::BitReader, invert::Bool = false, ::Type{T}=UInt8) where {T<:Unsigned}

Read a unary coding from the bitreader.

@param w::BitReader: the bitreader
@param invert::Bool: whether to invert the bits
@param T::Type: the type to return (default: UInt8)
@return::T: the decoded value
"""
function read_unary_coding(w::BitReader, invert::Bool = false, ::Type{T}=UInt8) where {T<:Unsigned}
    v = 0
    while read_bit(w) == invert
        v += 1
    end
    return convert(T, v)
end

"""
    write_truncated_binary_coding(w::BitWriter, v::T, n::Int) where {T<:Unsigned}

Write a truncated binary code to the bitwriter.

Truncated binary coding for integers in range [0, n-1]:
- If n is a power of 2, use standard binary with log2(n) bits
- Otherwise, use k = floor(log2(n)) bits for values < 2^(k+1) - n
  and k+1 bits for remaining values

@param w::BitWriter: the bitwriter
@param v::T: the value to write (must be in range [0, n-1])
@param n::T: the range size (number of possible values)
"""
function write_truncated_binary_coding(w::BitWriter, v::T, n::Int) where {T<:Unsigned}
    n == 0 && throw(ArgumentError("Range size n must be > 0"))
    v >= n && throw(ArgumentError("Value v=$v must be < n=$n"))
    
    # special case: if n is 1, no bits needed
    if n == 1
        return
    end
    
    # calculate k = floor(log2(n))
    k = floor(Int, log2(n))
    
    # calculate threshold u = 2^(k+1) - n using UInt64 to avoid overflow
    u_u64 = (UInt64(1) << (k + 1)) - UInt64(n)
    u = T(u_u64)
    
    if v < u
        # use k bits for values in [0, u-1]
        write_value(w, v, k)
    else
        # use k+1 bits for values in [u, n-1]
        # encode as v + u in k+1 bits
        write_value(w, v + u, k + 1)
    end
end

"""
    read_truncated_binary_coding(w::BitReader, n::Int, ::Type{T}=UInt8) where {T<:Unsigned}

Read a truncated binary code from the bitreader.

@param w::BitReader: the bitreader
@param n::Int: the range size (number of possible values)
@param T::Type: the type to return (default: UInt8)
@return::T: the decoded value
"""
function read_truncated_binary_coding(w::BitReader, n::Int, ::Type{T}=UInt8) where {T<:Unsigned}
    n == 0 && throw(ArgumentError("Range size n must be > 0"))
    
    # Special case: if n is 1, no bits needed
    if n == 1
        return T(0)
    end
    
    # Calculate k = floor(log2(n))
    k = floor(Int, log2(n))
    
    # Calculate threshold u = 2^(k+1) - n using UInt64 to avoid overflow
    u_u64 = (UInt64(1) << (k + 1)) - UInt64(n)
    u = Int(u_u64)
    
    # Read k bits first
    v = read_value(w, k, T)
    
    if v < u
        # Use k bits for values in [0, u-1]
        return v
    else
        # Read one more bit and reconstruct the value
        extra_bit = read_bit(w) ? 1 : 0
        return T(u + (v - u) * 2 + extra_bit)
    end
end

################################################################################
# Dispatching functions
################################################################################

"""
    write_encoded_value(w::BitWriter, value::T, compression::Symbol) where {T<:Unsigned}

Write a value to the bitwriter using the specified compression code.

@param w::BitWriter: the bitwriter
@param value::T: the value to write
@param compression::Symbol: the compression coding to use (:elias_gamma, :elias_delta, :golomb, :fibonacci, :zeta, :fed)
"""
function write_encoded_value(w::BitWriter, value::T, compression::Symbol) where {T<:Unsigned}
    if compression == :elias_gamma
        write_elias_gamma(w, value)
    elseif compression == :elias_delta
        write_elias_delta(w, value)
    elseif compression == :golomb
        write_golomb(w, value, GOLOMB_BASE)
    elseif compression == :fibonacci
        write_fibonacci(w, value)
    elseif compression == :zeta
        write_zeta(w, value, ZETA_BASE)
    elseif compression == :fed
        write_fed(w, value, FED_BLOCK_SIZE)
    else
        throw(ArgumentError("Invalid compression code: $compression"))
    end
end

"""
    read_encoded_value(r::BitReader, compression::Symbol, ::Type{T}=UInt8) where {T<:Unsigned}

Read a value from the bitreader using the specified compression code.

@param r::BitReader: the bitreader
@param compression::Symbol: the compression coding to use (:elias_gamma, :elias_delta, :golomb, :fibonacci, :zeta, :fed)
@param T::Type: the type to return (default: UInt8)
@return::T: the decoded value
"""
function read_encoded_value(r::BitReader, compression::Symbol, ::Type{T}=UInt8) where {T<:Unsigned}
    if compression == :elias_gamma
        return read_elias_gamma(r, T)
    elseif compression == :elias_delta
        return read_elias_delta(r, T)
    elseif compression == :golomb
        return read_golomb(r, GOLOMB_BASE, T)
    elseif compression == :fibonacci
        return read_fibonacci(r, T)
    elseif compression == :zeta
        return read_zeta(r, ZETA_BASE, T)
    elseif compression == :fed
        return read_fed(r, T, FED_BLOCK_SIZE)
    else
        throw(ArgumentError("Invalid compression code: $compression"))
    end
end

"""
    get_encoded_value(value::T, compression::Symbol, ::Type{T}=UInt8) where {T<:Unsigned}

Get what a value would be after encoding and then decoding through a compression scheme.
This is useful for comparing stop values with encoded stream values.

@param value::T: the value to encode/decode
@param compression::Symbol: the compression coding to use (:elias_gamma, :elias_delta, :golomb, :fibonacci, :zeta, :fed)
@param T::Type: the type to return (default: UInt8)
@return::T: the value as it would appear after encoding/decoding
"""
function get_encoded_value(value::T, compression::Symbol, ::Type{T}=UInt8) where {T<:Unsigned}
    # Create temporary buffer to encode and decode
    temp_io = IOBuffer()
    temp_writer = BitWriter(temp_io)
    
    # Encode the value
    write_encoded_value(temp_writer, value, compression)
    flush_bitwriter(temp_writer; flush_last_bits=true)
    
    # Read it back
    seekstart(temp_io)
    temp_reader = BitReader(temp_io)
    return read_encoded_value(temp_reader, compression, T)
end

################################################################################
# Huffman encoding
################################################################################

"""
    _sorted_pairs(frequencies::Dict{T, T})::Vector{Pair{T,T}} where {T<:Unsigned}

Helper – turn the Dict into a sorted Vector of (sym,weight) pairs

@param freqs::Dict{T, T}: The dictionary of frequencies
@return::Vector{Pair{T,T}}: The sorted vector of (sym,weight) pairs
"""
function _sorted_pairs(freqs::Dict{T, T})::Vector{Pair{T,T}} where {T<:Unsigned}
    # collect the frequencies into a vector of pairs
    pairs = collect(freqs)
    # sort the pairs by weight and then by symbol
    sort!(pairs, by = x -> (x.second, x.first))
    # return the sorted pairs
    return pairs
end

"""
    huffman_encoding(data::Vector{T})::Node{T} where {T<:Unsigned}

Build a Huffman tree from an array of symbols (any Unsigned type)

@param data::Vector{T}: The vector of symbols to encode
@return::Node{T}: The root of the Huffman tree
"""
function huffman_encoding(data::Vector{T})::Node{T} where {T<:Unsigned}
    # if the data is empty, throw an error
    isempty(data) && throw(ArgumentError("Input array must not be empty"))
    # build the frequency dictionary
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

* **`use_heap = true` (default)** - use an incremental binary min-heap.
  Good for medium-to-huge tables because each push / pop is *O(log n)* and no
  full vector re-sorts are needed.

* **`use_heap = false`** - fall back to the fully sorted-vector algorithm you
  already had.  Usually fastest for very small dictionaries because the heap
  setup costs more than one quick sort.

The function handles all corner cases (0, 1, 2 symbols) just like the other
overloads.
"""
function huffman_encoding(freqs::Dict{T,T}; use_heap::Bool = true) where {T<:Unsigned}
    # if the frequency dictionary is empty, throw an error
    isempty(freqs) && throw(ArgumentError("Frequency table must contain at least one symbol"))

    # handle trivial cases first
    if length(freqs) == 1
        (sym, _) = first(freqs)
        # return the single node
        return Node{T}(sym, EmptyNode, EmptyNode)
    # handle the case of two symbols
    elseif length(freqs) == 2
        # collect the frequencies into a vector of pairs
        (sym1, w1), (sym2, w2) = collect(freqs)
        # determine the left and right nodes
        left, right = w1 < w2 || (w1 == w2 && sym1 < sym2) ?
                      (Node{T}(sym1, EmptyNode, EmptyNode), Node{T}(sym2, EmptyNode, EmptyNode))  :
                      (Node{T}(sym2, EmptyNode, EmptyNode), Node{T}(sym1, EmptyNode, EmptyNode))
        return Node{T}(zero(T), left, right)
    end

    # fast path: using heap
    if use_heap
        # Each heap element is (weight, seq, Node).
        # - weight : the Huffman weight (frequency sum)
        # - seq    : a strictly increasing counter → guarantees uniqueness,
        #           so Node is never looked at when two weights are equal.
        #           That keeps the default tuple ordering well-defined.
        # - Node   : the actual tree node
        #
        # define the type of the heap elements
        HeapElem{T} = Tuple{T, Int, Node{T}}
        # create the heap
        heap        = BinaryMinHeap{HeapElem{T}}()      # default ForwardOrdering

        # initialize the sequence counter
        seq = 0
        # for each symbol in the frequency dictionary, push the weight and the node into the heap
        for (sym, w) in freqs
            # increment the sequence counter
            seq += 1
            # push the weight and the node into the heap
            push!(heap, (w, seq, Node{T}(sym, EmptyNode, EmptyNode)))
        end

        # while the heap has more than one element, pop the two smallest elements and create a new node
        while length(heap) > 1
            # pop the two smallest elements
            (w1, _, n1) = pop!(heap)
            (w2, _, n2) = pop!(heap)
            # determine the left and right nodes
            # deterministic left / right choice
            left, right = w1 < w2 || (w1 == w2 && n1.key < n2.key) ? (n1, n2) : (n2, n1)
            # increment the sequence counter
            seq += 1
            # push the new node into the heap
            push!(heap, (w1 + w2, seq, Node{T}(zero(T), left, right)))
        end

        # return the single remaining Node
        return pop!(heap)[3]
    else
        # slow-but-simple path: fully sort once, then vector queue
        # sorted vector of (symbol, weight) pairs
        pairs = _sorted_pairs(freqs)
        # create the priority queue
        pq = [(p.second, Node{T}(p.first, EmptyNode, EmptyNode)) for p in pairs]   # priority queue

        while length(pq) > 1
            # pop the two smallest elements
            (w1, n1) = popfirst!(pq)
            (w2, n2) = popfirst!(pq)
            # determine the left and right nodes
            left, right = w1 < w2 || (w1 == w2 && n1.key < n2.key) ? (n1, n2) : (n2, n1)
            # create the new node
            new = (w1 + w2, Node{T}(zero(T), left, right))
            # insert the new node in the priority queue
            pos = searchsortedfirst(pq, new; by = x -> x[1])
            # insert the new node in the priority queue
            insert!(pq, pos, new)
        end
        # return the single remaining Node
        return pq[1][2]
    end
end

"""
    encode_huffman_tree!(root::AbstractNode, S::BitArray{1}, D::Vector{T}) where {T<:Unsigned}

Encode a binary tree into a bitarray and a vector of leaf node values

@param root: root of the tree
@param S: bitarray to store the encoded tree
@param D: vector to store the leaf node values
"""
function encode_huffman_tree!(root::AbstractNode, S::BitArray{1}, D::Vector{T}) where {T<:Unsigned}
    # if the root is empty, return
    root === EmptyNode && return
    # if the root is a leaf, push true and the key to the bitarray and the vector
    if root.left === EmptyNode && root.right === EmptyNode      # leaf
        push!(S, true);   push!(D, root.key)
    # if the root is an internal node, push false and encode the left and right children
    else
        # push false to the bitarray
        push!(S, false)
        # encode the left child
        encode_huffman_tree!(root.left,  S, D)
        # encode the right child
        encode_huffman_tree!(root.right, S, D)
    end
end

"""
    decode_huffman_tree!(S::BitArray{1}, D::Vector{T}) where {T<:Unsigned}

Decode a binary tree from a bitarray and a vector of leaf node values

@param S: bitarray to decode
@param D: vector to store the leaf node values
"""
function decode_huffman_tree!(S::BitArray{1}, D::Vector{T})::AbstractNode where {T<:Unsigned}
    # if the bitarray is empty, return
    isempty(S) && return EmptyNode
    # pop the first bit
    b = popfirst!(S)
    # if the bit is 1, return a leaf node
    if b == 1
        # return the leaf node
        return Node{T}(popfirst!(D), EmptyNode, EmptyNode)
    # if the bit is 0, decode the left and right children
    else
        # decode the left child
        left  = decode_huffman_tree!(S, D)
        # decode the right child
        right = decode_huffman_tree!(S, D)
        # return the new node
        return Node{T}(zero(T), left, right)
    end
end

"""
    get_huffman_codes!(root::AbstractNode, C::Dict{BitArray{1},T}, B::BitArray{1}) where {T<:Unsigned}

Produce the canonical "code → symbol" dictionary

@param root: root of the tree
@param C: dictionary (bitarray -> value::T)
@param B: current bitarray
"""
function get_huffman_codes!(root::AbstractNode, C::Dict{BitArray{1},T}, B::BitArray{1}) where {T<:Unsigned}
    # if the root is empty, return
    root === EmptyNode && return
    # if the root is a leaf, copy the bitarray and add the key to the dictionary
    if root.left === EmptyNode && root.right === EmptyNode
        # copy the bitarray and add the key to the dictionary
        C[copy(B)] = root.key
        return
    end
    # push false to the bitarray
    push!(B, false);
    # get the codes for the left child
    get_huffman_codes!(root.left,  C, B);
    # pop the last bit
    pop!(B)
    # push true to the bitarray
    push!(B, true);
    # get the codes for the right child
    get_huffman_codes!(root.right, C, B);
    # pop the last bit
    pop!(B)
end

"""
    decode_huffman_values(tree::Node{T}, bits::BitArray{1})::Vector{T} where {T<:Unsigned}

Decode a bit-stream of data with the tree

@param tree: root of the tree
@param bits: bitarray to decode
"""
function decode_huffman_values(tree::Node{T}, bits::BitArray{1})::Vector{T} where {T<:Unsigned}
    # initialize the output vector
    out  = T[]
    # initialize the node
    node = tree
    # for each bit in the bitarray
    for bit in bits
        # if the bit is 0, go to the left child
        node = (bit == 0) ? node.left : node.right
        # if the node is a leaf, push the key to the output vector and restart for next symbol
        if node.left === EmptyNode && node.right === EmptyNode
            # push the key to the output vector
            push!(out, node.key)
            # restart for next symbol
            node = tree
        end
    end
    # if the last symbol is still in 'node' and the stream ended on a leaf
    if node.left === EmptyNode && node.right === EmptyNode && node !== tree
        # push the key to the output vector
        push!(out, node.key)
    end
    # return the output vector
    return out
end

################################################################################
# Elias gamma code
################################################################################

"""
    delta_encode_vector(lst::Vector{T}, shifted::Bool = false)::Vector{T} where {T<:Unsigned}

Delta encode a list of integers

@param lst::Vector{T}: The list of integers to delta encode
@param shifted::Bool: Whether to shift the list by 1 to avoid 0
@return::Vector{T}: The delta encoded list
"""
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
    write_elias_gamma(w::BitWriter, v::T) where {T<:Unsigned}

Write an Elias gamma code to the bitwriter.

- Gamma code for `v` (with v >= 1) is made of:
    - `k` zeros, where `k = floor(log2(v))`
    - followed by the binary representation of `v` on `k + 1` bits

For example, v = 5 (0b101) -> 2 zeros + 101 -> bits = 0 0 1 0 1
"""
function write_elias_gamma(w::BitWriter, v::T) where {T<:Unsigned}
    v == 0 && throw(ArgumentError("Elias gamma is undefined for 0"))

    # Determine the number of bits needed to represent v
    bits_len = sizeof(T) * 8 - leading_zeros(v)

    # Write prefix of (bits_len - 1) zeros
    for _ in 1:(bits_len - 1)
        write_bit(w, false)
    end

    # Write the actual binary representation of v
    write_value(w, v, bits_len)
end

"""
    write_elias_delta(w::BitWriter, v::T) where {T<:Unsigned}

Write an Elias delta code to the bitwriter.

Elias delta encoding of n ≥ 1:
- Let l = floor(log2(n)) + 1 = number of bits in binary representation of n
- Write l using Elias gamma coding
- Then write the (l - 1) lower bits of n (i.e., remove the leading 1)
"""
function write_elias_delta(w::BitWriter, v::T) where {T<:Unsigned}
    v == 0 && throw(ArgumentError("Elias delta is undefined for 0"))

    # l = number of bits needed to represent n
    bits_len = sizeof(T) * 8 - leading_zeros(v)

    # Step 1: write bits_len using Elias gamma
    write_elias_gamma(w, T(bits_len))

    # Step 2: write the binary of v without the leading 1
    mask = (T(1) << (bits_len - 1)) - 1
    tail_bits = convert(T, v & mask)
    write_value(w, tail_bits, bits_len - 1)
end

"""
    write_elias_coding(w::BitWriter, v::T, coding::Symbol) where {T<:Unsigned}

Write an Elias coding to the bitwriter.

@param w::BitWriter: The bit writer to write to
@param v::T: The value to write
@param coding::Symbol: The coding to use (:elias_gamma or :elias_delta)
"""
function write_elias_coding(w::BitWriter, v::T, coding::Symbol) where {T<:Unsigned}
    if coding == :elias_gamma
        write_elias_gamma(w, v)
    elseif coding == :elias_delta
        write_elias_delta(w, v)
    else
        throw(ArgumentError("Unsupported coding scheme: $coding. Supported schemes are :elias_gamma and :elias_delta"))
    end
end

"""
    read_elias_gamma(r::BitReader, ::Type{T}=UInt8) where {T<:Unsigned}

Read an Elias gamma encoded value from the bit reader.

Gamma code format:
- `k` leading zeros
- followed by `k + 1` bits representing the value

Returns the decoded value as type `T`.
"""
function read_elias_gamma(r::BitReader, ::Type{T}=UInt8) where {T<:Unsigned}
    # Step 1: Count leading zeros
    zeros = 0
    while !read_bit(r)
        zeros += 1
    end

    # Step 2: Read the remaining `zeros` bits and combine with implicit leading 1
    tail = read_value(r, zeros, T)
    return (T(1) << zeros) | tail
end

"""
    read_elias_delta(r::BitReader, ::Type{T}=UInt8) where {T<:Unsigned}

Read an Elias delta encoded value from the bit reader.

Delta code format:
1. Read l = read_elias_gamma(r)
2. Read (l - 1) bits -> the value without the leading 1
3. Return value = 1 << (l - 1) | tail
"""
function read_elias_delta(r::BitReader, ::Type{T}=UInt8) where {T<:Unsigned}
    # Step 1: Read the length of the value
    len = Int(read_elias_gamma(r, T))

    len < 1 && throw(ArgumentError("Invalid Elias delta encoding: length must be >= 1"))

    # Step 2: Read the (len - 1) lower bits
    tail = read_value(r, len-1, T)

    # Step 3: Combine with leading 1
    return (T(1) << (len-1)) | tail
end

"""
    read_elias_coding(r::BitReader, coding::Symbol, ::Type{T}=UInt8) where {T<:Unsigned}

Read an Elias coding from the bit reader

@param r::BitReader: The bit reader to read from
@param coding::Symbol: The coding to use (:elias_gamma or :elias_delta)
@param T::Type: The type to return (default: UInt8)
@return::T: The value read
"""
function read_elias_coding(r::BitReader, coding::Symbol, ::Type{T}=UInt8) where {T<:Unsigned}
    if coding == :elias_gamma
        return read_elias_gamma(r, T)
    elseif coding == :elias_delta
        return read_elias_delta(r, T)
    else
        throw(ArgumentError("Unsupported coding scheme: $coding. Supported schemes are :elias_gamma and :elias_delta"))
    end
end

################################################################################
# Elias-Fano encoding
################################################################################

"""
    write_elias_fano(w::BitWriter, values::Vector{T}) where {T<:Unsigned}

Write a sequence of monotonic values using Elias-Fano encoding.

Elias-Fano encoding is optimal for monotonic sequences and consists of:
1. High bits: Unary encoded gaps between values shifted by log₂(universe_size/n)
2. Low bits: Lower log₂(universe_size/n) bits of each value

@param w::BitWriter: The bit writer to write to
@param values::Vector{T}: Monotonic sequence of values to encode (must be sorted)
"""
function write_elias_fano(w::BitWriter, values::Vector{T}) where {T<:Unsigned}
    n = length(values)
    
    # Write header: number of values (add 1 since Elias gamma doesn't support 0)
    write_encoded_value(w, T(n + 1), :elias_gamma)
    
    # Handle empty sequence
    if n == 0
        return
    end
    
    # Universe size (highest value + 1)
    universe_size = Int(values[end]) + 1
    
    # Number of low bits per value
    low_bits = n == 1 ? 0 : max(0, floor(Int, log2(universe_size / n)))
    
    # Write low_bits in header (add 1 since Elias gamma doesn't support 0)
    write_encoded_value(w, T(low_bits + 1), :elias_gamma)
    
    # Extract high and low parts
    low_mask = (T(1) << low_bits) - T(1)
    
    # Write low bits first
    if low_bits > 0
        for value in values
            low_part = T(value & low_mask)
            write_value(w, low_part, low_bits)
        end
    end
    
    # Write high bits using unary gaps
    # Standard Elias-Fano: encode gaps between consecutive high parts
    prev_high = T(0) 
    for (i, value) in enumerate(values)
        high_part = value >> low_bits
        
        if i == 1
            # First value: write gap from 0 to high_part
            gap = high_part
        else
            # Subsequent values: gap from previous high_part  
            gap = high_part - prev_high
        end
        
        # Write gap as unary: gap zeros followed by one 1
        for _ in 1:gap
            write_bit(w, false)
        end
        write_bit(w, true)
        
        prev_high = high_part
    end
end

"""
    read_elias_fano(r::BitReader, ::Type{T}=UInt32) where {T<:Unsigned}

Read a sequence encoded with Elias-Fano encoding.

@param r::BitReader: The bit reader to read from
@param T::Type: The type to return (default: UInt32)
@returns Vector{T}: The decoded monotonic sequence
"""
function read_elias_fano(r::BitReader, ::Type{T}=UInt32) where {T<:Unsigned}
    # Read header (subtract 1 since we added 1 during encoding)
    n = Int(read_encoded_value(r, :elias_gamma, T)) - 1
    
    # Handle empty sequence
    if n == 0
        return T[]
    end
    
    # Read low_bits (subtract 1 since we added 1 during encoding)
    low_bits = Int(read_encoded_value(r, :elias_gamma, T)) - 1
    
    # Read low bits
    low_parts = Vector{T}(undef, n)
    low_mask = (T(1) << low_bits) - 1
    
    if low_bits > 0
        for i in 1:n
            low_parts[i] = T(read_value(r, low_bits, T)) & low_mask
        end
    else
        # No low bits to read, all values are 0
        for i in 1:n
            low_parts[i] = T(0)
        end
    end
    
    # Read high bits
    values = Vector{T}(undef, n)
    current_high = T(0)
    
    for i in 1:n
        # Count zeros before the next 1
        gap = 0
        while !read_bit(r)
            gap += 1
        end
        
        # Update current high part
        current_high += T(gap)
        
        # Reconstruct value
        high_part = current_high << low_bits
        values[i] = high_part | low_parts[i]
        
        # Note: don't increment current_high here - it stays at the actual high part
        # for the next iteration to calculate the correct gap
    end
    
    return values
end

################################################################################
# Golomb code
################################################################################

"""
    write_golomb(w::BitWriter, n::T, b::Int) where {T<:Unsigned}

Write a Golomb code for value `n` with base `b`.

Encoding steps:
1. Quotient q = div(n, b), encoded in unary: q zeros followed by one `1`
2. Remainder r = n % b, encoded in binary with ⌈log₂(b)⌉ bits
"""
function write_golomb(w::BitWriter, n::T, b::Int) where {T<:Unsigned}
    b <= 0 && throw(ArgumentError("Golomb base must be >= 1"))

    # Step 1: Unary representation of quotient
    q = div(n, b)
    for _ in 1:q
        write_bit(w, false)
    end
    write_bit(w, true)  # stop bit

    # Step 2: Binary representation of remainder
    r = n % b
    m = ceil(Int, log2(b))
    write_value(w, r, m)
end

"""
    read_golomb(r::BitReader, b::Int, ::Type{T}=UInt8) where {T<:Unsigned}

Read a Golomb-encoded value from the bit reader.

Golomb code:
- Unary-encoded quotient `q`: `q` zeros followed by a `1`
- Binary-encoded remainder `r` in `ceil(Int, log2(b))` bits

Returns `q * b + r`
"""
function read_golomb(r::BitReader, b::Int, ::Type{T}=UInt8) where {T<:Unsigned}
    b <= 0 && throw(ArgumentError("Golomb base must be >= 1"))

    # Step 1: Read unary prefix to get quotient
    q = 0
    while !read_bit(r)
        q += 1
    end

    # Step 2: Read fixed-width binary remainder
    m = ceil(Int, log2(b))
    rmod = read_value(r, m, T)

    # Step 3: Compute final value
    return T(q * b + rmod)
end

################################################################################
# Fibonacci code
################################################################################

"""
    write_fibonacci(w::BitWriter, n::T) where {T<:Unsigned}

Write `n` using Fibonacci coding (with '1' stop bit), using precomputed FIB_NUMBERS.
"""
function write_fibonacci(w::BitWriter, n::T) where {T<:Unsigned}
    n == 0 && throw(ArgumentError("Fibonacci code undefined for 0"))

    # find which Fibonacci numbers we need
    selected_indices = Int[]
    remaining = n
    i = length(FIB_NUMBERS)

    while i >= 1 && remaining > 0
        if FIB_NUMBERS[i] <= remaining
            push!(selected_indices, i)
            remaining -= FIB_NUMBERS[i]
            i -= 1  # skip next to avoid consecutive ones
        end
        i -= 1
    end

    # write bits in ascending order (index 1 first)
    reverse!(selected_indices)
    current_selected = 1
    
    for bit_pos in 1:length(FIB_NUMBERS)
        if current_selected <= length(selected_indices) && selected_indices[current_selected] == bit_pos
            write_bit(w, true)
            current_selected += 1
        else
            write_bit(w, false)
        end
        
        # if we've written all selected bits and the last bit was 1, we can add the stop bit
        if current_selected > length(selected_indices) && !isempty(selected_indices) && selected_indices[end] == bit_pos
            write_bit(w, true)  # stop bit
            break
        end
    end
end

"""
    read_fibonacci(r::BitReader, ::Type{T}=UInt8) where {T<:Unsigned}

Decode a Fibonacci code using precomputed FIB_NUMBERS.
"""
function read_fibonacci(r::BitReader, ::Type{T}=UInt8) where {T<:Unsigned}
    value = zero(T)
    prev = false
    i = 1

    while true
        b = read_bit(r)
        # if the bit is 1 and the previous bit is 1, we have reached the stop bit
        if b && prev
            break  # stop bit
        end
        # if the bit is 1, add the corresponding Fibonacci number to the value
        if b
            value += FIB_NUMBERS[i]
        end
        # update the previous bit
        prev = b
        # increment the index
        i += 1
        # if the index is greater than the length of the Fibonacci numbers, we have reached the end of the Fibonacci numbers
        if i > length(FIB_NUMBERS)
            throw(ArgumentError("Fibonacci code exceeds supported range (>= 2^40)"))
        end
    end

    return value
end

################################################################################
# Zeta code
################################################################################

"""
    write_zeta(w::BitWriter, v::T, k::Int) where {T<:Unsigned}

Write a zeta code for value `v` with parameter `k`.

Zeta coding of a natural number `v` depends on parameter `k` and is done in two steps:
1. Encode exponent h = floor(floor(log2(v))/k) using unary coding
2. Encode difference v - (2^k)^h using truncated binary coding
   with alphabet size (2^k)^(h+1) - (2^k)^h

@param w::BitWriter: the bitwriter
@param v::T: the value to encode (must be > 0)
@param k::Int: the zeta parameter (must be > 0)
"""
function write_zeta(w::BitWriter, v::T, k::Int) where {T<:Unsigned}
    if k <= 0
        throw(ArgumentError("k must be > 0"))
    end
    if v == 0
        throw(ArgumentError("Zeta coding is undefined for 0"))
    end
    
    # Special case: when k=1, zeta coding is identical to Elias gamma
    if k == 1
        write_elias_gamma(w, v)
        return
    end
    
    # Fast path: use precomputed bounds for common k values
    if k <= 8 && haskey(ZETA_H_BOUNDS, k)
        bounds = ZETA_H_BOUNDS[k]
        
        # Find h: bounds[i] = 2^(k*(i+1)) is first value requiring h = i+1
        # We want the number of bounds that v is >= to
        # This is equivalent to: count how many thresholds v exceeds
        h = 0
        for threshold in bounds
            if v >= threshold
                h += 1
            else
                break
            end
        end
        
        # Use precomputed power_base if available
        if h + 1 <= length(ZETA_POWER_BASES[k])
            power_base = T(ZETA_POWER_BASES[k][h + 1])  # +1 for 1-based indexing
        else
            # Fallback for very large h values
            power_base = T(1) << (k * h)
        end
    else
        # Fallback path: optimized h calculation using bit operations
        # Find h such that v is in range [2^(k*h), 2^(k*(h+1)))
        log2_v = sizeof(T) * 8 - leading_zeros(v) - 1  # floor(log2(v))
        h = div(log2_v, k)  # floor(log2_v / k)
        
        # Verify and adjust if needed (for edge cases)
        k_times_h = k * h
        power_base = T(1) << k_times_h
        if power_base > v
            h -= 1
            power_base >>= k 
        end
    end
    
    # Encode h using unary coding
    write_unary_coding(w, T(h))
    
    # Calculate remainder v - 2^(k*h)
    remainder = v - power_base
    
    # Calculate alphabet size = 2^(k*(h+1)) - 2^(k*h) = power_base * (2^k - 1)
    # Use UInt64 to avoid overflow
    alphabet_size_u64 = UInt64(power_base) * (UInt64(1 << k) - 1)
    
    # Check for overflow
    if alphabet_size_u64 > typemax(Int)
        throw(ArgumentError("Zeta coding error: alphabet_size overflow for v=$v, k=$k, h=$h, power_base=$power_base"))
    end
    
    alphabet_size = Int(alphabet_size_u64)
    
    # Encode remainder using truncated binary coding
    if alphabet_size > 1
        # Safety check: ensure remainder is in valid range
        if remainder >= alphabet_size
            throw(ArgumentError("Zeta coding error: remainder=$remainder >= alphabet_size=$alphabet_size for v=$v, k=$k, h=$h, power_base=$power_base"))
        end
        write_truncated_binary_coding(w, remainder, alphabet_size)
    end
end

"""
    read_zeta(r::BitReader, k::Int, ::Type{T}=UInt8) where {T<:Unsigned}

Read a zeta code from the bit reader.

@param r::BitReader: the bitreader
@param k::Int: the zeta parameter (must be > 0)
@param T::Type: the type to return (default: UInt8)
@return::T: the decoded value
"""
function read_zeta(r::BitReader, k::Int, ::Type{T}=UInt8) where {T<:Unsigned}
    h = read_unary_coding(r, false, UInt64)  # Read h as UInt64 to avoid overflow in calculations
    
    # Calculate power base: 2^(k*h)
    # Use UInt64 for all intermediate calculations to avoid overflow
    k_times_h = k * Int(h)
    power_base_u64 = UInt64(1) << k_times_h
    
    # Calculate alphabet size = 2^(k*(h+1)) - 2^(k*h) = power_base * (2^k - 1)
    # Use UInt64 to avoid overflow
    alphabet_size_u64 = power_base_u64 * (UInt64(1 << k) - 1)
    
    # Safety check for valid alphabet size
    if alphabet_size_u64 <= 0 || alphabet_size_u64 > typemax(Int)
        throw(ArgumentError("Zeta decoding error: invalid alphabet_size=$alphabet_size_u64 for h=$h, k=$k"))
    end
    
    alphabet_size = Int(alphabet_size_u64)
    
    # Read remainder using UInt64 to avoid overflow during reading
    remainder = read_truncated_binary_coding(r, alphabet_size, UInt64)
    
    # Final result: power_base + remainder
    result_u64 = power_base_u64 + UInt64(remainder)
    
    # Check if result fits in original type before converting
    if result_u64 > typemax(T)
        throw(ArgumentError("Zeta decoding error: result $result_u64 exceeds maximum value $(typemax(T)) for type $T. Debug info: h=$h, k=$k, power_base=$power_base_u64, remainder=$remainder"))
    end
    
    return T(result_u64)
end

################################################################################
# Fibonacci+Elias Delta Hybrid
################################################################################

"""
    write_fed(w::BitWriter, v::T, block_size::Int=FED_BLOCK_SIZE) where {T<:Unsigned}

Write a Fibonacci+Elias Delta hybrid code for value `v`.

Block-Based Fibonacci with Elias Delta Prefix (Hybrid Scheme):
- For a number N, compute q = floor(N/B) (block index) and r = N mod B (offset within block)
- Encode q using Elias delta code (efficient for large values)
- Encode r using Fibonacci code (fixed upper bound, short codes for small values)

@param w::BitWriter: the bitwriter
@param v::T: the value to encode (must be > 0)
@param block_size::Int: the block size B (default: FED_BLOCK_SIZE)
"""
function write_fed(w::BitWriter, v::T, block_size::Int=FED_BLOCK_SIZE) where {T<:Unsigned}
    v == 0 && throw(ArgumentError("FED encoding is undefined for 0"))
    block_size <= 0 && throw(ArgumentError("Block size must be > 0"))
    
    # Compute block index and offset
    q = div(v, block_size)  # block index
    r = v % block_size      # offset within block
    
    # Encode block index q using Elias delta
    # Note: If q == 0, we still need to encode it (Elias delta handles values ≥ 1)
    write_elias_delta(w, T(q + 1))  # shift by 1 since Elias delta requires v ≥ 1
    
    # Encode offset r using Fibonacci
    # Note: If r == 0, we need to handle it (Fibonacci requires v ≥ 1)
    write_fibonacci(w, T(r + 1))    # shift by 1 since Fibonacci requires v ≥ 1
end

"""
    read_fed(r::BitReader, ::Type{T}=UInt8, block_size::Int=FED_BLOCK_SIZE) where {T<:Unsigned}

Read a Fibonacci+Elias Delta hybrid code from the bit reader.

Block-Based Fibonacci with Elias Delta Prefix (Hybrid Scheme):
- First decode block index q from Elias delta code
- Then decode offset r from Fibonacci code  
- Return N = q * B + r where B is the block size

@param r::BitReader: the bit reader
@param T::Type: the type to return (default: UInt8)
@param block_size::Int: the block size B (default: FED_BLOCK_SIZE)
@return::T: the decoded value
"""
function read_fed(r::BitReader, ::Type{T}=UInt8, block_size::Int=FED_BLOCK_SIZE) where {T<:Unsigned}
    block_size <= 0 && throw(ArgumentError("Block size must be > 0"))
    
    # Decode block index q from Elias delta
    q = read_elias_delta(r, T) - T(1)  # unshift (we added 1 during encoding)
    
    # Decode offset r from Fibonacci
    r = read_fibonacci(r, T) - T(1)    # unshift (we added 1 during encoding)
    
    # Reconstruct original value: N = q * block_size + r
    # Use explicit type conversion to avoid promotion issues with custom types like UInt24
    result = T(UInt64(q) * UInt64(block_size) + UInt64(r))
    
    return result
end

################################################################################
# Interval compression (WebGraph-style)
################################################################################

"""
    compress_intervals(neighbors::Vector{T}, min_interval_length::Int=4) where {T<:Unsigned}

Compress consecutive neighbor ranges as intervals (WebGraph-style).
Returns (intervals, residuals) where intervals are consecutive ranges ≥ min_interval_length.

@param neighbors::Vector{T}: sorted list of neighbors
@param min_interval_length::Int: minimum length for interval compression (default: 4)
@return: (intervals, residuals) - intervals as [(start,length)], residuals as remaining values
"""
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

"""
    write_intervals_and_residuals(w::BitWriter, neighbors::Vector{T}, encoding::Symbol, min_interval_length::Int=4) where {T<:Unsigned}

Write neighbor list using interval compression + residual encoding.

@param w::BitWriter: the bitwriter
@param neighbors::Vector{T}: sorted list of neighbors  
@param encoding::Symbol: encoding for residuals
@param min_interval_length::Int: minimum interval length
"""
function write_intervals_and_residuals(w::BitWriter, neighbors::Vector{T}, encoding::Symbol, min_interval_length::Int=4) where {T<:Unsigned}
    isempty(neighbors) && return
    
    intervals, residuals = compress_intervals(neighbors, min_interval_length)
    
    # Write number of intervals (add 1 to avoid 0 with elias_gamma)
    write_encoded_value(w, T(length(intervals)) + T(1), :elias_gamma)
    
    # Write intervals: (start, length) pairs with delta encoding
    if !isempty(intervals)
        prev_start = T(0)
        for (start, length) in intervals
            # Delta encode start positions
            write_encoded_value(w, start - prev_start, encoding)
            # Encode interval length (already >= min_interval_length, add 1 to avoid 0)
            write_encoded_value(w, length - T(min_interval_length) + T(1), encoding)
            prev_start = start
        end
    end
    
    # Write residuals with delta encoding
    if !isempty(residuals)
        write_delta(w, residuals, encoding)
    end
end

"""
    read_intervals_and_residuals(r::BitReader, encoding::Symbol, min_interval_length::Int=4, ::Type{T}=UInt8) where {T<:Unsigned}

Read neighbor list from interval compression + residual encoding.

@param r::BitReader: the bitreader
@param encoding::Symbol: encoding used for residuals
@param min_interval_length::Int: minimum interval length used
@param T::Type: the type to return
@return::Vector{T}: reconstructed neighbor list
"""
function read_intervals_and_residuals(r::BitReader, encoding::Symbol, min_interval_length::Int=4, ::Type{T}=UInt8) where {T<:Unsigned}
    # Read number of intervals (subtract 1)
    num_intervals = read_encoded_value(r, :elias_gamma, T) - T(1)
    
    neighbors = T[]
    
    # Read intervals
    if num_intervals > 0
        prev_start = T(0)
        for _ in 1:num_intervals
            # Read delta-encoded start
            start_delta = read_encoded_value(r, encoding, T)
            start = prev_start + start_delta
            # Read length (subtract 1 and add back min_interval_length)
            length = read_encoded_value(r, encoding, T) - T(1) + T(min_interval_length)
            
            # Reconstruct interval
            for j in 0:(Int(length)-1)
                push!(neighbors, start + T(j))
            end
            prev_start = start
        end
    end
    
    # Read residuals
    residuals = read_delta(r, encoding, T)
    append!(neighbors, residuals)
    
    # Sort the combined result
    sort!(neighbors)
    return neighbors
end

################################################################################
# Run-length + delta code
################################################################################

"""
    write_delta(w::BitWriter, lst::Vector{T}, encoding::Symbol) where {T<:Unsigned}

Write a list of values to the bitwriter, using delta encoding.

@param w::BitWriter: the bitwriter
@param lst::Vector{T}: the list of values to write
@param encoding::Symbol: the compression coding to use (:elias_gamma, :elias_delta, :golomb, :fibonacci, :zeta)
"""
function write_delta(w::BitWriter, lst::Vector{T}, encoding::Symbol) where {T<:Unsigned}
    # if the list is empty, return
    isempty(lst) && return
    
    # delta encoding
    delta_lst = delta_encode_vector(lst)

    # write the first value (not delta encoded)
    # NB: we assume that the first value is not 0
    write_encoded_value(w, delta_lst[1], encoding)

    # write the rest of the values
    for i in 2:length(delta_lst)
        # NB: we shift by 1 to avoid zeros
        write_encoded_value(w, delta_lst[i] + T(1), encoding)
    end
end

"""
    read_delta(r::BitReader, encoding::Symbol, ::Type{T}=UInt8) where {T<:Unsigned}

Read a list of values from the bit reader, using delta encoding.

@param r::BitReader: the bit reader
@param encoding::Symbol: the compression coding to use (:elias_gamma, :elias_delta, :golomb, :fibonacci, :zeta)
@param T::Type: the type to return (default: UInt8)
@return::Vector{T}: the decoded list
"""
function read_delta(r::BitReader, encoding::Symbol, ::Type{T}=UInt8; max_elements::Union{Nothing,Int}=nothing, stop_value::Union{Nothing,T}=nothing) where {T<:Unsigned}
    lst = T[]

    # read the first value (not delta encoded)
    # NB: we assume that the first value is not 0
    try
        first_value = read_encoded_value(r, encoding, T)
        push!(lst, first_value)
    catch e
        # Empty list case: if we can't even read the first value
        if isa(e, EOFError) || isa(e, ErrorException)
            return T[]  # Return empty list
        else
            rethrow(e)
        end
    end

    # read the rest of the values
    while true
        # Check termination conditions
        if max_elements !== nothing && length(lst) >= max_elements
            break
        end

        # read the next value
        local value
        try
            value = read_encoded_value(r, encoding, T)
        catch e
            # End of stream - this is expected
            if isa(e, EOFError) || isa(e, ErrorException)
                break
            else
                rethrow(e)
            end
        end

        # Check for stop value before adding
        if stop_value !== nothing && value == stop_value
            break
        end

        # NB: shift back by 1
        push!(lst, lst[end] + value - T(1))
    end
    return lst
end

"""
    write_run_length_delta(w::BitWriter, compression::Symbol, lst::Vector{T}) where {T<:Unsigned}

Write a run-length + delta code to the bitwriter.

- [flag 0: 1 bit][delta: varint]  // delta
- [flag 1: 1 bit][value: varint][length: varint]   // run

@param w::BitWriter: the bitwriter
@param encoding::Symbol: the compression coding to use (:elias_gamma, :elias_delta, :golomb, :fibonacci, :zeta)
@param lst::Vector{T}: the list of vertex IDs to write
"""
function write_run_length_delta(w::BitWriter, encoding::Symbol, lst::Vector{T}) where {T<:Unsigned}
    # if the list is empty, return
    isempty(lst) && return

    # delta encoding
    delta_lst = delta_encode_vector(lst)

    # write the first value (not delta encoded)
    # NB: we assume that the first value is not 0
    write_encoded_value(w, delta_lst[1], encoding)

    # write the rest of the values
    i = 2
    while i <= length(delta_lst)
        current_value = delta_lst[i]

        # check for consecutive equal values (run)
        run_length = 1
        while i + run_length <= length(delta_lst) && delta_lst[i + run_length] == current_value
            run_length += 1
        end

        if run_length >= 3  # only use run-length for 3+ consecutive values
            # flag 1: run-length encoding
            write_bit(w, true)
            # encode the run length
            write_encoded_value(w, T(run_length), encoding)
            # encode the value shifted by 1 to avoid zeros
            write_encoded_value(w, current_value + T(1), encoding)
            i += run_length
        else
            # flag 0: delta encoding
            write_bit(w, false)
            # encode the delta value shifted by 1 to avoid zeros
            write_encoded_value(w, current_value + T(1), encoding)
            i += 1
        end
    end
end

"""
    read_run_length_delta(r::BitReader, encoding::Symbol, ::Type{T}=UInt8; max_elements::Union{Nothing,Int}=nothing, stop_value::Union{Nothing,T}=nothing) where {T<:Unsigned}

Read a run-length + delta code from the bit reader.

- [flag 0: 1 bit][delta: varint]  // delta
- [flag 1: 1 bit][value: varint][length: varint]   // run

@param r::BitReader: the bit reader
@param encoding::Symbol: the compression coding to use (:elias_gamma, :elias_delta, :golomb, :fibonacci, :zeta)
@param T::Type: the type to return (default: UInt8)
@param max_elements::Union{Nothing,Int}: maximum number of elements to read (for index mode)
@param stop_value::Union{Nothing,T}: value that signals end of list (for children mode)
@return::Vector{T}: the decoded list
"""
function read_run_length_delta(r::BitReader, encoding::Symbol, ::Type{T}=UInt8; max_elements::Union{Nothing,Int}=nothing, stop_value::Union{Nothing,T}=nothing) where {T<:Unsigned}
    # if the list is empty (immediate end of stream), return empty
    delta_lst = T[]
    
    # Pre-compute encoded stop value for comparison
    encoded_stop_value = nothing
    if stop_value !== nothing
        encoded_stop_value = get_encoded_value(stop_value, encoding, T)
    end
    
    try
        # read the first value (not delta encoded)
        # NB: we assume that the first value is not 0
        first_value = read_encoded_value(r, encoding, T)

        # if the first value is the stop value, return empty list
        if encoded_stop_value !== nothing && first_value == encoded_stop_value
            return T[]
        end

        push!(delta_lst, first_value)

        # read the rest of the values
        while true
            # Check termination conditions
            if max_elements !== nothing && length(delta_lst) >= max_elements
                break
            end
            
            # read the value flag (0: delta, 1: run-length)
            flag = read_bit(r)
            
            if flag
                # value flag 1: run-length encoding
                # Format: [run_length: varint][value: varint]
                run_length = Int(read_encoded_value(r, encoding, T))
                value = read_encoded_value(r, encoding, T)
                
                # add the run to the delta list (unshift the data value)
                for _ in 1:run_length
                    # NB: no need to check for max_elements here => we copy the whole run-length
                    push!(delta_lst, value - T(1))
                end
            else
                # value flag 0: delta encoding
                # Format: [delta_value: varint] (shifted by 1 to avoid zeros)
                # NB: first value is assumed to be not 0
                value = read_encoded_value(r, encoding, T)
                
                # Check for stop value (stop values are written as-is, not shifted)
                if encoded_stop_value !== nothing && value == encoded_stop_value
                    break
                end

                push!(delta_lst, value - T(1))
            end
        end
    catch e
        # end of stream reached, this is expected
        if !isa(e, EOFError) && !isa(e, ErrorException)
            rethrow(e)
        end
    end
    
    # reconstruct original values from delta encoding
    isempty(delta_lst) && return T[]
    
    original_lst = T[delta_lst[1]]
    for i in 2:length(delta_lst)
        push!(original_lst, original_lst[end] + delta_lst[i])
    end
    
    return original_lst
end

################################################################################
# Compressed graph data encoding
################################################################################

"""
    write_mix_encoded_list(w::BitWriter, delta_list::Vector{T}, encoding::Symbol, use_mix_mode::Bool=true) where {T<:Unsigned}

Write a list using mix mode (delta + run-length) following FORMAT_SPECS.md exactly.
This function implements the core mix mode format used by both index and children modes.

@param w::BitWriter: the bitwriter
@param delta_list::Vector{T}: the list of values to write
@param encoding::Symbol: the compression coding to use for values (:elias_gamma, :elias_delta, :golomb, :fibonacci, :zeta)
@param use_mix_mode::Bool: whether to use mix mode (default: true)
"""
function write_mix_encoded_list(w::BitWriter, delta_list::Vector{T}, encoding::Symbol, use_mix_mode::Bool=true) where {T<:Unsigned}
    # Handle empty list case first - don't call this function for empty lists!
    if isempty(delta_list)
        error("write_mix_encoded_list should not be called with empty lists - handle empty lists at caller level")
    end
    
    # write mix mode flag
    write_bit(w, use_mix_mode)
    
    # write the first value
    write_encoded_value(w, delta_list[1], encoding)
    
    if use_mix_mode
        # mix mode: write vertex flags and values
        i = 2
        while i <= length(delta_list)
            current_value = delta_list[i]
            
            # count consecutive occurrences of current_value
            run_length = 0
            j = i
            while j <= length(delta_list) && delta_list[j] == current_value
                run_length += 1
                j += 1
            end
            
            # run-length encoding is used for 3+ consecutive values
            if run_length >= 3
                # vertex flag = 1 (run-length encoding)
                write_bit(w, true)
                # run length
                write_encoded_value(w, T(run_length), encoding)
                # value
                write_encoded_value(w, current_value, encoding)
                i += run_length
            else
                # vertex flag = 0 (delta encoding)
                write_bit(w, false)
                # value
                write_encoded_value(w, current_value, encoding)
                i += 1
            end
        end
    else
        # delta-only: write remaining values directly
        for i in 2:length(delta_list)
            write_encoded_value(w, delta_list[i], encoding)
        end
    end
end

"""
    write_hybrid_mix_encoded_list(w::BitWriter, delta_list::Vector{T}, encoding::Symbol, use_run_length_and_interval::Bool=true, min_interval_length::Int=MIN_INTERVAL_LENGTH, is_children_mode::Bool=false) where {T<:Unsigned}

Write a delta-encoded list using hybrid mix mode encoding that adaptively combines:
- Delta encoding for irregular patterns
- Run-length encoding for repeated delta values  
- Interval encoding for consecutive sequences (reconstructed from deltas)

The encoding format uses bit flags to indicate section types:
- 0: Delta section (count + encoded delta values)
- 1: Section type flag follows
  - 10: Run-length section (count + value/length pairs)
  - 11: Interval section (count + start/length pairs)

@param w::BitWriter: the bitwriter
@param delta_list::Vector{T}: the delta-encoded neighbor list (like write_mix_encoded_list)
@param encoding::Symbol: the compression coding to use (:elias_gamma, :elias_delta, :golomb, :fibonacci, :zeta)
@param use_run_length_and_interval::Bool: enable hybrid features (default: true)
@param min_interval_length::Int: minimum length for interval compression (default: MIN_INTERVAL_LENGTH)
"""
function write_hybrid_mix_encoded_list(w::BitWriter, delta_list::Vector{T}, encoding::Symbol, use_run_length_and_interval::Bool=true, min_interval_length::Int=MIN_INTERVAL_LENGTH, is_children_mode::Bool=false) where {T<:Unsigned}
    if isempty(delta_list)
        error("write_hybrid_mix_encoded_list should not be called with empty lists")
    end
    
    # Decide whether to use hybrid mode for this list
    hybrid_active = use_run_length_and_interval && length(delta_list) > 1
    
    @debug "write_hybrid_mix_encoded_list: list_length=$(length(delta_list)), use_run_length_and_interval=$use_run_length_and_interval, hybrid_active=$hybrid_active"

    # Write hybrid mode flag
    write_bit(w, hybrid_active)
    
    # Write the first value (same as write_mix_encoded_list)
    write_encoded_value(w, delta_list[1], encoding)
    
    if !hybrid_active
        # Simple delta mode - write remaining values directly like write_mix_encoded_list
        for i in 2:length(delta_list)
            write_encoded_value(w, delta_list[i], encoding)
        end
        return
    end
    
    # length(delta_list) >= 2 here
    
    # Hybrid mode - analyze patterns in remaining delta values (2:end)
    remaining_deltas = delta_list[2:end]
    
    # Reconstruct original neighbors to find intervals. In children mode the
    # input deltas are shifted by +1, so unshift before reconstruction.
    local original_neighbors
    if is_children_mode
        if length(delta_list) == 0
            original_neighbors = T[]
        else
            unshifted = similar(delta_list)
            @inbounds begin
                unshifted[1] = delta_list[1] - one(T)
                for i in 2:length(delta_list)
                    unshifted[i] = delta_list[i] - one(T)
                end
            end
            original_neighbors = reconstruct_from_delta(unshifted)
        end
    else
        original_neighbors = reconstruct_from_delta(delta_list)
    end
    sections = analyze_delta_patterns_hybrid(remaining_deltas, original_neighbors[2:end], min_interval_length)
    
    # Count section types for debugging
    delta_count = count(s -> s.type == :delta, sections)
    run_length_count = count(s -> s.type == :run_length, sections) 
    interval_count = count(s -> s.type == :interval, sections)
    @debug "  Section summary: total=$(length(sections)), delta=$delta_count, run_length=$run_length_count, interval=$interval_count"
    
    # Write number of sections
    write_encoded_value(w, T(length(sections)), encoding)
    
    # Write each section with its type flag
    for section in sections
        if section.type == :delta
            # Delta section: flag=0, count, values
            write_bit(w, false)
            write_encoded_value(w, T(length(section.data)), encoding)
            for val in section.data
                write_encoded_value(w, val, encoding)
            end
        elseif section.type == :run_length
            # Run-length section: flag=1,0, count, value/length pairs
            write_bit(w, true)
            write_bit(w, false)
            write_encoded_value(w, T(length(section.data) ÷ 2), encoding)  # number of pairs
            for val in section.data
                write_encoded_value(w, val, encoding)
            end
        elseif section.type == :interval
            # Interval section: flag=1,1, count, start/length pairs
            write_bit(w, true)
            write_bit(w, true)
            write_encoded_value(w, T(length(section.data) ÷ 2), encoding)  # number of pairs
            for val in section.data
                write_encoded_value(w, val, encoding)
            end
        end
    end
end

"""
    analyze_delta_patterns_hybrid(delta_values::Vector{T}, original_neighbors::Vector{T}, min_interval_length::Int=MIN_INTERVAL_LENGTH) where {T<:Unsigned}

Analyze delta values (from position 2 onwards) to determine optimal hybrid encoding sections.
Uses original neighbors to detect intervals, then creates appropriate sections.
Returns a vector of encoding sections, each with a type (:delta, :run_length, :interval) and data.
"""
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
            # Create interval section: start and length
            @debug "  Found interval: start=$(original_neighbors[orig_i]), length=$interval_len at orig_pos=$orig_i"
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

"""
    find_consecutive_length(neighbors::Vector{T}, start::Int) where {T<:Unsigned}

Find the length of consecutive sequence starting at given position.
Returns the length of the consecutive sequence (1 if no sequence).
"""
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

"""
    find_run_length_patterns(delta_values::Vector{T}) where {T<:Unsigned}

Find run-length patterns in delta values and create appropriate sections.
Returns a vector of sections (delta or run_length).
"""
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
            @debug "  Found run-length: value=$current_val, length=$run_len at delta_pos=$i"
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

"""
    count_consecutive(values::Vector{T}, start::Int) where {T<:Unsigned}

Count consecutive occurrences of the same value starting at given position.
"""
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

"""
    read_hybrid_mix_encoded_list(r::BitReader, encoding::Symbol, mode::Symbol, ::Type{T}=UInt8; max_elements::Union{Int,Nothing}=nothing, stop_value::Union{T,Nothing}=nothing) where {T<:Unsigned}

Read a hybrid mix-encoded delta list that matches the write_hybrid_mix_encoded_list format.
Returns the reconstructed delta-encoded values (like read_mix_encoded_list).

@param r::BitReader: the bitreader
@param coding_scheme::Symbol: the coding scheme (:children or :index)
@param integer_encoding::Symbol: the integer encoding used (:elias_gamma, :elias_delta, :golomb, :fibonacci, :zeta)
@param T::Type: the unsigned integer type to use
@param max_elements: maximum elements to read (for index mode)
@param stop_value: stop value to watch for (for children mode)
@return Vector{T}: the reconstructed delta-encoded list
"""
function read_hybrid_mix_encoded_list(r::BitReader, coding_scheme::Symbol, integer_encoding::Symbol, ::Type{T}=UInt8; max_elements::Union{Int,Nothing}=nothing, stop_value::Union{T,Nothing}=nothing) where {T<:Unsigned}
    # Read hybrid mode flag
    use_run_length_and_interval = read_bit(r)
    
    # Read the first value (same as read_mix_encoded_list)
    first_value = read_encoded_value(r, integer_encoding, T)
    
    # Check if first value is the stop value for children mode
    if stop_value !== nothing && first_value == stop_value
        return T[]  # Empty list
    end
    
    result = T[first_value]
    elements_read = 1
    stop_consumed = false
    # Track last reconstructed original value on the fly to correctly
    # expand interval sections and keep bit alignment consistent.
    shift = coding_scheme == :children ? one(T) : zero(T)
    shift_int = Int(shift)
    last_original = coding_scheme == :children ? (Int(first_value) - shift_int) : Int(first_value)
    
    if !use_run_length_and_interval
        # Simple delta mode - read remaining values directly (same as read_mix_encoded_list)
        while true
            if max_elements !== nothing && elements_read >= max_elements
                break
            end
            
            # In children mode, a trailing stop may be prefixed by a vertex-flag 0.
            # If we detect it, consume and stop without trying to parse an encoded value.
            if stop_value !== nothing
                try
                    if !peek_bit(r)
                        _consume_children_trailing_stop(r, integer_encoding, T)
                        stop_consumed = true
                        break
                    end
                catch e
                    if isa(e, EOFError) || isa(e, ErrorException)
                        break
                    else
                        rethrow(e)
                    end
                end
            end

            value = read_encoded_value(r, integer_encoding, T)
            
            if stop_value !== nothing && value == stop_value
                stop_consumed = true
                break
            end
            
            push!(result, value)
            elements_read += 1
        end
        # In children mode, if not yet consumed above, consume trailing stop now.
        if stop_value !== nothing && !stop_consumed
            _consume_children_trailing_stop(r, integer_encoding, T)
            stop_consumed = true
        end
        # fall through to common post-processing (unshift + reconstruct)
    else
        # Hybrid mode - read sections
        num_sections = read_encoded_value(r, integer_encoding, T)
    
    for i in 1:num_sections
        section_flag = read_bit(r)
        
        if !section_flag
            # Delta section: flag=0, count, values
            count = read_encoded_value(r, integer_encoding, T)
            for j in 1:count
                if max_elements !== nothing && elements_read >= max_elements
                    break
                end
                
                value = read_encoded_value(r, integer_encoding, T)
                push!(result, value)
                elements_read += 1
                # update last original using decoded delta (accounting for shift)
                last_original += (Int(value) - shift_int)
            end
            
        else
            # Read second flag bit
            second_flag = read_bit(r)
            
            if !second_flag
                # Run-length section: flag=1,0, count, value/length pairs
                num_pairs = read_encoded_value(r, integer_encoding, T)
                
                for j in 1:num_pairs
                    value = read_encoded_value(r, integer_encoding, T)
                    length = read_encoded_value(r, integer_encoding, T)
                    
                    # Expand run-length back to individual values and update last_original
                    for k in 1:length
                        if max_elements !== nothing && elements_read >= max_elements
                            break
                        end
                        push!(result, value)
                        elements_read += 1
                        last_original += (Int(value) - shift_int)
                    end
                end
                
            else
                # Interval section: flag=1,1, count, start/length pairs
                num_pairs = read_encoded_value(r, integer_encoding, T)
                
                for j in 1:num_pairs
                    start = read_encoded_value(r, integer_encoding, T)
                    length = read_encoded_value(r, integer_encoding, T)
                    
                    # First encoded delta to reach the interval start from last_original
                    if max_elements !== nothing && elements_read >= max_elements
                        break
                    end
                    raw_delta = Int(start) - last_original
                    if raw_delta < 0
                        @warn "Interval start before last_original" last_original start elements_read
                        error("Invalid interval ordering: start < previous original")
                    end
                    first_delta = raw_delta + shift_int
                    push!(result, T(first_delta))
                    elements_read += 1
                    last_original = Int(start)
                    
                    # Add (length-1) consecutive 1 deltas (encoded with shift in children mode)
                    for k in 1:(Int(length) - 1)
                        if max_elements !== nothing && elements_read >= max_elements
                            break
                        end
                        push!(result, T(Int(one(T)) + shift_int))
                        elements_read += 1
                        last_original += 1
                    end
                end
            end
        end
        
        if max_elements !== nothing && elements_read >= max_elements
            break
        end
        end
    end
    
    # Consume trailing stop marker in children mode, matching writer behavior
    # In children mode, the writer may precede the encoded stop value with a
    # vertex-flag bit depending on the encoding path. Use the generic consumer
    # to remain aligned regardless of mix/hybrid mode details.
    if stop_value !== nothing && !stop_consumed
        _consume_children_trailing_stop(r, integer_encoding, T)
    end
    
    # Apply mode-specific unshifting (reverse of what write function does)
    if coding_scheme == :children && !isempty(result)
        # the write function applies delta_neighbors .+ T(1) in children mode, 
        # so we need to reverse this: result .- T(1)
        result = result .- T(1)
    end
    
    # Reconstruct original values from delta encoding
    if length(result) <= 1
        return result
    end
    
    original_list = T[result[1]]
    for i in 2:length(result)
        push!(original_list, original_list[end] + result[i])
    end
    
    return original_list
end

"""
    reconstruct_from_delta(delta_values::Vector{T}) where {T<:Unsigned}

Reconstruct original values from delta-encoded values.
"""
function reconstruct_from_delta(delta_values::Vector{T}) where {T<:Unsigned}
    if isempty(delta_values)
        return T[]
    end
    
    result = T[delta_values[1]]
    for i in 2:length(delta_values)
        push!(result, result[end] + delta_values[i])
    end
    
    return result
end

"""
    read_mix_encoded_list(r::BitReader, encoding::Symbol, mode::Symbol, ::Type{T}=UInt8; stop_value::Union{T,Nothing}=nothing, max_elements::Union{Int,Nothing}=nothing) where {T<:Unsigned}

Read a mix-encoded list that exactly matches the write_mix_encoded_list format.
This handles proper termination and value reconstruction according to the new format specifications.

@param r::BitReader: the bitreader
@param coding_scheme::Symbol: the coding scheme (:children or :index)
@param integer_encoding::Symbol: the integer encoding used (:elias_gamma, :elias_delta, :golomb, :fibonacci, :zeta)
@param T::Type: the unsigned integer type to use
@param max_elements: If provided, stop reading after this many elements have been reconstructed (>= 1 for index mode)
@param stop_value: If provided, stop reading when this encoded value is encountered in the stream (for children mode)
"""
function read_mix_encoded_list(r::BitReader, coding_scheme::Symbol, integer_encoding::Symbol, ::Type{T}=UInt8; max_elements::Union{Int,Nothing}=nothing, stop_value::Union{T,Nothing}=nothing) where {T<:Unsigned}
    # read mix mode flag
    use_mix_mode = read_bit(r)
    
    # read the first value from stream
    first_value = read_encoded_value(r, integer_encoding, T)
    
    # check if first value is the stop value (encoded)
    if stop_value !== nothing
        if first_value == stop_value
            # this is a stop value, indicating empty list
            return T[]
        end
    end
    
    delta_list = T[first_value]
    
    if use_mix_mode
        # mix mode: read vertex flags and values
        while true
            # check termination conditions
            if max_elements !== nothing && length(delta_list) >= max_elements
                break
            end
            
            # try to read vertex flag
            try
                vertex_flag = read_bit(r)
               
                # run-length encoding
                if vertex_flag 
                    # read the run length and the value
                    run_length = Int(read_encoded_value(r, integer_encoding, T))
                    value = read_encoded_value(r, integer_encoding, T)
                    
                    # add the run to the delta list
                    for _ in 1:run_length
                        if max_elements !== nothing && length(delta_list) >= max_elements
                            break
                        end
                        push!(delta_list, value)
                    end
                # delta-only mode
                else  
                    # read the delta value
                    delta_value = read_encoded_value(r, integer_encoding, T)
                    
                    # check for stop value (stop values are written as-is, not shifted)
                    if stop_value !== nothing
                        if delta_value == stop_value
                            # hit stop value, terminate reading
                            break
                        end
                    end
                    
                    push!(delta_list, delta_value)
                end
            catch e
                # end of available data - this is normal termination
                if isa(e, ErrorException) || isa(e, EOFError)
                    break
                else
                    rethrow(e)
                end
            end
        end
    # delta-only mode
    else
        # read remaining values
        while true
            # check termination conditions
            if max_elements !== nothing && length(delta_list) >= max_elements
                break
            end
            
            # try to read next value
            try
                # NB: no vertex flag is written in delta-only mode
                delta_value = read_encoded_value(r, integer_encoding, T)
                
                # check for stop value (stop values are written as-is, not shifted)
                if stop_value !== nothing
                    if delta_value == stop_value
                        # hit stop value, terminate reading
                        break
                    end
                end
                
                push!(delta_list, delta_value)
            catch e
                # end of available data - this is normal termination
                if isa(e, ErrorException) || isa(e, EOFError)
                    break
                else
                    rethrow(e)
                end
            end
        end
    end
    
    # apply mode-specific unshifting (reverse of what write function does)
    if coding_scheme == :children && !isempty(delta_list)
        # the write function applies delta_neighbors .+ T(1) in children mode, 
        # so we need to reverse this: delta_list .- T(1)
        delta_list = delta_list .- T(1)
    end
    
    # reconstruct original values from delta encoding
    if length(delta_list) <= 1
        return delta_list
    end
    
    original_list = T[delta_list[1]]
    for i in 2:length(delta_list)
        push!(original_list, original_list[end] + delta_list[i])
    end
    
    return original_list
end

"""
    _consume_children_trailing_stop(r::BitReader, encoding::Symbol, ::Type{T}=UInt8) where {T<:Unsigned}

Consume exactly one trailing stop marker in children mode, regardless of mix-mode.

End-of-vertex marker in children mode:
- if mix-mode is enabled, encoder writes: [vertex_flag=0: 1 bit] then [encoded stop value]
- if mix-mode is disabled, encoder writes: just [encoded stop value]

We do not know mix-mode here; instead, peek the next bit:
- if it is 0, consume it as the vertex flag, then read and discard one encoded value.
- if it is 1, directly read and discard one encoded value (the stop value begins with 1 for all supported codes).

@param r::BitReader: the bitreader
@param integer_encoding::Symbol: the integer encoding used (:elias_gamma, :elias_delta, :golomb, :fibonacci, :zeta)
@param T::Type: the unsigned integer type to use
"""
function _consume_children_trailing_stop(r::BitReader, integer_encoding::Symbol, ::Type{T}=UInt8) where {T<:Unsigned}
    # safely peek next bit; if at EOF, nothing to consume
    local next_is_zero
    try
        next_is_zero = !peek_bit(r)
    catch e
        if isa(e, EOFError) || isa(e, ErrorException)
            return
        else
            rethrow(e)
        end
    end

    # if next bit is 0, it's the vertex flag for mix-mode; consume it
    if next_is_zero
        _ = read_bit(r)
    end
    # read and discard exactly one encoded value (the stop value)
    _ = read_encoded_value(r, integer_encoding, T)
    return
end

"""
    estimate_reference_encoding_cost(ref_id::T, copy_bitmap::Vector{Bool}, residuals::Vector{T}, encoding::Symbol, mode::Symbol) where {T<:Unsigned}

Estimate the bit cost of encoding a vertex using reference encoding.
Returns the estimated number of bits required.

@param ref_id::T: the reference ID
@param copy_bitmap::Vector{Bool}: the copy bitmap
@param residuals::Vector{T}: the residuals
@param coding_scheme::Symbol: the coding scheme (:children or :index)
@param integer_encoding::Symbol: the integer encoding used (:elias_gamma, :elias_delta, :golomb, :fibonacci, :zeta)
@return::Int: the estimated cost in bits
"""
function estimate_reference_encoding_cost(ref_id::T, copy_bitmap::Vector{Bool}, residuals::Vector{T}, coding_scheme::Symbol, integer_encoding::Symbol) where {T<:Unsigned}
    cost = 0
    
    # ref_id cost
    cost += estimate_encoded_value_cost(T(ref_id), integer_encoding)
    
    # bitmap_len cost  
    cost += estimate_encoded_value_cost(T(length(copy_bitmap)), integer_encoding)
    
    # bitmap bits
    cost += length(copy_bitmap)
    
    # residuals_flag
    cost += 1
    
    if !isempty(residuals)
        # residuals_len for index mode
        if coding_scheme == :index
            cost += estimate_encoded_value_cost(T(length(residuals)), integer_encoding)
        end
        
        # residuals data (using traditional mix encoding)
        cost += estimate_mix_encoding_cost(residuals, integer_encoding)
    end
    
    return cost
end

"""
    estimate_hybrid_mix_encoding_cost(delta_neighbors::Vector{T}, encoding::Symbol, use_mix_mode::Bool) where {T<:Unsigned}

Estimate the bit cost of encoding a vertex using hybrid mix encoding.
Returns the estimated number of bits required.

@param delta_neighbors::Vector{T}: the delta-encoded neighbor list
@param integer_encoding::Symbol: the integer encoding used (:elias_gamma, :elias_delta, :golomb, :fibonacci, :zeta)
@param use_mix_mode::Bool: whether to use mix mode (default: true)
@return::Int: the estimated cost in bits
"""
function estimate_hybrid_mix_encoding_cost(delta_neighbors::Vector{T}, integer_encoding::Symbol, use_mix_mode::Bool) where {T<:Unsigned}
    if isempty(delta_neighbors)
        return 0
    end
    
    cost = 0
    hybrid_active = use_mix_mode && length(delta_neighbors) > 1
    
    # Hybrid mode flag
    cost += 1
    
    # First value
    cost += estimate_encoded_value_cost(delta_neighbors[1], integer_encoding)
    
    if !hybrid_active
        # Simple delta mode - remaining values
        for i in 2:length(delta_neighbors)
            cost += estimate_encoded_value_cost(delta_neighbors[i], integer_encoding)
        end
    else
        # Estimate hybrid sections cost (simplified estimation)
        remaining_deltas = delta_neighbors[2:end]
        
        # Rough estimation: assume mostly delta encoding with some patterns
        # This could be more sophisticated by actually analyzing patterns
        
        # Number of sections (estimated)
        estimated_sections = max(1, length(remaining_deltas) ÷ 4)  # Rough estimate
        cost += estimate_encoded_value_cost(T(estimated_sections), integer_encoding)
        
        # Section costs (simplified - assume mostly delta sections)
        for delta in remaining_deltas
            cost += 1  # section flag
            cost += estimate_encoded_value_cost(T(1), integer_encoding)  # count
            cost += estimate_encoded_value_cost(delta, integer_encoding)  # value
        end
    end
    
    return cost
end

"""
    estimate_mix_encoding_cost(values::Vector{T}, encoding::Symbol) where {T<:Unsigned}

Estimate the bit cost of traditional mix encoding.

@param values::Vector{T}: the values to encode
@param integer_encoding::Symbol: the integer encoding used (:elias_gamma, :elias_delta, :golomb, :fibonacci, :zeta)
@return::Int: the estimated cost in bits
"""
function estimate_mix_encoding_cost(values::Vector{T}, integer_encoding::Symbol) where {T<:Unsigned}
    if isempty(values)
        return 1  # mix_mode_flag
    end
    
    cost = 1  # mix_mode_flag
    cost += estimate_encoded_value_cost(values[1], integer_encoding)  # first value
    
    # Simplified: assume delta-only mode for estimation
    for i in 2:length(values)
        cost += 1  # vertex_flag
        cost += estimate_encoded_value_cost(values[i], integer_encoding)
    end
    
    return cost
end

"""
    estimate_encoded_value_cost(value::T, encoding::Symbol) where {T<:Unsigned}

Estimate the bit cost of encoding a single value with the given encoding scheme.
This is a simplified estimation - actual cost depends on bit-level details.

@param value::T: the value to encode
@param integer_encoding::Symbol: the integer encoding used (:elias_gamma, :elias_delta, :golomb, :fibonacci, :zeta)
@return::Int: the estimated cost in bits
"""
function estimate_encoded_value_cost(value::T, integer_encoding::Symbol) where {T<:Unsigned}
    if value == 0
        return 1  # Special case
    end
    
    # Rough bit cost estimates for different encodings
    if integer_encoding == :fibonacci
        # Fibonacci encoding roughly log_phi(n) + log_phi(n)/phi bits
        return ceil(Int, log(2, max(1, value)) * 1.44) + 2
    elseif integer_encoding == :elias_gamma
        # Elias gamma: 2⌊log2(n)⌋ + 1 bits  
        return 2 * floor(Int, log(2, max(1, value))) + 1
    elseif integer_encoding == :elias_delta
        # Elias delta: roughly log2(n) + 2log2(log2(n)) bits
        log_val = max(1, log(2, max(1, value)))
        return ceil(Int, log_val + 2 * log(2, max(1, log_val)))
    elseif integer_encoding == :zeta
        # Zeta coding with k=4: roughly log4(n) + k bits
        return ceil(Int, log(4, max(1, value))) + 4
    else
        # Default fallback
        return ceil(Int, log(2, max(1, value))) + 2
    end
end

"""
    write_compressed_graph_data(w::BitWriter, neighbor_lists::Dict{T,Vector{T}}, encoding::Symbol, mode::Symbol=:children, use_mix_mode::Bool=true, reference_enabled::Bool=true) where {T<:Unsigned}

Write compressed graph data with hybrid mix encoding (delta + run-length + intervals) and optional reference encoding.

# Format specifications: see FORMAT_SPECS.md

@param w::BitWriter: the bitwriter
@param neighbor_lists::Dict{T,Vector{T}}: neighbor lists for each vertex
@param coding_scheme::Symbol: the coding scheme to use (:children, :index)
@param integer_encoding::Symbol: the integer encoding to use (:elias_gamma, :elias_delta, :golomb, :fibonacci, :zeta)
@param use_mix_mode::Bool: whether to use mix mode (delta + run-length) (default: true)
@param reference_enabled::Bool: whether to enable reference encoding (default: true)
@param recursive_reference::Bool: whether to enable recursive reference encoding (default: true)
"""
function write_compressed_graph_data(w::BitWriter, neighbor_lists::Dict{T,Vector{T}}, coding_scheme::Symbol=:children, integer_encoding::Symbol=:elias_delta, use_mix_mode::Bool=true, reference_enabled::Bool=true, recursive_reference::Bool=true) where {T<:Unsigned}
    # get the number of vertices
    # NB: we assume that the vertices are numbered from 1 to vs
    vs = length(keys(neighbor_lists))

    # track which vertices have been encoded (can serve as references)
    # Use sliding window approach like WebGraph for better locality
    reference_window = T[]  # Circular buffer of recent vertices
    potential_references = Set{T}()  # For compatibility with existing functions
    
    # Build StreamingInvertedIndex for efficient reference candidate lookup
    ref_index = nothing
    vertex_to_stream_id = Dict{T,T}()
    stream_id_to_vertex = Dict{T,T}()
    # Dense vector mapping for faster lookup (index is stream_id)
    stream_id_to_vertex_vec = T[]
    ref_workspace = nothing
 
    if reference_enabled
        @debug "Building StreamingInvertedIndex for reference encoding..."
        index_build_start = time()
        
        # Find maximum neighbor ID to size the index appropriately
        @debug "  Finding maximum neighbor ID..."
        max_neighbor_start = time()
        max_neighbor_id = T(0)
        for neighbors in values(neighbor_lists)
            if !isempty(neighbors)
                max_neighbor_id = max(max_neighbor_id, maximum(neighbors))
            end
        end
        max_neighbor_time = time() - max_neighbor_start
        @debug "  Max neighbor ID: $max_neighbor_id (found in $(round(max_neighbor_time, digits=3))s)"
        
        # Create index with sufficient column space
        @debug "  Creating StreamingInvertedIndex..."
        index_create_start = time()
        ref_index = Index.StreamInvertedIndex{T}(max_neighbor_id)
        index_create_time = time() - index_create_start
        @debug "  Index created in $(round(index_create_time, digits=3))s"
        
        # Add candidates with degree >= REF_V_MIN_DEGREE to the index
        @debug "  Adding candidates to index (degree >= $REF_V_MIN_DEGREE)..."
        candidates_start = time()
        candidate_count = 0
        progress_interval = max(1, length(neighbor_lists) ÷ 20)  # Progress every 5%
        
        for (i, (vertex_id, neighbors)) in enumerate(neighbor_lists)
            if length(neighbors) >= REF_V_MIN_DEGREE
                stream_id = Index.add_candidate!(ref_index, neighbors)
                vertex_to_stream_id[vertex_id] = stream_id
                stream_id_to_vertex[stream_id] = vertex_id  # Pre-build reverse mapping
                # Ensure dense vector mapping capacity and set
                if stream_id > length(stream_id_to_vertex_vec)
                    resize!(stream_id_to_vertex_vec, Int(stream_id))
                end
                stream_id_to_vertex_vec[Int(stream_id)] = vertex_id
                candidate_count += 1
            end
            
            # progress logging
            if i % progress_interval == 0
                elapsed = time() - candidates_start
                @debug "    Progress: $i/$(length(neighbor_lists)) vertices ($(round(100*i/length(neighbor_lists), digits=1))%), $(candidate_count) candidates, $(round(elapsed, digits=3))s elapsed"
            end
        end
        
        candidates_time = time() - candidates_start
        
        # create reusable workspace for overlap computations
        @debug "  Creating reusable workspace..."
        workspace_start = time()
        ref_workspace = Index.OverlapWorkspace(ref_index)
        workspace_time = time() - workspace_start
        
        index_build_time = time() - index_build_start
        
        @debug "  Candidates added: $candidate_count/$(length(neighbor_lists)) ($(round(100*candidate_count/length(neighbor_lists), digits=1))%)"
        @debug "  Candidate addition time: $(round(candidates_time, digits=3))s"
        @debug "  Workspace creation time: $(round(workspace_time, digits=3))s"
        @debug "  Total index build time: $(round(index_build_time, digits=3))s"
    end

    # write reference flag
    write_bit(w, reference_enabled)

    # index section (if :index mode)
    if coding_scheme == :index
        for v in 1:vs
            # write the out-degree shifted by 1 using the specified encoding
            out_degree = T(length(neighbor_lists[T(v)]) + 1)  # shift by 1
            write_encoded_value(w, out_degree, integer_encoding)
        end
    end

    # data section
    @debug "Starting data section encoding..."
    data_section_start = time()
    progress_interval = max(1, vs ÷ 100)  # Progress every 1%
    reference_queries = 0
    references_found = 0
    # Stats counters
    reference_chosen_count = 0        # vertices encoded using reference (after cost compare)
    hybrid_chosen_count = 0           # vertices encoded using hybrid (no ref or ref lost)
    # Fast/slow lookup counters (top-level + recursive)
    fast_hit_count_ref = Ref(0)
    slow_hit_count_ref = Ref(0)
    # Recursive reference usage counters
    recursive_vertices_count = 0      # number of vertices where recursive residual encoding was used at least once
    recursive_event_total_ref = Ref(0) # total number of recursive selections across all levels
    
    # Dense mask of available references (by vertex id)
    available_mask = falses(Int(vs) + 1)
    # Workspace to build reference data without allocations
    ref_build_work = RefBuildWorkspace{T}()
    
    # Helper function to add vertex to reference window (WebGraph-style sliding window)
    function add_to_reference_window!(vertex::T)
        push!(reference_window, vertex)
        push!(potential_references, vertex)
        available_mask[vertex] = true
        
        # Maintain sliding window of size REF_WINDOW_SIZE
        if length(reference_window) > REF_WINDOW_SIZE
            # Remove oldest vertex from window
            old_vertex = popfirst!(reference_window)
            delete!(potential_references, old_vertex)
            available_mask[old_vertex] = false
        end
    end

    for v in 1:vs
        current_neighbors = neighbor_lists[T(v)]

        # handle empty lists
        if isempty(current_neighbors)
            if coding_scheme == :index
                # in index mode, we don't need to write anything for empty lists
                # since we already know from the out-degrees that there are no neighbors
                continue
            # children mode
            else
                # in children mode, empty lists must participate in reference encoding flow
                if reference_enabled
                    # children_flag = 0 (no existing reference for empty list)
                    write_bit(w, false)
                end
                
                # write mix mode flag
                # mix mode flag = 0 (delta-only)
                write_bit(w, true)
                
                # NB: if the list is empty, the stop value is the first value
                # and there is no vertex flag

                # stop value
                write_encoded_value(w, T(1), integer_encoding)  
                
                # continue to the next vertex
                continue
            end
        end
        
        # sort neighbors in place
        sorted_neighbors = sort!(current_neighbors)
        
        # compute delta encoding
        delta_neighbors = delta_encode_vector(sorted_neighbors)
        
        # apply mode-specific value shifting
        if coding_scheme == :children
            # shift all delta values by 1 to avoid zeros
            delta_neighbors = delta_neighbors .+ T(1)
        end
        if reference_enabled
            # check if we should use reference encoding
            use_reference = false
            best_ref = nothing
            copy_bitmap = Bool[]
            residuals = T[]
            
            if !isempty(potential_references)
                reference_queries += 1
                # try fast path when we have dense mappings
                used_fast = false
                best_ref = find_best_reference_fast(sorted_neighbors, ref_index, stream_id_to_vertex_vec, ref_workspace, available_mask)
                @debug "find_best_reference_fast result" best_ref=best_ref fast_lookup_available=(ref_index !== nothing)
                # if fast path fails (e.g., before any available set), fall back to original Set-based filter
                # TODO: to be removed
                if best_ref === nothing
                    @debug "fast path failed, falling back to original Set-based filter"
                    best_ref = find_best_reference(sorted_neighbors, ref_index, stream_id_to_vertex, ref_workspace, potential_references)
                else
                    used_fast = true
                end
                if best_ref !== nothing
                    # Record whether fast or slow hit (top-level)
                    if used_fast
                        fast_hit_count_ref[] += 1
                    else
                        slow_hit_count_ref[] += 1
                    end
                    references_found += 1
                    ref_neighbors = sort(neighbor_lists[best_ref])
                    # build copy bitmap + residuals using reusable workspace
                    create_reference_data!(ref_build_work, sorted_neighbors, ref_neighbors)
                    copy_bitmap = ref_build_work.copy_bitmap
                    residuals = ref_build_work.residuals
                    # reference selection is now threshold-based in find_best_reference
                    use_reference = true
                end
            end

            if use_reference
                # Compare reference encoding cost vs hybrid mix encoding cost
                ref_cost = estimate_reference_encoding_cost(best_ref, copy_bitmap, residuals, coding_scheme, integer_encoding)
                hybrid_cost = estimate_hybrid_mix_encoding_cost(delta_neighbors, integer_encoding, use_mix_mode)
                
                # Choose the more efficient encoding
                use_reference_final = ref_cost <= hybrid_cost
                @debug "Encoding comparison: vertex=$v, ref_cost=$ref_cost, hybrid_cost=$hybrid_cost, chosen=$(use_reference_final ? "reference" : "hybrid"), savings=$(abs(ref_cost - hybrid_cost)) bits"
                
                if use_reference_final
                    reference_chosen_count += 1
                    @debug "CHOICE: Reference encoding chosen: vertex=$v, ref_id=$best_ref, copy_bitmap_len=$(length(copy_bitmap)), residuals_len=$(length(residuals))"
                    # children_flag = 1 (reference mode)
                    write_bit(w, true)
                else
                    # Use hybrid encoding instead of reference
                    hybrid_chosen_count += 1
                    @debug "CHOICE: Hybrid mix encoding chosen over reference: vertex=$v, neighbors_count=$(length(sorted_neighbors))"
                    # children_flag = 0 (hybrid mix mode)
                    write_bit(w, false)
                    # Use hybrid encoding
                    write_hybrid_mix_encoded_list(w, delta_neighbors, integer_encoding, use_mix_mode)
                    add_to_reference_window!(T(v))
                end
                
                if use_reference_final
                    # write reference data
                    write_encoded_value(w, T(best_ref), integer_encoding)  # ref_id
                    write_encoded_value(w, T(length(copy_bitmap)), integer_encoding)  # bitmap_len
                    for bit in copy_bitmap
                        write_bit(w, bit)
                    end
                    
                    # write residuals
                    if !isempty(residuals)
                        write_bit(w, true)  # residuals_flag = 1

                        # residuals_len need to be written only for index mode
                        if coding_scheme == :index
                            write_encoded_value(w, T(length(residuals)), integer_encoding)  # residuals_len
                        end
                    
                        # encode residuals with same format as delta encoding
                        residual_deltas = delta_encode_vector(residuals)

                        # shift residuals by 1 to avoid zeros in children mode
                        if coding_scheme == :children
                            residual_deltas = residual_deltas .+ T(1)
                        end
                        
                        # write residuals using recursive reference encoding or standard mix mode
                        if recursive_reference
                            @debug "writing recursive reference residuals with use_mix_mode=$use_mix_mode, coding_scheme=$coding_scheme"
                            # Per-vertex recursive event counter
                            local_recursive_events = Ref(0)
                            write_recursive_reference_residuals(w, residual_deltas, coding_scheme, integer_encoding, use_mix_mode, neighbor_lists, ref_index, stream_id_to_vertex_vec, ref_workspace, available_mask, ref_build_work, potential_references, add_to_reference_window!, fast_hit_count_ref, slow_hit_count_ref, local_recursive_events)
                            # Aggregate recursive stats
                            if local_recursive_events[] > 0
                                recursive_vertices_count += 1
                                recursive_event_total_ref[] += local_recursive_events[]
                            end
                        else
                            @debug "writing standard mix mode residuals (recursive_reference=false)"
                            # Write recursive_flag = 0 (use standard mix mode)
                            write_bit(w, false)
                            if coding_scheme == :children
                                write_hybrid_mix_encoded_list(w, residual_deltas, integer_encoding, use_mix_mode, MIN_INTERVAL_LENGTH, true)
                            else
                                # For index mode, write the residuals directly
                                for residual in residual_deltas
                                    write_encoded_value(w, residual, integer_encoding)
                                end
                            end
                        end
                    else
                        # residuals_flag = 0 (no residuals)
                        write_bit(w, false)  
                    end
                    
                    add_to_reference_window!(T(v))
                end
            # no relevant reference: use mix mode
            else
                hybrid_chosen_count += 1
                @debug "CHOICE: Hybrid mix encoding (no reference found): vertex=$v, neighbors_count=$(length(sorted_neighbors))"
                # children_flag = 0 (hybrid mix mode)
                write_bit(w, false)
                # delta encode the neighbors
                delta_neighbors = delta_encode_vector(sorted_neighbors)
                if coding_scheme == :children
                    delta_neighbors = delta_neighbors .+ T(1)  # shift by 1 to avoid zeros
                end
                @debug "writing hybrid mix encoded list with use_mix_mode=$use_mix_mode, MIN_INTERVAL_LENGTH=$MIN_INTERVAL_LENGTH, coding_scheme=$coding_scheme"
                # use_mix_mode enables hybrid encoding (run-length + interval)
                write_hybrid_mix_encoded_list(w, delta_neighbors, integer_encoding, use_mix_mode, MIN_INTERVAL_LENGTH, coding_scheme == :children)
                add_to_reference_window!(T(v))
            end
        # reference disabled: use hybrid compression method
        else
            # delta encode the neighbors
            delta_neighbors = delta_encode_vector(sorted_neighbors)
            if coding_scheme == :children
                delta_neighbors = delta_neighbors .+ T(1)  # shift by 1 to avoid zeros
            end
            @debug "writing hybrid mix encoded list with use_mix_mode=$use_mix_mode, MIN_INTERVAL_LENGTH=$MIN_INTERVAL_LENGTH, coding_scheme=$coding_scheme"
            # use_mix_mode enables hybrid encoding (run-length + interval)
            write_hybrid_mix_encoded_list(w, delta_neighbors, integer_encoding, use_mix_mode, MIN_INTERVAL_LENGTH, coding_scheme == :children)
        end
        
        # progress logging
        if v % progress_interval == 0
            elapsed = time() - data_section_start
            @debug "  Progress: $v/$vs vertices ($(round(100*v/vs, digits=1))%), $(length(reference_window))/$(REF_WINDOW_SIZE) window refs, $reference_queries queries, $references_found found, $(round(elapsed, digits=3))s elapsed"
        end
        
        # write stop value after each vertex list in children mode
        if coding_scheme == :children
            # NB: the list is not empty, so the stop value is not the first value
            # and a vertex flag is needed
            if use_mix_mode
                # vertex flag = 0 (delta-only)
                write_bit(w, false)  
            end
            # stop value
            write_encoded_value(w, T(1), integer_encoding)  
        end
    end
    
    data_section_time = time() - data_section_start
    @debug "Data section completed!"
    @debug "  Total vertices processed: $vs"
    @debug "  Final potential references: $(length(potential_references))"
    @debug "  Reference queries made: $reference_queries"
    @debug "  References successfully found: $references_found ($(round(references_found > 0 ? 100*references_found/reference_queries : 0.0, digits=1))%)"
    @debug "  Data section encoding time: $(round(data_section_time, digits=3))s"
    @debug "  Average time per vertex: $(round(1000*data_section_time/vs, digits=3))ms"

    # Info-level statistics summary
    total_vertices = vs
    total_encoded_hybrid = hybrid_chosen_count
    total_encoded_reference = reference_chosen_count
    total_lookup_hits = fast_hit_count_ref[] + slow_hit_count_ref[]
    fast_ratio = total_lookup_hits > 0 ? (100 * fast_hit_count_ref[] / total_lookup_hits) : 0.0
    slow_ratio = total_lookup_hits > 0 ? (100 * slow_hit_count_ref[] / total_lookup_hits) : 0.0
    ref_ratio = total_vertices > 0 ? (100 * total_encoded_reference / total_vertices) : 0.0
    hyb_ratio = total_vertices > 0 ? (100 * total_encoded_hybrid / total_vertices) : 0.0
    recursive_vertices_ratio = total_encoded_reference > 0 ? (100 * recursive_vertices_count / total_encoded_reference) : 0.0

    @info "Encoding summary: reference chosen=$total_encoded_reference ($(round(ref_ratio, digits=1))%), hybrid chosen=$total_encoded_hybrid ($(round(hyb_ratio, digits=1))%)"
    @info "Lookup summary: fast hits=$(fast_hit_count_ref[]), slow hits=$(slow_hit_count_ref[]), total hits=$total_lookup_hits (fast=$(round(fast_ratio, digits=1))%, slow=$(round(slow_ratio, digits=1))%)"
    @info "Recursive summary: vertices with recursive residuals=$recursive_vertices_count ($(round(recursive_vertices_ratio, digits=1))% of ref-chosen), total recursive selections=$(recursive_event_total_ref[])"
end

"""
    read_compressed_graph_data(r::BitReader, vs::T, coding_scheme::Symbol=:children, integer_encoding::Symbol=:elias_delta, ::Type{T}=UInt8) where {T<:Unsigned}

Read compressed graph data with hybrid mix mode (delta + run-length + intervals) and optional reference mode.

# Format specifications: see FORMAT_SPECS.md

@param r::BitReader: the bitreader
@param vs::T: the number of vertices in the graph
@param coding_scheme::Symbol: the coding scheme to use (:children, :index)
@param integer_encoding::Symbol: the integer encoding to use (:elias_gamma, :elias_delta, :golomb, :fibonacci, :zeta)
@param T::Type: the type to return (default: UInt8)
@return::Dict{T,Vector{T}}: the decoded neighbor lists
"""
function read_compressed_graph_data(r::BitReader, vs::T, coding_scheme::Symbol=:children, integer_encoding::Symbol=:elias_delta, ::Type{T}=UInt8) where {T<:Unsigned}
    neighbor_lists = Dict{T,Vector{T}}()
    
    
    # 1. read reference flag  
    reference_enabled = read_bit(r)
    
    # 2. index section (if :index mode)
    out_degrees = T[]
    if coding_scheme == :index
        out_degrees = Vector{T}(undef, vs)
        for v in 1:vs
            encoded_degree = read_encoded_value(r, integer_encoding, T)
            out_degrees[v] = encoded_degree - T(1)  # unshift
        end
    end
    
    # 3. data section
    if coding_scheme == :index
        # index mode: read each vertex based on out-degree
        for v in 1:vs
            if out_degrees[v] == 0
                # empty list - no data written for this vertex
                neighbor_lists[T(v)] = T[]
            else
                # read exactly out_degrees[v] neighbors
                expected_neighbors = Int(out_degrees[v])
                
                # reference mode
                if reference_enabled
                    # read children flag
                    children_flag = read_bit(r)
                    if children_flag  # reference mode
                        # read reference data
                        ref_id = read_encoded_value(r, integer_encoding, T)
                        bitmap_len = read_encoded_value(r, integer_encoding, T)
                        
                        copy_bitmap = Bool[]
                        for _ in 1:Int(bitmap_len)
                            push!(copy_bitmap, read_bit(r))
                        end
                        
                        residuals_flag = read_bit(r)
                        residuals = T[]
                        
                        if residuals_flag
                            # in index mode, we know the length of the residuals
                            residuals_len = read_encoded_value(r, integer_encoding, T)
                            if residuals_len > 0
                                # NB: the original values are reconstructed from the delta encoding
                                residuals = read_recursive_reference_residuals(r, coding_scheme, integer_encoding, T; max_elements=Int(residuals_len), neighbor_lists=neighbor_lists)
                            end
                        end
                        
                        # reconstruct from reference
                        if haskey(neighbor_lists, ref_id)
                            ref_neighbors = neighbor_lists[ref_id]
                            current_neighbors = reconstruct_from_reference(ref_neighbors, copy_bitmap, residuals)
                        else
                            error("Invalid reference ID: $ref_id")
                        end
                    else  # read the appropriate format
                        # For now, assume hybrid format since that's what we're using in the new version
                        current_neighbors = read_hybrid_mix_encoded_list(r, coding_scheme, integer_encoding, T; max_elements=expected_neighbors)
                    end
                else
                    # reference disabled - read the appropriate format
                    # For now, assume hybrid format since that's what we're using in the new version
                    current_neighbors = read_hybrid_mix_encoded_list(r, coding_scheme, integer_encoding, T; max_elements=expected_neighbors)
                end
                neighbor_lists[T(v)] = current_neighbors
            end
        end
    # children mode
    else        
        for v in 1:vs
            current_neighbors = T[]
            
            # first try to read the vertex data
            vertex_data_read = false
            try
                if reference_enabled
                    # read children flag
                    children_flag = read_bit(r)
                    
                    # reference mode: children mode
                    if children_flag  
                        # read reference data
                        ref_id = read_encoded_value(r, integer_encoding, T)
                        bitmap_len = read_encoded_value(r, integer_encoding, T)
                        
                        copy_bitmap = Bool[]
                        for _ in 1:Int(bitmap_len)
                            push!(copy_bitmap, read_bit(r))
                        end
                        
                        residuals_flag = read_bit(r)
                        residuals = T[]
                        
                        # in children mode, residuals are terminated by a stop value
                        if residuals_flag
                            # NB: 
                            # - the unshifting is done in the read_recursive_reference_residuals function
                            # - the original values are reconstructed from the delta encoding
                            residuals = read_recursive_reference_residuals(r, coding_scheme, integer_encoding, T; stop_value=T(1), neighbor_lists=neighbor_lists)
                        end
                        
                        # reconstruct from reference
                        if haskey(neighbor_lists, ref_id)
                            ref_neighbors = neighbor_lists[ref_id]
                            current_neighbors = reconstruct_from_reference(ref_neighbors, copy_bitmap, residuals)
                        else
                            error("Invalid reference ID: $ref_id")
                        end
                        
                        # children mode: if there were no residuals, consume the trailing stop value now
                        if !residuals_flag
                            # Consume exactly one trailing stop marker independent of mix-mode
                            _consume_children_trailing_stop(r, integer_encoding, T)
                        end
                    # reference mode: read the appropriate format
                    else 
                        # The writer used either write_mix_encoded_list or write_hybrid_mix_encoded_list
                        # Both start with a flag, but we need to determine which format to use
                        # For now, assume hybrid format since that's what we're using in the new version
                        current_neighbors = read_hybrid_mix_encoded_list(r, coding_scheme, integer_encoding, T; stop_value=T(1))
                    end
                # reference disabled: read the appropriate format  
                else
                    # The writer used either write_mix_encoded_list or write_hybrid_mix_encoded_list
                    # For now, assume hybrid format since that's what we're using in the new version
                    current_neighbors = read_hybrid_mix_encoded_list(r, coding_scheme, integer_encoding, T; stop_value=T(1))
                end
                
                # successfully read vertex data
                neighbor_lists[T(v)] = current_neighbors
                vertex_data_read = true
            catch e
                # EOF while reading vertex data: current vertex and all remaining vertices are empty
                if isa(e, EOFError) || isa(e, ErrorException)
                    for remaining_v in v:vs
                        neighbor_lists[T(remaining_v)] = T[]
                    end
                    break
                else
                    rethrow(e)
                end
            end 
            # NB: stop values are consumed by read_mix_encoded_list when stop_value parameter is provided
        end
    end
    
    return neighbor_lists
end

"""
    find_best_reference(target::Vector{T}, ref_index::Union{Index.StreamInvertedIndex{T}, Nothing}, 
                       stream_id_to_vertex::Dict{T,T}, work::Union{Index.OverlapWorkspace{T}, Nothing},
                       available::Set{T}) where {T<:Unsigned}

Find the best reference vertex for the target neighbor list using StreamingInvertedIndex with pre-built workspace.
This optimized version reuses workspace and reverse mapping to eliminate per-query overhead.

@param target::Vector{T}: the neighbor list to encode  
@param ref_index::Union{StreamInvertedIndex{T}, Nothing}: pre-built inverted index with reference candidates
@param stream_id_to_vertex::Dict{T,T}: pre-built mapping from stream IDs to vertex IDs
@param work::Union{OverlapWorkspace{T}, Nothing}: pre-allocated workspace for overlap computation
@param available::Set{T}: vertices that can serve as references (already encoded)
@return::Union{T, Nothing}: the best reference vertex ID or nothing if no good reference
"""
function find_best_reference(target::Vector{T}, ref_index::Union{Index.StreamInvertedIndex{T}, Nothing}, 
                           stream_id_to_vertex::Dict{T,T}, work::Union{Index.OverlapWorkspace{T}, Nothing},
                           available::Set{T}) where {T<:Unsigned}
    # skip if no index available, no workspace, no candidates available, or target too small
    if ref_index === nothing || work === nothing || isempty(available) || length(target) <= 2
        return nothing
    end
    
    # compute overlaps with all candidates in the index (reusing workspace)
    counts, touched = Index.overlap!(ref_index, target, work)
    
    if isempty(touched)
        return nothing
    end
    
    # find best overlap among available references only
    best_ref = nothing
    best_overlap = 0
    
    for stream_id in touched
        overlap_count = counts[stream_id]
        
        # look up corresponding vertex ID using pre-built mapping
        if haskey(stream_id_to_vertex, stream_id)
            vertex_id = stream_id_to_vertex[stream_id]
            
            # skip if this candidate is not available as reference or overlap too low
            if vertex_id in available && overlap_count > best_overlap
                best_overlap = overlap_count
                best_ref = vertex_id
            end
        end
    end
    
    # apply threshold: use reference only if overlap exceeds REF_ENCODING_TH
    if best_overlap >= REF_ENCODING_TH
        return best_ref
    end
    
    return nothing
end

"""
    find_best_reference_fast(target::Vector{T}, ref_index::Union{Index.StreamInvertedIndex{T}, Nothing},
                             stream_id_to_vertex_vec::Union{Vector{T},Nothing},
                             work::Union{Index.OverlapWorkspace{T}, Nothing},
                             available_mask::Union{Vector{Bool},Nothing}) where {T<:Unsigned}

optimized best-reference selection:
- uses a dense `stream_id -> vertex` vector mapping when available (avoids Dict lookups)
- uses a dense `available_mask` Bool vector (avoids Set membership checks)
- falls back to `nothing` parameters to skip optimization.
"""
function find_best_reference_fast(target::Vector{T}, ref_index::Union{Index.StreamInvertedIndex{T}, Nothing},
                                 stream_id_to_vertex_vec::Union{Vector{T},Nothing},
                                 work::Union{Index.OverlapWorkspace{T}, Nothing},
                                 available_mask::Union{AbstractVector{Bool},Nothing}) where {T<:Unsigned}
    # skip if no index available, no workspace, or target too small
    if ref_index === nothing || work === nothing || isempty(target)
        return nothing
    end

    # compute overlaps with all candidates in the index (reusing workspace)
    counts, touched = Index.overlap!(ref_index, target, work)
    if isempty(touched)
        return nothing
    end

    # find best overlap among available references only
    best_ref = nothing
    best_overlap = 0

    # if dense mappings are available, use them
    if stream_id_to_vertex_vec !== nothing && available_mask !== nothing
        @inbounds for stream_id in touched
            # map to vertex id; skip invalid zeros
            if stream_id <= length(stream_id_to_vertex_vec)
                v = stream_id_to_vertex_vec[stream_id]
                if v != zero(T) && v <= length(available_mask) && available_mask[Int(v)]
                    ov = counts[stream_id]
                    if ov > best_overlap
                        best_overlap = ov
                        best_ref = v
                    end
                end
            end
        end
    end

    return best_overlap >= REF_ENCODING_TH ? best_ref : nothing
end

"""
    create_reference_data(target::Vector{T}, reference::Vector{T}) where {T<:Unsigned}

Create copy bitmap and residuals for reference encoding.

@param target::Vector{T}: the target neighbor list
@param reference::Vector{T}: the reference neighbor list  
@return::Tuple{Vector{Bool}, Vector{T}}: (copy_bitmap, residuals)
"""
function create_reference_data(target::Vector{T}, reference::Vector{T}) where {T<:Unsigned}
    # backward-compatible wrapper allocating fresh buffers
    work = RefBuildWorkspace{T}()
    create_reference_data!(work, target, reference)
    return work.copy_bitmap, work.residuals
end

"""
    create_reference_data!(work::RefBuildWorkspace{T}, target::Vector{T}, reference::Vector{T}) where {T<:Unsigned}

Create copy bitmap and residuals for reference encoding assuming both inputs are sorted.
This version reuses buffers in `work` to avoid per-call allocations and avoids Sets.
"""
function create_reference_data!(work::RefBuildWorkspace{T}, target::Vector{T}, reference::Vector{T}) where {T<:Unsigned}
    # clear the workspace
    empty!(work.copy_bitmap)
    empty!(work.residuals)
    
    # pointers in reference and target
    i = 1  # pointer in reference
    j = 1  # pointer in target
    nref = length(reference)
    ntgt = length(target)

    # build copy bitmap by scanning both sorted arrays once
    while i <= nref && j <= ntgt
        rv = reference[i]
        tv = target[j]
        if rv == tv
            push!(work.copy_bitmap, true)
            i += 1; j += 1
        elseif rv < tv
            push!(work.copy_bitmap, false)
            i += 1
        else # tv < rv
            # target element not in reference → residual
            push!(work.residuals, tv)
            j += 1
        end
    end

    # remaining reference elements are not in target
    while i <= nref
        push!(work.copy_bitmap, false)
        i += 1
    end

    # remaining target elements are residuals
    while j <= ntgt
        push!(work.residuals, target[j])
        j += 1
    end

    return work
end

"""
    reconstruct_from_reference(reference::Vector{T}, copy_bitmap::Vector{Bool}, residuals::Vector{T}) where {T<:Unsigned}

Reconstruct neighbor list from reference data.

@param reference::Vector{T}: the reference neighbor list
@param copy_bitmap::Vector{Bool}: bitmap indicating which reference elements to copy
@param residuals::Vector{T}: additional neighbors not in reference
@return::Vector{T}: the reconstructed neighbor list
"""
function reconstruct_from_reference(reference::Vector{T}, copy_bitmap::Vector{Bool}, residuals::Vector{T}) where {T<:Unsigned}
    # assumes `reference` is sorted and `residuals` is sorted.
    # produces a sorted merge without an extra sort step.
    # first, stream selected reference entries, then merge with residuals.

    # collect selected items from reference in order
    selected = T[]
    @inbounds for i in 1:min(length(copy_bitmap), length(reference))
        if copy_bitmap[i]
            push!(selected, reference[i])
        end
    end

    # if one side is empty, return the other directly
    isempty(selected) && return copy(residuals)
    isempty(residuals) && return selected

    # merge two sorted lists
    result = Vector{T}(undef, length(selected) + length(residuals))
    i = 1; j = 1; k = 1
    n1 = length(selected); n2 = length(residuals)
    @inbounds while i <= n1 && j <= n2
        if selected[i] <= residuals[j]
            result[k] = selected[i]; i += 1
        else
            result[k] = residuals[j]; j += 1
        end
        k += 1
    end
    @inbounds while i <= n1
        result[k] = selected[i]; i += 1; k += 1
    end
    @inbounds while j <= n2
        result[k] = residuals[j]; j += 1; k += 1
    end
    return result
end

################################################################################
# Recursive Reference Encoding
################################################################################

"""
    find_best_reference_set(neighbors::Vector{T}, neighbor_lists::Dict{T,Vector{T}}, potential_references::Set{T}) where {T<:Unsigned}

Find the best reference for the given neighbors using a simple Set-based approach.
This is a fallback for recursive reference encoding when the fast path is not available.

@param neighbors::Vector{T}: sorted neighbors to find a reference for
@param neighbor_lists::Dict{T,Vector{T}}: existing neighbor lists
@param potential_references::Set{T}: set of potential reference vertices  
@return::Union{T,Nothing}: best reference vertex ID, or nothing if no good reference found
"""
function find_best_reference_set(neighbors::Vector{T}, neighbor_lists::Dict{T,Vector{T}}, potential_references::Set{T}) where {T<:Unsigned}
    isempty(neighbors) && return nothing
    
    neighbors_set = Set(neighbors)
    best_ref = nothing
    best_overlap = 0
    min_overlap_threshold = max(1, length(neighbors) ÷ 4)  # Require at least 25% overlap
    
    for ref_id in potential_references
        if haskey(neighbor_lists, ref_id)
            ref_neighbors = neighbor_lists[ref_id]
            # Calculate overlap
            overlap = 0
            for neighbor in ref_neighbors
                if neighbor in neighbors_set
                    overlap += 1
                end
            end
            
            # Update best reference if this one has better overlap
            if overlap > best_overlap && overlap >= min_overlap_threshold
                best_overlap = overlap
                best_ref = ref_id
            end
        end
    end
    
    return best_ref
end

"""
    write_recursive_reference_residuals(w::BitWriter, residuals::Vector{T}, coding_scheme::Symbol, integer_encoding::Symbol, use_mix_mode::Bool, neighbor_lists::Dict{T,Vector{T}}, ref_index::Union{Index.StreamInvertedIndex{T}, Nothing}, stream_id_to_vertex_vec::Union{Vector{T},Nothing}, ref_workspace::Union{Index.OverlapWorkspace{T}, Nothing}, available_mask::Union{AbstractVector{Bool},Nothing}, ref_build_work::RefBuildWorkspace{T}, potential_references::Set{T}, add_to_reference_window!::Function; fast_hit_count_ref::Base.RefValue{Int}=Ref(0), slow_hit_count_ref::Base.RefValue{Int}=Ref(0), recursive_event_count_ref::Base.RefValue{Int}=Ref(0)) where {T<:Unsigned}

Write residuals using recursive reference encoding. This function attempts to find a reference
for the residuals themselves, creating a recursive reference structure when beneficial.

@param w::BitWriter: the bitwriter
@param residuals::Vector{T}: residuals to encode (already delta-encoded and shifted if needed)
@param coding_scheme::Symbol: coding scheme to use (:children or :index)
@param integer_encoding::Symbol: integer encoding to use
@param use_mix_mode::Bool: whether to use hybrid mix mode
@param neighbor_lists::Dict{T,Vector{T}}: existing neighbor lists for reference lookup
@param ref_index::Union{Index.StreamInvertedIndex{T}, Nothing}: reference index for fast lookups
@param stream_id_to_vertex_vec::Union{Vector{T},Nothing}: mapping for fast reference lookups
@param ref_workspace::Union{Index.OverlapWorkspace{T}, Nothing}: workspace for reference operations
@param available_mask::Union{AbstractVector{Bool},Nothing}: mask of available references
@param ref_build_work::RefBuildWorkspace{T}: workspace for reference building
@param potential_references::Set{T}: set of potential reference vertices
@param add_to_reference_window!::Function: function to add vertices to reference window
"""
function write_recursive_reference_residuals(w::BitWriter, residuals::Vector{T}, coding_scheme::Symbol, integer_encoding::Symbol, use_mix_mode::Bool, neighbor_lists::Dict{T,Vector{T}}, ref_index::Union{Index.StreamInvertedIndex{T}, Nothing}, stream_id_to_vertex_vec::Union{Vector{T},Nothing}, ref_workspace::Union{Index.OverlapWorkspace{T}, Nothing}, available_mask::Union{AbstractVector{Bool},Nothing}, ref_build_work::RefBuildWorkspace{T}, potential_references::Set{T}, add_to_reference_window!::Function,
    fast_hit_count_ref::Base.RefValue{Int}=Ref(0), slow_hit_count_ref::Base.RefValue{Int}=Ref(0), recursive_event_count_ref::Base.RefValue{Int}=Ref(0)) where {T<:Unsigned}
    
    isempty(residuals) && return
    
    # Reconstruct original neighbors from delta-encoded residuals
    original_residuals = if coding_scheme == :children
        # Unshift the residuals (they were shifted by +1 in children mode)
        unshifted = residuals .- T(1)
        reconstruct_from_delta(unshifted)
    else
        reconstruct_from_delta(residuals)
    end
    
    # Try to find a reference for the residuals
    best_ref = nothing
    if !isempty(potential_references)
        # Track which path found the reference
        used_fast = false
        best_ref = find_best_reference_fast(original_residuals, ref_index, stream_id_to_vertex_vec, ref_workspace, available_mask)
        @debug "Recursive: find_best_reference_fast result" best_ref=best_ref fast_lookup_available=(ref_index !== nothing)
        if best_ref === nothing
            # Note: stream_id_to_vertex is not available in the recursive context, fallback to Set-based approach
            # This is acceptable as recursive references are less frequent
            @debug "Recursive: fast path failed, using slower Set-based approach"
            best_ref = find_best_reference_set(original_residuals, neighbor_lists, potential_references)
        else
            used_fast = true
        end
        # Update fast/slow counters for recursive context
        if best_ref !== nothing
            if used_fast
                fast_hit_count_ref[] += 1
            else
                slow_hit_count_ref[] += 1
            end
        end
    end
    
    use_recursive_reference = false
    if best_ref !== nothing
        # Calculate costs to determine if recursive reference is beneficial
        ref_neighbors = sort(neighbor_lists[best_ref])
        create_reference_data!(ref_build_work, original_residuals, ref_neighbors)
        copy_bitmap = ref_build_work.copy_bitmap
        recursive_residuals = ref_build_work.residuals
        
        # Estimate costs
        recursive_ref_cost = estimate_recursive_reference_cost(best_ref, copy_bitmap, recursive_residuals, coding_scheme, integer_encoding)
        hybrid_cost = estimate_hybrid_mix_encoding_cost(residuals, integer_encoding, use_mix_mode)
        
        # Change comparison to <= so recursion wins on ties
        use_recursive_reference = recursive_ref_cost <= hybrid_cost
        @debug "Recursive reference decision: cost_recursive=$recursive_ref_cost, cost_hybrid=$hybrid_cost, chosen=$(use_recursive_reference ? "recursive" : "hybrid")"
    end
    
    if use_recursive_reference
        # Count this recursive selection
        recursive_event_count_ref[] += 1
        # Write recursive_flag = 1
        write_bit(w, true)
        
        # Write recursive reference data (same format as regular reference)
        ref_neighbors = sort(neighbor_lists[best_ref])
        create_reference_data!(ref_build_work, original_residuals, ref_neighbors)
        copy_bitmap = ref_build_work.copy_bitmap
        recursive_residuals = ref_build_work.residuals

        write_encoded_value(w, T(best_ref), integer_encoding)  # recursive_ref_id
        write_encoded_value(w, T(length(copy_bitmap)), integer_encoding)  # bitmap_len
        for bit in copy_bitmap
            write_bit(w, bit)
        end
        
        # Write recursive residuals
        if !isempty(recursive_residuals)
            write_bit(w, true)  # recursive_residuals_flag = 1
            
            if coding_scheme == :index
                write_encoded_value(w, T(length(recursive_residuals)), integer_encoding)  # recursive_residuals_len
            end
            
            # Encode recursive residuals
            recursive_residual_deltas = delta_encode_vector(recursive_residuals)
            if coding_scheme == :children
                recursive_residual_deltas = recursive_residual_deltas .+ T(1)
            end
            
            # Recursively try reference encoding on these residuals
            write_recursive_reference_residuals(w, recursive_residual_deltas, coding_scheme, integer_encoding, use_mix_mode, neighbor_lists, ref_index, stream_id_to_vertex_vec, ref_workspace, available_mask, ref_build_work, potential_references, add_to_reference_window!, fast_hit_count_ref, slow_hit_count_ref, recursive_event_count_ref)
        else
            write_bit(w, false)  # recursive_residuals_flag = 0
        end
    else
        # Write recursive_flag = 0, then use standard hybrid mix encoding
        write_bit(w, false)
        
        if coding_scheme == :children
            write_hybrid_mix_encoded_list(w, residuals, integer_encoding, use_mix_mode, MIN_INTERVAL_LENGTH, true)
        else
            # For index mode, write the residuals directly
            for residual in residuals
                write_encoded_value(w, residual, integer_encoding)
            end
        end
    end
end

"""
    read_recursive_reference_residuals(r::BitReader, encoding::Symbol, mode::Symbol, ::Type{T}; max_elements::Union{Int,Nothing}=nothing, stop_value::Union{T,Nothing}=nothing, neighbor_lists::Dict{T,Vector{T}}) where {T<:Unsigned}

Read residuals that may have been recursively reference-encoded.

@param r::BitReader: the bitreader
@param coding_scheme::Symbol: coding scheme to use (:children or :index)
@param integer_encoding::Symbol: integer encoding to use
@param T::Type: the type to return
@param max_elements::Union{Int,Nothing}: maximum number of elements (for index mode)
@param stop_value::Union{T,Nothing}: stop value (for children mode)
@param neighbor_lists::Dict{T,Vector{T}}: neighbor lists for reference reconstruction
@return::Vector{T}: the decoded residuals
"""
function read_recursive_reference_residuals(r::BitReader, coding_scheme::Symbol, integer_encoding::Symbol, ::Type{T}; max_elements::Union{Int,Nothing}=nothing, stop_value::Union{T,Nothing}=nothing, neighbor_lists::Dict{T,Vector{T}}) where {T<:Unsigned}
    
    # Read recursive_flag
    recursive_flag = read_bit(r)
    
    if recursive_flag
        # Recursive reference encoding
        recursive_ref_id = read_encoded_value(r, integer_encoding, T)
        bitmap_len = read_encoded_value(r, integer_encoding, T)
        
        copy_bitmap = Bool[]
        for _ in 1:Int(bitmap_len)
            push!(copy_bitmap, read_bit(r))
        end
        
        recursive_residuals_flag = read_bit(r)
        recursive_residuals = T[]
        
        if recursive_residuals_flag
            if coding_scheme == :index
                recursive_residuals_len = read_encoded_value(r, integer_encoding, T)
                if recursive_residuals_len > 0
                    # Recursively read the residuals
                    recursive_residuals = read_recursive_reference_residuals(r, coding_scheme, integer_encoding, T; max_elements=Int(recursive_residuals_len), neighbor_lists=neighbor_lists)
                end
            else  # children mode
                # Recursively read the residuals
                recursive_residuals = read_recursive_reference_residuals(r, coding_scheme, integer_encoding, T; stop_value=stop_value, neighbor_lists=neighbor_lists)
            end
        end
        
        # Reconstruct from recursive reference
        if haskey(neighbor_lists, recursive_ref_id)
            ref_neighbors = neighbor_lists[recursive_ref_id]
            return reconstruct_from_reference(ref_neighbors, copy_bitmap, recursive_residuals)
        else
            error("Invalid recursive reference ID: $recursive_ref_id")
        end
    else
        # Standard hybrid mix encoding
        if coding_scheme == :children
            return read_hybrid_mix_encoded_list(r, coding_scheme, integer_encoding, T; stop_value=stop_value)
        else
            # Index mode - read residuals directly
            residuals = T[]
            if max_elements !== nothing
                for _ in 1:max_elements
                    push!(residuals, read_encoded_value(r, integer_encoding, T))
                end
            end
            return residuals
        end
    end
end

"""
    estimate_recursive_reference_cost(ref_id::T, copy_bitmap::Vector{Bool}, recursive_residuals::Vector{T}, encoding::Symbol, mode::Symbol) where {T<:Unsigned}

Estimate the bit cost of recursive reference encoding.

@param ref_id::T: reference ID
@param copy_bitmap::Vector{Bool}: copy bitmap
@param recursive_residuals::Vector{T}: residuals that would be recursively encoded
@param coding_scheme::Symbol: coding scheme to use (:children or :index)
@param integer_encoding::Symbol: integer encoding to use
@return::Int: estimated cost in bits
"""
function estimate_recursive_reference_cost(ref_id::T, copy_bitmap::Vector{Bool}, recursive_residuals::Vector{T}, coding_scheme::Symbol, integer_encoding::Symbol) where {T<:Unsigned}
    cost = 0
    
    # Recursive flag (1 bit)
    cost += 1
    
    # Reference ID cost
    cost += estimate_encoded_value_cost(ref_id, integer_encoding)
    
    # Bitmap length cost
    cost += estimate_encoded_value_cost(T(length(copy_bitmap)), integer_encoding)
    
    # Bitmap cost
    cost += length(copy_bitmap)
    
    # Recursive residuals flag (1 bit)
    cost += 1
    
    if !isempty(recursive_residuals)
        if coding_scheme == :index
            # Residuals length cost
            cost += estimate_encoded_value_cost(T(length(recursive_residuals)), integer_encoding)
        end
        
        # For simplicity, assume recursive residuals use hybrid encoding
        # (In practice, this would require recursive cost estimation)
        delta_residuals = delta_encode_vector(recursive_residuals)
        if coding_scheme == :children
            delta_residuals = delta_residuals .+ T(1)
        end
        cost += estimate_hybrid_mix_encoding_cost(delta_residuals, integer_encoding, true)
    end
    
    return cost
end

end # module Compression
