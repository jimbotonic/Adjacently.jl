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
using ..IO: BitWriter, BitReader, write_bit, write_bits, read_bit, read_bits, flush, write_value, read_value, flush_bitwriter
using ..Constants: FIB_NUMBERS, BUFFER_SIZE, ZETA_H_BOUNDS, ZETA_POWER_BASES, ZETA_BASE, GOLOMB_BASE

# Export the functions we want to make available
export write_unary_coding,
       write_truncated_binary_coding,
       write_encoded_value,
       read_encoded_value,
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
       write_golomb,
       read_golomb,
       write_fibonacci,
       read_fibonacci,
       write_zeta,
       read_zeta,
       write_run_length_delta,
       read_run_length_delta,
       write_reference_encoding,
       read_reference_encoding

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
    
    # calculate threshold u = 2^(k+1) - n
    u = convert(T, (T(1) << (k + 1)) - n)
    
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
    
    # Calculate threshold u = 2^(k+1) - n
    u = (1 << (k + 1)) - n
    
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
@param compression::Symbol: the compression coding to use (:elias_gamma, :elias_delta, :golomb, :fibonacci, :zeta)
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
    else
        throw(ArgumentError("Invalid compression code: $compression"))
    end
end

"""
    read_encoded_value(r::BitReader, compression::Symbol, ::Type{T}=UInt8) where {T<:Unsigned}

Read a value from the bitreader using the specified compression code.

@param r::BitReader: the bitreader
@param compression::Symbol: the compression coding to use (:elias_gamma, :elias_delta, :golomb, :fibonacci, :zeta)
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
    else
        throw(ArgumentError("Invalid compression code: $compression"))
    end
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
    
    # Calculate alphabet size = 2^(k*(h+1)) - 2^(k*h) = 2^(k*h) * (2^k - 1)
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
    h = read_unary_coding(r, false, T)
    power_base = T(1) << (k * h)
    remainder = read_truncated_binary_coding(r, Int(power_base << k) - Int(power_base), T)
    return power_base + remainder
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
    
    try
        # read the first value (not delta encoded)
        # NB: we assume that the first value is not 0
        first_value = read_encoded_value(r, encoding, T)

        # if the first value is the stop value, return empty list
        if stop_value !== nothing && first_value == stop_value
            return T[]
        end

        push!(delta_lst, first_value)

        # read the rest of the values
        while true
            # Check termination conditions
            if max_elements !== nothing && length(delta_lst) >= max_elements
                break
            end
            
            # read the flag bit (0: delta, 1: run-length)
            flag = read_bit(r)
            
            if flag
                # flag 1: run-length encoding
                # Format: [run_length: varint][value: varint]
                run_length = Int(read_encoded_value(r, encoding, T))
                value = read_encoded_value(r, encoding, T)
                
                # add the run to the delta list (unshift the data value)
                for _ in 1:run_length
                    # NB: no need to check for max_elements here => we copy the whole run-length
                    push!(delta_lst, value - T(1))
                end
            else
                # flag 0: delta encoding
                # Format: [delta_value: varint] (shifted by 1 to avoid zeros)
                # NB: first value is assumed to be not 0
                value = read_encoded_value(r, encoding, T)
                
                # Check for stop value (stop values are written as-is, not shifted)
                if stop_value !== nothing && value == stop_value
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
# Reference encoding
################################################################################

"""
    write_reference_encoding(w::BitWriter, neighbor_lists::Dict{T,Vector{T}}, encoding::Symbol, mode::Symbol=:children, enable_reference::Bool=true) where {T<:Unsigned}

Write reference-encoded neighbor lists with run-length and delta compression to the bitwriter.

# Modes:

## Index Mode (:index)
- First writes out-degrees using run-length+delta encoding
- Then writes adjacency lists using delta encoding (no value shifting, no stop values)
- If enable_reference=true: uses reference encoding when beneficial for adjacency lists

## Children Mode (:children) 
- Writes adjacency lists using delta encoding with values shifted by 1 (minimum value 1)
- Adds stop value T(1) at end of each list except the last
- If enable_reference=true: uses reference encoding when beneficial

# Encoding Format:

## Global Header:
```
[Global Main Flag: 1 bit] [mode-specific data...]
```

- **Global Main Flag = 1**: Reference encoding is enabled
- **Global Main Flag = 0**: Reference encoding is disabled (direct encoding only)

## When Reference Encoding is Enabled (Global Main Flag = 1):

For each vertex, the format is:
```
[Vertex Flag: 1 bit] [Encoding-specific data...]
```

### Vertex Flag Values:
- **0**: Direct delta encoding for this vertex
- **1**: Reference encoding for this vertex

### Format for Vertex Flag = 0 (Direct Encoding):
- **Index mode**: delta encoded values without shifts or stops
- **Children mode**: delta encoded values shifted by 1 with stop values

### Format for Vertex Flag = 1 (Reference Encoding):
```
[1] [ref_id: varint] [bitmap_length: varint] [bitmap: N bits] [residuals_flag: 1 bit] [residuals?]
```

Where:
- **ref_id**: ID of reference vertex (encoded with specified compression)
- **bitmap_length**: Number of bits in copy bitmap (encoded with specified compression)  
- **bitmap**: Raw bits indicating which reference elements to copy (1=copy, 0=skip)
- **residuals_flag**: 1 bit (1=has residuals, 0=no residuals)
- **residuals**: If residuals_flag=1, delta encoded residual elements (following mode-specific encoding)

## When Reference Encoding is Disabled (Global Main Flag = 0):

All vertices use direct delta encoding without any per-vertex flags:
- **Index mode**: delta encoded values without shifts or stops
- **Children mode**: delta encoded values shifted by 1 with stop values

@param w::BitWriter: the bitwriter
@param neighbor_lists::Dict{T,Vector{T}}: neighbor lists for each vertex
@param encoding::Symbol: the compression coding to use for values (:elias_gamma, :elias_delta, :golomb, :fibonacci, :zeta)
@param mode::Symbol: the mode to use (:children, :index)
@param enable_reference::Bool: whether to enable reference encoding (default: true)
"""
function write_reference_encoding(w::BitWriter, neighbor_lists::Dict{T,Vector{T}}, encoding::Symbol, mode::Symbol=:children, enable_reference::Bool=true) where {T<:Unsigned}
    # get the number of vertices
    # NB: we assume that the vertices are numbered from 1 to vs
    vs = length(keys(neighbor_lists))

    # track which vertices have been encoded (can serve as references)
    potential_references = Set{T}()

    # write the main flag
    if enable_reference
        write_bit(w, true)  # main flag = 1: reference encoding
    else
        write_bit(w, false)  # main flag = 0: direct encoding
    end

    # if mode is :index, we need to write the out-degrees of the vertices first
    if mode == :index
        out_degrees = Vector{T}(undef, vs)
        for v in 1:vs
            out_degrees[v] = T(length(neighbor_lists[v]))
        end
        # shift the out-degrees by 1 to avoid zeros
        out_degrees = out_degrees .+ T(1)
        for v in 1:vs
            write_encoded_value(w, out_degrees[v], encoding)
        end
    end

    # write adjacency lists one after the other
    for v in 1:vs
        current_neighbors = neighbor_lists[v]
        
        if isempty(current_neighbors)
            # for empty lists in children mode, write stop value if not the last vertex
            if mode == :children && v < vs
                # write direct encoding flag
                if enable_reference
                    write_bit(w, false)  # flag = 0: direct encoding
                end
                # write stop value
                write_encoded_value(w, T(1), encoding)  # stop value
            end
            continue
        end
        
        # sort neighbors 
        sorted_neighbors = sort(copy(current_neighbors))
       
        ########################################################################
        # Index mode: write delta encoded values without shifting, no stop values
        ########################################################################
        if mode == :index
            # track if the vertex has been written
            vertex_written = false

            # reference encoding
            if enable_reference && length(sorted_neighbors) > 3
                # find the best reference
                best_ref = find_best_reference(sorted_neighbors, neighbor_lists, potential_references)
                
                if best_ref !== nothing
                    # write reference encoding flag
                    write_bit(w, true)   # flag = 1: reference encoding
                            
                    # write reference vertex ID
                    write_encoded_value(w, T(best_ref), encoding)
                            
                    # create copy bitmap and residuals
                    ref_neighbors = neighbor_lists[best_ref]
                    copy_bitmap, residuals = create_reference_data(sorted_neighbors, ref_neighbors)
                            
                    # write copy bitmap length and the bitmap itself
                    write_encoded_value(w, T(length(copy_bitmap)), encoding)
                    for bit in copy_bitmap
                        write_bit(w, bit)
                    end
                            
                    # write residuals flag and residuals if any exist
                    if !isempty(residuals)
                        write_bit(w, true)   # residuals flag = 1: has residuals

                        # write residuals
                        write_run_length_delta(w, encoding, residuals) 
                    else
                        # no residuals: write flag 0
                        write_bit(w, false)
                    end

                    # add current vertex to the potential references
                    push!(potential_references, T(v))

                    # set the vertex as written
                    vertex_written = true
                end
            # direct encoding
            elseif !vertex_written
                if enable_reference
                    # write direct encoding flag
                    write_bit(w, false)   # flag = 0: direct encoding
                end

                # write run-length + delta encoding for direct encoding
                write_run_length_delta(w, encoding, sorted_neighbors)

                # add current vertex to potential references
                push!(potential_references, T(v))
            end # end direct encoding
        ########################################################################
        # Children mode: write delta encoded values shifted by 1, with stop values
        ########################################################################
        else
            # track if the vertex has been written
            vertex_written = false

            # reference encoding
            if enable_reference && length(sorted_neighbors) > 3
                # find the best reference
                best_ref = find_best_reference(sorted_neighbors, neighbor_lists, potential_references)
                
                if best_ref !== nothing
                    # write reference encoding flag
                    write_bit(w, true)   # flag = 1: reference encoding
                            
                    # write reference vertex ID
                    write_encoded_value(w, T(best_ref), encoding)
                            
                    # create copy bitmap and residuals
                    ref_neighbors = neighbor_lists[best_ref]
                    copy_bitmap, residuals = create_reference_data(sorted_neighbors, ref_neighbors)
                            
                    # write copy bitmap length and the bitmap itself
                    write_encoded_value(w, T(length(copy_bitmap)), encoding)
                    for bit in copy_bitmap
                        write_bit(w, bit)
                    end
                            
                    # write residuals flag and residuals if any exist
                    if !isempty(residuals)
                        write_bit(w, true)   # residuals flag = 1: has residuals
                            
                        # write residuals
                        write_run_length_delta(w, encoding, residuals) 
                    else
                        # no residuals: write flag 0
                        write_bit(w, false)
                    end

                    # add current vertex to the potential references
                    push!(potential_references, T(v))

                    # set the vertex as written
                    vertex_written = true

                    # add stop value T(1) at the end of each list except the last one
                    if v < vs
                        # write direct encoding flag
                        # NB: reference encoding is enabled
                        write_bit(w, false)  # flag = 0: direct encoding
                        write_encoded_value(w, T(1), encoding)  # stop value
                    end
                end
            # direct encoding
            elseif !vertex_written
                if enable_reference
                    # write direct encoding flag
                    write_bit(w, false)   # flag = 0: direct encoding
                end
                
                # Use run-length + delta encoding for direct encoding
                write_run_length_delta(w, encoding, sorted_neighbors)

                # add stop value T(1) at the end of each list except the last one
                if v < vs
                    if enable_reference
                        write_bit(w, false)  # flag = 0: direct encoding
                    end
                    write_encoded_value(w, T(1), encoding)  # stop value
                end

                # add current vertex to potential references
                push!(potential_references, T(v))
            end # end direct encoding
        end
    end
end

"""
    read_reference_encoding(r::BitReader, encoding::Symbol, mode::Symbol=:children, enable_reference::Bool=true, ::Type{T}=UInt8) where {T<:Unsigned}

Read reference-encoded neighbor lists with run-length and delta compression from the bitreader.

# Modes:

## Index Mode (:index)
- First reads out-degrees using run-length+delta decoding
- Then reads adjacency lists using delta decoding (no value shifting, no stop values)
- If enable_reference=true: handles reference encoding when present

## Children Mode (:children) 
- Reads adjacency lists using delta decoding with values shifted by 1 (minimum value 1)
- Handles stop value T(1) at end of each list except the last
- If enable_reference=true: handles reference encoding when present

# Decoding Format:

## Global Header:
```
[Global Main Flag: 1 bit] [mode-specific data...]
```

- **Global Main Flag = 1**: Reference encoding is enabled
- **Global Main Flag = 0**: Reference encoding is disabled (direct encoding only)

## When Reference Encoding is Enabled (Global Main Flag = 1):

For each vertex, reads:
```
[Vertex Flag: 1 bit] [Encoding-specific data...]
```

### Vertex Flag Values:
- **0**: Direct delta encoding for this vertex
- **1**: Reference encoding for this vertex

### Format for Vertex Flag = 0 (Direct Encoding):
- **Index mode**: delta encoded values without shifts or stops
- **Children mode**: delta encoded values shifted by 1 with stop values

### Format for Vertex Flag = 1 (Reference Encoding):
```
[1] [ref_id: varint] [bitmap_length: varint] [bitmap: N bits] [residuals_flag: 1 bit] [residuals?]
```

Where:
- **ref_id**: ID of reference vertex (decoded with specified compression)
- **bitmap_length**: Number of bits in copy bitmap (decoded with specified compression)
- **bitmap**: Raw bits indicating which reference elements to copy (1=copy, 0=skip)
- **residuals_flag**: 1 bit (1=has residuals, 0=no residuals)
- **residuals**: If residuals_flag=1, delta decoded residual elements (following mode-specific decoding)

## When Reference Encoding is Disabled (Global Main Flag = 0):

All vertices use direct delta encoding without any per-vertex flags:
- **Index mode**: delta encoded values without shifts or stops
- **Children mode**: delta encoded values shifted by 1 with stop values

@param r::BitReader: the bitreader
@param vs::T: the number of vertices in the graph
@param encoding::Symbol: the compression coding used for values (:elias_gamma, :elias_delta, :golomb, :fibonacci, :zeta)
@param mode::Symbol: the mode to use (:children, :index)
@param enable_reference::Bool: whether to enable reference encoding (default: true)
@param T::Type: the type to return (default: UInt8)
@return::Dict{T,Vector{T}}: the decoded neighbor lists
"""
function read_reference_encoding(r::BitReader, vs::T, encoding::Symbol, mode::Symbol=:children, enable_reference::Bool=true, ::Type{T}=UInt8) where {T<:Unsigned}
    # initialize the neighbor lists
    neighbor_lists = Dict{T,Vector{T}}() 

    # read the global main flag to determine if reference encoding is enabled
    global_main_flag = read_bit(r)
    # global_main_flag = true means reference encoding is enabled
    # global_main_flag = false means direct encoding only
    
    # if index mode, read out-degrees
    out_degrees = T[]
    if mode == :index
        out_degrees = Vector{T}(undef, vs)
        for v in 1:vs
            out_degrees[v] = read_encoded_value(r, encoding, T)
        end
        # unshift the out-degrees (they were shifted by 1 during encoding)
        out_degrees = out_degrees .- T(1)
    end
    
    vertex_id = T(1)  # track current vertex ID
    
    try
        while vertex_id <= vs
            current_neighbors = T[]
            
            if global_main_flag && enable_reference
                # Reference encoding is enabled: read vertex flag for each vertex
                vertex_flag = read_bit(r)
                
                if vertex_flag
                    # vertex flag = 1: reference encoding for this vertex
                    
                    # read reference vertex id
                    ref_id = T(read_encoded_value(r, encoding, T))
                    
                    # read copy bitmap
                    bitmap_length = Int(read_encoded_value(r, encoding, T))
                    copy_bitmap = Bool[]
                    for _ in 1:bitmap_length
                        push!(copy_bitmap, read_bit(r))
                    end
                    
                    # read residuals flag
                    residuals_flag = read_bit(r)
                    residuals = if residuals_flag
                        # residuals flag = 1: has residuals
                        read_run_length_delta(r, encoding, T)
                    else
                        # residuals flag = 0: no residuals
                        T[]
                    end
                    
                    # reconstruct neighbor list from reference
                    if haskey(neighbor_lists, ref_id)
                        ref_neighbors = neighbor_lists[ref_id]
                        current_neighbors = reconstruct_from_reference(ref_neighbors, copy_bitmap, residuals)
                    else
                        throw(ArgumentError("Invalid reference ID: $ref_id (not found in already decoded vertices)"))
                    end
                else
                    # vertex flag = 0: direct encoding for this vertex
                    degrees_param = mode == :index ? out_degrees : T[]
                    current_neighbors = read_direct_encoded_list(r, encoding, mode, vertex_id, degrees_param, T)
                end
            else
                # Reference encoding is disabled: direct encoding only (no vertex flags)
                degrees_param = mode == :index ? out_degrees : T[]
                current_neighbors = read_direct_encoded_list(r, encoding, mode, vertex_id, degrees_param, T)
            end
            
            # stop values are handled within read_delta_encoded_list for children mode
            
            neighbor_lists[vertex_id] = current_neighbors
            vertex_id += T(1)
            
            # for children mode, check if we've reached the end
            if mode == :children && vertex_id > vs
                break
            end
        end
    catch e
        # end of stream reached - this is expected behavior for children mode
        if !isa(e, EOFError) && !isa(e, ErrorException)
            rethrow(e)
        end
        # for children mode, this is normal termination
        if mode == :children
            # normal end of stream
        else
            # for index mode, we should have read exactly vs vertices
            if vertex_id <= vs
                @warn "Unexpected end of stream in index mode at vertex $vertex_id (expected $vs)"
            end
        end
    end
    
    return neighbor_lists
end

"""
    read_direct_encoded_list(r::BitReader, encoding::Symbol, mode::Symbol, vertex_id::T, out_degrees::Vector{T}, ::Type{T}=UInt8) where {T<:Unsigned}

Read a directly-encoded neighbor list according to mode-specific format.

@param r::BitReader: the bitreader
@param encoding::Symbol: the compression coding used
@param mode::Symbol: the mode (:index or :children)
@param vertex_id::T: current vertex ID being processed
@param out_degrees::Vector{T}: the out-degrees for each vertex (for index mode), empty for children mode
@param T::Type: the type to return
@return::Vector{T}: the decoded neighbor list
"""
function read_direct_encoded_list(r::BitReader, encoding::Symbol, mode::Symbol, vertex_id::T, out_degrees::Vector{T}, ::Type{T}=UInt8) where {T<:Unsigned}
    try
        if mode == :children
            # For children mode, read run-length+delta stream with stop value detection
            # The stop value T(1) is written separately after the run-length+delta data
            # We can use the stop_value parameter in read_run_length_delta
            neighbors = read_run_length_delta(r, encoding, T; stop_value=T(1))
            return neighbors
        else
            # Index mode: read exactly out_degrees[vertex_id] neighbors
            expected_neighbors = out_degrees[vertex_id]
            
            if expected_neighbors == 0
                return T[]  # Empty list
            end
            
            # Use max_elements to limit reading to expected number of neighbors
            neighbors = read_run_length_delta(r, encoding, T; max_elements=Int(expected_neighbors))
            return neighbors
        end
    catch e
        # Empty list case: if we can't even read the first value
        if isa(e, EOFError) || isa(e, ErrorException)
            return T[]
        else
            rethrow(e)
        end
    end
end

"""
    read_delta_encoded_list(r::BitReader, encoding::Symbol, mode::Symbol, ::Type{T}=UInt8) where {T<:Unsigned}

Read a delta-encoded list according to mode-specific format.

@param r::BitReader: the bitreader
@param encoding::Symbol: the compression coding used
@param mode::Symbol: the mode (:index or :children)
@param T::Type: the type to return
@return::Vector{T}: the decoded neighbor list
"""
function read_delta_encoded_list(r::BitReader, encoding::Symbol, mode::Symbol, ::Type{T}=UInt8) where {T<:Unsigned}
    neighbors = T[]
    
    try
        # Read until we can't read anymore (for single list) or hit stop value (children mode)
        first_val = read_encoded_value(r, encoding, T)
        push!(neighbors, first_val)
        
        # Read subsequent delta values
        while true
            try
                val = read_encoded_value(r, encoding, T)
                
                if mode == :children
                    # In children mode, check for stop value T(1)
                    if val == T(1)
                        # This is the stop value, we're done with this list
                        break
                    end
                    # Unshift the value (it was shifted by 1 during encoding)
                    delta_val = val - T(1)
                else
                    # Index mode: no shifting
                    delta_val = val
                end
                
                # Reconstruct the actual value from delta
                next_val = neighbors[end] + delta_val
                push!(neighbors, next_val)
            catch e
                # End of stream or list
                break
            end
        end
    catch e
        # Empty list case
        if !isa(e, EOFError) && !isa(e, ErrorException)
            rethrow(e)
        end
    end
    
    return neighbors
end

"""
    find_best_reference(target::Vector{T}, all_lists::Dict{T,Vector{T}}, available::Set{T}) where {T<:Unsigned}

Find the best reference vertex for the target neighbor list.

@param target::Vector{T}: the neighbor list to encode
@param all_lists::Vector{Vector{T}}: all neighbor lists
@param available::Set{T}: vertices that can serve as references (already encoded)
@return::Union{T, Nothing}: the best reference vertex ID or nothing if no good reference
"""
function find_best_reference(target::Vector{T}, all_lists::Dict{T,Vector{T}}, available::Set{T}) where {T<:Unsigned}
    if isempty(available)
        return nothing
    end
    
    best_ref = nothing
    best_savings = 0
    
    for ref_id in available
        ref_list = all_lists[ref_id]
        
        # calculate potential savings
        savings = calculate_reference_savings(target, ref_list)
        
        # require at least 50% overlap for reference to be worthwhile
        if savings > best_savings && savings > length(target) * 0.5
            best_savings = savings
            best_ref = ref_id
        end
    end
    
    return best_ref
end

"""
    calculate_reference_savings(target::Vector{T}, reference::Vector{T}) where {T<:Unsigned}

Calculate the compression savings from using reference encoding.

@param target::Vector{T}: the target neighbor list
@param reference::Vector{T}: the potential reference list
@return::Int: estimated bit savings (higher is better)
"""
function calculate_reference_savings(target::Vector{T}, reference::Vector{T}) where {T<:Unsigned}
    if isempty(reference)
        return 0
    end
    
    # count shared elements
    target_set = Set(target)
    reference_set = Set(reference)
    shared_count = length(intersect(target_set, reference_set))
    
    # estimate savings: shared elements cost 1 bit each in bitmap vs full encoding
    # rough heuristic: each shared element saves ~log2(max_vertex_id) - 1 bits
    avg_savings_per_shared = 10  # heuristic for typical graphs
    
    return shared_count * avg_savings_per_shared
end

"""
    create_reference_data(target::Vector{T}, reference::Vector{T}) where {T<:Unsigned}

Create copy bitmap and residuals for reference encoding.

@param target::Vector{T}: the target neighbor list
@param reference::Vector{T}: the reference neighbor list  
@return::Tuple{Vector{Bool}, Vector{T}}: (copy_bitmap, residuals)
"""
function create_reference_data(target::Vector{T}, reference::Vector{T}) where {T<:Unsigned}
    target_set = Set(target)
    copy_bitmap = Bool[]
    residuals = T[]
    
    # create copy bitmap for reference elements
    for ref_neighbor in reference
        if ref_neighbor in target_set
            push!(copy_bitmap, true)
        else
            push!(copy_bitmap, false)
        end
    end
    
    # find residuals (elements in target but not in reference)
    reference_set = Set(reference)
    for target_neighbor in target
        if !(target_neighbor in reference_set)
            push!(residuals, target_neighbor)
        end
    end
    
    # sort residuals for delta encoding
    sort!(residuals)
    
    return copy_bitmap, residuals
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
    result = T[]
    
    # copy selected elements from reference
    for (i, should_copy) in enumerate(copy_bitmap)
        if should_copy && i <= length(reference)
            push!(result, reference[i])
        end
    end
    
    # add residuals
    append!(result, residuals)
    
    # sort the final result
    sort!(result)
    
    return result
end

end # module Compression
