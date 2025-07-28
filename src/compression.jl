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
using ..IO: BitWriter, BitReader, write_bit, write_bits, read_bit, read_bits, flush, write_value, read_value
using ..Constants: FIB_NUMBERS, BUFFER_SIZE

# Export the functions we want to make available
export write_unary_coding,
       write_truncated_binary_coding,
       huffman_encoding,
       encode_huffman_tree!,
       decode_huffman_tree!,
       get_huffman_codes!,
       decode_huffman_values,
       delta_encode_vector,
       write_elias_gamma,
       read_elias_gamma,
       write_elias_coding,
       read_elias_coding,
       write_golomb,
       read_golomb,
       write_fibonacci_code,
       read_fibonacci_code

################################################################################
# Basic encoding / decoding
################################################################################

"""
    write_unary_coding(w::BitWriter, v::T) where {T<:Unsigned}

Write a unary coding to the bitwriter.

@param w::BitWriter: the bitwriter
@param v::T: the value to write
"""
function write_unary_coding(w::BitWriter, v::T) where {T<:Unsigned}
    for i in 1:v
        write_bit(w, true)
    end
    write_bit(w, false)
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
@return::Vector{T}: The delta encoded list
"""
function delta_encode_vector(lst::Vector{T}, shifted::Bool = false)::Vector{T} where {T<:Unsigned}
    if shifted
        # shift the list by 1 to avoid 0
        lst_shifted = lst .+ 1
    else
        lst_shifted = lst
    end
    # if the list is empty, return an empty list
    isempty(lst_shifted) && return T[]
    # initialize the differences with the first element
    diffs = [T(lst_shifted[firstindex(lst_shifted)])]
    # for each element in the list, compute the difference with the previous element
    for i in eachindex(lst_shifted)[2:end]
        push!(diffs, T(lst_shifted[i] - lst_shifted[i-1]))
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
    len = read_elias_gamma(r, T)

    len < 1 && throw(ArgumentError("Invalid Elias delta encoding: length must be >= 1"))

    # Step 2: Read the (len - 1) lower bits
    tail = read_value(r, len - 1, T)

    # Step 3: Combine with leading 1
    return (T(1) << (len - 1)) | tail
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
    write_fibonacci_code(w::BitWriter, n::T) where {T<:Unsigned}

Write `n` using Fibonacci coding (with '1' stop bit), using precomputed FIB_NUMBERS.
"""
function write_fibonacci_code(w::BitWriter, n::T) where {T<:Unsigned}
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
    read_fibonacci_code(r::BitReader, ::Type{T}=UInt8) where {T<:Unsigned}

Decode a Fibonacci code using precomputed FIB_NUMBERS.
"""
function read_fibonacci_code(r::BitReader, ::Type{T}=UInt8) where {T<:Unsigned}
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

end # module Compression
