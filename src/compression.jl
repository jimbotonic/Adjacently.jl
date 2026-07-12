#
# Adjacently: Julia Complex Directed Networks Library
# Copyright (C) 2016-2026 Anonymous (double-blind review)
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
flush, write_value, read_value, flush_bitwriter, write_bytes, reset_bitwriter!, bytes_written, get_bytes

using ..Constants: FIB_NUMBERS, BUFFER_SIZE, ZETA_H_BOUNDS, ZETA_POWER_BASES, ZETA_BASE,
GOLOMB_BASE, REF_ENCODING_TH, REF_V_MIN_DEGREE, FED_BLOCK_SIZE, MIN_INTERVAL_LENGTH, MIN_RUN_LENGTH, REF_WINDOW_SIZE

using ..Index
using ..CompressionUtils

# Streaming context-adaptive range coder for the :context_range integer backend.
include("range_coder.jl")



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
       write_hybrid_mix_encoded_list,
       read_hybrid_mix_encoded_list,
       analyze_delta_patterns_hybrid,
       reconstruct_from_delta,
       read_elias_coding,
       write_elias_fano,
       read_elias_fano,
       write_golomb,
       read_golomb,
       write_golomb_rice,
       read_golomb_rice,
       write_fibonacci,
       read_fibonacci,
       write_zeta,
       read_zeta,
       write_fed,
       read_fed,
       compress_intervals,
       write_intervals_and_residuals,
       read_intervals_and_residuals,
       write_compressed_graph_data,
       read_compressed_graph_data,
       write_bitpacked_bitmap,
       write_rle_ones_deltas,
       write_bitmap_rle_ones,
       read_bitmap_rle_ones,
       write_bitmap_adaptive,
       read_bitmap_adaptive,
       write_block_encoding,
       read_block_encoding,
       estimate_block_encoding_cost,
       write_small_count,
       estimate_encoded_value_cost,
       estimate_interval_runlength_encoding_cost,
       write_greedy_graph_data,
       read_greedy_graph_data,
       parse_bg_ctxrange_index,
       write_cmdstream_graph_data,
       read_cmdstream_graph_data,
       CommandStream,
       COST_MODEL_FULL,
       COST_MODEL_FAST,
       DEFAULT_COST_MODEL

# Lightweight workspace types for reference building
struct RefBuildWorkspace{T<:Unsigned}
    copy_bitmap::Vector{Bool}
    residuals::Vector{T}
end

RefBuildWorkspace{T}() where {T<:Unsigned} = RefBuildWorkspace{T}(Bool[], T[])

# Pre-allocated buffers for exact bit-counting in BG greedy reference search.
# Uses buffer-only BitWriters (no IO target) for fast reset.
mutable struct GreedyCostBuffer
    w1::BitWriter   # copy encoding trial
    w2::BitWriter   # stop-delta / residual trial
    w3::BitWriter   # intervals trial
end

function GreedyCostBuffer()
    GreedyCostBuffer(BitWriter(capacity=1024), BitWriter(capacity=1024), BitWriter(capacity=1024))
end

@inline function _greedy_reset!(w::BitWriter)
    reset_bitwriter!(w)
end

# Pre-allocated workspace for BG/CS vertex search.
# Eliminates per-vertex allocation of overlap_candidates, merge buffers, and snapshots.
mutable struct VertexSearchWorkspace{T<:Unsigned}
    overlap_candidates::Vector{Tuple{Int,Vector{T},T}}
    merge_bmp::Vector{Bool}
    merge_res::Vector{T}
    best_bmp::Vector{Bool}   # only written on improvement
    best_res::Vector{T}      # only written on improvement
end

function VertexSearchWorkspace{T}() where {T<:Unsigned}
    VertexSearchWorkspace{T}(
        Vector{Tuple{Int,Vector{T},T}}(),
        Bool[], T[],
        Bool[], T[]
    )
end

# Resize dst to match src length and copy contents.
@inline function _copyto_resize!(dst::Vector{T}, src::Vector{T}) where T
    resize!(dst, length(src))
    copyto!(dst, src)
end
@inline function _copyto_resize!(dst::Vector{Bool}, src::Vector{Bool})
    resize!(dst, length(src))
    copyto!(dst, src)
end

# Two-pointer overlap count on sorted vectors. O(|a|+|b|), zero allocation.
@inline function _sorted_overlap_count(a::Vector{T}, b::Vector{T})::Int where {T}
    count = 0
    i = 1; j = 1
    @inbounds while i <= length(a) && j <= length(b)
        if a[i] == b[j]
            count += 1; i += 1; j += 1
        elseif a[i] < b[j]
            i += 1
        else
            j += 1
        end
    end
    return count
end

# Single-pass sorted merge producing copy bitmap + residual vector.
# `bitmap` is resized to length(ref_nl), `residuals` is emptied and filled.
function _sorted_merge_ref!(bitmap::Vector{Bool}, residuals::Vector{T},
                            neighbors::Vector{T}, ref_nl::Vector{T}) where {T}
    ref_len = length(ref_nl)
    resize!(bitmap, ref_len)
    empty!(residuals)
    i = 1; j = 1
    @inbounds while i <= length(neighbors) && j <= ref_len
        if neighbors[i] == ref_nl[j]
            bitmap[j] = true; i += 1; j += 1
        elseif neighbors[i] < ref_nl[j]
            push!(residuals, neighbors[i]); i += 1
        else
            bitmap[j] = false; j += 1
        end
    end
    @inbounds while j <= ref_len; bitmap[j] = false; j += 1; end
    @inbounds while i <= length(neighbors); push!(residuals, neighbors[i]); i += 1; end
    return nothing
end

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
# Measurement tap: when set to a Vector, records every coded integer so we can
# compare alternative integer backends (e.g. context-range) offline. nothing = off.
# Concretely typed so the hot-path check in write_encoded_value is type-stable and
# allocation-free (no boxed Any load per coded integer).
const _VALUE_TAP = Ref{Union{Nothing,Vector{UInt64}}}(nothing)
# Structural taps (lever-2 sizing): reference distances + copy-bitmap bits.
const _REFDIST_TAP = Ref{Any}(nothing)
const _COPYBMP_TAP = Ref{Any}(nothing)
# Exact per-stream bit accounting (Dict{Symbol,Int} when active) + residual accumulator.
const _STREAM_BITS = Ref{Any}(nothing)
const _RESID_ACC = Ref{Any}(nothing)
# Per-vertex residual metadata for typed-context experiments:
# Vector of (degree, n_residual_values, ref_mode_code) in encode order. ref_mode:
# 0=none, 1=reference, 2=multi_ref.
const _RESID_META = Ref{Any}(nothing)
# Per-vertex residual neighbor-id list (diagnostic tap for copy-aware rank-gap
# experiments): when set to a Vector, _encode_vertex! pushes a copy of each
# vertex's final residual (non-copied) neighbor id list in encode order.
const _RESIDIDS_TAP = Ref{Any}(nothing)
# Production :context_range backend: when set to a CtxRangeEncoder / CtxRangeDecoder,
# write_encoded_value / read_encoded_value redirect the residual-region integer
# stream to the range coder instead of the bitstream. Gated to the residual region
# by _encode_vertex! / read_greedy_graph_data.
# Concretely typed (CtxRangeEncoder/Decoder are defined in range_coder.jl, include'd
# above at L33) so the per-call redirect check in write/read_encoded_value is
# type-stable and allocation-free — the boxed Ref{Any} load on this hot path was
# the root cause of the flaky GC segfault.
const _RESID_SINK = Ref{Union{Nothing,CtxRangeEncoder}}(nothing)
const _RESID_SOURCE = Ref{Union{Nothing,CtxRangeDecoder}}(nothing)

# Phase 2b: structural FLAG bits (stop-delta continuation bits, hybrid-RLE section
# flags) redirect to an adaptive binary range coder when set — same gating pattern
# as _RESID_SINK. Set by write/read_cmdstream_graph_data (children mode, ctx) around
# the vertex loop; nothing = legacy raw bits in the struct bitstream.
const _FLAG_SINK = Ref{Union{Nothing,BinRangeEncoder}}(nothing)
const _FLAG_SOURCE = Ref{Union{Nothing,BinRangeDecoder}}(nothing)

@inline function _struct_write_bit(w, b::Bool)
    s = _FLAG_SINK[]
    if s === nothing
        write_bit(w, b)
    else
        brc_encode_bit!(s::BinRangeEncoder, b)
    end
    return nothing
end

@inline function _struct_read_bit(r)::Bool
    s = _FLAG_SOURCE[]
    return s === nothing ? read_bit(r) : brc_decode_bit!(s::BinRangeDecoder)
end

# --- Copy-aware rank gaps (Task A) -------------------------------------------
# In the :context_range residual path, residual neighbour ids are recoded as
# ranks in the complement of the already-copied set C (target ids copied from
# the reference), instead of absolute-coordinate gaps. This lets consecutive
# residuals whose gap is "explained" by an intervening copied id become adjacent
# ranks (interval-detectable), shaving residual bits. The transform is applied
# to the sorted residual id list R -> Q immediately BEFORE it enters the residual
# range coder, and inverted Q -> R immediately AFTER the residual range decoder;
# the coder / hybrid tokens / stop / count / interval logic are untouched
# (intervals are simply detected on Q). Tied to :context_range being active — no
# header bit. The transform is universe-agnostic: it depends only on C (sorted,
# unique, disjoint from R), so BG/CS (id universe 1:n) and CG (cluster-local
# 1:s) share it.
#
# rank_C(r) = r - #{c in C : c < r}   (1-based position of r in complement of C)
#
# Guards (caller-guaranteed): C sorted-unique, R sorted-strictly-increasing,
# R ∩ C = ∅. Then Q is sorted-strictly-increasing and 1-based.
@inline function _rankgap_forward(R::AbstractVector{I}, C::AbstractVector{I}) where {I<:Integer}
    Q = Vector{I}(undef, length(R))
    ci = 1; nc = length(C); copied_before = 0
    @inbounds for k in eachindex(R)
        r = R[k]
        while ci <= nc && C[ci] < r
            copied_before += 1; ci += 1
        end
        Q[k] = r - I(copied_before)
    end
    return Q
end

@inline function _rankgap_inverse(Q::AbstractVector{I}, C::AbstractVector{I}) where {I<:Integer}
    R = Vector{I}(undef, length(Q))
    ci = 1; nc = length(C)
    @inbounds for k in eachindex(Q)
        r = Q[k] + I(ci - 1)
        while ci <= nc && C[ci] <= r
            ci += 1; r += one(I)
        end
        R[k] = r
    end
    return R
end

# Map a single coordinate x (e.g. the vertex id) into the same complement-rank
# space as _rankgap_forward, so the residual coder's vertex-relative first-value
# zigzag keeps its locality after the residuals are remapped. Encode and decode
# must call this with the identical C to stay lossless.
@inline function _rank_of(x::I, C::AbstractVector{I}) where {I<:Integer}
    cnt = 0
    @inbounds for c in C
        c < x ? (cnt += 1) : break
    end
    return x - I(cnt)
end

function write_encoded_value(w::BitWriter, value::T, compression::Symbol) where {T<:Unsigned}
    let t = _VALUE_TAP[]; t === nothing || push!(t::Vector{UInt64}, UInt64(value)); end
    let s = _RESID_SINK[]
        if s !== nothing
            rc_encode_value!(s::CtxRangeEncoder, UInt64(value))
            return
        end
    end
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
    let s = _RESID_SOURCE[]
        if s !== nothing
            return T(rc_decode_value!(s::CtxRangeDecoder))
        end
    end
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

    # bits_len = number of significant bits in v
    bits_len = sizeof(T) * 8 - leading_zeros(v)

    # Gamma code = (bits_len-1) zero bits + bits_len bits of v
    # Total = 2*bits_len - 1 bits. Since v < 2^bits_len, the leading
    # zeros come from writing v in a wider field.
    write_value(w, UInt64(v), 2 * bits_len - 1)
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

"""
    write_golomb_rice(w::BitWriter, n::T, k::Int) where {T<:Unsigned}

Write `n` using Golomb-Rice encoding with parameter `k`.

Golomb-Rice is a special case of Golomb coding where b = 2^k.
This makes encoding/decoding more efficient since:
- The remainder is exactly k bits (no log2 calculation needed)
- Division by 2^k is a simple right shift
- Modulo 2^k is a simple bitwise AND

Encoding steps:
1. Quotient q = n >> k (divide by 2^k), encoded in unary: q zeros followed by one `1`
2. Remainder r = n & ((1 << k) - 1), encoded in binary with exactly k bits

@param w::BitWriter: the bitwriter
@param n::T: the value to encode (must be non-negative)
@param k::Int: the Rice parameter (b = 2^k)
"""
function write_golomb_rice(w::BitWriter, n::T, k::Int) where {T<:Unsigned}
    k < 0 && throw(ArgumentError("Rice parameter k must be >= 0"))

    # Step 1: Unary representation of quotient (n >> k)
    q = n >> k
    for _ in 1:q
        write_bit(w, false)
    end
    write_bit(w, true)  # stop bit

    # Step 2: Binary representation of remainder (n & ((1 << k) - 1))
    if k > 0
        r = n & ((T(1) << k) - T(1))
        write_value(w, r, k)
    end
end

"""
    read_golomb_rice(r::BitReader, k::Int, ::Type{T}=UInt8) where {T<:Unsigned}

Read a Golomb-Rice encoded value from the bit reader.

Golomb-Rice code with parameter k:
- Unary-encoded quotient `q`: `q` zeros followed by a `1`
- Binary-encoded remainder `r` in exactly `k` bits
- Returns `(q << k) + r` which equals `q * 2^k + r`

@param r::BitReader: the bitreader
@param k::Int: the Rice parameter (b = 2^k)
@param T::Type: the type to return
@return::T: the decoded value
"""
function read_golomb_rice(r::BitReader, k::Int, ::Type{T}=UInt8) where {T<:Unsigned}
    k < 0 && throw(ArgumentError("Rice parameter k must be >= 0"))

    # Step 1: Read unary prefix to get quotient
    q = T(0)
    while !read_bit(r)
        q += T(1)
    end

    # Step 2: Read k-bit binary remainder
    if k > 0
        remainder = read_value(r, k, T)
        return (q << k) + remainder
    else
        return q
    end
end

################################################################################
# Fibonacci code
################################################################################

# Precomputed Fibonacci encode table: for values 1..FIB_LUT_SIZE, each entry
# packs (msb_code << 6) | code_len into a UInt64. Lookup replaces greedy
# decomposition + bit reversal + write_value loop with a single table lookup.
const FIB_LUT_SIZE = 8192

const _FIB_ENCODE_TABLE = let
    table = Vector{UInt64}(undef, FIB_LUT_SIZE)
    for n in 1:FIB_LUT_SIZE
        code = UInt64(0); max_idx = 0; remaining = UInt64(n)
        for i in length(FIB_NUMBERS):-1:1
            if FIB_NUMBERS[i] <= remaining
                remaining -= FIB_NUMBERS[i]
                code |= (UInt64(1) << (i - 1))
                max_idx = max(max_idx, i)
                remaining == 0 && break
            end
        end
        code_len = max_idx + 1
        msb_code = UInt64(0)
        for i in 0:(max_idx - 1)
            msb_code = (msb_code << 1) | ((code >> i) & 1)
        end
        msb_code = (msb_code << 1) | UInt64(1)  # stop bit
        table[n] = (msb_code << 6) | UInt64(code_len)
    end
    table
end

"""
    write_fibonacci(w::BitWriter, n::T) where {T<:Unsigned}

Write `n` using Fibonacci coding (with '1' stop bit).
Uses a precomputed lookup table for values ≤ $FIB_LUT_SIZE.
"""
@inline function write_fibonacci(w::BitWriter, n::T) where {T<:Unsigned}
    n == 0 && throw(ArgumentError("Fibonacci code undefined for 0"))
    if n <= FIB_LUT_SIZE
        packed = @inbounds _FIB_ENCODE_TABLE[Int(n)]
        write_value(w, packed >> 6, Int(packed & 0x3f))
    else
        _write_fibonacci_fallback(w, n)
    end
end

# Fallback for values > FIB_LUT_SIZE (rare in practice)
function _write_fibonacci_fallback(w::BitWriter, n::T) where {T<:Unsigned}
    code = UInt64(0)
    max_idx = 0
    remaining = UInt64(n)
    @inbounds for i in length(FIB_NUMBERS):-1:1
        if FIB_NUMBERS[i] <= remaining
            remaining -= FIB_NUMBERS[i]
            code |= (UInt64(1) << (i - 1))
            if i > max_idx; max_idx = i; end
            remaining == 0 && break
        end
    end

    code_len = max_idx + 1
    msb_code = UInt64(0)
    @inbounds for i in 0:(max_idx - 1)
        msb_code = (msb_code << 1) | ((code >> i) & 1)
    end
    msb_code = (msb_code << 1) | UInt64(1)

    write_value(w, msb_code, code_len)
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
            value += T(FIB_NUMBERS[i])
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


"""
    write_intervals_and_residuals(w::BitWriter, neighbors::Vector{T}, encoding::Symbol, min_interval_length::Int=4) where {T<:Unsigned}

Write neighbor list using interval + run-length + residual encoding.
This is an enhanced WebGraph-style encoding that:
1. Extracts intervals (consecutive sequences)
2. Analyzes remaining values for run-length patterns (repeated gaps)
3. Encodes final residuals with delta encoding

@param w::BitWriter: the bitwriter
@param neighbors::Vector{T}: sorted list of neighbors
@param encoding::Symbol: encoding for residuals
@param min_interval_length::Int: minimum interval length
@param min_run_length::Int: minimum repetitions for run-length encoding
"""
function write_intervals_and_residuals(w::BitWriter, neighbors::Vector{T}, encoding::Symbol, min_interval_length::Int=MIN_INTERVAL_LENGTH; vertex_id=nothing, tight_intervals::Bool=false) where {T<:Unsigned}
    # Handle empty neighbor lists by writing zero counts
    if isempty(neighbors)
        @debug "WRITE intervals: empty neighbor list, writing two 1s"
        # Write num_intervals = 0 (encoded as 0+1=1)
        write_encoded_value(w, T(1), encoding)
        @debug "WRITE intervals: wrote num_intervals=0 (encoded as 1)"
        # Write num_residuals = 0 (encoded as 0+1=1)
        write_encoded_value(w, T(1), encoding)
        @debug "WRITE intervals: wrote num_residuals=0 (encoded as 1)"
        return
    end

    # Extract intervals from neighbors
    intervals, residuals = compress_intervals(neighbors, min_interval_length)

    # Write format: intervals, then residuals

    # Write number of intervals (increment by 1 to avoid 0)
    write_encoded_value(w, T(length(intervals)) + T(1), encoding)

    # Write intervals: (start, length) pairs with delta encoding
    if !isempty(intervals)
        prev_ref = T(0)  # reference point for delta encoding
        for (idx, (start, length)) in enumerate(intervals)
            if idx == 1 && vertex_id !== nothing
                # First interval start: zigzag offset from vertex_id
                offset = Int64(start) - Int64(vertex_id)
                encoded_start = UInt64(_zigzag_encode(offset) + 1)
                write_encoded_value(w, encoded_start, encoding)
            else
                # Delta encode start positions
                write_encoded_value(w, start - prev_ref, encoding)
            end
            # Encode interval length (already >= min_interval_length, increment by 1 to avoid 0)
            write_encoded_value(w, length - T(min_interval_length) + T(1), encoding)
            # tight_intervals: delta from end of interval; default: delta from start
            prev_ref = tight_intervals ? (start + length) : start
        end
    end

    # Write number of residuals (increment by 1 to avoid 0)
    write_encoded_value(w, T(length(residuals)) + T(1), encoding)

    # Write residuals with delta encoding
    if !isempty(residuals)
        @debug "WRITE residuals: residuals=$residuals"
        write_delta(w, residuals, encoding; vertex_id=vertex_id)
    end
end

"""
    read_intervals_and_residuals(r::BitReader, encoding::Symbol, min_interval_length::Int=4, min_run_length::Int=3, ::Type{T}=UInt8) where {T<:Unsigned}

Read neighbor list from interval + run-length + residual encoding.
Reads the format: intervals, then run-length pairs, then final residuals.

@param r::BitReader: the bitreader
@param encoding::Symbol: encoding used for residuals
@param min_interval_length::Int: minimum interval length used
@param min_run_length::Int: minimum run-length count used
@param T::Type: the type to return
@return::Vector{T}: reconstructed neighbor list
"""
function read_intervals_and_residuals(r::BitReader, encoding::Symbol, min_interval_length::Int=MIN_INTERVAL_LENGTH, ::Type{T}=UInt8; vertex_id=nothing, tight_intervals::Bool=false) where {T<:Unsigned}
    # Read number of intervals (subtract 1)
    num_intervals_raw = read_encoded_value(r, encoding, T)
    num_intervals = num_intervals_raw - T(1)
    @debug "READ intervals: num_intervals_raw=$num_intervals_raw, num_intervals=$num_intervals"

    neighbors = T[]
    @debug "READ intervals: starting with empty neighbors"

    # Read intervals
    if num_intervals > 0
        prev_ref = T(0)  # reference point for delta decoding
        for idx in 1:num_intervals
            if idx == 1 && vertex_id !== nothing
                # First interval start: zigzag decode from vertex_id
                raw_start = read_encoded_value(r, encoding, UInt64)
                start = T(Int64(vertex_id) + _zigzag_decode(raw_start - 1))
            else
                # Read delta-encoded start
                start_delta = read_encoded_value(r, encoding, T)
                start = prev_ref + start_delta
            end
            # Read length (subtract 1 and add back min_interval_length)
            length = read_encoded_value(r, encoding, T) - T(1) + T(min_interval_length)

            # Reconstruct interval
            for j in 0:(Int(length)-1)
                push!(neighbors, start + T(j))
            end
            # tight_intervals: delta from end of interval; default: delta from start
            prev_ref = tight_intervals ? (start + length) : start
        end
        @debug "READ intervals: after reading intervals, neighbors=$neighbors"
    end

    # Read number of residuals (subtract 1)
    num_residuals_raw = read_encoded_value(r, encoding, T)
    num_residuals = num_residuals_raw - T(1)
    @debug "READ residuals: num_residuals_raw=$num_residuals_raw, num_residuals=$num_residuals"

    # Read residuals
    if num_residuals > 0
        residuals = read_delta(r, encoding, T; max_elements=Int(num_residuals), vertex_id=vertex_id)
        @debug "READ residuals: got residuals=$residuals"
        append!(neighbors, residuals)
        @debug "READ residuals: neighbors after append=$neighbors"
    end

    # Sort the combined result
    @debug "READ: before sort, neighbors=$neighbors"
    sort!(neighbors)
    @debug "READ: after sort, neighbors=$neighbors"
    @debug "READ: returning neighbors=$neighbors"
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
function write_delta(w::BitWriter, lst::Vector{T}, encoding::Symbol; vertex_id=nothing, positive_gaps::Bool=false) where {T<:Unsigned}
    isempty(lst) && return

    # Write the first value (absolute or zigzag-relative)
    if vertex_id !== nothing
        offset = Int64(lst[1]) - Int64(vertex_id)
        encoded_first = UInt64(_zigzag_encode(offset) + 1)
        write_encoded_value(w, encoded_first, encoding)
    else
        write_encoded_value(w, lst[1], encoding)
    end

    # Write remaining values as inline gaps (no delta vector allocation)
    if positive_gaps
        @inbounds for i in 2:length(lst)
            write_encoded_value(w, lst[i] - lst[i-1], encoding)
        end
    else
        @inbounds for i in 2:length(lst)
            write_encoded_value(w, (lst[i] - lst[i-1]) + T(1), encoding)
        end
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
function read_delta(r::BitReader, encoding::Symbol, ::Type{T}=UInt8; max_elements::Union{Nothing,Int}=nothing, stop_value::Union{Nothing,T}=nothing, vertex_id=nothing, positive_gaps::Bool=false) where {T<:Unsigned}
    lst = T[]
    @debug "READ delta: starting, max_elements=$max_elements, stop_value=$stop_value"

    # read the first value (not delta encoded)
    # NB: we assume that the first value is not 0
    try
        raw_first = vertex_id !== nothing ? read_encoded_value(r, encoding, UInt64) : read_encoded_value(r, encoding, T)
        if vertex_id !== nothing
            # Relative first-value: vertex_id + zigzag_decode(raw - 1)
            first_value = T(Int64(vertex_id) + _zigzag_decode(raw_first - 1))
            @debug "READ delta: read zigzag first raw=$raw_first decoded=$first_value"
        else
            first_value = raw_first
            @debug "READ delta: read first_value=$first_value"
        end
        push!(lst, first_value)
    catch e
        # Empty list case: if we can't even read the first value
        if isa(e, EOFError) || isa(e, ErrorException)
            @debug "READ delta: caught exception reading first value: $e - returning empty list"
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

        # Shift back by 1 unless positive_gaps mode (gaps already ≥ 1, no shift applied)
        push!(lst, lst[end] + value - (positive_gaps ? T(0) : T(1)))
    end
    return lst
end

################################################################################
# Compressed graph data encoding
################################################################################

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
function write_hybrid_mix_encoded_list(w::BitWriter, delta_list::Vector{T}, encoding::Symbol, use_run_length_and_interval::Bool=true, min_interval_length::Int=MIN_INTERVAL_LENGTH, is_children_mode::Bool=false; vertex_id=nothing) where {T<:Unsigned}
    if isempty(delta_list)
        error("write_hybrid_mix_encoded_list should not be called with empty lists")
    end

    # Decide whether to use hybrid mode for this list
    hybrid_active = use_run_length_and_interval && length(delta_list) > 1

    @debug "write_hybrid_mix_encoded_list: list_length=$(length(delta_list)), use_run_length_and_interval=$use_run_length_and_interval, hybrid_active=$hybrid_active"

    # Write hybrid mode flag
    _struct_write_bit(w, hybrid_active)

    # Write the first value (same as write_mix_encoded_list)
    if vertex_id !== nothing
        # Relative first-value: zigzag(v₁ - vertex_id) + 1
        # delta_list[1] is the first absolute value from delta_encode_vector
        offset = Int64(delta_list[1]) - Int64(vertex_id)
        encoded_first = UInt64(_zigzag_encode(offset) + 1)
        write_encoded_value(w, encoded_first, encoding)
    else
        write_encoded_value(w, delta_list[1], encoding)
    end

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
            _struct_write_bit(w, false)
            write_encoded_value(w, T(length(section.data)), encoding)
            for val in section.data
                write_encoded_value(w, val, encoding)
            end
        elseif section.type == :run_length
            # Run-length section: flag=1,0, count, value/length pairs
            _struct_write_bit(w, true)
            _struct_write_bit(w, false)
            write_encoded_value(w, T(length(section.data) ÷ 2), encoding)  # number of pairs
            for val in section.data
                write_encoded_value(w, val, encoding)
            end
        elseif section.type == :interval
            # Interval section: flag=1,1, count, start/length pairs
            _struct_write_bit(w, true)
            _struct_write_bit(w, true)
            write_encoded_value(w, T(length(section.data) ÷ 2), encoding)  # number of pairs
            for val in section.data
                write_encoded_value(w, val, encoding)
            end
        end
    end
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
function read_hybrid_mix_encoded_list(r::BitReader, coding_scheme::Symbol, integer_encoding::Symbol, ::Type{T}=UInt8; max_elements::Union{Int,Nothing}=nothing, stop_value::Union{T,Nothing}=nothing, vertex_id=nothing) where {T<:Unsigned}
    # Read hybrid mode flag
    use_run_length_and_interval = _struct_read_bit(r)

    # Read the first value (same as read_mix_encoded_list)
    raw_first = vertex_id !== nothing ? read_encoded_value(r, integer_encoding, UInt64) : read_encoded_value(r, integer_encoding, T)

    if vertex_id !== nothing
        # Relative first-value: vertex_id + zigzag_decode(raw - 1)
        first_value = T(Int64(vertex_id) + _zigzag_decode(raw_first - 1))
    else
        first_value = raw_first
    end

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
        section_flag = _struct_read_bit(r)
        
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
            second_flag = _struct_read_bit(r)
            
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

    # bitmap cost (adaptive: 1 flag bit + min(raw, block))
    raw_cost = estimate_encoded_value_cost(UInt32(length(copy_bitmap)) + UInt32(1), integer_encoding) + length(copy_bitmap)
    block_cost = estimate_block_encoding_cost(copy_bitmap, integer_encoding)
    cost += 1 + min(raw_cost, block_cost)

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
    estimate_interval_runlength_encoding_cost(neighbors::Vector{T}, encoding::Symbol, min_interval_length::Int, min_run_length::Int) where {T<:Unsigned}

Estimate the bit cost of encoding a vertex using interval + run-length + residual encoding.
Returns the estimated number of bits required.

@param neighbors::Vector{T}: the sorted neighbor list
@param integer_encoding::Symbol: the integer encoding used (:elias_gamma, :elias_delta, :golomb, :fibonacci, :zeta)
@param min_interval_length::Int: minimum interval length (default: 4)
@param min_run_length::Int: minimum run-length (default: 3)
@return::Int: the estimated cost in bits
"""
function estimate_interval_runlength_encoding_cost(neighbors::Vector{T}, integer_encoding::Symbol, min_interval_length::Int, min_run_length::Int; vertex_id=nothing, tight_intervals::Bool=false) where {T<:Unsigned}
    if isempty(neighbors)
        # Empty list: just the count (0+1=1)
        return estimate_encoded_value_cost(T(1), integer_encoding)
    end

    intervals, residuals = compress_intervals(neighbors, min_interval_length)

    cost = 0

    # 1. Number of intervals (+1 to avoid zero, matching write_intervals_and_residuals)
    cost += estimate_encoded_value_cost(T(length(intervals)) + T(1), integer_encoding)

    # 2. Interval data: (start, length) pairs with delta encoding
    prev_ref = T(0)
    for (idx, (start, len)) in enumerate(intervals)
        if idx == 1 && vertex_id !== nothing
            encoded_start = UInt64(_zigzag_encode(Int64(start) - Int64(vertex_id)) + 1)
            cost += estimate_encoded_value_cost(encoded_start, integer_encoding)
        else
            # Delta from previous reference point (start or end depending on tight_intervals)
            cost += estimate_encoded_value_cost(start - prev_ref, integer_encoding)
        end
        cost += estimate_encoded_value_cost(T(len - min_interval_length + 1), integer_encoding)
        prev_ref = tight_intervals ? (start + T(len)) : start
    end

    # 3. Number of residuals (+1 to avoid zero, matching write_intervals_and_residuals)
    cost += estimate_encoded_value_cost(T(length(residuals)) + T(1), integer_encoding)

    # 4. Residuals (delta-encoded, matching write_delta: gaps shifted +1)
    if !isempty(residuals)
        if vertex_id !== nothing
            encoded_first = UInt64(_zigzag_encode(Int64(residuals[1]) - Int64(vertex_id)) + 1)
            cost += estimate_encoded_value_cost(encoded_first, integer_encoding)
        else
            cost += estimate_encoded_value_cost(residuals[1], integer_encoding)
        end
        @inbounds for i in 2:length(residuals)
            cost += estimate_encoded_value_cost((residuals[i] - residuals[i-1]) + T(1), integer_encoding)
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


################################################################################
# START: read / write compressed graph data
################################################################################

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

@returns nothing
"""
function write_compressed_graph_data(w::BitWriter, neighbor_lists::Dict{T,Vector{T}}, coding_scheme::Symbol=:children, integer_encoding::Symbol=:elias_delta, use_mix_mode::Bool=true, reference_enabled::Bool=true, recursive_reference::Bool=true, ref_window_size::Int=REF_WINDOW_SIZE) where {T<:Unsigned}
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

    # stats counters
    reference_chosen_count = 0        # vertices encoded using reference (after cost compare)
    hybrid_chosen_count = 0           # vertices encoded using hybrid (no ref or ref lost)
    # fast/slow lookup counters (top-level + recursive)
    fast_hit_count_ref = Ref(0)
    slow_hit_count_ref = Ref(0)
    # recursive reference usage counters
    recursive_vertices_count = 0      # number of vertices where recursive residual encoding was used at least once
    recursive_event_total_ref = Ref(0) # total number of recursive selections across all levels
    # interval encoding counters
    interval_first_count = 0          # vertices using interval encoding
    total_interval_neighbors = 0      # total neighbors encoded in intervals
    total_run_length_neighbors = 0    # total neighbors encoded in run-length pairs
    total_residual_neighbors = 0      # total neighbors encoded as residuals
    
    # Dense mask of available references (by vertex id)
    available_mask = falses(Int(vs) + 1)
    # Workspace to build reference data without allocations
    ref_build_work = RefBuildWorkspace{T}()
    
    # Helper function to add vertex to reference window (WebGraph-style sliding window)
    function add_to_reference_window!(vertex::T)
        push!(reference_window, vertex)
        push!(potential_references, vertex)
        available_mask[vertex] = true
        
        # Maintain sliding window of size ref_window_size
        if length(reference_window) > ref_window_size
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

                # Write empty interval+RL data (three 1s for zero counts)
                write_intervals_and_residuals(w, current_neighbors, integer_encoding, MIN_INTERVAL_LENGTH)

                # Write trailing stop value (required in children mode)
                @debug "STOP VALUE: Writing stop value for empty list vertex $v"
                write_encoded_value(w, T(1), integer_encoding)

                # Update counters - empty lists use interval encoding
                interval_first_count += 1

                # continue to the next vertex
                continue
            end
        end
        
        # sort neighbors in place
        sorted_neighbors = sort!(current_neighbors)

        # Track whether this vertex used reference encoding (needed for stop value logic)
        # Note: use_reference_final is different from use_reference:
        #   - use_reference: true if a suitable reference vertex was found
        #   - use_reference_final: true if reference encoding was chosen after cost comparison
        #   Both are needed because even if a reference is found, it might be rejected if interval+RL is cheaper
        use_reference_final = false

        if reference_enabled
            # check if we should use reference encoding
            use_reference = false
            best_ref = nothing
            copy_bitmap = Bool[]
            residuals = T[]
            
            if !isempty(potential_references)
                reference_queries += 1
                # WebGraph-style greedy cost-based selection
                best_ref = find_best_reference_greedy_cost(
                    sorted_neighbors,
                    neighbor_lists,
                    potential_references,
                    ref_build_work,
                    coding_scheme,
                    integer_encoding
                )
                @debug "find_best_reference_greedy_cost result" best_ref=best_ref

                if best_ref !== nothing
                    references_found += 1
                    # copy_bitmap and residuals are already computed in ref_build_work
                    # by find_best_reference_greedy_cost
                    copy_bitmap = ref_build_work.copy_bitmap
                    residuals = ref_build_work.residuals
                    use_reference = true
                    fast_hit_count_ref[] += 1  # Count greedy selection as "fast"
                end
            end

            if use_reference
                # Compare reference encoding cost vs interval+RL encoding cost
                ref_cost = estimate_reference_encoding_cost(best_ref, copy_bitmap, residuals, coding_scheme, integer_encoding)
                interval_rl_cost = estimate_interval_runlength_encoding_cost(sorted_neighbors, integer_encoding, MIN_INTERVAL_LENGTH, MIN_RUN_LENGTH)

                # Choose the more efficient encoding
                use_reference_final = ref_cost <= interval_rl_cost
                @debug "Encoding comparison: vertex=$v, ref_cost=$ref_cost, interval_rl_cost=$interval_rl_cost, chosen=$(use_reference_final ? "reference" : "interval+RL"), savings=$(abs(ref_cost - interval_rl_cost)) bits"
                
                if use_reference_final
                    reference_chosen_count += 1
                    @debug "CHOICE: Reference encoding chosen: vertex=$v, ref_id=$best_ref, copy_bitmap_len=$(length(copy_bitmap)), residuals_len=$(length(residuals))"
                    # children_flag = 1 (reference mode)
                    write_bit(w, true)
                else
                    # HYBRID MIX ENCODING DISABLED - use interval+RL instead
                    hybrid_chosen_count += 1
                    @debug "CHOICE: Interval+RL encoding chosen over reference: vertex=$v, neighbors_count=$(length(sorted_neighbors))"

                    # Compute interval+RL encoding
                    intervals_alt, residuals_alt = compress_intervals(sorted_neighbors, MIN_INTERVAL_LENGTH)
                    # DISABLED: run_length_pairs_alt, final_residuals_alt = compress_run_length(residuals_alt, MIN_RUN_LENGTH)
                    run_length_pairs_alt = Tuple{T,T}[]  # No RL encoding
                    final_residuals_alt = residuals_alt   # Pass through unchanged

                    interval_neighbor_count = sum(len for (_, len) in intervals_alt; init=0)
                    run_length_neighbor_count = sum(count for (_, count) in run_length_pairs_alt; init=0)

                    interval_first_count += 1
                    total_interval_neighbors += interval_neighbor_count
                    total_run_length_neighbors += run_length_neighbor_count
                    total_residual_neighbors += length(final_residuals_alt)

                    # children_flag = 0 (not using reference)
                    write_bit(w, false)
                    # Write intervals, run-length pairs, and residuals
                    write_intervals_and_residuals(w, sorted_neighbors, integer_encoding, MIN_INTERVAL_LENGTH)
                    add_to_reference_window!(T(v))
                end
                
                if use_reference_final
                    # write reference data
                    write_encoded_value(w, T(best_ref), integer_encoding)  # ref_id
                    write_bitmap_adaptive(w, copy_bitmap, integer_encoding)  # bitmap (adaptive encoding)
                    
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
                            # Use hybrid mix encoding for better compression in both modes
                            # write_hybrid_mix_encoded_list(w, residual_deltas, integer_encoding, use_mix_mode, MIN_INTERVAL_LENGTH, coding_scheme == :children)
                            # Write residuals directly for simplicity
                            for residual in residual_deltas
                                write_encoded_value(w, residual, integer_encoding)
                            end
                            # In children mode, write stop value to terminate the residuals
                            if coding_scheme == :children
                                write_encoded_value(w, T(1), integer_encoding)
                            end
                        end
                    else
                        # residuals_flag = 0 (no residuals)
                        write_bit(w, false)

                        # Write trailing stop value for reference vertices with no residuals (children mode only)
                        if coding_scheme == :children
                            @debug "STOP VALUE: Writing stop value for reference vertex $v with no residuals"
                            write_encoded_value(w, T(1), integer_encoding)
                        end
                    end

                    add_to_reference_window!(T(v))
                end
            # no relevant reference: use interval + residuals encoding
            else
                hybrid_chosen_count += 1

                # Use direct interval + run-length + residual encoding
                # This avoids the roundtrip: neighbors → deltas → reconstruct → find intervals
                intervals, residuals = compress_intervals(sorted_neighbors, MIN_INTERVAL_LENGTH)
                # DISABLED: run_length_pairs, final_residuals = compress_run_length(residuals, MIN_RUN_LENGTH)
                run_length_pairs = Tuple{T,T}[]  # No RL encoding
                final_residuals = residuals       # Pass through unchanged

                interval_neighbor_count = sum(len for (_, len) in intervals; init=0)
                run_length_neighbor_count = sum(count for (_, count) in run_length_pairs; init=0)

                interval_first_count += 1
                total_interval_neighbors += interval_neighbor_count
                total_run_length_neighbors += run_length_neighbor_count
                total_residual_neighbors += length(final_residuals)

                @debug "CHOICE: Interval+RL encoding: vertex=$v, neighbors=$(length(sorted_neighbors)), intervals=$(length(intervals)) ($(interval_neighbor_count) edges), rl_pairs=$(length(run_length_pairs)) ($(run_length_neighbor_count) edges), residuals=$(length(final_residuals))"

                # children_flag = 0 (not using reference)
                write_bit(w, false)
                # Write intervals, run-length pairs, and residuals
                write_intervals_and_residuals(w, sorted_neighbors, integer_encoding, MIN_INTERVAL_LENGTH)
                add_to_reference_window!(T(v))
            end
        # reference disabled: use interval + run-length + residuals encoding
        else
            # Use direct interval + run-length + residual encoding
            intervals, residuals = compress_intervals(sorted_neighbors, MIN_INTERVAL_LENGTH)
            # DISABLED: run_length_pairs, final_residuals = compress_run_length(residuals, MIN_RUN_LENGTH)
            run_length_pairs = Tuple{T,T}[]  # No RL encoding
            final_residuals = residuals       # Pass through unchanged

            interval_neighbor_count = sum(len for (_, len) in intervals; init=0)
            run_length_neighbor_count = sum(count for (_, count) in run_length_pairs; init=0)

            interval_first_count += 1
            total_interval_neighbors += interval_neighbor_count
            total_run_length_neighbors += run_length_neighbor_count
            total_residual_neighbors += length(final_residuals)

            @debug "CHOICE: Interval+RL encoding (ref disabled): vertex=$v, neighbors=$(length(sorted_neighbors)), intervals=$(length(intervals)) ($(interval_neighbor_count) edges), rl_pairs=$(length(run_length_pairs)) ($(run_length_neighbor_count) edges), residuals=$(length(final_residuals))"

            # Write intervals, run-length pairs, and residuals
            write_intervals_and_residuals(w, sorted_neighbors, integer_encoding, MIN_INTERVAL_LENGTH)
        end
        
        # progress logging
        if v % progress_interval == 0
            elapsed = time() - data_section_start
            @debug "  Progress: $v/$vs vertices ($(round(100*v/vs, digits=1))%), $(length(reference_window))/$(REF_WINDOW_SIZE) window refs, $reference_queries queries, $references_found found, $(round(elapsed, digits=3))s elapsed"
        end

        # write stop value after each vertex list in children mode
        # IMPORTANT: only for non-reference vertices, because reference vertices
        # handle their own stop values (within the residuals encoding if present,
        # or explicitly consumed in the decoder at line 3060-3063 if not present)
        if coding_scheme == :children && !use_reference_final
            # NB: the list is not empty, so the stop value is not the first value
            # and a vertex flag is needed
            @debug "STOP VALUE: Writing common stop value for non-reference vertex $v (use_mix_mode=$use_mix_mode)"
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

    # Interval+Run-length encoding statistics
    interval_vertex_ratio = total_vertices > 0 ? (100.0 * interval_first_count / total_vertices) : 0.0
    total_encoded_edges = total_interval_neighbors + total_run_length_neighbors + total_residual_neighbors
    interval_edge_density = total_encoded_edges > 0 ? (100.0 * total_interval_neighbors / total_encoded_edges) : 0.0
    run_length_edge_density = total_encoded_edges > 0 ? (100.0 * total_run_length_neighbors / total_encoded_edges) : 0.0
    residual_edge_density = total_encoded_edges > 0 ? (100.0 * total_residual_neighbors / total_encoded_edges) : 0.0

    @info "Interval+Run-length encoding summary:"
    @info "  Vertices using interval+RL encoding: $interval_first_count ($(round(interval_vertex_ratio, digits=1))% of all)"
    @info "  Total edges encoded: $total_encoded_edges"
    @info "    Intervals: $total_interval_neighbors ($(round(interval_edge_density, digits=1))%)"
    @info "    Run-length: $total_run_length_neighbors ($(round(run_length_edge_density, digits=1))%)"
    @info "    Residuals: $total_residual_neighbors ($(round(residual_edge_density, digits=1))%)"
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
                        copy_bitmap = read_bitmap_adaptive(r, integer_encoding)
                        
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
                        # The writer uses write_intervals_and_residuals for non-reference vertices
                        current_neighbors = read_intervals_and_residuals(r, integer_encoding, MIN_INTERVAL_LENGTH, T)
                    end
                else
                    # reference disabled - read the appropriate format
                    # The writer uses write_intervals_and_residuals for all vertices when reference is disabled
                    current_neighbors = read_intervals_and_residuals(r, integer_encoding, MIN_INTERVAL_LENGTH, T)
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
                    @info "Vertex $v: children_flag=$children_flag"

                    # reference mode: children mode
                    if children_flag
                        @info "Vertex $v: Reading reference data (children_flag=true)"
                        # read reference data
                        ref_id = read_encoded_value(r, integer_encoding, T)
                        @info "Vertex $v: Read ref_id=$ref_id"
                        copy_bitmap = read_bitmap_adaptive(r, integer_encoding)
                        @info "Vertex $v: Read bitmap with $(length(copy_bitmap)) bits"

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
                    # non-reference mode: use interval+RL encoding (hybrid mix disabled)
                    else
                        @info "Vertex $v: Using interval encoding (children_flag=0)"
                        current_neighbors = read_intervals_and_residuals(r, integer_encoding, MIN_INTERVAL_LENGTH, T)
                        @info "Vertex $v: Read $(length(current_neighbors)) neighbors via intervals"
                        # In children mode, consume the trailing stop value after each vertex list
                        _consume_children_trailing_stop(r, integer_encoding, T)
                        @info "Vertex $v: Consumed stop value"
                    end
                # reference disabled: read the appropriate format
                else
                    # The writer uses write_intervals_and_residuals for all vertices when reference is disabled
                    current_neighbors = read_intervals_and_residuals(r, integer_encoding, MIN_INTERVAL_LENGTH, T)
                    # In children mode, consume the trailing stop value after each vertex list
                    if coding_scheme == :children
                        _consume_children_trailing_stop(r, integer_encoding, T)
                    end
                end
                
                # successfully read vertex data
                neighbor_lists[T(v)] = current_neighbors
                vertex_data_read = true
            catch e
                # EOF while reading vertex data: current vertex and all remaining vertices are empty
                @warn "Caught exception at vertex $v: $e"
                @warn "Exception type: $(typeof(e))"
                if isa(e, EOFError) || isa(e, ErrorException)
                    @warn "Treating as EOF - marking remaining $(vs-v+1) vertices as empty"
                    for remaining_v in v:vs
                        neighbor_lists[T(remaining_v)] = T[]
                    end
                    break
                else
                    @warn "Rethrowing exception"
                    rethrow(e)
                end
            end 
            # NB: stop values are consumed by read_mix_encoded_list when stop_value parameter is provided
        end
    end
    
    return neighbor_lists
end

################################################################################
# END: read / write compressed graph data
################################################################################

################################################################################
# START: Greedy compressed graph data
################################################################################

# 3-bit encoding tags for self-describing greedy-compressed streams
const ENC_TAG_FIBONACCI = UInt8(0b000)
const ENC_TAG_ZETA = UInt8(0b001)
const ENC_TAG_ELIAS_GAMMA = UInt8(0b010)
const ENC_TAG_ELIAS_DELTA = UInt8(0b011)

const ENCODING_TAGS = Dict{Symbol,UInt8}(
    :fibonacci => ENC_TAG_FIBONACCI,
    :zeta => ENC_TAG_ZETA,
    :elias_gamma => ENC_TAG_ELIAS_GAMMA,
    :elias_delta => ENC_TAG_ELIAS_DELTA
)

const TAG_ENCODINGS = Dict{UInt8,Symbol}(
    ENC_TAG_FIBONACCI => :fibonacci,
    ENC_TAG_ZETA => :zeta,
    ENC_TAG_ELIAS_GAMMA => :elias_gamma,
    ENC_TAG_ELIAS_DELTA => :elias_delta
)

# Reference mode tags (2 bits)
const REF_NONE = UInt8(0b00)
const REF_REFERENCE = UInt8(0b01)
const REF_RECURSIVE = UInt8(0b10)

const REF_MODE_TAGS = Dict{Symbol,UInt8}(
    :none => REF_NONE,
    :reference => REF_REFERENCE,
    :recursive => REF_RECURSIVE
)

const TAG_REF_MODES = Dict{UInt8,Symbol}(
    REF_NONE => :none,
    REF_REFERENCE => :reference,
    REF_RECURSIVE => :recursive
)

function _write_encoding_tag(w, tag::UInt8)
    write_bit(w, (tag >> 2) & 0x1 == 1)
    write_bit(w, (tag >> 1) & 0x1 == 1)
    write_bit(w, tag & 0x1 == 1)
end

function _read_encoding_tag(r)::UInt8
    b2 = read_bit(r) ? UInt8(1) : UInt8(0)
    b1 = read_bit(r) ? UInt8(1) : UInt8(0)
    b0 = read_bit(r) ? UInt8(1) : UInt8(0)
    return (b2 << 2) | (b1 << 1) | b0
end

function _write_ref_mode_tag(w, tag::UInt8)
    write_bit(w, (tag >> 1) & 0x1 == 1)
    write_bit(w, tag & 0x1 == 1)
end

function _read_ref_mode_tag(r)::UInt8
    b1 = read_bit(r) ? UInt8(1) : UInt8(0)
    b0 = read_bit(r) ? UInt8(1) : UInt8(0)
    return (b1 << 1) | b0
end

# Min interval length tag: 2 bits encoding index into [2, 3, 4, 5]
const MIL_OPTIONS = [2, 3, 4, 5]

# Default encoding for greedy-compressed streams
const _FIXED_ENCODING = :fibonacci

# ── Cost model levels ────────────────────────────────────────────────────
# COST_MODEL_FULL (0): evaluate all encoding options (interval×4 MILs + RLE×4 MILs + delta)
# COST_MODEL_FAST (1): skip RLE, single MIL, simpler copy cost estimation
const COST_MODEL_FULL = 0
const COST_MODEL_FAST = 1

# Default cost model — can be overridden via ADJACENTLY_COST_MODEL env var
const DEFAULT_COST_MODEL = let
    v = get(ENV, "ADJACENTLY_COST_MODEL", "0")
    parse(Int, v)
end

# Total bits written so far (flushed bytes + buffered bits in accumulator)
@inline function _total_bits(w::BitWriter)
    return w.bit_count
end

# Flattened Encoding Options (Type, MIL) — full model
const ENCODING_OPTIONS = vcat(
    [(:interval, mil) for mil in MIL_OPTIONS],
    [(:rle, mil) for mil in MIL_OPTIONS],
    [(:delta, 0)]
)

# Fast encoding options: only interval(mil=4) + delta — skips all RLE
const FAST_ENCODING_OPTIONS = [(:interval, 4), (:delta, 0)]

# Max candidates for phase 2 of two-phase reference search
# Full model: evaluate all candidates (low cap was causing BPE regression)
const MAX_REF_CANDIDATES_PHASE2 = 1024
# Fast model: limit evaluation for speed
const MAX_REF_CANDIDATES_PHASE2_FAST = 8

# Fixed MIL for adaptive header mode (intervals path)
const ADAPTIVE_MIL = 4

# Encoding options for adaptive header mode: only STOP-delta and intervals(ADAPTIVE_MIL)
const ADAPTIVE_ENCODING_OPTIONS = [(:delta, 0), (:interval, ADAPTIVE_MIL)]

function _write_enc_opt_tag(w, enc_type::Symbol, mil::Int)
    tuple = (enc_type, mil)
    if enc_type == :delta; tuple = (:delta, 0); end
    idx = findfirst(==(tuple), ENCODING_OPTIONS)
    if idx === nothing; idx = 3; end # Default to (:interval, 4)
    write_value(w, UInt8(idx - 1), 4) # 4 bits for 9 options
end

function _read_enc_opt_tag(r)
    idx = Int(read_value(r, 4, UInt8)) + 1
    return ENCODING_OPTIONS[idx]
end

function _write_mil_tag(w, mil::Int)
    idx = findfirst(==(mil), MIL_OPTIONS)
    if idx === nothing; idx = 3; end  # default to 4
    tag = UInt8(idx - 1)  # 0-based
    write_bit(w, (tag >> 1) & 0x1 == 1)
    write_bit(w, tag & 0x1 == 1)
end

function _read_mil_tag(r)::Int
    b1 = read_bit(r) ? UInt8(1) : UInt8(0)
    b0 = read_bit(r) ? UInt8(1) : UInt8(0)
    idx = Int((b1 << 1) | b0) + 1  # 1-based
    return MIL_OPTIONS[idx]
end

# Zigzag encoding for relative first-value encoding in greedy path.
# Maps signed offsets to positive integers: 0→0, -1→1, +1→2, -2→3, +2→4, ...
_zigzag_encode(n::Int64)::UInt64 = n >= 0 ? UInt64(2n) : UInt64(2*(-n) - 1)
_zigzag_decode(v::UInt64)::Int64 = iseven(v) ? Int64(v >> 1) : -Int64((v >> 1) + 1)

# ===========================================================================
# Variable-Length Coded (VLC) vertex headers
#
# VLC v1/v2: The top 4 action combos cover ~97% of vertices in well-ordered
# web graphs. Short prefix codes for them, escape (9-10 bits) for the rest.
# VLC v3: All interval mils get equal cost (no escape penalty for mil≠2).
#
#   === VLC v3 ===
#   Code         Bits  Meaning
#   ----         ----  -------
#   0             1    reference + delta
#   10            2    none      + delta
#   110  + mm     5    none      + interval (mm = mil index, 2 bits)
#   1110 + mm     6    reference + interval (mm = mil index, 2 bits)
#   11110 + mm    7    none      + rle      (mm = mil index, 2 bits)
#   11111 + mm    7    reference + rle      (mm = mil index, 2 bits)
#
#   === VLC v1/v2 ===
#
#   Code   Bits  Meaning
#   ----   ----  -------
#   00      2    none      + delta
#   01      2    none      + interval (mil=2)
#   10      2    reference + delta
#   110     3    reference + interval (mil=2)
#   111     3+6  escape → full 6-bit header (ref_mode 2b + enc_opt 4b)
#
# Expected average: ~2.3 bits/vertex (vs 6.0 fixed), a ~62% reduction.
# ===========================================================================

# VLC v3: all interval/rle mils have equal cost, no expensive escape.
# VLC vertex header prefix codes with merged empty (see FORMAT_SPECS_STD.md).
# Codes:
#   0            (1)  ref + delta
#   10           (2)  none + delta
#   110          (3)  empty (no neighbors)
#   1110+mm      (6)  none + interval
#   11110+mm     (7)  ref + interval
#   111110+mm    (8)  none + rle
#   1111110+mm   (9)  ref + rle
#   11111110     (8)  multi-ref + delta
#   111111110+mm (11) multi-ref + interval
#   111111111+mm (11) multi-ref + rle
function _write_vertex_header_vlc(w, ref_mode::Symbol, enc_type::Symbol, mil::Int; adaptive::Bool=false)
    if adaptive
        # Adaptive header: 2 bits — ref_flag + adaptive_flag (0=STOP-delta, 1=intervals)
        write_bit(w, ref_mode == :reference)
        write_bit(w, enc_type == :interval)
        return
    end
    if ref_mode == :reference && enc_type == :delta
        write_bit(w, false)                                                           # 0
    elseif ref_mode == :none && enc_type == :delta
        write_bit(w, true); write_bit(w, false)                                       # 10
    elseif ref_mode == :empty
        write_bit(w, true); write_bit(w, true); write_bit(w, false)                   # 110
    elseif ref_mode == :none && enc_type == :interval
        write_bit(w, true); write_bit(w, true); write_bit(w, true); write_bit(w, false)  # 1110
        _write_mil_tag(w, mil)
    elseif ref_mode == :reference && enc_type == :interval
        for _ in 1:4; write_bit(w, true); end; write_bit(w, false)                    # 11110
        _write_mil_tag(w, mil)
    elseif ref_mode == :none && enc_type == :rle
        for _ in 1:5; write_bit(w, true); end; write_bit(w, false)                    # 111110
        _write_mil_tag(w, mil)
    elseif ref_mode == :reference && enc_type == :rle
        for _ in 1:6; write_bit(w, true); end; write_bit(w, false)                    # 1111110
        _write_mil_tag(w, mil)
    elseif ref_mode == :multi_ref && enc_type == :delta
        for _ in 1:7; write_bit(w, true); end; write_bit(w, false)                    # 11111110
    elseif ref_mode == :multi_ref && enc_type == :interval
        for _ in 1:8; write_bit(w, true); end; write_bit(w, false)                    # 111111110
        _write_mil_tag(w, mil)
    else  # multi_ref + rle
        for _ in 1:8; write_bit(w, true); end; write_bit(w, true)                     # 111111111
        _write_mil_tag(w, mil)
    end
end

function _read_vertex_header_vlc(r; adaptive::Bool=false)
    if adaptive
        # Adaptive header: 2 bits — ref_flag + adaptive_flag
        ref_flag = read_bit(r)
        interval_flag = read_bit(r)
        ref_mode = ref_flag ? :reference : :none
        enc_type = interval_flag ? :interval : :delta
        return (ref_mode, enc_type, interval_flag ? ADAPTIVE_MIL : 0)
    end
    if !read_bit(r);                           return (:reference, :delta, 0); end      # 0
    if !read_bit(r);                           return (:none, :delta, 0); end            # 10
    if !read_bit(r);                           return (:empty, :delta, 0); end           # 110
    if !read_bit(r); mil = _read_mil_tag(r);   return (:none, :interval, mil); end       # 1110+mm
    if !read_bit(r); mil = _read_mil_tag(r);   return (:reference, :interval, mil); end  # 11110+mm
    if !read_bit(r); mil = _read_mil_tag(r);   return (:none, :rle, mil); end            # 111110+mm
    if !read_bit(r); mil = _read_mil_tag(r);   return (:reference, :rle, mil); end       # 1111110+mm
    if !read_bit(r);                           return (:multi_ref, :delta, 0); end       # 11111110
    if !read_bit(r); mil = _read_mil_tag(r);   return (:multi_ref, :interval, mil); end  # 111111110+mm
    mil = _read_mil_tag(r);                    return (:multi_ref, :rle, mil)             # 111111111+mm
end

# Phase 2b: BG vertex header as a single small-int symbol for the context-modeled
# command stream. Non-adaptive alphabet = the 28 merged-VLC actions; adaptive
# alphabet = the 4 (ref_flag × interval_flag) combos. All ids < 2^K → direct tokens.
const _BG_CMD_ACTIONS = begin
    acts = Tuple{Symbol,Symbol,Int}[
        (:empty, :delta, 0), (:reference, :delta, 0),
        (:none, :delta, 0), (:multi_ref, :delta, 0),
    ]
    for rm in (:none, :reference, :multi_ref), et in (:interval, :rle), m in MIL_OPTIONS
        push!(acts, (rm, et, m))
    end
    acts
end
const _BG_CMD_INDEX = Dict{Tuple{Symbol,Symbol,Int},Int}(a => i - 1 for (i, a) in enumerate(_BG_CMD_ACTIONS))

function _bg_cmd_id(ref_mode::Symbol, enc_type::Symbol, mil::Int; adaptive::Bool=false)::Int
    if adaptive
        return (ref_mode == :reference ? 1 : 0) | (enc_type == :interval ? 2 : 0)
    end
    ref_mode == :empty && return 0
    if enc_type == :delta
        return _BG_CMD_INDEX[(ref_mode, :delta, 0)]
    end
    # Mirror _write_mil_tag's fallback: unknown mil → MIL_OPTIONS[3]
    m = mil in MIL_OPTIONS ? mil : MIL_OPTIONS[3]
    return _BG_CMD_INDEX[(ref_mode, enc_type, m)]
end

function _bg_cmd_decode(id::Int; adaptive::Bool=false)
    if adaptive
        ref_mode = (id & 1) != 0 ? :reference : :none
        interval = (id & 2) != 0
        return (ref_mode, interval ? :interval : :delta, interval ? ADAPTIVE_MIL : 0)
    end
    return _BG_CMD_ACTIONS[id + 1]
end

@inline function _vlc_header_cost(ref_mode::Symbol, enc_type::Symbol, mil::Int; adaptive::Bool=false)::Int
    if adaptive; return 2; end  # Always 2 bits: ref_flag + adaptive_flag
    if ref_mode == :reference && enc_type == :delta;    return 1; end   # 0
    if ref_mode == :none && enc_type == :delta;         return 2; end   # 10
    if ref_mode == :empty;                              return 3; end   # 110
    if ref_mode == :none && enc_type == :interval;      return 6; end   # 1110+mm
    if ref_mode == :reference && enc_type == :interval;  return 7; end  # 11110+mm
    if ref_mode == :none && enc_type == :rle;            return 8; end  # 111110+mm
    if ref_mode == :reference && enc_type == :rle;       return 9; end  # 1111110+mm
    if ref_mode == :multi_ref && enc_type == :delta;     return 8; end  # 11111110
    if ref_mode == :multi_ref && enc_type == :interval;  return 11; end # 111111110+mm
    return 11  # multi_ref + rle: 111111111+mm
end

"""
    read_greedy_graph_data(r, vs, coding_scheme, T; integer_encoding)

Read greedy-compressed graph data. The stream is self-describing with per-vertex headers.

Format:
- 1 bit: coding_scheme flag (true = :index, false = :children)
- 3 bits: encoding tag
- 2 bits: default ref_mode tag
- 2 bits: default mil tag
- Index section (if :index mode)
- Data section: per-vertex encoding with optional reference

Parameters:
- r: BitReader
- vs: number of vertices
- coding_scheme: :children or :index (used for fallback/verification)
- T: vertex type
- integer_encoding: default encoding (overridden by stream header)

Returns: Dict{T, Vector{T}} of neighbor lists
"""
function write_greedy_graph_data(w, neighbor_lists::Dict{T,Vector{T}},
        coding_scheme::Symbol=:children,
        ref_window_size::Int=7;
        integer_encoding::Symbol=_FIXED_ENCODING,
        stats::Union{Dict,Nothing}=nothing,
        copy_blocks::Bool=false,
        adaptive_copy::Bool=false,
        fixwidth_ref::Bool=false,
        stop_deltas::Bool=false,
        adaptive_deltas::Bool=false,
        split_residual::Bool=false,
        bv_blocks::Bool=false,
        compact_copy::Bool=false,
        tight_intervals::Bool=false,
        exact_costing::Bool=false,
        lr_split::Bool=false,
        multi_ref::Bool=false,
        adaptive_header::Bool=false,
        cost_model::Int=DEFAULT_COST_MODEL,
        index_sample_k::Int=0) where {T<:Unsigned}

    vs = length(keys(neighbor_lists))
    # :context_range routes residual-region integers to a streaming range coder;
    # structural integers (refdist, copy blocks, counts) still use fibonacci in the
    # bitstream. Residuals are collected across all vertices, finalized, returned.
    use_ctx_range = integer_encoding == :context_range
    ie = use_ctx_range ? :fibonacci : integer_encoding
    resid_sink   = use_ctx_range ? CtxRangeEncoder() : nothing
    refdist_sink = use_ctx_range ? CtxRangeEncoder() : nothing   # reference distances
    copy_sink    = use_ctx_range ? BinRangeEncoder() : nothing   # copy bitmap bits
    # Phase 2b: vertex headers → cmd range coder; stop-delta / RLE / split-residual /
    # adaptive-delta flags → binary coder. In sampled-index (RA) mode both streams are
    # chunked alongside resid/refdist/copy (headers decode sequentially within a chunk,
    # so they don't need to stay in the seekable struct stream).
    cmd_sink  = use_ctx_range ? CtxRangeEncoder() : nothing
    flag_sink = use_ctx_range ? BinRangeEncoder() : nothing
    if use_ctx_range && coding_scheme == :index && index_sample_k <= 0
        # Full (per-vertex) index needs a continuous range blob → not seekable.
        error(":context_range index mode requires a sampled index (index_sample_k > 0)")
    end

    # Index mode: full index uses count-prefixed (no stop_deltas); sampled uses children-mode data
    idx_mode = coding_scheme == :index
    sampled_idx = idx_mode && index_sample_k > 0
    # Chunked range streams: sampled index + :context_range. Each k-vertex chunk
    # finalizes the range coders into an independently-decodable byte blob so a
    # reader can seek to the chunk containing v (granularity k) instead of decoding
    # the whole graph. Per-chunk blob byte-offset tables are stored in the struct
    # stream right after the sampled struct-offset table.
    ctx_chunked = use_ctx_range && sampled_idx
    if idx_mode && !sampled_idx
        # Full index: override stop_deltas (redundant with offset table)
        stop_deltas = false
        adaptive_deltas = false
    end
    # Sampled index: keep stop_deltas as-is (data is children-mode, only the offset table changes)

    # Create cost buffer for exact bit-counting when requested
    gcb = exact_costing ? GreedyCostBuffer() : nothing

    # Pre-allocated workspace for vertex search
    _ws = VertexSearchWorkspace{T}()

    reference_window = T[]
    function add_to_ref_window!(vertex::T)
        push!(reference_window, vertex)
        if length(reference_window) > ref_window_size
            popfirst!(reference_window)
        end
    end

    # Fixed-width ref delta: number of bits to encode (distance-1) in [0, window-1]
    ref_dist_bits = fixwidth_ref ? ceil(Int, log2(max(2, ref_window_size))) : 0

    # Header: Coding scheme (1) + Encoding (3)
    write_bit(w, idx_mode)
    _write_encoding_tag(w, get(ENCODING_TAGS, ie, ENC_TAG_FIBONACCI))

    # --- Helper: encode one vertex's data into a BitWriter ---
    function _encode_vertex!(wr, v::T, current_neighbors::Vector{T})
        # Empty vertex handling
        if isempty(current_neighbors)
            if !idx_mode || sampled_idx
                # Children mode (or sampled index): write explicit empty marker
                if adaptive_header
                    if cmd_sink !== nothing
                        rc_encode_value!(cmd_sink, UInt64(_bg_cmd_id(:none, :delta, 0; adaptive=true)))
                    else
                        _write_vertex_header_vlc(wr, :none, :delta, 0; adaptive=true)
                    end
                    _struct_write_bit(wr, false)  # STOP terminator (0 elements)
                else
                    if cmd_sink !== nothing
                        rc_encode_value!(cmd_sink, UInt64(_bg_cmd_id(:empty, :delta, 0)))
                    else
                        _write_vertex_header_vlc(wr, :empty, :delta, 0)
                    end
                end
            end
            # Full index mode: nothing written (empty = offset[v] == offset[v+1])
            add_to_ref_window!(v)
            return
        end

        # Safety: shouldn't reach here with empty neighbors, but guard anyway
        if idx_mode && isempty(current_neighbors)
            add_to_ref_window!(v)
            return
        end

        actual_ref_mode = :none
        ref_result = nothing
        mil = MIN_INTERVAL_LENGTH
        enc_type = :interval

        # Greedy search
        _, actual_ref_mode, mil, ref_result, enc_type, use_stop, res_enc_type, res_mil = _greedy_vertex_search(current_neighbors, neighbor_lists, reference_window, ie; vertex_id=v, copy_blocks=copy_blocks, adaptive_copy=adaptive_copy, fixwidth_ref=fixwidth_ref, ref_dist_bits=ref_dist_bits, stop_deltas=stop_deltas, adaptive_deltas=adaptive_deltas, split_residual=split_residual, bv_blocks=bv_blocks, compact_copy=compact_copy, tight_intervals=tight_intervals, buf=gcb, lr_split=lr_split, multi_ref=multi_ref, adaptive_header=adaptive_header, ws=_ws, cost_model=cost_model)

        # Write vertex header
        _hb0 = _STREAM_BITS[] === nothing ? 0 : _total_bits(wr)
        if cmd_sink !== nothing
            rc_encode_value!(cmd_sink, UInt64(_bg_cmd_id(actual_ref_mode, enc_type, mil; adaptive=adaptive_header)))
        else
            _write_vertex_header_vlc(wr, actual_ref_mode, enc_type, mil; adaptive=adaptive_header)
        end
        if _STREAM_BITS[] !== nothing
            _sbh = _STREAM_BITS[]::Dict{Symbol,Int}
            _sbh[:header] = get(_sbh, :header, 0) + (_total_bits(wr) - _hb0)
        end

        # Collect action statistics if requested
        if stats !== nothing
            key = (actual_ref_mode, enc_type, mil)
            stats[key] = get(stats, key, 0) + 1
        end

        target_list = current_neighbors
        write_enc = enc_type
        write_mil = mil

        # Helper: write one reference (distance + copy bitmap)
        function _write_one_ref!(wr_out, ref_dist, copy_bmp)
            let t = _REFDIST_TAP[]; t === nothing || push!(t::Vector{UInt64}, UInt64(ref_dist)); end
            let t = _COPYBMP_TAP[]
                if t !== nothing
                    for b in copy_bmp; push!(t::Vector{Bool}, b != 0); end
                end
            end
            local _sbr = _STREAM_BITS[]
            local _rp0 = _sbr === nothing ? 0 : _total_bits(wr_out)
            if refdist_sink !== nothing
                rc_encode_value!(refdist_sink, UInt64(ref_dist))          # -> refdist range coder
            elseif fixwidth_ref
                local d = Int(ref_dist) - 1
                for b in (ref_dist_bits-1):-1:0
                    write_bit(wr_out, ((d >> b) & 1) == 1)
                end
            else
                write_encoded_value(wr_out, T(ref_dist), ie)
            end
            local _rp1 = _sbr === nothing ? 0 : _total_bits(wr_out)
            if copy_sink !== nothing
                for b in copy_bmp; brc_encode_bit!(copy_sink, b != 0); end  # -> copy binary coder
            elseif adaptive_copy && copy_blocks
                _write_adaptive_copy(wr_out, copy_bmp, ie; bv_blocks=bv_blocks, compact_copy=compact_copy)
            elseif copy_blocks
                _write_copy_blocks(wr_out, copy_bmp, ie)
            else
                write_bitmap_adaptive(wr_out, copy_bmp, ie)
            end
            if _sbr !== nothing
                _sbrd = _sbr::Dict{Symbol,Int}
                _sbrd[:refdist] = get(_sbrd, :refdist, 0) + (_rp1 - _rp0)
                _sbrd[:copy] = get(_sbrd, :copy, 0) + (_total_bits(wr_out) - _rp1)
            end
        end

        if actual_ref_mode == :multi_ref && ref_result !== nothing
            # Multi-reference: ref_result = (dist1, bmp1, dist2, bmp2, residuals)
            ref_dist1, copy_bmp1, ref_dist2, copy_bmp2, residuals = ref_result
            _write_one_ref!(wr, ref_dist1, copy_bmp1)
            _write_one_ref!(wr, ref_dist2, copy_bmp2)
            target_list = residuals
        elseif actual_ref_mode == :reference && ref_result !== nothing
            ref_distance, copy_bitmap, residuals = ref_result
            _write_one_ref!(wr, ref_distance, copy_bitmap)
            # Split residual encoding: signal if residuals use a different enc_type
            if split_residual
                same_enc = (res_enc_type == enc_type && (res_mil == mil || res_enc_type == :delta))
                _struct_write_bit(wr, !same_enc)  # 0 = same, 1 = different
                if !same_enc
                    if cmd_sink !== nothing
                        opt = res_enc_type == :delta ? (:delta, 0) : (res_enc_type, res_mil)
                        oidx = findfirst(==(opt), ENCODING_OPTIONS)
                        oidx === nothing && (oidx = 3)
                        rc_encode_value!(cmd_sink, UInt64(oidx - 1))
                    else
                        _write_enc_opt_tag(wr, res_enc_type, res_mil)
                    end
                end
            end
            target_list = residuals
            write_enc = split_residual ? res_enc_type : enc_type
            write_mil = split_residual ? res_mil : mil
        end

        # Copy-aware rank gaps: when the residual range coder is active, recode the
        # residual ids as ranks in the complement of the copied set
        # C = current_neighbors \ residuals. Identity when C is empty (:none), so
        # only reference / multi_ref vertices are affected. The vertex id is mapped
        # into the same rank space so the coder's vertex-relative first value keeps
        # its locality.
        _resid_vid = v
        if resid_sink !== nothing && actual_ref_mode != :none
            _C = setdiff(current_neighbors, target_list)   # sorted, unique, disjoint
            if !isempty(_C)
                target_list = _rankgap_forward(target_list, _C)
                _resid_vid = _rank_of(v, _C)
            end
        end

        # Encode target_list
        _rs0 = _STREAM_BITS[] === nothing ? 0 : _total_bits(wr)
        _rsn0 = (_RESID_META[] === nothing || _RESID_ACC[] === nothing) ? 0 : length(_RESID_ACC[]::Vector{UInt64})
        if _STREAM_BITS[] !== nothing && _RESID_ACC[] !== nothing
            _VALUE_TAP[] = _RESID_ACC[]   # capture only residual integers
        end
        if resid_sink !== nothing
            rc_reset_region!(resid_sink)
            _RESID_SINK[] = resid_sink    # redirect residual-region integers to range coder
        end
        if write_enc == :interval
            if lr_split
                _write_intervals_lr(wr, target_list, ie, write_mil; vertex_id=_resid_vid, tight_intervals=tight_intervals)
            else
                write_intervals_and_residuals(wr, target_list, ie, write_mil; vertex_id=_resid_vid, tight_intervals=tight_intervals)
            end
        elseif write_enc == :delta
            # Determine delta format: adaptive (per-vertex) or global
            vertex_stop = adaptive_deltas ? use_stop : stop_deltas
            if adaptive_deltas
                _struct_write_bit(wr, vertex_stop)  # 1-bit flag: true=stop, false=count
            end
            if vertex_stop
                _write_stop_delta(wr, target_list, ie; vertex_id=_resid_vid)
            else
                # Self-delimiting delta: count + values
                write_encoded_value(wr, T(length(target_list) + 1), ie)
                if !isempty(target_list)
                    write_delta(wr, target_list, ie; vertex_id=_resid_vid)
                end
            end
        elseif write_enc == :rle
            # Self-delimiting RLE: count + hybrid
            write_encoded_value(wr, T(length(target_list) + 1), ie)
            if !isempty(target_list)
                deltas = delta_encode_vector(target_list)
                write_hybrid_mix_encoded_list(wr, deltas, ie, true, write_mil, false; vertex_id=_resid_vid)
            end
        end
        if resid_sink !== nothing
            _RESID_SINK[] = nothing       # end residual region
        end
        if _STREAM_BITS[] !== nothing
            _VALUE_TAP[] = nothing
            _sbres = _STREAM_BITS[]::Dict{Symbol,Int}
            _sbres[:residual] = get(_sbres, :residual, 0) + (_total_bits(wr) - _rs0)
        end
        if _RESID_META[] !== nothing && _RESID_ACC[] !== nothing
            local _rmc = length(_RESID_ACC[]::Vector{UInt64}) - _rsn0
            local _rmm = (actual_ref_mode == :multi_ref) ? 2 : (actual_ref_mode == :reference ? 1 : 0)
            push!(_RESID_META[]::Vector{NTuple{3,Int}}, (length(current_neighbors), _rmc, _rmm))
        end
        let t = _RESIDIDS_TAP[]
            t === nothing || push!(t::Vector{Any}, copy(target_list))
        end

        add_to_ref_window!(v)
    end

    # Per-chunk range blobs (ctx_chunked only): concatenated finalized chunks +
    # cumulative byte-offset tables (length nchunks+1, first entry 0, last = blob len).
    resid_blob = UInt8[]; refdist_blob = UInt8[]; copy_blob = UInt8[]
    cmd_blob = UInt8[]; flag_blob = UInt8[]
    resid_coff = Int[0]; refdist_coff = Int[0]; copy_coff = Int[0]
    cmd_coff = Int[0]; flag_coff = Int[0]
    function _finalize_chunk!()
        append!(resid_blob,   rc_finish_and_reset!(resid_sink));   push!(resid_coff,   length(resid_blob))
        append!(refdist_blob, rc_finish_and_reset!(refdist_sink)); push!(refdist_coff, length(refdist_blob))
        append!(copy_blob,    brc_finish_and_reset!(copy_sink));   push!(copy_coff,    length(copy_blob))
        append!(cmd_blob,     rc_finish_and_reset!(cmd_sink));     push!(cmd_coff,     length(cmd_blob))
        append!(flag_blob,    brc_finish_and_reset!(flag_sink));   push!(flag_coff,    length(flag_blob))
        return nothing
    end
    if idx_mode
        # Two-pass encoding: encode per-vertex data to buffer, then write offset table + data
        buf_bw = BitWriter()
        vertex_offsets = Vector{Int}(undef, vs + 1)

        flag_sink === nothing || (_FLAG_SINK[] = flag_sink)
        try
            for v_idx in 1:vs
                # At each chunk boundary (aligned to the sampled-index sample points),
                # finalize + reset the range coders so each chunk is self-contained.
                if ctx_chunked && v_idx > 1 && (v_idx - 1) % index_sample_k == 0
                    _finalize_chunk!()
                end
                v = T(v_idx)
                vertex_offsets[v_idx] = _total_bits(buf_bw)
                current_neighbors = sort(get(neighbor_lists, v, T[]))
                _encode_vertex!(buf_bw, v, current_neighbors)
            end
        finally
            _FLAG_SINK[] = nothing
        end
        if ctx_chunked
            # Finalize the last chunk.
            _finalize_chunk!()
        end
        flush_bitwriter(buf_bw; flush_last_bits=true)
        vertex_offsets[vs + 1] = bytes_written(buf_bw) * 8

        # Compute entry width from max offset
        max_offset = vertex_offsets[vs + 1]
        entry_width = max_offset > 0 ? max(Int(ceil(log2(max_offset + 1))), 1) : 1

        if index_sample_k > 0
            # Two-level sampled offsets:
            # [1 bit: sampled=1] [8 bits: sample_k] [6 bits: entry_width]
            # [ceil(n/k)+2 entries × entry_width bits]
            # Data is identical to full index mode.
            write_bit(w, true)  # sampled flag
            @assert index_sample_k >= 4 && index_sample_k % 4 == 0 "index_sample_k must be a multiple of 4, >= 4"
            if use_ctx_range
                # :context_range (new path): store raw k in 32 bits so chunk sizes
                # can exceed the 8-bit (k/4-1) field's 1024 ceiling.
                write_value(w, UInt64(index_sample_k), 32)
            else
                # Legacy fibonacci sampled index: k as (k/4-1) in 8 bits (k ≤ 1024).
                write_value(w, UInt64(index_sample_k ÷ 4 - 1), 8)
            end
            write_value(w, UInt64(entry_width), 6)

            # Write sampled offsets: one entry per k-th vertex + end marker
            for v_idx in 1:index_sample_k:vs
                write_value(w, UInt64(vertex_offsets[v_idx]), entry_width)
            end
            # End marker: total data size
            write_value(w, UInt64(vertex_offsets[vs + 1]), entry_width)

            # :context_range — per-chunk byte-offset tables for the five range
            # streams (resid, refdist, copy, cmd, flag), each length nchunks+1
            # (cumulative, last = blob length). Lets a reader slice + spin up fresh
            # decoders on the chunk containing v.
            if ctx_chunked
                max_boff = max(resid_coff[end], refdist_coff[end], copy_coff[end],
                               cmd_coff[end], flag_coff[end])
                coff_width = max_boff > 0 ? max(Int(ceil(log2(max_boff + 1))), 1) : 1
                write_value(w, UInt64(coff_width), 6)
                for offs in (resid_coff, refdist_coff, copy_coff, cmd_coff, flag_coff)
                    for o in offs
                        write_value(w, UInt64(o), coff_width)
                    end
                end
            end
        else
            # Full offset table (original behavior):
            # [1 bit: sampled=0] [6 bits: entry_width] [(N+1) entries × entry_width bits]
            write_bit(w, false)  # not sampled
            write_value(w, UInt64(entry_width), 6)
            for i in 1:(vs + 1)
                write_value(w, UInt64(vertex_offsets[i]), entry_width)
            end
        end

        # Write buffered encoded data
        buf_data = collect(get_bytes(buf_bw))
        write_bytes(w, buf_data)
    else
        # Children mode: legacy degree-count index table (not used) + sequential encoding
        # Unified Data section. Gate the flag sink for the whole loop (2b).
        flag_sink === nothing || (_FLAG_SINK[] = flag_sink)
        try
            for v_idx in 1:vs
                v = T(v_idx)
                current_neighbors = sort(get(neighbor_lists, v, T[]))
                _encode_vertex!(w, v, current_neighbors)
            end
        finally
            _FLAG_SINK[] = nothing
        end
    end
    # :context_range — finalize the five range streams (all empty otherwise)
    if ctx_chunked
        # Already finalized per chunk; return the concatenated blobs.
        return (resid_blob, refdist_blob, copy_blob, cmd_blob, flag_blob)
    end
    if use_ctx_range
        return (rc_finish!(resid_sink), rc_finish!(refdist_sink), brc_finish!(copy_sink),
                rc_finish!(cmd_sink), brc_finish!(flag_sink))
    end
    return (UInt8[], UInt8[], UInt8[], UInt8[], UInt8[])
end

function read_greedy_graph_data(r::BitReader, vs::T, coding_scheme::Symbol, ::Type{T};
        integer_encoding::Symbol=:fibonacci, copy_blocks::Bool=false, adaptive_copy::Bool=false,
        ref_window_size::Int=7, fixwidth_ref::Bool=false, stop_deltas::Bool=false,
        adaptive_deltas::Bool=false,
        split_residual::Bool=false, bv_blocks::Bool=false, compact_copy::Bool=false,
        tight_intervals::Bool=false,
        lr_split::Bool=false,
        multi_ref::Bool=false,
        adaptive_header::Bool=false,
        ctx_range::Bool=false,
        resid_bytes::Vector{UInt8}=UInt8[],
        refdist_bytes::Vector{UInt8}=UInt8[],
        copy_bytes::Vector{UInt8}=UInt8[],
        cmd_bytes::Union{Nothing,Vector{UInt8}}=nothing,
        flag_bytes::Union{Nothing,Vector{UInt8}}=nothing,
        decode_range::Union{Nothing,UnitRange{Int}}=nothing,
        seed_window::Vector{T}=T[],
        ref_resolver::Union{Nothing,Function}=nothing) where {T<:Unsigned}

    neighbor_lists = Dict{T,Vector{T}}()
    ref_dist_bits = fixwidth_ref ? ceil(Int, log2(max(2, ref_window_size))) : 0
    # :context_range — residual/refdist/copy streams are pulled from range decoders.
    # For the chunked (sampled-index) path these are (re)initialised per chunk on a
    # byte slice; here we set up the whole-blob decoders used by the non-index path.
    resid_source   = ctx_range ? CtxRangeDecoder(resid_bytes) : nothing
    refdist_source = ctx_range ? CtxRangeDecoder(refdist_bytes) : nothing
    copy_source    = ctx_range ? BinRangeDecoder(copy_bytes) : nothing
    # Phase 2b: command + flag streams. nothing = non-ctx / fibonacci files.
    cmd_source  = (ctx_range && cmd_bytes  !== nothing) ? CtxRangeDecoder(cmd_bytes)  : nothing
    flag_source = (ctx_range && flag_bytes !== nothing) ? BinRangeDecoder(flag_bytes) : nothing
    # Capture the caller's flag-source global NOW — before any _setup_chunk! call
    # mutates it — so nested RA chunk decodes (via ref_resolver) restore the OUTER
    # call's decoder on exit, not their own.
    _flag_prev = _FLAG_SOURCE[]
    # Read `n` copy bits from the binary range decoder into a Bool bitmap.
    _read_copy_bits(n::Int) = Bool[brc_decode_bit!(copy_source::BinRangeDecoder) for _ in 1:n]
    # Referenced-vertex adjacency lookup: local (already-decoded) first, then an
    # optional external resolver (random-access: materialise a referenced vertex
    # that lives in an earlier, not-yet-decoded chunk).
    _ref_nbs(ref_v::T) = haskey(neighbor_lists, ref_v) ? neighbor_lists[ref_v] :
        (ref_resolver !== nothing ? ref_resolver(ref_v)::Vector{T} : T[])

    # Read stream metadata
    is_index_mode = read_bit(r)
    enc_tag = _read_encoding_tag(r)
    ie = get(TAG_ENCODINGS, enc_tag, :fibonacci)

    # Read offset table or legacy degree array
    vertex_offsets = Int[]
    sampled_offsets = Int[]
    is_sampled_index = false
    samp_k = 0
    # Per-chunk range-blob byte offsets (ctx_range + sampled index).
    resid_coff = Int[]; refdist_coff = Int[]; copy_coff = Int[]
    cmd_coff = Int[]; flag_coff = Int[]
    data_start_bit = 0
    if is_index_mode
        is_sampled_index = read_bit(r)
        if is_sampled_index
            # Two-level sampled offsets: data uses children-mode encoding.
            # :context_range stores raw k in 32 bits; fibonacci uses (k/4-1) in 8 bits.
            local _sk = ctx_range ? Int(read_value(r, 32, UInt64)) :
                (Int(read_value(r, 8, UInt64)) + 1) * 4
            samp_k = _sk
            entry_width = Int(read_value(r, 6, UInt64))
            n_sampled = length(1:_sk:Int(vs))
            for _ in 1:(n_sampled + 1)
                push!(sampled_offsets, Int(read_value(r, entry_width, UInt64)))
            end
            # :context_range — per-chunk byte-offset tables for the five range
            # streams (each length nchunks+1). Enables seeking to the chunk of v.
            if ctx_range
                coff_width = Int(read_value(r, 6, UInt64))
                for tbl in (resid_coff, refdist_coff, copy_coff, cmd_coff, flag_coff)
                    for _ in 1:(n_sampled + 1)
                        push!(tbl, Int(read_value(r, coff_width, UInt64)))
                    end
                end
            end
            # After the offset tables, the struct reader sits at the first vertex's
            # data — record the absolute bit position so we can seek per chunk.
            data_start_bit = r.bit_count
            # Sampled mode: data is children-mode, keep stop_deltas as passed
        else
            # Full offset table: data uses count-prefixed encoding
            entry_width = Int(read_value(r, 6, UInt64))
            for _ in 1:(Int(vs) + 1)
                push!(vertex_offsets, Int(read_value(r, entry_width, UInt64)))
            end
            # Full index mode: override stop_deltas (count-prefixed, no STOP needed)
            stop_deltas = false
            adaptive_deltas = false
        end
    end

    reference_window = T[]
    function add_to_ref_window!(vertex::T)
        push!(reference_window, vertex)
        if length(reference_window) > ref_window_size
            popfirst!(reference_window)
        end
    end

    # Chunked :context_range (sampled index): the three range streams are stored as
    # a concatenation of independently-finalized chunks (aligned to sample points,
    # every samp_k vertices). Spin up fresh decoders on chunk c's byte slice with
    # fresh adaptive models — matching the encoder's per-chunk reset.
    ctx_chunked = ctx_range && is_sampled_index && !isempty(resid_coff)
    function _setup_chunk!(c::Int)
        resid_source   = CtxRangeDecoder(resid_bytes[resid_coff[c] + 1 : resid_coff[c + 1]])
        refdist_source = CtxRangeDecoder(refdist_bytes[refdist_coff[c] + 1 : refdist_coff[c + 1]])
        copy_source    = BinRangeDecoder(copy_bytes[copy_coff[c] + 1 : copy_coff[c + 1]])
        if cmd_bytes !== nothing && !isempty(cmd_coff)
            cmd_source  = CtxRangeDecoder(cmd_bytes[cmd_coff[c] + 1 : cmd_coff[c + 1]])
            flag_source = BinRangeDecoder(flag_bytes[flag_coff[c] + 1 : flag_coff[c + 1]])
            _FLAG_SOURCE[] = flag_source   # refresh the global gate to the new chunk decoder
        end
        return nothing
    end

    # Random-access seek mode: decode only `decode_range`, seeding the reference
    # window with the caller-supplied predecessor ids and jumping the struct reader
    # + range decoders to the chunk containing the range start.
    if decode_range !== nothing
        for u in seed_window
            push!(reference_window, u)
        end
        c0 = ((first(decode_range) - 1) ÷ samp_k) + 1
        bitpos = data_start_bit + sampled_offsets[c0]
        r.index = bitpos + 1
        r.bit_count = bitpos
        _setup_chunk!(c0)
    end

    # Helper to read encoded list based on type (using captured ie, tight_intervals, etc.)
    # For :context_range this is the residual region: gate the range decoder so
    # read_encoded_value pulls from it (mirror of the encoder's sink region).
    function _read_encoded_list_body(r_in::BitReader, v::T, enc_type::Symbol, mil::Int)
        if resid_source !== nothing
            rc_reset_region!(resid_source::CtxRangeDecoder)
            _RESID_SOURCE[] = resid_source
        end
        try
            if enc_type == :interval
                if lr_split
                    return _read_intervals_lr(r_in, ie, mil, T; vertex_id=v, tight_intervals=tight_intervals)
                else
                    return read_intervals_and_residuals(r_in, ie, mil, T; vertex_id=v, tight_intervals=tight_intervals)
                end
            elseif enc_type == :delta
                # Determine delta format: adaptive (per-vertex) or global
                vertex_stop = adaptive_deltas ? _struct_read_bit(r_in) : stop_deltas
                if vertex_stop
                    return _read_stop_delta(r_in, ie, T; vertex_id=v)
                else
                    count = read_encoded_value(r_in, ie, T) - T(1)
                    if count > 0
                        return read_delta(r_in, ie, T; max_elements=Int(count), vertex_id=v)
                    end
                    return T[]
                end
            elseif enc_type == :rle
                count = read_encoded_value(r_in, ie, T) - T(1)
                if count > 0
                    return read_hybrid_mix_encoded_list(r_in, :index, ie, T; max_elements=Int(count), vertex_id=v)
                end
                return T[]
            end
            return T[]
        finally
            resid_source === nothing || (_RESID_SOURCE[] = nothing)
        end
    end

    # Helper: read one reference distance + copy bitmap, return copied neighbors
    function _read_one_ref(r_in::BitReader)
        if refdist_source !== nothing
            distance = T(rc_decode_value!(refdist_source::CtxRangeDecoder))
        elseif fixwidth_ref
            local d = 0
            for _ in 1:ref_dist_bits
                d = (d << 1) | Int(read_bit(r_in))
            end
            distance = T(d + 1)
        else
            distance = read_encoded_value(r_in, ie, T)
        end
        ref_idx = length(reference_window) - Int(distance) + 1
        ref_nbs = T[]
        if 1 <= ref_idx <= length(reference_window)
            ref_v = reference_window[ref_idx]
            ref_nbs = _ref_nbs(ref_v)
        end
        if copy_source !== nothing
            copy_bitmap = _read_copy_bits(length(ref_nbs))
            copied = reconstruct_from_reference(ref_nbs, copy_bitmap, T[])
        elseif adaptive_copy && copy_blocks
            ref_len = length(ref_nbs)
            copy_bitmap = _read_adaptive_copy(r_in, ref_len, ie, T; bv_blocks=bv_blocks, compact_copy=compact_copy)
            copied = reconstruct_from_reference(ref_nbs, copy_bitmap, T[])
        elseif copy_blocks
            copy_positions = _read_copy_blocks(r_in, ie, T)
            copied = reconstruct_from_copy_positions(ref_nbs, copy_positions, T[])
        else
            copy_bitmap = read_bitmap_adaptive(r_in, ie)
            copied = reconstruct_from_reference(ref_nbs, copy_bitmap, T[])
        end
        return copied
    end

    # Split-residual override tag: from the cmd stream when 2b is active, else raw.
    function _read_split_enc_opt(r_in::BitReader)
        if cmd_source !== nothing
            return ENCODING_OPTIONS[Int(rc_decode_value!(cmd_source::CtxRangeDecoder)) + 1]
        end
        return _read_enc_opt_tag(r_in)
    end

    function _decode_vertex!(r_in::BitReader, v::T)
        # Read vertex header
        ref_mode, enc_type, mil = cmd_source !== nothing ?
            _bg_cmd_decode(Int(rc_decode_value!(cmd_source::CtxRangeDecoder)); adaptive=adaptive_header) :
            _read_vertex_header_vlc(r_in; adaptive=adaptive_header)

        # Empty vertex (merged VLC code)
        if ref_mode == :empty
            return T[]
        end

        current_neighbors = T[]

        if ref_mode == :multi_ref
            # Multi-reference: read two ref+copy pairs, then residuals
            copied1 = _read_one_ref(r_in)
            copied2 = _read_one_ref(r_in)
            # Copy-aware rank gaps: C = copied1 ∪ copied2; decode Q in rank space
            # (rank-mapped vertex id) then invert Q -> R.
            _mr_vid = v; _mr_C = T[]
            if resid_source !== nothing
                _mr_C = sort!(unique!(vcat(copied1, copied2)))
                isempty(_mr_C) || (_mr_vid = _rank_of(v, _mr_C))
            end
            residuals = _read_encoded_list_body(r_in, _mr_vid, enc_type, mil)
            if resid_source !== nothing && !isempty(_mr_C)
                residuals = _rankgap_inverse(residuals, _mr_C)
            end
            # Union all three sets
            current_neighbors = sort!(unique!(vcat(copied1, copied2, residuals)))
        elseif ref_mode == :reference
            if refdist_source !== nothing
                distance = T(rc_decode_value!(refdist_source::CtxRangeDecoder))
            elseif fixwidth_ref
                local d = 0
                for _ in 1:ref_dist_bits
                    d = (d << 1) | Int(read_bit(r_in))
                end
                distance = T(d + 1)
            else
                distance = read_encoded_value(r_in, ie, T)
            end
            if copy_source !== nothing
                ref_idx = length(reference_window) - Int(distance) + 1
                ref_nbs = T[]
                if 1 <= ref_idx <= length(reference_window)
                    ref_v = reference_window[ref_idx]
                    ref_nbs = _ref_nbs(ref_v)
                end
                copy_bitmap = _read_copy_bits(length(ref_nbs))
                res_et, res_m = enc_type, mil
                if split_residual
                    if _struct_read_bit(r_in)
                        res_et, res_m = _read_split_enc_opt(r_in)
                    end
                end
                # Copy-aware rank gaps: C = copied ids; decode Q in rank space then
                # invert Q -> R.
                _rg_vid = v; _rg_C = T[]
                if resid_source !== nothing
                    _rg_C = reconstruct_from_reference(ref_nbs, copy_bitmap, T[])
                    isempty(_rg_C) || (_rg_vid = _rank_of(v, _rg_C))
                end
                residuals = _read_encoded_list_body(r_in, _rg_vid, res_et, res_m)
                if resid_source !== nothing && !isempty(_rg_C)
                    residuals = _rankgap_inverse(residuals, _rg_C)
                end
                current_neighbors = reconstruct_from_reference(ref_nbs, copy_bitmap, residuals)
            elseif adaptive_copy && copy_blocks
                ref_idx = length(reference_window) - Int(distance) + 1
                ref_nbs = T[]
                if 1 <= ref_idx <= length(reference_window)
                    ref_v = reference_window[ref_idx]
                    ref_nbs = _ref_nbs(ref_v)
                end
                ref_len = length(ref_nbs)
                copy_bitmap = _read_adaptive_copy(r_in, ref_len, ie, T; bv_blocks=bv_blocks, compact_copy=compact_copy)
                res_et, res_m = enc_type, mil
                if split_residual
                    if _struct_read_bit(r_in)
                        res_et, res_m = _read_split_enc_opt(r_in)
                    end
                end
                residuals = _read_encoded_list_body(r_in, v, res_et, res_m)
                current_neighbors = reconstruct_from_reference(ref_nbs, copy_bitmap, residuals)
            elseif copy_blocks
                copy_positions = _read_copy_blocks(r_in, ie, T)
                res_et, res_m = enc_type, mil
                if split_residual
                    if _struct_read_bit(r_in)
                        res_et, res_m = _read_split_enc_opt(r_in)
                    end
                end
                residuals = _read_encoded_list_body(r_in, v, res_et, res_m)
                ref_idx = length(reference_window) - Int(distance) + 1
                if 1 <= ref_idx <= length(reference_window)
                    ref_v = reference_window[ref_idx]
                    ref_nbs = _ref_nbs(ref_v)
                    current_neighbors = reconstruct_from_copy_positions(ref_nbs, copy_positions, residuals)
                else
                    current_neighbors = residuals
                end
            else
                copy_bitmap = read_bitmap_adaptive(r_in, ie)
                res_et, res_m = enc_type, mil
                if split_residual
                    if _struct_read_bit(r_in)
                        res_et, res_m = _read_split_enc_opt(r_in)
                    end
                end
                residuals = _read_encoded_list_body(r_in, v, res_et, res_m)
                ref_idx = length(reference_window) - Int(distance) + 1
                if 1 <= ref_idx <= length(reference_window)
                    ref_v = reference_window[ref_idx]
                    ref_nbs = _ref_nbs(ref_v)
                    current_neighbors = reconstruct_from_reference(ref_nbs, copy_bitmap, residuals)
                else
                    current_neighbors = residuals
                end
            end
        else
            current_neighbors = _read_encoded_list_body(r_in, v, enc_type, mil)
        end

        return current_neighbors
    end

    # Unified Data loop. In seek mode only `decode_range` is visited; otherwise the
    # whole graph is decoded sequentially. (_flag_prev captured at function entry;
    # restored in the finally so nested RA decodes hand back the outer decoder.)
    loop_range = decode_range === nothing ? (1:Int(vs)) : decode_range
    flag_source === nothing || (_FLAG_SOURCE[] = flag_source)
    try
        for v_idx in loop_range
            v = T(v_idx)

            # Chunked :context_range full decode: swap to fresh per-chunk decoders at
            # each sample boundary. (In seek mode the chunk is already set up and the
            # range never crosses a boundary.)
            if ctx_chunked && decode_range === nothing && (v_idx - 1) % samp_k == 0
                _setup_chunk!(((v_idx - 1) ÷ samp_k) + 1)
            end

            if is_index_mode && !isempty(vertex_offsets)
                # Full index mode: empty vertex = offset[v] == offset[v+1]
                if vertex_offsets[v_idx] == vertex_offsets[v_idx + 1]
                    neighbor_lists[v] = T[]
                    add_to_ref_window!(v)
                    continue
                end
            end
            # Sampled index mode: decode normally (count-prefixed handles empty vertices via count=0+1=1 → 0 elements)

            current_neighbors = _decode_vertex!(r, v)

            neighbor_lists[v] = current_neighbors
            add_to_ref_window!(v)
        end
    finally
        _FLAG_SOURCE[] = _flag_prev
    end

    return neighbor_lists
end


# Parse the header + offset tables of a sampled-index :context_range BG struct
# stream, returning everything a random-access reader needs to seek to a chunk.
# Decodes NO vertex data. `vs` = number of vertices.
function parse_bg_ctxrange_index(struct_bytes::Vector{UInt8}, vs::Int)
    r = BitReader(IOBuffer(struct_bytes))
    is_index = read_bit(r)
    enc_tag = _read_encoding_tag(r)
    ie = get(TAG_ENCODINGS, enc_tag, :fibonacci)
    is_index || error("parse_bg_ctxrange_index: not an index-mode stream")
    is_sampled = read_bit(r)
    is_sampled || error("parse_bg_ctxrange_index: random access requires a sampled index")
    samp_k = Int(read_value(r, 32, UInt64))   # :context_range stores raw k in 32 bits
    entry_width = Int(read_value(r, 6, UInt64))
    n_sampled = length(1:samp_k:vs)
    sampled_offsets = Int[]
    for _ in 1:(n_sampled + 1)
        push!(sampled_offsets, Int(read_value(r, entry_width, UInt64)))
    end
    coff_width = Int(read_value(r, 6, UInt64))
    resid_coff = Int[]; refdist_coff = Int[]; copy_coff = Int[]
    cmd_coff = Int[]; flag_coff = Int[]
    for tbl in (resid_coff, refdist_coff, copy_coff, cmd_coff, flag_coff)
        for _ in 1:(n_sampled + 1)
            push!(tbl, Int(read_value(r, coff_width, UInt64)))
        end
    end
    data_start_bit = r.bit_count
    return (samp_k=samp_k, sampled_offsets=sampled_offsets, resid_coff=resid_coff,
        refdist_coff=refdist_coff, copy_coff=copy_coff, cmd_coff=cmd_coff,
        flag_coff=flag_coff, data_start_bit=data_start_bit,
        ie=ie, n_sampled=n_sampled)
end

function _estimate_adaptive_bitmap_cost(bitmap::Vector{Bool}, ie::Symbol)::Int
    # Same logic as write_bitmap_adaptive: min of raw vs block encoding + 1-bit flag
    raw_cost = estimate_encoded_value_cost(UInt32(length(bitmap)) + UInt32(1), ie) + length(bitmap)
    block_cost = estimate_block_encoding_cost(bitmap, ie)
    return 1 + min(raw_cost, block_cost)  # 1-bit format flag
end

# --- Copy-blocks for greedy path ---
# More compact than adaptive bitmap: no bitmap length, no flag bit, no skip block lengths.
# Format: small_count(num_copy_blocks), then for each block:
#   first: varint(start) + varint(length)
#   subsequent: varint(gap_from_prev_end) + varint(length)

function _bitmap_to_copy_blocks(bitmap::Vector{Bool})
    positions = Int[i for i in 1:length(bitmap) if bitmap[i]]
    if isempty(positions)
        return Tuple{Int,Int}[]
    end
    blocks = Tuple{Int,Int}[]
    i = 1
    while i <= length(positions)
        start = positions[i]
        len = 1
        while i + len <= length(positions) && positions[i + len] == start + len
            len += 1
        end
        push!(blocks, (start, len))
        i += len
    end
    return blocks
end

function _write_copy_blocks(w::BitWriter, bitmap::Vector{Bool}, encoding::Symbol)
    blocks = _bitmap_to_copy_blocks(bitmap)
    ncb = length(blocks)
    write_small_count(w, UInt32(ncb), encoding)
    ncb == 0 && return
    write_encoded_value(w, UInt32(blocks[1][1]), encoding)
    write_encoded_value(w, UInt32(blocks[1][2]), encoding)
    prev_end = blocks[1][1] + blocks[1][2]
    for i in 2:ncb
        gap = blocks[i][1] - prev_end
        write_encoded_value(w, UInt32(gap), encoding)
        write_encoded_value(w, UInt32(blocks[i][2]), encoding)
        prev_end = blocks[i][1] + blocks[i][2]
    end
end

function _read_copy_blocks(r::BitReader, encoding::Symbol, ::Type{T}) where {T<:Unsigned}
    ncb = Int(read_small_count(r, encoding, T))
    ncb == 0 && return Int[]
    positions = Int[]
    start = Int(read_encoded_value(r, encoding, T))
    len = Int(read_encoded_value(r, encoding, T))
    for j in 0:(len-1); push!(positions, start + j); end
    prev_end = start + len
    for i in 2:ncb
        gap = Int(read_encoded_value(r, encoding, T))
        start = prev_end + gap
        len = Int(read_encoded_value(r, encoding, T))
        for j in 0:(len-1); push!(positions, start + j); end
        prev_end = start + len
    end
    return positions
end

@inline function _small_count_cost(v::UInt32, encoding::Symbol)::Int
    v <= 2 && return 2
    return 2 + estimate_encoded_value_cost(v, encoding)
end

function _estimate_copy_blocks_cost(bitmap::Vector{Bool}, encoding::Symbol)::Int
    blocks = _bitmap_to_copy_blocks(bitmap)
    ncb = length(blocks)
    cost = _small_count_cost(UInt32(ncb), encoding)
    ncb == 0 && return cost
    cost += estimate_encoded_value_cost(UInt32(blocks[1][1]), encoding)
    cost += estimate_encoded_value_cost(UInt32(blocks[1][2]), encoding)
    prev_end = blocks[1][1] + blocks[1][2]
    for i in 2:ncb
        gap = blocks[i][1] - prev_end
        cost += estimate_encoded_value_cost(UInt32(gap), encoding)
        cost += estimate_encoded_value_cost(UInt32(blocks[i][2]), encoding)
        prev_end = blocks[i][1] + blocks[i][2]
    end
    return cost
end

# --- BV-style alternating copy/skip block encoding ---
# Format: block_count, then alternating block sizes: copy_len, skip_len, copy_len, ...
# First block is always "copy" (can be 0 if bitmap starts with skip).
# More compact than (start, length) copy-blocks because starts are implicit.

function _bitmap_to_bv_blocks(bitmap::Vector{Bool})::Vector{Int}
    isempty(bitmap) && return Int[0]
    blocks = Int[]
    i = 1
    expecting_copy = true
    while i <= length(bitmap)
        run = 0
        while i <= length(bitmap) && bitmap[i] == expecting_copy
            run += 1; i += 1
        end
        push!(blocks, run)
        expecting_copy = !expecting_copy
    end
    # First block must be "copy" — if bitmap starts with skip, insert 0
    if !bitmap[1]
        pushfirst!(blocks, 0)
    end
    return blocks
end

function _write_bv_blocks(w::BitWriter, bitmap::Vector{Bool}, encoding::Symbol)
    blocks = _bitmap_to_bv_blocks(bitmap)
    write_small_count(w, UInt32(length(blocks)), encoding)
    for b in blocks
        write_encoded_value(w, UInt32(b), encoding)
    end
end

function _read_bv_blocks(r::BitReader, encoding::Symbol, ::Type{T}, ref_len::Int) where {T<:Unsigned}
    nblocks = Int(read_small_count(r, encoding, T))
    nblocks == 0 && return Bool[false for _ in 1:ref_len]
    bitmap = Bool[]
    is_copy = true
    for bi in 1:nblocks
        blen = Int(read_encoded_value(r, encoding, T))
        for _ in 1:blen
            push!(bitmap, is_copy)
        end
        is_copy = !is_copy
    end
    # Pad to ref_len if needed (trailing elements are !is_copy at this point)
    while length(bitmap) < ref_len
        push!(bitmap, is_copy)
    end
    return bitmap[1:ref_len]
end

function _estimate_bv_blocks_cost(bitmap::Vector{Bool}, encoding::Symbol)::Int
    blocks = _bitmap_to_bv_blocks(bitmap)
    cost = _small_count_cost(UInt32(length(blocks)), encoding)
    for b in blocks
        cost += estimate_encoded_value_cost(UInt32(b), encoding)
    end
    return cost
end

# --- Adaptive copy: bitmap / copy-blocks or BV-blocks / complement ---
# Nested-bit scheme:
# Default prefix code (compact_copy=false):
#   outer=1          → raw bitmap (ref_len bits)
#   outer=0, inner=0 → copy-blocks or BV-blocks (runs of copied positions)
#   outer=0, inner=1 → complement copy-blocks or BV-blocks (runs of SKIPPED positions)
# Costs: bitmap=1+ref_len, blocks=2+cb_bits, complement=2+cc_bits
#
# Compact prefix code (compact_copy=true):
#   0       → complement blocks (1 bit overhead, best for high-overlap)
#   10      → copy-blocks (2 bits overhead)
#   11      → raw bitmap (2 bits overhead)

function _estimate_adaptive_copy_cost(bitmap::Vector{Bool}, ie::Symbol; bv_blocks::Bool=false, compact_copy::Bool=false)::Int
    ref_len = length(bitmap)
    complement = Bool[!b for b in bitmap]
    if bv_blocks
        cb_raw = _estimate_bv_blocks_cost(bitmap, ie)
        cc_raw = _estimate_bv_blocks_cost(complement, ie)
    else
        cb_raw = _estimate_copy_blocks_cost(bitmap, ie)
        cc_raw = _estimate_copy_blocks_cost(complement, ie)
    end
    if compact_copy
        # 0=complement(1b), 10=blocks(2b), 11=bitmap(2b)
        bm_cost = 2 + ref_len
        cb_cost = 2 + cb_raw
        cc_cost = 1 + cc_raw
    else
        # 1=bitmap(1b), 00=blocks(2b), 01=complement(2b)
        bm_cost = 1 + ref_len
        cb_cost = 2 + cb_raw
        cc_cost = 2 + cc_raw
    end
    return min(bm_cost, cb_cost, cc_cost)
end

# Fast copy cost: just bitmap length, no copy-blocks analysis (O(1) instead of O(n))
@inline function _estimate_copy_cost_fast(bitmap::Vector{Bool})::Int
    return 1 + length(bitmap)  # 1-bit mode flag + raw bitmap bits
end

function _write_adaptive_copy(w::BitWriter, bitmap::Vector{Bool}, ie::Symbol; bv_blocks::Bool=false, compact_copy::Bool=false)
    ref_len = length(bitmap)
    complement = Bool[!b for b in bitmap]
    if bv_blocks
        cb_raw = _estimate_bv_blocks_cost(bitmap, ie)
        cc_raw = _estimate_bv_blocks_cost(complement, ie)
    else
        cb_raw = _estimate_copy_blocks_cost(bitmap, ie)
        cc_raw = _estimate_copy_blocks_cost(complement, ie)
    end
    if compact_copy
        bm_cost = 2 + ref_len
        cb_cost = 2 + cb_raw
        cc_cost = 1 + cc_raw
    else
        bm_cost = 1 + ref_len
        cb_cost = 2 + cb_raw
        cc_cost = 2 + cc_raw
    end
    if compact_copy
        if cc_cost <= cb_cost && cc_cost <= bm_cost
            write_bit(w, false)                           # 0: complement
            if bv_blocks
                _write_bv_blocks(w, complement, ie)
            else
                _write_copy_blocks(w, complement, ie)
            end
        elseif cb_cost <= bm_cost
            write_bit(w, true); write_bit(w, false)       # 10: copy blocks
            if bv_blocks
                _write_bv_blocks(w, bitmap, ie)
            else
                _write_copy_blocks(w, bitmap, ie)
            end
        else
            write_bit(w, true); write_bit(w, true)        # 11: bitmap
            for b in bitmap; write_bit(w, b); end
        end
    else
        if bm_cost <= cb_cost && bm_cost <= cc_cost
            write_bit(w, true)                            # outer=1: bitmap
            for b in bitmap; write_bit(w, b); end
        elseif cb_cost <= cc_cost
            write_bit(w, false); write_bit(w, false)      # outer=0, inner=0: blocks
            if bv_blocks
                _write_bv_blocks(w, bitmap, ie)
            else
                _write_copy_blocks(w, bitmap, ie)
            end
        else
            write_bit(w, false); write_bit(w, true)       # outer=0, inner=1: complement
            if bv_blocks
                _write_bv_blocks(w, complement, ie)
            else
                _write_copy_blocks(w, complement, ie)
            end
        end
    end
end

function _read_adaptive_copy(r::BitReader, ref_len::Int, ie::Symbol, ::Type{T}; bv_blocks::Bool=false, compact_copy::Bool=false) where {T<:Unsigned}
    if compact_copy
        if !read_bit(r)  # 0: complement
            if bv_blocks
                compl = _read_bv_blocks(r, ie, T, ref_len)
                return Bool[!b for b in compl]
            else
                skip_positions = _read_copy_blocks(r, ie, T)
                skip_set = Set(skip_positions)
                return Bool[!(i in skip_set) for i in 1:ref_len]
            end
        else
            inner = read_bit(r)
            if !inner  # 10: copy blocks
                if bv_blocks
                    return _read_bv_blocks(r, ie, T, ref_len)
                else
                    positions = _read_copy_blocks(r, ie, T)
                    pos_set = Set(positions)
                    return Bool[i in pos_set for i in 1:ref_len]
                end
            else  # 11: bitmap
                return Bool[read_bit(r) for _ in 1:ref_len]
            end
        end
    else
        if read_bit(r)  # outer=1: bitmap
            return Bool[read_bit(r) for _ in 1:ref_len]
        else
            inner = read_bit(r)
            if !inner   # inner=0: blocks → copied positions
                if bv_blocks
                    return _read_bv_blocks(r, ie, T, ref_len)
                else
                    positions = _read_copy_blocks(r, ie, T)
                    pos_set = Set(positions)
                    return Bool[i in pos_set for i in 1:ref_len]
                end
            else        # inner=1: complement → skipped positions
                if bv_blocks
                    compl = _read_bv_blocks(r, ie, T, ref_len)
                    return Bool[!b for b in compl]
                else
                    skip_positions = _read_copy_blocks(r, ie, T)
                    skip_set = Set(skip_positions)
                    return Bool[!(i in skip_set) for i in 1:ref_len]
                end
            end
        end
    end
end

# --- STOP-terminated delta lists ---
# Format: for each value: bit '1' + varint(gap); then bit '0'.
# First value: zigzag(v₁ - vertex_id) + 1 when vertex_id is provided.
# Remaining values: gap = v[i] - v[i-1] (>= 1 for sorted unique).
# Eliminates the count prefix used by standard delta encoding.

function _write_stop_delta(w::BitWriter, sorted_vals::Vector{T}, ie::Symbol; vertex_id=nothing) where {T<:Unsigned}
    if isempty(sorted_vals)
        _struct_write_bit(w, false)  # STOP
        return
    end
    # First value
    _struct_write_bit(w, true)
    if vertex_id !== nothing
        offset = Int64(sorted_vals[1]) - Int64(vertex_id)
        write_encoded_value(w, UInt64(_zigzag_encode(offset) + 1), ie)
    else
        write_encoded_value(w, sorted_vals[1], ie)
    end
    # Remaining: gaps (>= 1 for sorted unique lists)
    for i in 2:length(sorted_vals)
        _struct_write_bit(w, true)
        write_encoded_value(w, sorted_vals[i] - sorted_vals[i-1], ie)
    end
    _struct_write_bit(w, false)  # STOP
end

function _read_stop_delta(r::BitReader, ie::Symbol, ::Type{T}; vertex_id=nothing) where {T<:Unsigned}
    result = T[]
    first = true
    prev = zero(T)
    while _struct_read_bit(r)  # '1' = more values, '0' = STOP
        if first
            if vertex_id !== nothing
                raw64 = read_encoded_value(r, ie, UInt64)
                val = T(Int64(vertex_id) + _zigzag_decode(raw64 - 1))
            else
                val = read_encoded_value(r, ie, T)
            end
            first = false
        else
            val = prev + read_encoded_value(r, ie, T)
        end
        push!(result, val)
        prev = val
    end
    return result
end

function _estimate_stop_delta_cost(sorted_vals::Vector{T}, ie::Symbol; vertex_id=nothing)::Int where {T<:Unsigned}
    isempty(sorted_vals) && return 1  # just the STOP bit
    cost = 1  # STOP bit
    # First value
    if vertex_id !== nothing
        offset = Int64(sorted_vals[1]) - Int64(vertex_id)
        cost += 1 + estimate_encoded_value_cost(UInt64(_zigzag_encode(offset) + 1), ie)
    else
        cost += 1 + estimate_encoded_value_cost(sorted_vals[1], ie)
    end
    # Remaining gaps
    for i in 2:length(sorted_vals)
        cost += 1 + estimate_encoded_value_cost(sorted_vals[i] - sorted_vals[i-1], ie)
    end
    return cost
end

# --- LR-split interval encoding ---
# Intervals encoded as standard; residuals split at vertex_id into left (< vid) and
# right (>= vid) halves, each encoded as ascending distances from vertex_id.
# Produces smaller first-values than zigzag, saving ~0.1–0.4 BPE on web graphs.

function _write_intervals_lr(w::BitWriter, neighbors::Vector{T}, encoding::Symbol,
                              mil::Int; vertex_id=nothing, tight_intervals::Bool=false) where {T<:Unsigned}
    if vertex_id === nothing || isempty(neighbors)
        write_intervals_and_residuals(w, neighbors, encoding, mil; vertex_id=vertex_id, tight_intervals=tight_intervals)
        return
    end
    vid = T(vertex_id)

    intervals, residuals = compress_intervals(neighbors, mil)

    # Write intervals (same format as standard)
    write_encoded_value(w, T(length(intervals)) + T(1), encoding)
    if !isempty(intervals)
        prev_ref = T(0)
        for (idx, (start, len)) in enumerate(intervals)
            if idx == 1
                offset = Int64(start) - Int64(vid)
                write_encoded_value(w, UInt64(_zigzag_encode(offset) + 1), encoding)
            else
                write_encoded_value(w, start - prev_ref, encoding)
            end
            write_encoded_value(w, len - T(mil) + T(1), encoding)
            prev_ref = tight_intervals ? (start + T(len)) : start
        end
    end

    # Write residual count
    write_encoded_value(w, T(length(residuals)) + T(1), encoding)
    if !isempty(residuals)
        # Split at vertex_id
        split_idx = 0
        @inbounds for i in 1:length(residuals)
            if residuals[i] < vid
                split_idx = i
            else
                break
            end
        end
        n_left = split_idx
        n_right = length(residuals) - n_left

        # Left count (decoder derives right = total - left)
        write_encoded_value(w, T(n_left) + T(1), encoding)

        # Left residuals: reversed to ascending distances from vid
        if n_left > 0
            left_dists = Vector{T}(undef, n_left)
            @inbounds for i in 1:n_left
                left_dists[i] = vid - residuals[n_left - i + 1]
            end
            write_delta(w, left_dists, encoding)
        end

        # Right residuals: distances from vid (+1 to handle val==vid)
        if n_right > 0
            right_dists = Vector{T}(undef, n_right)
            @inbounds for i in 1:n_right
                right_dists[i] = residuals[n_left + i] - vid + T(1)
            end
            write_delta(w, right_dists, encoding)
        end
    end
end

function _read_intervals_lr(r::BitReader, encoding::Symbol, mil::Int,
                             ::Type{T}; vertex_id=nothing, tight_intervals::Bool=false) where {T<:Unsigned}
    if vertex_id === nothing
        return read_intervals_and_residuals(r, encoding, mil, T; vertex_id=vertex_id, tight_intervals=tight_intervals)
    end
    vid = T(vertex_id)

    # Read intervals (same format as standard)
    num_intervals = Int(read_encoded_value(r, encoding, T)) - 1
    neighbors = T[]

    if num_intervals > 0
        prev_ref = T(0)
        for idx in 1:num_intervals
            if idx == 1
                raw_start = read_encoded_value(r, encoding, UInt64)
                start = T(Int64(vid) + _zigzag_decode(raw_start - 1))
            else
                start = prev_ref + read_encoded_value(r, encoding, T)
            end
            len = Int(read_encoded_value(r, encoding, T)) - 1 + mil
            for j in 0:(len-1)
                push!(neighbors, start + T(j))
            end
            prev_ref = tight_intervals ? (start + T(len)) : start
        end
    end

    # Read residuals
    num_residuals = Int(read_encoded_value(r, encoding, T)) - 1
    if num_residuals > 0
        n_left = Int(read_encoded_value(r, encoding, T)) - 1
        n_right = num_residuals - n_left

        if n_left > 0
            left_dists = read_delta(r, encoding, T; max_elements=n_left)
            for i in n_left:-1:1
                push!(neighbors, vid - left_dists[i])
            end
        end

        if n_right > 0
            right_dists = read_delta(r, encoding, T; max_elements=n_right)
            for d in right_dists
                push!(neighbors, vid + d - T(1))
            end
        end
    end

    sort!(neighbors)
    return neighbors
end

function _estimate_intervals_lr_cost(neighbors::Vector{T}, ie::Symbol, mil::Int;
                                      vertex_id=nothing, tight_intervals::Bool=false) where {T<:Unsigned}
    if vertex_id === nothing || isempty(neighbors)
        return estimate_interval_runlength_encoding_cost(neighbors, ie, mil, 3; vertex_id=vertex_id, tight_intervals=tight_intervals)
    end
    vid = T(vertex_id)
    intervals, residuals = compress_intervals(neighbors, mil)
    cost = 0

    # Interval count
    cost += estimate_encoded_value_cost(T(length(intervals)) + T(1), ie)

    # Intervals (same as standard)
    prev_ref = T(0)
    for (idx, (start, len)) in enumerate(intervals)
        if idx == 1
            encoded_start = UInt64(_zigzag_encode(Int64(start) - Int64(vid)) + 1)
            cost += estimate_encoded_value_cost(encoded_start, ie)
        else
            cost += estimate_encoded_value_cost(start - prev_ref, ie)
        end
        cost += estimate_encoded_value_cost(T(len - mil + 1), ie)
        prev_ref = tight_intervals ? (start + T(len)) : start
    end

    # Residual count
    cost += estimate_encoded_value_cost(T(length(residuals)) + T(1), ie)

    if !isempty(residuals)
        # Split at vid
        split_idx = 0
        @inbounds for i in 1:length(residuals)
            residuals[i] < vid ? (split_idx = i) : break
        end
        n_left = split_idx
        n_right = length(residuals) - n_left

        # Left count
        cost += estimate_encoded_value_cost(T(n_left) + T(1), ie)

        # Left distances (ascending from vid)
        if n_left > 0
            prev = T(0)
            for i in 1:n_left
                d = vid - residuals[n_left - i + 1]
                gap = i == 1 ? d : (d - prev + T(1))
                cost += estimate_encoded_value_cost(gap, ie)
                prev = d
            end
        end

        # Right distances (ascending from vid, +1 shift)
        if n_right > 0
            prev = T(0)
            for i in 1:n_right
                d = residuals[n_left + i] - vid + T(1)
                gap = i == 1 ? d : (d - prev + T(1))
                cost += estimate_encoded_value_cost(gap, ie)
                prev = d
            end
        end
    end

    return cost
end

function reconstruct_from_copy_positions(reference::Vector{T}, positions::Vector{Int}, residuals::Vector{T}) where {T<:Unsigned}
    selected = T[]
    @inbounds for p in positions
        if 1 <= p <= length(reference)
            push!(selected, reference[p])
        end
    end
    isempty(selected) && return copy(residuals)
    isempty(residuals) && return selected
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
    @inbounds while i <= n1; result[k] = selected[i]; i += 1; k += 1; end
    @inbounds while j <= n2; result[k] = residuals[j]; j += 1; k += 1; end
    return result
end

function _estimate_delta_cost(neighbors::Vector{T}, ie::Symbol; vertex_id=nothing) where {T<:Unsigned}
    isempty(neighbors) && return 0
    cost = 0
    if vertex_id !== nothing
        encoded_first = UInt64(_zigzag_encode(Int64(neighbors[1]) - Int64(vertex_id)) + 1)
        cost += estimate_encoded_value_cost(encoded_first, ie)
    else
        cost += estimate_encoded_value_cost(neighbors[1], ie)
    end
    @inbounds for i in 2:length(neighbors)
        cost += estimate_encoded_value_cost((neighbors[i] - neighbors[i-1]) + T(1), ie)
    end
    return cost
end

function _estimate_rle_cost(neighbors::Vector{T}, ie::Symbol, mil::Int; vertex_id=nothing) where {T<:Unsigned}
    isempty(neighbors) && return 0
    cost = 1 # hybrid flag
    if vertex_id !== nothing
        encoded_first = UInt64(_zigzag_encode(Int64(neighbors[1]) - Int64(vertex_id)) + 1)
        cost += estimate_encoded_value_cost(encoded_first, ie)
    else
        cost += estimate_encoded_value_cost(neighbors[1], ie)
    end

    n = length(neighbors)
    if n > 1
        # Compute remaining deltas inline (avoids delta_encode_vector allocation)
        rem_deltas = Vector{T}(undef, n - 1)
        @inbounds for i in 2:n
            rem_deltas[i-1] = neighbors[i] - neighbors[i-1]
        end
        sections = analyze_delta_patterns_hybrid(rem_deltas, @view(neighbors[2:end]), mil)
        cost += estimate_encoded_value_cost(T(length(sections)), ie)
        for s in sections
            if s.type == :delta
                cost += 1
                cost += estimate_encoded_value_cost(T(length(s.data)), ie)
                for val in s.data; cost += estimate_encoded_value_cost(val, ie); end
            elseif s.type == :run_length
                cost += 2
                cost += estimate_encoded_value_cost(T(length(s.data)÷2), ie)
                for val in s.data; cost += estimate_encoded_value_cost(val, ie); end
            elseif s.type == :interval
                cost += 2
                cost += estimate_encoded_value_cost(T(length(s.data)÷2), ie)
                for val in s.data; cost += estimate_encoded_value_cost(val, ie); end
            end
        end
    end
    return cost
end

function _estimate_base_cost(target::Vector{T}, ie::Symbol, mil::Int, enc_type::Symbol; vertex_id=nothing, stop_deltas::Bool=false, tight_intervals::Bool=false, lr_split::Bool=false) where {T<:Unsigned}
    if enc_type == :interval
        if lr_split
            return _estimate_intervals_lr_cost(target, ie, mil; vertex_id=vertex_id, tight_intervals=tight_intervals)
        end
        return estimate_interval_runlength_encoding_cost(target, ie, mil, 3; vertex_id=vertex_id, tight_intervals=tight_intervals)
    elseif enc_type == :delta
        if stop_deltas
            return _estimate_stop_delta_cost(target, ie; vertex_id=vertex_id)
        else
            cost = estimate_encoded_value_cost(T(length(target) + 1), ie)
            return cost + _estimate_delta_cost(target, ie; vertex_id=vertex_id)
        end
    elseif enc_type == :rle
        cost = estimate_encoded_value_cost(T(length(target) + 1), ie)
        return cost + _estimate_rle_cost(target, ie, mil; vertex_id=vertex_id)
    end
    return 0
end

function _greedy_vertex_search(neighbors::Vector{T}, neighbor_lists::Dict{T,Vector{T}},
                                   ref_window::Vector{T}, ie::Symbol; vertex_id=nothing, copy_blocks::Bool=false, adaptive_copy::Bool=false, fixwidth_ref::Bool=false, ref_dist_bits::Int=0, stop_deltas::Bool=false, adaptive_deltas::Bool=false, split_residual::Bool=false, bv_blocks::Bool=false, compact_copy::Bool=false, tight_intervals::Bool=false, buf::Union{GreedyCostBuffer,Nothing}=nothing, lr_split::Bool=false, multi_ref::Bool=false, adaptive_header::Bool=false, ws::Union{VertexSearchWorkspace{T},Nothing}=nothing, cost_model::Int=COST_MODEL_FULL) where {T<:Unsigned}
    if isempty(neighbors)
        return (ie, :none, MIL_OPTIONS[1], nothing, :interval, false, :interval, MIL_OPTIONS[1])
    end

    best_cost = typemax(Int)
    best_ref_mode = :none
    best_mil = MIN_INTERVAL_LENGTH
    best_enc_type = :interval
    best_ref_result = nothing
    best_use_stop = false
    best_res_enc_type = :interval
    best_res_mil = MIL_OPTIONS[1]

    # Helper: compute exact or estimated cost for a target list
    @inline function _cost_for_list(target::Vector{T}, mil::Int, enc_type::Symbol, use_stop::Bool)
        if buf !== nothing
            # Exact bit-counting via trial encoding
            _greedy_reset!(buf.w2)
            if enc_type == :interval
                if lr_split
                    _write_intervals_lr(buf.w2, target, ie, mil; vertex_id=vertex_id, tight_intervals=tight_intervals)
                else
                    write_intervals_and_residuals(buf.w2, target, ie, mil; vertex_id=vertex_id, tight_intervals=tight_intervals)
                end
            elseif enc_type == :delta
                if use_stop
                    _write_stop_delta(buf.w2, target, ie; vertex_id=vertex_id)
                else
                    write_encoded_value(buf.w2, T(length(target) + 1), ie)
                    if !isempty(target)
                        write_delta(buf.w2, target, ie; vertex_id=vertex_id)
                    end
                end
            elseif enc_type == :rle
                write_encoded_value(buf.w2, T(length(target) + 1), ie)
                if !isempty(target)
                    deltas = delta_encode_vector(target)
                    write_hybrid_mix_encoded_list(buf.w2, deltas, ie, true, mil, false; vertex_id=vertex_id)
                end
            end
            return _total_bits(buf.w2)
        else
            return _estimate_base_cost(target, ie, mil, enc_type; vertex_id=vertex_id, stop_deltas=use_stop, tight_intervals=tight_intervals, lr_split=lr_split)
        end
    end

    # Helper: compute copy cost for a bitmap
    @inline function _copy_cost(bmp::Vector{Bool})
        if cost_model == COST_MODEL_FAST
            return _estimate_copy_cost_fast(bmp)
        end
        if buf !== nothing
            _greedy_reset!(buf.w1)
            if adaptive_copy && copy_blocks
                _write_adaptive_copy(buf.w1, bmp, ie; bv_blocks=bv_blocks, compact_copy=compact_copy)
            elseif copy_blocks
                _write_copy_blocks(buf.w1, bmp, ie)
            else
                write_bitmap_adaptive(buf.w1, bmp, ie)
            end
            return _total_bits(buf.w1)
        else
            if adaptive_copy && copy_blocks
                return _estimate_adaptive_copy_cost(bmp, ie; bv_blocks=bv_blocks, compact_copy=compact_copy)
            elseif copy_blocks
                return _estimate_copy_blocks_cost(bmp, ie)
            else
                return _estimate_adaptive_bitmap_cost(bmp, ie)
            end
        end
    end

    # Standard encoding options loop (always runs)
    # In adaptive header mode, only evaluate STOP-delta and intervals(ADAPTIVE_MIL)
    # In fast cost model, use reduced encoding options (skip RLE, single MIL)
    _enc_options = adaptive_header ? ADAPTIVE_ENCODING_OPTIONS : (cost_model == COST_MODEL_FAST ? FAST_ENCODING_OPTIONS : ENCODING_OPTIONS)
    for (enc_type, mil) in _enc_options
            delta_variants = if adaptive_deltas && enc_type == :delta
                [(false, true), (true, true)]
            else
                [(stop_deltas, false)]
            end

            for (use_stop, is_adaptive) in delta_variants
                flag_cost = is_adaptive ? 1 : 0

                noref_header = _vlc_header_cost(:none, enc_type, mil; adaptive=adaptive_header)
                # In adaptive mode, force STOP-delta for delta encoding
                actual_stop = adaptive_header ? true : use_stop
                base_cost = noref_header + _cost_for_list(neighbors, mil, enc_type, actual_stop) + flag_cost

                if base_cost < best_cost
                    best_cost = base_cost
                    best_ref_mode = :none
                    best_mil = mil
                    best_enc_type = enc_type
                    best_ref_result = nothing
                    best_use_stop = use_stop
                    best_res_enc_type = enc_type
                    best_res_mil = mil
                end
            end
        end

    # ── Two-phase reference search ─────────────────────────────────────────
    if !isempty(ref_window) && length(neighbors) >= 1
        # Phase 1: cheap overlap screening (reuse workspace vectors)
        _oc = ws !== nothing ? ws.overlap_candidates : Vector{Tuple{Int,Vector{T},T}}()
        empty!(_oc)
        for (i, ref_v) in enumerate(ref_window)
            ref_nl = get(neighbor_lists, ref_v, T[])
            isempty(ref_nl) && continue
            ov = _sorted_overlap_count(neighbors, ref_nl)
            ov < 1 && continue
            distance = T(length(ref_window) - i + 1)
            push!(_oc, (ov, ref_nl, distance))
        end

        # Sort by overlap descending, take top-K for full evaluation
        sort!(_oc, by=x -> x[1], rev=true)
        _max_k = cost_model == COST_MODEL_FAST ? MAX_REF_CANDIDATES_PHASE2_FAST : MAX_REF_CANDIDATES_PHASE2
        phase2_k = min(length(_oc), _max_k)

        # Pre-allocated workspace for merge
        _merge_bmp = ws !== nothing ? ws.merge_bmp : Bool[]
        _merge_res = ws !== nothing ? ws.merge_res : T[]
        _best_bmp = ws !== nothing ? ws.best_bmp : Bool[]
        _best_res = ws !== nothing ? ws.best_res : T[]
        best_dist = zero(T)

        # Phase 2: full evaluation of top-K candidates
        for ci in 1:phase2_k
            _, ref_nl, distance = _oc[ci]
            _sorted_merge_ref!(_merge_bmp, _merge_res, neighbors, ref_nl)

            # Copy cost (computed once per candidate)
            dist_cost = fixwidth_ref ? ref_dist_bits : estimate_encoded_value_cost(distance, ie)
            copy_c = _copy_cost(_merge_bmp)

            for (enc_type, mil) in _enc_options
                delta_variants = if adaptive_deltas && enc_type == :delta
                    [(false, true), (true, true)]
                else
                    [(stop_deltas, false)]
                end

                for (use_stop, is_adaptive) in delta_variants
                    flag_cost = is_adaptive ? 1 : 0
                    ref_header = _vlc_header_cost(:reference, enc_type, mil; adaptive=adaptive_header)
                    actual_stop = adaptive_header ? true : use_stop
                    base_ref_c = ref_header + dist_cost + copy_c

                    res_options = split_residual ? ENCODING_OPTIONS : [(enc_type, mil)]
                    for (res_et, res_mil) in res_options
                        res_use_stop = res_et == :delta ? actual_stop : false
                        res_c = base_ref_c
                        if !isempty(_merge_res)
                            res_c += _cost_for_list(_merge_res, res_mil, res_et, res_use_stop)
                        end
                        if split_residual
                            if res_et == enc_type && (res_mil == mil || res_et == :delta)
                                res_c += 1
                            else
                                res_c += 5
                            end
                        end
                        res_c += flag_cost

                        if res_c < best_cost
                            best_cost = res_c
                            best_ref_mode = :reference
                            best_mil = mil
                            best_enc_type = enc_type
                            best_use_stop = use_stop
                            best_res_enc_type = res_et
                            best_res_mil = res_mil
                            # Deferred copy: only copy on improvement
                            _copyto_resize!(_best_bmp, _merge_bmp)
                            _copyto_resize!(_best_res, _merge_res)
                            best_dist = distance
                        end
                    end
                end
            end
        end

        if best_ref_mode == :reference
            best_ref_result = (best_dist, copy(_best_bmp), copy(_best_res))
        end

        # ── Multi-reference search ─────────────────────────────────────────
        if multi_ref && !adaptive_header && length(neighbors) >= 6
            # Reuse phase 1 overlap data for ref1 top-5 selection
            top_k1 = min(length(_oc), 5)

            # Build ref1 candidates with full merge data
            # (ref1 copies are necessary since they're iterated over for ref2 pairing)
            ref1_data = Vector{Tuple{T, Vector{Bool}, Vector{T}, Int}}()
            for k in 1:top_k1
                _, ref_nl, dist1 = _oc[k]
                _sorted_merge_ref!(_merge_bmp, _merge_res, neighbors, ref_nl)
                length(_merge_res) < 3 && continue
                bmp1 = copy(_merge_bmp)
                res1 = copy(_merge_res)
                d_cost = fixwidth_ref ? ref_dist_bits : estimate_encoded_value_cost(dist1, ie)
                ref1_base = d_cost + _copy_cost(bmp1)
                push!(ref1_data, (dist1, bmp1, res1, ref1_base))
            end

            best_multi_cost = best_cost
            best_multi_result = nothing
            best_multi_enc_type = best_enc_type
            best_multi_mil = best_mil
            best_multi_use_stop = best_use_stop

            _merge_bmp2 = Bool[]
            _merge_res2 = T[]

            for (dist1, bmp1, res1, ref1_base) in ref1_data
                for (i, ref_v) in enumerate(ref_window)
                    ref2_nl = get(neighbor_lists, ref_v, T[])
                    isempty(ref2_nl) && continue
                    ov2 = _sorted_overlap_count(res1, ref2_nl)
                    ov2 < 3 && continue

                    _sorted_merge_ref!(_merge_bmp2, _merge_res2, res1, ref2_nl)
                    dist2 = T(length(ref_window) - i + 1)

                    ref2_dist_cost = fixwidth_ref ? ref_dist_bits : estimate_encoded_value_cost(dist2, ie)
                    ref2_copy_cost = _copy_cost(_merge_bmp2)

                    for (enc_type, mil) in _enc_options
                        use_stop = (enc_type == :delta) ? stop_deltas : false
                        multi_header = _vlc_header_cost(:multi_ref, enc_type, mil)
                        total = multi_header + ref1_base + ref2_dist_cost + ref2_copy_cost
                        if !isempty(_merge_res2)
                            total += _cost_for_list(_merge_res2, mil, enc_type, use_stop)
                        end

                        if total < best_multi_cost
                            best_multi_cost = total
                            best_multi_result = (dist1, bmp1, dist2, copy(_merge_bmp2), copy(_merge_res2))
                            best_multi_enc_type = enc_type
                            best_multi_mil = mil
                            best_multi_use_stop = use_stop
                        end
                    end
                end
            end

            if best_multi_result !== nothing
                best_cost = best_multi_cost
                best_ref_mode = :multi_ref
                best_ref_result = best_multi_result
                best_enc_type = best_multi_enc_type
                best_mil = best_multi_mil
                best_use_stop = best_multi_use_stop
                best_res_enc_type = best_multi_enc_type
                best_res_mil = best_multi_mil
            end
        end
    end

    return (ie, best_ref_mode, best_mil, best_ref_result, best_enc_type, best_use_stop, best_res_enc_type, best_res_mil)
end

################################################################################
# END: Greedy compressed graph data
################################################################################

################################################################################
# Command Stream (CS) compressed graph data
#
# Frequency-optimized prefix codes for per-vertex headers.
# Merges empty check, reference flag, and encoding type into a single code tree.
#
#   Code     Bits  Meaning                      Frequency (CNR-2000)
#   0          1   reference + stop_delta        ~53%
#   10         2   no-ref + stop_delta           ~36%
#   1100       4   no-ref + interval(mil=2)      ~6%
#   1101       4   ref + interval(mil=2)         ~1.5%
#   1110       4   empty vertex                  ~3%
#   1111+5b    9   escape: 1b ref + 4b enc_opt   ~0.5%
#
# Weighted average: ~1.72 bits/vertex vs VLC v2 + empty_prefix: ~2.57 bits/vertex
################################################################################

@inline function _cs_header_cost(ref_mode::Symbol, enc_type::Symbol, mil::Int, is_empty::Bool)::Int
    is_empty && return 4                                                          # 1110
    ref_mode == :reference && enc_type == :delta   && return 1                    # 0
    ref_mode == :none      && enc_type == :delta   && return 2                    # 10
    ref_mode == :none      && enc_type == :interval && mil == 2 && return 4       # 1100
    ref_mode == :reference && enc_type == :interval && mil == 2 && return 4       # 1101
    return 9  # 1111 escape: 4 + 1 (ref bit) + 4 (enc_opt)
end

function _write_cs_header(w, ref_mode::Symbol, enc_type::Symbol, mil::Int, is_empty::Bool)
    if is_empty
        # 1110
        write_bit(w, true); write_bit(w, true); write_bit(w, true); write_bit(w, false)
    elseif ref_mode == :reference && enc_type == :delta
        # 0
        write_bit(w, false)
    elseif ref_mode == :none && enc_type == :delta
        # 10
        write_bit(w, true); write_bit(w, false)
    elseif ref_mode == :none && enc_type == :interval && mil == 2
        # 1100
        write_bit(w, true); write_bit(w, true); write_bit(w, false); write_bit(w, false)
    elseif ref_mode == :reference && enc_type == :interval && mil == 2
        # 1101
        write_bit(w, true); write_bit(w, true); write_bit(w, false); write_bit(w, true)
    else
        # 1111 escape + 1-bit ref + 4-bit enc_opt
        write_bit(w, true); write_bit(w, true); write_bit(w, true); write_bit(w, true)
        write_bit(w, ref_mode == :reference)
        _write_enc_opt_tag(w, enc_type, mil)
    end
end

function _read_cs_header(r)
    if !read_bit(r)
        # 0 → reference + stop_delta
        return (false, :reference, :delta, 0)
    end
    if !read_bit(r)
        # 10 → no-ref + stop_delta
        return (false, :none, :delta, 0)
    end
    if !read_bit(r)
        # 110x
        if !read_bit(r)
            # 1100 → no-ref + interval(mil=2)
            return (false, :none, :interval, 2)
        else
            # 1101 → ref + interval(mil=2)
            return (false, :reference, :interval, 2)
        end
    end
    if !read_bit(r)
        # 1110 → empty vertex
        return (true, :none, :delta, 0)
    end
    # 1111 → escape: 1-bit ref + 4-bit enc_opt
    is_ref = read_bit(r)
    enc_type, mil = _read_enc_opt_tag(r)
    ref_mode = is_ref ? :reference : :none
    return (false, ref_mode, enc_type, mil)
end

# Phase 2b: CS command as a single small integer for the context-modeled command
# stream. 0 = empty; else 1 + (ref ? |OPTS| : 0) + (opt_idx - 1) over
# ENCODING_OPTIONS (9 options) → ids 0..18, all direct hybrid-uint tokens (< 2^K).
function _cs_cmd_id(ref_mode::Symbol, enc_type::Symbol, mil::Int, is_empty::Bool)::Int
    is_empty && return 0
    opt = enc_type == :delta ? (:delta, 0) : (enc_type, mil)
    idx = findfirst(==(opt), ENCODING_OPTIONS)
    idx === nothing && (idx = 3)
    return 1 + (ref_mode == :reference ? length(ENCODING_OPTIONS) : 0) + (idx - 1)
end

# Inverse of _cs_cmd_id; returns (is_empty, ref_mode, enc_type, mil) matching the
# _read_cs_header tuple shape (mil = 0 for :delta, as in the raw header codes).
function _cs_cmd_decode(id::Int)
    id == 0 && return (true, :none, :delta, 0)
    t = id - 1
    nopt = length(ENCODING_OPTIONS)
    ref_mode = t >= nopt ? :reference : :none
    enc_type, mil = ENCODING_OPTIONS[(t % nopt) + 1]
    return (false, ref_mode, enc_type, mil)
end

function _cs_vertex_search(neighbors::Vector{T}, neighbor_lists::Dict{T,Vector{T}},
                           ref_window::Vector{T}, ie::Symbol;
                           vertex_id=nothing, compact_copy::Bool=true,
                           tight_intervals::Bool=true,
                           lr_split::Bool=false, fixwidth_ref::Bool=false,
                           ref_dist_bits::Int=0,
                           ws::Union{VertexSearchWorkspace{T},Nothing}=nothing,
                           cost_model::Int=COST_MODEL_FULL) where {T<:Unsigned}
    if isempty(neighbors)
        return (ie, :none, MIL_OPTIONS[1], nothing, :interval, false)
    end

    best_cost = typemax(Int)
    best_ref_mode = :none
    best_mil = MIN_INTERVAL_LENGTH
    best_enc_type = :interval
    best_ref_result = nothing
    best_use_stop = false

    # No-ref evaluation
    _cs_enc_options = cost_model == COST_MODEL_FAST ? FAST_ENCODING_OPTIONS : ENCODING_OPTIONS
    for (enc_type, mil) in _cs_enc_options
        use_stop = enc_type == :delta
        noref_header = _cs_header_cost(:none, enc_type, mil, false)
        base_cost = noref_header + _estimate_base_cost(neighbors, ie, mil, enc_type;
            vertex_id=vertex_id, stop_deltas=use_stop, tight_intervals=tight_intervals,
            lr_split=lr_split)

        if base_cost < best_cost
            best_cost = base_cost
            best_ref_mode = :none
            best_mil = mil
            best_enc_type = enc_type
            best_ref_result = nothing
            best_use_stop = use_stop
        end
    end

    # ── Two-phase reference search ─────────────────────────────────────
    if !isempty(ref_window) && length(neighbors) >= 1
        # Phase 1: cheap overlap screening (reuse workspace vectors)
        _oc = ws !== nothing ? ws.overlap_candidates : Vector{Tuple{Int,Vector{T},T}}()
        empty!(_oc)
        for (i, ref_v) in enumerate(ref_window)
            ref_nl = get(neighbor_lists, ref_v, T[])
            isempty(ref_nl) && continue
            ov = _sorted_overlap_count(neighbors, ref_nl)
            ov < 1 && continue
            distance = T(length(ref_window) - i + 1)
            push!(_oc, (ov, ref_nl, distance))
        end

        sort!(_oc, by=x -> x[1], rev=true)
        _max_k_cs = cost_model == COST_MODEL_FAST ? MAX_REF_CANDIDATES_PHASE2_FAST : MAX_REF_CANDIDATES_PHASE2
        phase2_k = min(length(_oc), _max_k_cs)

        _merge_bmp = ws !== nothing ? ws.merge_bmp : Bool[]
        _merge_res = ws !== nothing ? ws.merge_res : T[]
        _best_bmp = ws !== nothing ? ws.best_bmp : Bool[]
        _best_res = ws !== nothing ? ws.best_res : T[]
        best_dist = zero(T)

        # Phase 2: full evaluation of top-K candidates
        for ci in 1:phase2_k
            _, ref_nl, distance = _oc[ci]
            _sorted_merge_ref!(_merge_bmp, _merge_res, neighbors, ref_nl)

            dist_cost = fixwidth_ref ? ref_dist_bits : estimate_encoded_value_cost(distance, ie)
            copy_c = cost_model == COST_MODEL_FAST ? _estimate_copy_cost_fast(_merge_bmp) : _estimate_adaptive_copy_cost(_merge_bmp, ie; compact_copy=compact_copy)

            for (enc_type, mil) in _cs_enc_options
                use_stop = enc_type == :delta
                ref_header = _cs_header_cost(:reference, enc_type, mil, false)
                ref_c = ref_header + dist_cost + copy_c
                if !isempty(_merge_res)
                    ref_c += _estimate_base_cost(_merge_res, ie, mil, enc_type;
                        vertex_id=vertex_id, stop_deltas=use_stop, tight_intervals=tight_intervals,
                        lr_split=lr_split)
                end

                if ref_c < best_cost
                    best_cost = ref_c
                    best_ref_mode = :reference
                    best_mil = mil
                    best_enc_type = enc_type
                    best_use_stop = use_stop
                    # Deferred copy: only copy on improvement
                    _copyto_resize!(_best_bmp, _merge_bmp)
                    _copyto_resize!(_best_res, _merge_res)
                    best_dist = distance
                end
            end
        end

        if best_ref_mode == :reference
            best_ref_result = (best_dist, copy(_best_bmp), copy(_best_res))
        end
    end

    return (ie, best_ref_mode, best_mil, best_ref_result, best_enc_type, best_use_stop)
end

function write_cmdstream_graph_data(w, neighbor_lists::Dict{T,Vector{T}},
        coding_scheme::Symbol=:children,
        ref_window_size::Int=64;
        integer_encoding::Symbol=:fibonacci,
        compact_copy::Bool=true,
        tight_intervals::Bool=true,
        lr_split::Bool=false,
        cost_model::Int=DEFAULT_COST_MODEL,
        index_sample_k::Int=0) where {T<:Unsigned}

    vs = length(keys(neighbor_lists))
    use_ctx_range = integer_encoding == :context_range
    ie = use_ctx_range ? :fibonacci : integer_encoding
    resid_sink   = use_ctx_range ? CtxRangeEncoder() : nothing
    refdist_sink = use_ctx_range ? CtxRangeEncoder() : nothing
    copy_sink    = use_ctx_range ? BinRangeEncoder() : nothing
    idx_mode = coding_scheme == :index
    # Phase 2b: CS command headers → context range coder; stop-delta / RLE-section
    # flags → adaptive binary coder. In sampled-index (RA) mode both streams are
    # chunked alongside resid/refdist/copy (headers decode sequentially within a
    # chunk, so they don't need to stay in the seekable struct stream).
    cmd_sink  = use_ctx_range ? CtxRangeEncoder() : nothing
    flag_sink = use_ctx_range ? BinRangeEncoder() : nothing
    sampled_idx = idx_mode && index_sample_k > 0
    if use_ctx_range && idx_mode && !sampled_idx
        # Full (per-vertex) index needs a continuous range blob → not seekable.
        error(":context_range index mode requires a sampled index (index_sample_k > 0)")
    end
    # Chunked range streams: sampled index + :context_range. Each k-vertex chunk
    # finalizes the range coders into an independently-decodable byte blob so a
    # reader can seek to the chunk containing v. Mirror of the BG path.
    ctx_chunked = use_ctx_range && sampled_idx

    fixwidth_ref = lr_split  # fixed-width ref distances when lr_split active (same as BG)
    ref_dist_bits = fixwidth_ref ? max(Int(ceil(log2(ref_window_size + 1))), 1) : 0

    # Pre-allocated workspace for vertex search (avoids per-vertex allocations)
    _ws = VertexSearchWorkspace{T}()

    reference_window = T[]
    function add_to_ref_window!(vertex::T)
        push!(reference_window, vertex)
        if length(reference_window) > ref_window_size
            popfirst!(reference_window)
        end
    end

    # Header: Coding scheme (1) + Encoding (3)
    write_bit(w, idx_mode)
    _write_encoding_tag(w, get(ENCODING_TAGS, ie, ENC_TAG_FIBONACCI))

    # --- Helper: encode one vertex's data into a BitWriter ---
    function _encode_cs_vertex!(wr, v::T, current_neighbors::Vector{T})
        # Empty vertex handling
        if isempty(current_neighbors)
            if !idx_mode || sampled_idx
                # Children mode (or sampled index): write explicit empty CS header
                if cmd_sink !== nothing
                    rc_encode_value!(cmd_sink, UInt64(_cs_cmd_id(:none, :delta, 0, true)))
                else
                    _write_cs_header(wr, :none, :delta, 0, true)
                end
            end
            # Full index mode: nothing written (empty = offset[v] == offset[v+1])
            add_to_ref_window!(v)
            return
        end

        # Greedy search with CS headers
        _, ref_mode, mil, ref_result, enc_type, use_stop = _cs_vertex_search(
            current_neighbors, neighbor_lists, reference_window, ie;
            vertex_id=v, compact_copy=compact_copy, tight_intervals=tight_intervals,
            lr_split=lr_split, fixwidth_ref=fixwidth_ref, ref_dist_bits=ref_dist_bits,
            ws=_ws, cost_model=cost_model)

        # Write CS header
        if cmd_sink !== nothing
            rc_encode_value!(cmd_sink, UInt64(_cs_cmd_id(ref_mode, enc_type, mil, false)))
        else
            _write_cs_header(wr, ref_mode, enc_type, mil, false)
        end

        target_list = current_neighbors
        if ref_mode != :none && ref_result !== nothing
            ref_distance, copy_bitmap, residuals = ref_result
            if refdist_sink !== nothing
                rc_encode_value!(refdist_sink, UInt64(ref_distance))
            elseif fixwidth_ref
                write_value(wr, UInt64(ref_distance), ref_dist_bits)
            else
                write_encoded_value(wr, T(ref_distance), ie)
            end
            if copy_sink !== nothing
                for b in copy_bitmap; brc_encode_bit!(copy_sink, b != 0); end
            else
                _write_adaptive_copy(wr, copy_bitmap, ie; compact_copy=compact_copy)
            end
            target_list = residuals
        end

        # Copy-aware rank gaps: recode residual ids as ranks in the complement of
        # the copied set C = current_neighbors \ residuals (identity when :none).
        # Vertex id is mapped into rank space to preserve first-value locality.
        _resid_vid = v
        if resid_sink !== nothing && ref_mode != :none
            _C = setdiff(current_neighbors, target_list)   # sorted, unique, disjoint
            if !isempty(_C)
                target_list = _rankgap_forward(target_list, _C)
                _resid_vid = _rank_of(v, _C)
            end
        end

        # Encode target_list (residual region — gate the range coder)
        if resid_sink !== nothing
            rc_reset_region!(resid_sink); _RESID_SINK[] = resid_sink
        end
        if enc_type == :interval
            if lr_split
                _write_intervals_lr(wr, target_list, ie, mil; vertex_id=_resid_vid, tight_intervals=tight_intervals)
            else
                write_intervals_and_residuals(wr, target_list, ie, mil; vertex_id=_resid_vid, tight_intervals=tight_intervals)
            end
        elseif enc_type == :delta
            if idx_mode && !sampled_idx
                # Full index mode: count-prefixed delta (no stop terminators)
                write_encoded_value(wr, T(length(target_list) + 1), ie)
                if !isempty(target_list)
                    write_delta(wr, target_list, ie; vertex_id=_resid_vid)
                end
            else
                # Children mode (or sampled index): always uses stop_deltas for :delta
                _write_stop_delta(wr, target_list, ie; vertex_id=_resid_vid)
            end
        elseif enc_type == :rle
            write_encoded_value(wr, T(length(target_list) + 1), ie)
            if !isempty(target_list)
                deltas = delta_encode_vector(target_list)
                write_hybrid_mix_encoded_list(wr, deltas, ie, true, mil, false; vertex_id=_resid_vid)
            end
        end
        if resid_sink !== nothing
            _RESID_SINK[] = nothing
        end

        add_to_ref_window!(v)
    end

    # Per-chunk range blobs (ctx_chunked only): concatenated finalized chunks +
    # cumulative byte-offset tables (length nchunks+1, first entry 0, last = blob len).
    resid_blob = UInt8[]; refdist_blob = UInt8[]; copy_blob = UInt8[]
    cmd_blob = UInt8[]; flag_blob = UInt8[]
    resid_coff = Int[0]; refdist_coff = Int[0]; copy_coff = Int[0]
    cmd_coff = Int[0]; flag_coff = Int[0]
    function _finalize_chunk!()
        append!(resid_blob,   rc_finish_and_reset!(resid_sink));   push!(resid_coff,   length(resid_blob))
        append!(refdist_blob, rc_finish_and_reset!(refdist_sink)); push!(refdist_coff, length(refdist_blob))
        append!(copy_blob,    brc_finish_and_reset!(copy_sink));   push!(copy_coff,    length(copy_blob))
        append!(cmd_blob,     rc_finish_and_reset!(cmd_sink));     push!(cmd_coff,     length(cmd_blob))
        append!(flag_blob,    brc_finish_and_reset!(flag_sink));   push!(flag_coff,    length(flag_blob))
        return nothing
    end
    if idx_mode
        # Two-pass encoding: encode per-vertex data to buffer, then write offset table + data
        buf_bw = BitWriter()
        vertex_offsets = Vector{Int}(undef, vs + 1)

        flag_sink === nothing || (_FLAG_SINK[] = flag_sink)
        try
            for v_idx in 1:vs
                # At each chunk boundary (aligned to the sample points), finalize + reset
                # the range coders so each chunk is independently decodable.
                if ctx_chunked && v_idx > 1 && (v_idx - 1) % index_sample_k == 0
                    _finalize_chunk!()
                end
                v = T(v_idx)
                vertex_offsets[v_idx] = _total_bits(buf_bw)
                current_neighbors = sort(get(neighbor_lists, v, T[]))
                _encode_cs_vertex!(buf_bw, v, current_neighbors)
            end
        finally
            _FLAG_SINK[] = nothing
        end
        if ctx_chunked
            # Finalize the last chunk.
            _finalize_chunk!()
        end
        flush_bitwriter(buf_bw; flush_last_bits=true)
        vertex_offsets[vs + 1] = bytes_written(buf_bw) * 8

        # Compute entry width from max offset
        max_offset = vertex_offsets[vs + 1]
        entry_width = max_offset > 0 ? max(Int(ceil(log2(max_offset + 1))), 1) : 1

        if index_sample_k > 0
            # Two-level sampled offsets
            write_bit(w, true)  # sampled flag
            @assert index_sample_k >= 4 && index_sample_k % 4 == 0 "index_sample_k must be a multiple of 4, >= 4"
            if use_ctx_range
                # :context_range: store raw k in 32 bits (chunk sizes can exceed the
                # 8-bit (k/4-1) field's 1024 ceiling). Matches BG + parse_bg_ctxrange_index.
                write_value(w, UInt64(index_sample_k), 32)
            else
                write_value(w, UInt64(index_sample_k ÷ 4 - 1), 8)
            end
            write_value(w, UInt64(entry_width), 6)
            for v_idx in 1:index_sample_k:vs
                write_value(w, UInt64(vertex_offsets[v_idx]), entry_width)
            end
            write_value(w, UInt64(vertex_offsets[vs + 1]), entry_width)

            # :context_range — per-chunk byte-offset tables for the five range
            # streams (each length nchunks+1). Same layout BG writes / the shared
            # parse_bg_ctxrange_index reads.
            if ctx_chunked
                max_boff = max(resid_coff[end], refdist_coff[end], copy_coff[end],
                               cmd_coff[end], flag_coff[end])
                coff_width = max_boff > 0 ? max(Int(ceil(log2(max_boff + 1))), 1) : 1
                write_value(w, UInt64(coff_width), 6)
                for offs in (resid_coff, refdist_coff, copy_coff, cmd_coff, flag_coff)
                    for o in offs
                        write_value(w, UInt64(o), coff_width)
                    end
                end
            end
        else
            # Full offset table
            write_bit(w, false)
            write_value(w, UInt64(entry_width), 6)
            for i in 1:(vs + 1)
                write_value(w, UInt64(vertex_offsets[i]), entry_width)
            end
        end

        # Write buffered encoded data
        buf_data = collect(get_bytes(buf_bw))
        write_bytes(w, buf_data)
    else
        # Children mode: sequential encoding. Gate the flag sink for the whole
        # loop (stop-delta / RLE flags divert to the binary coder in ctx mode).
        flag_sink === nothing || (_FLAG_SINK[] = flag_sink)
        try
            for v_idx in 1:vs
                v = T(v_idx)
                current_neighbors = sort(get(neighbor_lists, v, T[]))
                _encode_cs_vertex!(w, v, current_neighbors)
            end
        finally
            _FLAG_SINK[] = nothing
        end
    end
    # :context_range — finalize the five range streams.
    if ctx_chunked
        # Already finalized per chunk; return the concatenated blobs.
        return (resid_blob, refdist_blob, copy_blob, cmd_blob, flag_blob)
    elseif use_ctx_range
        return (rc_finish!(resid_sink), rc_finish!(refdist_sink), brc_finish!(copy_sink),
                rc_finish!(cmd_sink), brc_finish!(flag_sink))
    end
    return (UInt8[], UInt8[], UInt8[], UInt8[], UInt8[])
end

function read_cmdstream_graph_data(r::BitReader, vs::T, coding_scheme::Symbol, ::Type{T};
        integer_encoding::Symbol=:fibonacci, compact_copy::Bool=true,
        tight_intervals::Bool=true, ref_window_size::Int=64,
        lr_split::Bool=false,
        ctx_range::Bool=false,
        resid_bytes::Vector{UInt8}=UInt8[],
        refdist_bytes::Vector{UInt8}=UInt8[],
        copy_bytes::Vector{UInt8}=UInt8[],
        cmd_bytes::Union{Nothing,Vector{UInt8}}=nothing,
        flag_bytes::Union{Nothing,Vector{UInt8}}=nothing,
        decode_range::Union{Nothing,UnitRange{Int}}=nothing,
        seed_window::Vector{T}=T[],
        ref_resolver::Union{Nothing,Function}=nothing) where {T<:Unsigned}

    neighbor_lists = Dict{T,Vector{T}}()
    # Reassignable per-chunk decoders (swapped by _setup_chunk! in the chunked path).
    resid_source   = ctx_range ? CtxRangeDecoder(resid_bytes) : nothing
    refdist_source = ctx_range ? CtxRangeDecoder(refdist_bytes) : nothing
    copy_source    = ctx_range ? BinRangeDecoder(copy_bytes) : nothing
    # Phase 2b: command + flag streams. nothing = non-ctx / fibonacci files.
    cmd_source  = (ctx_range && cmd_bytes  !== nothing) ? CtxRangeDecoder(cmd_bytes)  : nothing
    flag_source = (ctx_range && flag_bytes !== nothing) ? BinRangeDecoder(flag_bytes) : nothing
    # Capture the caller's flag-source global NOW — before any _setup_chunk! call
    # mutates it — so nested RA chunk decodes (via ref_resolver) restore the OUTER
    # call's decoder on exit, not their own.
    _flag_prev = _FLAG_SOURCE[]
    _read_copy_bits(n::Int) = Bool[brc_decode_bit!(copy_source::BinRangeDecoder) for _ in 1:n]
    # Referenced-vertex adjacency lookup: local (already-decoded) first, then an
    # optional external resolver (random-access: materialise a referenced vertex in
    # an earlier, not-yet-decoded chunk).
    _ref_nbs(ref_v::T) = haskey(neighbor_lists, ref_v) ? neighbor_lists[ref_v] :
        (ref_resolver !== nothing ? ref_resolver(ref_v)::Vector{T} : T[])

    # Read stream metadata
    is_index_mode = read_bit(r)
    enc_tag = _read_encoding_tag(r)
    ie = get(TAG_ENCODINGS, enc_tag, :fibonacci)

    # Read offset table (index mode)
    vertex_offsets = Int[]
    sampled_offsets = Int[]
    is_sampled_index = false
    samp_k = 0
    resid_coff = Int[]; refdist_coff = Int[]; copy_coff = Int[]
    cmd_coff = Int[]; flag_coff = Int[]
    data_start_bit = 0
    if is_index_mode
        is_sampled_index = read_bit(r)
        if is_sampled_index
            # Two-level sampled offsets. :context_range stores raw k in 32 bits;
            # fibonacci uses (k/4-1) in 8 bits.
            local _sk = ctx_range ? Int(read_value(r, 32, UInt64)) :
                (Int(read_value(r, 8, UInt64)) + 1) * 4
            samp_k = _sk
            entry_width = Int(read_value(r, 6, UInt64))
            n_sampled = length(1:_sk:Int(vs))
            for _ in 1:(n_sampled + 1)
                push!(sampled_offsets, Int(read_value(r, entry_width, UInt64)))
            end
            # :context_range — per-chunk byte-offset tables (each length nchunks+1).
            if ctx_range
                coff_width = Int(read_value(r, 6, UInt64))
                for tbl in (resid_coff, refdist_coff, copy_coff, cmd_coff, flag_coff)
                    for _ in 1:(n_sampled + 1)
                        push!(tbl, Int(read_value(r, coff_width, UInt64)))
                    end
                end
            end
            # Struct reader now sits at the first vertex's data — record the absolute
            # bit position so we can seek per chunk.
            data_start_bit = r.bit_count
        else
            # Full offset table
            entry_width = Int(read_value(r, 6, UInt64))
            for _ in 1:(Int(vs) + 1)
                push!(vertex_offsets, Int(read_value(r, entry_width, UInt64)))
            end
        end
    end

    reference_window = T[]
    function add_to_ref_window!(vertex::T)
        push!(reference_window, vertex)
        if length(reference_window) > ref_window_size
            popfirst!(reference_window)
        end
    end

    fixwidth_ref = lr_split
    ref_dist_bits = fixwidth_ref ? max(Int(ceil(log2(ref_window_size + 1))), 1) : 0

    # Chunked :context_range (sampled index): the three range streams are stored as
    # a concatenation of independently-finalized chunks (aligned to sample points).
    # Spin up fresh decoders on chunk c's byte slice with fresh adaptive models.
    ctx_chunked = ctx_range && is_sampled_index && !isempty(resid_coff)
    function _setup_chunk!(c::Int)
        resid_source   = CtxRangeDecoder(resid_bytes[resid_coff[c] + 1 : resid_coff[c + 1]])
        refdist_source = CtxRangeDecoder(refdist_bytes[refdist_coff[c] + 1 : refdist_coff[c + 1]])
        copy_source    = BinRangeDecoder(copy_bytes[copy_coff[c] + 1 : copy_coff[c + 1]])
        if cmd_bytes !== nothing && !isempty(cmd_coff)
            cmd_source  = CtxRangeDecoder(cmd_bytes[cmd_coff[c] + 1 : cmd_coff[c + 1]])
            flag_source = BinRangeDecoder(flag_bytes[flag_coff[c] + 1 : flag_coff[c + 1]])
            _FLAG_SOURCE[] = flag_source   # refresh the global gate to the new chunk decoder
        end
        return nothing
    end

    # Random-access seek mode: decode only `decode_range`, seeding the reference
    # window and jumping the struct reader + range decoders to the range's chunk.
    if decode_range !== nothing
        for u in seed_window
            push!(reference_window, u)
        end
        c0 = ((first(decode_range) - 1) ÷ samp_k) + 1
        bitpos = data_start_bit + sampled_offsets[c0]
        r.index = bitpos + 1
        r.bit_count = bitpos
        _setup_chunk!(c0)
    end

    # Helper to read encoded list body (residual region — gate the range decoder)
    function _read_cs_body(r_in::BitReader, v::T, enc_type::Symbol, mil::Int)
        if resid_source !== nothing
            rc_reset_region!(resid_source::CtxRangeDecoder)
            _RESID_SOURCE[] = resid_source
        end
        try
            if enc_type == :interval
                if lr_split
                    return _read_intervals_lr(r_in, ie, mil, T; vertex_id=v, tight_intervals=tight_intervals)
                end
                return read_intervals_and_residuals(r_in, ie, mil, T; vertex_id=v, tight_intervals=tight_intervals)
            elseif enc_type == :delta
                if is_index_mode && !is_sampled_index
                    # Full index mode: count-prefixed delta (no stop terminators)
                    count = read_encoded_value(r_in, ie, T) - T(1)
                    if count > 0
                        return read_delta(r_in, ie, T; max_elements=Int(count), vertex_id=v)
                    end
                    return T[]
                else
                    # Children mode (or sampled index): always uses stop_deltas
                    return _read_stop_delta(r_in, ie, T; vertex_id=v)
                end
            elseif enc_type == :rle
                count = read_encoded_value(r_in, ie, T) - T(1)
                if count > 0
                    return read_hybrid_mix_encoded_list(r_in, :index, ie, T; max_elements=Int(count), vertex_id=v)
                end
                return T[]
            end
            return T[]
        finally
            resid_source === nothing || (_RESID_SOURCE[] = nothing)
        end
    end

    # Decode one vertex's neighbor list (header already positioned at r_in).
    function _decode_cs_vertex!(r_in::BitReader, v::T)
        is_empty, ref_mode, enc_type, mil = cmd_source !== nothing ?
            _cs_cmd_decode(Int(rc_decode_value!(cmd_source::CtxRangeDecoder))) :
            _read_cs_header(r_in)
        is_empty && return T[]

        if ref_mode != :none
            distance = refdist_source !== nothing ? T(rc_decode_value!(refdist_source::CtxRangeDecoder)) :
                       (fixwidth_ref ? T(read_value(r_in, ref_dist_bits, UInt64)) : read_encoded_value(r_in, ie, T))
            ref_idx = length(reference_window) - Int(distance) + 1
            ref_nbs = T[]
            if 1 <= ref_idx <= length(reference_window)
                ref_nbs = _ref_nbs(reference_window[ref_idx])
            end
            ref_len = length(ref_nbs)
            copy_bitmap = copy_source !== nothing ? _read_copy_bits(ref_len) :
                          _read_adaptive_copy(r_in, ref_len, ie, T; compact_copy=compact_copy)
            # Copy-aware rank gaps: C = copied ids; decode Q in rank space then invert.
            _cs_vid = v; _cs_C = T[]
            if resid_source !== nothing
                _cs_C = reconstruct_from_reference(ref_nbs, copy_bitmap, T[])
                isempty(_cs_C) || (_cs_vid = _rank_of(v, _cs_C))
            end
            residuals = _read_cs_body(r_in, _cs_vid, enc_type, mil)
            if resid_source !== nothing && !isempty(_cs_C)
                residuals = _rankgap_inverse(residuals, _cs_C)
            end
            return reconstruct_from_reference(ref_nbs, copy_bitmap, residuals)
        else
            return _read_cs_body(r_in, v, enc_type, mil)
        end
    end

    # Unified data loop. In seek mode only `decode_range` is visited. (_flag_prev
    # captured at function entry; restored in the finally so nested RA decodes
    # hand back the outer decoder.)
    loop_range = decode_range === nothing ? (1:Int(vs)) : decode_range
    flag_source === nothing || (_FLAG_SOURCE[] = flag_source)
    try
        for v_idx in loop_range
            v = T(v_idx)

            # Chunked full decode: swap to fresh per-chunk decoders at each sample
            # boundary. (In seek mode the chunk is already set up.)
            if ctx_chunked && decode_range === nothing && (v_idx - 1) % samp_k == 0
                _setup_chunk!(((v_idx - 1) ÷ samp_k) + 1)
            end

            if is_index_mode && !is_sampled_index
                # Full index mode: empty vertex = offset[v] == offset[v+1]
                if vertex_offsets[v_idx] == vertex_offsets[v_idx + 1]
                    neighbor_lists[v] = T[]
                    add_to_ref_window!(v)
                    continue
                end
            end

            neighbor_lists[v] = _decode_cs_vertex!(r, v)
            add_to_ref_window!(v)
        end
    finally
        _FLAG_SOURCE[] = _flag_prev
    end

    return neighbor_lists
end

################################################################################
# END: Command Stream (CS) compressed graph data
################################################################################

"""
    find_best_reference_greedy_cost(target::Vector{T}, neighbor_lists::Dict{T,Vector{T}},
                                     available::Set{T}, ref_workspace::RefBuildWorkspace{T},
                                     coding_scheme::Symbol, integer_encoding::Symbol) where {T<:Unsigned}

WebGraph-style greedy cost-based reference selection.
Evaluates all candidates in the reference window and selects the one that minimizes
the total encoding cost (reference ID + copy bitmap + residuals).

This is the key difference from overlap-based selection: instead of choosing the reference
with maximum overlap, we choose the one with minimum encoding cost.

@param target::Vector{T}: the target neighbor list to encode
@param neighbor_lists::Dict{T,Vector{T}}: all neighbor lists (for looking up candidates)
@param available::Set{T}: set of vertices available as references (in the window)
@param ref_workspace::RefBuildWorkspace{T}: workspace for building reference data
@param coding_scheme::Symbol: coding scheme (:children or :index)
@param integer_encoding::Symbol: integer encoding to use for cost estimation

@return best_ref::Union{T,Nothing}: the reference with minimum cost, or nothing if no good reference
"""
function find_best_reference_greedy_cost(target::Vector{T}, neighbor_lists::Dict{T,Vector{T}},
                                        available::Set{T}, ref_workspace::RefBuildWorkspace{T},
                                        coding_scheme::Symbol, integer_encoding::Symbol) where {T<:Unsigned}
    # Skip if no candidates available or target too small
    if isempty(available) || length(target) <= 2
        return nothing
    end

    best_ref = nothing
    best_cost = typemax(Int)  # Initialize to maximum possible cost
    best_copy_bitmap = Bool[]
    best_residuals = T[]

    # Evaluate each candidate reference
    for candidate_ref in available
        # Skip if candidate has no neighbors (can't be a useful reference)
        if !haskey(neighbor_lists, candidate_ref)
            continue
        end

        ref_neighbors = neighbor_lists[candidate_ref]

        # Skip empty references
        if isempty(ref_neighbors)
            continue
        end

        # Build reference encoding: compute copy bitmap and residuals
        create_reference_data!(ref_workspace, target, ref_neighbors)

        # Quick filter: skip if overlap is too low (< REF_ENCODING_TH)
        overlap_count = count(ref_workspace.copy_bitmap)
        if overlap_count < REF_ENCODING_TH
            continue
        end

        # Note: Bandit analysis showed overlap_ratio >= 0.8 predicts reference utility
        # with 83.3% accuracy, but using it as a hard threshold hurts compression.
        # The cost-based selection already captures this pattern, so we rely on that.

        # Calculate actual encoding cost using the cost estimation function
        cost = estimate_reference_encoding_cost(
            candidate_ref,
            ref_workspace.copy_bitmap,
            ref_workspace.residuals,
            coding_scheme,
            integer_encoding
        )

        # Greedy selection: choose minimum cost and save the bitmap/residuals
        if cost < best_cost
            best_cost = cost
            best_ref = candidate_ref
            best_copy_bitmap = copy(ref_workspace.copy_bitmap)
            best_residuals = copy(ref_workspace.residuals)
        end
    end

    # Copy the best result back to the workspace for the caller
    if best_ref !== nothing
        empty!(ref_workspace.copy_bitmap)
        append!(ref_workspace.copy_bitmap, best_copy_bitmap)
        empty!(ref_workspace.residuals)
        append!(ref_workspace.residuals, best_residuals)
    end

    return best_ref
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
        write_bitmap_adaptive(w, copy_bitmap, integer_encoding)  # bitmap (adaptive encoding)
        
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

        # Use hybrid mix encoding for better compression in both modes
        # write_hybrid_mix_encoded_list(w, residuals, integer_encoding, use_mix_mode, MIN_INTERVAL_LENGTH, coding_scheme == :children)
        # Write residuals directly for simplicity
        for residual in residuals
            write_encoded_value(w, residual, integer_encoding)
        end
        # In children mode, write stop value to terminate the residuals
        if coding_scheme == :children
            write_encoded_value(w, T(1), integer_encoding)
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
        copy_bitmap = read_bitmap_adaptive(r, integer_encoding)
        
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
        # Read residuals directly for simplicity (commenting out hybrid mix)
        # return read_hybrid_mix_encoded_list(r, coding_scheme, integer_encoding, T; stop_value=stop_value)

        # Read residuals directly for both modes
        delta_encoded_residuals = T[]
        if coding_scheme == :children
            # In children mode, read until stop value
            while true
                val = read_encoded_value(r, integer_encoding, T)
                if stop_value !== nothing && val == stop_value
                    break
                end
                push!(delta_encoded_residuals, val)
            end
            # Unshift residuals (they were shifted by +1 during writing to avoid zeros)
            delta_encoded_residuals = delta_encoded_residuals .- T(1)
        else
            # Index mode - read residuals directly with max_elements
            if max_elements !== nothing
                for _ in 1:max_elements
                    push!(delta_encoded_residuals, read_encoded_value(r, integer_encoding, T))
                end
            end
        end
        # Reconstruct from delta encoding
        residuals = reconstruct_from_delta(delta_encoded_residuals)
        return residuals
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

    # Bitmap cost (adaptive: 1 flag bit + min(raw, block))
    raw_cost = estimate_encoded_value_cost(UInt32(length(copy_bitmap)) + UInt32(1), integer_encoding) + length(copy_bitmap)
    block_cost = estimate_block_encoding_cost(copy_bitmap, integer_encoding)
    cost += 1 + min(raw_cost, block_cost)

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

###################
# Bitmap encoding
###################

"""
    write_bitpacked_bitmap(w::BitWriter, bits::Vector{Bool})

Write a bitmap packed into bytes (MSB-first within each byte). Length is implied by context.
Pads the last byte with zeros if needed.
"""
function write_bitpacked_bitmap(w::BitWriter, bits::Vector{Bool})
    n = length(bits)
    nbytes = (n + 7) >>> 3
    out = Vector{UInt8}(undef, nbytes)
    ptr = 1
    for b in 1:nbytes
        byte = 0x00
        for j in 0:7
            idx = (b - 1) * 8 + j + 1
            val = (idx <= n) && bits[idx]
            byte |= UInt8(val ? 1 : 0) << (7 - j)
        end
        out[ptr] = UInt8(byte)
        ptr += 1
    end
    write_bytes(w, out)
end

"""
    read_bitpacked_bitmap(r::BitReader, n::Int)::Vector{Bool}

Read a bitmap packed into bytes (MSB-first within each byte), inverse of `write_bitpacked_bitmap`.
`n` is the number of bits (must be known from context).
"""
function read_bitpacked_bitmap(r::BitReader, n::Int)::Vector{Bool}
    nbytes = (n + 7) >>> 3
    bits = Vector{Bool}(undef, n)
    for b in 1:nbytes
        for j in 0:7
            idx = (b - 1) * 8 + j + 1
            bit = read_bit(r)
            if idx <= n
                bits[idx] = bit
            end
        end
    end
    return bits
end

"""
    write_rle_ones_deltas(w::BitWriter, deltas::Vector{T}, varint::Symbol=:fibonacci) where {T<:Unsigned}

Encode a sequence of positive deltas with a run-length scheme for ones:
- For each token: write 1-bit flag (1 = run of ones, 0 = literal delta)
- If run: write run length using `varint`
- If literal: write delta using `varint`
Also writes sequence length up front using `varint` so the decoder knows how many tokens to read.
"""
function write_rle_ones_deltas(w::BitWriter, deltas::Vector{T}, varint::Symbol=:fibonacci) where {T<:Unsigned}
    tokens = 0
    i = 1
    # First pass to count tokens
    while i <= length(deltas)
        if deltas[i] == one(T)
            # count a run of ones
            j = i
            while j <= length(deltas) && deltas[j] == one(T)
                j += 1
            end
            tokens += 1
            i = j
        else
            tokens += 1
            i += 1
        end
    end
    # Write token count
    write_encoded_value(w, T(tokens), varint)
    # Second pass to write tokens
    i = 1
    while i <= length(deltas)
        if deltas[i] == one(T)
            # run token
            write_bit(w, true)
            j = i
            while j <= length(deltas) && deltas[j] == one(T)
                j += 1
            end
            runlen = T(j - i)
            write_encoded_value(w, runlen, varint)
            i = j
        else
            # literal token
            write_bit(w, false)
            write_encoded_value(w, deltas[i], varint)
            i += 1
        end
    end
end

"""
    read_rle_ones_deltas(r::BitReader, varint::Symbol=:fibonacci, ::Type{T}=UInt32) where {T<:Unsigned}

Decode a sequence of deltas encoded with `write_rle_ones_deltas`.
Reads token count, then for each token: 1-bit flag (1=run of ones, 0=literal delta).
"""
function read_rle_ones_deltas(r::BitReader, varint::Symbol=:fibonacci, ::Type{T}=UInt32) where {T<:Unsigned}
    token_count = Int(read_encoded_value(r, varint, T))
    deltas = T[]
    for _ in 1:token_count
        flag = read_bit(r)
        if flag  # run of ones
            runlen = Int(read_encoded_value(r, varint, T))
            for _ in 1:runlen
                push!(deltas, one(T))
            end
        else  # literal
            push!(deltas, read_encoded_value(r, varint, T))
        end
    end
    return deltas
end

"""
    write_bitmap_rle_ones(w::BitWriter, bitmap::Vector{Bool}, varint::Symbol=:fibonacci)

Encode a bitmap using RLE ones-delta encoding optimized for high-density bitmaps.

This converts the bitmap to positions of 1s, computes deltas between positions,
and uses the existing `write_rle_ones_deltas` to compress them efficiently.

Perfect for copy bitmaps which typically have:
- High density (80-90% of bits are 1)
- Consecutive runs of 1s (positions differ by 1)

@param w::BitWriter: the bitwriter
@param bitmap::Vector{Bool}: the bitmap to encode
@param varint::Symbol: encoding for deltas (default: :fibonacci)
"""
function write_bitmap_rle_ones(w::BitWriter, bitmap::Vector{Bool}, varint::Symbol=:fibonacci)
    # Handle empty bitmap
    if isempty(bitmap)
        # Write length = 0 (encoded as 1 to avoid Fibonacci zero issue)
        write_encoded_value(w, UInt32(1), varint)
        return
    end

    # Write bitmap length (add 1 to avoid zero)
    write_encoded_value(w, UInt32(length(bitmap)) + UInt32(1), varint)

    # Find positions of 1s
    ones_positions = UInt32[]
    for (i, bit) in enumerate(bitmap)
        if bit
            push!(ones_positions, UInt32(i))
        end
    end

    # Handle all-zeros bitmap
    if isempty(ones_positions)
        # Write ones_count = 0 (encoded as 1 to avoid zero)
        write_encoded_value(w, UInt32(1), varint)
        return
    end

    # Write number of 1s (add 1 to avoid zero)
    write_encoded_value(w, UInt32(length(ones_positions)) + UInt32(1), varint)

    # Compute deltas between positions (1-based indexing)
    # First position is absolute, rest are deltas
    deltas = UInt32[]
    push!(deltas, ones_positions[1])  # First position (absolute)

    for i in 2:length(ones_positions)
        push!(deltas, ones_positions[i] - ones_positions[i-1])
    end

    # Use existing RLE ones-deltas encoder
    write_rle_ones_deltas(w, deltas, varint)
end

"""
    read_bitmap_rle_ones(r::BitReader, varint::Symbol=:fibonacci)::Vector{Bool}

Decode a bitmap encoded with `write_bitmap_rle_ones`.

@param r::BitReader: the bitreader
@param varint::Symbol: encoding used for deltas (default: :fibonacci)
@return::Vector{Bool}: the decoded bitmap
"""
function read_bitmap_rle_ones(r::BitReader, varint::Symbol=:fibonacci)::Vector{Bool}
    # Read bitmap length (subtract 1)
    length_raw = Int(read_encoded_value(r, varint, UInt32))
    length_val = length_raw - 1

    # Handle empty bitmap
    if length_val == 0
        return Bool[]
    end

    # Read number of 1s (subtract 1)
    ones_count_raw = Int(read_encoded_value(r, varint, UInt32))
    ones_count = ones_count_raw - 1

    # Handle all-zeros bitmap
    if ones_count == 0
        return fill(false, length_val)
    end

    # Read token count for RLE ones-deltas
    token_count = Int(read_encoded_value(r, varint, UInt32))

    # Decode deltas
    deltas = UInt32[]
    for _ in 1:token_count
        flag = read_bit(r)
        if flag
            # Run of ones
            runlen = read_encoded_value(r, varint, UInt32)
            for _ in 1:runlen
                push!(deltas, UInt32(1))
            end
        else
            # Literal delta
            push!(deltas, read_encoded_value(r, varint, UInt32))
        end
    end

    # Reconstruct positions from deltas
    positions = UInt32[]
    if !isempty(deltas)
        current_pos = deltas[1]  # First is absolute position
        push!(positions, current_pos)

        for i in 2:length(deltas)
            current_pos += deltas[i]
            push!(positions, current_pos)
        end
    end

    # Build bitmap from positions
    bitmap = fill(false, length_val)
    for pos in positions
        if 1 <= pos <= length_val
            bitmap[pos] = true
        end
    end

    return bitmap
end

"""
    write_bitmap_adaptive(w::BitWriter, bitmap::Vector{Bool}, varint::Symbol=:fibonacci)

Adaptively encode a bitmap using the most efficient method between two options:
- Block encoding (WebGraph style): best for sparse bitmaps with clustered 1s
- Raw bit encoding: best for dense bitmaps

Format:
  - 1 bit: encoding flag (0=block, 1=raw)
  - Followed by the encoded data in the chosen format

@param w::BitWriter: the bitwriter
@param bitmap::Vector{Bool}: the bitmap to encode
@param varint::Symbol: encoding for length/deltas (default: :fibonacci)
"""
function write_bitmap_adaptive(w::BitWriter, bitmap::Vector{Bool}, varint::Symbol=:fibonacci)
    # ALWAYS write the length first to ensure the stream is self-delimiting
    write_small_count(w, UInt32(length(bitmap)), varint)
    
    if isempty(bitmap)
        return
    end

    # Cost comparison for adaptive choice
    raw_cost = length(bitmap) + 1
    block_cost = estimate_block_encoding_cost(bitmap, varint) + 1

    if block_cost < raw_cost
        write_bit(w, false) # Block flag
        write_block_encoding(w, bitmap, varint)
    else
        write_bit(w, true)  # Raw flag
        for bit in bitmap
            write_bit(w, bit)
        end
    end
end

"""
    read_bitmap_adaptive(r::BitReader, varint::Symbol=:fibonacci)::Vector{Bool}

Decode a self-delimiting adaptive bitmap.
"""
function read_bitmap_adaptive(r::BitReader, varint::Symbol=:fibonacci)::Vector{Bool}
    # Read explicit length
    length_val = Int(read_small_count(r, varint, UInt32))

    if length_val == 0
        return Bool[]
    end

    # Read encoding flag
    use_raw = read_bit(r)

    if !use_raw
        # Block encoding
        bitmap = read_block_encoding(r, varint)
        # Force length to match length_val
        if length(bitmap) < length_val
            append!(bitmap, fill(false, length_val - length(bitmap)))
        elseif length(bitmap) > length_val
            bitmap = bitmap[1:length_val]
        end
        return bitmap
    else
        # Raw encoding
        bitmap = Bool[]
        sizehint!(bitmap, length_val)
        for _ in 1:length_val
            push!(bitmap, read_bit(r))
        end
        return bitmap
    end
end

"""
    write_small_count(w::BitWriter, v::T, varint::Symbol=:fibonacci) where {T<:Unsigned}

Write small nonnegative counts efficiently:
- For v in {0,1,2}: write two-bit code 00, 01, 10 respectively
- Otherwise: write two-bit escape 11 followed by varint of v
"""
function write_small_count(w::BitWriter, v::T, varint::Symbol=:fibonacci) where {T<:Unsigned}
    if v == 0
        write_bit(w, false); write_bit(w, false)
    elseif v == 1
        write_bit(w, false); write_bit(w, true)
    elseif v == 2
        write_bit(w, true); write_bit(w, false)
    else
        write_bit(w, true); write_bit(w, true)
        write_encoded_value(w, v, varint)
    end
end

"""
    read_small_count(r::BitReader, varint::Symbol=:fibonacci, ::Type{T}=UInt32) where {T<:Unsigned}

Read a small nonnegative count encoded with `write_small_count`.
Two-bit codes: 00→0, 01→1, 10→2, 11→varint(v).
"""
function read_small_count(r::BitReader, varint::Symbol=:fibonacci, ::Type{T}=UInt32) where {T<:Unsigned}
    b1 = read_bit(r)
    b2 = read_bit(r)
    if !b1 && !b2
        return T(0)
    elseif !b1 && b2
        return T(1)
    elseif b1 && !b2
        return T(2)
    else
        return read_encoded_value(r, varint, T)
    end
end

"""
    write_bitmap_rle_ones_deltas(w::BitWriter, bits::Vector{Bool}, varint::Symbol=:fibonacci)

Encode a dense bitmap with long 1-runs by alternating:
- Literal zeros-run (flag 0) with length
- Run of ones (flag 1) with length
Writes token count up front, then (flag,length) tokens.
This leverages high density and long runs for compactness.
"""
function write_bitmap_rle_ones_deltas(w::BitWriter, bits::Vector{Bool}, varint::Symbol=:fibonacci)
    n = length(bits)
    tokens = Tuple{Bool,Int}[]  # (is_one_run, length)
    i = 1
    while i <= n
        # zeros run
        z = 0
        while i <= n && !bits[i]
            z += 1; i += 1
        end
        if z > 0
            push!(tokens, (false, z))
        end
        # ones run
        o = 0
        while i <= n && bits[i]
            o += 1; i += 1
        end
        if o > 0
            push!(tokens, (true, o))
        end
    end
    # Write token count and tokens
    write_encoded_value(w, UInt(length(tokens)), varint)
    for (is_one, len) in tokens
        write_bit(w, is_one)
        write_encoded_value(w, UInt(len), varint)
    end
end

"""
    read_bitmap_rle_ones_deltas(r::BitReader, varint::Symbol=:fibonacci, ::Type{T}=UInt32) where {T<:Unsigned}

Decode a dense bitmap encoded with `write_bitmap_rle_ones_deltas`.
Reads token count, then (flag, length) tokens to reconstruct the bitmap.
"""
function read_bitmap_rle_ones_deltas(r::BitReader, varint::Symbol=:fibonacci, ::Type{T}=UInt32) where {T<:Unsigned}
    token_count = Int(read_encoded_value(r, varint, T))
    bits = Bool[]
    for _ in 1:token_count
        is_one = read_bit(r)
        runlen = Int(read_encoded_value(r, varint, T))
        for _ in 1:runlen
            push!(bits, is_one)
        end
    end
    return bits
end

"""
    write_block_encoding(w::BitWriter, copy_bitmap::Vector{Bool}, varint::Symbol=:fibonacci)

Encode a bitmap using WebGraph's block encoding approach.
Alternates between copy blocks (1s) and skip blocks (0s):
[copy B₁, skip B₂, copy B₃, skip B₄, ...]

Format:
- Number of blocks (encoded)
- For each block: length (encoded)
- If block count is even, remaining positions are implicitly copied

This is more efficient than bitmap encoding for sparse reference patterns.
"""
function write_block_encoding(w::BitWriter, copy_bitmap::Vector{Bool}, varint::Symbol=:fibonacci)
    if isempty(copy_bitmap)
        # Empty bitmap: 0 blocks
        write_encoded_value(w, UInt32(1), varint)  # 0+1 to avoid Fibonacci zero
        return
    end

    # Extract blocks: alternating copy (1s) and skip (0s)
    blocks = UInt32[]
    i = 1
    n = length(copy_bitmap)
    expecting_copy = true  # First block is always a copy block

    while i <= n
        block_len = 0

        if expecting_copy
            # Copy block: count consecutive 1s
            while i <= n && copy_bitmap[i]
                block_len += 1
                i += 1
            end
        else
            # Skip block: count consecutive 0s
            while i <= n && !copy_bitmap[i]
                block_len += 1
                i += 1
            end
        end

        if block_len > 0
            push!(blocks, UInt32(block_len))
            expecting_copy = !expecting_copy
        else
            # Switching between copy/skip modes without consuming elements
            # This happens when we start with 0s (skip block length 0)
            if expecting_copy && i <= n && !copy_bitmap[i]
                # Start with skip instead of copy
                push!(blocks, UInt32(0))  # Zero-length copy block
                expecting_copy = false
            else
                break
            end
        end
    end

    # Write number of blocks
    write_encoded_value(w, UInt32(length(blocks)) + UInt32(1), varint)  # +1 to avoid zero

    # Write block lengths
    for block_len in blocks
        write_encoded_value(w, block_len + UInt32(1), varint)  # +1 to avoid zero
    end
end

"""
    read_block_encoding(r::BitReader, varint::Symbol=:fibonacci)::Vector{Bool}

Decode a bitmap from WebGraph's block encoding format.
Reconstructs the bitmap from alternating copy/skip blocks.
"""
function read_block_encoding(r::BitReader, varint::Symbol=:fibonacci)::Vector{Bool}
    # Read number of blocks
    num_blocks_raw = read_encoded_value(r, varint, UInt32)
    num_blocks = Int(num_blocks_raw - UInt32(1))

    if num_blocks == 0
        return Bool[]
    end

    # Read block lengths
    block_lengths = UInt32[]
    for _ in 1:num_blocks
        len_raw = read_encoded_value(r, varint, UInt32)
        push!(block_lengths, len_raw - UInt32(1))
    end

    # Reconstruct bitmap
    bitmap = Bool[]
    expecting_copy = true  # First block is always a copy block

    for block_len in block_lengths
        if expecting_copy
            # Copy block: append 1s
            for _ in 1:block_len
                push!(bitmap, true)
            end
        else
            # Skip block: append 0s
            for _ in 1:block_len
                push!(bitmap, false)
            end
        end
        expecting_copy = !expecting_copy
    end

    return bitmap
end


# -----------------------------------------------------------------------------
# Submodules
# -----------------------------------------------------------------------------

# CG (Clustered Greedy)
include("compression/cge.jl")
using .CG
export CG

# InterEncoding
include("compression/cge/inter_encoding.jl")
using .InterEncoding
export InterEncoding

# Command-Driven Bitstream Compression
include("compression/command_stream.jl")
using .CommandStream
export CommandStream

# -----------------------------------------------------------------------------
end # module Compression
