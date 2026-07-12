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

module IO

using LightGraphs, DataStructures, HDF5, JLD
using ..CustomTypes: UInt24, UInt40
using ..NodeTypes: Node, EmptyNode
using ..CustomLightGraphs: SimpleDiGraph, SimpleGraph, SimpleEdge
using ..Util: infer_uint_custom_type
using ..Constants: BUFFER_SIZE

export BitWriter,
	BitReader,
	write_bit,
	write_bits,
	write_value,
	write_bytes,
	read_bit,
	read_bits,
	read_value,
	peek_bit,
	flush_bitwriter,
	reset_bitwriter!,
	bytes_written,
	get_bytes,
	write_to_io,
	load_jls_serialized,
	serialize_to_jls,
	load_jld_serialized,
	serialize_to_jld,
	load_adjacency_list_from_csv,
	load_adjacency_list,
	load_graph_from_pajek,
	load_triangles

################################################################################
# BitWriter — UInt64-accumulator with Vector{UInt8} buffer
#
# Bits are packed MSB-first into a 64-bit accumulator (right-aligned: the
# first bit written occupies the highest position, the most recent bit is
# at position 0).  When the accumulator fills to 64 bits it is flushed to
# an internal Vector{UInt8} buffer as 8 big-endian bytes.
#
# An optional IO reference can be stored for backward compatibility with
# tests that write through BitWriter then read back via BitReader.
# flush_bitwriter syncs the internal buffer to io when present.
################################################################################

"""
    BitWriter

High-performance bit writer using a UInt64 accumulator backed by a
`Vector{UInt8}` byte buffer.  Bits are written MSB-first (big-endian bit order).
"""
mutable struct BitWriter
    io::Union{Base.IO, Nothing}  # optional IO target (for test backward compat)
    buf::Vector{UInt8}           # internal byte buffer
    pos::Int                     # next write position in buf (1-based)
    accum::UInt64                # accumulator — bits packed right-aligned, MSB-first
    bits_in::Int                 # number of valid bits in accum (0–63)
    bit_count::Int64             # total bits written since creation / last reset
end

"""
    BitWriter(; capacity=32768)

Construct a buffer-only BitWriter (no IO target).  This is the fast path
used by production encoders.
"""
BitWriter(; capacity::Int=32768) = BitWriter(nothing, Vector{UInt8}(undef, max(capacity, 64)), 1, UInt64(0), 0, 0)

"""
    BitWriter(io::Base.IO; capacity=4096*8)

Construct a BitWriter with an IO target for backward compatibility.
Bits are buffered internally; flush_bitwriter syncs to `io`.
"""
function BitWriter(io::Base.IO; capacity=BUFFER_SIZE*8)
    BitWriter(io, Vector{UInt8}(undef, max(capacity, 32768)), 1, UInt64(0), 0, 0)
end

# Write 8 bytes of the accumulator (big-endian) to internal buffer.
@inline function _write_accum_to_buf(w::BitWriter)
    p = w.pos
    need = p + 7
    if need > length(w.buf)
        resize!(w.buf, max(2 * length(w.buf), need))
    end
    v = w.accum
    @inbounds begin
        w.buf[p]   = (v >> 56) % UInt8
        w.buf[p+1] = (v >> 48) % UInt8
        w.buf[p+2] = (v >> 40) % UInt8
        w.buf[p+3] = (v >> 32) % UInt8
        w.buf[p+4] = (v >> 24) % UInt8
        w.buf[p+5] = (v >> 16) % UInt8
        w.buf[p+6] = (v >> 8) % UInt8
        w.buf[p+7] = v % UInt8
    end
    w.pos = p + 8
end

# Flush all 64 accumulator bits to internal buffer.
@inline function _flush_accum(w::BitWriter)
    _write_accum_to_buf(w)
    w.accum = UInt64(0)
    w.bits_in = 0
end

"""
    write_bit(writer::BitWriter, bit::Bool)

Write a single bit to the writer.
"""
@inline function write_bit(w::BitWriter, bit::Bool)
    w.accum = (w.accum << 1) | UInt64(bit)
    w.bits_in += 1
    w.bit_count += 1
    if w.bits_in == 64
        _flush_accum(w)
    end
end

"""
    write_bits(writer::BitWriter, bits::Vector{Bool})

Write a vector of bits to the writer.
"""
function write_bits(w::BitWriter, bits::Vector{Bool})
    @inbounds for bit in bits
        write_bit(w, bit)
    end
end

"""
    write_value(writer::BitWriter, value::T, n::Int) where {T<:Unsigned}

Write the lowest `n` bits of `value` in MSB-first order.
"""
@inline function write_value(w::BitWriter, value::T, n::Int) where {T<:Unsigned}
    n == 0 && return
    v = UInt64(value)
    total = w.bits_in + n
    w.bit_count += n
    if total <= 64
        # Fast path: everything fits in the accumulator
        w.accum = (w.accum << n) | v
        w.bits_in = total
        if total == 64
            _flush_accum(w)
        end
    else
        # Split across accumulator boundary
        first_n = 64 - w.bits_in
        w.accum = (w.accum << first_n) | (v >> (n - first_n))
        _write_accum_to_buf(w)
        rest = n - first_n
        w.accum = v & ((UInt64(1) << rest) - 1)
        w.bits_in = rest
    end
end

"""
    write_bytes(writer::BitWriter, bytes::Vector{UInt8})

Write a vector of bytes to the writer (8 bits each, MSB-first).
"""
function write_bytes(w::BitWriter, bytes::Vector{UInt8})
    @inbounds for byte in bytes
        write_value(w, byte, 8)
    end
end

"""
    flush_bitwriter(writer::BitWriter; flush_last_bits::Bool = false)

Flush buffered bits to the internal byte buffer.
Full bytes are always written. If `flush_last_bits` is true, the final
partial byte is zero-padded on the right and written as well.
When the writer has an IO target, the buffer is synced to IO afterward.
"""
function flush_bitwriter(w::BitWriter; flush_last_bits::Bool = false)
    if w.bits_in > 0
        # Left-align the valid bits for byte extraction
        aligned = w.accum << (64 - w.bits_in)
        full_bytes = w.bits_in >> 3          # div by 8
        remaining_bits = w.bits_in & 7       # mod 8

        # Ensure buffer capacity
        extra = full_bytes + (flush_last_bits && remaining_bits > 0 ? 1 : 0)
        if w.pos + extra - 1 > length(w.buf)
            resize!(w.buf, max(2 * length(w.buf), w.pos + extra))
        end

        # Write full bytes to internal buffer
        @inbounds for i in 0:(full_bytes - 1)
            w.buf[w.pos] = ((aligned >> (56 - 8 * i)) & 0xff) % UInt8
            w.pos += 1
        end

        if flush_last_bits && remaining_bits > 0
            # Write zero-padded last byte
            @inbounds w.buf[w.pos] = ((aligned >> (56 - 8 * full_bytes)) & 0xff) % UInt8
            w.pos += 1
            w.accum = UInt64(0)
            w.bits_in = 0
        else
            # Keep remaining bits (right-aligned in accum)
            if remaining_bits > 0
                w.accum = w.accum & ((UInt64(1) << remaining_bits) - 1)
            else
                w.accum = UInt64(0)
            end
            w.bits_in = remaining_bits
        end
    end

    # Sync internal buffer to IO target (for backward compat with tests)
    if w.io isa IOBuffer
        truncate(w.io, 0)
        seekstart(w.io)
        nb = w.pos - 1
        nb > 0 && write(w.io, @view w.buf[1:nb])
    end
end

"""
    reset_bitwriter!(w::BitWriter)

Reset the writer for reuse (e.g., trial encoding in greedy search).
"""
@inline function reset_bitwriter!(w::BitWriter)
    w.pos = 1
    w.accum = UInt64(0)
    w.bits_in = 0
    w.bit_count = 0
end

"""
    bytes_written(w::BitWriter) -> Int

Number of complete bytes in the internal buffer (excludes accumulator bits).
"""
@inline bytes_written(w::BitWriter) = w.pos - 1

"""
    get_bytes(w::BitWriter) -> SubArray

View of the internal buffer bytes written so far.
"""
@inline get_bytes(w::BitWriter) = @view w.buf[1:w.pos - 1]

"""
    write_to_io(w::BitWriter, io::Base.IO)

Write the internal buffer contents to an IO stream.
"""
function write_to_io(w::BitWriter, io::Base.IO)
    nb = w.pos - 1
    nb > 0 && write(io, @view w.buf[1:nb])
end

################################################################################
# BitReader
################################################################################

mutable struct BitReader
    io::Base.IO
    buffer::Vector{Bool}
    index::Int
    length::Int
    bit_count::Int64 # total bits read
end

"""
    BitReader(io::Base.IO; capacity=BUFFER_SIZE*8)

constructor for BitReader

@param io::Base.IO: the io to read from
@param capacity::Int: the capacity of the buffer
"""
function BitReader(io::Base.IO; capacity=BUFFER_SIZE*8)
    bytes = read(io)
    bits = Vector{Bool}(undef, length(bytes) * 8)
    k = 1
    for byte in bytes
        for i in 7:-1:0
            bits[k] = (byte >> i) & 0x01 == 1
            k += 1
        end
    end
    BitReader(io, bits, 1, k - 1, 0)
end

"""
    read_bit(reader::BitReader)::Bool

read a bit from the reader

@param reader::BitReader: the bit reader to read from
@return::Bool: the bit read
"""
function read_bit(reader::BitReader)::Bool
    if reader.index > reader.length
        error("Attempt to read past end of buffer")
    end
    bit = reader.buffer[reader.index]
    reader.index += 1
    reader.bit_count += 1
    return bit
end

"""
    peek_bit(reader::BitReader)::Bool

Look at the next bit without advancing the reader.

@param reader::BitReader: the bit reader to peek from
@return::Bool: the next bit, or error if past end
"""
function peek_bit(reader::BitReader)::Bool
    if reader.index > reader.length
        error("Attempt to peek past end of buffer")
    end
    return reader.buffer[reader.index]
end

"""
    read_bits(reader::BitReader, n::Int)::Vector{Bool}

read bits from the reader

@param reader::BitReader: the bit reader to read from
@param n::Int: the number of bits to read
@return::Vector{Bool}: the bits read
"""
function read_bits(reader::BitReader, n::Int)::Vector{Bool}
    bits = Vector{Bool}(undef, n)
    for i in 1:n
        bits[i] = read_bit(reader)
    end
    return bits
end

"""
    read_value(reader::BitReader, n::Int, ::Type{T}=UInt) where {T<:Unsigned}

Read the next `n` bits from the bit reader and return them
as an unsigned integer of type `T`.

@param reader::BitReader: the source of bits
@param n::Int: the number of bits to read
@param T::Type: the unsigned return type (default: UInt)
@return::T: the reconstructed unsigned value
"""
function read_value(reader::BitReader, n::Int, ::Type{T}=UInt8) where {T<:Unsigned}
    value = zero(T)
    for _ in 1:n
        value <<= 1
        value |= T(read_bit(reader))
    end
    return value
end

################################################################################
# Serialization
################################################################################

""" 
    load_jls_serialized(filename::AbstractString)

load serialized JLS data
"""
function load_jls_serialized(filename::AbstractString)
	x = open(filename, "r") do file
		deserialize(file)
	end
	return x
end

""" 
    serialize_to_jls(x::Any, filename::AbstractString)

serialize data to JLS format
"""
function serialize_to_jls(x::Any, filename::AbstractString)
	open("$filename.jls", "w") do file
		serialize(file, x)
	end
end

""" 
    load_jld_serialized(name::AbstractString, filename::AbstractString)

load serialized JLD data

NB: to be favored for long term storage
"""
function load_jld_serialized(name::AbstractString, filename::AbstractString)
	x = jldopen(filename, "r") do file
    		read(file, name)
  	end
	return x
end

""" 
    serialize_to_jld(x::Any, name::AbstractString, filename::AbstractString)

serialize data to JLD format

NB: to be favored for long term storage
"""
function serialize_to_jld(x::Any, name::AbstractString, filename::AbstractString)
	jldopen("$filename.jld", "w") do file
		write(file, name, x)
	end
end

################################################################################
# Load graph from various formats
################################################################################

""" 
    load_adjacency_list_from_csv(filename::AbstractString, separator::Char=',', use_header::Bool=false)

load graph from CSV adjacency list
"""
function load_adjacency_list_from_csv(filename::AbstractString, separator::AbstractChar=',', use_header::Bool=false)
	f = open(filename,"r")
	oni = Dict{UInt64,UInt64}()
	edges = Array{Tuple{UInt64,UInt64},1}()
	counter = convert(UInt64,1)
	
	# Skip header line if use_header is true
	if use_header
		readline(f)  # Skip first line (header)
	end
	
	while !eof(f)
		line = strip(readline(f))
		if !startswith(line, "#") && !isempty(line)
			edge = split(line, separator)
			if length(edge) >= 2
				v1 = parse(UInt64,edge[1])
				v2 = parse(UInt64,edge[2])

				if !haskey(oni, v1)
					oni[v1] = counter
					counter += convert(UInt64,1)
				end
				if !haskey(oni, v2)
					oni[v2] = counter
					counter += convert(UInt64,1)
				end
				push!(edges, (oni[v1], oni[v2]))
			end
		end
	end
	close(f)

	gs = length(keys(oni))
	nbits = convert(UInt8, ceil(log(2, gs)))
	V = infer_uint_custom_type(nbits)

	g = SimpleDiGraph{V}()

	# add vertices
	add_vertices!(g, gs)
	
	# add edges
	for edge in edges
		add_edge!(g, convert(V, edge[1]), convert(V, edge[2]))	
	end

	return g::AbstractGraph{V}
end

""" 
    load_adjacency_list(adj_list::Array{Unsigned,2})

load graph from adjacency list

NB: adjcency list should be represented as a 2-dimensional array with 2 rows and 1 column per edge
"""
function load_adjacency_list(adj_list::Array{Unsigned,2})
	oni = Dict{UInt64,UInt64}()
	edges = Array{Tuple{UInt64,UInt64},1}()
	counter = convert(UInt64,1)
	nes = size(adj_list)[2]

	for i in 1:nes
		edge = adj_list[i]
		v1 = edge[1]
		v2 = edge[2]

		if !haskey(oni, v1)
			oni[v1] = counter
			counter += convert(UInt64,1)
		end
		if !haskey(oni, v2)
			oni[v2] = counter
			counter += convert(UInt64,1)
		end
		push!(edges, (oni[v1], oni[v2]))
	end

	gs = length(keys(oni))
	nbits = convert(UInt8, ceil(log(2, gs)))
	V = infer_uint_custom_type(nbits)

	g = SimpleDiGraph{V}()
	# add vertices
	add_vertices!(g,gs)
	
	# add edges
	for edge in edges
		add_edge!(g, convert(V, edge[1]), convert(V, edge[2]))	
	end

	return g::AbstractGraph{V}
end

""" 
    load_graph_from_pajek(filename::AbstractString)

Load net Pajek file
"""
function load_graph_from_pajek(filename::AbstractString)
	f = open(filename,"r")
	
	inside_vertices_section = false
	inside_edges_section = false
	
	vdict = Dict{UInt64, UInt64}()
	vcounter = convert(UInt64,1)
	edges = Array{Tuple{UInt64,UInt64},1}()

	while !eof(f)
		line = lowercase(strip(readline(f)))
		if !startswith(line, "%")
			if startswith(line,"*vertices")
				inside_vertices_section = true
				continue
			elseif startswith(line,"*arcs")
				inside_edges_section = true
				inside_vertices_section = false
				continue
			end
			# Example of vertices section:
			#
			# *Vertices 82670
			# 1 "entity"
			# 2 "thing"
			# 3 "anything"
			# 4 "something"
			# 5 "nothing"
			# 6 "whole"
			if inside_vertices_section
				sa = split(line, ' ')
				# Dictionary of vertices ids and their associated counter
				vdict[parse(UInt64, sa[1])] = vcounter
				vcounter += convert(UInt64, 1)	
			# Example of edges section:
			#
			# *Arcs
			# 1 2
			# 1 3
			# 1 4
			# 1 5
			elseif inside_edges_section
				sa = split(line, ' ')
				vs = vdict[parse(UInt64, sa[1])]
				vt = vdict[parse(UInt64, sa[2])]
				push!(edges, (vs, vt))
			end
		end
	end
	close(f)

	gs = length(keys(vdict))
	nbits = convert(UInt8, ceil(log(2, gs)))
	V = infer_uint_custom_type(nbits)

	g = SimpleDiGraph{V}()
	add_vertices!(g, gs)
	for edge in edges
		add_edge!(g, convert(V, edge[1]), convert(V, edge[2]))
	end

	return g::AbstractGraph{V}
end

""" 
    load_triangles(::Type{T},filename::AbstractString) where {T<:Unsigned}

load triangles list from CSV text-formatted file 

TODO: infer type T 
"""
function load_triangles(::Type{T},filename::AbstractString) where {T<:Unsigned}
	f = open(filename,"r")
	a = (T,T,T)[]
	while !eof(f)
		te = split(readline(f)[2:(end-2)],',')
		v1 = parse(T,te[1])
		v2 = parse(T,te[2])
		v3 = parse(T,te[3])
		push!(a,(v1,v2,v3))
	end
	close(f)
	return a
end

end
