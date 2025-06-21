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
	flush_bitwriter,
	load_jls_serialized, 
	serialize_to_jls, 
	load_jld_serialized, 
	serialize_to_jld, 
	load_adjacency_list_from_csv, 
	load_adjacency_list, 
	load_graph_from_pajek, 
	load_triangles

################################################################################
# BitVectorIO
################################################################################

# # custom IO type for BitVector
# mutable struct BitVectorIO
# 	# bitvector to write to
#     bits::BitVector
# 	# position in the bitvector
# 	position::Int
# end

# # constructor for BitVectorIO with existing BitVector (default position 0)
# BitVectorIO(bits::BitVector) = BitVectorIO(bits, 0)

# # default constructor for BitVectorIO
# BitVectorIO() = BitVectorIO(BitVector(), 0)

# # required IO interface methods
# Base.isopen(io::BitVectorIO) = true
# Base.close(io::BitVectorIO) = nothing
# Base.eof(io::BitVectorIO) = io.position >= length(io.bits)
# Base.position(io::BitVectorIO) = io.position
# Base.seek(io::BitVectorIO, pos::Integer) = (io.position = pos)
# Base.seekend(io::BitVectorIO) = (io.position = length(io.bits))
# Base.skip(io::BitVectorIO, n::Integer) = (io.position += n)

# """
#     write(io::BitVectorIO, byte::UInt8)

# Write a byte to the BitVectorIO

# @param io::BitVectorIO: The bit vector io to write to
# @param byte::UInt8: The byte to write
# @return::Int: The number of bytes written
# """
# function Base.write(io::BitVectorIO, byte::UInt8)
#     # ensure we have enough space
#     needed_bits = io.position + 8
#     if needed_bits > length(io.bits)
#         resize!(io.bits, needed_bits)
#     end
    
#     # write the byte bit by bit
#     for i in 7:-1:0
#         bit = (byte >> i) & 0x01
#         io.bits[io.position + 1] = bit != 0
#         io.position += 1
#     end
#     return 1
# end

# """
#     write(io::BitVectorIO, bytes::Vector{UInt8})

# Write a vector of bytes to the BitVectorIO

# @param io::BitVectorIO: The bit vector io to write to
# @param bytes::Vector{UInt8}: The bytes to write
# @return::Int: The number of bytes written
# """
# function Base.write(io::BitVectorIO, bytes::Vector{UInt8})
#     total_bytes_written = 0
#     for byte in bytes
#         total_bytes_written += write(io, byte)
#     end
#     return total_bytes_written
# end

# """
#     write(io::BitVectorIO, bytes::SubArray{UInt8, 1, Vector{UInt8}, Tuple{UnitRange{Int64}}, true})

# Write a view of bytes to the BitVectorIO

# @param io::BitVectorIO: The bit vector io to write to
# @param bytes::SubArray{UInt8, 1, Vector{UInt8}, Tuple{UnitRange{Int64}}, true}: The bytes to write
# @return::Int: The number of bytes written
# """
# function Base.write(io::BitVectorIO, bytes::SubArray{UInt8, 1, Vector{UInt8}, Tuple{UnitRange{Int64}}, true})
#     total_bytes_written = 0
#     for byte in bytes
#         total_bytes_written += write(io, byte)
#     end
#     return total_bytes_written
# end

# """
#     write_bits!(io::BitVectorIO, bits::BitVector)

# Write individual bits to the BitVectorIO

# @param io::BitVectorIO: The bit vector io to write to
# @param bits::BitVector: The bits to write
# @return::Int: The number of bits written
# """
# function write_bits!(io::BitVectorIO, bits::BitVector)
#     # ensure we have enough space
#     needed_bits = io.position + length(bits)
#     if needed_bits > length(io.bits)
#         resize!(io.bits, needed_bits)
#     end
    
#     # write the bits
#     for i in eachindex(bits)
#         io.bits[io.position + i] = bits[i]
#     end
#     io.position += length(bits)
# end

# """
#     write_bit!(io::BitVectorIO, bit::Bool)

# Write a single bit to the BitVectorIO

# @param io::BitVectorIO: The bit vector io to write to
# @param bit::Bool: The bit to write
# @return::Int: The number of bits written
# """
# function write_bit!(io::BitVectorIO, bit::Bool)
#     # Ensure we have enough space
#     needed_bits = io.position + 1
#     if needed_bits > length(io.bits)
#         resize!(io.bits, needed_bits)
#     end
    
#     # Write the bit
#     io.bits[io.position + 1] = bit
#     io.position += 1
# end

################################################################################
# BitWriter
################################################################################

# mutable struct BitWriter
#     io::Any
#     buffer::Vector{UInt8}
#     nbits::Int
# end

# # Constructors for BitWriter
# BitWriter(io::Base.IO) = BitWriter(io, zeros(UInt8, BUFFER_SIZE), 0)
# BitWriter(io::BitVectorIO) = BitWriter(io, zeros(UInt8, BUFFER_SIZE), 0)
# BitWriter(bits::BitVector) = (empty!(bits); BitWriter(BitVectorIO(bits, 0), zeros(UInt8, BUFFER_SIZE), 0))
# BitWriter() = BitWriter(BitVectorIO(), zeros(UInt8, BUFFER_SIZE), 0)

# """
#     get_bits(w::BitWriter)

# Returns the bits that were actually written to `w.io`.

# @param w::BitWriter: The bit writer to get the bits from
# @return::BitVector: The bits that were actually written
# """
# function get_bits(w::BitWriter)
#     flush!(w)
#     if w.io isa BitVectorIO
#         # Return only the bits that were actually written
#         return w.io.bits[1:w.io.position]
#     else
#         error("BitWriter is not writing to a BitVector")
#     end
# end

# """
#     get_buffer(w::BitWriter)

# Returns the buffer contents of `w.io`.

# @param w::BitWriter: The bit writer to get the buffer from
# @return::Vector{UInt8}: The buffer contents
# """
# function get_buffer(w::BitWriter)
#     flush!(w)
#     if w.io isa BitVectorIO
#         return take!(w.io.buffer)
#     else
#         error("BitWriter is not writing to a BitVectorIO")
#     end
# end

# """
#     write_bits!(w::BitWriter, bits::BitVector)

# Writes `bits` to `w.io` using a buffered approach.

# @param w::BitWriter: The bit writer to write to
# @param bits::BitVector: The bits to write
# """
# function write_bits!(w::BitWriter, bits::BitVector)
#     # Ensure we have enough space
#     needed_bits = w.nbits + length(bits)
#     if needed_bits > length(w.buffer)
#         resize!(w.buffer, needed_bits)
#     end

# 	# write the bits to the buffer
# 	for i in eachindex(bits)
# 		w.buffer[w.nbits + i] = bits[i]
# 	end
# 	w.nbits += length(bits)
# end

# """
#     write_bits!(w::BitWriter, value::T) where {T<:Unsigned}

# Writes `value` to `w.io` using a buffered approach.

# @param w::BitWriter: The bit writer to write to
# @param value::T: The value to write
# """
# function write_bits!(w::BitWriter, value::T) where {T<:Unsigned}
# 	# get the number of bits in the value
#     count = sizeof(T) * 8
# 	# for each bit in the value, write the bit
#     for i in reverse(0:count-1)
#         # get the bit at the current index
#         bit = (value >> i) & 1
#         # get the byte index and bit position
#         byte_idx = div(w.nbits, 8) + 1
#         bit_pos = 7 - (w.nbits % 8)
        
#         # if the buffer is full, flush it
#         if byte_idx > BUFFER_SIZE
#             # write the buffer to the io
#             write(w.io, w.buffer)
#             # reset the buffer
#             w.buffer .= 0
#             # reset the number of bits in the buffer
#             w.nbits = 0
#             # reset the byte index
#             byte_idx = 1
#         end
        
# 		# set the bit in the buffer at the correct position
#         w.buffer[byte_idx] |= UInt8(bit) << bit_pos
#         # increment the number of bits in the buffer
#         w.nbits += 1
#     end
# end

# """
#     write_bytes!(w::BitWriter, bytes::Vector{UInt8})

# Writes `bytes` to `w.io` using a buffered approach.

# @param w::BitWriter: The bit writer to write to
# @param bytes::Vector{UInt8}: The bytes to write
# """
# function write_bytes!(w::BitWriter, bytes::Vector{UInt8})
#     # for each byte in the vector, write the byte
#     for byte in bytes
#         # write the byte to the bitwriter
#         write_bits!(w, byte)
#     end
# end

# """
#     flush!(w::BitWriter)

# Flushes the buffer to `w.io` and resets the bitwriter.

# NB: Should be called when writing is complete.

# @param w::BitWriter: The bit writer to flush
# """
# function flush!(w::BitWriter)
#     if w.nbits > 0
#         if w.io isa BitVectorIO
#             # For BitVectorIO, write individual bits to avoid padding
#             for i in 1:w.nbits
#                 byte_idx = div(i-1, 8) + 1
#                 bit_pos = 7 - ((i-1) % 8)
#                 bit = (w.buffer[byte_idx] >> bit_pos) & 0x01
#                 write_bit!(w.io, bit != 0)
#             end
#             # resize the BitVector to the number of bits written
#             resize!(w.io.bits, w.io.position)
#         else
#             # For regular IO, write complete bytes
#             bytes_to_write = div(w.nbits + 7, 8)
#             write(w.io, view(w.buffer, 1:bytes_to_write))
#         end
#         # reset the buffer
#         w.buffer .= 0
#         # reset the number of bits in the buffer
#         w.nbits = 0
#     end
# end

# mutable struct BitReader
#     io::Any
#     buffer::Vector{UInt8}
#     nbits::Int
# end

# # Constructors for BitReader
# BitReader(io::Base.IO) = BitReader(io, zeros(UInt8, BUFFER_SIZE), 0)
# BitReader(io::BitVectorIO) = BitReader(io, zeros(UInt8, BUFFER_SIZE), 0)

# """
#     read_bit!(br::BitReader)::Bool

# Reads a single bit from `br.io`.

# @param br::BitReader: The bit reader to read from
# @return::Bool: The bit read
# """
# function read_bit!(br::BitReader)::Bool
#     # if the buffer is empty, read the next chunk
#     if br.nbits == 0
#         # Read BUFFER_SIZE bytes into the buffer
#         readbytes!(br.io, br.buffer, BUFFER_SIZE)
#         br.nbits = BUFFER_SIZE * 8
#     end
#     # decrement the number of bits in the buffer
#     br.nbits -= 1
#     # get the byte index and bit position
#     byte_idx = div(br.nbits, 8) + 1
#     bit_pos = 7 - (br.nbits % 8)
#     # return the bit
#     return (br.buffer[byte_idx] >> bit_pos) & 0x01 != 0
# end

# """
#     read_bits!(br::BitReader, n::Int)::BitArray{1}

# Reads `n` bits from `br.io` and returns them as a BitArray.

# @param br::BitReader: The bit reader to read from
# @param n::Int: The number of bits to read
# @return::BitVector: Array of bits read
# """
# function read_bits!(br::BitReader, n::Int)::BitVector
#     # create a bitarray of size n
#     bits = BitArray(undef, n)
#     # for each bit in the bitarray, read the bit
#     for i in eachindex(bits)
#         # read the bit and set the bitarray at the current index
#         bits[i] = read_bit!(br) == 1
#     end
#     # return the bitarray
#     return bits
# end

# """
#     read_bytes!(br::BitReader, n::Int)::Vector{UInt8}

# Reads `n` bytes from `br.io` and returns them as a Vector{UInt8}.

# @param br::BitReader: The bit reader to read from
# @param n::Int: The number of bytes to read
# @return::Vector{UInt8}: Array of bytes read
# """
# function read_bytes!(br::BitReader, n::Int)::Vector{UInt8}
#     # create a vector of size n
#     bytes = Vector{UInt8}(undef, n)
#     # read the bytes into the vector
#     readbytes!(br.io, bytes, n)
#     # return the vector
#     return bytes
# end

"""
    BitWriter

writer for bits
"""
mutable struct BitWriter
    io::Base.IO
    buffer::Vector{Bool}
    index::Int  # next bit to write in buffer
end

"""
    BitWriter(io::Base.IO; capacity=4096*8)

constructor for BitWriter

@param io::Base.IO: The io to write to
@param capacity::Int: The capacity of the buffer
"""
function BitWriter(io::Base.IO; capacity=BUFFER_SIZE*8)
    BitWriter(io, Vector{Bool}(undef, capacity), 1)
end

"""
    write_bit(writer::BitWriter, bit::Bool)

write a bit to the writer

@param writer::BitWriter: The bit writer to write to
@param bit::Bool: The bit to write
"""
function write_bit(writer::BitWriter, bit::Bool)
    writer.buffer[writer.index] = bit
    writer.index += 1
    if writer.index > length(writer.buffer)
        flush_bitwriter(writer)
    end
end

"""
    write_bits(writer::BitWriter, bits::Vector{Bool})

write bits to the writer

@param writer::BitWriter: The bit writer to write to
@param bits::Vector{Bool}: The bits to write
"""
function write_bits(writer::BitWriter, bits::Vector{Bool})
    for bit in bits
        write_bit(writer, bit)
    end
end

"""
    write_value(writer::BitWriter, value::T, n::Int) where {T<:Unsigned}

write a value to the writer

@param writer::BitWriter: The bit writer to write to
@param value::UInt: The value to write
@param n::Int: The number of bits to write
"""
function write_value(writer::BitWriter, value::T, n::Int) where {T<:Unsigned}
    for i in (n-1):-1:0
        write_bit(writer, ((value >> i) & 1) == 1)
    end
end

"""
    write_bytes(writer::BitWriter, bytes::Vector{UInt8})

write bytes to the writer

@param writer::BitWriter: The bit writer to write to
@param bytes::Vector{UInt8}: The bytes to write
"""
function write_bytes(writer::BitWriter, bytes::Vector{UInt8})
    for byte in bytes
        write_value(writer, byte, 8)
    end
end

"""
    flush_bitwriter(writer::BitWriter; flush_last_bits::Bool = false)

flush the writer to the io

@param writer::BitWriter: The bit writer to flush
@param flush_last_bits::Bool: Whether to flush the last padded byte
"""
function flush_bitwriter(writer::BitWriter; flush_last_bits::Bool = false)	
    n = writer.index - 1  # total valid bits
    full_bytes = div(n, 8)
    remaining_bits = n % 8

    # write full bytes
    for i in 1:8:(full_bytes * 8)
        byte = UInt8(0)
        for j in 0:7
            byte |= UInt8(writer.buffer[i + j] ? 1 : 0) << (7 - j)
        end
        write(writer.io, byte)
    end

    # optionally flush last padded byte
    if flush_last_bits && remaining_bits > 0
        byte = UInt8(0)
        offset = full_bytes * 8
        for j in 0:(remaining_bits - 1)
            byte |= UInt8(writer.buffer[offset + j + 1] ? 1 : 0) << (7 - j)
        end
        # remaining low bits stay 0 (padding)
        write(writer.io, byte)
    end

    writer.index = 1  # reset buffer
end

################################################################################
# BitReader
################################################################################

mutable struct BitReader
    io::Base.IO
    buffer::Vector{Bool}
    index::Int
    length::Int
end

"""
    BitReader(io::Base.IO; capacity=BUFFER_SIZE*8)

constructor for BitReader

@param io::Base.IO: The io to read from
@param capacity::Int: The capacity of the buffer
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
    BitReader(io, bits, 1, k - 1)
end

"""
    read_bit(reader::BitReader)::Bool

read a bit from the reader

@param reader::BitReader: The bit reader to read from
@return::Bool: The bit read
"""
function read_bit(reader::BitReader)::Bool
    if reader.index > reader.length
        error("Attempt to read past end of buffer")
    end
    bit = reader.buffer[reader.index]
    reader.index += 1
    return bit
end

"""
    read_bits(reader::BitReader, n::Int)::Vector{Bool}

read bits from the reader

@param reader::BitReader: The bit reader to read from
@param n::Int: The number of bits to read
@return::Vector{Bool}: The bits read
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

@param reader::BitReader: The source of bits
@param n::Int: Number of bits to read
@param T::Type: Unsigned return type (default: UInt)
@return::T: The reconstructed unsigned value
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
    load_adjacency_list_from_csv(filename::AbstractString,separator::Char=',')

load graph from CSV adjacency list
"""
function load_adjacency_list_from_csv(filename::AbstractString, separator::AbstractChar=',')
	f = open(filename,"r")
	oni = Dict{UInt64,UInt64}()
	edges = Array{Tuple{UInt64,UInt64},1}()
	counter = convert(UInt64,1)
	while !eof(f)
		line = strip(readline(f))
		if !startswith(line, "#")
			edge = split(line, separator)
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
