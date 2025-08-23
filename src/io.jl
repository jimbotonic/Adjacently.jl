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
	peek_bit,
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
# BitWriter
################################################################################

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

@param writer::BitWriter: the bit writer to write to
@param value::UInt: the value to write
@param n::Int: the number of bits to write
"""
function write_value(writer::BitWriter, value::T, n::Int) where {T<:Unsigned}
    for i in (n-1):-1:0
        write_bit(writer, ((value >> i) & 1) == 1)
    end
end

"""
    write_bytes(writer::BitWriter, bytes::Vector{UInt8})

write bytes to the writer

@param writer::BitWriter: the bit writer to write to
@param bytes::Vector{UInt8}: the bytes to write
"""
function write_bytes(writer::BitWriter, bytes::Vector{UInt8})
    for byte in bytes
        write_value(writer, byte, 8)
    end
end

"""
    flush_bitwriter(writer::BitWriter; flush_last_bits::Bool = false)

flush the writer to the io

@param writer::BitWriter: the bit writer to flush
@param flush_last_bits::Bool: whether to flush the last padded byte
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
    BitReader(io, bits, 1, k - 1)
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
