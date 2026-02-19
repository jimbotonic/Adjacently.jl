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

module CommandStream

using ...IO: BitWriter, BitReader, write_bit, read_bit, write_value, read_value
using ..Compression: write_encoded_value, read_encoded_value, write_bitmap_adaptive, read_bitmap_adaptive

export AdjacencyCommand,
       StopCmd,
       GapCmd,
       ReferenceCmd,
       IntervalCmd,
       RLECmd,
       BlockCmd,
       write_command,
       read_command,
       write_vertex_commands,
       read_vertex_commands,
       reconstruct_from_commands

"""
    Command-Driven Bitstream Compression

This module implements a highly compact, flexible graph compression format where each 
vertex adjacency list is represented as a sequence of "instructions" for the decoder.

### The Flag Map (Prefix Codes)
To allow the decoder to differentiate between commands while maintaining high density, 
we use a prefix-free variable-length flag system.

| Flag | Command | Description | Bit Cost |
| :--- | :--- | :--- | :--- |
| `0` | **STOP** | End of current vertex adjacency list. | 1 bit |
| `10` | **GAP** | A single neighbor (stored as a delta from the previous). | 2 bits |
| `110` | **REFERENCE**| Copy neighbors from a previous vertex in the window. | 3 bits |
| `1110` | **INTERVAL** | A range of consecutive neighbors. | 4 bits |
| `11110` | **RLE** | A repeated gap sequence. | 5 bits |
| `11111` | **BLOCK** | Alternating included/excluded blocks of neighbors. | 5 bits |

### Key Features
1. **Granular Recursion**: Unlike standard WebGraph, multiple REFERENCE commands can 
   be mixed with INTERVAL and GAP commands for a single vertex.
2. **Stop-Value Efficiency**: The 1-bit '0' flag makes terminating lists extremely cheap.
3. **No Residual Headers**: Since commands are streamed, we don't need to pre-declare 
   the number of residuals or intervals.
"""

abstract type AdjacencyCommand end

"""
    StopCmd

Terminal command representing the end of an adjacency list.
Flag: `0`
"""
struct StopCmd <: AdjacencyCommand end

"""
    GapCmd{T}

A single neighbor stored as a delta from the previous neighbor.
Flag: `10`
"""
struct GapCmd{T<:Unsigned} <: AdjacencyCommand
    delta::T
end

"""
    ReferenceCmd{T}

Instruction to copy neighbors from a previous vertex.
Flag: `110`
"""
struct ReferenceCmd{T<:Unsigned} <: AdjacencyCommand
    ref_id::T
    copy_bitmap::Vector{Bool}
end

"""
    IntervalCmd{T}

A range of consecutive neighbors.
Flag: `1110`
"""
struct IntervalCmd{T<:Unsigned} <: AdjacencyCommand
    start::T
    len::T
end

"""
    RLECmd{T}

A sequence of neighbors with a repeating gap.
Flag: `1111`
"""
struct RLECmd{T<:Unsigned} <: AdjacencyCommand
    gap::T
    count::T
end

"""
    BlockCmd{T}

Alternating included and excluded blocks of neighbors starting from `start`.
Flag: `11111`
"""
struct BlockCmd{T<:Unsigned} <: AdjacencyCommand
    start::T
    block_lengths::Vector{T} # [included_len1, excluded_len1, included_len2, ...]
end

################################################################################
# Writing Functions
################################################################################

"""
    write_command(w::BitWriter, cmd::AdjacencyCommand, encoding::Symbol)

Write a single adjacency command to the bitstream.
"""
function write_command(w::BitWriter, cmd::StopCmd, encoding::Symbol)
    write_bit(w, false) # Flag 0
end

function write_command(w::BitWriter, cmd::GapCmd{T}, encoding::Symbol) where {T<:Unsigned}
    write_bit(w, true)  # Flag 1
    write_bit(w, false) # Flag 10
    write_encoded_value(w, cmd.delta, encoding)
end

function write_command(w::BitWriter, cmd::ReferenceCmd{T}, encoding::Symbol) where {T<:Unsigned}
    write_bit(w, true)  # Flag 1
    write_bit(w, true)  # Flag 11
    write_bit(w, false) # Flag 110
    write_encoded_value(w, cmd.ref_id, encoding)
    write_bitmap_adaptive(w, cmd.copy_bitmap, encoding)
end

function write_command(w::BitWriter, cmd::IntervalCmd{T}, encoding::Symbol) where {T<:Unsigned}
    write_bit(w, true)  # Flag 1
    write_bit(w, true)  # Flag 11
    write_bit(w, true)  # Flag 111
    write_bit(w, false) # Flag 1110
    write_encoded_value(w, cmd.start, encoding)
    write_encoded_value(w, cmd.len, encoding)
end

function write_command(w::BitWriter, cmd::RLECmd{T}, encoding::Symbol) where {T<:Unsigned}
    write_bit(w, true)  # Flag 1
    write_bit(w, true)  # Flag 11
    write_bit(w, true)  # Flag 111
    write_bit(w, true)  # Flag 1111
    write_bit(w, false) # Flag 11110
    write_encoded_value(w, cmd.gap, encoding)
    write_encoded_value(w, cmd.count, encoding)
end

function write_command(w::BitWriter, cmd::BlockCmd{T}, encoding::Symbol) where {T<:Unsigned}
    write_bit(w, true)  # Flag 1
    write_bit(w, true)  # Flag 11
    write_bit(w, true)  # Flag 111
    write_bit(w, true)  # Flag 1111
    write_bit(w, true)  # Flag 11111
    write_encoded_value(w, cmd.start, encoding)
    write_encoded_value(w, T(length(cmd.block_lengths)), encoding)
    for len in cmd.block_lengths
        write_encoded_value(w, len, encoding)
    end
end

"""
    write_vertex_commands(w::BitWriter, commands::Vector{AdjacencyCommand}, encoding::Symbol)

Write a sequence of commands for a vertex, ensuring it ends with a StopCmd.
"""
function write_vertex_commands(w::BitWriter, commands::Vector{AdjacencyCommand}, encoding::Symbol)
    for cmd in commands
        write_command(w, cmd, encoding)
        if cmd isa StopCmd
            return
        end
    end
    # Force stop if not present
    write_command(w, StopCmd(), encoding)
end

################################################################################
# Reading Functions
################################################################################

"""
    read_command(r::BitReader, encoding::Symbol, ::Type{T}) where {T<:Unsigned}

Read the next command from the bitstream.
"""
function read_command(r::BitReader, encoding::Symbol, ::Type{T}) where {T<:Unsigned}
    if !read_bit(r)
        return StopCmd()
    end
    
    # Prefix 1...
    if !read_bit(r)
        # Prefix 10
        delta = read_encoded_value(r, encoding, T)
        return GapCmd{T}(delta)
    end
    
    # Prefix 11...
    if !read_bit(r)
        # Prefix 110
        ref_id = read_encoded_value(r, encoding, T)
        copy_bitmap = read_bitmap_adaptive(r, encoding)
        return ReferenceCmd{T}(ref_id, copy_bitmap)
    end
    
    # Prefix 111...
    if !read_bit(r)
        # Prefix 1110
        start = read_encoded_value(r, encoding, T)
        len = read_encoded_value(r, encoding, T)
        return IntervalCmd{T}(start, len)
    end
    
    # Prefix 1111...
    if !read_bit(r)
        # Prefix 11110 (RLE)
        gap = read_encoded_value(r, encoding, T)
        count = read_encoded_value(r, encoding, T)
        return RLECmd{T}(gap, count)
    end
    
    # Prefix 11111 (Block)
    start = read_encoded_value(r, encoding, T)
    num_blocks = Int(read_encoded_value(r, encoding, T))
    block_lengths = T[]
    for _ in 1:num_blocks
        push!(block_lengths, read_encoded_value(r, encoding, T))
    end
    return BlockCmd{T}(start, block_lengths)
end

"""
    read_vertex_commands(r::BitReader, encoding::Symbol, ::Type{T}) where {T<:Unsigned}

Read commands for a vertex until a StopCmd is encountered.
"""
function read_vertex_commands(r::BitReader, encoding::Symbol, ::Type{T}) where {T<:Unsigned}
    commands = AdjacencyCommand[]
    while true
        cmd = read_command(r, encoding, T)
        push!(commands, cmd)
        if cmd isa StopCmd
            break
        end
    end
    return commands
end

"""
    reconstruct_from_commands(commands::Vector{AdjacencyCommand}, neighbor_lists::Dict{T,Vector{T}}, ::Type{T}) where {T<:Unsigned}

Reconstruct the adjacency list for a vertex from its sequence of commands.
"""
function reconstruct_from_commands(commands::Vector{AdjacencyCommand}, neighbor_lists::Dict{T,Vector{T}}, ::Type{T}) where {T<:Unsigned}
    neighbors = T[]
    last_val = zero(T)

    for cmd in commands
        if cmd isa StopCmd
            break
        elseif cmd isa GapCmd
            val = (length(neighbors) == 0) ? cmd.delta : (last_val + cmd.delta)
            push!(neighbors, val)
            last_val = val
        elseif cmd isa IntervalCmd
            for i in 0:(Int(cmd.len) - 1)
                val = T(cmd.start + i)
                push!(neighbors, val)
                last_val = val
            end
        elseif cmd isa RLECmd
            for _ in 1:Int(cmd.count)
                val = last_val + cmd.gap
                push!(neighbors, val)
                last_val = val
            end
        elseif cmd isa BlockCmd
            curr = cmd.start
            including = true
            for blen in cmd.block_lengths
                if including
                    for i in 0:(Int(blen) - 1)
                        val = T(curr + i)
                        push!(neighbors, val)
                        last_val = val
                    end
                end
                curr += blen
                including = !including
            end
        elseif cmd isa ReferenceCmd
            ref_neighbors = get(neighbor_lists, cmd.ref_id, T[])
            for (i, bit) in enumerate(cmd.copy_bitmap)
                if bit && i <= length(ref_neighbors)
                    val = ref_neighbors[i]
                    push!(neighbors, val)
                    last_val = val
                end
            end
        end
    end
    
    # Adjacency lists should be sorted and unique
    unique!(sort!(neighbors))
    return neighbors
end

end # module CommandStream
