#
# Adjacently: Julia Complex Directed Networks Library
# Copyright (C) 2016-2026 Jimmy Dubuisson <jimmy@dubuisson.ch>
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

module MGS

using LightGraphs, DataStructures, Logging
using ..CustomTypes: UInt24, UInt40
using ..NodeTypes: Node, EmptyNode
using ..CustomLightGraphs: SimpleDiGraph, SimpleGraph, SimpleEdge
using ..Util: infer_uint_custom_type, to_bytes
using ..IO: BitWriter, write_bytes, flush_bitwriter, BitReader,
	write_value, read_value, write_to_io, bytes_written, get_bytes, reset_bitwriter!

using ..Compression: huffman_encoding, get_huffman_codes!, decode_huffman_values,
	delta_encode_vector, write_elias_coding, read_elias_coding,
	write_golomb, read_golomb, write_fibonacci, read_fibonacci,
	write_zeta, read_zeta, write_compressed_graph_data, read_compressed_graph_data,
	write_greedy_graph_data, read_greedy_graph_data,
	write_cmdstream_graph_data, read_cmdstream_graph_data,
	COST_MODEL_FULL, COST_MODEL_FAST, DEFAULT_COST_MODEL

using ..Compression.CG: encode_level, decode_level, CGParams, CGStats

using ..Graph: get_basic_stats, get_in_out_degrees, get_out_degrees
using ..Relabeling: relabel_vertices, relabel_graph
using ..Constants: GOLOMB_BASE, ZETA_BASE


# constants for new header format
# Header structure: 'MGS' (3 bytes) + major/minor version (2 bytes) + flags (2 bytes) + vertices (5 bytes)
# 
# New flag byte structure:
# Byte 1: graph_type (2 bits) | coding_scheme (2 bits) | integer_encoding (4 bits)  
# Byte 2: option_flags (8 bits)
#
# Graph type constants (2 bits)
const GRAPH_TYPE_DIRECTED = 0b00
const GRAPH_TYPE_UNDIRECTED = 0b01
#
# Coding scheme constants   (2 bits)
const CODING_SCHEME_CHILDREN = 0b00
const CODING_SCHEME_INDEX = 0b01
#
# Varint encoding constants (4 bits, values 1-6)
const INT_ENCODING_ELIAS_GAMMA = 0x1
const INT_ENCODING_ELIAS_DELTA = 0x2
const INT_ENCODING_GOLOMB = 0x3
const INT_ENCODING_FED = 0x4
const INT_ENCODING_ZETA = 0x5
const INT_ENCODING_FIBONACCI = 0x6
#
# ── Self-describing byte 2 layout ───────────────────────────────────────────
# Algorithm IDs (0x00–0x0F): no additional params needed
const ALG_LEGACY_MGS = 0x00
const ALG_HUFFMAN    = 0x01
const ALG_BG        = 0x02   # BG with recommended defaults
const ALG_CS         = 0x03   # CS with recommended defaults
const ALG_CG        = 0x04   # CG with recommended defaults
# 0x05–0x0F reserved

# Parameter ranges: algorithm + explicit params
const PARAM_BG_BASE = 0x10   # BG + params (64 slots: 0x10–0x4F)
const PARAM_BG_MAX  = 0x4F
const PARAM_CS_BASE  = 0x50   # CS + params  (32 slots: 0x50–0x6F)
const PARAM_CS_MAX   = 0x6F
const PARAM_CG_BASE = 0x70   # CG + params (144 slots: 0x70–0xFF)
const PARAM_CG_MAX  = 0xFF

# ── BG window lookup tables ─────────────────────────────────────────────────
const BG_WINDOW_SIZES = [8, 16, 32, 64]   # 2-bit index
const BG_WINDOW_TO_IDX = Dict(8 => 0, 16 => 1, 32 => 2, 64 => 3)

# ── CS window lookup tables ─────────────────────────────────────────────────
const CS_WINDOW_SIZES = [4, 8, 16, 32, 64, 128, 256, 512]  # 3-bit index
const CS_WINDOW_TO_IDX = Dict(4 => 0, 8 => 1, 16 => 2, 32 => 3, 64 => 4, 128 => 5, 256 => 6, 512 => 7)

# ── CG window + membership + interval lookup tables ─────────────────────────
const CG_WINDOW_SIZES = [8, 16, 32, 64, 128, 256]  # 6 values, index 0–5
const CG_WINDOW_TO_IDX = Dict(8 => 0, 16 => 1, 32 => 2, 64 => 3, 128 => 4, 256 => 5)
# interval_mode: 0=none, 1=intervals, 2=intervals+lr_split, 3=intervals+lr_split+tight_deltas

# ── CG parameter encoding ────────────────────────────────────────────────────
# Mixed-radix: membership(2) × window(6) × interval_mode(4) × mil(3) = 144
# varint is stored in the header's integer_encoding field (byte 1, 4 bits)
# membership: 0=:implicit_ranges, 1=:stop
const CG_V32_MEMBERSHIP_ORDER = [:implicit_ranges, :stop]
const CG_V32_MEMBERSHIP_TO_IDX = Dict(:implicit_ranges => 0, :stop => 1)
# mil: 0=3, 1=4, 2=5
const CG_MIL_SIZES = [3, 4, 5]
const CG_MIL_TO_IDX = Dict(3 => 0, 4 => 1, 5 => 2)

# maximum number of vertices
MGS_MAX_SIZE = 0xffffffffff

# Helper functions for new header format
function create_header_flags(graph_type::Symbol, coding_scheme::Symbol, integer_encoding::Symbol, option_flags::UInt8)
    """Create header flags bytes according to new format."""
    
    # Map symbols to constants
    gt = graph_type == :directed ? GRAPH_TYPE_DIRECTED : GRAPH_TYPE_UNDIRECTED
    cs = coding_scheme == :children ? CODING_SCHEME_CHILDREN : CODING_SCHEME_INDEX
    
    ie = if integer_encoding == :elias_gamma
        INT_ENCODING_ELIAS_GAMMA
    elseif integer_encoding == :elias_delta
        INT_ENCODING_ELIAS_DELTA
    elseif integer_encoding == :golomb
        INT_ENCODING_GOLOMB
    elseif integer_encoding == :fed
        INT_ENCODING_FED
    elseif integer_encoding == :zeta
        INT_ENCODING_ZETA
    elseif integer_encoding == :fibonacci
        INT_ENCODING_FIBONACCI
    else
        error("Unsupported integer encoding: $integer_encoding")
    end
    
    # Construct byte 1: graph_type (2 bits) | coding_scheme (2 bits) | integer_encoding (4 bits)
    byte1 = UInt8((gt << 6) | (cs << 4) | ie)
    
    # Byte 2 is the option flags
    byte2 = option_flags
    
    return (byte1, byte2)
end

function decode_header_flags(byte1::UInt8, byte2::UInt8)
    """Decode header flags bytes according to new format."""
    
    # Extract fields from byte 1
    graph_type_bits = (byte1 >> 6) & 0b11
    coding_scheme_bits = (byte1 >> 4) & 0b11
    integer_encoding_bits = byte1 & 0b1111
    
    # Map bits back to symbols
    graph_type = graph_type_bits == GRAPH_TYPE_DIRECTED ? :directed : :undirected
    coding_scheme = coding_scheme_bits == CODING_SCHEME_CHILDREN ? :children : :index
    
    integer_encoding = if integer_encoding_bits == INT_ENCODING_ELIAS_GAMMA
        :elias_gamma
    elseif integer_encoding_bits == INT_ENCODING_ELIAS_DELTA
        :elias_delta
    elseif integer_encoding_bits == INT_ENCODING_GOLOMB
        :golomb
    elseif integer_encoding_bits == INT_ENCODING_FED
        :fed
    elseif integer_encoding_bits == INT_ENCODING_ZETA
        :zeta
    elseif integer_encoding_bits == INT_ENCODING_FIBONACCI
        :fibonacci
    else
        error("Unknown integer encoding: $integer_encoding_bits")
    end
    
    # Decode option flags
    option_flags = byte2
    
    return (graph_type, coding_scheme, integer_encoding, option_flags)
end

# ── BG parameter encode/decode (v3.2) ────────────────────────────────────────

"""
    encode_bg_params(; ref_window_size, copy_blocks, stop_deltas,
                       lr_split, exact_costing) → UInt8

Encode BG params into byte2 value (range 0x10–0x4F).
Bit layout (6 bits, offset = byte2 − 0x10):
  Bits 5-4: window (2 bits: 8/16/32/64)
  Bit 3:    copy_blocks
  Bit 2:    stop_deltas
  Bit 1:    lr_split
  Bit 0:    adaptive_header
"""
function encode_bg_params(; ref_window_size::Int=64, copy_blocks::Bool=true,
		stop_deltas::Bool=true, lr_split::Bool=false, adaptive_header::Bool=false)
	w = ref_window_size <= 8 ? 8 : ref_window_size
	haskey(BG_WINDOW_TO_IDX, w) || error("Unsupported BG window size: $ref_window_size (must be 8/16/32/64)")
	widx = BG_WINDOW_TO_IDX[w]
	offset = (widx << 4) | (copy_blocks ? 0x08 : 0x00) |
	         (stop_deltas ? 0x04 : 0x00) |
	         (lr_split ? 0x02 : 0x00) |
	         (adaptive_header ? 0x01 : 0x00)
	return UInt8(PARAM_BG_BASE + offset)
end

"""
    decode_bg_params(byte2::UInt8) → NamedTuple

Decode BG params from byte2 value. Returns all params needed by the loader.
Hardcoded defaults: compact_copy=true, tight_intervals=true.
When lr_split=true: implies fixwidth_ref=true.
"""
function decode_bg_params(byte2::UInt8)
	offset = Int(byte2) - PARAM_BG_BASE
	widx = (offset >> 4) & 0x03
	ref_window_size = BG_WINDOW_SIZES[widx + 1]
	copy_blocks = (offset & 0x08) != 0
	stop_deltas = (offset & 0x04) != 0
	lr_split = (offset & 0x02) != 0
	adaptive_header = (offset & 0x01) != 0
	return (ref_window_size=ref_window_size, copy_blocks=copy_blocks,
	        adaptive_copy=copy_blocks, stop_deltas=stop_deltas,
	        fixwidth_ref=lr_split,
	        compact_copy=true, tight_intervals=true,
	        lr_split=lr_split, adaptive_header=adaptive_header)
end

"""
    encode_cs_params(; ref_window_size, compact_copy, tight_intervals, lr_split) → UInt8

Encode CS params into a v3.2 byte2 value (range 0x50–0x6F).
Bit layout (5 bits, offset = byte2 − 0x50):
  Bit 4-2: ref_window_size  (3 bits → 0=4, 1=8, …, 7=512)
  Bit 1:   lr_split
  Bit 0:   (reserved, must be 0)
Implied: compact_copy=true, tight_intervals=true (always active).
"""
function encode_cs_params(; ref_window_size::Int=64, compact_copy::Bool=true,
		tight_intervals::Bool=true, lr_split::Bool=false)
	haskey(CS_WINDOW_TO_IDX, ref_window_size) || error("Unsupported CS window size: $ref_window_size")
	widx = CS_WINDOW_TO_IDX[ref_window_size]
	offset = (widx << 2) | (lr_split ? 0x02 : 0x00)
	return UInt8(PARAM_CS_BASE + offset)
end

"""
    decode_cs_params(byte2::UInt8) → NamedTuple

Decode CS params from byte2 value.
v3.2: bit 1 = lr_split, compact_copy=true, tight_intervals=true (hardcoded).
"""
function decode_cs_params(byte2::UInt8)
	offset = Int(byte2) - PARAM_CS_BASE
	widx = (offset >> 2) & 0x07
	ref_window_size = CS_WINDOW_SIZES[widx + 1]
	lr_split = (offset & 0x02) != 0
	return (ref_window_size=ref_window_size, compact_copy=true,
	        tight_intervals=true, lr_split=lr_split)
end

"""
    encode_cg_params(params::CGParams) → UInt8

Encode CG params into byte2 value (range 0x70–0xFF).
Mixed-radix: membership(2) × window(6) × interval_mode(4) × mil(3) = 144.
varint is stored separately in the header's integer_encoding field.
"""
function encode_cg_params(params::CGParams)
	haskey(CG_V32_MEMBERSHIP_TO_IDX, params.membership) || error("Unsupported CG membership: $(params.membership) (only :implicit_ranges, :stop)")
	haskey(CG_WINDOW_TO_IDX, params.intra_ref_window) || error("Unsupported CG window: $(params.intra_ref_window)")
	mil = params.intra_adapt_mil  # encode adapt_mil (= intra_mil in practice)
	haskey(CG_MIL_TO_IDX, mil) || error("Unsupported CG mil: $mil (must be 3, 4, or 5)")

	midx = CG_V32_MEMBERSHIP_TO_IDX[params.membership]
	widx = CG_WINDOW_TO_IDX[params.intra_ref_window]
	milidx = CG_MIL_TO_IDX[mil]
	interval_mode = if !params.intra_intervals
		0
	elseif params.intra_lr_split && params.intra_tight_deltas
		3
	elseif params.intra_lr_split
		2
	else
		1
	end
	offset = midx * 72 + widx * 12 + interval_mode * 3 + milidx
	@assert offset <= 143 "CG param offset $offset exceeds max 143"
	return UInt8(PARAM_CG_BASE + offset)
end

"""
    decode_cg_params(byte2::UInt8; varint::Symbol=:fibonacci) → CGParams

Decode CG params from byte2 value. Mixed-radix:
membership(2) × window(6) × interval_mode(4) × mil(3) = 144.
varint is passed from the header's integer_encoding field.
"""
function decode_cg_params(byte2::UInt8; varint::Symbol=:fibonacci)
	offset = Int(byte2) - PARAM_CG_BASE

	membership_idx = div(offset, 72)
	rem1 = offset - membership_idx * 72
	window_idx = div(rem1, 12)
	rem2 = rem1 - window_idx * 12
	interval_mode = div(rem2, 3)
	mil_idx = rem2 - interval_mode * 3

	membership = CG_V32_MEMBERSHIP_ORDER[membership_idx + 1]
	intra_ref_window = CG_WINDOW_SIZES[window_idx + 1]
	intra_intervals = interval_mode >= 1
	intra_lr_split = interval_mode >= 2
	intra_tight_deltas = interval_mode == 3
	mil = CG_MIL_SIZES[mil_idx + 1]

	return CGParams(
		L=128, varint=varint, count_varint=:fibonacci, gap=:fibonacci,
		degree=:elias_delta, perm_strategy=:blockpos, undirected_pairs=false,
		membership=membership, inter_strategy=:lists,
		intra_ref_enabled=true, intra_ref_window=intra_ref_window,
		intra_ref_rle=false, intra_block_try=false,
		positions_mode=:delta, additions_mode=:auto,
		min_cluster_density=0.0,
		intra_intervals=intra_intervals, intra_mil=mil, intra_greedy_mil=false,
		intra_mgs=false, intra_zigzag=true, intra_stop_deltas=true,
		intra_copy_blocks=true, intra_copy_adaptive=true,
		intra_ref_fixwidth=true, intra_ref_vlc=false,
		intra_add_adaptive=true, intra_raw_adaptive=true,
		intra_adapt_mil=mil, intra_lr_split=intra_lr_split,
		intra_tight_deltas=intra_tight_deltas,
	)
end

"""
    _bg_default_params() → NamedTuple

Hardcoded recommended BG defaults (matching ALG_BG = 0x02).
Equivalent to: window=64, all features on.
"""
_bg_default_params() = decode_bg_params(UInt8(PARAM_BG_BASE + 0x3F))

"""
    _cs_default_params() → NamedTuple

Hardcoded recommended CS defaults (matching ALG_CS = 0x03).
Equivalent to: window=64, compact_copy=true, tight_intervals=true.
"""
_cs_default_params() = decode_cs_params(encode_cs_params(ref_window_size=64, compact_copy=true, tight_intervals=true, lr_split=false))

"""
    _cg_default_params() → CGParams

Hardcoded recommended CG defaults (matching ALG_CG = 0x04).
implicit_ranges, perm, window=16, no intervals, no lr_split.
"""
_cg_default_params() = CGParams(
	L=128, varint=:fibonacci, count_varint=:fibonacci, gap=:fibonacci,
	degree=:elias_delta, perm_strategy=:blockpos, undirected_pairs=false,
	membership=:implicit_ranges, inter_strategy=:perm,
	intra_ref_enabled=true, intra_ref_window=16,
	intra_ref_rle=false, intra_block_try=false,
	positions_mode=:delta, additions_mode=:auto,
	min_cluster_density=0.0,
	intra_intervals=false, intra_mil=4, intra_greedy_mil=false,
	intra_mgs=false, intra_zigzag=true, intra_stop_deltas=true,
	intra_copy_blocks=true, intra_copy_adaptive=true,
	intra_ref_fixwidth=true, intra_ref_vlc=false,
	intra_add_adaptive=true, intra_raw_adaptive=true,
	intra_adapt_mil=4, intra_lr_split=false,
)

# Export the functions we want to make available
export write_mgs3_graph,
       write_compressed_mgs3_graph,
       load_mgs3_graph,
       load_compressed_mgs3_graph,
       write_huffman_compressed_mgs3_graph,
	   load_huffman_compressed_mgs3_graph,
	   load_greedy_mgs3_graph,
	   write_bg_mgs3_graph,
	   write_cs_mgs3_graph,
	   write_cg_mgs3_graph

################################################################################
# Write uncompressed MGS v3 graph
################################################################################

"""
    write_mgs3_graph(g::AbstractGraph{T},filename::AbstractString) where {T<:Unsigned}

write graph in format MGS v3

Supported encoding schemes:
- :children
- :index
"""
function write_mgs3_graph(g::AbstractGraph{T}, filename::AbstractString, encoding::Symbol=:children) where {T<:Unsigned}
  	# Header 12 bytes: 
	# -> version: 'MGS' 3 bytes string
	# -> major + minor version: 2 bytes
	# -> flags: 2 bytes
	#	 * Byte 1 (2 bits + 6 bits):
	#    	- graph type (2 bits): 			0x0: directed graph | 0x1: undirected graph
	#	 	- compression scheme (6 bits): 	0x0: no compression
	#	 * Byte 2:
	#		- coding scheme: 		0x0: data section only | 0x1: index+data section with implicit numbering of vertices
	#	 	- reserved flags: 		0x0: reserved
	# -> # vertices: 5 bytes position 
	#
	# <'MGS' string 3 bytes> + <16 bits major|minor version> + <flags 2 bytes> + <# vertices 5 bytes>
	vs = vertices(g)
	# number of vertices
	gs = convert(UInt64, length(vs))

	# if the graph has more than 2^40-1 vertices, `T` should be `UInt64`
	if gs > MGS_MAX_SIZE
		error("Input graph cannot have more than 2^40-1 vertices")
	end
	
	# `n_size_t` is the number of bytes needed to represent the graph vertices
	# NB: `n_size_t` (the computed size) may be lower than `sizeof(T)` of the graph type in parameter
	n_size_t = convert(UInt8, ceil(log(2, gs)/8))
	# NB: `T` should be an unsigned integer of size 1,2,4,8 bytes
	p_size_t = sizeof(T)
	
	#  The difference of size should >= 0 
	diff_size = p_size_t - n_size_t

	# Create header using new format (uncompressed graph = no options)
	option_flags = UInt8(ALG_LEGACY_MGS)  # No compression options for uncompressed graph

	# We don't have integer encoding for uncompressed graphs, so use elias_delta as default
	flag_byte1, flag_byte2 = create_header_flags(:directed, encoding, :elias_delta, option_flags)

	# Construct 12-byte header
	header_bytes = UInt8[
		# 'MGS' signature (3 bytes)
		0x4d, 0x47, 0x53,
		# Major version = 3, Minor version = 2 (2 bytes)
		0x03, 0x02,
		# Flag bytes (2 bytes)
		flag_byte1, flag_byte2,
		# Number of vertices (5 bytes, little-endian UInt40)
		(gs & 0xff), ((gs >> 8) & 0xff), ((gs >> 16) & 0xff), ((gs >> 24) & 0xff), ((gs >> 32) & 0xff)
	]

	# create the output file (with extension .mgs)
	f = open(filename * ".mgs", "w")
	
	### write header (12 bytes total with new format)
	write(f, header_bytes)
	
	# coding scheme: data section only with implicit numbering of vertices
	if encoding == :children
		# array of 0x00 of size `n_size_t`
		stop_seq = [0x00 for i in 1:n_size_t]
		for v in vs
			ovs = outneighbors(g, v)
			for c in ovs
				# `c` is of type `T` (e.g. UInt8, UInt16, UInt32, UInt64)
				# `diff_size` is the difference in size between `T` and `n_size_t`
				# `n_size_t` is the number of bytes needed to represent the graph vertices
				# `p_size_t` is the size of the graph type in parameter
				#bytes = reverse(reinterpret(UInt8, [c]))[(diff_size+1):end]
				bytes = to_bytes(c)[(diff_size+1):end]
				write(f, bytes)
			end
			write(f, stop_seq)
		end
	# coding scheme: index+data sections with implicit numbering of vertices
	elseif encoding == :index
		# number of children for each vertex
		# NB: `T` should have a sufficient size to store the number of children
		ods = T[]
		# flattened list of children for all the vertices
		# NB: `T` should have a sufficient size to store the children indices
		children = T[]

		for v in vs
			ovs = outneighbors(g, v)
			push!(ods, convert(T, length(ovs)))
			for o in ovs
				push!(children, o)	
			end
		end

		### write index section
		for o in ods
			# `o` is of type `T` (e.g. UInt8, UInt16, UInt32, UInt64)
			# `diff_size` is the difference in size between `T` and `n_size_t`
			# `n_size_t` is the number of bytes needed to represent the graph vertices
			# `p_size_t` is the size of the graph type in parameter
			#bytes = reverse(reinterpret(UInt8, [o]))[(diff_size+1):end]
			bytes = to_bytes(o)[(diff_size+1):end]
			write(f, bytes)
		end
		### write data section
		for c in children
			# `c` is of type `T` (e.g. UInt8, UInt16, UInt32, UInt64)
			# `diff_size` is the difference in size between `T` and `n_size_t`
			# `n_size_t` is the number of bytes needed to represent the graph vertices
			# `p_size_t` is the size of the graph type in parameter
			#bytes = reverse(reinterpret(UInt8, [c]))[(diff_size+1):end]
			bytes = to_bytes(c)[(diff_size+1):end]
			write(f, bytes)
		end
	end
	close(f)
end

################################################################################
# Load uncompressed MGS v3 graph
################################################################################

""" 
    load_mgs3_graph(filename::AbstractString)

load graph in format MGS v3
"""
function load_mgs3_graph(filename::AbstractString)
	f = open(filename, "r")
	### read header
  	# 7 bytes: <3 bytes string 'MGS'> + <2 bytes major/minor> + <2 bytes flags>
	# 5 bytes (40 bits): number of vertices
	###
	# 7-bytes version
	# e.g. UInt8[0x4d, 0x47, 0x53, 0x03, 0x00, 0x00, 0x00]
	version = read(f,7)
	minor_version = version[5]

	# number of vertices — format depends on version
	if minor_version >= 0x01
		# v3.1+: 5 bytes little-endian UInt40
		raw = read(f, 5)
		gs = UInt64(raw[1]) | (UInt64(raw[2]) << 8) | (UInt64(raw[3]) << 16) | (UInt64(raw[4]) << 24) | (UInt64(raw[5]) << 32)
	else
		# v3.0: 5 bytes big-endian
		gs = reinterpret(UInt64, vcat(reverse(read(f,5)),[0x00,0x00,0x00]))[1]
	end

	# `n_size_v` is the number of bits needed to represent the graph vertices
	n_bits_v = convert(UInt8, ceil(log(2, gs)))
	# Get appropriate custom UInt type based on number of bits needed
	V = infer_uint_custom_type(n_bits_v)

	# graph type is 2 first bits of 6th byte of header
	# 0x0: directed graph | 0x1: undirected graph
	graph_type = version[6] >> 6

	# coding scheme — v3.1 stores it in byte 1 bits 5-4; v3.0 in byte 2 top nibble
	if minor_version >= 0x01
		encoding = ((version[6] >> 4) & 0b11) == CODING_SCHEME_CHILDREN ? :children : :index
	else
		encoding = version[7] >> 4 == 0x0 ? :children : :index
	end

	# intialize graph g according to graph type
	g = if graph_type == 0x0
		SimpleDiGraph{V}()
	else
		SimpleGraph{V}()
	end
	
	# vertex set
	vs = range(1, stop=gs)

	# add vertices to graph
    add_vertices!(g, gs)

	# coding scheme: data section only with implicit numbering of vertices
	# NB: each list of children is terminated with 0
	if encoding == :children
		# read data
		children = V[]
		while !eof(f)
			c = read(f, sizeof(V))
			# NB: reinterpret function reads a byte array in little endian format
			# NB: reverse function reverses the byte array in little endian format
			push!(children, reinterpret(V, reverse(c))[1])
		end
		
		# add edges
		pos = 1
		n_children = length(children)

		for i in 1:length(vs)
			if pos <= n_children
				source = convert(V, i)
				while children[pos] != 0x00
					target = children[pos]
					add_edge!(g, source, target)
					pos += 1
				end
				# skip 0x00
				pos += 1
			else
				# if we reached the last child, we are done
				break
			end
		end
	# coding scheme: index+data sections with implicit numbering of vertices
	elseif encoding == :index
		# read index
		# NB: 
		# - position indices are 1-based
		# - each position indicates the index of the first child of a vertex
		ods = V[]
		for i in 1:gs
			p = read(f, sizeof(V))
			push!(ods, reinterpret(V, reverse(p))[1])
		end
		# read data
		children = V[]
		while !eof(f)
			c = read(f, sizeof(V))
			push!(children, reinterpret(V, reverse(c))[1])
		end

		# add edges
		current_pos = 1
		for v in vs
			source = convert(V, v)
			# if we reached the last parent vertex
			if v == length(vs)
				pos1 = current_pos
				pos2 = length(children)
			else
				pos1 = current_pos
				# position of the last child of vertex v
				pos2 = current_pos + ods[v] - 1
				current_pos = pos2 + 1
			end
			# add edges for each child
			for p in pos1:pos2
				target = children[p]
				add_edge!(g, source, target)
			end
		end
	end
	close(f)
	return g::AbstractGraph{V}
end

################################################################################
# Write compressed MGS v3 graph
################################################################################

"""
    write_compressed_mgs3_graph(g::AbstractGraph{T}, filename::AbstractString, compression::Symbol=:huffman, coding_scheme::UInt8=0x00) where {T<:Unsigned}

Write graph in MGS v3 format with specified encoding and compression scheme

Parameters:
- g: Input graph
- filename: Output filename
- compression: Compression scheme to use (default: :huffman)
- encoding: Coding scheme (:children for children section only, :index for index+children sections)

Supported encoding schemes:
- :children
- :index

Supported compression types:
- :huffman - Huffman coding (default)

@returns nothing
"""
function write_compressed_mgs3_graph(g::AbstractGraph{T}, filename::AbstractString, coding_scheme::Symbol=:children, compression::Symbol=:huffman, integer_encoding::Symbol=:fibonacci, use_mix_mode::Bool=true, reference_enabled::Bool=true, recursive_reference::Bool=true, ref_window_size::Int=1024) where {T<:Unsigned}
	# supported compression
	supported_compressions = [:huffman]
	supported_integer_encodings = [:elias_gamma, :elias_delta, :golomb, :fibonacci, :zeta, :fed]

	if compression == :huffman
        write_huffman_compressed_mgs3_graph(g, filename, coding_scheme)
    else
        error("Unsupported compression scheme: $compression. Supported schemes are :huffman")
    end
end

################################################################################
# Huffman compressed MGS v3 graph
################################################################################

"""
    write_huffman_compressed_mgs3_graph(g::AbstractGraph{T}, filename::AbstractString, coding_scheme::UInt8=0x00) where {T<:Unsigned}

write graph in a compressed MGS v3 (Huffman compression scheme)

Parameters:
- g: Input graph
- filename: Output filename
- encoding: Coding scheme (:children for children section only, :index for index+children sections)

@returns nothing
"""
function write_huffman_compressed_mgs3_graph(g::AbstractGraph{T}, filename::AbstractString, coding_scheme::Symbol=:children) where {T<:Unsigned}
	# Header format: see MGS_HEADER.md
	vs = vertices(g)
	# number of vertices
	gs = convert(UInt64, length(vs))

	# if the graph has more than 2^40-1 vertices, `T` should be `UInt64`
	if gs > MGS_MAX_SIZE
		error("Input graph cannot have more than 2^40-1 vertices")
	end
	
	# `n_bits_v` is the number of bits needed to represent the graph vertices
	n_bits_v = convert(UInt8, ceil(log(2, gs)))
	# Get appropriate custom UInt type based on number of bits needed
	V = infer_uint_custom_type(n_bits_v)
	# `n_size_t` is the number of bytes needed to represent the graph vertices
	# NB: `n_size_t` (the computed size) may be lower than `sizeof(T)` of the graph type in parameter
	n_size_t = sizeof(V)
	# NB: `T` should be an unsigned integer of size 1,2,4,8 bytes
	p_size_t = sizeof(T)
	
	#  The difference of size should >= 0 
	diff_size = p_size_t - n_size_t

	# Create header using new format (Huffman uses special option flag)
	option_flags = UInt8(ALG_HUFFMAN)  # Huffman compression flag

	# For Huffman, we still need an integer encoding (doesn't matter since Huffman overrides)
	flag_byte1, flag_byte2 = create_header_flags(:directed, coding_scheme, :elias_delta, option_flags)

	# Construct 12-byte header
	header_bytes = UInt8[
		# 'MGS' signature (3 bytes)
		0x4d, 0x47, 0x53,
		# Major version = 3, Minor version = 2 (2 bytes)
		0x03, 0x02,
		# Flag bytes (2 bytes)
		flag_byte1, flag_byte2,
		# Number of vertices (5 bytes, little-endian UInt40)
		(gs & 0xff), ((gs >> 8) & 0xff), ((gs >> 16) & 0xff), ((gs >> 24) & 0xff), ((gs >> 32) & 0xff)
	]

	# create the output file (with extension .mgz)
	f = open(filename * ".mgz", "w")

	# flattened list of children for the whole graph
	children = V[]

	# frequencies of each vertex (in- and out- degrees)
	in_degrees, out_degrees = get_in_out_degrees(g)

	# convert the out-degrees to the custom type `V`
	out_degrees = Dict{V,V}(k => convert(V, v) for (k, v) in out_degrees)
	# convert in-degrees to the custom type `V`
	in_degrees = Dict{V,V}(k => convert(V, v) for (k, v) in in_degrees)
	
	# collect all children and in-degrees
	for v in vs
		ovs = outneighbors(g, v)
		for o in ovs
			push!(children, convert(V, o))
		end
	end

	# if coding scheme is :children, we add a special stop sequence whose frequency is equal to `gs` the number of vertices
	# NB: vertex frequencies are in the range [0, gs-1]
	if coding_scheme == :children
		# add the stop sequence 0 to the frequencies
		# NB: 
		# - vertex IDs are in the range [1, gs]
		# - 0 is the stop sequence and its frequency is set to the total number of vertices
		in_degrees[zero(V)] = convert(V, gs)
	end

	@info("generating Huffman tree")
	# get Huffman encoding tree
	tree = huffman_encoding(in_degrees)

	@info("getting Huffman codes")
	# get Huffman codes in C (C:code -> value)
	C = Dict{BitArray{1}, V}()
	get_huffman_codes!(tree, C, BitArray{1}())

	# reverse dictionary (R: value -> code)
	R = Dict{V, BitArray{1}}()
	[R[value] = key for (key, value) in C]

	@info("writing header section (new format)")
	### write header (12 bytes total with new format)
	write(f, header_bytes)
	
	@info("writing frequency section")
	### write frequency section
	for v in vs
		bytes = to_bytes(in_degrees[v])[(diff_size+1):end]
		write(f, bytes)
	end

	if coding_scheme == :children
		# write stop sequence in frequency section
		bytes = to_bytes(in_degrees[zero(V)])[(diff_size+1):end]
		write(f, bytes)
		
		@info("writing data section with stop sequence")
		### write data section
		cdata = BitArray{1}()
		# stop sequence is the code associated to 0
		# NB: R has a length of `gs+1` because of the stop sequence
		stop_seq_code = R[zero(V)] 
		
		for v in vs
			ovs = outneighbors(g, v)
			for c in ovs
				code = R[convert(V, c)]
				append!(cdata, code)
			end
			append!(cdata, stop_seq_code)
		end
	elseif coding_scheme == :index
		@info("writing index section")
		### write index section
		for v in vs
			bytes = to_bytes(out_degrees[v])[(diff_size+1):end]
			write(f, bytes)
		end
		
		@info("writing data section")
		### write data section
		cdata = BitArray{1}()
		for c in children
			code = R[convert(V, c)]
			append!(cdata, code)
		end
	end
	
	# write the final padding
	# add padding if necessary to make the data section a multiple of 8 bytes
	scd = convert(UInt32, length(cdata))
	sp = 8 - scd%8
	for i in 1:sp
		push!(cdata,0)
	end
	
	# number of bytes to write
	nb = round(Int, length(cdata)/8)

	for i in 1:nb
		# BitArray chunks type is UInt64
		# we only need to keep the last byte of the chunks
		byte = reinterpret(UInt8, cdata[(i-1)*8+1:(i-1)*8+8].chunks)[1]
		write(f, byte)
	end

	# last byte indicates by how many bytes cdata was padded
	b = 0xff
	write(f, b >> sp)

	close(f)
end

################################################################################
# BG (Baseline greedy) compressed MGS v3 graph
################################################################################

"""
    write_bg_mgs3_graph(g, filename; ...)

Write graph in MGS v3 format with BG (Baseline greedy) compression.
Produces a `.mgz` file with v3.1 header (params encoded in byte2).
"""
function write_bg_mgs3_graph(g::AbstractGraph{T}, filename::AbstractString;
		coding_scheme::Symbol=:children, integer_encoding::Symbol=:fibonacci,
		ref_window_size::Int=64, copy_blocks::Bool=true,
		stop_deltas::Bool=true,
		lr_split::Bool=false, exact_costing::Bool=false,
		multi_ref::Bool=false,
		adaptive_header::Bool=false,
		cost_model::Int=DEFAULT_COST_MODEL) where {T<:Unsigned}
	vs = vertices(g)
	gs = convert(UInt64, length(vs))

	if gs > MGS_MAX_SIZE
		error("Input graph cannot have more than 2^40-1 vertices")
	end

	n_bits_v = convert(UInt8, ceil(log(2, gs)))
	V = infer_uint_custom_type(n_bits_v)

	option_flags = encode_bg_params(ref_window_size=ref_window_size,
		copy_blocks=copy_blocks, stop_deltas=stop_deltas,
		lr_split=lr_split, adaptive_header=adaptive_header)
	flag_byte1, flag_byte2 = create_header_flags(:directed, coding_scheme, integer_encoding, option_flags)

	header_bytes = UInt8[
		0x4d, 0x47, 0x53,
		0x03, 0x02,  # v3.2
		flag_byte1, flag_byte2,
		(gs & 0xff), ((gs >> 8) & 0xff), ((gs >> 16) & 0xff), ((gs >> 24) & 0xff), ((gs >> 32) & 0xff)
	]

	bw = BitWriter()

	write_bytes(bw, header_bytes)

	# Build sorted neighbor lists from graph
	nls = Dict{V,Vector{V}}()
	for v in vs
		nls[V(v)] = sort([V(o) for o in outneighbors(g, v)])
	end

	@info("writing BG compressed graph data")
	write_greedy_graph_data(bw, nls, coding_scheme, ref_window_size;
		integer_encoding=integer_encoding, copy_blocks=copy_blocks,
		adaptive_copy=copy_blocks, stop_deltas=stop_deltas,
		compact_copy=true,
		tight_intervals=true, fixwidth_ref=lr_split,
		exact_costing=exact_costing, lr_split=lr_split,
		multi_ref=multi_ref, adaptive_header=adaptive_header,
		cost_model=cost_model)

	flush_bitwriter(bw; flush_last_bits=true)
	open(filename * ".mgz", "w") do f
		write_to_io(bw, f)
	end
end

################################################################################
# CS (Command Stream) compressed MGS v3 graph
################################################################################

"""
    write_cs_mgs3_graph(g, filename; ...)

Write graph in MGS v3 format with CS (Command Stream) compression.
Produces a `.mgz` file with v3.1 header (params encoded in byte2).
"""
function write_cs_mgs3_graph(g::AbstractGraph{T}, filename::AbstractString;
		coding_scheme::Symbol=:children, integer_encoding::Symbol=:fibonacci,
		ref_window_size::Int=64, compact_copy::Bool=true,
		tight_intervals::Bool=true, lr_split::Bool=false,
		cost_model::Int=DEFAULT_COST_MODEL) where {T<:Unsigned}
	vs = vertices(g)
	gs = convert(UInt64, length(vs))

	if gs > MGS_MAX_SIZE
		error("Input graph cannot have more than 2^40-1 vertices")
	end

	n_bits_v = convert(UInt8, ceil(log(2, gs)))
	V = infer_uint_custom_type(n_bits_v)

	# v3.2: encode params into byte2
	option_flags = encode_cs_params(ref_window_size=ref_window_size, compact_copy=compact_copy,
		tight_intervals=tight_intervals, lr_split=lr_split)
	flag_byte1, flag_byte2 = create_header_flags(:directed, coding_scheme, integer_encoding, option_flags)

	header_bytes = UInt8[
		0x4d, 0x47, 0x53,
		0x03, 0x02,             # v3.2
		flag_byte1, flag_byte2,
		(gs & 0xff), ((gs >> 8) & 0xff), ((gs >> 16) & 0xff), ((gs >> 24) & 0xff), ((gs >> 32) & 0xff)
	]

	bw = BitWriter()

	write_bytes(bw, header_bytes)

	# Build sorted neighbor lists from graph
	nls = Dict{V,Vector{V}}()
	for v in vs
		nls[V(v)] = sort([V(o) for o in outneighbors(g, v)])
	end

	@info("writing CS compressed graph data")
	write_cmdstream_graph_data(bw, nls, coding_scheme, ref_window_size;
		integer_encoding=integer_encoding, compact_copy=compact_copy,
		tight_intervals=tight_intervals, lr_split=lr_split,
		cost_model=cost_model)

	flush_bitwriter(bw; flush_last_bits=true)
	open(filename * ".mgz", "w") do f
		write_to_io(bw, f)
	end
end

################################################################################
# CG (Clustered Graph) compressed MGS v3 graph
################################################################################

"""
    write_cg_mgs3_graph(g, filename, clusters; ...)

Write graph in MGS v3 format with CG compression.
Produces a `.mgz` file with self-describing header (params in byte2, varint in integer_encoding).

Parameters:
- g: Input graph
- filename: Output filename (without .mgz extension)
- clusters: Partition as Vector{Vector{T}} of vertex clusters
- params: CGParams encoding parameters (varint stored in header integer_encoding field)
"""
function write_cg_mgs3_graph(g::AbstractGraph{T}, filename::AbstractString,
		clusters::Vector{Vector{T}};
		coding_scheme::Symbol=:children,
		params::CGParams=CGParams(),
		progress::Union{Nothing,Function}=nothing) where {T<:Unsigned}
	vs = vertices(g)
	gs = convert(UInt64, length(vs))

	if gs > MGS_MAX_SIZE
		error("Input graph cannot have more than 2^40-1 vertices")
	end
	option_flags = encode_cg_params(params)
	flag_byte1, flag_byte2 = create_header_flags(:directed, coding_scheme, params.varint, option_flags)

	header_bytes = UInt8[
		0x4d, 0x47, 0x53,
		0x03, 0x02,             # v3.2
		flag_byte1, flag_byte2,
		(gs & 0xff), ((gs >> 8) & 0xff), ((gs >> 16) & 0xff), ((gs >> 24) & 0xff), ((gs >> 32) & 0xff)
	]

	bw = BitWriter()

	write_bytes(bw, header_bytes)

	stats = CGStats()
	if coding_scheme == :index
		# Two-pass encoding: encode to buffer, then write offset table + data
		@info("writing CG compressed graph data (index mode, two-pass)")
		buf_bw = BitWriter()

		# Pass 1: encode to buffer, record cluster start bit offsets
		K = length(clusters)
		cg_offsets = zeros(Int, 2 * K + 1)  # K intra offsets + 1 inter start + K inter per-source offsets
		encode_level(buf_bw, g, clusters; params=params, stats=stats, progress=progress,
			cluster_offsets=cg_offsets)
		flush_bitwriter(buf_bw; flush_last_bits=true)

		# Compute entry width (must cover all offsets including inter per-source)
		max_offset = maximum(cg_offsets[1:(K+1)])  # intra + inter start are absolute
		# Inter per-source offsets are relative; find max among them too
		if K > 0
			max_inter_rel = maximum(cg_offsets[(K+2):(2*K+1)])
			max_offset = max(max_offset, max_inter_rel)
		end
		entry_width = max_offset > 0 ? max(Int(ceil(log2(max_offset + 1))), 1) : 1

		# Write offset table: 6-bit entry_width + 32-bit K + (2K+1) entries
		write_value(bw, UInt64(entry_width), 6)
		write_value(bw, UInt32(K), 32)
		for i in 1:(2 * K + 1)
			write_value(bw, UInt64(cg_offsets[i]), entry_width)
		end

		# Write buffered encoded data
		buf_data = get_bytes(buf_bw)
		write_bytes(bw, collect(buf_data))
	else
		@info("writing CG compressed graph data")
		encode_level(bw, g, clusters; params=params, stats=stats, progress=progress)
	end

	flush_bitwriter(bw; flush_last_bits=true)
	open(filename * ".mgz", "w") do f
		write_to_io(bw, f)
	end

	# Log bit budget breakdown
	m = ne(g)
	total_bits = stats.bits_membership + stats.bits_intra_headers + stats.bits_intra_copy +
	             stats.bits_intra_add + stats.bits_intra_raw +
	             stats.bits_inter_headers + stats.bits_inter_degrees +
	             stats.bits_inter_perms + stats.bits_inter_lists
	@info "CG bit budget ($(m) edges):"
	@info "  membership:     $(stats.bits_membership) bits ($(round(stats.bits_membership/m, digits=4)) BPE)"
	@info "  intra headers:  $(stats.bits_intra_headers) bits ($(round(stats.bits_intra_headers/m, digits=4)) BPE)"
	@info "  intra copy:     $(stats.bits_intra_copy) bits ($(round(stats.bits_intra_copy/m, digits=4)) BPE)"
	@info "  intra add:      $(stats.bits_intra_add) bits ($(round(stats.bits_intra_add/m, digits=4)) BPE)"
	@info "  intra raw:      $(stats.bits_intra_raw) bits ($(round(stats.bits_intra_raw/m, digits=4)) BPE)"
	@info "  inter headers:  $(stats.bits_inter_headers) bits ($(round(stats.bits_inter_headers/m, digits=4)) BPE)"
	@info "  inter degrees:  $(stats.bits_inter_degrees) bits ($(round(stats.bits_inter_degrees/m, digits=4)) BPE)"
	@info "  inter perms:    $(stats.bits_inter_perms) bits ($(round(stats.bits_inter_perms/m, digits=4)) BPE)"
	@info "  inter lists:    $(stats.bits_inter_lists) bits ($(round(stats.bits_inter_lists/m, digits=4)) BPE)"
	@info "  ref used/no_ref: $(stats.intra_ref_used)/$(stats.intra_no_ref)"
	@info "  copy modes: bitmap=$(stats.intra_copy_bitmap_count) blocks=$(stats.intra_copy_blocks_count) complement=$(stats.intra_copy_complement_count)"
	@info "  overlap histogram: $(stats.intra_overlap_hist)"
end

################################################################################
# Load compressed MGS v3 graph
################################################################################

"""
    load_compressed_mgs3_graph(filename::AbstractString)

Load graph in MGS v3 format.

Parameters:
- filename: Input filename

Supported compression schemes:
- :huffman
- :elias_gamma
- :elias_delta
- :golomb
- :fibonacci
- :zeta
- :fed

Returns a graph loaded with the compression scheme specified in the header.
"""
function load_compressed_mgs3_graph(filename::AbstractString)
	f = open(filename, "r")

	# read header (12 bytes)
	header_bytes = read(f, 12)

	# Verify MGS signature
	if header_bytes[1:3] != [0x4d, 0x47, 0x53]  # 'MGS'
		error("Invalid MGS file signature")
	end

	# Extract and decode flags
	flag_byte1 = header_bytes[6]
	flag_byte2 = header_bytes[7]
	graph_type, encoding, compression, byte2 = decode_header_flags(flag_byte1, flag_byte2)

	# Extract number of vertices (5 bytes, little-endian)
	gs_bytes = header_bytes[8:12]
	gs = UInt64(gs_bytes[1]) | (UInt64(gs_bytes[2]) << 8) | (UInt64(gs_bytes[3]) << 16) |
		 (UInt64(gs_bytes[4]) << 24) | (UInt64(gs_bytes[5]) << 32)

	# Self-describing dispatch on byte2
	# Algorithm IDs (0x00–0x0F)
	g = if byte2 == ALG_LEGACY_MGS
		error("ALG_LEGACY_MGS (0x00) should not appear in compressed .mgz files")
	elseif byte2 == ALG_HUFFMAN
		load_huffman_compressed_mgs3_graph(f, graph_type, encoding, gs)
	elseif byte2 == ALG_BG
		p = _bg_default_params()
		load_greedy_mgs3_graph(f, graph_type, encoding, gs, compression;
			copy_blocks=p.copy_blocks, adaptive_copy=p.adaptive_copy,
			ref_window_size=p.ref_window_size, fixwidth_ref=p.fixwidth_ref,
			stop_deltas=p.stop_deltas,
			compact_copy=p.compact_copy, tight_intervals=p.tight_intervals,
			lr_split=p.lr_split, adaptive_header=p.adaptive_header)
	elseif byte2 == ALG_CS
		p = _cs_default_params()
		load_cs_mgs3_graph(f, graph_type, encoding, gs, compression;
			compact_copy=p.compact_copy, tight_intervals=p.tight_intervals,
			ref_window_size=p.ref_window_size, lr_split=p.lr_split)
	elseif byte2 == ALG_CG
		load_cg_mgs3_graph(f, graph_type, gs; params=_cg_default_params(), coding_scheme=encoding)
	elseif byte2 <= 0x0F
		error("Reserved algorithm ID: 0x$(string(byte2, base=16, pad=2))")
	# Parameter ranges (0x10–0xFF)
	elseif PARAM_BG_BASE <= byte2 <= PARAM_BG_MAX
		p = decode_bg_params(byte2)
		load_greedy_mgs3_graph(f, graph_type, encoding, gs, compression;
			copy_blocks=p.copy_blocks, adaptive_copy=p.adaptive_copy,
			ref_window_size=p.ref_window_size, fixwidth_ref=p.fixwidth_ref,
			stop_deltas=p.stop_deltas,
			compact_copy=p.compact_copy, tight_intervals=p.tight_intervals,
			lr_split=p.lr_split, adaptive_header=p.adaptive_header)
	elseif PARAM_CS_BASE <= byte2 <= PARAM_CS_MAX
		p = decode_cs_params(byte2)
		load_cs_mgs3_graph(f, graph_type, encoding, gs, compression;
			compact_copy=p.compact_copy, tight_intervals=p.tight_intervals,
			ref_window_size=p.ref_window_size, lr_split=p.lr_split)
	elseif PARAM_CG_BASE <= byte2 <= PARAM_CG_MAX
		p = decode_cg_params(byte2; varint=compression)
		load_cg_mgs3_graph(f, graph_type, gs; params=p, coding_scheme=encoding)
	else
		error("Unknown byte2: 0x$(string(byte2, base=16, pad=2))")
	end

	close(f)
	return g
end

""" 
    load_huffman_compressed_mgs3_graph(f::IO, graph_type::Symbol, encoding::Symbol, gs::UInt64)

load graph in compressed MGS v3 format (Huffman compression scheme)

Parameters:
- f: Input stream
- graph_type: graph type (:directed or :undirected)
- encoding: coding scheme (:children or :index)
- gs: number of vertices
"""
function load_huffman_compressed_mgs3_graph(io::IO, graph_type::Symbol, encoding::Symbol, gs::UInt64)
	# `n_size_u` is the number of bits needed to represent the graph vertices
	n_bits_v = convert(UInt8, ceil(log(2, gs)))
	# Get appropriate unsigned int type based on number of bits needed
	V = infer_uint_custom_type(n_bits_v)
	
	# intialize graph g according to graph type
	g = graph_type == :directed ? SimpleDiGraph{V}() : SimpleGraph{V}()
	
	# vertex set
	vs = range(1, stop=gs)
	# NB: stop sequence is the code associated to 0
	stop_seq = zero(V)

	@info("generating graph")
	@info("adding vertices")
	# add vertices to graph
    add_vertices!(g, gs)
	
	# frequencies of each vertex (in-degree)
	in_degrees = Dict{V,V}()

	# read frequency section for Huffman decoding
	for v in vs
		p = read(io, sizeof(V))
		in_degrees[v] = reinterpret(V, reverse(p))[1]
	end

	if encoding == :children
		# read stop sequence from frequency section
		p = read(io, sizeof(V))
		stop_seq_value = reinterpret(V, reverse(p))[1]
		in_degrees[stop_seq] = stop_seq_value
	end

	out_degrees = Dict{V,V}()
	if encoding == :index
		@info("reading index section")
		# read index
		for v in vs
			p = read(io, sizeof(V))
			out_degrees[v] = reinterpret(V, reverse(p))[1]
		end
	end

	@info("reading data section")
	# read data section
	cdata = BitVector()
	while !eof(io)
		# read a byte
		b = read(io, 1)[1]
		# read each bit of the byte
		for j in 0:7
			if ((b >> j) & 0x01) == 1
				push!(cdata, true)
			else
				push!(cdata, false)
			end
		end
	end

	# get last byte number of 0s
	sp = 8 - sum(cdata[end-7:end])
	cdata = cdata[1:end-(7+sp)]
	
	@info("generating Huffman tree")
	tree = huffman_encoding(in_degrees)
	
	@info("decoding values")
	children = decode_huffman_values(tree, cdata)

	# number of children
	n_children = length(children)

	@info("adding edges")
	if encoding == :children
		pos = 1
		for v in vs
			source = convert(V, v)
			while pos <= n_children && children[pos] != stop_seq
				target = children[pos]
				add_edge!(g, source, target)
				pos += 1
			end
			# skip stop sequence
			pos += 1
		end
	elseif encoding == :index
		current_pos = 1
		for v in vs
			source = convert(V,v)
			# if we reached the last parent vertex
			if v == length(vs)
				pos1 = current_pos
				pos2 = n_children
			else
				pos1 = current_pos
				# position of the last child of vertex v
				pos2 = current_pos + out_degrees[v] - 1
				current_pos = pos2 + 1
			end
			# add edges for each child
			for p in pos1:pos2
				target = children[p]
				add_edge!(g, source, target)
			end
		end
	end
	close(io)
	return g
end

################################################################################
# Greedy compressed MGS v3 graph
################################################################################

"""
    load_greedy_mgs3_graph(io, graph_type, coding_scheme, gs, integer_encoding; ...)

Load graph from greedy-compressed MGS v3 format. The data stream is self-describing.
"""
function load_greedy_mgs3_graph(io::IO, graph_type::Symbol, coding_scheme::Symbol, gs::UInt64,
		integer_encoding::Symbol=:fibonacci; copy_blocks::Bool=false, adaptive_copy::Bool=false,
		ref_window_size::Int=7, fixwidth_ref::Bool=false, stop_deltas::Bool=false,
		adaptive_deltas::Bool=false, split_residual::Bool=false,
		compact_copy::Bool=false, tight_intervals::Bool=false,
		lr_split::Bool=false, multi_ref::Bool=false,
		adaptive_header::Bool=false)
	n_bits_v = convert(UInt8, ceil(log(2, gs)))
	V = infer_uint_custom_type(n_bits_v)

	g = graph_type == :directed ? SimpleDiGraph{V}() : SimpleGraph{V}()

	@info("generating graph (greedy mode)")
	@info("adding vertices")
	add_vertices!(g, gs)

	reader = BitReader(io)

	@info("reading greedy-compressed graph data")
	neighbor_lists = read_greedy_graph_data(reader, V(gs), coding_scheme, V;
		integer_encoding=integer_encoding, copy_blocks=copy_blocks,
		adaptive_copy=adaptive_copy, ref_window_size=ref_window_size,
		fixwidth_ref=fixwidth_ref, stop_deltas=stop_deltas,
		adaptive_deltas=adaptive_deltas,
		split_residual=split_residual, compact_copy=compact_copy,
		tight_intervals=tight_intervals,
		lr_split=lr_split, multi_ref=multi_ref,
		adaptive_header=adaptive_header)

	@info("building graph from neighbor lists")
	for (source_vertex, neighbors) in neighbor_lists
		for target_vertex in neighbors
			add_edge!(g, source_vertex, target_vertex)
		end
	end

	@info("graph construction completed: $(nv(g)) vertices, $(ne(g)) edges")
	return g
end

################################################################################
# CS (Command Stream) load
################################################################################

"""
    load_cs_mgs3_graph(io, graph_type, coding_scheme, gs, integer_encoding; ...)

Load graph from CS (Command Stream) compressed MGS v3 format.
"""
function load_cs_mgs3_graph(io::IO, graph_type::Symbol, coding_scheme::Symbol, gs::UInt64,
		integer_encoding::Symbol=:fibonacci;
		compact_copy::Bool=true, tight_intervals::Bool=true, ref_window_size::Int=64,
		lr_split::Bool=false)
	n_bits_v = convert(UInt8, ceil(log(2, gs)))
	V = infer_uint_custom_type(n_bits_v)

	g = graph_type == :directed ? SimpleDiGraph{V}() : SimpleGraph{V}()

	@info("generating graph (CS mode)")
	add_vertices!(g, gs)

	reader = BitReader(io)

	@info("reading CS-compressed graph data")
	neighbor_lists = read_cmdstream_graph_data(reader, V(gs), coding_scheme, V;
		integer_encoding=integer_encoding, compact_copy=compact_copy,
		tight_intervals=tight_intervals, ref_window_size=ref_window_size,
		lr_split=lr_split)

	@info("building graph from neighbor lists")
	for (source_vertex, neighbors) in neighbor_lists
		for target_vertex in neighbors
			add_edge!(g, source_vertex, target_vertex)
		end
	end

	@info("graph construction completed: $(nv(g)) vertices, $(ne(g)) edges")
	return g
end

################################################################################
# CG load
################################################################################

"""
    load_cg_mgs3_graph(io, graph_type, gs; params=CGParams())

Load graph from CG compressed MGS v3 format.
"""
function load_cg_mgs3_graph(io::IO, graph_type::Symbol, gs::UInt64;
		params::CGParams=CGParams(), coding_scheme::Symbol=:children)
	n_bits_v = convert(UInt8, ceil(log(2, gs)))
	V = infer_uint_custom_type(n_bits_v)

	directed = graph_type == :directed

	@info("generating graph (CG mode)")

	reader = BitReader(io)

	if coding_scheme == :index
		# Read cluster offset table: 6-bit entry_width + 32-bit K + (2K+1) entries
		entry_width = Int(read_value(reader, 6, UInt64))
		K = Int(read_value(reader, 32, UInt32))
		# Read offset entries: K intra + 1 inter start + K inter per-source
		_cg_offsets = [Int(read_value(reader, entry_width, UInt64)) for _ in 1:(2 * K + 1)]
		@info "CG index mode: K=$K, entry_width=$entry_width, intra_offsets=$(_cg_offsets[1:K+1]), inter_offsets=$(_cg_offsets[K+2:2K+1])"
		# Fall through to normal sequential decode
	end

	@info("reading CG-compressed graph data")
	neighbor_lists = decode_level(reader, params; T=V, directed=directed, coding_scheme=coding_scheme)

	# Build graph from decoded neighbor lists
	g = directed ? SimpleDiGraph{V}() : SimpleGraph{V}()
	add_vertices!(g, gs)

	@info("building graph from neighbor lists")
	for (source_vertex, neighbors) in neighbor_lists
		for target_vertex in neighbors
			add_edge!(g, source_vertex, target_vertex)
		end
	end

	@info("graph construction completed: $(nv(g)) vertices, $(ne(g)) edges")
	return g
end

end # module MGS
