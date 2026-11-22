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
import ..Compression
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
	write_greedy_graph_data, read_greedy_graph_data, parse_bg_ctxrange_index,
	write_cmdstream_graph_data, read_cmdstream_graph_data, IndexedRangeState,
	_recycle_indexed_state!,
	COST_MODEL_FULL, COST_MODEL_FAST, DEFAULT_COST_MODEL

using ..Compression.CG: encode_level, decode_level, CGParams, CGStats, _read_membership

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
const INT_ENCODING_CONTEXT_RANGE = 0x7
#
# ── Self-describing byte 2 layout ───────────────────────────────────────────
# Algorithm IDs (0x00–0x0F): no additional params needed
const ALG_LEGACY_MGS = 0x00
const ALG_HUFFMAN    = 0x01
const ALG_BG        = 0x02   # BG with recommended defaults
const ALG_CS         = 0x03   # CS with recommended defaults
const ALG_CG        = 0x04   # CG with recommended defaults
const ALG_CAND      = 0x05   # candidate-bit coder (params stored after the header)
# 0x06–0x0F reserved

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
    elseif integer_encoding == :context_range
        INT_ENCODING_CONTEXT_RANGE
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
    elseif integer_encoding_bits == INT_ENCODING_CONTEXT_RANGE
        :context_range
    else
        error("Unknown integer encoding: $integer_encoding_bits")
    end
    
    # Decode option flags
    option_flags = byte2
    
    return (graph_type, coding_scheme, integer_encoding, option_flags)
end

# ── Header extension (v3.x + extension flag) ─────────────────────────────────
# Files that need no extension keep the plain 12-byte header and are byte-identical
# to earlier releases. A file that needs one sets the high bit of the minor-version
# byte (the low 7 bits keep the stream-layout version) and is followed directly by
#   varint(total_len) | { tag:UInt8 | varint(len) | payload[len] }*
# Older readers reject such files (unexpected minor byte) instead of misreading them.
# Tags are mandatory: a reader that meets an unknown tag must refuse the file.
const MGS_MINOR_EXT_FLAG = 0x80
const MGS_EXT_CG_WINDOW  = 0x01   # payload: varint(window); overrides the window slot in byte2
const MGS_EXT_SYMMETRIC  = 0x02   # payload: empty; rows hold only targets >= source; mirror on decode
const MGS_EXT_RECIPROCAL = 0x03   # payload: range-coded reciprocity bits (see `_reciprocal_split`)
const MGS_EXT_TRANSPOSED = 0x04   # payload: empty; the stored graph is the transpose (rows = in-lists); transpose on decode
const MGS_EXT_KNOWN_TAGS = (MGS_EXT_CG_WINDOW, MGS_EXT_SYMMETRIC, MGS_EXT_RECIPROCAL, MGS_EXT_TRANSPOSED)

function _mgs_varint_bytes(x::Integer)
	x >= 0 || throw(ArgumentError("varint must be nonnegative"))
	out = UInt8[]; v = UInt64(x)
	while true
		b = UInt8(v & 0x7f); v >>= 7
		if v == 0; push!(out, b); return out; end
		push!(out, b | 0x80)
	end
end
function _mgs_read_varint(io::IO)
	v = UInt64(0); shift = 0
	while true
		b = read(io, UInt8)
		v |= UInt64(b & 0x7f) << shift
		(b & 0x80) == 0 && return v
		shift += 7
		shift <= 63 || throw(ArgumentError("MGS extension varint overflow"))
	end
end

"""
    _mgs_ext_block(fields) -> Vector{UInt8}

Serialise extension fields (`tag => payload` pairs, written in the given order).
"""
function _mgs_ext_block(fields::AbstractVector{<:Pair})
	body = UInt8[]
	for (tag, payload) in fields
		UInt8(tag) in MGS_EXT_KNOWN_TAGS || throw(ArgumentError("unknown MGS extension tag $tag"))
		push!(body, UInt8(tag)); append!(body, _mgs_varint_bytes(length(payload))); append!(body, payload)
	end
	return vcat(_mgs_varint_bytes(length(body)), body)
end

"""
    _read_mgs_header(io) -> (header::Vector{UInt8}, ext::Dict{UInt8,Vector{UInt8}})

Read the 12-byte header and, when the extension flag is set, the extension block.
`header[5]` is returned with the flag cleared (the stream-layout minor version), and
`io` is left at the first byte after the header (and extension).
"""
function _read_mgs_header(io::IO)
	h = read(io, 12)
	length(h) == 12 || throw(EOFError())
	ext = Dict{UInt8,Vector{UInt8}}()
	if (h[5] & MGS_MINOR_EXT_FLAG) != 0
		h[5] &= ~MGS_MINOR_EXT_FLAG
		total = Int(_mgs_read_varint(io))
		body = read(io, total)
		length(body) == total || throw(EOFError())
		b = IOBuffer(body)
		while !eof(b)
			tag = read(b, UInt8)
			tag in MGS_EXT_KNOWN_TAGS || throw(ArgumentError("Unsupported MGS header extension tag 0x$(string(tag, base=16, pad=2))"))
			len = Int(_mgs_read_varint(b))
			payload = read(b, len)
			length(payload) == len || throw(EOFError())
			haskey(ext, tag) && throw(ArgumentError("Duplicate MGS header extension tag"))
			ext[tag] = payload
		end
	end
	return h, ext
end

_mgs_ext_symmetric(ext) = haskey(ext, MGS_EXT_SYMMETRIC)
function _mgs_ext_cg_window(ext)
	haskey(ext, MGS_EXT_CG_WINDOW) || return nothing
	w = Int(_mgs_read_varint(IOBuffer(ext[MGS_EXT_CG_WINDOW])))
	w >= 1 || throw(ArgumentError("invalid CG window in header extension"))
	return w
end
# Apply a CG window override from the extension (no-op when absent).
_apply_cg_ext(p, ext) = (w = _mgs_ext_cg_window(ext); w === nothing ? p : _with_cg_window(p, w))
function _with_cg_window(p::CGParams, w::Int)
	vals = ntuple(i -> fieldname(CGParams, i) === :intra_ref_window ? w : getfield(p, i), fieldcount(CGParams))
	return CGParams(vals...)
end
_with_cg_window(p::NamedTuple, w::Int) = merge(p, (intra_ref_window=w,))

"""
    _finalize_mgs_ext!(path, fields)

Rewrite a freshly written `.mgz` so that it carries the given extension fields:
set the extension flag in the minor byte and insert the block after the header.
No-op for an empty field list, so ordinary files are never touched.
"""
function _finalize_mgs_ext!(path::AbstractString, fields::AbstractVector{<:Pair})
	isempty(fields) && return nothing
	data = read(path)
	length(data) >= 12 && data[1:3] == UInt8[0x4d, 0x47, 0x53] || error("not an MGS file: $path")
	(data[5] & MGS_MINOR_EXT_FLAG) == 0 || error("MGS header extension already present: $path")
	h = data[1:12]; h[5] |= MGS_MINOR_EXT_FLAG
	tmp = path * ".tmp_ext"
	open(tmp, "w") do f
		write(f, h); write(f, _mgs_ext_block(fields)); write(f, view(data, 13:length(data)))
	end
	mv(tmp, path; force=true)
	return nothing
end

# ── Symmetric mode ───────────────────────────────────────────────────────────
# A graph whose arc set is symmetric (u->v iff v->u) is stored as its upper
# triangle (targets >= source, so self-loops are kept once) with the
# MGS_EXT_SYMMETRIC extension; decoders mirror every off-diagonal arc.
function _is_symmetric_arcset(g::AbstractGraph)
	for u in vertices(g), v in outneighbors(g, u)
		u == v && continue
		has_edge(g, v, u) || return false
	end
	return true
end
function _resolve_symmetric(g::AbstractGraph, symmetric)
	symmetric === false || symmetric === :off ? false :
	symmetric === true || symmetric === :on ?
		(_is_symmetric_arcset(g) || throw(ArgumentError("symmetric=:on requires a symmetric arc set")); true) :
	symmetric === :auto ? _is_symmetric_arcset(g) :
	throw(ArgumentError("symmetric must be :off, :auto or :on"))
end
function _upper_triangle_graph(g::AbstractGraph{T}) where {T}
	h = SimpleDiGraph{T}(nv(g))
	for u in vertices(g), v in outneighbors(g, u)
		v >= u && add_edge!(h, u, v)
	end
	return h
end
# Mirror decoded upper-triangle rows into full symmetric rows (sorted, unique).
function _mirror_upper_rows!(rows::AbstractVector{<:AbstractVector{T}}) where {T}
	n = length(rows)
	extra = [T[] for _ in 1:n]
	for u in 1:n, v in rows[u]
		Int(v) > u && push!(extra[Int(v)], T(u))
	end
	for v in 1:n
		isempty(extra[v]) && continue
		r = vcat(rows[v], extra[v]); sort!(r); unique!(r)
		rows[v] = r
	end
	return rows
end
# ── Reciprocal mode ──────────────────────────────────────────────────────────
# Rows are decoded in id order, so for row u every earlier v with a forward arc
# v->u (v < u) is a known back-candidate. The stored graph G' is G minus every
# backward arc u->v (v < u) whose reverse v->u exists. Forward arcs are never
# removed, so the decoder rebuilds the candidate lists from G' alone and reads
# one adaptive bit per (u, v) candidate saying whether u->v is present.
const _RECIP_NCTX = 3 * 5 * 16
@inline _recip_lb(x::Integer) = x <= 0 ? 0 : min(15, 64 - leading_zeros(UInt64(x)))
@inline _recip_ctx(prev::Int, hits::Int, tested::Int, nbb::Int) =
	1 + prev + 3 * (tested == 0 ? 0 : 1 + min(3, (4 * hits) ÷ (tested + 1))) + 15 * nbb
@inline function _insorted(r::AbstractVector{T}, x::T) where {T}
	i = searchsortedfirst(r, x); i <= length(r) && r[i] == x
end
function _back_candidates(rows::AbstractVector{<:AbstractVector{T}}) where {T}
	n = length(rows); B = [T[] for _ in 1:n]
	for v in 1:n, u in rows[v]
		Int(u) > v && push!(B[Int(u)], T(v))
	end
	return B
end

"""
    _reciprocal_split(g) -> (gp, bytes, removed)

Split `g` into the graph that is actually stored (`gp`, without reciprocated
backward arcs) and the range-coded reciprocity bits.
"""
function _reciprocal_split(g::AbstractGraph{T}) where {T}
	n = nv(g)
	a = [sort!(collect(T, outneighbors(g, u))) for u in 1:n]
	gp = SimpleDiGraph{T}(n); removed = 0
	for u in 1:n, v in a[u]
		if Int(v) < u && _insorted(a[Int(v)], T(u))
			removed += 1
		else
			add_edge!(gp, T(u), v)
		end
	end
	B = _back_candidates(a)
	e = Compression.BinRangeEncoder()
	c0 = ones(UInt32, _RECIP_NCTX); c1 = ones(UInt32, _RECIP_NCTX)
	for u in 1:n
		Bu = B[u]; isempty(Bu) && continue
		prev = 0; hits = 0; tested = 0; nbb = _recip_lb(length(Bu))
		for v in Bu
			b = _insorted(a[u], v)
			Compression.brc_encode_ctx_bit!(e, c0, c1, _recip_ctx(prev, hits, tested, nbb), b)
			prev = b ? 2 : 1; tested += 1; hits += b
		end
	end
	return gp, Compression.brc_finish!(e), removed
end

# Restore reciprocated arcs into decoded rows (sorted unique on input and output).
function _reciprocal_restore!(rows::AbstractVector{<:AbstractVector{T}}, bytes::Vector{UInt8}) where {T}
	n = length(rows)
	B = _back_candidates(rows)
	d = Compression.BinRangeDecoder(bytes)
	c0 = ones(UInt32, _RECIP_NCTX); c1 = ones(UInt32, _RECIP_NCTX)
	for u in 1:n
		Bu = B[u]; isempty(Bu) && continue
		prev = 0; hits = 0; tested = 0; nbb = _recip_lb(length(Bu)); add = T[]
		for v in Bu
			b = Compression.brc_decode_ctx_bit!(d, c0, c1, _recip_ctx(prev, hits, tested, nbb))
			b && push!(add, v)
			prev = b ? 2 : 1; tested += 1; hits += b
		end
		if !isempty(add)
			r = vcat(rows[u], add); sort!(r); rows[u] = r
		end
	end
	return rows
end
function _reciprocal_restore_graph!(g::AbstractGraph{T}, bytes::Vector{UInt8}) where {T}
	rows = [sort!(collect(T, outneighbors(g, u))) for u in vertices(g)]
	lens = length.(rows)
	_reciprocal_restore!(rows, bytes)
	for u in 1:length(rows)
		length(rows[u]) == lens[u] && continue
		for v in rows[u]; has_edge(g, T(u), v) || add_edge!(g, T(u), v); end
	end
	return g
end

# ── Transposed storage ───────────────────────────────────────────────────────
function _transpose_graph(g::AbstractGraph{T}) where {T}
	h = SimpleDiGraph{T}(nv(g))
	for u in vertices(g), v in outneighbors(g, u); add_edge!(h, v, u); end
	return h
end
function _transpose_rows(rows::AbstractVector{<:AbstractVector{T}}) where {T}
	n = length(rows); t = [T[] for _ in 1:n]
	for u in 1:n, v in rows[u]; push!(t[Int(v)], T(u)); end   # u ascending => rows sorted
	return t
end

function _mirror_upper_graph!(g::AbstractGraph)
	arcs = [(u, v) for u in vertices(g) for v in outneighbors(g, u) if v > u]
	for (u, v) in arcs; add_edge!(g, v, u); end
	return g
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
	params.intra_ref_window >= 1 || error("Unsupported CG window: $(params.intra_ref_window)")
	mil = params.intra_adapt_mil  # encode adapt_mil (= intra_mil in practice)
	haskey(CG_MIL_TO_IDX, mil) || error("Unsupported CG mil: $mil (must be 3, 4, or 5)")

	midx = CG_V32_MEMBERSHIP_TO_IDX[params.membership]
	# A window outside the 6-slot table is stored in slot 256 and overridden by the
	# MGS_EXT_CG_WINDOW header extension (see `_cg_window_ext_fields`).
	widx = get(CG_WINDOW_TO_IDX, params.intra_ref_window, CG_WINDOW_TO_IDX[256])
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

# Extension fields a CG file needs: none when the window fits the byte2 table.
_cg_window_ext_fields(params::CGParams) = haskey(CG_WINDOW_TO_IDX, params.intra_ref_window) ? Pair{UInt8,Vector{UInt8}}[] :
	[MGS_EXT_CG_WINDOW => _mgs_varint_bytes(params.intra_ref_window)]

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
       load_adjacency_mgs3_graph,
       write_huffman_compressed_mgs3_graph,
	   load_huffman_compressed_mgs3_graph,
	   load_greedy_mgs3_graph,
	   write_bg_mgs3_graph,
	   write_cs_mgs3_graph,
	   write_cg_mgs3_graph, write_cand_mgs3_graph,
	   load_indexed_mgs3_graph

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

`chunk_local_refs=true` restricts references to the current sampled context-range
chunk. This optional encoding policy trades compression for bounded access work;
it changes no file syntax and old readers can decode the result. Default: false.

`parallel_search=true` parallelizes FULL/FAST cost search in children mode only,
then emits decisions sequentially in bounded batches. `search_workers` defaults
to the Julia default thread-pool size; `search_batch_size` defaults to 4096.
Launch Julia with `--threads=N`. Bytes and decision tie-breaking match serial
encoding. Indexed mode and `exact_costing=true` trial writes are not supported
by this opt-in path. Concurrent top-level encoder calls remain unsupported
because emission uses shared stream hooks.
"""
function _write_bg_mgs3_graph_core(g::AbstractGraph{T}, filename::AbstractString;
		coding_scheme::Symbol=:children, integer_encoding::Symbol=:fibonacci,
		ref_window_size::Int=64, copy_blocks::Bool=true,
		stop_deltas::Bool=true,
		lr_split::Bool=false, exact_costing::Bool=false,
		multi_ref::Bool=false,
		adaptive_header::Bool=false,
		cost_model::Int=DEFAULT_COST_MODEL,
		index_sample_k::Int=0, chunk_local_refs::Bool=false,
		parallel_search::Bool=false,
		search_workers::Int=Threads.nthreads(:default),
		search_batch_size::Int=4096) where {T<:Unsigned}
	Compression._check_parallel_search(parallel_search, coding_scheme, exact_costing, search_workers, search_batch_size)
	chunk_local_refs && !(coding_scheme == :index && integer_encoding == :context_range) &&
		throw(ArgumentError("chunk_local_refs requires context-range index mode"))
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

	# Build sorted neighbor lists from graph
	nls = Dict{V,Vector{V}}()
	for v in vs
		nls[V(v)] = sort([V(o) for o in outneighbors(g, v)])
	end

	@info("writing BG compressed graph data")
	if integer_encoding == :context_range
		# Two-pass: structural bitstream + range-coded side streams.
		# Children (2b, v3.3): [12 header][8-byte lens ×5: struct, refdist, copy,
		#   cmd, flag][struct][refdist][copy][cmd][flag][resid]
		# Index/RA (legacy):   [12 header][8-byte lens ×3][struct][refdist][copy][resid]
		sbw = BitWriter()
		resid_bytes, refdist_bytes, copy_bytes, cmd_bytes, flag_bytes =
			write_greedy_graph_data(sbw, nls, coding_scheme, ref_window_size;
				integer_encoding=:context_range, copy_blocks=copy_blocks,
				adaptive_copy=copy_blocks, stop_deltas=stop_deltas,
				compact_copy=true,
				tight_intervals=true, fixwidth_ref=lr_split,
				exact_costing=exact_costing, lr_split=lr_split,
				multi_ref=multi_ref, adaptive_header=adaptive_header,
				cost_model=cost_model, index_sample_k=index_sample_k,
				chunk_local_refs=chunk_local_refs, parallel_search=parallel_search,
				search_workers=search_workers, search_batch_size=search_batch_size)
		flush_bitwriter(sbw; flush_last_bits=true)
		struct_bytes = collect(get_bytes(sbw))
		header_bytes[5] = 0x03  # v3.3: all BG ctx files use the 5-stream layout
		open(filename * ".mgz", "w") do f
			write(f, UInt8.(header_bytes))
			for L in (UInt64(length(struct_bytes)), UInt64(length(refdist_bytes)), UInt64(length(copy_bytes)),
			          UInt64(length(cmd_bytes)), UInt64(length(flag_bytes)))
				for i in 0:7; write(f, UInt8((L >> (8*i)) & 0xff)); end
			end
			write(f, struct_bytes)
			write(f, refdist_bytes)
			write(f, copy_bytes)
			write(f, cmd_bytes)
			write(f, flag_bytes)
			write(f, resid_bytes)
		end
		return
	end

	bw = BitWriter()
	write_bytes(bw, header_bytes)
	write_greedy_graph_data(bw, nls, coding_scheme, ref_window_size;
		integer_encoding=integer_encoding, copy_blocks=copy_blocks,
		adaptive_copy=copy_blocks, stop_deltas=stop_deltas,
		compact_copy=true,
		tight_intervals=true, fixwidth_ref=lr_split,
		exact_costing=exact_costing, lr_split=lr_split,
		multi_ref=multi_ref, adaptive_header=adaptive_header,
		cost_model=cost_model, index_sample_k=index_sample_k,
		parallel_search=parallel_search, search_workers=search_workers,
		search_batch_size=search_batch_size)

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

`chunk_local_refs=true` restricts references to the current sampled context-range
chunk without changing the file syntax. It may increase compressed size and is
disabled by default. It requires `coding_scheme=:index, integer_encoding=:context_range`.

`parallel_search=true`, `search_workers`, and `search_batch_size` have the same
children-only, byte-preserving behavior as `write_bg_mgs3_graph`. Stream emission
is serial; this does not make concurrent top-level encoder calls safe.
"""
function _write_cs_mgs3_graph_core(g::AbstractGraph{T}, filename::AbstractString;
		coding_scheme::Symbol=:children, integer_encoding::Symbol=:fibonacci,
		ref_window_size::Int=64, compact_copy::Bool=true,
		tight_intervals::Bool=true, lr_split::Bool=false,
		cost_model::Int=DEFAULT_COST_MODEL,
		index_sample_k::Int=0, chunk_local_refs::Bool=false,
		parallel_search::Bool=false,
		search_workers::Int=Threads.nthreads(:default),
		search_batch_size::Int=4096) where {T<:Unsigned}
	Compression._check_parallel_search(parallel_search, coding_scheme, false, search_workers, search_batch_size)
	chunk_local_refs && !(coding_scheme == :index && integer_encoding == :context_range) &&
		throw(ArgumentError("chunk_local_refs requires context-range index mode"))
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

	# Build sorted neighbor lists from graph
	nls = Dict{V,Vector{V}}()
	for v in vs
		nls[V(v)] = sort([V(o) for o in outneighbors(g, v)])
	end

	@info("writing CS compressed graph data")
	if integer_encoding == :context_range
		# Two-pass: structural bitstream + range-coded side streams.
		# Children (2b, v3.3): [12 header][8-byte lens ×5: struct, refdist, copy,
		#   cmd, flag][struct][refdist][copy][cmd][flag][resid]
		# Index/RA (legacy):   [12 header][8-byte lens ×3][struct][refdist][copy][resid]
		sbw = BitWriter()
		resid_bytes, refdist_bytes, copy_bytes, cmd_bytes, flag_bytes =
			write_cmdstream_graph_data(sbw, nls, coding_scheme, ref_window_size;
				integer_encoding=:context_range, compact_copy=compact_copy,
				tight_intervals=tight_intervals, lr_split=lr_split,
				cost_model=cost_model, index_sample_k=index_sample_k,
				chunk_local_refs=chunk_local_refs, parallel_search=parallel_search,
				search_workers=search_workers, search_batch_size=search_batch_size)
		flush_bitwriter(sbw; flush_last_bits=true)
		struct_bytes = collect(get_bytes(sbw))
		header_bytes[5] = 0x03  # v3.3: all CS ctx files use the 5-stream layout
		open(filename * ".mgz", "w") do f
			write(f, UInt8.(header_bytes))
			for L in (UInt64(length(struct_bytes)), UInt64(length(refdist_bytes)), UInt64(length(copy_bytes)),
			          UInt64(length(cmd_bytes)), UInt64(length(flag_bytes)))
				for i in 0:7; write(f, UInt8((L >> (8*i)) & 0xff)); end
			end
			write(f, struct_bytes); write(f, refdist_bytes); write(f, copy_bytes)
			write(f, cmd_bytes); write(f, flag_bytes)
			write(f, resid_bytes)
		end
		return
	end

	bw = BitWriter()
	write_bytes(bw, header_bytes)
	write_cmdstream_graph_data(bw, nls, coding_scheme, ref_window_size;
		integer_encoding=integer_encoding, compact_copy=compact_copy,
		tight_intervals=tight_intervals, lr_split=lr_split,
		cost_model=cost_model, index_sample_k=index_sample_k,
		parallel_search=parallel_search, search_workers=search_workers,
		search_batch_size=search_batch_size)

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

`parallel_search=true` evaluates per-vertex analytical costs concurrently, even
with one cluster, then writes the same decisions in order. `search_workers` is
capped by Julia's default thread pool (`--threads=N`). FULL stays FULL. Both
children and indexed modes are supported; MGS/block trial paths are excluded.
Decision storage remains cluster-sized, as in serial CG. Progress callbacks run
in order on the caller after each cluster's parallel search, not live on workers.
Concurrent top-level encoding calls remain unsupported (shared emission hooks).

`params.optimize_costs=true` enables byte-identical allocation-free analytical
cost scans and empty K=1 preparation removal. It defaults to false, is not an
encoded parameter, and can be combined with `parallel_search=true`.
"""
function _write_cg_mgs3_graph_core(g::AbstractGraph{T}, filename::AbstractString,
		clusters::Vector{Vector{T}};
		coding_scheme::Symbol=:children,
		params::CGParams=CGParams(),
		integer_encoding::Symbol=params.varint,
		progress::Union{Nothing,Function}=nothing,
		parallel_search::Bool=false,
		search_workers::Int=Threads.nthreads(:default)) where {T<:Unsigned}
	Compression.CG._check_parallel_search(params, parallel_search, search_workers)
	vs = vertices(g)
	gs = convert(UInt64, length(vs))

	if gs > MGS_MAX_SIZE
		error("Input graph cannot have more than 2^40-1 vertices")
	end
	ctx_range = integer_encoding == :context_range
	option_flags = encode_cg_params(params)
	hdr_ie = ctx_range ? :context_range : params.varint
	flag_byte1, flag_byte2 = create_header_flags(:directed, coding_scheme, hdr_ie, option_flags)

	header_bytes = UInt8[
		0x4d, 0x47, 0x53,
		0x03, 0x02,             # v3.2
		flag_byte1, flag_byte2,
		(gs & 0xff), ((gs >> 8) & 0xff), ((gs >> 16) & 0xff), ((gs >> 24) & 0xff), ((gs >> 32) & 0xff)
	]

	if ctx_range
		# Two-pass: structural bitstream + three range-coded side streams.
		# Layout: [12 header][8-byte lens ×3: struct, refdist, copy][struct][refdist][copy][resid]
		# Index (CG-RA): struct = [cluster offset table][chunk-offset tables ×3][data];
		# the three range streams are per-intra-cluster chunked (cluster = chunk).
		sbw = BitWriter()
		stats = CGStats()
		local resid_bytes, refdist_bytes, copy_bytes
		if coding_scheme == :index
			@info("writing CG compressed graph data (context_range, index/RA mode)")
			K = length(clusters)
			cg_offsets = zeros(Int, 2 * K + 1)
			buf_bw = BitWriter()
			resid_bytes, refdist_bytes, copy_bytes, rc_off, rd_off, cp_off =
				encode_level(buf_bw, g, clusters; params=params, stats=stats, progress=progress,
					ctx_range=true, cluster_offsets=cg_offsets, parallel_search=parallel_search, search_workers=search_workers)
			flush_bitwriter(buf_bw; flush_last_bits=true)
			max_offset = maximum(cg_offsets[1:(K + 1)])
			if K > 0
				max_offset = max(max_offset, maximum(cg_offsets[(K + 2):(2 * K + 1)]))
			end
			entry_width = max_offset > 0 ? max(Int(ceil(log2(max_offset + 1))), 1) : 1
			# Cluster offset table (same shape as the fibonacci index writer)
			write_value(sbw, UInt64(entry_width), 6)
			write_value(sbw, UInt32(K), 32)
			for i in 1:(2 * K + 1)
				write_value(sbw, UInt64(cg_offsets[i]), entry_width)
			end
			# Per-cluster chunk byte-offset tables for the three range streams
			max_boff = max(rc_off[end], rd_off[end], cp_off[end])
			coff_width = max_boff > 0 ? max(Int(ceil(log2(max_boff + 1))), 1) : 1
			write_value(sbw, UInt64(coff_width), 6)
			for offs in (rc_off, rd_off, cp_off), o in offs
				write_value(sbw, UInt64(o), coff_width)
			end
			write_bytes(sbw, collect(get_bytes(buf_bw)))
		else
			@info("writing CG compressed graph data (context_range)")
			resid_bytes, refdist_bytes, copy_bytes = encode_level(sbw, g, clusters; params=params, stats=stats, progress=progress, ctx_range=true, parallel_search=parallel_search, search_workers=search_workers)
		end
		flush_bitwriter(sbw; flush_last_bits=true)
		struct_bytes = collect(get_bytes(sbw))
		open(filename * ".mgz", "w") do f
			write(f, UInt8.(header_bytes))
			for L in (UInt64(length(struct_bytes)), UInt64(length(refdist_bytes)), UInt64(length(copy_bytes)))
				for i in 0:7; write(f, UInt8((L >> (8*i)) & 0xff)); end
			end
			write(f, struct_bytes); write(f, refdist_bytes); write(f, copy_bytes); write(f, resid_bytes)
		end
		return
	end

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
			cluster_offsets=cg_offsets, parallel_search=parallel_search, search_workers=search_workers)
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
		encode_level(bw, g, clusters; params=params, stats=stats, progress=progress, parallel_search=parallel_search, search_workers=search_workers)
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


# ── Public writers: symmetric/reciprocal modes + header extension around the core writers ──
# `symmetric` = :off (default, unchanged behaviour), :auto (use it when the arc set is
# symmetric) or :on (require it).
# `reciprocal` = :off (default), :on, or :auto (encode with and without, keep the smaller
# file; doubles encode time). It is skipped when symmetric mode is used, since the
# upper triangle has no reciprocated arcs.
# Both modes are sequential (:children) only.
function _symmetric_input(g, symmetric, coding_scheme)
	sym = _resolve_symmetric(g, symmetric)
	sym && coding_scheme != :children && throw(ArgumentError("symmetric mode supports coding_scheme=:children only"))
	return sym, (sym ? _upper_triangle_graph(g) : g)
end
_sym_fields(sym) = sym ? Pair{UInt8,Vector{UInt8}}[MGS_EXT_SYMMETRIC => UInt8[]] : Pair{UInt8,Vector{UInt8}}[]

function _write_with_modes(core!, g, filename, coding_scheme, symmetric, reciprocal, base_fields; transpose=:off)
	transpose in (:off, :on, :auto, false, true) || throw(ArgumentError("transpose must be :off, :auto or :on"))
	tmode = transpose === true ? :on : transpose === false ? :off : transpose
	tmode == :off && return _write_with_modes_dir(core!, g, filename, coding_scheme, symmetric, reciprocal, base_fields)
	coding_scheme == :children || throw(ArgumentError("transposed storage supports coding_scheme=:children only"))
	tfields = vcat(base_fields, [MGS_EXT_TRANSPOSED => UInt8[]])
	if tmode == :on
		return _write_with_modes_dir(core!, _transpose_graph(g), filename, coding_scheme, symmetric, reciprocal, tfields)
	end
	# :auto — symmetric arc sets are their own transpose
	_is_symmetric_arcset(g) && return _write_with_modes_dir(core!, g, filename, coding_scheme, symmetric, reciprocal, base_fields)
	r = _write_with_modes_dir(core!, g, filename, coding_scheme, symmetric, reciprocal, base_fields)
	tmp = filename * ".transpose_tmp"
	r2 = _write_with_modes_dir(core!, _transpose_graph(g), tmp, coding_scheme, symmetric, reciprocal, tfields)
	if filesize(tmp * ".mgz") < filesize(filename * ".mgz")
		mv(tmp * ".mgz", filename * ".mgz"; force=true); return r2
	end
	rm(tmp * ".mgz"; force=true)
	return r
end

function _write_with_modes_dir(core!, g, filename, coding_scheme, symmetric, reciprocal, base_fields)
	reciprocal in (:off, :on, :auto, false, true) || throw(ArgumentError("reciprocal must be :off, :auto or :on"))
	rmode = reciprocal === true ? :on : reciprocal === false ? :off : reciprocal
	sym, gw = _symmetric_input(g, symmetric, coding_scheme)
	sym && rmode == :on && throw(ArgumentError("reciprocal=:on cannot be combined with symmetric mode"))
	rmode != :off && coding_scheme != :children && throw(ArgumentError("reciprocal mode supports coding_scheme=:children only"))
	sym && (rmode = :off)
	fields0 = vcat(base_fields, _sym_fields(sym))
	path = filename * ".mgz"
	if rmode != :on
		r = core!(gw, filename)
		_finalize_mgs_ext!(path, fields0)
		rmode == :off && return r
	end
	gp, bits, _ = _reciprocal_split(gw)
	if rmode == :on
		r = core!(gp, filename)
		_finalize_mgs_ext!(path, vcat(fields0, [MGS_EXT_RECIPROCAL => bits]))
		return r
	end
	tmp = filename * ".recip_tmp"
	r2 = core!(gp, tmp)
	_finalize_mgs_ext!(tmp * ".mgz", vcat(fields0, [MGS_EXT_RECIPROCAL => bits]))
	if filesize(tmp * ".mgz") < filesize(path)
		mv(tmp * ".mgz", path; force=true); return r2
	end
	rm(tmp * ".mgz"; force=true)
	return r
end

function write_bg_mgs3_graph(g::AbstractGraph, filename::AbstractString; symmetric=:off, reciprocal=:off, transpose=:off,
		coding_scheme::Symbol=:children, kw...)
	_write_with_modes((h, f) -> _write_bg_mgs3_graph_core(h, f; coding_scheme=coding_scheme, kw...),
		g, filename, coding_scheme, symmetric, reciprocal, Pair{UInt8,Vector{UInt8}}[]; transpose=transpose)
end
function write_cs_mgs3_graph(g::AbstractGraph, filename::AbstractString; symmetric=:off, reciprocal=:off, transpose=:off,
		coding_scheme::Symbol=:children, kw...)
	_write_with_modes((h, f) -> _write_cs_mgs3_graph_core(h, f; coding_scheme=coding_scheme, kw...),
		g, filename, coding_scheme, symmetric, reciprocal, Pair{UInt8,Vector{UInt8}}[]; transpose=transpose)
end
function write_cg_mgs3_graph(g::AbstractGraph, filename::AbstractString, clusters; symmetric=:off, reciprocal=:off,
		transpose=:off, coding_scheme::Symbol=:children, params::CGParams=CGParams(), kw...)
	_write_with_modes((h, f) -> _write_cg_mgs3_graph_core(h, f, clusters; coding_scheme=coding_scheme, params=params, kw...),
		g, filename, coding_scheme, symmetric, reciprocal, _cg_window_ext_fields(params); transpose=transpose)
end
"""
    _write_cand_mgs3_graph_core(g, filename; coding_scheme=:children, params=CandParams())

Candidate-bit coder (`Compression.cand_encode`), MGS algorithm id `ALG_CAND` (0x05). Layout: the 12-byte header
(context-range backend, children mode), the optional extension block, `varint(len)` + `cand_params_bytes(params)`,
then the coder payload. Children mode only; vertex ids must fit in UInt32.
"""
function _write_cand_mgs3_graph_core(g::AbstractGraph, filename::AbstractString;
		coding_scheme::Symbol=:children, params::Compression.CandParams=Compression.CandParams())
	coding_scheme == :children || throw(ArgumentError("the candidate coder supports coding_scheme=:children only"))
	gs = UInt64(nv(g))
	gs > MGS_MAX_SIZE && error("Input graph cannot have more than 2^40-1 vertices")
	gs < typemax(UInt32) || throw(ArgumentError("the candidate coder needs vertex ids below 2^32"))
	rows = [sort!(UInt32[UInt32(x) for x in outneighbors(g, v)]) for v in vertices(g)]
	pb = Compression.cand_params_bytes(params)
	body = Compression.cand_encode(rows; params=params)
	fb1, fb2 = create_header_flags(:directed, :children, :context_range, ALG_CAND)
	header = UInt8[0x4d, 0x47, 0x53, 0x03, 0x02, fb1, fb2,
		(gs & 0xff), ((gs >> 8) & 0xff), ((gs >> 16) & 0xff), ((gs >> 24) & 0xff), ((gs >> 32) & 0xff)]
	open(filename * ".mgz", "w") do f
		write(f, header); write(f, _mgs_varint_bytes(length(pb))); write(f, pb); write(f, body)
	end
	return
end
function _read_cand_rows(io::IO, n::Int)
	pb = read(io, _mgs_read_varint(io))
	params = Compression.cand_params_from_bytes(pb)
	return Compression.cand_decode(read(io), n; params=params)
end
"""
    write_cand_mgs3_graph(g, filename; symmetric=:off, reciprocal=:off, transpose=:off, params=CandParams())

Public candidate-coder writer; the storage modes behave as for [`write_cg_mgs3_graph`](@ref).
"""
function write_cand_mgs3_graph(g::AbstractGraph, filename::AbstractString; symmetric=:off, reciprocal=:off,
		transpose=:off, coding_scheme::Symbol=:children, params::Compression.CandParams=Compression.CandParams())
	_write_with_modes((h, f) -> _write_cand_mgs3_graph_core(h, f; coding_scheme=coding_scheme, params=params),
		g, filename, coding_scheme, symmetric, reciprocal, Pair{UInt8,Vector{UInt8}}[]; transpose=transpose)
end
@doc (@doc _write_bg_mgs3_graph_core) write_bg_mgs3_graph
@doc (@doc _write_cs_mgs3_graph_core) write_cs_mgs3_graph
@doc (@doc _write_cg_mgs3_graph_core) write_cg_mgs3_graph

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

CG-only options `compact_lists`, `fused_lr`, `reuse_identity`,
`sorted_fastpath`, and `reuse_scratch` all default to `false`. Enable the
recorded set together for compact local IDs, fused LR decoding, identity-vector
transfer, sorted finalization, and scratch reuse without changing the file format.
`parallel_decode=false` also defaults off; enabling it requires independent
indexed context-range CG clusters and positive `decode_workers`. K=1 remains
serial. CG options on non-CG files are rejected rather than silently ignored.

`parallel_decode=true` is opt-in for indexed context-range CG files only. It
decodes independent clusters with up to `decode_workers` default-pool threads;
graph construction remains serial. K=1 falls back to serial stream decoding.
"""
function load_compressed_mgs3_graph(filename::AbstractString;
		parallel_decode::Bool=false, decode_workers::Int=Threads.nthreads(:default),
		compact_lists::Bool=false, fused_lr::Bool=false, reuse_identity::Bool=false,
		sorted_fastpath::Bool=false, reuse_scratch::Bool=false)
	f = open(filename, "r")

	# read header (12 bytes, plus the optional extension block)
	header_bytes, ext = _read_mgs_header(f)

	# Verify MGS signature
	if header_bytes[1:3] != [0x4d, 0x47, 0x53]  # 'MGS'
		error("Invalid MGS file signature")
	end

	# Extract and decode flags
	flag_byte1 = header_bytes[6]
	flag_byte2 = header_bytes[7]
	graph_type, encoding, compression, byte2 = decode_header_flags(flag_byte1, flag_byte2)
	if (compact_lists || fused_lr || reuse_identity || sorted_fastpath || reuse_scratch) &&
			!(byte2 == ALG_CG || PARAM_CG_BASE <= byte2 <= PARAM_CG_MAX)
		close(f)
		throw(ArgumentError("CG decoder options require a CG file"))
	end
	if parallel_decode && (!(byte2 == ALG_CG || PARAM_CG_BASE <= byte2 <= PARAM_CG_MAX) ||
			encoding != :index || compression != :context_range || decode_workers <= 0)
		close(f)
		throw(ArgumentError("parallel_decode requires indexed context-range CG and positive decode_workers"))
	end

	if (_mgs_ext_symmetric(ext) || haskey(ext, MGS_EXT_RECIPROCAL) || haskey(ext, MGS_EXT_TRANSPOSED)) && encoding != :children
		close(f); throw(ArgumentError("symmetric/reciprocal MGS files must use coding_scheme=:children"))
	end
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
		load_cg_mgs3_graph(f, graph_type, gs; params=_apply_cg_ext(_cg_default_params(), ext), coding_scheme=encoding,
			integer_encoding=compression, parallel_decode=parallel_decode, decode_workers=decode_workers,
			compact_lists=compact_lists, fused_lr=fused_lr, reuse_identity=reuse_identity,
			sorted_fastpath=sorted_fastpath, reuse_scratch=reuse_scratch)
	elseif byte2 == ALG_CAND
		(encoding == :children && compression == :context_range && header_bytes[5] == 0x02) ||
			(close(f); throw(ArgumentError("malformed candidate-coder header")))
		rows = _read_cand_rows(f, Int(gs))
		T = infer_uint_custom_type(UInt8(ndigits(Int(gs); base=2)))
		h = SimpleDiGraph{T}(Int(gs))
		for u in 1:Int(gs), v in rows[u]; add_edge!(h, T(u), T(v)); end
		h
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
		# :context_range files store fibonacci-coded structural integers; the range
		# coder handles residuals separately (signalled by compression==:context_range).
		p = _apply_cg_ext(decode_cg_params(byte2; varint=(compression == :context_range ? :fibonacci : compression)), ext)
		load_cg_mgs3_graph(f, graph_type, gs; params=p, coding_scheme=encoding, integer_encoding=compression,
			parallel_decode=parallel_decode, decode_workers=decode_workers,
			compact_lists=compact_lists, fused_lr=fused_lr, reuse_identity=reuse_identity,
			sorted_fastpath=sorted_fastpath, reuse_scratch=reuse_scratch)
	else
		error("Unknown byte2: 0x$(string(byte2, base=16, pad=2))")
	end

	close(f)
	haskey(ext, MGS_EXT_RECIPROCAL) && _reciprocal_restore_graph!(g, ext[MGS_EXT_RECIPROCAL])
	_mgs_ext_symmetric(ext) && _mirror_upper_graph!(g)
	haskey(ext, MGS_EXT_TRANSPOSED) && (g = _transpose_graph(g))
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

	# :context_range (v3.3): [8-byte lens ×5: struct, refdist, copy, cmd, flag]
	# [struct][refdist][copy][cmd][flag][resid] — children and index/RA alike.
	ctx_range = integer_encoding == :context_range
	resid_bytes = UInt8[]; refdist_bytes = UInt8[]; copy_bytes = UInt8[]
	cmd_bytes = nothing; flag_bytes = nothing
	if ctx_range
		rd8() = (b = read(io, 8); L = UInt64(0); for i in 1:8; L |= UInt64(b[i]) << (8*(i-1)); end; L)
		slen = rd8(); rdlen = rd8(); cplen = rd8()
		cmdlen = Int(rd8()); flaglen = Int(rd8())
		struct_bytes = read(io, Int(slen))
		refdist_bytes = read(io, Int(rdlen))
		copy_bytes = read(io, Int(cplen))
		cmd_bytes = read(io, cmdlen)
		flag_bytes = read(io, flaglen)
		resid_bytes = read(io)                      # remaining bytes = residual blob
		reader = BitReader(struct_bytes)
	else
		reader = BitReader(io)
	end

	@info("reading greedy-compressed graph data")
	neighbor_lists = read_greedy_graph_data(reader, V(gs), coding_scheme, V;
		integer_encoding=integer_encoding, copy_blocks=copy_blocks,
		adaptive_copy=adaptive_copy, ref_window_size=ref_window_size,
		fixwidth_ref=fixwidth_ref, stop_deltas=stop_deltas,
		adaptive_deltas=adaptive_deltas,
		split_residual=split_residual, compact_copy=compact_copy,
		tight_intervals=tight_intervals,
		lr_split=lr_split, multi_ref=multi_ref,
		adaptive_header=adaptive_header,
		ctx_range=ctx_range, resid_bytes=resid_bytes,
		refdist_bytes=refdist_bytes, copy_bytes=copy_bytes,
		cmd_bytes=cmd_bytes, flag_bytes=flag_bytes)

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

	# :context_range (v3.3): [8-byte lens ×5: struct, refdist, copy, cmd, flag]
	# [struct][refdist][copy][cmd][flag][resid] — children and index/RA alike.
	ctx_range = integer_encoding == :context_range
	resid_bytes = UInt8[]; refdist_bytes = UInt8[]; copy_bytes = UInt8[]
	cmd_bytes = nothing; flag_bytes = nothing
	if ctx_range
		rd8() = (b = read(io, 8); L = UInt64(0); for i in 1:8; L |= UInt64(b[i]) << (8*(i-1)); end; L)
		slen = rd8(); rdlen = rd8(); cplen = rd8()
		cmdlen = Int(rd8()); flaglen = Int(rd8())
		struct_bytes = read(io, Int(slen))
		refdist_bytes = read(io, Int(rdlen))
		copy_bytes = read(io, Int(cplen))
		cmd_bytes = read(io, cmdlen)
		flag_bytes = read(io, flaglen)
		resid_bytes = read(io)
		reader = BitReader(struct_bytes)
	else
		reader = BitReader(io)
	end

	@info("reading CS-compressed graph data")
	neighbor_lists = read_cmdstream_graph_data(reader, V(gs), coding_scheme, V;
		integer_encoding=integer_encoding, compact_copy=compact_copy,
		tight_intervals=tight_intervals, ref_window_size=ref_window_size,
		lr_split=lr_split,
		ctx_range=ctx_range, resid_bytes=resid_bytes,
		refdist_bytes=refdist_bytes, copy_bytes=copy_bytes,
		cmd_bytes=cmd_bytes, flag_bytes=flag_bytes)

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
		params::CGParams=CGParams(), coding_scheme::Symbol=:children,
		integer_encoding::Symbol=:fibonacci, parallel_decode::Bool=false,
		decode_workers::Int=Threads.nthreads(:default), compact_lists::Bool=false,
		fused_lr::Bool=false, reuse_identity::Bool=false, sorted_fastpath::Bool=false,
		reuse_scratch::Bool=false)
	n_bits_v = convert(UInt8, ceil(log(2, gs)))
	V = infer_uint_custom_type(n_bits_v)

	directed = graph_type == :directed

	@info("generating graph (CG mode)")

	# :context_range — split [8-byte lens ×3][struct][refdist][copy][resid]
	ctx_range = integer_encoding == :context_range
	resid_bytes = UInt8[]; refdist_bytes = UInt8[]; copy_bytes = UInt8[]
	if ctx_range
		rd8() = (b = read(io, 8); L = UInt64(0); for i in 1:8; L |= UInt64(b[i]) << (8*(i-1)); end; L)
		slen = rd8(); rdlen = rd8(); cplen = rd8()
		struct_bytes = read(io, Int(slen))
		refdist_bytes = read(io, Int(rdlen))
		copy_bytes = read(io, Int(cplen))
		resid_bytes = read(io)
		reader = BitReader(struct_bytes)
	else
		reader = BitReader(io)
	end

	chunk_offsets = nothing
	_cg_offsets = nothing
	if coding_scheme == :index
		# Read cluster offset table: 6-bit entry_width + 32-bit K + (2K+1) entries
		entry_width = Int(read_value(reader, 6, UInt64))
		K = Int(read_value(reader, 32, UInt32))
		# Read offset entries: K intra + 1 inter start + K inter per-source
		_cg_offsets = [Int(read_value(reader, entry_width, UInt64)) for _ in 1:(2 * K + 1)]
		@info "CG index mode: K=$K, entry_width=$entry_width"
		if ctx_range
			# CG-RA: per-cluster chunk byte-offset tables for the three range streams
			coff_width = Int(read_value(reader, 6, UInt64))
			rc_off = [Int(read_value(reader, coff_width, UInt64)) for _ in 1:(K + 1)]
			rd_off = [Int(read_value(reader, coff_width, UInt64)) for _ in 1:(K + 1)]
			cp_off = [Int(read_value(reader, coff_width, UInt64)) for _ in 1:(K + 1)]
			chunk_offsets = (rc_off, rd_off, cp_off)
		end
		# Fall through to normal sequential decode
	end

	@info("reading CG-compressed graph data")
	neighbor_lists = decode_level(reader, params; T=V, directed=directed, coding_scheme=coding_scheme,
		ctx_range=ctx_range, resid_bytes=resid_bytes, refdist_bytes=refdist_bytes, copy_bytes=copy_bytes,
		chunk_offsets=chunk_offsets, cg_offsets=_cg_offsets, data_start_bit=reader.bit_count,
		parallel_clusters=parallel_decode, decode_workers=decode_workers,
		compact_lists=compact_lists, fused_lr=fused_lr, reuse_identity=reuse_identity,
		sorted_fastpath=sorted_fastpath, reuse_scratch=reuse_scratch)

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

include("mgs_cache.jl")

const _INDEXED_MGS3_DOC = """
    load_indexed_mgs3_graph(filename; decode_prefix=true, max_live_chunks=16,
                          cache_budget_bytes=nothing)

Load an MGS v3.2 compressed graph with sampled index and return an object
supporting indexed access where available, or a full-decode fallback otherwise.

Returns a named tuple `(n, m, neighbors, algorithm, ...)` where:
- `n`: number of vertices
- `m`: number of edges (-1 in the true-seek path, where it is not computed)
- `neighbors(v)`: function that returns sorted neighbor list of vertex v (1-indexed)
- `algorithm`: symbol (:bg, :cs, or :cg)

Access paths:
- **BG/CS + :context_range + sampled index** → indexed seek: `neighbors(v)`
  decodes through v, not the rest of its chunk. Earlier cross-chunk references
  are recursively resolved; only files encoded with chunk-local references
  guarantee one-chunk access. Work also depends on decoded degrees, not just k.
  Prefix adjacency and live coder state are memoised, allowing later vertices in
  the same chunk to resume without replay. Coder models are released when a chunk
  is complete. At most `max_live_chunks` (default 16) suspended workspaces are
  retained, in least-recently-used order. An evicted prefix can be replayed if a
  later query extends it; its already memoised answers remain available. This
  bound excludes active recursive calls and the existing adjacency memo.
  `decode_prefix=false` restores whole-chunk materialisation for
  comparisons or workloads that prefer it. This keyword does not change files.
  `stats` includes distinct `chunks_decoded`, `decode_calls`, `vertices_decoded`,
  `arcs_decoded`, `total_chunks`, `samp_k`, and `active_workspaces()`.
  `reset_fn()` clears all memoised adjacency, continuation positions, and counters.
  At most one additional idle workspace retains reusable allocation buffers, not
  answers. Every reuse resets all probability tables and registers inside the
  decode call. `stats.scratch_bytes()` reports this idle storage excluding shared
  compressed input; it does not include live continuation workspaces or adjacency.
  For an opt-in byte-budgeted whole-chunk cache, set `decode_prefix=false` and
  `cache_budget_bytes` to a nonnegative byte count. Decoded chunks are packed into
  offsets/targets arrays and evicted by LRU; oversized chunks bypass the cache.
  The budget includes nonempty cache entries and their hash-table storage, but
  excludes the empty manager, compressed input, active query scratch and output.
  Cache stats distinguish top-level hits/misses, retained/peak bytes, evictions,
  oversized bypasses and actual manager heap size. This mode retains no live
  decoder models between queries. A zero budget disables persistent caching.
- **CG + :context_range + index** → indexed whole-cluster seek.
- **Other formats** → the graph is fully decoded once and cached, giving O(1)
  access to any vertex thereafter (no per-vertex seek).

# Example
```julia
idx = load_indexed_mgs3_graph("graph.mgz")
println("Graph: \$(idx.n) vertices, \$(idx.m) edges")
nbrs = idx.neighbors(42)  # get neighbors of vertex 42
```
"""
# Decode only the missing prefix. Keep live entropy state until the chunk is
# complete. One idle model allocation may be reused after a full reset; completed
# adjacency stays memoised separately and is never exposed as mutable scratch.
function _materialize_indexed_prefix!(decode!, memo, states, pool, c, target, samp_k, vs,
		decode_prefix, max_live_chunks, chunks_decoded, decode_calls, vertices_decoded, arcs_decoded)
	if haskey(memo, c) && haskey(memo[c], target)
		return memo[c]
	end
	cs = (c - 1) * samp_k + 1
	ce = min(c * samp_k, vs)
	stop = decode_prefix ? Int(target) : ce
	T = typeof(target)
	state = get!(states, c) do
		isempty(pool) ? valtype(typeof(states))() : pop!(pool)
	end
	start = state.next_vertex == 0 ? cs : state.next_vertex
	is_new = !haskey(memo, c)
	local nls
	try
		nls = decode!(start:stop, state)
	catch
		# A failed decode may have partially advanced the live models. Never reuse it.
		delete!(states, c); delete!(memo, c)
		rethrow()
	end
	memo[c] = nls
	chunks_decoded[] += is_new
	decode_calls[] += 1
	vertices_decoded[] += stop - start + 1
	for v in start:stop
		arcs_decoded[] += length(nls[T(v)])
	end
	delete!(states, c)
	if stop < ce
		states[c] = state # most recently resumed workspace goes to the end
		while length(states) > max_live_chunks
			# An evicted ancestor may still be active on the recursive call stack.
			# Remove its cache entry, but never recycle its mutable models here.
			popfirst!(states)
		end
	elseif isempty(pool)
		push!(pool, _recycle_indexed_state!(state))
	end
	return nls
end

function _reset_indexed_reader!(memo, states, pool)
	if isempty(pool) && !isempty(states)
		push!(pool, _recycle_indexed_state!(first(values(states))))
	end
	empty!(memo); empty!(states)
	return nothing
end
_indexed_scratch_bytes(pool) = Base.summarysize(pool; exclude=Union{Module,DataType,Task,Vector{UInt8}})

# ── True random-access reader: BG + :context_range + sampled index ────────────
# Returns (n, m, neighbors, algorithm, stats) where `neighbors(v)` seeks to the
# chunk containing v, decodes through v, and recursively materialises required
# prefixes in earlier chunks. `stats.chunks_decoded` counts distinct touched
# chunks, not continuation calls; vertices/arcs counters count actual decode work.
function _load_bg_ctxrange_random_access(filename::AbstractString, graph_type::Symbol,
		gs::UInt64, params; decode_prefix::Bool=true, max_live_chunks::Int=16, cache_budget_bytes=nothing)
	# Re-read the file body (header already validated by caller).
	local struct_bytes, refdist_bytes, copy_bytes, cmd_bytes, flag_bytes, resid_bytes
	open(filename, "r") do io
		_read_mgs_header(io)  # skip header (and extension)
		rd8() = (b = read(io, 8); L = UInt64(0); for i in 1:8; L |= UInt64(b[i]) << (8*(i-1)); end; L)
		slen = rd8(); rdlen = rd8(); cplen = rd8(); cmdlen = rd8(); flaglen = rd8()
		struct_bytes  = read(io, Int(slen))
		refdist_bytes = read(io, Int(rdlen))
		copy_bytes    = read(io, Int(cplen))
		cmd_bytes     = read(io, Int(cmdlen))
		flag_bytes    = read(io, Int(flaglen))
		resid_bytes   = read(io)
	end

	n_bits_v = convert(UInt8, ceil(log(2, gs)))
	T = infer_uint_custom_type(n_bits_v)
	vs = Int(gs)
	meta = parse_bg_ctxrange_index(struct_bytes, vs)
	samp_k = meta.samp_k
	rws = params.ref_window_size

	prototype = IndexedRangeState(T; dense=true)
	memo = Dict{Int,typeof(prototype.neighbor_lists)}()
	states = OrderedDict{Int,typeof(prototype)}()
	pool = [prototype] # at most one idle, answer-free workspace
	chunks_decoded = Ref(0)
	decode_calls = Ref(0); vertices_decoded = Ref(0); arcs_decoded = Ref(0)
	chunk_of(u::Int) = ((u - 1) ÷ samp_k) + 1

	# Forward declaration for mutual recursion (chunk → resolver → chunk).
	local get_neighbors_v
	function decode_prefix!(range, state)
		cs = ((first(range) - 1) ÷ samp_k) * samp_k + 1
		seed_window = T[T(u) for u in max(1, cs - rws):(cs - 1)]
		reader = BitReader(struct_bytes)
		nls = read_greedy_graph_data(reader, T(gs), :index, T;
			integer_encoding=:context_range,
			copy_blocks=params.copy_blocks, adaptive_copy=params.adaptive_copy,
			ref_window_size=params.ref_window_size, fixwidth_ref=params.fixwidth_ref,
			stop_deltas=params.stop_deltas, compact_copy=params.compact_copy,
			tight_intervals=params.tight_intervals, lr_split=params.lr_split,
			adaptive_header=params.adaptive_header,
			ctx_range=true, resid_bytes=resid_bytes, refdist_bytes=refdist_bytes,
			copy_bytes=copy_bytes, cmd_bytes=cmd_bytes, flag_bytes=flag_bytes,
			decode_range=range, seed_window=seed_window, ref_resolver=get_neighbors_v,
			preparsed_index=meta, workspace=state)
		return nls
	end
	if cache_budget_bytes !== nothing
		bounded=_budgeted_chunk_reader(decode_prefix!,T,vs,samp_k,:bg,cache_budget_bytes)
		get_neighbors_v=bounded.resolve
		return bounded.api
	end
	get_neighbors_v = function(u)
		return _materialize_indexed_prefix!(decode_prefix!, memo, states, pool, chunk_of(Int(u)), T(u),
			samp_k, vs, decode_prefix, max_live_chunks, chunks_decoded, decode_calls, vertices_decoded, arcs_decoded)[T(u)]
	end

	function _neighbors(v::Integer)
		1 <= v <= vs || error("Vertex $v out of range [1, $vs]")
		return sort(get_neighbors_v(T(v)))
	end

	return (n=vs, m=-1, neighbors=_neighbors, algorithm=:bg,
		stats=(chunks_decoded=chunks_decoded, total_chunks=meta.n_sampled, samp_k=samp_k,
			decode_calls=decode_calls, vertices_decoded=vertices_decoded, arcs_decoded=arcs_decoded,
			active_workspaces=() -> length(states), scratch_bytes=() -> _indexed_scratch_bytes(pool)),
		reset_fn=() -> (_reset_indexed_reader!(memo, states, pool); chunks_decoded[] = 0;
			decode_calls[] = 0; vertices_decoded[] = 0; arcs_decoded[] = 0))
end

# ── True random-access reader: CS + :context_range + sampled index ────────────
# Structurally identical to the BG reader (the CS struct-header layout matches, so
# parse_bg_ctxrange_index is reused); the only difference is the per-chunk decode
# calls read_cmdstream_graph_data with CS params.
function _load_cs_ctxrange_random_access(filename::AbstractString, graph_type::Symbol,
		gs::UInt64, params; decode_prefix::Bool=true, max_live_chunks::Int=16, cache_budget_bytes=nothing)
	local struct_bytes, refdist_bytes, copy_bytes, cmd_bytes, flag_bytes, resid_bytes
	open(filename, "r") do io
		_read_mgs_header(io)  # skip header (and extension)
		rd8() = (b = read(io, 8); L = UInt64(0); for i in 1:8; L |= UInt64(b[i]) << (8*(i-1)); end; L)
		slen = rd8(); rdlen = rd8(); cplen = rd8(); cmdlen = rd8(); flaglen = rd8()
		struct_bytes  = read(io, Int(slen))
		refdist_bytes = read(io, Int(rdlen))
		copy_bytes    = read(io, Int(cplen))
		cmd_bytes     = read(io, Int(cmdlen))
		flag_bytes    = read(io, Int(flaglen))
		resid_bytes   = read(io)
	end

	n_bits_v = convert(UInt8, ceil(log(2, gs)))
	T = infer_uint_custom_type(n_bits_v)
	vs = Int(gs)
	meta = parse_bg_ctxrange_index(struct_bytes, vs)
	samp_k = meta.samp_k
	rws = params.ref_window_size

	prototype = IndexedRangeState(T; dense=true)
	memo = Dict{Int,typeof(prototype.neighbor_lists)}()
	states = OrderedDict{Int,typeof(prototype)}()
	pool = [prototype]
	chunks_decoded = Ref(0)
	decode_calls = Ref(0); vertices_decoded = Ref(0); arcs_decoded = Ref(0)
	chunk_of(u::Int) = ((u - 1) ÷ samp_k) + 1

	local get_neighbors_v
	function decode_prefix!(range, state)
		cs = ((first(range) - 1) ÷ samp_k) * samp_k + 1
		seed_window = T[T(u) for u in max(1, cs - rws):(cs - 1)]
		reader = BitReader(struct_bytes)
		nls = read_cmdstream_graph_data(reader, T(gs), :index, T;
			integer_encoding=:context_range,
			compact_copy=params.compact_copy, tight_intervals=params.tight_intervals,
			ref_window_size=params.ref_window_size, lr_split=params.lr_split,
			ctx_range=true, resid_bytes=resid_bytes, refdist_bytes=refdist_bytes,
			copy_bytes=copy_bytes, cmd_bytes=cmd_bytes, flag_bytes=flag_bytes,
			decode_range=range, seed_window=seed_window, ref_resolver=get_neighbors_v,
			preparsed_index=meta, workspace=state)
		return nls
	end
	if cache_budget_bytes !== nothing
		bounded=_budgeted_chunk_reader(decode_prefix!,T,vs,samp_k,:cs,cache_budget_bytes)
		get_neighbors_v=bounded.resolve
		return bounded.api
	end
	get_neighbors_v = function(u)
		return _materialize_indexed_prefix!(decode_prefix!, memo, states, pool, chunk_of(Int(u)), T(u),
			samp_k, vs, decode_prefix, max_live_chunks, chunks_decoded, decode_calls, vertices_decoded, arcs_decoded)[T(u)]
	end

	function _neighbors(v::Integer)
		1 <= v <= vs || error("Vertex $v out of range [1, $vs]")
		return sort(get_neighbors_v(T(v)))
	end

	return (n=vs, m=-1, neighbors=_neighbors, algorithm=:cs,
		stats=(chunks_decoded=chunks_decoded, total_chunks=meta.n_sampled, samp_k=samp_k,
			decode_calls=decode_calls, vertices_decoded=vertices_decoded, arcs_decoded=arcs_decoded,
			active_workspaces=() -> length(states), scratch_bytes=() -> _indexed_scratch_bytes(pool)),
		reset_fn=() -> (_reset_indexed_reader!(memo, states, pool); chunks_decoded[] = 0;
			decode_calls[] = 0; vertices_decoded[] = 0; arcs_decoded[] = 0))
end

# ── True random-access reader: CG + :context_range + index ────────────────────
# Cluster = chunk. The struct stream carries [cluster offset table][chunk-offset
# tables ×3][data]; neighbors(v) decodes only v's cluster (its intra section +
# its per-source inter section) on fresh per-cluster range decoders. Intra
# references never cross clusters, so no cross-chunk resolver is needed.
function _load_cg_ctxrange_random_access(filename::AbstractString, graph_type::Symbol,
		gs::UInt64, params::CGParams; compact_lists::Bool=false, fused_lr::Bool=false,
		reuse_identity::Bool=false, sorted_fastpath::Bool=false, reuse_scratch::Bool=false)
	local struct_bytes, refdist_bytes, copy_bytes, resid_bytes
	open(filename, "r") do io
		_read_mgs_header(io)  # skip header (and extension)
		rd8() = (b = read(io, 8); L = UInt64(0); for i in 1:8; L |= UInt64(b[i]) << (8*(i-1)); end; L)
		slen = rd8(); rdlen = rd8(); cplen = rd8()
		struct_bytes  = read(io, Int(slen))
		refdist_bytes = read(io, Int(rdlen))
		copy_bytes    = read(io, Int(cplen))
		resid_bytes   = read(io)
	end

	n_bits_v = convert(UInt8, ceil(log(2, gs)))
	T = infer_uint_custom_type(n_bits_v)
	vs = Int(gs)
	directed = graph_type == :directed

	# Parse offset tables + membership once.
	r = BitReader(struct_bytes)
	entry_width = Int(read_value(r, 6, UInt64))
	K = Int(read_value(r, 32, UInt32))
	cg_offsets = [Int(read_value(r, entry_width, UInt64)) for _ in 1:(2 * K + 1)]
	coff_width = Int(read_value(r, 6, UInt64))
	rc_off = [Int(read_value(r, coff_width, UInt64)) for _ in 1:(K + 1)]
	rd_off = [Int(read_value(r, coff_width, UInt64)) for _ in 1:(K + 1)]
	cp_off = [Int(read_value(r, coff_width, UInt64)) for _ in 1:(K + 1)]
	data_start_bit = r.bit_count
	clusters = _read_membership(r, params, T)

	# vertex → cluster index lookup
	cluster_of = zeros(Int32, vs)
	for (ci, C) in enumerate(clusters), u in C
		cluster_of[Int(u)] = Int32(ci)
	end

	memo = Dict{Int,Dict{T,Vector{T}}}()
	clusters_decoded = Ref(0)
	function materialize_cluster(ci::Int)
		haskey(memo, ci) && return memo[ci]
		rr = BitReader(struct_bytes)
		nls = decode_level(rr, params; T=T, directed=directed, coding_scheme=:index,
			ctx_range=true, resid_bytes=resid_bytes, refdist_bytes=refdist_bytes,
			copy_bytes=copy_bytes, cg_offsets=cg_offsets, data_start_bit=data_start_bit,
			chunk_offsets=(rc_off, rd_off, cp_off), only_cluster=ci,
			preparsed_clusters=clusters, local_type=(compact_lists ? T : Int), fused_lr=fused_lr,
			reuse_identity=reuse_identity, sorted_fastpath=sorted_fastpath, reuse_scratch=reuse_scratch)
		memo[ci] = nls
		clusters_decoded[] += 1
		return nls
	end

	function _neighbors(v::Integer)
		1 <= v <= vs || error("Vertex $v out of range [1, $vs]")
		ci = Int(cluster_of[Int(v)])
		ci == 0 && return T[]
		nls = materialize_cluster(ci)
		# Returned storage stays caller-owned; the decoded cluster remains memoized.
		return sort(get(nls, T(v), T[]))
	end

	return (n=vs, m=-1, neighbors=_neighbors, algorithm=:cg,
		stats=(chunks_decoded=clusters_decoded, total_chunks=K, samp_k=-1),
		reset_fn=() -> (empty!(memo); clusters_decoded[] = 0))
end

function load_indexed_mgs3_graph(filename::AbstractString; decode_prefix::Bool=true, max_live_chunks::Int=16,
		cache_budget_bytes::Union{Nothing,Int}=nothing, compact_lists::Bool=false,
		fused_lr::Bool=false, reuse_identity::Bool=false, sorted_fastpath::Bool=false,
		reuse_scratch::Bool=false)
	max_live_chunks >= 1 || throw(ArgumentError("max_live_chunks must be positive"))
	if cache_budget_bytes !== nothing
		cache_budget_bytes >= 0 || throw(ArgumentError("cache_budget_bytes must be nonnegative"))
		!decode_prefix || throw(ArgumentError("byte-budgeted cache requires decode_prefix=false"))
	end
	# Indexed range streams need no upfront full decode. Cross-chunk references
	# in older files can still amplify the seek beyond the target chunk.
	open(filename, "r") do f
		header, ext = _read_mgs_header(f)
		(_mgs_ext_symmetric(ext) || haskey(ext, MGS_EXT_RECIPROCAL) || haskey(ext, MGS_EXT_TRANSPOSED)) &&
			throw(ArgumentError("random access is not supported for symmetric, reciprocal or transposed MGS files"))
		if header[1:3] == [0x4d, 0x47, 0x53]
			flag_byte1 = header[6]; flag_byte2 = header[7]
			graph_type, encoding, compression, byte2 = decode_header_flags(flag_byte1, flag_byte2)
			if (compact_lists || fused_lr || reuse_identity || sorted_fastpath || reuse_scratch) &&
					!(byte2 == ALG_CG || PARAM_CG_BASE <= byte2 <= PARAM_CG_MAX)
				throw(ArgumentError("CG decoder options require a CG file"))
			end
			gs_bytes = header[8:12]
			gs = UInt64(gs_bytes[1]) | (UInt64(gs_bytes[2]) << 8) | (UInt64(gs_bytes[3]) << 16) |
				(UInt64(gs_bytes[4]) << 24) | (UInt64(gs_bytes[5]) << 32)
			is_bg = byte2 == ALG_BG || (PARAM_BG_BASE <= byte2 <= PARAM_BG_MAX)
			if is_bg && encoding == :index && compression == :context_range
				p = byte2 == ALG_BG ? _bg_default_params() : decode_bg_params(byte2)
				return _load_bg_ctxrange_random_access(filename, graph_type, gs, p;
					decode_prefix=decode_prefix, max_live_chunks=max_live_chunks,cache_budget_bytes=cache_budget_bytes)
			end
			is_cs = byte2 == ALG_CS || (PARAM_CS_BASE <= byte2 <= PARAM_CS_MAX)
			if is_cs && encoding == :index && compression == :context_range
				p = byte2 == ALG_CS ? _cs_default_params() : decode_cs_params(byte2)
				return _load_cs_ctxrange_random_access(filename, graph_type, gs, p;
					decode_prefix=decode_prefix, max_live_chunks=max_live_chunks,cache_budget_bytes=cache_budget_bytes)
			end
			is_cg = byte2 == ALG_CG || (PARAM_CG_BASE <= byte2 <= PARAM_CG_MAX)
			if is_cg && encoding == :index && compression == :context_range
				cache_budget_bytes === nothing || throw(ArgumentError("byte-budgeted cache currently supports BG/CS indexed context-range files only"))
				p = _apply_cg_ext(byte2 == ALG_CG ? CGParams() : decode_cg_params(byte2), ext)
				return _load_cg_ctxrange_random_access(filename, graph_type, gs, p;
					compact_lists=compact_lists, fused_lr=fused_lr, reuse_identity=reuse_identity,
					sorted_fastpath=sorted_fastpath, reuse_scratch=reuse_scratch)
			end
		end
		cache_budget_bytes === nothing || throw(ArgumentError("byte-budgeted cache currently supports BG/CS indexed context-range files only"))
		return nothing
	end |> ra -> ra === nothing ? _load_indexed_mgs3_graph_fulldecode(filename;
		compact_lists=compact_lists, fused_lr=fused_lr, reuse_identity=reuse_identity,
		sorted_fastpath=sorted_fastpath, reuse_scratch=reuse_scratch) : ra
end

@doc _INDEXED_MGS3_DOC load_indexed_mgs3_graph

function _load_indexed_mgs3_graph_fulldecode(filename::AbstractString; kw...)
	# Load full graph (uses sequential decode internally)
	g = load_compressed_mgs3_graph(filename; kw...)
	n = Int(nv(g))
	m = Int(ne(g))

	# Pre-build sorted neighbor lists for O(1) access
	T = eltype(g)
	adj = Vector{Vector{T}}(undef, n)
	for v in vertices(g)
		adj[Int(v)] = sort(collect(outneighbors(g, v)))
	end

	# Determine algorithm from header
	alg = :unknown
	open(filename, "r") do f
		header, _ = _read_mgs_header(f)
		if header[1:3] == [0x4d, 0x47, 0x53]
			byte2 = header[7]
			if byte2 < 0x50
				alg = :bg
			elseif byte2 < 0x70
				alg = :cs
			else
				alg = :cg
			end
		end
	end

	function _neighbors(v::Integer)
		1 <= v <= n || error("Vertex $v out of range [1, $n]")
		return adj[Int(v)]
	end

	return (n=n, m=m, neighbors=_neighbors, algorithm=alg)
end

include("mgs_adjacency.jl")

end # module MGS
