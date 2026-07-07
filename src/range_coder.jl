# Streaming context-adaptive range coder for the `:context_range` integer backend.
#
# LZMA-style 32-bit carry-handling range coder + hybrid-token model:
#   value v -> token = bitlength(v), coded with a per-context adaptive frequency
#   table; the (token-1) mantissa bits below the MSB are written raw (bypass).
# Context = (position-in-list bucket) × (previous-token bucket), with the
# previous-token bucket reset to a sentinel at the first residual of each vertex.
# Both features are reconstructible by the decoder (position = pull index within
# the vertex's residual region; previous token = last decoded token), so the
# coder is losslessly decodable in a single pass in lockstep with the structural
# decode. This is a faithful streaming port of the verified prototype in
# research/graph_compression/v3/prototype/context_coder.jl (ctx_encode_typed).
#
# Included into `module Compression` (no module of its own).

const _RC_TOP   = UInt32(1) << 24

# ---- Hybrid-uint (JPEG-XL / Zuckerli) token parameters ----
# value v -> (token, nbits_raw, raw): small values (v < 2^K) are direct tokens;
# larger values fold the exponent + M high mantissa bits + L low mantissa bits
# into the entropy-coded token, and only the middle `nbits_raw` bits stay raw.
# Require K >= M + L. Alphabet size NSYM = 2^K + (64-K)*2^(M+L) (covers UInt64).
const _RC_K = 5                      # split_exponent
const _RC_M = 0                      # msb_in_token
const _RC_L = 1                      # lsb_in_token
@assert _RC_K >= _RC_M + _RC_L "hybrid-uint requires K >= M + L"
const _RC_ML   = _RC_M + _RC_L
const _RC_SPLIT = 1 << _RC_K
const _RC_NSYM  = _RC_SPLIT + (64 - _RC_K) * (1 << _RC_ML)   # hybrid token alphabet

const _RC_NPOS  = 4                  # position buckets: 0,1,2,3+
const _RC_NPREV = 25                 # prev-magnitude buckets 0..23 + reset sentinel (24)
const _RC_NCTX  = _RC_NPOS * _RC_NPREV
const _RC_RESET = _RC_NPREV - 1      # prev sentinel used at a region boundary
const _RC_INC   = UInt32(32)

@inline _rc_bitlen(v::UInt64) = v == 0 ? 0 : (64 - leading_zeros(v))
@inline _rc_posb(p::Int)      = p < _RC_NPOS - 1 ? p : _RC_NPOS - 1
@inline _rc_ctx(posb::Int, prev::Int) = posb * _RC_NPREV + prev + 1   # 1-based

# Hybrid-uint encode: v -> (token::Int in 0..NSYM-1, nbits_raw::Int, raw::UInt64).
@inline function _rc_hybrid_encode(v::UInt64)
    if v < UInt64(_RC_SPLIT)
        return (Int(v), 0, UInt64(0))
    end
    logv = 63 - leading_zeros(v)                       # floor_log2(v)
    exp  = logv - _RC_K
    nbits_raw = logv - _RC_M - _RC_L
    lsb = _RC_L == 0 ? UInt64(0) : (v & ((UInt64(1) << _RC_L) - 1))
    raw = nbits_raw == 0 ? UInt64(0) : ((v >> _RC_L) & ((UInt64(1) << nbits_raw) - 1))
    msb = _RC_M == 0 ? UInt64(0) : ((v >> (nbits_raw + _RC_L)) & ((UInt64(1) << _RC_M) - 1))
    token = _RC_SPLIT + (exp << _RC_ML) + (Int(msb) << _RC_L) + Int(lsb)
    return (token, nbits_raw, raw)
end

# ---------------- Encoder ----------------
mutable struct CtxRangeEncoder
    low::UInt64
    range::UInt32
    cache::UInt8
    cachesize::Int
    out::Vector{UInt8}
    freq::Matrix{UInt32}     # _RC_NCTX × _RC_NSYM
    tot::Vector{UInt32}      # _RC_NCTX
    pos::Int                 # running position within current region
    prev::Int                # running previous-token bucket
end
CtxRangeEncoder() = CtxRangeEncoder(0, typemax(UInt32), 0, 1, UInt8[],
    ones(UInt32, _RC_NCTX, _RC_NSYM), fill(UInt32(_RC_NSYM), _RC_NCTX),
    0, _RC_RESET)

@inline function rc_reset_region!(e::CtxRangeEncoder)
    e.pos = 0; e.prev = _RC_RESET
    return nothing
end

@inline function _rc_shiftlow!(e::CtxRangeEncoder)
    if e.low < 0xFF000000 || e.low > 0xFFFFFFFF
        c = e.cache
        while true
            push!(e.out, UInt8((c + (e.low >> 32)) & 0xFF))
            c = 0xFF
            e.cachesize -= 1
            e.cachesize == 0 && break
        end
        e.cache = UInt8((e.low >> 24) & 0xFF)
    end
    e.cachesize += 1
    e.low = (e.low << 8) & 0xFFFFFFFF
    return nothing
end

@inline function _rc_enc_freq!(e::CtxRangeEncoder, cum::UInt32, f::UInt32, tot::UInt32)
    r = e.range ÷ tot
    e.low += UInt64(r) * UInt64(cum)
    e.range = r * f
    while e.range < _RC_TOP
        e.range <<= 8
        _rc_shiftlow!(e)
    end
    return nothing
end

@inline function _rc_enc_bits!(e::CtxRangeEncoder, val::UInt64, nbits::Int)
    for i in (nbits-1):-1:0
        e.range >>= 1
        b = (val >> i) & 1
        e.low += UInt64(e.range) * b
        while e.range < _RC_TOP
            e.range <<= 8
            _rc_shiftlow!(e)
        end
    end
    return nothing
end

@inline function _rc_update!(freq::Matrix{UInt32}, tot::Vector{UInt32}, ctx::Int, sym::Int)
    freq[ctx, sym] += _RC_INC
    tot[ctx] += _RC_INC
    if tot[ctx] >= _RC_TOP ÷ 2
        t = UInt32(0)
        for s in 1:_RC_NSYM
            freq[ctx, s] = (freq[ctx, s] >> 1) | 1; t += freq[ctx, s]
        end
        tot[ctx] = t
    end
    return nothing
end

# Encode one value, advancing the running context state.
function rc_encode_value!(e::CtxRangeEncoder, v::UInt64)
    token, nbits_raw, raw = _rc_hybrid_encode(v)
    sym = token + 1
    ctx = _rc_ctx(_rc_posb(e.pos), e.prev)
    cum = UInt32(0)
    for s in 1:sym-1; cum += e.freq[ctx, s]; end
    _rc_enc_freq!(e, cum, e.freq[ctx, sym], e.tot[ctx])
    _rc_update!(e.freq, e.tot, ctx, sym)
    if nbits_raw > 0
        _rc_enc_bits!(e, raw, nbits_raw)
    end
    e.pos += 1
    # Context state stays coarse: previous-MAGNITUDE bucket (bitlength of the raw
    # value), NOT the hybrid token id — keeps context sparsity low.
    e.prev = min(_rc_bitlen(v), _RC_NPREV - 2)   # 0..23; 24 reserved for reset
    return nothing
end

function rc_finish!(e::CtxRangeEncoder)
    for _ in 1:5; _rc_shiftlow!(e); end
    return e.out
end

# Finalize the current chunk (flush to a byte-aligned blob) and reset ALL state —
# range coder registers + adaptive freq tables + position/prev context — so the
# encoder can start a fresh, independently-decodable chunk. Used by the random-
# access (chunked) :context_range path. Returns the finalized chunk bytes.
function rc_finish_and_reset!(e::CtxRangeEncoder)
    bytes = rc_finish!(e)          # returns e.out (its own vector)
    e.low = UInt64(0)
    e.range = typemax(UInt32)
    e.cache = UInt8(0)
    e.cachesize = 1
    e.out = UInt8[]               # new vector; returned bytes are not aliased
    fill!(e.freq, UInt32(1))
    fill!(e.tot, UInt32(_RC_NSYM))
    e.pos = 0
    e.prev = _RC_RESET
    return bytes
end

# ---------------- Decoder ----------------
mutable struct CtxRangeDecoder
    code::UInt32
    range::UInt32
    inp::Vector{UInt8}
    ip::Int
    freq::Matrix{UInt32}
    tot::Vector{UInt32}
    pos::Int
    prev::Int
end
function CtxRangeDecoder(bytes::Vector{UInt8})
    d = CtxRangeDecoder(UInt32(0), typemax(UInt32), bytes, 1,
        ones(UInt32, _RC_NCTX, _RC_NSYM), fill(UInt32(_RC_NSYM), _RC_NCTX),
        0, _RC_RESET)
    d.ip += 1                              # skip first (cache) byte
    for _ in 1:4
        d.code = (d.code << 8) | UInt32(get(d.inp, d.ip, 0x00)); d.ip += 1
    end
    return d
end

@inline function rc_reset_region!(d::CtxRangeDecoder)
    d.pos = 0; d.prev = _RC_RESET
    return nothing
end

@inline function _rc_dec_update!(d::CtxRangeDecoder, cum::UInt32, f::UInt32)
    d.code -= cum * d.range
    d.range *= f
    while d.range < _RC_TOP
        d.code = (d.code << 8) | UInt32(get(d.inp, d.ip, 0x00)); d.ip += 1
        d.range <<= 8
    end
    return nothing
end

@inline function _rc_dec_bits!(d::CtxRangeDecoder, nbits::Int)
    v = UInt64(0)
    for _ in 1:nbits
        d.range >>= 1
        t = (d.code - d.range) >> 31          # 1 if code < range else 0
        b = UInt64(1 - t)
        d.code -= d.range & (t - UInt32(1))    # subtract range iff bit==1
        v = (v << 1) | b
        while d.range < _RC_TOP
            d.code = (d.code << 8) | UInt32(get(d.inp, d.ip, 0x00)); d.ip += 1
            d.range <<= 8
        end
    end
    return v
end

# Hybrid-uint decode: read raw tail bits (if any) from `d` and reconstruct v.
@inline function _rc_hybrid_decode(d::CtxRangeDecoder, token::Int)
    if token < _RC_SPLIT
        return UInt64(token)
    end
    t = token - _RC_SPLIT
    exp    = t >> _RC_ML
    packed = _RC_ML == 0 ? 0 : (t & ((1 << _RC_ML) - 1))
    msb = _RC_M == 0 ? 0 : (packed >> _RC_L)
    lsb = _RC_L == 0 ? 0 : (packed & ((1 << _RC_L) - 1))
    logv = _RC_K + exp
    nbits_raw = logv - _RC_M - _RC_L
    raw = nbits_raw == 0 ? UInt64(0) : _rc_dec_bits!(d, nbits_raw)
    return (UInt64(1) << logv) | (UInt64(msb) << (nbits_raw + _RC_L)) |
           (raw << _RC_L) | UInt64(lsb)
end

# Decode one value, advancing the running context state (mirror of encode).
function rc_decode_value!(d::CtxRangeDecoder)
    ctx = _rc_ctx(_rc_posb(d.pos), d.prev)
    d.range ÷= d.tot[ctx]
    f = min(d.code ÷ d.range, d.tot[ctx] - 1)
    cum = UInt32(0); sym = 1
    while cum + d.freq[ctx, sym] <= f; cum += d.freq[ctx, sym]; sym += 1; end
    _rc_dec_update!(d, cum, d.freq[ctx, sym])
    _rc_update!(d.freq, d.tot, ctx, sym)
    token = sym - 1
    v = _rc_hybrid_decode(d, token)
    d.pos += 1
    d.prev = min(_rc_bitlen(v), _RC_NPREV - 2)
    return v
end

# ================= Binary range coder (copy bitmap; context = previous bit) =====
const _BRC_MAXTOT = UInt32(1) << 16
const _BRC_INC    = UInt32(24)

mutable struct BinRangeEncoder
    low::UInt64
    range::UInt32
    cache::UInt8
    cachesize::Int
    out::Vector{UInt8}
    c0::Vector{UInt32}       # count of 0 per context (2 ctx: prev bit 0/1)
    c1::Vector{UInt32}
    prev::Int
end
BinRangeEncoder() = BinRangeEncoder(0, typemax(UInt32), 0, 1, UInt8[],
    ones(UInt32, 2), ones(UInt32, 2), 1)

@inline function _brc_shiftlow!(e::BinRangeEncoder)
    if e.low < 0xFF000000 || e.low > 0xFFFFFFFF
        c = e.cache
        while true
            push!(e.out, UInt8((c + (e.low >> 32)) & 0xFF))
            c = 0xFF
            e.cachesize -= 1
            e.cachesize == 0 && break
        end
        e.cache = UInt8((e.low >> 24) & 0xFF)
    end
    e.cachesize += 1
    e.low = (e.low << 8) & 0xFFFFFFFF
    return nothing
end

@inline function _brc_enc_freq!(e::BinRangeEncoder, cum::UInt32, f::UInt32, tot::UInt32)
    r = e.range ÷ tot
    e.low += UInt64(r) * UInt64(cum)
    e.range = r * f
    while e.range < _RC_TOP
        e.range <<= 8
        _brc_shiftlow!(e)
    end
    return nothing
end

@inline function _brc_rescale!(c0::Vector{UInt32}, c1::Vector{UInt32}, ctx::Int)
    if c0[ctx] + c1[ctx] >= _BRC_MAXTOT
        c0[ctx] = (c0[ctx] >> 1) | 1
        c1[ctx] = (c1[ctx] >> 1) | 1
    end
    return nothing
end

function brc_encode_bit!(e::BinRangeEncoder, b::Bool)
    ctx = e.prev + 1
    c0 = e.c0[ctx]; c1 = e.c1[ctx]; tot = c0 + c1
    if b
        _brc_enc_freq!(e, c0, c1, tot)          # symbol 1 -> [c0, c0+c1)
        e.c1[ctx] += _BRC_INC
    else
        _brc_enc_freq!(e, UInt32(0), c0, tot)   # symbol 0 -> [0, c0)
        e.c0[ctx] += _BRC_INC
    end
    _brc_rescale!(e.c0, e.c1, ctx)
    e.prev = b ? 1 : 0
    return nothing
end

function brc_finish!(e::BinRangeEncoder)
    for _ in 1:5; _brc_shiftlow!(e); end
    return e.out
end

# Finalize the current chunk and reset ALL state (registers + per-context bit
# counts + prev-bit context). Mirror of rc_finish_and_reset! for the copy stream.
function brc_finish_and_reset!(e::BinRangeEncoder)
    bytes = brc_finish!(e)
    e.low = UInt64(0)
    e.range = typemax(UInt32)
    e.cache = UInt8(0)
    e.cachesize = 1
    e.out = UInt8[]
    fill!(e.c0, UInt32(1))
    fill!(e.c1, UInt32(1))
    e.prev = 1
    return bytes
end

mutable struct BinRangeDecoder
    code::UInt32
    range::UInt32
    inp::Vector{UInt8}
    ip::Int
    c0::Vector{UInt32}
    c1::Vector{UInt32}
    prev::Int
end
function BinRangeDecoder(bytes::Vector{UInt8})
    d = BinRangeDecoder(UInt32(0), typemax(UInt32), bytes, 1,
        ones(UInt32, 2), ones(UInt32, 2), 1)
    d.ip += 1
    for _ in 1:4
        d.code = (d.code << 8) | UInt32(get(d.inp, d.ip, 0x00)); d.ip += 1
    end
    return d
end

@inline function _brc_dec_update!(d::BinRangeDecoder, cum::UInt32, f::UInt32)
    d.code -= cum * d.range
    d.range *= f
    while d.range < _RC_TOP
        d.code = (d.code << 8) | UInt32(get(d.inp, d.ip, 0x00)); d.ip += 1
        d.range <<= 8
    end
    return nothing
end

function brc_decode_bit!(d::BinRangeDecoder)
    ctx = d.prev + 1
    c0 = d.c0[ctx]; c1 = d.c1[ctx]; tot = c0 + c1
    d.range ÷= tot
    f = min(d.code ÷ d.range, tot - 1)
    if f < c0
        _brc_dec_update!(d, UInt32(0), c0); d.c0[ctx] += _BRC_INC; b = false
    else
        _brc_dec_update!(d, c0, c1); d.c1[ctx] += _BRC_INC; b = true
    end
    _brc_rescale!(d.c0, d.c1, ctx)
    d.prev = b ? 1 : 0
    return b
end
