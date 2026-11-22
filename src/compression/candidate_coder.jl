# Candidate-bit graph coder.
#
# Rows are coded in id order. For row v the out-degree is coded first; arcs are then found by
#   (1) an optional reference copy from one earlier row inside a window,
#   (2) candidate bits: targets proposed from already-decoded structure (reciprocal arcs, co-cited targets,
#       2-hop out-out / sibling targets) are tested one adaptive bit each, the probability coming from a
#       fixed-point logistic mix of several context models,
#   (3) remainder arcs, coded in ascending order as gaps (or as Bernoulli cell bits on small graphs); every
#       remainder arc optionally seeds new candidate rounds.
# A row stops as soon as its degree is reached. Encoder and decoder run the same function (`_cand_run!`), so
# every modelling decision is shared; only integer arithmetic drives probabilities.
#
# Included into `module Compression` after range_coder.jl (uses the BinRange registers).

Base.@kwdef struct CandParams
    ref_window::Int     = 64      # reference-copy window (0 = no reference copy)
    ref_gain::Int       = 4       # a reference row r is used if overlap*ref_gain - length(r) - ref_bias > 0 (best such r)
    ref_bias::Int       = 4
    gen_recip::Bool     = true    # candidates: earlier in-neighbours of v
    gen_cocite::Bool    = true    # candidates: out-neighbours of earlier in-neighbours of v
    gen_2hop::Bool      = true    # candidates: out-neighbours of found targets, and of their in-neighbours
    gen_innbr::Bool     = false   # candidates: in-neighbours of found targets
    cap_a::Int          = 256     # most recent in-neighbours of v that are expanded
    cap_in::Int         = 32      # most recent in-neighbours of a found target that are expanded
    cap_row::Int        = 512     # rows longer than this are never expanded
    budget_base::Int    = 2048    # generation budget per row = budget_base + budget_per_arc * degree (list elements visited)
    budget_per_arc::Int = 1024
    min_score::Int      = 1       # a non-reciprocal candidate needs this many proposals in its round
    seed::Bool          = true    # remainder arcs seed further candidate rounds
    mix_lr::Int         = 6       # mixer learning rate (weight step = stretch*err*mix_lr >> 16)
    mix_w0::Int         = 19661   # initial mixer weight (65536 = 1.0)
    cnt_limit::Int      = 1023    # adaptation limit of the context counters
    cnt_p0::Int         = 209715  # initial counter probability (2^22 = 1.0)
    remainder::Symbol   = :gap    # :gap, :popgap (popularity-weighted target inside the gap bucket), or :cell (O(n) work per row: small graphs only)
    cell_ds::Int        = 2       # :cell distance bucket width (bits)
    cell_kappa::Int     = 200     # :cell prior strength (pairs)
    # Fields below were added after format version 1 shipped. They are serialised only when one of them differs from
    # its default, so default parameters keep their original bytes.
    pop_split::Int      = 5       # :popgap: gaps below 2^pop_split keep an exact adaptive token; larger gaps use the popularity code
    pop_prior::Int      = 1       # :popgap: initial hits and exposure of every popularity bucket
    gap_deg_buckets::Int = 1      # :gap/:popgap: the gap-token model is conditioned on min(floor(log2(degree+1)), gap_deg_buckets-1)
end

const _CAND_VERSION = 0x01
const _CAND_FIELDS = (:ref_window, :ref_gain, :ref_bias, :gen_recip, :gen_cocite, :gen_2hop, :gen_innbr, :cap_a, :cap_in,
                      :cap_row, :budget_base, :budget_per_arc, :min_score, :seed, :mix_lr, :mix_w0, :cnt_limit, :cnt_p0,
                      :remainder, :cell_ds, :cell_kappa, :pop_split, :pop_prior, :gap_deg_buckets)
const _CAND_NFIELDS_V1 = 21          # fields present in the first shipped layout
const _CAND_REMAINDERS = (:gap, :cell, :popgap)
const _CAND_NPOP = 22                # popularity buckets: min(floor(log2(1 + in-degree so far)), 21)

function _cand_check(p::CandParams)
    p.remainder in _CAND_REMAINDERS || throw(ArgumentError("CandParams.remainder must be :gap, :popgap or :cell"))
    (0 <= p.pop_split <= 6 && 1 <= p.pop_prior <= (1 << 20) && 1 <= p.gap_deg_buckets <= 32) ||
        throw(ArgumentError("CandParams: popgap constant out of range"))
    (0 <= p.ref_window && p.cap_a >= 0 && p.cap_in >= 0 && p.cap_row >= 0 && p.budget_base >= 0 && p.budget_per_arc >= 0) ||
        throw(ArgumentError("CandParams: negative cap"))
    (1 <= p.mix_lr <= 4096 && 0 <= p.mix_w0 <= (1 << 20) && 1 <= p.cnt_limit <= 65535 && 1 <= p.cnt_p0 < (1 << 22)) ||
        throw(ArgumentError("CandParams: mixer/counter constant out of range"))
    (1 <= p.cell_ds <= 32 && 0 <= p.cell_kappa <= (1 << 30) && 0 <= p.min_score && abs(p.ref_gain) < (1 << 20) && abs(p.ref_bias) < (1 << 20)) ||
        throw(ArgumentError("CandParams: constant out of range"))
    return p
end

# Layout: version byte, field count byte, then one little-endian Int32 per field (Bool 0/1; remainder 0=:gap 1=:cell
# 2=:popgap). The field count is _CAND_NFIELDS_V1 when every later field has its default value, else the full count;
# a reader fills missing trailing fields with their defaults.
function cand_params_bytes(p::CandParams)::Vector{UInt8}
    _cand_check(p)
    d = CandParams()
    # The 21-field (v1) layout is kept only for parameters a v1 reader decodes correctly: a v1 reader maps any remainder
    # code other than 1 to :gap, so :popgap always uses the extended layout, which v1 readers reject (field count).
    v1ok = p.remainder !== :popgap && all(getfield(p, f) == getfield(d, f) for f in _CAND_FIELDS[_CAND_NFIELDS_V1+1:end])
    nf = v1ok ? _CAND_NFIELDS_V1 : length(_CAND_FIELDS)
    b = UInt8[_CAND_VERSION, UInt8(nf)]
    for f in _CAND_FIELDS[1:nf]
        x = getfield(p, f)
        v = x isa Symbol ? Int32(findfirst(==(x), _CAND_REMAINDERS) - 1) : Int32(x)
        for k in 0:3; push!(b, UInt8((reinterpret(UInt32, v) >> (8k)) & 0xFF)); end
    end
    return b
end

function cand_params_from_bytes(b::AbstractVector{UInt8})::CandParams
    length(b) >= 2 || throw(ArgumentError("candidate-coder parameters: truncated"))
    b[1] == _CAND_VERSION || throw(ArgumentError("candidate-coder parameters: unsupported version $(b[1])"))
    nf = Int(b[2])
    (_CAND_NFIELDS_V1 <= nf <= length(_CAND_FIELDS) && length(b) == 2 + 4nf) || throw(ArgumentError("candidate-coder parameters: bad length"))
    d = CandParams()
    vals = Any[getfield(d, f) for f in _CAND_FIELDS]
    for i in 1:nf
        f = _CAND_FIELDS[i]
        u = UInt32(0); for k in 0:3; u |= UInt32(b[2 + 4(i-1) + k + 1]) << (8k); end
        v = Int(reinterpret(Int32, u))
        T = fieldtype(CandParams, f)
        if T === Symbol
            0 <= v < length(_CAND_REMAINDERS) || throw(ArgumentError("candidate-coder parameters: unknown remainder mode"))
            vals[i] = _CAND_REMAINDERS[v+1]
        else
            vals[i] = T === Bool ? (v != 0) : v
        end
    end
    return _cand_check(CandParams(vals...))
end

# ---------------- fixed-point logistic tables (lpaq-style) ----------------
const _CAND_SQ = Int32[1, 2, 3, 6, 10, 16, 27, 45, 73, 120, 194, 310, 488, 747, 1101, 1546, 2047, 2549, 2994, 3348, 3607,
                       3785, 3901, 3975, 4022, 4050, 4068, 4079, 4085, 4089, 4092, 4093, 4094]
# squash(d) ~ 4096 / (1 + exp(-d/256)), d clamped to [-2047, 2047]
@inline function _cand_squash(d::Int)
    d > 2047 && return 4095
    d < -2047 && return 0
    w = d & 127
    i = (d >> 7) + 16
    @inbounds return (Int(_CAND_SQ[i+1]) * (128 - w) + Int(_CAND_SQ[i+2]) * w + 64) >> 7
end
const _CAND_ST = let t = zeros(Int16, 4096), pi = 0
    for x in -2047:2047
        v = _cand_squash(x)
        for i in pi:v; t[i+1] = Int16(x); end
        pi = v + 1
    end
    for i in pi:4095; t[i+1] = Int16(2047); end
    t
end

# ---------------- coder front-ends over the BinRange registers ----------------
struct _CandEnc; e::BinRangeEncoder; end
struct _CandDec; d::BinRangeDecoder; end

# Code one bit; `z` = frequency of 0 out of `tot` (0 < z < tot <= 2^16). Encoder returns `b`, decoder the decoded bit.
@inline function _cand_bit!(c::_CandEnc, z::UInt32, tot::UInt32, b::Bool)
    b ? _brc_enc_freq!(c.e, z, tot - z, tot) : _brc_enc_freq!(c.e, UInt32(0), z, tot)
    return b
end
@inline function _cand_bit!(c::_CandDec, z::UInt32, tot::UInt32, b::Bool)
    d = c.d
    d.range ÷= tot
    f = min(d.code ÷ d.range, tot - UInt32(1))
    if f < z
        _brc_dec_update!(d, UInt32(0), z); return false
    else
        _brc_dec_update!(d, z, tot - z); return true
    end
end
# Count-based adaptive bit (same update rule as brc_encode_ctx_bit!).
@inline function _cand_cbit!(c, c0::Vector{UInt32}, c1::Vector{UInt32}, ctx::Int, b::Bool)
    @inbounds z = c0[ctx]; @inbounds o = c1[ctx]
    r = _cand_bit!(c, z, z + o, b)
    @inbounds r ? (c1[ctx] += _BRC_INC) : (c0[ctx] += _BRC_INC)
    _brc_rescale!(c0, c1, ctx)
    return r
end

# Adaptive integer model: hybrid token (range_coder.jl, K=5 L=1) through an 8-level binary tree per context,
# remaining mantissa bits bypassed.
struct _CandVal; c0::Vector{UInt32}; c1::Vector{UInt32}; end
_CandVal(nctx::Int) = _CandVal(ones(UInt32, 256nctx), ones(UInt32, 256nctx))
function _cand_val!(c, m::_CandVal, ctx::Int, v::Int)
    enc = c isa _CandEnc
    token = 0; nraw = 0; raw = UInt64(0)
    if enc
        token, nraw, raw = _rc_hybrid_encode(UInt64(v))
    end
    node = 1; base = (ctx - 1) * 256
    for i in 7:-1:0
        b = _cand_cbit!(c, m.c0, m.c1, base + node, enc && ((token >> i) & 1) == 1)
        node = 2node + (b ? 1 : 0)
    end
    enc || (token = node - 256)
    if token < 32
        return token
    end
    token < _RC_NSYM || error("candidate coder: corrupt stream (token)")
    logv = ((token - 32) >> 1) + 5
    nraw = logv - 1
    r = UInt64(0)
    for i in nraw-1:-1:0
        b = _cand_bit!(c, UInt32(1), UInt32(2), enc && ((raw >> i) & 1) == 1)
        r = (r << 1) | (b ? 1 : 0)
    end
    return Int((UInt64(1) << logv) | (r << 1) | UInt64((token - 32) & 1))
end

# One token 0..255 through the same binary tree (no mantissa).
function _cand_tok!(c, m::_CandVal, ctx::Int, token::Int)
    enc = c isa _CandEnc
    node = 1; base = (ctx - 1) * 256
    for i in 7:-1:0
        b = _cand_cbit!(c, m.c0, m.c1, base + node, enc && ((token >> i) & 1) == 1)
        node = 2node + (b ? 1 : 0)
    end
    return node - 256
end
# Symbol s (1-based) of a computed frequency table, total `tot` <= 2^16; decoder returns the decoded symbol.
@inline function _cand_sym!(c::_CandEnc, f::Vector{UInt32}, ns::Int, tot::UInt32, s::Int)
    cum = UInt32(0)
    @inbounds for i in 1:s-1; cum += f[i]; end
    @inbounds _brc_enc_freq!(c.e, cum, f[s], tot)
    return s
end
@inline function _cand_sym!(c::_CandDec, f::Vector{UInt32}, ns::Int, tot::UInt32, s::Int)
    d = c.d
    d.range ÷= tot
    t = min(d.code ÷ d.range, tot - UInt32(1))
    cum = UInt32(0); k = 1
    @inbounds while k < ns && cum + f[k] <= t; cum += f[k]; k += 1; end
    @inbounds f[k] > 0 || error("candidate coder: corrupt stream (symbol)")
    @inbounds _brc_dec_update!(d, cum, f[k])
    return k
end
# Uniform integer k in 0..N-1 (N < 2^32).
function _cand_uniform!(c, N::Int, k::Int)
    enc = c isa _CandEnc
    if N <= 65536
        N == 1 && return 0
        if enc
            _brc_enc_freq!(c.e, UInt32(k), UInt32(1), UInt32(N)); return k
        else
            d = c.d; d.range ÷= UInt32(N)
            t = min(d.code ÷ d.range, UInt32(N - 1)); _brc_dec_update!(d, t, UInt32(1)); return Int(t)
        end
    end
    nh = ((N - 1) >> 16) + 1
    kh = _cand_uniform!(c, nh, k >> 16)
    nl = kh == nh - 1 ? ((N - 1) & 0xFFFF) + 1 : 65536
    kl = _cand_uniform!(c, nl, k & 0xFFFF)
    return (kh << 16) | kl
end
# Fenwick trees over ids (Int32 counts).
@inline function _cand_fadd!(t::Vector{Int32}, i::Int, d::Int32)
    n = length(t)
    @inbounds while i <= n; t[i] += d; i += i & -i; end
    return nothing
end
@inline function _cand_fsum(t::Vector{Int32}, i::Int)
    s = 0
    @inbounds while i > 0; s += t[i]; i -= i & -i; end
    return s
end
@inline _cand_fcnt(t::Vector{Int32}, lo::Int, hi::Int) = lo > hi ? 0 : _cand_fsum(t, hi) - _cand_fsum(t, lo - 1)
# smallest i with prefix sum >= k (k >= 1); returns length(t)+1 if none
@inline function _cand_fsel(t::Vector{Int32}, k::Int)
    n = length(t); pos = 0
    step = n == 0 ? 0 : (1 << (63 - leading_zeros(n)))
    @inbounds while step > 0
        q = pos + step
        if q <= n && t[q] < k; pos = q; k -= t[q]; end
        step >>= 1
    end
    return pos + 1
end
@inline _cand_pb(indeg::Int) = min(63 - leading_zeros(indeg + 1), _CAND_NPOP - 1)

@inline _cand_bl(x::Int) = x <= 0 ? 0 : 64 - leading_zeros(x)
@inline _cand_bk6(s::Int) = s <= 3 ? s : s <= 5 ? 4 : 5
@inline _cand_bk7(s::Int) = s <= 3 ? s : s <= 5 ? 4 : s <= 9 ? 5 : 6
@inline _cand_nbk(q::Int) = q <= 2 ? q - 1 : q <= 4 ? 2 : q <= 8 ? 3 : 4          # q >= 1
@inline _cand_hl(z::Int) = (b = _cand_bl(z); 2b + (b >= 2 ? ((z >> (b - 2)) & 1) : 0))   # ~2*log2(z), integer

const _CAND_NIN = 9
const _CAND_SIZES = (2*7*3, 2*216, 8*7*2, 5*8*2, 5*7*2, 8*7*5, 48*2*3, 5*8, 1)
const _CAND_NSETS = 14

# The shared encode/decode pass. Encoder: `rows` is the input. Decoder: `rows` is filled row by row.
function _cand_run!(c, rows::Vector{Vector{UInt32}}, n::Int, P::CandParams)
    enc = c isa _CandEnc
    L = _cand_bl(n)
    inl = [UInt32[] for _ in 1:n]                    # in-neighbours among already coded rows, ascending
    tested = zeros(Int32, n); found = zeros(Int32, n)
    mark = enc ? zeros(Int32, n) : Int32[]
    s1 = zeros(Int32, n); s2 = zeros(Int32, n); s3 = zeros(Int32, n); isR = fill(false, n)
    touched = UInt32[]; keys = UInt64[]; newK = UInt32[]
    # models
    mdeg = _CandVal(2); moff = _CandVal(1); mgap = _CandVal(100 * P.gap_deg_buckets)
    popg = P.remainder === :popgap
    fen = Vector{Vector{Int32}}()
    if popg
        t1 = zeros(Int32, n)                       # every id starts in bucket 0
        for i in 1:n; t1[i] += Int32(1); j = i + (i & -i); j <= n && (t1[j] += t1[i]); end
        push!(fen, t1); for _ in 2:_CAND_NPOP; push!(fen, zeros(Int32, n)); end
    end
    phits = fill(UInt64(P.pop_prior), _CAND_NPOP); pexpo = fill(UInt64(P.pop_prior), _CAND_NPOP)
    pNA = zeros(Int, _CAND_NPOP); pNB = zeros(Int, _CAND_NPOP); pfr = zeros(UInt32, _CAND_NPOP); pq = zeros(UInt128, _CAND_NPOP)
    maxpb = 0
    self0 = ones(UInt32, 1); self1 = ones(UInt32, 1)
    href0 = ones(UInt32, 2); href1 = ones(UInt32, 2); same0 = ones(UInt32, 2); same1 = ones(UInt32, 2)
    cop0 = ones(UInt32, 80); cop1 = ones(UInt32, 80)
    offs = zeros(Int, _CAND_NIN); for i in 2:_CAND_NIN; offs[i] = offs[i-1] + _CAND_SIZES[i-1]; end
    ncm = sum(_CAND_SIZES)
    cp = fill(Int32(P.cnt_p0), ncm); cn = zeros(UInt16, ncm)
    wts = fill(Int32(P.mix_w0), _CAND_NIN, _CAND_NSETS)
    idx = zeros(Int, _CAND_NIN); st = zeros(Int, _CAND_NIN)
    cell = P.remainder === :cell
    nD = (L + P.cell_ds) ÷ P.cell_ds + 1
    CE = cell ? zeros(UInt32, nD, 7, 8, 5) : zeros(UInt32, 0, 0, 0, 0)
    CP = cell ? zeros(UInt32, nD, 7, 8, 5) : zeros(UInt32, 0, 0, 0, 0)
    Etot = UInt64(1); Ptot = UInt64(1000)
    lim = UInt16(P.cnt_limit); lr = P.mix_lr
    prevd = 0; lastref = 0; prevhadref = false
    W = P.ref_window

    for v in 1:n
        d = _cand_val!(c, mdeg, prevd == 0 ? 1 : 2, enc ? length(rows[v]) : 0)
        prevd = d
        cur = enc ? rows[v] : UInt32[]
        if d == 0
            enc || (rows[v] = cur)
            prevhadref = false
            continue
        end
        d <= n || error("candidate coder: corrupt stream (degree)")
        if enc
            @inbounds for x in cur; mark[x] = v; end
        end
        self = _cand_cbit!(c, self0, self1, 1, enc && mark[v] == v)
        need = d - (self ? 1 : 0)
        @inbounds tested[v] = v
        if self
            @inbounds found[v] = v
            enc || push!(cur, UInt32(v))
        end
        empty!(newK)
        # ---- reference copy
        if W > 0 && need > 0
            r = 0
            if enc
                bs = 0
                @inbounds for q in v-1:-1:max(1, v - W)
                    rq = rows[q]; isempty(rq) && continue
                    o = 0
                    for x in rq; o += (mark[x] == v && x != v) ? 1 : 0; end
                    sc = o * P.ref_gain - length(rq) - P.ref_bias
                    if sc > bs; bs = sc; r = q; end
                end
            end
            has = _cand_cbit!(c, href0, href1, prevhadref ? 1 : 2, r != 0)
            if has
                same = _cand_cbit!(c, same0, same1, prevhadref ? 1 : 2, enc && r == lastref)
                if same
                    r = lastref
                else
                    r = v - _cand_val!(c, moff, 1, enc ? v - r : 0)
                end
                (1 <= r < v) || error("candidate coder: corrupt stream (reference)")
                lastref = r; prevhadref = true
                pb = 1
                @inbounds for x32 in rows[r]
                    x = Int(x32)
                    x == v && continue
                    need == 0 && break
                    tested[x] = v
                    pop = min(_cand_bl(length(inl[x])), 7)
                    b = _cand_cbit!(c, cop0, cop1, 1 + pb + 2 * (pop + 8 * _cand_nbk(need)), enc && mark[x] == v)
                    pb = b ? 1 : 0
                    if b
                        found[x] = v; push!(newK, x32); need -= 1
                        enc || push!(cur, x32)
                    end
                end
            else
                prevhadref = false
            end
        else
            prevhadref = false
        end
        # ---- candidate rounds interleaved with remainder arcs
        Lv = inl[v]; hiA = length(Lv)
        first = true; rnd = 1; ri = 1; nrem = 0; prevrem = 0; prevbl = min(_cand_bl(d + 1), 23); cpos = 0
        gen = 0; budget = P.budget_base + P.budget_per_arc * d
        anygen = P.gen_recip || P.gen_cocite || P.gen_2hop || P.gen_innbr
        while need > 0
            while need > 0 && anygen && (first || !isempty(newK))
                empty!(touched)
                if first
                    @inbounds for i in hiA:-1:max(1, hiA - P.cap_a + 1)
                        u = Int(Lv[i])
                        if P.gen_recip
                            (s1[u] == 0 && s2[u] == 0 && s3[u] == 0 && !isR[u]) && push!(touched, UInt32(u))
                            isR[u] = true
                        end
                        ru = rows[u]
                        if P.gen_cocite && gen < budget && length(ru) <= P.cap_row
                            gen += length(ru)
                            for x32 in ru
                                x = Int(x32)
                                (s1[x] == 0 && s2[x] == 0 && s3[x] == 0 && !isR[x]) && push!(touched, x32)
                                s1[x] += Int32(1)
                            end
                        end
                    end
                    first = false
                end
                if P.gen_2hop || P.gen_innbr
                    @inbounds for w32 in newK
                        w = Int(w32)
                        if P.gen_2hop && w < v
                            rw = rows[w]
                            if gen < budget && length(rw) <= P.cap_row
                                gen += length(rw)
                                for x32 in rw
                                    x = Int(x32)
                                    (s1[x] == 0 && s2[x] == 0 && s3[x] == 0 && !isR[x]) && push!(touched, x32)
                                    s2[x] += Int32(1)
                                end
                            end
                        end
                        Iw = inl[w]; hw = length(Iw)
                        for i in hw:-1:max(1, hw - P.cap_in + 1)
                            u32 = Iw[i]; u = Int(u32)
                            if P.gen_innbr
                                (s1[u] == 0 && s2[u] == 0 && s3[u] == 0 && !isR[u]) && push!(touched, u32)
                                s2[u] += Int32(1)
                            end
                            ru = rows[u]
                            if P.gen_2hop && gen < budget && length(ru) <= P.cap_row
                                gen += length(ru)
                                for x32 in ru
                                    x = Int(x32)
                                    (s1[x] == 0 && s2[x] == 0 && s3[x] == 0 && !isR[x]) && push!(touched, x32)
                                    s3[x] += Int32(1)
                                end
                            end
                        end
                    end
                end
                empty!(newK); empty!(keys)
                @inbounds for (j, x32) in enumerate(touched)
                    x = Int(x32)
                    tot = Int(s1[x]) + Int(s2[x]) + Int(s3[x])
                    if tested[x] == v || found[x] == v || (!isR[x] && tot < P.min_score)
                        s1[x] = 0; s2[x] = 0; s3[x] = 0; isR[x] = false
                        continue
                    end
                    score = (isR[x] ? (1 << 24) : 0) + min(tot, (1 << 24) - 1)
                    push!(keys, (UInt64((1 << 26) - score) << 32) | UInt64(j))
                end
                sort!(keys)
                @inbounds for k in keys
                    x32 = touched[Int(k & 0xFFFFFFFF)]; x = Int(x32)
                    r = isR[x]; a = Int(s1[x]); b = Int(s2[x]); cc = Int(s3[x])
                    s1[x] = 0; s2[x] = 0; s3[x] = 0; isR[x] = false
                    need == 0 && continue
                    tested[x] = v
                    tot = a + b + cc; t7 = _cand_bk7(tot); rr = r ? 1 : 0; rd = min(rnd, 3) - 1
                    seen = length(inl[x])
                    pop = min(_cand_bl(seen), 7); bl = _cand_bl(abs(x - v)); db = min(4, (5bl) ÷ (L + 1))
                    xo = x < v ? min(_cand_bl(length(rows[x])), 5) + 1 : 0
                    nb = _cand_nbk(need)
                    lift = clamp(_cand_hl(1 + tot) - _cand_hl(1 + seen) + 40, 0, 47)
                    idx[1] = offs[1] + 1 + rr + 2 * (t7 + 7 * rd)
                    idx[2] = offs[2] + 1 + rr + 2 * (_cand_bk6(a) + 6 * (_cand_bk6(b) + 6 * _cand_bk6(cc)))
                    idx[3] = offs[3] + 1 + pop + 8 * (xo + 7 * rr)
                    idx[4] = offs[4] + 1 + db + 5 * (pop + 8 * (tot > 1 ? 1 : 0))
                    idx[5] = offs[5] + 1 + nb + 5 * (t7 + 7 * rr)
                    idx[6] = offs[6] + 1 + pop + 8 * (t7 + 7 * db)
                    idx[7] = offs[7] + 1 + lift + 48 * (rr + 2 * rd)
                    idx[8] = offs[8] + 1 + nb + 5 * pop
                    idx[9] = offs[9] + 1
                    set = 1 + t7 + 7 * rr
                    dot = 0
                    for i in 1:_CAND_NIN
                        s = Int(_CAND_ST[clamp(Int(cp[idx[i]]) >> 10, 1, 4095) + 1]); st[i] = s
                        dot += Int(wts[i, set]) * s
                    end
                    p = clamp(_cand_squash(clamp(dot >> 16, -2047, 2047)), 1, 4095)
                    hit = _cand_bit!(c, UInt32(4096 - p), UInt32(4096), enc && mark[x] == v)
                    err = ((hit ? 4096 : 0) - p) * lr
                    tgt = hit ? (1 << 22) : 0
                    for i in 1:_CAND_NIN
                        wts[i, set] = Int32(clamp(Int(wts[i, set]) + ((st[i] * err + 0x8000) >> 16), -(1 << 24), 1 << 24))
                        j = idx[i]; nn = Int(cn[j])
                        cp[j] = Int32(clamp(Int(cp[j]) + div((tgt - Int(cp[j])) * 2, 2nn + 3), 1, (1 << 22) - 1))
                        cn[j] < lim && (cn[j] += 0x0001)
                    end
                    if hit
                        push!(newK, x32); found[x] = v; need -= 1
                        enc || push!(cur, x32)
                    end
                end
                rnd += 1
            end
            need == 0 && break
            # ---- next remainder arc = smallest arc not yet found
            x = 0
            if enc
                @inbounds while found[cur[ri]] == v; ri += 1; end
                @inbounds x = Int(cur[ri])
            end
            if cell
                nbk = _cand_nbk(need) + 1
                g16 = (Etot << 16) ÷ Ptot
                kap = UInt64(P.cell_kappa)
                fx = 0
                @inbounds for y in cpos+1:n
                    tested[y] == v && continue
                    i1 = (_cand_bl(abs(y - v)) + P.cell_ds - 1) ÷ P.cell_ds + 1
                    i2 = (y < v ? min(_cand_bl(length(rows[y])), 5) + 1 : 0) + 1
                    i3 = min(_cand_bl(length(inl[y])), 7) + 1
                    p16 = clamp(((UInt64(CE[i1, i2, i3, nbk]) << 16) + kap * g16) ÷ (UInt64(CP[i1, i2, i3, nbk]) + kap + 1), UInt64(1), UInt64(65535))
                    hit = _cand_bit!(c, UInt32(65536 - p16), UInt32(65536), enc && y == x)
                    if CP[i1, i2, i3, nbk] < 0xFFFFFFF0
                        CP[i1, i2, i3, nbk] += 0x00000001; Ptot += 1
                        if hit; CE[i1, i2, i3, nbk] += 0x00000001; Etot += 1; end
                    end
                    if hit; fx = y; break; end
                    tested[y] = v
                end
                fx != 0 || error("candidate coder: corrupt stream (cell remainder)")
                x = fx; cpos = x
            else
                g = 0
                if enc
                    g = nrem == 0 ? (x >= v ? 2 * (x - v) : 2 * (v - x) - 1) : x - prevrem
                end
                gctx = 100 * min(_cand_bl(d + 1) - 1, P.gap_deg_buckets - 1) + (nrem == 0 ? 1 : min(nrem + 1, 3)) * 25 + prevbl + 1
                if popg
                    sp = 1 << P.pop_split
                    tok = _cand_tok!(c, mgap, gctx, enc ? (g < sp ? g : sp + (_cand_bl(g) - 1 - P.pop_split)) : 0)
                    if tok < sp
                        g = tok
                        enc || (x = nrem == 0 ? (iseven(g) ? v + (g >> 1) : v - ((g + 1) >> 1)) : prevrem + g)
                    else
                        dd = tok - sp + P.pop_split
                        dd <= 40 || error("candidate coder: corrupt stream (gap bucket)")
                        lo = 1 << dd; hi = (1 << (dd + 1)) - 1
                        # id intervals of the bucket: A (ids above the base), B (ids below v, first arc only)
                        a1 = 0; a2 = -1; b1 = 0; b2 = -1
                        if nrem == 0
                            a1 = v + ((lo + 1) >> 1); a2 = min(n, v + (hi >> 1))
                            b1 = max(1, v - ((hi + 1) >> 1)); b2 = v - ((lo + 2) >> 1)
                        else
                            a1 = prevrem + lo; a2 = min(n, prevrem + hi)
                        end
                        Q = UInt128(0)
                        @inbounds for b in 1:maxpb+1
                            na = _cand_fcnt(fen[b], a1, a2); nb2 = _cand_fcnt(fen[b], b1, b2)
                            pNA[b] = na; pNB[b] = nb2
                            q = ((UInt128(na + nb2) * phits[b]) << 32) ÷ pexpo[b]
                            pq[b] = q; Q += q
                        end
                        Q > 0 || error("candidate coder: corrupt stream (empty gap bucket)")
                        tot = UInt32(0)
                        @inbounds for b in 1:maxpb+1
                            f = pq[b] == 0 ? UInt32(0) : UInt32(max(UInt128(1), (pq[b] * 32000) ÷ Q))
                            (f == 0 && pNA[b] + pNB[b] > 0) && (f = UInt32(1))
                            pfr[b] = f; tot += f
                        end
                        bx = _cand_sym!(c, pfr, maxpb + 1, tot, enc ? _cand_pb(length(inl[x])) + 1 : 0)
                        @inbounds na = pNA[bx]; @inbounds nt = na + pNB[bx]
                        k = 0
                        if enc
                            @inbounds k = x >= a1 && x <= a2 ? _cand_fcnt(fen[bx], a1, x) - 1 : na + _cand_fcnt(fen[bx], b1, x) - 1
                        end
                        k = _cand_uniform!(c, nt, k)
                        if !enc
                            @inbounds x = k < na ? _cand_fsel(fen[bx], _cand_fsum(fen[bx], a1 - 1) + k + 1) :
                                                   _cand_fsel(fen[bx], _cand_fsum(fen[bx], b1 - 1) + (k - na) + 1)
                            g = nrem == 0 ? (x >= v ? 2 * (x - v) : 2 * (v - x) - 1) : x - prevrem
                        end
                        @inbounds phits[bx] += 1
                        @inbounds for b in 1:maxpb+1
                            pexpo[b] += UInt64(pNA[b] + pNB[b])
                            if pexpo[b] > (UInt64(1) << 44)
                                pexpo[b] >>= 1; phits[b] = max(UInt64(1), phits[b] >> 1)
                            end
                        end
                    end
                else
                    g = _cand_val!(c, mgap, gctx, g)
                    enc || (x = nrem == 0 ? (iseven(g) ? v + (g >> 1) : v - ((g + 1) >> 1)) : prevrem + g)
                end
                if !enc
                    (1 <= x <= n && found[x] != v && x != v) || error("candidate coder: corrupt stream (remainder)")
                end
                prevbl = min(_cand_bl(g), 23)
            end
            @inbounds found[x] = v
            @inbounds tested[x] = v
            need -= 1; nrem += 1; prevrem = x
            enc || push!(cur, UInt32(x))
            empty!(newK); P.seed && push!(newK, UInt32(x))
        end
        if !enc
            sort!(cur); rows[v] = cur
        end
        if popg
            @inbounds for x32 in cur
                x = Int(x32); l = length(inl[x]); ob = _cand_pb(l); nb2 = _cand_pb(l + 1)
                if nb2 != ob
                    _cand_fadd!(fen[ob+1], x, Int32(-1)); _cand_fadd!(fen[nb2+1], x, Int32(1))
                    nb2 > maxpb && (maxpb = nb2)
                end
                push!(inl[x], UInt32(v))
            end
        else
            @inbounds for x in cur; push!(inl[x], UInt32(v)); end
        end
    end
    return rows
end

"""
    cand_encode(rows; params=CandParams()) -> Vector{UInt8}

Encode a directed graph given as 1-based, sorted, duplicate-free out-neighbour lists (`n = length(rows)`;
self-loops and empty rows allowed). `n` and `params` are not stored in the output.
"""
function cand_encode(rows::Vector{Vector{UInt32}}; params::CandParams=CandParams())::Vector{UInt8}
    _cand_check(params)
    n = length(rows)
    for r in rows
        prev = UInt32(0)
        for x in r
            (prev < x <= n) || throw(ArgumentError("cand_encode: rows must be sorted, unique and within 1..n"))
            prev = x
        end
    end
    e = BinRangeEncoder()
    _cand_run!(_CandEnc(e), rows, n, params)
    return brc_finish!(e)
end

"""
    cand_decode(bytes, n; params=CandParams()) -> Vector{Vector{UInt32}}

Inverse of [`cand_encode`](@ref); `n` and `params` must be the ones used for encoding.
"""
function cand_decode(bytes::Vector{UInt8}, n::Int; params::CandParams=CandParams())::Vector{Vector{UInt32}}
    _cand_check(params)
    rows = Vector{Vector{UInt32}}(undef, n)
    n == 0 && return rows
    d = BinRangeDecoder(bytes)
    _cand_run!(_CandDec(d), rows, n, params)
    return rows
end
