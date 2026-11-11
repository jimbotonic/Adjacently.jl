"""
    load_adjacency_mgs3_graph(filename; compact_output=false,
        compact_lists=false, fused_lr=false, reuse_identity=false,
        sorted_fastpath=false, reuse_scratch=false, optimize_decode=false,
        parallel_reconstruct=false, decode_workers=Threads.nthreads(:default))

Decode a directed, non-indexed BG/CS/CG file into `Vector{Vector{T}}`, indexed
by the file's 1-based vertex IDs. Lists are sorted, unique and independently
owned by the caller. No incoming adjacency or graph object is constructed.

This is a full-file preload, **not compressed random access**. Subsequent
`adj[v]` accesses use decompressed RAM; copy a list if the caller needs to
mutate it without changing the retained adjacency.

`compact_output=true` copies completed lists to release spare vector capacity,
at the cost of additional allocation and loading time. This differs from
`compact_lists`, which selects CG's internal local-ID representation. The five
CG decoder switches have the same meaning as in `load_compressed_mgs3_graph`
and are rejected for BG/CS. All switches default to false.

`optimize_decode=true` enables dense output, call-owned scratch and guarded
normalization removal. For CG it also enables the five existing decoder
switches above, falling back internally when the identity K=1 LR invariant
does not hold. BG/CS acceleration requires non-indexed LR context-range;
other supported files use the legacy decoder.

`parallel_reconstruct=true` is an experimental, higher-memory two-pass
decoder for CS LR context-range or CG identity K=1 implicit-range LR
context-range files. It parses entropy streams serially, then reconstructs
reference trees with at most `decode_workers` tasks. Unsupported options
throw an ArgumentError. This is distinct from indexed CG parallel-cluster
decoding and does not enable random access. Staging requires O(vertices +
arcs) extra memory (about 41–45 MB on the measured CNR files). More workers
do not parallelize entropy parsing. File bytes and BPE are unchanged.

Supported: children-mode Fibonacci v3.2 and context-range files (BG/CS v3.3,
CG v3.2), with parameterized or default algorithm headers and any CG cluster
count. An `implicit_ranges` file retains only contiguous cluster membership,
not a mapping back to vertex IDs before an ordering step. Indexed files,
undirected graphs, other backends and a zero-vertex
universe are explicitly rejected. Use `load_compressed_mgs3_graph` for those
formats. The codec's existing trusted-input assumptions still apply; this is
not a validator for untrusted compressed payloads.
The accelerated BG/CS LR path and staged CS/CG paths use call-owned entropy
sources. Other paths retain their underlying shared-hook restrictions;
enabling this option is not a general thread-safety guarantee for all formats.
"""
function load_adjacency_mgs3_graph(filename::AbstractString;
        compact_output::Bool=false, compact_lists::Bool=false,
        fused_lr::Bool=false, reuse_identity::Bool=false,
        sorted_fastpath::Bool=false, reuse_scratch::Bool=false,
        optimize_decode::Bool=false, parallel_reconstruct::Bool=false,
        decode_workers::Int=Threads.nthreads(:default))
    parallel_reconstruct && decode_workers < 1 &&
        throw(ArgumentError("decode_workers must be positive"))
    open(filename, "r") do io
        h = read(io, 12)
        length(h) == 12 || throw(EOFError())
        h[1:3] == UInt8[0x4d, 0x47, 0x53] || throw(ArgumentError("Invalid MGS signature"))
        h[4] == 0x03 || throw(ArgumentError("Adjacency loading requires MGS v3"))
        h[6] >> 4 == 0 || throw(ArgumentError("Adjacency loading requires directed children mode"))
        backend_code = h[6] & 0x0f
        backend_code in (INT_ENCODING_FIBONACCI, INT_ENCODING_CONTEXT_RANGE) ||
            throw(ArgumentError("Adjacency loading supports Fibonacci and context-range only"))
        _, _, backend, b2 = decode_header_flags(h[6], h[7])
        alg = if b2 == ALG_BG || PARAM_BG_BASE <= b2 <= PARAM_BG_MAX
            :bg
        elseif b2 == ALG_CS || PARAM_CS_BASE <= b2 <= PARAM_CS_MAX
            :cs
        elseif b2 == ALG_CG || PARAM_CG_BASE <= b2 <= PARAM_CG_MAX
            :cg
        else
            throw(ArgumentError("Adjacency loading supports BG/CS/CG only"))
        end
        ctx = backend == :context_range
        expected_minor = ctx && alg != :cg ? 0x03 : 0x02
        h[5] == expected_minor || throw(ArgumentError("Unsupported MGS stream-layout version"))
        if alg != :cg && (compact_lists || fused_lr || reuse_identity || sorted_fastpath || reuse_scratch)
            throw(ArgumentError("CG decoder options require a CG file"))
        end
        n = sum(Int(h[7+i]) << (8*(i-1)) for i in 1:5)
        n > 0 || throw(ArgumentError("Adjacency loading requires a nonempty vertex universe"))
        # IDs are 1-based: an exact power of two needs one additional bit.
        T = infer_uint_custom_type(UInt8(ndigits(n; base=2)))
        p = if alg == :bg
            b2 == ALG_BG ? _bg_default_params() : decode_bg_params(b2)
        elseif alg == :cs
            b2 == ALG_CS ? _cs_default_params() : decode_cs_params(b2)
        else
            b2 == ALG_CG ? _cg_default_params() : decode_cg_params(b2; varint=:fibonacci)
        end
        rd = UInt8[]; cp = UInt8[]; resid = UInt8[]
        cmd = nothing; flag = nothing
        if ctx
            count = alg == :cg ? 3 : 5
            lengths = [ltoh(read(io, UInt64)) for _ in 1:count]
            # Check lengths before converting to Int or allocating buffers.
            remaining = UInt64(filesize(io) - position(io))
            for len in lengths
                len <= remaining || throw(EOFError())
                remaining -= len
            end
            structural, rd, cp = [read(io, Int(lengths[i])) for i in 1:3]
            if alg != :cg
                cmd = read(io, Int(lengths[4])); flag = read(io, Int(lengths[5]))
            end
            resid = read(io)
        else
            structural = read(io)
        end
        r = BitReader(structural)
        kw = (; ctx_range=ctx, resid_bytes=resid, refdist_bytes=rd, copy_bytes=cp)
        if parallel_reconstruct
            (ctx && alg in (:cs, :cg)) ||
                throw(ArgumentError("Parallel reconstruction requires CS/CG context-range"))
            f = (; n, T, alg, p, structural, rd, cp, resid, cmd, flag)
            rows = _nonra_parallel_adjacency(f, decode_workers)
            return _own_mgs3_adjacency(rows, n, compact_output, true)
        end
        fast_bgcs = optimize_decode && alg != :cg && ctx && p.lr_split
        rows = if alg == :cg
            decode_level(r, p; T=T, directed=true, coding_scheme=:children, kw...,
                compact_lists=compact_lists || optimize_decode, fused_lr=fused_lr || optimize_decode,
                reuse_identity=reuse_identity || optimize_decode,
                sorted_fastpath=sorted_fastpath || optimize_decode,
                reuse_scratch=reuse_scratch || optimize_decode,
                skip_normalize=optimize_decode, direct_output=optimize_decode)
        else
            decoder = alg == :bg ? read_greedy_graph_data : read_cmdstream_graph_data
            decoder(r, T(n), :children, T; integer_encoding=backend, p..., kw...,
                cmd_bytes=cmd, flag_bytes=flag, optimize_decode=fast_bgcs)
        end
        _own_mgs3_adjacency(rows, n, compact_output,
            fast_bgcs || (alg == :cg && optimize_decode))
    end
end

# Function barrier keeps the per-vertex loop specialized on the decoded ID type.
# Ownership transfers only after decoding has consumed all back-references.
function _own_mgs3_adjacency(rows::Dict{T,Vector{T}}, n::Int, compact::Bool, normalized::Bool=false) where T
    adj = Vector{Vector{T}}(undef, n)
    for v in 1:n
        ns = get(rows, T(v), nothing)
        list = ns === nothing ? T[] : ns
        normalized || Compression.CG._cg_sorted_unique!(list)
        adj[v] = compact ? copy(list) : list
    end
    return adj
end

function _own_mgs3_adjacency(rows::Vector{Vector{T}}, n::Int, compact::Bool, normalized::Bool=false) where T
    length(rows) == n || throw(ArgumentError("Wrong output universe"))
    if !normalized
        for list in rows; Compression.CG._cg_sorted_unique!(list); end
    end
    compact ? copy.(rows) : rows
end

function _own_mgs3_adjacency(rows::Compression._PrefixLists{T}, n::Int,
        compact::Bool, normalized::Bool=false) where T
    rows.first_vertex == 1 || throw(ArgumentError("Incomplete output universe"))
    _own_mgs3_adjacency(rows.lists, n, compact, normalized)
end

include("mgs_nonra_parallel.jl")
