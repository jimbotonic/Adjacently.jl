# Explicit-source readers: CG tasks never select a process-global entropy source.
function _cg_sorted_unique!(values::Vector)
    # Sort only if necessary, then compact adjacent duplicates in place. Julia
    # already specializes unique! for sorted integers; this mainly avoids the
    # redundant sort dispatch/check, not a guaranteed hash-set allocation.
    issorted(values) || sort!(values)
    isempty(values) && return values
    out = 1
    @inbounds for i in 2:length(values)
        if values[i] != values[out]
            out += 1
            values[out] = values[i]
        end
    end
    resize!(values, out)
    return values
end

_cg_read_value(r, encoding, ::Type{T}; source=nothing) where {T<:Unsigned} =
    Compression._read_encoded_value(r, encoding, T, source)

function _cg_read_small_count(r, encoding, ::Type{T}; source=nothing) where {T<:Unsigned}
    tag = 2Int(read_bit(r)) + Int(read_bit(r))
    return tag < 3 ? T(tag) : _cg_read_value(r, encoding, T; source=source)
end

function _cg_read_delta(r, encoding, ::Type{T}; max_elements::Int,
        vertex_id=nothing, positive_gaps=false, source=nothing) where {T<:Unsigned}
    Compression._read_delta_counted(r, encoding, T, max_elements, vertex_id,
        positive_gaps, source)
end

function _cg_read_intervals(r, encoding, mil, ::Type{T}; vertex_id=nothing,
        source=nothing) where {T<:Unsigned}
    count = Int(_cg_read_value(r, encoding, T; source=source)) - 1
    values = T[]; previous = zero(T)
    for i in 1:count
        start = if i == 1 && vertex_id !== nothing
            raw = _cg_read_value(r, encoding, UInt64; source=source)
            T(Int64(vertex_id) + Compression._zigzag_decode(raw - 1))
        else
            previous + _cg_read_value(r, encoding, T; source=source)
        end
        len = Int(_cg_read_value(r, encoding, T; source=source)) - 1 + mil
        for j in 0:len-1; push!(values, start + T(j)); end
        previous = start
    end
    count = Int(_cg_read_value(r, encoding, T; source=source)) - 1
    count > 0 && append!(values, _cg_read_delta(r, encoding, T;
        max_elements=count, vertex_id=vertex_id, source=source))
    sort!(values)
    return values
end

function _decode_clusters_parallel(r, params, clusters, ::Type{T}, ::Type{I},
        resid_bytes, refdist_bytes, copy_bytes, cg_offsets, data_start_bit,
        chunk_offsets, fused_lr, workers, options=(false,false,false)) where {T<:Unsigned,I<:Integer}
    count = length(clusters)
    decoded = Vector{Dict{T,Vector{T}}}(undef, count)
    ends = Vector{Tuple{Int,Int64}}(undef, count)
    next = Threads.Atomic{Int}(1)
    @sync for _ in 1:min(workers, Threads.nthreads(:default), count)
        Threads.@spawn begin
            workspace = (Compression.CtxRangeDecoder(UInt8[]),
                Compression.CtxRangeDecoder(UInt8[]), Compression.BinRangeDecoder(UInt8[]))
            while true
                ci = Threads.atomic_add!(next, 1)
                ci > count && break
                reader = BitReader(r.buffer)
                reader.length = r.length
                decoded[ci] = decode_level(reader, params; T=T, directed=true,
                    coding_scheme=:index, ctx_range=true, resid_bytes=resid_bytes,
                    refdist_bytes=refdist_bytes, copy_bytes=copy_bytes,
                    cg_offsets=cg_offsets, data_start_bit=data_start_bit,
                    chunk_offsets=chunk_offsets, only_cluster=ci,
                    preparsed_clusters=clusters, local_type=I, fused_lr=fused_lr,
                    decoder_workspace=workspace, reuse_identity=options[1],
                    sorted_fastpath=options[2], reuse_scratch=options[3])
                ends[ci] = (reader.index, reader.bit_count)
            end
        end
    end
    # Each source vertex belongs to one cluster; no worker mutates a shared Dict.
    result = Dict{T,Vector{T}}()
    for nls in decoded; merge!(result, nls); end
    # Preserve the serial API's cursor position after the last inter section.
    r.index, r.bit_count = ends[end]
    return result
end
