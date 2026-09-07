# Explicit-source primitives for opt-in CG decoding. Legacy BG/CS readers are unchanged.
@inline function _rankgap_inverse!(Q::Vector{I}, C::AbstractVector{I}) where {I<:Integer}
    ci = 1; nc = length(C)
    @inbounds for k in eachindex(Q)
        r = Q[k] + I(ci - 1)
        while ci <= nc && C[ci] <= r
            ci += 1; r += one(I)
        end
        Q[k] = r
    end
    return Q
end

function _merge_sorted_owned!(a::Vector{T}, b::Vector{T}) where T
    a === b && throw(ArgumentError("merge inputs must not alias"))
    isempty(b) && return a
    na = length(a); nb = length(b)
    resize!(a, na + nb)
    i = na; j = nb; k = na + nb
    @inbounds while i > 0 && j > 0
        if a[i] > b[j]
            a[k] = a[i]; i -= 1
        else
            a[k] = b[j]; j -= 1
        end
        k -= 1
    end
    @inbounds while j > 0
        a[k] = b[j]; k -= 1; j -= 1
    end
    return a
end

@inline _read_encoded_value(r::BitReader, compression::Symbol, ::Type{T}, source::CtxRangeDecoder) where {T<:Unsigned} =
    T(rc_decode_value!(source))

function _read_encoded_value(r::BitReader, compression::Symbol, ::Type{T}, ::Nothing) where {T<:Unsigned}
    if compression == :elias_gamma
        return read_elias_gamma(r, T)
    elseif compression == :elias_delta
        return read_elias_delta(r, T)
    elseif compression == :golomb
        return read_golomb(r, GOLOMB_BASE, T)
    elseif compression == :fibonacci
        return read_fibonacci(r, T)
    elseif compression == :zeta
        return read_zeta(r, ZETA_BASE, T)
    elseif compression == :fed
        return read_fed(r, T, FED_BLOCK_SIZE)
    else
        throw(ArgumentError("Invalid compression code: $compression"))
    end
end

function _read_delta_counted(r, encoding, ::Type{T}, count::Int, vertex_id,
                             positive_gaps::Bool, source::S) where {T<:Unsigned,S}
    count >= 0 || throw(ArgumentError("delta count must be nonnegative"))
    lst = T[]
    count == 0 && return lst
    sizehint!(lst, min(count, 1024))
    try
        first_value = if vertex_id === nothing
            _read_encoded_value(r, encoding, T, source)
        else
            raw = _read_encoded_value(r, encoding, UInt64, source)
            T(Int64(vertex_id) + _zigzag_decode(raw - 1))
        end
        push!(lst, first_value)
        shift = positive_gaps ? zero(T) : one(T)
        previous = first_value
        for _ in 2:count
            gap = _read_encoded_value(r, encoding, T, source)
            previous = previous + gap - shift
            push!(lst, previous)
        end
    catch e
        (e isa EOFError || e isa ErrorException) || rethrow(e)
    end
    return lst
end

function _read_intervals_lr_source(r, encoding, mil, ::Type{T}, vid::T,
                                  tight_intervals, source::S, buffer=nothing, merge_scratch=nothing) where {T<:Unsigned,S}

    # Read intervals (same format as standard)
    num_intervals = Int(_read_encoded_value(r, encoding, T, source)) - 1
    neighbors = buffer === nothing ? T[] : empty!(buffer::Vector{T})

    if num_intervals > 0
        prev_ref = T(0)
        for idx in 1:num_intervals
            if idx == 1
                raw_start = _read_encoded_value(r, encoding, UInt64, source)
                start = T(Int64(vid) + _zigzag_decode(raw_start - 1))
            else
                start = prev_ref + _read_encoded_value(r, encoding, T, source)
            end
            len = Int(_read_encoded_value(r, encoding, T, source)) - 1 + mil
            for j in 0:(len-1)
                push!(neighbors, start + T(j))
            end
            prev_ref = tight_intervals ? (start + T(len)) : start
        end
    end

    # Intervals are the first sorted run; left/right residuals form the second.
    interval_count = length(neighbors)
    num_residuals = Int(_read_encoded_value(r, encoding, T, source)) - 1
    if num_residuals > 0
        n_left = Int(_read_encoded_value(r, encoding, T, source)) - 1
        n_right = num_residuals - n_left

        if source isa CtxRangeDecoder
            # Reserve both halves together; separate size hints can shrink or
            # grow the same backing allocation twice on every vertex.
            sizehint!(neighbors, interval_count + min(num_residuals, 1024); shrink=false)
            if n_left > 0
                _append_lr_range_deltas!(neighbors, source, n_left, vid, true)
                reverse!(neighbors, interval_count+1, length(neighbors))
            end
            n_right > 0 && _append_lr_range_deltas!(neighbors, source, n_right, vid, false)
        elseif n_left > 0
            left_dists = _read_delta_counted(r, encoding, T, n_left, nothing, false, source)
            for i in n_left:-1:1
                push!(neighbors, vid - left_dists[i])
            end
        end

        if !(source isa CtxRangeDecoder) && n_right > 0
            right_dists = _read_delta_counted(r, encoding, T, n_right, nothing, false, source)
            for d in right_dists
                push!(neighbors, vid + d - T(1))
            end
        end
    end

    if source isa CtxRangeDecoder
        _merge_lr_runs!(neighbors, interval_count, merge_scratch)
    else
        sort!(neighbors) # retain legacy non-range behavior
    end
    return neighbors
end

function _merge_lr_runs!(a::Vector{T}, split::Int, scratch=nothing) where T
    n=length(a)
    0 <= split <= n || throw(ArgumentError("invalid sorted-run split"))
    (split == 0 || split == n || a[split] <= a[split+1]) && return a
    # Copy only the shorter run. Forward/backward merging cannot overwrite the
    # still-unread portion of the other run. Duplicates are preserved.
    if split <= n-split
        saved = if scratch === nothing
            a[1:split]
        else
            scratch === a && throw(ArgumentError("merge scratch must not alias output"))
            resize!(scratch, split)
            copyto!(scratch, 1, a, 1, split)
        end
        i=1; j=split+1; out=1
        @inbounds while i <= split && j <= n
            if saved[i] <= a[j]
                a[out]=saved[i]; i+=1
            else
                a[out]=a[j]; j+=1
            end
            out+=1
        end
        @inbounds while i <= split; a[out]=saved[i]; i+=1; out+=1; end
    else
        saved = if scratch === nothing
            a[split+1:n]
        else
            scratch === a && throw(ArgumentError("merge scratch must not alias output"))
            resize!(scratch, n-split)
            copyto!(scratch, 1, a, split+1, n-split)
        end
        i=split; j=n-split; out=n
        @inbounds while i >= 1 && j >= 1
            if a[i] > saved[j]
                a[out]=a[i]; i-=1
            else
                a[out]=saved[j]; j-=1
            end
            out-=1
        end
        @inbounds while j >= 1; a[out]=saved[j]; j-=1; out-=1; end
    end
    return a
end

function _append_lr_range_deltas!(neighbors::Vector{T}, source::CtxRangeDecoder,
                                  count::Int, vid::T, left::Bool) where T
    distance = T(rc_decode_value!(source))
    push!(neighbors, left ? vid-distance : vid+distance-one(T))
    for _ in 2:count
        distance += T(rc_decode_value!(source)) - one(T)
        push!(neighbors, left ? vid-distance : vid+distance-one(T))
    end
    return neighbors
end
