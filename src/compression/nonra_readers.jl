# Call-owned entropy sources for the opt-in non-RA LR decoder.
# The hybrid helper is used only for counted (:index) residual-list grammar.
@inline _nonra_flag(r,source)=source===nothing ? read_bit(r) : brc_decode_bit!(source)
function _read_nonra_stop(r,ie,::Type{T},vertex_id,source,flag,buffer=nothing) where T
    result=buffer===nothing ? T[] : empty!(buffer)
    first=true;previous=zero(T)
    while _nonra_flag(r,flag)
        value=if first
            raw=_read_encoded_value(r,ie,UInt64,source)
            T(Int64(vertex_id)+_zigzag_decode(raw-1))
        else
            previous+_read_encoded_value(r,ie,T,source)
        end
        push!(result,value);previous=value;first=false
    end
    result
end
function _read_nonra_hybrid(r::BitReader, coding_scheme::Symbol, integer_encoding::Symbol, ::Type{T}, source, flag_source; max_elements::Union{Int, Nothing} = nothing, stop_value::Union{T, Nothing} = nothing, vertex_id = nothing) where T <: Unsigned
    use_run_length_and_interval = _nonra_flag(r, flag_source)
    raw_first = if vertex_id !== nothing
            _read_encoded_value(r, integer_encoding, UInt64, source)
        else
            _read_encoded_value(r, integer_encoding, T, source)
        end
    if vertex_id !== nothing
        first_value = T(Int64(vertex_id) + _zigzag_decode(raw_first - 1))
    else
        first_value = raw_first
    end
    if stop_value !== nothing && first_value == stop_value
        return T[]
    end
    result = T[first_value]
    elements_read = 1
    stop_consumed = false
    shift = if coding_scheme == :children
            one(T)
        else
            zero(T)
        end
    shift_int = Int(shift)
    last_original = if coding_scheme == :children
            Int(first_value) - shift_int
        else
            Int(first_value)
        end
    if !use_run_length_and_interval
        while true
            if max_elements !== nothing && elements_read >= max_elements
                break
            end
            if stop_value !== nothing
                try
                    if !(peek_bit(r))
                        _consume_children_trailing_stop(r, integer_encoding, T)
                        stop_consumed = true
                        break
                    end
                catch e
                    if e isa EOFError || e isa ErrorException
                        break
                    else
                        rethrow(e)
                    end
                end
            end
            value = _read_encoded_value(r, integer_encoding, T, source)
            if stop_value !== nothing && value == stop_value
                stop_consumed = true
                break
            end
            push!(result, value)
            elements_read += 1
        end
        if stop_value !== nothing && !stop_consumed
            _consume_children_trailing_stop(r, integer_encoding, T)
            stop_consumed = true
        end
    else
        num_sections = _read_encoded_value(r, integer_encoding, T, source)
        for i = 1:num_sections
            section_flag = _nonra_flag(r, flag_source)
            if !section_flag
                count = _read_encoded_value(r, integer_encoding, T, source)
                for j = 1:count
                    if max_elements !== nothing && elements_read >= max_elements
                        break
                    end
                    value = _read_encoded_value(r, integer_encoding, T, source)
                    push!(result, value)
                    elements_read += 1
                    last_original += Int(value) - shift_int
                end
            else
                second_flag = _nonra_flag(r, flag_source)
                if !second_flag
                    num_pairs = _read_encoded_value(r, integer_encoding, T, source)
                    for j = 1:num_pairs
                        value = _read_encoded_value(r, integer_encoding, T, source)
                        length = _read_encoded_value(r, integer_encoding, T, source)
                        for k = 1:length
                            if max_elements !== nothing && elements_read >= max_elements
                                break
                            end
                            push!(result, value)
                            elements_read += 1
                            last_original += Int(value) - shift_int
                        end
                    end
                else
                    num_pairs = _read_encoded_value(r, integer_encoding, T, source)
                    for j = 1:num_pairs
                        start = _read_encoded_value(r, integer_encoding, T, source)
                        length = _read_encoded_value(r, integer_encoding, T, source)
                        if max_elements !== nothing && elements_read >= max_elements
                            break
                        end
                        raw_delta = Int(start) - last_original
                        if raw_delta < 0
                            #= none:129 =# @warn "Interval start before last_original" last_original start elements_read
                            error("Invalid interval ordering: start < previous original")
                        end
                        first_delta = raw_delta + shift_int
                        push!(result, T(first_delta))
                        elements_read += 1
                        last_original = Int(start)
                        for k = 1:Int(length) - 1
                            if max_elements !== nothing && elements_read >= max_elements
                                break
                            end
                            push!(result, T(Int(one(T)) + shift_int))
                            elements_read += 1
                            last_original += 1
                        end
                    end
                end
            end
            if max_elements !== nothing && elements_read >= max_elements
                break
            end
        end
    end
    if stop_value !== nothing && !stop_consumed
        _consume_children_trailing_stop(r, integer_encoding, T)
    end
    if coding_scheme == :children && !(isempty(result))
        result = result .- T(1)
    end
    if length(result) <= 1
        return result
    end
    original_list = T[result[1]]
    for i = 2:length(result)
        push!(original_list, original_list[end] + result[i])
    end
    return original_list
end
