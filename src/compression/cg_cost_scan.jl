
# Exact analytical LR cost without constructing interval/residual vectors.
# Sorted, unique input is the same precondition as the existing estimator.
function _estimate_ir_lr_scan(neighbors::Vector{T}, encoding::Symbol, mil::Int,
        vertex_id; tight_deltas::Bool=false)::Int where {T<:Unsigned}
    if vertex_id === nothing || isempty(neighbors)
        return estimate_interval_runlength_encoding_cost(neighbors,encoding,mil,3;vertex_id)
    end
    vid=T(vertex_id); nc=0; nr=0; nl=0; cost=0
    previous_start=zero(T); previous_left=zero(T); previous_right=zero(T)
    i=1
    @inbounds while i<=length(neighbors)
        len=1
        while i+len<=length(neighbors) && neighbors[i+len]==neighbors[i]+T(len)
            len+=1
        end
        if len>=mil
            start=neighbors[i]
            cost+=estimate_encoded_value_cost(nc==0 ? UInt64(Compression._zigzag_encode(Int64(start)-Int64(vid))+1) : UInt64(start-previous_start),encoding)
            cost+=estimate_encoded_value_cost(T(len-mil+1),encoding)
            previous_start=start;nc+=1
        else
            for j in i:i+len-1
                value=neighbors[j]
                if value<vid
                    if nl>0
                        gap=value-previous_left
                        cost+=estimate_encoded_value_cost(tight_deltas ? gap : gap+T(1),encoding)
                    end
                    previous_left=value;nl+=1
                else
                    d=value-vid+T(1)
                    cost+=estimate_encoded_value_cost(nr==0 ? d : (tight_deltas ? d-previous_right : d-previous_right+T(1)),encoding)
                    previous_right=d;nr+=1
                end
            end
        end
        i+=len
    end
    cost+=estimate_encoded_value_cost(T(nc)+T(1),encoding)
    cost+=estimate_encoded_value_cost(T(nl+nr)+T(1),encoding)
    if nl+nr>0
        cost+=estimate_encoded_value_cost(T(nl)+T(1),encoding)
        nl>0 && (cost+=estimate_encoded_value_cost(vid-previous_left,encoding))
    end
    cost
end
