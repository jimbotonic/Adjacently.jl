#
# Compression environment for RL policy learning.
# Wraps the graph compression pipeline as an RL environment where:
#   - State: per-vertex features (degree, interval density, gap size, reference overlap, window density)
#   - Action: encoding choices (integer encoding, reference mode, encoding type, min interval length)
#   - Reward: negative bits used to encode the vertex (minimize bits = maximize reward)
#   - Episode: one pass through all vertices in the graph
#

# Fixed integer encoding (fibonacci is near-optimal for RCM-ordered web graphs)
const FIXED_INTEGER_ENCODING = :fibonacci
const REFERENCE_OPTIONS = [:none, :reference, :recursive]
const MIN_INTERVAL_OPTIONS = [2, 3, 4, 5]

# Encoding types: Interval (default), RLE (Run-Length), Delta (Pure Delta)
const ENCODING_TYPES = [:interval, :rle, :delta]

# Flattened Encoding Options (Type, MIL)
# Interval: 4 MILs
# RLE: 4 MILs
# Delta: 1 (MIL ignored)
const ENCODING_OPTIONS = vcat(
    [(:interval, mil) for mil in MIN_INTERVAL_OPTIONS],
    [(:rle, mil) for mil in MIN_INTERVAL_OPTIONS],
    [(:delta, 0)]
)
const NUM_ENCODING_OPT = length(ENCODING_OPTIONS) # 9

const NUM_REF_OPT = length(REFERENCE_OPTIONS)       # 3
const NUM_ACTIONS = NUM_REF_OPT * NUM_ENCODING_OPT  # 27

using LightGraphs: AbstractGraph, vertices, outneighbors
using ..CompressionUtils

struct Action
    integer_encoding::Symbol
    reference_mode::Symbol       # :none, :reference, :recursive
    encoding_type::Symbol        # :interval, :rle, :delta
    min_interval_length::Int
end

function action_to_index(a::Action)::Int
    ref_idx = findfirst(==(a.reference_mode), REFERENCE_OPTIONS)
    
    enc_tuple = (a.encoding_type, a.min_interval_length)
    if a.encoding_type == :delta
        enc_tuple = (:delta, 0)
    end
    enc_idx = findfirst(==(enc_tuple), ENCODING_OPTIONS)
    
    return (ref_idx - 1) * NUM_ENCODING_OPT + enc_idx
end

function action_from_index(idx::Int)::Action
    idx -= 1
    enc_idx = idx % NUM_ENCODING_OPT + 1
    ref_idx = idx ÷ NUM_ENCODING_OPT + 1
    
    enc_type, mil = ENCODING_OPTIONS[enc_idx]
    
    return Action(
        FIXED_INTEGER_ENCODING,
        REFERENCE_OPTIONS[ref_idx],
        enc_type,
        mil
    )
end

mutable struct CompressionEnv{T<:Unsigned}
    neighbor_lists::Dict{T,Vector{T}}
    vertex_order::Vector{T}
    current_idx::Int
    ref_window::Vector{T}
    ref_window_size::Int
    total_bits::Int
    total_edges::Int
    done::Bool
end

function CompressionEnv(neighbor_lists::Dict{T,Vector{T}};
                        ref_window_size::Int=1024) where {T<:Unsigned}
    n = length(neighbor_lists)
    vertex_order = sort(collect(keys(neighbor_lists)))
    total_edges = sum(length(nl) for nl in values(neighbor_lists))
    return CompressionEnv{T}(
        neighbor_lists,
        vertex_order,
        0,
        T[],
        ref_window_size,
        0,
        total_edges,
        true
    )
end

"""
    CompressionEnv(g::AbstractGraph{T}; ref_window_size=1024)

Build an environment from a graph by extracting sorted neighbor lists.
"""
function CompressionEnv(g::AbstractGraph{T}; ref_window_size::Int=1024) where {T<:Unsigned}
    neighbor_lists = Dict{T,Vector{T}}()
    for v in vertices(g)
        neighbor_lists[v] = sort(T.(collect(outneighbors(g, v))))
    end
    return CompressionEnv(neighbor_lists; ref_window_size=ref_window_size)
end

"""
    current_vertex(env) -> Union{T,Nothing}

Return the current vertex being processed, or nothing if episode is done or not started.
"""
function current_vertex(env::CompressionEnv{T})::Union{T,Nothing} where {T<:Unsigned}
    if env.done || env.current_idx < 1 || env.current_idx > length(env.vertex_order)
        return nothing
    end
    return env.vertex_order[env.current_idx]
end

"""
    reset!(env) -> VertexFeatures

Reset the environment to the beginning of the graph. Returns features for the first vertex.
"""
function reset!(env::CompressionEnv{T})::VertexFeatures where {T<:Unsigned}
    env.current_idx = 1
    env.total_bits = 0
    env.done = false
    empty!(env.ref_window)

    v = env.vertex_order[1]
    neighbors = get(env.neighbor_lists, v, T[])
    return extract_features(neighbors, env.ref_window, env.neighbor_lists)
end

"""
    step!(env, action) -> (next_features, reward, done)

Encode the current vertex with the given action, advance to the next vertex.
Returns the next state features, reward (negative bits), and whether the episode is done.
"""
function step!(env::CompressionEnv{T}, action::Action) where {T<:Unsigned}
    if env.done
        dummy = VertexFeatures(1, 1, 1, 1, 1)
        return (dummy, 0.0, true)
    end

    v = env.vertex_order[env.current_idx]
    neighbors = sort(get(env.neighbor_lists, v, T[]))

    # Encode vertex and measure bits
    bits = encode_vertex(env, v, neighbors, action)
    env.total_bits += bits

    # Update reference window (all vertices, matching deployment writer)
    push!(env.ref_window, v)
    if length(env.ref_window) > env.ref_window_size
        popfirst!(env.ref_window)
    end

    # Advance
    env.current_idx += 1
    reward = -Float64(bits)

    if env.current_idx > length(env.vertex_order)
        env.done = true
        dummy = VertexFeatures(1, 1, 1, 1, 1)
        return (dummy, reward, true)
    end

    # Features for next vertex
    next_v = env.vertex_order[env.current_idx]
    next_neighbors = get(env.neighbor_lists, next_v, T[])
    next_features = extract_features(next_neighbors, env.ref_window, env.neighbor_lists)

    return (next_features, reward, false)
end

"""
    encode_vertex(env, vertex, neighbors, action) -> Int

Estimate the bit cost of encoding a vertex with the given action parameters.
Uses the production cost estimation functions from Compression module.
"""
function encode_vertex(env::CompressionEnv{T}, vertex::T, neighbors::Vector{T},
                       action::Action)::Int where {T<:Unsigned}
    ie = action.integer_encoding
    mil = action.min_interval_length
    enc_type = action.encoding_type

    if isempty(neighbors)
        # Empty list: just interval header (two 1s) for Interval/RLE, or 0 for Delta?
        # Actually, write_intervals_and_residuals writes two 1s.
        # write_delta simply returns if empty? No, need to write length.
        # Let's assume standard overhead for empty lists.
        return 2 * _est_value_cost(T(1), ie)
    end

    if action.reference_mode != :none && !isempty(env.ref_window)
        # Try reference encoding: find best reference in window
        recursive = action.reference_mode == :recursive
        ref_bits = _try_reference_encoding(env, neighbors, ie, mil, enc_type, recursive)
        if ref_bits !== nothing
            return ref_bits
        end
    end

    # Base encoding (no reference)
    return _estimate_base_encoding_cost(neighbors, ie, mil, enc_type)
end

"""
Estimate cost of base encoding (Interval, RLE, or Delta) for a neighbor list.
"""
function _estimate_base_encoding_cost(neighbors::Vector{T}, ie::Symbol, mil::Int, enc_type::Symbol)::Int where {T<:Unsigned}
    if enc_type == :interval
        return _estimate_interval_residual_cost(neighbors, ie, mil)
    elseif enc_type == :delta
        return _estimate_delta_cost(neighbors, ie)
    elseif enc_type == :rle
        return _estimate_rle_cost(neighbors, ie, mil)
    else
        # Fallback
        return _estimate_interval_residual_cost(neighbors, ie, mil)
    end
end

"""
Estimate cost of interval + residual encoding.
"""
function _estimate_interval_residual_cost(neighbors::Vector{T}, ie::Symbol,
                                          mil::Int)::Int where {T<:Unsigned}
    intervals, residuals = compress_intervals(neighbors, mil)
    cost = 0

    # Number of intervals (+1 to avoid 0)
    cost += _est_value_cost(T(length(intervals) + 1), ie)

    # Interval data
    prev_start = T(0)
    for (start, len) in intervals
        cost += _est_value_cost(start - prev_start, ie)
        cost += _est_value_cost(T(len - mil + 1), ie)
        prev_start = start
    end

    # Number of residuals (+1 to avoid 0)
    cost += _est_value_cost(T(length(residuals) + 1), ie)

    # Residuals (delta-encoded)
    if !isempty(residuals)
        deltas = delta_encode_vector(residuals)
        cost += _est_value_cost(deltas[1], ie)
        for i in 2:length(deltas)
            cost += _est_value_cost(deltas[i] + T(1), ie)
        end
    end

    return cost
end

"""
Estimate cost of pure delta encoding.
"""
function _estimate_delta_cost(neighbors::Vector{T}, ie::Symbol)::Int where {T<:Unsigned}
    isempty(neighbors) && return 0
    deltas = delta_encode_vector(neighbors)
    cost = 0
    cost += _est_value_cost(deltas[1], ie)
    for i in 2:length(deltas)
        cost += _est_value_cost(deltas[i] + T(1), ie)
    end
    return cost
end

"""
Estimate cost of Run-Length Encoding (via Hybrid Mix).
"""
function _estimate_rle_cost(neighbors::Vector{T}, ie::Symbol, mil::Int)::Int where {T<:Unsigned}
    isempty(neighbors) && return 0
    
    # RLE works on delta list
    deltas = delta_encode_vector(neighbors)
    
    # Hybrid mode cost
    cost = 1 # hybrid flag
    cost += _est_value_cost(deltas[1], ie) # first value
    
    if length(deltas) > 1
        remaining_deltas = deltas[2:end]
        
        # We need original neighbors to detect intervals in hybrid mode
        # neighbors is already the original list.
        # But analyze_delta_patterns_hybrid expects the *shifted* delta list if dealing with residuals?
        # Actually, neighbors are raw values here.
        # analyze_delta_patterns_hybrid expects deltas and original values.
        # It's used for residuals in `write_intervals_and_residuals`.
        # Here we are using it for the *whole* list.
        
        # `write_hybrid_mix_encoded_list` logic:
        # It analyzes `remaining_deltas`.
        
        sections = analyze_delta_patterns_hybrid(remaining_deltas, neighbors[2:end], mil)
        
        cost += _est_value_cost(T(length(sections)), ie)
        
        for section in sections
            if section.type == :delta
                cost += 1 # flag 0
                cost += _est_value_cost(T(length(section.data)), ie)
                for val in section.data
                    cost += _est_value_cost(val, ie)
                end
            elseif section.type == :run_length
                cost += 2 # flag 10
                cost += _est_value_cost(T(length(section.data) ÷ 2), ie)
                for val in section.data
                    cost += _est_value_cost(val, ie)
                end
            elseif section.type == :interval
                cost += 2 # flag 11
                cost += _est_value_cost(T(length(section.data) ÷ 2), ie)
                for val in section.data
                    cost += _est_value_cost(val, ie)
                end
            end
        end
    end
    return cost
end

"""
Try reference encoding against the window. Returns bits if a good reference found, nothing otherwise.
"""
function _try_reference_encoding(env::CompressionEnv{T}, neighbors::Vector{T},
                                 ie::Symbol, mil::Int, enc_type::Symbol, recursive::Bool=false)::Union{Int,Nothing} where {T<:Unsigned}
    ns = Set(neighbors)
    best_cost = nothing
    best_ref = nothing

    # Check recent references in window
    check_count = min(length(env.ref_window), 200)
    start_idx = max(1, length(env.ref_window) - check_count + 1)

    for i in start_idx:length(env.ref_window)
        ref_v = env.ref_window[i]
        ref_nl = get(env.neighbor_lists, ref_v, T[])
        isempty(ref_nl) && continue

        # Compute overlap
        ref_set = Set(ref_nl)
        overlap = length(intersect(ns, ref_set))
        overlap < 3 && continue  # REF_ENCODING_TH

        # Build copy bitmap and residuals
        copy_bitmap = Bool[n in ns for n in ref_nl]
        residuals = T[n for n in neighbors if !(n in ref_set)]

        # Estimate reference cost (use window distance)
        distance = UInt32(length(env.ref_window) - i + 1)
        ref_cost = _est_value_cost(distance, ie)
        ref_cost += _estimate_adaptive_bitmap_cost(copy_bitmap, ie)

        if !isempty(residuals)
            ref_cost += _estimate_base_encoding_cost(sort(residuals), ie, mil, enc_type)
        end

        if best_cost === nothing || ref_cost < best_cost
            best_cost = ref_cost
            best_ref = ref_v
        end
    end

    if best_cost === nothing
        return nothing
    end

    # Compare against non-reference encoding
    noref_cost = _estimate_base_encoding_cost(neighbors, ie, mil, enc_type)

    if best_cost < noref_cost
        return best_cost
    end

    return nothing
end

"""
Estimate the bit cost of encoding a single value using the production cost estimator.
Delegates to Compression.estimate_encoded_value_cost when available.
"""
function _est_value_cost(value::T, ie::Symbol)::Int where {T<:Unsigned}
    return estimate_encoded_value_cost(value, ie)
end

"""
Estimate adaptive bitmap cost (min of raw vs block encoding + 1-bit flag).
Matches the cost model of write_bitmap_adaptive in compression.jl.
"""
function _estimate_adaptive_bitmap_cost(bitmap::Vector{Bool}, ie::Symbol)::Int
    raw_cost = estimate_encoded_value_cost(UInt32(length(bitmap)) + UInt32(1), ie) + length(bitmap)
    block_cost = estimate_block_encoding_cost(bitmap, ie)
    return 1 + min(raw_cost, block_cost)
end

"""
    get_bits_per_edge(env) -> Float64

Return the bits per edge for the current episode.
"""
function get_bits_per_edge(env::CompressionEnv)::Float64
    env.total_edges == 0 && return 0.0
    return Float64(env.total_bits) / Float64(env.total_edges)
end
