#
# Compression environment for RL policy learning.
# Wraps the graph compression pipeline as an RL environment where:
#   - State: per-vertex features (degree, interval density, gap size, reference overlap, window density)
#   - Action: encoding choices (integer encoding, reference mode, min interval length)
#   - Reward: negative bits used to encode the vertex (minimize bits = maximize reward)
#   - Episode: one pass through all vertices in the graph
#

# Fixed integer encoding (fibonacci is near-optimal for RCM-ordered web graphs)
const FIXED_INTEGER_ENCODING = :fibonacci
const REFERENCE_OPTIONS = [:none, :reference, :recursive]
const MIN_INTERVAL_OPTIONS = [2, 3, 4, 5]

const NUM_REF_OPT = length(REFERENCE_OPTIONS)       # 3
const NUM_INTV_OPT = length(MIN_INTERVAL_OPTIONS)   # 4
const NUM_ACTIONS = NUM_REF_OPT * NUM_INTV_OPT      # 12

struct Action
    integer_encoding::Symbol
    reference_mode::Symbol       # :none, :reference, :recursive
    min_interval_length::Int
end

function action_to_index(a::Action)::Int
    ref_idx = findfirst(==(a.reference_mode), REFERENCE_OPTIONS)
    mil_idx = findfirst(==(a.min_interval_length), MIN_INTERVAL_OPTIONS)
    return (ref_idx - 1) * NUM_INTV_OPT + mil_idx
end

function action_from_index(idx::Int)::Action
    idx -= 1
    mil_idx = idx % NUM_INTV_OPT + 1
    ref_idx = idx ÷ NUM_INTV_OPT + 1
    return Action(
        FIXED_INTEGER_ENCODING,
        REFERENCE_OPTIONS[ref_idx],
        MIN_INTERVAL_OPTIONS[mil_idx]
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

    if isempty(neighbors)
        # Empty list: just interval header (two 1s)
        return 2 * _est_value_cost(T(1), ie)
    end

    if action.reference_mode != :none && !isempty(env.ref_window)
        # Try reference encoding: find best reference in window
        recursive = action.reference_mode == :recursive
        ref_bits = _try_reference_encoding(env, neighbors, ie, mil, recursive)
        if ref_bits !== nothing
            return ref_bits
        end
    end

    # Interval + residual encoding (estimated cost)
    return _estimate_interval_residual_cost(neighbors, ie, mil)
end

"""
Estimate cost of interval + residual encoding for a neighbor list.
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
Try reference encoding against the window. Returns bits if a good reference found, nothing otherwise.
When recursive=true, applies recursive reference estimation for residuals.
"""
function _try_reference_encoding(env::CompressionEnv{T}, neighbors::Vector{T},
                                 ie::Symbol, mil::Int, recursive::Bool=false)::Union{Int,Nothing} where {T<:Unsigned}
    ns = Set(neighbors)
    best_cost = nothing
    best_ref = nothing

    # Check recent references in window (balance speed vs accuracy)
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

        # Estimate reference cost (use vertex ID, matching deployment writer)
        ref_cost = _est_value_cost(ref_v, ie)
        ref_cost += _estimate_adaptive_bitmap_cost(copy_bitmap, ie)

        if !isempty(residuals)
            ref_cost += _estimate_interval_residual_cost(sort(residuals), ie, mil)
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
    noref_cost = _estimate_interval_residual_cost(neighbors, ie, mil)

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
