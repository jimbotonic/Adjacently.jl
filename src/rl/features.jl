#
# Feature extraction for RL-based compression policy learning.
# Converts per-vertex graph data into discretized feature vectors
# suitable for tabular Q-learning.
#

# Degree bins: [0, 1-3, 4-10, 11-50, 51-200, 201+]
const DEGREE_BINS = 6
# Interval density bins: [0-25%, 25-50%, 50-75%, 75-100%]
const INTERVAL_DENSITY_BINS = 4
# Max gap bins: [1-2, 3-10, 11-100, 101-1000, 1001+]
const MAX_GAP_BINS = 5
# Reference overlap bins: [0%, 1-25%, 25-50%, 50-75%, 75-100%]
const REF_OVERLAP_BINS = 5
# Reference window density bins: [0-25%, 25-50%, 50-75%, 75-100%]
const REF_WINDOW_DENSITY_BINS = 4

const NUM_STATES = DEGREE_BINS * INTERVAL_DENSITY_BINS * MAX_GAP_BINS * REF_OVERLAP_BINS * REF_WINDOW_DENSITY_BINS  # 2400

struct VertexFeatures
    degree_bin::Int              # 1-6
    interval_density::Int        # 1-4
    max_gap_bin::Int             # 1-5
    ref_overlap_bin::Int         # 1-5
    ref_window_density_bin::Int  # 1-4
end

function degree_to_bin(degree::Int)::Int
    degree == 0 && return 1
    degree <= 3 && return 2
    degree <= 10 && return 3
    degree <= 50 && return 4
    degree <= 200 && return 5
    return 6
end

function gap_to_bin(max_gap::Int)::Int
    max_gap <= 2 && return 1
    max_gap <= 10 && return 2
    max_gap <= 100 && return 3
    max_gap <= 1000 && return 4
    return 5
end

function density_to_bin(density::Float64)::Int
    density <= 0.25 && return 1
    density <= 0.50 && return 2
    density <= 0.75 && return 3
    return 4
end

function overlap_to_bin(overlap_ratio::Float64)::Int
    overlap_ratio <= 0.0 && return 1
    overlap_ratio <= 0.25 && return 2
    overlap_ratio <= 0.50 && return 3
    overlap_ratio <= 0.75 && return 4
    return 5
end

"""
    feature_index(f::VertexFeatures) -> Int

Map discretized features to a flat Q-table index (1-based).
"""
function feature_index(f::VertexFeatures)::Int
    return (f.degree_bin - 1) * (INTERVAL_DENSITY_BINS * MAX_GAP_BINS * REF_OVERLAP_BINS * REF_WINDOW_DENSITY_BINS) +
           (f.interval_density - 1) * (MAX_GAP_BINS * REF_OVERLAP_BINS * REF_WINDOW_DENSITY_BINS) +
           (f.max_gap_bin - 1) * (REF_OVERLAP_BINS * REF_WINDOW_DENSITY_BINS) +
           (f.ref_overlap_bin - 1) * REF_WINDOW_DENSITY_BINS +
           f.ref_window_density_bin
end

"""
    extract_features(neighbors, ref_candidates, neighbor_lists) -> VertexFeatures

Extract discretized features for a vertex given its sorted neighbor list,
available reference candidates, and the full neighbor list dictionary.

The ref_window_density is the fraction of ref_candidates that have degree >= 4
(the REF_V_MIN_DEGREE threshold), indicating how rich the reference window is.
"""
function extract_features(neighbors::Vector{T}, ref_candidates::Vector{T},
                          neighbor_lists::Dict{T,Vector{T}}) where {T<:Unsigned}
    degree = length(neighbors)

    # Degree bin
    d_bin = degree_to_bin(degree)

    # Interval density: fraction of neighbors in consecutive intervals
    if degree <= 1
        iv_density = 0.0
        max_gap = 1
    else
        intervals, residuals = compress_intervals(neighbors, 2)  # Use min_length=2 for feature extraction
        interval_count = sum(len for (_, len) in intervals; init=0)
        iv_density = interval_count / degree
        # Max gap from deltas
        deltas = delta_encode_vector(neighbors)
        max_gap = length(deltas) > 1 ? Int(maximum(deltas[2:end])) : 1
    end

    iv_bin = density_to_bin(iv_density)
    gap_bin = gap_to_bin(max_gap)

    # Best reference overlap
    best_overlap = 0.0
    if degree >= 2 && !isempty(ref_candidates)
        ns = Set(neighbors)
        # Check up to 20 candidates for speed
        check_limit = min(length(ref_candidates), 20)
        for i in 1:check_limit
            ref_v = ref_candidates[i]
            ref_nl = get(neighbor_lists, ref_v, T[])
            if !isempty(ref_nl)
                overlap = length(intersect(ns, Set(ref_nl)))
                ratio = overlap / max(degree, length(ref_nl))
                if ratio > best_overlap
                    best_overlap = ratio
                end
            end
        end
    end

    ov_bin = overlap_to_bin(best_overlap)

    # Reference window density: fraction of candidates with degree >= 4
    if isempty(ref_candidates)
        rw_density = 0.0
    else
        high_degree_count = count(ref_candidates) do ref_v
            length(get(neighbor_lists, ref_v, T[])) >= 4
        end
        rw_density = high_degree_count / length(ref_candidates)
    end
    rw_bin = density_to_bin(rw_density)

    return VertexFeatures(d_bin, iv_bin, gap_bin, ov_bin, rw_bin)
end
