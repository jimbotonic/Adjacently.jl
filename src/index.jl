#
# Adjacently: Julia Complex Directed Networks Library
# Copyright (C) 2016-2026 Jimmy Dubuisson <jimmy@dubuisson.ch>
#
# This program is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.
#
# This program is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.
#

module Index

export StreamInvertedIndex, OverlapWorkspace,
       ensure_ncols!, add_candidate!,
       overlap!, clear_workspace!, argmax_on_touched, topk_on_touched,
       jaccard_on_touched

# ──────────────────────────────────────────────────────────────────────────────
# Data structures
# ──────────────────────────────────────────────────────────────────────────────

"""
    StreamInvertedIndex{T<:Unsigned}

Incrementally maintain an inverted index (CSC-like) for streaming candidates.
- postings[u] :: Vector{T} - list of candidate IDs (rows) containing child u (column).
- rowdeg[i]  :: T - degree of candidate i.
- k          :: Ref{Int} - number of candidates currently stored.
- ncols      :: Ref{Int} - number of columns currently supported (size of postings).
"""
mutable struct StreamInvertedIndex{T<:Unsigned}
    postings :: Vector{Vector{T}}
    rowdeg   :: Vector{T}
    k        :: Int
    ncols    :: Int
end

"""
    OverlapWorkspace{T}

Reusable buffers for queries to avoid allocations.
- counts  :: Vector{T} - accumulator per candidate.
- touched :: Vector{T} - candidate IDs that received nonzero counts.
"""
mutable struct OverlapWorkspace{T<:Unsigned}
    counts  :: Vector{T}
    touched :: Vector{T}
end

# ──────────────────────────────────────────────────────────────────────────────
# Constructors and capacity management
# ──────────────────────────────────────────────────────────────────────────────

"""
    StreamInvertedIndex{T}(ncols::T; init_candidates::Int=0)

Create an empty streaming index for columns 1..ncols.
`init_candidates` optionally pre-sizes internal arrays for expected candidates.
"""
function StreamInvertedIndex{T}(ncols::Integer; init_candidates::Int=0) where {T<:Unsigned}
    postings = [Vector{T}() for _ in 1:ncols]
    rowdeg = Vector{T}()
    if init_candidates > 0
        sizehint!(rowdeg, init_candidates)
    end
    return StreamInvertedIndex{T}(postings, rowdeg, 0, Int(ncols))
end

"""
    ensure_ncols!(idx::StreamInvertedIndex, needed::T) where {T<:Unsigned}

Grow the postings array if you encounter a child ID > ncols.
Safe to call anytime; no-op if not needed.
"""
function ensure_ncols!(idx::StreamInvertedIndex{T}, needed::Integer) where {T<:Unsigned}
    if needed > idx.ncols
        extra = needed - idx.ncols
        append!(idx.postings, [Vector{T}() for _ in 1:extra])
        idx.ncols = needed
    end
    return idx
end

"""
    OverlapWorkspace(idx::StreamInvertedIndex; capacity_candidates::Int=0)

Create a workspace sized for the current number of candidates.
Optionally pre-reserve more space.
"""
function OverlapWorkspace(idx::StreamInvertedIndex{T}; capacity_candidates::Int=0) where {T<:Unsigned}
    k = max(idx.k, capacity_candidates)
    counts  = zeros(T, k > 0 ? k : 1)
    touched = Vector{T}()
    return OverlapWorkspace(counts, touched)
end

# ──────────────────────────────────────────────────────────────────────────────
# Adding candidates (streaming)
# ──────────────────────────────────────────────────────────────────────────────

"""
    add_candidate!(idx, children::AbstractVector{T}) where {T<:Unsigned} -> T

Add a candidate row with given children IDs (1-based). Returns the new candidate ID.
Per-candidate duplicates are removed (optional, but recommended).
"""
function add_candidate!(idx::StreamInvertedIndex{T},
                        children::AbstractVector) where {T<:Unsigned}
    # Ensure enough columns if streaming brought a new max child ID
    if !isempty(children)
        maxu = maximum(children)
        ensure_ncols!(idx, maxu)
    end

    # Deduplicate the row once to avoid double counting in queries
    tmp = T.(children)
    sort!(tmp)
    unique!(tmp)

    # Append row
    new_id = idx.k + 1
    push!(idx.rowdeg, T(length(tmp)))
    # Update postings
    @inbounds for u in tmp
        push!(idx.postings[u], T(new_id))
    end
    idx.k = new_id
    return T(new_id)
end

# Graphs.jl dependency removed - users can manually extract neighbors and use add_candidate!

# ──────────────────────────────────────────────────────────────────────────────
# Query: overlap counts y = C * v (SpMSpV via postings)
# ──────────────────────────────────────────────────────────────────────────────

"""
    clear_workspace!(work::OverlapWorkspace)

Zero out only the touched entries in `counts` and clear touched list.
"""
function clear_workspace!(work::OverlapWorkspace{T}) where {T<:Unsigned}
    @inbounds for i in work.touched
        work.counts[i] = 0
    end
    empty!(work.touched)
    return work
end

"""
    overlap!(idx::StreamInvertedIndex{T},
             target_children::AbstractVector{T},
             work::OverlapWorkspace)
        -> (counts::Vector{T}, touched::Vector{T})

Compute counts per candidate using the inverted index (postings).
Only increments candidates that share at least one child with the target.
"""
function overlap!(idx::StreamInvertedIndex{T},
                  target_children::AbstractVector{T},
                  work::OverlapWorkspace{T}) where {T<:Unsigned}
    # Make sure workspace can hold all candidates
    if length(work.counts) < idx.k
        resize!(work.counts, idx.k)
    end
    clear_workspace!(work)

    @inbounds for uT in target_children
        u = T(uT)
        if 1 <= u <= idx.ncols
            rows = idx.postings[u]
            for r in rows
                c = work.counts[r] + 1
                work.counts[r] = c
                if c == 1
                    push!(work.touched, r)
                end
            end
        end
    end
    return work.counts, work.touched
end

# ──────────────────────────────────────────────────────────────────────────────
# Utilities: argmax / top-k / Jaccard on touched
# ──────────────────────────────────────────────────────────────────────────────

"""
    argmax_on_touched(counts, touched) -> (i_max, val_max)

Argmax restricted to candidates that were touched. Returns (0, 0) if none.
"""
function argmax_on_touched(counts::AbstractVector{T},
                           touched::AbstractVector{T}) where {T<:Unsigned}
    if isempty(touched)
        return T(0), T(0)
    end
    best_i = touched[1]
    best_v = counts[best_i]
    @inbounds for j in 2:length(touched)
        i = touched[j]
        v = counts[i]
        if v > best_v
            best_v = v
            best_i = i
        end
    end
    return best_i, best_v
end

"""
    topk_on_touched(counts, touched, K::T) -> Vector{T}

Return the candidate IDs of the top-K counts among touched (descending).
"""
function topk_on_touched(counts::AbstractVector{T},
                         touched::AbstractVector{T}, K::Integer) where {T<:Unsigned}
    k_actual = min(K, length(touched))
    if k_actual <= 0
        return Tuple{T,T}[]
    end
    # partial sort by value, descending
    perm = sortperm(touched; by = i -> counts[i], rev = true)
    result = Tuple{T,T}[]
    for i in 1:k_actual
        idx = touched[perm[i]]
        push!(result, (idx, counts[idx]))
    end
    return result
end

"""
    jaccard_on_touched(counts, touched, rowdeg, d_t) -> Dict{T,Float32}

Compute Jaccard(i,t) for touched candidates only.
"""
function jaccard_on_touched(counts::AbstractVector{T},
                            touched::AbstractVector{T},
                            rowdeg::AbstractVector{T},
                            d_t::T) where {T<:Unsigned}
    scores = Dict{T,Float32}()
    @inbounds for i in touched
        inter = counts[i]
        uni = rowdeg[i] + d_t - inter
        scores[i] = inter == 0 ? 0f0 : Float32(inter) / Float32(uni)
    end
    return scores
end

end # module
