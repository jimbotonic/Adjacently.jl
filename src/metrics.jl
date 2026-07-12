#
# Adjacently: Julia Complex Directed Networks Library
# Copyright (C) 2016-2026 Anonymous (double-blind review)
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

"""
    Metrics

Unified module for graph analysis and ordering quality metrics.

Re-exports metric functions from Graph, Relabeling, and Distribution modules,
and adds new distributional metrics for compression analysis.

## Existing metrics (re-exported):

From `Graph`: `get_basic_stats`, `display_basic_stats`, `get_out_degree_stats`,
`get_sinks`, `get_sources`, `get_in_degrees`, `get_out_degrees`, `get_in_out_degrees`,
`get_clustering_coefficients`, `get_colink_coefficients`, `compute_volume`,
`compute_cut`, `compute_conductance`, `compute_modularity_gain`, `compute_support_score`

From `Relabeling`: `ordering_quality_metrics`, `print_ordering_metrics`,
`compare_ordering_metrics`

From `Distribution`: `get_graph_entropy`, `get_degree_entropy`, `get_entropy`

## New metrics:

- `compute_consecutive_metrics`: consecutive vertex similarity and best-ref distance
- `compute_gap_distribution_metrics`: successor/residual gap distributions
"""
module Metrics

using LightGraphs: AbstractGraph, nv, ne, outneighbors, vertices, eltype

# ── Re-exports from existing modules ────────────────────────────────────────

using ..Graph: get_basic_stats, display_basic_stats, get_out_degree_stats,
    get_sinks, get_sources, get_in_degrees, get_out_degrees, get_in_out_degrees,
    get_clustering_coefficients, get_colink_coefficients,
    compute_volume, compute_cut, compute_conductance,
    compute_modularity_gain, compute_support_score

using ..Relabeling: ordering_quality_metrics, print_ordering_metrics,
    compare_ordering_metrics

using ..Distribution: get_graph_entropy, get_degree_entropy, get_entropy

export # Re-exported from Graph
       get_basic_stats, display_basic_stats, get_out_degree_stats,
       get_sinks, get_sources, get_in_degrees, get_out_degrees, get_in_out_degrees,
       get_clustering_coefficients, get_colink_coefficients,
       compute_volume, compute_cut, compute_conductance,
       compute_modularity_gain, compute_support_score,
       # Re-exported from Relabeling
       ordering_quality_metrics, print_ordering_metrics, compare_ordering_metrics,
       # Re-exported from Distribution
       get_graph_entropy, get_degree_entropy, get_entropy,
       # New metrics
       compute_consecutive_metrics, compute_gap_distribution_metrics,
       compressibility_score

# ── New metric functions ────────────────────────────────────────────────────

"""
    compute_consecutive_metrics(g; window=64)

Compute consecutive vertex similarity and best-reference distance metrics.
These metrics quantify the "hierarchical locality" of a vertex ordering —
how well consecutive vertices share neighbors (enabling reference copy)
and how close the best reference is.

Returns a `Dict{Symbol,Any}` with keys:
- `:avg_consec_jaccard` — average Jaccard(v, v-1) across all consecutive pairs
- `:avg_consec_overlap` — average |N(v) ∩ N(v-1)| (absolute overlap count)
- `:pct_best_ref_d1` — fraction of vertices whose best reference is at distance 1
- `:avg_best_ref_dist` — average distance to the best reference within the window

High `avg_consec_jaccard` (>0.5) indicates web-crawl-quality ordering.
High `pct_best_ref_d1` (>50%) indicates strong sequential reference locality.
"""
function compute_consecutive_metrics(g::AbstractGraph; window::Int=64)
    n = Int(nv(g))

    # Build adjacency sets
    adj = Vector{Set{Int}}(undef, n)
    for v in vertices(g)
        adj[Int(v)] = Set{Int}(Int.(outneighbors(g, v)))
    end

    sum_jaccard = 0.0
    sum_overlap = 0.0
    n_consec = 0
    n_best_d1 = 0
    sum_best_dist = 0.0
    n_with_ref = 0

    for v in 2:n
        nv_set = adj[v]
        isempty(nv_set) && continue

        # Consecutive Jaccard with v-1
        nv1_set = adj[v-1]
        if !isempty(nv1_set)
            isect = length(intersect(nv_set, nv1_set))
            union_size = length(union(nv_set, nv1_set))
            sum_jaccard += union_size > 0 ? isect / union_size : 0.0
            sum_overlap += isect
            n_consec += 1
        end

        # Best reference in window: which distance d has max overlap?
        best_d = 0
        best_ov = 0
        for d in 1:min(window, v - 1)
            ref_set = adj[v - d]
            isempty(ref_set) && continue
            ov = length(intersect(nv_set, ref_set))
            if ov > best_ov
                best_ov = ov
                best_d = d
            end
        end
        if best_d > 0
            n_with_ref += 1
            sum_best_dist += best_d
            if best_d == 1
                n_best_d1 += 1
            end
        end
    end

    return Dict{Symbol,Any}(
        :avg_consec_jaccard => n_consec > 0 ? sum_jaccard / n_consec : 0.0,
        :avg_consec_overlap => n_consec > 0 ? sum_overlap / n_consec : 0.0,
        :pct_best_ref_d1    => n_with_ref > 0 ? n_best_d1 / n_with_ref : 0.0,
        :avg_best_ref_dist  => n_with_ref > 0 ? sum_best_dist / n_with_ref : 0.0,
    )
end

"""
    compute_gap_distribution_metrics(g; window=64)

Compute distributional metrics for successor gaps and residual gaps after
best-reference copy. These metrics reveal the "second-order locality" that
determines compression quality — not just whether neighbors are close, but
whether the *residuals after reference copy* are also close.

Returns a `Dict{Symbol,Any}` with keys:

Successor gap statistics:
- `:gap1_frac` — fraction of consecutive-neighbor gaps equal to 1
- `:gap_le4_frac` — fraction of gaps ≤ 4
- `:gap_le16_frac` — fraction of gaps ≤ 16
- `:gap_entropy` — entropy of gap distribution (bits)
- `:gap_mean` — mean gap value
- `:gap_p50` — median gap
- `:gap_p90` — 90th percentile gap

Structure metrics:
- `:deg_autocorrelation` — autocorrelation of the degree sequence (high = similar-degree vertices are grouped)

Residual metrics (after best-reference copy):
- `:residuals_per_vertex` — average residuals per non-empty vertex after best-ref copy
- `:residual_gap1_frac` — fraction of residual gaps equal to 1 (key discriminator!)
- `:residual_gap_mean` — mean residual gap

The `:residual_gap1_frac` metric is the key discriminator between crawl-ordered
web graphs (~66%) and Leiden+LLP reordered graphs (~46%). High values indicate
that not only are neighbors close to each other, but the neighbors that *can't*
be copied from the reference are *also* close to each other.
"""
function compute_gap_distribution_metrics(g::AbstractGraph; window::Int=64)
    n = Int(nv(g))

    # Build sorted neighbor lists
    adj = Vector{Vector{Int}}(undef, n)
    for v in vertices(g)
        adj[Int(v)] = sort(Int.(outneighbors(g, v)))
    end

    # ── Successor gap distribution ──
    all_gaps = Int[]
    total_edges = sum(length, adj)
    sizehint!(all_gaps, total_edges)
    for vi in 1:n
        nl = adj[vi]
        for j in 2:length(nl)
            push!(all_gaps, nl[j] - nl[j-1])
        end
    end

    n_gaps = length(all_gaps)
    if n_gaps == 0
        return Dict{Symbol,Any}(
            :gap1_frac => 0.0, :gap_le4_frac => 0.0, :gap_le16_frac => 0.0,
            :gap_entropy => 0.0, :gap_mean => 0.0, :gap_p50 => 0, :gap_p90 => 0,
            :deg_autocorrelation => 0.0,
            :residuals_per_vertex => 0.0, :residual_gap1_frac => 0.0, :residual_gap_mean => 0.0)
    end

    gap1_frac = count(==(1), all_gaps) / n_gaps
    gap_le4_frac = count(<=(4), all_gaps) / n_gaps
    gap_le16_frac = count(<=(16), all_gaps) / n_gaps
    gap_mean_val = sum(all_gaps) / n_gaps

    sort!(all_gaps)
    gap_p50 = all_gaps[max(1, round(Int, 0.50 * n_gaps))]
    gap_p90 = all_gaps[max(1, round(Int, 0.90 * n_gaps))]

    # Gap entropy (bits)
    gap_max_bin = min(all_gaps[end], 10000)
    gap_hist = zeros(Int, gap_max_bin + 1)
    for g_val in all_gaps
        idx = min(g_val, gap_max_bin) + 1
        gap_hist[idx] += 1
    end
    gap_entropy = 0.0
    for c in gap_hist
        c == 0 && continue
        p = c / n_gaps
        gap_entropy -= p * log2(p)
    end

    # ── Degree autocorrelation ──
    degrees = [length(adj[v]) for v in 1:n]
    deg_mean_val = sum(degrees) / n
    deg_var = sum((d - deg_mean_val)^2 for d in degrees) / n
    deg_autocorr = 0.0
    if deg_var > 0
        for v in 2:n
            deg_autocorr += (degrees[v] - deg_mean_val) * (degrees[v-1] - deg_mean_val)
        end
        deg_autocorr /= ((n - 1) * deg_var)
    end

    # ── Residual analysis after best-ref copy ──
    n_residuals_total = 0
    residual_gaps = Int[]
    n_nonempty = count(v -> !isempty(adj[v]), 1:n)

    for vi in 2:n
        nl = adj[vi]
        isempty(nl) && continue
        nl_set = Set(nl)

        # Find best ref in window
        best_ref = 0
        best_ov = 0
        for d in 1:min(window, vi - 1)
            ref_nl = adj[vi - d]
            isempty(ref_nl) && continue
            ov = length(intersect(nl_set, Set(ref_nl)))
            if ov > best_ov
                best_ov = ov
                best_ref = vi - d
            end
        end

        if best_ref > 0 && best_ov > 0
            ref_set = Set(adj[best_ref])
            residuals = sort(collect(setdiff(nl_set, ref_set)))
            n_residuals_total += length(residuals)
            for j in 2:length(residuals)
                push!(residual_gaps, residuals[j] - residuals[j-1])
            end
        else
            n_residuals_total += length(nl)
        end
    end

    avg_residuals = n_nonempty > 0 ? n_residuals_total / n_nonempty : 0.0
    res_gap_mean = isempty(residual_gaps) ? 0.0 : sum(residual_gaps) / length(residual_gaps)
    res_gap1_frac = isempty(residual_gaps) ? 0.0 : count(==(1), residual_gaps) / length(residual_gaps)

    return Dict{Symbol,Any}(
        :gap1_frac => gap1_frac,
        :gap_le4_frac => gap_le4_frac,
        :gap_le16_frac => gap_le16_frac,
        :gap_entropy => gap_entropy,
        :gap_mean => gap_mean_val,
        :gap_p50 => gap_p50,
        :gap_p90 => gap_p90,
        :deg_autocorrelation => deg_autocorr,
        :residuals_per_vertex => avg_residuals,
        :residual_gap1_frac => res_gap1_frac,
        :residual_gap_mean => res_gap_mean,
    )
end

"""Inverse sigmoid mapping: high input → low score, with midpoint and steepness control."""
_inv_sigmoid(x::Float64, midpoint::Float64, steepness::Float64) =
    clamp(1.0 / (1.0 + exp(steepness * (x - midpoint))), 0.0, 1.0)

"""
    compressibility_score(g; window=64)

Compute a compressibility score in [0, 1] that predicts how well a directed
graph will compress with reference-based encoders (BV, BG, CS, CG).

The score is derived from empirical correlation analysis (17 data points
across CNR-2000, LFR, Web, and Erdős-Rényi graphs with multiple orderings):

| Component          | Weight | Pearson r with BPE | Normalization                    |
|--------------------|--------|--------------------|---------------------------------|
| gap_entropy        | 0.30   | +0.990             | Sigmoid, midpoint=6.0, k=0.6    |
| avg_copy_fraction  | 0.20   | -0.982             | Direct [0, 1]                   |
| avg_log_gap        | 0.15   | +0.981             | Sigmoid, midpoint=5.0, k=0.7    |
| residual_gap1_frac | 0.20   | -0.839             | Direct [0, 1] (boosted weight)  |
| residuals/vertex   | 0.15   | +0.892             | Normalized by avg degree         |

The residual_gap1 weight is boosted (0.20 vs correlation-proportional 0.10) because
it is the key discriminator between orderings with similar gap metrics — e.g., CNR-2000
crawl order (res_gap1=0.66) vs Leiden+LLP (res_gap1=0.46) have nearly identical gap
entropy but 8% different BPE.

Calibration reference points (known BPE → score):
- CNR-2000 crawl order:   BPE ≈ 2.5  → score ≈ 0.92
- CNR-2000 Leiden+LLP:    BPE ≈ 2.7  → score ≈ 0.89
- LFR μ=0.05 Leiden+LLP:  BPE ≈ 8.0  → score ≈ 0.40
- Web deg=16 original:     BPE ≈ 9.0  → score ≈ 0.37
- LFR μ=0.50 Leiden+LLP:  BPE ≈ 11.3 → score ≈ 0.15
- ER p=0.005:             BPE ≈ 15.0 → score ≈ 0.09
- LFR random ordering:    BPE ≈ 14.0 → score ≈ 0.05

Score interpretation:
- **0.85–1.0**: Excellent — web crawl quality, BPE 2–4
- **0.55–0.85**: Good — strong community structure with good ordering, BPE 4–8
- **0.30–0.55**: Moderate — weak communities or moderate locality, BPE 8–11
- **0.10–0.30**: Poor — sparse locality, random-like gaps, BPE 11–14
- **0.00–0.10**: Minimal — no exploitable structure, BPE 14+

Returns a named tuple `(score, components)` where `components` is a Dict
with the raw metric values and individual normalized sub-scores.

# Example
```julia
s, c = compressibility_score(g)
println("Compressibility: \$(round(s.score, digits=3))")
println("Gap entropy: \$(round(c[:gap_entropy_raw], digits=2)) bits")
println("Copy fraction: \$(round(100*c[:copy_fraction_raw], digits=1))%")
```
"""
function compressibility_score(g::AbstractGraph; window::Int=64)
    n = Int(nv(g))
    m = Int(ne(g))

    m == 0 && return (score=0.0, components=Dict{Symbol,Float64}())

    # ── Compute needed metrics ──────────────────────────────────────────

    # Build sorted neighbor lists
    adj = Vector{Vector{Int}}(undef, n)
    for v in vertices(g)
        adj[Int(v)] = sort(Int.(outneighbors(g, v)))
    end

    # 1. Gap entropy and avg_log_gap
    all_gaps = Int[]
    sizehint!(all_gaps, m)
    for vi in 1:n
        nl = adj[vi]
        for j in 2:length(nl)
            push!(all_gaps, nl[j] - nl[j-1])
        end
    end
    n_gaps = length(all_gaps)

    gap_entropy = 0.0
    avg_log_gap = 0.0
    if n_gaps > 0
        gap_max_bin = min(maximum(all_gaps), 10000)
        gap_hist = zeros(Int, gap_max_bin + 1)
        for g_val in all_gaps
            gap_hist[min(g_val, gap_max_bin) + 1] += 1
            avg_log_gap += log2(max(1, g_val))
        end
        avg_log_gap /= n_gaps
        for c in gap_hist
            c == 0 && continue
            p = c / n_gaps
            gap_entropy -= p * log2(p)
        end
    end

    # 2. Reference overlap, copy fraction, residual stats (sampled)
    n_sample = min(n, 5000)
    sample_step = max(1, n ÷ n_sample)
    sum_copy_frac = 0.0
    n_with_nbrs = 0
    total_residuals = 0
    n_res_gap1 = 0
    n_res_gaps = 0

    for vi in 2:sample_step:n
        nl = adj[vi]
        isempty(nl) && continue
        nl_set = Set(nl)
        n_with_nbrs += 1

        best_ov = 0; best_ref = 0
        for d in 1:min(window, vi - 1)
            ref_nl = adj[vi - d]
            isempty(ref_nl) && continue
            ov = length(intersect(nl_set, Set(ref_nl)))
            if ov > best_ov; best_ov = ov; best_ref = vi - d; end
        end

        if best_ref > 0 && best_ov > 0
            sum_copy_frac += best_ov / length(nl)
            residuals = sort(collect(setdiff(nl_set, Set(adj[best_ref]))))
            total_residuals += length(residuals)
            for j in 2:length(residuals)
                n_res_gaps += 1
                if residuals[j] - residuals[j-1] == 1
                    n_res_gap1 += 1
                end
            end
        else
            total_residuals += length(nl)
        end
    end

    avg_copy_fraction = n_with_nbrs > 0 ? sum_copy_frac / n_with_nbrs : 0.0
    avg_residuals_per_v = n_with_nbrs > 0 ? total_residuals / n_with_nbrs : 0.0
    res_gap1_frac = n_res_gaps > 0 ? n_res_gap1 / n_res_gaps : 0.0
    avg_deg = m / n

    # ── Normalize components to [0, 1] ──────────────────────────────────
    # Calibrated using empirical range from correlation analysis:
    #   gap_entropy: 2.5 (CNR L+L) to 10.8 (LFR orig) → use [0, 12]
    #   avg_log_gap: 2.4 (CNR L+L) to 8.3 (LFR orig)  → use [0, 12]
    # Fixed range avoids graph-size dependence from log₂(n).

    # Empirical range from correlation analysis:
    #   gap_entropy: 2.5 (CNR best) to 10.8 (LFR random) — use sigmoid mapping
    #   avg_log_gap: 2.4 (CNR best) to 8.3 (LFR random) — use sigmoid mapping
    # Sigmoid centering: midpoint at gap_entropy ≈ 6.5 (moderate compressibility)
    # Steepness chosen so CNR (2.9) → ~0.92, LFR orig (10.8) → ~0.05

    s_gap_entropy = _inv_sigmoid(gap_entropy, 6.0, 0.6)
    s_log_gap     = _inv_sigmoid(avg_log_gap, 5.0, 0.7)
    s_copy_frac   = avg_copy_fraction          # already [0, 1]
    s_res_gap1    = res_gap1_frac              # already [0, 1]
    s_residuals   = avg_deg > 0 ? clamp(1.0 - avg_residuals_per_v / avg_deg, 0.0, 1.0) : 0.0

    # ── Weighted combination ────────────────────────────────────────────
    # Primary weights from correlation |r|, with residual_gap1 boosted
    # to discriminate orderings with similar gap metrics (e.g. CNR orig vs L+L)
    w_gap_entropy = 0.30
    w_copy_frac   = 0.20
    w_log_gap     = 0.15
    w_res_gap1    = 0.20   # boosted: key tiebreaker for top-quality orderings
    w_residuals   = 0.15

    raw = w_gap_entropy * s_gap_entropy +
          w_copy_frac   * s_copy_frac +
          w_log_gap     * s_log_gap +
          w_res_gap1    * s_res_gap1 +
          w_residuals   * s_residuals

    # Rescale from empirical [0.03, 0.82] range to [0.02, 0.95]
    score = clamp((raw - 0.03) / 0.79 * 0.93 + 0.02, 0.0, 1.0)

    components = Dict{Symbol,Float64}(
        :gap_entropy_raw      => gap_entropy,
        :gap_entropy_norm     => s_gap_entropy,
        :avg_log_gap_raw      => avg_log_gap,
        :avg_log_gap_norm     => s_log_gap,
        :copy_fraction_raw    => avg_copy_fraction,
        :copy_fraction_norm   => s_copy_frac,
        :residual_gap1_raw    => res_gap1_frac,
        :residual_gap1_norm   => s_res_gap1,
        :residuals_per_v_raw  => avg_residuals_per_v,
        :residuals_norm       => s_residuals,
    )

    return (score=score, components=components)
end

end # module Metrics
