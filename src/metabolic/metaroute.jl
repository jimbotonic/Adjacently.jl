#
# MetaRoute-style edge weights from BiGG chemical formulas.
#
# Original MetaRoute (Blum & Kohlbacher 2008) separates "main" substrate-product
# pairs (carbon-skeleton-preserving) from "cofactor" pairs using KEGG's RPAIR
# table. RPAIR was deprecated in 2016, so we approximate the main-pair labeling
# using only BiGG-derived atom overlap:
#
#   For each reaction R:
#     For each (s, p) ∈ substrates × products: compute atom_overlap(s, p).
#     Rank all pairs in R by overlap (descending).
#     Rank 1 = "main" pair (carbon skeleton preserved).
#
#   An edge (u, v) may appear in multiple reactions. Its final score is the
#   inverse of its average rank across reactions:
#     mr_score(u, v) = 1 / mean_rank(u, v)  ∈ [1/N, 1]
#   Main-pair edges → rank ≈ 1 → score ≈ 1.0.
#   Cofactor-like edges → rank ≈ N → score ≈ 1/N.
#

"""
    build_metaroute_weights(model_path, g, id_to_idx) -> Dict{Tuple{UInt32,UInt32},Float64}

Rank-based MetaRoute-style edge weights in `[0.01, 1.0]`. Requires
`parse_formula` and `atom_similarity` from `atom_weights.jl` — include order
matters when used outside this module.
"""
function build_metaroute_weights(model_path::String, g, id_to_idx)
    d = JSON.parsefile(model_path)

    formula_of = Dict{String, Dict{String,Int}}()
    for m in d["metabolites"]
        f = get(m, "formula", nothing)
        (f === nothing || f == "" || f == "null") && continue
        formula_of[m["id"]] = parse_formula(f)
    end

    pair_ranks = Dict{Tuple{String,String}, Vector{Int}}()

    for r in d["reactions"]
        mets = r["metabolites"]
        lb = get(r, "lower_bound", 0.0)

        substrates = String[]; products = String[]
        for (mid, coeff) in mets
            coeff < 0 && push!(substrates, mid)
            coeff > 0 && push!(products, mid)
        end
        (isempty(substrates) || isempty(products)) && continue

        scored = Tuple{String,String,Float64}[]
        for s in substrates, p in products
            s == p && continue
            (haskey(formula_of, s) && haskey(formula_of, p)) || continue
            push!(scored, (s, p, atom_similarity(formula_of[s], formula_of[p])))
        end
        isempty(scored) && continue

        sort!(scored, by=x -> -x[3])
        for (rank_idx, (s, p, _)) in enumerate(scored)
            key = (s, p)
            haskey(pair_ranks, key) || (pair_ranks[key] = Int[])
            push!(pair_ranks[key], rank_idx)
            if lb < 0
                rkey = (p, s)
                haskey(pair_ranks, rkey) || (pair_ranks[rkey] = Int[])
                push!(pair_ranks[rkey], rank_idx)
            end
        end
    end

    idx_to_bigg = Dict{Int, String}()
    for (bigg, idx) in id_to_idx
        idx_to_bigg[idx] = bigg
    end

    weights = Dict{Tuple{UInt32,UInt32}, Float64}()
    for e in edges(g)
        u, v = Int(src(e)), Int(dst(e))
        ub = get(idx_to_bigg, u, "")
        vb = get(idx_to_bigg, v, "")
        if isempty(ub) || isempty(vb)
            weights[(UInt32(u), UInt32(v))] = 0.01
            continue
        end
        ranks = get(pair_ranks, (ub, vb), Int[])
        if isempty(ranks)
            weights[(UInt32(u), UInt32(v))] = 0.01
        else
            mean_rank = sum(ranks) / length(ranks)
            weights[(UInt32(u), UInt32(v))] = 1.0 / mean_rank
        end
    end
    return weights
end
