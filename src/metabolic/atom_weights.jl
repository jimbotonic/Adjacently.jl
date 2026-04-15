#
# Chemistry-aware edge weights from BiGG chemical formulas.
#
# atom_similarity(u, v) = shared heavy atoms / min(|u|, |v|) over C/N/O/P/S.
# Min-denominator handles carrier molecules (acetyl-CoA, ACP esters) correctly
# so that pyr→accoa still scores ~1.0 despite the CoA carrier.
#
# This is a static graph property — no FBA needed, so atom weights don't have
# the zero-flux coverage problem that limits flux-weighted methods.
#

const ATOM_ELEMENTS = ("C", "N", "O", "P", "S")

"""
    parse_formula(formula) -> Dict{String, Int}

Parse a chemical formula string like "C6H12O6" into an element-count dict.
"""
function parse_formula(formula::String)
    counts = Dict{String, Int}()
    for m in eachmatch(r"([A-Z][a-z]?)(\d*)", formula)
        elem = String(m.captures[1])
        s = m.captures[2]
        n = (s === nothing || s == "") ? 1 : parse(Int, s)
        counts[elem] = get(counts, elem, 0) + n
    end
    return counts
end

heavy_atom_count(f::Dict{String,Int}) = sum(get(f, e, 0) for e in ATOM_ELEMENTS)

"""
    atom_similarity(f1, f2) -> Float64 ∈ [0, 1]

Fraction of the smaller molecule's heavy-atom skeleton shared with its partner.
Ignores hydrogen and trace elements. See file header for rationale.
"""
function atom_similarity(f1::Dict{String,Int}, f2::Dict{String,Int})
    shared = 0
    for e in ATOM_ELEMENTS
        shared += min(get(f1, e, 0), get(f2, e, 0))
    end
    n_min = min(heavy_atom_count(f1), heavy_atom_count(f2))
    return n_min == 0 ? 0.0 : shared / n_min
end

"""
    build_atom_weights(model_path, g, id_to_idx;
                       min_sim=0.01, transform=:linear) -> Dict

Build per-edge atom-similarity weights. `transform=:linear` gives similarity
directly (good for wPPR / multiplicative cost functions); `:neg_log` gives
`-log(max(sim, min_sim))` (good for Dijkstra-style costs).
"""
function build_atom_weights(model_path::String, g, id_to_idx;
                            min_sim::Float64=0.01,
                            transform::Symbol=:linear)
    d = JSON.parsefile(model_path)
    formula_of = Dict{String, Dict{String,Int}}()
    for m in d["metabolites"]
        f = get(m, "formula", nothing)
        (f === nothing || f == "" || f == "null") && continue
        formula_of[m["id"]] = parse_formula(f)
    end
    vertex_bigg = Dict{Int, String}()
    for (bigg, i) in id_to_idx
        vertex_bigg[i] = bigg
    end

    weights = Dict{Tuple{UInt32,UInt32}, Float64}()
    for e in edges(g)
        u, v = src(e), dst(e)
        bu = get(vertex_bigg, Int(u), nothing)
        bv = get(vertex_bigg, Int(v), nothing)
        sim = if bu !== nothing && bv !== nothing &&
                haskey(formula_of, bu) && haskey(formula_of, bv)
            atom_similarity(formula_of[bu], formula_of[bv])
        else
            0.0
        end
        sim = max(sim, min_sim)
        w = if transform == :neg_log
            -log(sim)
        else
            sim
        end
        weights[(UInt32(u), UInt32(v))] = w
    end
    return weights
end

"""
    atom_weight_fn(weights_dict)

Wrap a weight dict as a Yen-compatible weight function of signature
`(g, rg, u, v) -> weights_dict[(u, v)]`.
"""
function atom_weight_fn(weights::Dict{Tuple{UInt32,UInt32}, Float64})
    return (g, rg, u, v) -> get(weights, (u, v), 1.0)
end
