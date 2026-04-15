#
# BiGG metabolic model → directed graph.
# Vertices = metabolites. Edges = reaction connections (substrate → product,
# plus reverse for reversible reactions).
#

# ── Currency metabolites excluded by default (ubiquitous side-metabolites) ──
const CURRENCY_METABOLITES = Set([
    "h2o", "atp", "adp", "amp", "nad", "nadh", "nadp", "nadph",
    "coa", "co2", "h", "pi", "ppi", "o2", "nh4",
    "fad", "fadh2", "gtp", "gdp", "gmp",
    "udp", "ump", "utp", "cdp", "cmp", "ctp",
])

"""
    load_bigg_metabolic(model_path; compartment="c", keep_currency=false)

Parse a BiGG JSON model into (`met_ids`, `met_names`, `edge_list`, `n`).
`edge_list` is a vector of `(src_idx, tgt_idx)` tuples (1-based). Currency
metabolites (H₂O, ATP/ADP, NAD(P)(H), H⁺, Pi, O₂, CoA, …) are excluded by
default to avoid shortcut edges that dominate path scoring.
"""
function load_bigg_metabolic(model_path::String;
                              compartment::Union{String, Vector{String}}="c",
                              keep_currency::Bool=false)
    comp_set = compartment isa String ? Set([compartment]) : Set(compartment)
    @info "Loading BiGG model from $model_path (compartments=$(sort(collect(comp_set))), keep_currency=$keep_currency)"

    d = JSON.parsefile(model_path)

    met_ids = String[]
    met_names = String[]
    met_set = Set{String}()

    for m in d["metabolites"]
        m["compartment"] in comp_set || continue
        mid = m["id"]::String
        base_id = replace(mid, r"_[a-z]$" => "")
        if !keep_currency && base_id in CURRENCY_METABOLITES
            continue
        end
        push!(met_ids, mid)
        push!(met_names, get(m, "name", mid))
        push!(met_set, mid)
    end

    id_to_idx = Dict{String,Int}()
    for (i, mid) in enumerate(met_ids)
        id_to_idx[mid] = i
    end
    n = length(met_ids)
    @info "Metabolites in compartment '$compartment': $n ($(keep_currency ? "with" : "without") currency)"

    edges = Set{Tuple{Int,Int}}()
    n_reactions_used = 0

    for r in d["reactions"]
        mets = r["metabolites"]::Dict{String,Any}
        lb = get(r, "lower_bound", 0.0)

        substrates = Int[]
        products = Int[]
        for (mid, coeff) in mets
            haskey(id_to_idx, mid) || continue
            idx = id_to_idx[mid]
            if coeff < 0
                push!(substrates, idx)
            elseif coeff > 0
                push!(products, idx)
            end
        end

        isempty(substrates) && continue
        isempty(products) && continue
        n_reactions_used += 1

        for s in substrates, p in products
            s != p && push!(edges, (s, p))
        end

        if lb < 0
            for p in products, s in substrates
                s != p && push!(edges, (p, s))
            end
        end
    end

    edge_list = collect(edges)
    sort!(edge_list)

    @info "Reactions used: $n_reactions_used, Edges: $(length(edge_list))"

    return (met_ids, met_names, edge_list, n)
end

"""
    save_metabolic_graph(met_ids, met_names, edge_list, n; output_prefix)

Write `<output_prefix>.csv` (edge list) and `<output_prefix>_names.csv`
(vertex ID → BiGG ID + human-readable name).
"""
function save_metabolic_graph(met_ids, met_names, edge_list, n;
                               output_prefix::String="datasets/metabolic/iML1515_c")
    edge_path = output_prefix * ".csv"
    open(edge_path, "w") do f
        println(f, "source,target")
        for (s, t) in edge_list
            println(f, "$s,$t")
        end
    end
    @info "Edge list saved to $edge_path"

    names_path = output_prefix * "_names.csv"
    open(names_path, "w") do f
        println(f, "id,bigg_id,name")
        for (i, (mid, mname)) in enumerate(zip(met_ids, met_names))
            safe_name = replace(mname, "," => ";")
            println(f, "$i,$mid,$safe_name")
        end
    end
    @info "Vertex mapping saved to $names_path ($n metabolites)"

    return (edge_path, names_path)
end

"""
    build_metabolic_graph(model_path; compartment="c", keep_currency=false)
        -> (g, rg, met_ids, met_names, id_to_idx)

High-level loader. Returns forward graph `g`, reverse graph `rg` (for
bidirectional PPR / Yen's), plus identifier vectors + BiGG→vertex-index map.
"""
function build_metabolic_graph(model_path::String;
                                compartment::Union{String, Vector{String}}="c",
                                keep_currency::Bool=false)
    met_ids, met_names, edge_list, n = load_bigg_metabolic(
        model_path; compartment=compartment, keep_currency=keep_currency)

    g = SimpleDiGraph{UInt32}()
    add_vertices!(g, n)
    for (s, t) in edge_list
        add_edge!(g, UInt32(s), UInt32(t))
    end

    rg = SimpleDiGraph{UInt32}()
    add_vertices!(rg, n)
    for (s, t) in edge_list
        add_edge!(rg, UInt32(t), UInt32(s))
    end

    id_to_idx = Dict(mid => i for (i, mid) in enumerate(met_ids))

    @info "Built graph: $(nv(g)) vertices, $(ne(g)) edges"
    return (g, rg, met_ids, met_names, id_to_idx)
end
