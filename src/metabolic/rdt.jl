#
# Load RDT (Reaction Decoder Tool) atom-mapping-derived edge weights.
#
# Reads a TSV of (substrate_bigg, product_bigg, conservation_fraction) triples
# into a Dict{Tuple{UInt32,UInt32}, Float64} indexed by graph vertex pairs.
# Conservation ≈ fraction of heavy atoms traced between substrate and product
# by RDT v3.3.0. 1.0 = carbon-skeleton-preserving main pair; low = cofactor-like.
#

"""
    load_rdt_weights(path, id_to_idx) -> Dict{Tuple{UInt32,UInt32}, Float64}

Read the RDT TSV (produced by `scripts/run_rdt_batch.jl`) and index by graph
vertex IDs. Only entries where both BiGG IDs appear in `id_to_idx` are kept.
Missing file → empty dict (with a warning).
"""
function load_rdt_weights(path::String, id_to_idx::Dict{String, Int})
    weights = Dict{Tuple{UInt32,UInt32}, Float64}()
    isfile(path) || (@warn "RDT TSV not found: $path"; return weights)
    open(path) do io
        _ = readline(io)  # skip header
        for line in eachline(io)
            fields = split(line, '\t')
            length(fields) >= 3 || continue
            sub = String(fields[1]); prod = String(fields[2])
            (haskey(id_to_idx, sub) && haskey(id_to_idx, prod)) || continue
            c = tryparse(Float64, String(fields[3]))
            c === nothing && continue
            weights[(UInt32(id_to_idx[sub]), UInt32(id_to_idx[prod]))] = c
        end
    end
    return weights
end
