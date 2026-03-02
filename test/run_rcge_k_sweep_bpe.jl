#!/usr/bin/env julia

# RCGE K-sweep BPE study without heavy deps

using Logging
using LightGraphs: nv, ne, outneighbors
using Adjacently
using Adjacently.IO: load_adjacency_list_from_csv, BitWriter, flush_bitwriter
using Adjacently.Clustering: leiden_partition
using Adjacently.RCGE: encode_level, RCGEParams

function build_partition_for_K(g, part0::Vector{Int}, K::Int)
    n = nv(g)
    # Count labels
    counts = Dict{Int,Int}(); for c in part0; counts[c] = get(counts,c,0)+1; end
    labels_sorted = sort(collect(keys(counts)), by = c -> -counts[c])
    topK = min(K, length(labels_sorted))
    if topK == 0
        return Int[], Vector{Vector{eltype(vertices(g))}}()
    end
    top = labels_sorted[1:topK]
    top_index = Dict{Int,Int}(c => i for (i,c) in enumerate(top))
    partK = similar(part0)
    # Assign topK clusters; mark extras with 0
    for i in 1:n
        c = part0[i]
        partK[i] = haskey(top_index, c) ? top_index[c] : 0
    end
    # Reallocate extras by edge majority to topK clusters
    for i in 1:n
        if partK[i] == 0
            counts_c = zeros(Int, topK)
            for v in outneighbors(g, i)
                c0 = part0[Int(v)]
                if haskey(top_index, c0)
                    counts_c[top_index[c0]] += 1
                end
            end
            bestc = 1; bestw = -1
            for cidx in 1:topK
                if counts_c[cidx] > bestw
                    bestw = counts_c[cidx]; bestc = cidx
                end
            end
            partK[i] = bestc
        end
    end
    # Build clusters (compact, 1..keff)
    Vt = (typeof(g)).parameters[1]
    clusters = [Vt[] for _ in 1:topK]
    for i in 1:n
        push!(clusters[partK[i]], Vt(i))
    end
    clusters = filter(!isempty, clusters)
    return partK, clusters
end

function compute_bpe_stats(g, clusters, params)
    ncur = nv(g)
    # Count edges
    mcur = sum(length(outneighbors(g, v)) for v in 1:ncur)
    # Encode
    io = IOBuffer(); w = BitWriter(io)
    stats = Adjacently.RCGE.RCGEStats()
    encode_level(w, g, clusters; params=params, stats=stats)
    flush_bitwriter(w; flush_last_bits=true)
    bytes = take!(io)
    # Fractions intra/inter
    intra_edges = 0
    # Build quick label map
    lbl = Dict{Int,Int}()
    for (ci, C) in enumerate(clusters)
        for u in C
            lbl[Int(u)] = ci
        end
    end
    for i in 1:ncur
        li = lbl[i]
        for v in outneighbors(g, i)
            if lbl[Int(v)] == li
                intra_edges += 1
            end
        end
    end
    inter_edges = max(mcur - intra_edges, 0)
    total_bits = 8.0 * length(bytes)
    intra_bits = stats.bits_intra
    inter_bits = stats.bits_inter_headers + stats.bits_inter_lists
    avg_intra_bpe = intra_bits / max(intra_edges, 1)
    avg_inter_bpe = inter_bits / max(inter_edges, 1)
    frac_intra = intra_edges / max(mcur, 1)
    weighted_bpe = avg_intra_bpe * frac_intra + avg_inter_bpe * (1 - frac_intra)
    overall_bpe = total_bits / max(mcur, 1)
    return (; overall_bpe, weighted_bpe, avg_intra_bpe, avg_inter_bpe, frac_intra, mcur, total_bits, intra_bits, inter_bits)
end

function main()
    prev = current_logger()
    global_logger(ConsoleLogger(stderr, Logging.Info))
    try
        csv = get(ENV, "RCGE_DATASET", "datasets/webgraph/cnr-2000/cnr-2000.csv")
        if !isfile(csv)
            @error "Dataset not found" csv
            return
        end
        @info "Loading dataset..." csv
        g = load_adjacency_list_from_csv(csv, ',', true)
        n = nv(g)
        @info "Graph loaded" n=n

        # Initial partition once
        @info "Computing base Leiden partition..."
        part0 = leiden_partition(g; max_passes=8, max_levels=5)
        @info "Base partition labels" unique_labels=length(unique(part0))

        # Params
        params = RCGEParams(L=128, varint=:fibonacci, count_varint=:fibonacci, gap=:fibonacci,
                            degree=:elias_delta, undirected_pairs=false, perm_strategy=:blockpos,
                            membership=:elias_fano, inter_strategy=:lists, intra_ref_enabled=true,
                            intra_ref_window=32, intra_ref_rle=false,
                            intra_block_try=false, positions_mode=:delta, additions_mode=:delta)

        KSET = [2,4,8,16,32,64,128,256,512,1024,2048]
        best = nothing
        for K in KSET
            partK, clusters = build_partition_for_K(g, part0, K)
            keff = length(clusters)
            if keff == 0
                @info "K=$K produced no clusters; skipping"
                continue
            end
            stats = compute_bpe_stats(g, clusters, params)
            @info "K sweep" K=K keff=keff overall_bpe=round(stats.overall_bpe; digits=4) weighted_bpe=round(stats.weighted_bpe; digits=4) avg_intra_bpe=round(stats.avg_intra_bpe; digits=4) avg_inter_bpe=round(stats.avg_inter_bpe; digits=4) frac_intra=round(stats.frac_intra; digits=4)
            if best === nothing || stats.overall_bpe < best.overall_bpe
                best = (K=K, keff=keff, stats=stats)
            end
        end
        if best !== nothing
            s = best.stats
            @info "Best by overall bpe" K=best.K keff=best.keff overall_bpe=round(s.overall_bpe; digits=4) weighted_bpe=round(s.weighted_bpe; digits=4) avg_intra_bpe=round(s.avg_intra_bpe; digits=4) avg_inter_bpe=round(s.avg_inter_bpe; digits=4) frac_intra=round(s.frac_intra; digits=4)
        end
    finally
        global_logger(prev)
    end
end

main()

