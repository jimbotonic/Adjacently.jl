#!/usr/bin/env julia

# Study RCGE with pair-local inter lists across depths for CNR-2000

include("run_tests_main.jl")
using Logging
using Adjacently
using Adjacently.IO: load_adjacency_list_from_csv, BitWriter, flush_bitwriter
using Adjacently.Clustering: leiden_partition, aggregate_graph
using Adjacently.InterEncoding: cluster_density, should_stop_coarsening_by_density
using Adjacently.RCGE: encode_level, RCGEParams
using Adjacently.Graph: subgraph
using LightGraphs: nv, outneighbors, is_directed
import LightGraphs

@testset "RCGE local-pairs study (CNR-2000)" begin
    prev = current_logger()
    global_logger(ConsoleLogger(stderr, Logging.Info))
    try
        csv = "datasets/webgraph/cnr-2000/cnr-2000.csv"
        if !isfile(csv)
            @warn "CNR-2000 CSV not found at $csv; skipping"
            @test_skip "dataset unavailable"
            return
        end
        @info "Loading CNR-2000..."
        g = load_adjacency_list_from_csv(csv, ',', true)
        n = nv(g)
        @info "Graph loaded: n=$(n)"

        # Helper: count directed edges
        count_edges(h) = sum(length(outneighbors(h, v)) for v in 1:nv(h))

        # Convert a WeightedCoarseGraph to an adjacency TestGraph (repeat edges)
        function coarse_to_testgraph(Gc)
            nC = Gc.n
            Tt = (typeof(g)).parameters[1]
            adj = Dict{Int, Vector{Tt}}()
            for u in 1:nC
                lst = UInt32[]
                for (v, w) in Gc.out_w[u]
                    for _ in 1:Int(round(w))
                        push!(lst, Tt(v))
                    end
                end
                if !isempty(lst)
                    adj[u] = lst
                end
            end
            return (nC=nC, adj=adj, Tt=Tt)
        end

        # Timing helper
        ms_since(t0) = round(Int, (time_ns() - t0) / 1e6)

        # Compute per-cluster densities and a closer intra bpe estimate using the encoder on induced subgraphs
        function compute_cluster_stats(g, clusters, params)
            @info "Computing cluster stats for $(length(clusters)) clusters..."
            t_stats = time_ns()
            # densities
            dens = [cluster_density(g, C) for C in clusters]
            @info "Computed densities in $(ms_since(t_stats))ms"
            # intra bits per cluster and edges (measured via RCGE intra encoder on induced subgraph)
            per_cluster_bits = Float64[]
            per_cluster_edges = Int[]
            @info "Starting per-cluster encoding..."
            for (ci, C) in enumerate(clusters)
                t_c = time_ns()
                # Use existing subgraph function
                h, oni, noi = subgraph(g, C)
                s = length(C)
                Vt = (typeof(g)).parameters[1]
                mC = sum(length(outneighbors(h, v)) for v in 1:nv(h))
                # Encode one level with a single cluster [1..s] and measure intra bits as used by RCGE
                io = IOBuffer(); w = BitWriter(io)
                statsC = Adjacently.RCGE.RCGEStats()
                cluster_local = Vt[ Vt(i) for i in 1:s ]
                encode_level(w, h, [cluster_local]; params=params, stats=statsC)
                flush_bitwriter(w; flush_last_bits=true)
                push!(per_cluster_bits, Float64(statsC.bits_intra))
                push!(per_cluster_edges, mC)
                @info "Cluster $(ci)/$(length(clusters)) (size=$(s), edges=$(mC)) took $(ms_since(t_c))ms"
            end
            @info "All cluster stats took $(ms_since(t_stats))ms"
            return dens, per_cluster_bits, per_cluster_edges
        end

        # Partition helper (available for all levels)
        function build_partition_for_K(g2, part0::Vector{Int}, K::Int)
            n2 = nv(g2)
            counts = Dict{Int,Int}(); for c in part0; counts[c] = get(counts,c,0)+1; end
            labels_sorted2 = sort(collect(keys(counts)), by = c -> -counts[c])
            topK = min(K, length(labels_sorted2))
            if topK == 0
                return Int[], Vector{Vector{typeof(g2).parameters[1]}}()
            end
            top2 = labels_sorted2[1:topK]
            top_index2 = Dict{Int,Int}(c => i for (i,c) in enumerate(top2))
            partK2 = similar(part0)
            for i2 in 1:n2
                c = part0[i2]
                partK2[i2] = haskey(top_index2, c) ? top_index2[c] : 0
            end
            for i2 in 1:n2
                if partK2[i2] == 0
                    counts_c = zeros(Int, topK)
                    for v2 in outneighbors(g2, i2)
                        c0 = part0[Int(v2)]
                        if haskey(top_index2, c0)
                            counts_c[top_index2[c0]] += 1
                        end
                    end
                    bestc = 1; bestw = -1
                    for cidx in 1:topK
                        if counts_c[cidx] > bestw
                            bestw = counts_c[cidx]; bestc = cidx
                        end
                    end
                    partK2[i2] = clamp(bestc, 1, topK)
                end
            end
            Vt2 = (typeof(g2)).parameters[1]
            clusters2 = [Vt2[] for _ in 1:topK]
            for i2 in 1:n2
                idx = partK2[i2]
                if 1 <= idx <= topK
                    push!(clusters2[idx], Vt2(i2))
                else
                    push!(clusters2[1], Vt2(i2))
                end
            end
            clusters2 = filter(!isempty, clusters2)
            return partK2, clusters2
        end

        # Depth sweep
        MAX_LEVELS = try parse(Int, get(ENV, "RCGE_MAX_LEVELS", "3")) catch; 3 end
        KMAX = try parse(Int, get(ENV, "RCGE_K_MAX", "32")) catch; 32 end
        K_SELECT = get(ENV, "RCGE_K_SELECT", "modularity")  # or "ncut"
        min_clusters = try parse(Int, get(ENV, "RCGE_MIN_CLUSTERS", "2")) catch; 2 end
        min_density = try parse(Float64, get(ENV, "RCGE_MIN_CLUSTER_DENSITY", "0.01")) catch; 0.01 end
        min_cl_size = try parse(Int, get(ENV, "RCGE_MIN_CLUSTER_SIZE", "16")) catch; 16 end

        # Encoder params (used for selection and encoding at each level)
        params = RCGEParams(L=128, varint=:fibonacci, count_varint=:fibonacci, gap=:fibonacci, degree=:elias_delta, undirected_pairs=false, perm_strategy=:blockpos, membership=:elias_fano, inter_strategy=:lists, intra_ref_enabled=true, intra_ref_window=32, intra_ref_rle=false, intra_block_try=false, positions_mode=:delta, additions_mode=:delta)

        cur = g
        level = 1
        # Keep current clusters across levels; initialize on level 1 via K=2 selection
        clusters = Vector{Vector{(typeof(g)).parameters[1]}}()
        clusters_initialized = false
        prev_n = nv(cur)
        while level <= MAX_LEVELS
            level_start = time_ns()
            ncur = nv(cur); mcur = count_edges(cur)
            @info "Level $level: n=$ncur m=$mcur clusters=$(clusters_initialized ? string(length(clusters)) : "-")"

            # Initialize clusters only once at level 1; later levels reuse and split them
            if !clusters_initialized
                @info "Starting leiden_partition..."
                t0 = time_ns()
                part0 = leiden_partition(cur; max_passes=8, max_levels=5)
                @info "leiden_partition took $(ms_since(t0))ms"

            # (removed K-sweep helpers and logs; building K=2 directly)
            # Directly build K=2 partition for binary tree (no K sweep logging)
            @info "Building K=2 partition for binary tree..."
            partK2, clusters_init = build_partition_for_K(cur, part0, 2)
            clusters = clusters_init
            clusters_initialized = true
            else
                @info "Skipping K selection (using recursive clusters from previous level)"
            end
            # Compute stats

            # Encode and get overall bpe + intra/inter bpe
            @info "Encoding level with K=2..."
            t_encode = time_ns()
            io = IOBuffer(); w = BitWriter(io)
            stats = Adjacently.RCGE.RCGEStats()
            encode_level(w, cur, clusters; params=params, stats=stats)
            flush_bitwriter(w; flush_last_bits=true)
            bytes = take!(io)
            @info "Encoding took $(ms_since(t_encode))ms"
            # Fractions of intra/inter edges under best partition
            intra_edges = 0
            # Build label map from current clusters for intra/inter counts
            labels = Dict{Int,Int}()
            for (ci, C) in enumerate(clusters)
                for u in C
                    labels[Int(u)] = ci
                end
            end
            for i in 1:ncur
                li = get(labels, i, 0)
                for v in outneighbors(cur, i)
                    if get(labels, Int(v), -1) == li && li != 0
                        intra_edges += 1
                    end
                end
            end
            inter_edges = max(mcur - intra_edges, 0)
            # Bits used
            total_bits = 8.0 * length(bytes)
            intra_bits = stats.bits_intra
            inter_bits = stats.bits_inter_headers + stats.bits_inter_lists
            # Averages with correct denominators
            avg_intra_bpe = intra_bits / max(intra_edges, 1)
            avg_inter_bpe = inter_bits / max(inter_edges, 1)
            frac_intra = intra_edges / max(mcur, 1)
            weighted_bpe = avg_intra_bpe * frac_intra + avg_inter_bpe * (1 - frac_intra)
            bpe = total_bits / max(mcur, 1)
            @info "Stats: overall bpe=$(round(bpe,digits=4)), weighted_bpe=$(round(weighted_bpe,digits=4)), avg_intra_bpe=$(round(avg_intra_bpe,digits=4)), avg_inter_bpe=$(round(avg_inter_bpe,digits=4)), frac_intra=$(round(frac_intra,digits=4))"

            # (Removed duplicate helper definitions and global K sweep)

            # No global K sweep beyond level selection; we recurse to next level

            # Per-AB metrics
            ab_entries = getfield(stats, :ab_metrics)
            if !isempty(ab_entries)
                # sort by m_AB (x[3]) descending
                sorted_ab = sort(ab_entries, by = x -> -x[2])
                sorted_ab = sort(ab_entries, by = x -> -x[3])
                topN = min(5, length(sorted_ab))
                for i in 1:topN
                    ia, ib, mAB, bitsAB, refc, rawc = sorted_ab[i]
                    bp = bitsAB / max(mAB, 1)
                    @info "AB (ia=$(ia), ib=$(ib)): m=$(mAB) bits=$(bitsAB) inter_bpe=$(round(bp,digits=4)) refs=$(refc) raw=$(rawc)"
                end
            end

            dens, per_cluster_bits, per_cluster_edges = compute_cluster_stats(cur, clusters, params)
            @info "Cluster densities: min=$(minimum(dens)), max=$(maximum(dens)), avg=$(round(sum(dens)/length(dens),digits=4))"
            total_bits_approx = sum(per_cluster_bits)
            total_edges_approx = max(sum(per_cluster_edges), 1)
            avg_cluster_bpe_approx = total_bits_approx / total_edges_approx
            @info "Weighted approx intra bpe per cluster: $(round(avg_cluster_bpe_approx,digits=4))"


            # Decide recursion; if stopping by density/size, log and stop
            recursing = true
            stop_by_density = should_stop_coarsening_by_density(cur, clusters, min_density)
            stop_by_size = minimum(length.(clusters)) <= min_cl_size
            if stop_by_density || stop_by_size
                recursing = false
                @info "Stopping: recursing=false because density>= $(min_density) for all? $(stop_by_density) or min_cluster_size<= $(min_cl_size)? $(stop_by_size)"
                break
            else
                @info "Recursing: true (density/size OK)"
            end
            # Binary refinement: split each cluster into 2 children using subgraph partition
            function split_cluster_into_two(g, C::Vector)
                s = length(C)
                if s <= 2
                    return [C]
                end
                # Use existing subgraph function
                Vt = (typeof(g)).parameters[1]
                h, oni, noi = subgraph(g, C)
                # Partition and cap to K=2
                p0 = leiden_partition(h; max_passes=6, max_levels=4)
                partK, kids = build_partition_for_K(h, p0, 2)
                if length(kids) < 2
                    # Fallback: split by order
                    mid = cld(s,2)
                    return [C[1:mid], C[mid+1:end]]
                end
                # Map back to global ids: kids are local ids in h; use noi to recover originals
                child1 = Vt[ noi[i] for i in kids[1] ]
                child2 = Vt[ noi[i] for i in kids[2] ]
                return [child1, child2]
            end
            @info "Starting binary cluster splitting..."
            t_split = time_ns()
            next_clusters = eltype(clusters)[]
            for (ci, C) in enumerate(clusters)
                @info "Splitting cluster $(ci)/$(length(clusters)) (size=$(length(C)))..."
                for Kchild in split_cluster_into_two(cur, C)
                    if !isempty(Kchild)
                        push!(next_clusters, Kchild)
                    end
                end
            end
            @info "Binary splitting took $(ms_since(t_split))ms, produced $(length(next_clusters)) clusters"
            clusters = next_clusters
            @info "Level $(level) total time: $(ms_since(level_start))ms"
            level += 1
        end
    finally
        global_logger(prev)
    end
end
