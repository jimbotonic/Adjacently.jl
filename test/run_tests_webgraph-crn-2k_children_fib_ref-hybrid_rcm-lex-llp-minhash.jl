#!/usr/bin/env julia

# Test CNR-2000 with different vertex relabeling algorithms:
# - Original order (no relabeling)
# - Lexicographic order (sorted by out-degree)  
# - RCM (Reverse Cuthill-McKee) ordering
# - LLP (Local Layout Permutation) ordering
# - MinHash-based ordering
# 
# All using: children mode, Fibonacci encoding, reference + hybrid compression

using Pkg
Pkg.activate(".")

using Test
using Printf
using LightGraphs: nv, ne, outneighbors
using Adjacently.IO: load_adjacency_list_from_csv
using Adjacently.Relabeling: relabel_vertices_rcm, relabel_graph
using Adjacently.MGS: write_compressed_mgs3_graph, load_compressed_mgs3_graph
using Statistics: mean, median

"""
    create_lexicographic_relabeling(g)

Create lexicographic ordering based on out-degree (ascending order).
"""
function create_lexicographic_relabeling(g)
    T = eltype(g)
    n = nv(g)
    
    # Get out-degrees for all vertices
    degrees = [(T(v), length(outneighbors(g, v))) for v in 1:n]
    
    # Sort by out-degree (lexicographic order)
    sort!(degrees, by = x -> x[2])
    
    # Create mapping: old_vertex -> new_vertex
    relabel_map = Dict{T,T}()
    for (new_v, (old_v, _)) in enumerate(degrees)
        relabel_map[old_v] = T(new_v)
    end
    
    return relabel_graph(g, relabel_map)
end

"""
    create_llp_relabeling(g)

Create Local Layout Permutation ordering - a simplified version focusing on 
local neighborhood clustering.
"""
function create_llp_relabeling(g)
    T = eltype(g)
    n = nv(g)
    visited = falses(n)
    new_order = Int[]
    
    # Start with highest degree vertex
    degrees = [length(outneighbors(g, v)) for v in 1:n]
    start_vertex = argmax(degrees)
    
    # BFS-like traversal prioritizing high-degree unvisited neighbors
    queue = [start_vertex]
    visited[start_vertex] = true
    push!(new_order, start_vertex)
    
    while length(new_order) < n
        if isempty(queue)
            # Find next unvisited vertex with highest degree
            remaining = findall(!visited)
            if isempty(remaining)
                break
            end
            next_vertex = remaining[argmax([degrees[v] for v in remaining])]
            queue = [next_vertex]
        end
        
        current = popfirst!(queue)
        
        # Add unvisited neighbors, sorted by degree (descending)
        neighbors = collect(outneighbors(g, current))
        unvisited_neighbors = filter(v -> !visited[v], neighbors)
        sort!(unvisited_neighbors, by = v -> degrees[v], rev = true)
        
        for neighbor in unvisited_neighbors
            if !visited[neighbor]
                visited[neighbor] = true
                push!(new_order, neighbor)
                push!(queue, neighbor)
            end
        end
    end
    
    # Handle any remaining unvisited vertices
    for v in 1:n
        if !visited[v]
            push!(new_order, v)
        end
    end
    
    # Create relabeling map
    relabel_map = Dict{T,T}()
    for (new_v, old_v) in enumerate(new_order)
        relabel_map[T(old_v)] = T(new_v)
    end
    
    return relabel_graph(g, relabel_map)
end

"""
    create_minhash_relabeling(g)

Create MinHash-based ordering using Jaccard similarity of neighborhoods.
Groups similar vertices together.
"""
function create_minhash_relabeling(g)
    T = eltype(g)
    n = nv(g)
    
    # Simplified MinHash: use first k neighbors as "signature"
    k = 5  # signature size
    signatures = []
    
    for v in 1:n
        neighbors = sort(collect(outneighbors(g, v)))
        # Take first k neighbors as signature (pad with zeros if needed)
        signature = zeros(Int, k)
        for i in 1:min(k, length(neighbors))
            signature[i] = neighbors[i]
        end
        push!(signatures, (T(v), signature))
    end
    
    # Sort by signature (lexicographic comparison)
    sort!(signatures, by = x -> x[2])
    
    # Create relabeling map
    relabel_map = Dict{T,T}()
    for (new_v, (old_v, _)) in enumerate(signatures)
        relabel_map[T(old_v)] = T(new_v)
    end
    
    return relabel_graph(g, relabel_map)
end

"""
    calculate_locality_metrics(g)

Calculate locality metrics for the graph ordering.
"""
function calculate_locality_metrics(g)
    n = nv(g)
    gaps = Float64[]
    
    for v in 1:n
        neighbors = collect(outneighbors(g, v))
        if !isempty(neighbors)
            # Calculate gaps between vertex and its neighbors
            vertex_gaps = [abs(neighbor - v) for neighbor in neighbors]
            append!(gaps, vertex_gaps)
        end
    end
    
    avg_gap = isempty(gaps) ? 0.0 : mean(gaps)
    median_gap = isempty(gaps) ? 0.0 : median(gaps)
    
    # Locality score: lower is better (normalized by graph size)
    locality_score = avg_gap / n
    
    return (
        avg_neighbor_gap = avg_gap,
        median_neighbor_gap = median_gap,
        locality_score = locality_score
    )
end

@testset "CNR-2000 Vertex Relabeling Comparison Test" begin
    # Load CNR-2000 dataset
    cnr_csv = "datasets/webgraph/cnr-2000/cnr-2000.csv"
    
    if !isfile(cnr_csv)
        @warn "CNR-2000 dataset not found at $cnr_csv, skipping test"
        return
    end
    
    println("=== CNR-2000 Vertex Relabeling Comparison Test ===")
    
    # Load graph
    println("Loading CNR-2000...")
    g_original = load_adjacency_list_from_csv(cnr_csv, ',', true)
    n, m = nv(g_original), ne(g_original)
    println("Graph loaded: n=$n, m=$m")
    
    # Create test_data directory if it doesn't exist
    mkpath("test_data")
    
    # Storage for results
    results = []
    
    # Test configurations
    relabeling_methods = [
        ("Original", g_original, "No relabeling - original vertex order"),
        ("Lexicographic", nothing, "Lexicographic order by out-degree"),
        ("RCM-Out", nothing, "Reverse Cuthill-McKee with out-degree priority"),
        ("RCM-Sym", nothing, "Reverse Cuthill-McKee with symmetric adjacency"),
        ("LLP", nothing, "Local Layout Permutation ordering"),
        ("MinHash", nothing, "MinHash-based similarity ordering")
    ]
    
    for (method_name, graph, description) in relabeling_methods
        println("\n--- Testing $method_name Relabeling ---")
        println("Description: $description")
        
        # Apply relabeling
        if method_name == "Original"
            g_relabeled = graph
        elseif method_name == "Lexicographic"
            g_relabeled = create_lexicographic_relabeling(g_original)
        elseif method_name == "RCM-Out"
            rcm_map = relabel_vertices_rcm(g_original, :out)
            g_relabeled = relabel_graph(g_original, rcm_map)
        elseif method_name == "RCM-Sym"
            rcm_map = relabel_vertices_rcm(g_original, :sym)
            g_relabeled = relabel_graph(g_original, rcm_map)
        elseif method_name == "LLP"
            g_relabeled = create_llp_relabeling(g_original)
        elseif method_name == "MinHash"
            g_relabeled = create_minhash_relabeling(g_original)
        end
        
        # Verify graph integrity
        @test nv(g_relabeled) == n
        @test ne(g_relabeled) == m
        
        # Test compression
        output_file = "test_data/cnr2k_$(lowercase(replace(method_name, "-" => "_")))"
        
        println("  Compressing with children+fibonacci+reference+hybrid...")
        
        compression_start = time()
        @time write_compressed_mgs3_graph(
            g_relabeled,
            output_file,
            :children,    # children mode
            :fibonacci,   # fibonacci encoding
            true,         # mix encoding enabled (hybrid)
            true          # reference encoding enabled
        )
        compression_time = time() - compression_start
        
        # Check file was created
        mgz_file = "$(output_file).mgz"
        @test isfile(mgz_file)
        
        # Get compression statistics
        filesize = stat(mgz_file).size
        bpv = (filesize * 8) / n  # bits per vertex
        bpe = (filesize * 8) / m  # bits per edge
        
        println("  File size: $(round(filesize/1024/1024; digits=3)) MB")
        println("  Bits per vertex: $(round(bpv; digits=3))")
        println("  Bits per edge: $(round(bpe; digits=3))")
        println("  Compression time: $(round(compression_time; digits=2))s")
        
        # Test decompression
        println("  Testing decompression...")
        decompression_start = time()
        @time g_decoded = load_compressed_mgs3_graph(mgz_file)
        decompression_time = time() - decompression_start
        
        # Verify decompression
        @test nv(g_decoded) == n
        @test ne(g_decoded) == m
        
        # Quick neighbor verification (sample)
        sample_vertices = min(1000, n)
        mismatch_count = 0
        
        for v in 1:sample_vertices
            original_neighbors = sort(collect(outneighbors(g_relabeled, v)))
            decoded_neighbors = sort(collect(outneighbors(g_decoded, v)))
            
            if original_neighbors != decoded_neighbors
                mismatch_count += 1
                if mismatch_count <= 3  # Show first 3 mismatches
                    println("    Mismatch at vertex $v: expected $(length(original_neighbors)), got $(length(decoded_neighbors))")
                end
            end
        end
        
        @test mismatch_count == 0
        
        # Calculate locality metrics
        locality_metrics = calculate_locality_metrics(g_relabeled)
        
        # Store results
        push!(results, (
            method = method_name,
            description = description,
            filesize_bytes = filesize,
            filesize_mb = round(filesize/1024/1024; digits=3),
            bits_per_vertex = round(bpv; digits=3),
            bits_per_edge = round(bpe; digits=3),
            compression_time = round(compression_time; digits=2),
            decompression_time = round(decompression_time; digits=2),
            mismatch_count = mismatch_count,
            avg_neighbor_gap = locality_metrics.avg_neighbor_gap,
            median_neighbor_gap = locality_metrics.median_neighbor_gap,
            locality_score = locality_metrics.locality_score
        ))
        
        println("  Decompression time: $(round(decompression_time; digits=2))s")
        println("  ✓ Verification: $(mismatch_count == 0 ? "PASSED" : "FAILED ($mismatch_count mismatches)")")
        println("  Avg neighbor gap: $(round(locality_metrics.avg_neighbor_gap; digits=1))")
        println("  Locality score: $(round(locality_metrics.locality_score; digits=3))")
    end
    
    # Summary comparison
    println("\n=== COMPRESSION COMPARISON SUMMARY ===")
    println("Method           | Size (MB) | BPV    | BPE   | Comp(s) | Decomp(s) | Locality")
    println("-----------------|-----------|--------|-------|---------|-----------|----------")
    
    # Sort by file size for comparison
    sorted_results = sort(results, by = r -> r.filesize_bytes)
    
    for result in sorted_results
        @printf("%-15s | %8.3f | %6.3f | %5.3f | %7.2f | %9.2f | %8.3f\n",
            result.method,
            result.filesize_mb,
            result.bits_per_vertex,
            result.bits_per_edge,
            result.compression_time,
            result.decompression_time,
            result.locality_score
        )
    end
    
    # Best method analysis
    best_compression = sorted_results[1]
    worst_compression = sorted_results[end]
    
    improvement = (worst_compression.filesize_bytes - best_compression.filesize_bytes) / worst_compression.filesize_bytes * 100
    
    println("\n=== ANALYSIS ===")
    println("Best compression: $(best_compression.method) ($(best_compression.filesize_mb) MB)")
    println("Worst compression: $(worst_compression.method) ($(worst_compression.filesize_mb) MB)")
    println("Improvement: $(round(improvement; digits=1))% reduction from best to worst")
    
    # WebGraph comparison (CNR-2000 reference)
    webgraph_size = 1164843  # bytes
    webgraph_bpe = 2.897     # bits per edge
    
    println("\nWebGraph comparison:")
    for result in sorted_results
        gap = (result.filesize_bytes / webgraph_size - 1) * 100
        println("  $(result.method): $(round(gap; digits=1))% $(gap > 0 ? "larger" : "smaller") than WebGraph")
    end
    
    println("\n=== LOCALITY ANALYSIS ===")
    best_locality = minimum(r -> r.locality_score, results)
    println("Best locality: $(best_locality)")
    
    for result in results
        println("$(result.method): avg gap = $(round(result.avg_neighbor_gap; digits=1)), locality = $(round(result.locality_score; digits=3))")
    end
    
    println("\n✓ All relabeling methods tested successfully!")
end

"""
    create_lexicographic_relabeling(g)

Create lexicographic ordering based on out-degree (ascending order).
"""
function create_lexicographic_relabeling(g)
    T = eltype(g)
    n = nv(g)
    
    # Get out-degrees for all vertices
    degrees = [(T(v), length(outneighbors(g, v))) for v in 1:n]
    
    # Sort by out-degree (lexicographic order)
    sort!(degrees, by = x -> x[2])
    
    # Create mapping: old_vertex -> new_vertex
    relabel_map = Dict{T,T}()
    for (new_v, (old_v, _)) in enumerate(degrees)
        relabel_map[old_v] = T(new_v)
    end
    
    return relabel_graph(g, relabel_map)
end

"""
    create_llp_relabeling(g)

Create Local Layout Permutation ordering - a simplified version focusing on 
local neighborhood clustering.
"""
function create_llp_relabeling(g)
    T = eltype(g)
    n = nv(g)
    visited = falses(n)
    new_order = Int[]
    
    # Start with highest degree vertex
    degrees = [length(outneighbors(g, v)) for v in 1:n]
    start_vertex = argmax(degrees)
    
    # BFS-like traversal prioritizing high-degree unvisited neighbors
    queue = [start_vertex]
    visited[start_vertex] = true
    push!(new_order, start_vertex)
    
    while length(new_order) < n
        if isempty(queue)
            # Find next unvisited vertex with highest degree
            remaining = findall(!visited)
            if isempty(remaining)
                break
            end
            next_vertex = remaining[argmax([degrees[v] for v in remaining])]
            queue = [next_vertex]
        end
        
        current = popfirst!(queue)
        
        # Add unvisited neighbors, sorted by degree (descending)
        neighbors = collect(outneighbors(g, current))
        unvisited_neighbors = filter(v -> !visited[v], neighbors)
        sort!(unvisited_neighbors, by = v -> degrees[v], rev = true)
        
        for neighbor in unvisited_neighbors
            if !visited[neighbor]
                visited[neighbor] = true
                push!(new_order, neighbor)
                push!(queue, neighbor)
            end
        end
    end
    
    # Handle any remaining unvisited vertices
    for v in 1:n
        if !visited[v]
            push!(new_order, v)
        end
    end
    
    # Create relabeling map
    relabel_map = Dict{T,T}()
    for (new_v, old_v) in enumerate(new_order)
        relabel_map[T(old_v)] = T(new_v)
    end
    
    return relabel_graph(g, relabel_map)
end

"""
    create_minhash_relabeling(g)

Create MinHash-based ordering using Jaccard similarity of neighborhoods.
Groups similar vertices together.
"""
function create_minhash_relabeling(g)
    T = eltype(g)
    n = nv(g)
    
    # Simplified MinHash: use first k neighbors as "signature"
    k = 5  # signature size
    signatures = []
    
    for v in 1:n
        neighbors = sort(collect(outneighbors(g, v)))
        # Take first k neighbors as signature (pad with zeros if needed)
        signature = zeros(Int, k)
        for i in 1:min(k, length(neighbors))
            signature[i] = neighbors[i]
        end
        push!(signatures, (T(v), signature))
    end
    
    # Sort by signature (lexicographic comparison)
    sort!(signatures, by = x -> x[2])
    
    # Create relabeling map
    relabel_map = Dict{T,T}()
    for (new_v, (old_v, _)) in enumerate(signatures)
        relabel_map[T(old_v)] = T(new_v)
    end
    
    return relabel_graph(g, relabel_map)
end

"""
    calculate_locality_metrics(g)

Calculate locality metrics for the graph ordering.
"""
function calculate_locality_metrics(g)
    n = nv(g)
    gaps = Float64[]
    
    for v in 1:n
        neighbors = collect(outneighbors(g, v))
        if !isempty(neighbors)
            # Calculate gaps between vertex and its neighbors
            vertex_gaps = [abs(neighbor - v) for neighbor in neighbors]
            append!(gaps, vertex_gaps)
        end
    end
    
    avg_gap = isempty(gaps) ? 0.0 : mean(gaps)
    median_gap = isempty(gaps) ? 0.0 : median(gaps)
    
    # Locality score: lower is better (normalized by graph size)
    locality_score = avg_gap / n
    
    return (
        avg_neighbor_gap = avg_gap,
        median_neighbor_gap = median_gap,
        locality_score = locality_score
    )
end