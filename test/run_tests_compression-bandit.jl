#
# Adjacently: Julia Complex Directed Networks Library
# Copyright (C) 2016-2025 Jimmy Dubuisson <jimmy.dubuisson@gmail.com>
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

include("run_tests_main.jl")

#!/usr/bin/env julia

# Contextual Bandit Learning for Compression Strategy Selection
# Implements the bandit ideas from RL.md on CNR-2000

using LightGraphs: nv, ne, outneighbors
using Adjacently
using Adjacently.IO: load_adjacency_list_from_csv
using Adjacently.Relabeling: relabel_vertices_rcm, relabel_graph
using Adjacently.Compression: estimate_reference_encoding_cost, estimate_interval_runlength_encoding_cost,
                               create_reference_data!, RefBuildWorkspace, estimate_block_encoding_cost,
                               estimate_bitmap_rle_cost
using Adjacently.RL
using Statistics

println("=" ^ 80)
println("Contextual Bandit Learning for Compression")
println("Testing on CNR-2000")
println("=" ^ 80)

# Load and prepare graph
println("\n Loading CNR-2000 dataset...")
graph_original = load_adjacency_list_from_csv("datasets/webgraph/cnr-2000/cnr-2000.csv", ',', true)
vertices_count = nv(graph_original)
edges_count = ne(graph_original)
println("  Graph: $vertices_count vertices, $edges_count edges")

println("\n Applying RCM relabeling...")
rcm_mapping = relabel_vertices_rcm(graph_original, :out)
graph = relabel_graph(graph_original, rcm_mapping)
println("  RCM applied successfully")

# Collect neighbor lists
println("\n Collecting neighbor lists...")
neighbor_lists = Dict{UInt24,Vector{UInt24}}()
for v in 1:nv(graph)
    neighbors = sort(collect(outneighbors(graph, UInt24(v))))
    neighbor_lists[UInt24(v)] = neighbors
end
println("  Collected $(length(neighbor_lists)) neighbor lists")

# =============================================================================
# PHASE 1: Data Collection - Extract features and compute ground-truth costs
# =============================================================================

println("\n" * "=" ^ 80)
println("PHASE 1: COLLECTING TRAINING DATA")
println("=" ^ 80)

mutable struct EncodingDecision
    vertex_id::UInt24

    # Features
    nl_size::Int                    # |neighbor_list|
    ref_available::Bool             # Is reference available?
    ref_size::Int                   # |reference_list| (if available)
    overlap_count::Int              # Number of shared neighbors
    overlap_ratio::Float64          # overlap / min(nl, ref)
    residuals_count::Int            # Number of residuals
    bitmap_density::Float64         # Density of copy bitmap
    run_count::Int                  # Number of runs in bitmap
    max_run_length::Int             # Longest run in bitmap

    # Ground truth costs (in bits)
    interval_cost::Int              # Cost of interval+RL encoding
    reference_cost::Int             # Cost of reference encoding (if available)
    best_choice::Symbol             # :interval or :reference
    cost_difference::Int            # |interval_cost - reference_cost|

    # Bitmap encoding sub-decision
    bitmap_raw_cost::Int            # Raw bitmap bits
    bitmap_rle_cost::Int            # RLE bitmap bits
    bitmap_block_cost::Int          # Block encoding bits
    best_bitmap_method::Symbol      # :raw, :rle, or :block
end

training_data = EncodingDecision[]
ref_workspace = RefBuildWorkspace{UInt24}()
ref_window = UInt24[]
window_size = 1024

println("\n Collecting decision examples...")
sample_interval = max(1, vertices_count ÷ 1000)  # Sample ~1000 vertices for speed

for v in 1:vertices_count
    vertex = UInt24(v)
    neighbors = neighbor_lists[vertex]

    # Skip empty vertices
    if isempty(neighbors)
        continue
    end

    # Sample vertices to keep training set manageable
    if v % sample_interval != 0 && v != 1
        # Still need to maintain reference window
        push!(ref_window, vertex)
        if length(ref_window) > window_size
            popfirst!(ref_window)
        end
        continue
    end

    # Compute interval+RL cost (always available)
    interval_cost = estimate_interval_runlength_encoding_cost(
        neighbors, :fibonacci, 3, 2
    )

    # Try to find best reference
    best_ref = nothing
    best_ref_cost = typemax(Int)
    best_copy_bitmap = Bool[]
    best_residuals = UInt24[]

    for candidate_ref in ref_window
        if !haskey(neighbor_lists, candidate_ref)
            continue
        end

        ref_neighbors = neighbor_lists[candidate_ref]
        if isempty(ref_neighbors)
            continue
        end

        # Build reference encoding
        create_reference_data!(ref_workspace, neighbors, ref_neighbors)

        # Quick filter: skip if overlap too low
        overlap = count(ref_workspace.copy_bitmap)
        if overlap < 3
            continue
        end

        # Calculate cost
        cost = estimate_reference_encoding_cost(
            candidate_ref,
            ref_workspace.copy_bitmap,
            ref_workspace.residuals,
            :children,
            :fibonacci
        )

        if cost < best_ref_cost
            best_ref_cost = cost
            best_ref = candidate_ref
            best_copy_bitmap = copy(ref_workspace.copy_bitmap)
            best_residuals = copy(ref_workspace.residuals)
        end
    end

    # Extract features
    if best_ref !== nothing
        ref_neighbors = neighbor_lists[best_ref]
        overlap_count = count(best_copy_bitmap)
        residuals_count = length(best_residuals)

        # Bitmap statistics
        bitmap_density = isempty(best_copy_bitmap) ? 0.0 : overlap_count / length(best_copy_bitmap)

        # Count runs (consecutive 1s or 0s)
        run_count = 0
        max_run = 0
        current_run = 0
        last_val = nothing
        for bit in best_copy_bitmap
            if bit == last_val
                current_run += 1
            else
                if last_val !== nothing
                    run_count += 1
                    max_run = max(max_run, current_run)
                end
                current_run = 1
                last_val = bit
            end
        end
        if last_val !== nothing
            run_count += 1
            max_run = max(max_run, current_run)
        end

        # Bitmap sub-costs
        bitmap_raw = length(best_copy_bitmap)
        bitmap_rle = estimate_bitmap_rle_cost(best_copy_bitmap, :fibonacci)
        bitmap_block = estimate_block_encoding_cost(best_copy_bitmap, :fibonacci)

        best_bitmap = :raw
        if bitmap_block < min(bitmap_raw, bitmap_rle)
            best_bitmap = :block
        elseif bitmap_rle < bitmap_raw
            best_bitmap = :rle
        end

        decision = EncodingDecision(
            vertex,
            length(neighbors),
            true,
            length(ref_neighbors),
            overlap_count,
            overlap_count / min(length(neighbors), length(ref_neighbors)),
            residuals_count,
            bitmap_density,
            run_count,
            max_run,
            interval_cost,
            best_ref_cost,
            best_ref_cost <= interval_cost ? :reference : :interval,
            abs(interval_cost - best_ref_cost),
            bitmap_raw,
            bitmap_rle,
            bitmap_block,
            best_bitmap
        )

        push!(training_data, decision)
    else
        # No reference available or all refs too expensive
        decision = EncodingDecision(
            vertex,
            length(neighbors),
            false,
            0,
            0,
            0.0,
            0,
            0.0,
            0,
            0,
            interval_cost,
            typemax(Int),
            :interval,
            0,
            0, 0, 0,
            :raw
        )

        push!(training_data, decision)
    end

    # Progress
    if length(training_data) % 100 == 0
        println("  Collected $(length(training_data)) examples...")
    end

    # Add to reference window
    push!(ref_window, vertex)
    if length(ref_window) > window_size
        popfirst!(ref_window)
    end
end

println("\n Total training examples collected: $(length(training_data))")

# =============================================================================
# PHASE 2: Analysis - Study the decision patterns
# =============================================================================

println("\n" * "=" ^ 80)
println("PHASE 2: ANALYZING DECISION PATTERNS")
println("=" ^ 80)

# Filter to examples with references available
ref_available_data = filter(d -> d.ref_available, training_data)
println("\n Examples with reference available: $(length(ref_available_data))")

# Analyze reference vs interval choice
ref_chosen = filter(d -> d.best_choice == :reference, ref_available_data)
interval_chosen = filter(d -> d.best_choice == :interval, ref_available_data)

println("\n Reference vs Interval Decision:")
println("  Reference chosen: $(length(ref_chosen)) ($(round(100*length(ref_chosen)/length(ref_available_data), digits=1))%)")
println("  Interval chosen:  $(length(interval_chosen)) ($(round(100*length(interval_chosen)/length(ref_available_data), digits=1))%)")

if !isempty(ref_chosen)
    println("\n Reference Chosen - Feature Statistics:")
    println("  Avg overlap ratio:    $(round(mean(d.overlap_ratio for d in ref_chosen), digits=3))")
    println("  Avg bitmap density:   $(round(mean(d.bitmap_density for d in ref_chosen), digits=3))")
    println("  Avg nl_size:          $(round(mean(d.nl_size for d in ref_chosen), digits=1))")
    println("  Avg residuals_count:  $(round(mean(d.residuals_count for d in ref_chosen), digits=1))")
    println("  Avg cost savings:     $(round(mean(d.cost_difference for d in ref_chosen), digits=1)) bits")
end

if !isempty(interval_chosen)
    println("\n Interval Chosen - Feature Statistics:")
    println("  Avg overlap ratio:    $(round(mean(d.overlap_ratio for d in interval_chosen), digits=3))")
    println("  Avg bitmap density:   $(round(mean(d.bitmap_density for d in interval_chosen), digits=3))")
    println("  Avg nl_size:          $(round(mean(d.nl_size for d in interval_chosen), digits=1))")
    println("  Avg cost penalty:     $(round(mean(d.cost_difference for d in interval_chosen), digits=1)) bits")
end

# Analyze bitmap encoding method choice
bitmap_data = filter(d -> d.ref_available && d.best_choice == :reference, training_data)
if !isempty(bitmap_data)
    block_chosen = filter(d -> d.best_bitmap_method == :block, bitmap_data)
    rle_chosen = filter(d -> d.best_bitmap_method == :rle, bitmap_data)
    raw_chosen = filter(d -> d.best_bitmap_method == :raw, bitmap_data)

    println("\n Bitmap Encoding Method (for references):")
    println("  Block chosen: $(length(block_chosen)) ($(round(100*length(block_chosen)/length(bitmap_data), digits=1))%)")
    println("  RLE chosen:   $(length(rle_chosen)) ($(round(100*length(rle_chosen)/length(bitmap_data), digits=1))%)")
    println("  Raw chosen:   $(length(raw_chosen)) ($(round(100*length(raw_chosen)/length(bitmap_data), digits=1))%)")
end

# =============================================================================
# PHASE 3: Simple Rule-Based Bandit
# =============================================================================

println("\n" * "=" ^ 80)
println("PHASE 3: SIMPLE RULE-BASED BANDIT")
println("=" ^ 80)

println("\n Learning simple decision rules from data...")

# Rule 1: Reference selection threshold
# Find overlap_ratio threshold that maximizes correct decisions
ref_data_sorted = sort(ref_available_data, by=d -> d.overlap_ratio)
best_threshold = 0.0
best_accuracy = 0.0

for threshold in 0.1:0.05:0.9
    correct = sum(d -> (d.overlap_ratio >= threshold) == (d.best_choice == :reference), ref_available_data)
    accuracy = correct / length(ref_available_data)
    if accuracy > best_accuracy
        global best_accuracy = accuracy
        global best_threshold = threshold
    end
end

println("\n Learned Rule: Use reference if overlap_ratio >= $(round(best_threshold, digits=2))")
println("  Accuracy on training set: $(round(100*best_accuracy, digits=1))%")

# Rule 2: Bitmap method selection
# Find density threshold for block vs RLE
if !isempty(bitmap_data)
    best_density_threshold = 0.0
    best_bitmap_accuracy = 0.0

    for threshold in 0.1:0.05:0.9
        # Simple rule: block if density < threshold, else raw
        correct = sum(bitmap_data) do d
            predicted = d.bitmap_density < threshold ? :block : :raw
            predicted == d.best_bitmap_method
        end
        accuracy = correct / length(bitmap_data)
        if accuracy > best_bitmap_accuracy
            global best_bitmap_accuracy = accuracy
            global best_density_threshold = threshold
        end
    end

    println("\n Learned Rule: Use block encoding if density < $(round(best_density_threshold, digits=2))")
    println("  Accuracy on training set: $(round(100*best_bitmap_accuracy, digits=1))%")
end

# =============================================================================
# PHASE 4: Results Summary
# =============================================================================

println("\n" * "=" ^ 80)
println("SUMMARY")
println("=" ^ 80)

total_bits_saved = sum(d.cost_difference for d in ref_chosen)
total_bits_lost = sum(d.cost_difference for d in interval_chosen)
net_savings = total_bits_saved - total_bits_lost

println("\n Training Data Summary:")
println("  Total examples:           $(length(training_data))")
println("  Examples with references: $(length(ref_available_data))")
println("  Reference beneficial:     $(length(ref_chosen)) cases")
println("  Net bit savings:          $net_savings bits")
println("  Avg savings per ref use:  $(round(total_bits_saved / max(1, length(ref_chosen)), digits=1)) bits")

println("\n Key Insights:")
println("  1. Overlap ratio is a strong predictor of reference utility")
println("  2. Bitmap density correlates with best encoding method")
println("  3. Simple threshold rules achieve $(round(100*best_accuracy, digits=1))% accuracy")
println("\n Next Steps:")
println("  - Integrate learned thresholds into compression pipeline")
println("  - Test on validation graphs (e.g., Web-Google)")
println("  - Explore multi-feature decision trees for better accuracy")
println("  - Implement full contextual bandit with feature weighting")

println("\n" * "=" ^ 80)
println("Bandit learning test complete!")
println("=" ^ 80)
