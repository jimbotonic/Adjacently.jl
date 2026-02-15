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

"""
    ASTRA (Adaptive Streaming Adjacency) Compression

This module provides documentation and organization for the ASTRA compression algorithm.
All ASTRA functions are implemented in the parent Compression module and re-exported here.

ASTRA is a graph compression method inspired by WebGraph that uses:
- Greedy cost-based reference selection
- Adaptive bitmap encoding (block vs raw)
- Recursive reference compression
- Variable-length integer encoding (Fibonacci, Elias-Gamma/Delta, etc.)
- Mix encoding (intervals + run-length) for residuals

Key improvements over WebGraph:
- Greedy cost minimization instead of overlap heuristics
- Adaptive bitmap selection based on actual bit costs
- Simplified format (removed RLE, using 1-bit flag)

Performance on CNR-2000:
- 5.108 bits per edge
- 1.57 edges per byte compression ratio
- Perfect round-trip fidelity

## ASTRA-Specific Functions (defined in parent Compression module)

### Main Compression/Decompression:
- `write_compressed_graph_data()` - Main compression entry point with ASTRA algorithm
- `read_compressed_graph_data()` - Main decompression entry point

### Reference Selection (WebGraph-inspired):
- `find_best_reference_greedy_cost()` - Greedy cost-based reference selection (main ASTRA algorithm)
- `find_best_reference()` - Original WebGraph overlap-based selection
- `find_best_reference_fast()` - Fast reference lookup using index
- `find_best_reference_set()` - Reference selection from candidate set

### Reference Data Building:
- `create_reference_data()` - Build reference encoding (copy bitmap + residuals)
- `create_reference_data!()` - In-place version using workspace
- `reconstruct_from_reference()` - Reconstruct neighbors from reference
- `RefBuildWorkspace` - Workspace type for efficient reference building

### Recursive References:
- `write_recursive_reference_residuals()` - Recursive residual encoding
- `read_recursive_reference_residuals()` - Recursive residual decoding

### Adaptive Bitmap Encoding:
- `write_bitmap_adaptive()` - Adaptive bitmap compression (chooses block vs raw)
- `read_bitmap_adaptive()` - Adaptive bitmap decompression
- `write_block_encoding()` - Block encoding for sparse bitmaps
- `read_block_encoding()` - Block decoding
- `estimate_block_encoding_cost()` - Cost estimation for block encoding

### Mix Encoding (Intervals + Run-Length):
- `write_mix_encoded_list()` - Mix-mode encoding
- `read_mix_encoded_list()` - Mix-mode decoding
- `write_hybrid_mix_encoded_list()` - Advanced hybrid mix encoding
- `read_hybrid_mix_encoded_list()` - Hybrid mix decoding
- `write_intervals_and_residuals()` - Interval compression with residuals
- `read_intervals_and_residuals()` - Interval decompression

### Pattern Analysis:
- `analyze_delta_patterns_hybrid()` - Analyze delta patterns for hybrid encoding
- `find_run_length_patterns()` - Find run-length patterns
- `find_consecutive_length()` - Find consecutive sequence length
- `count_consecutive()` - Count consecutive values

### Cost Estimation:
- `estimate_reference_encoding_cost()` - Estimate reference encoding cost
- `estimate_recursive_reference_cost()` - Estimate recursive reference cost
- `estimate_hybrid_mix_encoding_cost()` - Estimate hybrid mix cost
- `estimate_interval_runlength_encoding_cost()` - Estimate interval+RL cost
- `estimate_mix_encoding_cost()` - Estimate mix encoding cost
- `estimate_bitmap_rle_cost()` - Estimate bitmap RLE cost

### Helper Functions:
- `reconstruct_from_delta()` - Reconstruct values from deltas
- `_consume_children_trailing_stop()` - Consume trailing stop marker

## Format Details

- Uses `OPTION_ASTRA` flag in MGS format
- Delta encoding for sorted adjacency lists
- Mix encoding (intervals + run-length) for residuals
- Greedy reference selection with adaptive bitmap compression
- Recursive reference support for multi-level chains

## See Also

- `ASTRALayered` - Level-based variant using BFS decomposition
- `Compression` - Parent module containing all implementations
"""
module ASTRA

# This is a documentation module.
# All ASTRA functions are implemented in the parent Compression module
# and accessible via Compression.function_name

end # module ASTRA
