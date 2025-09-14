# MGS Graph Compression Format Specification

## Overview

The MGS (Modified Graph Storage) format compresses directed graphs using multiple encoding techniques:
- **Delta encoding**: Compress sorted neighbor lists using differences
- **Hybrid mix-mode encoding**: Adaptively combine delta, run-length, and interval encoding
- **Reference encoding**: Reuse similar neighbor lists with bitmaps for differences
- **Recursive reference encoding**: Hierarchical reference compression with cost-based optimization
- **Multiple encodings**: Elias Gamma/Delta, Golomb, Fibonacci, Zeta, FED (Fibonacci+Elias Delta Hybrid)

## File Structure

```
[Header Section] [Data Section]
```

### Header Section
- Version identifier
- Graph metadata (vertices, edges)
- Encoding parameters

### Data Section Layout

#### 1. Configuration Flags
- **Reference flag** (1 bit): 0=hybrid mix only, 1=reference+hybrid mix enabled

#### 2. Index Section (`:index` mode only)
For each vertex:
- Out-degree encoded with specified encoding (shifted by +1 to avoid 0)

#### 3. Neighbor Lists

**For each vertex:**

##### With Reference Encoding Enabled:
```
children_flag (1 bit): 0=mix, 1=reference

If hybrid mix encoding (children_flag=0):
  [Hybrid mix-mode encoded neighbors]

If reference encoding (children_flag=1):
  ref_id (encoded): Reference vertex ID
  bitmap_len (encoded): Length of difference bitmap  
  bitmap (N bits): Bit vector of differences from reference
  residuals_flag (1 bit): 0=no residuals, 1=residuals follow
  
  If residuals_flag=1:
    recursive_flag (1 bit): 0=standard hybrid encoding, 1=recursive reference
    
    If recursive_flag=0:
      [Traditional mix-mode encoded residuals]
    
    If recursive_flag=1:
      [Recursive reference encoding - format repeats from children_flag]
```

##### Without Reference Encoding:
```
[Hybrid mix-mode encoded neighbors]
```

#### 4. Hybrid Mix-Mode Encoding

The hybrid approach adaptively combines delta, run-length, and interval encoding for optimal compression:

```
use_run_length_and_interval (1 bit): 0=delta only, 1=hybrid sections
first_value (encoding): First delta value

If delta only (flag=0):
  For each remaining value:
    delta_value (encoding): Delta-encoded gap
  
If hybrid sections (flag=1):
  num_sections (encoding): Number of encoding sections
  
  For each section:
    section_flag (1 bit): 0=delta section, 1=advanced section
    
    If delta section (section_flag=0):
      count (encoding): Number of delta values
      For each value:
        delta_value (encoding): Delta-encoded gap
    
    If advanced section (section_flag=1):
      second_flag (1 bit): 0=run-length, 1=interval
      
      If run-length (second_flag=0):
        count (encoding): Number of run-length pairs
        For each pair:
          gap_value (encoding): Delta gap
          run_length (encoding): Number of repetitions
      
      If interval (second_flag=1):
        count (encoding): Number of intervals
        For each interval:
          start_gap (encoding): Gap to interval start
          length (encoding): Interval length (adjusted by MIN_INTERVAL_LENGTH)
```

**Traditional Mix-Mode (for residuals):**
```
mix_mode_flag (1 bit): 0=delta only, 1=delta+run-length
first_neighbor (encoding): First neighbor value

If delta+run-length (mix_mode_flag=1):
  For each remaining neighbor:
    vertex_flag (1 bit): 0=delta, 1=run-length
    If delta: gap_value (encoding)
    If run-length: [gap_value (encoding)] [repeat_count (encoding)]
  
If delta only (mix_mode_flag=0):
  For each remaining neighbor:
    gap_value (encoding)
```

#### 5. Mode-Specific Details

**Index Mode (`:index`)**:
- No value shifting required
- No stop values needed
- Empty neighbor lists encoded implicitly via index

**Children Mode (`:children`)**:
- All gap values shifted by +1 before encoding (to avoid 0)
- Stop value written after each neighbor list
- Empty neighbor lists get stop value only

## Encoding Methods

| Method | Description | Use Case |
|--------|-------------|----------|
| `elias_gamma` | γ(n) = unary(⌊log₂(n)⌋+1) + binary(n-2^⌊log₂(n)⌋) | Small integers |
| `elias_delta` | δ(n) = γ(⌊log₂(n)⌋+1) + binary(n-2^⌊log₂(n)⌋) | Medium integers |
| `golomb` | Base-64 Golomb coding | Geometric distributions |
| `fibonacci` | Fibonacci sequence representation | Varied distributions |
| `zeta` | ζ_k(n) with parameter k=4 | Power-law distributions |
| `fed` | Block-based Fibonacci+Elias Delta hybrid | Large value ranges |

## Constants

```julia
REF_WINDOW_SIZE = 7           # Reference window size
MAX_REF_COUNT = 3             # Max references per vertex
MIN_INTERVAL_LENGTH = 4       # Min consecutive neighbors for intervals
FED_BLOCK_SIZE = 64           # FED encoding block size
REF_ENCODING_TH = 3           # Minimum overlap threshold for reference encoding
REF_V_MIN_DEGREE = 4          # Minimum degree for reference candidate vertices
```

## Hybrid Encoding Algorithm

The hybrid mix-mode encoding automatically selects the best encoding method for each section of delta values:

1. **Pattern Analysis**: Analyze delta sequences for:
   - Consecutive runs of identical values (run-length opportunities)
   - Consecutive integer sequences of length ≥ MIN_INTERVAL_LENGTH (interval opportunities)
   - Remaining irregular deltas (delta encoding)

2. **Adaptive Sectioning**: Group similar patterns into sections:
   - Delta sections: Irregular gap values
   - Run-length sections: Repeated gap values
   - Interval sections: Consecutive neighbor ranges

3. **Optimal Selection**: Choose encoding per section based on:
   - Delta: Standard gap encoding
   - Run-length: Gap + repetition count (efficient for social graphs)
   - Interval: Start + length (efficient for grid/mesh graphs)

This adaptive approach maximizes compression for diverse graph topologies while maintaining format compatibility.

## Recursive Reference Encoding

The recursive reference encoding system provides hierarchical compression by applying reference encoding recursively to residual neighbor lists. This approach is particularly effective for graphs with complex hierarchical similarity patterns.

### Algorithm Overview

1. **Primary Reference**: Encode vertex neighbors using standard reference encoding with a similar vertex
2. **Residual Analysis**: For remaining neighbors (residuals) after bitmap copy, perform cost analysis:
   - Cost of recursive reference encoding
   - Cost of standard hybrid mix encoding
3. **Recursive Decision**: If recursive reference is more efficient, apply the same reference encoding process to residuals
4. **Termination**: Recursion naturally terminates when no beneficial references exist

### Format Details

The recursive reference format extends the standard reference format by adding a `recursive_flag` after the `residuals_flag`:

```
Extended Reference Format:
children_flag (1 bit): 1 = reference mode
ref_id (encoded): Reference vertex ID  
bitmap_len (encoded): Length of copy bitmap
copy_bitmap (N bits): Bit vector indicating which neighbors to copy
residuals_flag (1 bit): 1 = residuals follow
IF residuals_flag = 1:
  recursive_flag (1 bit): 1 = recursive encoding, 0 = standard hybrid encoding
  IF recursive_flag = 0:
    [standard hybrid mix encoding of residuals]
  IF recursive_flag = 1:
    [recursive reference encoding of residuals - same format repeats]
```

### Cost-Based Optimization

The system uses intelligent cost estimation to decide between recursive and standard encoding:

```julia
function estimate_recursive_reference_cost(residuals, workspace)
    # Find best reference for residuals
    best_ref_vertex, overlap_count = find_best_reference_set(residuals, workspace)
    
    if overlap_count < REF_ENCODING_TH
        return typemax(Int)  # No viable reference
    end
    
    # Estimate costs
    ref_bits = estimate_reference_bits(residuals, best_ref_vertex, overlap_count)
    remaining_residuals = estimate_remaining_residuals(residuals, best_ref_vertex)
    recursive_cost = ref_bits + estimate_hybrid_cost(remaining_residuals)
    
    return recursive_cost
end
```

### Performance Characteristics

- **Compression Improvement**: 5-15% additional compression on graphs with hierarchical similarity
- **Processing Overhead**: Minimal due to cost-based early termination
- **Memory Efficiency**: Reuses existing workspace structures
- **Backward Compatibility**: Standard reference format remains unchanged

## Examples

### Example 1: Standard Reference Encoding

**Vertex 5 neighbors:** `[10, 15, 20, 25, 30]`
**Reference vertex 3 neighbors:** `[10, 15, 25, 35]`

**Encoding:**
```
children_flag: 1 (reference mode)
ref_id: encode(3)
bitmap_len: encode(4)
copy_bitmap: 1101 (copy positions 0,1,3 from reference)
residuals_flag: 1 (has residuals)
recursive_flag: 0 (use standard encoding)
[encode residuals [20, 30] with hybrid mix]
```

### Example 2: Recursive Reference Encoding

**Vertex 8 neighbors:** `[5, 10, 15, 20, 25, 30, 35, 40, 45]`
**Reference vertex 6 neighbors:** `[5, 10, 25, 30, 40]`
**Residuals after bitmap copy:** `[15, 20, 35, 45]`
**Secondary reference vertex 7 neighbors:** `[15, 20, 35, 50]`

**Encoding:**
```
children_flag: 1 (reference mode)
ref_id: encode(6)
bitmap_len: encode(5) 
copy_bitmap: 11011 (copy positions 0,1,2,4 from reference)
residuals_flag: 1 (has residuals)
recursive_flag: 1 (use recursive encoding)

  // Recursive encoding of residuals [15, 20, 35, 45]
  children_flag: 1 (reference mode)
  ref_id: encode(7)
  bitmap_len: encode(4)
  copy_bitmap: 1110 (copy positions 0,1,2 from secondary reference)
  residuals_flag: 1 (has final residuals)
  recursive_flag: 0 (use standard encoding for final residuals)
  [encode final residuals [45] with hybrid mix]
```

### Example 3: Cost-Based Decision Making

**Debug output showing cost-based optimization:**
```
┌ Debug: Recursive reference decision: cost_recursive=19, cost_hybrid=10, chosen=hybrid
└ @ Main.Adjacently.Compression ~/Documents/projects/Adjacently.jl/src/compression.jl:3167
```

In this case, the system calculated that standard hybrid encoding (10 bits) would be more efficient than recursive reference (19 bits), so it chose the standard approach automatically.
