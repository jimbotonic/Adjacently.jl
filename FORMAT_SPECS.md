# MGS Graph Compression Format Specification

## Overview

The MGS (Modified Graph Storage) format compresses directed graphs using multiple encoding techniques:
- **Delta encoding**: Compress sorted neighbor lists using differences
- **Hybrid mix-mode encoding**: Adaptively combine delta, run-length, and interval encoding
- **Reference encoding**: Reuse similar neighbor lists with bitmaps for differences
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
  [Traditional mix-mode encoded residuals] (if residuals_flag=1)
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
