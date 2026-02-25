# MGS Graph Compression Format Specification V2

## Overview

The MGS format compresses directed graphs using a **hierarchical cascade** of encoding techniques, from most to least efficient. Each vertex's neighbor list is encoded using the most efficient applicable method, selected by cost estimation.

## Encoding Hierarchy

### Level 1: Reference Encoding (Most Efficient)
- **When to use**: Similar neighbor lists exist within the reference window
- **Efficiency**: ~2-5 bits/edge
- **Components**: Reference vertex ID, copy bitmap, residuals

### Level 2: Recursive Reference Encoding
- **When to use**: Residuals from Level 1 can themselves be compressed via reference
- **Efficiency**: Additional 5-15% over standard reference
- **Components**: Primary reference encoding (Level 1) applied recursively on residuals; cost-based decision (only used if cheaper than Level 3)

### Level 3: Interval + Residual Encoding
- **When to use**: Neighbor list not using reference encoding
- **Efficiency**: ~2.5-4 bits/edge for intervals, delta for residuals
- **Components**: Intervals extracted directly from sorted neighbors as (start, length) pairs; remaining residuals delta-encoded

### Level 4: Hybrid Mix-Mode Encoding (Legacy)
- **When to use**: Residuals from reference encoding (when not using recursive reference)
- **Efficiency**: ~4-6 bits/edge
- **Components**: Delta sections (irregular gaps), run-length sections (repeated gaps), interval sections (consecutive sequences in delta space)

## Encoding Decision Tree

```
For each vertex's neighbor list:

1. Is there a good reference within the window?
   YES -> Use Reference Encoding (Level 1)
     |
     Are there residuals?
     YES -> Compare cost: Recursive Reference (Level 2) vs Interval+Residual (Level 3)
            Choose the cheaper option
   NO |

2. No reference available:
   -> Use Interval + Residual Encoding (Level 3)
      - Extract intervals (consecutive sequences >= MIN_INTERVAL_LENGTH)
      - Delta-encode remaining residuals
```

## File Structure

```
[Header Section] [Data Section]
```

### Header Section (12 bytes)

| Offset | Size | Field |
|--------|------|-------|
| 0 | 3 bytes | Signature: `MGS` (0x4d 0x47 0x53) |
| 3 | 1 byte | Major version (3) |
| 4 | 1 byte | Minor version (0) |
| 5 | 1 byte | Flags byte 1: graph_type (2 bits) &#124; coding_scheme (2 bits) &#124; integer_encoding (4 bits) |
| 6 | 1 byte | Flags byte 2: option_flags (8 bits) |
| 7 | 5 bytes | Number of vertices (little-endian UInt40) |

**Graph type**: 0b00 = directed, 0b01 = undirected

**Coding scheme**: 0b00 = children (stop-delimited), 0b01 = index (degree array)

**Integer encoding**: 0x1 = Elias gamma, 0x2 = Elias delta, 0x3 = Golomb, 0x4 = FED, 0x5 = Zeta, 0x6 = Fibonacci

**Option flags**:

| Value | Mode |
|-------|------|
| 0x00 | None (standard encoding) |
| 0x0F | ASTRA (full adaptive encoding) |
| 0x10–0x1F | RL policy mode (policy_id = value - 0x10 + 1, range 1–16) |
| 0x80 | Huffman (deprecated) |

### Data Section

#### Per-Vertex Encoding Format

```
encoding_mode (2 bits):
  00 = Pure Delta (Level 5)
  01 = Hybrid Mix (Level 4)
  10 = Interval+Residual (Level 3)
  11 = Reference (Level 1 or 2)

IF encoding_mode = 00 (Pure Delta):
  [delta-encoded neighbor list]

IF encoding_mode = 01 (Hybrid Mix):
  [hybrid mix-mode encoded list]

IF encoding_mode = 10 (Interval+Residual):
  num_intervals (encoded, +1 to avoid 0)
  For each interval:
    start_delta (encoded): delta from previous interval start
    length (encoded): interval length - MIN_INTERVAL_LENGTH + 1
  num_residuals (encoded, +1 to avoid 0)
  [delta-encoded residuals]

IF encoding_mode = 11 (Reference):
  ref_id (encoded): reference vertex ID
  bitmap_len (encoded): copy bitmap length
  copy_bitmap (N bits): which positions to copy from reference
  residuals_flag (1 bit): 0 = no residuals, 1 = has residuals

  IF residuals_flag = 1:
    residual_mode (2 bits):
      00 = Pure Delta
      01 = Hybrid Mix
      10 = Interval+Residual
      11 = Recursive Reference (repeat Reference format)
```

## Encoding Methods Performance

| Method | Typical Bits/Edge | Best For | Overhead |
|--------|-------------------|----------|----------|
| **Reference** | 2-5 | Similar neighbor patterns | Low (bitmap) |
| **Recursive Ref** | 2-4 | Hierarchical similarity | Medium (nested refs) |
| **Interval+Residual** | 2.5-4 | Consecutive node IDs (web graphs) | Low (count + pairs) |
| **Hybrid Mix** | 4-6 | Mixed patterns | Medium (section headers) |
| **Pure Delta** | 6-8 | Random/sparse graphs | Low (no overhead) |

## Parameters

| Parameter | Default | Description |
|-----------|---------|-------------|
| REF_WINDOW_SIZE | 1024 | Reference window size (sliding window of recent vertices) |
| MAX_REF_COUNT | 3 | Maximum references per vertex |
| REF_ENCODING_TH | 3 | Minimum overlap for reference encoding |
| REF_V_MIN_DEGREE | 4 | Minimum degree for a vertex to be a reference candidate |
| MIN_INTERVAL_LENGTH | 4 | Minimum consecutive neighbors to form an interval |
| ZETA_BASE | 3 | Zeta coding base parameter |
| FED_BLOCK_SIZE | 64 | FED block size |
| GOLOMB_BASE | 64 | Golomb coding base |

## Benchmark: CNR-2000

### WebGraph Reference
- **2.90 bits/edge**
- Intervals: 13.8% of arcs at 2.547 bits/interval-arc
- References: 68% of arcs copied
- Window size: 7, Zeta-3 encoding

### MGS V2 (Standard ASTRA)
- **5.108 bits/edge** (legacy ASTRA, Fibonacci encoding)

### MGS V2 + RL Greedy + LLP
- **5.093 bits/edge** (RL greedy mode with LLP vertex ordering)

### Key Gaps vs WebGraph
- **Ordering**: RCM vs LLP (WebGraph uses LLP for better locality)
- **Integer encoding**: Fibonacci vs Zeta-3 (WebGraph uses Zeta-3)
- **Copy bitmap compression**: RLE not yet applied to bitmaps
- **Window size**: 1024 vs 7 (larger window compensates for weaker ordering)

## RL-Compressed Stream Format

When `option_flags` is in the range 0x10–0x1F, the data section uses the RL (Reinforcement Learning) compressed format. The policy_id (1–16) is encoded as `option_flags - 0x10 + 1`.

### Stream Header (3 bits + 2 bits + 2 bits)

```
encoding_tag (3 bits):
  000 = Fibonacci
  001 = Zeta
  010 = Elias gamma
  011 = Elias delta

default_ref_mode (2 bits):
  00 = none (no reference encoding)
  01 = reference (standard reference)
  10 = recursive (recursive reference)

default_mil_tag (2 bits):
  00 = MIN_INTERVAL_LENGTH = 2
  01 = MIN_INTERVAL_LENGTH = 3
  10 = MIN_INTERVAL_LENGTH = 4
  11 = MIN_INTERVAL_LENGTH = 5
```

### Per-Vertex Record

Each vertex is encoded with a per-vertex header followed by its neighbor data:

```
per_vertex_header (4 bits):
  ref_mode (2 bits): override reference mode for this vertex
  mil_tag (2 bits): override MIN_INTERVAL_LENGTH for this vertex

neighbor_data:
  (same as standard V2 per-vertex format, using the per-vertex ref_mode and mil)
```

### Greedy Mode

When no trained policy file is provided, the encoder uses **greedy mode**: for each vertex, it exhaustively evaluates all combinations of reference mode and MIN_INTERVAL_LENGTH, then selects the combination yielding the fewest bits. The per-vertex header records the chosen settings so the decoder can reconstruct without the policy.

### Decoding

The RL format is self-describing — the policy is **not** needed for decompression. The per-vertex headers contain all encoding decisions. The `policy_id` in the file header is informational only.

## Backward Compatibility

Three strategies are supported:
1. **Version flag**: Format version in header (V1 = legacy, V2 = new hierarchy)
2. **Mode detection**: Decoder checks encoding_mode bits, falls back to hybrid mix
3. **Separate extension**: `.mgz` = current, `.mgz2` = new format
