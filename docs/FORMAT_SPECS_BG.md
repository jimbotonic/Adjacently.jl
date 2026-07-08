# MGS Graph Compression Format Specification V3

## Overview

> **Context-Range backend:** in addition to the prefix/varint integer codes below,
> BG supports `integer_encoding=:context_range` (header `0x7`) — a context-adaptive
> range coder with a 3-stream layout, copy-aware rank gaps, and chunked random
> access. This document specifies the classic (prefix-code) format; the shared
> context-range format is in [FORMAT_SPECS_CONTEXT_RANGE.md](FORMAT_SPECS_CONTEXT_RANGE.md).

MGS v3 extends v2 with two changes:

1. **Zeta-3 integer encoding** (z3), matching WebGraph's encoding choice.
2. **VLC vertex headers** — variable-length prefix codes replace the fixed 6-bit per-vertex headers, saving ~3.7 bits/vertex (~0.37 BPE on CNR-2000).

### Key changes from V2

| Aspect | V2 | V3 |
|--------|----|----|
| Integer encoding | Fibonacci (~1.44*log2(n) + 2 bits) | Zeta-3 (~4 bits for values 1-7) |
| Stream header tag | `000` (Fibonacci) | `001` (Zeta) |
| Per-vertex header | Fixed 6 bits (ref_mode 2b + enc_opt 4b) | VLC: 1-7 bits (avg ~1.5 bits) |

## Integer Encoding: Zeta-3

Zeta coding with parameter k=3 encodes a positive integer v as:

1. Compute h = floor(log2(v) / k)
2. Write h+1 in **unary** (h zeros followed by a 1): h+1 bits
3. Write the remainder in **truncated binary**: k*(h+1) bits

**Total bits**: (h+1) + k*(h+1) = (h+1)*(k+1)

| Value range | h | Bits |
|-------------|---|------|
| 1-7 | 0 | 4 |
| 8-63 | 1 | 8 |
| 64-511 | 2 | 12 |
| 512-4095 | 3 | 16 |

For web graph gaps (typically 1-100), most values cost 4-8 bits vs Fibonacci's 5-10 bits.

## File Structure

```
[Header Section (12 bytes)] [Bitstream Data Section]
```

### Header Section (12 bytes)

| Offset | Size | Field |
|--------|------|-------|
| 0 | 3 bytes | Signature: `MGS` (0x4d 0x47 0x53) |
| 3 | 1 byte | Major version (3) |
| 4 | 1 byte | Minor version (2 = self-describing v3.2) |
| 5 | 1 byte | Flags byte 1: graph_type (2 bits) &#124; coding_scheme (2 bits) &#124; integer_encoding (4 bits) |
| 6 | 1 byte | Flags byte 2: algorithm + encoded params (8 bits) — see below |
| 7 | 5 bytes | Number of vertices (little-endian UInt40) |

**Graph type**: 0b00 = directed, 0b01 = undirected

**Coding scheme**: 0b00 = children (stop-delimited), 0b01 = index (degree array / CG offset table)

**Integer encoding (4 bits)**:

| Code | Encoding |
|------|----------|
| 0x1 | Elias gamma |
| 0x2 | Elias delta |
| 0x3 | Golomb |
| 0x4 | FED |
| 0x5 | **Zeta** (k=3) |
| 0x6 | Fibonacci |

**Flags byte 2 — algorithm + params (v3.2)**:

| Range | Mode |
|-------|------|
| 0x00-0x0F | Algorithm IDs (hardcoded defaults) |
| 0x10-0x4F | BG + encoded params (6 bits) |
| 0x50-0x6F | CS + encoded params (5 bits) |
| 0x70-0xFF | CG + encoded params (mixed-radix, 144 slots) |

See `MGS_HEADER.md` for full byte 2 encoding details for CS and CG modes.

### BG Byte 2 Encoding (0x10-0x4F)

Offset = byte2 - 0x10. Bit layout (6 bits):

| Bits | Field | Values |
|------|-------|--------|
| 5-4 | ref_window_size | 0=8, 1=16, 2=32, 3=64 |
| 3 | copy_blocks | 0=false, 1=true |
| 2 | stop_deltas | 0=false, 1=true |
| 1 | lr_split | 0=false, 1=true |
| 0 | adaptive_header | 0=false, 1=true (per-vertex adaptive enc_type selection) |

**Implied parameters** (not encoded, always active):
- `adaptive_copy = copy_blocks` — 3-way adaptive copy is always used when copy_blocks is enabled
- `compact_copy = true` — compact 3-way prefix code always used
- `tight_intervals = true` — tight interval encoding always used
- `fixwidth_ref = lr_split` — fixed-width reference distances when LR-split is active
- `multi_ref` — not encoded in header; VLC codes are self-describing (decoder reads multi-ref VLC prefixes directly)

**Example**: Best config (window=64, copy_blocks, stop_deltas, lr_split, adaptive_header=false) encodes as:
offset = (3 << 4) | 0x08 | 0x04 | 0x02 | 0x00 = 0x3E, byte2 = 0x10 + 0x3E = 0x4E.

## Greedy Compressed Data Section

### Stream Header (4 bits)

```
coding_scheme_flag (1 bit):
  0 = :children (stop-delimited)
  1 = :index (degree array precedes data)

encoding_tag (3 bits):
  000 = Fibonacci
  001 = Zeta
  010 = Elias gamma
  011 = Elias delta
```

The encoding tag makes the stream **self-describing** — the decoder reads it and uses the corresponding integer encoding for all subsequent values. No policy file is needed for decompression.

### Optional Index Section (Per-Vertex Offset Table)

When `coding_scheme_flag = 1` (`:index` mode), a per-vertex offset table replaces
the old degree array, enabling O(1) random access to any vertex's compressed data:

```
entry_width (6 bits): number of bits per offset entry
N+1 entries (each entry_width bits):
  entry[0..N-1]: bit offset to vertex i's data (relative to start of per-vertex section)
  entry[N]:      total bit size of per-vertex section (end marker)
```

Empty vertex detection: if `entry[v] == entry[v+1]`, vertex v has no neighbors.

**Index mode encoding changes** (BG):
- `stop_deltas` forced to `false` — redundant with offset table boundaries
- `adaptive_deltas` forced to `false`
- Empty vertices are detected via `offset[v] == offset[v+1]` (no VLC empty code needed)
- All neighbor lists use count-prefixed delta encoding

**Two-pass encoding**: The encoder first writes all per-vertex data to a temporary
buffer, recording each vertex's start bit offset. Then `entry_width = ceil(log2(total_bits + 1))`
is computed, and the offset table + buffered data are written to the main stream.

**Random access**: To read vertex v, seek to `offset[v]` in the per-vertex section.
If vertex v references vertex v-d, recursively decode v-d (reference chains bounded
by `ref_window_size`, worst case: decode `ref_window_size` consecutive vertices).

### Per-Vertex Record

Each vertex is encoded with a **VLC (variable-length coded) header** followed by payload.

#### VLC Header (Merged)

The VLC header encodes three fields per vertex: `ref_mode` (none, reference, multi_ref, or empty), `enc_type` (delta, interval, or rle), and `mil` (minimum interval length, for interval/rle modes).

Empty vertices (degree-0) are encoded directly as a 3-bit VLC code (`110`), eliminating the need for a separate 1-bit empty-prefix flag. This saves 1 bit for the ~97% of non-empty vertices that previously paid for the empty-prefix flag, at the cost of 2 extra bits for the ~3% empty vertices. Net savings: ~0.08 BPE on CNR-2000.

All 28 action combinations (including empty) are covered without escape codes. The two most common actions (ref+delta, none+delta) receive the shortest codes. Interval and RLE actions include a 2-bit MIL tag, giving all MIL values equal cost. Multi-ref codes use longer prefixes since they are rare but high-value.

```
  0            (1 bit)   = reference + delta
  10           (2 bits)  = none      + delta
  110          (3 bits)  = empty (no neighbors)
  1110+mm      (6 bits)  = none      + interval, MIL from mm
  11110+mm     (7 bits)  = reference + interval, MIL from mm
  111110+mm    (8 bits)  = none      + rle, MIL from mm
  1111110+mm   (9 bits)  = reference + rle, MIL from mm
  11111110     (8 bits)  = multi-ref + delta
  111111110+mm (11 bits) = multi-ref + interval, MIL from mm
  111111111+mm (11 bits) = multi-ref + rle, MIL from mm
```

**MIL tag (mm, 2 bits)**:
```
  00 = MIL 2
  01 = MIL 3
  10 = MIL 4
  11 = MIL 5
```

**VLC header bit costs**:

| Action | Bits | Code |
|--------|------|------|
| reference + delta | 1 | `0` |
| none + delta | 2 | `10` |
| empty | 3 | `110` |
| none + interval | 6 | `1110` + 2-bit MIL |
| reference + interval | 7 | `11110` + 2-bit MIL |
| none + rle | 8 | `111110` + 2-bit MIL |
| reference + rle | 9 | `1111110` + 2-bit MIL |
| multi-ref + delta | 8 | `11111110` |
| multi-ref + interval | 11 | `111111110` + 2-bit MIL |
| multi-ref + rle | 11 | `111111111` + 2-bit MIL |

**VLC header cost in greedy search**: The `_vlc_header_cost` function returns the exact bit cost for each `(ref_mode, enc_type, mil)` combination. The greedy search includes this cost in all comparisons, ensuring the encoder accounts for the variable header overhead when selecting encoding strategies.

#### Historical VLC v1 and VLC v2

Earlier versions used two VLC schemes (v1 and v2) that only had shortcodes for MIL=2 actions and required a 9-10 bit escape path for all other combinations:

**VLC v1** (original):
```
  00   (2 bits) = none      + delta
  01   (2 bits) = none      + interval, MIL=2
  10   (2 bits) = reference + delta
  110  (3 bits) = reference + interval, MIL=2
  111  (3 bits) = escape -> full 6-bit header follows:
                    ref_mode (2 bits) + enc_opt (4 bits)
```

**VLC v2** (optimized for reference-heavy streams):
```
  0    (1 bit)  = reference + delta
  10   (2 bits) = none      + delta
  110  (3 bits) = none      + interval, MIL=2
  1110 (4 bits) = reference + interval, MIL=2
  1111 (4 bits) = escape -> full 6-bit header follows:
                    ref_mode (2 bits) + enc_opt (4 bits)
```

The current merged VLC scheme supersedes both by eliminating the expensive escape path (9-10 bits), integrating the empty vertex code directly into the prefix tree (removing the separate 1-bit empty-prefix flag), and giving all MIL values equal cost (6-11 bits for interval/rle). This makes the greedy search unbiased with respect to MIL choice and improves decoding speed (no conditional escape parsing).

### Per-Vertex Payload

#### When ref_mode = none (00)

The full neighbor list is encoded using the selected encoding option:

**Interval+Residual** (standard, `lr_split=false`):
```
varint(num_intervals + 1)
For each interval:
  varint(start_delta)               // zigzag(start - vertex_id) + 1 for first; delta from prev for subsequent
  varint(length - MIL + 1)          // length shifted by MIL
varint(num_residuals + 1)
delta_encoded_residuals             // first value: zigzag(v - vertex_id) + 1, then gaps+1
```

**Interval+LR-split Residual** (`lr_split=true`):
```
varint(num_intervals + 1)
For each interval:
  varint(start_delta)               // zigzag(start - vertex_id) + 1 for first; delta from prev for subsequent
  varint(length - MIL + 1)          // length shifted by MIL
varint(num_residuals + 1)
IF num_residuals > 0:
  varint(n_left + 1)               // count of residuals < vertex_id (decoder derives n_right = total - n_left)
  delta_encoded(left_distances)    // vid - residual[n_left-i+1] for i=1..n_left (ascending distances from vid)
  delta_encoded(right_distances)   // residual[i] - vid + 1 for i=1..n_right (ascending distances from vid)
```

LR-split separates residuals at `vertex_id` into left (< vid) and right (>= vid) halves. Each half is transformed to ascending distances from `vertex_id` and delta-encoded. This produces smaller first values than zigzag encoding, saving ~0.04 BPE on CNR-2000 with intervals enabled.

When `tight_intervals=true` (always active in v3.2), subsequent interval starts are delta-encoded from the **end** of the previous interval (start + length) rather than from the start.

**RLE/Hybrid Mix**:
```
varint(count + 1)                   // neighbor count
hybrid_mix_encoded_deltas           // delta-encoded then hybrid-compressed
  1 bit: hybrid_active flag
  varint(first_delta)
  if hybrid:
    varint(num_sections)
    per section:
      0 + varint(count) + values            // delta section
      1,0 + varint(num_pairs) + pairs       // run-length section
      1,1 + varint(num_pairs) + pairs       // interval section
```

**Delta**:

Two sub-formats depending on `stop_deltas` configuration:

*Count-prefixed delta* (`stop_deltas=false`):
```
varint(count + 1)                   // neighbor count
varint(zigzag(first - vertex_id) + 1)   // first neighbor
varint(gap2 + 1) ... varint(gapN + 1)  // subsequent gaps shifted by +1
```

*STOP-terminated delta* (`stop_deltas=true`):
```
For each value in sorted order:
  bit '1'                           // continuation flag
  varint(gap)                       // zigzag(first - vertex_id) + 1 for first; gap for rest
bit '0'                             // STOP terminator
```

STOP-terminated deltas eliminate the count prefix (2-6 bits), which saves bits for low-degree vertices and costs slightly more for high-degree vertices. Net savings on CNR-2000: ~0.08 BPE.

#### When ref_mode = reference (single reference)

```
varint(ref_distance)                // distance back in reference window
copy_encoding                       // which ref neighbors to copy (see below)
residual_payload                    // residuals encoded per enc_opt above
```

#### When ref_mode = multi_ref (two references)

When `multi_ref=true` is encoded in the header, the greedy encoder may select two references for a single vertex. The encoder picks the top-K (K=5) reference candidates by overlap count and tries each as ref1 paired with all window candidates as ref2, selecting the pair that minimizes total cost.

```
varint(ref1_distance)               // distance to first reference
copy_encoding_1                     // which ref1 neighbors to copy
varint(ref2_distance)               // distance to second reference
copy_encoding_2                     // which ref2 neighbors to copy (only neighbors not in ref1's copied set)
residual_payload                    // remaining neighbors encoded per enc_opt
```

The decoder reconstructs neighbors as the union of copied1, copied2, and residuals. The encoder ensures copy_encoding_2 only marks neighbors in N(v) that were not already copied from ref1, so no duplicates occur.

Multi-ref saves ~0.02 BPE on CNR-2000 by reducing residual lists for vertices whose neighbors span two different reference neighborhoods.

### Copy Position Encoding

Three modes available, controlled by `copy_blocks` and `adaptive_copy` parameters:

#### Adaptive Bitmap (default, `copy_blocks=false`)

```
small_count(bitmap_length):
  00 = 0 elements
  01 = 1 element
  10 = 2 elements
  11 + varint(length) = 3+ elements

if length > 0:
  format_flag (1 bit):
    0 = block encoding
    1 = raw bits

Block encoding:
  varint(num_blocks + 1)
  For each block:
    varint(block_length + 1)

Raw encoding:
  N raw bits (one per bitmap entry)
```

#### Copy-Blocks (`copy_blocks=true`, `adaptive_copy=false`)

WebGraph-style (start, length) run encoding:

```
small_count(num_copy_blocks)       // number of contiguous runs to copy
IF num_copy_blocks > 0:
  varint(first_start)              // 1-based index of first copied position
  varint(first_length)             // length of first run
  For each subsequent block:
    varint(gap)                    // positions skipped since end of last block
    varint(length)                 // length of this run
```

#### 3-Way Adaptive Copy (`adaptive_copy=true`, `copy_blocks=true`)

Per-reference vertex, picks cheapest of three modes using nested bit flags.

**Compact prefix code** (`compact_copy=true`, always active in v3.2):
```
  0  -> complement copy-blocks (1 bit overhead)
  10 -> copy-blocks (2 bits overhead)
  11 -> raw bitmap (2 bits overhead)
```

The compact prefix assigns the cheapest code (1 bit) to complement mode, which is the most common choice for high-overlap references. This saves ~0.01 BPE on CNR-2000.

The complement mode is highly effective for high-overlap references (>=90% of neighbors copied) — encoding 1-3 skipped positions is dramatically cheaper than encoding 15-20 copied positions.

## Greedy Search Cost Model

The greedy encoder evaluates all encoding options for each vertex and picks the one with minimum estimated bit cost. The cost model includes:

1. **VLC header cost** (1-11 bits depending on ref_mode and enc_type; 3 bits for empty)
2. **Reference distance cost** (Fibonacci varint of distance; x2 for multi-ref)
3. **Copy position cost** (adaptive: min of bitmap, copy-blocks, complement; x2 for multi-ref)
4. **Encoded list cost** (interval+residual, RLE, delta, or stop-delta)

### Cost Estimation Accuracy

**Fibonacci bit-length estimation** uses an exact precomputed Fibonacci lookup table (`_FIB_TABLE`) rather than the approximation formula `ceil(log2(n) * 1.44) + 2`. The exact computation eliminates ~1 bit systematic overestimation for small values (n <= 20), which is critical since most encoded values in web graphs are small.

**Interval encoding cost** estimation matches the actual writer exactly: interval/residual counts use `n+1` shift, subsequent interval starts use delta from previous end (not absolute), and residual gap values include the `+1` shift. A legacy phantom run-length pair count was also removed.

## Encoding Options (Action Space)

The greedy encoder evaluates all 28 combinations (no escape codes needed):

| ref_mode | enc_type | MIL | VLC bits | Description |
|----------|----------|-----|----------|-------------|
| reference | delta | - | 1 | Copy from ref + delta residuals |
| none | delta | - | 2 | Direct delta encoding |
| empty | - | - | 3 | No neighbors (degree-0 vertex) |
| none | interval | 2-5 | 6 | Direct interval+residual encoding |
| reference | interval | 2-5 | 7 | Copy from ref + interval residuals |
| none | rle | 2-5 | 8 | Direct hybrid mix encoding |
| reference | rle | 2-5 | 9 | Copy from ref + hybrid residuals |
| multi_ref | delta | - | 8 | Copy from 2 refs + delta residuals |
| multi_ref | interval | 2-5 | 11 | Copy from 2 refs + interval residuals |
| multi_ref | rle | 2-5 | 11 | Copy from 2 refs + hybrid residuals |

For each vertex, the combination with the lowest estimated bit cost is selected. The VLC per-vertex header records the choice, making the format fully self-describing. Multi-ref options are only evaluated when `multi_ref=true` is set in the header.

## Parameters

| Parameter | Default | Best config | Description |
|-----------|---------|-------------|-------------|
| integer_encoding | :zeta | :fibonacci | Global integer encoding for all varints |
| ref_window_size | 7 | 64 | Reference window size (recent vertices) |
| MIN_INTERVAL_LENGTH | 2-5 | 2-5 | Per-vertex, selected by greedy search |
| coding_scheme | :children | :children | Stop-delimited vertex encoding |
| copy_blocks | false | true | Copy-blocks for reference bitmap |
| stop_deltas | false | true | STOP-terminated delta lists |
| lr_split | false | true | LR-split residual encoding for intervals |
| multi_ref | false | true | Multi-reference: copy from two references per vertex |

**Implied parameters** (always active in v3.2, not configurable):

| Parameter | Value | Description |
|-----------|-------|-------------|
| adaptive_copy | =copy_blocks | 3-way adaptive copy when copy_blocks enabled |
| compact_copy | true | Compact prefix code for 3-way adaptive copy |
| tight_intervals | true | Interval start delta from end of prev interval |
| min_degree_for_ref | 1 | Minimum vertex degree to attempt reference search |
| min_overlap_for_ref | 1 | Minimum neighbor overlap to consider a reference candidate |

### Features Tested But Not Beneficial

| Feature | BPE Impact | Why It Didn't Help |
|---------|-----------|-------------------|
| `adaptive_deltas` | +0.06 BPE | Stop-deltas are almost universally better for web graphs; the 1-bit per-vertex flag overhead isn't recovered by the few vertices where count-prefix wins |
| `split_residual` | +0.04 BPE | The 1+4 bit encoding tag for residuals costs more than the savings from choosing a different encoding for the residual portion |
| `bv_blocks` | +0.12 BPE | BV-style alternating copy/skip blocks are less compact than existing (start, length) copy-blocks with gap encoding |
| `adaptive_encoding` | +0.15 BPE | Per-vertex 1-bit flag choosing intervals vs stop-delta is redundant with greedy search, which already evaluates all options |
| Recursive references | — | Shown not to work well in prior experiments |

## Vertex Ordering

Vertex ordering is the single largest factor in greedy MGS compression quality. The greedy encoder's reference window looks back at the previous `ref_window_size` vertices — if those vertices share many neighbors with the current vertex, compression is effective.

### Best Ordering: Leiden K=1 + Per-Group LLP

Two-step ordering: Leiden community detection (K=1, producing 2 groups) followed by per-group LLP on induced subgraphs. This creates tight sequential locality within each group while keeping the encoding simple (single-pass global window).

### Global LLP (slightly worse)

Global LLP with 10 passes and `:sym` mode across all vertices. Produces slightly worse locality than Leiden+LLP for the greedy encoder because it ignores community structure when ordering vertices.

### Key Finding: Greedy Beats BV and BV-HC

With Leiden+LLP ordering, BG (2.326 BPE) surpasses both WebGraph BV (2.898 BPE) and BV-HC (2.448 BPE) through a combination of per-vertex adaptive encoding, larger reference windows, multi-reference copying, aggressive low-degree reference search, and micro-optimizations (compact copy prefix, VLC, tight intervals). CS (2.304 BPE) is the overall best single-pass method at w=256. The ~0.003 BPE gap from BG to CG K=2 (2.329 BPE) is minimal:
- CG's **per-cluster local vertex indexing** (vertices renumbered 1...|C| within each cluster) enables much tighter intra-cluster reference encoding
- The greedy approach uses global vertex IDs throughout, limiting the locality gains achievable

The key breakthrough was lowering the minimum degree threshold for reference search from 3 to 1 (and overlap threshold from 3 to 1). On CNR-2000, 66K vertices have degree 1-2 (20% of non-empty vertices). Previously these were forced into raw encoding; now they can benefit from reference copying, saving ~0.30 BPE.

## Benchmark: CNR-2000

CNR-2000: 325,557 vertices, 3,216,152 edges (directed web graph).

### WebGraph BV Bit Budget (Reference: 2.898 BPE, LLP ordering)

From `cnr-2000.properties`:

| Component | Bits | BPE | Share |
|-----------|------|-----|-------|
| Outdegrees | 1,660,205 | 0.516 | 17.8% |
| References | 781,540 | 0.243 | 8.4% |
| Copy blocks | 1,353,080 | 0.421 | 14.5% |
| Intervals | 829,187 | 0.258 | 8.9% |
| Residuals | 4,694,729 | 1.460 | 50.4% |
| **Total** | **9,318,741** | **2.898** | 100% |

Key BV stats: `copiedarcs=2,195,145` (68.2%), `intervalisedarcs=443,657` (13.8%), `residualarcs=577,350` (17.9%), `avgdist=1.64`, `windowsize=7`, `minintervallength=4`, `zetak=3`.

### Optimization History

| Step | Change | BPE | Delta BPE |
|------|--------|-----|-----------|
| Baseline | LLP + Zeta-3, w=7, adaptive bitmap | 3.844 | — |
| Copy-blocks | Replace adaptive bitmap with copy-blocks | 3.528 | -0.316 |
| Larger window | w=7 -> w=64 | 3.355 | -0.173 |
| Fibonacci | Zeta-3 -> Fibonacci encoding | 3.268 | -0.087 |
| Leiden+LLP ordering | Global LLP -> Leiden K=1 + per-group LLP | 3.350 | +0.082 (ordering change) |
| 3-way adaptive copy | bitmap + copy-blocks + complement | 3.236 | -0.114 |
| STOP-terminated deltas | Eliminate count prefix on delta lists | 3.076 | -0.160 |
| Empty-prefix | 1-bit flag for degree-0 vertices | 3.031 | -0.044 |
| **Fix 1: VLC header cost** | Include header cost in greedy search | 2.994 | -0.037 |
| **Fix 2: Cost estimator** | Exact Fibonacci lengths + remove phantom count | 2.950 | -0.044 |
| **Fix 3: Interval cost** | Fix 4 discrepancies in interval/residual cost estimation | 2.932 | -0.018 |
| Empty-prefix (re-measured) | After all fixes | 2.956->2.912 | -0.044 |
| **Compact copy prefix** | 1-bit complement, 2-bit bitmap/blocks | 2.942->2.920 | -0.014 |
| **VLC v2 header** | 1-bit ref+delta code for reference-heavy streams | 2.946->2.900 | -0.010 |
| **Tight intervals** | Interval start delta from end of prev interval | 2.950->2.944 | -0.006 |
| **All winners combined** | empty + compact + vlc2 + tight | **2.881** | **-0.075** |
| **VLC v3** | Escape-free VLC with equal-cost MIL tags | ~2.881 | ~0 (faster decoding) |
| **LR-split residuals** | Split residuals at vertex_id, encode as distances | 2.817 | -0.048 |
| **Multi-reference** | Copy from two refs per vertex (top-K search) | 2.802 | -0.015 |
| **Merged VLC** | Empty code in VLC tree, remove 1-bit empty prefix | 2.792 | -0.010 |
| **Low-degree ref search** | Min degree 3→1, min overlap 3→1 for ref candidates | **2.493** | **-0.299** |

Fix 3 corrected four discrepancies between `estimate_interval_runlength_encoding_cost` and the actual `write_intervals_and_residuals`:
1. Interval count encoded as `n` instead of `n+1` (matching the writer's shift)
2. Residual count encoded as `n` instead of `n+1`
3. Subsequent interval starts used absolute position instead of delta from previous end
4. Residual gap values missing the `+1` shift applied by the writer

After all three fixes, the cost estimator matches actual written bits within 6 bits across 9.4M total bits.

Note: Fix 1, 2, and 3 improved all configurations simultaneously. The numbers above show Leiden+LLP ordering with cumulative improvements. With all fixes, the baseline (CB only) also improved from 3.349 to 3.154 BPE.

### Results Summary (Current Best)

| Config | Ordering | Features | BPE |
|--------|----------|----------|-----|
| **CS best (w=256)** | **Leiden+LLP** | **stop_deltas+adaptive_copy+low_deg_ref** | **2.304** |
| **BG best (w=64)** | **Leiden+LLP** | **copy_blocks+stop_deltas+multi_ref+low_deg_ref** | **2.326** |
| CG K=2 best | Original | all adaptive, w=64 | 2.329 |
| Pre-low-deg | Leiden+LLP | copy_blocks+stop_deltas+lr_split+multi_ref+merged_vlc | 2.792 |
| Best v3.1 (legacy) | Leiden+LLP | empty+compact+vlc2+tight | 2.881 |
| WebGraph BV | LLP | zeta-3, w=7 | 2.898 |
| WebGraph BV-HC | high compression (seq. access) | zeta-3, w=7 | 2.448 |

**CS beats WebGraph BV-HC**: 2.304 vs 2.448 BPE (-0.144 BPE)
**BG beats WebGraph BV**: 2.326 vs 2.898 BPE (-0.572 BPE)
**BG beats CG K=1**: 2.326 vs 2.638 BPE (-0.312 BPE) — HC uses high-compression sequential-access settings

### Multi-Dataset Benchmark

Best BG BPE across all tested datasets:

| Dataset | BG BPE | CG BPE | CS BPE | BV BPE | BG Config |
|---------|---------|---------|--------|--------|------------|
| cnr-2000 (no reorder) | 2.4929 | **2.3286** | 2.4348 | 2.898 | w=64, lr+mr |
| cnr-2000 (Leiden+LLP) | 2.3259 | 2.5652 | **2.3043** | 3.2335 | w=64, no-lr, mr |
| in-2004 | 1.895 | **1.7513** | 1.7839 | 2.172 | w=64, lr+mr |
| enwiki-2013 (Leiden+LLP) | **12.161** | 12.485 | 12.222 | 13.114 | w=64, lr+mr |
| enwiki-2013 (no reorder) | — | 15.718 | — | 13.114 | — |
| web-google core (Leiden+LLP) | 4.0735 | 4.3296 | **4.0288** | 5.0081 | w=64, no-lr, mr |
| web-google rcore (Leiden+LLP) | 3.7626 | 3.9359 | **3.7337** | 4.1751 | w=64, no-lr, mr |
| eat-core | 10.9191 | **10.5768** | 10.8679 | 10.705 | w=64, lr+mr |
| eat-core (Leiden+LLP) | 9.767 | **9.558** | 9.773 | 9.729 (HC) | w=64, lr+mr |
| eat-rcore | 9.5726 | **9.3192** | 9.5674 | 9.391 | w=64, lr+mr |
| arxiv-hep-ph core | 10.0910 | **9.8187** | 10.0767 | 10.262 | w=64, lr+mr |
| arxiv-hep-ph core (Leiden+LLP) | **7.252** | 7.162 | 7.269 | 7.706 (HC) | w=64, lr+mr |
| arxiv-hep-ph rcore | 9.0946 | **8.9406** | 8.9855 | 9.684 | w=64, lr+mr |
| amazon-0601 core | 12.4444 | **12.1853** | 12.5235 | 13.001 | w=8, lr+mr |
| amazon-0601 core (Leiden+LLP) | 8.058 | **7.903** | 8.095 | 8.722 (HC) | w=64, lr+mr |
| amazon-0601 rcore | 11.7310 | **11.7069** | 11.8497 | 12.064 | w=8, lr+mr |

**Key findings**:
- **enwiki-2013: BG beats CG** — BG (12.161, Leiden+LLP) outperforms CG (12.485, LLP) by 0.324 BPE. Multi-ref is valuable on the high-degree Wikipedia graph (avg degree 24.1). CG without reordering degrades to 15.718 BPE.
- Leiden+LLP ordering improves BG on Web-Google by ~0.3 BPE (core: 4.427→4.074, rcore: 4.080→3.763)
- BG still beats CG K=1 on Web-Google (core: 4.074 vs 4.330, rcore: 3.763 vs 3.936) and enwiki-2013, but CG K=1 beats BG on most other datasets
- LLP reordering is critical for SNAP datasets with poor original ordering (Web-Google: 1.2-1.4 BPE gain)
- lr_split helps on most datasets but hurts on cnr-2000 and Web-Google (crawl-ordered graphs)

### Synthetic Web Graph Benchmark (N=10000, original ordering)

BG parameter sweep against BV (w=64) on `random_web_digraph` graphs.
Sweep: 4 windows × 2 lr_split × 2 multi_ref × 2 encodings = 32 configs per degree.

| avg_deg | BV (w=64) | BG best | Config | Delta |
|---------|-----------|---------|--------|-------|
| 12 | 9.944 | **9.691** | fibonacci, w=8, lr=true, mr=true | -0.254 |
| 24 | 9.417 | **9.306** | zeta, w=8, lr=true, mr=false | -0.111 |
| 32 | 9.375 | **9.288** | zeta, w=8, lr=true, mr=false | -0.088 |
| 64 | 9.445 | **9.390** | zeta, w=64, lr=true, mr=false | -0.054 |

**Key findings**:
- **lr_split=true is required to beat BV** — same pattern as CS and CG; no lr_split=false config beats BV
- **Fibonacci wins at low degree (deg=12)**, zeta wins at medium-high degree (deg 24-64) — the crossover point is around deg~16
- **multi_ref has negligible impact** (<0.001 BPE difference) on synthetic web graphs — the tight locality means the best ref is always the immediate predecessor
- **Small window (w=8) is optimal for deg 12-32** — same pattern as CS and CG
- BG beats BV by less than CS (-0.25 vs -0.38 at deg=12) and much less than CG grid (-0.25 vs -0.63)
- **BG ranking**: worst of the three Adjacently methods on synthetic web graphs (BG < CS < CG grid)
- **Caveat**: tuned BV (zeta-5, i=2, m=-1) beats BG at all degrees, and beats all Adjacently methods at deg=64 (9.201 vs BG 9.390). BV's zeta-5 encoding is particularly effective for web-like gap distributions at higher density.

### LFR Benchmark (N=10000, avg_degree=15, tau1=2.5, tau2=1.5)

BG parameter sweep on LFR modular graphs with original and Leiden+LLP ordering.

| μ | Ordering | BV | BG best | Config | Delta |
|------|----------|------|---------|--------|-------|
| 0.05 | Original | 13.89 | **13.84** | zeta, w=64, lr, mr | -0.05 |
| 0.05 | Leiden+LLP | 8.44 | **8.36** | fib, w=64, lr, mr | -0.09 |
| 0.10 | Leiden+LLP | 9.13 | **9.06** | fib, w=64, lr, mr | -0.07 |
| 0.20 | Leiden+LLP | 9.95 | **9.91** | zeta, w=64, lr, mr | -0.04 |
| 0.30 | Leiden+LLP | 10.61 | **10.58** | zeta, w=64, lr, mr | -0.03 |
| 0.50 | Leiden+LLP | 11.60 | **11.56** | zeta, w=64, lr, mr | -0.04 |

**Key findings**:
- BG beats BV with Leiden+LLP at every μ, but margins are small (0.03-0.09 BPE)
- BG is the weakest Adjacently method on LFR (BG < CS < CG at every μ)
- Best config: zeta + w=64 + lr + mr (original) or fib + w=64 + lr + mr (Leiden+LLP at low μ)
- Fibonacci wins at low μ with Leiden+LLP (tight locality); zeta wins otherwise

### Best Greedy Config (v3.2)

```
Ordering:        Leiden K=1 + per-group LLP, :sym mode, 5 passes
Encoding:        Fibonacci varint
Window:          ref_window_size = 64
copy_blocks:     true   (implies adaptive_copy=true, compact_copy=true)
stop_deltas:     true
lr_split:        true   (LR-split residual encoding)
multi_ref:       true   (multi-reference per vertex, top-K search)
tight_intervals: true   (always active)
min_degree_ref:  1      (search refs for degree >= 1 vertices)
min_overlap_ref: 1      (consider refs with overlap >= 1)
BPE:             2.493 (varies slightly due to LLP randomness)
```

### Greedy Encoder Bit Budget (approximate)

Bit budget from per-vertex analysis (best v3.2 config, CNR-2000, **2.493 BPE**):

| Component | Est. BPE | Share |
|-----------|----------|-------|
| VLC headers (incl. empty) | ~0.46 | 18.5% |
| Copy positions (compact) | ~0.50 | 20.0% |
| Additions (ref residuals) | ~1.10 | 44.1% |
| Raw (non-ref vertices) | ~0.43 | 17.4% |
| Stream header | 0.00 | 0.0% |
| **Total** | **~2.493** | 100% |

**Vertex reference rate**: ~90%+ of non-empty vertices use a reference (single or multi), up from ~84% before the low-degree threshold change. The 66K degree 1-2 vertices previously forced into raw encoding now benefit from reference copying.

### Action Distribution (Actual, CNR-2000)

| Action | % | VLC bits |
|--------|---|----------|
| reference + delta | ~52% | 1 |
| none + delta | ~35% | 2 |
| none + interval(2) | ~6% | 5 |
| multi-ref + delta | ~2% | 7 |
| reference + interval(2) | ~1.5% | 6 |
| other (interval/rle, MIL 3-5) | ~3.5% | 5-10 |

Shannon entropy of action distribution: **~1.7 bits** (theoretical minimum). The VLC assigns 1 bit to the dominant action (ref+delta, ~52%), approaching the entropy limit.

### Comparison with WebGraph BV

With Leiden+LLP ordering, BG **beats** both WebGraph BV and BV-HC (2.326 vs 2.898 vs 2.448 BPE). The key structural differences:

| Aspect | Greedy | WebGraph BV |
|--------|--------|-------------|
| Outdegrees | Implicit (stop-delimited / VLC empty code) | Explicit (0.516 BPE) |
| Per-vertex header | VLC (1-11 bits, avg ~1.7 bits, incl. empty) | None (fixed pipeline) |
| Encoding selection | Per-vertex greedy cost search (27 options) | Fixed: intervals -> residuals |
| Reference count | 1-2 refs per vertex (multi-ref) | 1 ref per vertex |
| Reference window | w=64 | w=7 |
| Integer encoding | Fibonacci | Zeta-3 |
| Copy positions | 3-way adaptive (compact prefix) | Alternating block lengths |
| Residuals | LR-split (distances from vertex_id) | Zigzag from vertex_id |

The greedy encoder's advantage comes from eliminating explicit outdegree encoding (saves ~0.5 BPE via merged VLC empty code), using a larger reference window (w=64 vs w=7), multi-reference copying, and aggressive low-degree reference search (min degree 1, min overlap 1). These savings more than compensate for the per-vertex VLC header overhead (~0.46 BPE) and the cost of adaptive encoding selection.

### Potential Further Improvements

BG (2.326) and CS (2.304) have surpassed both WebGraph BV (2.898) and BV-HC (2.448), and are within 0.003-0.025 BPE of CG K=2 (2.329). The three methods have effectively converged on CNR-2000 with Leiden+LLP ordering:

| Opportunity | Est. savings | Feasibility |
|-------------|-------------|-------------|
| Per-cluster local IDs | 0.005-0.01 BPE | Requires architectural change (CG-style) |
| Arithmetic coding of headers | 0.005-0.01 BPE | Complexity vs marginal gain |
| Window size tuning per-cluster | 0.005 BPE | Adaptive window based on cluster density |

The near-convergence of all three methods suggests that the Leiden+LLP ordering is the dominant factor at this compression level, with diminishing returns from encoder improvements.
