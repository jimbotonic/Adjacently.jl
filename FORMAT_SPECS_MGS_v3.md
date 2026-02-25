# MGS Graph Compression Format Specification V3

## Overview

MGS v3 extends v2 with two changes:

1. **Zeta-3 integer encoding** (ζ₃), matching WebGraph's encoding choice.
2. **VLC vertex headers** — variable-length prefix codes replace the fixed 6-bit per-vertex headers, saving ~3.7 bits/vertex (~0.37 BPE on CNR-2000).

### Key changes from V2

| Aspect | V2 | V3 |
|--------|----|----|
| Integer encoding | Fibonacci (~1.44·log₂(n) + 2 bits) | Zeta-3 (~4 bits for values 1–7) |
| Stream header tag | `000` (Fibonacci) | `001` (Zeta) |
| Per-vertex header | Fixed 6 bits (ref_mode 2b + enc_opt 4b) | VLC: 2–9 bits (avg ~2.3 bits) |

## Integer Encoding: Zeta-3

Zeta coding with parameter k=3 encodes a positive integer v as:

1. Compute h = ⌊log₂(v) / k⌋
2. Write h+1 in **unary** (h zeros followed by a 1): h+1 bits
3. Write the remainder in **truncated binary**: k·(h+1) bits

**Total bits**: (h+1) + k·(h+1) = (h+1)·(k+1)

| Value range | h | Bits |
|-------------|---|------|
| 1–7 | 0 | 4 |
| 8–63 | 1 | 8 |
| 64–511 | 2 | 12 |
| 512–4095 | 3 | 16 |

For web graph gaps (typically 1–100), most values cost 4–8 bits vs Fibonacci's 5–10 bits.

## File Structure

```
[Header Section (12 bytes)] [Bitstream Data Section]
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

**Integer encoding (4 bits)**:

| Code | Encoding |
|------|----------|
| 0x1 | Elias gamma |
| 0x2 | Elias delta |
| 0x3 | Golomb |
| 0x4 | FED |
| 0x5 | **Zeta** (k=3) ← V3 default |
| 0x6 | Fibonacci |

**Option flags**:

| Value | Mode |
|-------|------|
| 0x00 | None (standard encoding) |
| 0x0F | ASTRA (full adaptive encoding) |
| 0x10–0x8F | Greedy mode (option_id = value - 0x10 + 1) |
| 0xFF | Huffman (deprecated) |

## Greedy Compressed Data Section

### Stream Header (4 bits)

```
coding_scheme_flag (1 bit):
  0 = :children (stop-delimited)
  1 = :index (degree array precedes data)

encoding_tag (3 bits):
  000 = Fibonacci
  001 = Zeta        ← V3 uses this
  010 = Elias gamma
  011 = Elias delta
```

The encoding tag makes the stream **self-describing** — the decoder reads it and uses the corresponding integer encoding for all subsequent values. No policy file is needed for decompression.

### Optional Index Section

When `coding_scheme_flag = 1` (`:index` mode):
```
For each vertex v = 1..N:
  varint(degree + 1)    // +1 to avoid zero
```

### Per-Vertex Record

Each vertex is encoded with an optional **empty-prefix flag**, then a **VLC (variable-length coded) header** followed by payload.

#### Empty-Prefix Optimization

When `empty_prefix=true`, each vertex begins with a 1-bit flag:
```
0 = vertex has no neighbors (done, skip to next vertex)
1 = vertex has neighbors (proceed with VLC header + payload)
```

This saves encoding overhead for degree-0 vertices. On CNR-2000 (~3% degree-0), this contributes ~0.04 BPE savings.

#### VLC Header

The top 4 action combinations cover ~97% of vertices in well-ordered web graphs. They receive short prefix codes; all other combinations use an escape prefix followed by the full 6-bit fixed encoding. Two VLC variants are available:

**VLC v1** (`vlc2=false`, original):
```
  00   (2 bits) = none      + delta
  01   (2 bits) = none      + interval, MIL=2
  10   (2 bits) = reference + delta
  110  (3 bits) = reference + interval, MIL=2
  111  (3 bits) = escape → full 6-bit header follows:
                    ref_mode (2 bits) + enc_opt (4 bits)
```

**VLC v2** (`vlc2=true`, optimized for reference-heavy streams):
```
  0    (1 bit)  = reference + delta
  10   (2 bits) = none      + delta
  110  (3 bits) = none      + interval, MIL=2
  1110 (4 bits) = reference + interval, MIL=2
  1111 (4 bits) = escape → full 6-bit header follows:
                    ref_mode (2 bits) + enc_opt (4 bits)
```

VLC v2 assigns a 1-bit code to the most common action (ref+delta, ~53% of vertices), saving ~0.5 bits/vertex on reference-heavy graphs at the cost of 1 extra bit for no-ref actions. Net savings on CNR-2000: ~0.01 BPE.

**Escape payload — ref_mode (2 bits):**
```
  00 = none
  01 = reference
  10 = recursive (unused in practice)
```

**Escape payload — enc_opt (4 bits):**
```
  0000 = Interval+Residual, MIL=2
  0001 = Interval+Residual, MIL=3
  0010 = Interval+Residual, MIL=4
  0011 = Interval+Residual, MIL=5
  0100 = RLE/Hybrid Mix,    MIL=2
  0101 = RLE/Hybrid Mix,    MIL=3
  0110 = RLE/Hybrid Mix,    MIL=4
  0111 = RLE/Hybrid Mix,    MIL=5
  1000 = Delta
  1001–1111 = reserved
```

**Average header cost** (CNR-2000, LLP + Fibonacci): ~2.3 bits/vertex vs 6.0 bits fixed.

**VLC header cost in greedy search**: The `_vlc_header_cost` function returns the exact bit cost (2, 3, or 9) for each `(ref_mode, enc_type, mil)` combination. The greedy search includes this cost in all comparisons, ensuring the encoder accounts for the variable header overhead when selecting encoding strategies. This is critical — without it, the search systematically undervalues cheap-header encodings (none+delta, ref+delta) vs expensive-header ones (escape path).

### Per-Vertex Payload

#### When ref_mode = none (00)

The full neighbor list is encoded using the selected encoding option:

**Interval+Residual** (enc_opt 0000–0011):
```
varint(num_intervals + 1)
For each interval:
  varint(start_delta)               // zigzag(start - vertex_id) + 1 for first; delta from prev for subsequent
  varint(length - MIL + 1)          // length shifted by MIL
varint(num_residuals + 1)
delta_encoded_residuals             // first value: zigzag(v - vertex_id) + 1, then gaps+1
```

When `tight_intervals=true`, subsequent interval starts are delta-encoded from the **end** of the previous interval (start + length) rather than from the start. This produces smaller deltas when intervals are close together, saving ~0.006 BPE on CNR-2000.

**RLE/Hybrid Mix** (enc_opt 0100–0111):
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

**Delta** (enc_opt 1000):

Two sub-formats depending on `stop_deltas` configuration:

*Count-prefixed delta* (`stop_deltas=false`):
```
varint(count + 1)                   // neighbor count
varint(zigzag(first - vertex_id) + 1)   // first neighbor
varint(gap₂ + 1) ... varint(gapₙ + 1)  // subsequent gaps shifted by +1
```

*STOP-terminated delta* (`stop_deltas=true`):
```
For each value in sorted order:
  bit '1'                           // continuation flag
  varint(gap)                       // zigzag(first - vertex_id) + 1 for first; gap for rest
bit '0'                             // STOP terminator
```

STOP-terminated deltas eliminate the count prefix (2–6 bits), which saves bits for low-degree vertices and costs slightly more for high-degree vertices. Net savings on CNR-2000: ~0.08 BPE.

#### When ref_mode ≠ none (01 or 10)

```
varint(ref_distance)                // distance back in reference window
copy_encoding                       // which ref neighbors to copy (see below)
residual_payload                    // residuals encoded per enc_opt above
```

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

**Default prefix code** (`compact_copy=false`):
```
outer_bit:
  1 → raw bitmap (ref_len bits, 1 bit overhead)
  0 → inner_bit:
    0 → copy-blocks encoding of COPIED positions (2 bits overhead)
    1 → copy-blocks encoding of SKIPPED positions (complement, 2 bits overhead)
```

**Compact prefix code** (`compact_copy=true`):
```
  0 → complement copy-blocks (1 bit overhead)
  10 → copy-blocks (2 bits overhead)
  11 → raw bitmap (2 bits overhead)
```

The compact prefix assigns the cheapest code (1 bit) to complement mode, which is the most common choice for high-overlap references. This saves ~0.01 BPE on CNR-2000.

The complement mode is highly effective for high-overlap references (≥90% of neighbors copied) — encoding 1–3 skipped positions is dramatically cheaper than encoding 15–20 copied positions.

## Greedy Search Cost Model

The greedy encoder evaluates all encoding options for each vertex and picks the one with minimum estimated bit cost. The cost model includes:

1. **VLC header cost** (2, 3, or 9 bits) — included since Fix 1
2. **Reference distance cost** (Fibonacci varint of distance)
3. **Copy position cost** (adaptive: min of bitmap, copy-blocks, complement)
4. **Encoded list cost** (interval+residual, RLE, delta, or stop-delta)
5. **Empty-prefix flag** (1 bit when enabled)

### Cost Estimation Accuracy

**Fibonacci bit-length estimation** uses an exact precomputed Fibonacci lookup table (`_FIB_TABLE`) rather than the approximation formula `ceil(log₂(n) × 1.44) + 2`. The exact computation eliminates ~1 bit systematic overestimation for small values (n ≤ 20), which is critical since most encoded values in web graphs are small.

**Interval encoding cost** estimation matches the actual writer exactly: interval/residual counts use `n+1` shift, subsequent interval starts use delta from previous end (not absolute), and residual gap values include the `+1` shift. A legacy phantom run-length pair count was also removed.

## Encoding Options (Action Space)

The greedy encoder evaluates all 27 combinations:

| ref_mode | enc_opt | Description |
|----------|---------|-------------|
| none | interval, MIL=2..5 | Direct interval+residual encoding |
| none | rle, MIL=2..5 | Direct hybrid mix encoding |
| none | delta | Direct delta encoding |
| reference | interval, MIL=2..5 | Copy from ref + interval residuals |
| reference | rle, MIL=2..5 | Copy from ref + hybrid residuals |
| reference | delta | Copy from ref + delta residuals |
| recursive | interval, MIL=2..5 | Recursive ref + interval residuals |
| recursive | rle, MIL=2..5 | Recursive ref + hybrid residuals |
| recursive | delta | Recursive ref + delta residuals |

For each vertex, the combination with the lowest estimated bit cost is selected. The VLC per-vertex header records the choice, making the format fully self-describing.

## Parameters

| Parameter | Default | Best greedy config | Description |
|-----------|---------|-------------------|-------------|
| integer_encoding | :zeta | :fibonacci | Global integer encoding for all varints |
| ZETA_BASE | 3 | — | Zeta coding parameter k |
| ref_window_size | 7 | 64 | Reference window size (recent vertices) |
| MIN_INTERVAL_LENGTH | 2–5 | 2–5 | Per-vertex, selected by greedy search |
| coding_scheme | :children | :children | Stop-delimited vertex encoding |
| copy_blocks | false | true | Copy-blocks for reference bitmap |
| adaptive_copy | false | true | 3-way adaptive copy (bitmap/blocks/complement) |
| stop_deltas | false | true | STOP-terminated delta lists |
| empty_prefix | false | true | 1-bit per-vertex empty flag |
| compact_copy | false | true | Compact prefix code for 3-way adaptive copy |
| vlc2 | false | true | VLC v2 header (1-bit ref+delta code) |
| tight_intervals | false | true | Interval start delta from end of prev interval |
| adaptive_deltas | false | false | Per-vertex stop vs count-prefix choice (not beneficial) |
| split_residual | false | false | Independent residual encoding (not beneficial) |
| bv_blocks | false | false | BV-style alternating blocks (not beneficial) |

### Features Tested But Not Beneficial

| Feature | BPE Impact | Why It Didn't Help |
|---------|-----------|-------------------|
| `adaptive_deltas` | +0.06 BPE | Stop-deltas are almost universally better for web graphs; the 1-bit per-vertex flag overhead isn't recovered by the few vertices where count-prefix wins |
| `split_residual` | +0.04 BPE | The 1+4 bit encoding tag for residuals costs more than the savings from choosing a different encoding for the residual portion |
| `bv_blocks` | +0.12 BPE | BV-style alternating copy/skip blocks are less compact than existing (start, length) copy-blocks with gap encoding. The explicit start positions in copy-blocks allow better compression when copied positions form sparse runs |
| Recursive references | Not tested recently | Shown not to work well in prior experiments |

## Vertex Ordering

Vertex ordering is the single largest factor in greedy MGS compression quality. The greedy encoder's reference window looks back at the previous `ref_window_size` vertices — if those vertices share many neighbors with the current vertex, compression is effective.

### Best Ordering: Leiden K=1 + Per-Group LLP

Two-step ordering: Leiden community detection (K=1, producing 2 groups) followed by per-group LLP on induced subgraphs. This creates tight sequential locality within each group while keeping the encoding simple (single-pass global window).

### Global LLP (slightly worse)

Global LLP with 10 passes and `:sym` mode across all vertices. Produces slightly worse locality than Leiden+LLP for the greedy encoder because it ignores community structure when ordering vertices.

### Key Finding: Greedy Beats BV, RCGE Remains Superior

The greedy approach now surpasses WebGraph BV (2.881 vs 2.897 BPE) through a combination of per-vertex adaptive encoding, larger reference windows, and micro-optimizations (compact copy prefix, VLC v2, tight intervals). However, the 0.45 BPE gap to RCGE (2.434 BPE) is structural:
- RCGE's **per-cluster local vertex indexing** (vertices renumbered 1…|C| within each cluster) enables much tighter intra-cluster reference encoding
- The greedy approach uses global vertex IDs throughout, limiting the locality gains achievable

## Benchmark: CNR-2000

CNR-2000: 325,557 vertices, 3,216,152 edges (directed web graph).

### WebGraph BV Bit Budget (Reference: 2.897 BPE)

From `cnr-2000.properties`:

| Component | Bits | BPE | Share |
|-----------|------|-----|-------|
| Outdegrees | 1,660,205 | 0.516 | 17.8% |
| References | 781,540 | 0.243 | 8.4% |
| Copy blocks | 1,353,080 | 0.421 | 14.5% |
| Intervals | 829,187 | 0.258 | 8.9% |
| Residuals | 4,694,729 | 1.460 | 50.4% |
| **Total** | **9,318,741** | **2.897** | 100% |

Key BV stats: `copiedarcs=2,195,145` (68.2%), `intervalisedarcs=443,657` (13.8%), `residualarcs=577,350` (17.9%), `avgdist=1.64`, `windowsize=7`, `minintervallength=4`, `zetak=3`.

### Optimization History

| Step | Change | BPE | Δ BPE |
|------|--------|-----|-------|
| Baseline | LLP + Zeta-3, w=7, adaptive bitmap | 3.844 | — |
| Copy-blocks | Replace adaptive bitmap with copy-blocks | 3.528 | −0.316 |
| Larger window | w=7 → w=64 | 3.355 | −0.173 |
| Fibonacci | Zeta-3 → Fibonacci encoding | 3.268 | −0.087 |
| Leiden+LLP ordering | Global LLP → Leiden K=1 + per-group LLP | 3.350 | +0.082 (ordering change) |
| 3-way adaptive copy | bitmap + copy-blocks + complement | 3.236 | −0.114 |
| STOP-terminated deltas | Eliminate count prefix on delta lists | 3.076 | −0.160 |
| Empty-prefix | 1-bit flag for degree-0 vertices | 3.031 | −0.044 |
| **Fix 1: VLC header cost** | Include header cost in greedy search | 2.994 | −0.037 |
| **Fix 2: Cost estimator** | Exact Fibonacci lengths + remove phantom count | 2.950 | −0.044 |
| **Fix 3: Interval cost** | Fix 4 discrepancies in interval/residual cost estimation | 2.932 | −0.018 |
| Empty-prefix (re-measured) | After all fixes | 2.956→2.912 | −0.044 |
| **Compact copy prefix** | 1-bit complement, 2-bit bitmap/blocks | 2.942→2.920 | −0.014 |
| **VLC v2 header** | 1-bit ref+delta code for reference-heavy streams | 2.946→2.900 | −0.010 |
| **Tight intervals** | Interval start delta from end of previous interval | 2.950→2.944 | −0.006 |
| **All winners combined** | empty + compact + vlc2 + tight | **2.881** | **−0.075** |

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
| CB only baseline | Leiden+LLP | copy_blocks | 3.154 |
| + adaptive | Leiden+LLP | + adaptive_copy | 3.039 |
| + stop_d | Leiden+LLP | + stop_deltas | 2.956 |
| + empty | Leiden+LLP | + empty_prefix | 2.912 |
| + compact | Leiden+LLP | + compact_copy | 2.920 |
| + vlc2 | Leiden+LLP | + vlc2 | 2.900 |
| + tight | Leiden+LLP | + tight_intervals | 2.886 |
| **Best combo** | **Leiden+LLP** | **empty+compact+vlc2+tight** | **2.881** |
| WebGraph BV | LLP | zeta-3, w=7 | 2.897 |
| RCGE best | Leiden+LLP | all adaptive | **2.434** |

**Greedy beats WebGraph BV**: 2.881 vs 2.897 BPE (−0.016 BPE, −51 Kbit, ~6.4 KB smaller)

### Best Greedy Config

```
Ordering:        Leiden K=1 + per-group LLP, :sym mode, 5 passes
Encoding:        Fibonacci varint
Window:          ref_window_size = 64
copy_blocks:     true
adaptive_copy:   true
stop_deltas:     true
empty_prefix:    true
compact_copy:    true
vlc2:            true
tight_intervals: true
BPE:             ~2.88 (varies slightly due to LLP randomness)
```

### Greedy Encoder Bit Budget (approximate)

Bit budget from per-vertex analysis (best config: empty+compact+vlc2+tight, CNR-2000, **2.881 BPE**, 1,158,204 bytes):

| Component | Est. Bits | Est. BPE | Share |
|-----------|-----------|----------|-------|
| Empty-prefix flags | 325,557 | 0.101 | 3.5% |
| VLC v2 headers | ~480,000 | ~0.149 | 5.2% |
| Reference distances | ~543,000 | ~0.169 | 5.9% |
| Copy positions (compact) | ~1,400,000 | ~0.435 | 15.1% |
| Body (enc data) | ~6,517,000 | ~2.027 | 70.3% |
| Stream header | 4 | 0.000 | 0.0% |
| **Total** | **~9,265,632** | **2.881** | 100% |

**Vertex reference rate**: ~84% of non-empty vertices use a reference. This is significantly higher than BV's 68.2% `copiedarcs` ratio, enabled by the larger window (w=64 vs w=7).

### Action Distribution (Actual, CNR-2000)

| Action | % | VLC v1 bits | VLC v2 bits |
|--------|---|------------|------------|
| reference + delta | 53.3% | 2 | **1** |
| none + delta | 35.9% | 2 | 2 |
| none + interval(2) | 6.3% | 2 | 3 |
| reference + interval(2) | ~1.5% | 3 | 4 |
| other (escape path) | ~3.0% | 9 | 10 |

Shannon entropy of action distribution: **~1.6 bits** (theoretical minimum). VLC v2 assigns 1 bit to the dominant action (ref+delta, 53%), approaching the entropy limit.

### Comparison with WebGraph BV

The greedy encoder now **beats** WebGraph BV (2.881 vs 2.897 BPE). The key structural differences:

| Aspect | Greedy | WebGraph BV |
|--------|--------|-------------|
| Outdegrees | Implicit (stop-delimited / empty-prefix) | Explicit (0.516 BPE) |
| Per-vertex header | VLC v2 (1–10 bits, avg ~2 bits) | None (fixed pipeline) |
| Encoding selection | Per-vertex greedy cost search (9 options) | Fixed: intervals → residuals |
| Reference window | w=64 | w=7 |
| Integer encoding | Fibonacci | Zeta-3 |
| Copy positions | 3-way adaptive (compact prefix) | Alternating block lengths |

The greedy encoder's advantage comes from eliminating explicit outdegree encoding (saves ~0.5 BPE) and using a larger reference window (w=64 vs w=7). These savings more than compensate for the per-vertex VLC header overhead (~0.16 BPE) and the cost of adaptive encoding selection.

### Potential Further Improvements

The greedy encoder has surpassed WebGraph BV. Remaining opportunities to close the gap to RCGE (2.434 BPE):

| Opportunity | Est. savings | Feasibility |
|-------------|-------------|-------------|
| Per-cluster local IDs | 0.10–0.15 BPE | Requires architectural change (RCGE-style) |
| Per-cluster bit-matrix | 0.10–0.20 BPE | Requires architectural change (RCGE-style) |
| Arithmetic coding of headers | 0.005–0.01 BPE | Complexity vs marginal gain |
| Window size tuning per-cluster | 0.005 BPE | Adaptive window based on cluster density |

The remaining 0.45 BPE gap to RCGE is largely structural — RCGE's per-cluster local vertex indexing and hierarchical encoding enable fundamentally tighter compression that the single-pass greedy approach cannot match.
