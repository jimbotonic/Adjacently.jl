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
| 0x10–0x8F | RL policy mode (policy_id = value - 0x10 + 1) |
| 0xFF | Huffman (deprecated) |

## RL-Compressed Data Section

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

Each vertex is encoded with a **VLC (variable-length coded) header** followed by payload.

The top 4 action combinations cover ~97% of vertices in well-ordered web graphs. They receive short prefix codes; all other combinations use a 3-bit escape prefix followed by the full 6-bit fixed encoding.

```
VLC header (prefix-free):
  00   (2 bits) = none      + delta
  01   (2 bits) = none      + interval, MIL=2
  10   (2 bits) = reference + delta
  110  (3 bits) = reference + interval, MIL=2
  111  (3 bits) = escape → full 6-bit header follows:
                    ref_mode (2 bits) + enc_opt (4 bits)
```

**Escape payload — ref_mode (2 bits):**
```
  00 = none
  01 = reference
  10 = recursive
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

**Average header cost** (CNR-2000, LLP + Zeta-3): ~2.3 bits/vertex vs 6.0 bits fixed.

### Per-Vertex Payload

#### When ref_mode = none (00)

The full neighbor list is encoded using the selected encoding option:

**Interval+Residual** (enc_opt 0000–0011):
```
varint(num_intervals + 1)
For each interval:
  varint(start_delta)               // delta from prev start (or absolute for first)
  varint(length - MIL + 1)          // length shifted by MIL
varint(num_residuals + 1)
delta_encoded_residuals             // first value absolute, then gaps+1
```

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
```
varint(count + 1)                   // neighbor count
varint(first_value)                 // first neighbor (absolute)
varint(gap₂ + 1) ... varint(gapₙ + 1)  // subsequent gaps shifted by +1
```

#### When ref_mode ≠ none (01 or 10)

```
varint(ref_distance)                // distance back in reference window
adaptive_bitmap                     // copy bitmap (which ref neighbors to copy)
residual_payload                    // residuals encoded per enc_opt above
```

### Adaptive Bitmap Format

```
small_count(bitmap_length):
  00 = 0 elements
  01 = 1 element
  10 = 2 elements
  11 + varint(length) = 3+ elements

if length > 0:
  format_flag (1 bit):
    0 = block encoding (WebGraph-style alternating blocks)
    1 = raw bits

Block encoding:
  varint(num_blocks + 1)
  For each block:
    varint(block_length + 1)
  Blocks alternate: copy, skip, copy, skip, ...
  If bitmap starts with 0s, first block has length 0

Raw encoding:
  N raw bits (one per bitmap entry)
```

The encoder chooses block vs raw by estimating both costs and picking the cheaper option.

## Copy-Blocks Reference Encoding

When `copy_blocks=true`, the reference bitmap (which neighbor positions to copy) is replaced by a run-length encoding of contiguous copy runs. This is the same format as WebGraph BV.

```
Copy-blocks format:
  small_count(num_copy_blocks)       // number of contiguous runs to copy
  IF num_copy_blocks > 0:
    varint(first_start)              // 1-based index of first copied position (≥ 1)
    For each subsequent block (2nd onward):
      varint(skip_length + 1)        // positions skipped since end of last block
      // run length is implicit: consumes one position per step until next skip
```

The encoder picks the cheaper of adaptive bitmap and copy-blocks for each vertex. Copy-blocks wins when copied positions form long contiguous runs — common under LLP ordering. It saves ~0.32 BPE vs adaptive bitmap on CNR-2000 at window=64.

**Note**: the `copy_blocks` flag is currently passed externally (not stored in the file header). The decoder must receive it out-of-band (e.g., from a companion configuration).

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
| copy_blocks | false | true | WebGraph-style copy-blocks for reference bitmap |

## Vertex Ordering

Vertex ordering is the single largest factor in greedy MGS compression quality. The greedy encoder's reference window looks back at the previous `ref_window_size` vertices — if those vertices share many neighbors with the current vertex, compression is effective.

### Global LLP (recommended)

Layered Label Propagation (LLP) with 10 passes and `:sym` mode is the best ordering for the greedy approach. LLP runs across all 325,557 vertices jointly, optimizing the global sequential locality. It considers both in- and out-edges, which matches the symmetric reference window behavior.

### Two-step Leiden + per-group LLP (not beneficial for greedy)

The RCGE format uses a two-step ordering: Leiden community detection → per-group LLP on induced subgraphs. This dramatically improves RCGE because RCGE's reference window is **cluster-local** (only looks within the same encoding cluster).

For the greedy approach, both the two-step ordering and cluster-local window resets were tested (Leiden K=1 → 2 groups of 18K + 307K vertices, LLP on each induced subgraph, optional window reset at group boundary):

| Ordering | Window mode | Codec | w | BPE |
|----------|-------------|-------|---|-----|
| Global LLP | global | Fibonacci | 64 | **3.268** |
| Leiden K=1 + per-group LLP | global | Fibonacci | 64 | 3.326 |
| Leiden K=1 + per-group LLP | **cluster-local** | Fibonacci | 64 | 3.326 |

The cluster-local window reset (resetting the reference window at the group boundary) produces **byte-for-byte identical output** to the global window. The reason: the greedy encoder already never picks cross-group references at the K=1 boundary — they are universally non-beneficial since group 1 (18K vertices, one large Leiden community) and group 2 (307K vertices, all others) have completely disjoint local neighborhoods. The greedy search already avoids them without needing an explicit reset.

The deeper reason per-group LLP still underperforms global LLP is that it ignores cross-group edges when computing the ordering, slightly degrading end-to-end sequential locality. Global LLP jointly optimizes all edges and produces better BPE for the greedy global reference window.

**Key finding**: the RCGE advantage over greedy is **architectural**, not replicable by vertex ordering or window management changes:
- RCGE uses **per-cluster local vertex indexing** (vertices renumbered 1…|C| within each cluster), enabling much tighter intra-cluster reference encoding
- RCGE's reference window is structurally cluster-local by design — it indexes into a local vertex space, not a global one
- The greedy approach uses global vertex IDs throughout; resetting the window boundary cannot reproduce RCGE's per-cluster locality gains

## Benchmark: CNR-2000

CNR-2000: 325,557 vertices, 3,216,152 edges (directed web graph).

### Results summary

| Config | Relabeling | Codec | Window | Local window | BPE |
|--------|-----------|-------|--------|--------------|-----|
| [1] LLP baseline | Global LLP | Zeta-3 | 7 | No | 3.844 |
| [2] LLP + CB | Global LLP | Zeta-3 | 7 | No | 3.528 |
| [3] LLP + CB w=64 | Global LLP | Zeta-3 | 64 | No | 3.355 |
| [4] **LLP + Fib + CB w=64** | **Global LLP** | **Fibonacci** | **64** | **No** | **3.268** |
| [5] Leiden+LLP, global w | Leiden K=1 + LLP | Fibonacci | 7 | No | 3.514 |
| [6] Leiden+LLP, global w | Leiden K=1 + LLP | Fibonacci | 64 | No | 3.326 |
| [7] Leiden+LLP, local w | Leiden K=1 + LLP | Fibonacci | 7 | Yes | 3.514 |
| [8] Leiden+LLP, local w | Leiden K=1 + LLP | Fibonacci | 64 | Yes | 3.326 |
| WebGraph BV | LLP | Zeta-3 | 7 | No | 2.897 |
| RCGE FW64 K=1 | Leiden + LLP | Fibonacci | 64 | Cluster-local | **2.887** |

Configs 7 and 8 are byte-for-byte identical to configs 5 and 6 — the cluster-local window reset has no effect (see Vertex Ordering section).

### Best greedy config (config 4)

```
Relabeling:  global LLP, :sym mode, 10 passes
Encoding:    Fibonacci varint
Window:      ref_window_size = 64
Copy-blocks: enabled
BPE:         ~3.27 (varies slightly due to LLP randomness)
```

### Optimization history

| Step | Change | BPE |
|------|--------|-----|
| Baseline | LLP + Zeta-3, w=7, adaptive bitmap | 3.844 |
| Copy-blocks | Replace adaptive bitmap with copy-blocks | 3.528 (−0.316) |
| Larger window | w=7 → w=64 | 3.355 (−0.173) |
| Fibonacci | Zeta-3 → Fibonacci encoding | **3.268** (−0.087) |

### Remaining gap vs WebGraph (0.37 BPE)

The greedy approach cannot close to WebGraph through vertex ordering, window management, or codec changes:

- **Two-step Leiden+LLP ordering**: slightly worse than global LLP for the greedy global window (3.326 vs 3.268 BPE)
- **Cluster-local window reset**: zero effect — the greedy encoder already avoids cross-cluster references
- **Architecture**: the gap is structural — RCGE's per-cluster local indexing gives ~0.4 BPE advantage that cannot be replicated in the greedy framework

See FORMAT_SPECS_RCGE.md for the RCGE approach which achieves 2.887 BPE.

### Action distribution (LLP + Zeta-3, w=7, CNR-2000)

| Action | % | VLC bits |
|--------|---|----------|
| none + delta | 31.0% | 2 |
| none + interval(2) | 29.1% | 2 |
| reference + delta | 25.0% | 2 |
| reference + interval(2) | 11.9% | 3 |
| other (14 combos) | 3.0% | 9 |

Shannon entropy of action distribution: **2.13 bits** (theoretical minimum).
