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

| Parameter | Default | Description |
|-----------|---------|-------------|
| integer_encoding | :zeta (V3) | Global integer encoding for all varints |
| ZETA_BASE | 3 | Zeta coding parameter k |
| ref_window_size | 7 | Reference window size (recent vertices) |
| MIN_INTERVAL_LENGTH | 2–5 | Per-vertex, selected by greedy search |
| coding_scheme | :children | Stop-delimited vertex encoding |

## Benchmark: CNR-2000

### WebGraph Reference
- **2.90 bits/edge** (BPE)
- Zeta-3 encoding, LLP ordering, window=7

### MGS V3 (Zeta-3 + VLC headers)
- **~4.68 BPE** (greedy, LLP ordering, window=7)

### Action distribution (LLP + Zeta-3, CNR-2000)

| Action | % | VLC bits |
|--------|---|----------|
| none + delta | 31.0% | 2 |
| none + interval(2) | 29.1% | 2 |
| reference + delta | 25.0% | 2 |
| reference + interval(2) | 11.9% | 3 |
| other (14 combos) | 3.0% | 9 |

Shannon entropy of action distribution: **2.13 bits** (theoretical minimum).

### Remaining gaps vs WebGraph
- **Vertex ordering**: LLP quality (number of passes, resolution levels)
- **Reference window**: WebGraph uses window=7 with highly optimized LLP; larger windows may help with weaker orderings
- **Cost estimation accuracy**: WebGraph uses exact bit counting; MGS uses estimates
