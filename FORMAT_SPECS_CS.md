# Command Stream (CS) Format Specification

## Overview

Command Stream (CS) is a single-pass graph compression format that replaces the greedy encoder's separate empty-prefix flag and VLC vertex header with a unified frequency-optimized prefix code tree. All proven greedy optimizations are preserved: 3-way adaptive copy, STOP-terminated deltas, tight intervals, Fibonacci encoding, and w=64 reference windows. Only the per-vertex header layer changes.

Key properties:
- Unified prefix code merges empty vertex flag, reference mode, and encoding type into a single code tree
- Frequency-optimal bit assignment based on CNR-2000 action distribution after Leiden+LLP ordering
- Hardcodes proven-best options: `copy_blocks=true`, `adaptive_copy=true`, `stop_deltas=true` (for delta), `compact_copy=true`, `tight_intervals=true`
- No `split_residual`, `bv_blocks`, `fixwidth_ref`, or `adaptive_deltas` — only features with demonstrated benefit
- Raw bitstream output (same stream structure as greedy, compatible with MGS v3 container)

**Result on CNR-2000**: **2.8621 BPE** — beats greedy best (2.8628 BPE) by 0.0007 BPE and WebGraph BV (2.897 BPE) by 0.035 BPE. All roundtrips verified.

---

## Motivation: Header Overhead Reduction

The greedy encoder's best configuration uses two separate per-vertex mechanisms:

1. **Empty-prefix flag** (1 bit per vertex): `0` = empty, `1` = non-empty
2. **VLC v2 header** (1–10 bits): encodes ref_mode + enc_type + mil for non-empty vertices

Combined cost: 1 bit (empty flag) + conditional VLC header. For the ~97% non-empty vertices on CNR-2000, the average is 1 + ~1.6 = ~2.6 bits/vertex.

CS merges these into a single prefix tree optimized for the actual action distribution:

| Code | Bits | Greedy equivalent | CS savings |
|------|------|-------------------|------------|
| `0` | 1 | empty(1) + VLC v2 ref+delta(1) = 2 bits | 1 bit |
| `10` | 2 | empty(1) + VLC v2 noref+delta(2) = 3 bits | 1 bit |
| `1100` | 4 | empty(1) + VLC v2 noref+int2(3) = 4 bits | 0 bits |
| `1101` | 4 | empty(1) + VLC v2 ref+int2(4) = 5 bits | 1 bit |
| `1110` | 4 | empty(0) = 1 bit | −3 bits |
| `1111`+5b | 9 | empty(1) + VLC v2 escape(10) = 11 bits | 2 bits |

The empty vertex code (`1110`, 4 bits) costs 3 more bits than the greedy's `0` flag (1 bit). However, empty vertices are only ~3% of all vertices, while the 1-bit savings on the ~53% ref+delta vertices more than compensates. The net effect is a modest reduction in average header cost.

---

## Frequency-Optimized Prefix Codes

Based on CNR-2000 vertex action distribution after Leiden+LLP ordering:

```
Code     Bits  Meaning                      Frequency   Weighted bits
0          1   reference + stop_delta        ~53%        0.53
10         2   no-ref + stop_delta           ~36%        0.72
1100       4   no-ref + interval(mil=2)      ~6%         0.24
1101       4   ref + interval(mil=2)         ~1.5%       0.06
1110       4   empty vertex                  ~3%         0.12
1111+5b    9   escape: 1b ref + 4b enc_opt   ~0.5%       0.05
                                              ───         ────
                                              100%        ~1.72 bits/vertex
```

**Weighted average**: 1.72 bits/vertex vs greedy VLC v2 + empty_prefix: ~2.57 bits/vertex.

The prefix code tree structure:

```
         ┌─ 0: ref + stop_delta (53%)
    root─┤
         └─ 1─┬─ 0: noref + stop_delta (36%)
               └─ 1─┬─ 0─┬─ 0: noref + interval(2) (6%)
                     │    └─ 1: ref + interval(2) (1.5%)
                     └─ 1─┬─ 0: empty vertex (3%)
                           └─ 1: escape (0.5%) → 1b ref + 4b enc_opt
```

### Escape Payload

When the prefix is `1111` (escape), 5 additional bits follow:

```
bit: is_reference (0 = none, 1 = reference)
4 bits: enc_opt tag (same as greedy):
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

This covers all 18 non-default `(ref_mode, enc_type, mil)` combinations that fall outside the 5 fast-path codes.

---

## Bitstream Structure

```
[Stream Header (4 bits)] [Optional Index Section] [Per-Vertex Records...]
```

### Stream Header (4 bits)

Identical to greedy:

```
coding_scheme_flag (1 bit):
  0 = :children (stop-delimited)
  1 = :index (degree array precedes data)

encoding_tag (3 bits):
  000 = Fibonacci     ← default
  001 = Zeta
  010 = Elias gamma
  011 = Elias delta
```

### Optional Index Section (Per-Vertex Offset Table)

When `coding_scheme_flag = 1` (`:index` mode), a per-vertex offset table enables
O(1) random access to any vertex's compressed data:

```
entry_width (6 bits): number of bits per offset entry
N+1 entries (each entry_width bits):
  entry[0..N-1]: bit offset to vertex i's data (relative to start of per-vertex section)
  entry[N]:      total bit size of per-vertex section (end marker)
```

Empty vertex detection: if `entry[v] == entry[v+1]`, vertex v has no neighbors — no
CS header is written for that vertex (unlike children mode where `1110` marks empty).

**Index mode encoding changes** (CS):
- Empty vertices are skipped entirely (nothing between their offsets)
- STOP-terminated delta is replaced by count-prefixed delta:
  `varint(count+1)` followed by delta-encoded values
- The CS prefix code tree is unchanged for non-empty vertices

**Two-pass encoding**: All per-vertex data is first written to a temporary buffer,
recording each vertex's start bit offset. Then the offset table + buffered data
are written to the main stream.

### Per-Vertex Record

Each vertex begins with a **CS header** (1–9 bits) that fully determines the vertex's encoding. No separate empty-prefix flag is needed — the empty case is encoded within the prefix tree.

```
CS_header:
  Decodes to (is_empty, ref_mode, enc_type, mil)

IF is_empty:
  Done. Vertex has no neighbors.

IF ref_mode = :reference:
  varint(ref_distance)              // lookback distance in reference window
  adaptive_copy(copy_bitmap)        // 3-way: complement / copy-blocks / bitmap
  [body: residuals encoded per enc_type]

ELSE (ref_mode = :none):
  [body: full neighbor list encoded per enc_type]
```

### Body Encoding

Identical to greedy. The `enc_type` from the CS header determines the encoding:

**Interval+Residual** (`:interval`):
```
varint(num_intervals + 1)
For each interval:
  varint(start_delta)               // zigzag(start - vertex_id) + 1 for first;
                                    // delta from end of prev interval (tight_intervals)
  varint(length - MIL + 1)          // length shifted by MIL
varint(num_residuals + 1)
delta_encoded_residuals             // first: zigzag(v - vertex_id) + 1, then gaps+1
```

**STOP-Terminated Delta** (`:delta`):
```
For each value in sorted order:
  bit '1'                           // continuation flag
  varint(gap)                       // zigzag(first - vertex_id) + 1 for first; gap for rest
bit '0'                             // STOP terminator
```

CS always uses STOP-terminated deltas for `:delta` encoding. Count-prefixed deltas are not used.

**RLE/Hybrid Mix** (`:rle`):
```
varint(count + 1)                   // neighbor count
hybrid_mix_encoded_deltas           // delta-encoded then hybrid-compressed
```

### Copy Position Encoding (3-Way Adaptive)

Compact prefix code (`compact_copy=true`):
```
  0 → complement copy-blocks (1 bit overhead)
  10 → copy-blocks (2 bits overhead)
  11 → raw bitmap (2 bits overhead)
```

The encoder evaluates all three modes per reference vertex and picks the cheapest.

---

## Greedy Search with CS Headers

The `_cs_vertex_search` function evaluates all encoding options for each vertex, using `_cs_header_cost` to include the exact prefix code bit cost in every comparison. This ensures the encoder accounts for the variable header overhead when selecting encoding strategies.

Search space: 9 `(enc_type, mil)` combinations from `ENCODING_OPTIONS`:
- Interval MIL=2,3,4,5 (4 options)
- RLE MIL=2,3,4,5 (4 options)
- Delta (1 option, always with stop_deltas=true)

For each option, both no-reference and reference modes are evaluated, yielding up to 18 candidates per vertex. The candidate with the lowest total bit cost (header + reference distance + copy positions + body) wins.

### Simplifications vs Greedy

| Feature | Greedy | CS |
|---------|--------|----|
| Header | VLC v2 (1–10b) + empty_prefix (1b) | CS prefix code (1–9b, empty integrated) |
| Copy mode | Configurable (bitmap / blocks / adaptive) | Always 3-way adaptive + compact |
| Delta format | Configurable (stop / count / adaptive) | Always stop_deltas for `:delta` |
| Split residual | Optional per-residual enc_opt | Not supported (not beneficial) |
| BV blocks | Optional alternating blocks | Not supported (not beneficial) |
| Fixed-width ref | Optional fixed-width distance | Not supported (varint always) |
| Adaptive deltas | Optional per-vertex stop/count | Not supported (stop always wins) |

---

## Parameters

| Parameter | Default | Description |
|-----------|---------|-------------|
| `integer_encoding` | `:fibonacci` | Integer encoding for all varints |
| `ref_window_size` | 64 | Reference window size (recent vertices) |
| `coding_scheme` | `:children` | Stop-delimited vertex encoding |
| `compact_copy` | `true` | Compact prefix code for 3-way adaptive copy |
| `tight_intervals` | `true` | Interval start delta from end of prev interval |

### Hardcoded Options (Always On)

| Feature | Value | Rationale |
|---------|-------|-----------|
| `copy_blocks` | `true` | Proven best for reference encoding |
| `adaptive_copy` | `true` | 3-way: bitmap/copy-blocks/complement |
| `stop_deltas` | `true` (for `:delta`) | Eliminates count prefix, saves ~0.08 BPE |
| Empty vertex handling | Integrated in CS header | No separate flag needed |

---

## API

### Encoding

```julia
using Adjacently.Compression: write_cmdstream_graph_data

io = IOBuffer()
bw = BitWriter(io)
write_cmdstream_graph_data(bw, neighbor_lists, :children, 64;
    integer_encoding=:fibonacci, compact_copy=true, tight_intervals=true)
flush_bitwriter(bw; flush_last_bits=true)
bytes = take!(io)
```

### Decoding

```julia
using Adjacently.Compression: read_cmdstream_graph_data

r = BitReader(IOBuffer(bytes))
decoded = read_cmdstream_graph_data(r, T(num_vertices), :children, T;
    integer_encoding=:fibonacci, compact_copy=true, tight_intervals=true,
    ref_window_size=64)
```

---

## Benchmark: CNR-2000

**Graph**: 325,557 vertices, 3,216,152 edges (directed web graph, Italian web crawl).
**Ordering**: Leiden K=1 + per-group LLP (`:sym`, 5 passes).

### Results

| Config | BPE | Size (bytes) | Roundtrip |
|--------|-----|-------------|-----------|
| **Command Stream** | **2.8621** | 1,150,608 | OK |
| Greedy best (vlc2+empty+compact+tight) | 2.8628 | 1,150,904 | OK |
| WebGraph BV | 2.897 | 1,164,843 | — |
| CGE (K=1, 3-way adaptive) | 2.4341 | 978,552 | OK |

### Header Savings Analysis

CS saves 296 bytes (0.0007 BPE) over greedy best. The theoretical savings of ~0.85 bits/vertex (277 Kbit = ~34 KB) from the prefix code improvement are partially offset by:

1. **Empty vertex overhead**: CS uses 4 bits for empty vertices vs greedy's 1 bit (empty_prefix flag). With ~9,767 empty vertices (3%), this costs +29 Kbit (+3.6 KB).
2. **Escape path savings less than expected**: The escape code is 9 bits (CS) vs 11 bits (greedy: 1 empty_flag + 10 VLC escape), saving 2 bits. But escapes are only ~0.5% of vertices.
3. **Vertex search cost model alignment**: The CS header cost differences are small (1–3 bits) and occasionally cause the encoder to make different ref/no-ref decisions than greedy, which can go either way.

The net result is a very modest improvement, indicating that VLC v2 + empty_prefix was already close to optimal for this action distribution.

### Comparison with Greedy Best

| Aspect | Command Stream | Greedy Best |
|--------|---------------|-------------|
| Per-vertex header | 1–9 bits (CS prefix code) | 1 bit (empty) + 1–10 bits (VLC v2) |
| Avg header bits/vertex | ~1.72 | ~2.57 |
| Empty handling | Integrated (4 bits) | Separate flag (1 bit) |
| Best action (ref+delta) | 1 bit | 1 bit (empty) + 1 bit (VLC) = 2 bits |
| Copy mode | Always 3-way adaptive compact | Configurable |
| Delta mode | Always stop_deltas | Configurable |
| Parameters | 5 (simplified) | 14 (full flexibility) |
| **BPE** | **2.8621** | **2.8628** |

### Comparison with WebGraph BV

| Aspect | Command Stream | WebGraph BV |
|--------|---------------|-------------|
| Outdegrees | Implicit (stop-delimited / CS header) | Explicit (0.516 BPE) |
| Per-vertex header | CS prefix (1–9 bits, avg ~1.7) | None (fixed pipeline) |
| Encoding selection | Per-vertex greedy cost search (9 options) | Fixed: intervals → residuals |
| Reference window | w=64 | w=7 |
| Integer encoding | Fibonacci | Zeta-3 |
| Copy positions | 3-way adaptive (compact prefix) | Alternating block lengths |
| **BPE** | **2.8621** | **2.897** |

---

## Key Findings

1. **VLC v2 + empty_prefix is near-optimal for this distribution.** The CS prefix code achieves only 0.0007 BPE improvement over the two-mechanism approach. The dominant action (ref+delta at 53%) already gets a 1-bit code in VLC v2; CS matches this. The main CS advantage — saving 1 bit on noref+delta (36% of vertices) — is partially offset by the higher cost for empty vertices (4 bits vs 1 bit).

2. **Merging empty flag into the prefix tree has mixed value.** Empty vertices are rare (~3%) but cheap to encode with a dedicated 1-bit flag. Giving them a 4-bit prefix code in the unified tree costs 3 extra bits × ~10K vertices = ~30 Kbit, which is close to the savings from shorter codes on the common actions.

3. **Hardcoding best-known options is cleaner.** CS eliminates 9 boolean parameters (`copy_blocks`, `adaptive_copy`, `stop_deltas`, `empty_prefix`, `vlc2`, `adaptive_deltas`, `split_residual`, `bv_blocks`, `fixwidth_ref`) by hardcoding the proven-best settings. The simpler API is easier to use and maintains equivalent compression quality.

4. **The 0.43 BPE gap to CGE is structural, not header-related.** Improving headers from ~2.57 to ~1.72 bits/vertex saves only ~0.085 BPE at most. The remaining gap comes from CGE's per-cluster local vertex indexing and hierarchical encoding, which create fundamentally tighter compression.

5. **Shannon entropy of the action distribution is ~1.6 bits.** The CS prefix code achieves ~1.72 bits/vertex, within 0.12 bits of the theoretical minimum. Further header improvements (e.g., arithmetic coding) could save at most ~0.004 BPE — negligible.

---

## Implementation

### Source Files

| File | Description |
|------|-------------|
| `src/compression.jl` | `_cs_header_cost`, `_write_cs_header`, `_read_cs_header`, `_cs_vertex_search`, `write_cmdstream_graph_data`, `read_cmdstream_graph_data` |
| `scripts/cmdstream_benchmark_cnr2000.jl` | Benchmark script comparing CS vs greedy best |

### Functions

| Function | Description |
|----------|-------------|
| `_cs_header_cost(ref_mode, enc_type, mil, is_empty)` | Returns exact bit cost for CS prefix code (1, 2, 4, or 9 bits) |
| `_write_cs_header(w, ref_mode, enc_type, mil, is_empty)` | Writes CS prefix bits to bitstream |
| `_read_cs_header(r)` | Reads CS prefix, returns `(is_empty, ref_mode, enc_type, mil)` |
| `_cs_vertex_search(neighbors, neighbor_lists, ref_window, ie; ...)` | Greedy vertex search using CS header costs |
| `write_cmdstream_graph_data(w, neighbor_lists, coding_scheme, ref_window_size; ...)` | Encode graph with CS headers |
| `read_cmdstream_graph_data(r, vs, coding_scheme, T; ...)` | Decode CS-encoded graph |

### Reused Functions (Unchanged)

All body encoding, reference search, and copy position functions are shared with the greedy encoder:
- `_try_find_reference`, `_estimate_base_cost`, `_estimate_adaptive_copy_cost`
- `_write_adaptive_copy`, `_read_adaptive_copy`
- `_write_stop_delta`, `_read_stop_delta`, `_estimate_stop_delta_cost`
- `write_intervals_and_residuals`, `read_intervals_and_residuals`
- `write_encoded_value`, `read_encoded_value`, `estimate_encoded_value_cost`
- `reconstruct_from_reference`, `delta_encode_vector`, `write_hybrid_mix_encoded_list`
- `_write_encoding_tag`, `_read_encoding_tag`, `_write_enc_opt_tag`, `_read_enc_opt_tag`
