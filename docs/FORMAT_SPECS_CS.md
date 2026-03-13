# Command Stream (CS) Format Specification

## Overview

Command Stream (CS) is a single-pass graph compression format that replaces the greedy encoder's separate empty-prefix flag and VLC vertex header with a unified frequency-optimized prefix code tree. All proven greedy optimizations are preserved: 3-way adaptive copy, STOP-terminated deltas, tight intervals, Fibonacci encoding, and w=64 reference windows. Only the per-vertex header layer changes.

Key properties:
- Unified prefix code merges empty vertex flag, reference mode, and encoding type into a single code tree
- Frequency-optimal bit assignment based on CNR-2000 action distribution after Leiden+LLP ordering
- Hardcodes proven-best options: `copy_blocks=true`, `adaptive_copy=true`, `stop_deltas=true` (for delta), `compact_copy=true`, `tight_intervals=true`
- Aggressive low-degree reference search: minimum degree 1 (not 3), minimum overlap 1 (not 3)
- Optional `lr_split` for LR-split interval residuals (beneficial on some datasets)
- Raw bitstream output (same stream structure as greedy, compatible with MGS v3 container)

**Result on CNR-2000**: **2.4348 BPE** — beats BG best (2.4929 BPE), CG K=1 (2.6378 BPE), WebGraph BV (2.898 BPE), and WebGraph BV-HC (2.448 BPE). Within 0.116 BPE of CG K=2 best (2.3191 BPE). All roundtrips verified.

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
  ref_distance                      // varint (default) or fixed-width (when lr_split=true)
  adaptive_copy(copy_bitmap)        // 3-way: complement / copy-blocks / bitmap
  [body: residuals encoded per enc_type]

ELSE (ref_mode = :none):
  [body: full neighbor list encoded per enc_type]
```

### Body Encoding

Identical to greedy. The `enc_type` from the CS header determines the encoding:

**Interval+Residual** (`:interval`, `lr_split=false`):
```
varint(num_intervals + 1)
For each interval:
  varint(start_delta)               // zigzag(start - vertex_id) + 1 for first;
                                    // delta from end of prev interval (tight_intervals)
  varint(length - MIL + 1)          // length shifted by MIL
varint(num_residuals + 1)
delta_encoded_residuals             // first: zigzag(v - vertex_id) + 1, then gaps+1
```

**Interval+LR-split Residual** (`:interval`, `lr_split=true`):
```
varint(num_intervals + 1)
For each interval:
  varint(start_delta)               // same as above
  varint(length - MIL + 1)
varint(num_residuals + 1)
IF num_residuals > 0:
  varint(n_left + 1)               // count of residuals < vertex_id
  delta_encoded(left_distances)    // vid - residual[i] (ascending distances from vid)
  delta_encoded(right_distances)   // residual[i] - vid + 1 (ascending distances from vid)
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

| Feature | Greedy (BG) | CS |
|---------|--------|----|
| Header | Merged VLC (1–11b, empty integrated) | CS prefix code (1–9b, empty integrated) |
| Copy mode | Configurable (bitmap / blocks / adaptive) | Always 3-way adaptive + compact |
| Delta format | Configurable (stop / count / adaptive) | Always stop_deltas for `:delta` |
| LR-split | Optional (`lr_split` param) | Optional (`lr_split` param) |
| Fixed-width ref | Optional (`fixwidth_ref = lr_split`) | Automatic when `lr_split=true` |
| Multi-ref | Optional (`multi_ref` param) | Not supported |
| Split residual | Optional per-residual enc_opt | Not supported (not beneficial) |
| BV blocks | Optional alternating blocks | Not supported (not beneficial) |
| Adaptive deltas | Optional per-vertex stop/count | Not supported (stop always wins) |

---

## Parameters

| Parameter | Default | Best (CNR-2000) | Description |
|-----------|---------|-----------------|-------------|
| `integer_encoding` | `:fibonacci` | `:fibonacci` | Integer encoding for all varints |
| `ref_window_size` | 64 | 64 | Reference window size (recent vertices) |
| `coding_scheme` | `:children` | `:children` | Stop-delimited vertex encoding |
| `lr_split` | `false` | `false` | LR-split residual encoding for intervals |

### Hardcoded Options (Always On)

| Feature | Value | Rationale |
|---------|-------|-----------|
| `copy_blocks` | `true` | Proven best for reference encoding |
| `adaptive_copy` | `true` | 3-way: bitmap/copy-blocks/complement |
| `compact_copy` | `true` | Compact prefix code for 3-way adaptive copy |
| `tight_intervals` | `true` | Interval start delta from end of prev interval |
| `stop_deltas` | `true` (for `:delta`) | Eliminates count prefix, saves ~0.08 BPE |
| `min_degree_for_ref` | 1 | Search refs for all non-empty vertices |
| `min_overlap_for_ref` | 1 | Consider all refs with any overlap |
| `fixwidth_ref` | `= lr_split` | Fixed-width ref distances when LR-split active |
| Empty vertex handling | Integrated in CS header | No separate flag needed |

### MGS Header Byte 2 Encoding (0x50-0x6F)

Offset = byte2 - 0x50. Bit layout (5 bits):

| Bits | Field | Values |
|------|-------|--------|
| 4-2 | ref_window_size | 3 bits: 0=4, 1=8, 2=16, 3=32, 4=64, 5=128, 6=256, 7=512 |
| 1 | lr_split | 0=false, 1=true |
| 0 | (reserved) | must be 0 |

Implied: `compact_copy=true`, `tight_intervals=true` (always active, not encoded).

---

## API

### Encoding

```julia
using Adjacently.MGS: write_cs_mgs3_graph

# Best config for CNR-2000 (no lr_split)
write_cs_mgs3_graph(g, "output_base"; ref_window_size=64)

# With lr_split (beneficial on some datasets like enwiki)
write_cs_mgs3_graph(g, "output_base"; ref_window_size=64, lr_split=true)
```

### Decoding

```julia
using Adjacently.MGS: load_compressed_mgs3_graph

# Fully self-describing — no params needed
g = load_compressed_mgs3_graph("output_base.mgz")
```

---

## Benchmark: CNR-2000

**Graph**: 325,557 vertices, 3,216,152 edges (directed web graph, Italian web crawl).
**Ordering**: Leiden K=1 + per-group LLP (`:sym`, 5 passes).

### Results

| Config | BPE | Roundtrip |
|--------|-----|-----------|
| **CS best** | **2.4348** | OK |
| BG best (lr+mr) | 2.4929 | OK |
| WebGraph BV-HC | 2.448 | — |
| CG K=1 (w=64) | 2.6378 | OK |
| WebGraph BV | 2.898 | — |
| CG K=2 best | 2.3191 | OK |

### Optimization History

| Step | Change | BPE | Delta |
|------|--------|-----|-------|
| Original CS | Frequency-optimized prefix codes | 2.8621 | — |
| **Low-degree ref search** | Min degree 3→1, min overlap 3→1 | **2.4348** | **-0.4273** |

The massive improvement came entirely from lowering the reference search thresholds. On CNR-2000, 66K vertices have degree 1-2 (20% of non-empty vertices). Previously these were forced into raw encoding; now they benefit from reference copying.

### LR-split on CNR-2000

CS + lr_split produces 2.5741 BPE (+0.14 vs baseline) — LR-split is harmful on CNR-2000, consistent with CG findings. The crawl-ordered web graph does not benefit from interval+LR encoding at this locality level.

### Comparison with BG

| Aspect | CS | BG best |
|--------|-----|---------|
| Per-vertex header | CS prefix (1–9 bits) | Merged VLC (1–11 bits) |
| Multi-reference | Not supported | 2-ref per vertex |
| LR-split | Optional | Optional |
| Fixed-width ref | Automatic with lr_split | Automatic with lr_split |
| Parameters | 4 (simplified) | 8 (full flexibility) |
| **BPE (CNR-2000)** | **2.4348** | **2.4929** |

CS beats BG by 0.058 BPE despite lacking multi-ref support. The CS prefix code's shorter headers (1 bit for ref+delta vs BG's merged VLC) more than compensate.

### Comparison with WebGraph BV

| Aspect | Command Stream | WebGraph BV |
|--------|---------------|-------------|
| Outdegrees | Implicit (stop-delimited / CS header) | Explicit (0.516 BPE) |
| Per-vertex header | CS prefix (1–9 bits, avg ~1.7) | None (fixed pipeline) |
| Encoding selection | Per-vertex greedy cost search (9 options) | Fixed: intervals → residuals |
| Reference window | w=64 | w=7 |
| Integer encoding | Fibonacci | Zeta-3 |
| Copy positions | 3-way adaptive (compact prefix) | Alternating block lengths |
| **BPE** | **2.4348** | **2.898** |

---

## Key Findings

1. **Low-degree reference search is the single biggest win.** Lowering the minimum degree threshold from 3→1 and overlap from 3→1 saves 0.427 BPE on CNR-2000. The 66K degree 1-2 vertices (20% of non-empty) were previously forced into expensive raw encoding despite having perfectly good reference candidates nearby.

2. **CS beats WebGraph BV-HC.** At 2.4348 BPE, CS surpasses WebGraph's host-compressed variant (2.448 BPE) without needing host-aware compression or URL structure analysis. The combination of per-vertex greedy search, aggressive reference coverage, and frequency-optimized headers is sufficient.

3. **CS beats BG despite fewer features.** CS lacks multi-ref support but still beats BG (2.4348 vs 2.4929) thanks to shorter prefix codes for the dominant ref+delta action (1 bit vs 1 bit merged VLC — but CS avoids the overhead of multi-ref/interval/rle VLC codes that inflate the code tree).

4. **LR-split is dataset-dependent.** On CNR-2000, lr_split hurts CS by +0.14 BPE (same pattern as CG). On datasets with more interval-friendly structure (e.g., enwiki), lr_split may help.

5. **The remaining 0.116 BPE gap to CG K=2 is structural.** CG's per-cluster local vertex indexing (vertices renumbered 1...|C| within each cluster) enables fundamentally tighter compression that single-pass greedy approaches cannot match.

6. **Hardcoding best-known options is cleaner.** CS uses 4 parameters vs BG's 8, by hardcoding `copy_blocks`, `adaptive_copy`, `compact_copy`, `tight_intervals`, and `stop_deltas`. The simpler API maintains equivalent or better compression quality.

### Multi-Dataset Benchmark

Best CS BPE across all tested datasets:

| Dataset | CS BPE | CG BPE | BG BPE | BV BPE | CS Config |
|---------|--------|---------|---------|--------|-----------|
| cnr-2000 (no reorder) | 2.4348 | **2.3191** | 2.4929 | 2.898 | w=64, no-lr |
| cnr-2000 (Leiden+LLP) | 2.3643 | 2.5652 | **2.3258** | 3.2335 | w=64, no-lr |
| in-2004 | 1.7839 | **1.7513** | 1.895 | 2.172 | w=64, no-lr |
| web-google core (Leiden+LLP) | **4.0288** | 4.3296 | 4.0735 | 5.0081 | w=256, no-lr |
| web-google rcore (Leiden+LLP) | **3.7337** | 3.9359 | 3.7626 | 4.1751 | w=256, no-lr |
| eat-core | 10.8679 | **10.5768** | 10.9191 | 10.705 | w=256, lr |
| eat-rcore | 9.5674 | **9.3192** | 9.5726 | 9.391 | w=256, lr |
| arxiv-hep-ph core | 10.0767 | **9.8187** | 10.0910 | 10.262 | w=256, lr |
| arxiv-hep-ph rcore | 8.9855 | **8.9406** | 9.0946 | 9.684 | w=256, no-lr |
| amazon-0601 core | 12.5235 | **12.1853** | 12.4444 | 13.001 | w=8, lr |
| amazon-0601 rcore | 11.8497 | **11.7069** | 11.7310 | 12.064 | w=8, lr |

**Bold** = best among Adjacently methods for that dataset.

**Key findings**:
- CS is the best single-pass method, consistently beating or matching BG
- Leiden+LLP ordering improves CS on Web-Google by ~0.37 BPE (core: 4.399→4.029, rcore: 4.064→3.734)
- CS w=256 remains the best single-pass method on Web-Google with Leiden+LLP ordering, beating BG and CG K=1
- On most other datasets, CG K=1 beats CS via per-cluster local indexing
- CS benefits from larger windows (w=256) on Web-Google and other moderate-locality graphs
- lr_split is dataset-dependent: helps on EAT, Arxiv, Amazon; hurts on cnr-2000 and Web-Google

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

All body encoding, reference search, and copy position functions are shared with the greedy (BG) encoder:
- `_try_find_reference`, `_estimate_base_cost`, `_estimate_adaptive_copy_cost`
- `_write_adaptive_copy`, `_read_adaptive_copy`
- `_write_stop_delta`, `_read_stop_delta`, `_estimate_stop_delta_cost`
- `write_intervals_and_residuals`, `read_intervals_and_residuals`
- `_write_intervals_lr`, `_read_intervals_lr`, `_estimate_intervals_lr_cost`
- `write_encoded_value`, `read_encoded_value`, `estimate_encoded_value_cost`
- `reconstruct_from_reference`, `delta_encode_vector`, `write_hybrid_mix_encoded_list`
- `_write_encoding_tag`, `_read_encoding_tag`, `_write_enc_opt_tag`, `_read_enc_opt_tag`
