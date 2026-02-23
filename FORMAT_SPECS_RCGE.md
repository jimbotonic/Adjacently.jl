# RCGE Format Specification

## Overview

RCGE (Reversible Coarsening Graph Encoding) compresses directed graphs by combining two-level community-based vertex ordering with per-cluster reference-based edge encoding. The key insight is that applying Leiden community detection followed by per-cluster LLP reordering creates tight sequential locality — consecutive vertices in encoding order share 80–90% of neighbors — which makes reference compression dramatically effective.

Key properties:
- Two-step vertex ordering: Leiden fine clustering → per-cluster LLP reordering on induced subgraphs
- Intra-cluster edges compressed with reference encoding (copy-blocks + STOP-terminated additions)
- Fixed-width reference delta encoding (no per-cluster bitmap or varint overhead)
- Cluster membership encoded with Elias-Fano, or eliminated entirely via implicit-ranges pre-relabeling
- Inter-cluster edges trivially small at K=1 (241 bytes on CNR-2000)
- Raw bitstream output (not wrapped in MGZ container)

**Best result on CNR-2000**: **2.550 BPE** — beats WebGraph BV (2.897 BPE) by 0.347 BPE with all roundtrips verified.

The 2.550 BPE result combines **implicit membership** encoding (7 bytes vs 88,764 bytes for Elias-Fano, saves 0.221 BPE) with **adaptive additions** (per-vertex choice between STOP-delta and interval encoding, saves 0.067 BPE) and **adaptive raw** (same per-vertex adaptive choice for non-referenced vertices, saves 0.049 BPE). Each optimization is independent and their savings are additive.

---

## Algorithm Workflow

```
Input: directed graph G, parameter K (default K=1)

1. Apply Leiden community detection to G
   → fine clusters F_1..F_n (typically ~34K clusters of ~9 vertices for CNR-2000)

2. Assign fine clusters to K top-level groups + overflow bin
   → top-K groups by Leiden cluster count; all remaining → group K+1
   → total effective groups: K+1

3. For each top-level group g = 1..K+1:
   a. Extract induced subgraph G_g (only intra-group edges)
   b. Apply LLP reordering on G_g with :sym criterion
   → vertices within each group are ordered so that Leiden fine-cluster
     members are consecutive and locally LLP-optimal

4. Encode all groups as one RCGE level:
   a. Write cluster membership (K+1 sorted vertex lists, Elias-Fano)
   b. Write intra-cluster induced edges with reference compression
   c. Write inter-cluster edge bundles (K+1 groups × K+1 groups)

5. Coarsen: build K+1-node aggregate graph
   → stopping condition always met for K ≤ 31 (K+1 ≤ 32 = min_clusters)

6. Concatenate bitstream (single level in practice for K=1)
```

### Why K=1 is Optimal

With K=1, there are 2 effective groups. The membership cost with Elias-Fano is ~0.221 BPE on CNR-2000. With implicit-ranges encoding (pre-relabeled graph, see below), membership becomes 7 bytes ≈ 0 BPE. Inter-cluster is negligible at K=1 (241 bytes), and intra encoding quality is unchanged. Going from K=8 to K=1 saves 0.104 BPE with Elias-Fano; combined with implicit ranges this gives a total 0.325 BPE improvement over K=8.

### Why Per-Cluster LLP Dominates

The two-step ordering (Leiden + per-cluster LLP on induced subgraphs) is the single most important factor. When LLP runs on a 9-vertex Leiden fine cluster's induced subgraph, consecutive vertices share 80–90% of neighbors. By contrast, LLP applied to the full graph yields much weaker locality. Benchmark evidence:

| Approach | Window | BPE |
|----------|--------|-----|
| No reordering (original IDs) | 7 | ~4.3 |
| Global LLP (whole graph) | 7 | 3.561 |
| Global LLP (whole graph) | 64 | 3.319 |
| **RCGE (Leiden + per-cluster LLP)** | **64** | **2.887** |

The 0.432 BPE advantage of RCGE over global LLP + w=64 + Fibonacci is entirely from ordering quality, not codec choice.

---

## Parameters

`RCGEParams` controls all encoding decisions:

| Parameter | Default | Description |
|-----------|---------|-------------|
| `L` | 32 | Cluster size threshold: ≤ L and undirected → upper-triangular bitset; else per-vertex lists |
| `varint` | `:fibonacci` | Integer encoding for sizes/lengths (positive-only fields) |
| `count_varint` | `:fibonacci` | Encoding for counts that may be zero |
| `gap` | `:fibonacci` | Encoding for gap-coded sorted lists |
| `degree` | `:golomb` | Encoding for degree vectors |
| `perm_strategy` | `:lehmer` | Permutation encoding strategy (`:lehmer`, `:raw`, `:blockpos`) |
| `undirected_pairs` | `true` | If true and graph is undirected, only encode cluster pairs A < B |
| `membership` | `:stop` | Membership encoding: `:elias_fano`, `:delta`, `:stop`, or `:implicit_ranges` |
| `inter_strategy` | `:lists` | Inter-cluster encoding: `:perm` or `:lists` |
| `intra_ref_enabled` | `true` | Enable reference encoding for intra-cluster adjacency |
| `intra_ref_window` | 16 | Lookback window size for reference candidates (best: 64) |
| `intra_ref_min_overlap` | 0.3 | Minimum overlap fraction for reference to be considered |
| `intra_ref_rle` | `true` | Use RLE ones-deltas for reference delta sequences (legacy; `intra_ref_fixwidth` is better) |
| `intra_ref_fixwidth` | `false` | Fixed-width per-vertex ref delta encoding (1 flag bit + ceil(log₂(window)) bits) |
| `intra_ref_vlc` | `false` | VLC (Fibonacci) for ref delta instead of fixed-width; harmful in practice (+0.093 BPE) |
| `intra_add_adaptive` | `false` | Per-vertex adaptive additions: pick cheaper of STOP-delta vs interval+residuals (requires `intra_stop_deltas=true`) |
| `intra_raw_adaptive` | `false` | Per-vertex adaptive raw: pick cheaper of STOP-delta vs interval+residuals for non-referenced vertices |
| `intra_adapt_mil` | `2` | Minimum interval length for adaptive interval encoding; MIL=2 is most effective |
| `intra_block_try` | `false` | Try block encoding (full MGS-style) for each cluster |
| `positions_mode` | `:delta` | Reference positions encoding: `:delta`, `:bitmap`, `:auto` |
| `additions_mode` | `:auto` | Reference additions encoding: `:auto`, `:delta`, `:intervals` |
| `min_cluster_density` | 0.0 | Minimum density threshold for clusters |
| `intra_intervals` | `false` | Use MGS-style interval+residual encoding for neighbor lists |
| `intra_mil` | 4 | Minimum interval length for interval detection |
| `intra_greedy_mil` | `false` | Per-vertex greedy search over mil values {2,3,4,5} |
| `intra_mgs` | `false` | Use full MGS encoder (reference+interval+recursive) per cluster |
| `intra_zigzag` | `false` | Zigzag relative first-value encoding (offset from local vertex index) |
| `intra_stop_deltas` | `false` | STOP-terminated delta lists (eliminates per-vertex count prefix) |
| `intra_copy_blocks` | `false` | WebGraph-style copy-blocks for reference positions (vs RLE bitmap) |

**Recommended best config** (2.550 BPE on CNR-2000 with K=1 + implicit ranges + adaptive encoding):
```julia
# Step 1: Leiden K=1 partition → 2 groups, then sort each by vertex ID
clusters_by_id = [sort(group1_vertices), sort(group2_vertices)]
# Step 2: Relabel vertices by vertex-ID rank within each group
g_rel = relabel_graph(g, vertex_map)  # vertex_map: v → rank in ID-sorted group
clusters_contiguous = [TV.(1:S1), TV.(S1+1:n)]
# Step 3: Encode with implicit_ranges + adaptive intervals
RCGEParams(
    L=128, membership=:implicit_ranges, inter_strategy=:perm,
    intra_ref_window=64, intra_ref_rle=false, intra_ref_fixwidth=true,
    intra_zigzag=true, intra_stop_deltas=true, intra_copy_blocks=true,
    intra_add_adaptive=true, intra_raw_adaptive=true, intra_adapt_mil=2,
    varint=:fibonacci, degree=:elias_delta, perm_strategy=:blockpos,
    undirected_pairs=false,
)
```

**Intermediate config without adaptive** (2.666 BPE, simpler — no per-vertex mode bits):
```julia
RCGEParams(
    L=128, membership=:implicit_ranges, inter_strategy=:perm,
    intra_ref_window=64, intra_ref_rle=false, intra_ref_fixwidth=true,
    intra_zigzag=true, intra_stop_deltas=true, intra_copy_blocks=true,
    varint=:fibonacci, degree=:elias_delta, perm_strategy=:blockpos,
    undirected_pairs=false,
)
```

**Standard config** (2.887 BPE, simpler — no pre-relabeling needed):
```julia
RCGEParams(
    L=128, membership=:elias_fano, inter_strategy=:perm,
    intra_ref_window=64, intra_ref_rle=false, intra_ref_fixwidth=true,
    intra_zigzag=true, intra_stop_deltas=true, intra_copy_blocks=true,
    varint=:fibonacci, degree=:elias_delta, perm_strategy=:blockpos,
    undirected_pairs=false,
)
# with K=1 in the encoding loop
```

---

## Bitstream Structure

Each coarsening level produces a bitstream with three sections:

```
┌──────────────────────────────────────────────┐
│ SECTION 1: Cluster Membership                │
│   K (varint): number of clusters             │
│   For each cluster c = 1..K:                 │
│     Sorted vertex list (Elias-Fano preferred)│
├──────────────────────────────────────────────┤
│ SECTION 2: Intra-Cluster Edges               │
│   For each cluster c = 1..K:                 │
│     IF |c| <= L AND undirected:              │
│       Upper-triangular bitset                │
│     ELSE:                                    │
│       [Per-vertex ref header]                │
│       Per-vertex adjacency lists             │
├──────────────────────────────────────────────┤
│ SECTION 3: Inter-Cluster Edges               │
│   For each source cluster A = 1..K:          │
│     STOP-terminated target cluster list      │
│     For each target cluster B:               │
│       Per-source-vertex neighbor records     │
│       STOP bit (0)                           │
└──────────────────────────────────────────────┘
```

### Section 1: Cluster Membership

The number of clusters K is written first as a varint. Then for each cluster, the sorted list of global vertex IDs is written.

**Elias-Fano** (`membership=:elias_fano`, recommended for no pre-relabeling): Self-describing — encodes its own length. Most compact for large sorted integer lists. At K=1 (2 groups), cost ≈ 0.221 BPE on CNR-2000.

**Delta** (`membership=:delta`): Explicit length varint followed by gap-coded vertex IDs.

**STOP-terminated** (`membership=:stop`): No explicit length. Each value preceded by `1` bit and delta-encoded; terminated by `0` bit.

**Implicit ranges** (`membership=:implicit_ranges`, best for pre-relabeled graphs): Requires vertices to be pre-relabeled so each cluster occupies a contiguous range of vertex IDs. Writes only the cluster sizes as varints (e.g., [S₁, S₂] = 7 bytes total for K=1). The decoder reconstructs cluster i as vertices `(offset+1)..(offset+S_i)` where offset = sum of previous sizes.

**Pre-relabeling requirement for implicit ranges**: The graph must be relabeled so that cluster i occupies IDs `(S₁+…+S_{i-1}+1)..(S₁+…+S_i)`. Critically, the relabeling must use **vertex-ID rank** within each cluster (not LLP rank), because `encode_level` re-sorts clusters by vertex ID internally. This preserves bit-for-bit identical intra encoding compared to Elias-Fano, changing only the membership encoding.

```julia
# Correct pre-relabeling for implicit ranges:
clusters_sorted = [sort(group1), sort(group2)]  # ID-sorted, not LLP-sorted
vertex_map = Dict(v => TV(rank) for (rank, v) in enumerate(Iterators.flatten(clusters_sorted)))
g_rel = relabel_graph(g, vertex_map)
```

### Section 2: Intra-Cluster Edges

For each cluster, the induced subgraph edges are encoded. Two modes:

#### Small Undirected Clusters (|C| ≤ L, undirected graph)

Upper-triangular bitset: for all pairs (i, j) where 1 ≤ i < j ≤ |C|, emit one bit indicating whether edge (i,j) exists. Total bits: |C| × (|C|−1) / 2.

#### Large or Directed Clusters

Per-vertex adjacency lists with optional reference encoding.

**Cluster-level reference header** (one of two formats):

*Fixed-width format* (`intra_ref_fixwidth=true`, recommended): Per vertex in order:
```
For each vertex v = 1..|C|:
  bit: has_ref flag (1 = uses reference, 0 = raw encoding)
  If has_ref:
    ceil(log₂(window)) bits: (delta − 1) as big-endian fixed-width integer
    where delta ∈ [1, window] is the lookback distance to the reference vertex
```
No byte-padding, no count prefix, no varint overhead. Saves ~0.095–0.126 BPE over legacy format.

*Legacy bitmap+delta format* (`intra_ref_fixwidth=false`):
1. **Reference bitmap**: byte-padded packed boolean vector, 1 bit/vertex
2. **Reference deltas**: count varint + individual Fibonacci varints (if `intra_ref_rle=false`), or RLE ones-deltas (if `intra_ref_rle=true`)

**Per-vertex payload** (in cluster vertex order):

For **referenced vertices** (`has_ref=true`):
```
positions:
  IF intra_copy_blocks:
    copy-blocks encoding (see Copy-Blocks below)
  ELSE:
    RLE ones-deltas bitmap over the reference vertex's neighbor list

additions:
  IF intra_add_adaptive AND intra_stop_deltas:
    bit: mode flag (1 = intervals, 0 = STOP-delta)
    IF mode = intervals:
      intervals_and_residuals encoding [†]
    ELSE:
      STOP-terminated delta list [†]
  ELIF intra_stop_deltas:
    STOP-terminated delta list [†]
  ELIF intra_intervals OR intra_greedy_mil:
    intervals_and_residuals encoding [†]
  ELSE (additions_mode = :delta):
    small_count(num_additions)
    delta-encoded addition list [†]
```

For **non-referenced vertices** (`has_ref=false`):
```
IF intra_raw_adaptive AND intra_stop_deltas:
  bit: mode flag (1 = intervals, 0 = STOP-delta)
  IF mode = intervals:
    intervals_and_residuals encoding [†]
  ELSE:
    STOP-terminated delta list [†]
ELIF intra_stop_deltas:
  STOP-terminated delta list [†]
ELIF intra_intervals OR intra_greedy_mil:
  intervals_and_residuals encoding [†]
ELSE:
  small_count(degree_in_cluster)
  delta-encoded sorted local neighbor IDs [†]
```

`[†]` When `intra_zigzag=true`, the first value is encoded as a zigzag offset from the current vertex's local index (see below). Does not apply to reference position lists.

#### Reference Encoding Decision

For each vertex, the encoder measures:
1. **Raw cost**: delta-coded full neighbor list
2. **Reference cost** (for each candidate in lookback window): copy-blocks positions + additions

The cheapest option wins. Two-pointer merge computes positions (shared elements) and additions (new elements). The lookback window contains at most `intra_ref_window` previous vertices in cluster order.

### Section 3: Inter-Cluster Edges

Organized by source cluster. For each source cluster A:

1. **Target list**: STOP-terminated delta-coded list of target cluster indices B with edges from A.

2. **Per-(A,B) group**: For each target cluster B, active source vertices are written as records:
```
For each source vertex u in A with neighbors in B:
  bit '1'                          # active marker
  varint(local_index_of_u_in_A)    # 1-based position in cluster A
  STOP-terminated delta list of B-local neighbor indices
bit '0'                            # end of AB group
```

At K=1 (2 groups), inter-cluster cost is trivially small — only 241 bytes on CNR-2000 (0.001 BPE).

---

## Encoding Details

### Copy-Blocks Positions Encoding

When `intra_copy_blocks=true`, reference positions are contiguous runs:

```
Copy-blocks format:
  small_count(num_copy_blocks)       # number of contiguous runs
  IF num_copy_blocks > 0:
    varint(first_start)              # 1-based start of first run (≥ 1)
    varint(first_length)             # length of first run (≥ 1)
    For each subsequent block i = 2..num_copy_blocks:
      varint(gap)                    # gap from end of previous block (≥ 1)
      varint(length)                 # length of this run (≥ 1)
```

**Example**: Reference has 10 neighbors; current vertex shares positions {1,2,3,6,7}.
- Copy blocks: [(start=1, len=3), (start=6, len=2)], gap = 6−(1+3) = 2
- Encoded: `small_count(2), fib(1), fib(3), fib(2), fib(2)` = 2+2+3+2+3 = 12 bits

Advantages over RLE ones-deltas bitmap:
1. Skip runs implicit (encoded as gaps) — no per-run type flag
2. No explicit bitmap length (decoder knows ref_len from reference list)
3. All values ≥ 1 — no +1 shifts needed
4. Saves ~0.398 BPE on CNR-2000 (vs old RLE bitmap)

### Fixed-Width Reference Delta Encoding

When `intra_ref_fixwidth=true`, the per-cluster reference header is a single per-vertex stream of fixed-width codes:

```
nbits = ceil(log₂(window))    # e.g., 6 bits for window=64

For each vertex v = 1..|C|:
  write_bit(has_ref[v])        # 1 bit
  if has_ref[v]:
    write bits: (delta[v] − 1) in nbits bits, big-endian
               where delta[v] ∈ [1, window]
```

**Savings over legacy bitmap+varint**:
- Eliminates byte-padding of bitmap (~3.5 bits/cluster × 34K clusters = ~4.5 KB)
- Eliminates delta count prefix (~8 bits/cluster × 34K clusters = ~34 KB)
- Fibonacci of delta values averages 6.5–7.5 bits; fixed 6 bits saves ~1.5 bits per ref vertex
- Total: ~38 KB saved = 0.095–0.126 BPE

At window=64 (nbits=6): non-ref cost = 1 bit, ref cost = 7 bits. Expected cost: 0.38×1 + 0.62×7 = 4.72 bits/vertex × 325K = 193 KB on CNR-2000.

### STOP-Terminated Delta Lists

When `intra_stop_deltas=true`, neighbor lists (raw and additions) use:

```
For each value v in sorted order:
  bit '1'          # more values
  varint(gap)      # v − prev (first value: absolute or zigzag-encoded)
bit '0'            # STOP terminator
```

Eliminates the `small_count` prefix (2–6 bits). Net saving: ~0.051 BPE on CNR-2000. Wins for vertex degree ≤ 5, loses for degree ≥ 7.

### Zigzag Relative First-Value Encoding

When `intra_zigzag=true`, the first value in neighbor lists is encoded as a signed offset from the vertex's local cluster index:

```
offset = v1 − local_idx
zigzag(offset) = offset ≥ 0 ? 2×offset : 2×(−offset) − 1
encoded = zigzag(offset) + 1     # +1 for Fibonacci compatibility
```

Decoding: `v1 = local_idx + zigzag_decode(raw − 1)`

Exploits LLP locality: neighbors cluster near the vertex in encoding order → small offsets → short Fibonacci codes. Saves ~0.741 BPE on CNR-2000 (single largest optimization).

### Adaptive Per-Vertex Codec Selection

When `intra_add_adaptive=true` (for additions) or `intra_raw_adaptive=true` (for raw/non-referenced vertices), the encoder independently decides per vertex which encoding is cheaper, and writes a 1-bit mode flag inline with the payload:

```
bit: mode (1 = interval+residuals, 0 = STOP-delta)
[payload: interval+residuals OR STOP-delta, depending on mode]
```

The encoder tries both encodings on a temporary buffer, compares byte counts (including the mode flag), and chooses the winner. `intra_adapt_mil` sets the minimum interval length threshold (MIL=2 is most aggressive; longer intervals required for MIL=3,4).

**Why adaptive works**: The additions section contains a mix of vertex neighborhoods. Some vertices have neighbors in arithmetic progressions (intervals encode cheaply), others have scattered neighbors (STOP-delta is cheaper). Fixed-mode encoding must pick one globally; adaptive picks optimally per vertex at the cost of 1 extra bit per vertex.

**Implementation note**: Requires `intra_stop_deltas=true` (the STOP-delta branch). The interval branch uses `write_intervals_and_residuals` with the same `vertex_id` zigzag first-value as the delta branch.

### RLE Ones-Deltas Bitmap (Legacy)

When `intra_copy_blocks=false`, positions use dense boolean bitmap:
- `varint(token_count)` alternating runs
- Each token: 1-bit flag (is_ones_run) + `varint(run_length)`

Superseded by copy-blocks, retained for backward compatibility.

---

## Integer Encodings

| Encoding | Zero Support | Description |
|----------|-------------|-------------|
| Fibonacci | No (≥ 1) | Zeckendorf representation terminated by `11` |
| Elias Delta | No (≥ 1) | Elias gamma prefix + binary suffix |
| Elias Fano | Yes | Monotone integer sequence: high/low bit split |
| Golomb | Yes (≥ 0) | Base-b quotient/remainder coding |
| Small-count | Yes | 2-bit escape: 00/01/10 for 0/1/2; 11 + varint for ≥ 3 |

For encodings without zero support (Fibonacci, Elias Delta), zero-valued fields are stored shifted by +1.

---

## Statistics Tracking

`RCGEStats` accumulates per-section bit counts:

| Field | Description |
|-------|-------------|
| `bits_membership` | Bits for Section 1 (cluster membership lists) |
| `bits_intra` | Total bits for Section 2 (intra-cluster edges) |
| `bits_intra_headers` | Bits for per-cluster ref headers (fixed-width or bitmap+deltas) |
| `bits_intra_copy` | Bits for copy-blocks positions |
| `bits_intra_add` | Bits for additions lists |
| `bits_intra_raw` | Bits for non-referenced vertex adjacency data |
| `bits_intra_ref_small_headers` | Bits for small reference headers |
| `intra_ref_used` | Count of vertices using reference encoding |
| `intra_no_ref` | Count of vertices using raw encoding |
| `bits_inter_headers` | Bits for inter-cluster target lists |
| `bits_inter_degrees` | Bits for inter-cluster degree vectors |
| `bits_inter_perms` | Bits for inter-cluster permutation data |
| `bits_inter_lists` | Bits for inter-cluster neighbor lists |

---

## Performance on CNR-2000

**Graph**: 325,557 vertices, 3,216,152 edges (Italian web crawl).
**WebGraph reference**: 2.897 BPE, 1,164,843 bytes.

### Optimization History

Ordering: Leiden fine clustering + per-cluster LLP, K=8 unless noted.

| Configuration | K | BPE | Δ BPE | File Size |
|---------------|---|-----|-------|-----------|
| Baseline (delta, window=32) | 8 | 4.307 | — | 1,691 KB |
| + Zigzag first-value | 8 | 3.566 | −0.741 | 1,400 KB |
| + STOP-terminated deltas | 8 | 3.515 | −0.051 | 1,380 KB |
| + Copy-blocks positions | 8 | 3.117 | −0.398 | 1,224 KB |
| + Fixed-width deltas (window=32) | 8 | 3.023 | −0.094 | 1,187 KB |
| + Fixed-width deltas (window=64) | 8 | 2.991 | −0.032 | 1,174 KB |
| Reduce K=4 | 4 | 2.954 | −0.037 | 1,160 KB |
| Reduce K=3 | 3 | 2.934 | −0.020 | 1,152 KB |
| Reduce K=2 | 2 | 2.914 | −0.020 | 1,144 KB |
| Reduce K=1 | 1 | 2.887 | −0.027 | 1,134 KB |
| WebGraph BV (reference) | — | 2.897 | — | 1,138 KB |
| + Implicit ranges membership | 1 | 2.666 | −0.221 | 1,047 KB |
| + Adaptive additions (MIL=2) | 1 | 2.599 | −0.067 | 1,022 KB |
| + Adaptive raw (MIL=2) | 1 | 2.617 | −0.049 | 1,028 KB |
| **+ Adaptive adds+raw (MIL=2)** | **1** | **2.550** | **−0.116** | **1,001 KB** |

**Total improvement**: 4.307 → 2.550 BPE (−40.8% over 9 independent optimizations).

Note: adaptive adds (−0.067 BPE) and adaptive raw (−0.049 BPE) are independent; combined they give −0.116 BPE because each targets a different section of the bitstream. All roundtrips verified OK.

### Best Config Bit Budget (FW64, K=1, Implicit Ranges + Adaptive — 2.550 BPE)

| Component | Bytes | BPE | Share |
|-----------|-------|-----|-------|
| Membership (2 sizes, implicit) | 7 | 0.000 | 0.0% |
| Headers (fixed-width ref info) | 193,406 | 0.481 | 19.1% |
| Copy positions (copy-blocks) | 276,433 | 0.688 | 27.5% |
| Additions (adaptive: STOP-delta or intervals) | 422,987 | 1.052 | 42.1% |
| Raw (adaptive: STOP-delta or intervals) | 132,121 | 0.329 | 13.2% |
| Inter-cluster edges | 241 | 0.001 | 0.0% |
| **Total** | **1,025,195** | **2.550** | 100% |

62.5% of vertices (203,615 / 325,557) use reference encoding. Adaptive savings: additions −27,015 bytes, raw −19,596 bytes.

### Intermediate Config Bit Budget (FW64, K=1, Implicit Ranges — 2.666 BPE)

| Component | Bytes | BPE | Share |
|-----------|-------|-----|-------|
| Membership (2 sizes, implicit) | 7 | 0.000 | 0.0% |
| Headers (fixed-width ref info) | 193,406 | 0.481 | 18.1% |
| Copy positions (copy-blocks) | 276,433 | 0.688 | 25.9% |
| Additions (STOP delta residuals) | 450,002 | 1.119 | 41.9% |
| Raw (non-referenced vertices) | 151,717 | 0.377 | 14.1% |
| Inter-cluster edges | 241 | 0.001 | 0.0% |
| **Total** | **1,071,806** | **2.666** | 100% |

62.5% of vertices (203,615 / 325,557) use reference encoding.

### Standard Config Bit Budget (FW64, K=1, Elias-Fano — 2.887 BPE)

| Component | Bytes | BPE | Share |
|-----------|-------|-----|-------|
| Membership (2 groups, Elias-Fano) | 88,764 | 0.221 | 7.6% |
| Headers (fixed-width ref info) | 193,406 | 0.481 | 16.7% |
| Copy positions (copy-blocks) | 276,433 | 0.688 | 23.8% |
| Additions (STOP delta residuals) | 450,002 | 1.119 | 38.8% |
| Raw (non-referenced vertices) | 151,717 | 0.377 | 13.1% |
| Inter-cluster edges | 241 | 0.001 | 0.0% |
| **Total** | **1,160,563** | **2.887** | 100% |

62.5% of vertices (203,615 / 325,557) use reference encoding.

### Comparison with WebGraph BV

| Component | RCGE K=1 + Adaptive | RCGE K=1 + Implicit | RCGE K=1 + EF | WebGraph BV |
|-----------|---------------------|---------------------|---------------|-------------|
| Membership / structure | 0.000 BPE | 0.000 BPE | 0.221 BPE | ~0.078 BPE |
| Headers | 0.481 BPE | 0.481 BPE | 0.481 BPE | ~0.126 BPE |
| Copy positions | 0.688 BPE | 0.688 BPE | 0.688 BPE | ~0.556 BPE |
| Residuals/additions | 1.052 BPE | 1.119 BPE | 1.119 BPE | ~1.583 BPE |
| Raw / non-referenced | 0.329 BPE | 0.377 BPE | 0.377 BPE | ~0.554 BPE |
| Inter-cluster | 0.001 BPE | 0.001 BPE | 0.001 BPE | — |
| **Total** | **2.550 BPE** | **2.666 BPE** | **2.887 BPE** | **2.897 BPE** |

RCGE + implicit ranges + adaptive encoding beats WebGraph BV by **0.347 BPE**. The reference encoding (copy + headers) is much more effective in RCGE because the two-step Leiden + per-cluster LLP ordering creates dense locality — 62.5% of vertices can reference a recent neighbor, vs ~56% for WebGraph's bisection ordering. The adaptive per-vertex codec choice (STOP-delta vs interval+residuals) further reduces both the additions and raw sections.

### Greedy Approach Comparison

| Approach | Ordering | Window | BPE |
|----------|----------|--------|-----|
| Greedy RL | None (crawl order) | 7 | ~4.3 |
| Greedy RL | Global LLP | 7 | 3.561 |
| Greedy RL | Global LLP | 64 | 3.319 |
| RCGE (K=1, EF) | Leiden + per-cluster LLP | 64 | 2.887 |
| WebGraph BV | Greedy bisection | 7 | 2.897 |
| RCGE (K=1, implicit ranges) | Leiden + per-cluster LLP | 64 | 2.666 |
| **RCGE (K=1, implicit + adaptive)** | **Leiden + per-cluster LLP** | **64** | **2.550** |

The 0.432 BPE advantage of RCGE over best-case greedy (global LLP + w=64 + Fibonacci) is entirely from the two-step ordering quality. Per-cluster LLP on Leiden-induced subgraphs creates locality an order of magnitude tighter than global LLP.

---

## Key Findings

1. **Vertex ordering is the dominant factor** (+1.43 BPE advantage over unordered greedy). The two-step approach — Leiden community detection creates ~34K tight fine clusters of ~9 vertices each, then per-cluster LLP ordering makes consecutive vertices share 80–90% of neighbors — is what enables RCGE to beat WebGraph.

2. **K=1 is the optimal cluster count** (−0.104 BPE from K=8 to K=1). At K=1, there are 2 effective top-level groups. Membership cost drops from 0.318 to 0.221 BPE (Elias-Fano) or to 0.000 BPE (implicit ranges). Inter-cluster becomes negligible (241 bytes). Intra encoding quality is unaffected — the reference window operates locally within each group's LLP-ordered sequence regardless of group count.

8b. **Implicit ranges membership** (−0.221 BPE). When vertices are pre-relabeled so the two K=1 groups occupy contiguous ID ranges [1..S₁] and [S₁+1..N], membership can be encoded as just two size varints (7 bytes) instead of Elias-Fano sorted lists (88,764 bytes). The key implementation insight: `encode_level` always re-sorts cluster arrays by vertex ID internally, so the pre-relabeling must use **vertex-ID rank** within each group (not LLP rank) to preserve bit-for-bit identical intra encoding. With this approach, the intra section is completely unchanged and only the membership section is eliminated.

3. **Copy-blocks is the largest codec optimization** (−0.398 BPE). WebGraph-style copy-blocks replace the RLE bitmap: skip runs are implicit (gaps between copy blocks), no per-run type flag, no bitmap length prefix. Additionally, the more accurate cost model causes 5K extra vertices to switch to reference encoding.

4. **Zigzag first-value encoding** (−0.741 BPE) exploits the per-cluster LLP locality. First neighbors are typically close to the vertex's own local index; the signed offset is small and compresses well with Fibonacci.

5. **Fixed-width ref deltas** (−0.094 to −0.126 BPE) replace the byte-padded bitmap + Fibonacci varint approach. With window=64, each ref vertex costs exactly 7 bits (1 flag + 6 data) vs ~7.5 bits average with Fibonacci + bitmap padding.

6. **Window=64 outperforms window=32** (−0.032 BPE). More reference candidates → 1,519 additional ref vertices, fewer additions, at the cost of +1 bit per ref delta (6 vs 5 fixed bits). Net positive.

7. **Zeta-3 encoding is harmful for RCGE** (+0.061 to +0.099 BPE vs Fibonacci). Despite being optimal for WebGraph's residuals, zeta-3 degrades copy-blocks positions encoding — it encodes copy-block starts/gaps/lengths less efficiently than Fibonacci for the small values typical in dense clusters.

8. **Global LLP + fixed interval encoding are both harmful**: Global LLP pre-ordering (−0.33 BPE from Leiden quality), fixed interval encoding (+0.119 BPE overhead not amortized on short intra-cluster lists).

9. **Adaptive per-vertex codec selection** (−0.116 BPE combined). With a single inline mode bit per vertex, the encoder picks the cheaper of STOP-terminated delta-list vs interval+residuals encoding. This recovers the benefit of interval encoding for vertices with long arithmetic progressions in their neighbor lists, while avoiding overhead for vertices where delta encoding is cheaper. At MIL=2 (aggressive interval detection): additions −0.067 BPE (450 KB → 423 KB), raw −0.049 BPE (152 KB → 132 KB), total −0.116 BPE.

10. **VLC reference-delta encoding is harmful** (+0.093 BPE). Using Fibonacci variable-length coding for ref deltas instead of fixed 6-bit encoding was expected to help but hurts in practice — ref deltas are distributed uniformly over [1,64] rather than concentrated near 1, so fixed-width (6 bits average) beats Fibonacci (~7–8 bits average for uniform input). This demonstrates that VLC only helps when the distribution is heavily skewed toward small values.
