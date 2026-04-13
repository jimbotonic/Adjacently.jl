# CG Format Specification

## Overview

CG (Clustered Graph Encoding) compresses directed graphs by combining two-level community-based vertex ordering with per-cluster reference-based edge encoding. The key insight is that applying Leiden community detection followed by per-cluster LLP reordering creates tight sequential locality — consecutive vertices in encoding order share 80–90% of neighbors — which makes reference compression dramatically effective.

Key properties:
- Two-step vertex ordering: Leiden fine clustering → per-cluster LLP reordering on induced subgraphs
- Intra-cluster edges compressed with reference encoding (copy-blocks + STOP-terminated additions)
- Fixed-width reference delta encoding (no per-cluster bitmap or varint overhead)
- Left/right residual split after interval extraction for smaller first-value encodings
- Cluster membership encoded with Elias-Fano, or eliminated entirely via implicit-ranges pre-relabeling
- Inter-cluster edges trivially small at K=1 (241 bytes on CNR-2000)
- Raw bitstream output (not wrapped in MGZ container)
- Analytical cost estimation with two-phase overlap-based candidate pruning (no IOBuffer trial encoding)
- Pre-built neighbor lists and double-buffer swap for zero-copy reference search

**Best results**:
- **CNR-2000**: **2.3286 BPE** (K=2, w=64, no LLP, adaptive STOP-delta) — beats WebGraph BV (2.898 BPE) by 0.569 BPE and BV-HC (2.448 BPE) by 0.119 BPE
- **in-2004**: **1.7513 BPE** (K=1, w=8, no clustering) — beats WebGraph BV (2.172 BPE) by 0.421 BPE and BV-HC (1.767 BPE) by 0.016 BPE
- **enwiki-2013**: **12.161 BPE** (BG w=64, lr+mr, Leiden+LLP) — beats BV-HC (12.639) by 3.8%. CG 12.485 (LLP), CS 12.222 (Leiden+LLP). CG original ordering: 15.718
- **Web-Google core**: **4.3296 BPE** (K=1, w=8, Leiden+LLP ordering, no intervals) — CS 4.029, BG 4.074 also benefit from same ordering
- **Web-Google rcore**: **3.9359 BPE** (K=1, w=8, Leiden+LLP ordering, no intervals) — CS 3.734, BG 3.763 also benefit from same ordering
- **EAT core**: **10.5768 BPE** (K=1, w=256, intervals + LR + tight_deltas, gamma)
- **EAT rcore**: **9.3192 BPE** (K=1, w=256, intervals + LR + tight_deltas, gamma)

All roundtrips verified.

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

4. Encode all groups as one CG level:
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
| **CG (Leiden + per-cluster LLP)** | **64** | **2.887** |

The 0.432 BPE advantage of CG over global LLP + w=64 + Fibonacci is entirely from ordering quality, not codec choice.

---

## Parameters

`CGParams` controls all encoding decisions:

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
| `intra_ref_rle` | `true` | Use RLE ones-deltas for reference delta sequences (legacy; `intra_ref_fixwidth` is better) |
| `intra_ref_fixwidth` | `false` | Fixed-width per-vertex ref delta encoding (1 flag bit + ceil(log₂(window)) bits) |
| `intra_copy_adaptive` | `false` | Per-ref-vertex: 3-way adaptive copy positions — pick cheapest of bitmap, copy-blocks, or complement copy-blocks (requires `intra_copy_blocks=true`). Nested 2-bit mode scheme: outer=1→bitmap (1 bit overhead), outer=0+inner=0→copy-blocks (2 bits), outer=0+inner=1→complement (2 bits). |
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
| `intra_lr_split` | `false` | Left/right residual split: after interval extraction, residuals are split at vertex_id into left/right halves and encoded as ascending distances (requires `intra_intervals=true` and `intra_zigzag=true`) |
| `intra_tight_deltas` | `false` | Remove redundant +1 shift on LR residual gaps that are always ≥1 (sorted unique lists). Uses `positive_gaps=true` in delta encoding. Requires `intra_lr_split=true`. Encoded in MGS header as `interval_mode=3`. Saves 0.159 BPE on EAT rcore, 0.064 BPE on EAT core. |

**Recommended config for CNR-2000** (2.3286 BPE, K=2, no LLP, adaptive STOP-delta):
```julia
# K=2 Leiden partition → 2 clusters, pre-relabel for implicit ranges, no LLP (neither global nor within-cluster)
CGParams(
    L=128, membership=:implicit_ranges, inter_strategy=:perm,
    intra_ref_window=64, intra_ref_rle=false, intra_ref_fixwidth=true,
    intra_zigzag=true, intra_stop_deltas=true, intra_copy_blocks=true,
    intra_copy_adaptive=true,
    intra_add_adaptive=true, intra_raw_adaptive=true, intra_adapt_mil=2,
    varint=:fibonacci, degree=:elias_delta, perm_strategy=:blockpos,
    undirected_pairs=false,
)
```

**Recommended config for enwiki-2013** (12.4854 BPE, K=1, LLP, intervals + LR split):
```julia
# K=1 with global LLP reordering, interval+residual encoding, LR residual split
CGParams(
    L=128, membership=:implicit_ranges, inter_strategy=:perm,
    intra_ref_window=64, intra_ref_rle=false, intra_ref_fixwidth=true,
    intra_zigzag=true, intra_stop_deltas=true, intra_copy_blocks=true,
    intra_copy_adaptive=true,
    intra_intervals=true, intra_mil=4,
    intra_add_adaptive=true, intra_raw_adaptive=true, intra_adapt_mil=2,
    intra_lr_split=true,
    varint=:fibonacci, degree=:elias_delta, perm_strategy=:blockpos,
    undirected_pairs=false,
)
```

**CNR-2000 K=1 config without LR split** (2.4341 BPE, adaptive STOP-delta):
```julia
CGParams(
    L=128, membership=:implicit_ranges, inter_strategy=:perm,
    intra_ref_window=64, intra_ref_rle=false, intra_ref_fixwidth=true,
    intra_zigzag=true, intra_stop_deltas=true, intra_copy_blocks=true,
    intra_copy_adaptive=true,
    intra_add_adaptive=true, intra_raw_adaptive=true, intra_adapt_mil=2,
    varint=:fibonacci, degree=:elias_delta, perm_strategy=:blockpos,
    undirected_pairs=false,
)
```

---

## Bitstream Structure

Each coarsening level produces a bitstream with three sections:

**Children mode** (`coding_scheme=:children`, default):
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

**Index mode** (`coding_scheme=:index`): Adds offset tables for random access.

```
┌──────────────────────────────────────────────┐
│ CLUSTER OFFSET TABLE (2K+1 entries)          │
│   entry_width (6 bits)                       │
│   K (32 bits)                                │
│   entries[1..K]: intra-cluster bit offsets    │
│   entry[K+1]: inter-section start offset     │
│   entries[K+2..2K+1]: per-source-cluster     │
│     inter offsets (relative to inter start)   │
├──────────────────────────────────────────────┤
│ SECTION 1: Cluster Membership (unchanged)    │
├──────────────────────────────────────────────┤
│ SECTION 2: Intra-Cluster Edges               │
│   For each cluster c = 1..K:                 │
│     [Per-vertex ref header]                  │
│     [Per-vertex mil tags (if greedy_mil)]    │
│     [INTRA-VERTEX OFFSET TABLE]  ← NEW      │
│     [Per-vertex payloads]                    │
├──────────────────────────────────────────────┤
│ SECTION 3: Inter-Cluster Edges (unchanged)   │
└──────────────────────────────────────────────┘
```

The cluster offset table is written between the MGS header and the compressed data
(see MGS_HEADER.md for placement). The intra-vertex offset tables are embedded
within each cluster's data.

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

**Intra-vertex offset table** (index mode only, after ref header and mil tags):

When `coding_scheme=:index`, a per-vertex offset table is written within each cluster
to enable O(1) random access to any vertex's payload:

```
vtx_entry_width (6 bits): number of bits per entry
s+1 entries (each vtx_entry_width bits):
  entry[0..s-1]: bit offset to local vertex i's payload (relative to payload start)
  entry[s]:      total payload bit size
```

The two-pass encoding first writes all per-vertex payloads to a temporary buffer,
recording bit offsets, then writes the offset table followed by bit-precise payload
data to the main stream.

**Per-vertex payload** (in cluster vertex order):

For **referenced vertices** (`has_ref=true`):
```
positions:
  IF intra_copy_adaptive AND intra_copy_blocks:
    bit: outer mode (1 = bitmap, 0 = copy-blocks or complement)
    IF outer = 1 (bitmap):
      ref_len bits: one bit per position in reference's neighbor list (1 = copied)
    ELSE:
      bit: inner mode (0 = copy-blocks, 1 = complement)
      IF inner = 0 (copy-blocks):
        copy-blocks encoding of copied positions (see Copy-Blocks below)
      ELSE (complement):
        copy-blocks encoding of SKIPPED positions (positions NOT copied)
        decoder: include all ref positions not in the decoded skip set
  ELIF intra_copy_blocks:
    copy-blocks encoding (see Copy-Blocks below)
  ELSE:
    RLE ones-deltas bitmap over the reference vertex's neighbor list

additions:
  IF (intra_intervals OR intra_greedy_mil) AND intra_lr_split:
    LR-split intervals+residuals encoding [‡]
  ELIF intra_intervals OR intra_greedy_mil:
    intervals_and_residuals encoding [†]
  ELIF intra_add_adaptive AND intra_stop_deltas:
    bit: mode flag (1 = intervals, 0 = STOP-delta)
    IF mode = intervals:
      intervals_and_residuals encoding [†]
    ELSE:
      STOP-terminated delta list [†]
  ELIF intra_stop_deltas:
    STOP-terminated delta list [†]
  ELSE (additions_mode = :delta):
    small_count(num_additions)
    delta-encoded addition list [†]
```

For **non-referenced vertices** (`has_ref=false`):
```
IF (intra_intervals OR intra_greedy_mil) AND intra_lr_split:
  LR-split intervals+residuals encoding [‡]
ELIF intra_intervals OR intra_greedy_mil:
  intervals_and_residuals encoding [†]
ELIF intra_raw_adaptive AND intra_stop_deltas:
  bit: mode flag (1 = intervals, 0 = STOP-delta)
  IF mode = intervals:
    intervals_and_residuals encoding [†]
  ELSE:
    STOP-terminated delta list [†]
ELIF intra_stop_deltas:
  STOP-terminated delta list [†]
ELSE:
  small_count(degree_in_cluster)
  delta-encoded sorted local neighbor IDs [†]
```

`[†]` When `intra_zigzag=true`, the first value is encoded as a zigzag offset from the current vertex's local index (see below). Does not apply to reference position lists.

`[‡]` Left/right residual split encoding (see below). Requires `intra_zigzag=true`.

#### Reference Encoding Decision

For each vertex, the encoder measures:
1. **Raw cost**: delta-coded full neighbor list
2. **Reference cost** (for each candidate in lookback window): copy-blocks positions + additions

The cheapest option wins. Two-pointer merge computes positions (shared elements) and additions (new elements). The lookback window contains at most `intra_ref_window` previous vertices in cluster order.

### Section 3: Inter-Cluster Edges

**Index mode per-source-cluster offsets**: In index mode, the extended cluster offset
table contains `2K+1` entries. Entries `K+2..2K+1` store per-source-cluster bit offsets
within the inter-cluster section (relative to inter section start). This enables O(1)
access to any source cluster's inter-cluster edge data.

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

**Complement variant** (when `intra_copy_adaptive=true`): For high-overlap references, it is cheaper to encode the SKIPPED positions using the same copy-blocks format. Example: reference has 20 neighbors; current vertex shares {1..18} (skips {19,20}).
- Copy-blocks would encode 18 copied positions: many runs, high cost
- Complement encodes 2 skipped positions: `small_count(1), fib(19), fib(2)` = ~12 bits
- Bitmap would use 20 bits

The encoder tries all three modes, picks the cheapest, and writes the appropriate mode bits.

Advantages over RLE ones-deltas bitmap:
1. Skip runs implicit (encoded as gaps) — no per-run type flag
2. No explicit bitmap length (decoder knows ref_len from reference list)
3. All values ≥ 1 — no +1 shifts needed
4. Saves ~0.398 BPE on CNR-2000 (vs old RLE bitmap); +3-way adaptive saves additional −0.116 BPE

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

### Left/Right Residual Split

When `intra_lr_split=true` (requires `intra_intervals=true` and `intra_zigzag=true`), neighbor lists are encoded using `_write_ir_lr` which applies a WebGraph-inspired left/right split to the residuals after interval extraction. This produces smaller first-value encodings by ensuring each half's residuals start near vertex_id.

**Encoding format** (`_write_ir_lr`):
```
1. Extract intervals from the full sorted neighbor list (standard algorithm)
2. Write intervals section (identical to standard encoding):
   varint(num_intervals + 1)
   For each interval (start, length):
     First interval: varint(zigzag(start - vertex_id) + 1)
     Subsequent: varint(start - prev_start)
     varint(length - mil + 1)

3. Write residuals with left/right split:
   varint(num_residuals + 1)
   IF num_residuals > 0:
     Split residuals at vertex_id:
       left  = residuals where value < vertex_id
       right = residuals where value > vertex_id
     varint(num_left + 1)

     Left distances (reversed to ascending):
       left_dists[i] = vertex_id - left[num_left - i + 1]    (i = 1..num_left)
       delta-encoded: first value + gaps (all values ≥ 1)

     Right distances (ascending):
       right_dists[i] = right[i] - vertex_id                 (i = 1..num_right)
       delta-encoded: first value + gaps (all values ≥ 1)
```

**Decoding** (`_read_ir_lr`):
```
1. Read intervals (identical to standard)
2. Read num_residuals; if > 0:
   Read num_left; num_right = num_residuals - num_left
   Read left_dists (delta-decode num_left values)
   Read right_dists (delta-decode num_right values)
   Reconstruct: left[i] = vertex_id - left_dists[num_left - i + 1]
                right[i] = vertex_id + right_dists[i]
3. Merge intervals + left + right, sort ascending
```

**Why LR split helps**: Standard zigzag encoding of residuals writes the first value as `zigzag(first_residual - vertex_id) + 1`. If the first residual is far from vertex_id (e.g., vertex_id=100, first_residual=50), the zigzag value is large: `zigzag(-50) = 99 → ~7 bits`. With LR split, the closest left residual (e.g., 95) has distance 5 (`~3 bits`), and the closest right residual (e.g., 102) has distance 2 (`~2 bits`). Additionally, the gap that crosses vertex_id in the standard encoding is eliminated entirely.

**Overhead**: One extra Fibonacci-encoded count (left_count) per list, typically 2–4 bits. For empty halves (all neighbors on one side of vertex_id), the overhead is just the left_count value.

**Results**:
- enwiki-2013: **−0.442 BPE** (12.927 → 12.485 BPE), crossing below WebGraph BV-HC (12.639)
- cnr-2000 (K=1, intervals mode): −0.113 BPE (2.684 → 2.570 BPE)

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

`CGStats` accumulates per-section bit counts:

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
| `intra_copy_bitmap_count` | Number of ref vertices that chose bitmap mode (when `intra_copy_adaptive`) |
| `intra_copy_blocks_count` | Number of ref vertices that chose copy-blocks mode (when `intra_copy_adaptive`) |
| `intra_copy_complement_count` | Number of ref vertices that chose complement mode (when `intra_copy_adaptive`) |
| `intra_overlap_hist` | 10-bucket histogram of overlap fractions [0%,10%)…[90%,100%) for ref vertices |
| `intra_add_count_total` | Total number of additions across all ref vertices |

---

## Performance

### CNR-2000

**Graph**: 325,557 vertices, 3,216,152 edges (Italian web crawl).
**WebGraph reference**: BV = 2.898 BPE (1,164,843 bytes), BV-HC = 2.448 BPE (984,261 bytes).

#### Best Result: 2.3286 BPE (K=2, no LLP)

The best CNR-2000 result uses K=2 Leiden clusters without global LLP pre-ordering, window=64, and adaptive STOP-delta encoding. The K=2 split produces two clusters (18,292 and 307,265 vertices) where the per-cluster LLP naturally orders vertices within each Leiden group. No global LLP reordering is applied (LLP_PASSES=0).

| Component | Bits | BPE | Share |
|-----------|------|-----|-------|
| Membership (implicit ranges) | 72 | 0.000 | 0.0% |
| Headers (fixed-width ref info) | 488,093 | 0.152 | 6.5% |
| Copy positions (3-way adaptive) | 1,606,737 | 0.500 | 21.4% |
| Additions (adaptive: STOP-delta or intervals) | 3,243,903 | 1.009 | 43.3% |
| Raw (adaptive: STOP-delta or intervals) | 1,098,678 | 0.342 | 14.7% |
| Inter-cluster edges | 1,918 | 0.001 | 0.0% |
| **Total** | **7,489,040** | **2.3286** | 100% |

Ref used: 202,015 / 325,557 (62.1%). Copy modes: bitmap=102,899 blocks=43,267 complement=55,849.

#### Index Mode BPE (CNR-2000)

Index mode adds per-vertex offset tables for random access at the cost of additional
storage. The overhead is consistent across all algorithms (~2.3 BPE):

| Algorithm | Children BPE | Index BPE | Delta |
|-----------|-------------|-----------|-------|
| BG | 2.5352 | 4.8297 | +2.2945 |
| CS | 2.5347 | 4.8351 | +2.3004 |
| CG (K=2) | 2.3286 | 4.6246 | +2.3055 |

The ~2.3 BPE overhead comes from N+1 offset entries (~25 bits each) for 325K vertices
over 3.2M edges. For CG, additional overhead includes the extended 2K+1 cluster offset
table and per-cluster intra-vertex offset tables.

### Multi-Dataset Benchmark

Best CG BPE across all tested datasets, compared with BG, CS, and WebGraph BV:

| Dataset | Vertices | Edges | CG BPE | BG BPE | CS BPE | BV BPE | CG Config |
|---------|----------|-------|---------|---------|--------|--------|------------|
| cnr-2000 (no reorder) | 325K | 3.2M | **2.3286** | 2.4929 | 2.4348 | 2.898 | K=2, w=64, no LLP |
| cnr-2000 (Leiden+LLP) | 325K | 3.2M | 2.5652 | 2.3259 | **2.3043** | 3.2335 | K=1, w=64; CS w=256 best |
| in-2004 | 1.38M | 16.9M | **1.7513** | 1.895 | 1.7839 | 2.172 | K=1, w=8 |
| enwiki-2013 (LLP) | 4.2M | 101M | 12.4854 | — | — | 13.114 | K=1, w=64, LLP, iv+LR |
| enwiki-2013 (Leiden+LLP) | 4.2M | 101M | — | **12.161** | 12.222 | 13.114 | BG w=64 lr+mr; CS w=128 lr |
| enwiki-2013 (no reorder) | 4.2M | 101M | 15.7183 | — | — | 13.114 | K=1, w=64, iv+LR |
| web-google core | 434K | 3.4M | 4.3296 | 4.0735 | **4.0288** | 5.0081 | K=1, w=8, Leiden+LLP, no-iv m4 |
| web-google rcore | 434K | 3.4M | 3.9359 | 3.7626 | **3.7337** | 4.1751 | K=1, w=8, Leiden+LLP, no-iv m4 |
| eat-core | 7.8K | 247K | **10.5768** | 10.9191 | 10.8679 | 10.705 | K=1, w=256, iv+LR+td, gamma |
| eat-core (Leiden+LLP) | 7.8K | 247K | **9.558** | 9.767 | 9.773 | 9.729 (HC) | K=1, w=256, iv+LR |
| eat-rcore | 7.8K | 247K | **9.3192** | 9.5726 | 9.5674 | 9.391 | K=1, w=256, iv+LR+td, gamma |
| arxiv-hep-ph core | 34.5K | 422K | **9.8187** | 10.0910 | 10.0767 | 10.262 | K=1, w=256, iv+LR+td |
| arxiv-hep-ph core (Leiden+LLP) | 34.5K | 422K | **7.162** | 7.252 | 7.269 | 7.706 (HC) | K=1, w=256, iv+LR |
| arxiv-hep-ph rcore | 12.7K | 140K | **8.9406** | 9.0946 | 8.9855 | 9.684 | K=1, w=256, iv+LR+td |
| amazon-0601 core | 395K | 3.3M | **12.1853** | 12.4444 | 12.5235 | 13.001 | K=1, w=8, iv+LR+td |
| amazon-0601 core (Leiden+LLP) | 395K | 3.3M | **7.903** | 8.058 | 8.095 | 8.722 (HC) | K=1, w=8, iv+LR |
| amazon-0601 rcore | 395K | 3.3M | **11.7069** | 11.7310 | 11.8497 | 12.064 | K=1, w=8, iv+LR+td |

**Bold** = best among Adjacently methods for that dataset.

**Key findings**:
- CG K=1 beats BG and CS on most datasets thanks to per-cluster local vertex indexing and adaptive encoding
- **Web-Google exception**: CS and BG still beat CG K=1 even with Leiden+LLP ordering. CG's fixwidth ref header overhead (1+ceil(log2(w)) bits per vertex) outweighs its advantages. All methods benefit equally (~0.3 BPE) from Leiden+LLP ordering.
- **enwiki-2013 exception**: BG (12.161, Leiden+LLP) beats CG (12.485, LLP). Multi-ref is valuable on the high-degree Wikipedia graph (avg degree 24.1). CG without reordering degrades to 15.718 BPE.
- **Leiden+LLP ordering** improves all methods: LLP→Leiden fine clusters→per-cluster LLP→concatenate. Improves CG by 0.35 BPE, CS by 0.37 BPE, BG by 0.35 BPE, BV by 0.34 BPE over plain LLP.
- CG w=8 beats w=64 on Web-Google: the +3 bits/ref-vertex for w=64 fixwidth outweighs payload savings
- `no-iv m4` mode (no intervals, mil=4) beats intervals+LR at high locality (Leiden+LLP produces very tight clusters)
- `tight_deltas` saves 0.06-0.16 BPE on datasets using intervals+LR-split

#### Optimization History (K=1 progression)

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
| WebGraph BV (reference) | — | 2.898 | — | 1,138 KB |
| + Implicit ranges membership | 1 | 2.666 | −0.221 | 1,047 KB |
| + Adaptive additions (MIL=2) | 1 | 2.599 | −0.067 | 1,022 KB |
| + Adaptive raw (MIL=2) | 1 | 2.617 | −0.049 | 1,028 KB |
| + Adaptive adds+raw (MIL=2) | 1 | 2.550 | −0.116 | 1,001 KB |
| + 2-way adaptive copy (bitmap vs copy-blocks) | 1 | 2.497 | −0.053 | 981 KB |
| + 3-way adaptive copy (+ complement) | 1 | 2.4341 | −0.063 | 956 KB |
| **K=2, no LLP** | **2** | **2.3286** | **−0.115** | **932 KB** |

**Total improvement**: 4.307 → 2.3286 BPE (−46.2% over 12 independent optimizations). The final K=2 step uses Leiden K=2 partitioning without global LLP, which produces better per-cluster locality than K=1 with the same LLP ordering.

Note: each adaptive optimization targets a different section (adds, raw, copy positions), and their savings are additive. Window size w=64 is optimal: w=128 adds +1 bit/ref-delta in headers (+26 KB) but saves only 23 KB in adds+raw — net negative. All roundtrips verified OK.

### Best Config Bit Budget (FW64, K=1, 3-Way Adaptive — 2.4341 BPE)

| Component | Bytes | BPE | Share |
|-----------|-------|-----|-------|
| Membership (2 sizes, implicit) | 7 | 0.000 | 0.0% |
| Headers (fixed-width ref info) | 193,406 | 0.481 | 20.3% |
| Copy positions (3-way adaptive: bitmap/copy-blocks/complement) | 229,790 | 0.572 | 24.1% |
| Additions (adaptive: STOP-delta or intervals) | 422,987 | 1.052 | 44.3% |
| Raw (adaptive: STOP-delta or intervals) | 132,121 | 0.329 | 13.9% |
| Inter-cluster edges | 241 | 0.001 | 0.0% |
| **Total** | **978,552** | **2.4341** | 100% |

62.5% of vertices (203,615 / 325,557) use reference encoding. Copy mode distribution: 61.3% bitmap (124,864), 22.6% copy-blocks (45,982), 16.1% complement (32,769). Overlap histogram peak: 49% of ref vertices have ≥90% of reference's neighbors copied — the dominant case for complement.

Adaptive savings vs fixed-encoding baseline (2.666 BPE): copy −46,643 bytes, additions −27,015 bytes, raw −19,596 bytes.

### Previous Best Config Bit Budget (FW64, K=1, 2-Way Adaptive — 2.497 BPE)

| Component | Bytes | BPE | Share |
|-----------|-------|-----|-------|
| Membership (2 sizes, implicit) | 7 | 0.000 | 0.0% |
| Headers (fixed-width ref info) | 193,406 | 0.481 | 19.3% |
| Copy positions (adaptive: copy-blocks or bitmap) | 255,092 | 0.635 | 25.4% |
| Additions (adaptive: STOP-delta or intervals) | 422,987 | 1.052 | 42.2% |
| Raw (adaptive: STOP-delta or intervals) | 132,121 | 0.329 | 13.2% |
| Inter-cluster edges | 241 | 0.001 | 0.0% |
| **Total** | **1,003,854** | **2.497** | 100% |

62.5% of vertices (203,615 / 325,557) use reference encoding.
Adaptive savings vs fixed-encoding baseline (2.666 BPE): copy −21,341 bytes, additions −27,015 bytes, raw −19,596 bytes.

### Previous Config Bit Budget (FW64, K=1, Adaptive Adds+Raw — 2.550 BPE)

| Component | Bytes | BPE | Share |
|-----------|-------|-----|-------|
| Membership (2 sizes, implicit) | 7 | 0.000 | 0.0% |
| Headers (fixed-width ref info) | 193,406 | 0.481 | 18.9% |
| Copy positions (copy-blocks, fixed) | 276,433 | 0.688 | 27.0% |
| Additions (adaptive: STOP-delta or intervals) | 422,987 | 1.052 | 41.3% |
| Raw (adaptive: STOP-delta or intervals) | 132,121 | 0.329 | 12.9% |
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

### Comparison with WebGraph BV (CNR-2000)

| Component | CG K=2 (best) | CG K=1 + 3-Way Adaptive | CG K=1 + EF | WebGraph BV |
|-----------|-----------------|---------------------------|---------------|-------------|
| Membership / structure | 0.000 BPE | 0.000 BPE | 0.221 BPE | ~0.078 BPE |
| Headers | 0.478 BPE | 0.481 BPE | 0.481 BPE | ~0.126 BPE |
| Copy positions | 0.499 BPE | 0.572 BPE | 0.688 BPE | ~0.556 BPE |
| Residuals/additions | 1.003 BPE | 1.052 BPE | 1.119 BPE | ~1.583 BPE |
| Raw / non-referenced | 0.338 BPE | 0.329 BPE | 0.377 BPE | ~0.554 BPE |
| Inter-cluster | 0.001 BPE | 0.001 BPE | 0.001 BPE | — |
| **Total** | **2.3286 BPE** | **2.4341 BPE** | **2.887 BPE** | **2.898 BPE** |

CG K=2 beats WebGraph BV by **0.579 BPE** (20.0% reduction) and BV-HC by **0.129 BPE** (5.3% reduction). The K=2 split uses **no within-cluster LLP** — the crawl-order locality within each Leiden community is already excellent for CNR-2000, and within-cluster LLP actually hurts by +0.290 BPE (see finding #17).

### Greedy Approach Comparison (CNR-2000)

| Approach | Ordering | Window | BPE |
|----------|----------|--------|-----|
| Greedy (baseline) | None (crawl order) | 7 | ~4.3 |
| Greedy | Global LLP | 7 | 3.561 |
| Greedy | Global LLP | 64 | 3.319 |
| **Greedy (best)** | **Leiden + per-group LLP** | **64** | **2.881** |
| WebGraph BV | Greedy bisection | 7 | 2.897 |
| CG (K=1, EF) | Leiden + per-cluster LLP | 64 | 2.887 |
| CG (K=1, implicit ranges) | Leiden + per-cluster LLP | 64 | 2.666 |
| CG (K=1, 3-way adaptive copy) | Leiden + per-cluster LLP | 64 | 2.4341 |
| **CG (K=2, no LLP)** | **Leiden (crawl order within clusters)** | **64** | **2.3286** |

### enwiki-2013

**Graph**: 4,203,325 vertices (4,206,785 with isolates), 101,355,853 edges (English Wikipedia 2013 link graph).
**WebGraph reference**: BV standard = 13.114 BPE, BV-HC (high compression) = 12.639 BPE.

#### Best Result: 12.161 BPE (BG, w=64, lr+mr, Leiden+LLP)

BG with Leiden+LLP ordering beats CG on enwiki-2013 — a reversal from most other datasets. BG's multi-reference feature is valuable on the high-degree Wikipedia graph (avg degree 24.1).

| Configuration | BPE | Bytes | vs BV-HC |
|---------------|-----|-------|----------|
| CG K=1, iv+LR, original ordering | 15.718 | 199,143,062 | +3.079 |
| WebGraph BV standard (w=7, zeta-3) | 13.114 | 166,148,696 | +0.475 |
| WebGraph BV-HC (high compression) | 12.639 | — | — |
| CG K=1, intervals + LR split (LLP) | 12.485 | 158,184,071 | −0.154 |
| CS w=128, lr (Leiden+LLP) | 12.222 | — | −0.417 |
| **BG w=64, lr+mr (Leiden+LLP)** | **12.161** | — | **−0.478** |

LLP ordering saves **3.233 BPE** (20.6%) over original ordering for CG K=1 with intervals+LR-split. Without LLP, CG (15.718) is worse than even BV standard (13.114), showing that enwiki-2013's original ordering has very weak sequential locality.

With Leiden+LLP, BG (12.161) beats CG (12.485) by 0.324 BPE. CS w=128 (12.222) is also ahead. CS w=256 is estimated at ~12.15 BPE.

The LR split saves **0.442 BPE** on enwiki-2013 (3.4% reduction), pushing CG below the WebGraph BV-HC baseline. This is a larger saving than on CNR-2000 because enwiki-2013 has higher average degree (24.1 vs 9.9) and more residuals per vertex after interval extraction, amplifying the benefit of better first-value encoding.

Encoding time: ~5 min (CG original ordering), ~38 min (CG LLP, Julia `--threads=auto`).
Ref used (CG original): 3,526,081 / 4,203,325 (83.9%). Ref used (CG LLP): 3,560,563 / 4,203,325 (84.7%).

---

## Synthetic Web Graph Benchmark (N=10000, original ordering)

CG parameter sweep against BV (w=64) on `random_web_digraph` graphs.
Four strategies compared: hand-tuned baseline, GNN-predicted params, auto_select_K, and exhaustive grid search (512 combos).

| avg_deg | BV (w=64) | CG baseline | CG GNN | CG K=auto | CG grid | Grid config |
|---------|-----------|-------------|--------|-----------|---------|-------------|
| 12 | 9.944 | 10.130 | 10.371 | 10.130 (K=1) | **9.318** (-0.627) | w=8, iv+lr, mil=3, zz, no-sd |
| 24 | 9.417 | 9.515 | 9.084 | 9.515 (K=1) | **9.074** (-0.343) | w=8, iv+lr, mil=4, zz, no-sd |
| 32 | 9.375 | 9.469 | 9.143 | 9.469 (K=1) | **9.132** (-0.244) | w=8, iv+lr, mil=5, zz, no-sd |
| 64 | 9.445 | 9.579 | 9.405 | 9.579 (K=1) | **9.390** (-0.055) | w=64, iv+lr, mil=5, zz, no-sd |

**Key findings**:
- **CG grid search beats BV by 0.05-0.63 BPE** across all densities — the largest gains of any method
- **intervals=true + lr_split=true is the critical combination** — the hand-tuned baseline (no intervals, no lr_split) loses to BV, but grid search with iv+lr wins decisively
- **GNN prediction is near-optimal at deg 24-64** (within 0.01-0.02 BPE of grid), correctly predicting iv+lr+w=8; overshoots at deg=12 (intervals=false)
- **auto_select_K always returns K=1** — synthetic web graphs lack community structure (coarse granularity <0.01), so Leiden clustering cannot help
- **Small window (w=8) is optimal for deg 12-32** — tight sequential locality in web graphs; w=64 wins only at deg=64 where locality is diluted
- **stop_deltas=false is consistently best** — the STOP-delta overhead isn't amortized when intervals+lr_split handles the encoding more efficiently
- CG grid search (-0.63 BPE at deg=12) outperforms CS zeta+lr (-0.38 BPE) — CG's fixwidth ref + adaptive copy is more effective than CS's prefix codes when properly tuned
- **Caveat**: tuned BV (zeta-5, i=2, m=-1) beats all Adjacently methods at deg=64 (9.201 vs CG 9.390) and nearly ties CG at deg=32 (9.152 vs 9.132). CG's advantage is largest at low degree where BV's fixed pipeline cannot adapt.

## LFR Benchmark (N=10000, avg_degree=15, tau1=2.5, tau2=1.5)

CG compression on LFR graphs with planted community structure, sweeping mixing parameter μ.

### Without reordering — CG K-sweep

| μ | BV | CG K=1 | CG K=2 | CG K=4 | CG K=8 | CG K=16 | Best |
|------|------|--------|--------|--------|--------|---------|------|
| 0.05 | 13.89 | 14.03 | 13.49 | 12.54 | 10.86 | **9.67** | K=16 (-4.21) |
| 0.10 | 13.81 | 13.97 | 13.84 | 13.39 | 12.10 | **10.99** | K=16 (-2.82) |
| 0.20 | 13.66 | 13.84 | 14.11 | 14.23 | 13.70 | **12.77** | K=16 (-0.89) |
| 0.30 | 13.53 | **13.73** | 14.26 | 14.87 | 14.77 | 14.09 | K=1 (+0.20) |
| 0.50 | **13.28** | 13.50 | 14.36 | 15.55 | 16.26 | 16.31 | BV |

### With Leiden+LLP — all methods comparison

| μ | BV | BG best | CS best | CG K=1 grid | CG config |
|------|------|---------|---------|-------------|-----------|
| 0.05 | 8.44 | 8.36 | 8.30 | **8.08** | w=32, iv+lr, mil=4, zz |
| 0.10 | 9.13 | 9.06 | 9.00 | **8.77** | w=32, iv+lr, mil=3, zz |
| 0.20 | 9.95 | 9.91 | 9.83 | **9.63** | w=32, iv+lr, mil=3, zz |
| 0.30 | 10.61 | 10.58 | 10.51 | **10.35** | w=64, iv+lr, mil=4, zz, sd |
| 0.50 | 11.60 | 11.56 | 11.49 | **11.42** | w=64, iv+lr, mil=4, zz |

### With Leiden+LLP — CG K-sweep (K=1 always best)

| μ | CG K=1 | CG K=2 | CG K=4 | CG K=8 | CG K=16 |
|------|--------|--------|--------|--------|---------|
| 0.05 | **7.94** | 8.07 | 8.25 | 8.36 | 8.33 |

**Key LFR findings**:
- **Leiden+LLP is the dominant factor**: 3.5-5.4 BPE gain at low μ, dwarfing all method differences
- **Without reordering, CG K=16 at μ=0.05 saves 4.2 BPE over BV** — per-cluster encoding exploits community structure directly
- **With Leiden+LLP, K=1 is always best** — reordering handles communities globally, making multi-cluster encoding counterproductive
- **CG K=1 grid beats all methods with Leiden+LLP** by 0.07-0.37 BPE at every μ
- **auto_select_K always returns K=1** (10K < 100K threshold) — misses the K=16 opportunity on unordered graphs

## Key Findings

1. **Vertex ordering is the dominant factor** (+1.43 BPE advantage over unordered greedy). The two-step approach — Leiden community detection creates ~34K tight fine clusters of ~9 vertices each, then per-cluster LLP ordering makes consecutive vertices share 80–90% of neighbors — is what enables CG to beat WebGraph.

2. **K=2 is the optimal cluster count for CNR-2000** (2.3286 BPE vs 2.4341 for K=1). The two-level Leiden split produces two clusters where the original crawl order already provides excellent locality — no within-cluster LLP is needed (LLP actually hurts, see finding #17). Membership cost is 0.000 BPE with implicit ranges (two contiguous ID ranges). For datasets where K=2 is infeasible or unhelpful (e.g., enwiki-2013 at 4.2M vertices), K=1 with global LLP remains effective. The general principle: reducing from K=8 to K=1 saves −0.104 BPE (membership 0.318 → 0.000 BPE with implicit ranges, inter-cluster negligible at 241 bytes).

8b. **Implicit ranges membership** (−0.221 BPE). When vertices are pre-relabeled so the two K=1 groups occupy contiguous ID ranges [1..S₁] and [S₁+1..N], membership can be encoded as just two size varints (7 bytes) instead of Elias-Fano sorted lists (88,764 bytes). The key implementation insight: `encode_level` always re-sorts cluster arrays by vertex ID internally, so the pre-relabeling must use **vertex-ID rank** within each group (not LLP rank) to preserve bit-for-bit identical intra encoding. With this approach, the intra section is completely unchanged and only the membership section is eliminated.

3. **Copy-blocks is the largest codec optimization** (−0.398 BPE). WebGraph-style copy-blocks replace the RLE bitmap: skip runs are implicit (gaps between copy blocks), no per-run type flag, no bitmap length prefix. Additionally, the more accurate cost model causes 5K extra vertices to switch to reference encoding.

4. **Zigzag first-value encoding** (−0.741 BPE) exploits the per-cluster LLP locality. First neighbors are typically close to the vertex's own local index; the signed offset is small and compresses well with Fibonacci.

5. **Fixed-width ref deltas** (−0.094 to −0.126 BPE) replace the byte-padded bitmap + Fibonacci varint approach. With window=64, each ref vertex costs exactly 7 bits (1 flag + 6 data) vs ~7.5 bits average with Fibonacci + bitmap padding.

6. **Window=64 outperforms window=32** (−0.032 BPE). More reference candidates → 1,519 additional ref vertices, fewer additions, at the cost of +1 bit per ref delta (6 vs 5 fixed bits). Net positive.

7. **Zeta-3 encoding is harmful for CG** (+0.061 to +0.099 BPE vs Fibonacci). Despite being optimal for WebGraph's residuals, zeta-3 degrades copy-blocks positions encoding — it encodes copy-block starts/gaps/lengths less efficiently than Fibonacci for the small values typical in dense clusters.

8. **Global LLP + fixed interval encoding are both harmful**: Global LLP pre-ordering (−0.33 BPE from Leiden quality), fixed interval encoding (+0.119 BPE overhead not amortized on short intra-cluster lists).

9. **Adaptive per-vertex codec selection** (−0.116 BPE combined). With a single inline mode bit per vertex, the encoder picks the cheaper of STOP-terminated delta-list vs interval+residuals encoding. This recovers the benefit of interval encoding for vertices with long arithmetic progressions in their neighbor lists, while avoiding overhead for vertices where delta encoding is cheaper. At MIL=2 (aggressive interval detection): additions −0.067 BPE (450 KB → 423 KB), raw −0.049 BPE (152 KB → 132 KB), total −0.116 BPE.

10. **VLC reference-delta encoding is harmful** (+0.093 BPE). Using Fibonacci variable-length coding for ref deltas instead of fixed 6-bit encoding was expected to help but hurts in practice — ref deltas are distributed uniformly over [1,64] rather than concentrated near 1, so fixed-width (6 bits average) beats Fibonacci (~7–8 bits average for uniform input). This demonstrates that VLC only helps when the distribution is heavily skewed toward small values.

11. **2-way adaptive copy-position encoding** (−0.053 BPE). Per-ref-vertex choice between copy-blocks and a dense bitmap over the reference's neighbor list (one bit per reference neighbor, 1=copied). Bitmap costs exactly `ref_len` bits, which beats copy-blocks when the reference list is short (≤ ~8 neighbors) or when copied positions are scattered. Mode flag overhead: 1 bit × 203,615 ref vertices = 25 KB; gross copy savings: 47 KB; net: −21 KB = −0.053 BPE.

12. **Window size w=64 is optimal** (larger windows are harmful). Increasing to w=128 adds 1 extra bit per ref-delta (+26 KB headers) but saves only 23 KB in additions+raw (0.4% more ref vertices found). Net: −3 KB (slightly worse). w=256 even more so: +53 KB headers, −44 KB adds+raw, net −9 KB. The additional references found beyond the 64-vertex window are too few to justify the wider fixed-width delta encoding.

13. **Complement copy-blocks (3-way adaptive)** (−0.063 BPE, copy 255 KB → 230 KB). A third copy mode encodes the SKIPPED positions instead of the copied ones, using the same copy-blocks format. For high-overlap references (≥90% copied), complement is dramatically cheaper: encoding 1–2 skipped positions costs ~4 bits vs encoding 18–20 copied positions in a run. The 3-way mode uses nested bits: outer=1→bitmap (1 bit overhead), outer=0+inner=0→copy-blocks (2 bits), outer=0+inner=1→complement (2 bits). Despite the 1-bit additional overhead for non-bitmap modes, complement dominates for the 16.1% of ref vertices (32,769 total) where it's cheapest. The overlap histogram shows why: 49% of ref vertices (99,760) have ≥90% overlap fraction — the main beneficiaries. Net savings: 255 KB − 230 KB = 25 KB beyond the 2-way baseline.

14. **Smart ref-delta encoding (delta=1 shorthand) is harmful** (+0.040 BPE). The intuition was that delta=1 (reference immediately preceding vertex) would dominate and could be encoded as just 2 bits instead of 7. In practice, **only 6% of ref vertices use delta=1**. The actual distribution is near-uniform over [1,64]: delta=49–64 is the largest bucket at 35.8%, and small deltas (1–4) account for only 17% combined. This confirms finding #10: the LLP-ordered large top-level group has useful references spread throughout the entire 64-vertex window, not concentrated at delta=1. Encoding delta=1 in 2 bits saves 5 bits for 6% of refs but costs 1 extra bit for the remaining 94%, netting +0.64 bits/ref × 203K refs = +16 KB ≈ +0.040 BPE.

15. **Elias-Fano for copy positions is harmful** (+0.034 BPE). EF was expected to help for large `ref_len` with scattered copied positions (e.g., k=10 scattered copies in ref_len=64 costs ~51 bits in EF vs 64 bits bitmap vs ~95 bits copy-blocks). In practice, EF is chosen for only **0.4% of ref vertices** (753 out of 203,615) — far too few to compensate for the mode-bit overhead. The flat 2-bit mode scheme (to accommodate a 4th mode) adds 1 extra bit for all 124,864 bitmap vertices (+15.6 KB), while 753 EF-mode vertices save only ~2 KB. The root cause: most ref vertices fall into the three cases already handled well by the existing 3-way adaptive — small ref_len (bitmap wins), high overlap (complement wins), or run-structured positions (copy-blocks wins). The large-ref_len + scattered regime where EF would excel is rare in this dataset.

16. **Left/right residual split** (−0.113 BPE on CNR-2000 K=1, −0.442 BPE on enwiki-2013 K=1). After interval extraction, remaining residuals are split at the vertex's own ID into left (< vertex_id) and right (> vertex_id) halves. Each half is transformed to ascending distances from vertex_id (left: `vertex_id − val`, reversed; right: `val − vertex_id`) and delta-encoded independently. This replaces zigzag encoding of the first residual value: instead of one potentially large zigzag-encoded signed offset, LR split produces two streams where the first value in each is a small positive distance (typically 1–10). Overhead is one Fibonacci-encoded `left_count` per list (~2–4 bits). The benefit scales with average degree — enwiki-2013 (avg degree 24.1) has more residuals per vertex than CNR-2000 (avg degree 9.9), so the saving is proportionally larger. On enwiki-2013, LR split pushes CG to **12.485 BPE**, surpassing the WebGraph BV-HC baseline (12.639 BPE) by 0.154 BPE.

17. **Within-cluster LLP is harmful for CNR-2000 K=2** (+0.290 BPE). With the Leiden K=2 split, applying LLP within each cluster degrades compression from 2.3286 to 2.6087 BPE. The crawl order within each Leiden community already provides strong locality for CNR-2000 — consecutive vertices in crawl order share many neighbors because the web crawler traverses links locally. LLP disrupts this natural ordering by rearranging vertices based on label propagation, which optimizes for a different notion of locality. This is dataset-specific: for enwiki-2013 (K=1), global LLP is essential because Wikipedia articles are not crawl-ordered. The key insight is that the optimal vertex ordering depends on both the dataset's inherent structure and the clustering strategy.

18. **Intervals + LR split are harmful for CNR-2000 K=2** (+0.167 BPE). Even with the LR split improvement, enabling interval+residual encoding on the best K=2 config (no within-cluster LLP) increases BPE from 2.3286 to 2.4865. The full comparison on the same Leiden partition: K=2 baseline 2.3286, K=2 + intervals 2.4865 (+0.167), K=2 + intervals + LLP 2.8093 (+0.490), K=2 + intervals + LR + LLP 2.7460 (+0.427). The STOP-delta encoding with zigzag first values remains superior for CNR-2000's short intra-cluster neighbor lists (avg degree 9.9), where interval detection overhead is not amortized. LR split helps on datasets with higher average degree (enwiki-2013, avg 24.1) where more residuals exist after interval extraction.

### Encoding Speed Optimizations

CG K=2 encoding on CNR-2000 was optimized from 1.558e-05 to **2.164e-06 sec/edge** (7.2× faster), making CG the fastest Adjacently encoder — 2.9× faster than BG (6.30e-06) and 1.8× faster than CS (3.97e-06). Five optimizations were applied:

| Optimization | Effect |
|-------------|--------|
| **Analytical cost estimation** | Replaced IOBuffer trial encoding with pure arithmetic (`estimate_encoded_value_cost`, `_fibonacci_bit_length`). Eliminates ~6 IOBuffer allocations per vertex. |
| **Two-phase candidate pruning** | Phase 1: cheap `_sorted_overlap_count` screens all 64 window candidates. Phase 2: full analytical evaluation on top 16 (`MAX_REF_CANDIDATES_PHASE2`). |
| **Early termination** | When positions cost alone exceeds current best, skip additions estimation. |
| **Pre-built neighbor lists** | All cluster neighbor lists built upfront before per-vertex reference search. |
| **Double-buffer swap** | Two position/adds vector pairs swapped on improvement, eliminating both `copy()` and re-merge overhead. |

BPE impact of 2-phase pruning (CNR-2000, CG K=2, no LLP):

| MAX_PHASE2 | BPE | SPE (sec/edge) | Speedup |
|------------|------|----------------|---------|
| 64 (no pruning) | 2.3295 | 3.595e-06 | 4.3× |
| 16 | 2.3360 | 2.164e-06 | 7.2× |
| 8 | 2.3424 | 1.555e-06 | 10.0× |

MAX=16 provides the best quality/speed tradeoff: 7.2× faster with only +0.007 BPE regression from the analytical-only baseline. The analytical estimators produce identical BPE to the original IOBuffer-based evaluation (2.3295 ≈ 2.3287, difference from Leiden randomness).

Cross-dataset speed results with the optimized encoder:

| Dataset | Edges | SPE (sec/edge) | Config |
|---------|-------|----------------|--------|
| cnr-2000 K=2 | 3.2M | 2.164e-06 | w=64, no LLP |
| in-2004 K=1 | 16.9M | 6.285e-07 | w=8, no LLP |
| enwiki-2013 K=1 (LLP) | 101M | 2.436e-06 | w=64, LLP, iv+LR |
| enwiki-2013 K=1 (original) | 101M | 3.072e-06 | w=64, iv+LR, no reorder |
