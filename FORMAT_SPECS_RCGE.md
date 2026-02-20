# RCGE Format Specification

## Overview

RCGE (Reversible Coarsening Graph Encoding) compresses directed graphs through multi-level hierarchical clustering combined with per-level bitstream encoding. At each coarsening level, the graph is partitioned into clusters, and three types of information are encoded: cluster membership, intra-cluster edges, and inter-cluster edges. The coarsened graph is then recursively partitioned until a stopping condition is met.

Key properties:
- Multi-level coarsening via Leiden community detection
- Cluster membership encoded with Elias-Fano, delta, or STOP-terminated lists
- Intra-cluster edges use optional reference encoding (greedy lookback with copy-blocks positions + additions)
- Inter-cluster edges organized by (source cluster, target cluster) pairs with STOP-terminated neighbor lists
- Raw bitstream output (not wrapped in MGZ container)

## Algorithm Workflow

```
1. Partition graph G into clusters C_1..C_K via Leiden
2. Cap clusters to top-K by size; assign remaining to overflow bin
3. Reorder vertices within each cluster (LLP, RCM, or MinHash)
4. Encode level:
   a. Write cluster membership (sorted vertex lists per cluster)
   b. Write intra-cluster induced edges
   c. Write inter-cluster edge bundles
5. Coarsen: build aggregate graph G' where clusters become super-vertices
6. If |G'.n| > min_clusters and G'.n changed: goto 1 with G = G'
7. Concatenate all level bitstreams
```

## Parameters

`RCGEParams` controls all encoding decisions:

| Parameter | Default | Description |
|-----------|---------|-------------|
| `L` | 32 | Cluster size threshold: clusters with size <= L and undirected use upper-triangular bitset; otherwise per-vertex adjacency lists |
| `varint` | `:fibonacci` | Integer encoding for sizes/lengths (positive-only fields) |
| `count_varint` | `:fibonacci` | Encoding for counts that may be zero |
| `gap` | `:fibonacci` | Encoding for gap-coded sorted lists |
| `degree` | `:golomb` | Encoding for degree vectors |
| `perm_strategy` | `:lehmer` | Permutation encoding strategy (`:lehmer`, `:raw`, `:blockpos`) |
| `undirected_pairs` | `true` | If true and graph is undirected, only encode cluster pairs A < B |
| `membership` | `:stop` | Membership encoding: `:elias_fano`, `:delta`, or `:stop` |
| `inter_strategy` | `:lists` | Inter-cluster encoding: `:perm` or `:lists` |
| `intra_ref_enabled` | `true` | Enable reference encoding for intra-cluster adjacency |
| `intra_ref_window` | 16 | Lookback window size for reference encoding |
| `intra_ref_min_overlap` | 0.3 | Minimum overlap fraction for reference to be considered |
| `intra_ref_rle` | `true` | Use RLE ones-deltas for reference delta vectors |
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

## Bitstream Structure

Each coarsening level produces a bitstream with three sections, written contiguously:

```
┌──────────────────────────────────────────────┐
│ SECTION 1: Cluster Membership                │
│   K (varint): number of clusters             │
│   For each cluster c = 1..K:                 │
│     Sorted vertex list (Elias-Fano / delta   │
│     / STOP-terminated)                       │
├──────────────────────────────────────────────┤
│ SECTION 2: Intra-Cluster Edges               │
│   For each cluster c = 1..K:                 │
│     IF |c| <= L AND undirected:              │
│       Upper-triangular bitset                │
│     ELSE:                                    │
│       [Optional ref bitmap + ref deltas]     │
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

The number of clusters K is written first as a varint. Then for each cluster, the sorted list of global vertex IDs is written using one of three encodings:

**Elias-Fano** (`membership=:elias_fano`): Self-describing — encodes its own length. Most compact for large sorted integer lists.

**Delta** (`membership=:delta`): Explicit length varint, followed by gap-coded (delta-encoded) vertex IDs.

**STOP-terminated** (`membership=:stop`): No explicit length. Each value is preceded by a `1` bit and delta-encoded; the list terminates with a `0` bit.

```
STOP-terminated format:
  For each vertex v in sorted order:
    bit '1'
    varint(v - prev)    # delta from previous value
  bit '0'               # terminator
```

### Section 2: Intra-Cluster Edges

For each cluster, the induced subgraph edges are encoded. Two modes:

#### Small Undirected Clusters (|C| <= L, undirected graph)

Upper-triangular bitset: for all pairs (i, j) where 1 <= i < j <= |C|, emit one bit indicating whether edge (i,j) exists. Total bits: |C| * (|C| - 1) / 2.

#### Large or Directed Clusters

Per-vertex adjacency lists with optional reference encoding. If reference encoding is enabled and beneficial:

**Cluster-level header** (emitted only when at least one vertex uses a reference):
1. **Reference bitmap**: packed boolean vector of length |C|, one bit per vertex indicating whether it uses reference encoding
2. **Reference deltas**: for each vertex with `use_ref=true`, the delta (distance in vertex order) to its reference. Written as:
   - If `intra_ref_rle=true`: RLE ones-deltas encoding
   - Otherwise: count varint + individual deltas

**Per-vertex payload** (in cluster vertex order):

For **referenced vertices** (use_ref=true):
```
positions:
  IF intra_copy_blocks:
    copy-blocks encoding (see Copy-Blocks below)
  ELSE:
    RLE ones-deltas bitmap over the reference vertex's neighbor list
additions:
  IF intra_intervals OR intra_greedy_mil:
    intervals_and_residuals encoding [†]
  ELIF additions_mode = :intervals:
    small_count(num_runs)
    For each run: (start, length) varints
    small_count(num_singles)
    delta-encoded singles [†]
  ELIF intra_stop_deltas:
    STOP-terminated delta list [†]
  ELSE (additions_mode = :delta):
    small_count(num_additions)
    delta-encoded addition list [†]
```

For **non-referenced vertices** (use_ref=false):
```
IF intra_intervals OR intra_greedy_mil:
  intervals_and_residuals encoding [†]
ELIF intra_stop_deltas:
  STOP-terminated delta list [†]
ELSE:
  small_count(degree_in_cluster)
  delta-encoded sorted local neighbor IDs [†]
```

`[†]` When `intra_zigzag=true`, the first value in delta-encoded and intervals_and_residuals lists is encoded as a zigzag offset from the current vertex's local index (see Zigzag Relative First-Value Encoding below), rather than as an absolute value. This does not apply to reference position lists, which are indices into the reference neighbor list.

#### Reference Encoding Decision

For each vertex, the encoder measures the bitstream cost of:
1. **Raw encoding**: delta-coded neighbor list (with STOP-terminated or small_count prefix)
2. **Reference encoding** (for each candidate in lookback window): copy-blocks positions + additions

The cheapest option is selected per vertex. Two-pointer merge is used to compute positions (shared elements with the reference) and additions (elements not in the reference).

### Section 3: Inter-Cluster Edges

Organized by source cluster. For each source cluster A (in order 1..K):

1. **Target list**: STOP-terminated delta-coded list of target cluster indices B that have edges from A.

2. **Per-(A,B) group**: For each target cluster B, active source vertices in A are written as records:
```
For each source vertex u in A with neighbors in B:
  bit '1'                          # active marker
  varint(local_index_of_u_in_A)    # 1-based position in cluster A
  STOP-terminated delta list of B-local neighbor indices
bit '0'                            # end of AB group
```

The B-local neighbor indices are 1-based positions within cluster B's vertex list.

## Reference Encoding Details

The intra-cluster reference encoding exploits similarity between successive adjacency lists within a cluster:

1. **Lookback window**: For vertex at position `idx`, consider reference candidates from position `max(1, idx - window)` to `idx - 1`
2. **Two-pointer merge**: Compute shared positions (elements present in both current and reference lists) and additions (elements only in current list)
3. **Cost estimation**: Measure actual bits for positions (copy-blocks or bitmap) + additions vs raw delta encoding
4. **Greedy selection**: Pick the reference yielding the fewest bits, or use raw encoding if no reference helps

### Copy-Blocks Positions Encoding

When `intra_copy_blocks=true`, reference positions are encoded as WebGraph-style copy-blocks — contiguous runs of copied positions:

```
Copy-blocks format:
  small_count(num_copy_blocks)     # number of contiguous runs
  IF num_copy_blocks > 0:
    varint(first_start)            # start position of first run (1-based, >= 1)
    varint(first_length)           # length of first run (>= 1)
    For each subsequent block i = 2..num_copy_blocks:
      varint(gap)                  # gap from end of previous block (>= 1)
      varint(length)               # length of this run (>= 1)
```

**Example**: Reference list has 10 neighbors. Current vertex shares neighbors at positions 1,2,3,6,7.
- Copy blocks: [(start=1, len=3), (start=6, len=2)]
- Encoded: small_count(2), fib(1), fib(3), fib(2), fib(2)
- gap = 6 - (1+3) = 2

This is more compact than the RLE ones-deltas bitmap because:
1. Skip runs are implicit (encoded by gaps between copy blocks)
2. No per-run type flag (alternation between copy/skip is implicit)
3. No explicit bitmap length needed (decoder knows ref_len from reference list)
4. All values >= 1, avoiding +1 shifts needed for zero-tolerant block encoding

### RLE Ones-Deltas Bitmap (Legacy)

When `intra_copy_blocks=false`, positions are encoded as a dense boolean bitmap over the reference list:
- `varint(token_count)`: number of alternating runs
- For each token: 1-bit flag (is_ones_run) + `varint(run_length)`
- The bitmap has length equal to the reference's neighbor count
- `true` at position j means the j-th neighbor of the reference is also a neighbor of the current vertex

### Additions

Neighbors of the current vertex that are not in the reference list. Encoded as:
- **STOP-terminated** (`intra_stop_deltas=true`): Each value preceded by `1` bit, delta-coded, terminated by `0` bit. No explicit count prefix.
- **Delta mode** (`additions_mode=:delta`): small_count + gap-coded sorted list (Fibonacci)
- **Intervals mode** (`intra_intervals` or `intra_greedy_mil`): intervals_and_residuals encoding (interval start/length pairs + delta-coded residuals)
- **Custom intervals mode** (`additions_mode=:intervals`): detected consecutive runs as (start, length) pairs + remaining singles delta-coded

### STOP-Terminated Delta Lists

When `intra_stop_deltas=true`, neighbor lists (both raw and additions) use STOP-terminated encoding instead of small_count + delta:

```
STOP-terminated delta list:
  For each value v in sorted order:
    bit '1'                        # more values
    varint(gap)                    # v - prev (first value: absolute or zigzag-encoded)
  bit '0'                          # STOP terminator
```

This eliminates the per-list count prefix. Cost: 1 bit per value + 1 stop bit, vs small_count (2 bits for degree 0-2, 2+varint for degree >= 3). STOP wins for degree 0-5, ties at 6, loses at 7+. Net saving: ~0.05 BPE on CNR-2000.

### Zigzag Relative First-Value Encoding

When `intra_zigzag=true`, the first value in neighbor lists and residual lists is encoded as a signed offset from the current vertex's local index within the cluster, rather than as an absolute value. This exploits the locality created by LLP reordering: neighbors tend to cluster near the vertex in the ordering, so the offset is typically small.

**Encoding**: For vertex at local index `v` with first neighbor value `v1`:
```
offset = v1 - v
zigzag(offset) = offset >= 0 ? 2*offset : 2*(-offset) - 1
encoded = zigzag(offset) + 1     # +1 shift for Fibonacci/Zeta compatibility (must be >= 1)
```

**Decoding**: Read raw value, then recover:
```
offset = zigzag_decode(raw - 1)
v1 = v + offset
```

The zigzag mapping converts signed offsets to positive integers: 0->0, -1->1, +1->2, -2->3, +2->4, ...

**Scope**: Zigzag applies to:
- First value in delta-encoded neighbor lists (raw and additions)
- First interval start in intervals_and_residuals encoding
- First value in residual sublists within intervals_and_residuals encoding

Zigzag does **not** apply to reference position lists (these are indices into the reference's neighbor list, not local vertex IDs).

## Integer Encodings

| Encoding | Zero Support | Description |
|----------|-------------|-------------|
| Fibonacci | No (values >= 1) | Zeckendorf representation terminated by `11` |
| Elias Delta | No (values >= 1) | Elias gamma prefix + binary suffix |
| Elias Fano | Yes | Monotone integer sequence encoding with high/low bit split |
| Golomb | Yes (values >= 0) | Base-b quotient/remainder coding |
| Small-count | Yes | 2-bit escape: 00/01/10 for 0/1/2, 11 + varint for >= 3 |

For encodings that do not support zero (Fibonacci, Elias Delta), zero-valued fields are stored shifted by +1.

## Multi-Level Coarsening

The multi-level loop proceeds as follows:

1. **Level 1**: Partition original graph via Leiden (max_passes=8, max_levels=5). Cap to K clusters (default K=8).
2. **Deeper levels**: Re-partition with Leiden (max_passes=5, max_levels=5). K adapted: `K = max(16, min(K, ceil(nclusters/2)))`.
3. **Coarsening**: `aggregate_graph` produces a weighted graph where each cluster becomes a super-vertex and edge weights count multi-edges.
4. **Stopping conditions**:
   - Coarse graph size unchanged from previous level
   - Coarse graph has <= `min_clusters` (default 32) super-vertices
   - Maximum level count reached (default 5)

The coarse weighted graph is converted to a `TestGraph` (unweighted with repeated edges matching weights) for the next encoding level.

## Statistics Tracking

`RCGEStats` accumulates per-section bit counts during encoding:

| Field | Description |
|-------|-------------|
| `bits_membership` | Bits for Section 1 (cluster membership lists) |
| `bits_intra` | Total bits for Section 2 (intra-cluster edges) |
| `bits_intra_headers` | Bits for intra-cluster structural headers (ref bitmaps, deltas, raw counts) |
| `bits_intra_copy` | Bits for reference copy-blocks / positions bitmaps |
| `bits_intra_add` | Bits for reference additions |
| `bits_intra_raw` | Bits for non-referenced vertex adjacency data |
| `bits_intra_ref_small_headers` | Bits for small reference headers |
| `intra_ref_used` | Count of vertices using reference encoding |
| `intra_no_ref` | Count of vertices using raw encoding |
| `bits_inter_headers` | Bits for inter-cluster target lists |
| `bits_inter_degrees` | Bits for inter-cluster degree vectors |
| `bits_inter_perms` | Bits for inter-cluster permutation data |
| `bits_inter_lists` | Bits for inter-cluster neighbor lists |

## Performance

### CNR-2000 Benchmark

Graph: 325,557 vertices, 3,216,152 edges. Partitioned via Leiden with K=8 clusters, LLP reordering (5 passes).

| Configuration | Bits/Edge | File Size |
|---------------|-----------|-----------|
| Baseline (ref_window=32, delta) | 4.307 | 1,691 KB |
| Zigzag | 3.566 | 1,400 KB |
| Zigzag + STOP deltas | 3.515 | 1,380 KB |
| **Zigzag + STOP + CopyBlocks** | **3.117** | **1,224 KB** |
| GreedyMIL (intervals, no zigzag) | 4.998 | 1,962 KB |
| GreedyMIL+Zigzag | 4.097 | 1,609 KB |
| Intervals (mil=2) | 5.062 | 1,987 KB |
| Intervals (mil=3) | 5.011 | 1,967 KB |
| Full MGS per cluster | 5.863 | 2,302 KB |

**Reference**: WebGraph (Zeta-3, LLP) achieves 2.90 BPE on this dataset.

### Intra-Cluster Breakdown

**Best config (Zigzag + STOP + CopyBlocks)**:

| Component | Size | BPE | Share |
|-----------|------|-----|-------|
| Membership | 125 KB | 0.318 | 10.2% |
| Headers (ref bitmaps, deltas) | 200 KB | 0.510 | 16.3% |
| Copy (copy-blocks positions) | 268 KB | 0.683 | 21.9% |
| Additions (STOP delta residuals) | 460 KB | 1.170 | 37.5% |
| Raw (non-referenced vertices) | 159 KB | 0.406 | 13.0% |
| Inter-cluster | 12 KB | 0.031 | 1.0% |

62.0% of vertices (201,902/325,557) use reference encoding.

**Zigzag + STOP baseline (before copy-blocks)**:

| Component | Size | BPE | Share |
|-----------|------|-----|-------|
| Membership | 125 KB | 0.318 | 9.1% |
| Headers (ref bitmaps, deltas) | 203 KB | 0.516 | 14.7% |
| Copy (RLE bitmap positions) | 343 KB | 0.873 | 24.8% |
| Additions (STOP delta residuals) | 485 KB | 1.236 | 35.1% |
| Raw (non-referenced vertices) | 213 KB | 0.542 | 15.4% |
| Inter-cluster | 12 KB | 0.031 | 0.9% |

60.4% of vertices (196,654/325,557) use reference encoding.

### Optimization History

| Change | BPE | Delta | Cumulative |
|--------|-----|-------|------------|
| Baseline (ref+delta) | 4.307 | — | — |
| + Zigzag first-value | 3.566 | -0.741 | -17.2% |
| + STOP-terminated deltas | 3.515 | -0.051 | -18.4% |
| + Copy-blocks positions | 3.117 | -0.398 | -27.6% |
| Gap to WebGraph (2.90) | | 0.217 | |

### Key Findings

- **Copy-blocks is the second most impactful optimization**: 3.117 BPE vs 3.515 with STOP (-0.398 BPE, -11.3%). By encoding reference positions as contiguous copy runs instead of RLE bitmaps, the format eliminates per-run type flags, implicit skip lengths, and explicit bitmap lengths. Copy positions shrink 21.7% (343 KB -> 268 KB). Additionally, 5,248 more vertices switch to reference encoding (improving raw by 25%), because the more accurate copy-blocks cost estimation reveals references that were marginally rejected under the old bitmap cost model.
- **Zigzag is the single most impactful optimization**: 3.566 BPE vs 4.307 baseline (-17.2%). By encoding the first neighbor value as a signed offset from the vertex's local index, the encoder exploits LLP locality. Raw encoding shrinks by 40% (372 KB -> 221 KB), additions by 22% (626 KB -> 489 KB).
- **STOP-terminated deltas provide modest savings**: 3.515 vs 3.566 (-0.051 BPE, -1.4%). Eliminates small_count prefix (2-6 bits) at the cost of 1 bit per value + 1 stop bit. Wins for degree 0-5 vertices, which are common in intra-cluster lists.
- **GreedyMIL (interval encoding) hurts on short lists**: 4.998 BPE, 16% worse than baseline. The intervals_and_residuals format has per-list overhead (num_intervals + num_residuals counts) that isn't amortized on RCGE's typically short intra-cluster lists.
- **Global LLP before Leiden hurts RCGE**: +0.33 BPE. Leiden on LLP-relabeled graph produces non-contiguous clusters; inter-cluster edges explode 5x (12 KB -> 64 KB).
- **Zeta-3 slightly worse than Fibonacci for RCGE**: Moderate-sized local IDs (1-36K) within clusters favor Fibonacci's fixed-width codewords.
- **K=8 is optimal**: Fewer clusters increase intra-cluster encoding cost; more clusters increase membership + inter-cluster overhead.
- **Recommended config**: `intra_zigzag=true, intra_stop_deltas=true, intra_copy_blocks=true` with moderate reference window (`intra_ref_window=32`).

## Example

### Small Graph

Consider a directed graph with 8 vertices partitioned into 3 clusters:
- C1 = {1, 2, 3} (vertices within a web domain)
- C2 = {4, 5, 6} (another domain)
- C3 = {7, 8} (small cluster)

Edges:
- Intra C1: 1->2, 2->3, 3->1
- Intra C2: 4->5, 5->6
- Inter C1->C2: 2->4, 3->5
- Inter C2->C3: 6->7
- Inter C3->C1: 8->1

### Encoded Bitstream (Sketch)

```
SECTION 1: Membership
  K = 3
  C1: [1, 2, 3]  (Elias-Fano or delta-coded)
  C2: [4, 5, 6]
  C3: [7, 8]

SECTION 2: Intra-Cluster
  C1 (directed, 3 vertices):
    v1 (vertex 1): count=1, neighbors=[2]  (local ID delta-coded)
    v2 (vertex 2): count=1, neighbors=[3]
    v3 (vertex 3): count=1, neighbors=[1]
  C2:
    v1 (vertex 4): count=1, neighbors=[2]
    v2 (vertex 5): count=1, neighbors=[3]
    v3 (vertex 6): count=0
  C3:
    v1 (vertex 7): count=0
    v2 (vertex 8): count=0

SECTION 3: Inter-Cluster
  A=C1: targets=[C2]  STOP
    (C1,C2) group:
      '1' u_local=2 neighbors_in_C2=[1]  STOP  # vertex 2 -> vertex 4 (C2 local 1)
      '1' u_local=3 neighbors_in_C2=[2]  STOP  # vertex 3 -> vertex 5 (C2 local 2)
      '0'  # end group

  A=C2: targets=[C3]  STOP
    (C2,C3) group:
      '1' u_local=3 neighbors_in_C3=[1]  STOP  # vertex 6 -> vertex 7 (C3 local 1)
      '0'

  A=C3: targets=[C1]  STOP
    (C3,C1) group:
      '1' u_local=2 neighbors_in_C1=[1]  STOP  # vertex 8 -> vertex 1 (C1 local 1)
      '0'
```

After encoding, the coarse graph has 3 super-vertices with edges C1->C2, C2->C3, C3->C1. If this is above `min_clusters`, coarsening continues recursively.
