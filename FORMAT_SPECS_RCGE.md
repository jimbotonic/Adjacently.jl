# RCGE Format Specification

## Overview

RCGE (Reversible Coarsening Graph Encoding) compresses directed graphs through multi-level hierarchical clustering combined with per-level bitstream encoding. At each coarsening level, the graph is partitioned into clusters, and three types of information are encoded: cluster membership, intra-cluster edges, and inter-cluster edges. The coarsened graph is then recursively partitioned until a stopping condition is met.

Key properties:
- Multi-level coarsening via Leiden community detection
- Cluster membership encoded with Elias-Fano, delta, or STOP-terminated lists
- Intra-cluster edges use optional reference encoding (greedy lookback with bitmap positions + additions)
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
positions_bitmap:
  RLE ones-deltas bitmap over the reference vertex's neighbor list
  indicating which positions are copied
additions:
  IF additions_mode = :delta:
    small_count(num_additions)
    delta-encoded addition list
  IF additions_mode = :intervals:
    small_count(num_runs)
    For each run: (start, length) varints
    small_count(num_singles)
    delta-encoded singles
```

For **non-referenced vertices** (use_ref=false):
```
small_count(degree_in_cluster)
delta-encoded sorted local neighbor IDs
```

#### Reference Encoding Decision

For each vertex, the encoder measures the bitstream cost of:
1. **Raw encoding**: small_count + delta-coded neighbor list
2. **Reference encoding** (for each candidate in lookback window): positions bitmap + additions

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
3. **Cost estimation**: Measure actual bits for positions bitmap + additions vs raw delta encoding
4. **Greedy selection**: Pick the reference yielding the fewest bits, or use raw encoding if no reference helps

### Positions Bitmap

The positions where the current vertex's neighbors match the reference's neighbors are encoded as a dense boolean bitmap over the reference list, using RLE ones-deltas:
- The bitmap has length equal to the reference's neighbor count
- `true` at position j means the j-th neighbor of the reference is also a neighbor of the current vertex

### Additions

Neighbors of the current vertex that are not in the reference list. Encoded as:
- **Delta mode**: small_count + gap-coded sorted list (Fibonacci)
- **Intervals mode**: detected consecutive runs as (start, length) pairs + remaining singles delta-coded

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
| `bits_intra_copy` | Bits for reference positions bitmaps |
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
| Baseline (ref_window=32) | 4.753 | 1,866 KB |
| ref_window=64 | 4.716 | 1,852 KB |
| **ref_window=128 (best)** | **4.693** | **1,842 KB** |
| ref_window=128, LLP 20 passes | 4.693 | 1,842 KB |
| K=4, ref_window=128 | 4.718 | 1,852 KB |
| K=12, ref_window=128 | 4.829 | 1,896 KB |
| K=16, ref_window=128 | 4.844 | 1,902 KB |
| Intervals (mil=2) | 5.062 | 1,987 KB |
| Intervals (mil=3) | 5.011 | 1,967 KB |
| Full MGS per cluster | 5.863 | 2,302 KB |

**Reference**: WebGraph (Zeta-3, LLP) achieves 2.90 BPE on this dataset.

### Intra-Cluster Breakdown (best config, ref_window=128)

| Component | Size | Share |
|-----------|------|-------|
| Headers (ref bitmaps, deltas) | 233 KB | 13.4% |
| Copy (positions bitmaps) | 370 KB | 21.2% |
| Additions (delta residuals) | 776 KB | 44.5% |
| Raw (non-referenced vertices) | 366 KB | 21.0% |

61.6% of vertices (200,408/325,557) use reference encoding with window=128.

### Key Findings

- **Wider reference windows help**: 32→128 saves 1.27% (diminishing returns beyond 128).
- **K=8 is optimal**: Fewer clusters increase intra-cluster encoding cost; more clusters increase membership + inter-cluster overhead.
- **Interval detection hurts**: Two count fields per vertex (num_intervals + num_residuals) add overhead that outweighs savings on cluster-local IDs which lack long consecutive runs.
- **Full MGS per cluster hurts**: Stop values (2 bits/vertex) and reference flags (1 bit/vertex) add ~121KB overhead without compensating gains.
- **Additions dominate**: 44.5% of intra-cluster bytes — the primary optimization target for future work.

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
