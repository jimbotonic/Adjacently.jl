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

**Best result on CNR-2000**: **2.4341 BPE** — beats WebGraph BV (2.897 BPE) by 0.463 BPE with all roundtrips verified.

The 2.4341 BPE result combines **implicit membership** (saves 0.221 BPE), **adaptive additions** (saves 0.067 BPE), **adaptive raw** (saves 0.049 BPE), and **3-way adaptive copy positions** (bitmap / copy-blocks / complement per ref vertex, saves 0.116 BPE total). All optimizations are independent and their savings are additive.

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

**Recommended best config** (2.4341 BPE on CNR-2000 with K=1 + implicit ranges + all adaptive):
```julia
# Step 1: Leiden K=1 partition → 2 groups, then sort each by vertex ID
clusters_by_id = [sort(group1_vertices), sort(group2_vertices)]
# Step 2: Relabel vertices by vertex-ID rank within each group
g_rel = relabel_graph(g, vertex_map)  # vertex_map: v → rank in ID-sorted group
clusters_contiguous = [TV.(1:S1), TV.(S1+1:n)]
# Step 3: Encode with implicit_ranges + all adaptive per-vertex codec choices
RCGEParams(
    L=128, membership=:implicit_ranges, inter_strategy=:perm,
    intra_ref_window=64, intra_ref_rle=false, intra_ref_fixwidth=true,
    intra_zigzag=true, intra_stop_deltas=true, intra_copy_blocks=true,
    intra_copy_adaptive=true,
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
| `intra_copy_bitmap_count` | Number of ref vertices that chose bitmap mode (when `intra_copy_adaptive`) |
| `intra_copy_blocks_count` | Number of ref vertices that chose copy-blocks mode (when `intra_copy_adaptive`) |
| `intra_copy_complement_count` | Number of ref vertices that chose complement mode (when `intra_copy_adaptive`) |
| `intra_overlap_hist` | 10-bucket histogram of overlap fractions [0%,10%)…[90%,100%) for ref vertices |
| `intra_add_count_total` | Total number of additions across all ref vertices |

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
| + Adaptive adds+raw (MIL=2) | 1 | 2.550 | −0.116 | 1,001 KB |
| + 2-way adaptive copy (bitmap vs copy-blocks) | 1 | 2.497 | −0.053 | 981 KB |
| **+ 3-way adaptive copy (+ complement)** | **1** | **2.4341** | **−0.063** | **956 KB** |

**Total improvement**: 4.307 → 2.4341 BPE (−43.5% over 11 independent optimizations).

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

### Comparison with WebGraph BV

| Component | RCGE K=1 + 3-Way Adaptive | RCGE K=1 + 2-Way Adaptive | RCGE K=1 + Adds/Raw Adaptive | RCGE K=1 + EF | WebGraph BV |
|-----------|---------------------------|---------------------------|------------------------------|---------------|-------------|
| Membership / structure | 0.000 BPE | 0.000 BPE | 0.000 BPE | 0.221 BPE | ~0.078 BPE |
| Headers | 0.481 BPE | 0.481 BPE | 0.481 BPE | 0.481 BPE | ~0.126 BPE |
| Copy positions | 0.572 BPE | 0.635 BPE | 0.688 BPE | 0.688 BPE | ~0.556 BPE |
| Residuals/additions | 1.052 BPE | 1.052 BPE | 1.052 BPE | 1.119 BPE | ~1.583 BPE |
| Raw / non-referenced | 0.329 BPE | 0.329 BPE | 0.329 BPE | 0.377 BPE | ~0.554 BPE |
| Inter-cluster | 0.001 BPE | 0.001 BPE | 0.001 BPE | 0.001 BPE | — |
| **Total** | **2.4341 BPE** | **2.497 BPE** | **2.550 BPE** | **2.887 BPE** | **2.897 BPE** |

RCGE + 3-way adaptive encoding beats WebGraph BV by **0.463 BPE**. The reference encoding (copy + headers) is much more effective in RCGE because the two-step Leiden + per-cluster LLP ordering creates dense locality — 62.5% of vertices can reference a recent neighbor, vs ~56% for WebGraph's bisection ordering. The 3-way adaptive copy positions (bitmap / copy-blocks / complement) reduces the copy section from 0.688 to 0.572 BPE: complement encoding dominates for the 49% of ref vertices with ≥90% overlap, where encoding 2–3 skipped positions costs far less than a full copy-blocks list of 18–20 entries.

### Greedy Approach Comparison

| Approach | Ordering | Window | BPE |
|----------|----------|--------|-----|
| Greedy (baseline) | None (crawl order) | 7 | ~4.3 |
| Greedy | Global LLP | 7 | 3.561 |
| Greedy | Global LLP | 64 | 3.319 |
| **Greedy (best)** | **Leiden + per-group LLP** | **64** | **2.881** |
| WebGraph BV | Greedy bisection | 7 | 2.897 |
| RCGE (K=1, EF) | Leiden + per-cluster LLP | 64 | 2.887 |
| RCGE (K=1, implicit ranges) | Leiden + per-cluster LLP | 64 | 2.666 |
| RCGE (K=1, adaptive adds+raw) | Leiden + per-cluster LLP | 64 | 2.550 |
| RCGE (K=1, 2-way adaptive copy) | Leiden + per-cluster LLP | 64 | 2.497 |
| **RCGE (K=1, 3-way adaptive copy)** | **Leiden + per-cluster LLP** | **64** | **2.4341** |

The greedy encoder now beats WebGraph BV (2.881 vs 2.897 BPE) with its best configuration (adaptive_copy + stop_deltas + empty_prefix + compact_copy + vlc2 + tight_intervals). The remaining 0.45 BPE gap to RCGE is from RCGE's per-cluster local vertex indexing and hierarchical encoding, which create fundamentally tighter locality than the single-pass greedy approach.

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

11. **2-way adaptive copy-position encoding** (−0.053 BPE). Per-ref-vertex choice between copy-blocks and a dense bitmap over the reference's neighbor list (one bit per reference neighbor, 1=copied). Bitmap costs exactly `ref_len` bits, which beats copy-blocks when the reference list is short (≤ ~8 neighbors) or when copied positions are scattered. Mode flag overhead: 1 bit × 203,615 ref vertices = 25 KB; gross copy savings: 47 KB; net: −21 KB = −0.053 BPE.

12. **Window size w=64 is optimal** (larger windows are harmful). Increasing to w=128 adds 1 extra bit per ref-delta (+26 KB headers) but saves only 23 KB in additions+raw (0.4% more ref vertices found). Net: −3 KB (slightly worse). w=256 even more so: +53 KB headers, −44 KB adds+raw, net −9 KB. The additional references found beyond the 64-vertex window are too few to justify the wider fixed-width delta encoding.

13. **Complement copy-blocks (3-way adaptive)** (−0.063 BPE, copy 255 KB → 230 KB). A third copy mode encodes the SKIPPED positions instead of the copied ones, using the same copy-blocks format. For high-overlap references (≥90% copied), complement is dramatically cheaper: encoding 1–2 skipped positions costs ~4 bits vs encoding 18–20 copied positions in a run. The 3-way mode uses nested bits: outer=1→bitmap (1 bit overhead), outer=0+inner=0→copy-blocks (2 bits), outer=0+inner=1→complement (2 bits). Despite the 1-bit additional overhead for non-bitmap modes, complement dominates for the 16.1% of ref vertices (32,769 total) where it's cheapest. The overlap histogram shows why: 49% of ref vertices (99,760) have ≥90% overlap fraction — the main beneficiaries. Net savings: 255 KB − 230 KB = 25 KB beyond the 2-way baseline.

14. **Smart ref-delta encoding (delta=1 shorthand) is harmful** (+0.040 BPE). The intuition was that delta=1 (reference immediately preceding vertex) would dominate and could be encoded as just 2 bits instead of 7. In practice, **only 6% of ref vertices use delta=1**. The actual distribution is near-uniform over [1,64]: delta=49–64 is the largest bucket at 35.8%, and small deltas (1–4) account for only 17% combined. This confirms finding #10: the LLP-ordered large top-level group has useful references spread throughout the entire 64-vertex window, not concentrated at delta=1. Encoding delta=1 in 2 bits saves 5 bits for 6% of refs but costs 1 extra bit for the remaining 94%, netting +0.64 bits/ref × 203K refs = +16 KB ≈ +0.040 BPE.

15. **Elias-Fano for copy positions is harmful** (+0.034 BPE). EF was expected to help for large `ref_len` with scattered copied positions (e.g., k=10 scattered copies in ref_len=64 costs ~51 bits in EF vs 64 bits bitmap vs ~95 bits copy-blocks). In practice, EF is chosen for only **0.4% of ref vertices** (753 out of 203,615) — far too few to compensate for the mode-bit overhead. The flat 2-bit mode scheme (to accommodate a 4th mode) adds 1 extra bit for all 124,864 bitmap vertices (+15.6 KB), while 753 EF-mode vertices save only ~2 KB. The root cause: most ref vertices fall into the three cases already handled well by the existing 3-way adaptive — small ref_len (bitmap wins), high overlap (complement wins), or run-structured positions (copy-blocks wins). The large-ref_len + scattered regime where EF would excel is rare in this dataset.
