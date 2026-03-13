# ASTRA-L File Format Specification

## Overview

ASTRA-L (Adaptive Streaming Adjacency — Layered) encodes a directed graph by partitioning vertices into levels derived from PageRank-guided radius-k ball explorations. Within each level, vertices are locally renumbered and intra-level edges are compressed using full MGS encoding (reference + interval + recursive reference). Cross-level edges are collected separately and encoded with delta-coded source/target lists.

Key properties:
- **Full-ball levels**: Each level contains all unassigned vertices from a radius-k BFS ball (not just the SCC). This produces fewer, larger levels with high intra-level edge density.
- **MGS compression for intra-level edges**: Each level's adjacency is compressed using the full MGS encoding pipeline — reference encoding with sliding window, interval detection, and recursive reference — applied to local IDs.
- **Compact cross-level encoding**: Cross-level edges are grouped by source vertex and delta-encoded, eliminating the expensive reconciliation dictionary.
- **Single-pass encoding**: No backward pass on the reverse graph (eliminated due to negligible benefit).
- **Small level merging**: Balls with fewer than `min_level_size` unassigned vertices are merged into the previous level.
- Seeds are picked by PageRank; frontier vertices beyond the ball are added to the priority queue for better coverage.

## Parameters

| Parameter | Default | Description |
|-----------|---------|-------------|
| `radius` (k) | 2 | Radius of BFS ball around each seed vertex |
| `damping` | 0.85 | PageRank damping factor |
| `epsilon` | 1e-6 | PageRank convergence threshold |
| `integer_encoding` | `:fibonacci` | Integer codec for all encoded values (`:fibonacci`, `:elias_delta`, etc.) |
| `ref_window_size` | 7 | Reference window size for intra-level MGS compression |
| `use_scc` | `false` | If true, extract SCC from ball (legacy); if false, use full ball |
| `min_level_size` | 1 | Minimum level size; smaller balls are merged into previous level |
| `log_every` | 1000 | Print progress every N levels (0 = silent) |

## Preliminaries

- Graph: directed, vertices have immutable global IDs in [1..n].
- PageRank: compute PageRank scores for all vertices upfront; maintain PR as a static ranking used for seed selection and tie-breaking.
- Radius k: user- or codec-chosen small integer (typically k = 2 or 3) used to explore balls around seeds.

## High-Level Workflow

1. Compute PageRank for all vertices; build `pagerank[global_id]`.
2. Iteratively build levels until all vertices are assigned:
   - Pick highest-PR unassigned vertex as seed (from frontier PQ or global PR order).
   - Explore radius-k ball B_k(seed) in G.
   - Filter to unassigned vertices in the ball; if too few, merge into previous level.
   - Assign contiguous local IDs 1..|L| sorted by descending PR.
   - Partition edges into intra-level (both endpoints in L) and cross-level (target outside L).
   - Add frontier vertices (1-hop beyond ball) to the seed priority queue.
3. Write output:
   - For each level: seed ID, level size, then intra-level adjacency compressed via full MGS encoding (reference + interval + recursive reference) on local IDs.
   - Cross-level edges: grouped by source vertex, delta-encoded source IDs and target lists.

## Data Tracked During Encoding

- `map_gl_to_local`: global_id → (level_id, local_id) for vertices assigned to current level.
- `map_local_to_gl`: (level_id, local_id) → global_id.
- `encoded_edges_set`: set of edges in global IDs already emitted as local adjacency (to avoid duplication across passes).
- `reconciliation`: list of tuples `(global_id, (level_x, local_x), (level_y, local_y))` for cross-level edges where the target is in a previous level.
- `remaining_cross_edges`: list of edges kept in global IDs when they cannot be encoded with local IDs yet (e.g., target level not finalized) or by design.
- `explored`: set of vertices seen during radius-k explorations; used to select next seeds.

## Level Construction (Forward Pass)

Given a current seed `s`:
1. Explore radius-k ball B_k(s) in G (e.g., BFS/DFS limited to k steps from `s`).
2. Restrict to the induced subgraph on B_k(s) and extract its SCC decomposition.
3. Select the SCC containing `s` as the current level L. If the induced subgraph is a tree (including the degenerate case of a single sink), merge its vertices into the previous level using the previous level’s current local IDs (no new level ID). Otherwise, create a new level for this SCC.
4. Assign local IDs 1..|L| and produce `map_gl_to_local` and `map_local_to_gl` for L.
5. Encode intra-level edges: for each `u ∈ L`, emit adjacency to `v ∈ L` using local IDs. Add `(u,v)` to `encoded_edges_set` in global IDs.
6. Record seed parents from the previous level: list local IDs in the previous level that have edges to `s`. Emit this list together with the level payload.
7. Cross-level edges from L to earlier levels: for any edge `(u,v)` with `v` in a previous level, add a reconciliation entry `(global(v), (level(v), local(v)), (level(u), local(u)))`. Do not re-encode the edge in local IDs.
8. Edges from L to not-yet-materialized levels: keep `(global(u), global(v))` in `remaining_cross_edges` (to be handled later if they never become intra-level in any pass).
9. Update `explored` with vertices in B_k(s) and L.

### Choosing the Next Seed (Forward Pass)

Among vertices in `explored` that are not part of the current level L, choose the vertex with the highest PageRank as the next seed `s_next`. Before finalizing level L, encode the parents of `s_next` that belong to L using L’s local IDs; this helps cheap cross-level connectivity at decode time. Then restart the radius-k exploration from `s_next` to build the next level.

### Tree Merge Heuristic

When the subgraph induced by the exploration is a tree (including a single sink), merge all vertices of the tentative new level into the previous level, re-using the previous level’s local IDs for already-present vertices and assigning local IDs for new vertices as needed, but without introducing a new level ID. This reduces per-level overhead on tree-like areas.

## Backward Pass (Reverse Graph)

Run the same procedure on the reversed graph G^R:
- Recompute levels by SCC-in-ball around PR-guided seeds (same k).
- Encode intra-level edges in local IDs exactly as for the forward pass.
- Skip edges already in `encoded_edges_set` (in global IDs) to avoid duplicates across passes.
- Any local edges encoded in the backward pass are reversed during decoding to match G.
- Maintain `reconciliation` and `remaining_cross_edges` with the same rules as in the forward pass (when applicable for edges landing in previous levels of the reverse layering).

## Final Assembly (What the File Encodes)

The ASTRA-L payload (inside the MGS container) encodes, in order:
1. Forward pass levels:
   - For each level: seed identifier, parent-local-IDs for the next seed (if any), and intra-level edges in local IDs.
2. Backward pass levels:
   - For each level: seed identifier and intra-level edges in local IDs not already encoded; these edges are marked as “reverse” and will be flipped during decoding.
3. Reconciliation dictionary:
   - Entries of the form `(global_id, (level_x, local_x), (level_y, local_y))` linking cross-level endpoints when an edge of the current level points into a previous level.
4. Remaining cross-level edges:
   - Any edges that could not be encoded with local IDs in either pass are listed in global IDs. Optionally, they can be bucketed by distance and direction (e.g., forward distance 1 from level t to t+1, backward distance 1 from t to t-1, distance 2, etc.).

All sections use STOP-delimited lists as in other ASTRA variants. Integer encoding and delta/interval schemes are identical to MGS/ASTRA defaults unless stated otherwise.

## File Structure

```
┌──────────────────────────────────────────────────────────────┐
│ LEVEL COUNT (varint)                                          │
├──────────────────────────────────────────────────────────────┤
│ LEVELS (repeated level_count times)                           │
│  For each level:                                              │
│    - Seed ID (varint, global ID)                              │
│    - Level size (varint, number of vertices)                  │
│    - Intra-level adjacency (full MGS compressed data)         │
│      • reference flag (1 bit)                                 │
│      • per-vertex: reference/interval/hybrid encoding         │
│      • STOP-delimited (:children mode)                        │
├──────────────────────────────────────────────────────────────┤
│ CROSS-LEVEL EDGES                                             │
│  - Active source count (varint)                               │
│  - For each active source (sorted by global ID):              │
│    - Delta source ID (varint, delta from previous source)     │
│    - Target count (varint)                                    │
│    - Delta-encoded sorted target list (hybrid mix)            │
│  - STOP (end marker)                                          │
└──────────────────────────────────────────────────────────────┘
```

Notes:
- Seed ID is always a global vertex ID.
- The intra-level adjacency block uses the same encoding as standard MGS files: reference encoding with copy bitmaps, interval detection, and recursive reference — all applied to local IDs (1..|L|).
- Cross-level edges use delta-coded source IDs to avoid encoding large absolute global IDs.

## Encoding Details

### Local ID Assignment
- Within each level, assign contiguous local IDs 1..|L| in a deterministic order (e.g., by decreasing PageRank, then by global ID). Determinism ensures stable delta encoding and reproducibility.

### Adjacency Lists
- For each local source ID, list local target IDs within the same level only.
- Sort targets by local ID, delta-encode consecutive differences, then apply the integer codec (e.g., Fibonacci) as in ASTRA.
- Use STOP to delimit each source’s adjacency and to end the level’s adjacency block.

### Integer and Difference Encoding
- Integers use Fibonacci encoding by default (self-delimiting).
- Target lists are sorted and delta-encoded.
- Optional interval/run-length compression may be used within a level consistent with existing ASTRA variants.

## Reconciliation Dictionary

Record an entry whenever an edge of the current level points to a vertex that already belongs to a previous level. Each entry captures the vertex identity and the involved local coordinates:

```
(global_id, (level_x, local_x), (level_y, local_y))
```

- `global_id`: the global ID of the target (or source) vertex being reconciled.
- `(level_x, local_x)`: the coordinates of the endpoint in the previous level.
- `(level_y, local_y)`: the coordinates of the endpoint in the current level.

The decoder uses these entries to stitch cross-level edges without needing to re-encode them in local form.

## Remaining Cross-Level Edges (Global IDs)

Some edges cannot be encoded in local IDs in either pass (e.g., edges whose endpoints never co-occur within the same level/SCC, or edges deferred whose levels never materialize). These are serialized in `(src_global, dst_global)` form after the reconciliation section. If desired, they can be additionally partitioned into forward/backward and by graph-theoretic distance buckets to aid specialized post-processing.

## Decoding Overview

1. Read and reconstruct the forward pass levels:
   - Recover local-to-global mappings per level from the level headers or the encoded vertex lists.
   - Build intra-level edges directly from local adjacency lists.
   - Attach parent-local-IDs information for the next seed (helps reconstruct short cross-level links if the implementation uses it for hints).
2. Read and reconstruct the backward pass levels:
   - Add intra-level edges (reverse-marked) after reversing them to match G.
   - Skip edges already present to avoid duplicates.
3. Apply reconciliation entries to create cross-level edges that point into previous levels.
4. Append remaining cross-level edges listed in global IDs.

The final graph is the union of:
- Encodable intra-level edges from the forward pass (including merged tree components) plus parents-of-seed edges captured as local references at level boundaries.
- Encodable intra-level edges from the backward pass that were not already encoded (reversed back to G’s orientation).
- All remaining cross-level edges in global IDs.

## Example (Sketch)

Suppose k=2 and the process builds Level 0 as the SCC around seed s0, Level 1 around seed s1, etc. Within each level, vertices are locally numbered 1..n. For Level t we emit:
- Seed s_t (global).
- Parents of s_{t+1} that lie in Level t as local IDs.
- Intra-level edges of Level t using local IDs.
If an edge from Level t points into Level t-1, we add a reconciliation tuple for the two endpoints. Any edge from Level t to future levels stays in the global edge list until it either becomes intra-level in a later pass or remains as a global edge.

## Performance

### CNR-2000 Benchmark

| Metric | Value |
|--------|-------|
| Graph | 325,557 vertices, 3,216,152 edges |
| Radius | k=2, full-ball mode, min_level_size=50 |
| Integer encoding | Fibonacci, ref_window=7 |
| Bits per edge | 7.549 |
| Output format | `.astral` (standalone bitstream, not MGZ container) |

### Comparison

| Method | Bits/Edge |
|--------|-----------|
| WebGraph (Zeta-3, LLP) | 2.90 |
| CG (Leiden + Elias-Fano) | 4.75 |
| Legacy ASTRA (MGS V2) | 5.11 |
| ASTRA-L (SCC-in-ball, old) | 450 (broken — see below) |
| ASTRA-L (full-ball + MGS, new) | 7.55 |

### Legacy SCC-in-Ball Issues

The original ASTRA-L implementation used SCC extraction from radius-2 balls, producing 18,500+ tiny levels (avg 17 vertices). This caused:
1. **Reconciliation explosion**: 13.7M entries at 5 integers each dominated the file
2. **No reference encoding**: Intra-level adjacency used only basic hybrid-mix
3. **Backward pass waste**: Full reverse-graph decomposition added ~60s and megabytes for near-zero new edges

### Decoder Status

`read_astra_layered_graph` is currently a stub (returns `nothing`). Roundtrip decoding is not yet implemented.

## Implementation Notes

- PageRank is computed once and cached for the entire run; ties may be broken by global ID.
- Radius-k exploration can be BFS layered by distance; SCC extraction is applied to the induced subgraph of the explored ball.
- The tree-merge heuristic significantly reduces overhead on acyclic regions; merging reuses previous level’s local IDs.
- Maintain `encoded_edges_set` in global IDs to avoid duplicate emission between forward and backward passes.
- Deterministic local ordering improves delta coding and reproducibility.

## Targets and Expectations

- Forward/backward split increases local encodability and shrinks the remainder list.
- Reconciliation entries are sparse relative to total edges and carry only coordinates, not full adjacency.
- Remaining cross-level edges may be further stratified to enable domain-specific optimizations, but this is optional.

## Concrete Example (k=2)

Graph (global IDs 1..7):
- Edges within SCC A: 1→2, 2→3, 3→1
- Cross edges from A: 2→4, 3→6
- SCC B: 4→5, 5→4
- Edges around tree: 5→6, 6→7
- Back edge to previous level: 5→1

Assumptions:
- PageRank order (for seed picking/ties): 3 > 1 > 2 > 5 > 4 > 6 > 7
- k = 2 (radius-2 ball)
- Local IDs within a level assigned by ascending global ID (deterministic for readability)

Forward pass (G):
- Level 0 from seed s0=1
  - B_2(1) includes {1,2,3,4}. SCC containing 1 is {1,2,3} → L0 = {1,2,3}
  - Local IDs: 1→1, 2→2, 3→3
  - Intra-level edges (local): 1→2, 2→3, 3→1
  - Explored (candidates for next seeds): {1,2,3,4}
  - Next seed s1 is 4 (highest PR among explored not in L0)
  - Parents of s1 in L0 (local IDs): {2} because 2→4 exists
  - Reconciliation from L0: none (no L0 edge points to a previous level)

- Level 1 from seed s1=4
  - B_2(4) includes {4,5,1,6,7} via 4→5, 5→4, 5→1, 5→6, 6→7
  - SCC containing 4 is {4,5} initially
  - Tentative next (remaining) region {6,7} is acyclic (a tree), so merge {6,7} into L1 (tree-merge heuristic)
  - Final L1 vertices: {4,5,6,7}
  - Local IDs: 4→1, 5→2, 6→3, 7→4
  - Intra-level edges (local): 1→2 (4→5), 2→1 (5→4), 2→3 (5→6), 3→4 (6→7)
  - Edge 5→1 points to previous level (L0), so add reconciliation:
    - (global_id=1, (level_0, local=1), (level_1, local=2))
  - No remaining vertices; forward pass ends

Backward pass (G^R):
- L0 reversed intra-edges would be 2→1, 3→2, 1→3 but these pairs are already represented by L0 forward edges; skip as duplicates
- L1 reversed intra-edges would be 2→1, 1→2, 3→2, 4→3 which reverse back to edges already encoded in forward pass; skip
- Backward pass contributes no new local edges in this example

Final assembly (what the file effectively contains):
- Forward local edges
  - Level 0 (L0, local IDs {1,2,3} for {1,2,3}): 1→2, 2→3, 3→1
  - Level 1 (L1, local IDs {1,2,3,4} for {4,5,6,7}): 1→2, 2→1, 2→3, 3→4
- Parents-of-next-seed lists
  - For L0: parents of s1=4 within L0 → {2}
  - For L1: no next seed (end of pass)
- Reconciliation entries
  - (global_id=1, (level_0, local=1), (level_1, local=2)) for edge 5→1
- Remaining cross-level edges (global IDs)
  - None (the tree-merge allowed 5→6 to be encoded locally in L1; 3→6 is not captured within k=2 from s0 and remains a cross-level relation only if never covered; here it becomes intra-level via L1 merge path 5→6)

Wire-format sketch (STOP-delimited, abstracted):
```
FORWARD PASS
  LEVEL 0
    SEED: 1
    PARENTS_OF_NEXT_SEED: [2]
    STOP
    ADJ (local):
      1: [2]   STOP
      2: [3]   STOP
      3: [1]   STOP
    STOP
    STOP

  LEVEL 1
    SEED: 4
    PARENTS_OF_NEXT_SEED: []
    STOP
    ADJ (local):
      1: [2]       # 4→5
      2: [1,3]     # 5→4, 5→6
      3: [4]       # 6→7
      4: []
      STOP per row
    STOP
    STOP

BACKWARD PASS
  (no additional intra-level edges in this example)

RECONCILIATION
  (1, (0,1), (1,2))
  STOP

REMAINING CROSS-LEVEL EDGES
  (none)
  STOP

STOP STOP
```
