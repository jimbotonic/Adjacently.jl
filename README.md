# Adjacently.jl

Adjacently.jl is a Julia library for the analysis and compression of large complex directed networks.

## MGS Graph Compression

Adjacently provides three compression methods for large-scale directed graphs, all stored in the **MGS** (Memory-Optimized Graph Store) v3.2 binary format. Each method writes a `.mgz` file with a self-describing 12-byte header. A single loader (`load_compressed_mgs3_graph`) reads the header and dispatches to the correct decoder — no parameters need to be passed at load time.

The three methods share a common foundation inspired by the [WebGraph framework](http://webgraph.di.unimi.it/): sorted adjacency lists, delta encoding, reference copy from a sliding window, interval encoding for runs of consecutive neighbors, and variable-length integer codes (Fibonacci, Elias-Gamma/Delta, Golomb, Zeta, FED).

### BG (Baseline Greedy)

Per-vertex greedy cost-based encoding. For each vertex the encoder evaluates every candidate reference in a sliding window and picks the combination of reference copy, intervals, and delta residuals with the lowest bit cost.

Key features:
- Greedy cost-based reference selection (minimum encoding cost, not just maximum overlap)
- Adaptive bitmap encoding: per-reference choice between copy-blocks, raw bitmap, and complement
- Merged VLC vertex header with 28 action codes including empty vertices
- STOP-terminated delta lists (no per-vertex count prefix)
- Multi-reference support (two-ref encoding for additional compression)
- Left/right residual split for interval mode

```julia
write_bg_mgs3_graph(g, "output"; ref_window_size=64, copy_blocks=true,
    stop_deltas=true, lr_split=true, multi_ref=true)
```

### CS (Command Stream)

Unified command-stream encoding with frequency-optimized prefix codes (1–9 bits). Uses the same building blocks as BG but with a simpler (non-greedy) strategy. Internally forces stop-terminated deltas and adaptive copy-blocks.

```julia
write_cs_mgs3_graph(g, "output"; ref_window_size=64, lr_split=false)
```

### CG (Clustered Greedy)

Two-level coarsening approach that exploits community structure. The graph is partitioned into clusters (e.g. via Leiden); intra-cluster edges use reference encoding with adaptive copy-blocks, while inter-cluster edges use list-based encoding.

```julia
write_cg_mgs3_graph(g, "output", clusters; params=CGParams(
    intra_ref_window=64, membership=:implicit_ranges, ...))
```

### Performance

On the CNR-2000 web graph (325,557 vertices, 3,216,152 edges):

| Method | Ordering | Bits per edge | sec/edge |
|--------|----------|--------------|----------|
| **BG** (w=64, multi-ref) | Leiden+LLP | **2.3258** | 6.30e-06 |
| CG K=2 (w=64) | Original | 2.3287 | 1.558e-05 |
| CS (w=64) | Leiden+LLP | 2.3643 | 3.97e-06 |
| WebGraph BV (w=7) | Original | 2.898 | 2.198e-07 |
| WebGraph BV-HC | Host-compressed | 2.448 | — |
| *Leiden+LLP ordering* | — | — | *8.741e-08* |

All methods achieve perfect round-trip fidelity. BG with Leiden+LLP ordering beats WebGraph BV by 20% and BV-HC by 5% on this dataset. BG is 2.5× faster and CS is 3.9× faster than CG K=2, thanks to two-phase sorted merge reference search. The Leiden+LLP ordering cost (load + relabel) is negligible compared to compression.

### MGS File Format (v3.2)

The `.mgz` binary format uses a 12-byte self-describing header:

```
Bytes 0–2:   'MGS' magic (3 bytes)
Bytes 3–4:   version — major (1 byte) + minor (1 byte)
Byte  5:     flags byte 1 — graph_type (2 bits) + coding_scheme (2 bits) + integer_encoding (4 bits)
Byte  6:     flags byte 2 — algorithm ID + encoded parameters (see below)
Bytes 7–11:  vertex count (5 bytes, little-endian UInt40)
```

**Byte 6 (option_flags)** fully identifies the algorithm and all its parameters:

| Range | Algorithm | Description |
|-------|-----------|-------------|
| `0x00` | Legacy MGS | Uncompressed |
| `0x02` | BG defaults | window=64, all features on |
| `0x03` | CS defaults | window=64, compact_copy, tight_intervals |
| `0x04` | CG defaults | implicit_ranges, window=16, no intervals |
| `0x10`–`0x4F` | BG + params | 64 slots: window(2b) + copy_blocks(1b) + stop_deltas(1b) + compact(1b) + vlc2(1b) |
| `0x50`–`0x6F` | CS + params | 32 slots: window(3b) + lr_split(1b) + reserved(1b) |
| `0x70`–`0xFF` | CG + params | 144 slots: membership(2) × window(6) × interval_mode(4) × mil(3) |

### Loading

All `.mgz` files are loaded through a single entry point. The header is fully self-describing — no parameters need to be passed:

```julia
g = load_compressed_mgs3_graph("graph.mgz")
```

### References

- [WebGraph: A Framework for Graph Compression](http://webgraph.di.unimi.it/) - Paolo Boldi and Sebastiano Vigna
- [The WebGraph Framework](https://dl.acm.org/doi/10.1145/988672.988752) - WWW 2004

## Tests

The test suite is organized into individual test files. Run from the project root:

```bash
# Run a specific test set
julia --project test/run_tests_{test-name}.jl

# Run the CNR-2000 compression roundtrip (BG, CS, CG)
julia --project test/run_tests_webgraph_best_compression.jl
```

## Notebooks

Interactive Jupyter notebooks are in `notebooks/`:

| Notebook | Description |
|----------|-------------|
| `cnr-2000-compression.ipynb` | CNR-2000 compression roundtrip with BG, CS, and CG (parameter documentation) |
| `Pagerank.ipynb` | PageRank computation on the Arxiv HEP-PH citation network |
| `shortest_paths.ipynb` | Shortest path algorithms with diffusion-based exploration |
| `shortest_paths2.ipynb` | Extended shortest path analysis |

```julia
julia> using IJulia
julia> notebook()
```

## Development

### Dependencies management

```julia
pkg> activate .
pkg> add {package-name}
pkg> update
```

### Datasets

- [Laboratory for Web Algorithmics - Datasets](https://law.di.unimi.it/datasets.php)
- [WebGraph](https://webgraph.di.unimi.it/)
- [WebGraph - Github](https://github.com/vigna/webgraph)
- [Stanford WebBase Project](http://diglib.stanford.edu:8091/~testbed/doc2/WebBase/)
