# Adjacently.jl

Adjacently.jl is a Julia library for the analysis and compression of large complex directed networks.

## MGS Graph Compression

Adjacently provides three compression methods for large-scale directed graphs, all stored in the **MGS** (Modified Graph Structure) v3 binary format. Each method writes a `.mgz` file with a 12-byte header that encodes graph type, coding scheme, integer encoding, and an option flag that identifies the compression approach. A single loader (`load_compressed_mgs3_graph`) reads the header and dispatches to the correct decoder.

The three methods share a common foundation inspired by the [WebGraph framework](http://webgraph.di.unimi.it/): sorted adjacency lists, delta encoding, reference copy from a sliding window, interval encoding for runs of consecutive neighbors, and variable-length integer codes (Fibonacci, Elias-Gamma/Delta, Golomb, Zeta, FED).

### STD (Standard Greedy)

Per-vertex greedy cost-based encoding. For each vertex the encoder evaluates every candidate reference in a sliding window and picks the combination of reference copy, intervals, and delta residuals with the lowest bit cost.

Key features:
- Greedy cost-based reference selection (minimum encoding cost, not just maximum overlap)
- Adaptive bitmap encoding: per-reference choice between copy-blocks, raw bitmap, and complement
- VLC v2 vertex header tags optimized for the most common encoding patterns
- STOP-terminated delta lists (no per-vertex count prefix)
- 1-bit empty-vertex prefix to skip isolated vertices

```julia
write_std_mgs3_graph(g, "output"; ref_window_size=64, copy_blocks=true,
    adaptive_copy=true, stop_deltas=true, empty_prefix=true,
    compact_copy=true, tight_intervals=true, vlc2=true)
```

### CS (Command Stream)

Unified command-stream encoding. The encoder writes a sequence of commands per vertex using the same building blocks as STD but with a simpler (non-greedy) strategy. Internally forces stop-terminated deltas and adaptive copy-blocks.

```julia
write_cs_mgs3_graph(g, "output"; ref_window_size=64,
    compact_copy=true, tight_intervals=true)
```

### CGE (Clustered Graph Encoding)

Two-level coarsening approach that exploits community structure. The graph is partitioned into clusters (e.g. via Leiden); intra-cluster edges use reference encoding with adaptive copy-blocks, while inter-cluster edges use permutation-based encoding of the bipartite structure.

```julia
write_cge_mgs3_graph(g, "output", clusters; params=CGEParams(L=128, ...))
```

### Performance

On the CNR-2000 web graph (325,557 vertices, 3,216,152 edges):

| Method | Bits per edge | File size |
|--------|--------------|-----------|
| STD    | 2.88         | 1,157,508 bytes |
| CS     | 2.88         | 1,157,224 bytes |
| CGE   | 2.43         | 978,562 bytes   |

All methods achieve perfect round-trip fidelity.

### MGS File Format

The `.mgz` binary format uses a 12-byte header:

```
'MGS' (3 bytes) + version (2 bytes) + flags (2 bytes) + vertex count (5 bytes)
```

The option flags byte identifies the compression method:
- `0x00`: uncompressed
- `0x0F`: ASTRA (legacy)
- `0x10`-`0x8F`: STD
- `0x90`-`0x9F`: CS
- `0xA0`-`0xAF`: CGE
- `0xFF`: Huffman (deprecated)

### Loading

All `.mgz` files are loaded through a single entry point:

```julia
# STD / CS — pass the same encoding flags used at write time
g = load_compressed_mgs3_graph("graph.mgz";
    copy_blocks=true, adaptive_copy=true, ref_window_size=64,
    stop_deltas=true, empty_prefix=true, compact_copy=true,
    tight_intervals=true, vlc2=true)

# CGE — pass the same CGEParams used at write time
g = load_compressed_mgs3_graph("graph.mgz"; cge_params=params)
```

### References

- [WebGraph: A Framework for Graph Compression](http://webgraph.di.unimi.it/) - Paolo Boldi and Sebastiano Vigna
- [The WebGraph Framework](https://dl.acm.org/doi/10.1145/988672.988752) - WWW 2004

## Tests

The test suite is organized into individual test files. Run from the project root:

```bash
# Run a specific test set
julia --project test/run_tests_{test-name}.jl

# Run the CNR-2000 compression roundtrip (STD, CS, CGE)
julia --project test/cnr_2000_best_known_compression.jl
```

## Notebooks

Interactive Jupyter notebooks are in `notebooks/`:

| Notebook | Description |
|----------|-------------|
| `cnr-2000-compression.ipynb` | CNR-2000 compression roundtrip with STD, CS, and CGE (parameter documentation) |
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
