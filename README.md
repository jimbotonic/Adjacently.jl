# Adjacently.jl

Adjacently.jl is a Julia library for the analysis and compression of large complex
directed networks. It started as a graph-compression toolkit (the **MGS** format and
its BG/CS/CG encoders) and now also hosts the research modules built on top of it:
neural diffusion fingerprints, metabolic pathway discovery, GNN-based vertex ordering,
and an agent-based governance simulation.

| Module | Purpose |
|--------|---------|
| `Compression`, `MGS` | BG / CS / CG encoders and the self-describing `.mgz` container |
| `Clustering`, `Relabeling` | Leiden / LLP partitioning and vertex-ordering pipelines |
| `Graph`, `Algo`, `Paths`, `PageRank`, `RandomWalks`, `Metrics` | Core graph analysis |
| `GNN` | 2-layer GNN / GAT scoring for compression-aware vertex ordering |
| `Fingerprints` | Neural diffusion fingerprints (NDF) and graph-based text models |
| `Metabolic` | BiGG model loading, FBA, and chemistry-aware pathway search |
| `MycelialPolis` | Agent-based simulation of latent societies in reactive hosts |

## Installation

Adjacently.jl is not in the general registry — clone it and instantiate the project
environment. Julia **1.12** or newer is required.

```bash
git clone https://github.com/jimbotonic/Adjacently.jl.git
cd Adjacently.jl
julia --project=. -e 'using Pkg; Pkg.instantiate()'
```

```julia
using Adjacently
using Adjacently.MGS
using Adjacently.Compression
```

Notes:

- The environment pulls Flux and CUDA (used by `Fingerprints` training). Everything
  outside the neural modules runs CPU-only; CUDA artifacts are downloaded at
  instantiate time but no GPU is required to use the compression stack.
- Run scripts and tests with `--project=.` from the repository root so relative
  dataset paths resolve.

## Datasets

Graphs used by the test suite and benchmark drivers ship **pre-compressed** so the
repository stays small. Raw edge lists and WebGraph binaries are gitignored.

| File | Vertices | Edges | Source |
|------|----------|-------|--------|
| `datasets/webgraph/cnr-2000/cnr-2000.mgz` | 325,557 | 3,216,152 | LAW (pre-reordered, not crawl order) |
| `datasets/Web_Google/web-Google.mgz` | 875,713 | 5,105,039 | SNAP web-Google |
| `datasets/Amazon_0601/Amazon0601.mgz` | 403,394 | 3,387,388 | SNAP Amazon-0601 |
| `datasets/Arxiv_HEP-PH/Cit-HepPh.mgz` | 34,546 | 421,578 | SNAP Arxiv HEP-PH citations |
| `datasets/EAT/EATnew.net` | 23,219 | 325,593 | Edinburgh Associative Thesaurus (Pajek) |

```julia
using Adjacently.MGS: load_compressed_mgs3_graph
using Adjacently.IO: load_graph_from_pajek

g = load_compressed_mgs3_graph("datasets/Web_Google/web-Google.mgz")
h = load_graph_from_pajek("datasets/EAT/EATnew.net")
```

The two largest LAW graphs (`in-2004`, `enwiki-2013`) and the true crawl-order
`cnr-2000` are too large to commit. Fetch them with:

```bash
bench/graph_compression/fetch_datasets.sh
```

(converting BVGraph → CSV needs the WebGraph jars on `WEBGRAPH_CP`).

External dataset sources:

- [Laboratory for Web Algorithmics — Datasets](https://law.di.unimi.it/datasets.php)
- [WebGraph](https://webgraph.di.unimi.it/) · [WebGraph on GitHub](https://github.com/vigna/webgraph)
- [SNAP](https://snap.stanford.edu/data/)

## Graph Compression

Adjacently provides three compression methods for large-scale directed graphs, all
stored in the **MGS** v3.x binary format. Each method writes a `.mgz` file with a
self-describing 12-byte header. A single loader (`load_compressed_mgs3_graph`) reads
the header and dispatches to the correct decoder — no parameters need to be passed at
load time.

The three methods share a common foundation inspired by the
[WebGraph framework](http://webgraph.di.unimi.it/): sorted adjacency lists, delta
encoding, reference copy from a sliding window, interval encoding for runs of
consecutive neighbors, and variable-length integer codes (Fibonacci, Elias-Gamma/Delta,
Golomb, Zeta, FED).

### BG (Baseline Greedy)

Per-vertex greedy cost-based encoding. For each vertex the encoder evaluates every
candidate reference in a sliding window and picks the combination of reference copy,
intervals, and delta residuals with the lowest bit cost.

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

Unified command-stream encoding with frequency-optimized prefix codes (1–9 bits). Uses
the same building blocks as BG but with a simpler (non-greedy) strategy. Internally
forces stop-terminated deltas and adaptive copy-blocks.

```julia
write_cs_mgs3_graph(g, "output"; ref_window_size=64, lr_split=false)
```

### CG (Clustered Greedy)

Two-level coarsening approach that exploits community structure. The graph is
partitioned into clusters (e.g. via Leiden); intra-cluster edges use reference encoding
with adaptive copy-blocks, while inter-cluster edges use list-based encoding.

```julia
write_cg_mgs3_graph(g, "output", clusters; params=CGParams(
    intra_ref_window=64, membership=:implicit_ranges, ...))
```

### Vertex ordering

Compression quality depends heavily on vertex order. The **Leiden-as-ordering**
pipeline — global LLP → Leiden fine clusters → per-cluster intra-LLP → concatenate,
then encode with K=1 — is the strongest ordering across the datasets above and is what
the benchmark drivers reproduce. `Relabeling` and `Clustering` expose the pieces.

### Integer coding backends

Besides the classic per-value codes, the encoders support `integer_encoding=:context_range`:
a context-adaptive range-coding backend shared by BG, CS, and CG. It splits a record
into separate symbol streams (residuals, reference distances, copy bitmaps, commands,
flags for BG/CS; three streams for CG) so adaptive statistics never mix, and it adds a
random-access mode built on independently decodable chunks.

```julia
write_bg_mgs3_graph(g, "output"; integer_encoding=:context_range)
```

Files written this way carry minor version `0x03` (v3.3 body layout) and remain fully
self-describing. See [docs/FORMAT_SPECS_CONTEXT_RANGE.md](docs/FORMAT_SPECS_CONTEXT_RANGE.md).

### Random access

`load_indexed_mgs3_graph` returns an object with a `neighbors(v)` accessor. With
**BG + `:context_range` + a sampled index** this is true O(k) random access: it seeks to
the chunk holding `v`, decodes only that chunk (memoised), and recursively materialises
referenced vertices from earlier chunks. All other formats fall back to a single full
decode plus O(1) lookups thereafter.

```julia
idx = load_indexed_mgs3_graph("graph.mgz")
nbrs = idx.neighbors(42)
```

### Performance

On the CNR-2000 web graph (325,557 vertices, 3,216,152 edges):

| Method | Ordering | Bits per edge | sec/edge |
|--------|----------|--------------|----------|
| **CS** (w=256) | Leiden+LLP | **2.3043** | 5.63e-05 |
| BG (w=64, multi-ref) | Leiden+LLP | 2.3259 | 1.95e-05 |
| CG K=2 (w=64) | Original | 2.3286 | 5.63e-06 |
| WebGraph BV (w=7) | Original | 2.898 | 2.198e-07 |
| WebGraph BV-HC | Host-compressed | 2.448 | — |

All methods achieve perfect round-trip fidelity. CS with Leiden+LLP ordering beats
WebGraph BV by 21% and BV-HC by 6% on this dataset. CG K=2 is the fastest Adjacently
encoder at 5.63e-06 sec/edge thanks to analytical cost estimation. A dual cost model is
available: full model (`cost_model=0`) maximizes compression quality, fast model
(`cost_model=1`) trades modest BPE for significant speedup (e.g. CS: 19× faster at
+0.28 BPE).

### References

- [WebGraph: A Framework for Graph Compression](http://webgraph.di.unimi.it/) — Paolo Boldi and Sebastiano Vigna
- [The WebGraph Framework](https://dl.acm.org/doi/10.1145/988672.988752) — WWW 2004

## MGS File Format

The `.mgz` / `.mgs` binary format uses a 12-byte self-describing header:

```
Bytes 0–2:   'MGS' magic (3 bytes)
Bytes 3–4:   version — major (1 byte) + minor (1 byte)
Byte  5:     flags byte 1 — graph_type (2 bits) + coding_scheme (2 bits) + integer_encoding (4 bits)
Byte  6:     flags byte 2 — algorithm ID + encoded parameters (see below)
Bytes 7–11:  vertex count (5 bytes, little-endian UInt40)
```

Minor version `0x02` is the self-describing v3.2 layout; `0x03` marks the v3.3
five-stream body used by BG/CS `:context_range` files. `integer_encoding` selects the
value code (`0x1` Elias gamma, `0x2` Elias delta, `0x3` Golomb, `0x4` FED, `0x5` Zeta,
`0x6` Fibonacci, `0x7` context-range).

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

All `.mgz` files are loaded through a single entry point. The header is fully
self-describing — no parameters need to be passed:

```julia
g = load_compressed_mgs3_graph("graph.mgz")
```

### Format specifications

| Document | Contents |
|----------|----------|
| [docs/MGS_HEADER.md](docs/MGS_HEADER.md) | Header layout and byte-2 encoding across versions |
| [docs/FORMAT_SPECS_MGS.md](docs/FORMAT_SPECS_MGS.md) | Uncompressed `.mgs` container |
| [docs/FORMAT_SPECS_BG.md](docs/FORMAT_SPECS_BG.md) | BG bitstream |
| [docs/FORMAT_SPECS_CS.md](docs/FORMAT_SPECS_CS.md) | CS bitstream |
| [docs/FORMAT_SPECS_CG.md](docs/FORMAT_SPECS_CG.md) | CG bitstream and cluster offset table |
| [docs/FORMAT_SPECS_CONTEXT_RANGE.md](docs/FORMAT_SPECS_CONTEXT_RANGE.md) | Context-range backend and chunked random access |

## Research Modules

These modules are the code side of ongoing research projects. Papers, experiment
scripts, and results live in a separate research repository and are not part of this
package.

### Fingerprints — neural diffusion fingerprints

APPNP-style generalization of the linear *Diffusion Fingerprints* construction
([arXiv:1408.4966](https://arxiv.org/abs/1408.4966)), which is recovered as the linear,
feature-free, K→∞ limit. Includes several propagation variants (fixed-α APPNP, GPR-style
learned coefficients, non-linear), directed PPR adjacency, and seed-gated pooling.

On top of the diffusion core the module carries graph-based text classification:
corpus readers, vocabulary and PPMI / collocation graph construction (symmetric or
directed), per-document graph-of-words models (`PerDocGNN`) with optional GloVe
initialization, hypergraph models (`HyperGAT`), character n-gram tokenization, and
BERT/word feature loading.

### Metabolic — pathway discovery

- BiGG JSON model loading → directed metabolite graph
- Flux Balance Analysis (JuMP + HiGHS) with condition presets
- Chemistry-aware edge weighting (atom similarity from BiGG formulas)
- MetaRoute-style reaction-context ranking, RDT atom-tracing weights

Pathway search itself uses the `Paths` module: bidirectional personalized-PageRank
corridors (`ppr_path_scores`, `wppr`) and fingerprint-guided branching
(`discover_paths_fp`, `discover_paths_bottleneck`).

### MycelialPolis — agent-based governance simulation

Simulation of non-hierarchical latent virtual societies inside reactive host platforms:
multiplex layers over configurable topologies (including scale-free stress tests),
adversarial host strategies (degree, betweenness, legibility, infiltration, attrition,
and an ε-greedy `AdaptiveHost`), adoption dynamics, repression and infrastructure
models, basin estimation, ablations, and distributed constitutional sensing. Order
parameters include Φ, Ψ_T, Λ, Γ, hierarchy scores, and `H_informal_excess`
(degree-preserving-null-corrected informal hierarchy).

### GNN — compression-aware vertex ordering

Lightweight 2-layer GNN and GAT producing per-vertex scores, with manual gradients and
both proxy and REINFORCE training (`train_gnn_proxy!`, `train_gnn_reinforce!`), used to
derive vertex orderings via `relabel_vertices_gnn`.

## Benchmarks

`bench/graph_compression/` holds reproduction drivers for the bits-per-edge tables of
the community-aware vertex-ordering paper. They are **not** part of the test suite and
can take a long time on the larger graphs:

```bash
julia --project=. bench/graph_compression/ord_ablation.jl   # orderings × encoders
julia --project=. bench/graph_compression/transfer.jl       # Leiden+LLP vs LLP gain
```

See [bench/graph_compression/README.md](bench/graph_compression/README.md) for the
table → driver map and for which cells cannot be reproduced from committed data.

## Tests

The test suite is organized into individual test files, each self-contained. Most
include `test/run_tests_main.jl`, which activates the project and loads the shared
imports and dataset helpers. Run from the project root:

```bash
# Any single test set
julia --project test/run_tests_{test-name}.jl

# Compression roundtrips on the committed fixtures
julia --project test/run_tests_compression.jl
julia --project test/run_tests_bg_adaptive.jl
julia --project test/run_tests_cg.jl
julia --project test/run_tests_cg_vs_bg.jl

# Index / random-access modes
julia --project test/run_tests_index_mode.jl

# Research modules
julia --project test/run_tests_mp_smoke.jl          # MycelialPolis (run_tests_mp_*.jl)
julia --project test/run_tests_ndf_linear_limit.jl  # Fingerprints (run_tests_ndf_*.jl)
```

## Development

### Dependencies management

```julia
pkg> activate .
pkg> add {package-name}
pkg> update
```

## License

GNU General Public License v3.0 — see [LICENSE](LICENSE).
