# Adjacently.jl

Adjacently.jl is a Julia library for the analysis and compression of large complex
directed networks. It started as a graph-compression toolkit (the **MGS** format and
its BG/CS/CG encoders) and now also hosts the research modules built on top of it:
neural diffusion fingerprints, metabolic pathway discovery, and an agent-based
governance simulation.

| Module | Purpose |
|--------|---------|
| `Compression`, `MGS` | BG / CS / CG encoders and the self-describing `.mgz` container |
| `Clustering`, `Relabeling` | Leiden / LLP partitioning and vertex-ordering pipelines |
| `Graph`, `Algo`, `Paths`, `PageRank`, `RandomWalks`, `Metrics` | Core graph analysis |
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

### Compression ratio

Bits per edge (BPE = `8 × filesize / m`) is a property of the encoded bitstream and is
therefore machine-independent, so every figure below is reproducible on any machine.
Encoder throughput is not — it depends on the CPU, the JIT, and the I/O path — so no
timings are published here; see the [companion paper](https://arxiv.org/abs/2605.21510)
for the speed discussion.

The tables report the `:context_range` backend at K=1 against the strongest
reference-based baselines: WebGraph **BV-HC** (high-compression) and **Zuckerli**
(max-compression mode), on seven datasets in both their native and Leiden+LLP orderings.
The best of the three Adjacently encoders is in bold; the last column is its relative
reduction over Zuckerli.

**Whole-graph compression (BPE)**

| Dataset | Ordering | BV-HC | Zuckerli | BG | CS | CG | vs Zuck. |
|---|---|---|---|---|---|---|---|
| amazon-0601 | native | 12.853 | 10.254 | 10.155 | 10.380 | **9.903** | +3.4% |
| amazon-0601 | leiden | 8.196 | 6.886 | 6.557 | 6.629 | **6.381** | +7.3% |
| arxiv-hep-ph | native | 10.132 | 9.379 | 8.954 | **8.823** | 8.899 | +5.9% |
| arxiv-hep-ph | leiden | 7.710 | 7.382 | 6.704 | 6.637 | **6.572** | +11.0% |
| cnr-2000 | native | 2.565 | 1.886 | 1.773 | **1.702** | 1.946 | +9.8% |
| cnr-2000 | leiden | 2.714 | 2.063 | 1.937 | **1.870** | 2.121 | +9.4% |
| EAT | native | 10.725 | 9.703 | 9.102 | **9.002** | 9.061 | +7.2% |
| EAT | leiden | 9.725 | 9.148 | 8.462 | 8.546 | **8.407** | +8.1% |
| enwiki-2013 | native | 15.625 | 13.299 | 13.007 | **12.945** | 13.121 | +2.7% |
| enwiki-2013 | leiden | 12.412 | 10.934 | **10.561** | 10.668 | 10.615 | +3.4% |
| in-2004 | native | 1.839 | 1.319 | 1.282 | **1.245** | 1.483 | +5.6% |
| in-2004 | leiden | 1.923 | 1.417 | 1.367 | **1.344** | 1.576 | +5.2% |
| web-google | native | 6.165 | 4.957 | 4.807 | **4.790** | 4.884 | +3.4% |
| web-google | leiden | 4.095 | 3.408 | 3.188 | **3.142** | 3.349 | +7.8% |

**Random access (BPE)**, each system in its own random-access mode: BV at w=7 with
bounded reference chains (m=3), Zuckerli with `--allow_random_access`, ours with the
embedded sampled index.

| Dataset | Ordering | BV w=7 | Zuck.-RA | BG | CS | CG | vs Zuck. |
|---|---|---|---|---|---|---|---|
| amazon-0601 | native | 13.001 | 10.514 | **10.260** | 10.485 | 20.629 | +2.4% |
| amazon-0601 | leiden | 8.343 | 7.055 | **6.695** | 6.800 | 8.451 | +5.1% |
| arxiv-hep-ph | native | 10.263 | 9.130 | 9.149 | **9.102** | 14.545 | +0.3% |
| arxiv-hep-ph | leiden | 7.964 | 7.176 | 6.948 | **6.918** | 8.437 | +3.6% |
| cnr-2000 | native | 3.178 | 2.443 | 1.749 | **1.698** | 5.437 | +30.5% |
| cnr-2000 | leiden | 3.225 | 2.497 | 1.890 | **1.845** | 3.050 | +26.1% |
| EAT | native | 10.754 | 9.322 | 9.130 | **9.081** | 12.608 | +2.6% |
| EAT | leiden | 9.769 | 8.756 | **8.508** | 8.598 | 10.411 | +2.8% |
| enwiki-2013 | native | 16.195 | 13.722 | **13.001** | 13.017 | 25.799 | +5.3% |
| enwiki-2013 | leiden | 12.867 | 11.304 | **10.505** | 10.672 | 13.775 | +7.1% |
| in-2004 | native | 2.388 | 1.932 | 1.278 | **1.251** | 4.185 | +35.2% |
| in-2004 | leiden | 2.352 | 1.883 | 1.354 | **1.337** | 2.354 | +29.0% |
| web-google | native | 6.717 | 5.526 | 4.920 | **4.917** | 11.772 | +11.0% |
| web-google | leiden | 4.690 | 3.900 | 3.328 | **3.299** | 4.272 | +15.4% |

Reading the tables:

- The best of the three encoders beats Zuckerli in all 28 cells, by +0.3% to +35%. The
  margin is largest on the most compressible graphs (in-2004, cnr-2000), where
  entropy-coding the control bits — previously 0.2–0.3 BPE of raw overhead — matters most.
- **BG-RA and CS-RA are the random-access winners.** CG-RA uses the K>1 cluster layer
  directly (cluster = seek chunk), so its cross-cluster edges cannot use reference or
  interval compression and pay the expensive inter-cluster stub encoding.
- BV-HC (w=16, unbounded reference chains) has no practical random access, which is why
  the random-access table uses BV w=7, m=3 instead.
- WebGraph's `.offsets` sidecar is excluded from every BV figure, conservatively against
  us. Zuckerli's random-access file can be *smaller* than its max-compression file on
  small graphs (per-block coder reset amortizes table overhead differently); each system
  is therefore compared in its own mode rather than across modes.
- All compressed files are verified by roundtrip decompression.

The `leiden` rows use the Leiden+LLP ordering with the small-cluster merge refinement
(`relabel_graph_leiden_llp(g; merge_clusters=:auto)`), so these absolute values are not
directly comparable to the ordering-ablation numbers reproduced by
`bench/graph_compression/ord_ablation.jl`, which use the Fibonacci backend and the
unmerged ordering. With that classic backend, the best cnr-2000 results are CS 2.304 and
BG 2.326 (Leiden+LLP) and CG K=2 2.329 (LAW order), against BV 2.898 and BV-HC 2.448.

> **Reproducibility of these two tables.** The shipped drivers reproduce the paper's
> *ordering-transfer* result, not the SOTA grids above. Reproducing these cells needs
> three things this repository does not currently provide end-to-end: the two largest
> LAW graphs (fetch them, see [Datasets](#datasets)), the WebGraph and Zuckerli binaries
> for the baseline columns, and a context-range SOTA driver — `sota_wholegraph.jl` /
> `sota_ra.jl` are not written yet. The committed encoder is the latest one; on the
> dataset swept end-to-end (EAT) it matches or beats the published CS/CG cells, so the
> residual is per-cell configuration, not an encoder regression. See
> [bench/graph_compression/README.md](bench/graph_compression/README.md) for the exact
> per-cell comparison.

A dual cost model is available for the greedy encoders: the full model (`cost_model=0`)
explores all encoding options and candidates for the best ratio, while the fast model
(`cost_model=1`) trades a few tenths of a BPE for a large reduction in encoding work.

### References

- [WebGraph: A Framework for Graph Compression](http://webgraph.di.unimi.it/) — Paolo Boldi and Sebastiano Vigna
- [The WebGraph Framework](https://dl.acm.org/doi/10.1145/988672.988752) — WWW 2004
- [Zuckerli: A New Compressed Representation for Graphs](https://arxiv.org/abs/2009.01353) — Versari et al., IEEE Access 2020
- [Community-Aware Vertex Ordering for Reference-Based Graph Compression: A Cross-Encoder Empirical Study](https://arxiv.org/abs/2605.21510) — the companion paper behind the tables above and the `bench/graph_compression` drivers

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

## Benchmarks and reproduction

`bench/graph_compression/` holds reproduction drivers for the bits-per-edge tables of
the [companion paper](https://arxiv.org/abs/2605.21510). They are **not** part of the
test suite and can take a long time on the larger graphs.

### Step 1 — environment

```bash
julia --project=. -e 'using Pkg; Pkg.instantiate()'
```

BPE is machine-independent, so results do not depend on the box. The ordering pipeline
is stochastic only through LLP's random visit order; both drivers fix the seeds, so runs
are deterministic. Defaults are overridable by environment variable —
`SEEDS=0,1,2`, `COST_MODEL=0`, `VERIFY=true`, and `BACKENDS=fib,ctx` for `transfer.jl`.

### Step 2 — datasets

The four SNAP graphs and a pre-reordered cnr-2000 are committed, so the ordering drivers
run with no downloads. For the LAW graphs (`in-2004`, `enwiki-2013`, and the true
crawl-order `cnr-2000`):

```bash
# downloads + md5-verifies, then converts BVGraph -> CSV
WEBGRAPH_CP='/path/to/webgraph-3.6.12.jar:/path/to/deps/*' \
  bash bench/graph_compression/fetch_datasets.sh
```

Without `WEBGRAPH_CP` the files are still downloaded and verified, and the script prints
the exact `ArcListASCIIGraph` command to finish the conversion by hand. Load a fetched
graph with `load_adjacency_list_from_csv("datasets/webgraph/<name>/<name>.csv", ',', true)`.

### Step 3 — run the drivers

```bash
julia --project=. bench/graph_compression/ord_ablation.jl   # orderings × encoders
julia --project=. bench/graph_compression/transfer.jl       # Leiden+LLP vs LLP gain
```

Both reorder in memory from the committed graphs — no intermediate ordering files — and
verify every encode by roundtrip decompression.

### What reproduces, and what does not

| Result | Status |
|--------|--------|
| Ordering ablation (orderings × encoders, Fibonacci) | reproduces from committed data (verified on EAT: 8/9 cells; see the Orig-CG caveat in the bench README) |
| Ordering transfer (3 seeds, 2 backends) | reproduces from committed data (verified on EAT, exact) |
| cnr-2000 *native* rows | needs the fetched LAW graph — the committed `.mgz` is pre-reordered |
| Context-range SOTA grids (whole-graph and random access) | **no driver yet** — `sota_wholegraph.jl` / `sota_ra.jl` are unwritten |
| Encoder feature ablation | blocked — uses encoder parameters that no longer exist |

Baseline columns need external tools, which this repository does not vendor: WebGraph
3.6.12 (Java 21) for the BV and BV-HC rows, and Zuckerli for the Zuckerli rows.

Note the ordering used: the drivers above call
`relabel_graph_leiden_llp(g; merge_clusters=nothing)`, matching the paper's ablation
tables. The SOTA grids instead use the small-cluster merge (`merge_clusters=:auto`),
which is implemented and public (`src/relabeling.jl`) but not driven by any script here.

See [bench/graph_compression/README.md](bench/graph_compression/README.md) for the full
table → driver map and the per-cell comparison against the published numbers.

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
