# Graph-compression benchmark drivers

Reproduction drivers for the bits-per-edge (BPE) tables in the
*Community-Aware Vertex Ordering for Reference-Based Graph Compression* paper.

These are **not** part of the test suite (`test/` holds correctness/roundtrip
checks). They regenerate paper tables and can take a long time on the larger
graphs. Run each from the repo root with the project activated:

```bash
~/.juliaup/bin/julia --project=. bench/graph_compression/<driver>.jl
```

## Table → driver map

| Paper table | Driver | Datasets | Status |
|---|---|---|---|
| `tab:ord-ablation` (orderings × encoders, Fibonacci) | `ord_ablation.jl` | web-google, amazon-0601, arxiv-hep-ph, EAT | **done** (4 SNAP; verified on EAT, 8/9 cells) |
| `tab:transfer` (Leiden+LLP-vs-LLP gain, 3 seeds, 2 backends) | `transfer.jl` | web-google, amazon-0601, arxiv-hep-ph, EAT | **done** (verified on EAT, exact) |
| `tab:ablation` (encoder feature history) | `feature_ablation.jl` | cnr-2000 | not started — needs reconstruction (uses removed encoder params; cnr-2000 ships only as `.mgz`, not the `.csv` the old script used) |
| cnr-2000 rows of `tab:ord-ablation` | — | cnr-2000 | not started — K=2 (Orig) vs K=1 (Leiden+LLP) split, medium-confidence params |
| `tab:encoder_comparison` (whole-graph, context-range) | `sota_wholegraph.jl` | all 7 | not started — needs LAW fetch |
| `tab:ra-sota` (random access, context-range) | `sota_ra.jl` | all 7 | not started — needs LAW fetch |

## Datasets

This first subset runs only on the **five datasets committed to the repo**
under `datasets/`:

- `datasets/webgraph/cnr-2000/cnr-2000.mgz`
- `datasets/Web_Google/web-Google.mgz`
- `datasets/Amazon_0601/Amazon0601.mgz`
- `datasets/Arxiv_HEP-PH/Cit-HepPh.mgz`
- `datasets/EAT/EATnew.net`

The two largest LAW graphs (`in-2004`, `enwiki-2013`) are too big to commit and
will be added via a `fetch_datasets.jl` downloader in a later pass, together
with the context-range whole-graph / random-access drivers.

## Reproducibility

BPE is a property of the encoded bitstream and is machine-independent. The
Leiden+LLP ordering is stochastic only through LLP's random visit order; drivers
fix the seed(s) so runs are deterministic.

`ord_ablation.jl` holds the encoder configuration fixed across the three
orderings. It reproduces every LLP / Leiden+LLP cell and the Orig. BG/CS cells to
three decimals. The Orig. CG cell in the paper table came from the per-dataset
best-config CG run (tight-deltas + gamma gap coding), so this fixed-config driver
reads it ~0.08 bpe higher; the ordering-transfer result is unaffected.
