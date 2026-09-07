# Opt-in codec acceleration

All acceleration flags default to **false**. Worker-count arguments do not enable
parallelism by themselves. The default `cluster_seed=nothing` retains the shared
RNG policy. No existing encoder, ordering, or published size table is rewritten.

The versioned settings are in `configs/codec_acceleration_v1.json`; verify the
matching source manifest from the repository root before a timing run:

Use a clean checkout of the promotion commit. Uncommitted reader/cache experiments
intentionally fail the source guard; do not regenerate the manifest to hide drift.

```sh
sha256sum -c bench/graph_compression/configs/codec_acceleration_v1.sha256
JULIA_NUM_THREADS=4 julia --project -e 'using Pkg; Pkg.test()'
```

Tests require only the public repository, including the tracked EAT input and
`test/fixtures/codec_legacy_v1.json`. They never silently skip exact-output checks
because the companion research checkout is absent.

## Encoding

Start Julia with `--threads=8`, then explicitly pass:

```julia
write_bg_mgs3_graph(g, path; parallel_search=true, search_workers=8,
    search_batch_size=4096, cost_model=0)
write_cs_mgs3_graph(g, path; parallel_search=true, search_workers=8,
    search_batch_size=4096, cost_model=0)
write_cg_mgs3_graph(g, path, clusters; params=full_params,
    parallel_search=true, search_workers=8)
```

BG/CS parallel search supports children mode and analytical FULL/FAST cost models;
exact-cost trial writing and index mode are rejected. CG supports children/index,
including K=1, but not MGS/block trial encoders. Workers evaluate read-only costs;
ordered stream emission stays serial. Candidate tie-breaking, bytes, and BPE are
unchanged. Concurrent top-level encoder calls remain unsupported because stream
emission uses shared hooks. Encoder parameters are separately pinned by
`parallel_encoding_v1.json` and `cg_speed_v1.json`; do not infer them from timing.

## CG decoding

Both full loading and indexed querying accept the same default-off options:

```julia
using JSON
cfg = JSON.parsefile("bench/graph_compression/configs/codec_acceleration_v1.json")
opts = (; (Symbol(k) => v for (k, v) in cfg["cg_decode"]["recommended"])...)
g = load_compressed_mgs3_graph(path; opts...)
idx = load_indexed_mgs3_graph(path; opts...)
```

`compact_lists` enables compact local IDs; `fused_lr` enables the compatible
range LR reader. `reuse_identity` avoids the second adjacency representation
when K=1, IDs are identity-mapped, and local/output types match. `sorted_fastpath`
avoids repeat sorts where ordering is guaranteed and uses adjacent deduplication.
`reuse_scratch` reuses additions and merge buffers in the fused reader. Enable
the full recorded option set for the validated fast path. Tight-gap and legacy
formats retain compatible fallbacks; returned lists never alias reusable scratch.

For independent indexed context-range clusters, additionally use
`parallel_decode=true, decode_workers=8` on the full loader. K=1 falls back to
serial decode. Options are forwarded into every cluster worker. Graph construction
is serial; indexed-reader caches are not thread-safe. Each cold K=1 query still
decodes the entire cluster. Report mean including misses, not cache-hit latency.

The promoted decoder requires packed-reader and bounded model-reset primitives.
Unrelated CDF-search and RA-cache/continuation experiments are not included in
this promotion. Historical timing artifacts retain their original source pins;
their baseline enabled compact/fused reading, while this public API defaults
those options off. Retiming this exact source/configuration is required before
assigning a new headline speedup to the promotion commit.

## Per-cluster LLP: changes ordering

This remains explicitly opt-in:

```julia
Random.seed!(53) # preceding Leiden partition
relabel_graph_leiden_llp(g; parallel_clusters=true, cluster_seed=53,
    cluster_workers=8, llp_mode=:sym, llp_passes=3)
```

The historical EAT scaling run used 4 workers; CNR used 8. The new per-cluster
seed policy changes RNG consumption and can change ordering and compressed BPE
relative to historical shared-RNG runs. **Re-run and re-record every affected
table.** Only seeded serial versus identically seeded parallel runs are promised
identical. Do not mix this policy into historical tables without recording it.
