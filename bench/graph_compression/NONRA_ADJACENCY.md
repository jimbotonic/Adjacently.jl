# Non-indexed outgoing-adjacency loading

`load_adjacency_mgs3_graph` is an opt-in, full-file preload. It avoids building
a graph and its incoming adjacency, then copying the outgoing adjacency again.
Existing loaders and encoder defaults are unchanged. Streaming remains experimental.

```julia
using Adjacently.MGS: load_adjacency_mgs3_graph
adj = load_adjacency_mgs3_graph("graph.mgz")  # BG, CS or CG, children mode
neighbors = adj[vertex]                     # borrowed decompressed RAM list
owned_neighbors = copy(adj[vertex])
```

Each returned vertex list is independent of the others. Changing `adj[v]` changes
the retained answer for that vertex. Disk bytes and BPE are unchanged.

For CG, the recommended measured configuration enables the previously promoted
decoder flags explicitly:

```julia
adj = load_adjacency_mgs3_graph("cg.mgz";
    compact_lists=true, fused_lr=true, reuse_identity=true,
    sorted_fastpath=true, reuse_scratch=true, compact_output=false)
```

Set `compact_output=true` to release spare list capacity after decoding. This
adds copies and allocation; it is a retained-memory tradeoff, not a latency win.
It is independent of CG's `compact_lists` local-ID option. All flags default off.
The recommended configurations, worker counts and unchanged EAT/CNR file hashes
are pinned in [nonra_adjacency_v1.json](configs/nonra_adjacency_v1.json).

To reproduce the loader comparison on a pinned file, use a clean matching
checkout (the driver checks both file and source hashes):

```bash
taskset -c 2 julia --project=. --threads=1 bench/graph_compression/nonra_decode.jl cnr_cg path/to/cnr_cg.mgz
```

It emits seven randomized paired rounds, allocation counts and retained-output
estimates as JSON. JIT and page cache are warm; full decoding is included in
each operation. Keep other benchmark jobs off the machine during timing.

The separately tuned CNR profile is [cg_cnr_nonra_v2.json](configs/cg_cnr_nonra_v2.json):
K=1, window 128, LR splitting, minimum interval length 5, FULL, children/context-range.
Use the key `cnr_cg_best` with that profile's file. Its 1.91736 BPE is the best of
25 tested configurations, still 1.35% larger than the measured non-RA Zuckerli
file. It does not replace the six fixed-file loader-ablation profiles above.

Supported containers are directed children-mode Fibonacci v3.2, BG/CS
context-range v3.3 and CG context-range v3.2. Both default algorithm IDs and
parameterized headers are supported. CG supports multiple clusters; explicit
membership is needed to preserve arbitrary original IDs. `implicit_ranges`
assumes contiguous clusters and does not preserve a pre-ordering permutation.
Unsupported layouts are rejected, not silently routed to another API. The
existing decoders assume trusted payloads; header/length checks do not make this
a hardened file validator.

For comparisons, report preload time and retained RAM separately from subsequent
in-RAM queries. Do not compare these warm queries against BVGraph or Zuckerli
decoding a compressed list. Non-RA cold queries still require sequential work;
this API does not add restart points, an index, or compressed random access.

The precursor experiment is archived at `research/graph_compression/nonra_v1/`;
its results and frozen source are not modified by this promotion. CNR parameter
exploration is a separate experiment, not a replacement for those fixed-file
timings. Tests exercise exact lists, ownership, unchanged bytes, default-off
switches, multi-cluster membership, and rejected containers.
