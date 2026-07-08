# Context-Range Backend — Format Specification

The `:context_range` integer encoding (header `integer_encoding = 0x7`) is a
context-adaptive range-coding backend shared by the **BG**, **CS**, and **CG**
encoders. It replaces the per-value prefix/varint codes (Fibonacci, Elias, Zeta)
of the classic backends with three adaptive range coders, and adds a random-access
mode built on independently-decodable chunks.

It is selected per-encoder with `integer_encoding=:context_range` and is fully
self-describing: `load_compressed_mgs3_graph(path)` needs no parameters.

- Code: `src/range_coder.jl` (coders), `src/compression.jl` (BG/CS wiring),
  `src/compression/cge.jl` (CG wiring), `src/mgs.jl` (file layout + loaders).
- Per-encoder specifics: [FORMAT_SPECS_BG.md](FORMAT_SPECS_BG.md),
  [FORMAT_SPECS_CS.md](FORMAT_SPECS_CS.md), [FORMAT_SPECS_CG.md](FORMAT_SPECS_CG.md).

---

## 1. Three-stream design

A reference-compressed adjacency record splits into three kinds of integers, each
routed to its own coder so the adaptive statistics never mix:

| Stream | Coder | Carries |
|--------|-------|---------|
| `resid` | `CtxRangeEncoder` (hybrid-uint) | residual neighbor ids (the non-copied part of each list) |
| `refdist` | `CtxRangeEncoder` (hybrid-uint) | reference distances (window back-pointers) |
| `copy` | `BinRangeEncoder` (binary) | the copy bitmap bits (which reference neighbors are reused) |

Structural bits that are **not** value integers — vertex/VLC headers, adaptive-stop
flags, interval-mode continuation bits, split-residual tags — stay in the classic
bitstream (`struct`), written with Fibonacci as before. Only value integers move to
the range coders. (This is why a naïve "range-code the whole residual region"
estimate over-counts: ~0.2–0.3 BPE of that region is structural flag bits that never
route through a coder.)

## 2. File layout (multi-section)

Whole-graph `:context_range` files use a four-section body after the 12-byte header:

```
[12-byte header]
[8-byte LE length: struct] [8-byte LE length: refdist] [8-byte LE length: copy]
[struct bytes] [refdist bytes] [copy bytes] [resid bytes]
```

The three explicit lengths let the loader slice the four blobs; `resid` runs to EOF.
`struct` is the classic bitstream (headers, offset tables, flags); the other three
are finalized range-coder byte streams. Identical layout for BG, CS, and CG.

## 3. Hybrid-uint token

Both `CtxRangeEncoder` streams code integers as **hybrid-uint tokens** (JPEG-XL /
Zuckerli style), config `K=5, M=0, L=1` (`_RC_K/_RC_M/_RC_L` in `range_coder.jl`):

- Values `< 2^K` (= 32) are their own token (direct).
- Larger values emit a token encoding `(exponent, top-M bits, bottom-L bits)`; the
  remaining `exponent − M − L` mantissa bits are written raw (bypass, uniform).
- Alphabet size `NSYM = 2^K + (64−K)·2^(M+L) = 150`.

The token is range-coded against the adaptive context model; only the token carries
learnable structure, the raw mantissa bits are incompressible by design.

## 4. Context model

Each hybrid token is coded under a context = **position bucket × previous-token
magnitude bucket**, reset at every vertex (region):

- position buckets `_RC_NPOS = 4`: list position `0, 1, 2, 3+`
- previous-magnitude buckets `_RC_NPREV = 25`: bucket of the previous token's
  magnitude `0..23`, plus a reset sentinel `24` at each region start
- total `_RC_NCTX = 100` contexts per stream

Models are **global across the graph** — only the position/prev-token *state* resets
per vertex (`rc_reset_region!`), not the learned frequencies. The binary copy coder
uses a 1-bit context (previous copy bit).

## 5. Copy-aware rank gaps

When a vertex references a predecessor, its residual ids are recoded as **ranks in
the complement of the copied set** before entering the `resid` coder, and inverted on
decode:

```
C = current_neighbors \ residuals            # the copied ids (sorted, unique)
rank_C(r) = r − #{c ∈ C : c < r}             # forward: strip copied slots
```

This removes the "holes" the copied neighbors would otherwise leave in the residual
gap sequence, tightening the deltas. The vertex id itself is rank-mapped into the
same space so the coder's vertex-relative first value keeps its locality. Identity
transform when nothing is copied (so `:none`-mode vertices are unaffected). No header
bit — it is implied by `:context_range` + reference mode. For multi-reference (BG),
`C` is the union of both copied sets. Helpers: `_rankgap_forward` / `_rankgap_inverse`
/ `_rank_of`.

---

## 6. Random access (sampled-index + chunked streams)

`:context_range` supports true seekable random access via **chunking**, available for
**BG** and **CS** (not CG — see §7). A full (per-vertex) index is *not* seekable under
range coding because the coder state is continuous; random access therefore **requires
a sampled index** (`index_sample_k > 0`, a multiple of 4). Requesting `coding_scheme=
:index` with `:context_range` and `index_sample_k ≤ 0` is an error.

### Chunking

Every `index_sample_k` vertices, the three coders are finalized to bytes and reset to
fresh adaptive models (`rc_finish_and_reset!` / `brc_finish_and_reset!`), producing a
sequence of independently-decodable chunks aligned to the sample points. The
concatenated chunk bytes are the `resid` / `refdist` / `copy` blobs.

### On-disk additions (in the `struct` stream)

Alongside the classic sampled per-vertex bit-offset table, the `struct` stream stores,
for `:context_range` only:

```
[1 bit: sampled = 1]
[32 bits: index_sample_k]           # raw k (not the legacy (k/4−1) 8-bit field)
[6 bits: entry_width]
[ (⌈n/k⌉+1) × entry_width ]         # sampled per-vertex bit offsets (jump to chunk's first vertex)
[6 bits: coff_width]
[ 3 × (⌈n/k⌉+1) × coff_width ]      # per-chunk cumulative BYTE offsets: resid, refdist, copy
```

The three per-chunk byte-offset tables let a reader slice out the byte range of any
chunk. `k` is stored in 32 bits (vs the classic 8-bit `(k/4−1)`, capped at 1024) so
large chunk sizes are expressible.

### Seek + cross-chunk references

To read vertex `v`, the loader finds its chunk `c = ⌊(v−1)/k⌋`, seeks the `struct`
reader to `data_start_bit + sampled_offset[c]`, spins up fresh range decoders on chunk
`c`'s byte slices, and decodes the whole chunk. A reference may point to a vertex in an
**earlier** chunk not yet materialized; a `ref_resolver` callback recursively
materializes that chunk (chunks are memoised). The reference window is seeded with the
`k`-vertex predecessors of the chunk's first vertex. This is **O(k)** per cold seek (the
chunk plus any referenced earlier chunks), not O(1); references cascade backward, so a
cold seek to chunk `c` may touch up to `c` chunks (memoised thereafter).

Loader: `MGS.load_indexed_mgs3_graph(path)` returns `(n, m, neighbors, algorithm,
stats, reset_fn)`; `stats.chunks_decoded` reports true seek cost. Parser:
`parse_bg_ctxrange_index` (reused verbatim by CS — the `struct`-header layout is
identical). Overhead is small: eat `k=1024` ≈ +0.07–0.08 BPE (~1%).

---

## 7. Encoder support matrix

| Encoder | Whole-graph `:context_range` | Random access (`:index` + sampled) |
|---------|:---:|:---:|
| **BG** | ✅ 3-stream, multi-ref, rank gaps | ✅ chunked (`_load_bg_ctxrange_random_access`) |
| **CS** | ✅ 3-stream, rank gaps | ✅ chunked (`_load_cs_ctxrange_random_access`) |
| **CG** | ✅ children mode (CG-1/2/3), rank gaps | ❌ not implemented — the ctx path falls back to Fibonacci when `cluster_offsets` is set; index-mode ctx is future work |

BG and CS share the same `struct`-header layout, so `parse_bg_ctxrange_index` and the
chunked reader logic are common. CG routes residual/copy/refdist through the same three
coders in children mode via a function barrier (`encode_level` → `@noinline
_encode_level_body!`) that keeps the union-typed sink locals out of the large body
(avoids a Julia-1.12 codegen miscompilation).

---

## 8. Notes / invariants

- `_RC_NSYM` is sized for full 64-bit values; zigzag intermediates in CG can reach the
  high buckets, so the alphabet must cover them (an undersized table segfaults).
- `_RESID_SINK` / `_RESID_SOURCE` are process globals gating the residual redirect;
  they are reset at each encode/decode entry so a prior run cannot leak into the next
  graph's Fibonacci path.
- Lossless and byte-exact: the range coders reproduce the prototype `ctx_encode_typed`
  output bit-for-bit; whole-graph and random-access decodes are verified against the
  original adjacency.
