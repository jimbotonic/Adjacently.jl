# MGS File Header

Header has 12 bytes:
- version: 'MGS' 3 bytes string
- major + minor version: 2 bytes
- flags: 2 bytes
  * Byte 1 (2 bits + 2 bits + 4 bits):
	- graph type (2 bits): 0b00: directed graph | 0b01: undirected graph
	- coding scheme (2 bits): 0b00: children | 0b01: index
	- integer encoding (4 bits): 0x1: Elias gamma | 0x2: Elias delta | 0x3: Golomb | 0x4: FED | 0x5: Zeta | 0x6: Fibonacci
  * Byte 2: option_flags — interpretation depends on version (see below)
- # vertices: 5 bytes (little-endian UInt40)

The header format is thus:
`<'MGS' string 3 bytes> + <16 bits major|minor version> + <flags 2 bytes> + <# vertices 5 bytes>`

---

## Version 3.0 (legacy) — Byte 2 Interpretation

| Value | Mode |
|-------|------|
| 0x00 | No option (uncompressed) |
| 0x0F | ASTRA (legacy adaptive streaming) |
| 0x10–0x8F | STD (Standard greedy compressed mode) |
| 0x90–0x9F | CS (Command Stream compressed mode) |
| 0xA0–0xAF | CGE (Clustered Graph Encoding) |
| 0xFF | Huffman (deprecated) |

In v3.0, byte 2 identifies the algorithm but **not** its encoding parameters.
The decoder must receive matching params out-of-band.

---

## Version 3.1 — Self-Describing Byte 2

In v3.1 (minor version = 0x01), byte 2 fully encodes the algorithm **and** all
variable parameters. The decoder reconstructs all needed params from byte 2 alone.

### Algorithm IDs (0x00–0x0F): no additional params needed

| Value | Algorithm | Params |
|-------|-----------|--------|
| 0x00 | Legacy MGS (uncompressed) | None |
| 0x01 | Huffman (deprecated) | None |
| 0x02 | STD with recommended defaults | window=64, all features on |
| 0x03 | CS with recommended defaults | compact_copy=true, tight_intervals=true, window=64 |
| 0x04 | CGE with recommended defaults | implicit_ranges, perm, window=16, no intervals |
| 0x05 | ASTRA (legacy adaptive) | None |
| 0x06–0x0F | Reserved | — |

### Parameter Ranges (0x10–0xFF): algorithm + explicit params

| Range | Algorithm | Slots | Param bits |
|-------|-----------|-------|------------|
| 0x10–0x4F | STD + params | 64 | 6 |
| 0x50–0x6F | CS + params | 32 | 5 |
| 0x70–0xFF | CGE + params | 144 | ~7 |

#### STD Parameter Encoding (6 bits, offset = byte2 − 0x10)

```
Bit 5-4: ref_window_size  (2 bits → 0=8, 1=16, 2=32, 3=64)
Bit 3:   copy_blocks      (1=on; adaptive_copy derived = same value)
Bit 2:   stop_deltas      (1=on; empty_prefix derived = same value)
Bit 1:   compact_copy + tight_intervals  (1=both on)
Bit 0:   vlc2             (1=on)
```

Derived (not stored): `adaptive_copy = copy_blocks`, `empty_prefix = stop_deltas`, `fixwidth_ref = false`.

#### CS Parameter Encoding (5 bits, offset = byte2 − 0x50)

```
Bit 4-2: ref_window_size  (3 bits → 0=4, 1=8, 2=16, 3=32, 4=64, 5=128, 6=256, 7=512)
Bit 1:   compact_copy     (1=on)
Bit 0:   tight_intervals  (1=on)
```

#### CGE Parameter Encoding — v3.2 (mixed-radix, offset = byte2 − 0x70, max 143)

```
offset = membership × 72 + window_idx × 12 + interval_mode × 4 + mil_idx
```

| Field | Values | Radix |
|-------|--------|-------|
| membership | 0=:implicit_ranges, 1=:stop | 2 |
| window_idx | 0→8, 1→16, 2→32, 3→64, 4→128, 5→256 | 6 |
| interval_mode | 0=none, 1=intervals, 2=intervals+lr_split | 3 |
| mil_idx | 0→2, 1→3, 2→4, 3→5 | 4 |

`inter_strategy` is not encoded (hardcoded to `:lists`).
Both `intra_mil` and `intra_adapt_mil` are set to the decoded mil value.
All other CGE fields are hardcoded to best-known defaults (see `decode_cge_params` in `src/mgs.jl`).

#### CGE Parameter Encoding — v3.1 (legacy, mixed-radix)

```
offset = membership × 36 + inter_strategy × 18 + window_idx × 3 + interval_mode
```

| Field | Values | Radix |
|-------|--------|-------|
| membership | 0=:stop, 1=:implicit_ranges, 2=:elias_fano, 3=:delta | 4 |
| inter_strategy | 0=:perm, 1=:lists | 2 |
| window_idx | 0→8, 1→16, 2→32, 3→64, 4→128, 5→256 | 6 |
| interval_mode | 0=none, 1=intervals, 2=intervals+lr_split | 3 |

mil not encoded: `intra_mil=4`, `intra_adapt_mil=2` (hardcoded).

### Index Mode (coding_scheme = :index)

When `coding_scheme = :index` (byte 1 bits 5-4 = 0b01), per-vertex offset tables
enable random access to individual vertex data without sequential decoding.

#### STD/CS Index Mode

A per-vertex offset table is written after the 4-bit VLC header:

```
[12-byte header] [4-bit VLC header] [per-vertex offset table] [per-vertex compressed data]
```

**Per-vertex offset table format:**
1. `entry_width` (6 bits): number of bits per offset entry
2. `N+1` entries, each `entry_width` bits wide:
   - entry[0..N-1]: bit offset to vertex i's data (relative to start of per-vertex section)
   - entry[N]: total bit size of per-vertex section (end marker)

Empty vertex detection: `offset[v] == offset[v+1]` → vertex has no neighbors.

In index mode, STD forces `stop_deltas=false` and `empty_prefix=false` (redundant with
offset table). CS replaces STOP-terminated delta with count-prefixed delta encoding:
`varint(count+1)` followed by delta-encoded values.

#### CGE Index Mode

A cluster offset table with `2K+1` entries is written between the header and the
compressed data, enabling O(1) random access to both intra-cluster and inter-cluster
sections:

```
[12-byte header] [cluster offset table] [compressed data]
```

**Cluster offset table format:**
1. `entry_width` (6 bits): number of bits per offset entry
2. `K` (32 bits): number of clusters
3. `2K+1` entries, each `entry_width` bits wide:
   - entries[1..K]: bit offset to cluster i's intra-edge section
   - entry[K+1]: bit offset to inter-cluster edge section start
   - entries[K+2..2K+1]: per-source-cluster inter offset (relative to inter section start)

**Intra-vertex offset table** (within each cluster, after cluster-wide headers):
1. `vtx_entry_width` (6 bits): number of bits per entry
2. `s+1` entries, each `vtx_entry_width` bits wide:
   - entry[0..s-1]: bit offset to local vertex i's payload (relative to payload start)
   - entry[s]: total payload bit size

This enables O(1) random access to any cluster and any vertex within a cluster.

---

## Backward Compatibility

The loader checks the version (bytes 3-4):
- **Version 3.0**: byte 2 interpreted with legacy range logic; caller must pass params
- **Version 3.1**: byte 2 self-describing; CGE uses v3.1 mixed-radix (`decode_cge_params_v31`)
- **Version 3.2**: byte 2 self-describing; CGE uses v3.2 mixed-radix with mil (`decode_cge_params`)

CGE writer now produces v3.2. STD/CS writers still produce v3.1.
Old `.mgz` files (v3.0, v3.1) continue to load unchanged.
