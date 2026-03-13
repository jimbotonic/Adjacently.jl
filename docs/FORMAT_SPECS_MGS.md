# MGS Uncompressed Graph Format Specification

## Overview

The MGS (Modified Graph Storage) uncompressed format stores directed or undirected
graphs as raw fixed-width adjacency lists. No variable-length coding, delta encoding,
or reference compression is applied. Vertex IDs are stored as fixed-width big-endian
integers using the minimum number of bytes needed to represent all vertices.

File extension: `.mgs`

## File Structure

```
[Header (12 bytes)] [Graph Data]
```

### Header (12 bytes)

| Offset | Size | Field |
|--------|------|-------|
| 0 | 3 bytes | Signature: `MGS` (0x4D 0x47 0x53) |
| 3 | 1 byte | Major version: `0x03` |
| 4 | 1 byte | Minor version: `0x02` (self-describing v3.2) |
| 5 | 1 byte | Flags byte 1: graph_type (2b) | coding_scheme (2b) | integer_encoding (4b) |
| 6 | 1 byte | Flags byte 2: `0x00` (ALG_LEGACY_MGS) |
| 7 | 5 bytes | Number of vertices N (little-endian UInt40) |

**Flags byte 1 layout** (MSB to LSB):

```
Bits 7-6: graph_type       0b00 = directed, 0b01 = undirected
Bits 5-4: coding_scheme    0b00 = :children, 0b01 = :index
Bits 3-0: integer_encoding 0x2 = Elias delta (default; unused by uncompressed format)
```

**Flags byte 2**: Always `0x00` for uncompressed graphs (ALG_LEGACY_MGS).

**Vertex count**: 40-bit unsigned integer, little-endian. Maximum N = 2^40 - 1.
Decoded as: `v[0] | (v[1] << 8) | (v[2] << 16) | (v[3] << 24) | (v[4] << 32)`.

The integer encoding field in byte 1 is set to a default value (`0x2`, Elias delta) but
is **not used** for encoding or decoding uncompressed graph data. It exists for header
format consistency with compressed MGS variants.

---

## Vertex Width

All vertex IDs and degree counts are stored as `n_size_t`-byte big-endian unsigned
integers, where:

```
n_size_t = ceil(log2(N) / 8)
```

| N (vertices) | n_size_t | Range |
|--------------|----------|-------|
| 1-255 | 1 byte | UInt8 |
| 256-65535 | 2 bytes | UInt16 |
| 65536-16M | 3 bytes | UInt24 |
| 16M-4G | 4 bytes | UInt32 |
| 4G-1T | 5 bytes | UInt40 |

Each value is written as `n_size_t` bytes in **big-endian** order (most significant
byte first). When the graph's internal type T is wider than `n_size_t`, the leading
`sizeof(T) - n_size_t` zero bytes are stripped.

---

## Graph Data: Children Mode (coding_scheme = 0b00)

In children mode, each vertex's adjacency list is written sequentially, terminated
by a zero stop sequence. Vertices are implicitly numbered 1 to N in order.

```
[12-byte header]
For vertex v = 1 to N:
  neighbor_1    (n_size_t bytes, big-endian)
  neighbor_2    (n_size_t bytes, big-endian)
  ...
  neighbor_k    (n_size_t bytes, big-endian)
  stop_sequence (n_size_t zero bytes: 0x00...00)
```

- Each neighbor is the raw vertex ID (1-based)
- The stop sequence is `n_size_t` bytes of `0x00`, marking the end of the list
- Vertices with no outgoing edges are represented by the stop sequence alone
- Neighbor order is preserved as stored in the graph (typically sorted ascending)

### Decoding

Read `n_size_t` bytes at a time. If all bytes are zero, advance to the next vertex.
Otherwise, interpret as a big-endian vertex ID and add edge (current_vertex -> ID).

### Example

4-vertex directed graph (n_size_t = 1 byte):
- Vertex 1 -> {2, 3}
- Vertex 2 -> {1}
- Vertex 3 -> {4}
- Vertex 4 -> {} (no outgoing edges)

```
Header: 4D 47 53 03 02 02 00 04 00 00 00 00
Data:   02 03 00    (v1: neighbors 2, 3; stop)
        01 00       (v2: neighbor 1; stop)
        04 00       (v3: neighbor 4; stop)
        00          (v4: no neighbors; stop)
```

Total data: 8 bytes. Total file: 20 bytes.

---

## Graph Data: Index Mode (coding_scheme = 0b01)

In index mode, an out-degree array precedes the flattened neighbor data. This
enables computing the position of any vertex's neighbors without scanning.

```
[12-byte header]
[Index section: N degree entries]
  degree[1]     (n_size_t bytes, big-endian)
  degree[2]     (n_size_t bytes, big-endian)
  ...
  degree[N]     (n_size_t bytes, big-endian)
[Data section: flattened neighbor lists]
  neighbors of vertex 1 (degree[1] entries)
  neighbors of vertex 2 (degree[2] entries)
  ...
  neighbors of vertex N (degree[N] entries)
```

- Index section: N entries, each the out-degree of vertex v (number of outgoing edges)
- Data section: all neighbor IDs concatenated without delimiters
- Each value (degree or neighbor ID) is `n_size_t` bytes, big-endian

### Random Access

To read vertex v's neighbors:
1. Read `degree[v]` from index section at byte offset `12 + (v-1) * n_size_t`
2. Compute data offset: `12 + N * n_size_t + sum(degree[1..v-1]) * n_size_t`
3. Read `degree[v]` consecutive `n_size_t`-byte entries

### Example

Same 4-vertex graph (n_size_t = 1 byte):

```
Header: 4D 47 53 03 02 12 00 04 00 00 00 00
                       ^^ coding_scheme = 0b01 (index mode)
Index:  02 01 01 00   (degrees: v1=2, v2=1, v3=1, v4=0)
Data:   02 03 01 04   (v1: 2,3; v2: 1; v3: 4; v4: nothing)
```

Total data: 8 bytes. Total file: 20 bytes.

---

## Size Analysis

### Children Mode

For a graph with N vertices and E edges:

```
File size = 12 + (E + N) * n_size_t bytes
```

Each edge contributes `n_size_t` bytes for the neighbor ID, and each vertex contributes
`n_size_t` bytes for its stop sequence.

### Index Mode

```
File size = 12 + (N + E) * n_size_t bytes
```

Each vertex contributes `n_size_t` bytes for its degree entry, and each edge contributes
`n_size_t` bytes for the neighbor ID. Both modes produce files of the same size.

### Bits Per Edge

```
BPE = 8 * file_size / E = 8 * (12 + (N + E) * n_size_t) / E
```

For CNR-2000 (N=325,557, E=3,216,152, n_size_t=3):

```
BPE = 8 * (12 + (325,557 + 3,216,152) * 3) / 3,216,152 = 8 * 10,625,139 / 3,216,152 = 26.43
```

This is the baseline against which compressed formats (BG, CS, CG) are measured.

---

## Relationship to Compressed Formats

The uncompressed MGS format shares only the 12-byte header with compressed formats.
After the header:

| Format | After header |
|--------|-------------|
| **Uncompressed** | Raw fixed-width adjacency lists (this spec) |
| **BG/CS** | 4-bit VLC stream header + bitstream-encoded vertex records |
| **CG** | Cluster offset table (index mode) or direct bitstream + clustered encoding |

The option_flags byte (byte 2) distinguishes the format:
- `0x00`: Uncompressed (ALG_LEGACY_MGS)
- `0x02`: BG with defaults
- `0x03`: CS with defaults
- `0x04`: CG with defaults
- `0x10-0xFF`: BG/CS/CG with explicit parameters

See `MGS_HEADER.md` for full byte 2 encoding, and `FORMAT_SPECS_BG.md`,
`FORMAT_SPECS_CS.md`, `FORMAT_SPECS_CG.md` for compressed format specifications.
