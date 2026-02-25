# MGS File Header

Header has 12 bytes: 
- version: 'MGS' 3 bytes string
- major + minor version: 2 bytes
- flags: 2 bytes
  * Byte 1 (2 bits + 2 bits + 4 bits):
	- graph type (2 bits): 0b00: directed graph | 0b01: undirected graph
	- coding scheme (2 bits): 0b00: children | 0b01: index
	- integer encoding (4 bits): 0x1: Elias gamma | 0x2: Elias delta | 0x3: Golomb | 0x4: FED | 0x5: Zeta | 0x6: Fibonacci
  * Byte 2:
    - option flags:
      - 0x00: no option
      - 0x01: complex [delta only] (deprecated tag; prefer ASTRA)
      - 0x03: complex [delta + mix (run-length + interval)] (deprecated tag)
      - 0x07: complex [delta + mix + reference] (deprecated tag)
      - 0x0F: ASTRA [delta + mix + recursive reference]
      - 0x10–0x8F: Greedy compressed mode (128 slots, reserved for future use)
        - 0x10: greedy mode 1 (default), up to 0x8F: greedy mode 128
        - Per-vertex encoding decisions via greedy cost-based search
        - Stream is self-describing (encoding tag in bitstream header)
      - 0xFF: Huffman (deprecated)
- # vertices: 5 bytes position 

The header format is thus:
<'MGS' string 3 bytes> + <16 bits major|minor version> + <flags 2 bytes> + <# vertices 5 bytes>
