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
      - 0x0F: ASTRA (legacy adaptive streaming)
      - 0x10–0x8F: STD (Standard greedy compressed mode)
      - 0x90–0x9F: CS (Command Stream compressed mode)
      - 0xA0–0xAF: RCGE (Reversible Coarsening Graph Encoding)
      - 0xFF: Huffman (deprecated)
- # vertices: 5 bytes position 

The header format is thus:
<'MGS' string 3 bytes> + <16 bits major|minor version> + <flags 2 bytes> + <# vertices 5 bytes>
