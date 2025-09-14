# MGS File Header

Header has 12 bytes: 
- version: 'MGS' 3 bytes string
- major + minor version: 2 bytes
- flags: 2 bytes
  * Byte 1 (2 bits + 2 bits + 6 bits):
	- graph type (2 bits): 0b00: directed graph | 0b01: undirected graph
	- coding scheme (2 bits): 0b00: children | 0b01: index
	- integer encoding (4 bits): 0x1: Elias gamma | 0x2: Elias delta | 0x3: Golomb | 0x4: FED | 0x5: Zeta | 0x6: Fibonacci
  * Byte 2:
	- option flags: 		
    	- 0000 0000: no option
    	- 0000 0001: complex [delta encoding only] - so called "delta"
    	- 0000 0011: complex [delta + mix encoding (run-length + interval)] - so called "mix"
    	- 0000 0111: complex [delta + mix + reference only] - so called "hybrid"
    	- 0000 1111: complex [delta + mix + recursive reference] - so called "hybrid+"
    	- 1000 0000: Huffman - so called "huffman"
- # vertices: 5 bytes position 

The header format is thus:
<'MGS' string 3 bytes> + <16 bits major|minor version> + <flags 2 bytes> + <# vertices 5 bytes>
