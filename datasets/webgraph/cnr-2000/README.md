# CNR-2000 WebGraph Dataset

## Dataset Information

**Source:** Laboratory for Web Algorithmics (LAW), University of Milan  
**URL:** https://law.di.unimi.it/webdata/cnr-2000/  
**Downloaded:** August 14, 2025

### Graph Properties

- **Nodes:** 325,557
- **Arcs (Edges):** 3,216,152
- **Graph Type:** Directed web graph
- **Domain:** Italian CNR (Consiglio Nazionale delle Ricerche) websites
- **Average Degree:** 9.879 (3,216,152 / 325,557)

### Compression Statistics

From `cnr-2000.properties`:
- **Compression Ratio:** 0.176 (17.6% of original size)
- **Bits per Link:** 2.897
- **Bits per Node:** 28.624
- **Window Size:** 7 (reference locality parameter)
- **Zeta K:** 3 (compression parameter)
- **Average Reference Distance:** 1.64

### File Contents

**Core Files:**
- `cnr-2000.graph` (1.16MB) - Compressed graph in BVGraph format
- `cnr-2000.properties` (982 bytes) - Graph metadata and compression parameters
- `cnr-2000.md5sums` (1.1KB) - Checksums for integrity verification

**Generated Files:**
- `cnr-2000.offsets` (318KB) - Random access index
- `cnr-2000.obl` (282KB) - Offset big list
- `cnr-2000.csv` (40.8MB) - Edge list in CSV format (3,216,153 rows)
- `cnr-2000-edges` (40.8MB) - Tab-separated edge list
- `cnr-2000-comprehensive-stats.txt` - Detailed graph analysis
- `cnr-2000.indegree` / `cnr-2000.outdegree` - Degree distribution files

**File Integrity:** ✅ Verified (checksums match)
- `cnr-2000.graph`: `24539022b181013f7930a2f45010f6de`
- `cnr-2000.properties`: `8b78a939fe399fb8784bcfe580c4e943`

## Usage

### Prerequisites
- WebGraph framework (see `$WEBGRAPH_FOLDER/BUILD.md`)
- Java 9+ with sufficient memory

### Generate Offset Files
```bash
cd $PROJECT_FOLDER/datasets/webgraph/cnr-2000
java -Xmx2G -cp "$WEBGRAPH_FOLDER/webgraph-3.6.12.jar:$WEBGRAPH_FOLDER/jars/runtime/*" \
  it.unimi.dsi.webgraph.BVGraph -o -O -L cnr-2000
```

**Output:**
```
15:17:59.229 [main] INFO it.unimi.dsi.webgraph.BVGraph - Writing offsets...
15:17:59.464 [main] INFO it.unimi.dsi.webgraph.BVGraph - Completed.
15:17:59.493 [main] INFO it.unimi.dsi.webgraph.BVGraph - Elapsed: 229ms [325,557 items, 1,427,881.58 items/s, 700.34 ns/item]
```

**Files created:** `cnr-2000.offsets` (318KB), `cnr-2000.obl` (282KB)

### Load in Java
```java
import it.unimi.dsi.webgraph.*;

ImmutableGraph graph = BVGraph.loadMapped("cnr-2000");
System.out.println("Nodes: " + graph.numNodes());
System.out.println("Arcs: " + graph.numArcs());
```

### Export to CSV Format

**Export to edge list (tab-separated):**
```bash
java -Xmx2G -cp "$WEBGRAPH_FOLDER/webgraph-3.6.12.jar:$WEBGRAPH_FOLDER/jars/runtime/*" \
  it.unimi.dsi.webgraph.ArcListASCIIGraph -g BVGraph cnr-2000 cnr-2000-edges
```

**Convert to CSV with header:**
```bash
echo "source,target" > cnr-2000.csv
sed 's/\t/,/g' cnr-2000-edges >> cnr-2000.csv
```

**Files created:** 
- `cnr-2000-edges` (40.8MB) - Tab-separated: 3,216,152 edges
- `cnr-2000.csv` (40.8MB) - CSV format: 3,216,153 rows (header + edges)

### Generate Statistics

**Built-in WebGraph statistics:**
```bash
java -Xmx2G -cp "$WEBGRAPH_FOLDER/webgraph-3.6.12.jar:$WEBGRAPH_FOLDER/jars/runtime/*" \
  it.unimi.dsi.webgraph.Stats -s cnr-2000 cnr-2000-detailed
```

**Comprehensive custom statistics:**
```bash
# Compile custom stats program
javac -cp "$WEBGRAPH_FOLDER/webgraph-3.6.12.jar:$WEBGRAPH_FOLDER/jars/runtime/*" CNRStats.java

# Run comprehensive analysis
java -Xmx2G -cp "$WEBGRAPH_FOLDER/webgraph-3.6.12.jar:$WEBGRAPH_FOLDER/jars/runtime/*:." \
  CNRStats > cnr-2000-comprehensive-stats.txt
```

**Key Statistics Generated:**
- **Nodes:** 325,557
- **Edges:** 3,216,152
- **Maximum outdegree:** 2,716
- **Maximum indegree:** 18,235
- **Dangling nodes:** 78,056 (23.98%)
- **Connected components:** 7
- **Largest component:** 325,239 nodes (99.90%)
- **Total compressed size:** 1.7MB
- **Bits per arc:** 4.426

### Additional Export Formats
See detailed instructions in `$WEBGRAPH_FOLDER/DATASETS.md`

## Research Context

This dataset represents a web crawl of Italian CNR (National Research Council) domains from 2000. It's commonly used for:

- **Algorithm testing** - Manageable size for debugging
- **Compression benchmarking** - Good example of web graph properties
- **Graph analysis research** - Real-world web link structure

### Citation

When using this dataset, please cite:

> Paolo Boldi and Sebastiano Vigna. The WebGraph framework I: Compression techniques. In Proc. of the 13th International World Wide Web Conference (WWW 2004), pages 595–601, Manhattan, USA, 2004. ACM Press.

## Technical Details

**Graph Format:** BVGraph (Boldi-Vigna compressed format)
**Compression Technique:** 
- Reference compression with locality
- Interval compression for consecutive nodes
- Zeta coding for gap encoding
- Block-based organization

**Storage Efficiency:**
- Original uncompressed: ~25MB (estimated)
- BVGraph compressed: 1.16MB (graph file only)
- Total with indices: 1.7MB (graph + offsets + obl)
- Compression ratio: 17.6%
- Space per edge: 2.897 bits (graph only), 4.426 bits (with indices)

## Processing Workflow Summary

This dataset has been fully processed with the following workflow:

1. **Download** - Core files downloaded and verified (August 14, 2025)
2. **Index Generation** - Offset files created for random access (229ms)
3. **CSV Export** - Edge list exported to CSV format (3.2M edges)
4. **Statistics** - Comprehensive analysis completed

**Total Processing Time:** ~5 minutes  
**Total Storage Used:** ~85MB (including all generated files)

## Additional Information

This dataset is part of a larger collection of web graphs available at:
https://law.di.unimi.it/datasets.php

For more datasets and tools, see the parent directory structure in this repository.

## Files Reference

All commands in this README assume you are in the CNR-2000 dataset directory:
```bash
cd $PROJECT_FOLDER/datasets/webgraph/cnr-2000
```

For building the WebGraph framework, see: `$WEBGRAPH_FOLDER/BUILD.md`  
For additional dataset instructions, see: `$WEBGRAPH_FOLDER/DATASETS.md`