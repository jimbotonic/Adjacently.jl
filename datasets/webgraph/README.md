# WebGraph Datasets

**Source:** Laboratory for Web Algorithmics (LAW), University of Milan
**URL:** https://law.di.unimi.it/datasets.php

This directory contains web graph datasets in BVGraph format, converted to CSV for use with Adjacently.jl.

## Datasets

### CNR-2000

**URL:** https://law.di.unimi.it/webdata/cnr-2000/
**Downloaded:** August 14, 2025

| Property | Value |
|----------|-------|
| Nodes | 325,557 |
| Arcs (Edges) | 3,216,152 |
| Domain | Italian CNR (Consiglio Nazionale delle Ricerche) websites |
| Average Degree | 9.879 |
| BVGraph Bits per Link | 2.897 |
| BVGraph Compression Ratio | 0.176 |
| Window Size | 7 |
| Zeta K | 3 |

**Core files:**
- `cnr-2000/cnr-2000.graph` (1.16MB) — BVGraph compressed format
- `cnr-2000/cnr-2000.properties` (982B) — graph metadata
- `cnr-2000/cnr-2000.md5sums` (1.1KB) — checksums (verified)

**Generated files:**
- `cnr-2000/cnr-2000.offsets` (318KB) — random access index
- `cnr-2000/cnr-2000.obl` (282KB) — offset big list
- `cnr-2000/cnr-2000.csv` (40.8MB) — CSV edge list (3,216,153 rows: header + edges)
- `cnr-2000/cnr-2000-edges` (40.8MB) — tab-separated edge list

**File integrity verified:**
- `cnr-2000.graph`: `24539022b181013f7930a2f45010f6de`
- `cnr-2000.properties`: `8b78a939fe399fb8784bcfe580c4e943`

---

### IN-2004

**URL:** https://law.di.unimi.it/webdata/in-2004/
**Downloaded:** February 25, 2026

| Property | Value |
|----------|-------|
| Nodes | 1,382,908 |
| Arcs (Edges) | 16,917,053 |
| Domain | Indian web domains (.in) |
| Average Degree | 12.232 |
| BVGraph Bits per Link | 2.172 |
| BVGraph Compression Ratio | 0.119 |
| Window Size | 7 |
| Zeta K | 3 |

**Core files:**
- `in-2004/in-2004.graph` (4.4MB) — BVGraph compressed format
- `in-2004/in-2004.properties` (1.1KB) — graph metadata
- `in-2004/in-2004.md5sums` (1.3KB) — checksums (verified)

**Generated files:**
- `in-2004/in-2004.offsets` (1.4MB) — random access index
- `in-2004/in-2004.obl` (1.2MB) — offset big list
- `in-2004/in-2004.csv` (236MB) — CSV edge list (16,917,054 rows: header + edges)
- `in-2004/in-2004-edges` (236MB) — tab-separated edge list

**File integrity verified:**
- `in-2004.graph`: `40ddfce0454afb7e85bdae197c18527b`
- `in-2004.properties`: `50f76c0c70d00ba0f6c0d1dd14c92c67`

---

### enwiki-2013

**URL:** https://law.di.unimi.it/webdata/enwiki-2013/
**Downloaded:** February 28, 2026

| Property | Value |
|----------|-------|
| Nodes | 4,206,785 |
| Arcs (Edges) | 101,355,853 |
| Domain | English Wikipedia (February 2013 snapshot) |
| Average Degree | 24.093 |
| BVGraph Bits per Link | 12.639 |
| BVGraph Compression Ratio | 0.695 |
| Window Size | 7 |
| Zeta K | 3 |

**Core files:**
- `enwiki-2013/enwiki-2013.graph` (159MB) — BVGraph compressed format
- `enwiki-2013/enwiki-2013.properties` (1.2KB) — graph metadata
- `enwiki-2013/enwiki-2013.md5sums` (1.4KB) — checksums (verified)

**Generated files:**
- `enwiki-2013/enwiki-2013.offsets` (7.5MB) — random access index
- `enwiki-2013/enwiki-2013.obl` (5.3MB) — offset big list
- `enwiki-2013/enwiki-2013.csv` (1.5GB) — CSV edge list (101,355,854 rows: header + edges)
- `enwiki-2013/enwiki-2013-edges` (1.5GB) — tab-separated edge list

**File integrity verified:**
- `enwiki-2013.graph`: `b8a88215f86895bb6ddb7f9ed175b7f1`
- `enwiki-2013.properties`: `0c99b4e15e2eb24a1fb69f9e496e54ca`

---

## Processing a New Dataset

The steps below show how to download and convert any LAW BVGraph dataset to CSV. Replace `DATASET` with the dataset name (e.g. `cnr-2000`, `in-2004`).

### Prerequisites

- WebGraph framework at `$WEBGRAPH_FOLDER` (see `apps/webgraph/BUILD.md`)
- Java 9+ with sufficient memory

### 1. Download core files

```bash
mkdir -p $PROJECT_FOLDER/datasets/webgraph/DATASET
cd $PROJECT_FOLDER/datasets/webgraph/DATASET

wget http://data.law.di.unimi.it/webdata/DATASET/DATASET.graph
wget http://data.law.di.unimi.it/webdata/DATASET/DATASET.properties
wget http://data.law.di.unimi.it/webdata/DATASET/DATASET.md5sums
```

### 2. Verify checksums

```bash
md5sum -c <(grep -E '(DATASET\.graph|DATASET\.properties)$' DATASET.md5sums)
```

### 3. Generate offset files

```bash
java -Xmx4G -cp "$WEBGRAPH_FOLDER/webgraph-3.6.12.jar:$WEBGRAPH_FOLDER/jars/runtime/*" \
  it.unimi.dsi.webgraph.BVGraph -o -O -L DATASET
```

Creates `DATASET.offsets` and `DATASET.obl`.

### 4. Export to CSV

```bash
# Export to tab-separated edge list
java -Xmx4G -cp "$WEBGRAPH_FOLDER/webgraph-3.6.12.jar:$WEBGRAPH_FOLDER/jars/runtime/*" \
  it.unimi.dsi.webgraph.ArcListASCIIGraph -g BVGraph DATASET DATASET-edges

# Convert to CSV with header
echo "source,target" > DATASET.csv
sed 's/\t/,/g' DATASET-edges >> DATASET.csv
```

### 5. Verify edge count

The CSV should have one header row plus the number of arcs listed in `DATASET.properties`:

```bash
wc -l DATASET.csv
# Should equal arcs + 1
```

### 6. Load in Java (optional)

```java
import it.unimi.dsi.webgraph.*;

ImmutableGraph graph = BVGraph.loadMapped("DATASET");
System.out.println("Nodes: " + graph.numNodes());
System.out.println("Arcs: " + graph.numArcs());
```

### 7. Generate statistics (optional)

```bash
java -Xmx4G -cp "$WEBGRAPH_FOLDER/webgraph-3.6.12.jar:$WEBGRAPH_FOLDER/jars/runtime/*" \
  it.unimi.dsi.webgraph.Stats -s DATASET DATASET-detailed
```

## Citation

When using these datasets, please cite:

> Paolo Boldi and Sebastiano Vigna. The WebGraph framework I: Compression techniques. In Proc. of the 13th International World Wide Web Conference (WWW 2004), pages 595-601, Manhattan, USA, 2004. ACM Press.

## References

- [LAW Datasets](https://law.di.unimi.it/datasets.php) — full list of available web graphs
- [WebGraph framework](https://webgraph.di.unimi.it/)
- [WebGraph on GitHub](https://github.com/vigna/webgraph)
