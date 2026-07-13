#!/usr/bin/env bash
#
# fetch_datasets.sh -- download the LAW web graphs that are NOT committed to the
# repository and convert them to the comma-separated edge-list CSV that the
# Julia drivers load (via Adjacently.IO.load_adjacency_list_from_csv).
#
# These are needed to unblock the reproduction tables that the committed
# datasets cannot cover (see README.md):
#   - cnr-2000   : the committed cnr-2000.mgz is a pre-reordered graph, NOT the
#                  LAW crawl order, so it cannot reproduce the cnr Orig./native
#                  cells. This fetches the true LAW cnr-2000 (its default,
#                  LLP-optimized order -- the paper's "Orig.(LLP)").
#   - in-2004    : ~16.9M edges, too large to commit.
#   - enwiki-2013: ~101M edges (166 MB .graph), too large to commit.
#
# The LAW graphs are WebGraph BVGraph binaries; converting them to an edge list
# needs the WebGraph Java library (the paper used WebGraph 3.6.12, Java 21).
# Point WEBGRAPH_CP at a classpath containing webgraph.jar and its dependencies
# (dsiutils, fastutil, sux4j, slf4j, jsap, commons-*). Without it, files are
# downloaded and verified but the conversion step prints the manual command.
#
# Usage:
#   bash bench/graph_compression/fetch_datasets.sh [dataset ...]   # default: all
#   WEBGRAPH_CP=/path/to/webgraph/'*' bash bench/graph_compression/fetch_datasets.sh cnr-2000
#
set -euo pipefail

REPO_ROOT="$(cd "$(dirname "${BASH_SOURCE[0]}")/../.." && pwd)"
HOST="http://data.law.di.unimi.it/webdata"     # NB: the data. subdomain, not law.di.unimi.it
ALL=(cnr-2000 in-2004 enwiki-2013)
if [ "$#" -gt 0 ]; then DATASETS=("$@"); else DATASETS=("${ALL[@]}"); fi

for name in "${DATASETS[@]}"; do
    dir="$REPO_ROOT/datasets/webgraph/$name"
    mkdir -p "$dir"
    echo "=== $name -> datasets/webgraph/$name ==="

    # 1. download .graph + .properties + .md5sums (skip if already present)
    for ext in graph properties md5sums; do
        out="$dir/$name.$ext"
        if [ -f "$out" ]; then
            echo "  have $name.$ext"
        else
            echo "  downloading $name.$ext ..."
            curl -fL --retry 3 --retry-delay 2 -o "$out" "$HOST/$name/$name.$ext"
        fi
    done

    # 2. verify md5 of the two files we use (LAW ships a .md5sums for all variants)
    ( cd "$dir" && grep -E "  $name\.(graph|properties)\$" "$name.md5sums" | md5sum -c - ) \
        && echo "  md5 OK" \
        || { echo "  !! md5 verification FAILED for $name"; exit 1; }

    # 3. convert BVGraph -> arc list -> src,dst CSV (sequential; no .offsets needed)
    csv="$dir/$name.csv"
    if [ -f "$csv" ]; then
        echo "  have $name.csv"
    elif [ -n "${WEBGRAPH_CP:-}" ]; then
        echo "  converting BVGraph -> $name.csv ..."
        tmp="$(mktemp)"
        java -Xmx8G -cp "$WEBGRAPH_CP" it.unimi.dsi.webgraph.ArcListASCIIGraph \
            -g it.unimi.dsi.webgraph.BVGraph "$dir/$name" "$tmp"
        { echo "src,dst"; awk 'NF>=2 {print $1","$2}' "$tmp"; } > "$csv"
        rm -f "$tmp"
        echo "  wrote $name.csv ($(wc -l < "$csv") lines incl. header)"
    else
        echo "  WEBGRAPH_CP not set -- downloaded only. To finish the conversion:"
        echo "    java -cp '<webgraph-3.6.12 + deps>' it.unimi.dsi.webgraph.ArcListASCIIGraph \\"
        echo "         -g it.unimi.dsi.webgraph.BVGraph '$dir/$name' /tmp/$name.arcs"
        echo "    { echo src,dst; awk 'NF>=2{print \$1\",\"\$2}' /tmp/$name.arcs; } > '$csv'"
    fi
done

echo
echo "Done. Load a fetched graph in a driver with:"
echo "  Adjacently.IO.load_adjacency_list_from_csv(\"datasets/webgraph/<name>/<name>.csv\", ',', true)"
echo "Note: cnr-2000.csv is the true LAW order; use it (not the committed .mgz) for the cnr Orig rows."
