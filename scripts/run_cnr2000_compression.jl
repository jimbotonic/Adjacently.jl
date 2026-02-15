#
# Adjacently: Julia Complex Directed Networks Library
# Copyright (C) 2016-2025 Jimmy Dubuisson <jimmy.dubuisson@gmail.com>
#
# This program is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.
#
# This program is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.
#

#!/usr/bin/env julia

# Run legacy compression on CNR-2000 and measure bits per edge.
# Expected result: ~5.5 bits/edge with RCM + Fibonacci + ASTRA encoding.

using Pkg
Pkg.activate(normpath(joinpath(@__DIR__, "..")))

using LightGraphs: nv, ne
using Adjacently
using Adjacently.IO: load_adjacency_list_from_csv
using Adjacently.Relabeling: relabel_vertices_rcm, relabel_graph
using Adjacently.MGS: write_compressed_mgs3_graph

const PROJECT_ROOT = normpath(joinpath(@__DIR__, ".."))
const CNR_CSV = joinpath(PROJECT_ROOT, "datasets", "webgraph", "cnr-2000", "cnr-2000.csv")
const OUTPUT_DIR = joinpath(PROJECT_ROOT, "test_data")
const OUTPUT_BASE = joinpath(OUTPUT_DIR, "cnr2000_benchmark")

function main()
    println("=" ^ 70)
    println("CNR-2000 Legacy Compression Benchmark")
    println("=" ^ 70)

    # 1. Load graph
    println("\nLoading CNR-2000 dataset...")
    if !isfile(CNR_CSV)
        error("CNR-2000 CSV not found at $CNR_CSV")
    end
    t0 = time()
    g = load_adjacency_list_from_csv(CNR_CSV, ',', true)
    load_time = time() - t0
    n, m = nv(g), ne(g)
    println("  Vertices: $n")
    println("  Edges:    $m")
    println("  Load time: $(round(load_time, digits=2))s")

    # 2. RCM relabeling
    println("\nApplying RCM relabeling...")
    t1 = time()
    rcm_mapping = relabel_vertices_rcm(g, :out)
    g_rcm = relabel_graph(g, rcm_mapping)
    relabel_time = time() - t1
    println("  Relabel time: $(round(relabel_time, digits=2))s")

    # 3. Compress
    mkpath(OUTPUT_DIR)
    println("\nCompressing with ASTRA encoding...")
    println("  Integer encoding: fibonacci")
    println("  Reference enabled: true (recursive)")
    println("  Window size: 1024")
    t2 = time()
    write_compressed_mgs3_graph(
        g_rcm,
        OUTPUT_BASE,
        :children,      # coding_scheme
        :complex,        # compression
        :fibonacci,      # integer_encoding
        true,            # use_mix_mode
        true,            # reference_enabled
        true,            # recursive_reference
        1024             # ref_window_size
    )
    compress_time = time() - t2

    # 4. Measure compression ratio
    output_file = OUTPUT_BASE * ".mgz"
    file_size = filesize(output_file)
    bits_per_edge = 8.0 * file_size / m

    println("\n" * "=" ^ 70)
    println("RESULTS")
    println("=" ^ 70)
    println("  File size:        $(round(file_size / (1024*1024), digits=3)) MB")
    println("  Bits per edge:    $(round(bits_per_edge, digits=3))")
    println("  Compression time: $(round(compress_time, digits=2))s")
    println("  Edges per byte:   $(round(m / file_size, digits=2))")
    println("")
    println("  WebGraph ref:     2.90 bits/edge")
    println("  Gap:              $(round(bits_per_edge - 2.90, digits=3)) bits/edge")
    println("=" ^ 70)

    # 5. Cleanup
    if isfile(output_file)
        rm(output_file)
        println("\nCleaned up $output_file")
    end
end

main()
