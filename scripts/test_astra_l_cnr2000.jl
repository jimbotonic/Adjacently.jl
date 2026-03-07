#
# Adjacently: Julia Complex Directed Networks Library
# Copyright (C) 2016-2026 Jimmy Dubuisson <jimmy@dubuisson.ch>
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

#
# Test ASTRA-L (Layered) compression on CNR-2000 and measure BPE
#

using Pkg
Pkg.activate(normpath(joinpath(@__DIR__, "..")))

using Logging
global_logger(ConsoleLogger(stderr, Logging.Info))

using LightGraphs: nv, ne, outneighbors
using Adjacently
using Adjacently.IO: load_adjacency_list_from_csv
using Adjacently.Graph: get_reverse_graph
using Adjacently.Compression.ASTRALayered: create_level_decomposition, write_astra_layered_graph

function main()
    println("=" ^ 70)
    println("ASTRA-L CNR-2000 Compression Test")
    println("=" ^ 70)

    # 1. Load CNR-2000 from CSV
    cnr_csv = joinpath(@__DIR__, "..", "datasets", "webgraph", "cnr-2000", "cnr-2000.csv")
    cnr_csv = normpath(cnr_csv)
    if !isfile(cnr_csv)
        error("CNR-2000 CSV not found at $cnr_csv")
    end

    println("\nLoading CNR-2000...")
    t0 = time()
    g = load_adjacency_list_from_csv(cnr_csv, ',', true)
    n = nv(g)
    m = ne(g)
    println("  Loaded: $n vertices, $m edges ($(round(time()-t0, digits=2))s)")

    # 2. Compute reverse graph (needed for merge incoming-edge capture)
    println("\nComputing reverse graph...")
    t1 = time()
    rg = get_reverse_graph(g)
    println("  Done ($(round(time()-t1, digits=2))s)")

    # 3. Build decomposition
    #    use_scc=false  → use full radius-k ball (not just SCC) for much larger levels
    #    min_level_size=50 → merge tiny balls into previous level
    println("\nBuilding level decomposition (radius=2, full-ball, min_level_size=50)...")
    t2 = time()
    decomp = create_level_decomposition(g, rg;
        radius=2, log_every=5000,
        use_scc=false, min_level_size=50)
    t_decomp = time() - t2

    # Statistics
    total_intra = sum(length(l.intra_edges) for l in decomp.levels)
    total_cross = length(decomp.remaining_cross_edges)
    level_sizes = [length(l.local_to_global) for l in decomp.levels]

    println("  Levels: $(length(decomp.levels))")
    println("  Level sizes: min=$(minimum(level_sizes)), max=$(maximum(level_sizes)), mean=$(round(sum(level_sizes)/length(level_sizes), digits=1))")
    println("  Intra-level edges: $total_intra ($(round(100*total_intra/m, digits=1))%)")
    println("  Cross-level edges: $total_cross ($(round(100*total_cross/m, digits=1))%)")
    println("  Decomposition time: $(round(t_decomp, digits=2))s")

    # 4. Write compressed file (with reference + interval encoding, no backward pass)
    output_file = joinpath(@__DIR__, "..", "datasets", "webgraph", "cnr-2000", "cnr2000_astra_l.astral")
    output_file = normpath(output_file)
    mkpath(dirname(output_file))

    println("\nWriting ASTRA-L file (Fibonacci, ref_window=7)...")
    t3 = time()
    write_astra_layered_graph(decomp, output_file, g;
        integer_encoding=:fibonacci, ref_window_size=7, log_every=5000)
    t_write = time() - t3
    println("  Write time: $(round(t_write, digits=2))s")

    # 5. Measure file size and BPE
    fsize = filesize(output_file)
    bpe = 8.0 * fsize / m

    println("\n" * "=" ^ 70)
    println("RESULTS")
    println("=" ^ 70)
    println("  Output file: $output_file")
    println("  File size:   $(fsize) bytes ($(round(fsize / 1024.0, digits=1)) KB)")
    println("  Edges:       $m")
    println("  Bits/edge:   $(round(bpe, digits=3))")
    println("  Comparison:  legacy ASTRA 5.108, CGE 4.753")
    println("=" ^ 70)
end

main()
