# Adjancently.jl

Adjancently.jl is Julia library for the analysis of large complex directed networks.

## ASTRA Compression

**ASTRA** (**A**daptive **STR**eaming **A**djacency) is a graph compression method for large-scale directed networks, inspired by the [WebGraph framework](http://webgraph.di.unimi.it/).

### Key Features

- **Greedy Cost-Based Reference Selection**: Unlike WebGraph's overlap-based heuristics, ASTRA evaluates all candidate references in a sliding window and selects the one with minimum encoding cost
- **Adaptive Bitmap Encoding**: Automatically chooses between block encoding (for sparse patterns) and raw encoding (for dense patterns) based on actual bit costs
- **Recursive Reference Compression**: Supports multi-level reference chains for graphs with strong locality
- **Variable-Length Integer Encoding**: Supports multiple schemes (Fibonacci, Elias-Gamma/Delta, Golomb, Zeta, FED)
- **Mix Encoding**: Combines interval encoding and run-length encoding for residual neighbors

### Performance

On the CNR-2000 web graph (325,557 vertices, 3,216,152 edges):
- **5.108 bits per edge**
- **1.57 edges per byte** compression ratio
- Perfect round-trip fidelity

### Format

ASTRA graphs are stored in the MGS (Modified Graph Structure) format with the `OPTION_ASTRA` flag, which enables:
1. Delta encoding for sorted adjacency lists
2. Mix encoding (intervals + run-length) for residuals
3. Greedy reference selection with adaptive bitmap compression
4. Recursive reference support

### References

- [WebGraph: A Framework for Graph Compression](http://webgraph.di.unimi.it/) - Paolo Boldi and Sebastiano Vigna
- [The WebGraph Framework](https://dl.acm.org/doi/10.1145/988672.988752) - WWW 2004

## Tests

The test suite is organized into several distinct test sets that can be run individually or all together. Use the following commands from the project root directory:

### Run tests
```
julia --project=. test/run_tests_{test-name}.jl

```

### Run tests with verbose output

Add the `-v` flag for verbose output:
```
julia -v run_tests_{test-name}.jl
```

### Run tests in parallel
Use the `-p` flag to run tests in parallel (requires multiple processors):
```
julia -p auto run_tests_{test-name}.jl
```

### List available test sets

To see all available test sets, run:

```
julia --project -e 'using Test; include("test/runtests.jl"); println("\nAvailable test sets:"); for ts in Test.get_testset_string() println("- ", ts) end'
```

## Development

### Launch notebooks
```
julia> using IJulia
julia> notebook()
```

### Dependencies management
```
pkg> activate .
pkg> add {package-name}
pkg> update
```

### Datasets

- [Laboratory for Web Algorithmics - Datasets](https://law.di.unimi.it/datasets.php)
- [WebGraph](https://webgraph.di.unimi.it/)
- [WebGraph - Github](https://github.com/vigna/webgraph)
- [Stanford WebBase Project](http://diglib.stanford.edu:8091/~testbed/doc2/WebBase/)
