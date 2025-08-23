# Adjancently.jl 

Adjancently.jl is Julia library for the analysis of large complex directed networks.

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
