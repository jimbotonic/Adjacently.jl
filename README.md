# Adjancently.jl 

Adjancently.jl is Julia library for the analysis of large complex directed networks.

## Tests

The test suite is organized into several distinct test sets that can be run individually or all together. Use the following commands from the project root directory:

### Run all tests
```
julia test/runtests.jl
```

### Run specific test sets
You can run individual test sets by using the `--testset` flag:

```
# Test sorting algorithms
julia test/runtests.jl --testset "Sorting Algorithms"

# Test binary search implementations
julia test/runtests.jl --testset "Binary Search"

# Test Huffman encoding
julia test/runtests.jl --testset "Huffman Encoding"

# Test graph serialization
julia test/runtests.jl --testset "Graph Serialization"

# Test Pajek graph format
julia test/runtests.jl --testset "Pajek Graph Format"
```

### Run tests with verbose output

Add the `-v` flag for verbose output:
```
julia -v test/runtests.jl
```

### Run tests in parallel
Use the `-p` flag to run tests in parallel (requires multiple processors):
```
julia -p auto test/runtests.jl
```

### List available test sets
To see all available test sets, run:

```
julia --project -e 'using Test; include("test/runtests.jl"); println("\nAvailable test sets:"); for ts in Test.get_testset_string() println("- ", ts) end'
```

This will output something like:

```
Available test sets:
- Sorting Algorithms
- Binary Search
- Huffman Encoding
- Graph Serialization
- Pajek Graph Format
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

### Programming notes

```julia
###
# Write in big endian
###

# reinterpret pos in an array of bytes
bytes = reinterpret(UInt8, [p]) 
#
# See for example:
# c = UInt16(1)
# bytes = reinterpret(UInt8, [c])
# 2-element reinterpret(UInt8, ::Vector{UInt16}):
# 0x01
# 0x00

# write bytes in big endian
write(f, reverse(bytes))

### 
# Read in big endian
###

# See for example:
# bytes = read(f, 2)
# pos = reinterpret(UInt16, reverse(bytes))

# read sizeof(T) bytes from file
child = read(f, UInt8, sizeof(T))

# add value to children list
# NB: `reinterpet` takes as input an array of bytes in little endian
# NB: `reinterpret` returns an array of T value -> select the first and only one
push!(children, reinterpret(T, reverse(child))[1])

```
