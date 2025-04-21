# Adjancently.jl 

Adjancently.jl is Julia library for the analysis of large complex directed networks.


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

### Notes
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
