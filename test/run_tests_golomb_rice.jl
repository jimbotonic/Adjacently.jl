#
# Adjacently: Julia Complex Directed Networks Library
# Copyright (C) 2016-2026 Anonymous (double-blind review)
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

# Test Golomb-Rice encoding/decoding

using Pkg
Pkg.activate(normpath(joinpath(@__DIR__, "..")))

using Adjacently.Compression: write_golomb_rice, read_golomb_rice
using Adjacently.IO: BitWriter, BitReader, flush_bitwriter

println("=" ^ 80)
println("Testing Golomb-Rice Encoding/Decoding")
println("=" ^ 80)

# Test function
function test_golomb_rice(values::Vector{T}, k::Int) where {T<:Unsigned}
    println("\nTesting k=$k with values: $values")

    # Write
    buffer = IOBuffer()
    w = BitWriter(buffer)

    for val in values
        write_golomb_rice(w, val, k)
    end
    flush_bitwriter(w, flush_last_bits=true)

    # Get encoded size
    encoded_bytes = position(buffer)
    encoded_bits = encoded_bytes * 8
    println("  Encoded size: $encoded_bytes bytes ($encoded_bits bits)")
    println("  Bits per value: $(round(encoded_bits / length(values), digits=2))")

    # Read back
    seekstart(buffer)
    r = BitReader(buffer)

    decoded = T[]
    for _ in 1:length(values)
        push!(decoded, read_golomb_rice(r, k, T))
    end

    # Verify
    if decoded == values
        println("  ✅ PASS: All values decoded correctly")
        return true
    else
        println("  ❌ FAIL: Mismatch!")
        println("    Original: $values")
        println("    Decoded:  $decoded")
        for i in 1:length(values)
            if values[i] != decoded[i]
                println("    Index $i: $(values[i]) != $(decoded[i])")
            end
        end
        return false
    end
end

# Test cases
all_passed = true

println("\n" * "=" ^ 80)
println("Test 1: Small values with k=0 (equivalent to unary)")
println("=" ^ 80)
all_passed &= test_golomb_rice(UInt8[0, 1, 2, 3, 4, 5], 0)

println("\n" * "=" ^ 80)
println("Test 2: Small values with k=2 (b=4)")
println("=" ^ 80)
all_passed &= test_golomb_rice(UInt8[0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10], 2)

println("\n" * "=" ^ 80)
println("Test 3: Small values with k=4 (b=16)")
println("=" ^ 80)
all_passed &= test_golomb_rice(UInt8[0, 5, 10, 15, 20, 25, 30], 4)

println("\n" * "=" ^ 80)
println("Test 4: Larger values with k=6 (b=64)")
println("=" ^ 80)
all_passed &= test_golomb_rice(UInt16[0, 10, 50, 100, 200, 500, 1000], 6)

println("\n" * "=" ^ 80)
println("Test 5: Edge cases")
println("=" ^ 80)
all_passed &= test_golomb_rice(UInt8[0], 3)  # Zero
all_passed &= test_golomb_rice(UInt8[255], 7)  # Max UInt8
all_passed &= test_golomb_rice(UInt16[0, 1, 255, 256, 65535], 8)  # Various sizes

println("\n" * "=" ^ 80)
println("Test 6: Power-of-2 boundaries")
println("=" ^ 80)
# Test values around powers of 2
all_passed &= test_golomb_rice(UInt8[15, 16, 17, 31, 32, 33, 63, 64, 65], 4)

println("\n" * "=" ^ 80)
println("Test 7: Comparison with regular Golomb (k=6 vs b=64)")
println("=" ^ 80)

using Adjacently.Compression: write_golomb, read_golomb

values = UInt16[0, 10, 50, 100, 200, 500, 1000]
k = 6
b = 64  # 2^6

# Golomb-Rice encoding
buffer_rice = IOBuffer()
w_rice = BitWriter(buffer_rice)
for val in values
    write_golomb_rice(w_rice, val, k)
end
flush_bitwriter(w_rice, flush_last_bits=true)
rice_bytes = position(buffer_rice)

# Regular Golomb encoding
buffer_golomb = IOBuffer()
w_golomb = BitWriter(buffer_golomb)
for val in values
    write_golomb(w_golomb, val, b)
end
flush_bitwriter(w_golomb, flush_last_bits=true)
golomb_bytes = position(buffer_golomb)

println("  Values: $values")
println("  Golomb-Rice (k=$k): $rice_bytes bytes")
println("  Golomb (b=$b): $golomb_bytes bytes")
println("  $(rice_bytes == golomb_bytes ? "✅ PASS: Sizes match" : "❌ FAIL: Sizes differ")")

# Verify they decode to the same values
seekstart(buffer_rice)
r_rice = BitReader(buffer_rice)
decoded_rice = [read_golomb_rice(r_rice, k, UInt16) for _ in 1:length(values)]

seekstart(buffer_golomb)
r_golomb = BitReader(buffer_golomb)
decoded_golomb = [read_golomb(r_golomb, b, UInt16) for _ in 1:length(values)]

if decoded_rice == values && decoded_golomb == values
    println("  ✅ PASS: Both decode correctly")
    all_passed &= true
else
    println("  ❌ FAIL: Decoding mismatch")
    all_passed &= false
end

# Final result
println("\n" * "=" ^ 80)
if all_passed
    println("🎉 ALL TESTS PASSED!")
else
    println("❌ SOME TESTS FAILED")
    exit(1)
end
println("=" ^ 80)
