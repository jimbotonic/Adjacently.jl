#
# Adjacently: Julia Complex Directed Networks Library
# Copyright (C) 2016-2025 Jimmy Dubuisson <jimmy.dubuisson@gmail.com>
#
# This program is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.
#
# This program is distributed in the hope that it should be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.
#

module CustomTypes

# Export the custom types
export UInt24, UInt40

# Trailing zeros
import Base: trailing_zeros, convert, reinterpret, >>, <<, Float64, Float32, 
    lastindex, firstindex, promote_rule, convert, AbstractFloat, Int64, UInt32, UInt64

# All the existing type definitions and implementations
primitive type UInt24 <: Unsigned 24 end
primitive type UInt40 <: Unsigned 40 end

# Implement convert from Int64 to UInt24
function Base.convert(::Type{UInt24}, x::Int64)
    if x < 0 || x > 0xFFFFFF
        throw(InexactError(:convert, UInt24, x))
    end
    return reinterpret(UInt24, UInt32(x) & 0xFFFFFF)
end

# Implement convert from Int64 to UInt40
function Base.convert(::Type{UInt40}, x::Int64)
    if x < 0 || x > 0xFFFFFFFFFF
        throw(InexactError(:convert, UInt40, x))
    end
    return reinterpret(UInt40, UInt64(x) & 0xFFFFFFFFFF)
end

# STD -> Custom conversion
function Base.reinterpret(::Type{UInt24}, x::UInt32)
    if x > 0xFFFFFF
        throw(InexactError(:reinterpret, UInt24, x))
    end
    return unsafe_load(Ptr{UInt24}(pointer_from_objref(Ref(x & 0xFFFFFF))))
end

function Base.reinterpret(::Type{UInt40}, x::UInt64)
    if x > 0xFFFFFFFFFF
        throw(InexactError(:reinterpret, UInt40, x))
    end
    return unsafe_load(Ptr{UInt40}(pointer_from_objref(Ref(x & 0xFFFFFFFFFF))))
end

# Custom -> STD conversion
function Base.reinterpret(::Type{UInt32}, x::UInt24)
    return unsafe_load(Ptr{UInt32}(pointer_from_objref(Ref(x))))
end

function Base.reinterpret(::Type{UInt64}, x::UInt40)
    return unsafe_load(Ptr{UInt64}(pointer_from_objref(Ref(x))))
end

# Basic arithmetic constants
Base.zero(::Type{UInt24}) = reinterpret(UInt24, 0x000000)
Base.one(::Type{UInt24}) = reinterpret(UInt24, 0x000001)
Base.zero(::Type{UInt40}) = reinterpret(UInt40, 0x0000000000)
Base.one(::Type{UInt40}) = reinterpret(UInt40, 0x0000000001)

# Comparison operators
Base.:<(x::UInt24, y::UInt24) = reinterpret(UInt32, x) < reinterpret(UInt32, y)
Base.:<=(x::UInt24, y::UInt24) = reinterpret(UInt32, x) <= reinterpret(UInt32, y)
Base.:>(x::UInt24, y::UInt24) = reinterpret(UInt32, x) > reinterpret(UInt32, y)
Base.:>=(x::UInt24, y::UInt24) = reinterpret(UInt32, x) >= reinterpret(UInt32, y)

Base.:<(x::UInt40, y::UInt40) = reinterpret(UInt64, x) < reinterpret(UInt64, y)
Base.:<=(x::UInt40, y::UInt40) = reinterpret(UInt64, x) <= reinterpret(UInt64, y)
Base.:>(x::UInt40, y::UInt40) = reinterpret(UInt64, x) > reinterpret(UInt64, y)
Base.:>=(x::UInt40, y::UInt40) = reinterpret(UInt64, x) >= reinterpret(UInt64, y)

# Display and utility functions
Base.show(io::IO, x::UInt24) = print(io, "UInt24(", reinterpret(UInt32, x), ")")
Base.show(io::IO, x::UInt40) = print(io, "UInt40(", reinterpret(UInt64, x), ")")
Base.typemax(::Type{UInt24}) = reinterpret(UInt24, 0xFFFFFF)
Base.typemax(::Type{UInt40}) = reinterpret(UInt40, 0xFFFFFFFFFF)

# Implement arithmetic operations
Base.:+(x::UInt24, y::UInt24) = reinterpret(UInt24, (reinterpret(UInt32, x) + reinterpret(UInt32, y)) & 0xFFFFFF)
Base.:-(x::UInt24, y::UInt24) = reinterpret(UInt24, (reinterpret(UInt32, x) - reinterpret(UInt32, y)) & 0xFFFFFF)

Base.:+(x::UInt40, y::UInt40) = reinterpret(UInt40, (reinterpret(UInt64, x) + reinterpret(UInt64, y)) & 0xFFFFFFFFFF)
Base.:-(x::UInt40, y::UInt40) = reinterpret(UInt40, (reinterpret(UInt64, x) - reinterpret(UInt64, y)) & 0xFFFFFFFFFF)

function >>(x::UInt24, n::Int)
    result = reinterpret(UInt32, x) >> n
    return reinterpret(UInt24, result & 0xFFFFFF)
end

function <<(x::UInt24, n::Int)
    result = reinterpret(UInt32, x) << n
    return reinterpret(UInt24, result & 0xFFFFFF)
end

function >>(x::UInt40, n::Int)
    result = reinterpret(UInt64, x) >> n
    return reinterpret(UInt40, result & 0xFFFFFFFFFF)
end

function <<(x::UInt40, n::Int)
    result = reinterpret(UInt64, x) << n
    return reinterpret(UInt40, result & 0xFFFFFFFFFF)
end

# Float conversions
Base.Float64(x::UInt24) = Float64(reinterpret(UInt32, x))
Base.Float32(x::UInt24) = Float32(reinterpret(UInt32, x))
Base.Float64(x::UInt40) = Float64(reinterpret(UInt64, x))
Base.Float32(x::UInt40) = Float32(reinterpret(UInt64, x))

# Dictionary indexing
function lastindex(d::Dict{UInt24, T}) where T
    return convert(UInt24, length(d))
end

function lastindex(d::Dict{UInt40, T}) where T
    return convert(UInt40, length(d))
end

function firstindex(d::Dict{UInt24, T}) where T
    return one(UInt24)
end

function firstindex(d::Dict{UInt40, T}) where T
    return one(UInt40)
end

# Additional conversions
Base.convert(::Type{UInt24}, x::UInt8) = reinterpret(UInt24, UInt32(x))
Base.convert(::Type{UInt24}, x::UInt16) = reinterpret(UInt24, UInt32(x))
Base.convert(::Type{UInt40}, x::UInt8) = reinterpret(UInt40, UInt64(x))
Base.convert(::Type{UInt40}, x::UInt16) = reinterpret(UInt40, UInt64(x))
Base.convert(::Type{UInt40}, x::UInt32) = reinterpret(UInt40, UInt64(x))

function trailing_zeros(x::UInt24)
    # Convert to UInt32 and mask the relevant bits
    return trailing_zeros(reinterpret(UInt32, x) & 0xFFFFFF)
end

function trailing_zeros(x::UInt40)
    # Convert to UInt64 and mask the relevant bits
    return trailing_zeros(reinterpret(UInt64, x) & 0xFFFFFFFFFF)
end

# Implement right shift for UInt24
function >>(x::UInt24, n::Unsigned)
    result = reinterpret(UInt32, x) >> n
    return reinterpret(UInt24, result & 0xFFFFFF)
end

# Implement left shift for UInt24
function <<(x::UInt24, n::Unsigned)
    result = reinterpret(UInt32, x) << n
    return reinterpret(UInt24, result & 0xFFFFFF)
end

# Implement right shift for UInt40
function >>(x::UInt40, n::Unsigned)
    result = reinterpret(UInt64, x) >> n
    return reinterpret(UInt40, result & 0xFFFFFFFFFF)
end

# Implement left shift for UInt40
function <<(x::UInt40, n::Unsigned)
    result = reinterpret(UInt64, x) << n
    return reinterpret(UInt40, result & 0xFFFFFFFFFF)
end

# Also implement versions that accept Integer for the shift amount
>>(x::UInt24, n::Integer) = >>(x, convert(UInt, n))
<<(x::UInt24, n::Integer) = <<(x, convert(UInt, n))
>>(x::UInt40, n::Integer) = >>(x, convert(UInt, n))
<<(x::UInt40, n::Integer) = <<(x, convert(UInt, n))

# Implement AbstractFloat conversion for UInt24 and UInt40
function AbstractFloat(x::UInt24)
    return AbstractFloat(reinterpret(UInt32, x))
end

function AbstractFloat(x::UInt40)
    return AbstractFloat(reinterpret(UInt64, x))
end

# Implement promotion rules for UInt24
promote_rule(::Type{UInt24}, ::Type{UInt8}) = UInt24
promote_rule(::Type{UInt24}, ::Type{UInt16}) = UInt24
promote_rule(::Type{UInt24}, ::Type{UInt32}) = UInt32
promote_rule(::Type{UInt24}, ::Type{UInt64}) = UInt64

# Implement promotion rules for UInt40
promote_rule(::Type{UInt40}, ::Type{UInt8}) = UInt40
promote_rule(::Type{UInt40}, ::Type{UInt16}) = UInt40
promote_rule(::Type{UInt40}, ::Type{UInt32}) = UInt40
promote_rule(::Type{UInt40}, ::Type{UInt64}) = UInt64

# Add conversion methods if needed
convert(::Type{UInt64}, x::UInt24) = UInt64(reinterpret(UInt32, x))
convert(::Type{UInt64}, x::UInt40) = UInt64(reinterpret(UInt64, x))

# Add promotion rules for Int64 and other integer types
Base.promote_rule(::Type{UInt24}, ::Type{Int64}) = Int64
Base.promote_rule(::Type{UInt24}, ::Type{Int32}) = Int64
Base.promote_rule(::Type{UInt40}, ::Type{Int32}) = Int64
Base.promote_rule(::Type{UInt40}, ::Type{Int64}) = Int64

# Add reverse conversion methods
Base.convert(::Type{Int64}, x::UInt24) = Int64(reinterpret(UInt32, x))
Base.convert(::Type{Int64}, x::UInt40) = Int64(reinterpret(UInt64, x))

# Add conversion methods for Int64 to standard unsigned types
function Base.convert(::Type{UInt32}, x::Int64)
    if x < 0 || x > typemax(UInt32)
        throw(InexactError(:convert, UInt32, x))
    end
    return UInt32(x)
end

function Base.convert(::Type{UInt64}, x::Int64)
    if x < 0
        throw(InexactError(:convert, UInt64, x))
    end
    return UInt64(x)
end


##########

# Add constructor methods for UInt24 and UInt40
function (::Type{UInt24})(x::UInt32)
    if x > 0xFFFFFF
        throw(InexactError(:UInt24, UInt24, x))
    end
    return reinterpret(UInt24, x & 0xFFFFFF)
end

function (::Type{UInt40})(x::UInt64)
    if x > 0xFFFFFFFFFF
        throw(InexactError(:UInt40, UInt40, x))
    end
    return reinterpret(UInt40, x & 0xFFFFFFFFFF)
end

function Int64(x::UInt24)
    return Int64(reinterpret(UInt32, x))
end

function Int64(x::UInt40)
    return Int64(reinterpret(UInt64, x))
end

function UInt32(x::UInt24)
    return reinterpret(UInt32, x)
end

# Add UInt64 constructor for UInt40 while we're at it
function UInt64(x::UInt40)
    return reinterpret(UInt64, x)
end

function (::Type{UInt24})(x::UInt64)
    if x > 0xFFFFFF
        throw(InexactError(:UInt24, UInt24, x))
    end
    return reinterpret(UInt24, UInt32(x) & 0xFFFFFF)
end

#####

end # module CustomTypes