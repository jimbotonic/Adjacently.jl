# Decoder-only containers for consecutive prefixes; sparse research workspaces
# and whole-graph decoders retain ordinary Dict/Vector storage.
mutable struct _PrefixLists{T<:Unsigned} <: AbstractDict{T,Vector{T}}
    first_vertex::Int
    lists::Vector{Vector{T}}
end
_PrefixLists{T}() where T = _PrefixLists{T}(0, Vector{T}[])
Base.length(p::_PrefixLists) = length(p.lists)
Base.haskey(p::_PrefixLists, v) = 1 <= Int(v)-p.first_vertex+1 <= length(p.lists)
function Base.getindex(p::_PrefixLists, v)
    haskey(p,v) || throw(KeyError(v))
    return p.lists[Int(v)-p.first_vertex+1]
end
function Base.setindex!(p::_PrefixLists{T}, ns::Vector{T}, v::T) where T
    isempty(p.lists) && (p.first_vertex=Int(v))
    i=Int(v)-p.first_vertex+1
    if i == length(p.lists)+1
        push!(p.lists,ns)
    elseif 1 <= i <= length(p.lists)
        p.lists[i]=ns
    else
        throw(ArgumentError("prefix lists must be appended consecutively"))
    end
    return p
end
function Base.iterate(p::_PrefixLists{T}, i::Int=1) where T
    i > length(p.lists) && return nothing
    return (T(p.first_vertex+i-1)=>p.lists[i], i+1)
end
Base.empty!(p::_PrefixLists) = (empty!(p.lists); p.first_vertex=0; p)

mutable struct _ReferenceRing{T} <: AbstractVector{T}
    data::Vector{T}
    head::Int
    count::Int
end
_ReferenceRing{T}() where T = _ReferenceRing{T}(T[],1,0)
Base.size(r::_ReferenceRing) = (r.count,)
Base.IndexStyle(::Type{<:_ReferenceRing}) = IndexLinear()
@inline function Base.getindex(r::_ReferenceRing, i::Int)
    @boundscheck checkbounds(r,i)
    j=r.head+i-1
    j > length(r.data) && (j-=length(r.data))
    return @inbounds r.data[j]
end
function Base.push!(r::_ReferenceRing{T}, x::T) where T
    if r.count == length(r.data)
        data=Vector{T}(undef,max(8,2length(r.data)))
        for i in 1:r.count; data[i]=r[i]; end
        r.data=data; r.head=1
    end
    j=r.head+r.count
    j > length(r.data) && (j-=length(r.data))
    r.data[j]=x; r.count+=1
    return r
end
function Base.popfirst!(r::_ReferenceRing)
    isempty(r) && throw(ArgumentError("empty reference ring"))
    x=r[1]; r.head=r.head == length(r.data) ? 1 : r.head+1; r.count-=1
    return x
end
Base.empty!(r::_ReferenceRing) = (r.head=1; r.count=0; r)
