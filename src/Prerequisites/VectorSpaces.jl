module VectorSpaces

using ..Traits: rawtype, efficientoperations, getcontent

import ..Traits: contentnames
import ...Interfaces: dimension, rank

export VectorSpace, NamedVectorSpace
export VectorSpaceStyle, VectorSpaceEnumerative, VectorSpaceCartesian, VectorSpaceDirectSummed, VectorSpaceDirectProducted, VectorSpaceZipped
# NamedVectorSpace, ZipNamedVectorSpace, ProductNamedVectorSpace, shape

"""
    VectorSpace{B} <: AbstractVector{B}

Abstract vector space.
"""
abstract type VectorSpace{B} <: AbstractVector{B} end
@inline Base.:(==)(vs₁::VectorSpace, vs₂::VectorSpace) = ==(efficientoperations, vs₁, vs₂)
@inline Base.isequal(vs₁::VectorSpace, vs₂::VectorSpace) = isequal(efficientoperations, vs₁, vs₂)
@inline Base.size(vs::VectorSpace) = (dimension(vs),)

"""
    VectorSpaceStyle

The style of a concrete type of vector space.
"""
abstract type VectorSpaceStyle end
@inline VectorSpaceStyle(vs::VectorSpace) = VectorSpaceStyle(typeof(vs))

@inline dimension(vs::VectorSpace) = dimension(VectorSpaceStyle(vs), vs)
@inline Base.getindex(vs::VectorSpace, i) = getindex(VectorSpaceStyle(vs), vs, i)
@inline Base.issorted(vs::VectorSpace) = issorted(VectorSpaceStyle(vs), vs)
@inline Base.findfirst(basis::B, vs::VectorSpace{B}) where B = findfirst(VectorSpaceStyle(vs), basis, vs)
@inline Base.searchsortedfirst(vs::VectorSpace{B}, basis::B) where B = searchsortedfirst(VectorSpaceStyle(vs), vs, basis)
@inline Base.in(basis::B, vs::VectorSpace{B}) where B = in(VectorSpaceStyle(vs), basis, vs)

@inline Base.issorted(::VectorSpaceStyle, vs::VectorSpace) = false
@inline function Base.findfirst(::VectorSpaceStyle, basis::B, vs::VectorSpace{B}) where B
    if issorted(vs)
        index = invoke(searchsortedfirst, Tuple{Any, Any}, vs, basis)
        return (0<index<=length(vs) && vs[index]==basis) ? index : nothing
    end
    for i = 1:length(vs)
        vs[i]==basis && return i
    end
end
@inline function Base.searchsortedfirst(::VectorSpaceStyle, vs::VectorSpace{B}, basis::B) where B
    issorted(vs) && return invoke(searchsortedfirst, Tuple{Any, Any}, vs, basis)
    for i = 1:length(vs) 
        vs[i]==basis && return i
    end
    return length(vs)+1
end
@inline Base.in(::VectorSpaceStyle, basis::B, vs::VectorSpace{B}) where B = isa(findfirst(basis, vs), Integer)

"""
    VectorSpaceEnumerative <: VectorSpaceStyle

Enumerative vector space style, which indicates that the vector space has a predefined content named `:table` that contains all its bases.
"""
struct VectorSpaceEnumerative <: VectorSpaceStyle end
@inline dimension(::VectorSpaceEnumerative, vs::VectorSpace) = length(getcontent(vs, :table))
@inline Base.getindex(::VectorSpaceEnumerative, vs::VectorSpace, i::Int) = getcontent(vs, :table)[i]

"""
    VectorSpaceCartesian <: VectorSpaceStyle

Cartesian vector space style, which indicates that every basis in it could be represented by a Cartesian index.
"""
struct VectorSpaceCartesian <: VectorSpaceStyle end
function shape end
@inline dimension(::VectorSpaceCartesian, vs::VectorSpace) = mapreduce(length, *, shape(vs))
@inline Base.getindex(::VectorSpaceCartesian, vs::VectorSpace, i::Integer) = getindex(VectorSpaceCartesian(), vs, CartesianIndices(shape(vs))[i])
@inline Base.getindex(::VectorSpaceCartesian, vs::VectorSpace, i::CartesianIndex) = rawtype(eltype(vs))(i, vs)
@inline Base.issorted(::VectorSpaceCartesian, vs::VectorSpace) = true
@inline function Base.findfirst(::VectorSpaceCartesian, basis::B, vs::VectorSpace{B}) where B
    index, idx = CartesianIndex(basis, vs), CartesianIndices(shape(vs))
    index∉idx && return nothing
    return LinearIndices(shape(vs))[index-first(idx)+oneunit(typeof(index))]
end
@inline function Base.searchsortedfirst(::VectorSpaceCartesian, vs::VectorSpace{B}, basis::B) where B
    index, idx = CartesianIndex(basis, vs), CartesianIndices(shape(vs))
    index∈idx && return LinearIndices(shape(vs))[index-first(idx)+oneunit(typeof(index))]
    return searchsortedfirst(idx, index)
end
@inline Base.in(::VectorSpaceCartesian, basis::B, vs::VectorSpace{B}) where B = CartesianIndex(basis, vs)∈CartesianIndices(shape(vs))

"""
    VectorSpaceDirectSummed <: VectorSpaceStyle

Vector space style which indicates that a vector space is the direct sum of its sub-components.
"""
struct VectorSpaceDirectSummed <: VectorSpaceStyle end
@inline dimension(::VectorSpaceDirectSummed, vs::VectorSpace) = mapreduce(dimension, +, getcontent(vs, :contents))
@inline function Base.getindex(::VectorSpaceDirectSummed, vs::VectorSpace, i::Integer)
    dimsum = cumsum(map(dimension, getcontent(vs, :contents)))
    m = searchsortedfirst(dimsum, i)
    n = i-dimsum[m]
    return getcontent(vs, :contents)[m][n]
end

"""
    VectorSpaceDirectProducted <: VectorSpaceStyle

Vector space style which indicates that a vector space is the direct product of its sub-components.
"""
struct VectorSpaceDirectProducted <: VectorSpaceStyle end
@inline dimension(::VectorSpaceDirectProducted, vs::VectorSpace) = mapreduce(dimension, *, getcontent(vs, :contents))
@inline function Base.getindex(::VectorSpaceDirectProducted, vs::VectorSpace, i::Integer)
    index = CartesianIndices(map(dimension, getcontent(vs, :contents)))[i]
    return rawtype(eltype(vs))(map(getindex, getcontent(vs, :contents), Tuple(index)), vs)
end

"""
    VectorSpaceZipped <: VectorSpaceStyle

Vector space style which indicates that a vector space is the zip of its sub-components.
"""
struct VectorSpaceZipped <: VectorSpaceStyle end
@inline dimension(::VectorSpaceZipped, vs::VectorSpace) = mapreduce(dimension, min, getcontent(vs, :contents))
@inline function Base.getindex(::VectorSpaceZipped, vs::VectorSpace, i::Integer)
    return rawtype(eltype(vs))(map(content->getindex(content, i), getcontent(vs, :contents)), vs)
end

"""
    NamedVectorSpace{N, T} <: VectorSpace{NamedTuple}

Abstract named vector space.

The subtypes must have the following predefined content:
* `:content::T`: the content of the vector space.
"""
abstract type NamedVectorSpace{N, T} <: VectorSpace{NamedTuple} end
@inline contentnames(::Type{<:NamedVectorSpace}) = (:content,)
@inline Base.eltype(::Type{<:NamedVectorSpace{N, T}}) where {N, T} = NamedTuple{(N,), Tuple{eltype(T)}}
@inline Base.:(==)(vs₁::NamedVectorSpace{N₁}, vs₂::NamedVectorSpace{N₂}) where {N₁, N₂} = N₁==N₂ && ==(efficientoperations, vs₁, vs₂)
@inline Base.isequal(vs₁::NamedVectorSpace{N₁}, vs₂::NamedVectorSpace{N₂}) where {N₁, N₂} = isequal(N₁, N₂) && isequal(efficientoperations, vs₁, vs₂)
@inline dimension(vs::NamedVectorSpace) = dimension(getcontent(vs, :content))
@inline Base.getindex(vs::NamedVectorSpace{N}, i) where N = NamedTuple{(N,)}((getcontent(vs, :content)[i],))

"""
    ZippedNamedVectorSpace{T<:Tuple{Vararg{NamedVectorSpace}}} <: VectorSpace{NamedTuple}

"""
struct ZippedNamedVectorSpace{T<:Tuple{Vararg{NamedVectorSpace}}} <: VectorSpace{NamedTuple}
    contents::T
end
@inline ZippedNamedVectorSpace(contents::NamedVectorSpace...) = ZippedNamedVectorSpace(contents)
@inline VectorSpaceStyle(::Type{<:ZippedNamedVectorSpace}) = VectorSpaceZipped()
@inline @generated function Base.eltype(::Type{<:ZippedNamedVectorSpace{T}}) where {T<:Tuple{Vararg{NamedVectorSpace}}}

end
# @inline function Base.CartesianIndex(basis::NamedTuple{NS}, nvs::NamedVectorSpace{:⊕, NS}) where NS
#     index = findfirst(isequal(basis[1]), getcontent(nvs, :contents)[1])::Int
#     for i = 2:length(NS)
#         @assert getcontent(nvs, :contents)[i][index] == basis[i] "CartesianIndex error: input basis out of range."
#     end
#     return CartesianIndex(index)
# end
# @inline function Base.CartesianIndex(basis::NamedTuple{NS, BS}, nvs::NamedVectorSpace{:⊗, NS, BS}) where {NS, BS<:Tuple}
#     return CartesianIndex(map((b, content)->findfirst(isequal(b), content), values(basis), getcontent(nvs, :contents)))
# end
# @inline function Base.NamedTuple(index::CartesianIndex{1}, nvs::NamedVectorSpace{:⊕, NS}) where NS
#     return NamedTuple{NS}(map(content->content[index[1]], getcontent(nvs, :contents)))
# end
# @inline function Base.NamedTuple(index::CartesianIndex{N}, nvs::NamedVectorSpace{:⊗, NS}) where {N, NS}
#     @assert N == length(NS) "NamedTuple error: mismatched input index and rank of named vector space."
#     return NamedTuple{NS}(map((i, content)->content[index[i]], ntuple(i->i, Val(N)), getcontent(nvs, :contents)))
# end

# """
#     keys(nvs::NamedVectorSpace) -> Tuple{Vararg{Symbol}}
#     keys(::Type{<:NamedVectorSpace{M, NS} where M}) where NS -> Tuple{Vararg{Symbol}}

# Get the names of a named vector space.
# """
# @inline Base.keys(nvs::NamedVectorSpace) = keys(typeof(nvs))
# @inline Base.keys(::Type{<:NamedVectorSpace{M, NS} where M}) where {NS} = NS

# """
#     values(nvs::NamedVectorSpace) -> Tuple{Vararg{AbstractVector}}

# Get the contents of a named vector space.
# """
# @inline Base.values(nvs::NamedVectorSpace) = getcontent(nvs, :contents)

# """
#     pairs(nvs::NamedVectorSpace) -> Base.Iterators.Pairs

# Get the name-value pairs of a named vector space.
# """
# @inline Base.pairs(nvs::NamedVectorSpace) = Base.Generator(=>, nvs|>keys, nvs|>values)

# """
#     eltype(nvs::NamedVectorSpace, i::Int)
#     eltype(::Type{<:NamedVectorSpace{M, NS, BS} where {M, NS}}, i::Int) where BS

# Get the eltype of the ith indexable dimensions of a named vector space.
# """
# @inline Base.eltype(nvs::NamedVectorSpace, i::Int) = eltype(nvs|>typeof, i)
# @inline Base.eltype(::Type{<:NamedVectorSpace{M, NS, BS} where {M, NS}}, i::Int) where BS = fieldtype(BS, i)

# """
#     ZipNamedVectorSpace{NS, BS<:Tuple, CS<:Tuple} <: NamedVectorSpace{:⊕, NS, BS}

# Zipped named vector space.
# """
# struct ZipNamedVectorSpace{NS, BS<:Tuple, CS<:Tuple} <: NamedVectorSpace{:⊕, NS, BS}
#     contents::CS
# end

# """
#     ZipNamedVectorSpace{NS}(contents...) where NS

# Construct a zipped named vector space.
# """
# function ZipNamedVectorSpace{NS}(contents...) where NS
#     @assert length(NS)==length(contents) && isa(NS, Tuple{Vararg{Symbol}}) "ZipNamedVectorSpace error: mismatched names and contents."
#     @assert mapreduce(length, ==, contents) "ZipNamedVectorSpace error: mismatched length of contents."
#     return ZipNamedVectorSpace{NS, Tuple{map(eltype, contents)...}, typeof(contents)}(contents)
# end

# """
#     ProductNamedVectorSpace{NS, BS<:Tuple, CS<:Tuple} <: NamedVectorSpace{:⊕, NS, BS}

# Product named vector space.
# """
# struct ProductNamedVectorSpace{NS, BS<:Tuple, CS<:Tuple} <: NamedVectorSpace{:⊗, NS, BS}
#     contents::CS
# end

# """
#     ProductNamedVectorSpace{NS}(contents...) where NS

# Construct a product named vector space.
# """
# function ProductNamedVectorSpace{NS}(contents...) where NS
#     @assert length(NS)==length(contents) && isa(NS, Tuple{Vararg{Symbol}}) "ProductNamedVectorSpace error: mismatched names and contents."
#     return ProductNamedVectorSpace{NS, Tuple{map(eltype, contents)...}, typeof(contents)}(contents)
# end

end # module
