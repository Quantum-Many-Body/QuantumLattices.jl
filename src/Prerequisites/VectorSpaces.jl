module VectorSpaces

using ..Traits: rawtype, efficientoperations, getcontent

import ..Traits: contentnames
import ...Interfaces: dimension, rank, ⊕, ⊗

export VectorSpace, NamedVectorSpace, SimpleNamedVectorSpace, ParameterSpace, ZippedNamedVectorSpace, DirectProductedNamedVectorSpace
export VectorSpaceStyle, VectorSpaceEnumerative, VectorSpaceCartesian, VectorSpaceDirectSummed, VectorSpaceDirectProducted, VectorSpaceZipped
export shape

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
    contents = getcontent(vs, :contents)
    dimsum = cumsum(map(dimension, contents))
    m = searchsortedfirst(dimsum, i)
    n = i-(m>1 ? dimsum[m-1] : 0)
    return contents[m][n]
end

"""
    VectorSpaceDirectProducted <: VectorSpaceStyle

Vector space style which indicates that a vector space is the direct product of its sub-components.
"""
struct VectorSpaceDirectProducted <: VectorSpaceStyle end
@inline dimension(::VectorSpaceDirectProducted, vs::VectorSpace) = mapreduce(dimension, *, getcontent(vs, :contents))
@inline function Base.getindex(::VectorSpaceDirectProducted, vs::VectorSpace, i::Integer)
    contents = getcontent(vs, :contents)
    index = CartesianIndices(map(dimension, contents))[i]
    return rawtype(eltype(vs))(map(getindex, contents, Tuple(index)), vs)
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
    NamedVectorSpace{B} <: VectorSpace{B}

Abstract named vector space.
"""
abstract type NamedVectorSpace{B} <: VectorSpace{B} end
@inline Base.keys(vs::NamedVectorSpace) = keys(typeof(vs))
@inline Base.:(==)(vs₁::NamedVectorSpace, vs₂::NamedVectorSpace) = keys(vs₁)==keys(vs₂) && invoke(==, Tuple{VectorSpace, VectorSpace}, vs₁, vs₂)
@inline Base.isequal(vs₁::NamedVectorSpace, vs₂::NamedVectorSpace) = isequal(keys(vs₁), keys(vs₂)) && invoke(isequal, Tuple{VectorSpace, VectorSpace}, vs₁, vs₂)
@inline Base.pairs(vs::NamedVectorSpace) = NamedVectorSpacePairIteration(vs)
struct NamedVectorSpacePairIteration{VS<:NamedVectorSpace} <: AbstractVector{NamedTuple}
    vs::VS
end
@inline Base.size(ps::NamedVectorSpacePairIteration) = size(ps.vs)
@inline Base.eltype(::Type{<:NamedVectorSpacePairIteration{VS}}) where {VS<:NamedVectorSpace} = NamedTuple{keys(VS), eltype(VS)}
@inline Base.getindex(ps::NamedVectorSpacePairIteration{VS}, i::Integer) where {VS<:NamedVectorSpace} = NamedTuple{keys(VS)}(ps.vs[i])

"""
    SimpleNamedVectorSpace{N, B} <: NamedVectorSpace{B}

Abstract simple named vector space.

The subtypes must have the following predefined content:
* `:content::T`: the content of the named vector space.
"""
abstract type SimpleNamedVectorSpace{N, B} <: NamedVectorSpace{B} end
@inline contentnames(::Type{<:SimpleNamedVectorSpace}) = (:content,)
@inline dimension(vs::SimpleNamedVectorSpace) = length(getcontent(vs, :content))
@inline Base.getindex(vs::SimpleNamedVectorSpace, i) = getcontent(vs, :content)[i]
@inline Base.in(basis::B, vs::SimpleNamedVectorSpace{N, B}) where {N, B} = in(basis, getcontent(vs, :content))
@inline Base.keys(::Type{<:SimpleNamedVectorSpace{N}}) where N = (N,)
@inline Base.eltype(::Type{<:NamedVectorSpacePairIteration{VS}}) where {VS<:SimpleNamedVectorSpace} = NamedTuple{keys(VS), Tuple{eltype(VS)}}
@inline Base.getindex(ps::NamedVectorSpacePairIteration{VS}, i::Integer) where {VS<:SimpleNamedVectorSpace} = NamedTuple{keys(VS)}((ps.vs[i],))

"""
    ParameterSpace{N, T, B} <: SimpleNamedVectorSpace{N, B}

Parameter space.
"""
struct ParameterSpace{N, T, B} <: SimpleNamedVectorSpace{N, B}
    content::T
    function ParameterSpace{N}(content) where N
        @assert isa(N, Symbol) "ParameterSpace error: $(N) is not a Symbol."
        new{N, typeof(content), eltype(content)}(content)
    end
end

"""
    ZippedNamedVectorSpace{T<:Tuple{Vararg{NamedVectorSpace}}} <: NamedVectorSpace{Any}

Zipped named vector space.
"""
struct ZippedNamedVectorSpace{T<:Tuple{Vararg{NamedVectorSpace}}} <: NamedVectorSpace{Any}
    contents::T
end
@inline ZippedNamedVectorSpace(contents::NamedVectorSpace...) = ZippedNamedVectorSpace(contents)
@inline VectorSpaceStyle(::Type{<:ZippedNamedVectorSpace}) = VectorSpaceZipped()
@inline @generated Base.eltype(::Type{<:ZippedNamedVectorSpace{T}}) where {T<:Tuple{Vararg{NamedVectorSpace}}} = Tuple{map(eltype, fieldtypes(T))...}
@inline @generated Base.keys(::Type{<:ZippedNamedVectorSpace{T}}) where {T<:Tuple{Vararg{NamedVectorSpace}}} = map(first, map(keys, fieldtypes(T)))
@inline Tuple(basis::Tuple, ::ZippedNamedVectorSpace) = basis

"""
    DirectProductedNamedVectorSpace{T<:Tuple{Vararg{NamedVectorSpace}}} <: NamedVectorSpace{Any}

Zipped named vector space.
"""
struct DirectProductedNamedVectorSpace{T<:Tuple{Vararg{NamedVectorSpace}}} <: NamedVectorSpace{Any}
    contents::T
end
@inline DirectProductedNamedVectorSpace(contents::NamedVectorSpace...) = DirectProductedNamedVectorSpace(contents)
@inline VectorSpaceStyle(::Type{<:DirectProductedNamedVectorSpace}) = VectorSpaceDirectProducted()
@inline @generated Base.eltype(::Type{<:DirectProductedNamedVectorSpace{T}}) where {T<:Tuple{Vararg{NamedVectorSpace}}} = Tuple{map(eltype, fieldtypes(T))...}
@inline @generated Base.keys(::Type{<:DirectProductedNamedVectorSpace{T}}) where {T<:Tuple{Vararg{NamedVectorSpace}}} = map(first, map(keys, fieldtypes(T)))
@inline Tuple(basis::Tuple, ::DirectProductedNamedVectorSpace) = basis

"""
    ⊕(vs₁::SimpleNamedVectorSpace, vs₂::SimpleNamedVectorSpace) -> ZipNamedVectorSpace
    ⊕(vs₁::SimpleNamedVectorSpace, vs₂::ZippedNamedVectorSpace) -> ZipNamedVectorSpace
    ⊕(vs₁::ZippedNamedVectorSpace, vs₂::SimpleNamedVectorSpace) -> ZipNamedVectorSpace
    ⊕(vs₁::ZippedNamedVectorSpace, vs₂::ZippedNamedVectorSpace) -> ZipNamedVectorSpace

The zip of named vector spaces.
"""
@inline ⊕(vs₁::SimpleNamedVectorSpace, vs₂::SimpleNamedVectorSpace) = ZippedNamedVectorSpace(vs₁, vs₂)
@inline ⊕(vs₁::SimpleNamedVectorSpace, vs₂::ZippedNamedVectorSpace) = ZippedNamedVectorSpace(vs₁, vs₂.contents...)
@inline ⊕(vs₁::ZippedNamedVectorSpace, vs₂::SimpleNamedVectorSpace) = ZippedNamedVectorSpace(vs₁.contents..., vs₂)
@inline ⊕(vs₁::ZippedNamedVectorSpace, vs₂::ZippedNamedVectorSpace) = ZippedNamedVectorSpace(vs₁.contents..., vs₂.contents...)

"""
    ⊗(vs₁::SimpleNamedVectorSpace, vs₂::SimpleNamedVectorSpace) -> DirectProductedNamedVectorSpace
    ⊗(vs₁::SimpleNamedVectorSpace, vs₂::DirectProductedNamedVectorSpace) -> DirectProductedNamedVectorSpace
    ⊗(vs₁::DirectProductedNamedVectorSpace, vs₂::SimpleNamedVectorSpace) -> DirectProductedNamedVectorSpace
    ⊗(vs₁::DirectProductedNamedVectorSpace, vs₂::DirectProductedNamedVectorSpace) -> DirectProductedNamedVectorSpace

The direct product of named vector spaces.
"""
@inline ⊗(vs₁::SimpleNamedVectorSpace, vs₂::SimpleNamedVectorSpace) = DirectProductedNamedVectorSpace(vs₁, vs₂)
@inline ⊗(vs₁::SimpleNamedVectorSpace, vs₂::DirectProductedNamedVectorSpace) = DirectProductedNamedVectorSpace(vs₁, vs₂.contents...)
@inline ⊗(vs₁::DirectProductedNamedVectorSpace, vs₂::SimpleNamedVectorSpace) = DirectProductedNamedVectorSpace(vs₁.contents..., vs₂)
@inline ⊗(vs₁::DirectProductedNamedVectorSpace, vs₂::DirectProductedNamedVectorSpace) = DirectProductedNamedVectorSpace(vs₁.contents..., vs₂.contents...)

end # module
