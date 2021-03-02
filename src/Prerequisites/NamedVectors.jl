module NamedVectors

using Printf: @printf
using ..TypeTraits: efficientoperations

export NamedVector, HomoNamedVector

"""
    NamedVector

Abstract type for all named vectors.
"""
abstract type NamedVector end

"""
    getindex(nv::NamedVector, index::Int)

Get the value by the `[]` syntax.
"""
Base.getindex(nv::NamedVector, index::Int) = getfield(nv, index)

"""
    setindex!(nv::NamedVector, value, index::Int)

Set the value by the `[]` syntax if mutable.
"""
Base.setindex!(nv::NamedVector, value, index::Int) = setfield!(nv, index, value)

"""
    ==(nv1::NamedVector, nv2::NamedVector) -> Bool

Overloaded equivalent operator. Two named vector are equal to each other if and only if their keys as well as their values are equal to each other.
!!! note
    It is not necessary for two named vectors to be of the same concrete type to be equal to each other.
"""
Base.:(==)(nv1::NamedVector, nv2::NamedVector) = (keys(nv1) == keys(nv2)) && (values(nv1 ) == values(nv2))

"""
    isequal(nv1::NamedVector, nv2::NamedVector) -> Bool

Overloaded equivalent operator. Two named vector are equal to each other if and only if their keys as well as their values are equal to each other.
!!! note
    It is not necessary for two named vectors to be of the same concrete type to be equal to each other.
"""
Base.isequal(nv1::NamedVector, nv2::NamedVector) = isequal(keys(nv1), keys(nv2)) && isequal(values(nv1), values(nv2))

"""
    <(nv1::NamedVector, nv2::NamedVector) -> Bool

Compare two named vectors and judge whether the first is less than the second.
"""
Base.:<(nv1::NamedVector, nv2::NamedVector) = <(efficientoperations, nv1, nv2)

"""
    isless(nv1::NamedVector, nv2::NamedVector) -> Bool

Compare two named vectors and judge whether the first is less than the second.
"""
Base.isless(nv1::NamedVector, nv2::NamedVector) = isless(efficientoperations, nv1, nv2)

"""
    show(io::IO, nv::NamedVector)

Show a concrete `NamedVector`.
"""
Base.show(io::IO, nv::NamedVector) = @printf io "%s(%s)" nameof(typeof(nv)) join(repr.(values(nv)), ", ")

"""
    hash(nv::NamedVector, h::UInt)

Hash a concrete `NamedVector`.
"""
Base.hash(nv::NamedVector, h::UInt) = hash(values(nv), h)

"""
    convert(::Type{Tuple}, nv::NamedVector) -> Tuple
    convert(::Type{NV}, nv::Tuple) where {NV<:NamedVector} -> NV

Convert a named vector to tuple and vice versa.
"""
@generated Base.convert(::Type{Tuple}, nv::NamedVector) = Expr(:tuple, (:(getfield(nv, $i)) for i = 1:fieldcount(nv))...)
function Base.convert(::Type{NV}, nv::Tuple) where NV<:NamedVector
    @assert fieldcount(NV) == length(nv) "convert error: dismatched length between $NV($(fieldcount(NV))) and input tuple($(length(nv)))."
    return NV(nv...)
end

"""
    length(::Type{NV}) where NV<:NamedVector -> Int
    length(nv::NamedVector) -> Int

Get the length of a concrete `NamedVector`.
"""
Base.length(::Type{NV}) where {NV<:NamedVector} = fieldcount(NV)
Base.length(nv::NamedVector) = length(typeof(nv))

"""
    zero(::Type{NV}) where NV<:NamedVector -> NV
    zero(nv::NamedVector) -> typeof(nv)

Get a concrete `NamedVector` with all values being zero.
"""
@generated function Base.zero(::Type{NV}) where NV<:NamedVector
    zeros = (zero(fieldtype(NV, i)) for i = 1:fieldcount(NV))
    return :(NV($(zeros...)))
end
Base.zero(nv::NamedVector) = zero(typeof(nv))

"""
    iterate(nv::NamedVector, state=1)
    iterate(rv::Iterators.Reverse{<:NamedVector}, state=length(rv.itr))

Iterate or reversely iterate over the values of a concrete `NamedVector`.
"""
Base.iterate(nv::NamedVector, state=1) = (state > length(nv)) ? nothing : (getfield(nv, state), state+1)
Base.iterate(rv::Iterators.Reverse{<:NamedVector}, state=length(rv.itr)) = (state < 1) ? nothing : (getfield(rv.itr, state), state-1)

"""
    keys(nv::NamedVector) -> NTuple{length(nv), Symbol}

Iterate over the names.
"""
Base.keys(nv::NamedVector) = fieldnames(typeof(nv))

"""
    values(nv::NamedVector) -> NTuple{length(nv)}

Iterate over the values.
"""
Base.values(nv::NamedVector) = convert(Tuple, nv)

"""
    pairs(nv::NamedVector) -> Base.Generator

Iterate over the name-value pairs.
"""
Base.pairs(nv::NamedVector) = Base.Generator(=>, keys(nv), values(nv))

"""
    replace(nv::NamedVector; kwargs...) -> typeof(nv)

Return a copy of a concrete `NamedVector` with some of the field values replaced by the keyword arguments.
"""
Base.replace(nv::NamedVector; kwargs...) = replace(efficientoperations, nv; kwargs...)

"""
    map(f, nvs::NV...) where NV<:NamedVector -> NV

Apply function `f` elementwise on the input named vectors.
"""
@generated function Base.map(f, nvs::NV...) where NV<:NamedVector
    exprs = Vector{Expr}(undef, fieldcount(NV))
    for i = 1:length(exprs)
        tmp = [:(nvs[$j][$i]) for j = 1:length(nvs)]
        exprs[i] = :(f($(tmp...)))
    end
    return :(($NV)($(exprs...)))
end

"""
    HomoNamedVector{T}

Abstract type for all homogeneous named vectors.
"""
abstract type HomoNamedVector{T} <: NamedVector end

"""
    eltype(::Type{<:HomoNamedVector{T}}) where T
    eltype(nv::HomoNamedVector)

Get the type parameter of a concrete `HomoNamedVector`.
"""
Base.eltype(::Type{<:HomoNamedVector{T}}) where T = T
Base.eltype(nv::HomoNamedVector) = eltype(typeof(nv))

end #module
