module NamedVectors

using Printf: @printf
using ..Traits: efficientoperations

export HomoNamedVector, NamedVector

"""
    NamedVector

Abstract type for all named vectors.
"""
abstract type NamedVector end
@inline Base.:(==)(nv₁::NamedVector, nv₂::NamedVector) = keys(nv₁)==keys(nv₂) && values(nv₁)==values(nv₂)
@inline Base.isequal(nv₁::NamedVector, nv₂::NamedVector) = isequal(keys(nv₁), keys(nv₂)) && isequal(values(nv₁), values(nv₂))
@inline Base.getindex(nv::NamedVector, index::Int) = getfield(nv, index)
@inline Base.setindex!(nv::NamedVector, value, index::Int) = setfield!(nv, index, value)
@inline Base.:<(nv₁::NamedVector, nv₂::NamedVector) = <(efficientoperations, nv₁, nv₂)
@inline Base.isless(nv₁::NamedVector, nv₂::NamedVector) = isless(efficientoperations, nv₁, nv₂)
@inline Base.hash(nv::NamedVector, h::UInt) = hash(values(nv), h)
@inline Base.length(::Type{NV}) where {NV<:NamedVector} = fieldcount(NV)
@inline Base.length(nv::NamedVector) = length(typeof(nv))
@inline Base.iterate(nv::NamedVector, state=1) = (state > length(nv)) ? nothing : (getfield(nv, state), state+1)
@inline Base.iterate(rv::Iterators.Reverse{<:NamedVector}, state=length(rv.itr)) = (state < 1) ? nothing : (getfield(rv.itr, state), state-1)
@inline Base.keys(nv::NamedVector) = fieldnames(typeof(nv))
@inline Base.values(nv::NamedVector) = convert(Tuple, nv)
@inline Base.pairs(nv::NamedVector) = Base.Generator(=>, keys(nv), values(nv))
Base.show(io::IO, nv::NamedVector) = @printf io "%s(%s)" nameof(typeof(nv)) join(repr.(values(nv)), ", ")

"""
    convert(::Type{Tuple}, nv::NamedVector) -> Tuple
    convert(::Type{NV}, nv::Tuple) where {NV<:NamedVector} -> NV

Convert a named vector to tuple and vice versa.
"""
@generated Base.convert(::Type{Tuple}, nv::NamedVector) = Expr(:tuple, (:(getfield(nv, $i)) for i = 1:fieldcount(nv))...)
function Base.convert(::Type{NV}, nv::Tuple) where NV<:NamedVector
    @assert fieldcount(NV) == length(nv) "convert error: mismatched length between $NV($(fieldcount(NV))) and input tuple($(length(nv)))."
    return NV(nv...)
end

"""
    zero(::Type{NV}) where NV<:NamedVector -> NV
    zero(nv::NamedVector) -> typeof(nv)

Get a concrete `NamedVector` with all values being zero.
"""
@inline @generated function Base.zero(::Type{NV}) where NV<:NamedVector
    zeros = (zero(fieldtype(NV, i)) for i = 1:fieldcount(NV))
    return :(NV($(zeros...)))
end
@inline Base.zero(nv::NamedVector) = zero(typeof(nv))

"""
    replace(nv::NamedVector; kwargs...) -> typeof(nv)

Return a copy of a concrete `NamedVector` with some of the field values replaced by the keyword arguments.
"""
@inline Base.replace(nv::NamedVector; kwargs...) = replace(efficientoperations, nv; kwargs...)

"""
    map(f, nvs::NV...) where NV<:NamedVector -> NV

Apply function `f` element-wisely on the input named vectors.
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
@inline Base.eltype(::Type{<:HomoNamedVector{T}}) where T = T
@inline Base.eltype(nv::HomoNamedVector) = eltype(typeof(nv))

end #module
