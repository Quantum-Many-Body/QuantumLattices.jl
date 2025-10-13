module DegreesOfFreedom

using DataStructures: OrderedDict
using Printf: @printf, @sprintf
using SparseArrays: SparseMatrixCSC, nnz
using StaticArrays: SVector
using ..QuantumLattices: OneAtLeast, OneOrMore, ZeroAtLeast, add!, decompose, str
using ..QuantumOperators: LinearTransformation, Operator, OperatorIndex, OperatorPack, Operators, QuantumOperator, scalartype, valuetolatextext
using ..Spatials: Bond, Point
using ..Toolkit: atol, efficientoperations, rtol, CompositeDict, Float, VectorSpace, VectorSpaceDirectProducted, VectorSpaceDirectSummed, concatenate, fulltype, parametertype, rawtype, reparameter

import LaTeXStrings: latexstring
import ..QuantumLattices: ⊕, ⊗, dimension, expand, expand!, id, ishermitian, kind, permute, rank, reset!, shape, update!, value
import ..QuantumOperators: idtype, operatortype, script
import ..Spatials: icoordinate, nneighbor, rcoordinate
import ..Toolkit: VectorSpaceStyle, contentnames, getcontent, parameternames

export CompositeInternal, CompositeIndex, CoordinatedIndex, Hilbert, Index, Internal, InternalIndex, InternalProd, InternalSum, SimpleInternal
export Boundary, Coupling, MatrixCoupling, MatrixCouplingComponent, MatrixCouplingProd, MatrixCouplingSum, Metric, OperatorIndexToTuple, Ordinal, Pattern, Table, Term, TermAmplitude, TermCoupling
export ˢᵗ, ⁿᵈ, ʳᵈ, ᵗʰ, plain, coordinatedindextype, diagonalfields, indextype, internalindextype, isdefinite, isdiagonal, partition, patternrule, showablefields, statistics, @pattern

# InternalIndex and Internal
"""
    InternalIndex <: OperatorIndex

Internal index of an internal degree of freedom.
"""
abstract type InternalIndex <: OperatorIndex end
@inline Base.show(io::IO, index::InternalIndex) = @printf io "%s(%s)" OperatorIndex[index] join(map(field->str(getfield(index, field)), showablefields(index)), ", ")

"""
    showablefields(index::InternalIndex) -> ZeroAtLeast{Symbol}
    showablefields(::Type{I}) where {I<:InternalIndex} -> ZeroAtLeast{Symbol}

Get the showable fields of an internal index or a type of internal index.
"""
@inline showablefields(index::InternalIndex) = showablefields(typeof(index))
@inline showablefields(::Type{I}) where {I<:InternalIndex} = fieldnames(I)

"""
    isdefinite(index::InternalIndex) -> Bool
    isdefinite(::Type{<:InternalIndex}) -> Bool

Judge whether an internal index or a type of an internal index denotes a definite internal degree of freedom.
"""
@inline isdefinite(index::InternalIndex) = isdefinite(typeof(index))
@inline isdefinite(::Type{<:InternalIndex}) = false

"""
    statistics(index::InternalIndex) -> Symbol

Get the statistics of an internal index.
"""
@inline statistics(index::InternalIndex) = statistics(typeof(index))

"""
    InternalIndex(index::OperatorIndex)

Get the `InternalIndex` part of an `OperatorIndex`.
"""
@inline InternalIndex(index::InternalIndex) = index

"""
    internalindextype(index::OperatorIndex)
    internalindextype(::Type{I}) where {I<:OperatorIndex}

Get the type of the `InternalIndex` part of an `OperatorIndex`.
"""
@inline internalindextype(index::OperatorIndex) = internalindextype(typeof(index))
@inline internalindextype(::Type{I}) where {I<:InternalIndex} = I

"""
    isdefinite(indexes::OneAtLeast{InternalIndex}) -> Bool
    isdefinite(::Type{T}) where {T<:OneAtLeast{InternalIndex}} -> Bool

Judge whether all of a set of internal indexes denote definite internal degrees of freedom.
"""
@inline isdefinite(indexes::OneAtLeast{InternalIndex}) = isdefinite(typeof(indexes))
@inline isdefinite(::Type{T}) where {T<:OneAtLeast{InternalIndex}} = all(map(isdefinite, fieldtypes(T)))

"""
    Internal{I<:OneOrMore{InternalIndex}} <: VectorSpace{I}

Internal space at a single point.
"""
abstract type Internal{I<:OneOrMore{InternalIndex}} <: VectorSpace{I} end

"""
    dimension(internal::Internal) -> Int

Get the dimension of an internal space at a single point.
"""
@inline dimension(internal::Internal) = length(internal)

"""
    SimpleInternal{I<:InternalIndex} <: Internal{I}

Simple internal space at a single point.
"""
abstract type SimpleInternal{I<:InternalIndex} <: Internal{I} end
@inline VectorSpaceStyle(::Type{<:SimpleInternal}) = VectorSpaceDirectProducted(:forward)
@inline Base.show(io::IO, internal::SimpleInternal) = @printf io "%s(%s)" typeof(internal) join(("$name=$(getfield(internal, name))" for name in fieldnames(typeof(internal))), ", ")

"""
    statistics(internal::SimpleInternal) -> Symbol
    statistics(::Type{<:SimpleInternal{I}}) where {I<:InternalIndex} -> Symbol

Get the statistics of a simple internal space.
"""
@inline statistics(internal::SimpleInternal) = statistics(typeof(internal))
@inline statistics(::Type{<:SimpleInternal{I}}) where {I<:InternalIndex} = statistics(I)

"""
    match(index::InternalIndex, internal::SimpleInternal) -> Bool
    match(::Type{I}, internal::SimpleInternal) where {I<:InternalIndex} -> Bool
    match(index::InternalIndex, ::Type{SI}) where {SI<:SimpleInternal} -> Bool
    match(::Type{I}, ::Type{SI}) where {I<:InternalIndex, SI<:SimpleInternal} -> Bool

Judge whether a simple internal space or the type of a simple internal space matches an internal index or the type of an internal index.

Here, "match" means that the `eltype` of the simple internal space is consistent with the type of the internal index, which usually means that they share the same type name.
"""
@inline Base.match(index::InternalIndex, internal::SimpleInternal) = match(typeof(index), typeof(internal))
@inline Base.match(::Type{I}, internal::SimpleInternal) where {I<:InternalIndex} = match(I, typeof(internal))
@inline Base.match(index::InternalIndex, ::Type{SI}) where {SI<:SimpleInternal} = match(typeof(index), SI)
@inline @generated Base.match(::Type{I}, ::Type{SI}) where {I<:InternalIndex, SI<:SimpleInternal} = nameof(I)==nameof(eltype(SI))

"""
    filter(index::InternalIndex, internal::SimpleInternal) -> Union{Nothing, typeof(internal)}
    filter(::Type{I}, internal::SimpleInternal) where {I<:InternalIndex} -> Union{Nothing, typeof(internal)}

Filter a simple internal space with respect to the input internal index or the type of an internal index.
"""
@inline Base.filter(index::InternalIndex, internal::SimpleInternal) = filter(typeof(index), internal)
@inline Base.filter(::Type{I}, internal::SimpleInternal) where {I<:InternalIndex} = filterhelper(Val(match(I, typeof(internal))), internal)
@inline @generated filterhelper(::Val{B}, internal) where B = B ? :(internal) : nothing

"""
    filter(index::InternalIndex, ::Type{SI}) where {SI<:SimpleInternal}
    filter(::Type{I}, ::Type{SI}) where {I<:InternalIndex, SI<:SimpleInternal}

Filter the type of a simple internal space with respect to the input internal index or the type of an internal index.
"""
@inline Base.filter(index::InternalIndex, ::Type{SI}) where {SI<:SimpleInternal} = filter(typeof(index), SI)
@inline Base.filter(::Type{I}, ::Type{SI}) where {I<:InternalIndex, SI<:SimpleInternal} = filterhelper(Val(match(I, SI)), SI)

"""
    shape(internal::SimpleInternal, index::InternalIndex) -> OrdinalRange{Int, Int}

Get the shape of a simple internal space when a labeled internal index is considered.

It can be overloaded to restrict the shape of a simple internal space based on the input internal index to significantly improve efficiency, but this is not necessary.
"""
@inline shape(internal::SimpleInternal, ::InternalIndex) = shape(internal)

"""
    CompositeInternal{T<:OneAtLeast{SimpleInternal}, I<:OneOrMore{InternalIndex}} <: Internal{I}

Abstract type of the composition (i.e., direct sum or direct product) of several simple internal spaces.
"""
abstract type CompositeInternal{T<:OneAtLeast{SimpleInternal}, I<:OneOrMore{InternalIndex}} <: Internal{I} end

"""
    rank(internal::CompositeInternal) -> Int
    rank(::Type{<:CompositeInternal{T}}) where {T<:OneAtLeast{SimpleInternal}} -> Int

Get the number of simple internal spaces in a composite internal space.
"""
@inline rank(internal::CompositeInternal) = rank(typeof(internal))
@inline rank(::Type{<:CompositeInternal{T}}) where {T<:OneAtLeast{SimpleInternal}} = fieldcount(T)

"""
    filter(index::InternalIndex, internal::CompositeInternal) -> Union{Nothing, SimpleInternal, CompositeInternal}
    filter(::Type{I}, internal::CompositeInternal) where {I<:InternalIndex} -> Union{Nothing, SimpleInternal, CompositeInternal}

Filter the composite internal space and select those that match the input internal index or the type of an internal index.
"""
@inline Base.filter(index::InternalIndex, internal::CompositeInternal) = filter(typeof(index), internal)
@inline Base.filter(::Type{I}, internal::CompositeInternal{T}) where {I<:InternalIndex, T<:OneAtLeast{SimpleInternal}} = filterhelper(Val(map(SI->match(I, SI), fieldtypes(T))), I, internal)
@inline @generated function filterhelper(::Val{BS}, ::Type{I}, internal::CompositeInternal) where {BS, I<:InternalIndex}
    exprs = []
    for (i, B) in enumerate(BS)
        B && push!(exprs, :(filter(I, internal.contents[$i])))
    end
    length(exprs)==0 && return
    length(exprs)==1 && return only(exprs)
    return :(rawtype(typeof(internal))($(exprs...)))
end

"""
    filter(index::InternalIndex, ::Type{CI}) where {CI<:CompositeInternal}
    filter(::Type{I}, ::Type{CI}) where {I<:InternalIndex, CI<:CompositeInternal}

Filter the type of a composite internal space and select those that match the input internal index or the type of an internal index.
"""
@inline Base.filter(index::InternalIndex, ::Type{CI}) where {CI<:CompositeInternal} = filter(typeof(index), CI)
@inline Base.filter(::Type{I}, ::Type{CI}) where {I<:InternalIndex, T<:OneAtLeast{SimpleInternal}, CI<:CompositeInternal{T}} = filterhelper(Val(map(SI->match(I, SI), fieldtypes(T))), I, CI)
@inline @generated function filterhelper(::Val{BS}, ::Type{I}, ::Type{CI}) where {BS, I<:InternalIndex, T<:OneAtLeast{SimpleInternal}, CI<:CompositeInternal{T}}
    exprs = []
    for (i, B) in enumerate(BS)
        B && push!(exprs, :(filter(I, fieldtype(T, $i))))
    end
    length(exprs)==0 && return
    length(exprs)==1 && return first(exprs)
    return Expr(:curly, nameof(CI), Expr(:curly, :Tuple, exprs...), :(eltype(rawtype(CI), $(exprs...))))
end

"""
    InternalSum{T<:OneAtLeast{SimpleInternal}, I<:InternalIndex} <: CompositeInternal{T, I}

Direct sum of several single internal spaces.
"""
struct InternalSum{T<:OneAtLeast{SimpleInternal}, I<:InternalIndex} <: CompositeInternal{T, I}
    contents::T
    InternalSum(contents::OneAtLeast{SimpleInternal}) = new{typeof(contents), eltype(InternalSum, fieldtypes(typeof(contents))...)}(contents)
end
@inline VectorSpaceStyle(::Type{<:InternalSum}) = VectorSpaceDirectSummed()
@inline Base.eltype(::Type{InternalSum}, types...) = mapreduce(eltype, typejoin, types; init=Union{})
@inline Base.show(io::IO, internal::InternalSum) = @printf io "%s" join((string(internal.contents[i]) for i = 1:rank(internal)), " ⊕ ")
@inline InternalSum(contents::SimpleInternal...) = InternalSum(contents)

"""
    ⊕(internal::SimpleInternal, internals::SimpleInternal...) -> InternalSum
    ⊕(internal::SimpleInternal, internals::InternalSum) -> InternalSum
    ⊕(internals::InternalSum, internal::SimpleInternal) -> InternalSum
    ⊕(internals₁::InternalSum, internals₂::InternalSum) -> InternalSum

Direct sum between simple internal spaces and composite internal spaces.
"""
@inline ⊕(internal::SimpleInternal, internals::SimpleInternal...) = InternalSum(internal, internals...)
@inline ⊕(internal::SimpleInternal, internals::InternalSum) = InternalSum(internal, internals.contents...)
@inline ⊕(internals::InternalSum, internal::SimpleInternal) = InternalSum(internals.contents..., internal)
@inline ⊕(internals₁::InternalSum, internals₂::InternalSum) = InternalSum(internals₁.contents..., internals₂.contents...)

"""
    InternalProd{T<:OneAtLeast{SimpleInternal}, I<:OneAtLeast{InternalIndex}} <: CompositeInternal{T, I}

Direct product of several single internal spaces.
"""
struct InternalProd{T<:OneAtLeast{SimpleInternal}, I<:OneAtLeast{InternalIndex}} <: CompositeInternal{T, I}
    contents::T
    InternalProd(contents::OneAtLeast{SimpleInternal}) = new{typeof(contents), eltype(InternalProd, fieldtypes(typeof(contents))...)}(contents)
end
@inline VectorSpaceStyle(::Type{<:InternalProd}) = VectorSpaceDirectProducted(:forward)
@inline Base.eltype(::Type{InternalProd}, types...) = Tuple{map(eltype, types)...}
@inline Base.convert(::Type{<:OneAtLeast{InternalIndex}}, index::CartesianIndex, internal::InternalProd) = map(getindex, internal.contents, index.I)
@inline Base.show(io::IO, internal::InternalProd) = @printf io "%s" join((string(internal.contents[i]) for i = 1:rank(internal)), " ⊗ ")
@inline InternalProd(contents::SimpleInternal...) = InternalProd(contents)

"""
    ⊗(internal::SimpleInternal, internals::SimpleInternal...) -> InternalProd
    ⊗(internal::SimpleInternal, internals::InternalProd) -> InternalProd
    ⊗(internals::InternalProd, internal::SimpleInternal) -> InternalProd
    ⊗(internals₁::InternalProd, internals₂::InternalProd) -> InternalProd

Direct product between simple internal spaces and composite internal spaces.
"""
@inline ⊗(internal::SimpleInternal, internals::SimpleInternal...) = InternalProd(internal, internals...)
@inline ⊗(internal::SimpleInternal, internals::InternalProd) = InternalProd(internal, internals.contents...)
@inline ⊗(internals::InternalProd, internal::SimpleInternal) = InternalProd(internals.contents..., internal)
@inline ⊗(internals₁::InternalProd, internals₂::InternalProd) = InternalProd(internals₁.contents..., internals₂.contents...)

# Index
"""
    Ordinal

Ordinal of an `Int`.
"""
struct Ordinal
    n::Int
end
@inline Base.show(io::IO, nth::Ordinal) = @printf io "%s" nth.n==1 ? "1ˢᵗ" : nth.n==2 ? "2ⁿᵈ" : nth.n==3 ? "3ʳᵈ" : "$(nth.n)ᵗʰ"
@inline Base.:*(i::Int, nth::Ordinal) = Ordinal(i*nth.n)
@inline Base.getindex(v::Tuple, i::Ordinal) = v[i.n]
@inline Base.getindex(v::Bond, i::Ordinal) = v[i.n]

"""
    const ˢᵗ = ⁿᵈ = ʳᵈ = ᵗʰ = Ordinal(1)

Constant ordinals.
"""
const ˢᵗ = ⁿᵈ = ʳᵈ = ᵗʰ = Ordinal(1)

"""
    Index(site::Union{Int, Ordinal, Colon}, internal::InternalIndex)

Index of a degree of freedom, which consist of the spatial part (i.e., the site index) and the internal part (i.e., the internal index).
"""
struct Index{I<:InternalIndex, S<:Union{Int, Ordinal, Colon}} <: OperatorIndex
    site::S
    internal::I
end
@inline parameternames(::Type{<:Index}) = (:internal, :site)
function Base.show(io::IO, index::Index)
    internal = String[]
    for field in showablefields(internalindextype(index))
        value = getfield(index.internal, field)
        push!(internal, str(value))
    end
    internal = join(internal, ", ")
    @printf io "%s(%s%s%s)" OperatorIndex[index.internal] str(index.site) (length(internal)>0 ? ", " : "") internal
end
@inline InternalIndex(index::Index) = index.internal
@inline internalindextype(::Type{I}) where {I<:Index} = parametertype(I, :internal)

"""
    Index(index::OperatorIndex)

Get the `Index` part of an `OperatorIndex`.
"""
@inline Index(index::Index) = index

"""
    indextype(index::OperatorIndex)
    indextype(::Type{I}) where {I<:OperatorIndex}

Get the type of the `Index` part of an `OperatorIndex`.
"""
@inline indextype(index::OperatorIndex) = indextype(typeof(index))
@inline indextype(::Type{I}) where {I<:Index} = I

"""
    isdefinite(index::Index) -> Bool
    isdefinite(::Type{<:Index{I}}) where {I<:InternalIndex} -> Bool

Determine whether an index denotes a definite degree of freedom.
"""
@inline isdefinite(index::Index) = isdefinite(typeof(index))
@inline isdefinite(::Type{<:Index{I}}) where {I<:InternalIndex} = isdefinite(I)

"""
    isdefinite(indexes::OneAtLeast{Index}) -> Bool
    isdefinite(::Type{T}) where {T<:OneAtLeast{Index}} -> Bool

Determine whether a tuple of indexes denotes a definite degree of freedom.
"""
@inline isdefinite(indexes::OneAtLeast{Index}) = isdefinite(typeof(indexes))
@inline isdefinite(::Type{T}) where {T<:OneAtLeast{Index}} = all(map(isdefinite, fieldtypes(T)))

"""
    statistics(index::Index) -> Symbol
    statistics(::Type{I}) where {I<:Index} -> Symbol

Get the statistics of an index.
"""
@inline statistics(index::Index) = statistics(typeof(index))
@inline statistics(::Type{I}) where {I<:Index} = statistics(internalindextype(I))

"""
    adjoint(index::Index) -> typeof(index)

Get the adjoint of an index.
"""
@inline Base.adjoint(index::Index) = rawtype(typeof(index))(index.site, adjoint(index.internal))

"""
    script(index::Index, ::Val{:site}; kwargs...) -> String
    script(index::Index, attr::Val; kwargs...) -> String

Get the required script of an index.
"""
@inline script(index::Index, ::Val{:site}; kwargs...) = str(index.site)
@inline script(index::Index, attr::Val; kwargs...) = script(index.internal, attr; kwargs...)

"""
    permute(index₁::Index, index₂::Index) -> ZeroAtLeast{Operator}

Get the permutation of two indexes.
"""
function permute(index₁::Index, index₂::Index)
    if index₁.site ≠ index₂.site
        return (Operator(statistics(index₁)==:f && statistics(index₂)==:f ? -1 : 1, index₂, index₁),)
    else
        return map(op->Operator(value(op), map(internal->Index(index₁.site, internal), id(op))), permute(index₁.internal, index₂.internal))
    end
end

"""
    indextype(I::Type{<:SimpleInternal})

Get the compatible type of the index based on the type of an internal space.
"""
@inline indextype(I::Type{<:SimpleInternal}) = fulltype(Index, NamedTuple{(:internal, :site), Tuple{eltype(I), Int}})

# CompositeIndex and CoordinatedIndex
"""
    CompositeIndex{I<:Index} <: OperatorIndex

Abstract type of a composite index.
"""
abstract type CompositeIndex{I<:Index} <: OperatorIndex end
@inline contentnames(::Type{<:CompositeIndex}) = (:index,)
@inline parameternames(::Type{<:CompositeIndex}) = (:index,)
@inline InternalIndex(index::CompositeIndex) = InternalIndex(Index(index))
@inline internalindextype(::Type{<:CompositeIndex{I}}) where {I<:Index} = internalindextype(I)
@inline Index(index::CompositeIndex) = getcontent(index, :index)
@inline @generated indextype(::Type{I}) where {I<:CompositeIndex} = parametertype(supertype(I, :CompositeIndex), :index)

"""
    statistics(index::CompositeIndex) -> Symbol
    statistics(::Type{I}) where {I<:CompositeIndex} -> Symbol

Get the statistics of a composite index.
"""
@inline statistics(index::CompositeIndex) = statistics(typeof(index))
@inline statistics(::Type{I}) where {I<:CompositeIndex} = statistics(indextype(I))

"""
    script(index::CompositeIndex, ::Val{attr}; kwargs...) where attr

Get the `attr` script of a composite index.
"""
@inline script(index::CompositeIndex, ::Val{attr}; kwargs...) where attr = script(getcontent(index, :index), Val(attr); kwargs...)

"""
    CoordinatedIndex{I<:Index, V<:SVector} <: CompositeIndex{I}

Coordinated index, i.e., index with coordinates in the unitcell.
"""
struct CoordinatedIndex{I<:Index, V<:SVector} <: CompositeIndex{I}
    index::I
    rcoordinate::V
    icoordinate::V
    CoordinatedIndex(index::Index, rcoordinate::V, icoordinate::V) where {V<:SVector} = new{typeof(index), V}(index, compositeindexcoordinate(rcoordinate), compositeindexcoordinate(icoordinate))
end
@inline contentnames(::Type{<:CoordinatedIndex}) = (:index, :rcoordinate, :icoordinate)
@inline parameternames(::Type{<:CoordinatedIndex}) = (:index, :coordination)
@inline Base.hash(index::CoordinatedIndex, h::UInt) = hash((index.index, Tuple(index.rcoordinate)), h)
@inline Base.propertynames(::OneAtLeast{CoordinatedIndex}) = (:indexes, :rcoordinates, :icoordinates)
function Base.show(io::IO, index::CoordinatedIndex)
    internal = String[]
    for field in showablefields(internalindextype(index))
        value = getfield(InternalIndex(index), field)
        push!(internal, str(value))
    end
    internal = join(internal, ", ")
    @printf io "%s(%s%s%s, %s, %s)" OperatorIndex[index.index.internal] str(index.index.site) (length(internal)>0 ? ", " : "") internal index.rcoordinate index.icoordinate
end
@inline compositeindexcoordinate(vector::SVector) = vector
@inline compositeindexcoordinate(vector::SVector{N, Float}) where N = SVector(ntuple(i->vector[i]===-0.0 ? 0.0 : vector[i], Val(N)))

"""
    CoordinatedIndex(index::Index, rcoordinate, icoordinate)
    CoordinatedIndex(index::Index; rcoordinate, icoordinate)

Construct a coordinated index.
"""
@inline CoordinatedIndex(index::Index, rcoordinate, icoordinate) = CoordinatedIndex(index, SVector{length(rcoordinate)}(rcoordinate), SVector{length(icoordinate)}(icoordinate))
@inline CoordinatedIndex(index::Index; rcoordinate, icoordinate) = CoordinatedIndex(index, rcoordinate, icoordinate)

"""
    adjoint(index::CoordinatedIndex) -> typeof(index)

Get the adjoint of a coordinated index.
"""
@inline Base.adjoint(index::CoordinatedIndex) = CoordinatedIndex(index.index', index.rcoordinate, index.icoordinate)

"""
    script(index::CoordinatedIndex, ::Val{:rcoordinate}; kwargs...) -> String
    script(index::CoordinatedIndex, ::Val{:icoordinate}; kwargs...) -> String

Get the rcoordinate/icoordinate script of a coordinated index.
"""
@inline script(index::CoordinatedIndex, ::Val{:rcoordinate}; kwargs...) = @sprintf "[%s]" join(valuetolatextext.(index.rcoordinate), ", ")
@inline script(index::CoordinatedIndex, ::Val{:icoordinate}; kwargs...) = @sprintf "[%s]" join(valuetolatextext.(index.icoordinate), ", ")

"""
    script(index::CoordinatedIndex, ::Val{:integercoordinate}; vectors, kwargs...)

Get the integer coordinate script of a coordinated index.
"""
function script(index::CoordinatedIndex, ::Val{:integercoordinate}; vectors, kwargs...)
    rcoeff = decompose(index.icoordinate, vectors...)
    icoeff = Int.(round.(rcoeff))
    @assert isapprox(efficientoperations, rcoeff, icoeff) "script error: mismatched icoordinate of the input coordinated index and vectors."
    return @sprintf "[%s]" join(icoeff, ", ")
end

"""
    permute(index₁::CoordinatedIndex, index₂::CoordinatedIndex) -> ZeroAtLeast{Operator}

Get the permutation of two coordinated indexes.
"""
function permute(index₁::CoordinatedIndex, index₂::CoordinatedIndex)
    if index₁.rcoordinate ≠ index₂.rcoordinate || index₁.icoordinate ≠ index₂.icoordinate
        return (Operator(statistics(index₁)==:f && statistics(index₂)==:f ? -1 : 1, index₂, index₁),)
    else
        rcoordinate, icoordinate = index₁.rcoordinate, index₁.icoordinate
        return map(op->Operator(value(op), map(index->CoordinatedIndex(index, rcoordinate, icoordinate), id(op))), permute(index₁.index, index₂.index))
    end
end

"""
    coordinatedindextype(I::Type{<:SimpleInternal}, P::Type{<:Point})

Get the compatible type of the coordinated index based on the type of an internal space and the type of a point.
"""
@inline coordinatedindextype(I::Type{<:SimpleInternal}, P::Type{<:Point}) = fulltype(CoordinatedIndex, NamedTuple{(:index, :coordination), Tuple{indextype(I), SVector{dimension(P), scalartype(P)}}})

"""
    rcoordinate(opt::Operator{<:Number, <:ZeroAtLeast{CoordinatedIndex}}) -> SVector

Get the whole rcoordinate of an operator.
"""
@inline function rcoordinate(opt::Operator{<:Number, <:ZeroAtLeast{CoordinatedIndex}})
    rank(opt)==1 && return id(opt)[1].rcoordinate
    rank(opt)==2 && return id(opt)[2].rcoordinate-id(opt)[1].rcoordinate
    error("rcoordinate error: not supported rank($(rank(opt))) of $(nameof(opt)).")
end

"""
    icoordinate(opt::Operator{<:Number, <:ZeroAtLeast{CoordinatedIndex}}) -> SVector

Get the whole icoordinate of an operator.
"""
@inline function icoordinate(opt::Operator{<:Number, <:ZeroAtLeast{CoordinatedIndex}})
    rank(opt)==1 && return id(opt)[1].icoordinate
    rank(opt)==2 && return id(opt)[2].icoordinate-id(opt)[1].icoordinate
    error("icoordinate error: not supported rank($(rank(opt))) of $(nameof(opt)).")
end

# Hilbert
"""
    Hilbert{I<:Internal} <: CompositeDict{Int, I}

Hilbert space at a lattice.
"""
struct Hilbert{I<:Internal} <: CompositeDict{Int, I}
    contents::OrderedDict{Int, I}
    Hilbert(contents::OrderedDict{Int, <:Internal}) = new{valtype(contents)}(contents)
end

"""
    Hilbert(ps::Pair...)
    Hilbert(kv)

Construct a Hilbert space the same way as an `OrderedDict`.
"""
@inline Hilbert(p::Pair, ps::Pair...) = Hilbert((p, ps...))
@inline Hilbert(kv) = Hilbert(OrderedDict(kv))

"""
    Hilbert(internals::Internal...)
    Hilbert(internals::OneAtLeast{Internal})
    Hilbert(internals::AbstractVector{<:Internal})

Construct a Hilbert space with the given internal spaces.
"""
@inline Hilbert(internal::Internal, internals::Internal...) = Hilbert((internal, internals...))
@inline Hilbert(internals::OneAtLeast{Internal}) = Hilbert(i=>internal for (i, internal) in enumerate(internals))
@inline Hilbert(internals::AbstractVector{<:Internal}) = Hilbert(i=>internal for (i, internal) in enumerate(internals))

"""
    Hilbert(internal::Internal, num::Integer)

Construct a Hilbert space with all same internal spaces.
"""
@inline Hilbert(internal::Internal, num::Integer) = Hilbert(i=>internal for i in 1:num)

# Pattern
"""
    diagonalfields(::Type{I}) where {I<:InternalIndex} -> ZeroAtLeast{Symbol}

Get the field names that can be subject to all-equal constraint based on the type of an internal index.
"""
@inline diagonalfields(::Type{I}) where {I<:InternalIndex} = fieldnames(I)

"""
    isdiagonal(indexes::OneAtLeast{InternalIndex}, ::Val{fields}) where fields -> Bool

Judge whether a set of homogenous internal indexes is subject to the "diagonal" constraint.
"""
@inline isdiagonal(::Type{I}, indexes::OneAtLeast{InternalIndex}) where {I<:InternalIndex} = _isdiagonal_(I, I|>diagonalfields|>Val, indexes)
@generated function _isdiagonal_(::Type{I}, ::Val{candidates}, indexes::OneAtLeast{InternalIndex}) where {I<:InternalIndex, candidates}
    fields = Symbol[]
    for field in candidates
        F = fieldtype(I, field)
        F==Symbol && error("isdiagonal error: $I is too complicated for the \"diagonal\" constraint since it host the `:$field` field represented by a `Symbol` other than a `Colon`.")
        F==Colon && push!(fields, field)
    end
    length(fields)==0 && return true
    exprs = []
    for field in QuoteNode.(fields)
        push!(exprs, Expr(:call, allequal, Expr(:tuple, [:(getfield(indexes[$i], $field)) for i=1:fieldcount(indexes)]...)))
    end
    return Expr(:call, all, Expr(:tuple, exprs...))
end

"""
    Pattern{I, P, N, C<:OneAtLeast{Function}} <: QuantumOperator

Coupling pattern.
"""
struct Pattern{I, P, N, C<:OneAtLeast{Function}} <: QuantumOperator
    indexes::I
    constraints::C
    representations::NTuple{N, String}
    function Pattern{P}(indexes::OneAtLeast{Index}, constraints::NTuple{N, Function}, representations::NTuple{N, String}=map(string, constraints))  where {P, N}
        @assert N>0 "Pattern error: non-positive N($N)."
        @assert isa(P, NTuple{N, Int}) "Pattern error: partition ($P) is not $N-tuple of integers."
        @assert isa(indexes.sites, OneAtLeast{Union{Ordinal, Colon}}) "Pattern error: for each index, the site must be an `Ordinal` or the `:`."
        @assert sum(P)==length(indexes) "Pattern error: mismatched sum of partition ($(sum(P))) and number of indexes ($(length(indexes)))."
        new{typeof(indexes), P, N, typeof(constraints)}(indexes, constraints, representations)
    end
end
@inline parameternames(::Type{<:Pattern}) = (:indexes, :partition, :npartition, :constraints)
@inline Base.:(==)(pattern₁::Pattern, pattern₂::Pattern) = partition(pattern₁)==partition(pattern₂) && pattern₁.indexes==pattern₂.indexes && pattern₁.representations==pattern₂.representations
@inline Base.:isequal(pattern₁::Pattern, pattern₂::Pattern) = isequal(partition(pattern₁), partition(pattern₂)) && isequal(pattern₁.indexes, pattern₂.indexes) && isequal(pattern₁.representations, pattern₂.representations)
@inline Base.hash(pattern::Pattern, h::UInt) = hash((partition(pattern)..., pattern.indexes..., pattern.representations...), h)
function Base.show(io::IO, pattern::Pattern)
    len, count = length(pattern.representations), 1
    for i = 1:len
        r = rank(pattern, i)
        start, stop = count, count+r-1
        indexes = pattern[start:stop]
        if i==1
            @printf io "%s" (isdefinite(indexes) ? "" : "∑")
            (len>1 || !isdefinite(indexes)) && @printf io "%s" "["
        else
            @printf io " ⊗ %s[" (isdefinite(indexes) ? "" : "∑")
        end
        @printf io "%s" join(indexes, " ")
        (len>1 || !isdefinite(indexes)) && @printf io "%s" "]"
        representation = pattern.representations[i]
        length(representation)==0 || occursin("isdiagonal", representation) || @printf io "(%s)" representation
        count = stop+1
    end
end

"""
    Pattern(indexes::OneAtLeast{Index}, constraint::Function, representation::String=string(constraint))
    Pattern{P}(indexes::OneAtLeast{Index}, constraints::NTuple{N, Function}, representations::NTuple{N, String}=map(string, constraints))  where {P, N}

Construct a coupling pattern with 1) only one constraint function, and 2) several constraint functions.
"""
@inline Pattern(indexes::OneAtLeast{Index}, constraint::Function, representation::String=string(constraint)) = Pattern{(rank(indexes),)}(indexes, (constraint,), (representation,))

"""
    Pattern(indexes::OneAtLeast{Index})
    Pattern(index::Index, indexes::Index...)

Construct a coupling pattern subject to the "diagonal" constraint for a homogeneous set of indexes.
"""
@inline Pattern(index::Index, indexes::Index...) = Pattern((index, indexes...))
@inline Pattern(indexes::OneAtLeast{Index}) = Pattern(indexes, isdiagonal)

"""
    getindex(pattern::Pattern, slice)

Get the indexes specified by `slice` in a coupling pattern.
"""
@inline Base.getindex(pattern::Pattern, slice) = pattern.indexes[slice]

"""
    partition(pattern::Pattern) -> OneAtLeast{Int}
    partition(::Type{<:Pattern{I, P} where I}) where P -> P

Get the partition of the coupling pattern.
"""
@inline partition(pattern::Pattern) = partition(typeof(pattern))
@inline partition(::Type{<:Pattern{I, P} where I}) where P = P

"""
    rank(pattern::Pattern) -> Int
    rank(::Type{P}) where {P<:Pattern} -> Int

Get the total rank of the coupling pattern.
"""
@inline rank(pattern::Pattern) = rank(typeof(pattern))
@inline @generated rank(::Type{P}) where {P<:Pattern} = sum(partition(P))

"""
    rank(pattern::Pattern, i::Integer) -> Int
    rank(::Type{P}, i::Integer) where {P<:Pattern} -> Int

Get the ith rank of the coupling pattern where the ith constraint can apply.
"""
@inline rank(pattern::Pattern, i::Integer) = rank(typeof(pattern), i)
@inline rank(::Type{P}, i::Integer) where {P<:Pattern} = partition(P)[i]

"""
    latexstring(pattern::Pattern) -> String

Convert a coupling pattern to the latex format.
"""
function latexstring(pattern::Pattern)
    result = String[]
    len, count = length(pattern.representations), 1
    for i = 1:len
        r = rank(pattern, i)
        start, stop = count, count+r-1
        indexes = pattern[start:stop]
        representation = pattern.representations[i]
        summation = occursin("isdiagonal", representation) ? "" : replace(replace(representation, "&&"=>"\\,\\text{and}\\,"), "||"=>"\\,\\text{or}\\,")
        summation=="" || (summation = join(push!(symbols(indexes.internals, representation), summation), ",\\,"))
        context = isdefinite(indexes) ? "" : "\\sum_{$summation}" 
        for content in indexes
            context = @sprintf "%s%s%s" context (length(context)>0 ? " " : "") latexstring(content)
        end
        push!(result, context)
        count = stop+1
    end
    return join(result, " \\cdot ")
end
function symbols(indexes::OneAtLeast{InternalIndex}, constraint::String)
    result = String[]
    for index in indexes
        for i = 1:fieldcount(typeof(index))
            attr = getfield(index, i)
            if isa(attr, Symbol)
                value = string(attr)
                occursin(value, constraint) || push!(result, value)
            end
        end
    end
    return sort!(unique!(result))
end

"""
    match(pattern::Pattern, indexes::OneAtLeast{InternalIndex}) -> Bool

Judge whether a set of internal indexes satisfies a coupling pattern.
"""
@generated function Base.match(pattern::Pattern, indexes::OneAtLeast{InternalIndex})
    @assert rank(pattern)==rank(indexes) "match error: mismatched ranks of the coupling pattern and the input indexes."
    exprs, count = [], 1
    for (i, r) in enumerate(partition(pattern))
        start, stop = count, count+r-1
        segment = Expr(:tuple, [:(indexes[$pos]) for pos=start:stop]...)
        if fieldtype(fieldtype(pattern, :constraints), i) == typeof(isdiagonal)
            I = Expr(:call, eltype, Expr(:tuple, [:(pattern.indexes[$pos].internal) for pos=start:stop]...))
            push!(exprs, :(pattern.constraints[$i]($I, $segment)::Bool))
        else
            push!(exprs, :(pattern.constraints[$i]($segment)::Bool))
        end
        count = stop+1
    end
    return Expr(:call, :all, Expr(:tuple, exprs...))
end

"""
    ⊗(pattern₁::Pattern, pattern₂::Pattern) -> Pattern

Get the combination of two coupling patterns.
"""
@inline function ⊗(pattern₁::Pattern, pattern₂::Pattern)
    return Pattern{(partition(pattern₁)..., partition(pattern₂)...)}((pattern₁.indexes..., pattern₂.indexes...), (pattern₁.constraints..., pattern₂.constraints...), (pattern₁.representations..., pattern₂.representations...))
end

"""
    @pattern index₁ index₂ ...
    @pattern(index₁, index₂, ...; constraint=...)

Construct a coupling pattern according to the pattern of the input indexes and an optional constraint.
"""
macro pattern(exprs...)
    if exprs[1].head==:parameters
        @assert(
            length(exprs[1].args)==1 && exprs[1].args[1].head==:kw && exprs[1].args[1].args[1]==:constraint,
            "@pattern error: constraint must be specified by the keyword argument `constraint`."
        )
        constraint = exprs[1].args[1].args[2]
        patterns = exprs[2:end]
    else
        constraint = true
        patterns = exprs
    end
    indexes = []
    blocks = []
    for (i, pattern) in enumerate(patterns)
        if pattern.head==:call && length(pattern.args)==3 && isa(pattern.args[3], Expr) && pattern.args[3].head==:call && pattern.args[3].args[1]≠://
            candidates = pattern.args[3].args[2:end]
            flag = 1
        elseif pattern.head==:call && length(pattern.args)>=3
            candidates = pattern.args[3:end]
            flag = 2
        else
            error("@pattern error: wrong pattern.")
        end
        attrs = []
        for (j, attr) in enumerate(candidates)
            isa(attr, Expr) && !(attr.head==:call && length(attr.args)==3 && attr.args[1]==:// && isa(attr.args[2], Int) && isa(attr.args[3], Int)) && (attr = Symbol(attr))
            if isa(attr, Symbol)
                @assert attr≠Symbol(":") "@pattern error: Colon(:) cannot be used for internal indexes in @pattern."
                @assert Base.isidentifier(attr) "@pattern error: wrong pattern."
                if !isdefined(__module__, attr)
                    push!(blocks, quote
                        local $attr
                        if @isdefined($attr)
                            ($attr)!=getfield(indexes[$i], $j) && return false
                        else
                            $attr = getfield(indexes[$i], $j)
                        end
                    end)
                    attr = QuoteNode(attr)
                end
            else
                push!(blocks, quote
                    (getfield(indexes[$i], $j)!=$attr) && return false
                end)
            end
            push!(attrs, attr)
        end
        if flag == 1
            index = Expr(:call, pattern.args[1], pattern.args[2], Expr(:call, pattern.args[3].args[1], attrs...))
        elseif flag == 2
            index = Expr(:call, pattern.args[1], pattern.args[2], attrs...)
        end
        push!(indexes, index)
    end
    N = length(indexes)
    indexes = Expr(:tuple, indexes...)
    blocks = Expr(:block, vcat([block.args for block in blocks]...)...)
    iname, fname = gensym("indexes"), gensym("constraint")
    representation = constraint==true ? "" : string(constraint)
    return quote
        function ($fname)(indexes::NTuple{$N, InternalIndex})
            $blocks
            return $constraint
        end
        local $iname = $(esc(indexes))
        Pattern(map(Index, $iname.sites, $iname.internals), $fname, $representation)
    end
end

# patternrule
"""
    patternrule(value, ::Val{Name}, args...; kwargs...) where Name

By overriding this function, a named rule for the value of some attributes in a pattern can be defined so that the value can be transformed into the desired one.

In most cases, such overridden functions define the default rules for the attributes in a pattern when they take on the default value `:`. 
"""
function patternrule end

"""
    patternrule(value, ::Val, args...; kwargs...)-> typeof(value)

Default pattern rule unless otherwise specified.
"""
@inline patternrule(value, ::Val, args...; kwargs...) = value

"""
    patternrule(sites::OneAtLeast{Colon}, ::Val, bondlength::Integer) -> OneAtLeast{Ordinal}

Define the default rule for the sites of a set of indexes in a coupling pattern.
"""
function patternrule(sites::OneAtLeast{Colon}, ::Val, bondlength::Integer)
    N = fieldcount(typeof(sites))
    bondlength==1 && return ntuple(i->1ˢᵗ, Val(N))
    bondlength==2 && begin
        N==2 && return (1ˢᵗ, 2ⁿᵈ)
        N==4 && return (1ˢᵗ, 1ˢᵗ, 2ⁿᵈ, 2ⁿᵈ)
        error("patternrule error: not implemented.")
    end
    error("patternrule error: not implemented.")
end

"""
    patternrule(indexes::OneAtLeast{InternalIndex}, ::Val{Name}) where Name -> OneAtLeast{InternalIndex}

Define the default rule for the internal index in a coupling pattern.
"""
@inline @generated function patternrule(indexes::OneAtLeast{InternalIndex}, ::Val{Name}) where Name
    allequal(map(nameof, fieldtypes(indexes))) || return :(indexes)
    exprs = []
    for name in QuoteNode.(fieldnames(eltype(indexes)))
        vs = Expr(:tuple, [:(getfield(indexes[$i], $name)) for i = 1:fieldcount(indexes)]...)
        push!(exprs, :(patternrule($vs, Val($(QuoteNode(Name))), eltype(indexes), Val($name))))
    end
    return :(map(apply, fieldtypes(typeof(indexes)), $(exprs...)))
end
@inline apply(::Type{F}, args...) where F = F(args...)

# Coupling
"""
    Coupling{V, P<:Pattern} <: OperatorPack{V, P}

Coupling among internal degrees of freedom at different or same lattice points.
"""
struct Coupling{V, P<:Pattern} <: OperatorPack{V, P}
    value::V
    pattern::P
end
@inline getcontent(coupling::Coupling, ::Val{:id}) = coupling.pattern
@inline Base.length(coupling::Coupling) = length(typeof(coupling))
@inline Base.length(::Type{<:Coupling}) = 1
@inline Base.eltype(coupling::Coupling) = eltype(typeof(coupling))
@inline Base.eltype(C::Type{<:Coupling}) = C
@inline Base.iterate(coupling::Coupling) = (coupling, nothing)
@inline Base.iterate(::Coupling, ::Nothing) = nothing
@inline Base.show(io::IO, coupling::Coupling) = @printf io "%s%s" (coupling.value≈1 ? "" : coupling.value≈-1 ? "- " : string(str(coupling.value), " ")) coupling.pattern
@inline Base.summary(io::IO, couplings::Vector{<:Coupling}) = @printf io "%s-element Vector{Coupling}" length(couplings)

"""
    Coupling(value, pattern::Pattern)
    Coupling(pattern::Pattern)

Construct a coupling.
"""
@inline Coupling(pattern::Pattern) = Coupling(1, pattern)

"""
    Coupling(index::Index, indexes::Index...)
    Coupling(value::Number, index::Index, indexes::Index...)
    Coupling(value::Number, indexes::OneAtLeast{Index})

Construct a coupling with the input indexes as the pattern.
"""
@inline Coupling(index::Index, indexes::Index...) = Coupling(1, index, indexes...)
@inline Coupling(value::Number, index::Index, indexes::Index...) = Coupling(value, (index, indexes...))
@inline Coupling(value::Number, indexes::OneAtLeast{Index}) = Coupling(value, Pattern(indexes))

"""
    Coupling{I}(sites::Union{NTuple{N, Ordinal}, Colon}, fields::Union{NTuple{N}, Colon}...) where {N, I<:InternalIndex}
    Coupling{I}(value::Number, sites::Union{NTuple{N, Ordinal}, Colon}, fields::Union{NTuple{N}, Colon}...) where {N, I<:InternalIndex}

Construct a `Coupling` with the input sites and the fields of a kind of internal index.
"""
@inline Coupling{I}(sites::Union{NTuple{N, Ordinal}, Colon}, fields::Union{NTuple{N}, Colon}...) where {N, I<:InternalIndex} = Coupling{I}(1, sites, fields...)
@inline function Coupling{I}(value::Number, sites::Union{NTuple{N, Ordinal}, Colon}, fields::Union{NTuple{N}, Colon}...) where {N, I<:InternalIndex}
    return Coupling(value, map(Index, default(sites, Val(N)), map(I, map(field->default(field, Val(N)), fields)...)))
end
@inline default(fields, ::Val) = fields
@inline default(::Colon, N::Val) = ntuple(i->:, N)

"""
    Coupling{F}(sites::Union{NTuple{N, Ordinal}, Colon}, fields::Union{NTuple{N}, Colon}...) where {N, F}
    Coupling{F}(value::Number, sites::Union{NTuple{N, Ordinal}, Colon}, fields::Union{NTuple{N}, Colon}...) where {N, F}

Construct a `Coupling` by a function that can construct an `Index` with the input sites and the fields of a kind of internal index.
"""
@inline Coupling{F}(sites::Union{NTuple{N, Ordinal}, Colon}, fields::Union{NTuple{N}, Colon}...) where {N, F} = Coupling{F}(1, sites, fields...)
@inline function Coupling{F}(value::Number, sites::Union{NTuple{N, Ordinal}, Colon}, fields::Union{NTuple{N}, Colon}...) where {N, F}
    return Coupling(value, map(F, default(sites, Val(N)), map(field->default(field, Val(N)), fields)...))
end

"""
    rank(coupling::Coupling) -> Int
    rank(::Type{M}) where {M<:Coupling} -> Int

Get the rank of a coupling.
"""
@inline rank(coupling::Coupling) = rank(typeof(coupling))
@inline rank(::Type{<:Coupling{V, P} where V}) where {P<:Pattern} = rank(P)

"""
    *(cp₁::Coupling, cp₂::Coupling) -> Coupling

Get the multiplication between two coupling.
"""
@inline Base.:*(cp₁::Coupling, cp₂::Coupling) = Coupling(cp₁.value*cp₂.value, cp₁.pattern⊗cp₂.pattern)

"""
    latexstring(coupling::Coupling) -> String

Convert a coupling to the latex format.
"""
@inline latexstring(coupling::Coupling) = @sprintf "%s%s" (coupling.value≈1 ? "" : coupling.value≈-1 ? "- " : string(str(coupling.value), "\\,")) latexstring(coupling.pattern)

"""
    expand(coupling::Coupling, ::Val{Rule}, bond::Bond, hilbert::Hilbert) where Rule

Expand a coupling with the given bond and Hilbert space under a given named pattern rule.
"""
function expand(coupling::Coupling, ::Val{Rule}, bond::Bond, hilbert::Hilbert) where Rule
    sites = patternrule(coupling.pattern.indexes.sites, Val(Rule), length(bond))
    indexes = patternrule(coupling.pattern.indexes.internals, Val(Rule))
    points = ntuple(i->bond[sites[i]], Val(rank(coupling)))
    internal = InternalProd(map((index, internal)->filter(index, internal), indexes, ntuple(i->hilbert[points[i].site], Val(rank(coupling)))))
    pattern = Pattern{partition(coupling.pattern)}(map(Index, sites, indexes), coupling.pattern.constraints, coupling.pattern.representations)
    return CouplingExpand(coupling.value, points, internal, pattern)
end
struct CouplingExpand{V, N, SV<:SVector, I<:InternalProd, P<:Pattern}
    value::V
    sites::NTuple{N, Int}
    rcoordinates::NTuple{N, SV}
    icoordinates::NTuple{N, SV}
    internal::I
    pattern::P
end
function CouplingExpand(value, points::NTuple{N, Point}, internal::InternalProd, pattern::Pattern) where N
    sites = ntuple(i->points[i].site, Val(N))
    rcoordinates = ntuple(i->points[i].rcoordinate, Val(N))
    icoordinates = ntuple(i->points[i].icoordinate, Val(N))
    return CouplingExpand(value, sites, rcoordinates, icoordinates, internal, pattern)
end
@inline Base.eltype(ex::CouplingExpand) = eltype(typeof(ex))
@inline @generated function Base.eltype(::Type{<:CouplingExpand{V, N, SV, I}}) where {V, N, SV<:SVector, I<:InternalProd}
    return Operator{V, Tuple{map(I->CoordinatedIndex{Index{I, Int}, SV}, fieldtypes(eltype(I)))...}}
end
@inline Base.IteratorSize(::Type{<:CouplingExpand}) = Base.SizeUnknown()
function Base.iterate(ex::CouplingExpand, state=(space=InternalIndexSpace(ex.pattern.indexes.internals, ex.internal); (space, iterate(space))))
    space, state = state
    while !isnothing(state)
        indexes, state = state
        state = iterate(space, state)
        if match(ex.pattern, indexes)
            return Operator(ex.value, map(CoordinatedIndex, map(Index, ex.sites, indexes), ex.rcoordinates, ex.icoordinates)), (space, state)
        end
    end
    return
end
struct InternalIndexSpace{I<:OneOrMore{InternalIndex}, E<:OneOrMore{InternalIndex}, V<:Internal{E}} <: VectorSpace{E}
    indexes::I
    internal::V
end
@inline VectorSpaceStyle(::Type{<:InternalIndexSpace{<:OneAtLeast{InternalIndex}, <:OneAtLeast{InternalIndex}, <:CompositeInternal}}) = VectorSpaceDirectProducted(:forward)
@inline function getcontent(space::InternalIndexSpace{<:OneAtLeast{InternalIndex}, <:OneAtLeast{InternalIndex}, <:CompositeInternal}, ::Val{:contents})
    return map((internal, index)->InternalIndexSpace(index, internal), space.internal.contents, space.indexes)
end
@inline function Base.convert(::Type{<:OneAtLeast{InternalIndex}}, index::CartesianIndex, space::InternalIndexSpace{<:OneAtLeast{InternalIndex}, <:OneAtLeast{InternalIndex}, <:CompositeInternal})
    return map(getindex, getcontent(space, :contents), index.I)
end
@inline VectorSpaceStyle(::Type{<:InternalIndexSpace{<:InternalIndex, <:InternalIndex, <:SimpleInternal}}) = VectorSpaceDirectProducted(:forward)
@inline shape(space::InternalIndexSpace{<:InternalIndex, <:InternalIndex, <:SimpleInternal}) = shape(space.internal, space.indexes)
@inline function Base.convert(::Type{<:CartesianIndex}, index::E, space::InternalIndexSpace{<:InternalIndex, E, <:SimpleInternal}) where {E<:InternalIndex}
    return convert(CartesianIndex, index, space.internal)
end
@inline function Base.convert(::Type{E}, index::CartesianIndex, space::InternalIndexSpace{<:InternalIndex, E, <:SimpleInternal}) where {E<:InternalIndex}
    return convert(E, index, space.internal)
end

# MatrixCoupling
"""
    MatrixCouplingComponent{T₁, T₂, V<:AbstractVector{T₁}} <: VectorSpace{Tuple{T₁, T₁, T₂}}

A component of a matrix coupling, i.e., a matrix acting on a separated internal space.
"""
struct MatrixCouplingComponent{T₁, T₂, V<:AbstractVector{T₁}} <: VectorSpace{Tuple{T₁, T₁, T₂}}
    left::V
    right::V
    matrix::SparseMatrixCSC{T₂, Int}
    function MatrixCouplingComponent(left::V, right::V, matrix::AbstractMatrix{T}) where {V<:AbstractVector, T}
        @assert length(left)==size(matrix)[1] && length(right)==size(matrix)[2] "MatrixCouplingComponent error: mismatched inputs."
        new{eltype(V), T, V}(left, right, convert(SparseMatrixCSC{T, Int}, matrix))
    end
end
@inline parameternames(::Type{<:MatrixCouplingComponent}) = (:basistype, :datatype, :basis)
@inline Base.length(component::MatrixCouplingComponent) = nnz(component.matrix)
@inline function Base.getindex(component::MatrixCouplingComponent, i::Integer)
    row = component.matrix.rowval[i]
    col = searchsortedlast(component.matrix.colptr, i)
    return (component.left[row], component.right[col], component.matrix.nzval[i])
end
@inline Base.promote_rule(::Type{MatrixCouplingComponent{T, T₁, V}}, ::Type{MatrixCouplingComponent{T, T₂, V}}) where {T, T₁, T₂, V<:AbstractVector{T}} = MatrixCouplingComponent{T, promote_type(T₁, T₂), V}
@inline function Base.convert(::Type{MatrixCouplingComponent{T, T₁, V}}, component::MatrixCouplingComponent{T, T₂, V}) where {T, T₁, T₂, V<:AbstractVector{T}}
    return MatrixCouplingComponent(component.left, component.right, convert(SparseMatrixCSC{T₁, Int}, component.matrix))
end

"""
    MatrixCoupling{I<:InternalIndex, S<:Union{Ordinal, Colon}, C<:OneAtLeast{MatrixCouplingComponent}, E<:Coupling} <: VectorSpace{E}

Matrix coupling, i.e., a set of couplings whose coefficients are specified by matrices acting on separated internal spaces.
"""
struct MatrixCoupling{I<:InternalIndex, S<:Union{Ordinal, Colon}, C<:OneAtLeast{MatrixCouplingComponent}, E<:Coupling} <: VectorSpace{E}
    sites::Tuple{S, S}
    contents::C
    function MatrixCoupling{I}(sites::Union{NTuple{2, Ordinal}, NTuple{2, Colon}, Colon}, contents::OneAtLeast{MatrixCouplingComponent}) where {I<:InternalIndex}
        @assert fieldcount(I)==length(contents) "MatrixCoupling error: mismatched type of internal index ($nameof(I)) and components (len=$length(contents))."
        sites = default(sites, Val(2))
        new{I, eltype(sites), typeof(contents), _eltype_(I, eltype(sites), typeof(contents))}(sites, contents)
    end
end
@inline parameternames(::Type{<:MatrixCoupling}) = (:internal, :site, :components)
@inline @generated function _eltype_(::Type{I}, ::Type{S}, ::Type{C}) where {I<:InternalIndex, S<:Union{Ordinal, Colon}, C<:OneAtLeast{MatrixCouplingComponent}}
    types = fieldtypes(C)
    V = Expr(:call, :promote_type, [:(parametertype($C, :datatype)) for C in types]...)
    I = Expr(:call, :internalindextype, I, [:(parametertype($C, :basistype)) for C in types]...)
    return :(Coupling{$V, Pattern{NTuple{2, Index{$I, $S}}, (2,), 1, Tuple{typeof(isdiagonal)}}})
end
@inline VectorSpaceStyle(::Type{<:MatrixCoupling}) = VectorSpaceDirectProducted(:forward)
function Base.convert(::Type{<:Coupling}, index::CartesianIndex, mc::MatrixCoupling{I}) where {I<:InternalIndex}
    contents = map(getindex, getcontent(mc, :contents), index.I)
    value = mapreduce(content->content[3], *, contents, init=1)
    index₁ = Index(mc.sites[1], I(map(content->content[1], contents)...))
    index₂ = Index(mc.sites[2], I(map(content->content[2], contents)...))
    return Coupling(value, index₁, index₂)
end
@inline function Base.promote_rule(::Type{<:MatrixCoupling{I, S, C₁}}, ::Type{<:MatrixCoupling{I, S, C₂}}) where {I<:InternalIndex, S<:Union{Ordinal, Colon}, C₁<:OneAtLeast{MatrixCouplingComponent}, C₂<:OneAtLeast{MatrixCouplingComponent}}
    C = _promote_type_(C₁, C₂)
    return MatrixCoupling{I, S, C, _eltype_(I, S, C)}
end
@inline Base.convert(::Type{<:MatrixCoupling{I, S, C}}, mc::MatrixCoupling{I, S}) where {I<:InternalIndex, S<:Union{Ordinal, Colon}, C<:OneAtLeast{MatrixCouplingComponent}} = MatrixCoupling{I}(mc.sites, convert(C, mc.contents))
@inline @generated function _promote_type_(::Type{T₁}, ::Type{T₂}) where {N, T₁<:NTuple{N, Any}, T₂<:NTuple{N, Any}}
    return Expr(:curly, Tuple, [:(promote_type(fieldtype(T₁, $i), fieldtype(T₂, $i))) for i=1:N]...)
end

"""
    MatrixCoupling{I}(sites::Union{NTuple{2, Ordinal}, NTuple{2, Colon}, Colon}, contents::MatrixCouplingComponent...) where {I<:InternalIndex}
    MatrixCoupling{I}(sites::Union{NTuple{2, Ordinal}, NTuple{2, Colon}, Colon}, contents::OneAtLeast{MatrixCouplingComponent}) where {I<:InternalIndex}

Construct a `MatrixCoupling`.
"""
@inline MatrixCoupling{I}(sites::Union{NTuple{2, Ordinal}, NTuple{2, Colon}, Colon}, contents::MatrixCouplingComponent...) where {I<:InternalIndex} = MatrixCoupling{I}(sites, contents)

"""
    MatrixCouplingProd{V<:Number, C<:OneAtLeast{MatrixCoupling}, E<:Coupling} <: VectorSpace{E}

Product of matrix couplings together with an overall coefficient.
"""
struct MatrixCouplingProd{V<:Number, C<:OneAtLeast{MatrixCoupling}, E<:Coupling} <: VectorSpace{E}
    value::V
    contents::C
    function MatrixCouplingProd(value::Number, contents::OneAtLeast{MatrixCoupling})
        new{typeof(value), typeof(contents), _eltype_(typeof(value), _eltypes_(typeof(contents)))}(value, contents)
    end
end
@inline @generated _eltypes_(::Type{TS}) where {TS<:OneAtLeast{MatrixCoupling}} = Expr(:curly, :Tuple, [:(eltype($C)) for C in fieldtypes(TS)]...)
@inline @generated function _eltype_(::Type{V}, ::Type{TS}) where {V<:Number, TS<:Tuple}
    types = fieldtypes(TS)
    MV = promote_type(V, map(valtype, types)...)
    N = length(types)
    RS = ntuple(i->2, Val(N))
    IS = Tuple{concatenate(map(C->fieldtypes(parametertype(fieldtype(C, :pattern), :indexes)), types)...)...}
    FS = Tuple{concatenate(map(C->fieldtypes(parametertype(fieldtype(C, :pattern), :constraints)), types)...)...}
    return Coupling{MV, Pattern{IS, RS, N, FS}}
end
@inline VectorSpaceStyle(::Type{<:MatrixCouplingProd}) = VectorSpaceDirectProducted(:forward)
@inline function Base.convert(::Type{<:Coupling}, index::CartesianIndex, mcp::MatrixCouplingProd)
    contents = map(getindex, getcontent(mcp, :contents), index.I)
    return prod(contents; init=mcp.value)
end
@inline function Base.promote_rule(::Type{<:MatrixCouplingProd{V₁, C₁}}, ::Type{<:MatrixCouplingProd{V₂, C₂}}) where {V₁<:Number, C₁<:OneAtLeast{MatrixCoupling}, V₂<:Number, C₂<:OneAtLeast{MatrixCoupling}}
    V = promote_type(V₁, V₂)
    C = _promote_type_(C₁, C₂)
    return MatrixCouplingProd{V, C, _eltype_(V, _eltypes_(C))}
end
@inline Base.convert(::Type{<:MatrixCouplingProd{V, C}}, mcp::MatrixCouplingProd) where {V<:Number, C<:OneAtLeast{MatrixCoupling}} = MatrixCouplingProd(convert(V, mcp.value), convert(C, mcp.contents))

"""
    MatrixCouplingProd(value::Number, contents::OneAtLeast{MatrixCoupling})
    MatrixCouplingProd(content::MatrixCoupling, contents::MatrixCoupling...)
    MatrixCouplingProd(value::Number, content::MatrixCoupling, contents::MatrixCoupling...)

Construct a `MatrixCouplingProd`.
"""
@inline MatrixCouplingProd(content::MatrixCoupling, contents::MatrixCoupling...) = MatrixCouplingProd(1, content, contents...)
@inline MatrixCouplingProd(value::Number, content::MatrixCoupling, contents::MatrixCoupling...) = MatrixCouplingProd(value, (content, contents...))

"""
    MatrixCouplingSum{C<:MatrixCouplingProd, N, E<:Coupling} <: VectorSpace{E}

Sum of the products of matrix couplings.
"""
struct MatrixCouplingSum{C<:MatrixCouplingProd, N, E<:Coupling} <: VectorSpace{E}
    contents::NTuple{N, C}
    function MatrixCouplingSum(contents::OneAtLeast{MatrixCouplingProd})
        contents = promote(contents...)
        new{eltype(contents), fieldcount(typeof(contents)), eltype(eltype(contents))}(contents)
    end
end
@inline VectorSpaceStyle(::Type{<:MatrixCouplingSum}) = VectorSpaceDirectSummed()

"""
    MatrixCouplingSum(contents::OneAtLeast{Union{MatrixCoupling, MatrixCouplingProd}})
    MatrixCouplingSum(content::Union{MatrixCoupling, MatrixCouplingProd}, contents::Union{MatrixCoupling, MatrixCouplingProd}...)

Construct a `MatrixCouplingSum`.
"""
@inline MatrixCouplingSum(content::Union{MatrixCoupling, MatrixCouplingProd}, contents::Union{MatrixCoupling, MatrixCouplingProd}...) = MatrixCouplingSum((content, contents...))
@inline MatrixCouplingSum(contents::OneAtLeast{Union{MatrixCoupling, MatrixCouplingProd}}) = MatrixCouplingSum(promote(map(m->1*m, contents)...))

"""
    *(mc₁::MatrixCoupling, mc₂::MatrixCoupling) -> MatrixCouplingProd
    *(mc::MatrixCoupling, mcp::MatrixCouplingProd) -> MatrixCouplingProd
    *(mcp::MatrixCouplingProd, mc::MatrixCoupling) -> MatrixCouplingProd
    *(mcp₁::MatrixCouplingProd, mcp₂::MatrixCouplingProd) -> MatrixCouplingProd
    *(mc::MatrixCoupling, factor::Number) -> MatrixCouplingProd
    *(factor::Number, mc::MatrixCoupling) -> MatrixCouplingProd
    *(factor::Number, mcp::MatrixCouplingProd) -> MatrixCouplingProd
    *(mcp::MatrixCouplingProd, factor::Number) -> MatrixCouplingProd
    *(mcs::MatrixCouplingSum, element::Union{Number, MatrixCoupling, MatrixCouplingProd}) -> MatrixCouplingSum
    *(element::Union{Number, MatrixCoupling, MatrixCouplingProd}, mcs::MatrixCouplingSum) -> MatrixCouplingSum
    *(mcs₁::MatrixCouplingSum, mcs₂::MatrixCouplingSum) -> MatrixCouplingSum

Product between `MatrixCoupling`s and `MatrixCouplingProd`s.
"""
@inline Base.:*(mc₁::MatrixCoupling, mc₂::MatrixCoupling) = MatrixCouplingProd(mc₁, mc₂)
@inline Base.:*(mc::MatrixCoupling, mcp::MatrixCouplingProd) = MatrixCouplingProd(mcp.value, mc, mcp.contents...)
@inline Base.:*(mcp::MatrixCouplingProd, mc::MatrixCoupling) = MatrixCouplingProd(mcp.value, mcp.contents..., mc)
@inline Base.:*(mcp₁::MatrixCouplingProd, mcp₂::MatrixCouplingProd) = MatrixCouplingProd(mcp₁.value*mcp₂.value, mcp₁.contents..., mcp₂.contents...)
@inline Base.:*(factor::Number, mc::MatrixCoupling) = MatrixCouplingProd(factor, mc)
@inline Base.:*(mc::MatrixCoupling, factor::Number) = MatrixCouplingProd(factor, mc)
@inline Base.:*(factor::Number, mcp::MatrixCouplingProd) = MatrixCouplingProd(factor*mcp.value, mcp.contents)
@inline Base.:*(mcp::MatrixCouplingProd, factor::Number) = MatrixCouplingProd(factor*mcp.value, mcp.contents)
@inline Base.:*(mcs::MatrixCouplingSum, factor::Number) = MatrixCouplingSum(map(m->m*factor, mcs.contents))
@inline Base.:*(factor::Number, mcs::MatrixCouplingSum) = MatrixCouplingSum(map(m->m*factor, mcs.contents))
@inline Base.:*(mcs::MatrixCouplingSum, element::Union{MatrixCoupling, MatrixCouplingProd}) = MatrixCouplingSum(map(m->m*element, mcs.contents))
@inline Base.:*(element::Union{MatrixCoupling, MatrixCouplingProd}, mcs::MatrixCouplingSum) = MatrixCouplingSum(map(m->element*m, mcs.contents))
@inline Base.:*(mcs₁::MatrixCouplingSum, mcs₂::MatrixCouplingSum) = MatrixCouplingSum(concatenate(map(m₁->map(m₂->m₁*m₂, mcs₂.contents), mcs₁.contents)...))

"""
    +(mc₁::Union{MatrixCoupling, MatrixCouplingProd}, mc₂::Union{MatrixCoupling, MatrixCouplingProd}) -> MatrixCouplingSum
    +(mc::Union{MatrixCoupling, MatrixCouplingProd}, mcs::MatrixCouplingSum) -> MatrixCouplingSum
    +(mcs::MatrixCouplingSum, mc::Union{MatrixCoupling, MatrixCouplingProd}) -> MatrixCouplingSum
    +(mcs₁::MatrixCouplingSum, mcs₂::MatrixCouplingSum) -> MatrixCouplingSum

Addition between `MatrixCoupling`s and `MatrixCouplingProd`s.
"""
@inline Base.:+(mc₁::Union{MatrixCoupling, MatrixCouplingProd}, mc₂::Union{MatrixCoupling, MatrixCouplingProd}) = MatrixCouplingSum(mc₁, mc₂)
@inline Base.:+(mc::Union{MatrixCoupling, MatrixCouplingProd}, mcs::MatrixCouplingSum) = MatrixCouplingSum(mc, mcs.contents...)
@inline Base.:+(mcs::MatrixCouplingSum, mc::Union{MatrixCoupling, MatrixCouplingProd}) = MatrixCouplingSum(mcs.contents..., mc)
@inline Base.:+(mcs₁::MatrixCouplingSum, mcs₂::MatrixCouplingSum) = MatrixCouplingSum(mcs₁.contents..., mcs₂.contents...)

"""
    ^(mc::Union{MatrixCoupling, MatrixCouplingProd, MatrixCouplingSum}, n::Integer) -> Union{MatrixCoupling, MatrixCouplingProd, MatrixCouplingSum}

Get the nth power of a `MatrixCoupling`/`MatrixCouplingProd`/`MatrixCouplingSum`.
"""
@inline Base.:^(mc::Union{MatrixCoupling, MatrixCouplingProd, MatrixCouplingSum}, n::Integer) = prod(ntuple(i->mc, Val(n)); init=1)

"""
    /(mcp::MatrixCouplingProd, factor::Number) -> MatrixCouplingProd
    /(mcs::MatrixCouplingSum, factor::Number) -> MatrixCouplingSum
    /(mc::MatrixCoupling, factor::Number) -> MatrixCouplingProd
    //(mcp::MatrixCouplingProd, factor::Number) -> MatrixCouplingProd
    //(mcs::MatrixCouplingSum, factor::Number) -> MatrixCouplingSum
    //(mc::MatrixCoupling, factor::Number) -> MatrixCouplingProd
    -(mc::MatrixCoupling) -> MatrixCouplingProd
    -(mcp::MatrixCouplingProd) -> MatrixCouplingProd
    -(mcs::MatrixCouplingSum) -> MatrixCouplingSum
    -(mc₁::Union{MatrixCoupling, MatrixCouplingProd}, mc₂::Union{MatrixCoupling, MatrixCouplingProd}) -> MatrixCouplingSum
    -(mc::Union{MatrixCoupling, MatrixCouplingProd}, mcs::MatrixCouplingSum) -> MatrixCouplingSum
    -(mcs::MatrixCouplingSum, mc::Union{MatrixCoupling, MatrixCouplingProd}) -> MatrixCouplingSum
    -(mcs₁::MatrixCouplingSum, mcs₂::MatrixCouplingSum) -> MatrixCouplingSum

Define right-division, minus and subtraction operator for a `MatrixCoupling`/`MatrixCouplingProd`/`MatrixCouplingSum`.
"""
@inline Base.:/(mcp::MatrixCouplingProd, factor::Number) = MatrixCouplingProd(mcp.value/factor, mcp.contents)
@inline Base.:/(mcs::MatrixCouplingSum, factor::Number) = MatrixCouplingSum(map(m->m/factor, mcs.contents))
@inline Base.:/(mc::MatrixCoupling, factor::Number) = MatrixCouplingProd(1/factor, mc)
@inline Base.://(mcp::MatrixCouplingProd, factor::Number) = MatrixCouplingProd(mcp.value//factor, mcp.contents)
@inline Base.://(mcs::MatrixCouplingSum, factor::Number) = MatrixCouplingSum(map(m->m//factor, mcs.contents))
@inline Base.://(mc::MatrixCoupling, factor::Number) = MatrixCouplingProd(1//factor, mc)
@inline Base.:-(mc::MatrixCoupling) = MatrixCouplingProd(-1, mc)
@inline Base.:-(mcp::MatrixCouplingProd) = MatrixCouplingProd(-1*mcp.value, mcp.contents)
@inline Base.:-(mcs::MatrixCouplingSum) = MatrixCouplingSum(map(m->-m, mcs.contents))
@inline Base.:-(mc₁::Union{MatrixCoupling, MatrixCouplingProd}, mc₂::Union{MatrixCoupling, MatrixCouplingProd}) = MatrixCouplingSum(mc₁, -mc₂)
@inline Base.:-(mc::Union{MatrixCoupling, MatrixCouplingProd}, mcs::MatrixCouplingSum) = MatrixCouplingSum(mc, map(m->-m, mcs.contents)...)
@inline Base.:-(mcs::MatrixCouplingSum, mc::Union{MatrixCoupling, MatrixCouplingProd}) = MatrixCouplingSum(mcs.contents..., -mc)
@inline Base.:-(mcs₁::MatrixCouplingSum, mcs₂::MatrixCouplingSum) = MatrixCouplingSum(mcs₁.contents..., map(m->-m, mcs₂.contents)...)

# Term
"""
    TermAmplitude{F} <: Function

Function for the amplitude of a term.
"""
struct TermAmplitude{F} <: Function
    TermAmplitude(amplitude::Union{Function, Nothing}) = new{amplitude}()
end
@inline Base.:(==)(::TermAmplitude{F₁}, ::TermAmplitude{F₂}) where {F₁, F₂} = F₁==F₂
@inline Base.isequal(::TermAmplitude{F₁}, ::TermAmplitude{F₂}) where {F₁, F₂} = isequal(F₁, F₂)
@inline (::TermAmplitude{nothing})(::Bond) = 1
@inline (::TermAmplitude{F})(bond::Bond) where F = F(bond)
@inline Base.valtype(termamplitude::TermAmplitude, bond::Bond) = valtype(typeof(termamplitude), typeof(bond))
@inline Base.valtype(::Type{TermAmplitude{nothing}}, ::Type{<:Bond}) = Int
@inline Base.valtype(::Type{TermAmplitude{F}}, ::Type{B}) where {F, B<:Bond} = Core.Compiler.return_type(F, Tuple{B})

"""
    TermCoupling{C<:Coupling, F} <: Function

Function for the coupling of a term.
"""
struct TermCoupling{C<:Coupling, F} <: Function
    coupling::F
    TermCoupling(coupling) = new{eltype(coupling), typeof(coupling)}(coupling)
    TermCoupling(coupling::Function) = new{eltype(Core.Compiler.return_type(coupling, Tuple{Bond})), typeof(coupling)}(coupling)
    TermCoupling{C}(coupling::Function) where {C<:Coupling} = new{C, typeof(coupling)}(coupling)
end
@inline Base.:(==)(termcoupling₁::TermCoupling, termcoupling₂::TermCoupling) = ==(termcoupling₁.coupling, termcoupling₂.coupling)
@inline Base.isequal(termcoupling₁::TermCoupling, termcoupling₂::TermCoupling) = isequal(termcoupling₁.coupling, termcoupling₂.coupling)
@inline Base.valtype(termcoupling::TermCoupling) = valtype(typeof(termcoupling))
@inline Base.valtype(::Type{<:TermCoupling{C}}) where {C<:Coupling} = C
@inline (termcoupling::TermCoupling)(::Bond) = termcoupling.coupling
@inline (termcoupling::TermCoupling{<:Coupling, <:Function})(bond::Bond) = termcoupling.coupling(bond)

"""
    Term{K, I, V, B, C<:TermCoupling, A<:TermAmplitude}

Term of a quantum lattice system.
"""
mutable struct Term{K, I, V, B, C<:TermCoupling, A<:TermAmplitude}
    value::V
    const bondkind::B
    const coupling::C
    const amplitude::A
    const ishermitian::Bool
    const ismodulatable::Bool
    const factor::V
    function Term{K, I}(value, bondkind, coupling::TermCoupling, amplitude::TermAmplitude, ishermitian::Bool, ismodulatable::Bool, factor) where {K, I}
        @assert isa(K, Symbol) "Term error: kind must be a Symbol."
        @assert isa(I, Symbol) "Term error: id must be a Symbol."
        (isa(bondkind, Number) && iszero(bondkind) && !ishermitian) || @assert(
            value==value',
            "Term error: real value required. For an Hermitian term, the value must be real. So is the case for a term beyond onsite even in the on-Hermitian situation because a complex value always has the positive direction and it should be specified by the amplitude function.")
        new{K, I, typeof(value), typeof(bondkind), typeof(coupling), typeof(amplitude)}(value, bondkind, coupling, amplitude, ishermitian, ismodulatable, factor)
    end
end
@inline Base.:(==)(term₁::Term, term₂::Term) = ==(efficientoperations, term₁, term₂)
@inline Base.isequal(term₁::Term, term₂::Term) = isequal(efficientoperations, term₁, term₂)

"""
    Term{K}(id::Symbol, value, bondkind, coupling, ishermitian::Bool; amplitude::Union{Function, Nothing}=nothing, ismodulatable::Bool=true) where K

Construct a term.
"""
@inline function Term{K}(id::Symbol, value, bondkind, coupling, ishermitian::Bool; amplitude::Union{Function, Nothing}=nothing, ismodulatable::Bool=true) where K
    return Term{K, id}(value, bondkind, TermCoupling(coupling), TermAmplitude(amplitude), ishermitian, ismodulatable, 1)
end

"""
    kind(term::Term) -> Symbol
    kind(::Type{<:Term) -> Symbol

Get the kind of a term.
"""
@inline kind(term::Term) = kind(typeof(term))
@inline kind(::Type{<:Term{K}}) where K = K

"""
    id(term::Term) -> Symbol
    id(::Type{<:Term) -> Symbol

Get the id of a term.
"""
@inline id(term::Term) = id(typeof(term))
@inline id(::Type{<:Term{K, I} where K}) where I = I

"""
    value(term::Term) -> valtype(term)

Get the value of a term.
"""
@inline value(term::Term) = term.value

"""
    valtype(term::Term)
    valtype(::Type{<:Term)

Get the value type of a term.
"""
@inline Base.valtype(term::Term) = valtype(typeof(term))
@inline Base.valtype(::Type{<:Term{K, I, V} where {K, I}}) where V = V

"""
    valtype(terms::OneAtLeast{Term})
    valtype(::Type{<:T}) where {T<:OneAtLeast{Term}}

Get the common value type of a set of terms.
"""
@inline Base.valtype(terms::OneAtLeast{Term}) = valtype(typeof(terms))
@inline @generated Base.valtype(::Type{<:T}) where {T<:OneAtLeast{Term}} = mapreduce(valtype, promote_type, fieldtypes(T))

"""
    rank(term::Term) -> Int
    rank(::Type{<:Term) -> Int

Get the rank of a term.
"""
@inline rank(term::Term) = rank(typeof(term))
@inline rank(::Type{<:Term{K, I, V, B, C} where {K, I, V, B}}) where {C<:TermCoupling} = rank(valtype(C))

"""
    string(term::Term, bond::Bond, hilbert::Hilbert) -> String

Get the string representation of a term on a bond with a given Hilbert space.
"""
function Base.string(term::Term, bond::Bond, hilbert::Hilbert)
    cache = String[]
    if term.bondkind == bond.kind
        value = term.value * term.amplitude(bond) * term.factor
        if !isapprox(value, 0.0, atol=atol, rtol=rtol)
            for coupling in term.coupling(bond)
                if !isnothing(iterate(expand(coupling, Val(kind(term)), bond, hilbert)))
                    representation = string(value * coupling)
                    term.ishermitian || (representation = string(representation, " + h.c."))
                    push!(cache, @sprintf "%s" representation)
                end
            end
        end
    end
    return join(cache, "\n")
end

"""
    replace(term::Term, value) -> Term

Replace the value of a term.
"""
@inline function Base.replace(term::Term, value)
    return Term{kind(term), id(term)}(value, term.bondkind, term.coupling, term.amplitude, term.ishermitian, term.ismodulatable, convert(typeof(value), term.factor))
end

"""
    one(term::Term) -> Term

Get a unit term.
"""
@inline Base.one(term::Term) = replace(term, one(term.value))

"""
    zero(term::Term) -> Term

Get a zero term.
"""
@inline Base.zero(term::Term) = replace(term, zero(term.value))

"""
    update!(term::Term, args...; kwargs...) -> Term

Update the value of a term if it is ismodulatable.
"""
function update!(term::Term, args...; kwargs...)
    @assert term.ismodulatable "update! error: not modulatable term."
    term.value = get(kwargs, id(term), term.value)
    return term
end

"""
    operatortype(::Type{B}, ::Type{H}, ::Type{T}) where {B<:Bond, H<:Hilbert, T<:Term}

Get the compatible `Operator` type from the type of a term, a Hilbert space and a bond.
"""
@inline function operatortype(::Type{B}, ::Type{H}, ::Type{T}) where {B<:Bond, H<:Hilbert, T<:Term}
    C = valtype(fieldtype(T, :coupling))
    @assert C<:Coupling "operatortype error: not supported."
    V, V′, V′′ = valtype(T), valtype(C), valtype(fieldtype(T, :amplitude), B)
    isconcretetype(V′) && (V = promote_type(V, V′))
    isconcretetype(V′′) && (V = promote_type(V, V′′))
    coordinatedindextypes = ntuple(i->coordinatedindextype(filter(internalindextype(fieldtype(fieldtype(fieldtype(C, :pattern), :indexes), i)), valtype(H)), eltype(B)), Val(rank(C)))
    return fulltype(Operator, NamedTuple{(:value, :id), Tuple{V, Tuple{coordinatedindextypes...}}})
end

"""
    expand!(operators::Operators, term::Term, bond::Bond, hilbert::Hilbert; half::Bool=false) -> Operators
    expand!(operators::Operators, term::Term, bonds, hilbert::Hilbert; half::Bool=false) -> Operators

Expand the operators of a term on a bond/set-of-bonds with a given Hilbert space.

The `half` parameter determines the behavior of generating operators, which falls into the following two categories
* `true`: "Hermitian half" of the generated operators
* `false`: "Hermitian whole" of the generated operators
"""
function expand!(operators::Operators, term::Term, bond::Bond, hilbert::Hilbert; half::Bool=false)
    if term.bondkind == bond.kind
        value = term.value * term.amplitude(bond) * term.factor
        if !isapprox(value, 0, atol=atol, rtol=rtol)
            for coupling in term.coupling(bond)
                for opt in expand(coupling, Val(kind(term)), bond, hilbert)
                    isapprox(opt.value, 0, atol=atol, rtol=rtol) && continue
                    if half
                        add!(operators, opt*scalartype(operators)(value/(term.ishermitian ? 2 : 1)))
                    else
                        opt = opt*scalartype(operators)(value)
                        add!(operators, opt)
                        term.ishermitian || add!(operators, opt')
                    end
                end
            end
        end
    end
    return operators
end
@inline function expand!(operators::Operators, term::Term, bonds, hilbert::Hilbert; half::Bool=false)
    for bond in bonds
        expand!(operators, term, bond, hilbert; half=half)
    end
    return operators
end

"""
    expand(term::Term, bond::Bond, hilbert::Hilbert; half::Bool=false) -> Operators
    expand(term::Term, bonds, hilbert::Hilbert; half::Bool=false) -> Operators

Expand the operators of a term on a bond/set-of-bonds with a given Hilbert space.
"""
@inline function expand(term::Term, bond::Bond, hilbert::Hilbert; half::Bool=false)
    M = operatortype(bond|>typeof, hilbert|>typeof, term|>typeof)
    expand!(Operators{M}(), term, bond, hilbert; half=half)
end
@inline function expand(term::Term, bonds, hilbert::Hilbert; half::Bool=false)
    M = operatortype(bonds|>eltype, hilbert|>typeof, term|>typeof)
    expand!(Operators{M}(), term, bonds, hilbert; half=half)
end

"""
    nneighbor(term::Term) -> Int
    nneighbor(terms::OneAtLeast{Term}) -> Int

Get the
1) order of neighbor in a single term;
2) highest order of neighbors in a set of terms.
"""
@inline nneighbor(term::Term) = term.bondkind::Int
@inline nneighbor(terms::OneAtLeast{Term}) = maximum(nneighbor, terms)::Int

# Metric and Table
"""
    Metric <: Function

Rules for measuring an operator index so that different operator indexes can be compared.

As a function, every instance should accept only one positional argument, i.e. the operator index to be measured.
"""
abstract type Metric <: Function end
@inline Base.:(==)(m₁::T, m₂::T) where {T<:Metric} = ==(efficientoperations, m₁, m₂)
@inline Base.isequal(m₁::T, m₂::T) where {T<:Metric} = isequal(efficientoperations, m₁, m₂)
@inline (M::Type{<:Metric})(::Type{I}) where {I<:CompositeIndex} = M(indextype(I))
@inline (M::Type{<:Metric})(::Type{H}) where {H<:Hilbert} = M(Index{H|>valtype|>eltype, Int})
@inline (metric::Metric)(index::CompositeIndex) = metric(getcontent(index, :index))
@inline Base.valtype(::Type{M}, ::Type{I}) where {M<:Metric, I<:CompositeIndex} = valtype(M, indextype(I))

"""
    OperatorIndexToTuple{Fields} <: Metric

A rule that converts an operator index into a tuple based on the specified type parameter `Fields`.

Here, `Fields` must be a tuple of `Union{Symbol, Function}`, which determines the elements of the converted tuple on an element-by-element basis.

For the ith element of `Fields`:

1) If it is a Symbol, it represents the name of a single index of an `OperatorIndex`, and its value will become the corresponding element in the converted tuple.
2) If it is a Function, it should be a trait function of an `OperatorIndex`, and its return value will become the corresponding element in the converted tuple.
"""
struct OperatorIndexToTuple{Fields} <: Metric
    OperatorIndexToTuple(fields::ZeroAtLeast{Union{Symbol, Function}}) = new{fields}()
end
@inline OperatorIndexToTuple(fields::Union{Symbol, Function}...) = OperatorIndexToTuple(fields)

"""
    keys(::OperatorIndexToTuple{Fields}) where Fields -> Fields
    keys(::Type{<:OperatorIndexToTuple{Fields}}) where Fields -> Fields

Get the values of the type parameter `Fields`.
"""
@inline Base.keys(::OperatorIndexToTuple{Fields}) where Fields = Fields
@inline Base.keys(::Type{<:OperatorIndexToTuple{Fields}}) where Fields = Fields

"""
    OperatorIndexToTuple(::Type{I}) where {I<:Index}

Construct the metric rule from the information of the `Index` type.
"""
@inline OperatorIndexToTuple(::Type{I}) where {I<:Index} = OperatorIndexToTuple(:site, (fieldnames(internalindextype(I)))...)

"""
    valtype(::Type{<:OperatorIndexToTuple}, ::Type{<:Index})

Get the valtype of applying an `OperatorIndexToTuple` rule to an `Index`.
"""
@inline @generated function Base.valtype(::Type{M}, ::Type{I}) where {M<:OperatorIndexToTuple, I<:Index}
    types = []
    for field in keys(M)
        if isa(field, Function)
            push!(types, :(Core.Compiler.return_type($field, Tuple{Type{I}})))
        elseif field==:site
            push!(types, Int)
        elseif hasfield(internalindextype(I), field)
            push!(types, fieldtype(internalindextype(I), field))
        end
    end
    return  Expr(:curly, :Tuple, types...)
end

"""
    (operatorunittotuple::OperatorIndexToTuple)(index::Index) -> Tuple

Convert an index to a tuple.
"""
@inline @generated function (operatorunittotuple::OperatorIndexToTuple)(index::Index)
    exprs = []
    for name in keys(operatorunittotuple)
        if isa(name, Function)
            push!(exprs, :(($name)(typeof(index))))
        elseif name==:site
            push!(exprs, :(index.site))
        elseif hasfield(internalindextype(index), name)
            push!(exprs, :(getfield(InternalIndex(index), $(QuoteNode(name)))))
        end
    end
    return Expr(:tuple, exprs...)
end

"""
    Table{I, B<:Metric} <: CompositeDict{I, Int}

Table of operator index v.s. sequence pairs.
"""
struct Table{I, B<:Metric} <: CompositeDict{I, Int}
    by::B
    contents::OrderedDict{I, Int}
end
@inline contentnames(::Type{<:Table}) = (:by, :contents)
@inline Table{I}(by::Metric) where {I<:OperatorIndex} = Table(by, OrderedDict{valtype(typeof(by), I), Int}())
@inline vec2dict(vs::AbstractVector) = OrderedDict{eltype(vs), Int}(v=>i for (i, v) in enumerate(vs))

"""
    getindex(table::Table, index::OperatorIndex) -> Int

Inquiry the sequence of an operator index.
"""
@inline Base.getindex(table::Table, index::OperatorIndex) = table[convert(keytype(table), table.by(index))]

"""
    haskey(table::Table, index::OperatorIndex) -> Bool
    haskey(table::Table, indexes::ZeroAtLeast{OperatorIndex}) -> ZeroAtLeast{Bool}

Judge whether a single operator index or a set of operator indexes have been assigned with sequences in table.
"""
@inline Base.haskey(table::Table, index::OperatorIndex) = haskey(table, table.by(index))
@inline Base.haskey(table::Table, indexes::ZeroAtLeast{OperatorIndex}) = map(index->haskey(table, index), indexes)

"""
    Table(indexes::AbstractVector{<:OperatorIndex}, by::Metric=OperatorIndexToTuple(eltype(indexes)))

Convert a set of operator units to the corresponding table of operator index vs. sequence pairs.

The input operator units are measured by the input `by` function with the duplicates removed. The resulting unique values are sorted, which determines the sequence of the input `indexes`. Note that two operator units have the same sequence if their converted values are equal to each other.
"""
@inline Table(indexes::AbstractVector{<:OperatorIndex}, by::Metric=OperatorIndexToTuple(eltype(indexes))) = Table(by, [by(index) for index in indexes]|>unique!|>sort!|>vec2dict)

"""
    Table(hilbert::Hilbert, by::Metric=OperatorIndexToTuple(typeof(hilbert))) -> Table

Get the index-sequence table of a Hilbert space.
"""
function Table(hilbert::Hilbert, by::Metric=OperatorIndexToTuple(typeof(hilbert)))
    result = fulltype(Index, Tuple{hilbert|>valtype|>eltype, Int})[]
    for (site, internal) in hilbert
        for index in internal
            push!(result, Index(site, index))
        end
    end
    return Table(result, by)
end

"""
    union(tables::Table...) -> Table

Unite several operator index vs. sequence tables.
"""
function Base.union(tables::Table...)
    @assert allequal(map(table->table.by, tables)) "union error: all input tables should have the same `by` attribute."
    indices = (tables|>eltype|>keytype)[]
    for table in tables
        for index in keys(table)
            push!(indices, index)
        end
    end
    return Table(tables[1].by, indices|>unique!|>sort!|>vec2dict)
end

"""
    reset!(table::Table, indexes::AbstractVector{<:OperatorIndex}) -> Table

Reset a table by a new set of indexes.
"""
function reset!(table::Table, indexes::AbstractVector{<:OperatorIndex})
    empty!(table)
    for (i, id) in enumerate([table.by(index) for index in indexes]|>unique!|>sort!)
        table[id] = i
    end
    return table
end

"""
    reset!(table::Table, hilbert::Hilbert) -> Table

Reset a table by a Hilbert space.
"""
function reset!(table::Table, hilbert::Hilbert)
    indices = Index{hilbert|>valtype|>eltype, Int}[]
    for (site, internal) in pairs(hilbert)
        for internalindex in internal
            push!(indices, (indices|>eltype)(site, internalindex))
        end
    end
    return reset!(table, indices)
end

"""
    findall(select::Function, hilbert::Hilbert, table::Table) -> Vector{Int}

Find all the sequences of indexes contained in a Hilbert space according to a table and a select function.
"""
function Base.findall(select::Function, hilbert::Hilbert, table::Table)
    result = Int[]
    for (site, internal) in pairs(hilbert)
        for internalindex in internal
            index = Index(site, internalindex)
            if select(index)
                push!(result, table[index])
            end
        end
    end
    return sort!(unique!(result))
end

# Boundary
"""
    Boundary{Names}(values::AbstractVector{<:Number}, vectors::AbstractVector{<:AbstractVector{<:Number}}) where Names

Boundary twist of operators.
"""
struct Boundary{Names, D<:Number, V<:AbstractVector} <: LinearTransformation
    values::Vector{D}
    vectors::Vector{V}
    function Boundary{Names}(values::AbstractVector{<:Number}, vectors::AbstractVector{<:AbstractVector{<:Number}}) where Names
        @assert length(Names)==length(values)==length(vectors) "Boundary error: mismatched names, values and vectors."
        datatype = promote_type(eltype(values), Float)
        new{Names, datatype, eltype(vectors)}(convert(Vector{datatype}, values), vectors)
    end
end
@inline Base.:(==)(bound₁::Boundary, bound₂::Boundary) = keys(bound₁)==keys(bound₂) && ==(efficientoperations, bound₁, bound₂)
@inline Base.isequal(bound₁::Boundary, bound₂::Boundary) = isequal(keys(bound₁), keys(bound₂)) && isequal(efficientoperations, bound₁, bound₂)
@inline Base.valtype(::Type{<:Boundary}, M::Type{<:Operator}) = reparameter(M, :value, promote_type(Complex{Int}, scalartype(M)))
@inline Base.valtype(B::Type{<:Boundary}, MS::Type{<:Operators}) = (M = valtype(B, eltype(MS)); Operators{M, idtype(M)})

"""
    keys(bound::Boundary) -> ZeroAtLeast{Symbol}
    keys(::Type{<:Boundary{Names}}) where Names -> Names

Get the names of the boundary parameters.
"""
@inline Base.keys(bound::Boundary) = keys(typeof(bound))
@inline Base.keys(::Type{<:Boundary{Names}}) where Names = Names

"""
    (bound::Boundary)(operator::Operator; origin::Union{AbstractVector, Nothing}=nothing) -> Operator

Get the boundary twisted operator.
"""
@inline function (bound::Boundary)(operator::Operator; origin::Union{AbstractVector, Nothing}=nothing)
    values = isnothing(origin) ? bound.values : bound.values-origin
    return replace(operator, operator.value*exp(1im*mapreduce(u->angle(u, bound.vectors, values), +, id(operator))))
end

"""
    update!(bound::Boundary; parameters...) -> Boundary

Update the values of the boundary twisted phase.
"""
@inline @generated function update!(bound::Boundary; parameters...)
    exprs = []
    for (i, name) in enumerate(QuoteNode.(keys(bound)))
        push!(exprs, :(bound.values[$i] = get(parameters, $name, bound.values[$i])))
    end
    return Expr(:block, exprs..., :(return bound))
end

"""
    merge!(bound::Boundary, another::Boundary) -> typeof(bound)

Merge the values and vectors of the twisted boundary condition from another one.
"""
@inline function Base.merge!(bound::Boundary, another::Boundary)
    @assert keys(bound)==keys(another) "merge! error: mismatched names of boundary parameters."
    bound.values .= another.values
    bound.vectors .= another.vectors
    return bound
end

"""
    replace(bound::Boundary; values=bound.values, vectors=bound.vectors) -> Boundary

Replace the values or vectors of a twisted boundary condition and get the new one.

!!! note
    The plain boundary condition keeps plain even when replaced with new values or new vectors.
"""
@inline Base.replace(bound::Boundary; values=bound.values, vectors=bound.vectors) = Boundary{keys(bound)}(values, vectors)

"""
    plain

Plain boundary condition without any twist.
"""
const plain = Boundary{()}(Float[], SVector{0, Float}[])
@inline Base.valtype(::Type{typeof(plain)}, M::Type{<:Operator}) = M
@inline Base.valtype(::Type{typeof(plain)}, M::Type{<:Operators}) = M
@inline (::typeof(plain))(operator::Operator; kwargs...) = operator
@inline Base.replace(::(typeof(plain)); kwargs...) = plain

end #module
