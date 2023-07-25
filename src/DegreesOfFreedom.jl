module DegreesOfFreedom

using Printf: @printf, @sprintf
using SparseArrays: SparseMatrixCSC, nnz
using StaticArrays: SVector
using ..QuantumLattices: add!, decompose, dimension, dtype
using ..QuantumOperators: ID, LinearTransformation, Operator, OperatorPack, Operators, OperatorSum, OperatorUnit, valuetolatextext, valuetostr
using ..Spatials: Bond, Point
using ..Toolkit: atol, efficientoperations, rtol, CompositeDict, CompositeTuple, Float, NamedContainer, VectorSpace, VectorSpaceCartesian, VectorSpaceDirectProducted, VectorSpaceDirectSummed, VectorSpaceStyle, commontype, concatenate, decimaltostr, fulltype, parametertype, rawtype, reparameter

import LaTeXStrings: latexstring
import ..QuantumLattices: ⊕, ⊗, expand, expand!, kind, id, ishermitian, rank, reset!, update!, value
import ..QuantumOperators: idtype, optype, script
import ..Spatials: icoordinate, rcoordinate
import ..Toolkit: contentnames, getcontent, isparameterbound, parameternames, shape

export CompositeIID, CompositeInternal, IID, Internal, SimpleIID, SimpleInternal
export AbstractCompositeIndex, CompositeIndex, Hilbert, Index, iidtype, indextype, statistics
export wildcard, Component, Constraint, Coupling, Diagonal, IIDSpace, MatrixCoupling, MatrixCouplingProd, MatrixCouplingSum, Pattern, diagonalizablefields, isdefinite, @indexes
export Term, TermAmplitude, TermCoupling, TermFunction, TermModulate, ismodulatable, sitestructure
export Metric, OperatorUnitToTuple, Table
export Boundary, plain

# IID and Internal
"""
    IID <: OperatorUnit

The id of an internal degree of freedom.
"""
abstract type IID <: OperatorUnit end

"""
    SimpleIID <: IID

The id of a simple internal degree of freedom.
"""
abstract type SimpleIID <: IID end
@inline statistics(iid::SimpleIID) = statistics(typeof(iid))
@inline isdefinite(iid::SimpleIID) = isdefinite(typeof(iid))
@inline isdefinite(::Type{<:SimpleIID}) = false

"""
    CompositeIID{T<:Tuple{Vararg{SimpleIID}}} <: IID

The composition of several single internal ids.
"""
struct CompositeIID{T<:Tuple{Vararg{SimpleIID}}} <: IID
    contents::T
end
Base.show(io::IO, ciid::CompositeIID) = @printf io "%s" join((string(ciid[i]) for i = 1:rank(ciid)), " ⊗ ")
@inline Base.length(ciid::CompositeIID) = length(typeof(ciid))
@inline Base.length(::Type{<:CompositeIID{T}}) where {T<:Tuple{Vararg{SimpleIID}}} = fieldcount(T)
@inline Base.getindex(ciid::CompositeIID, i::Int) = ciid.contents[i]
@inline Base.getproperty(ciid::CompositeIID, name::Symbol) = ciidgetproperty(ciid, Val(name))
@inline ciidgetproperty(ciid::CompositeIID, ::Val{:contents}) = getfield(ciid, :contents)
@inline ciidgetproperty(ciid::CompositeIID, ::Val{name}) where name = getproperty(getfield(ciid, :contents), name)

"""
    CompositeIID(contents::SimpleIID...)

Construct a composite iid from a set of simple iids.
"""
@inline CompositeIID(contents::SimpleIID...) = CompositeIID(contents)

"""
    rank(ciid::CompositeIID) -> Int
    rank(::Type{<:CompositeIID{T}}) where {T<:Tuple{Vararg{SimpleIID}}} -> Int

Get the number of simple iids in a composite iid.
"""
@inline rank(ciid::CompositeIID) = rank(typeof(ciid))
@inline rank(::Type{<:CompositeIID{T}}) where {T<:Tuple{Vararg{SimpleIID}}} = fieldcount(T)

"""
    iidtype(ciid::CompositeIID, i::Integer)
    iidtype(::Type{<:CompositeIID{T}}, i::Integer) where {T<:Tuple{Vararg{SimpleIID}}}

Get the type of the ith simple iid in a composite iid.
"""
@inline iidtype(ciid::CompositeIID, i::Integer) = iidtype(typeof(ciid), i)
@inline iidtype(::Type{<:CompositeIID{T}}, i::Integer) where {T<:Tuple{Vararg{SimpleIID}}} = fieldtype(T, i)

"""
    ⊗(iid₁::SimpleIID, iid₂::SimpleIID) -> CompositeIID
    ⊗(iid::SimpleIID, ciid::CompositeIID) -> CompositeIID
    ⊗(ciid::CompositeIID, iid::SimpleIID) -> CompositeIID
    ⊗(ciid₁::CompositeIID, ciid₂::CompositeIID) -> CompositeIID

Direct product between simple iids and composite iids.
"""
@inline ⊗(iid₁::SimpleIID, iid₂::SimpleIID) = CompositeIID(iid₁, iid₂)
@inline ⊗(iid::SimpleIID, ciid::CompositeIID) = CompositeIID(iid, ciid.contents...)
@inline ⊗(ciid::CompositeIID, iid::SimpleIID) = CompositeIID(ciid.contents..., iid)
@inline ⊗(ciid₁::CompositeIID, ciid₂::CompositeIID) = CompositeIID(ciid₁.contents..., ciid₂.contents...)

"""
    Internal{I<:IID} <: VectorSpace{I}

The whole internal degrees of freedom at a single point.
"""
abstract type Internal{I<:IID} <: VectorSpace{I} end

"""
    SimpleInternal{I<:SimpleIID} <: Internal{I}

The simple internal degrees of freedom at a single point.
"""
abstract type SimpleInternal{I<:SimpleIID} <: Internal{I} end
@inline VectorSpaceStyle(::Type{<:SimpleInternal}) = VectorSpaceCartesian()
@inline Base.show(io::IO, i::SimpleInternal) = @printf io "%s(%s)" i|>typeof|>nameof join(("$name=$(getfield(i, name))" for name in i|>typeof|>fieldnames), ", ")

"""
    statistics(i::SimpleInternal) -> Symbol
    statistics(::Type{<:SimpleInternal{I}}) where {I<:SimpleIID} -> Symbol

Get the statistics of a simple internal space.
"""
@inline statistics(i::SimpleInternal) = statistics(typeof(i))
@inline statistics(::Type{<:SimpleInternal{I}}) where {I<:SimpleIID} = statistics(I)

"""
    match(iid::SimpleIID, i::SimpleInternal) -> Bool
    match(::Type{I}, ::Type{SI}) where {I<:SimpleIID, SI<:SimpleInternal}

Judge whether a simple iid or a simple iid type matches a simple internal space or a simple internal space type.

Here, "match" means that the eltype of the simple internal space has the same type name with the simple iid.
"""
@inline Base.match(iid::SimpleIID, i::SimpleInternal) = match(typeof(iid), typeof(i))
@inline Base.match(::Type{I}, ::Type{SI}) where {I<:SimpleIID, SI<:SimpleInternal} = nameof(I)==nameof(eltype(SI))

"""
    filter(iid::SimpleIID, i::SimpleInternal) -> Union{Nothing, typeof(i)}
    filter(::Type{I}, i::SimpleInternal) where {I<:SimpleIID} -> Union{Nothing, typeof(i)}

Filter a simple internal space with respect to the input `iid` or type `I`.
"""
@inline Base.filter(iid::SimpleIID, i::SimpleInternal) = filter(typeof(iid), i)
@inline Base.filter(::Type{I}, i::SimpleInternal) where {I<:SimpleIID} = match(I, typeof(i)) ? i : nothing

"""
    filter(iid::SimpleIID, ::Type{T}) where {T<:SimpleInternal}
    filter(::Type{I}, ::Type{T}) where {I<:SimpleIID, T<:SimpleInternal}

Filter the type of a simple internal space with respect to the input `iid` or type `I`.
"""
@inline Base.filter(iid::SimpleIID, ::Type{T}) where {T<:SimpleInternal} = filter(typeof(iid), T)
@inline Base.filter(::Type{I}, ::Type{T}) where {I<:SimpleIID, T<:SimpleInternal} = match(I, T) ? T : nothing

"""
    CompositeInternal{K, T<:Tuple{Vararg{SimpleInternal}}} <: Internal{IID}

The composition of several single internal spaces.
"""
struct CompositeInternal{K, T<:Tuple{Vararg{SimpleInternal}}} <: Internal{IID}
    contents::T
    function CompositeInternal{K}(contents::Tuple{Vararg{SimpleInternal}}) where K
        @assert K∈(:⊕, :⊗) "CompositeInternal error: kind must be either :⊕ for direct sum or :⊗ for direct product."
        new{K, typeof(contents)}(contents)
    end
end
@inline kind(ci::CompositeInternal) = kind(typeof(ci))
@inline kind(::Type{<:CompositeInternal{K}}) where K = K
@inline Base.show(io::IO, ci::CompositeInternal) = @printf io "%s" join((string(ci.contents[i]) for i = 1:rank(ci)), @sprintf(" %s ", kind(ci)))
@inline VectorSpaceStyle(::Type{<:CompositeInternal{:⊕}}) = VectorSpaceDirectSummed()
@inline VectorSpaceStyle(::Type{<:CompositeInternal{:⊗}}) = VectorSpaceDirectProducted()
@inline @generated Base.eltype(::Type{<:CompositeInternal{:⊕, T}}) where {T<:Tuple{Vararg{SimpleInternal}}} = mapreduce(eltype, typejoin, fieldtypes(T), init=Union{})
@inline @generated Base.eltype(::Type{<:CompositeInternal{:⊗, T}}) where {T<:Tuple{Vararg{SimpleInternal}}} = CompositeIID{Tuple{map(eltype, fieldtypes(T))...}}
@inline SimpleIID(iid::SimpleIID, ::CompositeInternal) = iid
@inline CompositeIID(iid::Tuple, ::CompositeInternal) = CompositeIID(iid)

"""
    CompositeInternal{K}(contents::SimpleInternal...) where K

Construct a composite internal space from a set of simple internal spaces.
"""
@inline CompositeInternal{K}(contents::SimpleInternal...) where K = CompositeInternal{K}(contents)

"""
    rank(ci::CompositeInternal) -> Int
    rank(::Type{<:CompositeInternal{K, T}}) where {K, T<:Tuple{Vararg{SimpleInternal}}} -> Int

Get the number of simple internal spaces in a composite internal space.
"""
@inline rank(ci::CompositeInternal) = rank(typeof(ci))
@inline rank(::Type{<:CompositeInternal{K, T}}) where {K, T<:Tuple{Vararg{SimpleInternal}}} = fieldcount(T)

"""
    ⊗(i₁::SimpleInternal, i₂::SimpleInternal) -> CompositeInternal{:⊗}
    ⊗(i::SimpleInternal, ci::CompositeInternal{:⊗}) -> CompositeInternal{:⊗}
    ⊗(ci::CompositeInternal{:⊗}, i::SimpleInternal) -> CompositeInternal{:⊗}
    ⊗(ci₁::CompositeInternal{:⊗}, ci₂::CompositeInternal{:⊗}) -> CompositeInternal{:⊗}

Direct product between simple internal spaces and composite internal spaces.
"""
@inline ⊗(i₁::SimpleInternal, i₂::SimpleInternal) = CompositeInternal{:⊗}(i₁, i₂)
@inline ⊗(i::SimpleInternal, ci::CompositeInternal{:⊗}) = CompositeInternal{:⊗}(i, ci.contents...)
@inline ⊗(ci::CompositeInternal{:⊗}, i::SimpleInternal) = CompositeInternal{:⊗}(ci.contents..., i)
@inline ⊗(ci₁::CompositeInternal{:⊗}, ci₂::CompositeInternal{:⊗}) = CompositeInternal{:⊗}(ci₁.contents..., ci₂.contents...)

"""
    ⊕(i₁::SimpleInternal, i₂::SimpleInternal) -> CompositeInternal{:⊕}
    ⊕(i::SimpleInternal, ci::CompositeInternal{:⊕}) -> CompositeInternal{:⊕}
    ⊕(ci::CompositeInternal{:⊕}, i::SimpleInternal) -> CompositeInternal{:⊕}
    ⊕(ci₁::CompositeInternal{:⊕}, ci₂::CompositeInternal{:⊕}) -> CompositeInternal{:⊕}

Direct product between simple internal spaces and composite internal spaces.
"""
@inline ⊕(i₁::SimpleInternal, i₂::SimpleInternal) = CompositeInternal{:⊕}(i₁, i₂)
@inline ⊕(i::SimpleInternal, ci::CompositeInternal{:⊕}) = CompositeInternal{:⊕}(i, ci.contents...)
@inline ⊕(ci::CompositeInternal{:⊕}, i::SimpleInternal) = CompositeInternal{:⊕}(ci.contents..., i)
@inline ⊕(ci₁::CompositeInternal{:⊕}, ci₂::CompositeInternal{:⊕}) = CompositeInternal{:⊕}(ci₁.contents..., ci₂.contents...)

"""
    filter(iid::SimpleIID, ci::CompositeInternal) -> Union{Nothing, SimpleInternal, CompositeInternal}
    filter(::Type{I}, ci::CompositeInternal) where {I<:SimpleIID} -> Union{Nothing, SimpleInternal, CompositeInternal}

Filter the composite internal space and select those that matches `I` or the type of `iid`.
"""
@inline Base.filter(iid::SimpleIID, ci::CompositeInternal) = filter(typeof(iid), ci)
@inline Base.filter(::Type{I}, ci::CompositeInternal) where {I<:SimpleIID} = filterhelper₁(I, ci, filtermatches(I, typeof(ci))|>Val)
@inline @generated function filtermatches(::Type{I}, ::Type{CompositeInternal{K, C}}) where {I<:SimpleIID, K, C<:Tuple{Vararg{SimpleInternal}}}
    exprs = []
    for i = 1:fieldcount(C)
        T = fieldtype(C, i)
        push!(exprs, :(match(I, $T)))
    end
    return Expr(:tuple, exprs...)
end
@inline @generated function filterhelper₁(::Type{I}, ci::CompositeInternal, ::Val{BS}) where {I<:SimpleIID, BS}
    exprs = []
    for (i, B) in enumerate(BS)
        B && push!(exprs, :(ci.contents[$i]))
    end
    length(exprs)==0 && return
    length(exprs)==1 && return first(exprs)
    return Expr(:call, Expr(:curly, :CompositeInternal, QuoteNode(kind(ci))), exprs...)
end

"""
    filter(iid::SimpleIID, ::Type{C}) where {C<:CompositeInternal}
    filter(::Type{I}, ::Type{C}) where {I<:SimpleIID, C<:CompositeInternal}

Filter the type of a composite internal space and select those that matches `I` or the type of `iid`.
"""
@inline Base.filter(iid::SimpleIID, ::Type{C}) where {C<:CompositeInternal} = filter(typeof(iid), C)
@inline Base.filter(::Type{I}, ::Type{C}) where {I<:SimpleIID, C<:CompositeInternal} = filterhelper₂(I, C, filtermatches(I, C)|>Val)
@inline @generated function filterhelper₂(::Type{I}, ::Type{CompositeInternal{K, C}}, ::Val{BS}) where {I<:SimpleIID, K, C<:Tuple{Vararg{SimpleInternal}}, BS}
    exprs = []
    for (i, B) in enumerate(BS)
        B && push!(exprs, :(fieldtype(C, $i)))
    end
    length(exprs)==0 && return
    length(exprs)==1 && return first(exprs)
    return Expr(:curly, :CompositeInternal, QuoteNode(K), Expr(:curly, :Tuple, exprs...))
end

# Index and CompositeIndex
"""
    Index{S<:Union{Int, Colon}, I<:SimpleIID} <: OperatorUnit

The index of a degree of freedom, which consist of the spatial part and the internal part.
"""
struct Index{S<:Union{Int, Colon}, I<:SimpleIID} <: OperatorUnit
    site::S
    iid::I
end
@inline parameternames(::Type{<:Index}) = (:site, :iid)
@inline isparameterbound(::Type{<:Index}, ::Val{:site}, ::Type{S}) where {S<:Union{Int, Colon}} = !isconcretetype(S)
@inline isparameterbound(::Type{<:Index}, ::Val{:iid}, ::Type{I}) where {I<:SimpleIID} = !isconcretetype(I)
@inline isdefinite(index::Index) = isdefinite(typeof(index))
@inline isdefinite(::Type{<:Index{<:Union{Int, Colon}, I}}) where {I<:SimpleIID} = isdefinite(I)
@inline isdefinite(indexes::Tuple{Vararg{Index}}) = isdefinite(typeof(indexes))
@inline @generated isdefinite(::Type{T}) where {T<:Tuple{Vararg{Index}}} = Expr(:call, :all, Expr(:tuple, [:(isdefinite(fieldtype(T, $i))) for i=1:fieldcount(T)]...))
@inline Base.show(io::IO, index::Index{Colon}) = @printf io "Index(:, %s)" index.iid

"""
    iidtype(index::Index)
    iidtype(::Type{I}) where {I<:Index}

Get the type of the internal part of an index.
"""
@inline iidtype(index::Index) = iidtype(typeof(index))
@inline iidtype(::Type{I}) where {I<:Index} = parametertype(I, :iid)

"""
    statistics(index::Index) -> Symbol
    statistics(::Type{<:Index{I}}) where {I<:SimpleIID} -> Symbol

Get the statistics of an index.
"""
@inline statistics(index::Index) = statistics(typeof(index))
@inline statistics(::Type{<:Index{<:Union{Int, Colon}, I}}) where {I<:SimpleIID} = statistics(I)

"""
    adjoint(index::Index) -> typeof(index)

Get the adjoint of an index.
"""
@inline Base.adjoint(index::Index) = rawtype(typeof(index))(index.site, adjoint(index.iid))

"""
    AbstractCompositeIndex{I<:Index} <: OperatorUnit

The abstract type of a composite index.
"""
abstract type AbstractCompositeIndex{I<:Index} <: OperatorUnit end
@inline contentnames(::Type{<:AbstractCompositeIndex}) = (:index,)
@inline parameternames(::Type{<:AbstractCompositeIndex}) = (:index,)
@inline isparameterbound(::Type{<:AbstractCompositeIndex}, ::Val{:index}, ::Type{I}) where {I<:Index} = !isconcretetype(I)

"""
    indextype(::AbstractCompositeIndex)
    indextype(::Type{<:AbstractCompositeIndex})

Get the index type of a composite index.
"""
@inline indextype(index::AbstractCompositeIndex) = indextype(typeof(index))
@inline @generated indextype(::Type{I}) where {I<:AbstractCompositeIndex} = parametertype(supertype(I, :AbstractCompositeIndex), :index)

"""
    statistics(index::AbstractCompositeIndex) -> Symbol
    statistics(::Type{<:AbstractCompositeIndex{I}}) where {I<:Index} -> Symbol

Get the statistics of a composite operator id.
"""
@inline statistics(index::AbstractCompositeIndex) = statistics(typeof(index))
@inline statistics(::Type{<:AbstractCompositeIndex{I}}) where {I<:Index} = statistics(I)

"""
    CompositeIndex{I<:Index, V<:SVector} <: AbstractCompositeIndex{I}

Composite index of a quantum operator.
"""
struct CompositeIndex{I<:Index, V<:SVector} <: AbstractCompositeIndex{I}
    index::I
    rcoordinate::V
    icoordinate::V
    CompositeIndex(index::Index, rcoordinate::V, icoordinate::V) where {V<:SVector} = new{typeof(index), V}(index, compositeindexcoordinate(rcoordinate), compositeindexcoordinate(icoordinate))
end
@inline contentnames(::Type{<:CompositeIndex}) = (:index, :rcoordinate, :icoordinate)
@inline parameternames(::Type{<:CompositeIndex}) = (:index, :coordination)
@inline isparameterbound(::Type{<:CompositeIndex}, ::Val{:coordination}, ::Type{V}) where {V<:SVector} = !isconcretetype(V)
@inline Base.hash(index::CompositeIndex, h::UInt) = hash((index.index, Tuple(index.rcoordinate)), h)
@inline Base.propertynames(::ID{CompositeIndex}) = (:indexes, :rcoordinates, :icoordinates)
@inline Base.show(io::IO, index::CompositeIndex) = @printf io "CompositeIndex(%s, %s, %s)" index.index index.rcoordinate index.icoordinate
@inline compositeindexcoordinate(vector::SVector) = vector
@inline compositeindexcoordinate(vector::SVector{N, Float}) where N = SVector(ntuple(i->vector[i]===-0.0 ? 0.0 : vector[i], Val(N)))

"""
    CompositeIndex(index::Index, rcoordinate, icoordinate)
    CompositeIndex(index::Index; rcoordinate, icoordinate)

Construct an operator id.
"""
@inline CompositeIndex(index::Index, rcoordinate, icoordinate) = CompositeIndex(index, SVector{length(rcoordinate)}(rcoordinate), SVector{length(icoordinate)}(icoordinate))
@inline CompositeIndex(index::Index; rcoordinate, icoordinate) = CompositeIndex(index, rcoordinate, icoordinate)

"""
    adjoint(index::CompositeIndex) -> typeof(index)

Get the adjoint of an operator id.
"""
@inline Base.adjoint(index::CompositeIndex) = CompositeIndex(index.index', index.rcoordinate, index.icoordinate)

"""
    indextype(I::Type{<:SimpleInternal}, P::Type{<:Point}, ::Val)

Get the compatible composite index type based on the information of its internal part.
"""
@inline function indextype(I::Type{<:SimpleInternal}, P::Type{<:Point}, ::Val)
    return fulltype(CompositeIndex, NamedTuple{(:index, :coordination), Tuple{fulltype(Index, NamedTuple{(:site, :iid), Tuple{Int, eltype(I)}}), SVector{dimension(P), dtype(P)}}})
end

"""
    rcoordinate(opt::Operator{<:Number, <:ID{CompositeIndex}}) -> SVector

Get the whole rcoordinate of an operator.
"""
@inline function rcoordinate(opt::Operator{<:Number, <:ID{CompositeIndex}})
    rank(opt)==1 && return id(opt)[1].rcoordinate
    rank(opt)==2 && return id(opt)[2].rcoordinate-id(opt)[1].rcoordinate
    error("rcoordinate error: not supported rank($(rank(opt))) of $(nameof(opt)).")
end

"""
    icoordinate(opt::Operator{<:Number, <:ID{CompositeIndex}}) -> SVector

Get the whole icoordinate of an operator.
"""
@inline function icoordinate(opt::Operator{<:Number, <:ID{CompositeIndex}})
    rank(opt)==1 && return id(opt)[1].icoordinate
    rank(opt)==2 && return id(opt)[2].icoordinate-id(opt)[1].icoordinate
    error("icoordinate error: not supported rank($(rank(opt))) of $(nameof(opt)).")
end

"""
    script(::Val{:rcoordinate}, index::CompositeIndex; kwargs...) -> String
    script(::Val{:icoordinate}, index::CompositeIndex; kwargs...) -> String

Get the `rcoordinate/icoordinate` script of a composite index.
"""
@inline script(::Val{:rcoordinate}, index::CompositeIndex; kwargs...) = @sprintf "[%s]" join(valuetolatextext.(index.rcoordinate), ", ")
@inline script(::Val{:icoordinate}, index::CompositeIndex; kwargs...) = @sprintf "[%s]" join(valuetolatextext.(index.icoordinate), ", ")

"""
    script(::Val{:integercoordinate}, index::CompositeIndex; vectors, kwargs...)

Get the integral script of the icoordinate of an composite index.
"""
function script(::Val{:integercoordinate}, index::CompositeIndex; vectors, kwargs...)
    rcoeff = decompose(index.icoordinate, vectors...)
    icoeff = Int.(round.(rcoeff))
    @assert isapprox(efficientoperations, rcoeff, icoeff) "script error: mismatched icoordinate of the input composite index and vectors."
    return @sprintf "[%s]" join(icoeff, ", ")
end

"""
    script(::Val{attr}, index::CompositeIndex; kwargs...) where attr

Get the `attr` script of an index, which is contained in its index.
"""
@inline script(::Val{attr}, index::CompositeIndex; kwargs...) where attr = script(Val(attr), index.index; kwargs...)

"""
    script(::Val{:site}, index::Index; kwargs...) -> String
    script(attr::Val, index::Index; kwargs...) -> String

Get the required script of a spin index.
"""
@inline script(::Val{:site}, index::Index; kwargs...) = string(index.site)
@inline script(::Val{:site}, index::Index{Colon}; kwargs...) = ":"
@inline script(attr::Val, index::Index; kwargs...) = script(attr, index.iid; kwargs...)

# Hilbert
"""
    Hilbert{I<:Internal} <: CompositeDict{Int, I}

Hilbert space at a lattice.
"""
struct Hilbert{I<:Internal} <: CompositeDict{Int, I}
    contents::Dict{Int, I}
    Hilbert(contents::Dict{Int, <:Internal}) = new{valtype(contents)}(contents)
end

"""
    Hilbert(ps::Pair...)
    Hilbert(kv)

Construct a Hilbert space the same way as a Dict.
"""
@inline Hilbert(ps::Pair...) = Hilbert(ps)
@inline Hilbert(kv) = Hilbert(Dict(kv))

"""
    Hilbert(internals::Internal...)
    Hilbert(internals::Tuple{Vararg{Internal}})
    Hilbert(internals::AbstractVector{<:Internal})

Construct a Hilbert space with the given internals.
"""
@inline Hilbert(internals::Internal...) = Hilbert(internals)
@inline Hilbert(internals::Tuple{Vararg{Internal}}) = Hilbert(i=>internal for (i, internal) in enumerate(internals))
@inline Hilbert(internals::AbstractVector{<:Internal}) = Hilbert(i=>internal for (i, internal) in enumerate(internals))

"""
    Hilbert(internal::Internal, num::Int)

Construct a Hilbert space with all internal spaces the same.
"""
@inline Hilbert(internal::Internal, num::Int) = Hilbert(i=>internal for i in 1:num)

# Coupling
## IIDSpace
"""
    IIDSpace{I<:IID, V<:Internal} <: VectorSpace{IID}

The space expanded by a "labeled" iid.
"""
struct IIDSpace{I<:IID, V<:Internal} <: VectorSpace{IID}
    iid::I
    internal::V
end
@inline Base.eltype(::Type{<:IIDSpace{<:IID, V}}) where {V<:Internal} = eltype(V)

### IIDSpace of SimpleIID and SimpleInternal
@inline VectorSpaceStyle(::Type{<:IIDSpace{<:SimpleIID, <:SimpleInternal}}) = VectorSpaceCartesian()
@inline Base.CartesianIndex(iid::SimpleIID, iidspace::IIDSpace{<:SimpleIID, <:SimpleInternal}) = CartesianIndex(iid, iidspace.internal)
@inline function Base.getindex(::VectorSpaceCartesian, iidspace::IIDSpace{<:SimpleIID, <:SimpleInternal}, index::CartesianIndex)
    return rawtype(eltype(iidspace))(index, iidspace.internal)
end

### IIDSpace of CompositeIID and CompositeInternal
@inline VectorSpaceStyle(::Type{<:IIDSpace{<:CompositeIID, <:CompositeInternal}}) = VectorSpaceDirectProducted()
@inline function getcontent(iidspace::IIDSpace{<:CompositeIID, <:CompositeInternal}, ::Val{:contents})
    return map((iid, internal)->IIDSpace(iid, internal), iidspace.iid.contents, iidspace.internal.contents)
end
@inline CompositeIID(iids::Tuple{Vararg{SimpleIID}}, ::IIDSpace{<:CompositeIID, <:CompositeInternal}) = CompositeIID(iids)

"""
    expand(iid::SimpleIID, internal::SimpleInternal) -> IIDSpace
    expand(iids::NTuple{N, SimpleIID}, internals::NTuple{N, SimpleInternal}) where N -> IIDSpace

Get the space expanded by a set of "labeled" iids.
"""
@inline expand(iid::SimpleIID, internal::SimpleInternal) = expand((iid,), (internal,))
@inline expand(iids::NTuple{N, SimpleIID}, internals::NTuple{N, SimpleInternal}) where N = IIDSpace(CompositeIID(iids), CompositeInternal{:⊗}(internals))

## Constraint
const wildcard = Symbol("*")
"""
    Pattern <: Function

Construct a pattern for a set of `Index`es.
"""
struct Pattern <: Function
    pattern::Dict{Symbol, Vector{Tuple{Int, Symbol}}}
end
function Pattern(indexes::Tuple{Index, Vararg{Index}})
    pattern = Dict{Symbol, Vector{Tuple{Int, Symbol}}}()
    for (i, index) in enumerate(indexes)
        iid = index.iid
        for attr in fieldnames(typeof(iid))
            value = getfield(iid, attr)
            if isa(value, Symbol) && value≠wildcard
                haskey(pattern, value) || (pattern[value] = Tuple{Int, Symbol}[])
                push!(pattern[value], (i, attr))
            end
        end
    end
    return Pattern(pattern)
end
function (pattern::Pattern)(indexes::Tuple{Index, Vararg{Index}})
    length(pattern.pattern)==0 && return true
    for component in values(pattern.pattern)
        length(component)>1 && (allequal(map(pair::Tuple{Int, Symbol}->getfield(indexes[pair[1]].iid, pair[2]), component)) || return false)
    end
    return true
end

"""
    Diagonal{Fields} <: Function

Construct a pattern for a set of homogenous `Index`es that all the specified fields of their contained iids should be diagonal, respectively.
"""
struct Diagonal{Fields} <: Function
    Diagonal(fields::Tuple{Vararg{Symbol}}) = new{fields}()
end
@inline Diagonal(fields::Symbol...) = Diagonal(fields)
@inline Diagonal(::Type{I}) where {I<:Index} = _diagonal(I, I|>diagonalizablefields|>Val)
@generated function _diagonal(::Type{I}, ::Val{choices}) where {I<:Index, choices}
    fields = Symbol[]
    T = iidtype(I)
    for choice in choices
        fieldtype(T, choice)==Colon && push!(fields, choice)
    end
    return Diagonal(fields...)
end
@inline @generated diagonalizablefields(::Type{I}) where {I<:Index} = fieldnames(iidtype(I))
@generated function (diagonal::Diagonal{fields})(indexes::Tuple{I, Vararg{I, N}}) where {fields, I<:Index, N}
    exprs = []
    for field in QuoteNode.(fields)
        push!(exprs, Expr(:call, allequal, Expr(:tuple, [:(getfield(indexes[$i].iid, $field)) for i=1:(N+1)]...)))
    end
    return Expr(:call, all, Expr(:tuple, exprs...))
end

"""
    Constraint{RS, N, C<:NTuple{N, Function}}

The constraint of the indexes of internal degrees of freedom in a coupling.
"""
struct Constraint{RS, N, C<:NTuple{N, Function}}
    representations::NTuple{N, String}
    conditions::C
    function Constraint{RS}(representations::NTuple{N, String}, conditions::NTuple{N, Function})  where {RS, N}
        @assert isa(RS, NTuple) "Constraint error: ranks (`RS`) must be tuple of integers."
        @assert length(RS)==N "Constraint error: mismatched number of ranks ($RS), representations ($representations) and conditions ($conditions)."
        new{RS, N, typeof(conditions)}(representations, conditions)
    end
end
@inline Base.:(==)(constraint₁::Constraint{RS₁}, constraint₂::Constraint{RS₂}) where {RS₁, RS₂} = RS₁==RS₂ && constraint₁.representations==constraint₂.representations
@inline Base.:isequal(constraint₁::Constraint{RS₁}, constraint₂::Constraint{RS₂}) where {RS₁, RS₂} = isequal(RS₁, RS₂) && isequal(constraint₁.representations, constraint₂.representations)
@inline Base.hash(constraint::Constraint{RS}, h::UInt) where RS = hash((RS, constraint.representations), h)

"""
    Constraint{R}() where R
    Constraint{R}(condition::Union{Pattern, Diagonal}) where R
    Constraint{R}(representation::String, condition::Function) where R

Construct a constraint with only one condition.
"""
@inline Constraint{R}() where R = Constraint{R}(Diagonal())
@inline Constraint{R}(condition::Union{Pattern, Diagonal}) where R = Constraint{R}("pattern", condition)
@inline Constraint{R}(representation::String, condition::Function) where R = Constraint{(R,)}((representation,), (condition,))

"""
    Constraint(indexes::Index...)
    Constraint(indexes::NTuple{N, Index}) where N

Construct a constraint based on the pattern of the input indexes.
"""
@inline Constraint(index::Index, indexes::Index...) = Constraint((index, indexes...))
@inline Constraint(indexes::Tuple{Index, Vararg{Index}}) = Constraint{fieldcount(typeof(indexes))}("pattern", constraint(indexes, indexes|>typeof|>isdiagonalable|>Val))
@inline constraint(indexes::Tuple{Index, Vararg{Index}}, ::Val{true}) = Diagonal(eltype(indexes))
@inline constraint(indexes::Tuple{Index, Vararg{Index}}, ::Val{false}) = Pattern(indexes)
@generated function isdiagonalable(::Type{T}) where {T<:Tuple{Index, Vararg{Index}}}
    types = fieldtypes(T)
    allequal(types) || return false
    for type in map(iidtype, types)
        for i = 1:fieldcount(type)
            fieldtype(type, i)==Symbol && return false
        end
    end
    return true
end

"""
    rank(constraint::Constraint) -> Int
    rank(::Type{<:Constraint{RS}}) where RS -> Int

Get the rank of the coupling indexes that a constraint can apply.
"""
@inline rank(constraint::Constraint) = rank(typeof(constraint))
@inline @generated rank(::Type{<:Constraint{RS}}) where RS = sum(RS)

"""
    rank(constraint::Constraint, i::Integer) -> Int
    rank(::Type{<:Constraint{RS}}, i::Integer) where RS -> Int

Get the rank of the ith homogenous segment of the coupling indexes that a constraint can apply.
"""
@inline rank(constraint::Constraint, i::Integer) = rank(typeof(constraint), i)
@inline rank(::Type{<:Constraint{RS}}, i::Integer) where RS = RS[i]

"""
    match(constraint::Constraint, indexes::Tuple{Vararg{Index}}) -> Bool

Judge whether a composite iid fulfills a constraint.
"""
@generated function Base.match(constraint::Constraint{RS}, indexes::Tuple{Index, Vararg{Index}}) where RS
    @assert rank(constraint)==fieldcount(indexes) "match error: mismatched rank of indexes and constraint."
    exprs, count = [], 1
    for (i, r) in enumerate(RS)
        start, stop = count, count+r-1
        segment = Expr(:tuple, [:(indexes[$pos]) for pos=start:stop]...)
        push!(exprs, :(constraint.conditions[$i]($segment)::Bool))
        count = stop+1
    end
    return Expr(:call, :all, Expr(:tuple, exprs...))
end

"""
    *(constraint₁::Constraint, constraint₂::Constraint) -> Constraint

Get the combination of two sets of constraints.
"""
@inline function Base.:*(constraint₁::Constraint{RS₁}, constraint₂::Constraint{RS₂}) where {RS₁, RS₂}
    return Constraint{(RS₁..., RS₂...)}((constraint₁.representations..., constraint₂.representations...), (constraint₁.conditions..., constraint₂.conditions...))
end

"""
    @indexes index₁ index₂ ...
    @indexes(index₁, index₂, ...; constraint=...)

Construct an set of indexes and its constraint according to the input index pattern and an optional constraint.
"""
macro indexes(exprs...)
    if exprs[1].head==:parameters
        @assert(
            length(exprs[1].args)==1 && exprs[1].args[1].head==:kw && exprs[1].args[1].args[1]==:constraint,
            "@indexes error: constraint must be specified by the keyword argument `constraint`."
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
        @assert pattern.head==:call && pattern.args[end].head==:call "@indexes error: wrong pattern."
        attrs = []
        for (j, attr) in enumerate(pattern.args[end].args[2:end])
            isa(attr, Expr) && !(attr.head==:call && length(attr.args)==3 && attr.args[1]==:// && isa(attr.args[2], Int) && isa(attr.args[3], Int)) && (attr = Symbol(attr))
            if isa(attr, Symbol) && attr≠wildcard
                @assert Base.isidentifier(attr) "@indexes error: wrong pattern."
                push!(blocks, quote
                    local $attr
                    if @isdefined($attr)
                        ($attr)!=getfield(indexes[$i].iid, $j) && return false
                    else
                        $attr = getfield(indexes[$i].iid, $j)
                    end
                end)
                attr = QuoteNode(attr)
            end
            push!(attrs, attr)
        end
        push!(indexes, Expr(:call, pattern.args[1], pattern.args[2], Expr(:call, pattern.args[end].args[1], attrs...)))
    end
    N = length(indexes)
    indexes = Expr(:tuple, indexes...)
    blocks = Expr(:block, vcat([block.args for block in blocks]...)...)
    name = gensym("constraint")
    representation = constraint==true ? "pattern" : string(constraint)
    return quote
        function ($name)(indexes::NTuple{$N, Index})
            $blocks
            return $constraint
        end
        ($(esc(indexes)), Constraint{$N}($representation, $name))
    end
end

## Coupling
"""
    Coupling{V, I<:ID{Index}, C<:Constraint} <: OperatorPack{V, Tuple{I, C}}

The coupling intra/inter internal degrees of freedom at different lattice points.
"""
struct Coupling{V, I<:ID{Index}, C<:Constraint} <: OperatorPack{V, Tuple{I, C}}
    value::V
    indexes::I
    constraint::C
end
@inline parameternames(::Type{<:Coupling}) = (:value, :indexes, :constraint)
@inline isparameterbound(::Type{<:Coupling}, ::Val{:indexes}, ::Type{I}) where {I<:ID{Index}} = !isconcretetype(I)
@inline isparameterbound(::Type{<:Coupling}, ::Val{:constraint}, ::Type{C}) where {C<:Constraint} = !isconcretetype(C)
@inline idtype(M::Type{<:Coupling}) = Tuple{parametertype(M, :indexes), parametertype(M, :constraint)}
@inline getcontent(coupling::Coupling, ::Val{:id}) = (coupling.indexes, coupling.constraint)
@inline rank(M::Type{<:Coupling}) = fieldcount(parametertype(M, :indexes))
@inline Coupling(id::Tuple{ID{Index}, Constraint}) = Coupling(1, id)
@inline Coupling(value, id::Tuple{ID{Index}, Constraint}) = Coupling(value, id...)
@inline CompositeIID(coupling::Coupling, ::Val=Val(:term)) = CompositeIID(coupling.indexes.iids)
@inline Constraint(coupling::Coupling, ::Val=Val(:term)) = coupling.constraint
@inline Base.iterate(coupling::Coupling) = (coupling, nothing)
@inline Base.iterate(coupling::Coupling, ::Nothing) = nothing
@inline Base.eltype(coupling::Coupling) = eltype(typeof(coupling))
@inline Base.eltype(C::Type{<:Coupling}) = C
@inline Base.length(coupling::Coupling) = length(typeof(coupling))
@inline Base.length(::Type{<:Coupling}) = 1
@inline function Base.show(io::IO, coupling::Coupling)
    @printf io "%s" (coupling.value≈1 ? "" : coupling.value≈-1 ? "- " : string(decimaltostr(coupling.value), " "))
    len, count = length(coupling.constraint.representations), 1
    for i = 1:len
        r = rank(coupling.constraint, i)
        start, stop = count, count+r-1
        indexes = coupling.indexes[start:stop]
        if i==1
            @printf io "%s" (isdefinite(indexes) ? "" : "∑")
            (len>1 || !isdefinite(indexes)) && @printf io "%s" "["
        else
            @printf io " ⋅ %s[" (isdefinite(indexes) ? "" : "∑")
        end
        @printf io "%s" join(indexes, " ")
        (len>1 || !isdefinite(indexes)) && @printf io "%s" "]"
        coupling.constraint.representations[i]=="pattern" || @printf io "(%s)" coupling.constraint.representations[i]
        count = stop+1
    end
end
@inline Base.summary(io::IO, couplings::Vector{<:Coupling}) = @printf io "%s-element Vector{Coupling}" length(couplings)

"""
    Coupling(indexes::Index...)
    Coupling(value, indexes::Index...)
    Coupling(value, indexes::Tuple{Vararg{Index}})

Construct a `Coupling` with the input `indexes` as the pattern.
"""
@inline Coupling(index::Index, indexes::Index...) = Coupling(1, (index, indexes...))
@inline Coupling(value, index::Index, indexes::Index...) = Coupling(value, (index, indexes...))
@inline Coupling(value, indexes::Tuple{Index, Vararg{Index}}) = Coupling(value, indexes, Constraint(indexes))

"""
    Coupling(sites::Union{NTuple{N, Int}, Colon}, ::Type{I}, fields::Union{NTuple{N}, Colon}...) where {N, I<:SimpleIID}
    Coupling(value, sites::Union{NTuple{N, Int}, Colon}, ::Type{I}, fields::Union{NTuple{N}, Colon}...) where {N, I<:SimpleIID}
    Coupling{N}(sites::Union{NTuple{N, Int}, Colon}, ::Type{I}, fields::Union{NTuple{N}, Colon}...) where {N, I<:SimpleIID}
    Coupling{N}(value, sites::Union{NTuple{N, Int}, Colon}, ::Type{I}, fields::Union{NTuple{N}, Colon}...) where {N, I<:SimpleIID}

Construct a `Coupling` with the input sites and the fields of a kind of simple iid.
"""


@inline Coupling(sites::Union{NTuple{N, Int}, Colon}, ::Type{I}, fields::Union{NTuple{N}, Colon}...) where {N, I<:SimpleIID} = Coupling{N}(sites, I, fields...)
@inline Coupling(value, sites::Union{NTuple{N, Int}, Colon}, ::Type{I}, fields::Union{NTuple{N}, Colon}...) where {N, I<:SimpleIID} = Coupling{N}(value, sites, I, fields...)
@inline Coupling{N}(sites::Union{NTuple{N, Int}, Colon}, ::Type{I}, fields::Union{NTuple{N}, Colon}...) where {N, I<:SimpleIID} = Coupling{N}(1, sites, I, fields...)
@inline function Coupling{N}(value, sites::Union{NTuple{N, Int}, Colon}, ::Type{I}, fields::Union{NTuple{N}, Colon}...) where {N, I<:SimpleIID}
    return Coupling(value, ID(Index, default(sites, Val(N)), ID(I, map(field->default(field, Val(N)), fields)...)))
end
@inline default(fields, ::Val) = fields
@inline default(::Colon, N::Val) = ntuple(i->:, N)

"""
    *(cp₁::Coupling, cp₂::Coupling) -> Coupling

Get the multiplication between two coupling.
"""
@inline Base.:*(cp₁::Coupling, cp₂::Coupling) = Coupling(cp₁.value*cp₂.value, (cp₁.indexes..., cp₂.indexes...), cp₁.constraint*cp₂.constraint)

"""
    latexstring(coupling::Coupling) -> String

Convert a `Coupling` to the latex format.
"""
function latexstring(coupling::Coupling)
    result = String[]
    len, count = length(coupling.constraint.representations), 1
    for i = 1:len
        r = rank(coupling.constraint, i)
        start, stop = count, count+r-1
        indexes = coupling.indexes[start:stop]
        pattern = coupling.constraint.representations[i]
        summation = pattern=="pattern" ? "" : replace(replace("$(coupling.constraint.representations[i])", "&&"=>"\\,\\text{and}\\,"), "||"=>"\\,\\text{or}\\,")
        summation=="" || (summation = join(push!(symbols(indexes, pattern), summation), ",\\,"))
        temp = isdefinite(indexes) ? "" : @sprintf "\\sum_{%s}" summation
        for index in indexes
            temp = @sprintf "%s %s" temp latexstring(index)
        end
        push!(result, temp)
        count = stop+1
    end
    return @sprintf "%s%s" valuetostr(coupling.value) join(result, " \\cdot ")
end
function symbols(indexes::Tuple{Vararg{Index}}, constraint::String)
    result = String[]
    for index in indexes
        iid = index.iid
        for i = 1:fieldcount(typeof(iid))
            attr = getfield(iid, i)
            if isa(attr, Symbol)
                value = string(attr)
                occursin(value, constraint) || push!(result, value)
            end
        end
    end
    return sort!(unique!(result))
end

## MatrixCoupling
"""
    Component{T₁, T₂} <: VectorSpace{Tuple{T₁, T₁, T₂}}

A component of a `MatrixCoupling`, i.e., a matrix acting on a separated internal space.
"""
struct Component{T₁, T₂, V<:AbstractVector{T₁}} <: VectorSpace{Tuple{T₁, T₁, T₂}}
    left::V
    right::V
    matrix::SparseMatrixCSC{T₂, Int}
    function Component(left::V, right::V, matrix::AbstractMatrix{T}) where {V<:AbstractVector, T}
        @assert length(left)==size(matrix)[1] && length(right)==size(matrix)[2] "Component error: mismatched inputs."
        new{eltype(V), T, V}(left, right, convert(SparseMatrixCSC{T, Int}, matrix))
    end
end
@inline Base.length(component::Component) = nnz(component.matrix)
@inline function Base.getindex(component::Component, i::Int)
    row = component.matrix.rowval[i]
    col = searchsortedlast(component.matrix.colptr, i)
    return (component.left[row], component.right[col], component.matrix.nzval[i])
end
@inline Base.promote_rule(::Type{Component{T, T₁, V}}, ::Type{Component{T, T₂, V}}) where {T, T₁, T₂, V<:AbstractVector{T}} = Component{T, promote_type(T₁, T₂), V}
@inline function Base.convert(::Type{Component{T, T₁, V}}, component::Component{T, T₂, V}) where {T, T₁, T₂, V<:AbstractVector{T}}
    return Component(component.left, component.right, convert(SparseMatrixCSC{T₁, Int}, component.matrix))
end

"""
    MatrixCoupling{I<:SimpleIID, C<:Tuple{Vararg{Component}}} <: VectorSpace{Coupling}

A set of `Coupling`s whose coefficients are specified by matrices acting on separated internal spaces.
"""
struct MatrixCoupling{I<:SimpleIID, S<:Union{Int, Colon}, C<:Tuple{Vararg{Component}}} <: VectorSpace{Coupling}
    sites::Tuple{S, S}
    contents::C
    function MatrixCoupling{I}(sites::Union{NTuple{2, Int}, NTuple{2, Colon}}, contents::Tuple{Vararg{Component}}) where {I<:SimpleIID}
        @assert fieldcount(I)==length(contents) "MatrixCoupling error: mismatched iidtype ($nameof(I)) and components (len=$length(contents))."
        new{I, eltype(sites), typeof(contents)}(sites, contents)
    end
end
@inline MatrixCoupling(sites::Union{NTuple{2, Int}, Colon}, ::Type{I}, contents::Component...) where {I<:SimpleIID} = MatrixCoupling{I}(default(sites, Val(2)), contents)
@inline @generated function Base.eltype(::Type{MC}) where {MC<:MatrixCoupling}
    types = fieldtypes(fieldtype(MC, :contents))
    V = Expr(:call, :promote_type, [:(parametertype($C, 2)) for C in types]...)
    R = Expr(:curly, :Index, parametertype(MC, 2), Expr(:call, :iidtype, parametertype(MC, 1), [:(parametertype($C, 1)) for C in types]...))
    C = Expr(:curly, :Constraint, (2,), 1, Expr(:curly, Tuple, :(typeof(Diagonal($R)))))
    return :(Coupling{$V, Tuple{$R, $R}, $C})
end
@inline VectorSpaceStyle(::Type{<:MatrixCoupling}) = VectorSpaceDirectProducted()
function Coupling(contents::Tuple, mc::MatrixCoupling{I}) where {I<:SimpleIID}
    value = mapreduce(content->content[3], *, contents, init=1)
    index₁ = Index(mc.sites[1], I(map(content->content[1], contents)...))
    index₂ = Index(mc.sites[2], I(map(content->content[2], contents)...))
    return Coupling(value, index₁, index₂)
end
@inline function Base.promote_rule(::Type{MatrixCoupling{I, S, C₁}}, ::Type{MatrixCoupling{I, S, C₂}}) where {I<:SimpleIID, S<:Union{Int, Colon}, C₁<:Tuple{Vararg{Component}}, C₂<:Tuple{Vararg{Component}}}
    return MatrixCoupling{I, S, _promote_type(C₁, C₂)}
end
@inline Base.convert(::Type{MatrixCoupling{I, S, C}}, mc::MatrixCoupling{I, S}) where {I<:SimpleIID, S<:Union{Int, Colon}, C<:Tuple{Vararg{Component}}} = MatrixCoupling{I}(mc.sites, convert(C, mc.contents))
@inline @generated function _promote_type(::Type{T₁}, ::Type{T₂}) where {N, T₁<:NTuple{N, Any}, T₂<:NTuple{N, Any}}
    return Expr(:curly, Tuple, [:(promote_type(fieldtype(T₁, $i), fieldtype(T₂, $i))) for i=1:N]...)
end

"""
    MatrixCouplingProd{V<:Number, C<:Tuple{Vararg{MatrixCoupling}}} <: VectorSpace{Coupling}

The product of a set of `Coupling`s whose coefficients are specified by matrices.
"""
struct MatrixCouplingProd{V<:Number, C<:Tuple{Vararg{MatrixCoupling}}} <: VectorSpace{Coupling}
    value::V
    contents::C
end
@inline MatrixCouplingProd(contents::MatrixCoupling...) = MatrixCouplingProd(1, contents)
@inline MatrixCouplingProd(value::Number, contents::MatrixCoupling...) = MatrixCouplingProd(value, contents)
@inline Base.eltype(::Type{MCP}) where {MCP<:MatrixCouplingProd} = _eltype(fieldtype(MCP, :value), _eltypes(MCP))
@inline @generated _eltypes(::Type{MCP}) where {MCP<:MatrixCouplingProd} = Expr(:curly, :Tuple, [:(eltype($C)) for C in fieldtypes(fieldtype(MCP, :contents))]...)
@inline @generated function _eltype(::Type{V}, ::Type{TS}) where {V<:Number, TS<:Tuple}
    types = fieldtypes(TS)
    MV = promote_type(V, map(valtype, types)...)
    IS = Tuple{concatenate(map(C->fieldtypes(parametertype(C, :indexes)), types)...)...}
    FS = Tuple{concatenate(map(C->fieldtypes(fieldtype(parametertype(C, :constraint), :conditions)), types)...)...}
    N = length(types)
    RS = ntuple(i->2, Val(N))
    return Coupling{MV, IS, Constraint{RS, N, FS}}
end
@inline VectorSpaceStyle(::Type{<:MatrixCouplingProd}) = VectorSpaceDirectProducted()
@inline Coupling(contents::Tuple{Vararg{Coupling}}, mcp::MatrixCouplingProd) = prod(contents; init=mcp.value)
@inline function Base.promote_rule(::Type{MatrixCouplingProd{V₁, C₁}}, ::Type{MatrixCouplingProd{V₂, C₂}}) where {V₁<:Number, C₁<:Tuple{Vararg{MatrixCoupling}}, V₂<:Number, C₂<:Tuple{Vararg{MatrixCoupling}}}
    return MatrixCouplingProd{promote_type(V₁, V₂), _promote_type(C₁, C₂)}
end
@inline Base.convert(::Type{MatrixCouplingProd{V, C}}, mcp::MatrixCouplingProd) where {V<:Number, C<:Tuple{Vararg{MatrixCoupling}}} = MatrixCouplingProd(convert(V, mcp.value), convert(C, mcp.contents))

"""
    MatrixCouplingSum{C<:MatrixCouplingProd, N} <: VectorSpace{Coupling}

The sum of a set of `Coupling`s whose coefficients are specified by matrices.
"""
struct MatrixCouplingSum{C<:MatrixCouplingProd, N} <: VectorSpace{Coupling}
    contents::NTuple{N, C}
end
@inline MatrixCouplingSum(contents::Union{MatrixCoupling, MatrixCouplingProd}...) = MatrixCouplingSum(promote(map(m->1*m, contents)...))
@inline Base.eltype(::Type{<:MatrixCouplingSum{C}}) where {C<:MatrixCouplingProd} = eltype(eltype(C))
@inline VectorSpaceStyle(::Type{<:MatrixCouplingSum}) = VectorSpaceDirectSummed()

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

The product between `MatrixCoupling`s and `MatrixCouplingProd`s.
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

The addition between `MatrixCoupling`s and `MatrixCouplingProd`s.
"""
@inline Base.:+(mc₁::Union{MatrixCoupling, MatrixCouplingProd}, mc₂::Union{MatrixCoupling, MatrixCouplingProd}) = MatrixCouplingSum(mc₁, mc₂)
@inline Base.:+(mc::Union{MatrixCoupling, MatrixCouplingProd}, mcs::MatrixCouplingSum) = MatrixCouplingSum(mc, mcs.contents...)
@inline Base.:+(mcs::MatrixCouplingSum, mc::Union{MatrixCoupling, MatrixCouplingProd}) = MatrixCouplingSum(mcs.contents..., mc)
@inline Base.:+(mcs₁::MatrixCouplingSum, mcs₂::MatrixCouplingSum) = MatrixCouplingSum(mcs₁.contents..., mcs₂.contents...)

"""
    ^(mc::Union{MatrixCoupling, MatrixCouplingProd, MatrixCouplingSum}, n::Int) -> Union{MatrixCoupling, MatrixCouplingProd, MatrixCouplingSum}

Get the nth power of a `MatrixCoupling`/`MatrixCouplingProd`/`MatrixCouplingSum`.
"""
@inline Base.:^(mc::Union{MatrixCoupling, MatrixCouplingProd, MatrixCouplingSum}, n::Int) = prod(ntuple(i->mc, Val(n)); init=1)
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

Define right-division, minus and subtraction operator.
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
    TermFunction <: Function

Abstract type for concrete term functions.
"""
abstract type TermFunction <: Function end
@inline Base.:(==)(tf₁::TermFunction, tf₂::TermFunction) = ==(efficientoperations, tf₁, tf₂)
@inline Base.isequal(tf₁::TermFunction, tf₂::TermFunction) = isequal(efficientoperations, tf₁, tf₂)

"""
    TermAmplitude(amplitude::Union{Function, Nothing}=nothing)

The function for the amplitude of a term.
"""
struct TermAmplitude{A<:Union{Function, Nothing}} <: TermFunction
    amplitude::A
    TermAmplitude(amplitude::Union{Function, Nothing}=nothing) = new{typeof(amplitude)}(amplitude)
end
@inline (termamplitude::TermAmplitude{Nothing})(::Bond) = 1
@inline (termamplitude::TermAmplitude{<:Function})(bond::Bond) = termamplitude.amplitude(bond)

"""
    TermCoupling{E<:Coupling, C} <: TermFunction

The function for the coupling of a term.
"""
struct TermCoupling{E<:Coupling, C} <: TermFunction
    coupling::C
    TermCoupling(coupling) = new{eltype(coupling), typeof(coupling)}(coupling)
    TermCoupling(coupling::Function) = new{eltype(commontype(coupling, Tuple{Bond}, Any)), typeof(coupling)}(coupling)
    TermCoupling{E}(coupling::Function) where {E<:Coupling} = new{E, typeof(coupling)}(coupling)
end
@inline Base.valtype(termcoupling::TermCoupling) = valtype(typeof(termcoupling))
@inline Base.valtype(::Type{<:TermCoupling{E}}) where {E<:Coupling} = E
@inline (termcoupling::TermCoupling)(::Bond) = termcoupling.coupling
@inline (termcoupling::TermCoupling{<:Coupling, <:Function})(bond::Bond) = termcoupling.coupling(bond)

"""
    TermModulate(id::Symbol, modulate::Function)
    TermModulate(id::Symbol, modulate::Bool)

The function for the modulation of a term.
"""
struct TermModulate{M<:Union{Function, Val{true}, Val{false}}, id} <: TermFunction
    modulate::M
    TermModulate(id::Symbol, modulate::Function) = new{typeof(modulate), id}(modulate)
    TermModulate(id::Symbol, modulate::Bool=true) = new{Val{modulate}, id}(modulate|>Val)
end
@inline (termmodulate::TermModulate{Val{true}, id})(args...; kwargs...) where id = get(kwargs, id, nothing)
@inline (termmodulate::TermModulate{<:Function})(args...; kwargs...) = termmodulate.modulate(args...; kwargs...)
@inline ismodulatable(termmodulate::TermModulate) = ismodulatable(typeof(termmodulate))
@inline ismodulatable(::Type{<:TermModulate{Val{B}}}) where B = B
@inline ismodulatable(::Type{<:TermModulate{<:Function}}) = true

"""
    Term{K, I, V, B, C<:TermCoupling, A<:TermAmplitude, M<:TermModulate}

A term of a quantum lattice system.
"""
mutable struct Term{K, I, V, B, C<:TermCoupling, A<:TermAmplitude, M<:TermModulate}
    value::V
    const bondkind::B
    const coupling::C
    const amplitude::A
    const ishermitian::Bool
    const modulate::M
    const factor::V
    function Term{K, I}(value, bondkind, coupling::TermCoupling, amplitude::TermAmplitude, ishermitian::Bool, modulate::TermModulate, factor) where {K, I}
        @assert isa(K, Symbol) "Term error: kind must be a Symbol."
        @assert isa(I, Symbol) "Term error: id must be a Symbol."
        @assert value==value' "Term error: only real values are allowed. Complex values should be specified by the amplitude function."
        V, B, C, A, M = typeof(value), typeof(bondkind), typeof(coupling), typeof(amplitude), typeof(modulate)
        new{K, I, V, B, C, A, M}(value, bondkind, coupling, amplitude, ishermitian, modulate, factor)
    end
end
@inline Base.:(==)(term₁::Term, term₂::Term) = ==(efficientoperations, term₁, term₂)
@inline Base.isequal(term₁::Term, term₂::Term) = isequal(efficientoperations, term₁, term₂)

"""
    Term{K}(id::Symbol, value, bondkind, coupling, ishermitian::Bool; amplitude::Union{Function, Nothing}=nothing, modulate::Union{Function, Bool}=true) where K

Construct a term.
"""
@inline function Term{K}(id::Symbol, value, bondkind, coupling, ishermitian::Bool; amplitude::Union{Function, Nothing}=nothing, modulate::Union{Function, Bool}=true) where K
    return Term{K, id}(value, bondkind, TermCoupling(coupling), TermAmplitude(amplitude), ishermitian, TermModulate(id, modulate), 1)
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
    valtype(term::Term)
    valtype(::Type{<:Term)

Get the value type of a term.
"""
@inline Base.valtype(term::Term) = valtype(typeof(term))
@inline Base.valtype(::Type{<:Term{K, I, V} where {K, I}}) where V = V

"""
    rank(term::Term) -> Int
    rank(::Type{<:Term) -> Int

Get the rank of a term.
"""
@inline rank(term::Term) = rank(typeof(term))
@inline rank(::Type{<:Term{K, I, V, B, C} where {K, I, V, B}}) where {C<:TermCoupling} = rank(valtype(C))

"""
    ismodulatable(term::Term) -> Bool
    ismodulatable(::Type{<:Term}) -> Bool

Judge whether a term could be modulated by its modulate function.
"""
@inline ismodulatable(term::Term) = ismodulatable(typeof(term))
@inline ismodulatable(::Type{<:Term{K, I, V, B, <:TermCoupling, <:TermAmplitude, M} where {K, I, V, B}}) where M = ismodulatable(M)

"""
    repr(term::Term, bond::Bond, hilbert::Hilbert) -> String

Get the repr representation of a term on a bond with a given Hilbert space.
"""
function Base.repr(term::Term, bond::Bond, hilbert::Hilbert)
    cache = String[]
    if term.bondkind == bond.kind
        value = term.value * term.amplitude(bond) * term.factor
        if !isapprox(value, 0.0, atol=atol, rtol=rtol)
            for coupling in term.coupling(bond)
                if !isnothing(iterate(expand(term|>kind|>Val, coupling, bond, hilbert)))
                    representation = repr(value*coupling)
                    term.ishermitian || (representation = string(replace(representation, " * "=>"*"), " + h.c."))
                    push!(cache, @sprintf "%s" representation)
                end
            end
        end
    end
    return join(cache, "\n")
end

"""
    replace(term::Term; kwargs...) -> Term

Replace some attributes of a term with key word arguments.
"""
@inline @generated function Base.replace(term::Term; kwargs...)
    exprs = [:(get(kwargs, $name, getfield(term, $name))) for name in QuoteNode.(term|>fieldnames)]
    return :(Term{kind(term), id(term)}($(exprs...)))
end

"""
    one(term::Term) -> Term

Get a unit term.
"""
@inline Base.one(term::Term) = replace(term, value=one(term.value))

"""
    zero(term::Term) -> Term

Get a zero term.
"""
@inline Base.zero(term::Term) = replace(term, value=zero(term.value))

"""
    update!(term::Term, args...; kwargs...) -> Term

Update the value of a term by its `modulate` function.
"""
function update!(term::Term, args...; kwargs...)
    @assert ismodulatable(term) "update! error: not modulatable term."
    value = term.modulate(args...; kwargs...)
    isnothing(value) || (term.value = value)
    return term
end

"""
    optype(::Type{T}, ::Type{H}, ::Type{B}) where {T<:Term, H<:Hilbert, B<:Bond}

Get the compatible `Operator` type from the type of a term, a Hilbert space and a bond.
"""
@inline function optype(::Type{T}, ::Type{H}, ::Type{B}) where {T<:Term, H<:Hilbert, B<:Bond}
    C = valtype(fieldtype(T, :coupling))
    @assert C<:Coupling "optype error: not supported."
    indextypes = ntuple(i->indextype(filter(iidtype(fieldtype(parametertype(C, :indexes), i)) , valtype(H)), eltype(B), Val(kind(T))), Val(rank(C)))
    return fulltype(Operator, NamedTuple{(:value, :id), Tuple{valtype(T), Tuple{indextypes...}}})
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
        if !isapprox(value, 0.0, atol=atol, rtol=rtol)
            for coupling in term.coupling(bond)
                for opt in expand(term|>kind|>Val, coupling, bond, hilbert)
                    isapprox(opt.value, 0.0, atol=atol, rtol=rtol) && continue
                    if half
                        add!(operators, opt*valtype(eltype(operators))(value/(term.ishermitian ? 2 : 1)))
                    else
                        opt = opt*valtype(eltype(operators))(value)
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
    M = optype(term|>typeof, hilbert|>typeof, bond|>typeof)
    expand!(Operators{M}(), term, bond, hilbert; half=half)
end
@inline function expand(term::Term, bonds, hilbert::Hilbert; half::Bool=false)
    M = optype(term|>typeof, hilbert|>typeof, bonds|>eltype)
    expand!(Operators{M}(), term, bonds, hilbert; half=half)
end

"""
    sitestructure(::Val{termkind}, ::Val{termrank}, bondlength::Integer) where {termkind, termrank} -> NTuple{termrank, Int}

Get the site structure, i.e. the acting centers of the coupling on a bond, of a certain kind of term.
"""
function sitestructure(::Val, ::Val{termrank}, bondlength::Integer) where termrank
    bondlength==1 && return ntuple(i->1, Val(termrank))
    bondlength==2 && begin
        termrank==2 && return (1, 2)
        termrank==4 && return (1, 1, 2, 2)
        error("sitestructure error: not supported for rank-$termrank.")
    end
    error("sitestructure error: not supported for a generic bond containing $bondlength points.")
end

"""
    expand(::Val{termkind}, coupling::Coupling, bond::Bond, hilbert::Hilbert) where termkind

Expand a coupling with the given bond and Hilbert space of a certain kind of term.
"""
function expand(::Val{termkind}, coupling::Coupling, bond::Bond, hilbert::Hilbert) where termkind
    centers = isa(coupling.indexes.sites, Tuple{Vararg{Int}}) ? coupling.indexes.sites : sitestructure(Val(termkind), Val(rank(coupling)), length(bond))
    points = NTuple{rank(coupling), eltype(bond)}(bond[centers[i]] for i = 1:rank(coupling))
    internals = map((iid, internal)->filter(iid, internal), CompositeIID(coupling, Val(termkind)).contents, ntuple(i->hilbert[bond[centers[i]].site], Val(rank(coupling))))
    @assert rank(coupling)==length(points)==length(internals) "expand error: mismatched rank."
    return CExpand(value(coupling), points, IIDSpace(CompositeIID(coupling, Val(termkind)), CompositeInternal{:⊗}(internals)), Constraint(coupling, Val(termkind)))
end
struct CExpand{V, N, SV<:SVector, S<:IIDSpace, C<:Constraint}
    value::V
    sites::NTuple{N, Int}
    rcoordinates::NTuple{N, SV}
    icoordinates::NTuple{N, SV}
    iidspace::S
    constraint::C
end
function CExpand(value, points::NTuple{N, P}, iidspace::IIDSpace, constraint::Constraint) where {N, P<:Point}
    sites = NTuple{N, Int}(points[i].site for i = 1:N)
    rcoordinates = NTuple{N, SVector{dimension(P), dtype(P)}}(points[i].rcoordinate for i = 1:N)
    icoordinates = NTuple{N, SVector{dimension(P), dtype(P)}}(points[i].icoordinate for i = 1:N)
    return CExpand(value, sites, rcoordinates, icoordinates, iidspace, constraint)
end
@inline Base.eltype(ex::CExpand) = eltype(typeof(ex))
@inline @generated function Base.eltype(::Type{<:CExpand{V, N, SV, S}}) where {V, N, SV<:SVector, S<:IIDSpace}
    return Operator{V, Tuple{map(I->CompositeIndex{Index{Int, I}, SV}, fieldtypes(fieldtype(eltype(S), :contents)))...}}
end
@inline Base.IteratorSize(::Type{<:CExpand}) = Base.SizeUnknown()
function Base.iterate(ex::CExpand, state=iterate(ex.iidspace))
    result = nothing
    while !isnothing(state)
        ciid, state = state
        indexes = ID(Index, ex.sites, ciid.contents)
        if match(ex.constraint, indexes)
            result = Operator(ex.value, ID(CompositeIndex, indexes, ex.rcoordinates, ex.icoordinates)), iterate(ex.iidspace, state)
            break
        else
            state = iterate(ex.iidspace, state)
        end
    end
    return result
end

# Metric and Table
"""
    Metric <: Function

The rules for measuring an operator unit so that different operator units can be compared.

As a function, every instance should accept only one positional argument, i.e. the operator unit to be measured.
"""
abstract type Metric <: Function end
@inline Base.:(==)(m₁::T, m₂::T) where {T<:Metric} = ==(efficientoperations, m₁, m₂)
@inline Base.isequal(m₁::T, m₂::T) where {T<:Metric} = isequal(efficientoperations, m₁, m₂)
@inline (M::Type{<:Metric})(::Type{I}) where {I<:AbstractCompositeIndex} = M(indextype(I))
@inline (metric::Metric)(index::AbstractCompositeIndex) = metric(getcontent(index, :index))
@inline Base.valtype(::Type{M}, ::Type{I}) where {M<:Metric, I<:AbstractCompositeIndex} = valtype(M, indextype(I))
@inline (M::Type{<:Metric})(::Type{H}) where {H<:Hilbert} = M(Index{Int, H|>valtype|>eltype})

"""
    OperatorUnitToTuple{Fields} <: Metric

A rule that converts an operator unit to a tuple by iterating over a set of selected fields in a specific order.
"""
struct OperatorUnitToTuple{Fields} <: Metric
    OperatorUnitToTuple(fields::Tuple{Vararg{Symbol}}) = new{fields}()
end
@inline OperatorUnitToTuple(fields::Symbol...) = OperatorUnitToTuple(fields)

"""
    keys(::OperatorUnitToTuple{Fields}) where Fields -> Fields
    keys(::Type{<:OperatorUnitToTuple{Fields}}) where Fields -> Fields

Get the names of the selected fields.
"""
@inline Base.keys(::OperatorUnitToTuple{Fields}) where Fields = Fields
@inline Base.keys(::Type{<:OperatorUnitToTuple{Fields}}) where Fields = Fields

"""
    OperatorUnitToTuple(::Type{I}) where {I<:Index}

Construct the metric rule from the information of the `Index` type.
"""
@inline OperatorUnitToTuple(::Type{I}) where {I<:Index} = OperatorUnitToTuple(:site, (fieldnames(iidtype(I)))...)

"""
    valtype(::Type{<:OperatorUnitToTuple}, ::Type{<:Index})

Get the valtype of applying an `OperatorUnitToTuple` rule to an `Index`.
"""
@inline @generated function Base.valtype(::Type{M}, ::Type{I}) where {M<:OperatorUnitToTuple, I<:Index}
    types = []
    for field in keys(M)
        if field==:site
            push!(types, Int)
        elseif hasfield(iidtype(I), field)
            push!(types, fieldtype(iidtype(I), field))
        end
    end
    return  Expr(:curly, :Tuple, types...)
end

"""
    (operatorunittotuple::OperatorUnitToTuple)(index::Index) -> Tuple

Convert an index to a tuple.
"""
@inline @generated function (operatorunittotuple::OperatorUnitToTuple)(index::Index)
    exprs = []
    for name in keys(operatorunittotuple)
        field = QuoteNode(name)
        if name==:site
            push!(exprs, :(index.site))
        elseif hasfield(iidtype(index), name)
            push!(exprs, :(getfield(index.iid, $field)))
        end
    end
    return Expr(:tuple, exprs...)
end

"""
    Table{I, B<:Metric} <: CompositeDict{I, Int}

The table of operator unit vs. sequence pairs.
"""
struct Table{I, B<:Metric} <: CompositeDict{I, Int}
    by::B
    contents::Dict{I, Int}
end
@inline contentnames(::Type{<:Table}) = (:by, :contents)
@inline Table{I}(by::Metric) where {I<:OperatorUnit} = Table(by, Dict{valtype(typeof(by), I), Int}())
@inline vec2dict(vs::AbstractVector) = Dict{eltype(vs), Int}(v=>i for (i, v) in enumerate(vs))

"""
    getindex(table::Table, operatorunit::OperatorUnit) -> Int

Inquiry the sequence of an operator unit.
"""
@inline Base.getindex(table::Table, operatorunit::OperatorUnit) = table[table.by(operatorunit)]

"""
    haskey(table::Table, operatorunit::OperatorUnit) -> Bool
    haskey(table::Table, operatorunits::ID{OperatorUnit}) -> Tuple{Vararg{Bool}}

Judge whether a single operator unit or a set of operator units have been assigned with sequences in table.
"""
@inline Base.haskey(table::Table, operatorunit::OperatorUnit) = haskey(table, table.by(operatorunit))
@inline Base.haskey(table::Table, operatorunits::ID{OperatorUnit}) = map(operatorunit->haskey(table, operatorunit), operatorunits)

"""
    Table(operatorunits::AbstractVector{<:OperatorUnit}, by::Metric=OperatorUnitToTuple(eltype(operatorunits)))

Convert a set of operator units to the corresponding table of operator unit vs. sequence pairs.

The input operator units are measured by the input `by` function with the duplicates removed. The resulting unique values are sorted, which determines the sequence of the input `operatorunits`. Note that two operator units have the same sequence if their converted values are equal to each other.
"""
@inline Table(operatorunits::AbstractVector{<:OperatorUnit}, by::Metric=OperatorUnitToTuple(eltype(operatorunits))) = Table(by, [by(operatorunit) for operatorunit in operatorunits]|>unique!|>sort!|>vec2dict)

"""
    Table(hilbert::Hilbert, by::Metric=OperatorUnitToTuple(typeof(hilbert))) -> Table

Get the index-sequence table of a Hilbert space.
"""
function Table(hilbert::Hilbert, by::Metric=OperatorUnitToTuple(typeof(hilbert)))
    result = Index{Int, hilbert|>valtype|>eltype}[]
    for (site, internal) in hilbert
        for iid in internal
            push!(result, (result|>eltype)(site, iid))
        end
    end
    return Table(result, by)
end

"""
    union(tables::Table...) -> Table

Unite several operator unit vs. sequence tables.
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
    reset!(table::Table, operatorunits::AbstractVector{<:OperatorUnit}) -> Table

Reset a table by a new set of operatorunits.
"""
function reset!(table::Table, operatorunits::AbstractVector{<:OperatorUnit})
    empty!(table)
    for (i, id) in enumerate([table.by(operatorunit) for operatorunit in operatorunits]|>unique!|>sort!)
        table[id] = i
    end
    return table
end

"""
    reset!(table::Table, hilbert::Hilbert) -> Table

Reset a table by a Hilbert space.
"""
function reset!(table::Table, hilbert::Hilbert)
    indices = Index{Int, hilbert|>valtype|>eltype}[]
    for (site, internal) in pairs(hilbert)
        for iid in internal
            push!(indices, (indices|>eltype)(site, iid))
        end
    end
    return reset!(table, indices)
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
@inline Base.valtype(::Type{<:Boundary}, M::Type{<:Operator}) = reparameter(M, :value, promote_type(Complex{Int}, valtype(M)))
@inline Base.valtype(B::Type{<:Boundary}, MS::Type{<:Operators}) = (M = valtype(B, eltype(MS)); Operators{M, idtype(M)})

"""
    keys(bound::Boundary) -> Tuple{Vararg{Symbol}}
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
