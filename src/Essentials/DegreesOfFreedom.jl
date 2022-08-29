module DegreesOfFreedom

using MacroTools: @capture
using Printf: @printf, @sprintf
using StaticArrays: SVector
using ..QuantumOperators: ID, Operator, OperatorPack, Operators, OperatorSum, OperatorUnit, Transformation, idtype, valuetolatextext
using ..Spatials: Bond, Point
using ...Essentials: dtype
using ...Interfaces: add!, decompose, dimension
using ...Prerequisites: atol, rtol, Float, decimaltostr
using ...Prerequisites.Traits: efficientoperations, commontype, fulltype, parametertype, rawtype, reparameter
using ...Prerequisites.CompositeStructures: CompositeDict, CompositeTuple, NamedContainer
using ...Prerequisites.VectorSpaces: VectorSpace, VectorSpaceCartesian, VectorSpaceDirectProducted, VectorSpaceDirectSummed, VectorSpaceStyle

import LinearAlgebra: ishermitian
import ..QuantumOperators: optype, script
import ..Spatials: icoordinate, rcoordinate
import ...Essentials: kind, reset!, update!
import ...Interfaces: ⊕, ⊗, expand, expand!, id, rank, value
import ...Prerequisites.Traits: contentnames, getcontent, isparameterbound, parameternames

export CompositeIID, CompositeInternal, IID, Internal, SimpleIID, SimpleInternal
export AbstractCompositeIndex, CompositeIndex, Hilbert, Index, iidtype, indextype, ishermitian, statistics
export wildcard, IIDSpace, Subscript, Subscripts, diagonal, noconstrain, subscriptexpr, @subscript_str
export AbstractCoupling, Coupling, Couplings, couplingcenters, couplinginternals, couplingpoints, @couplings
export Metric, OperatorUnitToTuple, Table
export Term, TermAmplitude, TermCouplings, TermFunction, TermModulate, ismodulatable
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
    Index{I<:SimpleIID} <: OperatorUnit

The index of a degree of freedom, which consist of the spatial part and the internal part.
"""
struct Index{I<:SimpleIID} <: OperatorUnit
    site::Int
    iid::I
end
@inline parameternames(::Type{<:Index}) = (:iid,)
@inline isparameterbound(::Type{<:Index}, ::Val{:iid}, ::Type{I}) where {I<:SimpleIID} = !isconcretetype(I)

"""
    iidtype(index::Index)
    iidtype(::Type{I}) where {I<:Index}

Get the type of the internal part of an index.
"""
@inline iidtype(index::Index) = iidtype(typeof(index))
@inline iidtype(::Type{I}) where {I<:Index} = parametertype(I, 1)

"""
    statistics(index::Index) -> Symbol
    statistics(::Type{<:Index{I}}) where {I<:SimpleIID} -> Symbol

Get the statistics of an index.
"""
@inline statistics(index::Index) = statistics(typeof(index))
@inline statistics(::Type{<:Index{I}}) where {I<:SimpleIID} = statistics(I)

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
    return fulltype(CompositeIndex, NamedTuple{(:index, :coordination), Tuple{fulltype(Index, NamedTuple{(:iid,), Tuple{eltype(I)}}), SVector{dimension(P), dtype(P)}}})
end

"""
    rcoordinate(opt::Operator{<:Number, <:ID{CompositeIndex}}) -> SVector

Get the whole rcoordinate of an operator.
"""
@inline function rcoordinate(opt::Operator{<:Number, <:ID{CompositeIndex}})
    rank(opt)==1 && return id(opt)[1].rcoordinate
    rank(opt)==2 && return id(opt)[1].rcoordinate-id(opt)[2].rcoordinate
    error("rcoordinate error: not supported rank($(rank(opt))) of $(nameof(opt)).")
end

"""
    icoordinate(opt::Operator{<:Number, <:ID{CompositeIndex}}) -> SVector

Get the whole icoordinate of an operator.
"""
@inline function icoordinate(opt::Operator{<:Number, <:ID{CompositeIndex}})
    rank(opt)==1 && return id(opt)[1].icoordinate
    rank(opt)==2 && return id(opt)[1].icoordinate-id(opt)[2].icoordinate
    error("icoordinate error: not supported rank($(rank(opt))) of $(nameof(opt)).")
end

"""
    script(::Val{:rcoordinate}, index::CompositeIndex; kwargs...) -> String
    script(::Val{:icoordinate}, index::CompositeIndex; kwargs...) -> String

Get the `:rcoordinate/:icoordinate` script of a composite index.
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

# Hilbert
"""
    Hilbert{I<:Internal} <: CompositeDict{Int, I}

Hilbert space at a lattice.
"""
struct Hilbert{I<:Internal} <: CompositeDict{Int, I}
    map::Function
    contents::Dict{Int, I}
end
@inline contentnames(::Type{<:Hilbert}) = (:map, :contents)

"""
    Hilbert(ps::Pair...)
    Hilbert(kv)

Construct a Hilbert space the same way as a Dict.
"""
@inline Hilbert(ps::Pair...) = Hilbert(ps)
function Hilbert(kv)
    contents = Dict(kv)
    map = site->contents[site]
    return Hilbert(map, contents)
end

"""
    Hilbert(map::Function, sites::AbstractVector{Int})
    Hilbert{I}(map::Function, sites::AbstractVector{Int}) where {I<:Internal}

Construct a Hilbert space from a function and a set of sites.

Here, `map` maps an integer to an `Internal`.
"""
@inline function Hilbert(map::Function, sites::AbstractVector{Int})
    I = commontype(map, Tuple{Vararg{Any}}, Internal)
    return Hilbert{I}(map, sites)
end
function Hilbert{I}(map::Function, sites::AbstractVector{Int}) where {I<:Internal}
    contents = Dict{sites|>eltype, I}()
    for site in sites
        contents[site] = map(site)
    end
    return Hilbert(map, contents)
end

"""
    reset!(hilbert::Hilbert, sites::AbstractVector{Int}) -> Hilbert

Reset the Hilbert space with new sites.
"""
function reset!(hilbert::Hilbert, sites::AbstractVector{Int})
    empty!(hilbert)
    for site in sites
        hilbert[site] = hilbert.map(site)::valtype(hilbert)
    end
    return hilbert
end

# Coupling and Couplings
"""
    IIDSpace{I<:IID, V<:Internal, Kind} <: VectorSpace{IID}

The space expanded by a "labeled" iid.
"""
struct IIDSpace{I<:IID, V<:Internal, Kind} <: VectorSpace{IID}
    iid::I
    internal::V
    IIDSpace(iid::IID, internal::Internal, ::Val{Kind}=Val(:info)) where Kind = new{typeof(iid), typeof(internal), Kind}(iid, internal)
end
@inline kind(iidspace::IIDSpace) = kind(typeof(iidspace))
@inline kind(::Type{<:IIDSpace{<:IID, <:Internal, Kind}}) where Kind = Kind
@inline Base.eltype(::Type{<:IIDSpace{<:IID, V}}) where {V<:Internal} = eltype(V)

@inline VectorSpaceStyle(::Type{<:IIDSpace{<:SimpleIID, <:SimpleInternal}}) = VectorSpaceCartesian()
@inline Base.CartesianIndex(iid::SimpleIID, iidspace::IIDSpace{<:SimpleIID, <:SimpleInternal}) = CartesianIndex(iid, iidspace.internal)
@inline function Base.getindex(::VectorSpaceCartesian, iidspace::IIDSpace{<:SimpleIID, <:SimpleInternal}, index::CartesianIndex)
    return rawtype(eltype(iidspace))(index, iidspace.internal)
end

@inline VectorSpaceStyle(::Type{<:IIDSpace{<:CompositeIID, <:CompositeInternal}}) = VectorSpaceDirectProducted()
@inline function getcontent(iidspace::IIDSpace{<:CompositeIID, <:CompositeInternal}, ::Val{:contents})
    return map((iid, internal)->IIDSpace(iid, internal, Val(kind(iidspace))), CompositeIID(iidspace).contents, CompositeInternal(iidspace).contents)
end
@inline CompositeIID(iids::Tuple{Vararg{SimpleIID}}, ::IIDSpace{<:CompositeIID, <:CompositeInternal}) = CompositeIID(iids)
@inline CompositeIID(iidspace::IIDSpace{<:CompositeIID, <:CompositeInternal}) = iidspace.iid
@inline CompositeInternal(iidspace::IIDSpace{<:CompositeIID, <:CompositeInternal}) = iidspace.internal

"""
    expand(iid::SimpleIID, internal::SimpleInternal) -> IIDSpace
    expand(iids::NTuple{N, SimpleIID}, internals::NTuple{N, SimpleInternal}) where N -> IIDSpace

Get the space expanded by a set of "labeled" iids.
"""
@inline expand(iid::SimpleIID, internal::SimpleInternal) = expand((iid,), (internal,))
@inline expand(iids::NTuple{N, SimpleIID}, internals::NTuple{N, SimpleInternal}) where N = IIDSpace(CompositeIID(iids), CompositeInternal{:⊗}(internals))

@inline diagonal(xs...) = length(xs)<2 ? true : all(map(==(xs[1]), xs))
@inline noconstrain(_...) = true
const wildcard = Symbol("*")
"""
    Subscript{P<:Tuple} <: CompositeTuple{P}

The subscript representative of a certain internal degree of freedom.
"""
struct Subscript{P<:Tuple} <: CompositeTuple{P}
    pattern::P
    rep::String
    constraint::Function
end
@inline contentnames(::Type{<:Subscript}) = (:contents, :rep, :constraint)
@inline getcontent(subscript::Subscript, ::Val{:contents}) = subscript.pattern
@inline Base.:(==)(subs₁::Subscript, subs₂::Subscript) = subs₁.pattern==subs₂.pattern && subs₁.rep==subs₂.rep
@inline Base.:isequal(subs₁::Subscript, subs₂::Subscript) = isequal(subs₁.pattern, subs₂.pattern) && isequal(subs₁.rep, subs₂.rep)
function Base.show(io::IO, subscript::Subscript)
    if subscript.rep ∈ ("diagonal", "noconstrain", "constant")
        @printf io "[%s]" join(subscript.pattern, " ")
    else
        @printf io "%s" subscript.rep
    end
end

"""
    Subscript(N::Int)
    Subscript(pattern::Tuple, check_constant::Bool=false)

Construct the subscript representative of a certain internal degree of freedom.
"""
@inline Subscript(N::Int) = Subscript(Val(N))
@inline Subscript(::Val{N}) where N = Subscript(ntuple(i->wildcard, Val(N)), "diagonal", diagonal)
@inline Subscript(pattern::Tuple, check_constant::Bool=false) = Subscript(pattern, Val(check_constant))
@inline function Subscript(pattern::Tuple, ::Val{false})
    any(map(p->isa(p, Symbol), pattern)) && error("Subscript error: wrong constant pattern.")
    return Subscript(pattern, "noconstrain", noconstrain)
end
@inline function Subscript(pattern::Tuple, ::Val{true})
    any(map(p->isa(p, Symbol), pattern)) && error("Subscript error: wrong constant pattern.")
    return Subscript(pattern, "constant", (xs...)->xs==pattern)
end

"""
    rank(subscript::Subscript) -> Int
    rank(::Type{<:Subscript}) -> Int

Get the number of the whole variables of the subscript.
"""
@inline rank(subscript::Subscript) = rank(typeof(subscript))
@inline rank(::Type{T}) where {T<:Subscript} = length(T)

"""
    match(subscript::Subscript, values::Tuple) -> Bool

Judge whether a set of values matches the pattern specified by the subscript.
"""
@inline function Base.match(subscript::Subscript, values::Tuple)
    @assert length(subscript)==length(values) "match error: mismatched length of values."
    return subscript.constraint(values...)::Bool
end

"""
    subscript"..." -> Subscript

Construct the subscript from a literal string.
"""
macro subscript_str(str)
    expr = Meta.parse(str)
    expr.head==:toplevel || return subscriptexpr(expr)
    @assert length(expr.args)==2 && isa(expr.args[2], Bool) "@subscript_str error: wrong pattern."
    return subscriptexpr(expr.args[1], expr.args[2])
end
function subscriptexpr(expr::Expr, check_constant::Bool=false)
    if @capture(expr, op_(cp_))
        @assert op.head∈(:hcat, :vect) "subscriptexpr error: wrong pattern."
        pattern, condition = Tuple(op.args), cp
        rep = @sprintf "[%s](%s)" join(pattern, " ") condition
    else
        @assert expr.head∈(:hcat, :vect) "subscriptexpr error: wrong pattern."
        pattern, condition = Tuple(expr.args), true
        rep = @sprintf "[%s]" join(pattern, " ")
    end
    if !any(map(p->isa(p, Symbol), pattern))
        check_constant && return :(Subscript($pattern, Val(true)))
        return :(Subscript($pattern, Val(false)))
    end
    paramargs, groups = Symbol[], Dict{Symbol, Vector{Symbol}}()
    for sub in pattern
        isa(sub, Symbol) || begin
            paramarg = gensym("paramarg")
            push!(paramargs, paramarg)
            condition = Expr(Symbol("&&"), condition, Expr(:call, :(==), paramarg, sub))
            continue
        end
        if sub∉paramargs
            push!(paramargs, sub)
            groups[sub]=[sub]
        else
            paramarg = gensym("paramarg")
            push!(paramargs, paramarg)
            push!(groups[sub], :(==))
            push!(groups[sub], paramarg)
        end
    end
    for group in values(groups)
        length(group)==1 && continue
        condition = Expr(Symbol("&&"), condition, Expr(:comparison, group...))
    end
    name = gensym("subconstraint")
    constraint = :($name($(paramargs...)) = $condition)
    return Expr(:block, constraint, :(Subscript($pattern, $rep, $name)))
end

"""
    Subscripts{T<:Tuple{Vararg{NamedContainer{Subscript}}}, R<:Tuple{Vararg{Pair{UnitRange{Int}, <:Tuple{Vararg{String}}}}}} <: OperatorUnit

The complete set of subscripts of the internal degrees of freedom.
"""
struct Subscripts{T<:Tuple{Vararg{NamedContainer{Subscript}}}, R<:Tuple{Vararg{Pair{UnitRange{Int}, <:Tuple{Vararg{String}}}}}} <: OperatorUnit
    contents::T
    rep::R
end
@inline Base.:(==)(subs₁::Subscripts, subs₂::Subscripts) = subs₁.contents==subs₂.contents
@inline Base.hash(subscripts::Subscripts, h::UInt) = hash(subscripts.rep, h)
@inline Base.length(subscripts::Subscripts) = length(typeof(subscripts))
@inline Base.length(::Type{<:Subscripts{T}}) where {T<:Tuple{Vararg{NamedContainer{Subscript}}}} = fieldcount(T)
@inline Base.getindex(subscripts::Subscripts, i::Integer) = subscripts.contents[i]
@inline Base.iterate(subscripts::Subscripts, state=1) = state>length(subscripts) ? nothing : (subscripts[state], state+1)
function Base.show(io::IO, subscripts::Subscripts)
    for (i, segment) in enumerate(subscripts.contents)
        i>1 && @printf io "%s" " × "
        for (j, (field, subscript)) in enumerate(pairs(segment))
            j>1 && @printf io "%s" " ⊗ "
            @printf io "%s%s" field subscript
        end
    end
end
function Base.repr(subscripts::Subscripts, slice, field::Symbol)
    result = []
    for (i, segment) in enumerate(slice)
        i>1 && push!(result, "×")
        push!(result, @sprintf "%s" getfield(subscripts.contents[segment], field))
    end
    return join(result)
end

"""
    Subscripts(contents::NamedContainer{Subscript}...)

Construct the complete set of subscripts.
"""
function Subscripts(contents::NamedContainer{Subscript}...)
    for segment in contents
        length(segment)>1 && @assert mapreduce(length, ==, values(segment)) "Subscripts error: mismatched ranks."
    end
    reps = map(content->map(subscript->subscript.rep, values(content)), contents)
    counts = (0, cumsum(map(content->rank(first(content)), contents))...)
    return Subscripts(contents, map((i, rep)->(counts[i]+1:counts[i+1])=>rep, ntuple(i->i, Val(fieldcount(typeof(reps)))), reps))
end

"""
    rank(subscripts::Subscripts) -> Int
    rank(::Type{<:Subscripts{T}}) where {T<:Tuple{Vararg{NamedContainer{Subscript}}}} -> Int

Get the rank of the subscripts.
"""
@inline rank(subscripts::Subscripts) = rank(typeof(subscripts))
@inline @generated function rank(::Type{<:Subscripts{T}}) where {T<:Tuple{Vararg{NamedContainer{Subscript}}}}
    sum(rank(fieldtype(fieldtype(T, i), 1)) for i = 1:fieldcount(T))
end

"""
    rank(subscripts::Subscripts, i::Integer) -> Int
    rank(::Type{<:Subscripts{T}}, i::Integer) where {T<:Tuple{Vararg{NamedContainer{Subscript}}}} -> Int

Get the rank of the ith homogenous segment of the subscripts.
"""
@inline rank(subscripts::Subscripts, i::Integer) = rank(typeof(subscripts), i)
@inline rank(::Type{<:Subscripts{T}}, i::Integer) where {T<:Tuple{Vararg{NamedContainer{Subscript}}}} = rank(fieldtype(fieldtype(T, i), 1))

"""
    match(subscripts::Subscripts, iids::Tuple{Vararg{SimpleIID}}) -> Bool
    match(subscripts::Subscripts, ciid::CompositeIID) -> Bool

Judge whether a composite iid matches the patterns specified by the subscripts.
"""
@inline Base.match(subscripts::Subscripts, ciid::CompositeIID) = match(subscripts, ciid.contents)
@generated function Base.match(subscripts::Subscripts, iids::Tuple{Vararg{SimpleIID}})
    length(subscripts)==0 && return true
    @assert rank(subscripts)==fieldcount(iids) "match error: mismatched rank of iids and subscripts."
    exprs, count = [], 1
    for i = 1:length(subscripts)
        start, stop = count, count+rank(subscripts, i)-1
        for field in fieldnames(fieldtype(fieldtype(subscripts, :contents), i))
            field = QuoteNode(field)
            paramvalue = Expr(:tuple, [:(getfield(iids[$j], $field)) for j = start:stop]...)
            push!(exprs, :(match(getfield(subscripts[$i], $field), $paramvalue)))
        end
        count = stop+1
    end
    return Expr(:call, :all, Expr(:tuple, exprs...))
end

"""
    *(subscripts₁::Subscripts, subscripts₂::Subscripts) -> Subscripts

Get the combination of two sets of subscripts.
"""
@inline Base.:*(subscripts₁::Subscripts, subscripts₂::Subscripts) = Subscripts((subscripts₁.contents..., subscripts₂.contents...)...)

"""
    AbstractCoupling{V, I<:ID{OperatorUnit}} <: OperatorPack{V, I}

The abstract coupling intra/inter internal degrees of freedom at different lattice points.
"""
abstract type AbstractCoupling{V, I<:ID{OperatorUnit}} <: OperatorPack{V, I} end
@inline ID{SimpleIID}(coupling::AbstractCoupling) = id(coupling)
@inline Subscripts(coupling::AbstractCoupling) = Subscripts()

"""
    couplingcenters(coupling::AbstractCoupling, bond::Bond, info::Val) -> NTuple{rank(coupling), Int}

Get the acting centers of the coupling on a bond.
"""
function couplingcenters(coupling::AbstractCoupling, bond::Bond, info::Val)
    length(bond)==1 && return ntuple(i->1, Val(rank(coupling)))
    length(bond)==2 && begin
        rank(coupling)==2 && return (1, 2)
        rank(coupling)==4 && return (1, 1, 2, 2)
        error("couplingcenters error: not supported for a rank-$(rank(coupling)) coupling.")
    end
    error("couplingcenters error: not supported for a generic bond containing $(length(bond)) points.")
end

"""
    couplingpoints(coupling::AbstractCoupling, bond::Bond, info::Val) -> NTuple{rank(coupling), eltype(bond)}

Get the points where each order of the coupling acts on.
"""
@inline function couplingpoints(coupling::AbstractCoupling, bond::Bond, info::Val)
    centers = couplingcenters(coupling, bond, info)
    return NTuple{rank(coupling), eltype(bond)}(bond[centers[i]] for i = 1:rank(coupling))
end

"""
    couplinginternals(coupling::AbstractCoupling, bond::Bond, hilbert::Hilbert, info::Val) -> NTuple{rank(coupling), SimpleInternal}

Get the internal spaces where each order of the coupling acts on.
"""
@inline function couplinginternals(coupling::AbstractCoupling, bond::Bond, hilbert::Hilbert, info::Val)
    centers = couplingcenters(coupling, bond, info)
    internals = ntuple(i->hilbert[bond[centers[i]].site], Val(rank(coupling)))
    return map((iid, internal)->filter(iid, internal), ID{SimpleIID}(coupling), internals)
end

"""
    expand(coupling::AbstractCoupling, bond::Bond, hilbert::Hilbert, info::Val)

Expand a coupling with the given bond and Hilbert space.
"""
function expand(coupling::AbstractCoupling, bond::Bond, hilbert::Hilbert, info::Val)
    points = couplingpoints(coupling, bond, info)
    internals = couplinginternals(coupling, bond, hilbert, info)
    @assert rank(coupling)==length(points)==length(internals) "expand error: mismatched rank."
    return CExpand(value(coupling), points, IIDSpace(CompositeIID(ID{SimpleIID}(coupling)), CompositeInternal{:⊗}(internals), info), Subscripts(coupling))
end
struct CExpand{V, N, SV<:SVector, S<:IIDSpace, C<:Subscripts}
    value::V
    sites::NTuple{N, Int}
    rcoordinates::NTuple{N, SV}
    icoordinates::NTuple{N, SV}
    iidspace::S
    subscripts::C
end
function CExpand(value, points::NTuple{N, P}, iidspace::IIDSpace, subscripts::Subscripts) where {N, P<:Point}
    sites = NTuple{N, Int}(points[i].site for i = 1:N)
    rcoordinates = NTuple{N, SVector{dimension(P), dtype(P)}}(points[i].rcoordinate for i = 1:N)
    icoordinates = NTuple{N, SVector{dimension(P), dtype(P)}}(points[i].icoordinate for i = 1:N)
    return CExpand(value, sites, rcoordinates, icoordinates, iidspace, subscripts)
end
@inline Base.eltype(ex::CExpand) = eltype(typeof(ex))
@inline @generated function Base.eltype(::Type{<:CExpand{V, N, SV, S}}) where {V, N, SV<:SVector, S<:IIDSpace}
    return Operator{V, Tuple{map(I->CompositeIndex{Index{I}, SV}, fieldtypes(fieldtype(eltype(S), :contents)))...}}
end
@inline Base.IteratorSize(::Type{<:CExpand}) = Base.SizeUnknown()
function Base.iterate(ex::CExpand, state=iterate(ex.iidspace))
    result = nothing
    while !isnothing(state)
        ciid, state = state
        if match(ex.subscripts, ciid)
            result = Operator(ex.value, ID(CompositeIndex, ID(Index, ex.sites, ciid.contents), ex.rcoordinates, ex.icoordinates)), iterate(ex.iidspace, state)
            break
        else
            state = iterate(ex.iidspace, state)
        end
    end
    return result
end

"""
    Coupling{V, I<:ID{SimpleIID}, C<:Subscripts} <: AbstractCoupling{V, Tuple{CompositeIID{I}, C}}

The coupling intra/inter internal degrees of freedom at different lattice points.
"""
struct Coupling{V, I<:ID{SimpleIID}, C<:Subscripts} <: AbstractCoupling{V, Tuple{CompositeIID{I}, C}}
    value::V
    iids::I
    subscripts::C
    function Coupling(value::Number, iids::ID{SimpleIID}, subscripts::Subscripts=Subscripts())
        new{typeof(value), typeof(iids), typeof(subscripts)}(value, iids, subscripts)
    end
end
@inline parameternames(::Type{<:Coupling}) = (:value, :iids, :subscripts)
@inline isparameterbound(::Type{<:Coupling}, ::Val{:iids}, ::Type{I}) where {I<:ID{SimpleIID}} = !isconcretetype(I)
@inline isparameterbound(::Type{<:Coupling}, ::Val{:subscripts}, ::Type{C}) where {C<:Subscripts} = !isconcretetype(C)
@inline getcontent(coupling::Coupling, ::Val{:id}) = ID(CompositeIID(coupling.iids), coupling.subscripts)
@inline rank(::Type{<:Coupling{V, I} where V}) where {I<:ID{SimpleIID}} = fieldcount(I)
@inline Coupling(value::Number, id::Tuple{CompositeIID, Subscripts}) = Coupling(value, id[1].contents, id[2])
@inline ID{SimpleIID}(coupling::Coupling) = coupling.iids
@inline Subscripts(coupling::Coupling) = coupling.subscripts

"""
    *(cp₁::Coupling, cp₂::Coupling) -> Coupling

Get the multiplication between two couplings.
"""
@inline Base.:*(cp₁::Coupling, cp₂::Coupling) = Coupling(cp₁.value*cp₂.value, ID(cp₁.iids, cp₂.iids), cp₁.subscripts*cp₂.subscripts)

"""
    Couplings(cps::AbstractCoupling...)

A pack of couplings intra/inter internal degrees of freedom at different lattice points.

Alias for `OperatorSum{<:AbstractCoupling, <:ID{OperatorUnit}}`.
"""
const Couplings{C<:AbstractCoupling, I<:ID{OperatorUnit}} = OperatorSum{C, I}
@inline Couplings(cps::AbstractCoupling...) = OperatorSum(cps)
@inline Couplings(cps::Couplings) = cps

"""
    @couplings cps -> Couplings

Convert an expression/literal to a set of couplings.
"""
macro couplings(cps) :(Couplings($(esc(cps)))) end

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
@inline (M::Type{<:Metric})(::Type{H}) where {H<:Hilbert} = M(Index{H|>valtype|>eltype})

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
    result = Index{hilbert|>valtype|>eltype}[]
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
    @assert mapreduce(table->table.by, ==, tables) "union error: all input tables should have the same `by` attribute."
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
    indices = Index{hilbert|>valtype|>eltype}[]
    for (site, internal) in hilbert
        for iid in internal
            push!(indices, (indices|>eltype)(site, iid))
        end
    end
    return reset!(table, indices)
end

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
@inline (termamplitude::TermAmplitude{Nothing})(args...; kwargs...) = 1
@inline (termamplitude::TermAmplitude{<:Function})(args...; kwargs...) = termamplitude.amplitude(args...; kwargs...)

"""
    TermCouplings(couplings::Union{Couplings, Function})

The function for the couplings of a term.
"""
struct TermCouplings{C₁<:Union{Function, Couplings}, C₂<:Union{Couplings, Nothing}} <: TermFunction
    couplings::C₁
    TermCouplings(couplings::Couplings) = new{typeof(couplings), Nothing}(couplings)
    TermCouplings{C}(couplings::Function) where {C<:Couplings} = new{typeof(couplings), C}(couplings)
    TermCouplings(couplings::Function) = new{typeof(couplings), commontype(couplings, Tuple{Vararg{Any}}, Couplings)}(couplings)
end
@inline Base.valtype(termcouplings::TermCouplings) = valtype(typeof(termcouplings))
@inline Base.valtype(::Type{<:TermCouplings{C}}) where {C<:Couplings} = C
@inline Base.valtype(::Type{<:TermCouplings{<:Function, C}}) where {C<:Couplings} = C
@inline (termcouplings::TermCouplings{<:Couplings})(args...; kwargs...) = termcouplings.couplings
@inline (termcouplings::TermCouplings{<:Function})(args...; kwargs...) = termcouplings.couplings(args...; kwargs...)

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
    Term{K, I, V, B, C<:TermCouplings, A<:TermAmplitude, M<:TermModulate}

A term of a quantum lattice system.
"""
mutable struct Term{K, I, V, B, C<:TermCouplings, A<:TermAmplitude, M<:TermModulate}
    value::V
    bondkind::B
    couplings::C
    amplitude::A
    ishermitian::Bool
    modulate::M
    factor::V
    function Term{K, I}(value, bondkind, couplings::TermCouplings, amplitude::TermAmplitude, ishermitian::Bool, modulate::TermModulate, factor) where {K, I}
        @assert isa(K, Symbol) "Term error: kind must be a Symbol."
        @assert isa(I, Symbol) "Term error: id must be a Symbol."
        @assert value==value' "Term error: only real values are allowed. Complex values should be specified by the amplitude function."
        V, B, C, A, M = typeof(value), typeof(bondkind), typeof(couplings), typeof(amplitude), typeof(modulate)
        new{K, I, V, B, C, A, M}(value, bondkind, couplings, amplitude, ishermitian, modulate, factor)
    end
end
@inline Base.:(==)(term₁::Term, term₂::Term) = ==(efficientoperations, term₁, term₂)
@inline Base.isequal(term₁::Term, term₂::Term) = isequal(efficientoperations, term₁, term₂)

"""
    Term{K}(id::Symbol, value, bondkind; couplings::Union{Function, Couplings}, ishermitian::Bool, amplitude::Union{Function, Nothing}=nothing, modulate::Union{Function, Bool}=false) where K

Construct a term.
"""
@inline function Term{K}(id::Symbol, value, bondkind; couplings::Union{Function, Couplings}, ishermitian::Bool, amplitude::Union{Function, Nothing}=nothing, modulate::Union{Function, Bool}=false) where K
    return Term{K, id}(value, bondkind, TermCouplings(couplings), TermAmplitude(amplitude), ishermitian, TermModulate(id, modulate), 1)
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
@inline rank(::Type{<:Term{K, I, V, B, C} where {K, I, V, B}}) where {C<:TermCouplings} = rank(eltype(valtype(C)))

"""
    ismodulatable(term::Term) -> Bool
    ismodulatable(::Type{<:Term}) -> Bool

Judge whether a term could be modulated by its modulate function.
"""
@inline ismodulatable(term::Term) = ismodulatable(typeof(term))
@inline ismodulatable(::Type{<:Term{K, I, V, B, <:TermCouplings, <:TermAmplitude, M} where {K, I, V, B}}) where M = ismodulatable(M)

"""
    repr(term::Term, bond::Bond, hilbert::Hilbert) -> String

Get the repr representation of a term on a bond with a given Hilbert space.
"""
function Base.repr(term::Term, bond::Bond, hilbert::Hilbert)
    cache = String[]
    if term.bondkind == bond.kind
        value = term.value * term.amplitude(bond) * term.factor
        if !isapprox(value, 0.0, atol=atol, rtol=rtol)
            for coupling in term.couplings(bond)
                if !isnothing(iterate(expand(coupling, bond, hilbert, term|>kind|>Val)))
                    push!(cache, @sprintf "%s: %s%s" kind(term) repr(value*coupling) (term.ishermitian ? "" : " + h.c."))
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
    C = eltype(valtype(fieldtype(T, :couplings)))
    @assert C<:Coupling "optype error: not supported."
    indextypes = ntuple(i->indextype(filter(fieldtype(parametertype(C, :iids), i) , valtype(H)), eltype(B), Val(kind(T))), Val(rank(C)))
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
            for coupling in term.couplings(bond)
                for opt in expand(coupling, bond, hilbert, term|>kind|>Val)
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

# Boundary
"""
    Boundary{Names}(values::AbstractVector{<:Number}, vectors::AbstractVector{<:AbstractVector{<:Number}}) where Names

Boundary twist of operators.
"""
struct Boundary{Names, D<:Number, V<:AbstractVector} <: Transformation
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
