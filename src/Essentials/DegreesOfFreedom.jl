module DegreesOfFreedom

using MacroTools: @capture
using Printf: @printf, @sprintf
using StaticArrays: SVector
using LaTeXStrings: latexstring
using ..QuantumOperators: SingularID, ID, OperatorProd, OperatorSum
using ..Spatials: AbstractPID, AbstractBond, Point, Bonds, pidtype
using ...Essentials: dtype
using ...Interfaces: id, value, decompose, dimension, add!
using ...Prerequisites: Float, atol, rtol, decimaltostr
using ...Prerequisites.Traits: rawtype, efficientoperations, commontype, parametertype
using ...Prerequisites.CompositeStructures: CompositeDict, CompositeTuple, NamedContainer
using ...Prerequisites.VectorSpaces: CartesianVectorSpace

import LinearAlgebra: ishermitian
import ..QuantumOperators: sequence, idtype
import ..Spatials: pidtype, rcoord, icoord
import ...Essentials: kind, reset!, update!
import ...Interfaces: id, value, rank, expand, expand!, ⊗
import ...Prerequisites.Traits: parameternames, isparameterbound, contentnames, getcontent
import ...Prerequisites.VectorSpaces: shape, ndimshape

export IID, SimpleIID, CompositeIID, Internal, SimpleInternal, CompositeInternal
export AbstractOID, Index, CompositeOID, OID, Operator, Operators, IIDSpace, Hilbert
export statistics, iidtype, ishermitian, indextype, oidtype
export Subscript, Subscripts, SubscriptsID, @subscript_str, subscriptexpr, wildcard, diagonal, noconstrain
export AbstractCoupling, Coupling, Couplings, couplingcenters, couplingpoints, couplinginternals, @couplings
export Metric, OIDToTuple, Table
export TermFunction, TermAmplitude, TermCouplings, TermModulate, Term, ismodulatable, abbr, otype
export LaTeX, latexname, latexformat, superscript, subscript, script
export Boundary, twist, plain

"""
    IID <: SingularID

The id of an internal degree of freedom.
"""
abstract type IID <: SingularID end

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
iidtype(ciid::CompositeIID, i::Integer) = iidtype(typeof(ciid), i)
iidtype(::Type{<:CompositeIID{T}}, i::Integer) where {T<:Tuple{Vararg{SimpleIID}}} = fieldtype(T, i)

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
    Internal{I<:IID} <: CartesianVectorSpace{I}

The whole internal degrees of freedom at a single point.
"""
abstract type Internal{I<:IID} <: CartesianVectorSpace{I} end

"""
    SimpleInternal{I<:SimpleIID} <: Internal{I}

The simple internal degrees of freedom at a single point.
"""
abstract type SimpleInternal{I<:SimpleIID} <: Internal{I} end
Base.show(io::IO, i::SimpleInternal) = @printf io "%s(%s)" i|>typeof|>nameof join(("$name=$(getfield(i, name))" for name in i|>typeof|>fieldnames), ", ")

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
@inline @generated Base.match(::Type{I}, ::Type{SI}) where {I<:SimpleIID, SI<:SimpleInternal} = nameof(I)==nameof(eltype(SI))

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
    CompositeInternal{T<:Tuple{Vararg{SimpleInternal}}} <: Internal{CompositeIID}

The composition of several single internal spaces.
"""
struct CompositeInternal{T<:Tuple{Vararg{SimpleInternal}}} <: Internal{CompositeIID}
    contents::T
end
@inline Base.eltype(ci::CompositeInternal) = eltype(typeof(ci))
@inline @generated function Base.eltype(::Type{<:CompositeInternal{T}}) where {T<:Tuple{Vararg{SimpleInternal}}}
    return CompositeIID{Tuple{[eltype(fieldtype(T, i)) for i = 1:fieldcount(T)]...}}
end
Base.show(io::IO, ci::CompositeInternal) = @printf io "%s" join((string(ci.contents[i]) for i = 1:rank(ci)), " ⊗ ")

"""
    CompositeInternal(contents::SimpleInternal...)

Construct a composite internal space from a set of simple internal spaces.
"""
@inline CompositeInternal(contents::SimpleInternal...) = CompositeInternal(contents)

"""
    rank(ci::CompositeInternal) -> Int
    rank(::Type{<:CompositeInternal{T}}) where {T<:Tuple{Vararg{SimpleInternal}}} -> Int

Get the number of simple internal spaces in a composite internal space.
"""
@inline rank(ci::CompositeInternal) = rank(typeof(ci))
@inline rank(::Type{<:CompositeInternal{T}}) where {T<:Tuple{Vararg{SimpleInternal}}} = fieldcount(T)

@inline @generated shape(ci::CompositeInternal) = Expr(:tuple, [:(shape(ci.contents[$i])...) for i = 1:rank(ci)]...)
@inline ndimshape(::Type{<:CompositeInternal{T}}) where {T<:Tuple{Vararg{SimpleInternal}}} = sum(ndimshape(fieldtype(T, i)) for i = 1:fieldcount(T))
@inline @generated function Base.CartesianIndex(ciid::CompositeIID, ci::CompositeInternal)
    exprs = []
    for i = 1:rank(ci)
        push!(exprs, :(Tuple(CartesianIndex(ciid[$i], ci.contents[$i]))...))
    end
    return Expr(:call, :CartesianIndex, exprs...)
end
@inline CompositeIID(index::CartesianIndex, ci::CompositeInternal) = compositeiid(index, ci, segment(ci)|>Val)
@inline @generated function segment(ci::CompositeInternal{T}) where {T<:Tuple{Vararg{SimpleInternal}}}
    Expr(:tuple, [:(ndimshape(fieldtype(T, $i))) for i = 1:fieldcount(T)]...)
end
@inline @generated function compositeiid(index::CartesianIndex, ci::CompositeInternal, ::Val{dims}) where dims
    count = 1
    exprs = []
    for (i, dim) in enumerate(dims)
        cartesianindex = Expr(:call, :CartesianIndex, [:(index[$j]) for j = count:(count+dim-1)]...)
        push!(exprs, :(ci.contents[$i][$cartesianindex]))
        count += dim
    end
    return Expr(:call, :CompositeIID, exprs...)
end

"""
    ⊗(i₁::SimpleInternal, i₂::SimpleInternal) -> CompositeInternal
    ⊗(i::SimpleInternal, ci::CompositeInternal) -> CompositeInternal
    ⊗(ci::CompositeInternal, i::SimpleInternal) -> CompositeInternal
    ⊗(ci₁::CompositeInternal, ci₂::CompositeInternal) -> CompositeInternal

Direct product between simple internal spaces and composite internal spaces.
"""
@inline ⊗(i₁::SimpleInternal, i₂::SimpleInternal) = CompositeInternal(i₁, i₂)
@inline ⊗(i::SimpleInternal, ci::CompositeInternal) = CompositeInternal(i, ci.contents...)
@inline ⊗(ci::CompositeInternal, i::SimpleInternal) = CompositeInternal(ci.contents..., i)
@inline ⊗(ci₁::CompositeInternal, ci₂::CompositeInternal) = CompositeInternal(ci₁.contents..., ci₂.contents...)

"""
    filter(iid::SimpleIID, ci::CompositeInternal) -> Union{Nothing, SimpleInternal, CompositeInternal}
    filter(::Type{I}, ci::CompositeInternal) where {I<:SimpleIID} -> Union{Nothing, SimpleInternal, CompositeInternal}

Filter the composite internal space and select those that matches `I` or the type of `iid`.
"""
@inline Base.filter(iid::SimpleIID, ci::CompositeInternal) = filter(typeof(iid), ci)
@inline Base.filter(::Type{I}, ci::CompositeInternal) where {I<:SimpleIID} = filterhelper₁(I, ci, filtermatches(I, typeof(ci))|>Val)
@inline @generated function filtermatches(::Type{I}, ::Type{CompositeInternal{C}}) where {I<:SimpleIID, C<:Tuple{Vararg{SimpleInternal}}}
    exprs = []
    for i = 1:fieldcount(C)
        T = fieldtype(C, i)
        push!(exprs, :(match(I, $T)))
    end
    Expr(:tuple, exprs...)
end
@inline @generated function filterhelper₁(::Type{I}, ci::CompositeInternal, ::Val{BS}) where {I<:SimpleIID, BS}
    exprs = []
    for (i, B) in enumerate(BS)
        B && push!(exprs, :(ci.contents[$i]))
    end
    length(exprs)==0 && return
    length(exprs)==1 && return first(exprs)
    return Expr(:call, :CompositeInternal, exprs...)
end

"""
    filter(iid::SimpleIID, ::Type{C}) where {C<:CompositeInternal}
    filter(::Type{I}, ::Type{C}) where {I<:SimpleIID, C<:CompositeInternal}

Filter the type of a composite internal space and select those that matches `I` or the type of `iid`.
"""
@inline Base.filter(iid::SimpleIID, ::Type{C}) where {C<:CompositeInternal} = filter(typeof(iid), C)
@inline Base.filter(::Type{I}, ::Type{C}) where {I<:SimpleIID, C<:CompositeInternal} = filterhelper₂(I, C, filtermatches(I, C)|>Val)
@inline @generated function filterhelper₂(::Type{I}, ::Type{CompositeInternal{C}}, ::Val{BS}) where {I<:SimpleIID, C<:Tuple{Vararg{SimpleInternal}}, BS}
    exprs = []
    for (i, B) in enumerate(BS)
        B && push!(exprs, :(fieldtype(C, $i)))
    end
    length(exprs)==0 && return
    length(exprs)==1 && return first(exprs)
    return Expr(:curly, :CompositeInternal, Expr(:curly, :Tuple, exprs...))
end

"""
    AbstractOID <: SingularID

Abstract type of operator id.
"""
abstract type AbstractOID <: SingularID end

"""
    ishermitian(id::ID{AbstractOID, N}) where N -> Bool

Judge whether an operator id is Hermitian.
"""
function ishermitian(id::ID{AbstractOID, N}) where N
    for i = 1:((N+1)÷2)
        id[i]'==id[N+1-i] || return false
    end
    return true
end

"""
    Index{P<:AbstractPID, I<:SimpleIID} <: AbstractOID

The index of a degree of freedom, which consist of the spatial part and the internal part.
"""
struct Index{P<:AbstractPID, I<:SimpleIID} <: AbstractOID
    pid::P
    iid::I
end

"""
    pidtype(index::Index)
    pidtype(::Type{<:Index{P}}) where {P<:AbstractPID}

Get the type of the spatial part of an index.
"""
@inline pidtype(index::Index) = pidtype(typeof(index))
@inline pidtype(::Type{<:Index{P}}) where {P<:AbstractPID} = P

"""
    iidtype(index::Index)
    iidtype(::Type{<:Index{<:AbstractPID, I}}) where {I<:SimpleIID}

Get the type of the internal part of an index.
"""
@inline iidtype(index::Index) = iidtype(typeof(index))
@inline iidtype(::Type{<:Index{<:AbstractPID, I}}) where {I<:SimpleIID} = I

"""
    statistics(index::Index) -> Symbol
    statistics(::Type{<:Index{<:AbstractPID, I}}) where {I<:SimpleIID} -> Symbol

Get the statistics of an index.
"""
@inline statistics(index::Index) = statistics(typeof(index))
@inline statistics(::Type{<:Index{<:AbstractPID, I}}) where {I<:SimpleIID} = statistics(I)

"""
    adjoint(index::Index) -> typeof(index)

Get the adjoint of an index.
"""
@inline Base.adjoint(index::Index) = rawtype(typeof(index))(index.pid, adjoint(index.iid))

"""
    CompositeOID{I<:Index} <: AbstractOID

The abstract type of composite operator id.
"""
abstract type CompositeOID{I<:Index} <: AbstractOID end
@inline contentnames(::Type{<:CompositeOID}) = (:index,)

"""
    indextype(::CompositeOID)
    indextype(::Type{<:CompositeOID})

Get the index type of a composite operator id.
"""
@inline indextype(oid::CompositeOID) = indextype(typeof(oid))
@inline indextype(::Type{<:CompositeOID{I}}) where {I<:Index} = I

"""
    statistics(oid::CompositeOID) -> Symbol
    statistics(::Type{<:CompositeOID{I}}) where {I<:Index} -> Symbol

Get the statistics of a composite operator id.
"""
@inline statistics(oid::CompositeOID) = statistics(typeof(oid))
@inline statistics(::Type{<:CompositeOID{I}}) where {I<:Index} = statistics(I)

"""
    OID{I<:Index, V<:SVector} <: CompositeOID{I}

Operator id.
"""
struct OID{I<:Index, V<:SVector} <: CompositeOID{I}
    index::I
    rcoord::V
    icoord::V
    OID(index::Index, rcoord::V, icoord::V) where {V<:SVector} = new{typeof(index), V}(index, oidcoord(rcoord), oidcoord(icoord))
end
@inline contentnames(::Type{<:OID}) = (:index, :rcoord, :icoord)
@inline Base.hash(oid::OID, h::UInt) = hash((oid.index, Tuple(oid.rcoord)), h)
@inline Base.propertynames(::ID{OID}) = (:indexes, :rcoords, :icoords)
@inline Base.show(io::IO, oid::OID) = @printf io "OID(%s, %s, %s)" oid.index oid.rcoord oid.icoord
@inline oidcoord(vector::SVector) = vector
@generated function oidcoord(vector::SVector{N, Float}) where N
    exprs = [:((vector[$i] === -0.0) ? 0.0 : vector[$i]) for i = 1:N]
    return :(SVector($(exprs...)))
end

"""
    OID(index::Index, rcoord, icoord)
    OID(index::Index; rcoord, icoord)

Construct an operator id.
"""
@inline OID(index::Index, rcoord, icoord) = OID(index, SVector{length(rcoord)}(rcoord), SVector{length(icoord)}(icoord))
@inline OID(index::Index; rcoord, icoord) = OID(index, rcoord, icoord)

"""
    adjoint(oid::OID) -> typeof(oid)
    adjoint(oid::ID{AbstractOID}) -> ID{AbstractOID}

Get the adjoint of an operator id.
"""
@inline Base.adjoint(oid::OID) = OID(oid.index', oid.rcoord, oid.icoord)
@inline @generated Base.adjoint(oid::ID{AbstractOID}) = Expr(:call, :ID, [:(oid[$i]') for i = fieldcount(oid):-1:1]...)

"""
    oidtype(I::Type{<:SimpleInternal}, P::Type{<:Point}, ::Val)

Get the compatible oid type from the combination of the internal part and the spatial part.
"""
@inline oidtype(I::Type{<:SimpleInternal}, P::Type{<:Point}, ::Val) = OID{Index{P|>pidtype, I|>eltype}, SVector{P|>dimension, P|>dtype}}

"""
    Operator{V<:Number, I<:ID{AbstractOID}} <: OperatorProd{V, I}

Operator.
"""
struct Operator{V<:Number, I<:ID{AbstractOID}} <: OperatorProd{V, I}
    value::V
    id::I
end
@inline Operator(value::Number) = Operator(value, ())
Base.show(io::IO, opt::Operator) = @printf io "%s(%s, %s)" nameof(typeof(opt)) decimaltostr(value(opt)) id(opt)
Base.show(io::IO, ::MIME"text/plain", opt::Operator) = show(io, opt)

"""
    adjoint(opt::Operator) -> Operator

Get the adjoint of an operator.
"""
@inline Base.adjoint(opt::Operator) = rawtype(typeof(opt))(value(opt)', id(opt)')

"""
    ishermitian(opt::Operator) -> Bool

Judge whether an operator is Hermitian.
"""
@inline ishermitian(opt::Operator) = isa(value(opt), Real) && ishermitian(id(opt))

"""
    rcoord(opt::Operator) -> SVector

Get the whole rcoord of an operator.
"""
@inline function rcoord(opt::Operator{<:Number, <:ID{OID}})
    rank(opt)==1 && return id(opt)[1].rcoord
    rank(opt)==2 && return id(opt)[1].rcoord-id(opt)[2].rcoord
    error("rcoord error: not supported rank($(rank(opt))) of $(nameof(opt)).")
end

"""
    icoord(opt::Operator) -> SVector

Get the whole icoord of an operator.
"""
@inline function icoord(opt::Operator{<:Number, <:ID{OID}})
    rank(opt)==1 && return id(opt)[1].icoord
    rank(opt)==2 && return id(opt)[1].icoord-id(opt)[2].icoord
    error("icoord error: not supported rank($(rank(opt))) of $(nameof(opt)).")
end

"""
    Operators(opts::Operator...)

A set of operators.

Type alias for `OperatorSum{<:ID{AbstractOID}, <:Operator}`.
"""
const Operators{I<:ID{AbstractOID}, O<:Operator} = OperatorSum{I, O}
@inline Operators(opts::Operator...) = OperatorSum(opts...)

"""
    adjoint(opts::Operators) -> Operators

Get the adjoint of a set of operators.
"""
function Base.adjoint(opts::Operators)
    result = Operators{opts|>keytype, opts|>valtype}()
    for opt in values(opts)
        nopt = opt |> adjoint
        result[id(nopt)] = nopt
    end
    return result
end

"""
    summary(io::IO, opts::Operators)

Print a brief description of a set of operators to an io.
"""
Base.summary(io::IO, opts::Operators) = @printf io "Operators{%s}" valtype(opts)

"""
    ishermitian(opts::Operators) -> Bool

Judge whether a set of operators as a whole is Hermitian.
"""
@inline ishermitian(opts::Operators) = opts == opts'

"""
    IIDSpace{I<:IID, V<:Internal, Kind} <: CartesianVectorSpace{IID}

The space expanded by a "labeled" iid.
"""
struct IIDSpace{I<:IID, V<:Internal, Kind} <: CartesianVectorSpace{IID}
    iid::I
    internal::V
    IIDSpace(iid::IID, internal::Internal, ::Val{Kind}=Val(:info)) where Kind = new{typeof(iid), typeof(internal), Kind}(iid, internal)
end
Base.eltype(iidspace::IIDSpace) = eltype(typeof(iidspace))
Base.eltype(::Type{<:IIDSpace{<:IID, V}}) where {V<:Internal} = eltype(V)
kind(iidspace::IIDSpace) = kind(typeof(iidspace))
kind(::Type{<:IIDSpace{<:IID, <:Internal, Kind}}) where Kind = Kind
@generated function shape(iidspace::IIDSpace{I, V}) where {I<:CompositeIID, V<:CompositeInternal}
    @assert rank(I)==rank(V) "shape error: mismatched composite iid and composite internal space."
    Kind = Val(kind(iidspace))
    Expr(:tuple, [:(shape(IIDSpace(iidspace.iid[$i], iidspace.internal.contents[$i], $Kind))...) for i = 1:rank(I)]...)
end
ndimshape(::Type{<:IIDSpace{<:IID, V}}) where {V<:Internal} = ndimshape(V)
Base.CartesianIndex(iid::IID, iidspace::IIDSpace) = CartesianIndex(iid, iidspace.internal)
Base.getindex(iidspace::IIDSpace, index::CartesianIndex) = rawtype(eltype(iidspace))(index, iidspace.internal)

"""
    expand(iids::NTuple{N, IID}, internals::NTuple{N, Internal}) where N -> IIDSpace

Get the space expanded by a set of "labeled" iids.
"""
@inline expand(iids::NTuple{N, IID}, internals::NTuple{N, Internal}) where N = IIDSpace(CompositeIID(iids), CompositeInternal(internals))

"""
    Hilbert{I<:Internal, P<:AbstractPID, M<:Function} <: CompositeDict{P, I}

Hilbert space at a lattice.
"""
struct Hilbert{I<:Internal, P<:AbstractPID, M<:Function} <: CompositeDict{P, I}
    map::M
    contents::Dict{P, I}
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
    map = pid -> contents[pid]
    return Hilbert(map, contents)
end

"""
    Hilbert(map::Function, pids::AbstractVector{<:AbstractPID})

Construct a Hilbert space from a function and a set of point ids.

Here, `map` maps a `AbstractPID` to an `Internal`.
"""
@inline function Hilbert(map::Function, pids::AbstractVector{<:AbstractPID})
    I = commontype(map, Tuple{Vararg{Any}}, Internal)
    return Hilbert{I}(map, pids)
end

"""
    Hilbert{I}(map::Function, pids::AbstractVector{<:AbstractPID}) where {I<:Internal}

Construct a Hilbert space from a function and a set of point ids.

Here, `map` maps a `AbstractPID` to an `Internal`.
"""
function Hilbert{I}(map::Function, pids::AbstractVector{<:AbstractPID}) where {I<:Internal}
    contents = Dict{pids|>eltype, I}()
    for pid in pids
        contents[pid] = map(pid)
    end
    return Hilbert(map, contents)
end

"""
    reset!(hilbert::Hilbert, pids) -> Hilbert

Reset the Hilbert space with new pids.
"""
function reset!(hilbert::Hilbert, pids)
    empty!(hilbert)
    for pid in pids
        hilbert[pid] = hilbert.map(pid)
    end
    hilbert
end

@inline diagonal(xs...) = length(xs)<2 ? true : all(map(==(xs[1]), xs))
@inline noconstrain(_...) = true
const wildcard = Symbol("*")
"""
    Subscript{P<:Tuple, C<:Function} <: CompositeTuple{P}

The subscript representative of a certain internal degree of freedom.
"""
struct Subscript{P<:Tuple, C<:Function} <: CompositeTuple{P}
    pattern::P
    rep::String
    constraint::C
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
    return subscript.constraint(values...)
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
    Subscripts{T<:Tuple{Vararg{NamedContainer{Subscript}}}} <: CompositeTuple{T}

The complete set of subscripts of the internal degrees of freedom.
"""
struct Subscripts{T<:Tuple{Vararg{NamedContainer{Subscript}}}} <: CompositeTuple{T}
    contents::T
end
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
    return Subscripts(contents)
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
@inline Base.:*(subscripts₁::Subscripts, subscripts₂::Subscripts) = Subscripts((subscripts₁.contents..., subscripts₂.contents...))

"""
    idtype(subscripts::Subscripts)
    idtype(::Type{<:Subscripts{T}}) where {T<:Tuple{Vararg{NamedContainer{Subscript}}}}

Get the type of the id of the subscripts.
"""
@inline idtype(subscripts::Subscripts) = idtype(typeof(subscripts))
@generated function idtype(::Type{<:Subscripts{T}}) where {T<:Tuple{Vararg{NamedContainer{Subscript}}}}
    exprs = []
    for i = 1:fieldcount(T)
        push!(exprs, Pair{UnitRange{Int}, NTuple{fieldcount(fieldtype(T, 1)), String}})
    end
    return Expr(:curly, :SubscriptsID, Expr(:curly, :Tuple, exprs...))
end

"""
    SubscriptsID{T<:Tuple{Vararg{Pair{UnitRange{Int}, <:Tuple{Vararg{String}}}}}} <: SingularID

The id of the subscripts.
"""
struct SubscriptsID{T<:Tuple{Vararg{Pair{UnitRange{Int}, <:Tuple{Vararg{String}}}}}} <: SingularID
    rep::T
end

"""
    SubscriptsID(subscripts::Subscripts)

Construct the id of the subscripts.
"""
@generated function SubscriptsID(subscripts::Subscripts)
    exprs, count = [], 1
    for i = 1:length(subscripts)
        rep = []
        for field in fieldnames(fieldtype(fieldtype(subscripts, :contents), i))
            field = QuoteNode(field)
            push!(rep, :(getfield(getfield(subscripts[$i], $field), :rep)))
        end
        push!(exprs, Expr(:call, :(=>), count:(count+rank(subscripts, i)-1), Expr(:tuple, rep...)))
        count += rank(subscripts, i)
    end
    return Expr(:call, :SubscriptsID, Expr(:tuple, exprs...))
end

"""
    AbstractCoupling{V, I<:ID{SingularID}} <: OperatorProd{V, I}

The abstract coupling intra/inter internal degrees of freedom at different lattice points.
"""
abstract type AbstractCoupling{V, I<:ID{SingularID}} <: OperatorProd{V, I} end
ID{SimpleIID}(coupling::AbstractCoupling) = id(coupling)
Subscripts(coupling::AbstractCoupling) = Subscripts()

"""
    couplingcenters(coupling::AbstractCoupling, bond::AbstractBond, info::Val) -> NTuple{rank(coupling), Int}

Get the acting centers of the coupling on a bond.
"""
@inline couplingcenters(coupling::AbstractCoupling, point::Point, ::Val) = ntuple(i->1, Val(rank(coupling)))

"""
    couplingpoints(coupling::AbstractCoupling, bond::AbstractBond, info::Val) -> NTuple{rank(coupling), eltype(bond)}

Get the points where each order of the coupling acts on.
"""
@inline function couplingpoints(coupling::AbstractCoupling, bond::AbstractBond, info::Val)
    centers = couplingcenters(coupling, bond, info)
    return NTuple{rank(coupling), eltype(bond)}(bond[centers[i]] for i = 1:rank(coupling))
end

"""
    couplinginternals(coupling::AbstractCoupling, bond::AbstractBond, hilbert::Hilbert, info::Val) -> NTuple{rank(coupling), SimpleInternal}

Get the internal spaces where each order of the coupling acts on.
"""
@inline @generated function couplinginternals(coupling::AbstractCoupling, bond::AbstractBond, hilbert::Hilbert, info::Val)
    exprs = []
    for i = 1:rank(coupling)
        push!(exprs, :(filter(ID{SimpleIID}(coupling)[$i], hilbert[bond[centers[$i]].pid])))
    end
    return Expr(:block, :(centers = couplingcenters(coupling, bond, info)), Expr(:tuple, exprs...))
end

"""
    expand(coupling::AbstractCoupling, bond::AbstractBond, hilbert::Hilbert, info::Val)

Expand a coupling with the given bond and Hilbert space.
"""
function expand(coupling::AbstractCoupling, bond::AbstractBond, hilbert::Hilbert, info::Val)
    points = couplingpoints(coupling, bond, info)
    internals = couplinginternals(coupling, bond, hilbert, info)
    @assert rank(coupling)==length(points)==length(internals) "expand error: mismatched rank."
    return CExpand(value(coupling), points, IIDSpace(CompositeIID(ID{SimpleIID}(coupling)), CompositeInternal(internals), info), Subscripts(coupling))
end
struct CExpand{V, N, P<:AbstractPID, SV<:SVector, S<:IIDSpace, C<:Subscripts}
    value::V
    pids::NTuple{N, P}
    rcoords::NTuple{N, SV}
    icoords::NTuple{N, SV}
    iidspace::S
    subscripts::C
end
function CExpand(value, points::NTuple{N, P}, iidspace::IIDSpace, subscripts::Subscripts) where {N, P<:Point}
    pids = NTuple{N, pidtype(P)}(points[i].pid for i = 1:N)
    rcoords = NTuple{N, SVector{dimension(P), dtype(P)}}(points[i].rcoord for i = 1:N)
    icoords = NTuple{N, SVector{dimension(P), dtype(P)}}(points[i].icoord for i = 1:N)
    return CExpand(value, pids, rcoords, icoords, iidspace, subscripts)
end
@inline Base.eltype(ex::CExpand) = eltype(typeof(ex))
@inline @generated function Base.eltype(::Type{<:CExpand{V, N, P, SV, S}}) where {V, N, P<:AbstractPID, SV<:SVector, S<:IIDSpace}
    IIDS = fieldtypes(fieldtype(eltype(S), :contents))
    return Tuple{V, Tuple{[OID{Index{P, IIDS[i]}, SV} for i = 1:N]...}}
end
@inline Base.IteratorSize(::Type{<:CExpand}) = Base.SizeUnknown()
function Base.iterate(ex::CExpand, state=iterate(ex.iidspace))
    result = nothing
    while !isnothing(state)
        ciid, state = state
        if match(ex.subscripts, ciid)
            result = (ex.value, ID(OID, ID(Index, ex.pids, ciid.contents), ex.rcoords, ex.icoords)), iterate(ex.iidspace, state)
            break
        else
            state = iterate(ex.iidspace, state)
        end
    end
    return result
end

"""
    Coupling{V, I<:ID{SimpleIID}, C<:Subscripts, CI<:SubscriptsID} <: AbstractCoupling{V, Tuple{CompositeIID{I}, CI}}

The coupling intra/inter internal degrees of freedom at different lattice points.
"""
struct Coupling{V, I<:ID{SimpleIID}, C<:Subscripts, CI<:SubscriptsID} <: AbstractCoupling{V, Tuple{CompositeIID{I}, CI}}
    value::V
    cid::I
    subscripts::C
    function Coupling(value::Number, cid::ID{SimpleIID}, subscripts::Subscripts=Subscripts())
        new{typeof(value), typeof(cid), typeof(subscripts), idtype(subscripts)}(value, cid, subscripts)
    end
end
@inline parameternames(::Type{<:Coupling}) = (:value, :cid, :subscripts, :subscriptsid)
@inline isparameterbound(::Type{<:Coupling}, ::Val{:cid}, ::Type{I}) where {I<:ID{SimpleIID}} = !isconcretetype(I)
@inline isparameterbound(::Type{<:Coupling}, ::Val{:subscripts}, ::Type{C}) where {C<:Subscripts} = !isconcretetype(C)
@inline isparameterbound(::Type{<:Coupling}, ::Val{:subscriptsid}, ::Type{CI}) where {CI<:SubscriptsID} = !isconcretetype(CI)
@inline contentnames(::Type{<:Coupling}) = (:value, :id, :subscripts)
@inline getcontent(coupling::Coupling, ::Val{:id}) = ID(CompositeIID(coupling.cid), SubscriptsID(coupling.subscripts))
@inline rank(::Type{<:Coupling{V, I} where V}) where {I<:ID{SimpleIID}} = fieldcount(I)
@inline Coupling(value::Number, id::Tuple{CompositeIID, SubscriptsID}, subscripts::Subscripts) = Coupling(value, first(id).contents, subscripts)
@inline ID{SimpleIID}(coupling::Coupling) = coupling.cid
@inline Subscripts(coupling::Coupling) = coupling.subscripts

"""
    *(cp₁::Coupling, cp₂::Coupling) -> Coupling

Get the multiplication between two couplings.
"""
@inline Base.:*(cp₁::Coupling, cp₂::Coupling) = Coupling(cp₁.value*cp₂.value, cp₁.cid*cp₂.cid, cp₁.subscripts*cp₂.subscripts)

"""
    Couplings(cps::AbstractCoupling...)

A pack of couplings intra/inter internal degrees of freedom at different lattice points.

Alias for `OperatorSum{<:ID{SingularID}, <:AbstractCoupling}`.
"""
const Couplings{I<:ID{SingularID}, C<:AbstractCoupling} = OperatorSum{I, C}
@inline Couplings(cps::AbstractCoupling...) = OperatorSum(cps...)
@inline Couplings(cps::Couplings) = cps

"""
    @couplings cps -> Couplings

Convert an expression/literal to a set of couplings.
"""
macro couplings(cps) :(Couplings($(esc(cps)))) end

"""
    Metric <: Function

The rules for measuring a concrete oid so that oids can be compared.

As a function, every instance should accept only one positional argument, i.e. the concrete oid to be measured.
"""
abstract type Metric <: Function end
@inline Base.:(==)(m₁::T, m₂::T) where {T<:Metric} = ==(efficientoperations, m₁, m₂)
@inline Base.isequal(m₁::T, m₂::T) where {T<:Metric} = isequal(efficientoperations, m₁, m₂)
@inline (M::Type{<:Metric})(::Type{I}) where {I<:CompositeOID} = M(indextype(I))
@inline (metric::Metric)(oid::CompositeOID) = metric(getcontent(oid, :index))
@inline Base.valtype(::Type{M}, ::Type{I}) where {M<:Metric, I<:CompositeOID} = valtype(M, indextype(I))
@inline (M::Type{<:Metric})(::Type{H}) where {H<:Hilbert} = M(Index{H|>keytype, H|>valtype|>eltype})

"""
    OIDToTuple{Fields} <: Metric

A rule that converts an oid to a tuple by iterating over a set of selected fields in a specific order.
"""
struct OIDToTuple{Fields} <: Metric
    OIDToTuple(fields::Tuple{Vararg{Symbol}}) = new{fields}()
end
@inline OIDToTuple(fields::Symbol...) = OIDToTuple(fields)

"""
    keys(::OIDToTuple{Fields}) where Fields -> Fields
    keys(::Type{<:OIDToTuple{Fields}}) where Fields -> Fields

Get the names of the selected fields.
"""
@inline Base.keys(::OIDToTuple{Fields}) where Fields = Fields
@inline Base.keys(::Type{<:OIDToTuple{Fields}}) where Fields = Fields

"""
    filter(f::Function, oidtotuple::OIDToTuple) -> OIDToTuple

Filter the selected fields.
"""
@inline Base.filter(f::Function, oidtotuple::OIDToTuple) = OIDToTuple(Tuple(field for field in keys(oidtotuple) if f(field)))

"""
    OIDToTuple(::Type{I}) where {I<:Index}

Construct the conversion rule from the information of subtypes of `AbstractOID`.
"""
@inline @generated OIDToTuple(::Type{I}) where {I<:Index} = OIDToTuple(fieldnames(pidtype(I))..., (fieldnames(iidtype(I)))...)

"""
    valtype(::Type{<:OIDToTuple}, ::Type{<:Index})

Get the valtype of applying an `OIDToTuple` rule to a subtype of `AbstractOID`.
"""
@inline @generated function Base.valtype(::Type{M}, ::Type{I}) where {M<:OIDToTuple, I<:Index}
    types = []
    for field in keys(M)
        if hasfield(pidtype(I), field)
            push!(types, fieldtype(pidtype(I), field))
        elseif hasfield(iidtype(I), field)
            push!(types, fieldtype(iidtype(I), field))
        end
    end
    return  Expr(:curly, :Tuple, types...)
end

"""
    (oidtotuple::OIDToTuple)(index::Index) -> Tuple

Convert a concrete oid to a tuple.
"""
@inline @generated function (oidtotuple::OIDToTuple)(index::Index)
    exprs = []
    for name in keys(oidtotuple)
        field = QuoteNode(name)
        if hasfield(pidtype(index), name)
            push!(exprs, :(getfield(index.pid, $field)))
        elseif hasfield(iidtype(index), name)
            push!(exprs, :(getfield(index.iid, $field)))
        end
    end
    return Expr(:tuple, exprs...)
end

"""
    Table{I, B<:Metric} <: CompositeDict{I, Int}

The table of oid-sequence pairs.
"""
struct Table{I, B<:Metric} <: CompositeDict{I, Int}
    by::B
    contents::Dict{I, Int}
end
@inline contentnames(::Type{<:Table}) = (:by, :contents)
@inline Table{I}(by::Metric) where {I<:AbstractOID} = Table(by, Dict{valtype(typeof(by), I), Int}())
@inline vec2dict(vs::AbstractVector) = Dict{eltype(vs), Int}(v=>i for (i, v) in enumerate(vs))

"""
    getindex(table::Table, oid::AbstractOID) -> Int

Inquiry the sequence of an oid.
"""
@inline Base.getindex(table::Table, oid::AbstractOID) = table[table.by(oid)]

"""
    haskey(table::Table, oid::AbstractOID) -> Bool
    haskey(table::Table, id::ID{AbstractOID}) -> Tuple{Vararg{Bool}}

Judge whether a single oid or a set of oids have been assigned with sequences in table.
"""
@inline Base.haskey(table::Table, oid::AbstractOID) = haskey(table, table.by(oid))
@inline @generated function Base.haskey(table::Table, id::ID{AbstractOID})
    exprs = []
    for i = 1:fieldcount(id)
        push!(exprs, :(haskey(table, id[$i])))
    end
    return Expr(:tuple, exprs...)
end

"""
    Table(oids::AbstractVector{<:AbstractOID}, by::Metric=OIDToTuple(eltype(oids)))

Convert a set of concrete oids to the corresponding table of oid-sequence pairs.

The input oids are measured by the input `by` function with the duplicates removed. The resulting unique values are sorted, which determines the sequence of the input `oids`. Note that two oids have the same sequence if their converted values are equal to each other.
"""
@inline Table(oids::AbstractVector{<:AbstractOID}, by::Metric=OIDToTuple(eltype(oids))) = Table(by, [by(oid) for oid in oids]|>unique!|>sort!|>vec2dict)

"""
    Table(hilbert::Hilbert, by::Metric=OIDToTuple(typeof(hilbert))) -> Table

Get the oid-sequence table of a Hilbert space.
"""
function Table(hilbert::Hilbert, by::Metric=OIDToTuple(typeof(hilbert)))
    result = Index{hilbert|>keytype, hilbert|>valtype|>eltype}[]
    for (pid, internal) in hilbert
        for iid in internal
            push!(result, (result|>eltype)(pid, iid))
        end
    end
    return Table(result, by)
end

"""
    union(tables::Table...) -> Table

Unite several oid-sequence tables.
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
    reset!(table::Table, oids::AbstractVector{<:AbstractOID}) -> Table

Reset a table by a new set of oids.
"""
function reset!(table::Table, oids::AbstractVector{<:AbstractOID})
    empty!(table)
    for (i, id) in enumerate([table.by(oid) for oid in oids]|>unique!|>sort!)
        table[id] = i
    end
    return table
end

"""
    reset!(table::Table, hilbert::Hilbert) -> Table

Reset a table by a Hilbert space.
"""
function reset!(table::Table, hilbert::Hilbert)
    indices = Index{hilbert|>keytype, hilbert|>valtype|>eltype}[]
    for (pid, internal) in hilbert
        for iid in internal
            push!(indices, (indices|>eltype)(pid, iid))
        end
    end
    reset!(table, indices)
end

"""
    sequence(opt::Operator, table::AbstractDict) -> NTuple{rank(opt), Int}

Get the sequence of the oids of an operator according to a table.
"""
@inline @generated function sequence(opt::Operator, table::AbstractDict{<:Any,Int})
    return Expr(:tuple, [:(table[id(opt)[$i]]) for i = 1:rank(opt)]...)
end

"""
    TermFunction <: Function

Abstract type for concrete term functions.
"""
abstract type TermFunction <: Function end
@inline Base.:(==)(tf1::TermFunction, tf2::TermFunction) = ==(efficientoperations, tf1, tf2)
@inline Base.isequal(tf1::TermFunction, tf2::TermFunction) = isequal(efficientoperations, tf1, tf2)

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
Base.valtype(termcouplings::TermCouplings) = valtype(typeof(termcouplings))
Base.valtype(::Type{<:TermCouplings{C}}) where {C<:Couplings} = C
Base.valtype(::Type{<:TermCouplings{<:Function, C}}) where {C<:Couplings} = C
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
    modulate::M
    factor::V
    function Term{K, I}(value, bondkind, couplings::TermCouplings, amplitude::TermAmplitude, modulate::TermModulate, factor) where {K, I}
        @assert isa(K, Symbol) "Term error: kind must be a Symbol."
        @assert isa(I, Symbol) "Term error: id must be a Symbol."
        V, B, C, A, M = typeof(value), typeof(bondkind), typeof(couplings), typeof(amplitude), typeof(modulate)
        new{K, I, V, B, C, A, M}(value, bondkind, couplings, amplitude, modulate, factor)
    end
end
@inline Base.:(==)(term1::Term, term2::Term) = ==(efficientoperations, term1, term2)
@inline Base.isequal(term1::Term, term2::Term) = isequal(efficientoperations, term1, term2)
function Base.show(io::IO, term::Term)
    @printf io "%s{%s}(id=%s, value=%s, bondkind=%s, factor=%s)" kind(term) rank(term) id(term) decimaltostr(term.value) term.bondkind decimaltostr(term.factor)
end

"""
    Term{K}(id::Symbol, value, bondkind;
        couplings::Union{Function, Couplings},
        amplitude::Union{Function, Nothing}=nothing,
        modulate::Union{Function, Bool}=false
        ) where K

Construct a term.
"""
@inline function Term{K}(id::Symbol, value, bondkind;
        couplings::Union{Function, Couplings},
        amplitude::Union{Function, Nothing}=nothing,
        modulate::Union{Function, Bool}=false
        ) where K
    Term{K, id}(value, bondkind, TermCouplings(couplings), TermAmplitude(amplitude), TermModulate(id, modulate), 1)
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
@inline rank(::Type{<:Term{K, I, V, B, C} where {K, I, V, B}}) where {C<:TermCouplings} = rank(valtype(valtype(C)))

"""
    abbr(term::Term) -> Symbol
    abbr(::Type{<:Term}) -> Symbol

Get the abbreviation of the kind of a term.
"""
@inline abbr(term::Term) = abbr(typeof(term))
@inline abbr(::Type{<:Term}) = :tm

"""
    ismodulatable(term::Term) -> Bool
    ismodulatable(::Type{<:Term}) -> Bool

Judge whether a term could be modulated by its modulate function.
"""
@inline ismodulatable(term::Term) = ismodulatable(typeof(term))
@inline ismodulatable(::Type{<:Term{K, I, V, B, <:TermCouplings, <:TermAmplitude, M} where {K, I, V, B}}) where M = ismodulatable(M)

"""
    ishermitian(term::Term) -> Bool
    ishermitian(::Type{<:Term}) -> Bool
"""
@inline ishermitian(term::Term) = ishermitian(typeof(term))
@inline ishermitian(::Type{<:Term}) = error("ishermitian error: not implemented.")

"""
    repr(term::Term, bond::AbstractBond, hilbert::Hilbert) -> String

Get the repr representation of a term on a bond with a given Hilbert space.
"""
function Base.repr(term::Term, bond::AbstractBond, hilbert::Hilbert)
    cache = String[]
    if term.bondkind == bond|>kind
        value = term.value * term.amplitude(bond) * term.factor
        if abs(value) ≠ 0
            for coupling in values(term.couplings(bond))
                isnothing(iterate(expand(coupling, bond, hilbert, term|>kind|>Val))) || push!(cache, @sprintf "%s: %s" abbr(term) repr(value*coupling))
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
    +(term::Term) -> Term
    -(term::Term) -> Term
    *(term::Term, factor) -> Term
    *(factor, term::Term) -> Term
    /(term::Term, factor) -> Term

Allowed arithmetic operations for a term.
"""
@inline Base.:+(term::Term) = term
@inline Base.:-(term::Term) = term * (-1)
@inline Base.:*(term::Term, factor) = factor * term
@inline Base.:*(factor, term::Term) = replace(term, factor=factor*term.factor)
@inline Base.:/(term::Term, factor) = term * (one(term.value)/factor)

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
    otype(::Type{T}, ::Type{H}, ::Type{B}) where {T<:Term, H<:Hilbert, B<:AbstractBond}

Get the compatible operator type from the type of a term, a Hilbert space and a bond.
"""
@inline @generated function otype(::Type{T}, ::Type{H}, ::Type{B}) where {T<:Term, H<:Hilbert, B<:AbstractBond}
    C = valtype(valtype(fieldtype(T, :couplings)))
    @assert C<:Coupling "otype error: not supported."
    exprs = []
    for i = 1:rank(C)
        I = fieldtype(parametertype(C, :cid), i)
        push!(exprs, :(oidtype(filter($I, valtype(H)), eltype(B), Val(kind(T)))))
    end
    return Expr(:curly, Operator, valtype(T), Expr(:curly, Tuple, exprs...))
end

"""
    expand!(operators::Operators, term::Term, bond::AbstractBond, hilbert::Hilbert; half::Bool=false, table::Union{Nothing, Table}=nothing) -> Operators
    expand!(operators::Operators, term::Term, bonds::Bonds, hilbert::Hilbert; half::Bool=false, table::Union{Nothing, Table}=nothing) -> Operators

Expand the operators of a term on a bond/set-of-bonds with a given Hilbert space.

The `half` parameter determines the behavior of generating operators, which falls into the following two categories
* `true`: "Hermitian half" of the generated operators
* `false`: "Hermitian whole" of the generated operators
"""
function expand!(operators::Operators, term::Term, bond::AbstractBond, hilbert::Hilbert; half::Bool=false, table::Union{Nothing, Table}=nothing)
    if term.bondkind == bond|>kind
        value = term.value * term.amplitude(bond) * term.factor
        if abs(value) ≠ 0
            hermitian = ishermitian(term)
            optype = otype(term|>typeof, hilbert|>typeof, bond|>typeof)
            record = (isnothing(hermitian) && length(operators)>0) ? Set{optype|>idtype}() : nothing
            for coupling in values(term.couplings(bond))
                for (coeff, id) in expand(coupling, bond, hilbert, term|>kind|>Val)
                    isapprox(coeff, 0.0, atol=atol, rtol=rtol) && continue
                    !isnothing(table) && !all(haskey(table, id)) && continue
                    if hermitian == true
                        add!(operators, rawtype(optype)(convert(optype|>valtype, value*coeff/(half ? 2 : 1)), id))
                    elseif hermitian == false
                        opt = rawtype(optype)(convert(optype|>valtype, value*coeff), id)
                        add!(operators, opt)
                        half || add!(operators, opt')
                    else
                        if !(isnothing(record) ? haskey(operators, id') : id'∈record)
                            isnothing(record) || push!(record, id)
                            ovalue = valtype(optype)(value*coeff/termfactor(id, term|>kind|>Val))
                            opt = rawtype(optype)(ovalue, id)
                            add!(operators, opt)
                            half || add!(operators, opt')
                        end
                    end
                end
            end
        end
    end
    return operators
end
@generated function expand!(operators::Operators, term::Term, bonds::Bonds, hilbert::Hilbert; half::Bool=false, table::Union{Nothing, Table}=nothing)
    exprs = []
    for i = 1:rank(bonds)
        push!(exprs, :(for bond in bonds.bonds[$i] expand!(operators, term, bond, hilbert; half=half, table=table) end))
    end
    push!(exprs, :(return operators))
    return Expr(:block, exprs...)
end
@inline termfactor(id::ID{AbstractOID}, ::Val) = ishermitian(id) ? 2 : 1

"""
    expand(term::Term, bond::AbstractBond, hilbert::Hilbert; half::Bool=false, table::Union{Nothing, Table}=nothing) -> Operators
    expand(term::Term, bonds::Bonds, hilbert::Hilbert; half::Bool=false, table::Union{Nothing, Table}=nothing) -> Operators

Expand the operators of a term on a bond/set-of-bonds with a given Hilbert space.
"""
@inline function expand(term::Term, bond::AbstractBond, hilbert::Hilbert; half::Bool=false, table::Union{Nothing, Table}=nothing)
    optype = otype(term|>typeof, hilbert|>typeof, bond|>typeof)
    expand!(Operators{idtype(optype), optype}(), term, bond, hilbert; half=half, table=table)
end
@inline function expand(term::Term, bonds::Bonds, hilbert::Hilbert; half::Bool=false, table::Union{Nothing, Table}=nothing)
    optype = otype(term|>typeof, hilbert|>typeof, bonds|>eltype)
    expand!(Operators{idtype(optype), optype}(), term, bonds, hilbert; half=half, table=table)
end

"""
    LaTeX{SP, SB}(body, spdelimiter::String=", ", sbdelimiter::String=", "; options...) where {SP, SB}

LaTeX string representation.
"""
struct LaTeX{SP, SB, B, O}
    body::B
    spdelimiter::String
    sbdelimiter::String
    options::O
    function LaTeX{SP, SB}(body, spdelimiter::String=", ", sbdelimiter::String=", "; options...) where {SP, SB}
        @assert isa(SP, Tuple{Vararg{Symbol}}) && isa(SB, Tuple{Vararg{Symbol}}) "LaTeX error: SP and SB must be tuple of symbols."
        new{SP, SB, typeof(body), typeof(options)}(body, spdelimiter, sbdelimiter, options)
    end
end
@inline superscript(::Type{<:LaTeX{SP}}) where SP = SP
@inline subscript(::Type{<:LaTeX{SP, SB} where SP}) where SB = SB

"""
    latexname(T::Type{<:AbstractOID}) -> Symbol

Get the name of a type of `AbstractOID` in the latex format lookups.
"""
@inline latexname(T::Type{<:AbstractOID}) = nameof(T)

const latexformats = Dict{Symbol, LaTeX}()
"""
    latexformat(T::Type{<:AbstractOID}) -> LaTeX
    latexformat(T::Type{<:AbstractOID}, l::LaTeX) -> LaTeX

Get/Set the LaTeX format for a subtype of `AbstractOID`.
"""
@inline latexformat(T::Type{<:AbstractOID}) = latexformats[latexname(T)]
@inline latexformat(T::Type{<:AbstractOID}, l::LaTeX) = latexformats[latexname(T)] = l

"""
    script(::Val{:BD}, oid::AbstractOID, l::LaTeX) -> Any
    script(::Val{:SP}, oid::AbstractOID, l::LaTeX) -> Tuple
    script(::Val{:SB}, oid::AbstractOID, l::LaTeX) -> Tuple

Get the body/superscript/subscript of the LaTeX string representation of an oid.
"""
@inline script(::Val{:BD}, oid::AbstractOID, l::LaTeX) = l.body
@inline @generated script(::Val{:SP}, oid::AbstractOID, l::LaTeX) = Expr(:tuple, [:(script(Val($sup), oid; l.options...)) for sup in QuoteNode.(l|>superscript)]...)
@inline @generated script(::Val{:SB}, oid::AbstractOID, l::LaTeX) = Expr(:tuple, [:(script(Val($sub), oid; l.options...)) for sub in QuoteNode.(l|>subscript)]...)

"""
    script(::Val{:rcoord}, oid::OID; kwargs...) -> String
    script(::Val{:icoord}, oid::OID; kwargs...) -> String

Get the `:rcoord/:icoord` script of an oid.
"""
@inline script(::Val{:rcoord}, oid::OID; kwargs...) = @sprintf "[%s]" join(valuetolatextext.(oid.rcoord), ", ")
@inline script(::Val{:icoord}, oid::OID; kwargs...) = @sprintf "[%s]" join(valuetolatextext.(oid.icoord), ", ")

"""
    script(::Val{:integralicoord}, oid::OID; vectors, kwargs...)

Get the integral script of the icoord of an oid.
"""
function script(::Val{:integralicoord}, oid::OID; vectors, kwargs...)
    rcoeff = decompose(oid.icoord, vectors...)
    icoeff = Int.(round.(rcoeff))
    @assert isapprox(efficientoperations, rcoeff, icoeff) "script error: mismatched icoord of oid and input vectors."
    return @sprintf "[%s]" join(icoeff, ", ")
end

"""
    script(::Val{attr}, oid::OID; kwargs...) where attr

Get the `attr` script of an oid, which is contained in its index.
"""
@inline script(::Val{attr}, oid::OID; kwargs...) where attr = script(Val(attr), oid.index; kwargs...)

"""
    repr(oid::AbstractOID) -> String
    repr(oid::AbstractOID, l::LaTeX) -> String

LaTeX string representation of an oid.
"""
@inline Base.repr(oid::AbstractOID) = repr(oid, latexformat(typeof(oid)))
@inline function Base.repr(oid::AbstractOID, l::LaTeX)
    @sprintf "%s^{%s}_{%s}" script(Val(:BD), oid, l) join(script(Val(:SP), oid, l), l.spdelimiter) join(script(Val(:SB), oid, l), l.sbdelimiter)
end

"""
    repr(opt::Operator) -> String

Get the string representation of an operator in the LaTeX format.
"""
function Base.repr(opt::Operator)
    rank(opt)==0 && return replace(valuetolatextext(value(opt)), " "=>"")
    poses = Int[]
    push!(poses, 1)
    for i = 2:rank(opt)
        id(opt)[i]≠id(opt)[i-1] && push!(poses, i)
    end
    push!(poses, rank(opt)+1)
    result = valuetostr(value(opt))
    for i = 1:length(poses)-1
        order = poses[i+1] - poses[i]
        if order == 1
            result = @sprintf "%s%s" result repr(id(opt)[poses[i]])
        else
            result = @sprintf "%s(%s)^%s" result repr(id(opt)[poses[i]]) order
        end
    end
    return result
end
function valuetostr(v)
    v==+1 && return ""
    v==-1 && return "-"
    result = valuetolatextext(v)
    if occursin(" ", result) || (isa(v, Complex) && real(v)≠0 && imag(v)≠0)
        bra, ket = occursin("(", result) ? ("[", "]") : ("(", ")")
        result = @sprintf "%s%s%s" bra replace(result, " "=>"") ket
    end
    return result
end
@inline valuetolatextext(value::Union{Real, Complex}) = decimaltostr(value)
function valuetolatextext(value)
    io = IOBuffer()
    if showable(MIME"text/latex"(), value)
        show(IOContext(io, :limit=>false), MIME"text/latex"(), value)
    else
        show(IOContext(io, :limit=>false), MIME"text/plain"(), value)
    end
    return replace(replace(replace(replace(String(take!(io)), "\\begin{equation*}" => ""), "\\end{equation*}" => ""), "\$" => ""), "\n" => "")
end

"""
    repr(opts::Operators) -> String

Get the string representation of a set of operators in the LaTeX format.
"""
function Base.repr(opts::Operators)
    result = String[]
    for (i, opt) in enumerate(values(opts))
        rep = repr(opt)
        i>1 && rep[1]≠'-' && push!(result, "+")
        push!(result, rep)
    end
    return join(result, "")
end

"""
    show(io::IO, ::MIME"text/latex", opt::Operator)

Show an operator.
"""
Base.show(io::IO, ::MIME"text/latex", opt::Operator) = show(io, MIME"text/latex"(), latexstring(repr(opt)))

"""
    show(io::IO, ::MIME"text/latex", opts::Operators)

Show LaTeX formed operators.
"""
Base.show(io::IO, ::MIME"text/latex", opts::Operators) = show(io, MIME"text/latex"(), latexstring(repr(opts)))

"""
    Boundary{Names}(values::AbstractVector{<:Number}, vectors::AbstractVector{<:AbstractVector{<:Number}}) where Names

Boundary twist of operators.
"""
struct Boundary{Names, D<:Number, V<:AbstractVector} <: Function
    values::Vector{D}
    vectors::Vector{V}
    function Boundary{Names}(values::AbstractVector{<:Number}, vectors::AbstractVector{<:AbstractVector{<:Number}}) where Names
        @assert length(Names)==length(values)==length(vectors) "Boundary error: mismatched names, values and vectors."
        datatype = promote_type(eltype(values), Float)
        new{Names, datatype, eltype(vectors)}(convert(Vector{datatype}, values), vectors)
    end
end
@inline Base.:(==)(bound1::Boundary, bound2::Boundary) = ==(efficientoperations, bound1, bound2)
@inline Base.isequal(bound1::Boundary, bound2::Boundary) = isequal(efficientoperations, bound1, bound2)

"""
    keys(bound::Boundary) -> Tuple{Vararg{Symbol}}
    keys(::Type{<:Boundary{Names}}) where Names -> Names

Get the names of the boundary parameters.
"""
@inline Base.keys(bound::Boundary) = keys(typeof(bound))
@inline Base.keys(::Type{<:Boundary{Names}}) where Names = Names

"""
    (bound::Boundary)(operator::Operator) -> Operator

Get the boundary twisted operator.
"""
@inline (bound::Boundary)(operator::Operator) = twist(operator, bound.vectors, bound.values)

"""
    twist(operator::Operator, vectors::AbstractVector{<:AbstractVector{<:Number}}, values::AbstractVector{<:Number}) -> Operator

Twist an operator.
"""
@inline function twist(operator::Operator, vectors::AbstractVector{<:AbstractVector{<:Number}}, values::AbstractVector{<:Number})
    return replace(operator, operator.value*exp(1im*angle(operator.id, vectors, values)))
end

"""
    angle(bound::Boundary, operator::Operator) -> Number

Get the boundary twist phase of an operator.
"""
@inline Base.angle(bound::Boundary, operator::Operator) = angle(operator.id, bound.vectors, bound.values)

"""
    angle(id::ID{AbstractOID}, vectors::AbstractVector{<:AbstractVector{<:Number}}, values::AbstractVector{<:Number}) -> Number

Get the total twist phase of an id.
"""
@inline @generated function Base.angle(id::ID{AbstractOID}, vectors::AbstractVector{<:AbstractVector{<:Number}}, values::AbstractVector{<:Number})
    Expr(:call, :+, [:(angle(id[$i], vectors, values)) for i = 1:rank(id)]...)
end

"""
    update!(bound::Boundary, args...; kwargs...) -> Boundary

Update the values of the boundary twisted phase.
"""
@inline @generated function update!(bound::Boundary, args...; kwargs...)
    exprs = []
    for (i, name) in enumerate(QuoteNode.(keys(bound)))
        push!(exprs, :(bound.values[$i] = get(kwargs, $name, bound.values[$i])))
    end
    return Expr(:block, exprs..., :(return bound))
end

"""
    plain

Plain boundary condition without any twist.
"""
const plain = Boundary{()}(Float[], SVector{0, Float}[])
@inline Base.angle(::typeof(plain), operator::Operator) = 0
@inline (::typeof(plain))(operator::Operator) = operator

end #module
