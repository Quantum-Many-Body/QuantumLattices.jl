module DegreesOfFreedom

using Printf: @printf, @sprintf
using StaticArrays: SVector
using LaTeXStrings: latexstring
using ..QuantumAlgebras: SimpleID, ID, Element, Elements
using ..Spatials: AbstractPID, Point
using ...Essentials: dtype
using ...Interfaces: id, value, decompose, dimension
using ...Prerequisites: Float, decimaltostr
using ...Prerequisites.Traits: rawtype, efficientoperations
using ...Prerequisites.CompositeStructures: CompositeDict
using ...Prerequisites.VectorSpaces: CartesianVectorSpace

import ..QuantumAlgebras: sequence
import ..Spatials: pidtype, rcoord, icoord
import ...Essentials: reset!, update!
import ...Interfaces: rank, ⊗
import ...Prerequisites.Traits: contentnames, getcontent
import ...Prerequisites.VectorSpaces: shape, ndimshape

export IID, SimpleIID, CompositeIID, Internal, InternalIndex, SimpleInternal, CompositeInternal
export Hilbert, AbstractOID, Index, AbstractCompositeOID, OID, AbstractOperator, Operator, Operators, Metric, OIDToTuple, Table
export internaltype, iidtype, isHermitian, indextype, oidtype
export LaTeX, latexname, latexformat, superscript, subscript, script
export Boundary, twist, plain

"""
    IID <: SimpleID

The id of an internal degree of freedom.
"""
abstract type IID <: SimpleID end

"""
    SimpleIID <: IID

The id of a simple internal degree of freedom.
"""
abstract type SimpleIID <: IID end

"""
    CompositeIID{T<:Tuple{Vararg{SimpleIID}}} <: IID

The composition of several single internal ids.
"""
struct CompositeIID{T<:Tuple{Vararg{SimpleIID}}} <: IID
    content::T
end
Base.show(io::IO, ciid::CompositeIID) = @printf io "%s" join((string(ciid[i]) for i = 1:rank(ciid)), " ⊗ ")
@inline Base.length(ciid::CompositeIID) = length(typeof(ciid))
@inline Base.length(::Type{<:CompositeIID{T}}) where {T<:Tuple{Vararg{SimpleIID}}} = fieldcount(T)
@inline Base.getindex(ciid::CompositeIID, i::Int) = ciid.content[i]
@inline Base.getproperty(ciid::CompositeIID, name::Symbol) = ciidgetproperty(ciid, Val(name))
@inline ciidgetproperty(ciid::CompositeIID, ::Val{:content}) = getfield(ciid, :content)
@inline ciidgetproperty(ciid::CompositeIID, ::Val{name}) where name = getproperty(getfield(ciid, :content), name)

"""
    CompositeIID(content::SimpleIID...)

Construct a composite iid from a set of simple iids.
"""
@inline CompositeIID(content::SimpleIID...) = CompositeIID(content)

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
@inline ⊗(iid::SimpleIID, ciid::CompositeIID) = CompositeIID(iid, ciid.content...)
@inline ⊗(ciid::CompositeIID, iid::SimpleIID) = CompositeIID(ciid.content..., iid)
@inline ⊗(ciid₁::CompositeIID, ciid₂::CompositeIID) = CompositeIID(ciid₁.content..., ciid₂.content...)

"""
    Internal{I<:IID} <: CartesianVectorSpace{I}

The whole internal degrees of freedom at a single point.
"""
abstract type Internal{I<:IID} <: CartesianVectorSpace{I} end

"""
    InternalIndex{I}

Index of a simple internal space in the composite internal spaces.
"""
struct InternalIndex{I} end
InternalIndex(I::Int) = InternalIndex{I}()

"""
    SimpleInternal{I<:SimpleIID} <: Internal{I}

The simple internal degrees of freedom at a single point.
"""
abstract type SimpleInternal{I<:SimpleIID} <: Internal{I} end
Base.show(io::IO, i::SimpleInternal) = @printf io "%s(%s)" i|>typeof|>nameof join(("$name=$(getfield(i, name))" for name in i|>typeof|>fieldnames), ", ")

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
    CompositeInternal{T<:Tuple{Vararg{SimpleInternal}}} <: Internal{CompositeIID}

The composition of several single internal spaces.
"""
struct CompositeInternal{T<:Tuple{Vararg{SimpleInternal}}} <: Internal{CompositeIID}
    content::T
end
@inline Base.eltype(ci::CompositeInternal) = eltype(typeof(ci))
@inline @generated function Base.eltype(::Type{<:CompositeInternal{T}}) where {T<:Tuple{Vararg{SimpleInternal}}}
    return CompositeIID{Tuple{[eltype(fieldtype(T, i)) for i = 1:fieldcount(T)]...}}
end
Base.show(io::IO, ci::CompositeInternal) = @printf io "%s" join((string(ci[InternalIndex(i)]) for i = 1:rank(ci)), " ⊗ ")

"""
    CompositeInternal(content::SimpleInternal...)

Construct a composite internal space from a set of simple internal spaces.
"""
@inline CompositeInternal(content::SimpleInternal...) = CompositeInternal(content)

"""
    rank(ci::CompositeInternal) -> Int
    rank(::Type{<:CompositeInternal{T}}) where {T<:Tuple{Vararg{SimpleInternal}}} -> Int

Get the number of simple internal spaces in a composite internal space.
"""
@inline rank(ci::CompositeInternal) = rank(typeof(ci))
@inline rank(::Type{<:CompositeInternal{T}}) where {T<:Tuple{Vararg{SimpleInternal}}} = fieldcount(T)

"""
    internaltype(ci::CompositeInternal, index::Int)
    internaltype(::Type{<:CompositeInternal{T}}, index::Int) where {T<:Tuple{Vararg{SimpleInternal}}}

Get the type of the ith simple internal space in a composite internal space.
"""
@inline internaltype(ci::CompositeInternal, index::Int) = internaltype(typeof(ci), index)
@inline internaltype(::Type{<:CompositeInternal{T}}, index::Int) where {T<:Tuple{Vararg{SimpleInternal}}} = fieldtype(T, index)

"""
    getindex(ci::CompositeInternal, ::InternalIndex{I}) where I

Get the Ith simple internal space in a composite internal space.
"""
@inline Base.getindex(ci::CompositeInternal, ::InternalIndex{I}) where I = ci.content[I]

@inline @generated shape(ci::CompositeInternal) = Expr(:tuple, [:(shape(ci[InternalIndex($i)])...) for i = 1:rank(ci)]...)
@inline ndimshape(::Type{<:CompositeInternal{T}}) where {T<:Tuple{Vararg{SimpleInternal}}} = sum(ndimshape(fieldtype(T, i)) for i = 1:fieldcount(T))
@inline @generated function Base.CartesianIndex(ciid::CompositeIID, ci::CompositeInternal)
    exprs = []
    for i = 1:rank(ci)
        push!(exprs, :(Tuple(CartesianIndex(ciid[$i], ci[InternalIndex($i)]))...))
    end
    return Expr(:call, :CartesianIndex, exprs...)
end
@inline CompositeIID(index::CartesianIndex, ci::CompositeInternal) = compositeiid(index, ci, segment(ci)|>Val)
@inline @generated segment(ci::CompositeInternal) = Expr(:tuple, [:(ndimshape(internaltype(ci, $i))) for i = 1:rank(ci)]...)
@inline @generated function compositeiid(index::CartesianIndex, ci::CompositeInternal, ::Val{dims}) where dims
    count = 1
    exprs = []
    for (i, dim) in enumerate(dims)
        cartesianindex = Expr(:call, :CartesianIndex, [:(index[$j]) for j = count:(count+dim-1)]...)
        push!(exprs, :(ci[InternalIndex($i)][$cartesianindex]))
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
@inline ⊗(i::SimpleInternal, ci::CompositeInternal) = CompositeInternal(i, ci.content...)
@inline ⊗(ci::CompositeInternal, i::SimpleInternal) = CompositeInternal(ci.content..., i)
@inline ⊗(ci₁::CompositeInternal, ci₂::CompositeInternal) = CompositeInternal(ci₁.content..., ci₂.content...)

"""
    filter(iid::SimpleIID, ci::CompositeInternal) -> Union{Nothing, SimpleInternal, CompositeInternal}
    filter(::Type{I}, ci::CompositeInternal) where {I<:SimpleIID} -> Union{Nothing, SimpleInternal, CompositeInternal}

Filter the composite internal space and select those that matches `I` or the type of `iid`.
"""
@inline Base.filter(iid::SimpleIID, ci::CompositeInternal) = filter(typeof(iid), ci)
@inline Base.filter(::Type{I}, ci::CompositeInternal) where {I<:SimpleIID} = filterhelper(I, ci, filtermatches(I, ci)|>Val)
@inline @generated function filtermatches(::Type{I}, ci::CompositeInternal) where {I<:SimpleIID}
    exprs = []
    for i = 1:rank(ci)
        T = internaltype(ci, i)
        push!(exprs, :(match(I, $T)))
    end
    Expr(:tuple, exprs...)
end
@inline @generated function filterhelper(::Type{I}, ci::CompositeInternal, ::Val{BS}) where {I<:SimpleIID, BS}
    exprs = []
    for (i, B) in enumerate(BS)
        B && push!(exprs, :(ci[InternalIndex($i)]))
    end
    length(exprs)==0 && return
    length(exprs)==1 && return first(exprs)
    return Expr(:call, :CompositeInternal, exprs...)
end

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

"""
    AbstractOID <: SimpleID

Abstract type of operator id.
"""
abstract type AbstractOID <: SimpleID end

"""
    isHermitian(id::ID{AbstractOID, N}) where N -> Bool

Judge whether an operator id is Hermitian.
"""
function isHermitian(id::ID{AbstractOID, N}) where N
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
    adjoint(index::Index) -> typeof(index)

Get the adjoint of an index.
"""
@inline Base.adjoint(index::Index) = rawtype(typeof(index))(index.pid, adjoint(index.iid))

"""
    AbstractCompositeOID{I<:Index} <: AbstractOID

The abstract type of composite operator id.
"""
abstract type AbstractCompositeOID{I<:Index} <: AbstractOID end
@inline contentnames(::Type{<:AbstractCompositeOID}) = (:index,)

"""
    indextype(::AbstractCompositeOID)
    indextype(::Type{<:AbstractCompositeOID})

Get the index type of an composite operator id.
"""
@inline indextype(oid::AbstractCompositeOID) = indextype(typeof(oid))
@inline indextype(::Type{<:AbstractCompositeOID{I}}) where {I<:Index} = I

"""
    OID{I<:Index, V<:SVector} <: AbstractCompositeOID{I}

Operator id.
"""
struct OID{I<:Index, V<:SVector} <: AbstractCompositeOID{I}
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
    AbstractOperator{V<:Number, I<:ID{AbstractOID}} <: Element{V, I}

Abstract type for an operator.
"""
abstract type AbstractOperator{V<:Number, I<:ID{AbstractOID}} <: Element{V, I} end

"""
    Operator{V<:Number, I<:ID{AbstractOID}} <: AbstractOperator{V, I}

Operator.
"""
struct Operator{V<:Number, I<:ID{AbstractOID}} <: AbstractOperator{V, I}
    value::V
    id::I
end
@inline Operator(value::Number) = Operator(value, ())
Base.show(io::IO, opt::Operator) = @printf io "%s(%s, %s)" nameof(typeof(opt)) decimaltostr(value(opt)) id(opt)

"""
    adjoint(opt::Operator) -> Operator

Get the adjoint of an operator.
"""
@inline Base.adjoint(opt::Operator) = rawtype(typeof(opt))(value(opt)', id(opt)')

"""
    isHermitian(opt::Operator) -> Bool

Judge whether an operator is Hermitian.
"""
@inline isHermitian(opt::Operator) = isa(value(opt), Real) && isHermitian(id(opt))

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
    Operators(opts::AbstractOperator...)

A set of operators.

Type alias for `Elements{<:ID{AbstractOID}, <:AbstractOperator}`.
"""
const Operators{I<:ID{AbstractOID}, O<:AbstractOperator} = Elements{I, O}
@inline Operators(opts::AbstractOperator...) = Elements(opts...)

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
    isHermitian(opts::Operators) -> Bool

Judge whether a set of operators as a whole is Hermitian.
"""
@inline isHermitian(opts::Operators) = opts == opts'

"""
    Metric <: Function

The rules for measuring a concrete oid so that oids can be compared.

As a function, every instance should accept only one positional argument, i.e. the concrete oid to be measured.
"""
abstract type Metric <: Function end
@inline Base.:(==)(m₁::T, m₂::T) where {T<:Metric} = ==(efficientoperations, m₁, m₂)
@inline Base.isequal(m₁::T, m₂::T) where {T<:Metric} = isequal(efficientoperations, m₁, m₂)
@inline Base.valtype(::Type{<:Metric}, ::Type{I}) where {I<:AbstractOID} = I
@inline (::Metric)(oid::AbstractOID) = oid

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
    OIDToTuple(::Type{I}) where {I<:AbstractCompositeOID}

Construct the conversion rule from the information of subtypes of `AbstractOID`.
"""
@inline @generated OIDToTuple(::Type{I}) where {I<:Index} = OIDToTuple(fieldnames(pidtype(I))..., (fieldnames(iidtype(I)))...)
@inline OIDToTuple(::Type{I}) where {I<:AbstractCompositeOID} = OIDToTuple(indextype(I))

"""
    OIDToTuple(::Type{H}) where {H<:Hilbert}

Construct the conversion rule from the information of `Hilbert`.
"""
@inline OIDToTuple(::Type{H}) where {H<:Hilbert} = OIDToTuple(Index{H|>keytype, H|>valtype|>eltype})

"""
    valtype(::Type{<:OIDToTuple}, ::Type{<:Index})
    valtype(::Type{<:OIDToTuple}, ::Type{<:AbstractCompositeOID})

Get the valtype of applying an `OIDToTuple` rule to a subtype of `AbstractOID`.
"""
@inline @generated function Base.valtype(::Type{M}, ::Type{I}) where {M<:OIDToTuple, I<:Index}
    types = []
    for field in keys(M)
        hasfield(pidtype(I), field) && push!(types, fieldtype(pidtype(I), field))
        hasfield(iidtype(I), field) && push!(types, fieldtype(iidtype(I), field))
    end
    return  Expr(:curly, :Tuple, types...)
end
@inline @generated Base.valtype(::Type{M}, ::Type{I}) where {M<:OIDToTuple, I<:AbstractCompositeOID} = valtype(M, indextype(I))

"""
    (oidtotuple::OIDToTuple)(index::Index) -> Tuple
    (oidtotuple::OIDToTuple)(oid::AbstractCompositeOID) -> Tuple

Convert a concrete oid to a tuple.
"""
@inline @generated function (oidtotuple::OIDToTuple)(index::Index)
    exprs = []
    for name in keys(oidtotuple)
        field = QuoteNode(name)
        hasfield(pidtype(index), name) && push!(exprs, :(getfield(index.pid, $field)))
        hasfield(iidtype(index), name) && push!(exprs, :(getfield(index.iid, $field)))
    end
    return Expr(:tuple, exprs...)
end
@inline (oidtotuple::OIDToTuple)(oid::AbstractCompositeOID) = oidtotuple(getcontent(oid, :index))

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
@inline Base.haskey(table::Table, oid::AbstractCompositeOID) = haskey(table, getcontent(oid, :index))
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
    sequence(opt::AbstractOperator, table::AbstractDict) -> NTuple{rank(opt), Int}

Get the sequence of the oids of an operator according to a table.
"""
@inline @generated function sequence(opt::AbstractOperator, table::AbstractDict{<:Any,Int})
    return Expr(:tuple, [:(table[id(opt)[$i]]) for i = 1:rank(opt)]...)
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
    repr(opt::AbstractOperator) -> String

Get the string representation of an operator in the LaTeX format.
"""
function Base.repr(opt::AbstractOperator)
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
    show(io::IO, ::MIME"text/latex", opt::AbstractOperator)

Show an operator.
"""
Base.show(io::IO, ::MIME"text/latex", opt::AbstractOperator) = show(io, MIME"text/latex"(), latexstring(repr(opt)))

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
    (bound::Boundary)(operator::AbstractOperator) -> AbstractOperator

Get the boundary twisted operator.
"""
@inline (bound::Boundary)(operator::AbstractOperator) = twist(operator, bound.vectors, bound.values)

"""
    twist(operator::AbstractOperator, vectors::AbstractVector{<:AbstractVector{<:Number}}, values::AbstractVector{<:Number}) -> AbstractOperator

Twist an operator.
"""
@inline function twist(operator::AbstractOperator, vectors::AbstractVector{<:AbstractVector{<:Number}}, values::AbstractVector{<:Number})
    return replace(operator, operator.value*exp(1im*angle(operator.id, vectors, values)))
end

"""
    angle(bound::Boundary, operator::AbstractOperator) -> Number

Get the boundary twist phase of an operator.
"""
@inline Base.angle(bound::Boundary, operator::AbstractOperator) = angle(operator.id, bound.vectors, bound.values)

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
@inline Base.angle(::typeof(plain), operator::AbstractOperator) = 0
@inline (::typeof(plain))(operator::AbstractOperator) = operator

end #module
