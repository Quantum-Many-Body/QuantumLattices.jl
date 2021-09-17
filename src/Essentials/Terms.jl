module Terms

using MacroTools: @capture
using Printf: @printf, @sprintf
using StaticArrays: SVector
using ..QuantumAlgebras: SimpleID, ID, Element, Elements, Transformation
using ..Spatials: AbstractPID, AbstractBond, Point, Bonds, AbstractLattice, acrossbonds, isintracell, pidtype
using ..DegreesOfFreedom: IID, SimpleIID, CompositeIID, Internal, CompositeInternal, InternalIndex
using ..DegreesOfFreedom: Hilbert, AbstractOID, Index, OID, oidtype, Operator, Operators, Table, Boundary, plain
using ...Essentials: dtype
using ...Interfaces: add!, dimension
using ...Prerequisites: atol, rtol, decimaltostr
using ...Prerequisites.Traits: rawtype, efficientoperations, commontype, parametertype
using ...Prerequisites.CompositeStructures: CompositeTuple, NamedContainer
using ...Prerequisites.VectorSpaces: CartesianVectorSpace

import ..QuantumAlgebras: idtype
import ..DegreesOfFreedom: isHermitian
import ...Interfaces: id, value, rank, expand, expand!
import ...Essentials: kind, update!, reset!
import ...Prerequisites.Traits: parameternames, isparameterbound, contentnames, getcontent
import ...Prerequisites.VectorSpaces: shape, ndimshape

export Subscript, Subscripts, SubscriptsID, IIDSpace, @subscript_str, subscriptexpr, wildcard, diagonal, noconstrain
export AbstractCoupling, Coupling, Couplings, couplingcenters, couplingpoints, couplinginternals, @couplings
export TermFunction, TermAmplitude, TermCouplings, TermModulate, Term, abbr, ismodulatable, otype
export Parameters, GenOperators, AbstractGenerator, AbstractCompleteGenerator, Generator, AbstractSimplifiedGenerator, SimplifiedGenerator

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
    SubscriptsID{T<:Tuple{Vararg{Pair{UnitRange{Int}, <:Tuple{Vararg{String}}}}}} <: SimpleID

The id of the subscripts.
"""
struct SubscriptsID{T<:Tuple{Vararg{Pair{UnitRange{Int}, <:Tuple{Vararg{String}}}}}} <: SimpleID
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
    Expr(:tuple, [:(shape(IIDSpace(iidspace.iid[$i], iidspace.internal[InternalIndex($i)], $Kind))...) for i = 1:rank(I)]...)
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
    AbstractCoupling{V, I<:ID{SimpleID}} <: Element{V, I}

The abstract coupling intra/inter internal degrees of freedom at different lattice points.
"""
abstract type AbstractCoupling{V, I<:ID{SimpleID}} <: Element{V, I} end
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

Alias for `Elements{<:ID{SimpleID}, <:AbstractCoupling}`.
"""
const Couplings{I<:ID{SimpleID}, C<:AbstractCoupling} = Elements{I, C}
@inline Couplings(cps::AbstractCoupling...) = Elements(cps...)
@inline Couplings(cps::Couplings) = cps

"""
    @couplings cps -> Couplings

Convert an expression/literal to a set of couplings.
"""
macro couplings(cps) :(Couplings($(esc(cps)))) end

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
@inline (termmodulate::TermModulate{Val{true}, id})(; kwargs...) where id = get(kwargs, id, nothing)
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
    isHermitian(term::Term) -> Bool
    isHermitian(::Type{<:Term}) -> Bool
"""
@inline isHermitian(term::Term) = isHermitian(typeof(term))
@inline isHermitian(::Type{<:Term}) = error("isHermitian error: not implemented.")

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
    expand!(operators::Operators, term::Term, bond::AbstractBond, hilbert::Hilbert, half::Bool=false; table::Union{Nothing, Table}=nothing) -> Operators
    expand!(operators::Operators, term::Term, bonds::Bonds, hilbert::Hilbert, half::Bool=false; table::Union{Nothing, Table}=nothing) -> Operators

Expand the operators of a term on a bond/set-of-bonds with a given Hilbert space.

The `half` parameter determines the behavior of generating operators, which falls into the following two categories
* `true`: "Hermitian half" of the generated operators
* `false`: "Hermitian whole" of the generated operators
"""
function expand!(operators::Operators, term::Term, bond::AbstractBond, hilbert::Hilbert, half::Bool=false; table::Union{Nothing, Table}=nothing)
    if term.bondkind == bond|>kind
        value = term.value * term.amplitude(bond) * term.factor
        if abs(value) ≠ 0
            hermitian = isHermitian(term)
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
@generated function expand!(operators::Operators, term::Term, bonds::Bonds, hilbert::Hilbert, half::Bool=false; table::Union{Nothing, Table}=nothing)
    exprs = []
    for i = 1:rank(bonds)
        push!(exprs, :(for bond in bonds.bonds[$i] expand!(operators, term, bond, hilbert, half; table=table) end))
    end
    push!(exprs, :(return operators))
    return Expr(:block, exprs...)
end
@inline termfactor(id::ID{AbstractOID}, ::Val) = isHermitian(id) ? 2 : 1

"""
    expand(term::Term, bond::AbstractBond, hilbert::Hilbert, half::Bool=false; table::Union{Nothing, Table}=nothing) -> Operators
    expand(term::Term, bonds::Bonds, hilbert::Hilbert, half::Bool=false; table::Union{Nothing, Table}=nothing) -> Operators

Expand the operators of a term on a bond/set-of-bonds with a given Hilbert space.
"""
@inline function expand(term::Term, bond::AbstractBond, hilbert::Hilbert, half::Bool=false; table::Union{Nothing, Table}=nothing)
    optype = otype(term|>typeof, hilbert|>typeof, bond|>typeof)
    expand!(Operators{idtype(optype), optype}(), term, bond, hilbert, half; table=table)
end
@inline function expand(term::Term, bonds::Bonds, hilbert::Hilbert, half::Bool=false; table::Union{Nothing, Table}=nothing)
    optype = otype(term|>typeof, hilbert|>typeof, bonds|>eltype)
    expand!(Operators{idtype(optype), optype}(), term, bonds, hilbert, half, table=table)
end

"""
    Parameters{Names}(values::Number...) where Names

A NamedTuple that contain the key-value pairs.
"""
const Parameters{Names} = NamedContainer{Number, Names}
@inline Parameters{Names}(values::Number...) where {Names} = NamedContainer{Names}(values)

"""
    match(params₁::Parameters, params₂::Parameters, atol=atol, rtol=rtol) -> Bool

Judge whether the second set of parameters matches the first.
"""
@generated function Base.match(params₁::Parameters, params₂::Parameters, atol=atol, rtol=rtol)
    names = intersect(fieldnames(params₁), fieldnames(params₂))
    length(names)==0 && return true
    name = QuoteNode(names[1])
    expr = :(isapprox(getfield(params₁, $name), getfield(params₂, $name), atol=atol, rtol=rtol))
    for i = 2:length(names)
        name = QuoteNode(names[i])
        expr = Expr(:&&, expr, :(isapprox(getfield(params₁, $name), getfield(params₂, $name), atol=atol, rtol=rtol)))
    end
    return expr
end

"""
    GenOperators{C<:Operators, A<:NamedContainer{Operators}, B<:NamedContainer{Operators}}

A set of operators. This is the core of [`AbstractGenerator`](@ref).
"""
struct GenOperators{C<:Operators, A<:NamedContainer{Operators}, B<:NamedContainer{Operators}}
    constops::C
    alterops::A
    boundops::B
end
@inline Base.:(==)(genops₁::GenOperators, genops₂::GenOperators) = ==(efficientoperations, genops₁, genops₂)
@inline Base.isequal(genops₁::GenOperators, genops₂::GenOperators) = isequal(efficientoperations, genops₁, genops₂)

"""
    GenOperators(terms::Tuple{Vararg{Term}}, bonds::Bonds, hilbert::Hilbert, half::Bool; table::Union{Nothing, Table}=nothing)

Construct a set of operators.
"""
@generated function GenOperators(terms::Tuple{Vararg{Term}}, bonds::Bonds, hilbert::Hilbert, half::Bool; table::Union{Nothing, Table}=nothing)
    constterms, alterterms = [], []
    for term in fieldtypes(terms)
        ismodulatable(term) ? push!(alterterms, term) : push!(constterms, term)
    end
    names = NTuple{fieldcount(terms), Symbol}(id(term) for term in fieldtypes(terms))
    alternames = NTuple{length(alterterms), Symbol}(id(term) for term in alterterms)
    exprs, alterops, boundops = [], [], []
    push!(exprs, quote
        choosedterms = length($constterms)>0 ? $constterms : $alterterms
        constoptp = otype(choosedterms[1], hilbert|>typeof, bonds|>eltype)
        constidtp = constoptp |> idtype
        for i = 2:length(choosedterms)
            tempoptp = otype(choosedterms[i], hilbert|>typeof, bonds|>eltype)
            constoptp = promote_type(constoptp, tempoptp)
            constidtp = promote_type(constidtp, tempoptp|>idtype)
        end
        constops = Operators{constidtp, constoptp}()
        innerbonds = filter(acrossbonds, bonds, Val(:exclude))
        boundbonds = filter(acrossbonds, bonds, Val(:include))
    end)
    for i = 1:fieldcount(terms)
        push!(boundops, :(expand(one(terms[$i]), boundbonds, hilbert, half, table=table)))
        if ismodulatable(fieldtype(terms, i))
            push!(alterops, :(expand(one(terms[$i]), innerbonds, hilbert, half, table=table)))
        else
            push!(exprs, :(expand!(constops, terms[$i], innerbonds, hilbert, half, table=table)))
        end
    end
    alterops = Expr(:tuple, alterops...)
    boundops = Expr(:tuple, boundops...)
    push!(exprs, quote
        alterops = NamedContainer{$alternames}($alterops)
        boundops = NamedContainer{$names}($boundops)
        return GenOperators(constops, alterops, boundops)
    end)
    return Expr(:block, exprs...)
end

"""
    eltype(ops::GenOperators)
    eltype(::Type{<:GenOperators})

Get the eltype of a set of operators, which is defined to be the common operator type of all operators it contains.
"""
@inline Base.eltype(ops::GenOperators) = eltype(typeof(ops))
@inline @generated function Base.eltype(::Type{<:GenOperators{S, D, B}}) where {S<:Operators, D<:NamedContainer{Operators}, B<:NamedContainer{Operators}}
    optp = S |> valtype
    (fieldcount(D) > 0) && (optp = promote_type(optp, mapreduce(valtype, promote_type, fieldtypes(D))))
    (fieldcount(B) > 0) && (optp = promote_type(optp, mapreduce(valtype, promote_type, fieldtypes(B))))
    return optp
end

"""
    idtype(ops::GenOperators)
    idtype(::Type{<:GenOperators})

Get the idtype of a set of operators, which is defined to be the common idtype of all operators it contains.
"""
@inline idtype(ops::GenOperators) = idtype(typeof(ops))
@inline @generated function idtype(::Type{<:GenOperators{S, D, B}}) where {S<:Operators, D<:NamedContainer{Operators}, B<:NamedContainer{Operators}}
    idtp = S |> keytype
    (fieldcount(D) > 0) && (idtp = promote_type(idtp, mapreduce(keytype, promote_type, fieldtypes(D))))
    (fieldcount(B) > 0) && (idtp = promote_type(idtp, mapreduce(keytype, promote_type, fieldtypes(B))))
    return idtp
end

"""
    empty!(ops::GenOperators) -> GenOperators

Empty a set of operators.
"""
@generated function Base.empty!(ops::GenOperators)
    exprs = [:(empty!(ops.constops))]
    for i = 1:fieldcount(fieldtype(ops, :alterops)) push!(exprs, :(empty!(ops.alterops[$i]))) end
    for i = 1:fieldcount(fieldtype(ops, :boundops)) push!(exprs, :(empty!(ops.boundops[$i]))) end
    push!(exprs, :(return ops))
    return Expr(:block, exprs...)
end

"""
    empty(ops::GenOperators) -> GenOperators

Get an empty copy of a set of operators.
"""
@generated function Base.empty(ops::GenOperators)
    exprs = [:(constops = empty(ops.constops))]
    alterops, boundops = [], []
    for i = 1:fieldcount(fieldtype(ops, :alterops)) push!(alterops, :(empty(ops.alterops[$i]))) end
    for i = 1:fieldcount(fieldtype(ops, :boundops)) push!(boundops, :(empty(ops.boundops[$i]))) end
    alterops = Expr(:tuple, alterops...)
    boundops = Expr(:tuple, boundops...)
    push!(exprs, quote
        alterops = NamedContainer{fieldnames(fieldtype(ops|>typeof, :alterops))}($alterops)
        boundops = NamedContainer{fieldnames(fieldtype(ops|>typeof, :boundops))}($boundops)
        return GenOperators(constops, alterops, boundops)
    end)
    return Expr(:block, exprs...)
end

"""
    merge!(genops::GenOperators, another::GenOperators) -> GenOperators

Merge the operators from `another` to `genops`.
"""
@generated function Base.merge!(genops::GenOperators, another::GenOperators)
    exprs = [:(merge!(genops.constops, another.constops))]
    for name in QuoteNode.(fieldnames(fieldtype(genops, :alterops)))
        push!(exprs, :(merge!(getfield(genops.alterops, $name), getfield(another.alterops, $name))))
    end
    for name in QuoteNode.(fieldnames(fieldtype(genops, :boundops)))
        push!(exprs, :(merge!(getfield(genops.boundops, $name), getfield(another.boundops, $name))))
    end
    push!(exprs, :(return genops))
    return Expr(:block, exprs...)
end

"""
    merge(genops::GenOperators, another::GenOperators) -> GenOperators

Merge the operators from `another` and `genops`.
"""
@inline Base.merge(genops::GenOperators, another::GenOperators) = merge!(empty(genops), another)

"""
    expand!(operators::Operators, ops::GenOperators, boundary::Boundary; kwargs...) -> Operators

Expand the operators with the given boundary twist and term coefficients.
"""
@generated function expand!(operators::Operators, ops::GenOperators, boundary::Boundary; kwargs...)
    exprs = [:(add!(operators, ops.constops))]
    for name in QuoteNode.(fieldnames(fieldtype(ops, :alterops)))
        push!(exprs, :(value = get(kwargs, $name, nothing)))
        push!(exprs, :(for opt in values(getfield(ops.alterops, $name)) add!(operators, opt*value) end))
    end
    for name in QuoteNode.(fieldnames(fieldtype(ops, :boundops)))
        push!(exprs, :(value = get(kwargs, $name, nothing)))
        push!(exprs, :(for opt in values(getfield(ops.boundops, $name)) add!(operators, boundary(opt)*value) end))
    end
    push!(exprs, :(return operators))
    return Expr(:block, exprs...)
end

"""
    reset!(genops::GenOperators, terms::Tuple{Vararg{Term}}, bonds::Bonds, hilbert::Hilbert, half::Bool; table::Union{Nothing, Table}=nothing) -> GenOperators

Reset a set of operators by the new terms, bonds and Hilbert space.
"""
@generated function reset!(genops::GenOperators, terms::Tuple{Vararg{Term}}, bonds::Bonds, hilbert::Hilbert, half::Bool; table::Union{Nothing, Table}=nothing)
    exprs = []
    push!(exprs, quote
        empty!(genops)
        innerbonds = filter(acrossbonds, bonds, Val(:exclude))
        boundbonds = filter(acrossbonds, bonds, Val(:include))
    end)
    for (i, term) in enumerate(fieldtypes(terms))
        name = QuoteNode(term|>id)
        push!(exprs, :(expand!(getfield(genops.boundops, $name), one(terms[$i]), boundbonds, hilbert, half, table=table)))
        if ismodulatable(term)
            push!(exprs, :(expand!(getfield(genops.alterops, $name), one(terms[$i]), innerbonds, hilbert, half, table=table)))
        else
            push!(exprs, :(expand!(genops.constops, terms[$i], innerbonds, hilbert, half, table=table)))
        end
    end
    push!(exprs, :(return genops))
    return Expr(:block, exprs...)
end

"""
    (transformation::Transformation)(genops::GenOperators) -> GenOperators

Get the transformed set of operators.
"""
function (transformation::Transformation)(genops::GenOperators)
    constops = transformation(genops.constops)
    alterops = NamedContainer{keys(genops.alterops)}(map(transformation, values(genops.alterops)))
    boundops = NamedContainer{keys(genops.boundops)}(map(transformation, values(genops.boundops)))
    return GenOperators(constops, alterops, boundops)
end

"""
    AbstractGenerator{T<:Union{Table, Nothing}, B<:Boundary, OS<:GenOperators}

Abstract generator.

By protocol, a concrete generator should have the following predefined contents:
* `table::T`: the index-sequence table if it is not nothing
* `boundary::B`: boundary twist for the generated operators
* `operators::OS`: the generated operators
"""
abstract type AbstractGenerator{T<:Union{Table, Nothing}, B<:Boundary, OS<:GenOperators} end
@inline contentnames(::Type{<:AbstractGenerator}) = (:table, :boundary, :operators)
@inline Base.:(==)(gen₁::AbstractGenerator, gen₂::AbstractGenerator) = ==(efficientoperations, gen₁, gen₂)
@inline Base.isequal(gen₁::AbstractGenerator, gen₂::AbstractGenerator) = isequal(efficientoperations, gen₁, gen₂)

"""
    otype(gen::AbstractGenerator)
    otype(::Type{<:AbstractGenerator})

Get the operator type of the generated operators.
"""
@inline otype(gen::AbstractGenerator) = otype(typeof(gen))
@inline otype(::Type{<:AbstractGenerator{<:Union{Table, Nothing}, <:Boundary, OS}}) where {OS<:GenOperators} = eltype(OS)

"""
    expand!(operators::Operators, gen::AbstractGenerator) -> Operators

Expand the operators of a generator.
"""
@inline expand!(operators::Operators, gen::AbstractGenerator) = expand!(operators, getcontent(gen, :operators), getcontent(gen, :boundary); Parameters(gen)...)

"""
    expand(gen::AbstractGenerator) -> Operators

Expand the operators of a generator.
"""
@inline expand(gen::AbstractGenerator) = expand!(Operators{idtype(otype(gen)), otype(gen)}(), gen)

"""
    AbstractCompleteGenerator{
        TS<:NamedContainer{Term},
        BS<:Bonds,
        H<:Hilbert,
        T<:Union{Table, Nothing},
        B<:Boundary,
        OS<:GenOperators
        } <: AbstractGenerator{T, B, OS}

Abstract complete generator.

By protocol, a concrete complete generator should have the following predefined contents in addition to those contained in its supertype:
* `terms::TS`: the terms contained in a complete generator
* `bonds::BS`: the bonds on which the terms are defined
* `hilbert::H`: the total Hilbert space
* `half::Bool`: true for generating an Hermitian half of the operators and false for generating the whole
"""
abstract type AbstractCompleteGenerator{
    TS<:NamedContainer{Term},
    BS<:Bonds,
    H<:Hilbert,
    T<:Union{Table, Nothing},
    B<:Boundary,
    OS<:GenOperators
    } <: AbstractGenerator{T, B, OS}
end
@inline contentnames(::Type{<:AbstractCompleteGenerator}) = (:terms, :bonds, :hilbert, :half, :table, :boundary, :operators)

"""
    Parameters(gen::AbstractCompleteGenerator) -> Parameters

Get the parameters of the terms of a generator.
"""
@generated function Parameters(gen::AbstractCompleteGenerator{TS}) where {TS<:NamedContainer{Term}}
    names, values = fieldnames(TS), [:(getcontent(gen, :terms)[$i].value) for i = 1:fieldcount(TS)]
    return :(Parameters{$names}($(values...)))
end

"""
    expand(gen::AbstractCompleteGenerator, name::Symbol) -> Operators
    expand(gen::AbstractCompleteGenerator, i::Int) -> Operators
    expand(gen::AbstractCompleteGenerator, name::Symbol, i::Int) -> Operators

Expand the operators of a generator:
1) the operators of a specific term;
2) the operators on a specific bond;
3) the operators of a specific term on a specific bond.
"""
function expand(gen::AbstractCompleteGenerator, name::Symbol)
    bonds = getcontent(gen, :bonds)
    hilbert = getcontent(gen, :hilbert)
    ops = getcontent(gen, :operators)
    term = getfield(getcontent(gen, :terms), name)
    optp = otype(term|>typeof, hilbert|>typeof, bonds|>eltype)
    result = Operators{idtype(optp), optp}()
    if ismodulatable(term)
        for opt in getfield(ops.alterops, name)|>values add!(result, opt*term.value) end
    else
        expand!(result, term, filter(acrossbonds, bonds, Val(:exclude)), hilbert, getcontent(gen, :half), table=getcontent(gen, :table))
    end
    for opt in getfield(ops.boundops, name)|>values
        add!(result, getcontent(gen, :boundary)(opt)*term.value)
    end
    return result
end
@generated function expand(gen::AbstractCompleteGenerator{TS}, i::Int) where {TS<:NamedContainer{Term}}
    exprs = []
    push!(exprs, quote
        bond = getcontent(gen, :bonds)[i]
        result = Operators{idtype(otype(gen)), otype(gen)}()
    end)
    for i = 1:fieldcount(TS)
        push!(exprs, :(expand!(result, getcontent(gen, :terms)[$i], bond, getcontent(gen, :hilbert), getcontent(gen, :half), table=getcontent(gen, :table))))
    end
    push!(exprs, quote
        isintracell(bond) && for opt in values(result)
            result[id(opt)] = getcontent(gen, :boundary)(opt)
        end
        return result
    end)
    return Expr(:block, exprs...)
end
function expand(gen::AbstractCompleteGenerator, name::Symbol, i::Int)
    bond = getcontent(gen, :bonds)[i]
    result = expand(getfield(getcontent(gen, :terms), name), bond, getcontent(gen, :hilbert), getcontent(gen, :half), table=getcontent(gen, :table))
    isintracell(bond) && for opt in values(result)
        result[id(opt)] = getcontent(gen, :boundary)(opt)
    end
    return result
end

"""
    update!(gen::AbstractCompleteGenerator; kwargs...) -> typeof(gen)

Update the coefficients of the terms in a generator.
"""
@generated function update!(gen::AbstractCompleteGenerator{TS}; kwargs...) where {TS<:NamedContainer{Term}}
    exprs = [:(update!(getcontent(gen, :boundary); kwargs...))]
    for i = 1:fieldcount(TS)
        ismodulatable(fieldtype(TS, i)) && push!(exprs, :(update!(getcontent(gen, :terms)[$i]; kwargs...)))
    end
    push!(exprs, :(return gen))
    return Expr(:block, exprs...)
end

"""
    Generator{
        TS<:NamedContainer{Term},
        BS<:Bonds,
        H<:Hilbert,
        T<:Union{Table, Nothing},
        B<:Boundary,
        OS<:GenOperators
        } <: AbstractCompleteGenerator{TS, BS, H, T, B, OS}

A generator of operators based on terms, configuration of internal degrees of freedom, and boundary twist.
"""
struct Generator{
        TS<:NamedContainer{Term},
        BS<:Bonds,
        H<:Hilbert,
        T<:Union{Table, Nothing},
        B<:Boundary,
        OS<:GenOperators
        } <: AbstractCompleteGenerator{TS, BS, H, T, B, OS}
    terms::TS
    bonds::BS
    hilbert::H
    half::Bool
    table::T
    boundary::B
    operators::OS
end

"""
    Generator(terms::Tuple{Vararg{Term}}, bonds::Bonds, hilbert::Hilbert;
        half::Bool=false, table::Union{Table,Nothing}=nothing, boundary::Boundary=plain
        )

Construct a generator of operators.
"""
@inline function Generator(terms::Tuple{Vararg{Term}}, bonds::Bonds, hilbert::Hilbert;
        half::Bool=false, table::Union{Table,Nothing}=nothing, boundary::Boundary=plain
        )
    Generator(namedterms(terms), bonds, hilbert, half, table, boundary, GenOperators(terms, bonds, hilbert, half, table=table))
end
@generated function namedterms(terms::Tuple{Vararg{Term}})
    names = NTuple{fieldcount(terms), Symbol}(id(fieldtype(terms, i)) for i = 1:fieldcount(terms))
    return :(NamedContainer{$names}(terms))
end

"""
    empty!(gen::Generator) -> Generator

Empty the :bonds, :hilbert, :table and :operators of a generator.
"""
function Base.empty!(gen::Generator)
    empty!(gen.bonds)
    empty!(gen.hilbert)
    isnothing(gen.table) || empty!(gen.table)
    empty!(gen.operators)
    return gen
end

"""
    empty(gen::Generator) -> Generator

Get an empty copy of a generator.
"""
@inline function Base.empty(gen::Generator)
    Generator(gen.terms,
        empty(gen.bonds),
        empty(gen.hilbert),
        gen.half,
        isnothing(gen.table) ? nothing : empty(gen.table),
        gen.boundary,
        empty(gen.operators)
        )
end

"""
    reset!(gen::Generator, lattice::AbstractLattice) -> Generator

Reset a generator by a new lattice.
"""
function reset!(gen::Generator, lattice::AbstractLattice)
    reset!(gen.bonds, lattice)
    reset!(gen.hilbert, lattice.pids)
    isnothing(gen.table) || reset!(gen.table, gen.hilbert)
    reset!(gen.operators, Tuple(gen.terms), gen.bonds, gen.hilbert, gen.half, table=gen.table)
    return gen
end

"""
    AbstractSimplifiedGenerator{P<:Parameters, T<:Union{Table, Nothing}, B<:Boundary, OS<:GenOperators} <: AbstractGenerator{T, B, OS}

Abstract simplified generator.

By protocol, a concrete simplified generator should have the following predefined contents in addition to those contained in its supertype:
* `parameters::P`: the parameters of the generator
"""
abstract type AbstractSimplifiedGenerator{P<:Parameters, T<:Union{Table, Nothing}, B<:Boundary, OS<:GenOperators} <: AbstractGenerator{T, B, OS} end
@inline contentnames(::Type{<:AbstractSimplifiedGenerator}) = (:parameters, :table, :boundary, :operators)

"""
    Parameters(gen::AbstractSimplifiedGenerator) -> Parameters

Get the parameters of the terms of a generator.
"""
@inline Parameters(gen::AbstractSimplifiedGenerator) = getcontent(gen, :parameters)

"""
    SimplifiedGenerator{P<:Parameters, T<:Union{Table, Nothing}, B<:Boundary, OS<:GenOperators} <: AbstractSimplifiedGenerator{P, T, B, OS}

The simplified generator of operators.
"""
mutable struct SimplifiedGenerator{P<:Parameters, T<:Union{Table, Nothing}, B<:Boundary, OS<:GenOperators} <: AbstractSimplifiedGenerator{P, T, B, OS}
    parameters::P
    table::T
    boundary::B
    operators::OS
end

"""
    SimplifiedGenerator(parameters::Parameters, operators::GenOperators; table::Union{Table,Nothing}=nothing, boundary::Boundary=plain)

Construct an simplified generator of operators.
"""
@inline function SimplifiedGenerator(parameters::Parameters, operators::GenOperators; table::Union{Table,Nothing}=nothing, boundary::Boundary=plain)
    return SimplifiedGenerator(parameters, table, boundary, operators)
end

"""
    empty!(gen::SimplifiedGenerator) -> SimplifiedGenerator

Empty the operators of an simplified generator.
"""
@inline function Base.empty!(gen::SimplifiedGenerator)
    empty!(gen.operators)
    isnothing(gen.table) || empty!(gen.table)
    return gen
end

"""
    empty(gen::SimplifiedGenerator) -> SimplifiedGenerator

Get an empty copy of an simplified generator.
"""
@inline function Base.empty(gen::SimplifiedGenerator)
    return SimplifiedGenerator(gen.parameters, isnothing(gen.table) ? nothing : empty(gen.table), gen.boundary, empty(gen.operators))
end

"""
    update!(gen::SimplifiedGenerator; kwargs...) -> typeof(gen)

Update the parameters of a generator.
"""
@generated function update!(gen::SimplifiedGenerator; kwargs...)
    names = fieldnames(fieldtype(gen, :parameters))
    values = [:(get(kwargs, $name, getfield(gen.parameters, $name))) for name in QuoteNode.(names)]
    return Expr(:block, :(update!(gen.boundary; kwargs...)), :(gen.parameters = Parameters{$names}($(values...))), :(return gen))
end

"""
    reset!(gen::SimplifiedGenerator, operators::GenOperators; table::Union{Table, Nothing}=nothing) -> SimplifiedGenerator

Reset a generator by a new lattice.
"""
function reset!(gen::SimplifiedGenerator, operators::GenOperators; table::Union{Table, Nothing}=nothing)
    isnothing(gen.table) || isnothing(table) || merge!(empty!(gen.table), table)
    merge!(empty!(gen.operators), operators)
    return gen
end

"""
    (transformation::Transformation)(gen::AbstractGenerator;
        table::Union{Table, Nothing}=getcontent(gen, :table),
        boundary::Boundary=getcontent(gen, :boundary)
        )

Get the result of a transformation on the generated operators.
"""
function (transformation::Transformation)(gen::AbstractGenerator;
        table::Union{Table, Nothing}=getcontent(gen, :table),
        boundary::Boundary=getcontent(gen, :boundary)
        )
    return SimplifiedGenerator(Parameters(gen), transformation(getcontent(gen, :operators)), table=table, boundary=boundary)
end

end  # module
