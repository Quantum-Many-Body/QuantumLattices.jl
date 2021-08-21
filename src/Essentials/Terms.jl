module Terms

using Printf: @printf, @sprintf
using StaticArrays: SVector
using MacroTools: @capture
using ..Spatials: AbstractBond, Point, Bonds, AbstractLattice, acrossbonds, isintracell
using ..DegreesOfFreedom: Config, AbstractOID, oidtype, Operator, Operators, Table, Boundary, plain
using ...Interfaces: add!
using ...Prerequisites: atol, rtol, decimaltostr
using ...Prerequisites.Traits: rawtype, efficientoperations
using ...Prerequisites.CompositeStructures: CompositeTuple, NamedContainer
using ...Mathematics.AlgebraOverFields: SimpleID, ID, Element, Elements

import ..DegreesOfFreedom: isHermitian
import ...Interfaces: id, rank, expand, expand!, dimension
import ...Essentials: kind, update!, reset!
import ...Prerequisites.Traits: contentnames, getcontent
import ...Mathematics.AlgebraOverFields: idtype

export Subscripts, @subscripts, @subscripts_str, SubID
export Coupling, Couplings, @couplings, couplingcenters, couplingpoints, couplinginternals
export TermFunction, TermAmplitude, TermCouplings, TermModulate, Term, abbr, ismodulatable, otype
export Parameters, GenOperators, AbstractGenerator, Generator

const wildcard = Symbol("*")
@inline nonconstrain(_...) = true

"""
    Subscripts{DS, RS, OP<:Tuple, CP<:AbstractString, M<:Tuple{Vararg{Function}}, C<:Tuple{Vararg{Function}}} <: CompositeTuple{OP}

A subscript set of a certain internal degree of freedom.
"""
struct Subscripts{DS, RS, OP<:Tuple, CP<:Tuple, M<:Tuple{Vararg{Function}}, C<:Tuple{Vararg{Function}}} <: CompositeTuple{OP}
    opattern::OP
    cpattern::CP
    mapping::M
    constrain::C
    function Subscripts{DS, RS}(opattern::Tuple, cpattern::Tuple, mapping::Tuple{Vararg{Function}}, constrain::Tuple{Vararg{Function}}) where {DS, RS}
        @assert isa(DS, Tuple{Vararg{Int}}) && isa(RS, Tuple{Vararg{Int}}) "Subscripts error: wrong ranks and dimensions."
        @assert length(DS)==length(RS)==length(cpattern)==length(mapping)==length(constrain) "Subscripts error: dismatched number of segments."
        @assert sum(DS)==length(opattern) "Subscripts error: dismatched total dimension."
        new{DS, RS, typeof(opattern), typeof(cpattern), typeof(mapping), typeof(constrain)}(opattern, cpattern, mapping, constrain)
    end
end
@inline contentnames(::Type{<:Subscripts}) = (:contents, :cpattern, :mapping, :constrain)
@inline getcontent(subscripts::Subscripts, ::Val{:contents}) = subscripts.opattern

"""
    Subscripts(N::Int)
    Subscripts(opattern::Tuple{Vararg{Union{Integer, Symbol}}})
    Subscripts(opattern::Tuple{Vararg{Union{Integer, Symbol}}}, cpattern::AbstractString)

Construct a subscript set of a certain internal degree of freedom.
"""
@inline Subscripts(N::Int) = Subscripts(Val(N))
function Subscripts(::Val{N}) where N
    mapping = value::Int -> ntuple(i->value, Val(N))
    opattern = ntuple(i->wildcard, Val(N))
    return Subscripts{(N,), (1,)}(opattern, ("nonconstrain",), (mapping,), (nonconstrain,))
end
function Subscripts(opattern::Tuple{Vararg{Union{Integer, Symbol}}})
    ipattern = Tuple(subscriptsipattern(opattern))
    mapping = Meta.eval(Expr(:->, Expr(:tuple, ipattern...), Expr(:block, Expr(:tuple, opattern...))))
    cpattern = "nonconstrain"
    return Subscripts{(length(opattern),), (length(ipattern),)}(opattern, (cpattern,), (mapping,), (nonconstrain,))
end
function Subscripts(opattern::Tuple{Vararg{Union{Integer, Symbol}}}, cpattern::AbstractString)
    ipattern = Tuple(subscriptsipattern(opattern))
    mapping = Meta.eval(Expr(:->, Expr(:tuple, ipattern...), Expr(:block, Expr(:tuple, opattern...))))
    constrain = Meta.eval(Expr(:->, Expr(:tuple, ipattern...), Expr(:block, Meta.parse(cpattern))))
    return Subscripts{(length(opattern),), (length(ipattern),)}(opattern, (cpattern,), (mapping,), (constrain,))
end
function subscriptsipattern(opattern)
    result = Symbol[]
    for op in opattern
        isa(op, Symbol) && op∉result && push!(result, op)
    end
    return result
end

"""
    @subscripts expr::Expr -> Subscripts

Construct a subscript set from an opattern and optionally with a constrain.
"""
macro subscripts(expr::Expr)
    return subscriptsexpr(expr)
end
function subscriptsexpr(expr::Expr)
    if @capture(expr, [ops__;](cps__))
        @assert length(ops)==length(cps) "@subscripts error: dismatched opattern and constrain."
    elseif @capture(expr, op_(cp_))
        ops, cps = [op], [cp]
    elseif @capture(expr, [ops__;])
        cps = [true for i = 1:length(ops)]
    else
        @capture(expr, op_)
        ops, cps = [op], [true]
    end
    DS, RS, opattern, cpattern, mappings, constrains, fdefs = [], [], [], [], [], [], [], [], []
    for (op, cp) in zip(ops, cps)
        push!(DS, length(op.args))
        append!(opattern, op.args)
        ip = subscriptsipattern(op.args)
        push!(RS, length(ip))
        mapname = gensym("submap")
        op = Expr(:tuple, op.args...)
        push!(fdefs, :($mapname($(ip...)) = $op))
        push!(mappings, mapname)
        if isa(cp, Expr)
            push!(cpattern, String(Symbol(cp)))
            constrainname = gensym("subconstrain")
            push!(fdefs, :($constrainname($(ip...)) = $cp))
            push!(constrains, constrainname)
        else
            push!(cpattern, "nonconstrain")
            push!(constrains, :nonconstrain)
        end
    end
    DS, RS = Tuple(DS), Tuple(RS)
    opattern, cpattern = Tuple(opattern), Tuple(cpattern)
    mappings, constrains = Expr(:tuple, mappings...), Expr(:tuple, constrains...)
    return Expr(:block, fdefs..., :(Subscripts{$DS, $RS}($opattern, $cpattern, $mappings, $constrains)))
end

"""
    subscripts"..." -> Subscripts

Construct a subscript set from a literal string.
"""
macro subscripts_str(str)
    return subscriptsexpr(Meta.parse(str))
end

"""
    ==(subs₁::Subscripts, subs₂::Subscripts) -> Bool

Judge whether two subscript sets are equivalent to each other.
"""
@inline function Base.:(==)(subs₁::Subscripts{DS₁}, subs₂::Subscripts{DS₂}) where {DS₁, DS₂}
    return DS₁==DS₂ && subs₁.opattern==subs₂.opattern && subs₁.cpattern==subs₂.cpattern
end

"""
    isequal(subs₁::Subscripts, subs₂::Subscripts) -> Bool

Judge whether two subscript sets are equivalent to each other.
"""
@inline function Base.:isequal(subs₁::Subscripts{DS₁}, subs₂::Subscripts{DS₂}) where {DS₁, DS₂}
    return isequal(DS₁, DS₂) && isequal(subs₁.opattern, subs₂.opattern) && isequal(subs₁.cpattern, subs₂.cpattern)
end

"""
    dimension(subscripts::Subscripts) -> Int
    dimension(::Type{<:Subscripts}) -> Int

Get the total number of the whole variables of a subscript set.
"""
@inline dimension(subscripts::Subscripts) = dimension(typeof(subscripts))
@inline dimension(::Type{<:Subscripts{DS}}) where DS = sum(DS)

"""
    dimension(subscripts::Subscripts, i::Int) -> Int
    dimension(::Type{<:Subscripts}, i::Int) -> Int

Get the total number of the whole variables of the ith segment of a subscript set.
"""
@inline dimension(subscripts::Subscripts, i::Int) = dimension(typeof(subscripts), i)
@inline dimension(::Type{<:Subscripts{DS}}, i::Int) where DS = DS[i]

"""
    rank(subscripts::Subscripts) -> Int
    rank(::Type{<:Subscripts}) -> Int

Get the total number of the independent variables of a subscript set.
"""
@inline rank(subscripts::Subscripts) = rank(typeof(subscripts))
@inline rank(::Type{<:Subscripts{DS, RS} where DS}) where RS = sum(RS)

"""
    rank(subscripts::Subscripts, i::Int) -> Int
    rank(::Type{<:Subscripts}, i::Int) -> Int

Get the number of the independent variables of the ith segment of a subscript set.
"""
@inline rank(subscripts::Subscripts, i::Int) = rank(typeof(subscripts), i)
@inline rank(::Type{<:Subscripts{DS, RS} where DS}, i::Int) where RS = RS[i]

"""
    split(subscripts::Subscripts) -> Tuple{Vararg{Subscripts}}

Split a subscript set into individual independent segments.
"""
@generated function Base.split(subscripts::Subscripts{DS, RS}) where {DS, RS}
    exprs = []
    count = 0
    for i = 1:length(DS)
        NDS, NRS = (DS[i],), (RS[i],)
        opattern = Expr(:tuple, [:(subscripts[$j]) for j = (1+count):(DS[i]+count)]...)
        push!(exprs, :(Subscripts{$NDS, $NRS}($opattern, (subscripts.cpattern[$i],), (subscripts.mapping[$i],), (subscripts.constrain[$i],))))
        count += DS[i]
    end
    return Expr(:tuple, exprs...)
end

"""
    show(io::IO, subscripts::Subscripts)

Show a subscript set.
"""
function Base.show(io::IO, subscripts::Subscripts)
    segments = split(subscripts)
    @printf io "["
    for (i, segment) in enumerate(segments)
        @printf io "%s" join(segment.opattern, " ")
        i<length(segments) && @printf io "; "
    end
    @printf io "]"
    count = 1
    cpattern = String[]
    for segment in segments
        if segment.cpattern[1] == "nonconstrain"
            push!(cpattern, ":")
        else
            count += 1
            push!(cpattern, String(segment.cpattern[1]))
        end
    end
    count>1 && @printf io "(%s)" join(cpattern, ", ")
end

"""
    (subscripts::Subscripts)(values::Tuple{Vararg{Int}}) -> NTuple{dimension(subscripts), Int}

Construct a subscript set from an independent variable set.
"""
@generated function (subscripts::Subscripts)(values::Tuple{Vararg{Int}})
    @assert rank(subscripts)==fieldcount(values) "subscripts error: dismatched rank and input arguments."
    exprs = []
    count = 1
    for i = 1:fieldcount(fieldtype(subscripts, :mapping))
        vs = [:(values[$j]) for j = count:(count+rank(subscripts, i)-1)]
        push!(exprs, :(subscripts.mapping[$i]($(vs...))...))
        count += rank(subscripts, i)
    end
    return Expr(:tuple, exprs...)
end

"""
    isvalid(subscripts::Subscripts, values::Tuple{Vararg{Int}}) -> Bool

Judge whether an independent variable set are valid to construct a subscript set.
"""
@generated function Base.isvalid(subscripts::Subscripts, values::Tuple{Vararg{Int}})
    @assert rank(subscripts)==fieldcount(values) "subscripts error: dismatched rank and input arguments."
    rank(subscripts)==0 && return true
    count = rank(subscripts, 1)
    vs = [:(values[$j]) for j = 1:count]
    expr = :(subscripts.constrain[1]($(vs...)))
    for i = 2:fieldcount(fieldtype(subscripts, :constrain))
        vs = [:(values[$j]) for j = (count+1):(count+rank(subscripts, i))]
        expr = Expr(:(&&), expr, :(subscripts.constrain[$i]($(vs...))))
        count += rank(subscripts, i)
    end
    return expr
end

"""
    *(subs₁::Subscripts, subs₂::Subscripts) -> Subscripts

Get the multiplication between two subscript sets.
"""
@inline function Base.:*(subs₁::Subscripts{DS₁, RS₁}, subs₂::Subscripts{DS₂, RS₂}) where {DS₁, RS₁, DS₂, RS₂}
    DS, RS = (DS₁..., DS₂...), (RS₁..., RS₂...)
    opattern = (subs₁.opattern..., subs₂.opattern...)
    cpattern = (subs₁.cpattern..., subs₂.cpattern...)
    mapping = (subs₁.mapping..., subs₂.mapping...)
    constrain = (subs₁.constrain..., subs₂.constrain...)
    return Subscripts{DS, RS}(opattern, cpattern, mapping, constrain)
end

"""
    expand(subscripts::Subscripts, dimensions::Tuple{Vararg{Int}}) -> SbExpand

Expand a subscript set with the given variable ranges.
"""
function expand(subscripts::Subscripts, dimensions::Tuple{Vararg{Int}})
    @assert dimension(subscripts)==length(dimensions) "expand error: dismatched input dimensions $dimensions."
    dims, dcount, rcount = Vector{Int}(undef, rank(subscripts)), 0, 0
    for (i, segment) in enumerate(split(subscripts))
        opattern = segment.opattern
        ipattern = subscriptsipattern(opattern)
        for j = 1:length(opattern)
            isa(opattern[j], Int) && @assert 0<opattern[j]<=dimensions[dcount+j] "expand error: opattern($opattern) out of range."
        end
        for j = 1:length(ipattern)
            index, flag = 1, false
            while !isnothing((index = findnext(isequal(ipattern[j]), opattern, index); index))
                flag || (dims[rcount+j] = dimensions[dcount+index]; flag = true)
                flag && @assert dimensions[dcount+index]==dims[rcount+j] "expand error: dismatched input dimensions."
                index = index + 1
            end
        end
        dcount = dcount + length(opattern)
        rcount = rcount + length(ipattern)
    end
    return SbExpand(subscripts, NTuple{rank(subscripts), Int}(dims))
end
struct SbExpand{S<:Subscripts, D}
    subscripts::S
    candidates::CartesianIndices{D, NTuple{D, Base.OneTo{Int}}}
    function SbExpand(subscripts::Subscripts, dims::NTuple{D, Int}) where D
        @assert D==rank(subscripts) "SbExpand error: dismatched inputs."
        new{typeof(subscripts), D}(subscripts, CartesianIndices(dims))
    end
end
@inline Base.eltype(::Type{T}) where {T<:SbExpand} = NTuple{dimension(fieldtype(T, :subscripts)), Int}
@inline Base.IteratorSize(::Type{<:SbExpand}) = Base.SizeUnknown()
function Base.iterate(sbe::SbExpand, state::Int=1)
    while state <= length(sbe.candidates)
        values = Tuple(sbe.candidates[state])
        isvalid(sbe.subscripts, values) && return (sbe.subscripts(values), state+1)
        state = state + 1
    end
    return nothing
end

"""
    SubID{O<:Tuple, C<:Tuple, S<:Tuple{Vararg{Int}}} <: SimpleID

The id of a subscript set.
"""
struct SubID{O<:Tuple, C<:Tuple, S<:Tuple{Vararg{Int}}} <: SimpleID
    opattern::O
    cpattern::C
    segment::S
end
@inline SubID(subscripts::Subscripts) = id(subscripts)

"""
    id(subs::Subscripts) -> SubID

Get the id of a subscript set.
"""
@inline id(subscripts::Subscripts{DS}) where DS = SubID(subscripts.opattern, subscripts.cpattern, DS)

"""
    Coupling{V, I<:ID{SimpleID}} <: Element{V, I}

The coupling intra/inter interanl degrees of freedom at different lattice points.
"""
abstract type Coupling{V, I<:ID{SimpleID}} <: Element{V, I} end

"""
    couplingcenters(coupling::Coupling, bond::AbstractBond, info::Val) -> NTuple{rank(coupling), Int}

Get the centers of the coupling on a bond.
"""
@inline couplingcenters(coupling::Coupling, point::Point, ::Val) = ntuple(i->1, Val(rank(coupling)))

"""
    couplingpoints(coupling::Coupling, bond::AbstractBond, info::Val) -> NTuple{rank(coupling), eltype(bond)}

Get the points that correspond to where each order of the coupling acts on.
"""
@inline function couplingpoints(coupling::Coupling, bond::AbstractBond, info::Val)
    centers = couplingcenters(coupling, bond, info)
    return NTuple{rank(coupling), eltype(bond)}(bond[centers[i]] for i = 1:rank(coupling))
end

"""
    couplinginternals(coupling::Coupling, bond::AbstractBond, config::Config, info::Val) -> NTuple{rank(coupling), valtype(config)}

Get the interanl spaces that correspond to where each order of the coupling acts on.
"""
@inline function couplinginternals(coupling::Coupling, bond::AbstractBond, config::Config, info::Val)
    centers = couplingcenters(coupling, bond, info)
    return NTuple{rank(coupling), valtype(config)}(config[bond[centers[i]].pid] for i = 1:rank(coupling))
end

"""
    Couplings(cps::Coupling...)

A pack of couplings intra/inter interanl degrees of freedom at different lattice points.

Alias for `Elements{<:ID{SimpleID}, <:Coupling}`.
"""
const Couplings{I<:ID{SimpleID}, C<:Coupling} = Elements{I, C}
@inline Couplings(cps::Coupling...) = Elements(cps...)

"""
    @couplings cps -> Couplings

Convert an expression/literal to a set of couplings.
"""
macro couplings(cps)
    result = @eval(__module__, $cps)
    return isa(result, Coupling) ? Couplings(result) : isa(result, Couplings) ? result : error("@couplings error: inputs contain strangers that are not coupling/couplings.")
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
    TermCouplings(coupling::Coupling)
    TermCouplings(couplings::Couplings)
    TermCouplings(couplings::Function)

The function for the couplings of a term.
"""
struct TermCouplings{C<:Union{Function, Couplings}} <: TermFunction
    couplings::C
    TermCouplings(coupling::Coupling) = new{Couplings{coupling|>idtype, coupling|>typeof}}(Couplings(coupling))
    TermCouplings(couplings::Couplings) = new{couplings|>typeof}(couplings)
    TermCouplings(couplings::Function) = new{couplings|>typeof}(couplings)
end
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
    Term{K, R, I}(value, bondkind, couplings::TermCouplings, amplitude::TermAmplitude, modulate::TermModulate, factor) where {K, R, I}
    Term{K, R}(id::Symbol, value, bondkind;
        couplings::Union{Function, Coupling, Couplings, TermCouplings},
        amplitude::Union{Function, TermAmplitude, Nothing}=nothing,
        modulate::Union{Function, TermModulate, Bool}=false
        ) where {K, R}

A term of a quantum lattice system.
"""
mutable struct Term{K, R, I, V, B, C<:TermCouplings, A<:TermAmplitude, M<:TermModulate}
    value::V
    bondkind::B
    couplings::C
    amplitude::A
    modulate::M
    factor::V
    function Term{K, R, I}(value, bondkind, couplings::TermCouplings, amplitude::TermAmplitude, modulate::TermModulate, factor) where {K, R, I}
        @assert isa(K, Symbol) "Term error: kind must be a Symbol."
        @assert isa(I, Symbol) "Term error: id must be a Symbol."
        V, B, C, A, M = typeof(value), typeof(bondkind), typeof(couplings), typeof(amplitude), typeof(modulate)
        new{K, R, I, V, B, C, A, M}(value, bondkind, couplings, amplitude, modulate, factor)
    end
end
function Term{K, R}(id::Symbol, value, bondkind;
        couplings::Union{Function, Coupling, Couplings, TermCouplings},
        amplitude::Union{Function, TermAmplitude, Nothing}=nothing,
        modulate::Union{Function, TermModulate, Bool}=false
        ) where {K, R}
    isa(couplings, TermCouplings) || (couplings = TermCouplings(couplings))
    isa(amplitude, TermAmplitude) || (amplitude = TermAmplitude(amplitude))
    isa(modulate, TermModulate) || (modulate = TermModulate(id, modulate))
    Term{K, R, id}(value, bondkind, couplings, amplitude, modulate, 1)
end

"""
    kind(term::Term) -> Symbol
    kind(::Type{<:Term) -> Symbol

Get the kind of a term.
"""
@inline kind(term::Term) = kind(typeof(term))
@inline kind(::Type{<:Term{K}}) where K = K

"""
    rank(term::Term) -> Int
    rank(::Type{<:Term) -> Int

Get the rank of a term.
"""
@inline rank(term::Term) = rank(typeof(term))
@inline rank(::Type{<:Term{K, R} where K}) where R = R

"""
    id(term::Term) -> Symbol
    id(::Type{<:Term) -> Symbol

Get the id of a term.
"""
@inline id(term::Term) = id(typeof(term))
@inline id(::Type{<:Term{K, R, I} where {K, R}}) where I = I

"""
    valtype(term::Term)
    valtype(::Type{<:Term)

Get the value type of a term.
"""
@inline Base.valtype(term::Term) = valtype(typeof(term))
@inline Base.valtype(::Type{<:Term{K, R, I, V} where {K, R, I}}) where V = V

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
@inline ismodulatable(::Type{<:Term{K, R, I, V, B, <:TermCouplings, <:TermAmplitude, M} where {K, R, I, V, B}}) where M = ismodulatable(M)

"""
    isHermitian(term::Term) -> Bool
    isHermitian(::Type{<:Term}) -> Bool
"""
@inline isHermitian(term::Term) = isHermitian(typeof(term))
@inline isHermitian(::Type{<:Term}) = error("isHermitian error: not implemented.")

"""
    ==(term1::Term, term2::Term) -> Bool

Judge whether two terms are equivalent to each other.
"""
@inline Base.:(==)(term1::Term, term2::Term) = ==(efficientoperations, term1, term2)

"""
    isequal(term1::Term, term2::Term) -> Bool

Judge whether two terms are equivalent to each other.
"""
@inline Base.isequal(term1::Term, term2::Term) = isequal(efficientoperations, term1, term2)

"""
    show(io::IO, term::Term)

Show a term.
"""
function Base.show(io::IO, term::Term)
    @printf io "%s{%s}(id=%s, value=%s, bondkind=%s, factor=%s)" kind(term) rank(term) id(term) decimaltostr(term.value) term.bondkind decimaltostr(term.factor)
end

"""
    repr(term::Term, bond::AbstractBond, config::Config) -> String

Get the repr representation of a term on a bond with a given config.
"""
function Base.repr(term::Term, bond::AbstractBond, config::Config)
    cache = String[]
    if term.bondkind == bond|>kind
        value = term.value * term.amplitude(bond) * term.factor
        if abs(value) ≠ 0
            for coupling in values(term.couplings(bond))
                length(expand(coupling, bond, config, term|>kind|>Val))>0 &&  push!(cache, @sprintf "%s: %s" abbr(term) repr(value*coupling))
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
    return :(Term{kind(term), rank(term), id(term)}($(exprs...)))
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
@inline Base.:/(term::Term, factor) = term * (1/factor)

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
    otype(T::Type{<:Term}, C::Type{<:Config}, B::Type{<:AbstractBond})

Get the compatible operator type from the type of a term, a configuration of the internal degrees of freedom and a bond.
"""
otype(T::Type{<:Term}, C::Type{<:Config}, B::Type{<:AbstractBond}) = Operator{T|>valtype, ID{oidtype(C|>valtype, B|>eltype, Val(:info)), T|>rank}}

"""
    expand!(operators::Operators, term::Term, bond::AbstractBond, config::Config, half::Bool=false; table::Union{Nothing, Table}=nothing) -> Operators
    expand!(operators::Operators, term::Term, bonds::Bonds, config::Config, half::Bool=false; table::Union{Nothing, Table}=nothing) -> Operators

Expand the operators of a term on a bond/set-of-bonds with a given config.

The `half` parameter determines the behavior of generating operators, which falls into the following two categories
* `true`: "Hermitian half" of the generated operators
* `false`: "Hermitian whole" of the generated operators
"""
function expand!(operators::Operators, term::Term, bond::AbstractBond, config::Config, half::Bool=false; table::Union{Nothing, Table}=nothing)
    if term.bondkind == bond|>kind
        value = term.value * term.amplitude(bond) * term.factor
        if abs(value) ≠ 0
            hermitian = isHermitian(term)
            optype = otype(term|>typeof, config|>typeof, bond|>typeof)
            record = (isnothing(hermitian) && length(operators)>0) ? Set{optype|>idtype}() : nothing
            for coupling in values(term.couplings(bond))
                for (coeff, id) in expand(coupling, bond, config, term|>kind|>Val)
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
@generated function expand!(operators::Operators, term::Term, bonds::Bonds, config::Config, half::Bool=false; table::Union{Nothing, Table}=nothing)
    exprs = []
    for i = 1:rank(bonds)
        push!(exprs, :(for bond in bonds.bonds[$i] expand!(operators, term, bond, config, half; table=table) end))
    end
    push!(exprs, :(return operators))
    return Expr(:block, exprs...)
end
@inline termfactor(id::ID{AbstractOID}, ::Val) = isHermitian(id) ? 2 : 1

"""
    expand(term::Term, bond::AbstractBond, config::Config, half::Bool=false; table::Union{Nothing, Table}=nothing) -> Operators
    expand(term::Term, bonds::Bonds, config::Config, half::Bool=false; table::Union{Nothing, Table}=nothing) -> Operators

Expand the operators of a term on a bond/set-of-bonds with a given config.
"""
@inline function expand(term::Term, bond::AbstractBond, config::Config, half::Bool=false; table::Union{Nothing, Table}=nothing)
    optype = otype(term|>typeof, config|>typeof, bond|>typeof)
    expand!(Operators{idtype(optype), optype}(), term, bond, config, half; table=table)
end
@inline function expand(term::Term, bonds::Bonds, config::Config, half::Bool=false; table::Union{Nothing, Table}=nothing)
    optype = otype(term|>typeof, config|>typeof, bonds|>eltype)
    expand!(Operators{idtype(optype), optype}(), term, bonds, config, half, table=table)
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
    GenOperators(terms::Tuple{Vararg{Term}}, bonds::Bonds, config::Config, half::Bool; table::Union{Nothing, Table}=nothing)

A set of operators. This is the core of [`AbstractGenerator`](@ref).
"""
struct GenOperators{C<:Operators, A<:NamedContainer{Operators}, B<:NamedContainer{Operators}}
    constops::C
    alterops::A
    boundops::B
end
@generated function GenOperators(terms::Tuple{Vararg{Term}}, bonds::Bonds, config::Config, half::Bool; table::Union{Nothing, Table}=nothing)
    constterms, alterterms = [], []
    for term in fieldtypes(terms)
        ismodulatable(term) ? push!(alterterms, term) : push!(constterms, term)
    end
    names = NTuple{fieldcount(terms), Symbol}(id(term) for term in fieldtypes(terms))
    alternames = NTuple{length(alterterms), Symbol}(id(term) for term in alterterms)
    exprs, alterops, boundops = [], [], []
    push!(exprs, quote
        choosedterms = length($constterms)>0 ? $constterms : $alterterms
        constoptp = otype(choosedterms[1], config|>typeof, bonds|>eltype)
        constidtp = constoptp |> idtype
        for i = 2:length(choosedterms)
            tempoptp = otype(choosedterms[i], config|>typeof, bonds|>eltype)
            constoptp = promote_type(constoptp, tempoptp)
            constidtp = promote_type(constidtp, tempoptp|>idtype)
        end
        constops = Operators{constidtp, constoptp}()
        innerbonds = filter(acrossbonds, bonds, Val(:exclude))
        boundbonds = filter(acrossbonds, bonds, Val(:include))
    end)
    for i = 1:fieldcount(terms)
        push!(boundops, :(expand(one(terms[$i]), boundbonds, config, half, table=table)))
        if ismodulatable(fieldtype(terms, i))
            push!(alterops, :(expand(one(terms[$i]), innerbonds, config, half, table=table)))
        else
            push!(exprs, :(expand!(constops, terms[$i], innerbonds, config, half, table=table)))
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
    ==(genops1::GenOperators, genops2::GenOperators) -> Bool

Judge whether two sets of operators are equivalent to each other.
"""
@inline Base.:(==)(genops1::GenOperators, genops2::GenOperators) = ==(efficientoperations, genops1, genops2)

"""
    isequal(genops1::GenOperators, genops2::GenOperators) -> Bool

Judge whether two sets of operators are equivalent to each other.
"""
@inline Base.isequal(genops1::GenOperators, genops2::GenOperators) = isequal(efficientoperations, genops1, genops2)

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
    reset!(genops::GenOperators, terms::Tuple{Vararg{Term}}, bonds::Bonds, config::Config, half::Bool; table::Union{Nothing, Table}=nothing) -> GenOperators

Reset a set of operators by new terms, bonds, config, etc..
"""
@generated function reset!(genops::GenOperators, terms::Tuple{Vararg{Term}}, bonds::Bonds, config::Config, half::Bool; table::Union{Nothing, Table}=nothing)
    exprs = []
    push!(exprs, quote
        empty!(genops)
        innerbonds = filter(acrossbonds, bonds, Val(:exclude))
        boundbonds = filter(acrossbonds, bonds, Val(:include))
    end)
    for (i, term) in enumerate(fieldtypes(terms))
        name = QuoteNode(term|>id)
        push!(exprs, :(expand!(getfield(genops.boundops, $name), one(terms[$i]), boundbonds, config, half, table=table)))
        if ismodulatable(term)
            push!(exprs, :(expand!(getfield(genops.alterops, $name), one(terms[$i]), innerbonds, config, half, table=table)))
        else
            push!(exprs, :(expand!(genops.constops, terms[$i], innerbonds, config, half, table=table)))
        end
    end
    push!(exprs, :(return genops))
    return Expr(:block, exprs...)
end

"""
    AbstractGenerator{TS<:NamedContainer{Term}, BS<:Bonds, C<:Config, T<:Table, B<:Boundary, OS<:GenOperators}

Abstract generator.

By protocol, a concrete generator should have the following predefined contents:
* `terms::TS`: the terms contained in a generator
* `bonds::BS`: the bonds on which the terms are defined
* `config::C`: the configuration of the interanl degrees of freedom
* `half::Bool`: true for generating an Hermitian half of the operators and false for generating the whole
* `table::Table`: the index-sequence table
* `boundary::B`: boundary twist for the generated operators, `nothing` for no twist
* `operators::OS`: the generated operators
"""
abstract type AbstractGenerator{TS<:NamedContainer{Term}, BS<:Bonds, C<:Config, T<:Table, B<:Boundary, OS<:GenOperators} end
@inline contentnames(::Type{<:AbstractGenerator}) = (:terms, :bonds, :config, :half, :table, :boundary, :operators)

"""
    ==(gen₁::AbstractGenerator, gen₂::AbstractGenerator) -> Bool

Judge whether generators are equivalent to each other.
"""
@inline Base.:(==)(gen₁::AbstractGenerator, gen₂::AbstractGenerator) = ==(efficientoperations, gen₁, gen₂)

"""
    isequal(gen₁::AbstractGenerator, gen₂::AbstractGenerator) -> Bool

Judge whether generators are equivalent to each other.
"""
@inline Base.isequal(gen₁::AbstractGenerator, gen₂::AbstractGenerator) = isequal(efficientoperations, gen₁, gen₂)

"""
    otype(gen::AbstractGenerator)
    otype(::Type{<:AbstractGenerator})

Get the operator type of the generated opeators.
"""
@inline otype(gen::AbstractGenerator) = otype(typeof(gen))
@inline otype(::Type{<:AbstractGenerator{<:NamedContainer{Term}, <:Bonds, <:Config, <:Table, <:Boundary, OS}}) where {OS<:GenOperators} = eltype(OS)

"""
    Parameters(gen::AbstractGenerator) -> Parameters

Get the parameters of the terms of a generator.
"""
@generated function Parameters(gen::AbstractGenerator{TS}) where {TS<:NamedContainer{Term}}
    names, values = fieldnames(TS), [:(getcontent(gen, :terms)[$i].value) for i = 1:fieldcount(TS)]
    return :(Parameters{$names}($(values...)))
end

"""
    expand!(operators::Operators, gen::AbstractGenerator) -> Operators

Expand the operators of a generator.
"""
@inline expand!(operators::Operators, gen::AbstractGenerator) = expand!(operators, getcontent(gen, :operators), getcontent(gen, :boundary); Parameters(gen)...)

"""
    expand(gen::AbstractGenerator) -> Operators
    expand(gen::AbstractGenerator, name::Symbol) -> Operators
    expand(gen::AbstractGenerator, i::Int) -> Operators
    expand(gen::AbstractGenerator, name::Symbol, i::Int) -> Operators

Expand the operators of a generator:
1) the total operators;
2) the operators of a specific term;
3) the operators on a specific bond;
4) the operators of a specific term on a specific bond.
"""
@inline expand(gen::AbstractGenerator) = expand!(Operators{idtype(otype(gen)), otype(gen)}(), gen)
function expand(gen::AbstractGenerator, name::Symbol)
    bonds = getcontent(gen, :bonds)
    config = getcontent(gen, :config)
    ops = getcontent(gen, :operators)
    term = getfield(getcontent(gen, :terms), name)
    optp = otype(term|>typeof, config|>typeof, bonds|>eltype)
    result = Operators{idtype(optp), optp}()
    if ismodulatable(term)
        for opt in getfield(ops.alterops, name)|>values add!(result, opt*term.value) end
    else
        expand!(result, term, filter(acrossbonds, bonds, Val(:exclude)), config, getcontent(gen, :half), table=getcontent(gen, :table))
    end
    for opt in getfield(ops.boundops, name)|>values
        add!(result, getcontent(gen, :boundary)(opt)*term.value)
    end
    return result
end
@generated function expand(gen::AbstractGenerator{TS}, i::Int) where {TS<:NamedContainer{Term}}
    exprs = []
    push!(exprs, quote
        bond = getcontent(gen, :bonds)[i]
        result = Operators{idtype(otype(gen)), otype(gen)}()
    end)
    for i = 1:fieldcount(TS)
        push!(exprs, :(expand!(result, getcontent(gen, :terms)[$i], bond, getcontent(gen, :config), getcontent(gen, :half), table=getcontent(gen, :table))))
    end
    push!(exprs, quote
        isintracell(bond) && for opt in values(result)
            result[id(opt)] = getcontent(gen, :boundary)(opt)
        end
        return result
    end)
    return Expr(:block, exprs...)
end
function expand(gen::AbstractGenerator, name::Symbol, i::Int)
    bond = getcontent(gen, :bonds)[i]
    result = expand(getfield(getcontent(gen, :terms), name), bond, getcontent(gen, :config), getcontent(gen, :half), table=getcontent(gen, :table))
    isintracell(bond) && for opt in values(result)
        result[id(opt)] = getcontent(gen, :boundary)(opt)
    end
    return result
end

"""
    update!(gen::AbstractGenerator; kwargs...) -> typeof(gen)

Update the coefficients of the terms in a generator.
"""
@generated function update!(gen::AbstractGenerator{TS}; kwargs...) where {TS<:NamedContainer{Term}}
    exprs = [:(update!(getcontent(gen, :boundary); kwargs...))]
    for i = 1:fieldcount(TS)
        ismodulatable(fieldtype(TS, i)) && push!(exprs, :(update!(getcontent(gen, :terms)[$i]; kwargs...)))
    end
    push!(exprs, :(return gen))
    return Expr(:block, exprs...)
end

"""
    Generator(terms::Tuple{Vararg{Term}}, bonds::Bonds, config::Config;
        half::Bool=false, table::Union{Table,Nothing}=nothing, boundary::Boundary=plain
        )

A generator of operators based on terms, configuration of internal degrees of freedom, and boundary twist.
"""
struct Generator{TS<:NamedContainer{Term}, BS<:Bonds, C<:Config, T<:Table, B<:Boundary, OS<:GenOperators} <: AbstractGenerator{TS, BS, C, T, B, OS}
    terms::TS
    bonds::BS
    config::C
    half::Bool
    table::T
    boundary::B
    operators::OS
end
@inline function Generator(terms::Tuple{Vararg{Term}}, bonds::Bonds, config::Config;
        half::Bool=false, table::Union{Table,Nothing}=nothing, boundary::Boundary=plain
        )
    isnothing(table) && (table = Table(config))
    Generator(namedterms(terms), bonds, config, half, table, boundary, GenOperators(terms, bonds, config, half, table=table))
end
@generated function namedterms(terms::Tuple{Vararg{Term}})
    names = NTuple{fieldcount(terms), Symbol}(id(fieldtype(terms, i)) for i = 1:fieldcount(terms))
    return :(NamedContainer{$names}(terms))
end

"""
    empty!(gen::Generator) -> Generator

Empty the :bonds, :config, :table and :operators of a generator.
"""
function Base.empty!(gen::Generator)
    empty!(gen.bonds)
    empty!(gen.config)
    empty!(gen.table)
    empty!(gen.operators)
    return gen
end

"""
    empty(gen::Generator) -> Generator

Get an empty copy of a generator.
"""
@inline Base.empty(gen::Generator) = Generator(gen.terms, empty(gen.bonds), empty(gen.config), gen.half, empty(gen.table), gen.boundary, empty(gen.operators))

"""
    reset!(gen::Generator, lattice::AbstractLattice) -> Generator

Reset a generator by a new lattice.
"""
function reset!(gen::Generator, lattice::AbstractLattice)
    reset!(gen.bonds, lattice)
    reset!(gen.config, lattice.pids)
    reset!(gen.table, gen.config)
    reset!(gen.operators, Tuple(gen.terms), gen.bonds, gen.config, gen.half, table=gen.table)
    return gen
end

end  # module
