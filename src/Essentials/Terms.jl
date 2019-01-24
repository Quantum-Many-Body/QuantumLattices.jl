module Terms

using Printf: @printf,@sprintf
using StaticArrays: SVector
using ..Spatials: AbstractBond,neighbor,pidtype
using ..DegreesOfFreedom: Index,IDFConfig,Couplings
using ...Prerequisites: Float,atol
using ...Prerequisites.TypeTraits: efficientoperations
using ...Mathematics.AlgebraOverFields: SimpleID,ID,Element,Elements

import ...Prerequisites.Interfaces: rank,expand

export rank
export OID,Operator,Operators,isHermitian
export Term,statistics,species,abbr

@generated function propercoord(vector::SVector{N,Float}) where N
    exprs=[:(vector[$i]===-0.0 ? 0.0 : vector[$i]) for i=1:N]
    return :(SVector($(exprs...)))
end

"""
    OID(index::Index,rcoord::Union{Nothing,SVector{N,Float}},icoord::Union{Nothing,SVector{N,Float}},seq::Union{Nothing,Int}) where N
    OID(index::Index;rcoord::Union{Nothing,SVector}=nothing,icoord::Union{Nothing,SVector}=nothing,seq::Union{Nothing,Int}=nothing)

Operator id.
"""
struct OID{I<:Index,RC<:Union{Nothing,SVector},IC<:Union{Nothing,SVector},S<:Union{Nothing,Int}} <: SimpleID
    index::I
    rcoord::RC
    icoord::IC
    seq::S
    function OID(index::Index,rcoord::Union{Nothing,SVector{N,Float}},icoord::Union{Nothing,SVector{N,Float}},seq::Union{Nothing,Int}) where N
        isa(rcoord,SVector) && (rcoord=propercoord(rcoord))
        isa(icoord,SVector) && (icoord=propercoord(icoord))
        new{typeof(index),typeof(rcoord),typeof(icoord),typeof(seq)}(index,rcoord,icoord,seq)
    end
end
OID(index::Index;rcoord::Union{Nothing,SVector}=nothing,icoord::Union{Nothing,SVector}=nothing,seq::Union{Nothing,Int}=nothing)=OID(index,rcoord,icoord,seq)
Base.propertynames(::Type{<:ID{N,<:OID}},private::Bool=false) where N=private ? (:contents,:indexes,:rcoords,:icoords) : (:indexes,:rcoords,:icoords)

"""
    adjoint(oid::OID) -> OID

Get the adjoint of an operator id.
"""
Base.adjoint(oid::OID)=OID(oid.index',oid.rcoord,oid.icoord,oid.seq)

"""
    Operator{N,V<:Number,I<:ID} <: Element{V,I}

Abstract type for an operator.
"""
abstract type Operator{N,V<:Number,I<:ID} <: Element{V,I} end
function (O::Type{<:Operator})( value::Number,
                                indexes::NTuple{N,Index};
                                rcoords::Union{Nothing,NTuple{N,SVector{M,Float}}}=nothing,
                                icoords::Union{Nothing,NTuple{N,SVector{M,Float}}}=nothing,
                                seqs::Union{Nothing,NTuple{N,Int}}=nothing
                                ) where {N,M}
    rcoords===nothing && (rcoords=ntuple(i->nothing,N))
    icoords===nothing && (icoords=ntuple(i->nothing,N))
    seqs===nothing && (seqs=ntuple(i->nothing,N))
    return O(value,ID(OID,indexes,rcoords,icoords,seqs))
end

"""
    adjoint(opt::Operator{N}) where N -> Operator

Get the adjoint of an operator.
"""
@generated function Base.adjoint(opt::Operator{N}) where N
    exprs=[:(opt.id[$i]') for i=N:-1:1]
    return :(typeof(opt).name.wrapper(opt.value',ID($(exprs...))))
end

"""
    isHermitian(opt::Operator) -> Bool

Judge whether an operator is Hermitian.
"""
isHermitian(opt::Operator)=opt==opt'

"""
    Operators(opts::Operator...)

A set of operators.

Type alias of `Operators{I<:ID,O<:Operator}=Elements{I,O}`.
"""
const Operators{I<:ID,O<:Operator}=Elements{I,O}
Operators(opts::Operator...)=Elements(opts...)
function Base.adjoint(opts::Operators)
    result=Operators{opts|>keytype,opts|>valtype}()
    for opt in values(opts)
        nopt=opt|>adjoint
        result[nopt.id]=nopt
    end
    return result
end

"""
    isHermitian(opts::Operators) -> Bool

Judge whether a set of operators as a whole is Hermitian.
"""
isHermitian(opts::Operators)=opts==opts'

abstract type TermFunction <: Function end
Base.:(==)(tf1::TermFunction,tf2::TermFunction) = ==(efficientoperations,tf1,tf2)
Base.isequal(tf1::TermFunction,tf2::TermFunction)=isequal(efficientoperations,tf1,tf2)
struct TermCouplings{C<:Couplings} <: TermFunction
    couplings::C
end
(termcouplings::TermCouplings)(args...;kwargs...)=termcouplings.couplings
struct TermAmplitude <: TermFunction end
(termamplitude::TermAmplitude)(args...;kwargs...)=1
const termamplitude=TermAmplitude()
struct TermModulate <: TermFunction
    id::Symbol
end
(termmodulate::TermModulate)(;kwargs...)=get(kwargs,termmodulate.id,nothing)

"""
    Term{ST,SP,RK}(id::Symbol,value::Number,neighbor::Any,couplings::Function,amplitude::Function,modulate::Union{Function,Nothing},factor::Number) where {ST,SP,RK}
    Term{ST,SP,RK}( id::Symbol,value::Number,neighbor::Any;
                    couplings::Union{Function,Couplings},
                    amplitude::Union{Function,Nothing}=nothing,
                    modulate::Union{Function,Bool}=false,
                    factor::Number=1
                    ) where {ST,SP,RK}

A term of a quantum lattice system.
"""
struct Term{Statistics,Species,Rank,V<:Number,N<:Any,C<:Function,A<:Function,M<:Union{Function,Nothing},F<:Number}
    id::Symbol
    value::V
    neighbor::N
    couplings::C
    amplitude::A
    modulate::M
    factor::F
    function Term{ST,SP,RK}(id::Symbol,value::Number,neighbor::Any,couplings::Function,amplitude::Function,modulate::Union{Function,Nothing},factor::Number) where {ST,SP,RK}
        @assert ST in ('F','B') && isa(SP,Symbol) && isa(RK,Int) "Term error: not supported type parameter."
        new{ST,SP,RK,typeof(value),typeof(neighbor),typeof(couplings),typeof(amplitude),typeof(modulate),typeof(factor)}(id,value,neighbor,couplings,amplitude,modulate,factor)
    end
end
function Term{ST,SP,RK}(id::Symbol,value::Number,neighbor::Any;
                        couplings::Union{Function,Couplings},
                        amplitude::Union{Function,Nothing}=nothing,
                        modulate::Union{Function,Bool}=false,
                        factor::Number=1
                        ) where {ST,SP,RK}
    isa(couplings,Couplings) && (couplings=TermCouplings(couplings))
    amplitude===nothing && (amplitude=termamplitude)
    isa(modulate,Bool) && (modulate=modulate ? TermModulate(id) : nothing)
    Term{ST,SP,RK}(id,value,neighbor,couplings,amplitude,modulate,factor)
end

"""
    statistics(term::Term) -> Char
    statistics(::Type{<:Term{ST}}) where ST -> Char

Get the statistics of a term.
"""
statistics(term::Term)=term|>typeof|>statistics
statistics(::Type{<:Term{ST}}) where ST=ST

"""
    species(term::Term) -> Symbol
    species(::Type{<:Term{ST,SP}}) where {ST,SP} -> Symbol

Get the species of a term.
"""
species(term::Term)=term|>typeof|>species
species(::Type{<:Term{ST,SP}}) where {ST,SP}=SP

"""
    rank(term::Term) -> Int
    rank(::Type{<:Term{ST,SP,RK}}) where {ST,SP,RK}

Get the rank of a term.
"""
rank(term::Term)=term|>typeof|>rank
rank(::Type{<:Term{ST,SP,RK}}) where {ST,SP,RK}=RK

"""
    abbr(::Term) -> Symbol

Get the abbreviation of the species of a term.
"""
abbr(::Term)=:tm

"""
    ==(term1::Term,term2::Term) -> Bool
    isequal(term1::Term,term2::Term) -> Bool

Judge whether two terms are equivalent to each other.
"""
Base.:(==)(term1::Term,term2::Term) = ==(efficientoperations,term1,term2)
Base.isequal(term1::Term,term2::Term)=isequal(efficientoperations,term1,term2)

"""
    show(io::IO,term::Term)

Show a term.
"""
Base.show(io::IO,term::Term)=@printf io "%s{%s,%s}(id=%s,value=%s,neighbor=%s)" species(term) statistics(term) rank(term) term.id term.value term.neighbor

"""
    repr(term::Term,bond::AbstractBond,config::IDFConfig) -> String

Get the repr representation of a term on a bond with a given config.
"""
function Base.repr(term::Term,bond::AbstractBond,config::IDFConfig)
    cache=String[]
    if term.neighbor==bond|>neighbor
        pids=NTuple{length(bond),pidtype(bond)}(point.pid for point in bond)
        interanls=NTuple{length(bond),valtype(config)}(config[pid] for pid in pids)
        value=term.value*term.amplitude(bond)*term.factor
        if !isapprox(abs(value),0.0,atol=atol)
            for coupling in values(term.couplings(bond))
                length(expand(coupling,pids,interanls))>0 &&  push!(cache,@sprintf "%s: %s" abbr(term) repr(term.value*coupling))
            end
        end
    end
    return join(cache,"\n")
end

"""
    replace(term::Term{ST,SP,RK};kwargs...) where {ST,SP,RK} -> Term

Replace some attributes of a term with key word arguments.
"""
@generated function Base.replace(term::Term{ST,SP,RK};kwargs...) where {ST,SP,RK}
    exprs=[:(get(kwargs,$name,getfield(term,$name))) for name in QuoteNode.(term|>fieldnames)]
    return :(Term{ST,SP,RK}($(exprs...)))
end

"""
    +(term::Term) -> Term
    -(term::Term) -> Term
    *(term::Term,factor::Number) -> Term
    *(factor::Number,term::Term) -> Term
    /(term::Term,factor::Number) -> Term

Allowed arithmetic operations for a term.
"""
Base.:+(term::Term)=term
Base.:-(term::Term)=term*(-1)
Base.:*(term::Term,factor::Number)=factor*term
Base.:*(factor::Number,term::Term)=replace(term,factor=factor*term.factor)
Base.:/(term::Term,factor::Number)=term*(1/factor)

"""
    one(term::Term) -> Term

Get a unit term.
"""
Base.one(term::Term)=replace(term,value=one(term.value))

"""
    zero(term::Term) -> Term

Get a zero term.
"""
Base.zero(term::Term)=replace(term,value=zero(term.value))

end  # module
