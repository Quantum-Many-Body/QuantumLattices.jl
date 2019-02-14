module Terms

using Printf: @printf,@sprintf
using StaticArrays: SVector
using ..Spatials: AbstractBond,neighbor,pidtype
using ..DegreesOfFreedom: Index,IDFConfig,Couplings,Table,propercenters
using ...Prerequisites: Float,atol,decimaltostr
using ...Prerequisites.TypeTraits: efficientoperations
using ...Mathematics.AlgebraOverFields: SimpleID,ID,Element,Elements,idtype

import ...Prerequisites.Interfaces: rank,expand,dimension,add!,update!
import ..Spatials: rcoord,icoord

export rank,expand,rcood,icoord,update!
export OID,Operator,Operators,isHermitian,oidtype,otype
export TermFunction,TermAmplitude,TermCouplings,TermModulate
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
Base.hash(oid::OID,h::UInt)=hash(oid|>values,h)
Base.fieldnames(::Type{<:OID})=(:index,:rcoord,:icoord,:seq)
Base.propertynames(::Type{<:ID{<:NTuple{N,OID}}},private::Bool=false) where N=private ? (:contents,:indexes,:rcoords,:icoords) : (:indexes,:rcoords,:icoords)

"""
    show(io::IO,oid::OID)

Show an operator id.
"""
function Base.show(io::IO,oid::OID)
    @printf io "OID(%s" oid.index
    oid.rcoord===nothing ? (@printf io ",:") : (@printf io ",[%s]" join(oid.rcoord,","))
    oid.icoord===nothing ? (@printf io ",:") : (@printf io ",[%s]" join(oid.icoord,","))
    @printf io ",%s)" oid.seq===nothing ? ":" : oid.seq
end

"""
    adjoint(oid::OID) -> typeof(oid)
    adjoint(oid::ID{<:NTuple{N,OID}}) where N -> typeof(oid)

Get the adjoint of an operator id.
"""
Base.adjoint(oid::OID)=OID(oid.index',oid.rcoord,oid.icoord,oid.seq)
@generated Base.adjoint(oid::ID{<:NTuple{N,OID}}) where N=Expr(:call,:ID,[:(oid[$i]') for i=N:-1:1]...)

"""
    isHermitian(oid::ID{<:NTuple{N,OID}}) where N -> Bool

Judge whether an operator id is Hermitian.
"""
function isHermitian(oid::ID{<:NTuple{N,OID}}) where N
    for i=1:((N+1)÷2)
        oid[i]'==oid[N+1-i] || return false
    end
    return true
end

"""
    oidtype

Get the compatible oid type from a bond type and a table type.
"""
function oidtype end

"""
    Operator{N,V<:Number,I<:ID{<:NTuple{N,OID}}} <: Element{N,V,I}

Abstract type for an operator.
"""
abstract type Operator{N,V<:Number,I<:ID{<:NTuple{N,OID}}} <: Element{N,V,I} end
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
    show(io::IO,opt::Operator)

Show an operator.
"""
Base.show(io::IO,opt::Operator)=@printf io "%s(value=%s,id=%s)" nameof(typeof(opt)) decimaltostr(opt.value) opt.id

"""
    adjoint(opt::Operator{N}) where N -> Operator

Get the adjoint of an operator.
"""
Base.adjoint(opt::Operator)=typeof(opt).name.wrapper(opt.value',opt.id')

"""
    isHermitian(opt::Operator) -> Bool

Judge whether an operator is Hermitian.
"""
isHermitian(opt::Operator)=isa(opt.value,Real) && isHermitian(opt.id)

"""
    rcoord(opt::Operator{1}) -> SVector
    rcoord(opt::Operator{2}) -> SVector

Get the whole rcoord of an operator.
"""
rcoord(opt::Operator{1})=opt.id[1].rcoord
rcoord(opt::Operator{2})=opt.id[1].rcoord-opt.id[2].rcoord

"""
    icoord(opt::Operator{1}) -> SVector
    icoord(opt::Operator{2}) -> SVector

Get the whole icoord of an operator.
"""
icoord(opt::Operator{1})=opt.id[1].icoord
icoord(opt::Operator{2})=opt.id[1].icoord-opt.id[2].icoord

"""
    otype

Get the compatible operator type from a term type, a bond type and a table type.
"""
function otype end

"""
    Operators(opts::Operator...)

A set of operators.

Type alias of `Operators{I<:ID,O<:Operator}=Elements{I,O}`.
"""
const Operators{I<:ID,O<:Operator}=Elements{I,O}
Operators(opts::Operator...)=Elements(opts...)

"""
    adjoint(opts::Operators) -> Operators

Get the adjoint of a set of operators.
"""
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

"""
    TermFunction <: Function

Abstract type for concrete term functions.
"""
abstract type TermFunction <: Function end

"""
    ==(tf1::TermFunction,tf2::TermFunction) -> Bool
    isequal(tf1::TermFunction,tf2::TermFunction) -> Bool

Judge whether two concrete term functions are equivalent to each other.
"""
Base.:(==)(tf1::TermFunction,tf2::TermFunction) = ==(efficientoperations,tf1,tf2)
Base.isequal(tf1::TermFunction,tf2::TermFunction)=isequal(efficientoperations,tf1,tf2)

"""
    TermAmplitude <: TermFunction

The function for the amplitude of a term.
"""
struct TermAmplitude <: TermFunction end
(termamplitude::TermAmplitude)(args...;kwargs...)=1
const termamplitude=TermAmplitude()

"""
    TermCouplings(candidate::Couplings)
    TermCouplings(candidates::NTuple{N,<:Couplings},choice::Function) where N

The function for the couplings of a term.
"""
struct TermCouplings{T<:Tuple,C<:Union{Function,Nothing}} <: TermFunction
    candidates::T
    choice::C
    TermCouplings(candidate::Couplings)=new{Tuple{typeof(candidate)},Nothing}((candidate,),nothing)
    TermCouplings(candidates::NTuple{N,<:Couplings},choice::Function) where N=new{typeof(candidates),typeof(choice)}(candidates,choice)
end
(termcouplings::TermCouplings{<:Tuple,Nothing})(args...;kwargs...)=termcouplings.candidates[1]
(termcouplings::TermCouplings{<:Tuple,<:Function})(args...;kwargs...)=termcouplings.candidates[termcouplings.choice(args...;kwargs...)]

"""
    rank(tcs::TermCouplings) -> Int
    rank(TCS::Type{<:TermCouplings}) -> Int

Get the rank of the couplings it contained.
"""
rank(tcs::TermCouplings)=tcs|>typeof|>rank
rank(TCS::Type{<:TermCouplings})=fieldtype(TCS,:candidates)|>eltype|>valtype|>rank

"""
    TermModulate(id::Symbol)

The function for the modulation of a term.
"""
struct TermModulate <: TermFunction
    id::Symbol
end
(termmodulate::TermModulate)(;kwargs...)=get(kwargs,termmodulate.id,nothing)

"""
    Term{ST,SP}(id::Symbol,value::Number,neighbor::Any,couplings::TermCouplings,amplitude::Function,modulate::Union{Function,Nothing},factor::Number) where {ST,SP}
    Term{ST,SP}(id::Symbol,value::Number,neighbor::Any;
                couplings::Union{TermCouplings,Couplings},
                amplitude::Union{Function,Nothing}=nothing,
                modulate::Union{Function,Bool}=false,
                factor::Number=1
                ) where {ST,SP}

A term of a quantum lattice system.
"""
mutable struct Term{Statistics,Species,V<:Number,N<:Any,C<:TermCouplings,A<:Function,M<:Union{Function,Nothing}}
    id::Symbol
    value::V
    neighbor::N
    couplings::C
    amplitude::A
    modulate::M
    factor::V
    function Term{ST,SP}(id::Symbol,value::Number,neighbor::Any,couplings::TermCouplings,amplitude::Function,modulate::Union{Function,Nothing},factor::Number) where {ST,SP}
        @assert ST in ('F','B') && isa(SP,Symbol) "Term error: not supported type parameter."
        new{ST,SP,typeof(value),typeof(neighbor),typeof(couplings),typeof(amplitude),typeof(modulate)}(id,value,neighbor,couplings,amplitude,modulate,factor)
    end
end
function Term{ST,SP}(   id::Symbol,value::Number,neighbor::Any;
                        couplings::Union{TermCouplings,Couplings},
                        amplitude::Union{Function,Nothing}=nothing,
                        modulate::Union{Function,Bool}=false,
                        factor::Number=1
                        ) where {ST,SP}
    isa(couplings,Couplings) && (couplings=TermCouplings(couplings))
    amplitude===nothing && (amplitude=termamplitude)
    isa(modulate,Bool) && (modulate=modulate ? TermModulate(id) : nothing)
    Term{ST,SP}(id,value,neighbor,couplings,amplitude,modulate,factor)
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
    valtype(term::Term)
    valtype(::Type{<:Term{ST,SP,V}}) where {ST,SP,V<:Number}

Get the value type of a term.
"""
Base.valtype(term::Term)=term|>typeof|>valtype
Base.valtype(::Type{<:Term{ST,SP,V}}) where {ST,SP,V<:Number}=V

"""
    rank(term::Term) -> Int
    rank(::Type{T}) where T<:Term -> Int

Get the rank of a term.
"""
rank(term::Term)=term|>typeof|>rank
rank(::Type{T}) where T<:Term=fieldtype(T,:couplings)|>rank

"""
    abbr(term::Term) -> Symbol
    abbr(::Type{<:Term}) -> Symbol

Get the abbreviation of the species of a term.
"""
abbr(term::Term)=term|>typeof|>abbr
abbr(::Type{<:Term})=:tm

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
function Base.show(io::IO,term::Term)
    @printf io "%s{%s}(id=%s,value=%s,neighbor=%s,factor=%s)" species(term) statistics(term) term.id decimaltostr(term.value) term.neighbor decimaltostr(term.factor)
end

"""
    repr(term::Term,bond::AbstractBond,config::IDFConfig) -> String

Get the repr representation of a term on a bond with a given config.
"""
function Base.repr(term::Term,bond::AbstractBond,config::IDFConfig)
    cache=String[]
    if term.neighbor==bond|>neighbor
        value=term.value*term.amplitude(bond)*term.factor
        if !isapprox(abs(value),0.0,atol=atol)
            pids=NTuple{length(bond),pidtype(bond)}(point.pid for point in bond)
            interanls=NTuple{length(bond),valtype(config)}(config[pid] for pid in pids)
            for coupling in values(term.couplings(bond))
                length(expand(coupling,pids,interanls))>0 &&  push!(cache,@sprintf "%s: %s" abbr(term) repr(value*coupling))
            end
        end
    end
    return join(cache,"\n")
end

"""
    replace(term::Term{ST,SP};kwargs...) where {ST,SP} -> Term

Replace some attributes of a term with key word arguments.
"""
@generated function Base.replace(term::Term{ST,SP};kwargs...) where {ST,SP}
    exprs=[:(get(kwargs,$name,getfield(term,$name))) for name in QuoteNode.(term|>fieldnames)]
    return :(Term{ST,SP}($(exprs...)))
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

"""
    expand(otype::Type{<:Operator},term::Term,bond::AbstractBond,config::IDFConfig,table::Union{Table,Nothing}=nothing,half::Union{Bool,Nothing}=false)

Expand the operators of a term on a bond with a given config.

The `half` parameter determines the behavior of generating operators, which falls into the following three categories
* `false`: no extra operations on the generated operators
* `true`: an extra multiplication by 0.5 with the generated operators
* `nothing`: "Hermitian half" of the generated operators
"""
function expand(otype::Type{<:Operator},term::Term,bond::AbstractBond,config::IDFConfig,table::Union{Table,Nothing}=nothing,half::Union{Bool,Nothing}=false)
    result=nothing
    if term.neighbor==bond|>neighbor
        value=term.value*term.amplitude(bond)*term.factor
        if !isapprox(abs(value),0.0,atol=atol)
            result=Operators{idtype(otype),otype}()
            @assert (fieldtype(eltype(idtype(otype)),:seq)===Nothing)==(table===nothing) "expand error: `table` must be assigned if the sequences are required."
            @assert rank(otype)==rank(term) "expand error: dismatched ranks between operator and term."
            rtype,itype=fieldtype(eltype(idtype(otype)),:rcoord),fieldtype(eltype(idtype(otype)),:icoord)
            pids=NTuple{length(bond),pidtype(bond)}(point.pid for point in bond)
            rcoords=NTuple{length(bond),SVector{dimension(bond),Float}}(point.rcoord for point in bond)
            icoords=NTuple{length(bond),SVector{dimension(bond),Float}}(point.icoord for point in bond)
            interanls=NTuple{length(bond),valtype(config)}(config[pid] for pid in pids)
            for coupling in values(term.couplings(bond))
                perm=propercenters(typeof(coupling),coupling.id.centers,Val(rank(bond)))::NTuple{rank(otype),Int}
                orcoords=getcoords(rtype,rcoords,perm)
                oicoords=getcoords(itype,icoords,perm)
                for (coeff,oindexes) in expand(coupling,pids,interanls)
                    isa(table,Table) && any(NTuple{rank(otype),Bool}(!haskey(table,index) for index in oindexes)) && continue
                    id=ID(OID,oindexes,orcoords,oicoords,getseqs(table,oindexes))
                    if !(half===nothing && haskey(result,id'))
                        ovalue=valtype(otype)(value*coeff*getfactor(half,id))
                        add!(result,otype.name.wrapper(ovalue,id))
                    end
                end
            end
        end
    end
    return result
end
getcoords(::Type{<:Nothing},rcoords::NTuple{N,SVector{M,Float}},perm::NTuple{R,Int}) where {N,M,R}=NTuple{R,Nothing}(nothing for i=1:R)
getcoords(::Type{<:SVector},rcoords::NTuple{N,SVector{M,Float}},perm::NTuple{R,Int}) where {N,M,R}=NTuple{R,SVector{M,Float}}(rcoords[p] for p in perm)
getseqs(::Nothing,indexes::NTuple{N,<:OID}) where N=NTuple{N,Nothing}(nothing for i=1:N)
getseqs(table::Table{I},indexes::NTuple{N,I}) where {N,I<:Index}=NTuple{N,Int}(table[index] for index in indexes)
getfactor(half::Bool,::ID)=half ? 0.5 : 1.0
getfactor(::Nothing,id::ID)=isHermitian(id) ? 0.5 : 1.0

"""
    update!(term::Term,args...;kwargs...) -> Term

Update the value of a term by its `modulate` function.
"""
update!(term::Term,args...;kwargs...)=(term.value=term.modulate(args...;kwargs...)::valtype(term);term)

end  # module
