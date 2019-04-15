module Terms

using Printf: @printf,@sprintf
using StaticArrays: SVector
using ..Spatials: AbstractBond,pidtype
using ..DegreesOfFreedom: Index,IDFConfig,Table,OID,Operators,oidtype,otype
using ...Interfaces: add!
using ...Prerequisites: Float,atol,decimaltostr
using ...Prerequisites.TypeTraits: efficientoperations,indtosub,corder
using ...Prerequisites.CompositeStructures: CompositeTuple
using ...Mathematics.AlgebraOverFields: SimpleID,ID,Element,Elements,idtype

import ..DegreesOfFreedom: isHermitian
import ...Interfaces: id,rank,kind,expand,expand!,dimension,update!

export Subscript,Subscripts,@subscript
export Coupling,Couplings
export TermFunction,TermAmplitude,TermCouplings,TermModulate
export Term,statistics,abbr

const wildcard='*'
const constant=':'
"""
    Subscript(  ipattern::NTuple{N1,Any},
                opattern::NTuple{N2,Any},
                mapping::Union{Function,Nothing}=nothing,
                constrain::Union{Function,Nothing}=nothing,
                identifier::Union{Symbol,Char}=wildcard
                ) where {N1,N2}
    Subscript{N}() where N
    Subscript(opattern::NTuple{N,Int}) where N

The subscripts of some orbital/spin degrees of freedom.
"""
struct Subscript{N1,N2,I<:Tuple,O<:Tuple,M<:Union{Function,Nothing},C<:Union{Function,Nothing},D<:Union{Symbol,Char}}
    ipattern::I
    opattern::O
    mapping::M
    constrain::C
    identifier::D
    function Subscript( ipattern::NTuple{N1,Any},
                        opattern::NTuple{N2,Any},
                        mapping::Union{Function,Nothing}=nothing,
                        constrain::Union{Function,Nothing}=nothing,
                        identifier::Union{Symbol,Char}=wildcard
                        ) where {N1,N2}
        new{N1,N2,typeof(ipattern),typeof(opattern),typeof(mapping),typeof(constrain),typeof(identifier)}(ipattern,opattern,mapping,constrain,identifier)
    end
end
Subscript{N}() where N=Subscript((wildcard,),NTuple{N,Char}(wildcard for i=1:N),nothing,nothing,wildcard)
Subscript(opattern::NTuple{N,Int}) where N=Subscript((),opattern,nothing,nothing,constant)

"""
    ==(sub1::Subscript,sub2::Subscript) -> Bool
    isequal(sub1::Subscript,sub2::Subscript) -> Bool

Judge whether two subscripts are equivalent to each other.
"""
Base.:(==)(sub1::Subscript,sub2::Subscript)=sub1.ipattern==sub2.ipattern && sub1.opattern==sub2.opattern && sub1.identifier==sub2.identifier
Base.isequal(sub1::Subscript,sub2::Subscript)=isequal(sub1.ipattern,sub2.ipattern) && isequal(sub1.opattern,sub2.opattern) && isequal(sub1.identifier,sub2.identifier)

"""
    show(io::IO,subscript::Subscript)

Show a subscript.
"""
function Base.show(io::IO,subscript::Subscript)
    if subscript.identifier==constant || subscript.identifier==wildcard
        @printf io "(%s)" join(subscript.opattern,',')
    else
        @printf io "(%s) with %s" join(subscript.opattern,',') subscript.identifier
    end
end

"""
    rank(subscript::Subscript) -> Int
    rank(::Type{<:Subscript{N}}) where N -> Int

Get the number of the independent variables that are used to describe the subscripts of some orbital/spin degrees of freedom.
"""
rank(subscript::Subscript)=subscript|>typeof|>rank
rank(::Type{<:Subscript{N}}) where N=N

"""
    dimension(subscript::Subscript) -> Int
    dimension(::Type{<:(Subscript{N1,N2} where N1)}) where N2 -> Int

Get the number of the whole variables that are used to describe the subscripts of some orbital/spin degrees of freedom.
"""
dimension(subscript::Subscript)=subscript|>typeof|>dimension
dimension(::Type{<:(Subscript{N1,N2} where N1)}) where N2=N2

"""
    (subscript::Subscript{N})(::Val{'M'},values::Vararg{Int,N}) where N -> NTuple{dimension(subscript),Int}
    (subscript::Subscript{N})(::Val{'C'},values::Vararg{Int,N}) where N -> Bool

* Construct the subscripts from a set of independent variables.
* Judge whether a set of independent variables are valid to construct the subscripts.
"""
function (subscript::Subscript)() end
(subscript::Subscript{N})(::Val{'M'},values::Vararg{Int,N}) where N=subscript.mapping(values...)
(subscript::Subscript{0,N,<:Tuple,<:Tuple,Nothing,Nothing})(::Val{'M'}) where N=subscript.opattern
(subscript::Subscript{1,N,<:Tuple,<:Tuple,Nothing,Nothing})(::Val{'M'},value::Int) where N=NTuple{N,Int}(value for i=1:N)
(subscript::Subscript{N})(::Val{'C'},values::Vararg{Int,N}) where N=subscript.constrain(values...)
(subscript::Subscript{N1,N2,<:Tuple,<:Tuple,<:Union{Function,Nothing},Nothing})(::Val{'C'},values::Vararg{Int,N1}) where {N1,N2}=true

"""
    @subscript expr::Expr with constrain::Expr -> Subscript

Construct a subscript from a map and optionally with a constrain.
"""
macro subscript(expr::Expr,with::Symbol=:with,constrain::Union{Expr,Symbol}=:nothing)
    @assert with==:with "@subscript error: pattern and constrain must be separated by :with."
    subscriptexpr(expr,constrain)
end
function subscriptexpr(pattern::Expr,::Symbol)
    @assert pattern.head==:tuple "subscriptexpr error: wrong input pattern (not a tuple)."
    ip,op=Tuple(subscriptipattern(pattern.args)),Tuple(pattern.args)
    mapname,identifier=gensym(),QuoteNode(gensym())
    return quote
        $mapname($(ip...))=$pattern
        Subscript($ip,$op,$mapname,nothing,$identifier)
    end
end
function subscriptexpr(pattern::Expr,constrain::Expr)
    @assert pattern.head==:tuple "subscriptexpr error: wrong input pattern (not a tuple)."
    ip,op=Tuple(subscriptipattern(pattern.args)),Tuple(pattern.args)
    mapname,constrainname,identifier=gensym(),gensym(),QuoteNode(gensym())
    return quote
        $mapname($(ip...))=$pattern
        $constrainname($(ip...))=$(constrain)
        Subscript($ip,$op,$mapname,$constrainname,$identifier)
    end
end
function subscriptipattern(opattern)
    result=[]
    for op in opattern
        @assert op≠Symbol(wildcard) && op≠Symbol(constant) "subscriptipattern error: wrong symbols."
        isa(op,Symbol) && op ∉ result && push!(result,op)
    end
    return result
end

"""
    Subscripts(contents::Subscript...)

A complete set of all the independent subscripts of the orbital/spin degrees of freedom.
"""
struct Subscripts{T<:Tuple,R,D} <: CompositeTuple{T}
    contents::T
    function Subscripts(contents::T) where T<:Tuple
        R=sum(rank(fieldtype(T,i)) for i=1:fieldcount(T))
        D=sum(dimension(fieldtype(T,i)) for i=1:fieldcount(T))
        new{T,R,D}(contents)
    end
end
Subscripts(contents::Subscript...)=Subscripts(contents)

"""
    rank(subscripts::Subscripts) -> Int
    rank(::Type{S}) where S<:Subscripts -> Int

Get the total number of the independent variables of the complete subscript set.
"""
rank(subscripts::Subscripts)=subscripts|>typeof|>rank
rank(::Type{<:Subscripts{<:Tuple,R}}) where R=R

"""
    rank(subscripts::Subscripts,i::Int) -> Int
    rank(::Type{<:Subscripts{T}},i::Int) where T -> Int

Get the number of the independent variables of a component of the complete subscript set.
"""
rank(subscripts::Subscripts,i::Int)=rank(typeof(subscripts),i)
rank(::Type{<:Subscripts{T}},i::Int) where T=rank(fieldtype(T,i))

"""
    dimension(subscripts::Subscripts) -> Int
    dimension(::Type{<:(Subscripts{<:Tuple,R,D} where R)}) where D -> Int

Get the total number of the whole variables of the complete subscript set.
"""
dimension(subscripts::Subscripts)=subscripts|>typeof|>dimension
dimension(::Type{<:(Subscripts{<:Tuple,R,D} where R)}) where D=D

"""
    dimension(subscripts::Subscripts,i::Int) -> Int
    dimension(::Type{<:Subscripts{T}},i::Int) where T -> Int

Get the total number of the whole variables of a component of the complete subscript set.
"""
dimension(subscripts::Subscripts,i::Int)=dimension(typeof(subscripts),i)
dimension(::Type{<:Subscripts{T}},i::Int) where T=dimension(fieldtype(T,i))

"""
    (subscripts::Subscripts)(::Val{'M'},values::NTuple{N,Int}) where N -> NTuple{dimension(subscripts),Int}
    (subscripts::Subscripts)(::Val{'C'},values::NTuple{N,Int}) where N -> Bool

* Construct the complete set of subscripts from a complete set of independent variables.
* Judge whether a complete set of independent variables are valid to construct the complete subscripts.
"""
function (subscripts::Subscripts)() end
@generated function (subscripts::Subscripts)(::Val{'M'},values::NTuple{N,Int}) where N
    @assert rank(subscripts)==N
    exprs,count=[],1
    for i=1:length(subscripts)
        vs=Tuple(:(values[$j]) for j=count:(count+rank(subscripts,i)-1))
        push!(exprs,:(subscripts.contents[$i](Val('M'),$(vs...))...))
        count=count+rank(subscripts,i)
    end
    return Expr(:tuple,exprs...)
end
@generated function (subscripts::Subscripts)(::Val{'C'},values::NTuple{N,Int}) where N
    @assert rank(subscripts)==N
    length(subscripts)==0 && return :(true)
    count=rank(subscripts,1)
    vs=Tuple(:(values[$j]) for j=1:count)
    expr=:(subscripts.contents[1](Val('C'),$(vs...)))
    for i=2:length(subscripts)
        vs=Tuple(:(values[$j]) for j=(count+1):(count+rank(subscripts,i)))
        expr=Expr(:(&&),expr,:(subscripts.contents[$i](Val('C'),$(vs...))))
        count=count+rank(subscripts,i)
    end
    return expr
end

"""
    *(sub1::Subscript,sub2::Subscript) -> Subscripts
    *(subs::Subscripts,sub::Subscript) -> Subscripts
    *(sub::Subscript,subs::Subscripts) -> Subscripts
    *(subs1::Subscripts,subs2::Subscripts) -> Subscripts

Get the multiplication between subscripts or complete sets of subscripts.
"""
Base.:*(sub1::Subscript,sub2::Subscript)=Subscripts(sub1,sub2)
Base.:*(subs::Subscripts,sub::Subscript)=Subscripts(subs.contents...,sub)
Base.:*(sub::Subscript,subs::Subscripts)=Subscripts(sub,subs.contents...)
Base.:*(subs1::Subscripts,subs2::Subscripts)=Subscripts(subs1.contents...,subs2.contents...)

"""
    expand(subscripts::Subscripts,dimensions::NTuple{N,Int}) where N -> SbExpand

Expand a complete set of subscripts with a given set of variable ranges.
"""
function expand(subscripts::Subscripts,dimensions::NTuple{N,Int}) where N
    @assert dimension(subscripts)==N "expand error: dismatched input dimensions $dimensions."
    dims,dcount,rcount=Vector{Int}(undef,rank(subscripts)),0,0
    for i=1:length(subscripts)
        ipattern=subscripts[i].ipattern
        opattern=subscripts[i].opattern
        for j=1:dimension(subscripts,i)
            isa(opattern[j],Int) && @assert 0<opattern[j]<=dimensions[dcount+j] "expand error: opattern($opattern) out of range."
        end
        for j=1:rank(subscripts,i)
            index,flag=1,false
            while (index=findnext(isequal(ipattern[j]),opattern,index))≠nothing
                flag || (dims[rcount+j]=dimensions[dcount+index];flag=true)
                flag && @assert dimensions[dcount+index]==dims[rcount+j] "expand error: dismatched input dimensions."
                index=index+1
            end
        end
        dcount=dcount+dimension(subscripts,i)
        rcount=rcount+rank(subscripts,i)
    end
    return SbExpand(subscripts,NTuple{rank(subscripts),Int}(dims))
end
struct SbExpand{N,S<:Subscripts,D}
    subscripts::S
    dims::NTuple{D,Int}
    function SbExpand(subscripts::Subscripts,dims::NTuple{D,Int}) where D
        @assert D==rank(subscripts) "SbExpand error: dismatched inputs."
        new{dimension(subscripts),typeof(subscripts),D}(subscripts,dims)
    end
end
Base.eltype(::Type{<:SbExpand{N}}) where N=NTuple{N,Int}
Base.IteratorSize(::Type{<:SbExpand})=Base.SizeUnknown()
function Base.iterate(sbe::SbExpand,state::Int=1)
    while state<=prod(sbe.dims)
        inds=indtosub(sbe.dims,state,corder)
        sbe.subscripts(Val('C'),inds) && return (sbe.subscripts(Val('M'),inds),state+1)
        state=state+1
    end
    return nothing
end

"""
    Coupling{N,V<:Number,I<:ID{<:NTuple{N,SimpleID}}} <: Element{N,V,I}

The coupling intra/inter interanl degrees of freedom at different lattice points.
"""
abstract type Coupling{N,V<:Number,I<:ID{<:NTuple{N,SimpleID}}} <: Element{N,V,I} end

couplingcenter(::Type{<:Coupling},i::Int,n::Int,::Val{1})=1
couplingcenter(::Type{<:Coupling},i::Int,n::Int,::Val{R}) where R=error("couplingcenter error: no default center for a rank-$R bond.")
@generated function couplingcenters(::Type{C},centers::NTuple{N,Any},::Val{R}) where {C<:Coupling,N,R}
    exprs=[:(isa(centers[$i],Int) ? (0<centers[$i]<=R ? centers[$i] : error("couplingcenters error: center out of range.")) : couplingcenter(C,$i,N,Val(R))) for i=1:N]
    return Expr(:tuple,exprs...)
end

"""
    Couplings{I<:ID,C<:Coupling} <: AbstractDict{I,C}

A pack of couplings intra/inter interanl degrees of freedom at different lattice points.

Alias for `Elements{I,C}`.
"""
const Couplings{I<:ID,C<:Coupling}=Elements{I,C}
Couplings(cps::Coupling...)=Elements(cps...)

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
    TermAmplitude(amplitude::Union{Function,Nothing}=nothing)

The function for the amplitude of a term.
"""
struct TermAmplitude{A<:Union{Function,Nothing}} <: TermFunction
    amplitude::A
    TermAmplitude(amplitude::Union{Function,Nothing}=nothing)=new{typeof(amplitude)}(amplitude)
end
(termamplitude::TermAmplitude{Nothing})(args...;kwargs...)=1
(termamplitude::TermAmplitude{<:Function})(args...;kwargs...)=termamplitude.amplitude(args...;kwargs...)

"""
    TermCouplings(candidate::Coupling)
    TermCouplings(candidate::Couplings)
    TermCouplings(contents::Tuple{<:Tuple{Vararg{Couplings}},<:Function})
    TermCouplings(candidates::NTuple{N,<:Couplings},choice::Function) where N

The function for the couplings of a term.
"""
struct TermCouplings{C<:Union{Function,Couplings}} <: TermFunction
    couplings::C
    TermCouplings(coupling::Coupling)=new{Couplings{coupling.id|>typeof,coupling|>typeof}}(Couplings(coupling))
    TermCouplings(couplings::Couplings)=new{couplings|>typeof}(couplings)
    TermCouplings(couplings::Function)=new{couplings|>typeof}(couplings)
end
(termcouplings::TermCouplings{<:Couplings})(args...;kwargs...)=termcouplings.couplings
(termcouplings::TermCouplings{<:Function})(args...;kwargs...)=termcouplings.couplings(args...;kwargs...)

"""
    TermModulate(id::Symbol,modulate::Union{Function,Nothing}=nothing)
    TermModulate(id::Symbol,modulate::Bool)

The function for the modulation of a term.
"""
struct TermModulate{M<:Union{Function,Nothing},id} <: TermFunction
    modulate::M
    TermModulate(id::Symbol,modulate::Union{Function,Nothing}=nothing)=new{typeof(modulate),id}(modulate)
    TermModulate(id::Symbol,modulate::Bool)=(@assert modulate "TermModulate error: input `modulate` must be `true`.";new{Nothing,id}(nothing))
end
(termmodulate::TermModulate{Nothing,id})(;kwargs...) where id=get(kwargs,id,nothing)
(termmodulate::TermModulate{<:Function})(args...;kwargs...)=termmodulate.modulate(args...;kwargs...)

"""
    Term{ST,K,R,I}(value::Number,bondkind::Any,couplings::TermCouplings,amplitude::TermAmplitude,modulate::Union{TermModulate,Nothing},factor::Number) where {ST,K,R,I}
    Term{ST,K,R}(  id::Symbol,value::Number,bondkind::Any;
                    couplings::Union{Function,Coupling,Couplings,TermCouplings},
                    amplitude::Union{Function,Nothing}=nothing,
                    modulate::Union{Function,Bool}=false,
                    ) where {ST,K,R}

A term of a quantum lattice system.
"""
mutable struct Term{statistics,kind,R,id,V<:Number,B<:Any,C<:TermCouplings,A<:TermAmplitude,M<:Union{TermModulate,Nothing}}
    value::V
    bondkind::B
    couplings::C
    amplitude::A
    modulate::M
    factor::V
    function Term{ST,K,R,I}(value::Number,bondkind::Any,couplings::TermCouplings,amplitude::TermAmplitude,modulate::Union{TermModulate,Nothing},factor::Number) where {ST,K,R,I}
        @assert ST in ('F','B') "Term error: statistics must be 'F' or 'B'."
        @assert isa(K,Symbol) "Term error: kind must be a Symbol."
        @assert isa(I,Symbol) "Term error: id must be a Symbol."
        new{ST,K,R,I,typeof(value),typeof(bondkind),typeof(couplings),typeof(amplitude),typeof(modulate)}(value,bondkind,couplings,amplitude,modulate,factor)
    end
end
function Term{ST,K,R}(  id::Symbol,value::Number,bondkind::Any;
                        couplings::Union{Function,Coupling,Couplings,TermCouplings},
                        amplitude::Union{Function,Nothing}=nothing,
                        modulate::Union{Function,Bool}=false,
                        ) where {ST,K,R}
    isa(couplings,TermCouplings) || (couplings=TermCouplings(couplings))
    isa(amplitude,TermAmplitude) || (amplitude=TermAmplitude(amplitude))
    isa(modulate,TermModulate) || (modulate=modulate===false ? nothing : TermModulate(id,modulate))
    Term{ST,K,R,id}(value,bondkind,couplings,amplitude,modulate,1)
end

"""
    statistics(term::Term) -> Char
    statistics(::Type{<:Term{ST}}) where ST -> Char

Get the statistics of a term.
"""
statistics(term::Term)=term|>typeof|>statistics
statistics(::Type{<:Term{ST}}) where ST=ST

"""
    kind(term::Term) -> Symbol
    kind(::Type{<:Term{ST,K} where ST}) where K -> Symbol

Get the kind of a term.
"""
kind(term::Term)=term|>typeof|>kind
kind(::Type{<:Term{ST,K} where ST}) where K=K

"""
    rank(term::Term) -> Int
    rank(::Type{<:Term{ST,K,R} where {ST,K}}) where R -> Int

Get the rank of a term.
"""
rank(term::Term)=term|>typeof|>rank
rank(::Type{<:Term{ST,K,R} where {ST,K}}) where R=R

"""
    id(term::Term) -> Symbol
    id(::Type{<:Term{ST,K,R,I} where {ST,K,R}}) where I -> Symbol

Get the id of a term.
"""
id(term::Term)=term|>typeof|>id
id(::Type{<:Term{ST,K,R,I} where {ST,K,R}}) where I=I

"""
    valtype(term::Term)
    valtype(::Type{<:Term{ST,K,R,I,V} where {ST,K,R,I}}) where {V<:Number}

Get the value type of a term.
"""
Base.valtype(term::Term)=term|>typeof|>valtype
Base.valtype(::Type{<:Term{ST,K,R,I,V} where {ST,K,R,I}}) where {V<:Number}=V

"""
    abbr(term::Term) -> Symbol
    abbr(::Type{<:Term}) -> Symbol

Get the abbreviation of the kind of a term.
"""
abbr(term::Term)=term|>typeof|>abbr
abbr(::Type{<:Term})=:tm

"""
    isHermitian(term::Term) -> Bool
    isHermitian(::Type{<:Term}) -> Bool
"""
isHermitian(term::Term)=term|>typeof|>isHermitian
isHermitian(::Type{<:Term})=error("isHermitian error: not implemented.")

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
    @printf io "%s{%s%s}(id=%s,value=%s,bondkind=%s,factor=%s)" kind(term) rank(term) statistics(term) term|>id decimaltostr(term.value) term.bondkind decimaltostr(term.factor)
end

"""
    repr(term::Term,bond::AbstractBond,config::IDFConfig) -> String

Get the repr representation of a term on a bond with a given config.
"""
function Base.repr(term::Term,bond::AbstractBond,config::IDFConfig)
    cache=String[]
    if term.bondkind==bond|>kind
        value=term.value*term.amplitude(bond)*term.factor
        if !isapprox(abs(value),0.0,atol=atol)
            pids=NTuple{length(bond),pidtype(bond)}(point.pid for point in bond)
            interanls=NTuple{length(bond),valtype(config)}(config[pid] for pid in pids)
            for coupling in values(term.couplings(bond))
                length(expand(coupling,pids,interanls,term|>kind|>Val))>0 &&  push!(cache,@sprintf "%s: %s" abbr(term) repr(value*coupling))
            end
        end
    end
    return join(cache,"\n")
end

"""
    replace(term::Term{ST,K,R,I};kwargs...) where {ST,K,R,I} -> Term

Replace some attributes of a term with key word arguments.
"""
@generated function Base.replace(term::Term{ST,K,R,I};kwargs...) where {ST,K,R,I}
    exprs=[:(get(kwargs,$name,getfield(term,$name))) for name in QuoteNode.(term|>fieldnames)]
    return :(Term{ST,K,R,I}($(exprs...)))
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
    expand!(operators::Operators,term::Term,bond::AbstractBond,config::IDFConfig,table::Union{Table,Nothing}=nothing,half::Bool=false) -> Operators

Expand the operators of a term on a bond with a given config.

The `half` parameter determines the behavior of generating operators, which falls into the following two categories
* `true`: "Hermitian half" of the generated operators
* `false`: "Hermitian whole" of the generated operators
"""
function expand!(operators::Operators,term::Term,bond::AbstractBond,config::IDFConfig,table::Union{Table,Nothing}=nothing,half::Bool=false)
    if term.bondkind==bond|>kind
        value=term.value*term.amplitude(bond)*term.factor
        if !isapprox(abs(value),0.0,atol=atol)
            otype,ktype=operators|>typeof|>valtype,operators|>typeof|>keytype
            @assert (fieldtype(ktype|>eltype,:seq)===Nothing)==(table===nothing) "expand! error: `table` must be assigned if the sequences are required."
            @assert rank(otype)==rank(term) "expand! error: dismatched ranks between operator and term."
            apriori=isHermitian(term)
            record=apriori===nothing && length(operators)>0 ? Set{ktype}() : nothing
            rtype,itype=fieldtype(ktype|>eltype,:rcoord),fieldtype(ktype|>eltype,:icoord)
            pids=NTuple{length(bond),pidtype(bond)}(point.pid for point in bond)
            rcoords=NTuple{length(bond),SVector{dimension(bond),Float}}(point.rcoord for point in bond)
            icoords=NTuple{length(bond),SVector{dimension(bond),Float}}(point.icoord for point in bond)
            interanls=NTuple{length(bond),valtype(config)}(config[pid] for pid in pids)
            for coupling in values(term.couplings(bond))
                perm=couplingcenters(typeof(coupling),coupling.id.centers,Val(rank(bond)))::NTuple{rank(otype),Int}
                orcoords,oicoords=termcoords(rtype,rcoords,perm),termcoords(itype,icoords,perm)
                for (coeff,oindexes) in expand(coupling,pids,interanls,term|>kind|>Val) # needs improvement memory allocation 7 times for each fock coupling
                    isa(table,Table) && any(NTuple{rank(otype),Bool}(!haskey(table,index) for index in oindexes)) && continue
                    id=ID(OID,oindexes,orcoords,oicoords,termseqs(table,oindexes))
                    if apriori===true
                        add!(operators,otype.name.wrapper(convert(otype|>valtype,value*coeff*(half ? 0.5 : 1.0)),id))
                    elseif apriori==false
                        opt=otype.name.wrapper(convert(otype|>valtype,value*coeff),id)
                        add!(operators,opt)
                        half || add!(operators,opt')
                    else
                        if !(record===nothing ? haskey(operators,id') : in(id',record))
                            record===nothing || push!(record,id)
                            ovalue=valtype(otype)(value*coeff*termfactor(nothing,id,term|>kind|>Val)) # needs improvement memory allocation 2 times for each
                            opt=otype.name.wrapper(ovalue,id)
                            add!(operators,opt)
                            half || add!(operators,opt')
                        end
                    end
                end
            end
        end
    end
    return operators
end
termcoords(::Type{<:Nothing},rcoords::NTuple{N,SVector{M,Float}},perm::NTuple{R,Int}) where {N,M,R}=NTuple{R,Nothing}(nothing for i=1:R)
termcoords(::Type{<:SVector},rcoords::NTuple{N,SVector{M,Float}},perm::NTuple{R,Int}) where {N,M,R}=NTuple{R,SVector{M,Float}}(rcoords[p] for p in perm)
termseqs(::Nothing,indexes::NTuple{N,<:Index}) where N=NTuple{N,Nothing}(nothing for i=1:N)
@generated termseqs(table::Table{I},indexes::NTuple{N,I}) where {N,I<:Index}=Expr(:tuple,[:(table[indexes[$i]]) for i=1:N]...)
termfactor(::Nothing,id::ID{<:NTuple{N,OID}},::Val{K}) where {N,K}=isHermitian(id) ? 0.5 : 1.0

"""
    expand(term::Term,bond::AbstractBond,config::IDFConfig,table::Union{Table,Nothing}=nothing,half::Bool=false) -> Operators

Expand the operators of a term on a bond with a given config.
"""
function expand(term::Term,bond::AbstractBond,config::IDFConfig,table::Union{Table,Nothing}=nothing,half::Bool=false)
    optype=otype(term|>typeof,oidtype(config|>typeof|>valtype|>eltype,bond|>typeof,table|>typeof))
    expand!(Operators{idtype(optype),optype}(),term,bond,config,table,half)
end

"""
    update!(term::Term,args...;kwargs...) -> Term

Update the value of a term by its `modulate` function.
"""
update!(term::Term,args...;kwargs...)=(term.value=term.modulate(args...;kwargs...)::valtype(term);term)

end  # module
