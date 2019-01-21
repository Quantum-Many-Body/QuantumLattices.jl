using Printf: @printf,@sprintf
using ..Spatials: PID,AbstractBond,Point,Bond,pidtype
using ..DegreesOfFreedom: IID,Index,Internal,FilteredAttributes,Coupling,Couplings
using ...Prerequisites: delta,decimaltostr
using ...Prerequisites.TypeTraits: indtosub,corder
using ...Prerequisites.CompositeStructures: CompositeTuple
using ...Mathematics.AlgebraOverFields: SimpleID,ID
using ...Mathematics.VectorSpaces: AbstractVectorSpace,IsMultiIndexable,MultiIndexOrderStyle

import ...Prerequisites.Interfaces: dims,inds,rank,dimension,⊗,expand

export dims,inds,rank,dimension,⊗,expand
export ANNIHILATION,CREATION,FID,FIndex,Fock
export usualfockindextotuple,nambufockindextotuple
export Subscript,Subscripts,@subscript
export FCID,FockCoupling
export σ⁰,σˣ,σʸ,σᶻ,σ⁺,σ⁻

"""
    ANNIHILATION

Indicate that the nambu index is ANNIHILATION.
"""
const ANNIHILATION=1
"""
    CREATION

Indicate that the nambu index is CREATION.
"""
const CREATION=2

"""
    FID <: IID

The Fock id.
"""
struct FID <: IID
    orbital::Int
    spin::Int
    nambu::Int
end
Base.fieldnames(::Type{FID})=(:orbital,:spin,:nambu)

"""
    FID(;orbital::Int=1,spin::Int=1,nambu::Int=ANNIHILATION)

Create a Fock id.
"""
FID(;orbital::Int=1,spin::Int=1,nambu::Int=ANNIHILATION)=FID(orbital,spin,nambu)

"""
    adjoint(fid::FID) -> FID

Get the adjoint of a Fock id.
"""
Base.adjoint(fid::FID)=FID(fid.orbital,fid.spin,3-fid.nambu)

"""
    Fock <: Internal{FID}

The Fock interanl degrees of freedom.
"""
struct Fock <: Internal{FID}
    atom::Int
    norbital::Int
    nspin::Int
    nnambu::Int
end
IsMultiIndexable(::Type{Fock})=IsMultiIndexable(true)
MultiIndexOrderStyle(::Type{Fock})=MultiIndexOrderStyle('C')
dims(fock::Fock)=(fock.norbital,fock.nspin,fock.nnambu)
inds(fid::FID,::Fock)=(fid.orbital,fid.spin,fid.nambu)
FID(inds::NTuple{3,Int},::Fock)=FID(inds[1],inds[2],inds[3])

"""
    Fock(;atom::Int=1,norbital::Int=1,nspin::Int=2,nnambu::Int=2)

Construct a Fock degrees of freedom.
"""
Fock(;atom::Int=1,norbital::Int=1,nspin::Int=2,nnambu::Int=2)=Fock(atom,norbital,nspin,nnambu)

"""
    FIndex{S} <: Index{PID{S},FID}

The Fock index.
"""
struct FIndex{S} <: Index{PID{S},FID}
    scope::S
    site::Int
    orbital::Int
    spin::Int
    nambu::Int
end
Base.fieldnames(::Type{<:FIndex})=(:scope,:site,:orbital,:spin,:nambu)

"""
    union(::Type{P},::Type{FID})

Get the union type of `PID` and `FID`.
"""
Base.union(::Type{P},::Type{FID}) where P<:PID=FIndex{fieldtype(P,:scope)}

"""
    usualfockindextotuple

Indicate that the filtered attributes are `(:scope,:site,:orbital,:spin)` when converting a Fock index to tuple.
"""
const usualfockindextotuple=FilteredAttributes(:scope,:site,:orbital,:spin)
"""
    nambufockindextotuple

Indicate that the filtered attributes are `(:scope,:nambu,:site,:orbital,:spin)` when converting a Fock index to tuple.
"""
const nambufockindextotuple=FilteredAttributes(:scope,:nambu,:site,:orbital,:spin)

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
    if subscript.identifier==constant
        @printf io "(%s)" join(subscript.opattern,',')
    elseif subscript.identifier==wildcard
        @printf io "%s=>(%s)" join(subscript.ipattern,',') join(subscript.opattern,',')
    else
        @printf io "(%s)=>(%s) with %s" join(subscript.ipattern,',') join(subscript.opattern,',') subscript.identifier
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
    dimension(::Type{<:Subscript{N1,N2}}) where {N1,N2} -> Int

Get the number of the whole variables that are used to describe the subscripts of some orbital/spin degrees of freedom.
"""
dimension(subscript::Subscript)=subscript|>typeof|>dimension
dimension(::Type{<:Subscript{N1,N2}}) where {N1,N2}=N2

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
macro subscript(expr::Expr,with::Symbol=:with,constrain::Union{Expr,Symbol}=:nothing,gensym::Expr=:(gensym=false))
    @assert expr.head==:call && expr.args[1]==:(=>)
    @assert isa(expr.args[2],Expr) && expr.args[2].head==:tuple && all(expr.args[2].args.≠Symbol(wildcard)) && all(expr.args[2].args.≠Symbol(constant))
    @assert isa(expr.args[3],Expr) && expr.args[3].head==:tuple && all(expr.args[3].args.≠Symbol(wildcard)) && all(expr.args[3].args.≠Symbol(constant))
    @assert with==:with
    @assert isa(gensym,Expr) && gensym.head==:(=) && gensym.args[1]==:gensym && isa(gensym.args[2],Bool)
    ip=gensym.args[2] ? Tuple(Base.gensym(arg) for arg in expr.args[2].args) : Tuple(expr.args[2].args)
    op=gensym.args[2] ? Tuple(ip[findfirst(isequal(arg),expr.args[2].args)] for arg in expr.args[3].args) : Tuple(expr.args[3].args)
    mapname,identifier=Base.gensym(),QuoteNode(Base.gensym())
    if constrain===:nothing
        return quote
            $mapname($(expr.args[2].args...))=$(expr.args[3])
            Subscript($ip,$op,$mapname,nothing,$identifier)
        end
    else
        constrainname=Base.gensym()
        return quote
            $mapname($(expr.args[2].args...))=$(expr.args[3])
            $constrainname($(expr.args[2].args...))=$(constrain)
            Subscript($ip,$op,$mapname,$constrainname,$identifier)
        end
    end
end

"""
    Subscripts(contents::Subscript...)

A complete set of all the independent subscripts of the orbital/spin degrees of freedom.
"""
struct Subscripts{T<:Tuple} <: CompositeTuple{T}
    contents::T
end
Subscripts(contents::Subscript...)=Subscripts(contents)

"""
    rank(subscripts::Subscripts) -> Int
    rank(::Type{S}) where S<:Subscripts -> Int

Get the total number of the independent variables of the complete subscript set.
"""
rank(subscripts::Subscripts)=subscripts|>typeof|>rank
rank(::Type{S}) where S<:Subscripts=sum(rank(S,i) for i=1:length(S))

"""
    rank(subscripts::Subscripts,i::Int) -> Int
    rank(::Type{<:Subscripts{T}},i::Int) where T -> Int

Get the number of the independent variables of a component of the complete subscript set.
"""
rank(subscripts::Subscripts,i::Int)=rank(typeof(subscripts),i)
rank(::Type{<:Subscripts{T}},i::Int) where T=rank(fieldtype(T,i))

"""
    dimension(subscripts::Subscripts) -> Int
    dimension(::Type{S}) where S<:Subscripts -> Int

Get the total number of the whole variables of the complete subscript set.
"""
dimension(subscripts::Subscripts)=subscripts|>typeof|>dimension
dimension(::Type{S}) where S<:Subscripts=sum(dimension(S,i) for i=1:length(S))

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
    FCID(;center=wildcard,atom=wildcard,orbital=wildcard,spin=wildcard,nambu=wildcard,obsub=wildcard,spsub=wildcard)

The id of a Fock coupling.
"""
struct FCID{C,A,O,S,N,OS,SS} <: SimpleID
    center::C
    atom::A
    orbital::O
    spin::S
    nambu::N
    obsub::OS
    spsub::SS
end
Base.fieldnames(::Type{<:FCID})=(:center,:atom,:orbital,:spin,:nambu,:obsub,:spsub)
FCID(;center=wildcard,atom=wildcard,orbital=wildcard,spin=wildcard,nambu=wildcard,obsub=wildcard,spsub=wildcard)=FCID(center,atom,orbital,spin,nambu,obsub,spsub)

"""
    FockCoupling(value::Number,id::ID{N,I},obsubscripts::Subscripts,spsubscripts::Subscripts) where {N,I<:FCID}
    FockCoupling{N}(value::Number;
                    centers::Union{NTuple{N,Int},Nothing}=nothing,
                    atoms::Union{NTuple{N,Int},Nothing}=nothing,
                    orbitals::Union{NTuple{N,Int},Subscript,Nothing}=nothing,
                    spins::Union{NTuple{N,Int},Subscript,Nothing}=nothing,
                    nambus::Union{NTuple{N,Int},Nothing}=nothing) where N

Fock coupling.
"""
struct FockCoupling{N,V<:Number,I<:ID{N,<:FCID},OS<:Subscripts,SS<:Subscripts} <: Coupling{V,I}
    value::V
    id::I
    obsubscripts::OS
    spsubscripts::SS
    function FockCoupling(value::Number,id::ID{N,I},obsubscripts::Subscripts,spsubscripts::Subscripts) where {N,I<:FCID}
        new{N,value|>typeof,id|>typeof,obsubscripts|>typeof,spsubscripts|>typeof}(value,id,obsubscripts,spsubscripts)
    end
end
function FockCoupling{N}(   value::Number;
                            centers::Union{NTuple{N,Int},Nothing}=nothing,
                            atoms::Union{NTuple{N,Int},Nothing}=nothing,
                            orbitals::Union{NTuple{N,Int},Subscript,Nothing}=nothing,
                            spins::Union{NTuple{N,Int},Subscript,Nothing}=nothing,
                            nambus::Union{NTuple{N,Int},Nothing}=nothing) where N
    centers===nothing && (centers=ntuple(i->wildcard,N))
    atoms===nothing && (atoms=ntuple(i->wildcard,N))
    obsubscript=isa(orbitals,NTuple{N,Int}) ? Subscript(orbitals) : isa(orbitals,Nothing) ? Subscript{N}() : orbitals
    spsubscript=isa(spins,NTuple{N,Int}) ? Subscript(spins) : isa(spins,Nothing) ? Subscript{N}() : spins
    orbitals=obsubscript.opattern
    spins=spsubscript.opattern
    nambus==nothing && (nambus=ntuple(i->wildcard,N))
    obsubs=ntuple(i->obsubscript.identifier,N)
    spsubs=ntuple(i->spsubscript.identifier,N)
    return FockCoupling(value,ID(FCID,centers,atoms,orbitals,spins,nambus,obsubs,spsubs),Subscripts(obsubscript),Subscripts(spsubscript))
end

"""
    show(io::IO,fc::FockCoupling)

Show a Fock coupling.
"""
function Base.show(io::IO,fc::FockCoupling)
    @printf io "FockCoupling(value=%s" decimaltostr(fc.value)
    cache=[]
    for attrname in (:centers,:atoms,:orbitals,:spins,:nambus)
        any((attrvalue=getproperty(fc.id,attrname)).≠wildcard) && @printf io ",%s=(%s)" attrname join(attrvalue,",")
    end
    for attrname in (:obsubs,:spsubs)
        all((attrvalue=getproperty(fc.id,attrname)).==wildcard) || all((attrvalue.==constant)) || @printf io ",%s=(%s)" attrname join(attrvalue,",")
    end
    @printf io ")"
end

"""
    repr(fc::FockCoupling) -> String

Get the repr representation of a Fock coupling.
"""
function Base.repr(fc::FockCoupling)
    cache=[]
    for (attrname,abbr) in zip((:atoms,:orbitals,:spins,:nambus),("sl","ob","sp","ph"))
        any((attrvalue=getproperty(fc.id,attrname)).≠wildcard) && push!(cache,@sprintf "%s(%s)" abbr join(attrvalue,':'))
    end
    result=decimaltostr(fc.value)
    length(cache)>0 && (result=@sprintf "%s %s" result join(cache,"⊗"))
    any((centers=fc.id.centers).≠wildcard) && (result=@sprintf "%s@(%s)" result join(centers,'-'))
    obsubs,spsubs=fc.id.obsubs,fc.id.spsubs
    ((all(obsubs.==wildcard) || all(obsubs.==constant)) && (all(spsubs.==wildcard) || all(spsubs.==constant))) || (
            result=@sprintf "%s with %s" result join(((@sprintf "(%s,%s)" obsub spsub) for (obsub,spsub) in zip(obsubs,spsubs))," && ")
            )
    return result
end

"""
    *(fc1::FockCoupling,fc2::FockCoupling) -> FockCoupling

Get the multiplication between two Fock couplings.
"""
function Base.:*(fc1::FockCoupling,fc2::FockCoupling)
    FockCoupling(fc1.value*fc2.value,fc1.id*fc2.id,fc1.obsubscripts*fc2.obsubscripts,fc1.spsubscripts*fc2.spsubscripts)
end

"""
    ⊗(fc1::FockCoupling{N},fc2::FockCoupling{N}) where N -> FockCoupling

Get the direct product between two Fock couplings.
"""
function ⊗(fc1::FockCoupling{N},fc2::FockCoupling{N}) where N
    value,attrvalues,subscripts=fc1.value*fc2.value,[],[]
    wildcards=NTuple{N,Char}(wildcard for i=1:N)
    constants=NTuple{N,Char}(constant for i=1:N)
    v1,v2=fc1.id.centers,fc2.id.centers; t1,t2=v1==wildcards,v2==wildcards
    push!(attrvalues,t1 && t2 ? wildcards : t1 ? v2 : t2 ? v1 : v1==v2 ? v1 : error("⊗ error: dismatched centers ($v1 and $v2)."))
    for (i,attrname) in enumerate((:atoms,:orbitals,:spins,:nambus,:obsubs,:spsubs))
        v1,v2=getproperty(fc1.id,attrname),getproperty(fc2.id,attrname)
        if isa(v1,NTuple{2,Int}) && isa(v2,NTuple{2,Int})
            value=value*delta(v1[2],v2[1])
            push!(attrvalues,(v1[1],v2[2]))
        else
            t1,t2=v1==wildcards,v2==wildcards
            if attrname ∉ (:obsubs,:spsubs)
                @assert (t1 || t2) "⊗ error: dismatched $attrname ($v1 and $v2)."
            else
                f1,f2=v1==constants,v2==constants
                @assert t1 || t2 || (f1 && f2) "⊗ error: dismatched $attrname ($v1 and $v2)."
                attrname==:obsubs && push!(subscripts,t1 ? fc2.obsubscripts : fc1.obsubscripts)
                attrname==:spsubs && push!(subscripts,t1 ? fc2.spsubscripts : fc1.spsubscripts)
            end
            push!(attrvalues,t1 ? v2 : v1)
        end
    end
    return FockCoupling(value,ID(FCID,attrvalues...),subscripts...)
end

"""
    expand(fc::FockCoupling,pid::PID,fock::Fock) -> FCExpand
    expand(fc::FockCoupling,pids::NTuple{R,PID},focks::NTuple{R,Fock}) where R -> FCExpand

Expand a Fock coupling with the given set of point ids and Fock degrees of freedom.
"""
expand(fc::FockCoupling,pid::PID,fock::Fock)=expand(fc,(pid,),(fock,))
function expand(fc::FockCoupling,pids::NTuple{R,PID},focks::NTuple{R,Fock}) where R
    centers=propercenters(fc.id.centers,Val(R))
    rpids=NTuple{rank(fc),eltype(pids)}(pids[centers[i]] for i=1:rank(fc))
    rfocks=NTuple{rank(fc),eltype(focks)}(focks[centers[i]] for i=1:rank(fc))
    nambus=propernambus(fc.id.nambus,NTuple{rank(fc),Int}(rfocks[i].nnambu for i=1:rank(fc)))
    for (i,atom) in enumerate(fc.id.atoms)
        isa(atom,Int) && (atom≠rfocks[i].atom) && return FCExpand(fc.value,rpids,emptyexpands(rank(fc)),emptyexpands(rank(fc)),nambus)
    end
    obsbexpands=expand(fc.obsubscripts,NTuple{rank(fc),Int}(rfocks[i].norbital for i=1:rank(fc)))|>collect
    spsbexpands=expand(fc.spsubscripts,NTuple{rank(fc),Int}(rfocks[i].nspin for i=1:rank(fc)))|>collect
    return FCExpand(fc.value,rpids,obsbexpands,spsbexpands,nambus)
end

defaultcenter(i::Int,n::Int,::Val{1})=1
defaultcenter(i::Int,n::Int,::Val{2})=i<=n/2 ? 1 : 2
defaultcenter(i::Int,n::Int,::Val{R}) where R=error("defaultcenter error: not defined default center for rank-$R bond.")
@generated function propercenters(centers::NTuple{N,Any},::Val{R}) where {N,R}
    exprs=[:(isa(centers[$i],Int) ? (0<centers[$i]<=R ? centers[$i] : error("propercenters error: center out of range.")) : defaultcenter($i,N,Val(R))) for i=1:N]
    return Expr(:tuple,exprs...)
end

defaultnambu(i::Int,nnambu::Int)=i%2==1 ? min(CREATION,nnambu) : ANNIHILATION
@generated function propernambus(nambus::NTuple{N,Any},ranges::NTuple{N,Int}) where N
    exprs=[:(isa(nambus[$i],Int) ? (0<nambus[$i]<=ranges[$i] ? nambus[$i] : error("propernambus error: nambu out of range.")) : defaultnambu($i,ranges[$i])) for i=1:N]
    return Expr(:tuple,exprs...)
end

const empty2expands=Vector{NTuple{2,Int}}()
const empty4expands=Vector{NTuple{4,Int}}()
emptyexpands(n::Int)=n==2 ? empty2expands : n==4 ? empty4expands : Vector{NTuple{n,Int}}()

struct FCExpand{V<:Number,N,S} <: AbstractVectorSpace{Tuple{V,NTuple{N,FIndex{S}}}}
    value::V
    pids::NTuple{N,PID{S}}
    obsbexpands::Vector{NTuple{N,Int}}
    spsbexpands::Vector{NTuple{N,Int}}
    nambus::NTuple{N,Int}
end
IsMultiIndexable(::Type{<:FCExpand})=IsMultiIndexable(true)
MultiIndexOrderStyle(::Type{<:FCExpand})=MultiIndexOrderStyle('C')
dims(fce::FCExpand)=(length(fce.obsbexpands),length(fce.spsbexpands))
@generated function Tuple(inds::NTuple{2,Int},fce::FCExpand{<:Number,N}) where N
    exprs=[:(FIndex(fce.pids[$i],FID(fce.obsbexpands[inds[1]][$i],fce.spsbexpands[inds[2]][$i],fce.nambus[$i]))) for i=1:N]
    return Expr(:tuple,:(fce.value),Expr(:tuple,exprs...))
end

"""
    σ⁰(mode::String;centers::Union{NTuple{2,Int},Nothing}=nothing) -> Couplings{ID{2,FCID},FockCoupling{2,Int,ID{2,FCID}}}

The Pauli matrix σ⁰, which can act on the space of spins("sp"), orbitals("ob"), sublattices("sl") or particle-holes("ph").
"""
function σ⁰(mode::String;centers::Union{NTuple{2,Int},Nothing}=nothing)
    @assert (mode=="sp" || mode=="ob" || mode=="sl" || mode=="ph") "σ⁰ error: not supported mode($mode)."
    attrname=mode=="sp" ? :spins : mode=="ob" ? :orbitals : mode=="sl" ? :atoms : :nambus
    (attrval1,attrval2)=mode=="ph" ? ((ANNIHILATION,CREATION),(CREATION,ANNIHILATION)) : ((1,1),(2,2))
    FockCoupling{2}(1;attrname=>attrval1,:centers=>centers)+FockCoupling{2}(1;attrname=>attrval2,:centers=>centers)
end

"""
    σˣ(mode::String;centers::Union{NTuple{2,Int},Nothing}=nothing) -> Couplings{ID{2,FCID},FockCoupling{2,Int,ID{2,FCID}}}

The Pauli matrix σˣ, which can act on the space of spins("sp"), orbitals("ob"), sublattices("sl") or particle-holes("ph").
"""
function σˣ(mode::String;centers::Union{NTuple{2,Int},Nothing}=nothing)
    @assert (mode=="sp" || mode=="ob" || mode=="sl" || mode=="ph") "σˣ error: not supported mode($mode)."
    attrname=mode=="sp" ? :spins : mode=="ob" ? :orbitals : mode=="sl" ? :atoms : :nambus
    (attrval1,attrval2)=mode=="ph" ? ((ANNIHILATION,ANNIHILATION),(CREATION,CREATION)) : ((1,2),(2,1))
    FockCoupling{2}(1;attrname=>attrval1,:centers=>centers)+FockCoupling{2}(1;attrname=>attrval2,:centers=>centers)
end

"""
    σʸ(mode::String;centers::Union{NTuple{2,Int},Nothing}=nothing) -> Couplings{ID{2,FCID},FockCoupling{2,Int,ID{2,FCID}}}

The Pauli matrix σʸ, which can act on the space of spins("sp"), orbitals("ob"), sublattices("sl") or particle-holes("ph").
"""
function σʸ(mode::String;centers::Union{NTuple{2,Int},Nothing}=nothing)
    @assert (mode=="sp" || mode=="ob" || mode=="sl" || mode=="ph") "σʸ error: not supported mode($mode)."
    attrname=mode=="sp" ? :spins : mode=="ob" ? :orbitals : mode=="sl" ? :atoms : :nambus
    (attrval1,attrval2)=mode=="ph" ? ((ANNIHILATION,ANNIHILATION),(CREATION,CREATION)) : ((1,2),(2,1))
    FockCoupling{2}(1im;attrname=>attrval1,:centers=>centers)+FockCoupling{2}(-1im;attrname=>attrval2,:centers=>centers)
end

"""
    σᶻ(mode::String;centers::Union{NTuple{2,Int},Nothing}=nothing) -> Couplings{ID{2,FCID},FockCoupling{2,Int,ID{2,FCID}}}

The Pauli matrix σᶻ, which can act on the space of spins("sp"), orbitals("ob"), sublattices("sl") or particle-holes("ph").
"""
function σᶻ(mode::String;centers::Union{NTuple{2,Int},Nothing}=nothing)
    @assert (mode=="sp" || mode=="ob" || mode=="sl" || mode=="ph") "σᶻ error: not supported mode($mode)."
    attrname=mode=="sp" ? :spins : mode=="ob" ? :orbitals : mode=="sl" ? :atoms : :nambus
    (attrval1,attrval2)=mode=="ph" ? ((ANNIHILATION,CREATION),(CREATION,ANNIHILATION)) : ((1,1),(2,2))
    FockCoupling{2}(-1;attrname=>attrval1,:centers=>centers)+FockCoupling{2}(1;attrname=>attrval2,:centers=>centers)
end

"""
    σ⁺(mode::String;centers::Union{NTuple{2,Int},Nothing}=nothing) -> Couplings{ID{2,FCID},FockCoupling{2,Int,ID{2,FCID}}}

The Pauli matrix σ⁺, which can act on the space of spins("sp"), orbitals("ob"), sublattices("sl") or particle-holes("ph").
"""
function σ⁺(mode::String;centers::Union{NTuple{2,Int},Nothing}=nothing)
    @assert (mode=="sp" || mode=="ob" || mode=="sl" || mode=="ph") "σᶻ error: not supported mode($mode)."
    attrname=mode=="sp" ? :spins : mode=="ob" ? :orbitals : mode=="sl" ? :atoms : :nambus
    attrval=mode=="ph" ? (CREATION,CREATION) : (2,1)
    Couplings(FockCoupling{2}(1;attrname=>attrval,:centers=>centers))
end

"""
    σ⁻(mode::String;centers::Union{NTuple{2,Int},Nothing}=nothing) -> Couplings{ID{2,FCID},FockCoupling{2,Int,ID{2,FCID}}}

The Pauli matrix σ⁻, which can act on the space of spins("sp"), orbitals("ob"), sublattices("sl") or particle-holes("ph").
"""
function σ⁻(mode::String;centers::Union{NTuple{2,Int},Nothing}=nothing)
    @assert (mode=="sp" || mode=="ob" || mode=="sl" || mode=="ph") "σᶻ error: not supported mode($mode)."
    attrname=mode=="sp" ? :spins : mode=="ob" ? :orbitals : mode=="sl" ? :atoms : :nambus
    attrval=mode=="ph" ? (ANNIHILATION,ANNIHILATION) : (1,2)
    Couplings(FockCoupling{2}(1;attrname=>attrval,:centers=>centers))
end
