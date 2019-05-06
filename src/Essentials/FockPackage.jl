module FockPackage

using LinearAlgebra: dot
using Printf: @printf,@sprintf
using ..Spatials: PID,AbstractBond,Bond,decompose
using ..DegreesOfFreedom: IID,Index,Internal,FilteredAttributes,IDFConfig,Table,OID,Operator,Operators,LaTeX,latexformat,coordpresent
using ..Terms: wildcard,constant,Subscript,Subscripts,subscriptexpr,Coupling,Couplings,@subscript,couplingcenters,Term,TermCouplings,TermAmplitude,TermModulate
using ...Interfaces: id,rank,kind
using ...Prerequisites: Float,delta,decimaltostr
using ...Mathematics.AlgebraOverFields: SimpleID,ID,Element
using ...Mathematics.VectorSpaces: VectorSpace,IsMultiIndexable,MultiIndexOrderStyle

import ..DegreesOfFreedom: script,otype,isHermitian
import ..Terms: couplingcenter,statistics,abbr,termfactor
import ...Interfaces: dims,inds,⊗,⋅,expand,expand!,permute
import ...Mathematics.AlgebraOverFields: rawelement

export ANNIHILATION,CREATION,MAJORANA,FID,FIndex,Fock
export usualfockindextotuple,nambufockindextotuple
export FockOperator,FOperator,BOperator,isnormalordered
export foptdefaultlatex,boptdefaultlatex
export FCID,FockCoupling
export σ⁰,σˣ,σʸ,σᶻ,σ⁺,σ⁻
export @fc_str,@σ⁰_str,@σˣ_str,@σʸ_str,@σᶻ_str,@σ⁺_str,@σ⁻_str
export Onsite,Hopping,Pairing
export hubbard,Hubbard
export interorbitalinterspin,InterOrbitalInterSpin
export interorbitalintraspin,InterOrbitalIntraSpin
export spinflip,SpinFlip
export pairhopping,PairHopping
export Coulomb
export FFockTerm,BFockTerm

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
    MAJORANA

Indicate that the nambu index is MAJORANA.
"""
const MAJORANA=0

"""
    FID <: IID

The Fock id.
"""
struct FID <: IID
    orbital::Int
    spin::Int
    nambu::Int
    FID(orbital::Int,spin::Int,nambu::Int)=(@assert nambu ∈ (0,1,2) "FID error: wrong input nambu($nambu).";new(orbital,spin,nambu))
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
Base.adjoint(fid::FID)=FID(fid.orbital,fid.spin,fid.nambu==0 ? 0 : 3-fid.nambu)

"""
    Fock <: Internal{FID}

The Fock interanl degrees of freedom.
"""
struct Fock <: Internal{FID}
    atom::Int
    norbital::Int
    nspin::Int
    nnambu::Int
    Fock(atom::Int,norbital::Int,nspin::Int,nnambu::Int)=(@assert nnambu ∈ (1,2) "Fock error: wrong input nnambu($nnambu).";new(atom,norbital,nspin,nnambu))
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
    function FIndex(scope::S,site::Int,orbital::Int,spin::Int,nambu::Int) where S
        @assert nambu ∈ (0,1,2) "FIndex error: wrong input nambu($nambu)."
        new{S}(scope,site,orbital,spin,nambu)
    end
end
Base.fieldnames(::Type{<:FIndex})=(:scope,:site,:orbital,:spin,:nambu)

"""
    union(::Type{P},::Type{FID}) where P<:PID

Get the union type of `PID` and `FID`.
"""
Base.union(::Type{P},::Type{FID}) where P<:PID=FIndex{fieldtype(P,:scope)}

"""
    script(oid::OID{<:FIndex},::Val{:site}) -> Int
    script(oid::OID{<:FIndex},::Val{:orbital}) -> Int
    script(oid::OID{<:FIndex},::Val{:spinint}) -> Int
    script(oid::OID{<:FIndex},::Val{:spinsym}) -> String
    script(oid::OID{<:FIndex},::Val{:nambu}) -> String

Get the required script of an Fock oid.
"""
script(oid::OID{<:FIndex},::Val{:site})=oid.index.site
script(oid::OID{<:FIndex},::Val{:orbital})=oid.index.orbital
script(oid::OID{<:FIndex},::Val{:spinint})=oid.index.spin
script(oid::OID{<:FIndex},::Val{:spinsym})=oid.index.spin==1 ? "↓" : oid.index.spin==2 ? "↑" : error("script error: spin($oid.index.spin) not in (1,2).")
script(oid::OID{<:FIndex},::Val{:nambu})=oid.index.nambu==CREATION ? "\\dagger" : ""

"""
    angle(id::OID{<:FIndex},vectors::AbstractVector{<:AbstractVector{Float}},values::AbstractVector{Float}) -> Complex{Float}

Get the twist phase corresponding to a Fock oid.
"""
function Base.angle(id::OID{<:FIndex},vectors::AbstractVector{<:AbstractVector{Float}},values::AbstractVector{Float})
    phase=  length(vectors)==1 ? 2pi*dot(decompose(id.icoord,vectors[1]),values) :
            length(vectors)==2 ? 2pi*dot(decompose(id.icoord,vectors[1],vectors[2]),values) :
            length(vectors)==3 ? 2pi*dot(decompose(id.icoord,vectors[1],vectors[2],vectors[3]),values) :
            error("angle error: not supported number of input basis vectors.")
    id.index.nambu==ANNIHILATION ? phase : id.index.nambu==CREATION ? -phase : error("angle error: not supported Fock index.")
end

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

"""
    FockOperator{V,I<:ID} <: Operator{V,I}

Abstract type for all Fock operators.
"""
abstract type FockOperator{V,I<:ID} <: Operator{V,I} end

"""
    isnormalordered(opt::FockOperator) -> Bool

Judge whether a FockOperator is normal ordered.
"""
function isnormalordered(opt::FockOperator)
    flag=true
    for i=1:rank(opt)
        flag && opt.id[i].index.nambu==ANNIHILATION && (flag=false)
        flag || opt.id[i].index.nambu==CREATION && return false
    end
    return true
end

"""
    FOperator(value,id::ID{<:Tuple{Vararg{OID}}}=ID())

Fermionic Fock operator.
"""
struct FOperator{V,I<:ID} <: FockOperator{V,I}
    value::V
    id::I
    FOperator(value,id::ID{<:Tuple{Vararg{OID}}}=ID())=new{typeof(value),typeof(id)}(value,id)
end

"""
    foptdefaultlatex

The default LaTeX format of the oids of a fermionic operator.
"""
const foptdefaultlatex=LaTeX{(:nambu,),(:site,:orbital,:spinsym)}('c')
latexformat(FOperator,foptdefaultlatex)

"""
    statistics(opt::FOperator) -> Char
    statistics(::Type{<:FOperator}) -> Char

Get the statistics of FOperator.
"""
statistics(opt::FOperator)=opt|>typeof|>statistics
statistics(::Type{<:FOperator})='F'

"""
    rawelement(::Type{<:FOperator})

Get the raw name of a type of FOperator.
"""
rawelement(::Type{<:FOperator})=FOperator

"""
    permute(::Type{<:FOperator},id1::OID{<:FIndex},id2::OID{<:FIndex},::Any=nothing) -> Tuple{Vararg{FOperator}}

Permute two fermionic oid and get the result.
"""
function permute(::Type{<:FOperator},id1::OID{<:FIndex},id2::OID{<:FIndex},::Any=nothing)
    @assert id1.index≠id2.index || id1.rcoord≠id2.rcoord || id1.icoord≠id2.icoord "permute error: permuted ids should not be equal to each other."
    if id1.index'==id2.index && id1.rcoord==id2.rcoord && id1.icoord==id2.icoord
        return (FOperator(1,ID()),FOperator(-1,ID(id2,id1)))
    else
        return (FOperator(-1,ID(id2,id1)),)
    end
end

"""
    *(f1::FOperator,f2::FOperator) -> Union{Nothing,FOperator}

Get the multiplication of two fermionic Fock operators.
"""
function Base.:*(f1::FOperator,f2::FOperator)
    rank(f1)>0 && rank(f2)>0 && f1.id[end]==f2.id[1] && return nothing
    return invoke(*,Tuple{Element,Element},f1,f2)
end

"""
    BOperator(value,id::ID{<:Tuple{Vararg{OID}}}=ID())

Bosonic Fock operator.
"""
struct BOperator{V,I<:ID} <: FockOperator{V,I}
    value::V
    id::I
    BOperator(value,id::ID{<:Tuple{Vararg{OID}}}=ID())=new{typeof(value),typeof(id)}(value,id)
end

"""
    boptdefaultlatex

The default LaTeX format of the oids of a bosonic operator.
"""
const boptdefaultlatex=LaTeX{(:nambu,),(:site,:orbital,:spinsym)}('b')
latexformat(BOperator,boptdefaultlatex)

"""
    statistics(opt::BOperator)
    statistics(::Type{<:BOperator})

Get the statistics of BOperator.
"""
statistics(opt::BOperator)=opt|>typeof|>statistics
statistics(::Type{<:BOperator})='B'

"""
    rawelement(::Type{<:BOperator})

Get the raw name of a type of BOperator.
"""
rawelement(::Type{<:BOperator})=BOperator

"""
    permute(::Type{<:BOperator},id1::OID{<:FIndex},id2::OID{<:FIndex},::Any=nothing) -> Tuple{Vararg{BOperator}}

Permute two bosonic oid and get the result.
"""
function permute(::Type{<:BOperator},id1::OID{<:FIndex},id2::OID{<:FIndex},::Any=nothing)
    @assert id1.index≠id2.index || id1.rcoord≠id2.rcoord || id1.icoord≠id2.icoord "permute error: permuted ids should not be equal to each other."
    if id1.index'==id2.index && id1.rcoord==id2.rcoord && id1.icoord==id2.icoord
        if id1.index.nambu==CREATION
            return (BOperator(1,ID()),BOperator(1,ID(id2,id1)))
        else
            return (BOperator(-1,ID()),BOperator(1,ID(id2,id1)))
        end
    else
        return (BOperator(1,ID(id2,id1)),)
    end
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
    FockCoupling(value,id::ID{<:Tuple{Vararg{FCID}}},obsubscripts::Subscripts,spsubscripts::Subscripts)
    FockCoupling{N}(value1;
                    centers::Union{NTuple{N,Int},Nothing}=nothing,
                    atoms::Union{NTuple{N,Int},Nothing}=nothing,
                    orbitals::Union{NTuple{N,Int},Subscript,Nothing}=nothing,
                    spins::Union{NTuple{N,Int},Subscript,Nothing}=nothing,
                    nambus::Union{NTuple{N,Int},Nothing}=nothing) where N

Fock coupling.
"""
struct FockCoupling{V,I<:ID,OS<:Subscripts,SS<:Subscripts} <: Coupling{V,I}
    value::V
    id::I
    obsubscripts::OS
    spsubscripts::SS
    function FockCoupling(value,id::ID{<:Tuple{Vararg{FCID}}},obsubscripts::Subscripts,spsubscripts::Subscripts)
        new{typeof(value),typeof(id),typeof(obsubscripts),typeof(spsubscripts)}(value,id,obsubscripts,spsubscripts)
    end
end
function FockCoupling{N}(   value=1;
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
    @printf io "FockCoupling{%s}(value=%s" rank(fc) decimaltostr(fc.value)
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
    flag=true
    for (attrname,abbr) in zip((:atoms,:orbitals,:spins,:nambus),("sl","ob","sp","ph"))
        any((attrvalue=getproperty(fc.id,attrname)).≠wildcard) && push!(cache,@sprintf "%s(%s)" abbr join(attrvalue,','))
    end
    result=decimaltostr(fc.value)
    length(cache)>0 && (flag=false;result=@sprintf "%s %s" result join(cache,"⊗"))
    any((centers=fc.id.centers).≠wildcard) && (flag=false;result=@sprintf "%s%s@(%s)" result (length(cache)>0 ? "" : " ") join(centers,','))
    obsubs,spsubs=fc.id.obsubs,fc.id.spsubs
    ((all(obsubs.==wildcard) || all(obsubs.==constant)) && (all(spsubs.==wildcard) || all(spsubs.==constant))) || (
            flag=false;
            result=@sprintf "%s with %s" result join(((@sprintf "(%s,%s)" obsub spsub) for (obsub,spsub) in zip(obsubs,spsubs))," && ")
            )
    flag && (result=@sprintf "%s {%s}" result rank(fc))
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
    ⊗(fc1::FockCoupling,fc2::FockCoupling) -> FockCoupling

Get the direct product between two Fock couplings.
"""
function ⊗(fc1::FockCoupling,fc2::FockCoupling)
    @assert fc1|>rank==fc2|>rank "⊗ error: dismatched rank."
    value,attrvalues,subscripts,N=fc1.value*fc2.value,[],[],rank(fc1)
    wildcards=NTuple{N,Char}(wildcard for i=1:N)
    v1,v2=fc1.id.centers,fc2.id.centers; t1,t2=v1==wildcards,v2==wildcards
    push!(attrvalues,t1 && t2 ? wildcards : t1 ? v2 : t2 ? v1 : v1==v2 ? v1 : error("⊗ error: dismatched centers ($v1 and $v2)."))
    for (i,attrname) in enumerate((:atoms,:orbitals,:spins,:nambus,:obsubs,:spsubs))
        v1,v2=getproperty(fc1.id,attrname),getproperty(fc2.id,attrname)
        t1,t2=v1==wildcards,v2==wildcards
        if attrname ∉ (:obsubs,:spsubs)
            @assert (t1 || t2) "⊗ error: dismatched $attrname ($v1 and $v2)."
        else
            @assert t1 || t2 "⊗ error: dismatched $attrname ($v1 and $v2)."
            attrname==:obsubs && push!(subscripts,t1 ? fc2.obsubscripts : fc1.obsubscripts)
            attrname==:spsubs && push!(subscripts,t1 ? fc2.spsubscripts : fc1.spsubscripts)
        end
        push!(attrvalues,t1 ? v2 : v1)
    end
    return FockCoupling(value,ID(FCID,attrvalues...),subscripts...)
end

"""
    ⋅(fc1::FockCoupling,fc2::FockCoupling) -> FockCoupling

Get the dot product of two rank-2 Fock couplings.

A rank-2 FockCoupling can be considered as a matrix acting on the sublattice, orbital, spin and nambu spaces.
The dot product here is defined as the multiplication between such matrices.
"""
function ⋅(fc1::FockCoupling,fc2::FockCoupling)
    @assert rank(fc1)==rank(fc2)==2 "⋅ error: input fockcouplings must be rank-2."
    attrpairs=[]
    value=fc1.value*fc2.value
    wildcards=(wildcard,wildcard)
    v1,v2=fc1.id.centers,fc2.id.centers
    t1,t2=v1==wildcards,v2==wildcards
    push!(attrpairs,:centers=>(t1 && t2 ? nothing : t1 ? v2 : t2 ? v1 : v1==v2 ? v1 : error("⋅ error: dismatched centers ($v1 and $v2).")))
    for (i,attrname) in enumerate((:atoms,:orbitals,:spins,:nambus))
        v1,v2=getproperty(fc1.id,attrname),getproperty(fc2.id,attrname)
        if isa(v1,NTuple{2,Int}) && isa(v2,NTuple{2,Int})
            value=value*delta(v1[2],v2[1])
            push!(attrpairs,attrname=>(v1[1],v2[2]))
        else
            @assert v1==wildcards && v2==wildcards "⋅ error: dismatched $attrname ($v1 and $v2)."
        end
    end
    return FockCoupling{2}(value;attrpairs...)
end

"""
    rawelement(::Type{<:FockCoupling})

Get the raw name of a type of FockCoupling.
"""
rawelement(::Type{<:FockCoupling})=FockCoupling

"""
    expand(fc::FockCoupling,pid::PID,fock::Fock,kind::Union{Val{K},Nothing}=nothing) where K -> Union{FCExpand,Tuple{}}
    expand(fc::FockCoupling,pids::NTuple{R,PID},focks::NTuple{R,Fock},kind::Union{Val{K},Nothing}=nothing) where {R,K} -> Union{FCExpand,Tuple{}}

Expand a Fock coupling with the given set of point ids and Fock degrees of freedom.
"""
expand(fc::FockCoupling,pid::PID,fock::Fock,kind::Union{Val{K},Nothing}=nothing) where K=expand(fc,(pid,),(fock,),kind)
function expand(fc::FockCoupling,pids::NTuple{R,PID},focks::NTuple{R,Fock},kind::Union{Val{K},Nothing}=nothing) where {R,K}
    centers=couplingcenters(typeof(fc),fc.id.centers,Val(R))
    rpids=NTuple{rank(fc),eltype(pids)}(pids[centers[i]] for i=1:rank(fc))
    rfocks=NTuple{rank(fc),eltype(focks)}(focks[centers[i]] for i=1:rank(fc))
    nambus=fockcouplingnambus(kind,fc.id.nambus,NTuple{rank(fc),Int}(rfocks[i].nnambu for i=1:rank(fc)))
    for (i,atom) in enumerate(fc.id.atoms)
       isa(atom,Int) && (atom≠rfocks[i].atom) && return ()
    end
    obsbexpands=expand(fc.obsubscripts,NTuple{rank(fc),Int}(rfocks[i].norbital for i=1:rank(fc)))|>collect
    spsbexpands=expand(fc.spsubscripts,NTuple{rank(fc),Int}(rfocks[i].nspin for i=1:rank(fc)))|>collect
    return FCExpand(fc.value,rpids,obsbexpands,spsbexpands,nambus)
end
couplingcenter(::Type{<:FockCoupling},i::Int,n::Int,::Val{2})=i<=n/2 ? 1 : 2
@generated function fockcouplingnambus(::Union{Val{K},Nothing},nambus::NTuple{N,Any},ranges::NTuple{N,Int}) where {N,K}
    exprs=[:(isa(nambus[$i],Int) ?
            ((ranges[$i]==1 && nambus[$i]==0) || (ranges[$i]==2 && 0<nambus[$i]<=2) ? nambus[$i] : error("fockcouplingnambus error: nambu out of range.")) :
            (ranges[$i]==1 ? 0 : ($i)%2==1 ? CREATION : ANNIHILATION)
            ) for i=1:N]
    return Expr(:tuple,exprs...)
end
struct FCExpand{V,N,S} <: VectorSpace{Tuple{V,NTuple{N,FIndex{S}}}}
    value::V
    pids::NTuple{N,PID{S}}
    obsbexpands::Vector{NTuple{N,Int}}
    spsbexpands::Vector{NTuple{N,Int}}
    nambus::NTuple{N,Int}
end
IsMultiIndexable(::Type{<:FCExpand})=IsMultiIndexable(true)
MultiIndexOrderStyle(::Type{<:FCExpand})=MultiIndexOrderStyle('C')
dims(fce::FCExpand)=(length(fce.obsbexpands),length(fce.spsbexpands))
@generated function Tuple(inds::NTuple{2,Int},fce::FCExpand{V,N}) where {V,N}
    exprs=[:(FIndex(fce.pids[$i],FID(fce.obsbexpands[inds[1]][$i],fce.spsbexpands[inds[2]][$i],fce.nambus[$i]))) for i=1:N]
    return Expr(:tuple,:(fce.value),Expr(:tuple,exprs...))
end

"""
    σ⁰(mode::String;centers::Union{NTuple{2,Int},Nothing}=nothing) -> Couplings{ID,FockCoupling{Int,ID{<:NTuple{2,FCID}}}}

The Pauli matrix σ⁰, which can act on the space of spins("sp"), orbitals("ob"), sublattices("sl") or particle-holes("ph").
"""
function σ⁰(mode::String;centers::Union{NTuple{2,Int},Nothing}=nothing)
    @assert (mode=="sp" || mode=="ob" || mode=="sl" || mode=="ph") "σ⁰ error: not supported mode($mode)."
    attrname=mode=="sp" ? :spins : mode=="ob" ? :orbitals : mode=="sl" ? :atoms : :nambus
    (attrval1,attrval2)=mode=="ph" ? ((ANNIHILATION,CREATION),(CREATION,ANNIHILATION)) : ((1,1),(2,2))
    FockCoupling{2}(1;attrname=>attrval1,:centers=>centers)+FockCoupling{2}(1;attrname=>attrval2,:centers=>centers)
end

"""
    σˣ(mode::String;centers::Union{NTuple{2,Int},Nothing}=nothing) -> Couplings{ID,FockCoupling{Int,ID{<:NTuple{2,FCID}}}}

The Pauli matrix σˣ, which can act on the space of spins("sp"), orbitals("ob"), sublattices("sl") or particle-holes("ph").
"""
function σˣ(mode::String;centers::Union{NTuple{2,Int},Nothing}=nothing)
    @assert (mode=="sp" || mode=="ob" || mode=="sl" || mode=="ph") "σˣ error: not supported mode($mode)."
    attrname=mode=="sp" ? :spins : mode=="ob" ? :orbitals : mode=="sl" ? :atoms : :nambus
    (attrval1,attrval2)=mode=="ph" ? ((ANNIHILATION,ANNIHILATION),(CREATION,CREATION)) : ((1,2),(2,1))
    FockCoupling{2}(1;attrname=>attrval1,:centers=>centers)+FockCoupling{2}(1;attrname=>attrval2,:centers=>centers)
end

"""
    σʸ(mode::String;centers::Union{NTuple{2,Int},Nothing}=nothing) -> Couplings{ID,FockCoupling{Int,ID{<:NTuple{2,FCID}}}}

The Pauli matrix σʸ, which can act on the space of spins("sp"), orbitals("ob"), sublattices("sl") or particle-holes("ph").
"""
function σʸ(mode::String;centers::Union{NTuple{2,Int},Nothing}=nothing)
    @assert (mode=="sp" || mode=="ob" || mode=="sl" || mode=="ph") "σʸ error: not supported mode($mode)."
    attrname=mode=="sp" ? :spins : mode=="ob" ? :orbitals : mode=="sl" ? :atoms : :nambus
    (attrval1,attrval2)=mode=="ph" ? ((ANNIHILATION,ANNIHILATION),(CREATION,CREATION)) : ((1,2),(2,1))
    FockCoupling{2}(1im;attrname=>attrval1,:centers=>centers)+FockCoupling{2}(-1im;attrname=>attrval2,:centers=>centers)
end

"""
    σᶻ(mode::String;centers::Union{NTuple{2,Int},Nothing}=nothing) -> Couplings{ID,FockCoupling{Int,ID{<:NTuple{2,FCID}}}}

The Pauli matrix σᶻ, which can act on the space of spins("sp"), orbitals("ob"), sublattices("sl") or particle-holes("ph").
"""
function σᶻ(mode::String;centers::Union{NTuple{2,Int},Nothing}=nothing)
    @assert (mode=="sp" || mode=="ob" || mode=="sl" || mode=="ph") "σᶻ error: not supported mode($mode)."
    attrname=mode=="sp" ? :spins : mode=="ob" ? :orbitals : mode=="sl" ? :atoms : :nambus
    (attrval1,attrval2)=mode=="ph" ? ((ANNIHILATION,CREATION),(CREATION,ANNIHILATION)) : ((1,1),(2,2))
    FockCoupling{2}(-1;attrname=>attrval1,:centers=>centers)+FockCoupling{2}(1;attrname=>attrval2,:centers=>centers)
end

"""
    σ⁺(mode::String;centers::Union{NTuple{2,Int},Nothing}=nothing) -> Couplings{ID,FockCoupling{Int,ID{<:NTuple{2,FCID}}}}

The Pauli matrix σ⁺, which can act on the space of spins("sp"), orbitals("ob"), sublattices("sl") or particle-holes("ph").
"""
function σ⁺(mode::String;centers::Union{NTuple{2,Int},Nothing}=nothing)
    @assert (mode=="sp" || mode=="ob" || mode=="sl" || mode=="ph") "σ⁺ error: not supported mode($mode)."
    attrname=mode=="sp" ? :spins : mode=="ob" ? :orbitals : mode=="sl" ? :atoms : :nambus
    attrval=mode=="ph" ? (CREATION,CREATION) : (2,1)
    Couplings(FockCoupling{2}(1;attrname=>attrval,:centers=>centers))
end

"""
    σ⁻(mode::String;centers::Union{NTuple{2,Int},Nothing}=nothing) -> Couplings{ID,FockCoupling{Int,ID{<:NTuple{2,FCID}}}}

The Pauli matrix σ⁻, which can act on the space of spins("sp"), orbitals("ob"), sublattices("sl") or particle-holes("ph").
"""
function σ⁻(mode::String;centers::Union{NTuple{2,Int},Nothing}=nothing)
    @assert (mode=="sp" || mode=="ob" || mode=="sl" || mode=="ph") "σ⁻ error: not supported mode($mode)."
    attrname=mode=="sp" ? :spins : mode=="ob" ? :orbitals : mode=="sl" ? :atoms : :nambus
    attrval=mode=="ph" ? (ANNIHILATION,ANNIHILATION) : (1,2)
    Couplings(FockCoupling{2}(1;attrname=>attrval,:centers=>centers))
end

"""
    fc"..." -> FockCoupling

Construct a FockCoupling from a literal string.
"""
macro fc_str(str)
    ps=split(str," with ")
    conditions=length(ps)==2 ? fcconditions(ps[2]) : nothing
    ps=split(ps[1],' ')
    coeff=eval(Meta.parse(ps[1]))
    if ps[2][1]=='{' && ps[2][3]=='}'
        N=parse(Int,ps[2][2])
        return FockCoupling{N}(coeff)
    else
        attrpairs=[]
        components=split(ps[2],'@')
        centers=length(components)==2 ? couplingcenters(components[2]) : nothing
        push!(attrpairs,:centers=>centers)
        N=centers===nothing ? nothing : length(centers)
        count=0
        if length(components[1])>0
            for component in split(components[1],'⊗')
                attrname,attrvalue=fccomponent(component)
                N===nothing && (N=length(attrvalue))
                if isa(attrvalue,Expr)
                    @assert attrname==:orbitals || attrname==:spins "@fc_str error: wrong input pattern."
                    @assert N==length(attrvalue.args) "@fc_str error: dismatched ranks."
                    count=count+1
                    condition=conditions===nothing ? :nothing : conditions[count]
                    push!(attrpairs,attrname=>eval(subscriptexpr(attrvalue,condition)))
                else
                    @assert N==length(attrvalue) "@fc_str error: dismatched ranks."
                    push!(attrpairs,attrname=>attrvalue)
                end
            end
        end
        return FockCoupling{N}(coeff;attrpairs...)
    end
end
function fcconditions(str::AbstractString)
    conditions=Meta.parse(str)
    conditions=conditions.head==:tuple ? conditions.args : [conditions]
    return [condition=="*" ? :nothing : condition for condition in conditions]
end
function fccomponent(str::AbstractString)
    @assert str[3]=='(' && str[end]==')' "fccomponent error: wrong input pattern."
    attrname=str[1:2]=="sl" ? :atoms : str[1:2]=="ob" ? :orbitals : str[1:2]=="sp" ? :spins : str[1:2]=="ph" ? :nambus : error("fccomponent error: wrong input pattern.")
    expr=Meta.parse(str[3:end])
    attrvalue=isa(expr,Expr) ? (all(isa(arg,Int) for arg in expr.args) ? Tuple(expr.args) : expr) : (expr,)
    return attrname=>attrvalue
end

σᵅsplit(str::AbstractString)=(ps=split(str,'@'); length(ps)==1 ? (ps[1],nothing) : length(ps)==2 ? (ps[1],couplingcenters(ps[2])) : "σᵅsplit error: wrong input pattern.")
σᵅname(mode::AbstractString)=mode=="sp" ? :spins : mode=="ob" ? :orbitals : mode=="sl" ? :atoms : mode=="ph" ? :nambus : error("σᵅname error: wrong input mode.")

"""
    σ⁰"sp"/σ⁰"sp@(c₁,c₂)" -> Couplings
    σ⁰"ob"/σ⁰"ob@(c₁,c₂)" -> Couplings
    σ⁰"sl"/σ⁰"sl@(c₁,c₂)" -> Couplings
    σ⁰"ph"/σ⁰"ph@(c₁,c₂)" -> Couplings

The Pauli matrix σ⁰, which can act on the space of spins("sp"), orbitals("ob"), sublattices("sl") or particle-holes("ph").
"""
macro σ⁰_str(str::String)
    mode,centers=σᵅsplit(str)
    attrname=σᵅname(mode)
    (attrval1,attrval2)=mode=="ph" ? ((ANNIHILATION,CREATION),(CREATION,ANNIHILATION)) : ((1,1),(2,2))
    FockCoupling{2}(1;attrname=>attrval1,:centers=>centers)+FockCoupling{2}(1;attrname=>attrval2,:centers=>centers)
end

"""
    σˣ"sp"/σˣ"sp@(c₁,c₂)" -> Couplings
    σˣ"ob"/σˣ"ob@(c₁,c₂)" -> Couplings
    σˣ"sl"/σˣ"sl@(c₁,c₂)" -> Couplings
    σˣ"ph"/σˣ"ph@(c₁,c₂)" -> Couplings

The Pauli matrix σˣ, which can act on the space of spins("sp"), orbitals("ob"), sublattices("sl") or particle-holes("ph").
"""
macro σˣ_str(str::String)
    mode,centers=σᵅsplit(str)
    attrname=σᵅname(mode)
    (attrval1,attrval2)=mode=="ph" ? ((ANNIHILATION,ANNIHILATION),(CREATION,CREATION)) : ((1,2),(2,1))
    FockCoupling{2}(1;attrname=>attrval1,:centers=>centers)+FockCoupling{2}(1;attrname=>attrval2,:centers=>centers)
end

"""
    σʸ"sp"/σʸ"sp@(c₁,c₂)" -> Couplings
    σʸ"ob"/σʸ"ob@(c₁,c₂)" -> Couplings
    σʸ"sl"/σʸ"sl@(c₁,c₂)" -> Couplings
    σʸ"ph"/σʸ"ph@(c₁,c₂)" -> Couplings

The Pauli matrix σʸ, which can act on the space of spins("sp"), orbitals("ob"), sublattices("sl") or particle-holes("ph").
"""
macro σʸ_str(str::String)
    mode,centers=σᵅsplit(str)
    attrname=σᵅname(mode)
    (attrval1,attrval2)=mode=="ph" ? ((ANNIHILATION,ANNIHILATION),(CREATION,CREATION)) : ((1,2),(2,1))
    FockCoupling{2}(1im;attrname=>attrval1,:centers=>centers)+FockCoupling{2}(-1im;attrname=>attrval2,:centers=>centers)
end

"""
    σᶻ"sp"/σᶻ"sp@(c₁,c₂)" -> Couplings
    σᶻ"ob"/σᶻ"ob@(c₁,c₂)" -> Couplings
    σᶻ"sl"/σᶻ"sl@(c₁,c₂)" -> Couplings
    σᶻ"ph"/σᶻ"ph@(c₁,c₂)" -> Couplings

The Pauli matrix σᶻ, which can act on the space of spins("sp"), orbitals("ob"), sublattices("sl") or particle-holes("ph").
"""
macro σᶻ_str(str::String)
    mode,centers=σᵅsplit(str)
    attrname=σᵅname(mode)
    (attrval1,attrval2)=mode=="ph" ? ((ANNIHILATION,CREATION),(CREATION,ANNIHILATION)) : ((1,1),(2,2))
    FockCoupling{2}(-1;attrname=>attrval1,:centers=>centers)+FockCoupling{2}(1;attrname=>attrval2,:centers=>centers)
end

"""
    σ⁺"sp"/σ⁺"sp@(c₁,c₂)" -> Couplings
    σ⁺"ob"/σ⁺"ob@(c₁,c₂)" -> Couplings
    σ⁺"sl"/σ⁺"sl@(c₁,c₂)" -> Couplings
    σ⁺"ph"/σ⁺"ph@(c₁,c₂)" -> Couplings

The Pauli matrix σ⁺, which can act on the space of spins("sp"), orbitals("ob"), sublattices("sl") or particle-holes("ph").
"""
macro σ⁺_str(str::String)
    mode,centers=σᵅsplit(str)
    attrname=σᵅname(mode)
    attrval=mode=="ph" ? (CREATION,CREATION) : (2,1)
    Couplings(FockCoupling{2}(1;attrname=>attrval,:centers=>centers))
end

"""
    σ⁻"sp"/σ⁻"sp@(c₁,c₂)" -> Couplings
    σ⁻"ob"/σ⁻"ob@(c₁,c₂)" -> Couplings
    σ⁻"sl"/σ⁻"sl@(c₁,c₂)" -> Couplings
    σ⁻"ph"/σ⁻"ph@(c₁,c₂)" -> Couplings

The Pauli matrix σ⁻, which can act on the space of spins("sp"), orbitals("ob"), sublattices("sl") or particle-holes("ph").
"""
macro σ⁻_str(str::String)
    mode,centers=σᵅsplit(str)
    attrname=σᵅname(mode)
    attrval=mode=="ph" ? (ANNIHILATION,ANNIHILATION) : (1,2)
    Couplings(FockCoupling{2}(1;attrname=>attrval,:centers=>centers))
end

"""
    Onsite{ST}( id::Symbol,value::Any;
                couplings::Union{Function,Coupling,Couplings,Nothing}=nothing,
                amplitude::Union{Function,Nothing}=nothing,
                modulate::Union{Function,Bool}=false,
                ) where {ST}

Onsite term.

Type alias for `Term{statistics,:Onsite,2,id,V,Int,<:TermCouplings,<:TermAmplitude,<:Union{TermModulate,Nothing}}`.
"""
const Onsite{statistics,id,V,C<:TermCouplings,A<:TermAmplitude,M<:Union{TermModulate,Nothing}}=Term{statistics,:Onsite,2,id,V,Int,C,A,M}
function Onsite{ST}(id::Symbol,value::Any;
                    couplings::Union{Function,Coupling,Couplings,Nothing}=nothing,
                    amplitude::Union{Function,Nothing}=nothing,
                    modulate::Union{Function,Bool}=false,
                    ) where {ST}
    couplings=TermCouplings(couplings===nothing ? FockCoupling{2}() : couplings)
    Term{ST,:Onsite,2}(id,value,0,couplings=couplings,amplitude=amplitude,modulate=modulate)
end
abbr(::Type{<:Onsite})=:st
isHermitian(::Type{<:Onsite})=nothing

"""
    Hopping{ST}(id::Symbol,value::Any,bondkind::Int=1;
                couplings::Union{Function,Coupling,Couplings,Nothing}=nothing,
                amplitude::Union{Function,Nothing}=nothing,
                modulate::Union{Function,Bool}=false,
                ) where {ST}

Hopping term.

Type alias for `Term{statistics,:Hopping,2,id,V,Int,<:TermCouplings,<:TermAmplitude,<:Union{TermModulate,Nothing}}`.
"""
const Hopping{statistics,id,V,C<:TermCouplings,A<:TermAmplitude,M<:Union{TermModulate,Nothing}}=Term{statistics,:Hopping,2,id,V,Int,C,A,M}
function Hopping{ST}(id::Symbol,value::Any,bondkind::Int=1;
                    couplings::Union{Function,Coupling,Couplings,Nothing}=nothing,
                    amplitude::Union{Function,Nothing}=nothing,
                    modulate::Union{Function,Bool}=false,
                    ) where {ST}
    couplings=TermCouplings(couplings===nothing ? FockCoupling{2}() : couplings)
    @assert bondkind≠0 "Hopping error: input bondkind (neighbor) cannot be 0. Use `Onsite` instead."
    Term{ST,:Hopping,2}(id,value,bondkind,couplings=couplings,amplitude=amplitude,modulate=modulate)
end
abbr(::Type{<:Hopping})=:hp
isHermitian(::Type{<:Hopping})=false

"""
    Pairing{ST}(id::Symbol,value::Any,bondkind::Int=0;
                couplings::Union{Function,Coupling,Couplings,Nothing}=nothing,
                amplitude::Union{Function,Nothing}=nothing,
                modulate::Union{Function,Bool}=false,
                ) where {ST}

Pairing term.

Type alias for `Term{statistics,:Pairing,2,id,V,Int,<:TermCouplings,<:TermAmplitude,<:Union{TermModulate,Nothing}}`.
"""
const Pairing{statistics,id,V,C<:TermCouplings,A<:TermAmplitude,M<:Union{TermModulate,Nothing}}=Term{statistics,:Pairing,2,id,V,Int,C,A,M}
function Pairing{ST}(id::Symbol,value::Any,bondkind::Int=0;
                    couplings::Union{Function,Coupling,Couplings,Nothing}=nothing,
                    amplitude::Union{Function,Nothing}=nothing,
                    modulate::Union{Function,Bool}=false,
                    ) where {ST}
    couplings=TermCouplings(couplings===nothing ? FockCoupling{2}() : couplings)
    Term{ST,:Pairing,2}(id,value,bondkind,couplings=couplings,amplitude=amplitude,modulate=modulate)
end
abbr(::Type{<:Pairing})=:pr
isHermitian(::Type{<:Pairing})=false
function fockcouplingnambus(::Val{:Pairing},nambus::NTuple{2,Any},ranges::NTuple{2,Int})
    @assert ranges==(2,2) "fockcouplingnambus error: ranges for Pairing terms must be (2,2)."
    (   isa(nambus[1],Int) ? (0<nambus[1]<=2 ? nambus[1] : error("fockcouplingnambus error: nambu out of range.")) : ANNIHILATION,
        isa(nambus[2],Int) ? (0<nambus[2]<=2 ? nambus[2] : error("fockcouplingnambus error: nambu out of range.")) : ANNIHILATION
        )
end
function expand!(   operators::Operators,term::Pairing,bond::AbstractBond,config::IDFConfig,
                    table::Union{Table,Nothing}=nothing,
                    half::Bool=false,
                    coord::Union{Val{true},Val{false}}=coordpresent
                    )
    argtypes=Tuple{Operators,Term,AbstractBond,IDFConfig,Union{Table,Nothing},Bool,Union{Val{true},Val{false}}}
    invoke(expand!,argtypes,operators,term,bond,config,table,half,coord)
    isa(bond,Bond) && invoke(expand!,argtypes,operators,term,bond|>reverse,config,table,half,coord)
    return operators
end

const hubbard=FockCoupling{4}(spins=(2,2,1,1),nambus=(CREATION,ANNIHILATION,CREATION,ANNIHILATION))
"""
    Hubbard{ST}(id::Symbol,value::Any;
                amplitude::Union{Function,Nothing}=nothing,
                modulate::Union{Function,Bool}=false,
                ) where {ST}

Hubbard term.

Type alias for `Term{statistics,:Hubbard,4,id,V,Int,<:TermCouplings,<:TermAmplitude,<:Union{TermModulate,Nothing}}`.
"""
const Hubbard{statistics,id,V,C<:TermCouplings,A<:TermAmplitude,M<:Union{TermModulate,Nothing}}=Term{statistics,:Hubbard,4,id,V,Int,C,A,M}
function Hubbard{ST}(id::Symbol,value::Any;
                    amplitude::Union{Function,Nothing}=nothing,
                    modulate::Union{Function,Bool}=false,
                    ) where {ST}
    @assert value==value' "Hubbard error: only real values are allowed."
    Term{ST,:Hubbard,4}(id,value,0,couplings=hubbard,amplitude=amplitude,modulate=modulate)
end
abbr(::Type{<:Hubbard})=:hb
isHermitian(::Type{<:Hubbard})=true

const interorbitalinterspin=FockCoupling{4}(orbitals=(@subscript (α,α,β,β) with α<β),spins=(@subscript (σ₁,σ₁,σ₂,σ₂) with σ₁≠σ₂),nambus=(CREATION,ANNIHILATION,CREATION,ANNIHILATION))
"""
    InterOrbitalInterSpin{ST}(  id::Symbol,value::Any;
                                amplitude::Union{Function,Nothing}=nothing,
                                modulate::Union{Function,Bool}=false,
                                ) where {ST}

Interorbital-interspin term.

Type alias for `Term{statistics,:InterOrbitalInterSpin,4,id,V,Int,<:TermCouplings,<:TermAmplitude,<:Union{TermModulate,Nothing}}`.
"""
const InterOrbitalInterSpin{statistics,id,V,C<:TermCouplings,A<:TermAmplitude,M<:Union{TermModulate,Nothing}}=Term{statistics,:InterOrbitalInterSpin,4,id,V,Int,C,A,M}
function InterOrbitalInterSpin{ST}( id::Symbol,value::Any;
                                    amplitude::Union{Function,Nothing}=nothing,
                                    modulate::Union{Function,Bool}=false,
                                    ) where {ST}
    @assert value==value' "InterOrbitalInterSpin error: only real values are allowed."
    Term{ST,:InterOrbitalInterSpin,4}(id,value,0,couplings=interorbitalinterspin,amplitude=amplitude,modulate=modulate)
end
abbr(::Type{<:InterOrbitalInterSpin})=:nons
isHermitian(::Type{<:InterOrbitalInterSpin})=true

const interorbitalintraspin=FockCoupling{4}(orbitals=(@subscript (α,α,β,β) with α<β),nambus=(CREATION,ANNIHILATION,CREATION,ANNIHILATION))
"""
    InterOrbitalIntraSpin{ST}(  id::Symbol,value::Any;
                                amplitude::Union{Function,Nothing}=nothing,
                                modulate::Union{Function,Bool}=false,
                                ) where {ST}

Interorbital-intraspin term.

Type alias for `Term{statistics,:InterOrbitalIntraSpin,4,id,V,Int,<:TermCouplings,<:TermAmplitude,<:Union{TermModulate,Nothing}}`.
"""
const InterOrbitalIntraSpin{statistics,id,V,C<:TermCouplings,A<:TermAmplitude,M<:Union{TermModulate,Nothing}}=Term{statistics,:InterOrbitalIntraSpin,4,id,V,Int,C,A,M}
function InterOrbitalIntraSpin{ST}( id::Symbol,value::Any;
                                    amplitude::Union{Function,Nothing}=nothing,
                                    modulate::Union{Function,Bool}=false,
                                    ) where {ST}
    @assert value==value' "InterOrbitalIntraSpin error: only real values are allowed."
    Term{ST,:InterOrbitalIntraSpin,4}(id,value,0,couplings=interorbitalintraspin,amplitude=amplitude,modulate=modulate)
end
abbr(::Type{<:InterOrbitalIntraSpin})=:noes
isHermitian(::Type{<:InterOrbitalIntraSpin})=true

const spinflip=FockCoupling{4}(orbitals=(@subscript (α,β,α,β) with α<β),spins=(2,1,1,2),nambus=(CREATION,CREATION,ANNIHILATION,ANNIHILATION))
"""
    SpinFlip{ST}(   id::Symbol,value::Any;
                    amplitude::Union{Function,Nothing}=nothing,
                    modulate::Union{Function,Bool}=false,
                    ) where {ST}

Spin-flip term.

Type alias for `Term{statistics,:SpinFlip,4,id,V,Int,<:TermCouplings,<:TermAmplitude,<:Union{TermModulate,Nothing}}`.
"""
const SpinFlip{statistics,id,V,C<:TermCouplings,A<:TermAmplitude,M<:Union{TermModulate,Nothing}}=Term{statistics,:SpinFlip,4,id,V,Int,C,A,M}
function SpinFlip{ST}(  id::Symbol,value::Any;
                        amplitude::Union{Function,Nothing}=nothing,
                        modulate::Union{Function,Bool}=false,
                        ) where {ST}
    Term{ST,:SpinFlip,4}(id,value,0,couplings=spinflip,amplitude=amplitude,modulate=modulate)
end
abbr(::Type{<:SpinFlip})=:sf
isHermitian(::Type{<:SpinFlip})=false

const pairhopping=FockCoupling{4}(orbitals=(@subscript (α,α,β,β) with α<β),spins=(2,1,1,2),nambus=(CREATION,CREATION,ANNIHILATION,ANNIHILATION))
"""
    PairHopping{ST}(    id::Symbol,value::Any;
                        amplitude::Union{Function,Nothing}=nothing,
                        modulate::Union{Function,Bool}=false,
                        ) where {ST}

Pair-hopping term.

Type alias for `Term{statistics,:PairHopping,4,id,V,Int,<:TermCouplings,<:TermAmplitude,<:Union{TermModulate,Nothing}}`.
"""
const PairHopping{statistics,id,V,C<:TermCouplings,A<:TermAmplitude,M<:Union{TermModulate,Nothing}}=Term{statistics,:PairHopping,4,id,V,Int,C,A,M}
function PairHopping{ST}(   id::Symbol,value::Any;
                            amplitude::Union{Function,Nothing}=nothing,
                            modulate::Union{Function,Bool}=false,
                            ) where {ST}
    Term{ST,:PairHopping,4}(id,value,0,couplings=pairhopping,amplitude=amplitude,modulate=modulate)
end
abbr(::Type{<:PairHopping})=:ph
isHermitian(::Type{<:PairHopping})=false

"""
    Coulomb{ST}(    id::Symbol,value::Any,bondkind::Int=1;
                    couplings::Union{Function,Coupling,Couplings,Nothing}=nothing,
                    amplitude::Union{Function,Nothing}=nothing,
                    modulate::Union{Function,Bool}=false,
                    ) where {ST}

Coulomb term.

Type alias for `Term{statistics,:Coulomb,4,id,V,Int,<:TermCouplings,<:TermAmplitude,<:Union{TermModulate,Nothing}}`.
"""
const Coulomb{statistics,id,V,C<:TermCouplings,A<:TermAmplitude,M<:Union{TermModulate,Nothing}}=Term{statistics,:Coulomb,4,id,V,Int,C,A,M}
function Coulomb{ST}(   id::Symbol,value::Any,bondkind::Int=1;
                        couplings::Union{Function,Coupling,Couplings,Nothing}=nothing,
                        amplitude::Union{Function,Nothing}=nothing,
                        modulate::Union{Function,Bool}=false,
                        ) where {ST}
    couplings=TermCouplings(couplings===nothing ? FockCoupling{2}()*FockCoupling{2}() : couplings)
    @assert bondkind≠0 "Coulomb error: input bondkind (neighbor) cannot be 0. Use `Hubbard/InterOrbitalInterSpin/InterOrbitalIntraSpin/SpinFlip/PairHopping` instead."
    Term{ST,:Coulomb,4}(id,value,bondkind,couplings=couplings,amplitude=amplitude,modulate=modulate)
end
abbr(::Type{<:Coulomb})=:cl
isHermitian(::Type{<:Coulomb})=nothing
termfactor(::Nothing,id::ID{<:NTuple{4,OID}},::Val{:Coulomb})=id[2]'==id[1] && id[4]'==id[3] ? 2 : 1

"""
    FFockTerm

Fermionic Fock term types.
"""
const FFockTerm=Union{Onsite{'F'},Hopping{'F'},Pairing{'F'},Hubbard{'F'},InterOrbitalInterSpin{'F'},InterOrbitalIntraSpin{'F'},SpinFlip{'F'},PairHopping{'F'},Coulomb{'F'}}

"""
    BFockTerm

Bosonic Fock term types.
"""
const BFockTerm=Union{Onsite{'B'},Hopping{'B'},Pairing{'B'},Hubbard{'B'},InterOrbitalInterSpin{'B'},InterOrbitalIntraSpin{'B'},SpinFlip{'B'},PairHopping{'B'},Coulomb{'B'}}

"""
    otype(T::Type{<:FFockTerm},I::Type{<:OID})
    otype(T::Type{<:BFockTerm},I::Type{<:OID})

Get the operator type of a Fock term.
"""
otype(T::Type{<:FFockTerm},I::Type{<:OID})=FOperator{T|>valtype,ID{NTuple{T|>rank,I}}}
otype(T::Type{<:BFockTerm},I::Type{<:OID})=BOperator{T|>valtype,ID{NTuple{T|>rank,I}}}

end # module
