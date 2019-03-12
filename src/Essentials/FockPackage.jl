module FockPackage

using StaticArrays: SVector
using LinearAlgebra: dot
using Printf: @printf,@sprintf
using ..Spatials: PID,AbstractBond,Point,Bond,pidtype,decompose
using ..DegreesOfFreedom: IID,Index,Internal,FilteredAttributes,IDFConfig,Table,OID,Operator,Operators
using ..Terms: wildcard,constant,Subscript,Subscripts,Coupling,Couplings,@subscript,propercenters,Term,TermCouplings,TermAmplitude,TermModulate
using ...Interfaces: rank,dimension,add!
using ...Prerequisites: Float,delta,decimaltostr
using ...Mathematics.AlgebraOverFields: SimpleID,ID
using ...Mathematics.VectorSpaces: VectorSpace,IsMultiIndexable,MultiIndexOrderStyle

import ..DegreesOfFreedom: twist,oidtype,otype
import ..Terms: defaultcenter,statistics,abbr,properfactor
import ...Interfaces: dims,inds,⊗,expand

export ANNIHILATION,CREATION,FID,FIndex,Fock
export usualfockindextotuple,nambufockindextotuple
export FockOperator,FOperator,BOperator,isnormalordered
export FCID,FockCoupling
export σ⁰,σˣ,σʸ,σᶻ,σ⁺,σ⁻
export Onsite,Hopping,Pairing
export hubbard,Hubbard
export interorbitalinterspin,InterOrbitalInterSpin
export interorbitalintraspin,InterOrbitalIntraSpin
export spinflip,SpinFlip
export pairhopping,PairHopping
export Coulomb

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
    twist(id::OID{<:FIndex},vectors::AbstractVector{<:AbstractVector{Float}},values::AbstractVector{Float}) -> Complex{Float}

Get the twist phase corresponding to a Fock oid.
"""
function twist(id::OID{<:FIndex},vectors::AbstractVector{<:AbstractVector{Float}},values::AbstractVector{Float})
    phase=  length(vectors)==1 ? exp(2.0im*pi*dot(decompose(id.icoord,vectors[1]),values)) :
            length(vectors)==2 ? exp(2.0im*pi*dot(decompose(id.icoord,vectors[1],vectors[2]),values)) :
            length(vectors)==3 ? exp(2.0im*pi*dot(decompose(id.icoord,vectors[1],vectors[2],vectors[3]),values)) :
            error("twist error: not supported number of input basis vectors.")
    id.index.nambu==ANNIHILATION ? phase : id.index.nambu==CREATION ? conj(phase) : error("twist error: not supported Fock index.")
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
    oidtype(::Val{:Fock},B::Type{<:AbstractBond},::Type{Nothing})
    oidtype(::Val{:Fock},B::Type{<:AbstractBond},::Type{<:Table})

Get the compatible Fock OID type with an AbstractBond type and a Table/Nothing type.
"""
oidtype(::Val{:Fock},B::Type{<:AbstractBond},::Type{Nothing})=OID{FIndex{fieldtype(B|>pidtype,:scope)},SVector{B|>dimension,Float},SVector{B|>dimension,Float},Nothing}
oidtype(::Val{:Fock},B::Type{<:AbstractBond},::Type{<:Table})=OID{FIndex{fieldtype(B|>pidtype,:scope)},SVector{B|>dimension,Float},SVector{B|>dimension,Float},Int}

"""
    FockOperator{N,V<:Number,I<:ID{<:NTuple{N,OID}}} <: Operator{N,V,I}

Abstract type for all Fock operators.
"""
abstract type FockOperator{N,V<:Number,I<:ID{<:NTuple{N,OID}}} <: Operator{N,V,I} end

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
    FOperator(value::Number,id::ID{<:NTuple{N,OID}}) where N

Fermionic Fock operator.
"""
struct FOperator{N,V<:Number,I<:ID{<:NTuple{N,OID}}} <: FockOperator{N,V,I}
    value::V
    id::I
    FOperator(value::Number,id::ID{<:NTuple{N,OID}}) where N=new{N,typeof(value),typeof(id)}(value,id)
end

"""
    statistics(opt::FOperator) -> Char
    statistics(::Type{<:FOperator}) -> Char

Get the statistics of FOperator.
"""
statistics(opt::FOperator)=opt|>typeof|>statistics
statistics(::Type{<:FOperator})='F'

"""
    BOperator(value::Number,id::ID{<:NTuple{N,OID}}) where N

Bosonic Fock operator.
"""
struct BOperator{N,V<:Number,I<:ID{<:NTuple{N,OID}}} <: FockOperator{N,V,I}
    value::V
    id::I
    BOperator(value::Number,id::ID{<:NTuple{N,OID}}) where N=new{N,typeof(value),typeof(id)}(value,id)
end

"""
    statistics(opt::BOperator)
    statistics(::Type{<:BOperator})

Get the statistics of BOperator.
"""
statistics(opt::BOperator)=opt|>typeof|>statistics
statistics(::Type{<:BOperator})='B'

"""
    otype(::Val{:Fock},O::Type{<:Term{'F'}},B::Type{<:AbstractBond},T::Type{<:Union{Nothing,Table}})
    otype(::Val{:Fock},O::Type{<:Term{'B'}},B::Type{<:AbstractBond},T::Type{<:Union{Nothing,Table}})

Get the compatible Fock operator type with a Term type, an AbstractBond type and a Table/Nothing type.
"""
otype(::Val{:Fock},O::Type{<:Term{'F'}},B::Type{<:AbstractBond},T::Type{<:Union{Nothing,Table}})=FOperator{O|>rank,O|>valtype,ID{NTuple{O|>rank,oidtype(Val(:Fock),B,T)}}}
otype(::Val{:Fock},O::Type{<:Term{'B'}},B::Type{<:AbstractBond},T::Type{<:Union{Nothing,Table}})=BOperator{O|>rank,O|>valtype,ID{NTuple{O|>rank,oidtype(Val(:Fock),B,T)}}}

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
    FockCoupling(value::Number,id::ID{<:NTuple{N,FCID}},obsubscripts::Subscripts,spsubscripts::Subscripts) where N
    FockCoupling{N}(value::Number=1;
                    centers::Union{NTuple{N,Int},Nothing}=nothing,
                    atoms::Union{NTuple{N,Int},Nothing}=nothing,
                    orbitals::Union{NTuple{N,Int},Subscript,Nothing}=nothing,
                    spins::Union{NTuple{N,Int},Subscript,Nothing}=nothing,
                    nambus::Union{NTuple{N,Int},Nothing}=nothing) where N

Fock coupling.
"""
struct FockCoupling{N,V<:Number,I<:ID{<:NTuple{N,FCID}},OS<:Subscripts,SS<:Subscripts} <: Coupling{N,V,I}
    value::V
    id::I
    obsubscripts::OS
    spsubscripts::SS
    function FockCoupling(value::Number,id::ID{<:NTuple{N,FCID}},obsubscripts::Subscripts,spsubscripts::Subscripts) where N
        new{N,value|>typeof,id|>typeof,obsubscripts|>typeof,spsubscripts|>typeof}(value,id,obsubscripts,spsubscripts)
    end
end
function FockCoupling{N}(   value::Number=1;
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
    expand(fc::FockCoupling,pid::PID,fock::Fock,species::Union{Val{S},Nothing}=nothing) where S -> Union{FCExpand,Tuple{}}
    expand(fc::FockCoupling,pids::NTuple{R,PID},focks::NTuple{R,Fock},species::Union{Val{S},Nothing}=nothing) where {R,S} -> Union{FCExpand,Tuple{}}

Expand a Fock coupling with the given set of point ids and Fock degrees of freedom.
"""
expand(fc::FockCoupling,pid::PID,fock::Fock,species::Union{Val{S},Nothing}=nothing) where S=expand(fc,(pid,),(fock,),species)
function expand(fc::FockCoupling,pids::NTuple{R,PID},focks::NTuple{R,Fock},species::Union{Val{S},Nothing}=nothing) where {R,S}
    centers=propercenters(typeof(fc),fc.id.centers,Val(R))
    rpids=NTuple{rank(fc),eltype(pids)}(pids[centers[i]] for i=1:rank(fc))
    rfocks=NTuple{rank(fc),eltype(focks)}(focks[centers[i]] for i=1:rank(fc))
    nambus=propernambus(species,fc.id.nambus,NTuple{rank(fc),Int}(rfocks[i].nnambu for i=1:rank(fc)))
    for (i,atom) in enumerate(fc.id.atoms)
       isa(atom,Int) && (atom≠rfocks[i].atom) && return ()
    end
    obsbexpands=expand(fc.obsubscripts,NTuple{rank(fc),Int}(rfocks[i].norbital for i=1:rank(fc)))|>collect
    spsbexpands=expand(fc.spsubscripts,NTuple{rank(fc),Int}(rfocks[i].nspin for i=1:rank(fc)))|>collect
    return FCExpand(fc.value,rpids,obsbexpands,spsbexpands,nambus)
end
defaultcenter(::Type{<:FockCoupling},i::Int,n::Int,::Val{2})=i<=n/2 ? 1 : 2
@generated function propernambus(::Union{Val{S},Nothing},nambus::NTuple{N,Any},ranges::NTuple{N,Int}) where {N,S}
    exprs=[:(isa(nambus[$i],Int) ?
            ((ranges[$i]==1 && nambus[$i]==0) || (ranges[$i]==2 && 0<nambus[$i]<=2) ? nambus[$i] : error("propernambus error: nambu out of range.")) :
            (ranges[$i]==1 ? 0 : ($i)%2==1 ? CREATION : ANNIHILATION)
            ) for i=1:N]
    return Expr(:tuple,exprs...)
end
struct FCExpand{V<:Number,N,S} <: VectorSpace{Tuple{V,NTuple{N,FIndex{S}}}}
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
    σ⁰(mode::String;centers::Union{NTuple{2,Int},Nothing}=nothing) -> Couplings{ID{<:NTuple{2,FCID}},FockCoupling{2,Int,ID{<:NTuple{2,FCID}}}}

The Pauli matrix σ⁰, which can act on the space of spins("sp"), orbitals("ob"), sublattices("sl") or particle-holes("ph").
"""
function σ⁰(mode::String;centers::Union{NTuple{2,Int},Nothing}=nothing)
    @assert (mode=="sp" || mode=="ob" || mode=="sl" || mode=="ph") "σ⁰ error: not supported mode($mode)."
    attrname=mode=="sp" ? :spins : mode=="ob" ? :orbitals : mode=="sl" ? :atoms : :nambus
    (attrval1,attrval2)=mode=="ph" ? ((ANNIHILATION,CREATION),(CREATION,ANNIHILATION)) : ((1,1),(2,2))
    FockCoupling{2}(1;attrname=>attrval1,:centers=>centers)+FockCoupling{2}(1;attrname=>attrval2,:centers=>centers)
end

"""
    σˣ(mode::String;centers::Union{NTuple{2,Int},Nothing}=nothing) -> Couplings{ID{<:NTuple{2,FCID}},FockCoupling{2,Int,ID{<:NTuple{2,FCID}}}}

The Pauli matrix σˣ, which can act on the space of spins("sp"), orbitals("ob"), sublattices("sl") or particle-holes("ph").
"""
function σˣ(mode::String;centers::Union{NTuple{2,Int},Nothing}=nothing)
    @assert (mode=="sp" || mode=="ob" || mode=="sl" || mode=="ph") "σˣ error: not supported mode($mode)."
    attrname=mode=="sp" ? :spins : mode=="ob" ? :orbitals : mode=="sl" ? :atoms : :nambus
    (attrval1,attrval2)=mode=="ph" ? ((ANNIHILATION,ANNIHILATION),(CREATION,CREATION)) : ((1,2),(2,1))
    FockCoupling{2}(1;attrname=>attrval1,:centers=>centers)+FockCoupling{2}(1;attrname=>attrval2,:centers=>centers)
end

"""
    σʸ(mode::String;centers::Union{NTuple{2,Int},Nothing}=nothing) -> Couplings{ID{<:NTuple{2,FCID}},FockCoupling{2,Int,ID{<:NTuple{2,FCID}}}}

The Pauli matrix σʸ, which can act on the space of spins("sp"), orbitals("ob"), sublattices("sl") or particle-holes("ph").
"""
function σʸ(mode::String;centers::Union{NTuple{2,Int},Nothing}=nothing)
    @assert (mode=="sp" || mode=="ob" || mode=="sl" || mode=="ph") "σʸ error: not supported mode($mode)."
    attrname=mode=="sp" ? :spins : mode=="ob" ? :orbitals : mode=="sl" ? :atoms : :nambus
    (attrval1,attrval2)=mode=="ph" ? ((ANNIHILATION,ANNIHILATION),(CREATION,CREATION)) : ((1,2),(2,1))
    FockCoupling{2}(1im;attrname=>attrval1,:centers=>centers)+FockCoupling{2}(-1im;attrname=>attrval2,:centers=>centers)
end

"""
    σᶻ(mode::String;centers::Union{NTuple{2,Int},Nothing}=nothing) -> Couplings{ID{<:NTuple{2,FCID}},FockCoupling{2,Int,ID{<:NTuple{2,FCID}}}}

The Pauli matrix σᶻ, which can act on the space of spins("sp"), orbitals("ob"), sublattices("sl") or particle-holes("ph").
"""
function σᶻ(mode::String;centers::Union{NTuple{2,Int},Nothing}=nothing)
    @assert (mode=="sp" || mode=="ob" || mode=="sl" || mode=="ph") "σᶻ error: not supported mode($mode)."
    attrname=mode=="sp" ? :spins : mode=="ob" ? :orbitals : mode=="sl" ? :atoms : :nambus
    (attrval1,attrval2)=mode=="ph" ? ((ANNIHILATION,CREATION),(CREATION,ANNIHILATION)) : ((1,1),(2,2))
    FockCoupling{2}(-1;attrname=>attrval1,:centers=>centers)+FockCoupling{2}(1;attrname=>attrval2,:centers=>centers)
end

"""
    σ⁺(mode::String;centers::Union{NTuple{2,Int},Nothing}=nothing) -> Couplings{ID{<:NTuple{2,FCID}},FockCoupling{2,Int,ID{<:NTuple{2,FCID}}}}

The Pauli matrix σ⁺, which can act on the space of spins("sp"), orbitals("ob"), sublattices("sl") or particle-holes("ph").
"""
function σ⁺(mode::String;centers::Union{NTuple{2,Int},Nothing}=nothing)
    @assert (mode=="sp" || mode=="ob" || mode=="sl" || mode=="ph") "σᶻ error: not supported mode($mode)."
    attrname=mode=="sp" ? :spins : mode=="ob" ? :orbitals : mode=="sl" ? :atoms : :nambus
    attrval=mode=="ph" ? (CREATION,CREATION) : (2,1)
    Couplings(FockCoupling{2}(1;attrname=>attrval,:centers=>centers))
end

"""
    σ⁻(mode::String;centers::Union{NTuple{2,Int},Nothing}=nothing) -> Couplings{ID{<:NTuple{2,FCID}},FockCoupling{2,Int,ID{<:NTuple{2,FCID}}}}

The Pauli matrix σ⁻, which can act on the space of spins("sp"), orbitals("ob"), sublattices("sl") or particle-holes("ph").
"""
function σ⁻(mode::String;centers::Union{NTuple{2,Int},Nothing}=nothing)
    @assert (mode=="sp" || mode=="ob" || mode=="sl" || mode=="ph") "σᶻ error: not supported mode($mode)."
    attrname=mode=="sp" ? :spins : mode=="ob" ? :orbitals : mode=="sl" ? :atoms : :nambus
    attrval=mode=="ph" ? (ANNIHILATION,ANNIHILATION) : (1,2)
    Couplings(FockCoupling{2}(1;attrname=>attrval,:centers=>centers))
end

"""
    Onsite{ST}( id::Symbol,value::Number;
                couplings::Union{Tuple{<:Tuple{Vararg{Couplings}},<:Function},Coupling,Couplings,Nothing}=nothing,
                amplitude::Union{Function,Nothing}=nothing,
                modulate::Union{Function,Bool}=false,
                factor::Number=1
                ) where {ST}

Onsite term.

Type alias for `Term{Statistics,:Onsite,<:Number,Int,<:TermCouplings,<:TermAmplitude,<:Union{TermModulate,Nothing}}`.
"""
const Onsite{Statistics,V<:Number,C<:TermCouplings,A<:TermAmplitude,M<:Union{TermModulate,Nothing}}=Term{Statistics,:Onsite,V,Int,C,A,M}
function Onsite{ST}(id::Symbol,value::Number;
                    couplings::Union{Tuple{<:Tuple{Vararg{Couplings}},<:Function},Coupling,Couplings,Nothing}=nothing,
                    amplitude::Union{Function,Nothing}=nothing,
                    modulate::Union{Function,Bool}=false,
                    factor::Number=1
                    ) where {ST}
    couplings=TermCouplings(couplings===nothing ? FockCoupling{2}() : couplings)
    @assert rank(couplings)==2 "Onsite error: input couplings must be rank-2."
    Term{ST,:Onsite}(id,value,0,couplings=couplings,amplitude=amplitude,modulate=modulate,factor=factor)
end
abbr(::Type{<:Onsite})=:st
function expand(term::Onsite,bond::AbstractBond,config::IDFConfig,table::Union{Table,Nothing}=nothing,half::Bool=false)
    expand(otype(Val(:Fock),typeof(term),typeof(bond),typeof(table)),term,bond,config,table,half ? nothing : false)
end

"""
    Hopping{ST}(id::Symbol,value::Number;
                neighbor::Int=1,
                couplings::Union{Tuple{<:Tuple{Vararg{Couplings}},<:Function},Coupling,Couplings,Nothing}=nothing,
                amplitude::Union{Function,Nothing}=nothing,
                modulate::Union{Function,Bool}=false,
                factor::Number=1
                ) where {ST}

Hopping term.

Type alias for `Term{Statistics,:Hopping,<:Number,Int,<:TermCouplings,<:TermAmplitude,<:Union{TermModulate,Nothing}}`.
"""
const Hopping{Statistics,V<:Number,C<:TermCouplings,A<:TermAmplitude,M<:Union{TermModulate,Nothing}}=Term{Statistics,:Hopping,V,Int,C,A,M}
function Hopping{ST}(id::Symbol,value::Number;
                    neighbor::Int=1,
                    couplings::Union{Tuple{<:Tuple{Vararg{Couplings}},<:Function},Coupling,Couplings,Nothing}=nothing,
                    amplitude::Union{Function,Nothing}=nothing,
                    modulate::Union{Function,Bool}=false,
                    factor::Number=1
                    ) where {ST}
    couplings=TermCouplings(couplings===nothing ? FockCoupling{2}() : couplings)
    @assert rank(couplings)==2 "Hopping error: input couplings must be rank-2."
    @assert neighbor≠0 "Hopping error: input neighbor cannot be 0. Use `Onsite` instead."
    Term{ST,:Hopping}(id,value,neighbor,couplings=couplings,amplitude=amplitude,modulate=modulate,factor=factor)
end
abbr(::Type{<:Hopping})=:hp
function expand(term::Hopping,bond::AbstractBond,config::IDFConfig,table::Union{Table,Nothing}=nothing,half::Bool=false)
    result=expand(otype(Val(:Fock),typeof(term),typeof(bond),typeof(table)),term,bond,config,table,false)
    half || length(result)==0 || for opt in result|>values|>collect add!(result,opt') end
    return result
end

"""
    Pairing{ST}(id::Symbol,value::Number;
                neighbor::Int=0,
                couplings::Union{Tuple{<:Tuple{Vararg{Couplings}},<:Function},Coupling,Couplings,Nothing}=nothing,
                amplitude::Union{Function,Nothing}=nothing,
                modulate::Union{Function,Bool}=false,
                factor::Number=1
                ) where {ST}

Pairing term.

Type alias for `Term{Statistics,:Pairing,<:Number,Int,<:TermCouplings,<:TermAmplitude,<:Union{TermModulate,Nothing}}`.
"""
const Pairing{Statistics,V<:Number,C<:TermCouplings,A<:TermAmplitude,M<:Union{TermModulate,Nothing}}=Term{Statistics,:Pairing,V,Int,C,A,M}
function Pairing{ST}(id::Symbol,value::Number;
                    neighbor::Int=0,
                    couplings::Union{Tuple{<:Tuple{Vararg{Couplings}},<:Function},Coupling,Couplings,Nothing}=nothing,
                    amplitude::Union{Function,Nothing}=nothing,
                    modulate::Union{Function,Bool}=false,
                    factor::Number=1
                    ) where {ST}
    couplings=TermCouplings(couplings===nothing ? FockCoupling{2}() : couplings)
    @assert rank(couplings)==2 "Pairing error: input couplings must be rank-2."
    Term{ST,:Pairing}(id,value,neighbor,couplings=couplings,amplitude=amplitude,modulate=modulate,factor=factor)
end
abbr(::Type{<:Pairing})=:pr
function propernambus(::Val{:Pairing},nambus::NTuple{2,Any},ranges::NTuple{2,Int})
    @assert ranges==(2,2) "propernambus error: ranges for Pairing terms must be (2,2)."
    (   isa(nambus[1],Int) ? (0<nambus[1]<=2 ? nambus[1] : error("propernambus error: nambu out of range.")) : ANNIHILATION,
        isa(nambus[2],Int) ? (0<nambus[2]<=2 ? nambus[2] : error("propernambus error: nambu out of range.")) : ANNIHILATION
        )
end
function expand(term::Pairing,bond::AbstractBond,config::IDFConfig,table::Union{Table,Nothing}=nothing,half::Bool=false)
    result=expand(otype(Val(:Fock),typeof(term),typeof(bond),typeof(table)),term,bond,config,table,false)
    isa(result,Operators) && isa(bond,Bond) && add!(result,expand(otype(Val(:Fock),typeof(term),typeof(bond),typeof(table)),term,bond|>reverse,config,table,false))
    half || length(result)==0 || for opt in result|>values|>collect add!(result,opt') end
    return result
end

const hubbard=FockCoupling{4}(spins=(2,2,1,1),nambus=(CREATION,ANNIHILATION,CREATION,ANNIHILATION))
"""
    Hubbard{ST}(id::Symbol,value::Real;
                amplitude::Union{Function,Nothing}=nothing,
                modulate::Union{Function,Bool}=false,
                factor::Number=1
                ) where {ST}

Hubbard term.

Type alias for `Term{Statistics,:Hubbard,<:Number,Int,<:TermCouplings,<:TermAmplitude,<:Union{TermModulate,Nothing}}`.
"""
const Hubbard{Statistics,V<:Number,C<:TermCouplings,A<:TermAmplitude,M<:Union{TermModulate,Nothing}}=Term{Statistics,:Hubbard,V,Int,C,A,M}
function Hubbard{ST}(id::Symbol,value::Real;
                    amplitude::Union{Function,Nothing}=nothing,
                    modulate::Union{Function,Bool}=false,
                    factor::Number=1
                    ) where {ST}
    Term{ST,:Hubbard}(id,value,0,couplings=hubbard,amplitude=amplitude,modulate=modulate,factor=factor)
end
abbr(::Type{<:Hubbard})=:hb
function expand(term::Hubbard,bond::AbstractBond,config::IDFConfig,table::Union{Table,Nothing}=nothing,half::Bool=false)
    expand(otype(Val(:Fock),typeof(term),typeof(bond),typeof(table)),term,bond,config,table,half)
end

const interorbitalinterspin=FockCoupling{4}(orbitals=(@subscript (α,β)=>(α,α,β,β) with α<β),
                                            spins=(@subscript (σ₁,σ₂)=>(σ₁,σ₁,σ₂,σ₂) with σ₁≠σ₂),
                                            nambus=(CREATION,ANNIHILATION,CREATION,ANNIHILATION)
                                            )
"""
    InterOrbitalInterSpin{ST}(  id::Symbol,value::Real;
                                amplitude::Union{Function,Nothing}=nothing,
                                modulate::Union{Function,Bool}=false,
                                factor::Number=1
                                ) where {ST}

Interorbital-interspin term.

Type alias for `Term{Statistics,:InterOrbitalInterSpin,<:Number,Int,<:TermCouplings,<:TermAmplitude,<:Union{TermModulate,Nothing}}`.
"""
const InterOrbitalInterSpin{Statistics,V<:Number,C<:TermCouplings,A<:TermAmplitude,M<:Union{TermModulate,Nothing}}=Term{Statistics,:InterOrbitalInterSpin,V,Int,C,A,M}
function InterOrbitalInterSpin{ST}( id::Symbol,value::Real;
                                    amplitude::Union{Function,Nothing}=nothing,
                                    modulate::Union{Function,Bool}=false,
                                    factor::Number=1
                                    ) where {ST}
    Term{ST,:InterOrbitalInterSpin}(id,value,0,couplings=interorbitalinterspin,amplitude=amplitude,modulate=modulate,factor=factor)
end
abbr(::Type{<:InterOrbitalInterSpin})=:nons
function expand(term::InterOrbitalInterSpin,bond::AbstractBond,config::IDFConfig,table::Union{Table,Nothing}=nothing,half::Bool=false)
    expand(otype(Val(:Fock),typeof(term),typeof(bond),typeof(table)),term,bond,config,table,half)
end

const interorbitalintraspin=FockCoupling{4}(orbitals=(@subscript (α,β)=>(α,α,β,β) with α<β),
                                            nambus=(CREATION,ANNIHILATION,CREATION,ANNIHILATION)
                                            )
"""
    InterOrbitalIntraSpin{ST}(  id::Symbol,value::Real;
                                amplitude::Union{Function,Nothing}=nothing,
                                modulate::Union{Function,Bool}=false,
                                factor::Number=1
                                ) where {ST}

Interorbital-intraspin term.

Type alias for `Term{Statistics,:InterOrbitalIntraSpin,<:Number,Int,<:TermCouplings,<:TermAmplitude,<:Union{TermModulate,Nothing}}`.
"""
const InterOrbitalIntraSpin{Statistics,V<:Number,C<:TermCouplings,A<:TermAmplitude,M<:Union{TermModulate,Nothing}}=Term{Statistics,:InterOrbitalIntraSpin,V,Int,C,A,M}
function InterOrbitalIntraSpin{ST}( id::Symbol,value::Real;
                                    amplitude::Union{Function,Nothing}=nothing,
                                    modulate::Union{Function,Bool}=false,
                                    factor::Number=1
                                    ) where {ST}
    Term{ST,:InterOrbitalIntraSpin}(id,value,0,couplings=interorbitalintraspin,amplitude=amplitude,modulate=modulate,factor=factor)
end
abbr(::Type{<:InterOrbitalIntraSpin})=:noes
function expand(term::InterOrbitalIntraSpin,bond::AbstractBond,config::IDFConfig,table::Union{Table,Nothing}=nothing,half::Bool=false)
    expand(otype(Val(:Fock),typeof(term),typeof(bond),typeof(table)),term,bond,config,table,half)
end

const spinflip=FockCoupling{4}( orbitals=(@subscript (α,β)=>(α,β,α,β) with α<β),
                                spins=(2,1,1,2),
                                nambus=(CREATION,CREATION,ANNIHILATION,ANNIHILATION)
                                )
"""
    SpinFlip{ST}(   id::Symbol,value::Real;
                    amplitude::Union{Function,Nothing}=nothing,
                    modulate::Union{Function,Bool}=false,
                    factor::Number=1
                    ) where {ST}

Spin-flip term.

Type alias for `Term{Statistics,:SpinFlip,<:Number,Int,<:TermCouplings,<:TermAmplitude,<:Union{TermModulate,Nothing}}`.
"""
const SpinFlip{Statistics,V<:Number,C<:TermCouplings,A<:TermAmplitude,M<:Union{TermModulate,Nothing}}=Term{Statistics,:SpinFlip,V,Int,C,A,M}
function SpinFlip{ST}(  id::Symbol,value::Real;
                        amplitude::Union{Function,Nothing}=nothing,
                        modulate::Union{Function,Bool}=false,
                        factor::Number=1
                        ) where {ST}
    Term{ST,:SpinFlip}(id,value,0,couplings=spinflip,amplitude=amplitude,modulate=modulate,factor=factor)
end
abbr(::Type{<:SpinFlip})=:sf
function expand(term::SpinFlip,bond::AbstractBond,config::IDFConfig,table::Union{Table,Nothing}=nothing,half::Bool=false)
    result=expand(otype(Val(:Fock),typeof(term),typeof(bond),typeof(table)),term,bond,config,table,false)
    half || length(result)==0 || for opt in result|>values|>collect add!(result,opt') end
    return result
end

const pairhopping=FockCoupling{4}(  orbitals=(@subscript (α,β)=>(α,α,β,β) with α<β),
                                    spins=(2,1,1,2),
                                    nambus=(CREATION,CREATION,ANNIHILATION,ANNIHILATION)
                                    )
"""
    PairHopping{ST}(    id::Symbol,value::Real;
                        amplitude::Union{Function,Nothing}=nothing,
                        modulate::Union{Function,Bool}=false,
                        factor::Number=1
                        ) where {ST}

Pair-hopping term.

Type alias for `Term{Statistics,:PairHopping,<:Number,Int,<:TermCouplings,<:TermAmplitude,<:Union{TermModulate,Nothing}}`.
"""
const PairHopping{Statistics,V<:Number,C<:TermCouplings,A<:TermAmplitude,M<:Union{TermModulate,Nothing}}=Term{Statistics,:PairHopping,V,Int,C,A,M}
function PairHopping{ST}(   id::Symbol,value::Real;
                            amplitude::Union{Function,Nothing}=nothing,
                            modulate::Union{Function,Bool}=false,
                            factor::Number=1
                            ) where {ST}
    Term{ST,:PairHopping}(id,value,0,couplings=pairhopping,amplitude=amplitude,modulate=modulate,factor=factor)
end
abbr(::Type{<:PairHopping})=:ph
function expand(term::PairHopping,bond::AbstractBond,config::IDFConfig,table::Union{Table,Nothing}=nothing,half::Bool=false)
    result=expand(otype(Val(:Fock),typeof(term),typeof(bond),typeof(table)),term,bond,config,table,false)
    half || length(result)==0 || for opt in result|>values|>collect add!(result,opt') end
    return result
end

"""
    Coulomb{ST}(    id::Symbol,value::Number;
                    neighbor::Int=1,
                    couplings::Union{Tuple{<:Tuple{Vararg{Couplings}},<:Function},Coupling,Couplings,Nothing}=nothing,
                    amplitude::Union{Function,Nothing}=nothing,
                    modulate::Union{Function,Bool}=false,
                    factor::Number=1
                    ) where {ST}

Coulomb term.

Type alias for `Term{Statistics,:Coulomb,<:Number,Int,<:TermCouplings,<:TermAmplitude,<:Union{TermModulate,Nothing}}`.
"""
const Coulomb{Statistics,V<:Number,C<:TermCouplings,A<:TermAmplitude,M<:Union{TermModulate,Nothing}}=Term{Statistics,:Coulomb,V,Int,C,A,M}
function Coulomb{ST}(   id::Symbol,value::Number;
                        neighbor::Int=1,
                        couplings::Union{Tuple{<:Tuple{Vararg{Couplings}},<:Function},Coupling,Couplings,Nothing}=nothing,
                        amplitude::Union{Function,Nothing}=nothing,
                        modulate::Union{Function,Bool}=false,
                        factor::Number=1
                        ) where {ST}
    couplings=TermCouplings(couplings===nothing ? FockCoupling{2}()*FockCoupling{2}() : couplings)
    @assert rank(couplings)==4 "Coulomb error: input couplings must be rank-4."
    @assert neighbor≠0 "Coulomb error: input neighbor cannot be 0. Use `Hubbard/InterOrbitalInterSpin/InterOrbitalIntraSpin/SpinFlip/PairHopping` instead."
    Term{ST,:Coulomb}(id,value,neighbor,couplings=couplings,amplitude=amplitude,modulate=modulate,factor=factor)
end
abbr(::Type{<:Coulomb})=:cl
properfactor(::Nothing,id::ID{<:NTuple{4,OID}},::Val{:Coulomb})=id[2]'==id[1] && id[4]'==id[3] ? 0.5 : 1.0
function expand(term::Coulomb,bond::AbstractBond,config::IDFConfig,table::Union{Table,Nothing}=nothing,half::Bool=false)
    result=expand(otype(Val(:Fock),typeof(term),typeof(bond),typeof(table)),term,bond,config,table,nothing)
    half || length(result)==0 || for opt in result|>values|>collect add!(result,opt') end
    return result
end

end # module
