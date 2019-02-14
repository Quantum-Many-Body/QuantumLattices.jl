using Printf: @printf,@sprintf
using ..Spatials: PID,AbstractBond,Point,Bond,pidtype
using ..DegreesOfFreedom: IID,Index,Internal,FilteredAttributes,Subscript,Subscripts,Coupling,Couplings
using ...Prerequisites: delta,decimaltostr
using ...Mathematics.AlgebraOverFields: SimpleID,ID
using ...Mathematics.VectorSpaces: AbstractVectorSpace,IsMultiIndexable,MultiIndexOrderStyle

import ..DegreesOfFreedom: wildcard,constant,defaultcenter,propercenters
import ...Prerequisites.Interfaces: dims,inds,rank,⊗,expand

export dims,inds,⊗,expand
export ANNIHILATION,CREATION,FID,FIndex,Fock
export usualfockindextotuple,nambufockindextotuple
export FCID,FockCoupling
export σ⁰,σˣ,σʸ,σᶻ,σ⁺,σ⁻
export FockCouplings

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
    union(::Type{P},::Type{FID}) where P<:PID

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
    expand(fc::FockCoupling,pid::PID,fock::Fock) -> Union{FCExpand,Tuple{}}
    expand(fc::FockCoupling,pids::NTuple{R,PID},focks::NTuple{R,Fock}) where R -> Union{FCExpand,Tuple{}}

Expand a Fock coupling with the given set of point ids and Fock degrees of freedom.
"""
expand(fc::FockCoupling,pid::PID,fock::Fock)=expand(fc,(pid,),(fock,))
function expand(fc::FockCoupling,pids::NTuple{R,PID},focks::NTuple{R,Fock}) where R
    centers=propercenters(typeof(fc),fc.id.centers,Val(R))
    rpids=NTuple{rank(fc),eltype(pids)}(pids[centers[i]] for i=1:rank(fc))
    rfocks=NTuple{rank(fc),eltype(focks)}(focks[centers[i]] for i=1:rank(fc))
    nambus=propernambus(fc.id.nambus,NTuple{rank(fc),Int}(rfocks[i].nnambu for i=1:rank(fc)))
    for (i,atom) in enumerate(fc.id.atoms)
       isa(atom,Int) && (atom≠rfocks[i].atom) && return ()
    end
    obsbexpands=expand(fc.obsubscripts,NTuple{rank(fc),Int}(rfocks[i].norbital for i=1:rank(fc)))|>collect
    spsbexpands=expand(fc.spsubscripts,NTuple{rank(fc),Int}(rfocks[i].nspin for i=1:rank(fc)))|>collect
    return FCExpand(fc.value,rpids,obsbexpands,spsbexpands,nambus)
end
defaultcenter(::Type{<:FockCoupling},i::Int,n::Int,::Val{2})=i<=n/2 ? 1 : 2
@generated function propernambus(nambus::NTuple{N,Any},ranges::NTuple{N,Int}) where N
    exprs=[:(isa(nambus[$i],Int) ?
            (0<nambus[$i]<=ranges[$i] ? nambus[$i] : error(" error: nambu out of range.")) :
            (($i)%2==1 ? min(CREATION,ranges[$i]) : ANNIHILATION)
            ) for i=1:N]
    return Expr(:tuple,exprs...)
end
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
    FockCouplings(  ::Val{N},value::Number=1;
                    centers::Union{NTuple{N,Int},Nothing}=nothing,
                    atoms::Union{NTuple{N,Int},Nothing}=nothing,
                    orbitals::Union{NTuple{N,Int},Subscript,Nothing}=nothing,
                    spins::Union{NTuple{N,Int},Subscript,Nothing}=nothing,
                    nambus::Union{NTuple{N,Int},Nothing}=nothing) where N -> Couplings

Construct an instance of `Couplings` that contains only one element of `FockCoupling`.
"""
function FockCouplings( ::Val{N},value::Number=1;
                        centers::Union{NTuple{N,Int},Nothing}=nothing,
                        atoms::Union{NTuple{N,Int},Nothing}=nothing,
                        orbitals::Union{NTuple{N,Int},Subscript,Nothing}=nothing,
                        spins::Union{NTuple{N,Int},Subscript,Nothing}=nothing,
                        nambus::Union{NTuple{N,Int},Nothing}=nothing) where N
    return Couplings(FockCoupling{N}(value,centers=centers,atoms=atoms,orbitals=orbitals,spins=spins,nambus=nambus))
end
