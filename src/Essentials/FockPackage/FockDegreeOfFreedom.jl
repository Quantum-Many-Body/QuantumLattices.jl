using Printf: @printf,@sprintf
using ..Spatial: PID,AbstractBond,Point,Bond,pidtype
using ..DegreeOfFreedom: IID,Index,Internal,FilteredAttributes,Coupling,Couplings
using ...Utilities: delta,indtosub,corder,decimaltostr
using ...Utilities.AlgebraOverField: SimpleID,ID

import ...Utilities: expand

export ANNIHILATION,CREATION,FID,FIndex,Fock
export usualfockindextotuple,nambufockindextotuple
export FCID,FockCoupling,expand
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

"""
    FID(;orbital::Int=1,spin::Int=1,nambu::Int=ANNIHILATION)

Create a Fock id.
"""
FID(;orbital::Int=1,spin::Int=1,nambu::Int=ANNIHILATION)=FID(orbital,spin,nambu)

"""
    fieldnames(::Type{FID}) -> NTuple{3,Symbol}

Get the fieldnames of `FID`.
"""
Base.fieldnames(::Type{FID})=(:orbital,:spin,:nambu)

"""
    adjoint(fid::FID) -> FID

Get the adjoint of a Fock id.
"""
Base.adjoint(fid::FID)=FID(fid.orbital,fid.spin,3-fid.nambu)

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

"""
    fieldnames(::Type{<:FIndex}) -> NTuple{5,Symbol}

Get the fieldnames of a Fock index.
"""
Base.fieldnames(::Type{<:FIndex})=(:scope,:site,:orbital,:spin,:nambu)

"""
    union(::Type{P},::Type{FID})

Get the union type of `PID` and `FID`.
"""
Base.union(::Type{P},::Type{FID}) where P<:PID=FIndex{fieldtype(P,:scope)}

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

"""
    Fock(;atom::Int=1,norbital::Int=1,nspin::Int=2,nnambu::Int=2)

Construct a Fock degrees of freedom.
"""
Fock(;atom::Int=1,norbital::Int=1,nspin::Int=2,nnambu::Int=2)=Fock(atom,norbital,nspin,nnambu)

"""
    length(fock::Fock) -> Int

Get the dimension of a Fock degrees of freedom.
"""
Base.length(fock::Fock)=prod((fock.norbital,fock.nspin,fock.nnambu))

"""
    iterate(fock::Fock,state=1)

Iterate over a Fock degrees of freedom.
"""
Base.iterate(fock::Fock,state=1)=state>length(fock) ? nothing : (FID(indtosub((fock.norbital,fock.nspin,fock.nnambu),state,corder)...),state+1)

"""
    usualfockindextotuple

Indicate that the filtered attributes are (:scope,:site,:orbital,:spin) when converting a Fock index to tuple.
"""
const usualfockindextotuple=FilteredAttributes(:scope,:site,:orbital,:spin)
"""
    nambufockindextotuple

Indicate that the filtered attributes are (:scope,:nambu,:site,:orbital,:spin) when converting a Fock index to tuple.
"""
const nambufockindextotuple=FilteredAttributes(:scope,:nambu,:site,:orbital,:spin)

"""
    FCID{A,O,S,N} <: SimpleID

The id of a Fock coupling.
"""
struct FCID{A,O,S,N} <: SimpleID
    atom::A
    orbital::O
    spin::S
    nambu::N
end

"""
    FCID(;atom=nothing,orbital=nothing,spin=nothing,nambu=nothing)

Construct an id of a Fock coupling.
"""
FCID(;atom=nothing,orbital=nothing,spin=nothing,nambu=nothing)=FCID(atom,orbital,spin,nambu)

"""
    fieldnames(::Type{<:FCID}) -> NTuple{4,Symbol}

Get the fieldnames of `FCID`.
"""
Base.fieldnames(::Type{<:FCID})=(:atom,:orbital,:spin,:nambu)

"""
    FockCoupling{N,V,I<:ID{N,<:FCID}} <: Coupling{V,I,N}

A Fock coupling.
"""
struct FockCoupling{N,V,I<:ID{N,<:FCID}} <: Coupling{V,I,N}
    value::V
    id::I
    FockCoupling(value::Number,id::ID{N,I}) where {N,I<:FCID}=new{N,value|>typeof,id|>typeof}(value,id)
end

"""
    FockCoupling{N}(    value::Number;
                        atoms::Union{NTuple{N,Int},Nothing}=nothing,
                        orbitals::Union{NTuple{N,Int},Nothing}=nothing,
                        spins::Union{NTuple{N,Int},Nothing}=nothing,
                        nambus::Union{NTuple{N,Int},Nothing}=nothing) where N

Construct a Fock coupling.
"""
function FockCoupling{N}(   value::Number;
                            atoms::Union{NTuple{N,Int},Nothing}=nothing,
                            orbitals::Union{NTuple{N,Int},Nothing}=nothing,
                            spins::Union{NTuple{N,Int},Nothing}=nothing,
                            nambus::Union{NTuple{N,Int},Nothing}=nothing) where N
    atoms==nothing && (atoms=ntuple(i->nothing,N))
    orbitals==nothing && (orbitals=ntuple(i->nothing,N))
    spins==nothing && (spins=ntuple(i->nothing,N))
    nambus==nothing && (nambus=ntuple(i->nothing,N))
    return FockCoupling(value,ID(FCID,atoms,orbitals,spins,nambus))
end

"""
    show(io::IO,fc::FockCoupling)

Show a Fock coupling.
"""
function Base.show(io::IO,fc::FockCoupling)
    @printf io "FockCoupling(value=%s" decimaltostr(fc.value)
    cache=[]
    for attrname in (:atoms,:orbitals,:spins,:nambus)
        attrvalue=getproperty(fc.id,attrname)
        all(attrvalue.≠nothing) && @printf io ",%s=(%s)" attrname join(attrvalue,",")
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
        attrvalue=getproperty(fc.id,attrname)
        all(attrvalue.≠nothing) && push!(cache,@sprintf "%s%s" abbr join(attrvalue))
    end
    length(cache)>0 ? (@sprintf "%s %s" decimaltostr(fc.value) join(cache,"⊗")) : decimaltostr(fc.value)
end

"""
    *(fc1::FockCoupling{2},fc2::FockCoupling{2}) -> FockCoupling{2}

Get the multiplication between two rank-2 Fock couplings.
"""
function Base.:*(fc1::FockCoupling{2},fc2::FockCoupling{2})
    value=fc1.value*fc2.value
    attrpairs=[]
    for attrname in (:atoms,:orbitals,:spins,:nambus)
        attrvalue1=getproperty(fc1.id,attrname)
        attrvalue2=getproperty(fc2.id,attrname)
        if attrvalue1≠(nothing,nothing) && attrvalue2≠(nothing,nothing)
            value=value*delta(attrvalue1[2],attrvalue2[1])
            push!(attrpairs,attrname=>(attrvalue1[1],attrvalue2[2]))
        elseif attrvalue1≠(nothing,nothing)
            push!(attrpairs,attrname=>attrvalue1)
        elseif attrvalue2≠(nothing,nothing)
            push!(attrpairs,attrname=>attrvalue2)
        end
    end
    FockCoupling{2}(value;attrpairs...)
end

"""
    expand(fc::FockCoupling,bond::AbstractBond,focks::Fock...) -> FCExpand

Expand a Fock coupling on a bond with the given Fock degrees of freedom.
"""
expand(fc::FockCoupling,bond::AbstractBond,focks::Fock...)=FCExpand(fc,bond,focks...)
struct FCExpand{M,N,V<:Number,P<:PID,O<:Union{NTuple{N,Int},Int},S<:Union{NTuple{N,Int},Int}}
    length::Int
    value::V
    pids::NTuple{N,P}
    obs::O
    sps::S
    nbs::NTuple{N,Int}
end
function FCExpand(fc::FockCoupling{2},point::Point,fock::Fock)
    if fc.id.atoms ∉ ((nothing,nothing),(fock.atom,fock.atom))
        FCExpand{'E',2,fc|>valtype,point|>pidtype,Int,Int}(0,fc.value,(point.pid,point.pid),0,0,(ANNIHILATION,ANNIHILATION))
    else
        wildcard=(nothing,nothing)
        norbital,nspin,nnambu=fock.norbital,fock.nspin,fock.nnambu
        orbitals,spins,nambus=fc.id.orbitals,fc.id.spins,fc.id.nambus
        obs=orbitals==wildcard ? norbital : (@assert 0<orbitals[1]<=norbital && 0<orbitals[2]<=norbital "FCExpand error: orbitals($orbitals) out of range.";orbitals)
        sps=spins==wildcard ? nspin : (@assert 0<spins[1]<=nspin && 0<spins[2]<=nspin "FCExpand error: spins($spins) out of range.";spins)
        nbs=nambus==wildcard ? (min(CREATION,nnambu),ANNIHILATION) : (@assert 0<nambus[1]<=nnambu && 0<nambus[2]<=nnambu "FCExpand error: nambus($nambus) out of range.";nambus)
        (mode,len)=orbitals==wildcard ? (spins==wildcard ? ('F',norbital*nspin) : ('S',norbital)) : (spins==wildcard ? ('O',nspin) : ('U',1))
        FCExpand{mode,2,fc|>valtype,point|>pidtype,obs|>typeof,sps|>typeof}(len,fc.value,(point.pid,point.pid),obs,sps,nbs)
    end
end
function FCExpand(fc::FockCoupling{2},bond::Bond,sf::Fock,ef::Fock)
    if fc.id.atoms ∉ ((nothing,nothing),(ef.atom,sf.atom))
        FCExpand{'E',2,fc|>valtype,bond|>pidtype,Int,Int}(0,fc.value,(bond.epoint.pid,bond.spoint.pid),0,0,(ANNIHILATION,ANNIHILATION))
    else
        wildcard=(nothing,nothing)
        orbitals,spins,nambus=fc.id.orbitals,fc.id.spins,fc.id.nambus
        obs=orbitals==wildcard ?
            (@assert ef.norbital==sf.norbital "FCExpand error: not equal norbital for efock($(ef.norbital)) and sfock($(sf.norbital)).";ef.norbital) :
            (@assert 0<orbitals[1]<=ef.norbital && 0<orbitals[2]<=sf.norbital "FCExpand error: orbitals($orbitals) out of range.";orbitals)
        sps=spins==wildcard ?
            (@assert sf.nspin==ef.nspin "FCExpand error: not equal nspin for efock($(ef.nspin)) and sfock($(sf.nspin)).";ef.nspin) :
            (@assert 0<spins[1]<=ef.nspin && 0<spins[2]<=sf.nspin "FCExpand error: spins($spins) out of range.";spins)
        nbs=nambus==wildcard ?
            (min(CREATION,ef.nnambu),ANNIHILATION) :
            (@assert 0<nambus[1]<=ef.nnambu && 0<nambus[2]<=sf.nnambu "FCExpand error: nambus($nambus) out of range.";nambus)
        (mode,len)=orbitals==wildcard ? (spins==wildcard ? ('F',ef.norbital*ef.nspin) : ('S',ef.norbital)) : (spins==wildcard ? ('O',ef.nspin) : ('U',1))
        FCExpand{mode,2,fc|>valtype,bond|>pidtype,obs|>typeof,sps|>typeof}(len,fc.value,(bond.epoint.pid,bond.spoint.pid),obs,sps,nbs)
    end
end
Base.eltype(::Type{<:FCExpand{M,N,V,P,O,S}}) where {M,N,V,P,O,S}=Tuple{V,NTuple{N,FIndex{fieldtype(P,:scope)}}}
Base.length(fce::FCExpand)=fce.length
@generated function FCExpandIndex(pids::NTuple{N,<:PID},obs::NTuple{N,Int},sps::NTuple{N,Int},nbs::NTuple{N,Int}) where N
    exprs=[:(FIndex(pids[$i],FID(obs[$i],sps[$i],nbs[$i]))) for i=1:N]
    return Expr(:tuple,exprs...)
end
Base.iterate(fce::FCExpand{'E'},t::Int=1)=nothing
Base.iterate(fce::FCExpand{'U'},t::Int=1)=t>length(fce) ? nothing : ((fce.value,FCExpandIndex(fce.pids,fce.obs,fce.sps,fce.nbs)),t+1)
Base.iterate(fce::FCExpand{'O'},t::Int=1)=t>length(fce) ? nothing : ((fce.value,FCExpandIndex(fce.pids,fce.obs,(t,t),fce.nbs)),t+1)
Base.iterate(fce::FCExpand{'S'},t::Int=1)=t>length(fce) ? nothing : ((fce.value,FCExpandIndex(fce.pids,(t,t),fce.sps,fce.nbs)),t+1)
Base.iterate(fce::FCExpand{'F'},t::Int=1)=t>length(fce) ? nothing : ((o,p)=indtosub((fce.obs,fce.sps),t,corder);((fce.value,FCExpandIndex(fce.pids,(o,o),(p,p),fce.nbs)),t+1))

"""
    σ⁰(mode::String)

The Pauli matrix σ⁰, which can act on the space of spins("sp"), orbitals("ob"), sublattices("sl") or particle-holes("ph").
"""
function σ⁰(mode::String)
    @assert (mode=="sp" || mode=="ob" || mode=="sl" || mode=="ph") "σ⁰ error: not supported mode($mode)."
    attrname=mode=="sp" ? :spins : mode=="ob" ? :orbitals : mode=="sl" ? :atoms : :nambus
    (attrval1,attrval2)=mode=="ph" ? ((ANNIHILATION,CREATION),(CREATION,ANNIHILATION)) : ((1,1),(2,2))
    FockCoupling{2}(1;attrname=>attrval1)+FockCoupling{2}(1;attrname=>attrval2)
end

"""
    σˣ(mode::String)

The Pauli matrix σˣ, which can act on the space of spins("sp"), orbitals("ob"), sublattices("sl") or particle-holes("ph").
"""
function σˣ(mode::String)
    @assert (mode=="sp" || mode=="ob" || mode=="sl" || mode=="ph") "σˣ error: not supported mode($mode)."
    attrname=mode=="sp" ? :spins : mode=="ob" ? :orbitals : mode=="sl" ? :atoms : :nambus
    (attrval1,attrval2)=mode=="ph" ? ((ANNIHILATION,ANNIHILATION),(CREATION,CREATION)) : ((1,2),(2,1))
    FockCoupling{2}(1;attrname=>attrval1)+FockCoupling{2}(1;attrname=>attrval2)
end

"""
    σʸ(mode::String)

The Pauli matrix σʸ, which can act on the space of spins("sp"), orbitals("ob"), sublattices("sl") or particle-holes("ph").
"""
function σʸ(mode::String)
    @assert (mode=="sp" || mode=="ob" || mode=="sl" || mode=="ph") "σʸ error: not supported mode($mode)."
    attrname=mode=="sp" ? :spins : mode=="ob" ? :orbitals : mode=="sl" ? :atoms : :nambus
    (attrval1,attrval2)=mode=="ph" ? ((ANNIHILATION,ANNIHILATION),(CREATION,CREATION)) : ((1,2),(2,1))
    FockCoupling{2}(1im;attrname=>attrval1)+FockCoupling{2}(-1im;attrname=>attrval2)
end

"""
    σᶻ(mode::String)

The Pauli matrix σᶻ, which can act on the space of spins("sp"), orbitals("ob"), sublattices("sl") or particle-holes("ph").
"""
function σᶻ(mode::String)
    @assert (mode=="sp" || mode=="ob" || mode=="sl" || mode=="ph") "σᶻ error: not supported mode($mode)."
    attrname=mode=="sp" ? :spins : mode=="ob" ? :orbitals : mode=="sl" ? :atoms : :nambus
    (attrval1,attrval2)=mode=="ph" ? ((ANNIHILATION,CREATION),(CREATION,ANNIHILATION)) : ((1,1),(2,2))
    FockCoupling{2}(-1;attrname=>attrval1)+FockCoupling{2}(1;attrname=>attrval2)
end

"""
    σ⁺(mode::String)

The Pauli matrix σ⁺, which can act on the space of spins("sp"), orbitals("ob"), sublattices("sl") or particle-holes("ph").
"""
function σ⁺(mode::String)
    @assert (mode=="sp" || mode=="ob" || mode=="sl" || mode=="ph") "σᶻ error: not supported mode($mode)."
    attrname=mode=="sp" ? :spins : mode=="ob" ? :orbitals : mode=="sl" ? :atoms : :nambus
    attrval=mode=="ph" ? (CREATION,CREATION) : (2,1)
    Couplings(FockCoupling{2}(1;attrname=>attrval))
end

"""
    σ⁻(mode::String)

The Pauli matrix σ⁻, which can act on the space of spins("sp"), orbitals("ob"), sublattices("sl") or particle-holes("ph").
"""
function σ⁻(mode::String)
    @assert (mode=="sp" || mode=="ob" || mode=="sl" || mode=="ph") "σᶻ error: not supported mode($mode)."
    attrname=mode=="sp" ? :spins : mode=="ob" ? :orbitals : mode=="sl" ? :atoms : :nambus
    attrval=mode=="ph" ? (ANNIHILATION,ANNIHILATION) : (1,2)
    Couplings(FockCoupling{2}(1;attrname=>attrval))
end
