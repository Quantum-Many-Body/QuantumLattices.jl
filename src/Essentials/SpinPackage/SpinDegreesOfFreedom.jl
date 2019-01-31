using Printf: @printf,@sprintf
using ..Spatials: PID
using ..DegreesOfFreedom: IID,Internal,Index,FilteredAttributes,Subscript,Subscripts,Coupling,Couplings
using ...Prerequisites: Float,decimaltostr
using ...Mathematics.VectorSpaces: AbstractVectorSpace,IsMultiIndexable,MultiIndexOrderStyle
using ...Mathematics.AlgebraOverFields: SimpleID,ID

import ..DegreesOfFreedom: wildcard,constant,defaultcenter,propercenters
import ...Prerequisites.Interfaces: dims,inds,expand,rank

export dims,inds,expand
export SID,Spin,SIndex,usualspinindextotuple
export SCID,SpinCoupling
export Heisenberg,Ising,Sˣ,Sʸ,Sᶻ
export SpinCouplings

const sidtagmap=Dict(1=>'i',2=>'x',3=>'y',4=>'z',5=>'+',6=>'-')
const sidseqmap=Dict(v=>k for (k,v) in sidtagmap)
const sidajointmap=Dict('i'=>'i','x'=>'x','y'=>'y','z'=>'z','+'=>'-','-'=>'+')

"""
    SID <: IID

The spin id.
"""
struct SID <: IID
    orbital::Int
    spin::Float
    tag::Char
    function SID(orbital::Int,spin::Real,tag::Char)
        @assert tag in ('i','x','y','z','+','-') "SID error: not supported tag($tag)."
        new(orbital,spin,tag)
    end
end
Base.fieldnames(::Type{<:SID})=(:orbital,:spin,:tag)

"""
    SID(;orbital::Int=1,spin::Real=0.5,tag::Char='i')

Create a spin id.
"""
SID(;orbital::Int=1,spin::Real=0.5,tag::Char='i')=SID(orbital,spin,tag)

"""
    adjoint(sid::SID) -> SID

Get the adjoint of a spin id.
"""
Base.adjoint(sid::SID)=SID(sid.orbital,sid.spin,sidajointmap[sid.tag])

"""
    Spin <: Internal{SID}

The spin interanl degrees of freedom.
"""
struct Spin <: Internal{SID}
    atom::Int
    norbital::Int
    spin::Float
end
IsMultiIndexable(::Type{Spin})=IsMultiIndexable(true)
MultiIndexOrderStyle(::Type{Spin})=MultiIndexOrderStyle('C')
dims(sp::Spin)=(sp.norbital,length(sidtagmap))
inds(sid::SID,::Spin)=(sid.orbital,sidseqmap[sid.tag])
SID(inds::NTuple{2,Int},sp::Spin)=SID(inds[1],sp.spin,sidtagmap[inds[2]])

"""
    Spin(;atom::Int=1,norbital::Int=1,spin::Real=0.5)

Construct a spin degrees of freedom.
"""
Spin(;atom::Int=1,norbital::Int=1,spin::Real=0.5)=Spin(atom,norbital,spin)

"""
    SIndex{S} <: Index{PID{S},SID}

The spin index.
"""
struct SIndex{S} <: Index{PID{S},SID}
    scope::S
    site::Int
    orbital::Int
    spin::Float
    tag::Char
end
Base.fieldnames(::Type{<:SIndex})=(:scope,:site,:orbital,:spin,:tag)

"""
    union(::Type{P},::Type{SID}) where {P<:PID}

Get the union type of `PID` and `SID`.
"""
Base.union(::Type{P},::Type{SID}) where {P<:PID}=SIndex{fieldtype(P,:scope)}

"""
    usualspinindextotuple

Indicate that the filtered attributes are `(:scope,:site,:orbital)` when converting a spin index to tuple.
"""
const usualspinindextotuple=FilteredAttributes(:scope,:site,:orbital)

"""
    SCID(;center=wildcard,atom=wildcard,orbital=wildcard,tag='i',subscript=wildcard)=SCID(center,atom,orbital,tag,subscript)

The id of a spin coupling.
"""
struct SCID{C,A,O,S} <: SimpleID
    center::C
    atom::A
    orbital::O
    tag::Char
    subscript::S
    function SCID(center,atom,orbital,tag::Char,subscript)
        @assert tag in ('i','x','y','z','+','-') "SCID error: not supported tag($tag)."
        new{typeof(center),typeof(atom),typeof(orbital),typeof(subscript)}(center,atom,orbital,tag,subscript)
    end
end
SCID(;center=wildcard,atom=wildcard,orbital=wildcard,tag='i',subscript=wildcard)=SCID(center,atom,orbital,tag,subscript)
Base.fieldnames(::Type{<:SCID})=(:center,:atom,:orbital,:tag,:subscript)

"""
    SpinCoupling(value::Number,id::ID{N,I},subscripts::Subscripts) where {N,I<:SCID}
    SpinCoupling{N}(    value::Number=1;
                        tags::NTuple{N,Char},
                        centers::Union{NTuple{N,Int},Nothing}=nothing,
                        atoms::Union{NTuple{N,Int},Nothing}=nothing,
                        orbitals::Union{NTuple{N,Int},Subscript,Nothing}=nothing
                        ) where N

Spin coupling.
"""
struct SpinCoupling{N,V<:Number,I<:ID{N,<:SCID},S<:Subscripts} <: Coupling{N,V,I}
    value::V
    id::I
    subscripts::S
    SpinCoupling(value::Number,id::ID{N,I},subscripts::Subscripts) where {N,I<:SCID}=new{N,value|>typeof,id|>typeof,typeof(subscripts)}(value,id,subscripts)
end
function SpinCoupling{N}(   value::Number=1;
                            tags::NTuple{N,Char},
                            centers::Union{NTuple{N,Int},Nothing}=nothing,
                            atoms::Union{NTuple{N,Int},Nothing}=nothing,
                            orbitals::Union{NTuple{N,Int},Subscript,Nothing}=nothing
                            ) where N
    centers===nothing && (centers=ntuple(i->wildcard,N))
    atoms===nothing && (atoms=ntuple(i->wildcard,N))
    subscript=isa(orbitals,NTuple{N,Int}) ? Subscript(orbitals) : isa(orbitals,Nothing) ? Subscript{N}() : orbitals
    orbitals=subscript.opattern
    subscripts=ntuple(i->subscript.identifier,N)
    return SpinCoupling(value,ID(SCID,centers,atoms,orbitals,tags,subscripts),Subscripts(subscript))
end

"""
    show(io::IO,sc::SpinCoupling)

Show a spin coupling.
"""
function Base.show(io::IO,sc::SpinCoupling)
    @printf io "SpinCoupling(value=%s" decimaltostr(sc.value)
    cache=[]
    for attrname in (:centers,:atoms,:orbitals,:tags)
        attrvalue=getproperty(sc.id,attrname)
        all(attrvalue.≠wildcard) && @printf io ",%s=(%s)" attrname join(attrvalue,",")
    end
    all((subscripts=sc.id.subscripts).==wildcard) || all((subscripts.==constant)) || @printf io ",subscripts=(%s)" join(subscripts,",")
    @printf io ")"
end

"""
    repr(sc::SpinCoupling) -> String

Get the repr representation of a spin coupling.
"""
function Base.repr(sc::SpinCoupling)
    cache=[]
    for (attrname,abbr) in zip((:atoms,:orbitals),("sl","ob"))
        attrvalue=getproperty(sc.id,attrname)
        all(attrvalue.≠wildcard) && push!(cache,@sprintf "%s(%s)" abbr join(attrvalue,':'))
    end
    result=@sprintf "%s S%s" decimaltostr(sc.value) join(sc.id.tags,'S')
    length(cache)>0 && (result=@sprintf "%s %s" result join(cache,"⊗"))
    any((centers=sc.id.centers).≠wildcard) && (result=@sprintf "%s@(%s)" result join(centers,'-'))
    subscripts=sc.id.subscripts
    all(subscripts.==wildcard) || all(subscripts.==constant) || (result=@sprintf "%s with %s" result join(subscripts," && "))
    return result
end

"""
    *(sc1::SpinCoupling,sc2::SpinCoupling) -> SpinCoupling

Get the multiplication between two spin couplings.
"""
Base.:*(sc1::SpinCoupling,sc2::SpinCoupling)=SpinCoupling(sc1.value*sc2.value,sc1.id*sc2.id,sc1.subscripts*sc2.subscripts)

"""
    expand(sc::SpinCoupling,pid::PID,spin::Spin) -> SCExpand
    expand(sc::SpinCoupling,pids::NTuple{N,PID},spins::NTuple{N,Spin}) where N -> SCExpand

Expand a spin coupling with the given set of point ids and spin degrees of freedom.
"""
expand(sc::SpinCoupling,pid::PID,spin::Spin)=expand(sc,(pid,),(spin,))
function expand(sc::SpinCoupling,pids::NTuple{N,PID},spins::NTuple{N,Spin}) where N
    centers=propercenters(typeof(sc),sc.id.centers,Val(N))
    rpids=NTuple{rank(sc),eltype(pids)}(pids[centers[i]] for i=1:rank(sc))
    rspins=NTuple{rank(sc),eltype(spins)}(spins[centers[i]] for i=1:rank(sc))
    sps=NTuple{rank(sc),Float}(rspins[i].spin for i=1:rank(sc))
    for (i,atom) in enumerate(sc.id.atoms)
        isa(atom,Int) && (atom≠rspins[i].atom) && return SCExpand(sc.value,rpids,Vector{NTuple{rank(sc),Int}}(),sps,sc.id.tags)
    end
    sbexpands=expand(sc.subscripts,NTuple{rank(sc),Int}(rspins[i].norbital for i=1:rank(sc)))|>collect
    return SCExpand(sc.value,rpids,sbexpands,sps,sc.id.tags)
end
defaultcenter(::Type{<:SpinCoupling},i::Int,n::Int,::Val{2})=i%2==1 ? 1 : 2
struct SCExpand{V<:Number,N,S} <: AbstractVectorSpace{Tuple{V,NTuple{N,SIndex{S}}}}
    value::V
    pids::NTuple{N,PID{S}}
    sbexpands::Vector{NTuple{N,Int}}
    spins::NTuple{N,Float}
    tags::NTuple{N,Char}
end
IsMultiIndexable(::Type{<:SCExpand})=IsMultiIndexable(true)
MultiIndexOrderStyle(::Type{<:SCExpand})=MultiIndexOrderStyle('C')
dims(sce::SCExpand)=(length(sce.sbexpands),)
@generated function Tuple(index::Tuple{Int},sce::SCExpand{<:Number,N}) where N
    exprs=[:(SIndex(sce.pids[$i],SID(sce.sbexpands[index[1]][$i],sce.spins[$i],sce.tags[$i]))) for i=1:N]
    return Expr(:tuple,:(sce.value),Expr(:tuple,exprs...))
end

"""
    Heisenberg(;centers::Union{NTuple{2,Int},Nothing}=nothing,
                atoms::Union{NTuple{2,Int},Nothing}=nothing,
                orbitals::Union{NTuple{2,Int},Subscript,Nothing}=nothing
                ) -> Couplings{ID{2,SCID},SpinCoupling{2,Float,ID{2,SCID}}}

The Heisenberg couplings.
"""
function Heisenberg(;centers::Union{NTuple{2,Int},Nothing}=nothing,atoms::Union{NTuple{2,Int},Nothing}=nothing,orbitals::Union{NTuple{2,Int},Subscript,Nothing}=nothing)
    sc1=SpinCoupling{2}(0.5,centers=centers,atoms=atoms,orbitals=orbitals,tags=('+','-'))
    sc2=SpinCoupling{2}(0.5,centers=centers,atoms=atoms,orbitals=orbitals,tags=('-','+'))
    sc3=SpinCoupling{2}(1.0,centers=centers,atoms=atoms,orbitals=orbitals,tags=('z','z'))
    return Couplings(sc1,sc2,sc3)
end

"""
    Ising(  tag::Char;
            centers::Union{NTuple{2,Int},Nothing}=nothing,
            atoms::Union{NTuple{2,Int},Nothing}=nothing,
            orbitals::Union{NTuple{2,Int},Subscript,Nothing}=nothing
            ) -> Couplings{ID{2,SCID},SpinCoupling{2,Float,ID{2,SCID}}}

The Ising couplings.
"""
function Ising(tag::Char;centers::Union{NTuple{2,Int},Nothing}=nothing,atoms::Union{NTuple{2,Int},Nothing}=nothing,orbitals::Union{NTuple{2,Int},Subscript,Nothing}=nothing)
    @assert tag in ('x','y','z') "Ising error: not supported input tag($tag)."
    return Couplings(SpinCoupling{2}(1.0,centers=centers,atoms=atoms,orbitals=orbitals,tags=(tag,tag)))
end

"""
    Sˣ(;center::Union{Int,Nothing}=nothing,atom::Union{Int,Nothing}=nothing,orbital::Union{Int,Nothing}=nothing) -> Couplings{ID{1,SCID},SpinCoupling{1,Float,ID{1,SCID}}}

The single Sˣ coupling.
"""
function Sˣ(;center::Union{Int,Nothing}=nothing,atom::Union{Int,Nothing}=nothing,orbital::Union{Int,Nothing}=nothing)
    Couplings(SpinCoupling{1}(  1.0,tags=('x',),
                                centers=center===nothing ? nothing : (center,),
                                atoms=atom===nothing ? nothing : (atom,),
                                orbitals=orbital==nothing ? nothing : (orbital,)
                                ))
end

"""
    Sʸ(;center::Union{Int,Nothing}=nothing,atom::Union{Int,Nothing}=nothing,orbital::Union{Int,Nothing}=nothing) -> Couplings{ID{1,SCID},SpinCoupling{1,Float,ID{1,SCID}}}

The single Sʸ coupling.
"""
function Sʸ(;center::Union{Int,Nothing}=nothing,atom::Union{Int,Nothing}=nothing,orbital::Union{Int,Nothing}=nothing)
    Couplings(SpinCoupling{1}(  1.0,tags=('y',),
                                centers=center===nothing ? nothing : (center,),
                                atoms=atom===nothing ? nothing : (atom,),
                                orbitals=orbital==nothing ? nothing : (orbital,)
                                ))
end

"""
    Sᶻ(;center::Union{Int,Nothing}=nothing,atom::Union{Int,Nothing}=nothing,orbital::Union{Int,Nothing}=nothing) -> Couplings{ID{1,SCID},SpinCoupling{1,Float,ID{1,SCID}}}

The single Sᶻ coupling.
"""
function Sᶻ(;center::Union{Int,Nothing}=nothing,atom::Union{Int,Nothing}=nothing,orbital::Union{Int,Nothing}=nothing)
    Couplings(SpinCoupling{1}(  1.0,tags=('z',),
                                centers=center===nothing ? nothing : (center,),
                                atoms=atom===nothing ? nothing : (atom,),
                                orbitals=orbital==nothing ? nothing : (orbital,)
                                ))
end

"""
    SpinCouplings(  ::Val{N},value::Number=1;
                    tags::NTuple{N,Char},
                    centers::Union{NTuple{N,Int},Nothing}=nothing,
                    atoms::Union{NTuple{N,Int},Nothing}=nothing,
                    orbitals::Union{NTuple{N,Int},Subscript,Nothing}=nothing
                    ) where N -> Couplings

Construct an instance of `Couplings` that contains only one element of `SpinCoupling`.
"""
function SpinCouplings( ::Val{N},value::Number=1;
                        tags::NTuple{N,Char},
                        centers::Union{NTuple{N,Int},Nothing}=nothing,
                        atoms::Union{NTuple{N,Int},Nothing}=nothing,
                        orbitals::Union{NTuple{N,Int},Subscript,Nothing}=nothing
                        ) where N
    return Couplings(SpinCoupling{N}(value,tags=tags,centers=centers,atoms=atoms,orbitals=orbitals))
end
