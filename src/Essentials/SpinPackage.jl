module SpinPackage

using Printf: @printf,@sprintf
using ..Spatials: PID
using ..DegreesOfFreedom: IID,Internal,Index,FilteredAttributes,OID,Operator,LaTeX,latexformat
using ..Terms: wildcard,constant,Subscript,Subscripts,subscriptexpr,Coupling,Couplings,couplingcenters,Term,TermCouplings,TermAmplitude,TermModulate
using ...Interfaces: rank,kind
using ...Prerequisites: Float,decimaltostr,delta
using ...Mathematics.VectorSpaces: VectorSpace,IsMultiIndexable,MultiIndexOrderStyle
using ...Mathematics.AlgebraOverFields: SimpleID,ID

import ..DegreesOfFreedom: script,otype,isHermitian
import ..Terms: couplingcenter,statistics,abbr
import ...Interfaces: dims,inds,expand,matrix,permute

export SID,Spin,SIndex
export usualspinindextotuple
export SOperator,soptdefaultlatex
export SCID,SpinCoupling
export Heisenberg,Ising,Gamma,DM,Sˣ,Sʸ,Sᶻ
export @sc_str,@heisenberg_str,@ising_str,@gamma_str,@dm_str,@sˣ_str,@sʸ_str,@sᶻ_str
export SpinTerm

const sidtagmap=Dict(1=>'x',2=>'y',3=>'z',4=>'+',5=>'-')
const sidseqmap=Dict(v=>k for (k,v) in sidtagmap)
const sidajointmap=Dict('x'=>'x','y'=>'y','z'=>'z','+'=>'-','-'=>'+')

"""
    SID <: IID

The spin id.
"""
struct SID <: IID
    orbital::Int
    spin::Float
    tag::Char
    function SID(orbital::Int,spin::Real,tag::Char)
        @assert tag in ('x','y','z','+','-') "SID error: not supported tag($tag)."
        new(orbital,spin,tag)
    end
end
Base.fieldnames(::Type{<:SID})=(:orbital,:spin,:tag)

"""
    SID(;orbital::Int=1,spin::Real=0.5,tag::Char='z')

Create a spin id.
"""
SID(;orbital::Int=1,spin::Real=0.5,tag::Char='z')=SID(orbital,spin,tag)

"""
    adjoint(sid::SID) -> SID

Get the adjoint of a spin id.
"""
Base.adjoint(sid::SID)=SID(sid.orbital,sid.spin,sidajointmap[sid.tag])

"""
    matrix(sid::SID,dtype::Type{<:Number}=Complex{Float}) -> Matrix{dtype}

Get the matrix representation of a sid.
"""
function matrix(sid::SID,dtype::Type{<:Number}=Complex{Float})
    N=Int(2*sid.spin+1)
    result=zeros(dtype,(N,N))
    for i=1:N,j=1:N
        row,col=N+1-i,N+1-j
        m,n=sid.spin+1-i,sid.spin+1-j
        result[row,col]=sid.tag=='x' ? (delta(i+1,j)+delta(i,j+1))*sqrt(sid.spin*(sid.spin+1)-m*n)/2 :
                        sid.tag=='y' ? (delta(i+1,j)-delta(i,j+1))*sqrt(sid.spin*(sid.spin+1)-m*n)/2im :
                        sid.tag=='z' ? delta(i,j)*m :
                        sid.tag=='+' ? delta(i+1,j)*sqrt(sid.spin*(sid.spin+1)-m*n) :
                                       delta(i,j+1)*sqrt(sid.spin*(sid.spin+1)-m*n)
    end
    return result
end

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
    script(oid::OID{<:SIndex},::Val{:site}) -> Int
    script(oid::OID{<:SIndex},::Val{:orbital}) -> Int
    script(oid::OID{<:SIndex},::Val{:spin}) -> Float
    script(oid::OID{<:SIndex},::Val{:tag}) -> Char

Get the required script of a spin oid.
"""
script(oid::OID{<:SIndex},::Val{:site})=oid.index.site
script(oid::OID{<:SIndex},::Val{:orbital})=oid.index.orbital
script(oid::OID{<:SIndex},::Val{:spin})=oid.index.spin
script(oid::OID{<:SIndex},::Val{:tag})=oid.index.tag

"""
    usualspinindextotuple

Indicate that the filtered attributes are `(:scope,:site,:orbital)` when converting a spin index to tuple.
"""
const usualspinindextotuple=FilteredAttributes(:scope,:site,:orbital)

"""
    SOperator(value::Number,id::ID{OID}=ID())

Spin operator.
"""
struct SOperator{V<:Number,I<:ID{OID}} <: Operator{V,I}
    value::V
    id::I
    SOperator(value::Number)=new{typeof(vlalue),Tuple{}}(value,())
    SOperator(value::Number,id::ID{OID})=new{typeof(value),typeof(id)}(value,id)
end

"""
    soptdefaultlatex

The default LaTeX format of the oids of a spin operator.
"""
const soptdefaultlatex=LaTeX{(:tag,),(:site,:orbital)}('S')
latexformat(SOperator,soptdefaultlatex)

"""
    statistics(opt::SOperator) -> Char
    statistics(::Type{<:SOperator}) -> Char

Get the statistics of SOperator.
"""
statistics(opt::SOperator)=opt|>typeof|>statistics
statistics(::Type{<:SOperator})='B'

"""
    permute(::Type{<:SOperator},id1::OID{<:SIndex},id2::OID{<:SIndex},table) -> Tuple{Vararg{SOperator}}

Permute two fermionic oid and get the result.
"""
function permute(::Type{<:SOperator},id1::OID{<:SIndex},id2::OID{<:SIndex},table)
    @assert id1.index≠id2.index || id1.rcoord≠id2.rcoord || id1.icoord≠id2.icoord "permute error: permuted ids should not be equal to each other."
    if usualspinindextotuple(id1.index)==usualspinindextotuple(id2.index) && id1.rcoord==id2.rcoord && id1.icoord==id2.icoord
        @assert id1.index.spin==id2.index.spin "permute error: noncommutable ids should have the same spin attribute."
        if id1.index.tag=='x'
            id2.index.tag=='y' && return (SOperator(+1im,ID(permutesoid(id1,'z',table))),SOperator(1,ID(id2,id1)))
            id2.index.tag=='z' && return (SOperator(-1im,ID(permutesoid(id1,'y',table))),SOperator(1,ID(id2,id1)))
            id2.index.tag=='+' && return (SOperator(-1,ID(permutesoid(id1,'z',table))),SOperator(1,ID(id2,id1)))
            id2.index.tag=='-' && return (SOperator(+1,ID(permutesoid(id1,'z',table))),SOperator(1,ID(id2,id1)))
        elseif id1.index.tag=='y'
            id2.index.tag=='x' && return (SOperator(-1im,ID(permutesoid(id1,'z',table))),SOperator(1,ID(id2,id1)))
            id2.index.tag=='z' && return (SOperator(+1im,ID(permutesoid(id1,'x',table))),SOperator(1,ID(id2,id1)))
            id2.index.tag=='+' && return (SOperator(-1im,ID(permutesoid(id1,'z',table))),SOperator(1,ID(id2,id1)))
            id2.index.tag=='-' && return (SOperator(-1im,ID(permutesoid(id1,'z',table))),SOperator(1,ID(id2,id1)))
        elseif id1.index.tag=='z'
            id2.index.tag=='x' && return (SOperator(+1im,ID(permutesoid(id1,'y',table))),SOperator(1,ID(id2,id1)))
            id2.index.tag=='y' && return (SOperator(-1im,ID(permutesoid(id1,'x',table))),SOperator(1,ID(id2,id1)))
            id2.index.tag=='+' && return (SOperator(+1,ID(id2)),SOperator(1,ID(id2,id1)))
            id2.index.tag=='-' && return (SOperator(-1,ID(id2)),SOperator(1,ID(id2,id1)))
        elseif id1.index.tag=='+'
            id2.index.tag=='x' && return (SOperator(+1,ID(permutesoid(id1,'z',table))),SOperator(1,ID(id2,id1)))
            id2.index.tag=='y' && return (SOperator(+1im,ID(permutesoid(id1,'z',table))),SOperator(1,ID(id2,id1)))
            id2.index.tag=='z' && return (SOperator(-1,ID(id1)),SOperator(1,ID(id2,id1)))
            id2.index.tag=='-' && return (SOperator(+2,ID(permutesoid(id1,'z',table))),SOperator(1,ID(id2,id1)))
        elseif id1.index.tag=='-'
            id2.index.tag=='x' && return (SOperator(-1,ID(permutesoid(id1,'z',table))),SOperator(1,ID(id2,id1)))
            id2.index.tag=='y' && return (SOperator(1im,ID(permutesoid(id1,'z',table))),SOperator(1,ID(id2,id1)))
            id2.index.tag=='z' && return (SOperator(+1,ID(id1)),SOperator(1,ID(id2,id1)))
            id2.index.tag=='+' && return (SOperator(-2,ID(permutesoid(id1,'z',table))),SOperator(1,ID(id2,id1)))
        end
    else
        return (SOperator(1,ID(id2,id1)),)
    end
end
function permutesoid(id::OID{<:SIndex},tag::Char,table)
    index=replace(id.index,tag=tag)
    fieldtype(id|>typeof,:seq)<:Nothing && return replace(id,index=index)
    seq=get(table,index,nothing)
    seq===nothing && error("permutesoid error: incomplete index-sequence table with $index not found.")
    return replace(id,index=index,seq=seq)
end

"""
    SCID(;center=wildcard,atom=wildcard,orbital=wildcard,tag='z',subscript=wildcard)=SCID(center,atom,orbital,tag,subscript)

The id of a spin coupling.
"""
struct SCID{C,A,O,S} <: SimpleID
    center::C
    atom::A
    orbital::O
    tag::Char
    subscript::S
    function SCID(center,atom,orbital,tag::Char,subscript)
        @assert tag in ('x','y','z','+','-') "SCID error: not supported tag($tag)."
        new{typeof(center),typeof(atom),typeof(orbital),typeof(subscript)}(center,atom,orbital,tag,subscript)
    end
end
SCID(;center=wildcard,atom=wildcard,orbital=wildcard,tag='z',subscript=wildcard)=SCID(center,atom,orbital,tag,subscript)
Base.fieldnames(::Type{<:SCID})=(:center,:atom,:orbital,:tag,:subscript)

"""
    SpinCoupling(value::Number,id::ID{SCID},subscripts::Subscripts)
    SpinCoupling{N}(    value::Number=1;
                        tags::NTuple{N,Char},
                        centers::Union{NTuple{N,Int},Nothing}=nothing,
                        atoms::Union{NTuple{N,Int},Nothing}=nothing,
                        orbitals::Union{NTuple{N,Int},Subscript,Nothing}=nothing
                        ) where N

Spin coupling.
"""
struct SpinCoupling{V,I<:ID{SCID},S<:Subscripts} <: Coupling{V,I}
    value::V
    id::I
    subscripts::S
    function SpinCoupling(value::Number,id::ID{SCID},subscripts::Subscripts)
        new{typeof(value),typeof(id),typeof(subscripts)}(value,id,subscripts)
    end
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
        all(attrvalue.≠wildcard) && push!(cache,@sprintf "%s(%s)" abbr join(attrvalue,','))
    end
    result=@sprintf "%s S%s" decimaltostr(sc.value) join(sc.id.tags,'S')
    length(cache)>0 && (result=@sprintf "%s %s" result join(cache,"⊗"))
    any((centers=sc.id.centers).≠wildcard) && (result=@sprintf "%s%s@(%s)" result (length(cache)>0 ? "" : " ") join(centers,','))
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
    expand(sc::SpinCoupling,pid::PID,spin::Spin,kind::Union{Val{K},Nothing}=nothing) where K -> Union{SCExpand,Tuple{}}
    expand(sc::SpinCoupling,pids::NTuple{N,PID},spins::NTuple{N,Spin},kind::Union{Val{K},Nothing}=nothing) where {N,K} -> Union{SCExpand,Tuple{}}

Expand a spin coupling with the given set of point ids and spin degrees of freedom.
"""
expand(sc::SpinCoupling,pid::PID,spin::Spin,kind::Union{Val{K},Nothing}=nothing) where K=expand(sc,(pid,),(spin,),kind)
function expand(sc::SpinCoupling,pids::NTuple{N,PID},spins::NTuple{N,Spin},kind::Union{Val{K},Nothing}=nothing) where {N,K}
    centers=couplingcenters(typeof(sc),sc.id.centers,Val(N))
    rpids=NTuple{rank(sc),eltype(pids)}(pids[centers[i]] for i=1:rank(sc))
    rspins=NTuple{rank(sc),eltype(spins)}(spins[centers[i]] for i=1:rank(sc))
    sps=NTuple{rank(sc),Float}(rspins[i].spin for i=1:rank(sc))
    for (i,atom) in enumerate(sc.id.atoms)
        isa(atom,Int) && (atom≠rspins[i].atom) && return ()
    end
    sbexpands=expand(sc.subscripts,NTuple{rank(sc),Int}(rspins[i].norbital for i=1:rank(sc)))|>collect
    return SCExpand(sc.value,rpids,sbexpands,sps,sc.id.tags)
end
couplingcenter(::Type{<:SpinCoupling},i::Int,n::Int,::Val{2})=i%2==1 ? 1 : 2
struct SCExpand{V,N,S} <: VectorSpace{Tuple{V,NTuple{N,SIndex{S}}}}
    value::V
    pids::NTuple{N,PID{S}}
    sbexpands::Vector{NTuple{N,Int}}
    spins::NTuple{N,Float}
    tags::NTuple{N,Char}
end
IsMultiIndexable(::Type{<:SCExpand})=IsMultiIndexable(true)
MultiIndexOrderStyle(::Type{<:SCExpand})=MultiIndexOrderStyle('C')
dims(sce::SCExpand)=(length(sce.sbexpands),)
@generated function Tuple(index::Tuple{Int},sce::SCExpand{V,N}) where {V,N}
    exprs=[:(SIndex(sce.pids[$i],SID(sce.sbexpands[index[1]][$i],sce.spins[$i],sce.tags[$i]))) for i=1:N]
    return Expr(:tuple,:(sce.value),Expr(:tuple,exprs...))
end

"""
    Heisenberg( mode::String="+-z";
                centers::Union{NTuple{2,Int},Nothing}=nothing,
                atoms::Union{NTuple{2,Int},Nothing}=nothing,
                orbitals::Union{NTuple{2,Int},Subscript,Nothing}=nothing
                ) -> Couplings{ID{SCID,2},SpinCoupling{Rational{Int}/Int,<:ID{SCID,2}}}

The Heisenberg couplings.
"""
function Heisenberg(mode::String="+-z";
                    centers::Union{NTuple{2,Int},Nothing}=nothing,
                    atoms::Union{NTuple{2,Int},Nothing}=nothing,
                    orbitals::Union{NTuple{2,Int},Subscript,Nothing}=nothing
                    )
    @assert mode=="+-z" || mode=="xyz" "Heisenberg error: not supported mode($mode)."
    if mode=="+-z"
        sc1=SpinCoupling{2}(1//2,centers=centers,atoms=atoms,orbitals=orbitals,tags=('+','-'))
        sc2=SpinCoupling{2}(1//2,centers=centers,atoms=atoms,orbitals=orbitals,tags=('-','+'))
        sc3=SpinCoupling{2}(1//1,centers=centers,atoms=atoms,orbitals=orbitals,tags=('z','z'))
    else
        sc1=SpinCoupling{2}(1,centers=centers,atoms=atoms,orbitals=orbitals,tags=('x','x'))
        sc2=SpinCoupling{2}(1,centers=centers,atoms=atoms,orbitals=orbitals,tags=('y','y'))
        sc3=SpinCoupling{2}(1,centers=centers,atoms=atoms,orbitals=orbitals,tags=('z','z'))
    end
    return Couplings(sc1,sc2,sc3)
end

"""
    Ising(  tag::Char;
            centers::Union{NTuple{2,Int},Nothing}=nothing,
            atoms::Union{NTuple{2,Int},Nothing}=nothing,
            orbitals::Union{NTuple{2,Int},Subscript,Nothing}=nothing
            ) -> Couplings{ID{SCID,2},SpinCoupling{Rational{Int}/Int,<:ID{SCID,2}}}

The Ising couplings.
"""
function Ising(tag::Char;centers::Union{NTuple{2,Int},Nothing}=nothing,atoms::Union{NTuple{2,Int},Nothing}=nothing,orbitals::Union{NTuple{2,Int},Subscript,Nothing}=nothing)
    @assert tag in ('x','y','z') "Ising error: not supported input tag($tag)."
    return Couplings(SpinCoupling{2}(1,centers=centers,atoms=atoms,orbitals=orbitals,tags=(tag,tag)))
end

"""
    Gamma(  tag::Char;
            centers::Union{NTuple{2,Int},Nothing}=nothing,
            atoms::Union{NTuple{2,Int},Nothing}=nothing,
            orbitals::Union{NTuple{2,Int},Subscript,Nothing}=nothing
            ) -> Couplings{ID{SCID,2},SpinCoupling{Rational{Int}/Int,<:ID{SCID,2}}}

The Gamma couplings.
"""
function Gamma(tag::Char;centers::Union{NTuple{2,Int},Nothing}=nothing,atoms::Union{NTuple{2,Int},Nothing}=nothing,orbitals::Union{NTuple{2,Int},Subscript,Nothing}=nothing)
    @assert tag in ('x','y','z') "Gamma error: not supported input tag($tag)."
    t1,t2=tag=='x' ? ('y','z') : tag=='y' ? ('z','x') : ('x','y')
    sc1=SpinCoupling{2}(1,centers=centers,atoms=atoms,orbitals=orbitals,tags=(t1,t2))
    sc2=SpinCoupling{2}(1,centers=centers,atoms=atoms,orbitals=orbitals,tags=(t2,t1))
    return Couplings(sc1,sc2)
end

"""
    DM( tag::Char;
        centers::Union{NTuple{2,Int},Nothing}=nothing,
        atoms::Union{NTuple{2,Int},Nothing}=nothing,
        orbitals::Union{NTuple{2,Int},Subscript,Nothing}=nothing
        ) -> Couplings{ID{SCID,2},SpinCoupling{Rational{Int}/Int,<:ID{SCID,2}}}
"""
function DM(tag::Char;centers::Union{NTuple{2,Int},Nothing}=nothing,atoms::Union{NTuple{2,Int},Nothing}=nothing,orbitals::Union{NTuple{2,Int},Subscript,Nothing}=nothing)
    @assert tag in ('x','y','z') "DM error: not supported input tag($tag)."
    t1,t2=tag=='x' ? ('y','z') : tag=='y' ? ('z','x') : ('x','y')
    sc1=SpinCoupling{2}(1,centers=centers,atoms=atoms,orbitals=orbitals,tags=(t1,t2))
    sc2=SpinCoupling{2}(-1,centers=centers,atoms=atoms,orbitals=orbitals,tags=(t2,t1))
    return Couplings(sc1,sc2)
end

scsinglewrapper(::Nothing)=nothing
scsinglewrapper(value::Int)=(value,)

"""
    Sˣ(;atom::Union{Int,Nothing}=nothing,orbital::Union{Int,Nothing}=nothing) -> Couplings{ID{SCID,1},SpinCoupling{Rational{Int}/Int,<:ID{SCID,1}}}

The single Sˣ coupling.
"""
function Sˣ(;atom::Union{Int,Nothing}=nothing,orbital::Union{Int,Nothing}=nothing)
    Couplings(SpinCoupling{1}(1,tags=('x',),atoms=scsinglewrapper(atom),orbitals=scsinglewrapper(orbital)))
end

"""
    Sʸ(;atom::Union{Int,Nothing}=nothing,orbital::Union{Int,Nothing}=nothing) -> Couplings{ID{SCID,1},SpinCoupling{Rational{Int}/Int,<:ID{SCID,1}}}

The single Sʸ coupling.
"""
function Sʸ(;atom::Union{Int,Nothing}=nothing,orbital::Union{Int,Nothing}=nothing)
    Couplings(SpinCoupling{1}(1,tags=('y',),atoms=scsinglewrapper(atom),orbitals=scsinglewrapper(orbital)))
end

"""
    Sᶻ(;atom::Union{Int,Nothing}=nothing,orbital::Union{Int,Nothing}=nothing) -> Couplings{ID{SCID,1},SpinCoupling{Rational{Int}/Int,<:ID{SCID,1}}}

The single Sᶻ coupling.
"""
function Sᶻ(;atom::Union{Int,Nothing}=nothing,orbital::Union{Int,Nothing}=nothing)
    Couplings(SpinCoupling{1}(1,tags=('z',),atoms=scsinglewrapper(atom),orbitals=scsinglewrapper(orbital)))
end

"""
    sc"..." -> SpinCoupling

Construct a SpinCoupling from a literal string.
"""
macro sc_str(str::String)
    ps=split(str," with ")
    condition=length(ps)==2 ? Meta.parse(ps[2]) : :nothing
    ps=split(ps[1],' ')
    coeff=eval(Meta.parse(ps[1]))
    tags=Tuple(replace(ps[2],"S"=>""))
    attrpairs=Any[:tags=>tags]
    if length(ps)==3
        components=split(ps[3],'@')
        centers=length(components)==2 ? couplingcenters(components[2]) : nothing
        centers===nothing || @assert length(centers)==length(tags) "@sc_str error: dismatched ranks."
        push!(attrpairs,:centers=>centers)
        if length(components[1])>0
            for component in split(components[1],'⊗')
                attrname,attrvalue=sccomponent(component)
                if isa(attrvalue,Expr)
                    @assert attrname==:orbitals "@sc_str error: wrong input pattern."
                    @assert length(tags)==length(attrvalue.args) "@sc_str error: dismatched ranks."
                    push!(attrpairs,attrname=>eval(subscriptexpr(attrvalue,condition)))
                else
                    @assert length(tags)==length(attrvalue) "@sc_str error: dismatched ranks."
                    push!(attrpairs,attrname=>attrvalue)
                end
            end
        end
    end
    return SpinCoupling{length(tags)}(coeff;attrpairs...)
end
function sccomponent(str::AbstractString)
    @assert str[3]=='(' && str[end]==')' "sccomponent error: wrong input pattern."
    attrname=str[1:2]=="sl" ? :atoms : str[1:2]=="ob" ? :orbitals : error("sccomponent error: wrong input pattern.")
    expr=Meta.parse(str[3:end])
    attrvalue=isa(expr,Expr) ? (all(isa(arg,Int) for arg in expr.args) ? Tuple(expr.args) : expr) : (expr,)
    return attrname=>attrvalue
end

function scpairs(str::AbstractString,::Val{R}) where R
    attrpairs=Pair{Symbol,NTuple{R,Int}}[]
    if length(str)>0
        components=split(str,'@')
        length(components)==2 && push!(attrpairs,:centers=>couplingcenters(components[2]))
        if length(components[1])>0
            for component in split(components[1],'⊗')
                attrname,attrvalue=sccomponent(component)
                @assert isa(attrvalue,NTuple{R,Int}) "scpairs error: wrong input pattern."
                push!(attrpairs,attrname=>attrvalue)
            end
        end
    end
    return attrpairs
end

"""
    heisenberg"sl(a₁,a₂)⊗ob(o₁,o₂)@(c₁,c₂)" -> Couplings

The Heisenberg couplings.
"""
macro heisenberg_str(str::String)
    ps=split(str," ")
    @assert length(ps)∈(1,2) "heisenberg_str error: not supported pattern."
    mode,str=length(ps)==2 ? (ps[1],ps[2]) : ps[1]∈("+-z","xyz") ?  (ps[1],"") : ("+-z",ps[1])
    @assert mode=="+-z" || mode=="xyz" "heisenberg_str error: not supported mode($mode)."
    attrpairs=scpairs(str,Val(2))
    if mode=="+-z"
        sc1=SpinCoupling{2}(1//2;:tags=>('+','-'),attrpairs...)
        sc2=SpinCoupling{2}(1//2;:tags=>('-','+'),attrpairs...)
        sc3=SpinCoupling{2}(1//1;:tags=>('z','z'),attrpairs...)
    else
        sc1=SpinCoupling{2}(1;:tags=>('x','x'),attrpairs...)
        sc2=SpinCoupling{2}(1;:tags=>('y','y'),attrpairs...)
        sc3=SpinCoupling{2}(1;:tags=>('z','z'),attrpairs...)
    end
    return Couplings(sc1,sc2,sc3)
end

"""
    ising"x sl(a₁,a₂)⊗ob(o₁,o₂)@(c₁,c₂)" -> Couplings
    ising"y sl(a₁,a₂)⊗ob(o₁,o₂)@(c₁,c₂)" -> Couplings
    ising"z sl(a₁,a₂)⊗ob(o₁,o₂)@(c₁,c₂)" -> Couplings

The Ising couplings.
"""
macro ising_str(str::String)
    @assert str[1] in ('x','y','z') "@ising_str error: wrong input pattern."
    attrpairs=length(str)>1 ? (@assert str[2]==' ' "@ising_str error: wrong input pattern."; scpairs(str[3:end],Val(2))) : Pair{Symbol,NTuple{2,Int}}[]
    return Couplings(SpinCoupling{2}(1;:tags=>(str[1],str[1]),attrpairs...))
end

"""
    gamma"x sl(a₁,a₂)⊗ob(o₁,o₂)@(c₁,c₂)" -> Couplings
    gamma"y sl(a₁,a₂)⊗ob(o₁,o₂)@(c₁,c₂)" -> Couplings
    gamma"z sl(a₁,a₂)⊗ob(o₁,o₂)@(c₁,c₂)" -> Couplings

The Gamma couplings.
"""
macro gamma_str(str::String)
    @assert str[1] in ('x','y','z') "@gamma_str error: wrong input pattern."
    t1,t2=str[1]=='x' ? ('y','z') : str[1]=='y' ? ('z','x') : ('x','y')
    attrpairs=length(str)>1 ? (@assert str[2]==' ' "@gamma_str error: wrong input pattern."; scpairs(str[3:end],Val(2))) : Pair{Symbol,NTuple{2,Int}}[]
    sc1=SpinCoupling{2}(1;:tags=>(t1,t2),attrpairs...)
    sc2=SpinCoupling{2}(1;:tags=>(t2,t1),attrpairs...)
    return Couplings(sc1,sc2)
end

"""
    dm"x sl(a₁,a₂)⊗ob(o₁,o₂)@(c₁,c₂)" -> Couplings
    dm"y sl(a₁,a₂)⊗ob(o₁,o₂)@(c₁,c₂)" -> Couplings
    dm"z sl(a₁,a₂)⊗ob(o₁,o₂)@(c₁,c₂)" -> Couplings

The DM couplings.
"""
macro dm_str(str::String)
    @assert str[1] in ('x','y','z') "@dm_str error: wrong input pattern."
    t1,t2=str[1]=='x' ? ('y','z') : str[1]=='y' ? ('z','x') : ('x','y')
    attrpairs=length(str)>1 ? (@assert str[2]==' ' "@dm_str error: wrong input pattern."; scpairs(str[3:end],Val(2))) : Pair{Symbol,NTuple{2,Int}}[]
    sc1=SpinCoupling{2}(1;:tags=>(t1,t2),attrpairs...)
    sc2=SpinCoupling{2}(-1;:tags=>(t2,t1),attrpairs...)
    return Couplings(sc1,sc2)
end


"""
    sˣ"sl(a)⊗ob(o)" -> Couplings
    sʸ"sl(a)⊗ob(o)" -> Couplings
    sᶻ"sl(a)⊗ob(o)" -> Couplings

The single Sˣ/Sʸ/Sᶻ coupling.
"""
macro sˣ_str(str::String) Couplings(SpinCoupling{1}(1;:tags=>('x',),scpairs(str,Val(1))...)) end
macro sʸ_str(str::String) Couplings(SpinCoupling{1}(1;:tags=>('y',),scpairs(str,Val(1))...)) end
macro sᶻ_str(str::String) Couplings(SpinCoupling{1}(1;:tags=>('z',),scpairs(str,Val(1))...)) end

"""
    SpinTerm{R}(id::Symbol,value::Any,bondkind::Any;
                couplings::Union{Function,Coupling,Couplings},
                amplitude::Union{Function,Nothing}=nothing,
                modulate::Union{Function,Bool}=false,
                ) where R

Spin term.

Type alias for `Term{'B',:SpinTerm,R,id,V,<:Any,<:TermCouplings,<:TermAmplitude,<:Union{TermModulate,Nothing}}`.
"""
const SpinTerm{R,id,V,B<:Any,C<:TermCouplings,A<:TermAmplitude,M<:Union{TermModulate,Nothing}}=Term{'B',:SpinTerm,R,id,V,B,C,A,M}
function SpinTerm{R}(id::Symbol,value::Any,bondkind::Any;
                    couplings::Union{Function,Coupling,Couplings},
                    amplitude::Union{Function,Nothing}=nothing,
                    modulate::Union{Function,Bool}=false,
                    ) where R
    isa(couplings,TermCouplings) || (couplings=TermCouplings(couplings))
    isa(amplitude,TermAmplitude) || (amplitude=TermAmplitude(amplitude))
    isa(modulate,TermModulate) || (modulate=modulate===false ? nothing : TermModulate(id,modulate))
    Term{'B',:SpinTerm,R,id}(value,bondkind,couplings,amplitude,modulate,1)
end
abbr(::Type{<:SpinTerm})=:sp
isHermitian(::Type{<:SpinTerm})=true

"""
    otype(T::Type{<:SpinTerm},I::Type{<:OID})

Get the operator type of a spin term.
"""
otype(T::Type{<:SpinTerm},I::Type{<:OID})=SOperator{T|>valtype,ID{I,T|>rank}}

end # module
