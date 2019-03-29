using Test
using Printf: @sprintf
using LinearAlgebra: dot
using StaticArrays: SVector
using QuantumLattices.Essentials.Terms
using QuantumLattices.Essentials.Spatials: Point,PID,Bond,decompose
using QuantumLattices.Essentials.DegreesOfFreedom: IDFConfig,Table,IID,Index,Internal,FilteredAttributes,OID,Operator,Operators
using QuantumLattices.Interfaces: rank,update!,id
using QuantumLattices.Prerequisites: Float,decimaltostr
using QuantumLattices.Mathematics.AlgebraOverFields: ID,SimpleID
import QuantumLattices.Interfaces: dimension,expand
import QuantumLattices.Essentials.DegreesOfFreedom: twist,isHermitian,otype
import QuantumLattices.Essentials.Terms: couplingcenter,couplingcenters,abbr

struct TID <: IID nambu::Int end
Base.adjoint(sl::TID)=TID(3-sl.nambu)

struct TIndex{S} <: Index{PID{S},TID}
    scope::S
    site::Int
    nambu::Int
end
Base.fieldnames(::Type{<:TIndex})=(:scope,:site,:nambu)
Base.union(::Type{P},::Type{I}) where {P<:PID,I<:TID}=TIndex{fieldtype(P,:scope)}
function twist(id::OID{<:TIndex},vectors::AbstractVector{<:AbstractVector{Float}},values::AbstractVector{Float})
    phase=  length(vectors)==1 ? exp(2.0im*pi*dot(decompose(id.icoord,vectors[1]),values)) :
            length(vectors)==2 ? exp(2.0im*pi*dot(decompose(id.icoord,vectors[1],vectors[2]),values)) :
            length(vectors)==3 ? exp(2.0im*pi*dot(decompose(id.icoord,vectors[1],vectors[2],vectors[3]),values)) :
            error("twist error: not supported number of input basis vectors.")
    id.index.nambu==1 ? phase : conj(phase)
end

struct TFock <: Internal{TID} end
dimension(f::TFock)=2
Base.getindex(f::TFock,i::Int)=TID(i)

struct TCID <: SimpleID
    center::Int
    nambu::Int
end
Base.fieldnames(::Type{<:TCID})=(:center,:nambu)

struct TCoupling{N,V<:Number,I<:ID{<:NTuple{N,TCID}}} <: Coupling{N,V,I}
    value::V
    id::I
end
Base.repr(tc::TCoupling)=@sprintf "%s ph(%s)@(%s)" decimaltostr(tc.value) join(tc.id.nambus,',') join(tc.id.centers,',')
function expand(tc::TCoupling,pids::NTuple{R,PID},focks::NTuple{R,TFock},species::Union{Val{S},Nothing}=nothing) where {R,S}
    nambus=tc.id.nambus
    pids=NTuple{rank(tc),eltype(pids)}(pids[tc.id.centers[i]] for i=1:rank(tc))
    return ((tc.value,NTuple{rank(tc),TIndex{fieldtype(pids|>eltype,:scope)}}(TIndex(pids[i].scope,pids[i].site,nambus[i]) for i=1:rank(tc))),)
end

struct TOperator{N,V<:Number,I<:ID{<:NTuple{N,OID}}} <: Operator{N,V,I}
    value::V
    id::I
    TOperator(value::Number,id::ID{<:NTuple{N,OID}}) where N=new{N,typeof(value),typeof(id)}(value,id)
end

abbr(::Type{<:Term{ST,:TermMu}}) where ST=:tmu
isHermitian(::Type{<:Term{ST,:TermMu}}) where ST=true
otype(T::Type{<:Term{ST,:TermMu}},I::Type{<:OID}) where ST=TOperator{T|>rank,T|>valtype,ID{NTuple{T|>rank,I}}}

abbr(::Type{<:Term{ST,:TermHopping}}) where ST=:thp
isHermitian(::Type{<:Term{ST,:TermHopping}}) where ST=false
otype(T::Type{<:Term{ST,:TermHopping}},I::Type{<:OID}) where ST=TOperator{T|>rank,T|>valtype,ID{NTuple{T|>rank,I}}}

@testset "Subscript" begin
    sub=@subscript (x1,4,4,x2) with x1<x2
    @test sub==deepcopy(sub)
    @test isequal(sub,deepcopy(sub))
    @test rank(sub)==rank(typeof(sub))==2
    @test dimension(sub)==dimension(typeof(sub))==4
    @test string(sub)=="(x1,4,4,x2) with $(sub.identifier)"
    @test sub(Val('M'),1,2)==(1,4,4,2)
    @test sub(Val('C'),1,2)==true
    @test sub(Val('C'),2,1)==false

    sub=@subscript (x1,x2,x1,x2)
    @test sub(Val('M'),1,2)==(1,2,1,2)
    @test sub(Val('C'),1,2)==true

    sub=Subscript{4}()
    @test rank(sub)==1 && dimension(sub)==4
    @test string(sub)=="(*,*,*,*)"
    @test sub(Val('M'),2)==(2,2,2,2)
    @test sub(Val('C'),2)==true

    sub=Subscript((1,2,2,1))
    @test rank(sub)==0 && dimension(sub)==4
    @test string(sub)=="(1,2,2,1)"
    @test sub(Val('M'))==(1,2,2,1)
    @test sub(Val('C'))==true
end

@testset "Subscripts" begin
    sub1=@subscript (x1,2) with x1<2
    sub2=@subscript (y1,y2,y1,y2) with y1<y2
    subs=Subscripts(sub1,sub2)
    @test rank(subs)==rank(typeof(subs))==3
    @test rank(subs,1)==rank(typeof(subs),1)==1
    @test rank(subs,2)==rank(typeof(subs),2)==2
    @test dimension(subs)==dimension(typeof(subs))==6
    @test dimension(subs,1)==dimension(typeof(subs),1)==2
    @test dimension(subs,2)==dimension(typeof(subs),2)==4
    @test subs(Val('M'),(1,2,3))==(1,2,2,3,2,3)
    @test subs(Val('C'),(1,2,3))==true
    @test subs(Val('C'),(2,2,3))==false
    @test subs(Val('C'),(1,3,2))==false
    @test subs(Val('C'),(2,3,2))==false

    sub3=Subscript((6,6))
    @test subs==sub1*sub2
    @test subs*sub3==Subscripts(sub1,sub2,sub3)
    @test sub3*subs==Subscripts(sub3,sub1,sub2)
    @test subs*Subscripts(sub3)==Subscripts(sub1,sub2,sub3)
    @test expand(subs*sub3,(3,3,2,3,2,3,8,9))|>collect==[(1,2,1,2,1,2,6,6),(1,2,1,3,1,3,6,6),(1,2,2,3,2,3,6,6)]
end

@testset "centerrelated" begin
    @test couplingcenter(Coupling,2,2,Val(1))==1
    @test couplingcenter(Coupling,1,2,Val(1))==1
    @test couplingcenters(Coupling,('*','*'),Val(1))==(1,1)
    @test couplingcenters(Coupling,(1,2,3),Val(3))==(1,2,3)
end

@testset "TermFunction" begin
    ta=TermAmplitude()
    @test ta()==1

    ta=TermAmplitude(x->x+3.0)
    @test ta(1.0)==4.0

    tcs=TCoupling(1.0,ID(TCID(1,1),TCID(2,1)))+TCoupling(2.0,ID(TCID(1,2),TCID(2,2)))
    termcouplings=TermCouplings(tcs)
    @test termcouplings==deepcopy(TermCouplings(tcs))
    @test isequal(termcouplings,deepcopy(TermCouplings(tcs)))
    @test termcouplings()==tcs

    tcs1=TCoupling(1.0,ID(TCID(1,1),TCID(2,1)))+TCoupling(1.0,ID(TCID(1,2),TCID(2,2)))
    tcs2=TCoupling(1.0,ID(TCID(1,2),TCID(2,1)))+TCoupling(1.0,ID(TCID(2,1),TCID(1,2)))
    termcouplings=TermCouplings(i->(tcs1,tcs2)[(i-1)%2+1])
    @test termcouplings(1)==tcs1
    @test termcouplings(2)==tcs2

    termmodulate=TermModulate(:t)
    @test termmodulate(t=1)==1
    @test termmodulate(mu=1)==nothing
    @test termmodulate()==nothing

    termmodulate=TermModulate(:t,t->t*2.0)
    @test termmodulate(2)==4
end

@testset "Term" begin
    point=Point(PID(1,1),(0.0,0.0),(0.0,0.0))
    config=IDFConfig{TFock}(pid->TFock(),[PID(1,1)])
    term=Term{'F',:TermMu,2}(:mu,1.5,0,couplings=TCoupling(1.0,ID(TCID(1,2),TCID(1,1))),amplitude=(bond->3),modulate=false)
    @test term|>statistics==term|>typeof|>statistics=='F'
    @test term|>species==term|>typeof|>species==:TermMu
    @test term|>id==term|>typeof|>id==:mu
    @test term|>valtype==term|>typeof|>valtype==Float
    @test term|>rank==term|>typeof|>rank==2
    @test term|>abbr==term|>typeof|>abbr==:tmu
    @test term|>isHermitian==term|>typeof|>isHermitian==true
    @test term==deepcopy(term)
    @test isequal(term,deepcopy(term))
    @test string(term)=="TermMu{2F}(id=mu,value=1.5,neighbor=0,factor=1.0)"
    @test repr(term,point,config)=="tmu: 4.5 ph(2,1)@(1,1)"
    @test +term==term
    @test -term==term*(-1)==replace(term,factor=-term.factor)
    @test 2*term==term*2==replace(term,factor=term.factor*2)
    @test term/2==term*(1/2)==replace(term,factor=term.factor/2)

    p1=Point(PID(1,1),(0.0,0.0),(0.0,0.0))
    p2=Point(PID(1,2),(0.0,0.0),(0.0,0.0))
    config=IDFConfig{TFock}(pid->TFock(),[PID(1,1),PID(1,2)])
    tcs1=Couplings(TCoupling(1.0,ID(TCID(1,2),TCID(1,2))))
    tcs2=Couplings(TCoupling(1.0,ID(TCID(1,1),TCID(1,1))))
    term=Term{'F',:TermMu,2}(:mu,1.5,0,couplings=bond->(tcs1,tcs2)[bond.pid.site%2+1],amplitude=bond->3,modulate=true)
    @test repr(term,p1,config)=="tmu: 4.5 ph(1,1)@(1,1)"
    @test repr(term,p2,config)=="tmu: 4.5 ph(2,2)@(1,1)"
    @test one(term)==replace(term,value=1.0)
    @test zero(term)==replace(term,value=0.0)
    @test term.modulate(mu=4.0)==4.0
    @test term.modulate(t=1.0)==nothing
    @test update!(term,mu=4.25)==replace(term,value=4.25)
    @test term.value==4.25
end

@testset "expand" begin
    point=Point(PID(1,1),(0.0,0.0),(0.0,0.0))
    config=IDFConfig{TFock}(pid->TFock(),[PID(1,1)])
    table=Table(config,by=filter(attr->attr≠:nambu,FilteredAttributes(TIndex)))
    term=Term{'F',:TermMu,2}(:mu,1.5,0,couplings=TCoupling(1.0,ID(TCID(1,2),TCID(1,1))),amplitude=bond->3.0,modulate=true)
    operators=Operators(TOperator(+2.25,(TIndex(1,1,2),TIndex(1,1,1)),rcoords=(SVector(0.0,0.0),SVector(0.0,0.0)),icoords=(SVector(0.0,0.0),SVector(0.0,0.0)),seqs=(1,1)))
    @test expand(term,point,config,table,true)==operators
    @test expand(term,point,config,table,false)==operators*2

    bond=Bond(1,Point(PID('b',2),(1.5,1.5),(1.0,1.0)),Point(PID('a',1),(0.5,0.5),(0.0,0.0)))
    config=IDFConfig{TFock}(pid->TFock(),[PID('a',1),PID('b',2)])
    table=Table(config,by=filter(attr->attr≠:nambu,FilteredAttributes(TIndex)))
    term=Term{'F',:TermHopping,2}(:t,1.5,1,couplings=TCoupling(1.0,ID(TCID(1,2),TCID(2,1))),amplitude=bond->3.0,modulate=true)
    operators=Operators(TOperator(4.5,(TIndex('a',1,2),TIndex('b',2,1)),rcoords=(SVector(0.5,0.5),SVector(1.5,1.5)),icoords=(SVector(0.0,0.0),SVector(1.0,1.0)),seqs=(1,2)))
    @test expand(term,bond,config,table,true)==operators
    @test expand(term,bond,config,table,false)==operators+operators'
end

@testset "Boundary" begin
    opt=TOperator(4.5,(TIndex('a',1,2),TIndex('b',2,1)),rcoords=(SVector(0.5,0.5),SVector(1.5,1.5)),icoords=(SVector(0.0,0.0),SVector(1.0,1.0)),seqs=(1,2))
    bound=Boundary{(:θ₁,:θ₂)}([0.1,0.2],[[1.0,0.0],[0.0,1.0]],twist)
    @test bound(opt)≈replace(opt,value=4.5*exp(2im*pi*0.3))
    update!(bound,θ₁=0.3)
    @test bound(opt)≈replace(opt,value=4.5*exp(2im*pi*0.5))
end
