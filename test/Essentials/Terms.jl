using Printf: @sprintf
using StaticArrays: SVector
using Hamiltonian.Essentials.Terms
using Hamiltonian.Essentials.Spatials: Point,PID,Bond
using Hamiltonian.Essentials.DegreesOfFreedom: IDFConfig,Table,IID,Index,Internal,FilteredAttributes,Coupling,Couplings
using Hamiltonian.Prerequisites: Float,decimaltostr
using Hamiltonian.Mathematics.AlgebraOverFields: ID,SimpleID
import Hamiltonian.Prerequisites.Interfaces: dimension,expand

struct TID <: IID nambu::Int end
Base.adjoint(sl::TID)=TID(3-sl.nambu)

struct TIndex{S} <: Index{PID{S},TID}
    scope::S
    site::Int
    nambu::Int
end
Base.fieldnames(::Type{<:TIndex})=(:scope,:site,:nambu)
Base.union(::Type{P},::Type{I}) where {P<:PID,I<:TID}=TIndex{fieldtype(P,:scope)}

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
Base.repr(tc::TCoupling)=@sprintf "%s ph(%s)@(%s)" decimaltostr(tc.value) join(tc.id.nambus,':') join(tc.id.centers,'-')
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

@testset "OID" begin
    oid=OID(TIndex(1,1,1),rcoord=SVector(0.0,-0.0),icoord=nothing,seq=1)
    @test oid'==OID(TIndex(1,1,2),rcoord=SVector(0.0,0.0),icoord=nothing,seq=1)
    @test hash(oid,UInt(1))==hash(OID(TIndex(1,1,1),rcoord=SVector(0.0,0.0),icoord=nothing,seq=1),UInt(1))
    @test propertynames(ID{<:NTuple{2,OID}},true)==(:contents,:indexes,:rcoords,:icoords)
    @test propertynames(ID{<:NTuple{2,OID}},false)==(:indexes,:rcoords,:icoords)
    @test fieldnames(OID)==(:index,:rcoord,:icoord,:seq)
    @test string(oid)=="OID(TIndex(1,1,1),[0.0,0.0],:,1)"
    @test ID(oid',oid)'==ID(oid',oid)
    @test isHermitian(ID(oid',oid))==true
    @test isHermitian(ID(oid,oid))==false
end

@testset "Operator" begin
    opt=TOperator(1.0im,(TIndex(1,2,2),TIndex(1,1,1)),rcoords=(SVector(1.0,0.0),SVector(0.0,0.0)),seqs=(2,1))
    @test opt'==TOperator(-1.0im,(TIndex(1,1,2),TIndex(1,2,1)),rcoords=(SVector(0.0,0.0),SVector(1.0,0.0)),seqs=(1,2))
    @test isHermitian(opt)==false
    @test string(opt)=="TOperator(value=1.0im,id=ID(OID(TIndex(1,2,2),[1.0,0.0],:,2),OID(TIndex(1,1,1),[0.0,0.0],:,1)))"

    opt=TOperator(1.0,(TIndex(1,1,2),TIndex(1,1,1)),rcoords=(SVector(0.5,0.5),SVector(0.5,0.5)),icoords=(SVector(1.0,1.0),SVector(1.0,1.0)),seqs=(1,1))
    @test opt'==opt
    @test isHermitian(opt)==true

    opt=TOperator(1.0,(TIndex(1,1,2),TIndex(1,1,1)),rcoords=(SVector(0.5,0.5),SVector(0.0,0.5)),icoords=(SVector(1.0,1.0),SVector(0.0,1.0)),seqs=(1,1))
    @test rcoord(opt)==SVector(0.5,0.0)
    @test icoord(opt)==SVector(1.0,0.0)

    opt=TOperator(1.0,ID(OID(TIndex(1,1,2),SVector(0.5,0.0),SVector(1.0,0.0),1)))
    @test rcoord(opt)==SVector(0.5,0.0)
    @test icoord(opt)==SVector(1.0,0.0)
end

@testset "Operators" begin
    opt1=TOperator(1.0im,(TIndex(1,2,2),TIndex(1,1,1)),rcoords=(SVector(1.0,0.0),SVector(0.0,0.0)),seqs=(2,1))
    opt2=TOperator(1.0,(TIndex(1,1,2),TIndex(1,1,1)),rcoords=(SVector(0.0,0.0),SVector(0.0,0.0)),seqs=(1,1))
    opts=Operators(opt1,opt2)
    @test opts'==Operators(opt1',opt2')
    @test opts'+opts==Operators(opt1,opt1',opt2*2)
    @test isHermitian(opts)==false
    @test isHermitian(opts'+opts)==true
end

@testset "TermFunction" begin
    ta=TermAmplitude()
    @test ta()==1

    ta=TermAmplitude(x->x+3.0)
    @test ta(1.0)==4.0

    tcs=TCoupling(1.0,ID(TCID(1,1),TCID(2,1)))+TCoupling(2.0,ID(TCID(1,2),TCID(2,2)))
    termcouplings=TermCouplings(tcs)
    @test termcouplings()==tcs
    @test termcouplings|>rank==termcouplings|>typeof|>rank==2

    tcs1=TCoupling(1.0,ID(TCID(1,1),TCID(2,1)))+TCoupling(1.0,ID(TCID(1,2),TCID(2,2)))
    tcs2=TCoupling(1.0,ID(TCID(1,2),TCID(2,1)))+TCoupling(1.0,ID(TCID(2,1),TCID(1,2)))
    termcouplings=TermCouplings((tcs1,tcs2),i->(i-1)%2+1)
    @test termcouplings==TermCouplings((termcouplings.candidates,termcouplings.choice))
    @test isequal(termcouplings,TermCouplings((termcouplings.candidates,termcouplings.choice)))
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
    term=Term{'F',:Term}(:mu,1.5,0,couplings=TCoupling(1.0,ID(TCID(1,2),TCID(1,1))),amplitude=(bond->3),modulate=false)
    @test term|>statistics==term|>typeof|>statistics=='F'
    @test term|>species==term|>typeof|>species==:Term
    @test term|>valtype==term|>typeof|>valtype==Float
    @test term|>rank==term|>typeof|>rank==2
    @test term|>abbr==term|>typeof|>abbr==:tm
    @test term==deepcopy(term)
    @test isequal(term,deepcopy(term))
    @test string(term)=="Term{F}(id=mu,value=1.5,neighbor=0,factor=1.0)"
    @test repr(term,point,config)=="tm: 4.5 ph(2:1)@(1-1)"
    @test +term==term
    @test -term==term*(-1)==replace(term,factor=-term.factor)
    @test 2*term==term*2==replace(term,factor=term.factor*2)
    @test term/2==term*(1/2)==replace(term,factor=term.factor/2)

    p1=Point(PID(1,1),(0.0,0.0),(0.0,0.0))
    p2=Point(PID(1,2),(0.0,0.0),(0.0,0.0))
    config=IDFConfig{TFock}(pid->TFock(),[PID(1,1),PID(1,2)])
    tcs1=Couplings(TCoupling(1.0,ID(TCID(1,2),TCID(1,2))))
    tcs2=Couplings(TCoupling(1.0,ID(TCID(1,1),TCID(1,1))))
    term=Term{'F',:Term}(:mu,1.5,0,couplings=((tcs1,tcs2),bond->bond.pid.site%2+1),amplitude=bond->3,modulate=true)
    @test repr(term,p1,config)=="tm: 4.5 ph(1:1)@(1-1)"
    @test repr(term,p2,config)=="tm: 4.5 ph(2:2)@(1-1)"
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
    O=OID{TIndex{Int},SVector{2,Float},Nothing,Int}
    term=Term{'F',:Term}(:mu,1.5,0,couplings=TCoupling(1.0,ID(TCID(1,2),TCID(1,1))),amplitude=bond->3.0,modulate=true)
    operators=Operators(TOperator(+2.25,(TIndex(1,1,2),TIndex(1,1,1)),rcoords=(SVector(0.0,0.0),SVector(0.0,0.0)),seqs=(1,1)))
    @test expand(TOperator{2,Float,ID{NTuple{2,O}}},term,point,config,table,nothing)==operators
    @test expand(TOperator{2,Float,ID{NTuple{2,O}}},term,point,config,table,true)==operators
    @test expand(TOperator{2,Float,ID{NTuple{2,O}}},term,point,config,table,false)==operators*2

    bond=Bond(1,Point(PID('b',2),(1.5,1.5),(1.0,1.0)),Point(PID('a',1),(0.5,0.5),(0.0,0.0)))
    config=IDFConfig{TFock}(pid->TFock(),[PID('a',1),PID('b',2)])
    table=Table(config,by=filter(attr->attr≠:nambu,FilteredAttributes(TIndex)))
    O=OID{TIndex{Char},SVector{2,Float},SVector{2,Float},Int}
    term=Term{'F',:Term}(:t,1.5,1,couplings=TCoupling(1.0,ID(TCID(1,2),TCID(2,1))),amplitude=bond->3.0,modulate=true)
    operators=Operators(TOperator(4.5,(TIndex('a',1,2),TIndex('b',2,1)),rcoords=(SVector(0.5,0.5),SVector(1.5,1.5)),icoords=(SVector(0.0,0.0),SVector(1.0,1.0)),seqs=(1,2)))
    @test expand(TOperator{2,Float,ID{NTuple{2,O}}},term,bond,config,table,nothing)==operators
    @test expand(TOperator{2,Float,ID{NTuple{2,O}}},term,bond,config,table,true)==operators/2
    @test expand(TOperator{2,Float,ID{NTuple{2,O}}},term,bond,config,table,false)==operators
end
