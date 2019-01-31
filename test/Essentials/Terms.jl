using StaticArrays: SVector
using Hamiltonian.Essentials.Terms
using Hamiltonian.Essentials.Spatials: Point,PID,Bond
using Hamiltonian.Essentials.DegreesOfFreedom: IDFConfig,Table
using Hamiltonian.Essentials.FockPackage: FIndex,Fock,σˣ,σᶻ,usualfockindextotuple
using Hamiltonian.Prerequisites: Float
using Hamiltonian.Mathematics.AlgebraOverFields: ID
import Hamiltonian.Essentials.Terms: termamplitude

@testset "OID" begin
    oid=OID(FIndex(1,1,1,1,1),rcoord=SVector(0.0,-0.0),icoord=nothing,seq=1)
    @test oid'==OID(FIndex(1,1,1,1,2),rcoord=SVector(0.0,0.0),icoord=nothing,seq=1)
    @test hash(oid,UInt(1))==hash(OID(FIndex(1,1,1,1,1),rcoord=SVector(0.0,0.0),icoord=nothing,seq=1),UInt(1))
    @test propertynames(ID{2,OID},true)==(:contents,:indexes,:rcoords,:icoords)
    @test propertynames(ID{2,OID},false)==(:indexes,:rcoords,:icoords)
    @test fieldnames(OID)==(:index,:rcoord,:icoord,:seq)
    @test string(oid)=="OID(FIndex(1,1,1,1,1),[0.0,0.0],:,1)"
    @test ID(oid',oid)'==ID(oid',oid)
    @test isHermitian(ID(oid',oid))==true
    @test isHermitian(ID(oid,oid))==false
end

struct TOperator{N,V<:Number,I<:ID{N,<:OID}} <: Operator{N,V,I}
    value::V
    id::I
    TOperator(value::Number,id::ID{N,<:OID}) where N=new{N,typeof(value),typeof(id)}(value,id)
end

@testset "Operator" begin
    opt=TOperator(1.0im,(FIndex(1,1,1,2,2),FIndex(1,1,1,1,1)),rcoords=(SVector(1.0,0.0),SVector(0.0,0.0)),seqs=(2,1))
    @test opt'==TOperator(-1.0im,(FIndex(1,1,1,1,2),FIndex(1,1,1,2,1)),rcoords=(SVector(0.0,0.0),SVector(1.0,0.0)),seqs=(1,2))
    @test isHermitian(opt)==false
    @test string(opt)=="TOperator(value=1.0im,id=ID(OID(FIndex(1,1,1,2,2),[1.0,0.0],:,2),OID(FIndex(1,1,1,1,1),[0.0,0.0],:,1)))"

    opt=TOperator(1.0,(FIndex(1,1,1,1,2),FIndex(1,1,1,1,1)),rcoords=(SVector(0.5,0.5),SVector(0.5,0.5)),icoords=(SVector(1.0,1.0),SVector(1.0,1.0)),seqs=(1,1))
    @test opt'==opt
    @test isHermitian(opt)==true

    opt=TOperator(1.0,(FIndex(1,1,1,1,2),FIndex(1,1,1,1,1)),rcoords=(SVector(0.5,0.5),SVector(0.0,0.5)),icoords=(SVector(1.0,1.0),SVector(0.0,1.0)),seqs=(1,1))
    @test rcoord(opt)==SVector(0.5,0.0)
    @test icoord(opt)==SVector(1.0,0.0)

    opt=TOperator(1.0,ID(OID(FIndex(1,1,1,1,2),SVector(0.5,0.0),SVector(1.0,0.0),1)))
    @test rcoord(opt)==SVector(0.5,0.0)
    @test icoord(opt)==SVector(1.0,0.0)
end

@testset "Operators" begin
    opt1=TOperator(1.0im,(FIndex(1,1,1,2,2),FIndex(1,1,1,1,1)),rcoords=(SVector(1.0,0.0),SVector(0.0,0.0)),seqs=(2,1))
    opt2=TOperator(1.0,(FIndex(1,1,1,1,2),FIndex(1,1,1,1,1)),rcoords=(SVector(0.0,0.0),SVector(0.0,0.0)),seqs=(1,1))
    opts=Operators(opt1,opt2)
    @test opts'==Operators(opt1',opt2')
    @test opts'+opts==Operators(opt1,opt1',opt2*2)
    @test isHermitian(opts)==false
    @test isHermitian(opts'+opts)==true
end

@testset "TermFunction" begin
    @test termamplitude==termamplitude
    @test isequal(termamplitude,termamplitude)
    @test termamplitude()==1

    termcouplings=TermCouplings(σˣ("sp"))
    @test termcouplings()==σˣ("sp")
    @test termcouplings|>rank==termcouplings|>typeof|>rank==2

    termcouplings=TermCouplings((σˣ("sp"),σᶻ("sp")),i->(i-1)%2+1)
    @test termcouplings(1)==σˣ("sp")
    @test termcouplings(2)==σᶻ("sp")

    termmodulate=TermModulate(:t)
    @test termmodulate(t=1)==1
    @test termmodulate(mu=1)==nothing
    @test termmodulate()==nothing
end

@testset "Term" begin
    point=Point(PID(1,1),(0.0,0.0),(0.0,0.0))
    config=IDFConfig(pid->Fock(norbital=2,nspin=2,nnambu=2),Fock,[PID(1,1)])

    term=Term{'F',:Term}(:mu,1.5,0,couplings=σˣ("sp"),amplitude=(bond->3),modulate=false)
    @test term|>statistics==term|>typeof|>statistics=='F'
    @test term|>species==term|>typeof|>species==:Term
    @test term|>valtype==term|>typeof|>valtype==Float
    @test term|>rank==term|>typeof|>rank==2
    @test term|>abbr==term|>typeof|>abbr==:tm
    @test term==deepcopy(term)
    @test isequal(term,deepcopy(term))
    @test string(term)=="Term{F}(id=mu,value=1.5,neighbor=0,factor=1.0)"
    @test repr(term,point,config)=="tm: 4.5 sp(2:1)\ntm: 4.5 sp(1:2)"
    @test +term==term
    @test -term==term*(-1)==replace(term,factor=-term.factor)
    @test 2*term==term*2==replace(term,factor=term.factor*2)
    @test term/2==term*(1/2)==replace(term,factor=term.factor/2)

    p1=Point(PID(1,1),(0.0,0.0),(0.0,0.0))
    p2=Point(PID(1,2),(0.0,0.0),(0.0,0.0))
    config=IDFConfig(pid->Fock(norbital=2,nspin=2,nnambu=2),Fock,[PID(1,1),PID(1,2)])
    term=Term{'F',:Term}(:mu,1.5,0,couplings=TermCouplings((σᶻ("sp")⊗σˣ("ob"),σˣ("sp")⊗σˣ("ob")),bond->bond.pid.site%2+1),amplitude=bond->3,modulate=true)
    @test repr(term,p1,config)=="tm: 4.5 ob(2:1)⊗sp(1:2)\ntm: 4.5 ob(2:1)⊗sp(2:1)\ntm: 4.5 ob(1:2)⊗sp(2:1)\ntm: 4.5 ob(1:2)⊗sp(1:2)"
    @test repr(term,p2,config)=="tm: 4.5 ob(2:1)⊗sp(2:2)\ntm: -4.5 ob(1:2)⊗sp(1:1)\ntm: -4.5 ob(2:1)⊗sp(1:1)\ntm: 4.5 ob(1:2)⊗sp(2:2)"
    @test one(term)==replace(term,value=1.0)
    @test zero(term)==replace(term,value=0.0)
    @test term.modulate(mu=4.0)==4.0
    @test term.modulate(t=1.0)==nothing
    @test update!(term,mu=4.25)==replace(term,value=4.25)
    @test term.value==4.25
end

@testset "expand" begin
    point=Point(PID(1,1),(0.0,0.0),(0.0,0.0))
    config=IDFConfig(pid->Fock(norbital=2,nspin=2,nnambu=2),Fock,[PID(1,1)])
    table=Table(config,by=usualfockindextotuple)
    O=OID{FIndex{Int},SVector{2,Float},Nothing,Int}
    term=Term{'F',:Term}(:mu,1.5,0,couplings=σᶻ("sp")⊗σᶻ("ob"),amplitude=bond->3.0,modulate=true)
    operators=Operators(TOperator(+2.25,(FIndex(1,1,2,2,2),FIndex(1,1,2,2,1)),rcoords=(SVector(0.0,0.0),SVector(0.0,0.0)),seqs=(4,4)),
                        TOperator(-2.25,(FIndex(1,1,2,1,2),FIndex(1,1,2,1,1)),rcoords=(SVector(0.0,0.0),SVector(0.0,0.0)),seqs=(3,3)),
                        TOperator(-2.25,(FIndex(1,1,1,2,2),FIndex(1,1,1,2,1)),rcoords=(SVector(0.0,0.0),SVector(0.0,0.0)),seqs=(2,2)),
                        TOperator(+2.25,(FIndex(1,1,1,1,2),FIndex(1,1,1,1,1)),rcoords=(SVector(0.0,0.0),SVector(0.0,0.0)),seqs=(1,1))
    )
    @test expand(TOperator{2,Float,ID{2,O,NTuple{2,O}}},term,point,config,table,nothing)==operators
    @test expand(TOperator{2,Float,ID{2,O,NTuple{2,O}}},term,point,config,table,true)==operators
    @test expand(TOperator{2,Float,ID{2,O,NTuple{2,O}}},term,point,config,table,false)==operators*2

    bond=Bond(1,Point(PID('b',2),(1.5,1.5),(1.0,1.0)),Point(PID('a',1),(0.5,0.5),(0.0,0.0)))
    config=IDFConfig(pid->Fock(norbital=2,nspin=2,nnambu=2),Fock,[PID('a',1),PID('b',2)])
    table=Table(config,by=usualfockindextotuple)
    O=OID{FIndex{Char},SVector{2,Float},SVector{2,Float},Int}
    term=Term{'F',:Term}(:t,1.5,1,couplings=σᶻ("sp")⊗σᶻ("ob"),amplitude=bond->3.0,modulate=true)
    operators=Operators(TOperator(-4.5,(FIndex('a',1,1,2,2),FIndex('b',2,1,2,1)),rcoords=(SVector(0.5,0.5),SVector(1.5,1.5)),icoords=(SVector(0.0,0.0),SVector(1.0,1.0)),seqs=(2,6)),
                        TOperator(+4.5,(FIndex('a',1,1,1,2),FIndex('b',2,1,1,1)),rcoords=(SVector(0.5,0.5),SVector(1.5,1.5)),icoords=(SVector(0.0,0.0),SVector(1.0,1.0)),seqs=(1,5)),
                        TOperator(+4.5,(FIndex('a',1,2,2,2),FIndex('b',2,2,2,1)),rcoords=(SVector(0.5,0.5),SVector(1.5,1.5)),icoords=(SVector(0.0,0.0),SVector(1.0,1.0)),seqs=(4,8)),
                        TOperator(-4.5,(FIndex('a',1,2,1,2),FIndex('b',2,2,1,1)),rcoords=(SVector(0.5,0.5),SVector(1.5,1.5)),icoords=(SVector(0.0,0.0),SVector(1.0,1.0)),seqs=(3,7))
    )
    @test expand(TOperator{2,Float,ID{2,O,NTuple{2,O}}},term,bond,config,table,nothing)==operators
    @test expand(TOperator{2,Float,ID{2,O,NTuple{2,O}}},term,bond,config,table,true)==operators/2
    @test expand(TOperator{2,Float,ID{2,O,NTuple{2,O}}},term,bond,config,table,false)==operators
end
