using StaticArrays: SVector
using Hamiltonian.Essentials.Terms
using Hamiltonian.Essentials.Spatials: Point,PID
using Hamiltonian.Essentials.DegreesOfFreedom: IDFConfig
using Hamiltonian.Essentials.FockPackage: FIndex,Fock,σˣ,σᶻ
using Hamiltonian.Mathematics.AlgebraOverFields: ID

@testset "OID" begin
    oid=OID(FIndex(1,1,1,1,1),rcoord=SVector(0.0,-0.0),icoord=nothing,seq=1)
    @test oid'==OID(FIndex(1,1,1,1,2),rcoord=SVector(0.0,0.0),icoord=nothing,seq=1)
    @test hash(oid,UInt(1))==hash(OID(FIndex(1,1,1,1,1),rcoord=SVector(0.0,0.0),icoord=nothing,seq=1),UInt(1))
    @test propertynames(ID{2,OID},true)==(:contents,:indexes,:rcoords,:icoords)
    @test propertynames(ID{2,OID},false)==(:indexes,:rcoords,:icoords)
end

struct TOperator{N,V<:Number,I<:ID} <: Operator{N,V,I}
    value::V
    id::I
    TOperator(value::Number,id::ID{N,<:OID}) where N=new{N,typeof(value),typeof(id)}(value,id)
end

@testset "Operator" begin
    opt=TOperator(1.0im,(FIndex(1,1,1,2,2),FIndex(1,1,1,1,1)),rcoords=(SVector(1.0,0.0),SVector(0.0,0.0)),seqs=(2,1))
    @test opt'==TOperator(-1.0im,(FIndex(1,1,1,1,2),FIndex(1,1,1,2,1)),rcoords=(SVector(0.0,0.0),SVector(1.0,0.0)),seqs=(1,2))
    @test isHermitian(opt)==false
    opt=TOperator(1.0,(FIndex(1,1,1,1,2),FIndex(1,1,1,1,1)),rcoords=(SVector(0.0,0.0),SVector(0.0,0.0)),seqs=(1,1))
    @test opt'==opt
    @test isHermitian(opt)==true
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

@testset "Term" begin
    point=Point(PID(1,1),(0.0,0.0),(0.0,0.0))
    config=IDFConfig(pid->Fock(norbital=2,nspin=2,nnambu=2),Fock,[PID(1,1)])

    term=Term{'F',:Term,2}(:mu,2.0,0,couplings=σˣ("sp"),amplitude=(bond->2),modulate=false)
    @test term|>statistics==term|>typeof|>statistics=='F'
    @test term|>species==term|>typeof|>species==:Term
    @test term|>rank==term|>typeof|>rank==2
    @test term|>abbr==:tm
    @test term==deepcopy(term)
    @test isequal(term,deepcopy(term))
    @test string(term)=="Term{F,2}(id=mu,value=2.0,neighbor=0)"
    @test repr(term,point,config)=="tm: 2.0 sp(2:1)\ntm: 2.0 sp(1:2)"
    @test +term==term
    @test -term==term*(-1)==replace(term,factor=-term.factor)
    @test 2*term==term*2==replace(term,factor=term.factor*2)
    @test term/2==term*(1/2)==replace(term,factor=term.factor/2)

    term=Term{'F',:Term,2}(:mu,2.0,0,couplings=bond->σᶻ("sp")⊗σˣ("ob"),amplitude=bond->2,modulate=true)
    @test repr(term,point,config)=="tm: 2.0 ob(2:1)⊗sp(2:2)\ntm: -2.0 ob(1:2)⊗sp(1:1)\ntm: -2.0 ob(2:1)⊗sp(1:1)\ntm: 2.0 ob(1:2)⊗sp(2:2)"
    @test one(term)==replace(term,value=1.0)
    @test zero(term)==replace(term,value=0.0)
    @test term.modulate(mu=4.0)==4.0
    @test term.modulate(t=1.0)==nothing
end
