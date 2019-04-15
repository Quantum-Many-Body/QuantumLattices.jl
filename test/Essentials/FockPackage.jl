using Test
using StaticArrays: SVector
using QuantumLattices.Essentials.FockPackage
using QuantumLattices.Essentials.Spatials: Bond,Point,PID,rcoord,azimuthd
using QuantumLattices.Essentials.DegreesOfFreedom: Table,IDFConfig,OID,Operators,twist,oidtype,otype
using QuantumLattices.Essentials.Terms: Couplings,@subscript,statistics,abbr
using QuantumLattices.Interfaces: dims,inds,⊗,⋅,expand
using QuantumLattices.Prerequisites: Float
using QuantumLattices.Mathematics.VectorSpaces: IsMultiIndexable,MultiIndexOrderStyle
using QuantumLattices.Mathematics.AlgebraOverFields: ID,rawelement
import QuantumLattices.FockPackage: fockcouplingnambus

@testset "FID" begin
    fid=FID(orbital=1,spin=1)
    @test fid|>typeof|>fieldnames==(:orbital,:spin,:nambu)
    @test fid'==FID(1,1,2)
    @test fid''==FID(1,1,1)

    fid=FID(1,1,0)
    @test fid'==fid
end

@testset "Fock" begin
    fock=Fock(norbital=1,nspin=2,nnambu=2)
    @test IsMultiIndexable(Fock)==IsMultiIndexable(true)
    @test MultiIndexOrderStyle(Fock)==MultiIndexOrderStyle('C')
    @test dims(fock)==(1,2,2)
    @test inds(FID(1,1,1),fock)==(1,1,1)
    @test FID((1,1,1),fock)==FID(1,1,1)
    @test fock|>collect==[FID(1,1,1),FID(1,1,2),FID(1,2,1),FID(1,2,2)]
end

@testset "FIndex" begin
    @test FIndex|>fieldnames==(:scope,:site,:orbital,:spin,:nambu)
    @test union(PID{Int},FID)==FIndex{Int}
    @test twist(OID(FIndex(1,1,1,1,1),[0.0,0.0],[1.0,2.0],1),[[1.0,0.0],[0.0,1.0]],[0.1,0.0])≈exp(2im*pi*0.1)
    @test twist(OID(FIndex(1,1,1,1,2),[0.0,0.0],[1.0,2.0],1),[[1.0,0.0],[0.0,1.0]],[0.0,0.2])≈exp(-2im*pi*0.4)
end

@testset "oidtype" begin
    @test oidtype(FID,Point{2,PID{Int}},Nothing)==OID{FIndex{Int},SVector{2,Float},SVector{2,Float},Nothing}
    @test oidtype(FID,Point{2,PID{Int}},Table)==OID{FIndex{Int},SVector{2,Float},SVector{2,Float},Int}
end

@testset "FockOperator" begin
    opt=FOperator(1.0,(FIndex(1,1,1,1,2),FIndex(1,1,1,1,1)))
    @test opt|>isnormalordered
    @test opt|>statistics==opt|>typeof|>statistics=='F'

    @test rawelement(FOperator{N,<:Number,<:ID{<:NTuple{N,OID}}} where N)==FOperator
    opt=FOperator(1.0,(FIndex(1,2,1,1,2),FIndex(1,2,1,1,1),FIndex(1,1,1,1,2),FIndex(1,1,1,1,1)))
    @test opt|>isnormalordered==false

    @test rawelement(BOperator{N,<:Number,<:ID{<:NTuple{N,OID}}} where N)==BOperator
    opt=BOperator(1.0,(FIndex(1,1,1,1,2),FIndex(1,1,1,1,1)))
    @test opt|>statistics==opt|>typeof|>statistics=='B'
end

@testset "FCID" begin
    @test FCID(atom=1,orbital=1,spin=1,nambu=1)|>typeof|>fieldnames==(:center,:atom,:orbital,:spin,:nambu,:obsub,:spsub)
end

@testset "FockCoupling" begin
    @test FockCoupling{2}(1.0)|>string=="FockCoupling{2}(value=1.0)"
    @test FockCoupling{2}(1.0,atoms=(1,1))|>string=="FockCoupling{2}(value=1.0,atoms=(1,1))"
    @test FockCoupling{2}(1.0,atoms=(1,1),spins=(1,2))|>string=="FockCoupling{2}(value=1.0,atoms=(1,1),spins=(1,2))"
    @test FockCoupling{2}(2.0)|>repr=="2.0 {2}"

    fc1=FockCoupling{2}(1.5,atoms=(2,1),spins=(@subscript (x,1)),centers=(1,2))
    fc2=FockCoupling{2}(2.0,atoms=(1,2),orbitals=(@subscript (x,y) with x<y),centers=(1,2))
    @test fc1|>repr=="1.5 sl(2,1)⊗sp(x,1)@(1,2) with (*,$(fc1.id[1].spsub)) && (*,$(fc1.id[2].spsub))"
    @test fc2|>repr=="2.0 sl(1,2)⊗ob(x,y)@(1,2) with ($(fc2.id[1].obsub),*) && ($(fc2.id[2].obsub),*)"
    fc=fc1*fc2
    @test fc|>repr=="3.0 sl(2,1,1,2)⊗ob(*,*,x,y)⊗sp(x,1,*,*)@(1,2,1,2) with (*,$(fc1.id[1].spsub)) && (*,$(fc1.id[2].spsub)) && ($(fc2.id[1].obsub),*) && ($(fc2.id[2].obsub),*)"

    fc1=FockCoupling{2}(1.5,spins=(@subscript (x,1)),centers=(1,2))
    fc2=FockCoupling{2}(2.0,orbitals=(@subscript (x,y) with x<y),centers=(1,2))
    fc=fc1⊗fc2
    @test fc|>repr=="3.0 ob(x,y)⊗sp(x,1)@(1,2) with ($(fc2.id[1].obsub),$(fc1.id[2].spsub)) && ($(fc2.id[1].obsub),$(fc1.id[2].spsub))"

    fc1=FockCoupling{2}(1.5,atoms=(2,1),centers=(1,2))
    fc2=FockCoupling{2}(2.0,atoms=(1,2),centers=(1,2))
    fc=fc1⋅fc2
    @test fc|>repr=="3.0 sl(2,2)@(1,2)"

    ex=expand(FockCoupling{2}(2.0,atoms=(1,1)),PID(1,1),Fock(atom=2,norbital=2,nspin=2,nnambu=2))
    @test collect(ex)==[]

    ex=expand(FockCoupling{2}(2.0,atoms=(1,2),orbitals=(1,2)),(PID(1,1),PID(1,2)),(Fock(atom=1,norbital=2,nspin=2,nnambu=2),Fock(atom=2,norbital=2,nspin=2,nnambu=2)))
    @test IsMultiIndexable(typeof(ex))==IsMultiIndexable(true)
    @test MultiIndexOrderStyle(typeof(ex))==MultiIndexOrderStyle('C')
    @test dims(ex)==(1,2)
    @test collect(ex)==[(2.0,(FIndex(1,1,1,1,2),FIndex(1,2,2,1,1))),
                        (2.0,(FIndex(1,1,1,2,2),FIndex(1,2,2,2,1)))
    ]

    ex=expand(FockCoupling{4}(2.0,centers=(1,1,1,1),spins=(2,2,1,1)),PID(1,1),Fock(atom=1,norbital=2,nspin=2,nnambu=2))
    @test dims(ex)==(2,1)
    @test collect(ex)==[(2.0,(FIndex(1,1,1,2,2),FIndex(1,1,1,2,1),FIndex(1,1,1,1,2),FIndex(1,1,1,1,1))),
                        (2.0,(FIndex(1,1,2,2,2),FIndex(1,1,2,2,1),FIndex(1,1,2,1,2),FIndex(1,1,2,1,1)))
    ]

    ex=expand(FockCoupling{4}(2.0,orbitals=(@subscript (α,α,β,β) with α<β),spins=(2,1,1,2),nambus=(2,2,1,1)),PID(1,1),Fock(atom=1,norbital=3,nspin=2,nnambu=2))
    @test dims(ex)==(3,1)
    @test collect(ex)==[(2.0,(FIndex(1,1,1,2,2),FIndex(1,1,1,1,2),FIndex(1,1,2,1,1),FIndex(1,1,2,2,1))),
                        (2.0,(FIndex(1,1,1,2,2),FIndex(1,1,1,1,2),FIndex(1,1,3,1,1),FIndex(1,1,3,2,1))),
                        (2.0,(FIndex(1,1,2,2,2),FIndex(1,1,2,1,2),FIndex(1,1,3,1,1),FIndex(1,1,3,2,1)))
    ]

    fc1=FockCoupling{2}(1.0,spins=(2,2))
    fc2=FockCoupling{2}(-1.0,spins=(1,1))
    ex=expand(fc1*fc2,PID(1,1),Fock(atom=1,norbital=2,nspin=2,nnambu=2))
    @test dims(ex)==(4,1)
    @test collect(ex)==[(-1.0,(FIndex(1,1,1,2,2),FIndex(1,1,1,2,1),FIndex(1,1,1,1,2),FIndex(1,1,1,1,1))),
                        (-1.0,(FIndex(1,1,1,2,2),FIndex(1,1,1,2,1),FIndex(1,1,2,1,2),FIndex(1,1,2,1,1))),
                        (-1.0,(FIndex(1,1,2,2,2),FIndex(1,1,2,2,1),FIndex(1,1,1,1,2),FIndex(1,1,1,1,1))),
                        (-1.0,(FIndex(1,1,2,2,2),FIndex(1,1,2,2,1),FIndex(1,1,2,1,2),FIndex(1,1,2,1,1)))
    ]
end

@testset "fockcouplingnambus" begin
    @test fockcouplingnambus(nothing,('*','*','*','*'),(1,1,2,2))==(0,0,2,1)
    @test fockcouplingnambus(nothing,(0,0,1,1),(1,1,2,2))==(0,0,1,1)
    @test fockcouplingnambus(Val(:Pairing),('*','*'),(2,2))==(1,1)
    @test fockcouplingnambus(Val(:Pairing),(1,2),(2,2))==(1,2)
end

@testset "σ⁰" begin
    @test σ⁰("sp")==FockCoupling{2}(1,spins=(1,1))+FockCoupling{2}(1,spins=(2,2))
    @test σ⁰("ob")==FockCoupling{2}(1,orbitals=(1,1))+FockCoupling{2}(1,orbitals=(2,2))
    @test σ⁰("sl")==FockCoupling{2}(1,atoms=(1,1))+FockCoupling{2}(1,atoms=(2,2))
    @test σ⁰("ph")==FockCoupling{2}(1,nambus=(1,2))+FockCoupling{2}(1,nambus=(2,1))
end

@testset "σˣ" begin
    @test σˣ("sp")==FockCoupling{2}(1,spins=(1,2))+FockCoupling{2}(1,spins=(2,1))
    @test σˣ("ob")==FockCoupling{2}(1,orbitals=(1,2))+FockCoupling{2}(1,orbitals=(2,1))
    @test σˣ("sl")==FockCoupling{2}(1,atoms=(1,2))+FockCoupling{2}(1,atoms=(2,1))
    @test σˣ("ph")==FockCoupling{2}(1,nambus=(1,1))+FockCoupling{2}(1,nambus=(2,2))
end

@testset "σʸ" begin
    @test σʸ("sp")==FockCoupling{2}(1im,spins=(1,2))+FockCoupling{2}(-1im,spins=(2,1))
    @test σʸ("ob")==FockCoupling{2}(1im,orbitals=(1,2))+FockCoupling{2}(-1im,orbitals=(2,1))
    @test σʸ("sl")==FockCoupling{2}(1im,atoms=(1,2))+FockCoupling{2}(-1im,atoms=(2,1))
    @test σʸ("ph")==FockCoupling{2}(1im,nambus=(1,1))+FockCoupling{2}(-1im,nambus=(2,2))
end

@testset "σᶻ" begin
    @test σᶻ("sp")==FockCoupling{2}(-1,spins=(1,1))+FockCoupling{2}(1,spins=(2,2))
    @test σᶻ("ob")==FockCoupling{2}(-1,orbitals=(1,1))+FockCoupling{2}(1,orbitals=(2,2))
    @test σᶻ("sl")==FockCoupling{2}(-1,atoms=(1,1))+FockCoupling{2}(1,atoms=(2,2))
    @test σᶻ("ph")==FockCoupling{2}(-1,nambus=(1,2))+FockCoupling{2}(1,nambus=(2,1))
end

@testset "σ⁺" begin
    @test σ⁺("sp")==Couplings(FockCoupling{2}(1,spins=(2,1)))
    @test σ⁺("ob")==Couplings(FockCoupling{2}(1,orbitals=(2,1)))
    @test σ⁺("sl")==Couplings(FockCoupling{2}(1,atoms=(2,1)))
    @test σ⁺("ph")==Couplings(FockCoupling{2}(1,nambus=(2,2)))
end

@testset "σ⁻" begin
    @test σ⁻("sp")==Couplings(FockCoupling{2}(1,spins=(1,2)))
    @test σ⁻("ob")==Couplings(FockCoupling{2}(1,orbitals=(1,2)))
    @test σ⁻("sl")==Couplings(FockCoupling{2}(1,atoms=(1,2)))
    @test σ⁻("ph")==Couplings(FockCoupling{2}(1,nambus=(1,1)))
end

@testset "Onsite" begin
    term=Onsite{'F'}(:mu,1.5)
    @test term|>abbr==:st
    @test otype(term|>typeof,OID{FIndex{Int},SVector{2,Float},SVector{2,Float},Nothing})==FOperator{2,Float,ID{NTuple{2,OID{FIndex{Int},SVector{2,Float},SVector{2,Float},Nothing}}}}
    @test otype(term|>typeof,OID{FIndex{Int},SVector{2,Float},SVector{2,Float},Int})==FOperator{2,Float,ID{NTuple{2,OID{FIndex{Int},SVector{2,Float},SVector{2,Float},Int}}}}

    term=Onsite{'B'}(:mu,1.5,couplings=σˣ("sp")⊗σᶻ("ob"),modulate=true)
    @test term|>abbr==:st
    @test otype(term|>typeof,OID{FIndex{Int},SVector{2,Float},SVector{2,Float},Nothing})==BOperator{2,Float,ID{NTuple{2,OID{FIndex{Int},SVector{2,Float},SVector{2,Float},Nothing}}}}
    @test otype(term|>typeof,OID{FIndex{Int},SVector{2,Float},SVector{2,Float},Int})==BOperator{2,Float,ID{NTuple{2,OID{FIndex{Int},SVector{2,Float},SVector{2,Float},Int}}}}

    point=Point(PID('a',1),(0.5,0.5),(0.0,0.0))
    config=IDFConfig{Fock}(pid->Fock(atom=pid.site%2,norbital=2,nspin=2,nnambu=2),[point.pid])
    table=Table(config,by=usualfockindextotuple)

    term=Onsite{'F'}(:mu,1.5,couplings=σˣ("sp")⊗σᶻ("ob"),modulate=true)
    operators=Operators(FOperator(+1.5,ID(OID(FIndex('a',1,2,2,2),[0.5,0.5],[0.0,0.0],4),OID(FIndex('a',1,2,1,1),[0.5,0.5],[0.0,0.0],3))),
                        FOperator(-1.5,ID(OID(FIndex('a',1,1,2,2),[0.5,0.5],[0.0,0.0],2),OID(FIndex('a',1,1,1,1),[0.5,0.5],[0.0,0.0],1)))
    )
    @test expand(term,point,config,table,true)==operators
    @test expand(term,point,config,table,false)==operators+operators'
    operators=Operators(FOperator(+1.5,ID(OID(FIndex('a',1,2,2,2),[0.5,0.5],[0.0,0.0],nothing),OID(FIndex('a',1,2,1,1),[0.5,0.5],[0.0,0.0],nothing))),
                        FOperator(-1.5,ID(OID(FIndex('a',1,1,2,2),[0.5,0.5],[0.0,0.0],nothing),OID(FIndex('a',1,1,1,1),[0.5,0.5],[0.0,0.0],nothing)))
    )
    @test expand(term,point,config,nothing,true)==operators
    @test expand(term,point,config,nothing,false)==operators+operators'

    term=Onsite{'F'}(:mu,1.5,couplings=σᶻ("sp")⊗σᶻ("ob"),modulate=true)
    operators=Operators(FOperator(-0.75,ID(OID(FIndex('a',1,2,1,2),[0.5,0.5],[0.0,0.0],3),OID(FIndex('a',1,2,1,1),[0.5,0.5],[0.0,0.0],3))),
                        FOperator(-0.75,ID(OID(FIndex('a',1,1,2,2),[0.5,0.5],[0.0,0.0],2),OID(FIndex('a',1,1,2,1),[0.5,0.5],[0.0,0.0],2))),
                        FOperator(+0.75,ID(OID(FIndex('a',1,2,2,2),[0.5,0.5],[0.0,0.0],4),OID(FIndex('a',1,2,2,1),[0.5,0.5],[0.0,0.0],4))),
                        FOperator(+0.75,ID(OID(FIndex('a',1,1,1,2),[0.5,0.5],[0.0,0.0],1),OID(FIndex('a',1,1,1,1),[0.5,0.5],[0.0,0.0],1)))
    )
    @test expand(term,point,config,table,true)==operators
    @test expand(term,point,config,table,false)==operators+operators'
end

@testset "Hopping" begin
    bond=Bond(1,Point(PID('a',1),(0.5,0.5),(0.0,0.0)),Point(PID('b',2),(0.0,0.0),(0.0,0.0)))
    config=IDFConfig{Fock}(pid->Fock(atom=pid.site%2,norbital=2,nspin=2,nnambu=2),[bond.spoint.pid,bond.epoint.pid])
    table=Table(config,by=usualfockindextotuple)
    term=Hopping{'F'}(:t,1.5,1)
    operators=Operators(FOperator(1.5,ID(OID(FIndex('b',2,2,2,2),[0.0,0.0],[0.0,0.0],8),OID(FIndex('a',1,2,2,1),[0.5,0.5],[0.0,0.0],4))),
                        FOperator(1.5,ID(OID(FIndex('b',2,2,1,2),[0.0,0.0],[0.0,0.0],7),OID(FIndex('a',1,2,1,1),[0.5,0.5],[0.0,0.0],3))),
                        FOperator(1.5,ID(OID(FIndex('b',2,1,1,2),[0.0,0.0],[0.0,0.0],5),OID(FIndex('a',1,1,1,1),[0.5,0.5],[0.0,0.0],1))),
                        FOperator(1.5,ID(OID(FIndex('b',2,1,2,2),[0.0,0.0],[0.0,0.0],6),OID(FIndex('a',1,1,2,1),[0.5,0.5],[0.0,0.0],2)))
    )
    @test term|>abbr==:hp
    @test expand(term,bond,config,table,true)==operators
    @test expand(term,bond,config,table,false)==operators+operators'
end

@testset "Pairing" begin
    bond=Bond(1,Point(PID('a',1),(0.5,0.5),(0.0,0.0)),Point(PID('b',2),(0.0,0.0),(0.0,0.0)))
    config=IDFConfig{Fock}(pid->Fock(atom=pid.site%2,norbital=1,nspin=2,nnambu=2),[bond.spoint.pid,bond.epoint.pid])
    table=Table(config,by=nambufockindextotuple)
    term=Pairing{'F'}(:Δ,1.5,1,couplings=FockCoupling{2}(spins=(2,2)),amplitude=bond->bond|>rcoord|>azimuthd≈45 ? 1 : -1)
    operators=Operators(FOperator(-1.5,ID(OID(FIndex('b',2,1,2,1),[0.0,0.0],[0.0,0.0],6),OID(FIndex('a',1,1,2,1),[0.5,0.5],[0.0,0.0],2))),
                        FOperator(+1.5,ID(OID(FIndex('a',1,1,2,1),[0.5,0.5],[0.0,0.0],2),OID(FIndex('b',2,1,2,1),[0.0,0.0],[0.0,0.0],6)))
    )
    @test term|>abbr==:pr
    @test expand(term,bond,config,table,true)==operators
    @test expand(term,bond,config,table,false)==operators+operators'

    point=Point(PID('a',1),(0.5,0.5),(0.0,0.0))
    config=IDFConfig{Fock}(pid->Fock(atom=pid.site%2,norbital=1,nspin=2,nnambu=2),[point.pid])
    table=Table(config,by=nambufockindextotuple)
    term=Pairing{'F'}(:Δ,1.5,0,couplings=FockCoupling{2}(spins=(2,1))-FockCoupling{2}(spins=(1,2)))
    operators=Operators(FOperator(+1.5,ID(OID(FIndex('a',1,1,2,1),[0.5,0.5],[0.0,0.0],2),OID(FIndex('a',1,1,1,1),[0.5,0.5],[0.0,0.0],1))),
                        FOperator(-1.5,ID(OID(FIndex('a',1,1,1,1),[0.5,0.5],[0.0,0.0],1),OID(FIndex('a',1,1,2,1),[0.5,0.5],[0.0,0.0],2)))
    )
    @test term|>abbr==:pr
    @test expand(term,point,config,table,true)==operators
    @test expand(term,point,config,table,false)==operators+operators'
end

@testset "Hubbard" begin
    point=Point(PID('a',1),(0.5,0.5),(0.0,0.0))
    config=IDFConfig{Fock}(pid->Fock(atom=pid.site%2,norbital=2,nspin=2,nnambu=2),[point.pid])
    table=Table(config,by=usualfockindextotuple)
    term=Hubbard{'F'}(:H,2.5)
    operators=Operators(FOperator(1.25,ID(  OID(FIndex('a',1,1,2,2),[0.5,0.5],[0.0,0.0],2),OID(FIndex('a',1,1,2,1),[0.5,0.5],[0.0,0.0],2),
                                            OID(FIndex('a',1,1,1,2),[0.5,0.5],[0.0,0.0],1),OID(FIndex('a',1,1,1,1),[0.5,0.5],[0.0,0.0],1))),
                        FOperator(1.25,ID(  OID(FIndex('a',1,2,2,2),[0.5,0.5],[0.0,0.0],4),OID(FIndex('a',1,2,2,1),[0.5,0.5],[0.0,0.0],4),
                                            OID(FIndex('a',1,2,1,2),[0.5,0.5],[0.0,0.0],3),OID(FIndex('a',1,2,1,1),[0.5,0.5],[0.0,0.0],3)))
    )
    @test term|>abbr==:hb
    @test expand(term,point,config,table,true)==operators
    @test expand(term,point,config,table,false)==operators*2
end

@testset "InterOrbitalInterSpin" begin
    point=Point(PID('a',1),(0.5,0.5),(0.0,0.0))
    config=IDFConfig{Fock}(pid->Fock(atom=pid.site%2,norbital=2,nspin=2,nnambu=2),[point.pid])
    table=Table(config,by=usualfockindextotuple)
    term=InterOrbitalInterSpin{'F'}(:H,2.5)
    operators=Operators(FOperator(1.25,ID(  OID(FIndex('a',1,1,2,2),[0.5,0.5],[0.0,0.0],2),OID(FIndex('a',1,1,2,1),[0.5,0.5],[0.0,0.0],2),
                                            OID(FIndex('a',1,2,1,2),[0.5,0.5],[0.0,0.0],3),OID(FIndex('a',1,2,1,1),[0.5,0.5],[0.0,0.0],3))),
                        FOperator(1.25,ID(  OID(FIndex('a',1,1,1,2),[0.5,0.5],[0.0,0.0],1),OID(FIndex('a',1,1,1,1),[0.5,0.5],[0.0,0.0],1),
                                            OID(FIndex('a',1,2,2,2),[0.5,0.5],[0.0,0.0],4),OID(FIndex('a',1,2,2,1),[0.5,0.5],[0.0,0.0],4)))
    )
    @test term|>abbr==:nons
    @test expand(term,point,config,table,true)==operators
    @test expand(term,point,config,table,false)==operators*2
end

@testset "InterOrbitalIntraSpin" begin
    point=Point(PID('a',1),(0.5,0.5),(0.0,0.0))
    config=IDFConfig{Fock}(pid->Fock(atom=pid.site%2,norbital=2,nspin=2,nnambu=2),[point.pid])
    table=Table(config,by=usualfockindextotuple)
    term=InterOrbitalIntraSpin{'F'}(:H,2.5)
    operators=Operators(FOperator(1.25,ID(  OID(FIndex('a',1,1,1,2),[0.5,0.5],[0.0,0.0],1),OID(FIndex('a',1,1,1,1),[0.5,0.5],[0.0,0.0],1),
                                            OID(FIndex('a',1,2,1,2),[0.5,0.5],[0.0,0.0],3),OID(FIndex('a',1,2,1,1),[0.5,0.5],[0.0,0.0],3))),
                        FOperator(1.25,ID(  OID(FIndex('a',1,1,2,2),[0.5,0.5],[0.0,0.0],2),OID(FIndex('a',1,1,2,1),[0.5,0.5],[0.0,0.0],2),
                                            OID(FIndex('a',1,2,2,2),[0.5,0.5],[0.0,0.0],4),OID(FIndex('a',1,2,2,1),[0.5,0.5],[0.0,0.0],4)))
    )
    @test term|>abbr==:noes
    @test expand(term,point,config,table,true)==operators
    @test expand(term,point,config,table,false)==operators*2
end

@testset "SpinFlip" begin
    point=Point(PID('a',1),(0.5,0.5),(0.0,0.0))
    config=IDFConfig{Fock}(pid->Fock(atom=pid.site%2,norbital=2,nspin=2,nnambu=2),[point.pid])
    table=Table(config,by=usualfockindextotuple)
    term=SpinFlip{'F'}(:H,2.5)
    operators=Operators(FOperator(2.5,ID(   OID(FIndex('a',1,1,2,2),[0.5,0.5],[0.0,0.0],2),OID(FIndex('a',1,2,1,2),[0.5,0.5],[0.0,0.0],3),
                                            OID(FIndex('a',1,1,1,1),[0.5,0.5],[0.0,0.0],1),OID(FIndex('a',1,2,2,1),[0.5,0.5],[0.0,0.0],4)))
    )
    @test term|>abbr==:sf
    @test expand(term,point,config,table,true)==operators
    @test expand(term,point,config,table,false)==operators+operators'
end

@testset "PairHopping" begin
    point=Point(PID('a',1),(0.5,0.5),(0.0,0.0))
    config=IDFConfig{Fock}(pid->Fock(atom=pid.site%2,norbital=2,nspin=2,nnambu=2),[point.pid])
    table=Table(config,by=usualfockindextotuple)
    term=PairHopping{'F'}(:H,2.5)
    operators=Operators(FOperator(2.5,ID(   OID(FIndex('a',1,1,2,2),[0.5,0.5],[0.0,0.0],2),OID(FIndex('a',1,1,1,2),[0.5,0.5],[0.0,0.0],1),
                                            OID(FIndex('a',1,2,1,1),[0.5,0.5],[0.0,0.0],3),OID(FIndex('a',1,2,2,1),[0.5,0.5],[0.0,0.0],4)))
    )
    @test term|>abbr==:ph
    @test expand(term,point,config,table,true)==operators
    @test expand(term,point,config,table,false)==operators+operators'
end

@testset "Coulomb" begin
    bond=Bond(1,Point(PID('a',1),(0.5,0.5),(0.0,0.0)),Point(PID('b',2),(0.0,0.0),(0.0,0.0)))
    config=IDFConfig{Fock}(pid->Fock(atom=pid.site%2,norbital=1,nspin=2,nnambu=2),[bond.spoint.pid,bond.epoint.pid])
    table=Table(config,by=usualfockindextotuple)

    term=Coulomb{'F'}(:V,2.5,1,couplings=σᶻ("sp")*σᶻ("sp"))
    operators=Operators(FOperator(-1.25,ID( OID(FIndex('b',2,1,1,2),[0.0,0.0],[0.0,0.0],3),OID(FIndex('b',2,1,1,1),[0.0,0.0],[0.0,0.0],3),
                                            OID(FIndex('a',1,1,2,2),[0.5,0.5],[0.0,0.0],2),OID(FIndex('a',1,1,2,1),[0.5,0.5],[0.0,0.0],2))),
                        FOperator(+1.25,ID( OID(FIndex('b',2,1,1,2),[0.0,0.0],[0.0,0.0],3),OID(FIndex('b',2,1,1,1),[0.0,0.0],[0.0,0.0],3),
                                            OID(FIndex('a',1,1,1,2),[0.5,0.5],[0.0,0.0],1),OID(FIndex('a',1,1,1,1),[0.5,0.5],[0.0,0.0],1))),
                        FOperator(-1.25,ID( OID(FIndex('b',2,1,2,2),[0.0,0.0],[0.0,0.0],4),OID(FIndex('b',2,1,2,1),[0.0,0.0],[0.0,0.0],4),
                                            OID(FIndex('a',1,1,1,2),[0.5,0.5],[0.0,0.0],1),OID(FIndex('a',1,1,1,1),[0.5,0.5],[0.0,0.0],1))),
                        FOperator(+1.25,ID( OID(FIndex('b',2,1,2,2),[0.0,0.0],[0.0,0.0],4),OID(FIndex('b',2,1,2,1),[0.0,0.0],[0.0,0.0],4),
                                            OID(FIndex('a',1,1,2,2),[0.5,0.5],[0.0,0.0],2),OID(FIndex('a',1,1,2,1),[0.5,0.5],[0.0,0.0],2)))
    )
    @test term|>abbr==:cl
    @test expand(term,bond,config,table,true)==operators
    @test expand(term,bond,config,table,false)==operators+operators'

    term=Coulomb{'F'}(:V,2.5,1,couplings=σˣ("sp")*σᶻ("sp"))
    operators=Operators(FOperator(-2.5,ID(  OID(FIndex('b',2,1,2,2),[0.0,0.0],[0.0,0.0],4),OID(FIndex('b',2,1,1,1),[0.0,0.0],[0.0,0.0],3),
                                            OID(FIndex('a',1,1,1,2),[0.5,0.5],[0.0,0.0],1),OID(FIndex('a',1,1,1,1),[0.5,0.5],[0.0,0.0],1))),
                        FOperator(+2.5,ID(  OID(FIndex('b',2,1,1,2),[0.0,0.0],[0.0,0.0],3),OID(FIndex('b',2,1,2,1),[0.0,0.0],[0.0,0.0],4),
                                            OID(FIndex('a',1,1,2,2),[0.5,0.5],[0.0,0.0],2),OID(FIndex('a',1,1,2,1),[0.5,0.5],[0.0,0.0],2))),
                        FOperator(+2.5,ID(  OID(FIndex('b',2,1,2,2),[0.0,0.0],[0.0,0.0],4),OID(FIndex('b',2,1,1,1),[0.0,0.0],[0.0,0.0],3),
                                            OID(FIndex('a',1,1,2,2),[0.5,0.5],[0.0,0.0],2),OID(FIndex('a',1,1,2,1),[0.5,0.5],[0.0,0.0],2))),
                        FOperator(-2.5,ID(  OID(FIndex('b',2,1,1,2),[0.0,0.0],[0.0,0.0],3),OID(FIndex('b',2,1,2,1),[0.0,0.0],[0.0,0.0],4),
                                            OID(FIndex('a',1,1,1,2),[0.5,0.5],[0.0,0.0],1),OID(FIndex('a',1,1,1,1),[0.5,0.5],[0.0,0.0],1)))
    )
    @test term|>abbr==:cl
    @test expand(term,bond,config,table,true)==operators
    @test expand(term,bond,config,table,false)==operators+operators'
end
