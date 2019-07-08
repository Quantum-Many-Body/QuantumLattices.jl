using Test
using StaticArrays: SVector
using QuantumLattices.Essentials.FockPackage
using QuantumLattices.Essentials.Spatials: Bond,Point,PID,rcoord,azimuthd
using QuantumLattices.Essentials.DegreesOfFreedom: Table,IDFConfig,OID,Operators,oidtype,otype,script,latexformat
using QuantumLattices.Essentials.Terms: Couplings,@subscript,statistics,abbr
using QuantumLattices.Interfaces: dims,inds,⊗,⋅,expand,permute
using QuantumLattices.Prerequisites: Float
using QuantumLattices.Mathematics.VectorSpaces: IsMultiIndexable,MultiIndexOrderStyle
using QuantumLattices.Mathematics.AlgebraOverFields: ID
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
    @test angle(OID(FIndex(1,1,1,1,1),[0.0,0.0],[1.0,2.0],1),[[1.0,0.0],[0.0,1.0]],[0.1,0.0])≈2pi*0.1
    @test angle(OID(FIndex(1,1,1,1,2),[0.0,0.0],[1.0,2.0],1),[[1.0,0.0],[0.0,1.0]],[0.0,0.2])≈-2pi*0.4
end

@testset "script" begin
    @test script(OID(FIndex('c',1,2,1,1)),Val(:site))==1
    @test script(OID(FIndex('c',1,2,1,1)),Val(:orbital))==2
    @test script(OID(FIndex('c',1,2,3,1)),Val(:spinint))==3
    @test script(OID(FIndex('c',1,2,2,1)),Val(:spinsym))=="↑"
    @test script(OID(FIndex('c',1,2,1,1)),Val(:spinsym))=="↓"
    @test script(OID(FIndex('c',1,2,3,1)),Val(:nambu))==""
    @test script(OID(FIndex('c',1,2,3,2)),Val(:nambu))=="\\dagger"
end

@testset "oidtype" begin
    @test oidtype(FID,Point{2,PID{Int}},Nothing,Val(true))==OID{FIndex{Int},SVector{2,Float},SVector{2,Float},Nothing}
    @test oidtype(FID,Point{2,PID{Int}},Table,Val(true))==OID{FIndex{Int},SVector{2,Float},SVector{2,Float},Int}
    @test oidtype(FID,Point{2,PID{Int}},Nothing,Val(false))==OID{FIndex{Int},Nothing,Nothing,Nothing}
    @test oidtype(FID,Point{2,PID{Int}},Table,Val(false))==OID{FIndex{Int},Nothing,Nothing,Int}
end

@testset "FockOperator" begin
    @test latexformat(FOperator)==foptdefaultlatex

    opt=FOperator(1.0,(FIndex(1,1,1,1,2),FIndex(1,1,1,1,1)))
    @test opt|>isnormalordered
    @test opt|>statistics==opt|>typeof|>statistics=='F'

    opt=FOperator(1.0,(FIndex(1,2,1,1,2),FIndex(1,2,1,1,1),FIndex(1,1,1,2,2),FIndex(1,1,1,2,1)))
    @test opt|>isnormalordered==false
    @test repr(opt)=="c^{\\dagger}_{2,1,↓}c^{}_{2,1,↓}c^{\\dagger}_{1,1,↑}c^{}_{1,1,↑}"

    opt1=FOperator(1.5,(FIndex(1,2,1,1,2),FIndex(1,2,1,1,1)))
    opt2=FOperator(2.0,(FIndex(1,2,1,1,1),FIndex(1,2,1,1,2)))
    @test opt1*opt2==nothing

    opt1=FOperator(1.5,(FIndex(1,2,1,1,2),FIndex(1,2,1,1,1)))
    opt2=FOperator(2.0,(FIndex(1,2,1,1,2),FIndex(1,2,1,1,1)))
    @test opt1*opt2==FOperator(3.0,(FIndex(1,2,1,1,2),FIndex(1,2,1,1,1),FIndex(1,2,1,1,2),FIndex(1,2,1,1,1)))

    @test latexformat(BOperator)==boptdefaultlatex

    opt=BOperator(1.0,(FIndex(1,1,1,1,2),FIndex(1,1,1,1,1)))
    @test opt|>statistics==opt|>typeof|>statistics=='B'
    @test repr(opt)=="b^{\\dagger}_{1,1,↓}b^{}_{1,1,↓}"
end

@testset "permute" begin
    i1,i2=OID(FIndex(1,1,1,1,2)),OID(FIndex(1,1,1,1,1))
    @test permute(FOperator,i1,i2)==(FOperator(1),FOperator(-1,ID(i2,i1)))
    @test permute(FOperator,i2,i1)==(FOperator(1),FOperator(-1,ID(i1,i2)))
    @test permute(BOperator,i1,i2)==(BOperator(1),BOperator(1,ID(i2,i1)))
    @test permute(BOperator,i2,i1)==(BOperator(-1),BOperator(1,ID(i1,i2)))

    i1,i2=OID(FIndex(1,1,1,2,2)),OID(FIndex(1,1,1,1,1))
    @test permute(FOperator,i1,i2)==(FOperator(-1,ID(i2,i1)),)
    @test permute(FOperator,i2,i1)==(FOperator(-1,ID(i1,i2)),)
    @test permute(BOperator,i1,i2)==(BOperator(1,ID(i2,i1)),)
    @test permute(BOperator,i2,i1)==(BOperator(1,ID(i1,i2)),)
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

@testset "fockcoupling" begin
    fc=fc"1.0 sl(1,1,1,1)⊗ob(α,α,β,β)⊗sp(σ,γ,σ,γ)⊗ph(2,1,2,1)@(1,1,2,2) with (α<β,σ≠γ)"
    ob,sp=fc.obsubscripts[1].identifier,fc.spsubscripts[1].identifier
    @test repr(fc)=="1.0 sl(1,1,1,1)⊗ob(α,α,β,β)⊗sp(σ,γ,σ,γ)⊗ph(2,1,2,1)@(1,1,2,2) with ($ob,$sp) && ($ob,$sp) && ($ob,$sp) && ($ob,$sp)"

    fc=fc"1.0 sl(1,1,1,1)⊗ob(α,α,β,β)⊗sp(σ,γ,σ,γ)⊗ph(2,1,2,1)@(1,1,2,2)"
    ob,sp=fc.obsubscripts[1].identifier,fc.spsubscripts[1].identifier
    @test repr(fc)=="1.0 sl(1,1,1,1)⊗ob(α,α,β,β)⊗sp(σ,γ,σ,γ)⊗ph(2,1,2,1)@(1,1,2,2) with ($ob,$sp) && ($ob,$sp) && ($ob,$sp) && ($ob,$sp)"

    fc=fc"1.0 sl(1,1,1,1)⊗ob(α,α,β,β)⊗sp(σ,γ,σ,γ)⊗ph(2,1,2,1)"
    ob,sp=fc.obsubscripts[1].identifier,fc.spsubscripts[1].identifier
    @test repr(fc)=="1.0 sl(1,1,1,1)⊗ob(α,α,β,β)⊗sp(σ,γ,σ,γ)⊗ph(2,1,2,1) with ($ob,$sp) && ($ob,$sp) && ($ob,$sp) && ($ob,$sp)"

    fc=fc"1.0 sl(1,1,1,1)⊗ob(α,α,β,β)⊗sp(σ,γ,σ,γ)⊗ph(2,1,2,1) with (α<β,σ≠γ)"
    ob,sp=fc.obsubscripts[1].identifier,fc.spsubscripts[1].identifier
    @test repr(fc)=="1.0 sl(1,1,1,1)⊗ob(α,α,β,β)⊗sp(σ,γ,σ,γ)⊗ph(2,1,2,1) with ($ob,$sp) && ($ob,$sp) && ($ob,$sp) && ($ob,$sp)"

    fc=fc"1.0 sl(1,1,1,1)⊗ob(α,α,β,β)⊗sp(σ,γ,σ,γ)⊗ph(2,1,2,1) with α<β,σ≠γ"
    ob,sp=fc.obsubscripts[1].identifier,fc.spsubscripts[1].identifier
    @test repr(fc)=="1.0 sl(1,1,1,1)⊗ob(α,α,β,β)⊗sp(σ,γ,σ,γ)⊗ph(2,1,2,1) with ($ob,$sp) && ($ob,$sp) && ($ob,$sp) && ($ob,$sp)"

    fc=fc"1.0 sl(1,1,1,1)⊗ob(α,α,β,β)⊗sp(σ,γ,σ,γ)⊗ph(2,1,2,1) with (α<β,*)"
    ob,sp=fc.obsubscripts[1].identifier,fc.spsubscripts[1].identifier
    @test repr(fc)=="1.0 sl(1,1,1,1)⊗ob(α,α,β,β)⊗sp(σ,γ,σ,γ)⊗ph(2,1,2,1) with ($ob,$sp) && ($ob,$sp) && ($ob,$sp) && ($ob,$sp)"

    fc=fc"1.0 sl(1,1,1,1)⊗ob(α,α,β,β)⊗ph(2,1,2,1) with α<β"
    ob,sp=fc.obsubscripts[1].identifier,fc.spsubscripts[1].identifier
    @test repr(fc)=="1.0 sl(1,1,1,1)⊗ob(α,α,β,β)⊗ph(2,1,2,1) with ($ob,$sp) && ($ob,$sp) && ($ob,$sp) && ($ob,$sp)"

    fc=fc"1.0 sl(1,1,1,1)⊗ob(α,α,β,β)⊗sp(2,1,2,1) with α<β"
    ob,sp=fc.obsubscripts[1].identifier,fc.spsubscripts[1].identifier
    @test repr(fc)=="1.0 sl(1,1,1,1)⊗ob(α,α,β,β)⊗sp(2,1,2,1) with ($ob,$sp) && ($ob,$sp) && ($ob,$sp) && ($ob,$sp)"

    fc=fc"1.0 sl(1,1,1,1)⊗ob(1,1,1,1)⊗ph(2,1,2,1)"
    @test repr(fc)=="1.0 sl(1,1,1,1)⊗ob(1,1,1,1)⊗ph(2,1,2,1)"

    fc=fc"1.0 sl(1,1,1,1)⊗ph(2,1,2,1)"
    @test repr(fc)=="1.0 sl(1,1,1,1)⊗ph(2,1,2,1)"

    fc=fc"1.0 @(1,1,1,1)"
    @test repr(fc)=="1.0 @(1,1,1,1)"

    fc=fc"1.0im {2}"
    @test repr(fc)=="1.0im {2}"
end

@testset "fockcouplings" begin
    @test σ⁰"sp"==σ⁰("sp") && σ⁰"sp@(1,2)"==σ⁰("sp",centers=(1,2))
    @test σ⁰"ob"==σ⁰("ob") && σ⁰"ob@(1,2)"==σ⁰("ob",centers=(1,2))
    @test σ⁰"sl"==σ⁰("sl") && σ⁰"sl@(1,2)"==σ⁰("sl",centers=(1,2))
    @test σ⁰"ph"==σ⁰("ph") && σ⁰"ph@(1,2)"==σ⁰("ph",centers=(1,2))

    @test σˣ"sp"==σˣ("sp") && σˣ"sp@(1,2)"==σˣ("sp",centers=(1,2))
    @test σˣ"ob"==σˣ("ob") && σˣ"ob@(1,2)"==σˣ("ob",centers=(1,2))
    @test σˣ"sl"==σˣ("sl") && σˣ"sl@(1,2)"==σˣ("sl",centers=(1,2))
    @test σˣ"ph"==σˣ("ph") && σˣ"ph@(1,2)"==σˣ("ph",centers=(1,2))

    @test σʸ"sp"==σʸ("sp") && σʸ"sp@(1,2)"==σʸ("sp",centers=(1,2))
    @test σʸ"ob"==σʸ("ob") && σʸ"ob@(1,2)"==σʸ("ob",centers=(1,2))
    @test σʸ"sl"==σʸ("sl") && σʸ"sl@(1,2)"==σʸ("sl",centers=(1,2))
    @test σʸ"ph"==σʸ("ph") && σʸ"ph@(1,2)"==σʸ("ph",centers=(1,2))

    @test σᶻ"sp"==σᶻ("sp") && σᶻ"sp@(1,2)"==σᶻ("sp",centers=(1,2))
    @test σᶻ"ob"==σᶻ("ob") && σᶻ"ob@(1,2)"==σᶻ("ob",centers=(1,2))
    @test σᶻ"sl"==σᶻ("sl") && σᶻ"sl@(1,2)"==σᶻ("sl",centers=(1,2))
    @test σᶻ"ph"==σᶻ("ph") && σᶻ"ph@(1,2)"==σᶻ("ph",centers=(1,2))

    @test σ⁺"sp"==σ⁺("sp") && σ⁺"sp@(1,2)"==σ⁺("sp",centers=(1,2))
    @test σ⁺"ob"==σ⁺("ob") && σ⁺"ob@(1,2)"==σ⁺("ob",centers=(1,2))
    @test σ⁺"sl"==σ⁺("sl") && σ⁺"sl@(1,2)"==σ⁺("sl",centers=(1,2))
    @test σ⁺"ph"==σ⁺("ph") && σ⁺"ph@(1,2)"==σ⁺("ph",centers=(1,2))

    @test σ⁻"sp"==σ⁻("sp") && σ⁻"sp@(1,2)"==σ⁻("sp",centers=(1,2))
    @test σ⁻"ob"==σ⁻("ob") && σ⁻"ob@(1,2)"==σ⁻("ob",centers=(1,2))
    @test σ⁻"sl"==σ⁻("sl") && σ⁻"sl@(1,2)"==σ⁻("sl",centers=(1,2))
    @test σ⁻"ph"==σ⁻("ph") && σ⁻"ph@(1,2)"==σ⁻("ph",centers=(1,2))
end

@testset "Onsite" begin
    term=Onsite{'F'}(:mu,1.5)
    @test term|>abbr==:st
    @test otype(term|>typeof,OID{FIndex{Int},SVector{2,Float},SVector{2,Float},Nothing})==FOperator{Float,ID{OID{FIndex{Int},SVector{2,Float},SVector{2,Float},Nothing},2}}
    @test otype(term|>typeof,OID{FIndex{Int},SVector{2,Float},SVector{2,Float},Int})==FOperator{Float,ID{OID{FIndex{Int},SVector{2,Float},SVector{2,Float},Int},2}}

    term=Onsite{'B'}(:mu,1.5,couplings=σˣ("sp")⊗σᶻ("ob"),modulate=true)
    @test term|>abbr==:st
    @test otype(term|>typeof,OID{FIndex{Int},SVector{2,Float},SVector{2,Float},Nothing})==BOperator{Float,ID{OID{FIndex{Int},SVector{2,Float},SVector{2,Float},Nothing},2}}
    @test otype(term|>typeof,OID{FIndex{Int},SVector{2,Float},SVector{2,Float},Int})==BOperator{Float,ID{OID{FIndex{Int},SVector{2,Float},SVector{2,Float},Int},2}}

    point=Point(PID('a',1),(0.5,0.5),(0.0,0.0))
    config=IDFConfig{Fock}(pid->Fock(atom=pid.site%2,norbital=2,nspin=2,nnambu=2),[point.pid])
    table=Table(config,usualfockindextotuple)

    term=Onsite{'F'}(:mu,1.5,couplings=σˣ("sp")⊗σᶻ("ob"),modulate=true)
    operators=Operators(FOperator(+1.5,ID(OID(FIndex('a',1,2,2,2),[0.5,0.5],[0.0,0.0],4),OID(FIndex('a',1,2,1,1),[0.5,0.5],[0.0,0.0],3))),
                        FOperator(-1.5,ID(OID(FIndex('a',1,1,1,2),[0.5,0.5],[0.0,0.0],1),OID(FIndex('a',1,1,2,1),[0.5,0.5],[0.0,0.0],2)))
    )
    @test expand(term,point,config,table,true)==operators
    @test expand(term,point,config,table,false)==operators+operators'
    operators=Operators(FOperator(+1.5,ID(OID(FIndex('a',1,2,2,2),[0.5,0.5],[0.0,0.0],nothing),OID(FIndex('a',1,2,1,1),[0.5,0.5],[0.0,0.0],nothing))),
                        FOperator(-1.5,ID(OID(FIndex('a',1,1,1,2),[0.5,0.5],[0.0,0.0],nothing),OID(FIndex('a',1,1,2,1),[0.5,0.5],[0.0,0.0],nothing)))
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
    table=Table(config,usualfockindextotuple)
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
    table=Table(config,nambufockindextotuple)
    term=Pairing{'F'}(:Δ,1.5,1,couplings=FockCoupling{2}(spins=(2,2)),amplitude=bond->bond|>rcoord|>azimuthd≈45 ? 1 : -1)
    operators=Operators(FOperator(-1.5,ID(OID(FIndex('b',2,1,2,1),[0.0,0.0],[0.0,0.0],6),OID(FIndex('a',1,1,2,1),[0.5,0.5],[0.0,0.0],2))),
                        FOperator(+1.5,ID(OID(FIndex('a',1,1,2,1),[0.5,0.5],[0.0,0.0],2),OID(FIndex('b',2,1,2,1),[0.0,0.0],[0.0,0.0],6)))
    )
    @test term|>abbr==:pr
    @test expand(term,bond,config,table,true)==operators
    @test expand(term,bond,config,table,false)==operators+operators'

    point=Point(PID('a',1),(0.5,0.5),(0.0,0.0))
    config=IDFConfig{Fock}(pid->Fock(atom=pid.site%2,norbital=1,nspin=2,nnambu=2),[point.pid])
    table=Table(config,nambufockindextotuple)
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
    table=Table(config,usualfockindextotuple)
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
    table=Table(config,usualfockindextotuple)
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
    table=Table(config,usualfockindextotuple)
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
    table=Table(config,usualfockindextotuple)
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
    table=Table(config,usualfockindextotuple)
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
    table=Table(config,usualfockindextotuple)

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
