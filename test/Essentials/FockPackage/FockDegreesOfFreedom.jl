using Hamiltonian.Essentials.FockPackage
using Hamiltonian.Prerequisites: Float
using Hamiltonian.Essentials.Spatials: PID,Point,Bond
using Hamiltonian.Essentials.DegreesOfFreedom: Couplings
using Hamiltonian.Mathematics.VectorSpaces: IsMultiIndexable,MultiIndexOrderStyle

@testset "FID" begin
    fid=FID(orbital=1,spin=1)
    @test fid|>typeof|>fieldnames==(:orbital,:spin,:nambu)
    @test fid'==FID(1,1,2)
    @test fid''==FID(1,1,1)
end

@testset "FIndex" begin
    @test FIndex|>fieldnames==(:scope,:site,:orbital,:spin,:nambu)
    @test union(PID{Int},FID)==FIndex{Int}
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

@testset "Subscript" begin
    sub=@subscript (x1,x2)=>(x1,4,4,x2) with x1<x2
    @test sub==deepcopy(sub)
    @test isequal(sub,deepcopy(sub))
    @test rank(sub)==rank(typeof(sub))==2
    @test dimension(sub)==dimension(typeof(sub))==4
    @test string(sub)=="(x1,x2)=>(x1,4,4,x2) with $(sub.identifier)"
    @test sub(Val('M'),1,2)==(1,4,4,2)
    @test sub(Val('C'),1,2)==true
    @test sub(Val('C'),2,1)==false

    sub=@subscript (x1,x2)=>(x1,x2,x1,x2)
    @test sub(Val('M'),1,2)==(1,2,1,2)
    @test sub(Val('C'),1,2)==true

    sub=Subscript{4}()
    @test rank(sub)==1 && dimension(sub)==4
    @test string(sub)=="*=>(*,*,*,*)"
    @test sub(Val('M'),2)==(2,2,2,2)
    @test sub(Val('C'),2)==true

    sub=Subscript((1,2,2,1))
    @test rank(sub)==0 && dimension(sub)==4
    @test string(sub)=="(1,2,2,1)"
    @test sub(Val('M'))==(1,2,2,1)
    @test sub(Val('C'))==true
end

@testset "Subscripts" begin
    sub1=@subscript (x1,)=>(x1,2) with x1<2
    sub2=@subscript (y1,y2)=>(y1,y2,y1,y2) with y1<y2
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

@testset "FCID" begin
    @test FCID(atom=1,orbital=1,spin=1,nambu=1)|>typeof|>fieldnames==(:center,:atom,:orbital,:spin,:nambu,:obsub,:spsub)
end

@testset "FockCoupling" begin
    @test FockCoupling{2}(1.0)|>string=="FockCoupling(value=1.0)"
    @test FockCoupling{2}(1.0,atoms=(1,1))|>string=="FockCoupling(value=1.0,atoms=(1,1))"
    @test FockCoupling{2}(1.0,atoms=(1,1),spins=(1,2))|>string=="FockCoupling(value=1.0,atoms=(1,1),spins=(1,2))"
    @test FockCoupling{2}(2.0)|>repr=="2.0"

    fc1=FockCoupling{2}(1.5,atoms=(2,1),spins=(@subscript (x,)=>(x,1)),centers=(1,2))
    fc2=FockCoupling{2}(2.0,atoms=(1,2),orbitals=(@subscript (x,y)=>(x,y) with x<y),centers=(1,2))
    @test fc1|>repr=="1.5 sl(2:1)⊗sp(x:1)@(1-2) with (*,$(fc1.id[1].spsub)) && (*,$(fc1.id[1].spsub))"
    @test fc2|>repr=="2.0 sl(1:2)⊗ob(x:y)@(1-2) with ($(fc2.id[2].obsub),*) && ($(fc2.id[2].obsub),*)"
    fc=fc1*fc2
    @test fc|>repr=="3.0 sl(2:1:1:2)⊗ob(*:*:x:y)⊗sp(x:1:*:*)@(1-2-1-2) with (*,$(fc1.id[1].spsub)) && (*,$(fc1.id[2].spsub)) && ($(fc2.id[2].obsub),*) && ($(fc2.id[2].obsub),*)"
    fc=fc1⊗fc2
    @test fc|>repr=="3.0 sl(2:2)⊗ob(x:y)⊗sp(x:1)@(1-2) with ($(fc2.id[2].obsub),$(fc1.id[1].spsub)) && ($(fc2.id[2].obsub),$(fc1.id[1].spsub))"

    ex=expand(FockCoupling{2}(2.0,atoms=(1,1)),PID(1,1),Fock(atom=2,norbital=2,nspin=2,nnambu=2))
    @test IsMultiIndexable(typeof(ex))==IsMultiIndexable(true)
    @test MultiIndexOrderStyle(typeof(ex))==MultiIndexOrderStyle('C')
    @test dims(ex)==(0,0)
    @test collect(ex)==[]

    ex=expand(FockCoupling{2}(2.0,atoms=(1,2),orbitals=(1,2)),(PID(1,1),PID(1,2)),(Fock(atom=1,norbital=2,nspin=2,nnambu=2),Fock(atom=2,norbital=2,nspin=2,nnambu=2)))
    @test dims(ex)==(1,2)
    @test collect(ex)==[(2.0,(FIndex(1,1,1,1,2),FIndex(1,2,2,1,1))),
                        (2.0,(FIndex(1,1,1,2,2),FIndex(1,2,2,2,1)))
    ]

    ex=expand(FockCoupling{4}(2.0,centers=(1,1,1,1),spins=(2,2,1,1)),PID(1,1),Fock(atom=1,norbital=2,nspin=2,nnambu=2))
    @test dims(ex)==(2,1)
    @test collect(ex)==[(2.0,(FIndex(1,1,1,2,2),FIndex(1,1,1,2,1),FIndex(1,1,1,1,2),FIndex(1,1,1,1,1))),
                        (2.0,(FIndex(1,1,2,2,2),FIndex(1,1,2,2,1),FIndex(1,1,2,1,2),FIndex(1,1,2,1,1)))
    ]

    ex=expand(FockCoupling{4}(2.0,orbitals=(@subscript (α,β)=>(α,α,β,β) with α<β),spins=(2,1,1,2),nambus=(2,2,1,1)),PID(1,1),Fock(atom=1,norbital=3,nspin=2,nnambu=2))
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
