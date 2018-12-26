using Hamiltonian.Utilities: Float
using Hamiltonian.Essentials.FockPackage
using Hamiltonian.Essentials.Spatial: PID,Point,Bond
using Hamiltonian.Essentials.DegreeOfFreedom: Couplings

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
    @test fock|>length==4
    @test fock|>collect==[FID(1,1,1),FID(1,1,2),FID(1,2,1),FID(1,2,2)]
end

@testset "FCID" begin
    @test FCID(atom=1,orbital=1,spin=1,nambu=1)|>typeof|>fieldnames==(:atom,:orbital,:spin,:nambu)
end

@testset "FockCoupling" begin
    @test FockCoupling{2}(1.0)|>string=="FockCoupling(value=1.0)"
    @test FockCoupling{2}(1.0,atoms=(1,1))|>string=="FockCoupling(value=1.0,atoms=(1,1))"
    @test FockCoupling{2}(1.0,atoms=(1,1),spins=(1,2))|>string=="FockCoupling(value=1.0,atoms=(1,1),spins=(1,2))"
    @test FockCoupling{2}(2.0)|>repr=="2.0"
    @test FockCoupling{2}(2.0,atoms=(1,1))|>repr=="2.0 sl11"
    @test FockCoupling{2}(2.0,atoms=(1,1),orbitals=(2,2),spins=(1,2),nambus=(1,2))|>repr=="2.0 sl11⊗ob22⊗sp12⊗ph12"

    fc1=FockCoupling{2}(1.5,atoms=(1,2),orbitals=(1,1),spins=(1,1))
    fc2=FockCoupling{2}(2.0,atoms=(2,1),orbitals=(1,2),nambus=(1,2))
    @test fc1*fc2==FockCoupling{2}(3.0,atoms=(1,1),orbitals=(1,2),spins=(1,1),nambus=(1,2))
    fc1=FockCoupling{2}(1.5,atoms=(1,2))
    fc2=FockCoupling{2}(2.0,atoms=(1,2))
    @test fc1*fc2==FockCoupling{2}(0.0,atoms=(1,2))

    point,fock=Point(PID(1,1),(0.0,0.0)),Fock(atom=2,norbital=2,nspin=2,nnambu=2)
    ex=expand(FockCoupling{2}(2.0,atoms=(1,1)),point,fock)
    @test ex|>typeof|>eltype==Tuple{Float,NTuple{2,FIndex{Int}}}
    @test ex|>length==0
    @test ex|>collect==[]

    ex=expand(FockCoupling{2}(2.0,atoms=(2,2)),point,fock)
    @test ex|>length==4
    @test ex|>collect==[(2.0,(FIndex(1,1,1,1,2),FIndex(1,1,1,1,1))),
                        (2.0,(FIndex(1,1,1,2,2),FIndex(1,1,1,2,1))),
                        (2.0,(FIndex(1,1,2,1,2),FIndex(1,1,2,1,1))),
                        (2.0,(FIndex(1,1,2,2,2),FIndex(1,1,2,2,1)))
                        ]

    ex=expand(FockCoupling{2}(2.0,orbitals=(1,2)),point,fock)
    @test ex|>length==2
    @test ex|>collect==[(2.0,(FIndex(1,1,1,1,2),FIndex(1,1,2,1,1))),
                        (2.0,(FIndex(1,1,1,2,2),FIndex(1,1,2,2,1)))
                        ]

    ex=expand(FockCoupling{2}(2.0,spins=(1,2)),point,fock)
    @test ex|>length==2
    @test ex|>collect==[(2.0,(FIndex(1,1,1,1,2),FIndex(1,1,1,2,1))),
                        (2.0,(FIndex(1,1,2,1,2),FIndex(1,1,2,2,1)))
                        ]

    ex=expand(FockCoupling{2}(2.0,orbitals=(2,1),spins=(1,2)),point,fock)
    @test ex|>length==1
    @test ex|>collect==[(2.0,(FIndex(1,1,2,1,2),FIndex(1,1,1,2,1)))]

    bond=Bond(1,Point(PID(1,1),(0.0,0.0)),Point(PID(1,2),(0.0,0.0)))
    sfock,efock=Fock(atom=1,norbital=2,nspin=2,nnambu=2),Fock(atom=2,norbital=2,nspin=2,nnambu=2)
    ex=expand(FockCoupling{2}(2.0,atoms=(2,2)),bond,sfock,efock)
    @test ex|>length==0
    @test ex|>collect==[]

    ex=expand(FockCoupling{2}(2.0,atoms=(2,1)),bond,sfock,efock)
    @test ex|>length==4
    @test ex|>collect==[(2.0,(FIndex(1,2,1,1,2),FIndex(1,1,1,1,1))),
                        (2.0,(FIndex(1,2,1,2,2),FIndex(1,1,1,2,1))),
                        (2.0,(FIndex(1,2,2,1,2),FIndex(1,1,2,1,1))),
                        (2.0,(FIndex(1,2,2,2,2),FIndex(1,1,2,2,1)))
                        ]

    ex=expand(FockCoupling{2}(2.0,orbitals=(1,2)),bond,sfock,efock)
    @test ex|>length==2
    @test ex|>collect==[(2.0,(FIndex(1,2,1,1,2),FIndex(1,1,2,1,1))),
                        (2.0,(FIndex(1,2,1,2,2),FIndex(1,1,2,2,1)))
                        ]

    ex=expand(FockCoupling{2}(2.0,spins=(1,2)),bond,sfock,efock)
    @test ex|>length==2
    @test ex|>collect==[(2.0,(FIndex(1,2,1,1,2),FIndex(1,1,1,2,1))),
                        (2.0,(FIndex(1,2,2,1,2),FIndex(1,1,2,2,1)))
                        ]

    ex=expand(FockCoupling{2}(2.0,orbitals=(2,1),spins=(1,2)),bond,sfock,efock)
    @test ex|>length==1
    @test ex|>collect==[(2.0,(FIndex(1,2,2,1,2),FIndex(1,1,1,2,1)))]
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
