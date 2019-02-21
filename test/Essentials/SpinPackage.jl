using Hamiltonian.Essentials.SpinPackage
using Hamiltonian.Essentials.Spatials: PID
using Hamiltonian.Essentials.Terms: Couplings,@subscript
using Hamiltonian.Mathematics.VectorSpaces: IsMultiIndexable,MultiIndexOrderStyle

@testset "SID" begin
    @test SID|>fieldnames==(:orbital,:spin,:tag)
    @test SID(orbital=1,spin=0.5,tag='z')'==SID(orbital=1,spin=0.5,tag='z')
    @test SID(orbital=1,spin=0.5,tag='x')'==SID(orbital=1,spin=0.5,tag='x')
    @test SID(orbital=1,spin=0.5,tag='y')'==SID(orbital=1,spin=0.5,tag='y')
    @test SID(orbital=1,spin=0.5,tag='+')'==SID(orbital=1,spin=0.5,tag='-')
    @test SID(orbital=1,spin=0.5,tag='-')'==SID(orbital=1,spin=0.5,tag='+')
    @test SID(orbital=1,spin=0.5,tag='i')'==SID(orbital=1,spin=0.5,tag='i')
end

@testset "Spin" begin
    spin=Spin(atom=1,norbital=2,spin=1.0)
    @test IsMultiIndexable(Spin)==IsMultiIndexable(true)
    @test MultiIndexOrderStyle(Spin)==MultiIndexOrderStyle('C')
    @test dims(spin)==(2,6)
    @test inds(SID(1,1.0,'i'),spin)==(1,1)
    @test SID((1,2),spin)==SID(1,1.0,'x')
    @test spin|>collect==[  SID(1,1.0,'i'),SID(1,1.0,'x'),SID(1,1.0,'y'),SID(1,1.0,'z'),SID(1,1.0,'+'),SID(1,1.0,'-'),
                            SID(2,1.0,'i'),SID(2,1.0,'x'),SID(2,1.0,'y'),SID(2,1.0,'z'),SID(2,1.0,'+'),SID(2,1.0,'-')
                            ]
end

@testset "SIndex" begin
    @test SIndex|>fieldnames==(:scope,:site,:orbital,:spin,:tag)
    @test union(PID{Char},SID)==SIndex{Char}
end

@testset "SCID" begin
    @test SCID(center=1,atom=1,orbital=1,tag='z')|>typeof|>fieldnames==(:center,:atom,:orbital,:tag,:subscript)
end

@testset "SpinCoupling" begin
    @test SpinCoupling{2}(1.0,tags=('+','-'))|>string=="SpinCoupling(value=1.0,tags=(+,-))"
    @test SpinCoupling{2}(1.0,atoms=(1,1),tags=('z','z'))|>string=="SpinCoupling(value=1.0,atoms=(1,1),tags=(z,z))"
    @test SpinCoupling{2}(1.0,atoms=(1,1),orbitals=(1,2),tags=('-','+'))|>string=="SpinCoupling(value=1.0,atoms=(1,1),orbitals=(1,2),tags=(-,+))"
    @test SpinCoupling{2}(2.0,tags=('x','y'))|>repr=="2.0 SxSy"

    sc1=SpinCoupling{2}(1.5,tags=('+','-'),atoms=(1,2),orbitals=(@subscript (x,y)=>(x,y) with x>y),centers=(1,2))
    sc2=SpinCoupling{2}(2.0,tags=('+','-'),atoms=(1,2),orbitals=(@subscript (x,y)=>(x,y) with x<y),centers=(1,2))
    @test sc1|>repr=="1.5 S+S- sl(1:2)⊗ob(x:y)@(1-2) with $(sc1.id[1].subscript) && $(sc1.id[2].subscript)"
    @test sc2|>repr=="2.0 S+S- sl(1:2)⊗ob(x:y)@(1-2) with $(sc2.id[1].subscript) && $(sc2.id[2].subscript)"

    sc=sc1*sc2
    @test sc|>repr=="3.0 S+S-S+S- sl(1:2:1:2)⊗ob(x:y:x:y)@(1-2-1-2) with $(sc1.id[1].subscript) && $(sc1.id[2].subscript) && $(sc2.id[1].subscript) && $(sc2.id[2].subscript)"

    ex=expand(SpinCoupling{2}(2.0,tags=('+','-'),atoms=(1,1)),PID(1,1),Spin(atom=2,norbital=2,spin=1.0))
    @test collect(ex)==[]

    ex=expand(SpinCoupling{2}(2.0,tags=('+','-'),atoms=(1,2),orbitals=(1,2)),(PID(1,1),PID(1,2)),(Spin(atom=1,norbital=2,spin=1.0),Spin(atom=2,norbital=2,spin=1.0)))
    @test IsMultiIndexable(typeof(ex))==IsMultiIndexable(true)
    @test MultiIndexOrderStyle(typeof(ex))==MultiIndexOrderStyle('C')
    @test dims(ex)==(1,)
    @test collect(ex)==[(2.0,(SIndex(1,1,1,1.0,'+'),SIndex(1,2,2,1.0,'-')))]

    ex=expand(SpinCoupling{4}(2.0,tags=('+','-','+','-'),orbitals=(@subscript (α,β)=>(α,α,β,β) with α<β)),PID(1,1),Spin(norbital=3,spin=1.0))
    @test dims(ex)==(3,)
    @test collect(ex)==[(2.0,(SIndex(1,1,1,1.0,'+'),SIndex(1,1,1,1.0,'-'),SIndex(1,1,2,1.0,'+'),SIndex(1,1,2,1.0,'-'))),
                        (2.0,(SIndex(1,1,1,1.0,'+'),SIndex(1,1,1,1.0,'-'),SIndex(1,1,3,1.0,'+'),SIndex(1,1,3,1.0,'-'))),
                        (2.0,(SIndex(1,1,2,1.0,'+'),SIndex(1,1,2,1.0,'-'),SIndex(1,1,3,1.0,'+'),SIndex(1,1,3,1.0,'-')))
    ]
end

@testset "Heisenberg" begin
    @test Heisenberg(orbitals=(1,2))==Couplings(SpinCoupling{2}(1.0,tags=('z','z'),orbitals=(1,2)),
                                                SpinCoupling{2}(0.5,tags=('+','-'),orbitals=(1,2)),
                                                SpinCoupling{2}(0.5,tags=('-','+'),orbitals=(1,2))
                                                )
end

@testset "Ising" begin
    @test Ising('x',atoms=(1,2))==Couplings(SpinCoupling{2}(1.0,tags=('x','x'),atoms=(1,2)))
    @test Ising('y',atoms=(1,2))==Couplings(SpinCoupling{2}(1.0,tags=('y','y'),atoms=(1,2)))
    @test Ising('z',atoms=(1,2))==Couplings(SpinCoupling{2}(1.0,tags=('z','z'),atoms=(1,2)))
end

@testset "Sᵅ" begin
    @test Sˣ(center=1)==Couplings(SpinCoupling{1}(1.0,tags=('x',),centers=(1,)))
    @test Sʸ(atom=1)==Couplings(SpinCoupling{1}(1.0,tags=('y',),atoms=(1,)))
    @test Sᶻ(orbital=1)==Couplings(SpinCoupling{1}(1.0,tags=('z',),orbitals=(1,)))
end

@testset "matrix" begin
    @test isapprox(matrix(SID(1,0.5,'i')),[[1.0,0.0] [0.0,1.0]])
    @test isapprox(matrix(SID(1,0.5,'z')),[[-0.5,0.0] [0.0,0.5]])
    @test isapprox(matrix(SID(1,0.5,'x')),[[0.0,0.5] [0.5,0.0]])
    @test isapprox(matrix(SID(1,0.5,'y')),[[0.0,-0.5im] [0.5im,0.0]])
    @test isapprox(matrix(SID(1,0.5,'+')),[[0.0,1.0] [0.0,0.0]])
    @test isapprox(matrix(SID(1,0.5,'-')),[[0.0,0.0] [1.0,0.0]])

    @test isapprox(matrix(SID(1,1.0,'i')),[[1.0,0.0,0.0] [0.0,1.0,0.0] [0.0,0.0,1.0]])
    @test isapprox(matrix(SID(1,1.0,'z')),[[-1.0,0.0,0.0] [0.0,0.0,0.0] [0.0,0.0,1.0]])
    @test isapprox(matrix(SID(1,1.0,'x')),[[0.0,√2/2,0.0] [√2/2,0.0,√2/2] [0.0,√2/2,0.0]])
    @test isapprox(matrix(SID(1,1.0,'y')),[[0.0,-√2im/2,0.0] [√2im/2,0.0,-√2im/2] [0.0,√2im/2,0.0]])
    @test isapprox(matrix(SID(1,1.0,'+')),[[0.0,√2,0.0] [0.0,0.0,√2] [0.0,0.0,0.0]])
    @test isapprox(matrix(SID(1,1.0,'-')),[[0.0,0.0,0.0] [√2,0.0,0.0] [0.0,√2,0.0]])
end
