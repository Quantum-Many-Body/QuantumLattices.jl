using Test
using StaticArrays: SVector
using QuantumLattices.Essentials.SpinPackage
using QuantumLattices.Essentials.Spatials: PID,Point,Bond
using QuantumLattices.Essentials.DegreesOfFreedom: Table,OID,isHermitian,IDFConfig,Operators,oidtype,otype,script,iid
using QuantumLattices.Essentials.Terms: Couplings,@subscript,statistics,abbr
using QuantumLattices.Interfaces: dims,inds,expand,matrix,permute,rank
using QuantumLattices.Prerequisites: Float
using QuantumLattices.Mathematics.Combinatorics: Permutations
using QuantumLattices.Mathematics.AlgebraOverFields: ID,rawelement
using QuantumLattices.Mathematics.VectorSpaces: IsMultiIndexable,MultiIndexOrderStyle

@testset "SID" begin
    @test SID|>fieldnames==(:orbital,:spin,:tag)
    @test SID(orbital=1,spin=0.5,tag='z')'==SID(orbital=1,spin=0.5,tag='z')
    @test SID(orbital=1,spin=0.5,tag='x')'==SID(orbital=1,spin=0.5,tag='x')
    @test SID(orbital=1,spin=0.5,tag='y')'==SID(orbital=1,spin=0.5,tag='y')
    @test SID(orbital=1,spin=0.5,tag='+')'==SID(orbital=1,spin=0.5,tag='-')
    @test SID(orbital=1,spin=0.5,tag='-')'==SID(orbital=1,spin=0.5,tag='+')
end

@testset "matrix" begin
    @test isapprox(matrix(SID(1,0.5,'z')),[[-0.5,0.0] [0.0,0.5]])
    @test isapprox(matrix(SID(1,0.5,'x')),[[0.0,0.5] [0.5,0.0]])
    @test isapprox(matrix(SID(1,0.5,'y')),[[0.0,-0.5im] [0.5im,0.0]])
    @test isapprox(matrix(SID(1,0.5,'+')),[[0.0,1.0] [0.0,0.0]])
    @test isapprox(matrix(SID(1,0.5,'-')),[[0.0,0.0] [1.0,0.0]])

    @test isapprox(matrix(SID(1,1.0,'z')),[[-1.0,0.0,0.0] [0.0,0.0,0.0] [0.0,0.0,1.0]])
    @test isapprox(matrix(SID(1,1.0,'x')),[[0.0,√2/2,0.0] [√2/2,0.0,√2/2] [0.0,√2/2,0.0]])
    @test isapprox(matrix(SID(1,1.0,'y')),[[0.0,-√2im/2,0.0] [√2im/2,0.0,-√2im/2] [0.0,√2im/2,0.0]])
    @test isapprox(matrix(SID(1,1.0,'+')),[[0.0,√2,0.0] [0.0,0.0,√2] [0.0,0.0,0.0]])
    @test isapprox(matrix(SID(1,1.0,'-')),[[0.0,0.0,0.0] [√2,0.0,0.0] [0.0,√2,0.0]])
end

@testset "Spin" begin
    spin=Spin(atom=1,norbital=2,spin=1.0)
    @test IsMultiIndexable(Spin)==IsMultiIndexable(true)
    @test MultiIndexOrderStyle(Spin)==MultiIndexOrderStyle('C')
    @test dims(spin)==(2,5)
    @test inds(SID(1,1.0,'z'),spin)==(1,3)
    @test SID((1,1),spin)==SID(1,1.0,'x')
    @test spin|>collect==[  SID(1,1.0,'x'),SID(1,1.0,'y'),SID(1,1.0,'z'),SID(1,1.0,'+'),SID(1,1.0,'-'),
                            SID(2,1.0,'x'),SID(2,1.0,'y'),SID(2,1.0,'z'),SID(2,1.0,'+'),SID(2,1.0,'-')
                            ]
end

@testset "SIndex" begin
    @test SIndex|>fieldnames==(:scope,:site,:orbital,:spin,:tag)
    @test union(PID{Char},SID)==SIndex{Char}
end

@testset "script" begin
    @test script(OID(SIndex('S',1,2,0.5,'z')),Val(:site))==1
    @test script(OID(SIndex('S',1,2,0.5,'z')),Val(:orbital))==2
    @test script(OID(SIndex('S',1,2,0.5,'z')),Val(:spin))==0.5
    @test script(OID(SIndex('S',1,2,0.5,'z')),Val(:tag))=='z'
end

@testset "oidtype" begin
    @test oidtype(SID,Point{2,PID{Int}},Nothing,Val(true))==OID{SIndex{Int},SVector{2,Float},SVector{2,Float},Nothing}
    @test oidtype(SID,Point{2,PID{Int}},Table,Val(true))==OID{SIndex{Int},SVector{2,Float},SVector{2,Float},Int}
    @test oidtype(SID,Point{2,PID{Int}},Nothing,Val(false))==OID{SIndex{Int},Nothing,Nothing,Nothing}
    @test oidtype(SID,Point{2,PID{Int}},Table,Val(false))==OID{SIndex{Int},Nothing,Nothing,Int}
end

@testset "SOperator" begin
    @test rawelement(SOperator{V} where V)==SOperator
    @test latexformat(SOperator)==soptdefaultlatex
    opt=SOperator(1.0,(SIndex('a',1,1,0.5,'+'),SIndex('a',1,1,0.5,'-')))
    @test opt|>statistics==opt|>typeof|>statistics=='B'
    @test opt'==SOperator(1.0,(SIndex('a',1,1,0.5,'+'),SIndex('a',1,1,0.5,'-')))
    @test isHermitian(opt)
    @test repr(opt)=="S^{+}_{1,1}S^{-}_{1,1}"
end

@testset "permute" begin
    soptrep(opt::SOperator)=opt.value*prod([matrix(opt.id[i].index|>iid) for i=1:rank(opt)])
    for S in (0.5,1.0,1.5)
        oids=[OID(SIndex('S',1,2,S,tag),seq=i) for (i,tag) in enumerate(('x','y','z','+','-'))]
        table=Dict(oid.index=>oid.seq for oid in oids)
        for (id1,id2) in Permutations{2}(oids)
            left=soptrep(SOperator(1,ID(id1,id2)))
            right=sum([soptrep(opt) for opt in permute(SOperator,id1,id2,table)])
            @test isapprox(left,right)
        end
    end
    id1=OID(SIndex('S',1,2,0.5,'z'))
    id2=OID(SIndex('S',2,2,0.5,'z'))
    @test permute(SOperator,id1,id2,nothing)==(SOperator(1,ID(id2,id1)),)
end

@testset "SCID" begin
    @test SCID(center=1,atom=1,orbital=1,tag='z')|>typeof|>fieldnames==(:center,:atom,:orbital,:tag,:subscript)
end

@testset "SpinCoupling" begin
    @test rawelement(SpinCoupling{V} where V)==SpinCoupling

    @test SpinCoupling{2}(1.0,tags=('+','-'))|>string=="SpinCoupling(value=1.0,tags=(+,-))"
    @test SpinCoupling{2}(1.0,atoms=(1,1),tags=('z','z'))|>string=="SpinCoupling(value=1.0,atoms=(1,1),tags=(z,z))"
    @test SpinCoupling{2}(1.0,atoms=(1,1),orbitals=(1,2),tags=('-','+'))|>string=="SpinCoupling(value=1.0,atoms=(1,1),orbitals=(1,2),tags=(-,+))"
    @test SpinCoupling{2}(2.0,tags=('x','y'))|>repr=="2.0 SxSy"

    sc1=SpinCoupling{2}(1.5,tags=('+','-'),atoms=(1,2),orbitals=(@subscript (x,y) with x>y),centers=(1,2))
    sc2=SpinCoupling{2}(2.0,tags=('+','-'),atoms=(1,2),orbitals=(@subscript (x,y) with x<y),centers=(1,2))
    @test sc1|>repr=="1.5 S+S- sl(1,2)⊗ob(x,y)@(1,2) with $(sc1.id[1].subscript) && $(sc1.id[2].subscript)"
    @test sc2|>repr=="2.0 S+S- sl(1,2)⊗ob(x,y)@(1,2) with $(sc2.id[1].subscript) && $(sc2.id[2].subscript)"

    sc=sc1*sc2
    @test sc|>repr=="3.0 S+S-S+S- sl(1,2,1,2)⊗ob(x,y,x,y)@(1,2,1,2) with $(sc1.id[1].subscript) && $(sc1.id[2].subscript) && $(sc2.id[1].subscript) && $(sc2.id[2].subscript)"

    ex=expand(SpinCoupling{2}(2.0,tags=('+','-'),atoms=(1,1)),PID(1,1),Spin(atom=2,norbital=2,spin=1.0))
    @test collect(ex)==[]

    ex=expand(SpinCoupling{2}(2.0,tags=('+','-'),atoms=(1,2),orbitals=(1,2)),(PID(1,1),PID(1,2)),(Spin(atom=1,norbital=2,spin=1.0),Spin(atom=2,norbital=2,spin=1.0)))
    @test IsMultiIndexable(typeof(ex))==IsMultiIndexable(true)
    @test MultiIndexOrderStyle(typeof(ex))==MultiIndexOrderStyle('C')
    @test dims(ex)==(1,)
    @test collect(ex)==[(2.0,(SIndex(1,1,1,1.0,'+'),SIndex(1,2,2,1.0,'-')))]

    ex=expand(SpinCoupling{4}(2.0,tags=('+','-','+','-'),orbitals=(@subscript (α,α,β,β) with α<β)),PID(1,1),Spin(norbital=3,spin=1.0))
    @test dims(ex)==(3,)
    @test collect(ex)==[(2.0,(SIndex(1,1,1,1.0,'+'),SIndex(1,1,1,1.0,'-'),SIndex(1,1,2,1.0,'+'),SIndex(1,1,2,1.0,'-'))),
                        (2.0,(SIndex(1,1,1,1.0,'+'),SIndex(1,1,1,1.0,'-'),SIndex(1,1,3,1.0,'+'),SIndex(1,1,3,1.0,'-'))),
                        (2.0,(SIndex(1,1,2,1.0,'+'),SIndex(1,1,2,1.0,'-'),SIndex(1,1,3,1.0,'+'),SIndex(1,1,3,1.0,'-')))
    ]
end

@testset "Heisenberg" begin
    @test Heisenberg(orbitals=(1,2))==Couplings(SpinCoupling{2}(1//1,tags=('z','z'),orbitals=(1,2)),
                                                SpinCoupling{2}(1//2,tags=('+','-'),orbitals=(1,2)),
                                                SpinCoupling{2}(1//2,tags=('-','+'),orbitals=(1,2))
                                                )
    @test Heisenberg("xyz")==Couplings( SpinCoupling{2}(1,tags=('x','x')),
                                        SpinCoupling{2}(1,tags=('y','y')),
                                        SpinCoupling{2}(1,tags=('z','z'))
                                        )
end

@testset "Ising" begin
    @test Ising('x',atoms=(1,2))==Couplings(SpinCoupling{2}(1,tags=('x','x'),atoms=(1,2)))
    @test Ising('y',atoms=(1,2))==Couplings(SpinCoupling{2}(1,tags=('y','y'),atoms=(1,2)))
    @test Ising('z',atoms=(1,2))==Couplings(SpinCoupling{2}(1,tags=('z','z'),atoms=(1,2)))
end

@testset "Gamma" begin
    @test Gamma('x',orbitals=(1,1))==SpinCoupling{2}(1,tags=('y','z'),orbitals=(1,1))+SpinCoupling{2}(1,tags=('z','y'),orbitals=(1,1))
    @test Gamma('y',atoms=(1,2))==SpinCoupling{2}(1,tags=('z','x'),atoms=(1,2))+SpinCoupling{2}(1,tags=('x','z'),atoms=(1,2))
    @test Gamma('z')==SpinCoupling{2}(1,tags=('x','y'))+SpinCoupling{2}(1,tags=('y','x'))
end

@testset "DM" begin
    @test DM('x',orbitals=(1,1))==SpinCoupling{2}(1,tags=('y','z'),orbitals=(1,1))-SpinCoupling{2}(1,tags=('z','y'),orbitals=(1,1))
    @test DM('y',atoms=(1,2))==SpinCoupling{2}(1,tags=('z','x'),atoms=(1,2))-SpinCoupling{2}(1,tags=('x','z'),atoms=(1,2))
    @test DM('z')==SpinCoupling{2}(1,tags=('x','y'))-SpinCoupling{2}(1,tags=('y','x'))
end

@testset "Sᵅ" begin
    @test Sˣ(atom=1,orbital=1)==Couplings(SpinCoupling{1}(1,tags=('x',),atoms=(1,),orbitals=(1,)))
    @test Sʸ(atom=1)==Couplings(SpinCoupling{1}(1,tags=('y',),atoms=(1,)))
    @test Sᶻ(orbital=1)==Couplings(SpinCoupling{1}(1,tags=('z',),orbitals=(1,)))
end

@testset "spincoupling" begin
    sc=sc"1.0 S+S- sl(1,1)⊗ob(α,β)@(1,2) with α<β"
    ob=sc.subscripts[1].identifier
    @test repr(sc)=="1.0 S+S- sl(1,1)⊗ob(α,β)@(1,2) with $ob && $ob"

    sc=sc"1.0 S+S- sl(1,1)⊗ob(α,β)@(1,2)"
    ob=sc.subscripts[1].identifier
    @test repr(sc)=="1.0 S+S- sl(1,1)⊗ob(α,β)@(1,2) with $ob && $ob"

    sc=sc"1.0 S+S- sl(1,1)⊗ob(α,β)"
    ob=sc.subscripts[1].identifier
    @test repr(sc)=="1.0 S+S- sl(1,1)⊗ob(α,β) with $ob && $ob"

    sc=sc"1.0 S+S- sl(1,1)⊗ob(1,2)@(1,2)"
    @test repr(sc)=="1.0 S+S- sl(1,1)⊗ob(1,2)@(1,2)"

    sc=sc"1.0 S+S- sl(1,1)@(1,2)"
    @test repr(sc)=="1.0 S+S- sl(1,1)@(1,2)"

    sc=sc"1.0 S+S- @(1,2)"
    @test repr(sc)=="1.0 S+S- @(1,2)"

    sc=sc"1.0 S+S-"
    @test repr(sc)=="1.0 S+S-"
end

@testset "spincouplings" begin
    @test heisenberg"sl(1,1)⊗ob(1,3)@(1,2)"==Heisenberg(centers=(1,2),atoms=(1,1),orbitals=(1,3))
    @test heisenberg"sl(1,1)⊗ob(1,3)"==Heisenberg(atoms=(1,1),orbitals=(1,3))
    @test heisenberg"@(1,2)"==Heisenberg(centers=(1,2))
    @test heisenberg"ob(1,3)"==Heisenberg(orbitals=(1,3))
    @test heisenberg"sl(1,1)"==Heisenberg(atoms=(1,1))
    @test heisenberg""==heisenberg"+-z"==Heisenberg()
    @test heisenberg"xyz"==Heisenberg("xyz")

    @test ising"x"==Ising('x') && ising"x sl(1,1)⊗ob(1,3)@(1,2)"==Ising('x',centers=(1,2),atoms=(1,1),orbitals=(1,3))
    @test ising"y"==Ising('y') && ising"y sl(1,1)⊗ob(1,3)@(1,2)"==Ising('y',centers=(1,2),atoms=(1,1),orbitals=(1,3))
    @test ising"z"==Ising('z') && ising"z sl(1,1)⊗ob(1,3)@(1,2)"==Ising('z',centers=(1,2),atoms=(1,1),orbitals=(1,3))

    @test gamma"x"==Gamma('x') && gamma"x sl(1,1)⊗ob(1,3)@(1,2)"==Gamma('x',centers=(1,2),atoms=(1,1),orbitals=(1,3))
    @test gamma"y"==Gamma('y') && gamma"y sl(1,1)⊗ob(1,3)@(1,2)"==Gamma('y',centers=(1,2),atoms=(1,1),orbitals=(1,3))
    @test gamma"z"==Gamma('z') && gamma"z sl(1,1)⊗ob(1,3)@(1,2)"==Gamma('z',centers=(1,2),atoms=(1,1),orbitals=(1,3))

    @test dm"x"==DM('x') && dm"x sl(1,1)⊗ob(1,3)@(1,2)"==DM('x',centers=(1,2),atoms=(1,1),orbitals=(1,3))
    @test dm"y"==DM('y') && dm"y sl(1,1)⊗ob(1,3)@(1,2)"==DM('y',centers=(1,2),atoms=(1,1),orbitals=(1,3))
    @test dm"z"==DM('z') && dm"z sl(1,1)⊗ob(1,3)@(1,2)"==DM('z',centers=(1,2),atoms=(1,1),orbitals=(1,3))

    @test sˣ""==Sˣ() && sˣ"sl(1)⊗ob(2)"==Sˣ(atom=1,orbital=2)
    @test sʸ""==Sʸ() && sʸ"sl(1)⊗ob(2)"==Sʸ(atom=1,orbital=2)
    @test sᶻ""==Sᶻ() && sᶻ"sl(1)⊗ob(2)"==Sᶻ(atom=1,orbital=2)
end

@testset "SpinTerm" begin
    term=SpinTerm{1}(:h,1.5,0,couplings=Sᶻ())
    @test term|>abbr==:sp
    @test otype(term|>typeof,OID{SIndex{Int},SVector{2,Float},SVector{2,Float},Nothing})==SOperator{Float,ID{NTuple{1,OID{SIndex{Int},SVector{2,Float},SVector{2,Float},Nothing}}}}
    @test otype(term|>typeof,OID{SIndex{Int},SVector{2,Float},SVector{2,Float},Int})==SOperator{Float,ID{NTuple{1,OID{SIndex{Int},SVector{2,Float},SVector{2,Float},Int}}}}

    point=Point(PID('a',1),(0.5,0.5),(0.0,0.0))
    config=IDFConfig{Spin}(pid->Spin(atom=pid.site%2,norbital=2,spin=0.5),[point.pid])
    table=Table(config,usualspinindextotuple)
    term=SpinTerm{1}(:h,1.5,0,couplings=Sᶻ())
    operators=Operators(SOperator(1.5,ID(OID(SIndex('a',1,1,0.5,'z'),[0.5,0.5],[0.0,0.0],1))),
                        SOperator(1.5,ID(OID(SIndex('a',1,2,0.5,'z'),[0.5,0.5],[0.0,0.0],2)))
    )
    @test expand(term,point,config,table)==operators

    bond=Bond(1,Point(PID('a',1),(0.0,0.0),(0.0,0.0)),Point(PID('b',1),(0.5,0.5),(0.0,0.0)))
    config=IDFConfig{Spin}(pid->Spin(atom=pid.site%2,norbital=2,spin=0.5),[bond.spoint.pid,bond.epoint.pid])
    table=Table(config,usualspinindextotuple)
    term=SpinTerm{2}(:J,1.5,1,couplings=Heisenberg())
    operators=Operators(SOperator(1.50,ID(OID(SIndex('b',1,2,0.5,'z'),[0.5,0.5],[0.0,0.0],4),OID(SIndex('a',1,2,0.5,'z'),[0.0,0.0],[0.0,0.0],2))),
                        SOperator(0.75,ID(OID(SIndex('b',1,2,0.5,'-'),[0.5,0.5],[0.0,0.0],4),OID(SIndex('a',1,2,0.5,'+'),[0.0,0.0],[0.0,0.0],2))),
                        SOperator(0.75,ID(OID(SIndex('b',1,1,0.5,'-'),[0.5,0.5],[0.0,0.0],3),OID(SIndex('a',1,1,0.5,'+'),[0.0,0.0],[0.0,0.0],1))),
                        SOperator(0.75,ID(OID(SIndex('b',1,1,0.5,'+'),[0.5,0.5],[0.0,0.0],3),OID(SIndex('a',1,1,0.5,'-'),[0.0,0.0],[0.0,0.0],1))),
                        SOperator(1.50,ID(OID(SIndex('b',1,1,0.5,'z'),[0.5,0.5],[0.0,0.0],3),OID(SIndex('a',1,1,0.5,'z'),[0.0,0.0],[0.0,0.0],1))),
                        SOperator(0.75,ID(OID(SIndex('b',1,2,0.5,'+'),[0.5,0.5],[0.0,0.0],4),OID(SIndex('a',1,2,0.5,'-'),[0.0,0.0],[0.0,0.0],2)))
    )
    @test expand(term,bond,config,table)==operators
end
