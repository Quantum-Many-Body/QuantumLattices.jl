using Test
using StaticArrays: SVector
using QuantumLattices.Essentials.SpinPackage
using QuantumLattices.Essentials.Spatials: PID,Point,Bond
using QuantumLattices.Essentials.DegreesOfFreedom: Table,OID,isHermitian,IDFConfig,Operators,oidtype,otype,optdefaultlatex
using QuantumLattices.Essentials.Terms: Couplings,@subscript,statistics,abbr
using QuantumLattices.Interfaces: dims,inds,expand,matrix
using QuantumLattices.Prerequisites: Float
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
    @test oidtype(SID,Point{2,PID{Int}},Nothing)==OID{SIndex{Int},SVector{2,Float},SVector{2,Float},Nothing}
    @test oidtype(SID,Point{2,PID{Int}},Table)==OID{SIndex{Int},SVector{2,Float},SVector{2,Float},Int}
end

@testset "SOperator" begin
    @test rawelement(SOperator{N,<:Number,<:ID{<:NTuple{N,OID}}} where N)==SOperator
    @test optdefaultlatex(SOperator)==soptdefaultlatex
    opt=SOperator(1.0,(SIndex('a',1,1,0.5,'+'),SIndex('a',1,1,0.5,'-')))
    @test opt|>statistics==opt|>typeof|>statistics=='B'
    @test opt'==SOperator(1.0,(SIndex('a',1,1,0.5,'+'),SIndex('a',1,1,0.5,'-')))
    @test isHermitian(opt)
    @test repr(opt)=="1.0S^{+}_{1,1}S^{-}_{1,1}"
end

@testset "SCID" begin
    @test SCID(center=1,atom=1,orbital=1,tag='z')|>typeof|>fieldnames==(:center,:atom,:orbital,:tag,:subscript)
end

@testset "SpinCoupling" begin
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

@testset "Gamma" begin
    @test Gamma('x',orbitals=(1,1))==SpinCoupling{2}(1.0,tags=('y','z'),orbitals=(1,1))+SpinCoupling{2}(1.0,tags=('z','y'),orbitals=(1,1))
    @test Gamma('y',atoms=(1,2))==SpinCoupling{2}(1.0,tags=('z','x'),atoms=(1,2))+SpinCoupling{2}(1.0,tags=('x','z'),atoms=(1,2))
    @test Gamma('z')==SpinCoupling{2}(1.0,tags=('x','y'))+SpinCoupling{2}(1.0,tags=('y','x'))
end

@testset "Sᵅ" begin
    @test Sˣ(atom=1,orbital=1)==Couplings(SpinCoupling{1}(1.0,tags=('x',),atoms=(1,),orbitals=(1,)))
    @test Sʸ(atom=1)==Couplings(SpinCoupling{1}(1.0,tags=('y',),atoms=(1,)))
    @test Sᶻ(orbital=1)==Couplings(SpinCoupling{1}(1.0,tags=('z',),orbitals=(1,)))
end

@testset "SpinTerm" begin
    term=SpinTerm{1}(:h,1.5,0,couplings=Sᶻ())
    @test term|>abbr==:sp
    @test otype(term|>typeof,OID{SIndex{Int},SVector{2,Float},SVector{2,Float},Nothing})==SOperator{1,Float,ID{NTuple{1,OID{SIndex{Int},SVector{2,Float},SVector{2,Float},Nothing}}}}
    @test otype(term|>typeof,OID{SIndex{Int},SVector{2,Float},SVector{2,Float},Int})==SOperator{1,Float,ID{NTuple{1,OID{SIndex{Int},SVector{2,Float},SVector{2,Float},Int}}}}

    point=Point(PID('a',1),(0.5,0.5),(0.0,0.0))
    config=IDFConfig{Spin}(pid->Spin(atom=pid.site%2,norbital=2,spin=0.5),[point.pid])
    table=Table(config,by=usualspinindextotuple)
    term=SpinTerm{1}(:h,1.5,0,couplings=Sᶻ())
    operators=Operators(SOperator(1.5,ID(OID(SIndex('a',1,1,0.5,'z'),[0.5,0.5],[0.0,0.0],1))),
                        SOperator(1.5,ID(OID(SIndex('a',1,2,0.5,'z'),[0.5,0.5],[0.0,0.0],2)))
    )
    @test expand(term,point,config,table)==operators

    bond=Bond(1,Point(PID('a',1),(0.0,0.0),(0.0,0.0)),Point(PID('b',1),(0.5,0.5),(0.0,0.0)))
    config=IDFConfig{Spin}(pid->Spin(atom=pid.site%2,norbital=2,spin=0.5),[bond.spoint.pid,bond.epoint.pid])
    table=Table(config,by=usualspinindextotuple)
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
