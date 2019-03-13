using Test
using Random: seed!
using StaticArrays: SVector
using QuantumLattices.Essentials.Spatials
using QuantumLattices.Interfaces: decompose,rank,dimension

@testset "distance" begin
    @test distance([0.0,0.0],[1.0,1.0])≈sqrt(2.0)
end

@testset "azimuthd" begin
    @test azimuthd([+1.0])≈0.0
    @test azimuthd([-1.0])≈180.0
    @test azimuthd([+1.0,+1.0])≈45.0
    @test azimuthd([-1.0,-1.0])≈225.0
    @test azimuthd([+1.0,+1.0,-2.0])≈45.0
    @test azimuthd([+1.0,-1.0,+2.0])≈315.0
end

@testset "azimuth" begin
    @test azimuth([+1.0])≈0.0
    @test azimuth([-1.0])≈pi
    @test azimuth([+1.0,+1.0])≈1//4*pi
    @test azimuth([-1.0,-1.0])≈5//4*pi
    @test azimuth([+1.0,+1.0,-2.0])≈1//4*pi
    @test azimuth([+1.0,-1.0,+2.0])≈7//4*pi
end

@testset "polard" begin
    @test polard([1.0,0.0,1.0])≈45.0
end

@testset "polard" begin
    @test polar([1.0,0.0,1.0])≈1//4*pi
end

@testset "volume" begin
    @test volume([1,0],[0,1],[1,1])==0
    @test volume([1,0,0],[0,1,0],[0,0,1])==1
end

@testset "isparallel" begin
    @test isparallel([0.0,0.0],[1.0,1.0])==1
    @test isparallel([2.0,2.0],[1.0,1.0])==1
    @test isparallel([-2.0,-2.0],[1.0,1.0])==-1
    @test isparallel([1.0,-1.0],[1.0,1.0])==0
end

@testset "isonline" begin
    seed!()
    p1,p2=[0.0,0.0],[1.0,1.0]
    @test isonline([0.0,0.0],p1,p2,ends=(true,rand(Bool)))==true
    @test isonline([0.0,0.0],p1,p2,ends=(false,rand(Bool)))==false
    @test isonline([1.0,1.0],p1,p2,ends=(rand(Bool),true))==true
    @test isonline([1.0,1.0],p1,p2,ends=(rand(Bool),false))==false
    @test isonline([0.5,0.5],p1,p2,ends=(rand(Bool),rand(Bool)))==true
    @test isonline([1.1,1.1],p1,p2,ends=(rand(Bool),rand(Bool)))==false
end

@testset "decompose" begin
    a,c=randn(3),randn()
    @test all(decompose(c*a,a).≈(c,))

    a1,a2=randn(2),randn(2)
    c1,c2=randn(2)
    @test all(decompose(c1*a1+c2*a2,a1,a2).≈(c1,c2))

    a1,a2=randn(3),randn(3)
    c1,c2=randn(2)
    @test all(decompose(c1*a1+c2*a2,a1,a2).≈(c1,c2))

    a1,a2,a3=randn(3),randn(3),randn(3)
    c1,c2,c3=randn(3)
    @test all(decompose(c1*a1+c2*a2+c3*a3,a1,a2,a3).≈(c1,c2,c3))
end

@testset "isintratriangle" begin
    seed!()
    p1,p2,p3=[-1.0,0.0],[0.0,1.0],[1.0,0.0]
    for (i,p) in enumerate((p1,p2,p3))
        @test isintratriangle(p,p1,p2,p3,vertexes=Tuple(j==i ? true : rand(Bool) for j=1:3),edges=(rand(Bool),rand(Bool),rand(Bool)))==true
        @test isintratriangle(p,p1,p2,p3,vertexes=Tuple(j==i ? false : rand(Bool) for j=1:3),edges=(rand(Bool),rand(Bool),rand(Bool)))==false
    end
    for (i,p) in enumerate(((p1+p2)/2,(p2+p3)/2,(p3+p1)/2))
        @test isintratriangle(p,p1,p2,p3,vertexes=(rand(Bool),rand(Bool),rand(Bool)),edges=Tuple(j==i ? true : rand(Bool) for j=1:3))==true
        @test isintratriangle(p,p1,p2,p3,vertexes=(rand(Bool),rand(Bool),rand(Bool)),edges=Tuple(j==i ? false : rand(Bool) for j=1:3))==false
    end
    @assert isintratriangle([0.0,0.5],p1,p2,p3,vertexes=(rand(Bool),rand(Bool),rand(Bool)),edges=(rand(Bool),rand(Bool),rand(Bool)))==true
    @assert isintratriangle([0.0,1.5],p1,p2,p3,vertexes=(rand(Bool),rand(Bool),rand(Bool)),edges=(rand(Bool),rand(Bool),rand(Bool)))==false
end

@testset "issubordinate" begin
    @test issubordinate([2.0],[[1.0]])==true
    @test issubordinate([1.5],[[1.0]])==false
    @test issubordinate([0.5,3.5],[[0.5,0.5],[-0.5,0.5]])==true
    @test issubordinate([0.5,3.0],[[0.5,0.5],[-0.5,0.5]])==false
    @test issubordinate([4.0,4.0,4.0],[[1.0,0.0,0.0],[0.0,1.0,0.0],[0.0,0.0,1.0]])==true
    @test issubordinate([4.0,4.0,2.5],[[1.0,0.0,0.0],[0.0,1.0,0.0],[0.0,0.0,1.0]])==false
end

@testset "reciprocals" begin
    @test reciprocals([[1.0]])≈[[1.0]]*2pi
    @test reciprocals([[1.0,0.0],[0.0,1.0]])≈[[1.0,0.0],[0.0,1.0]]*2pi
    @test reciprocals([[1.0,0.0,0.0],[0.0,1.0,0.0],[0.0,0.0,1.0]])≈[[1.0,0.0,0.0],[0.0,1.0,0.0],[0.0,0.0,1.0]]*2pi
end

@testset "translate" begin
    @test translate([1.0 2.0; 3.0 4.0],[1.0,2.0])==[2.0 3.0; 5.0 6.0]
end

@testset "rotate" begin
    @test rotate([1.0 2.0; 1.0 2.0],pi/4,axis=([1.0,1.0],(0.0,0.0)))≈[1.0 1.0; 1.0 1.0+√2.0]
    @test rotate(reshape([2.0,2.0,2.0],3,1),pi/12,axis=([0.0,1.0,1.0],(pi/2,0.0)))≈reshape([2.0,1.0+√2/2,1.0+√6/2],3,1)
    @test rotate(reshape([2.0,2.0,2.0],3,1),pi/12,axis=([1.0,0.0,1.0],(pi/2,pi/2)))≈reshape([1.0+√6/2,2.0,1.0+√2/2],3,1)
    @test rotate(reshape([2.0,2.0,2.0],3,1),pi/12,axis=([1.0,1.0,0.0],(0.0,0.0)))≈reshape([1.0+√2/2,1.0+√6/2,2.0],3,1)
end

@testset "tile" begin
    @test tile(reshape([0.0,0.0],2,1),[[1.0,0.0],[0.0,1.0]],((0.5,0.5),(-0.5,-0.5)))==[0.5 -0.5; 0.5 -0.5]
end

@testset "minimumlengths" begin
    @test minimumlengths(reshape([0.0,0.0],2,1),[[1.0,0.0],[0.0,1.0]],7)≈[1.0,√2,2.0,√5,2*√2,3.0,√10]
end

@testset "intralinks" begin
    ps,a1,a2=[0.0 0.5; 0.0 0.5],[1.0,0.0],[0.0,1.0]
    neighbors=Dict(1=>1.0)
    links=[(1,1,1,SVector(0.0,-1.0)),(1,1,1,SVector(-1.0,0.0)),(1,2,2,SVector(0.0,-1.0)),(1,2,2,SVector(-1.0,0.0))]
    @test intralinks(ps,[a1,a2],neighbors)==links
end

@testset "interlinks" begin
    ps1=[0.0 1.0; 0.0 0.0]
    ps2=[0.0 1.0; 1.0 1.0]
    neighbors=Dict(0=>0.0,1=>1.0)
    links=[(1,1,1,SVector(0.0,0.0)),(1,2,2,SVector(0.0,0.0))]
    @test interlinks(ps1,ps2,neighbors)==links
end

@testset "PID" begin
    @test PID(scope="tz",site=1)|>string=="PID(\"tz\",1)"
end

@testset "Point" begin
    @test Point(PID(0,1),(0.0,0.0),(0.0,0.0))==Point(PID(0,1),[0.0,0.0],[0.0,0.0])
    @test isequal(Point(PID(0,1),(0.0,0.0),(0.0,0.0)),Point(PID(0,1),[0.0,0.0],[0.0,0.0]))
    @test Point{PID{Int},2}|>length==1
    @test Point{PID{Int},2}|>eltype==Point{PID{Int},2}
    @test Point{PID{Int},2}|>rank==1
    @test Point{PID{Int},2}|>pidtype==PID{Int}
    @test Point{PID{Int},2}|>dimension==2
    @test Point{PID{Int},2}|>neighbor==0
    @test Point(PID(0,1),(0.0,0.0))|>length==1
    @test Point(PID(0,1),(0.0,0.0))|>eltype==Point{PID{Int},2}
    @test Point(PID(0,1),(0.0,0.0))|>rank==1
    @test Point(PID(0,1),(0.0,0.0))|>pidtype==PID{Int}
    @test Point(PID(0,1),(0.0,0.0))|>dimension==2
    @test Point(PID(0,1),(0.0,0.0))|>neighbor==0
    @test Point(PID(0,1),(0.0,0.0),(0.0,0.0))|>string=="Point(PID(0,1),[0.0,0.0],[0.0,0.0])"
    @test Point(PID(0,1),(0.0,0.0))|>collect==[Point(PID(0,1),(0.0,0.0))]
end

@testset "Bond" begin
    bond=Bond(1,Point(PID(1,1),(0.0,0.0),(0.0,0.0)),Point(PID(1,2),(0.0,1.0),(0.0,1.0)))
    @test bond|>deepcopy==bond
    @test isequal(bond|>deepcopy,bond)
    @test bond|>string=="Bond(1,Point(PID(1,1),[0.0,0.0],[0.0,0.0]),Point(PID(1,2),[0.0,1.0],[0.0,1.0]))"
    @test bond|>reverse==Bond(1,Point(PID(1,2),(0.0,1.0),(0.0,1.0)),Point(PID(1,1),(0.0,0.0),(0.0,0.0)))
    @test bond|>length==2
    @test bond|>typeof|>length==2
    @test bond|>eltype==Point{PID{Int},2}
    @test bond|>typeof|>eltype==Point{PID{Int},2}
    @test bond|>rank==2
    @test bond|>typeof|>rank==2
    @test bond|>pidtype==PID{Int}
    @test bond|>neighbor==1
    @test bond|>typeof|>pidtype==PID{Int}
    @test bond|>dimension==2
    @test bond|>typeof|>dimension==2
    @test bond|>rcoord==[0.0,1.0]
    @test bond|>icoord==[0.0,1.0]
    @test bond|>isintracell==false
    @test bond|>collect==[Point(PID(1,2),(0.0,1.0),(0.0,1.0)),Point(PID(1,1),(0.0,0.0),(0.0,0.0))]
end

@testset "Lattice" begin
    lattice=Lattice("Tuanzi",[Point(PID(1,1),(0.5,0.5),(0.0,0.0))],vectors=[[1.0,0.0],[0.0,1.0]],neighbors=1)
    @test lattice|>deepcopy==lattice
    @test isequal(lattice|>deepcopy,lattice)
    @test lattice|>string=="Lattice(Tuanzi)"
    @test lattice|>length==1
    @test lattice|>dimension==2
    @test lattice|>keytype==PID{Int}
    @test lattice|>valtype==Point{PID{Int},2}
    @test lattice[LatticeIndex{'R'}(1)]==SVector(0.5,0.5)
    @test lattice[LatticeIndex{'R'}(PID(1,1))]==SVector(0.5,0.5)
    @test lattice[LatticeIndex{'I'}(1)]==SVector(0.0,0.0)
    @test lattice[LatticeIndex{'I'}(PID(1,1))]==SVector(0.0,0.0)
    @test lattice[LatticeIndex{'P'}(1)]==Point(PID(1,1),(0.5,0.5),(0.0,0.0))
    @test lattice[LatticeIndex{'P'}(PID(1,1))]==Point(PID(1,1),(0.5,0.5),(0.0,0.0))
    @test lattice|>nneighbor==1
    @test lattice|>bonds==[ Point(PID(1,1),(0.5,0.5),(0.0,0.0)),
                            Bond(1,Point(PID(1,1),[0.5,0.5],[0.0,0.0]),Point(PID(1,1),[0.5,-0.5],[0.0,-1.0])),
                            Bond(1,Point(PID(1,1),[0.5,0.5],[0.0,0.0]),Point(PID(1,1),[-0.5,0.5],[-1.0,0.0]))
                            ]
    lattice=Lattice(    "SuperTuanzi",
                        [   Lattice("Tuanzi1",[Point(PID(1,1),(0.0,0.0))],neighbors=0),
                            Lattice("Tuanzi2",[Point(PID(2,1),(0.5,0.5))],neighbors=0),
                            ],
                        vectors=    [[1.0,0.0],[0.0,1.0]],
                        neighbors=  2,
                            )
    @test setdiff(bonds(lattice),[bonds(lattice,zerothbonds);bonds(lattice,insidebonds);bonds(lattice,acrossbonds)])|>length==0
end

@testset "SuperLattice" begin
    lattice=SuperLattice(   "SuperTuanzi",
                            [   Lattice("TuanziSys",[Point(PID(1,1),(0.0,0.0)),Point(PID(1,2),(0.5,0.5))],neighbors=Dict(0=>0.0,1=>√2/2)),
                                Lattice("TuanziEnv",[Point(PID(2,1),(-0.05,-0.05)),Point(PID(2,2),(0.55,0.55))],neighbors=Dict(0=>0.0)),
                                ],
                            vectors=[[1.0,0.0],[0.0,1.0]],
                            neighbors=Dict(1=>√2/2,-1=>√2/20)
                        )
    @test setdiff(bonds(lattice),[bonds(lattice,zerothbonds);bonds(lattice,insidebonds);bonds(lattice,acrossbonds)])|>length==0
    @test setdiff(bonds(lattice,insidebonds),[bonds(lattice,intrabonds);bonds(lattice,interbonds)])|>length==0
end

@testset "Cylinder" begin
    cylinder=Cylinder{PID{String}}("Tuanzi",[0.0 0.0; 0.0 1.0],[1.0,0.0],vector=[0.0,2.0],neighbors=1)
    insert!(cylinder,"A","B")
    insert!(cylinder,"C3",cut=2,scopes=["C1","C1","C2","C2"])
    @test cylinder.pids==[PID("C1",1),PID("C3",2),PID("C3",1),PID("C1",2),PID("C2",1),PID("C2",2)]
    @test cylinder.rcoords==[-1.0 -1.0 0.0 0.0 1.0 1.0; 0.0 1.0 0.0 1.0 0.0 1.0]
    @test cylinder.icoords==[0.0 0.0 0.0 0.0 0.0 0.0; 0.0 0.0 0.0 0.0 0.0 0.0]
    lattice=cylinder("C1","C2","C3")
    @test lattice.pids==[PID("C1",1),PID("C1",2),PID("C2",1),PID("C2",2),PID("C3",1),PID("C3",2)]
    @test lattice.rcoords==[-1.0 -1.0 0.0 0.0 1.0 1.0; 0.0 1.0 0.0 1.0 0.0 1.0]
    @test lattice.icoords==[0.0 0.0 0.0 0.0 0.0 0.0; 0.0 0.0 0.0 0.0 0.0 0.0]
end
