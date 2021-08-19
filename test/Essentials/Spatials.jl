using Test
using Random: seed!
using StaticArrays: SVector
using QuantumLattices.Essentials.Spatials
using QuantumLattices.Essentials: dtype, kind, reset!
using QuantumLattices.Interfaces: decompose, rank, dimension, expand
using QuantumLattices.Prerequisites: Float
using QuantumLattices.Prerequisites.Traits: contentnames, getcontent
using QuantumLattices.Prerequisites.SimpleTrees: leaves
using QuantumLattices.Mathematics.QuantumNumbers: Momentum2D, AbelianNumbers

@testset "distance" begin
    @test distance([0.0, 0.0], [1.0, 1.0]) ≈ sqrt(2.0)
end

@testset "azimuthd" begin
    @test azimuthd([+1.0]) ≈ 0.0
    @test azimuthd([-1.0]) ≈ 180.0
    @test azimuthd([+1.0, +1.0]) ≈ 45.0
    @test azimuthd([-1.0, -1.0]) ≈ 225.0
    @test azimuthd([+1.0, +1.0, -2.0]) ≈ 45.0
    @test azimuthd([+1.0, -1.0, +2.0]) ≈ 315.0
end

@testset "azimuth" begin
    @test azimuth([+1.0]) ≈ 0.0
    @test azimuth([-1.0]) ≈ pi
    @test azimuth([+1.0, +1.0]) ≈ 1//4*pi
    @test azimuth([-1.0, -1.0]) ≈ 5//4*pi
    @test azimuth([+1.0, +1.0, -2.0]) ≈ 1//4*pi
    @test azimuth([+1.0, -1.0, +2.0]) ≈ 7//4*pi
end

@testset "polard" begin
    @test polard([1.0, 0.0, 1.0]) ≈ 45.0
end

@testset "polard" begin
    @test polar([1.0, 0.0, 1.0]) ≈ 1//4*pi
end

@testset "volume" begin
    @test volume([1, 0], [0, 1], [1, 1]) == 0
    @test volume([1, 0, 0], [0, 1, 0], [0, 0, 1]) == 1
end

@testset "isparallel" begin
    @test isparallel([0.0, 0.0], [1.0, 1.0]) == 1
    @test isparallel([2.0, 2.0], [1.0, 1.0]) == 1
    @test isparallel([-2.0, -2.0], [1.0, 1.0]) == -1
    @test isparallel([1.0, -1.0], [1.0, 1.0]) == 0
end

@testset "isonline" begin
    seed!()
    p1, p2 = [0.0, 0.0], [1.0, 1.0]
    @test isonline([0.0, 0.0], p1, p2, ends = (true, rand(Bool))) == true
    @test isonline([0.0, 0.0], p1, p2, ends = (false, rand(Bool))) == false
    @test isonline([1.0, 1.0], p1, p2, ends = (rand(Bool), true)) == true
    @test isonline([1.0, 1.0], p1, p2, ends = (rand(Bool), false)) == false
    @test isonline([0.5, 0.5], p1, p2, ends = (rand(Bool), rand(Bool))) == true
    @test isonline([1.1, 1.1], p1, p2, ends = (rand(Bool), rand(Bool))) == false
end

@testset "decompose" begin
    a, c = randn(3), randn()
    @test all(decompose(c*a, a) .≈ (c,))

    a1, a2 = randn(2), randn(2)
    c1, c2 = randn(2)
    @test all(decompose(c1*a1+c2*a2, a1, a2) .≈ (c1, c2))

    a1, a2 = randn(3), randn(3)
    c1, c2 = randn(2)
    @test all(decompose(c1*a1+c2*a2, a1, a2) .≈ (c1, c2))

    a1, a2, a3 = randn(3), randn(3), randn(3)
    c1, c2, c3 = randn(3)
    @test all(decompose(c1*a1+c2*a2+c3*a3, a1, a2, a3) .≈ (c1, c2, c3))
end

@testset "isintratriangle" begin
    seed!()
    p1, p2, p3 = [-1.0, 0.0], [0.0, 1.0], [1.0, 0.0]
    for (i, p) in enumerate((p1, p2, p3))
        @test isintratriangle(p, p1, p2, p3, vertexes=Tuple(j == i ? true : rand(Bool) for j = 1:3), edges=(rand(Bool), rand(Bool), rand(Bool))) == true
        @test isintratriangle(p, p1, p2, p3, vertexes=Tuple(j == i ? false : rand(Bool) for j = 1:3), edges=(rand(Bool), rand(Bool), rand(Bool))) == false
    end
    for (i, p) in enumerate(((p1+p2)/2, (p2+p3)/2, (p3+p1)/2))
        @test isintratriangle(p, p1, p2, p3, vertexes=(rand(Bool), rand(Bool), rand(Bool)), edges=Tuple(j == i ? true : rand(Bool) for j = 1:3)) == true
        @test isintratriangle(p, p1, p2, p3, vertexes=(rand(Bool), rand(Bool), rand(Bool)), edges=Tuple(j == i ? false : rand(Bool) for j = 1:3)) == false
    end
    @assert isintratriangle([0.0, 0.5], p1, p2, p3, vertexes=(rand(Bool), rand(Bool), rand(Bool)), edges=(rand(Bool), rand(Bool), rand(Bool))) == true
    @assert isintratriangle([0.0, 1.5], p1, p2, p3, vertexes=(rand(Bool), rand(Bool), rand(Bool)), edges=(rand(Bool), rand(Bool), rand(Bool))) == false
end

@testset "issubordinate" begin
    @test issubordinate([2.0], [[1.0]]) == true
    @test issubordinate([1.5], [[1.0]]) == false
    @test issubordinate([0.5, 3.5], [[0.5, 0.5], [-0.5, 0.5]]) == true
    @test issubordinate([0.5, 3.0], [[0.5, 0.5], [-0.5, 0.5]]) == false
    @test issubordinate([4.0, 4.0, 4.0], [[1.0, 0.0, 0.0], [0.0, 1.0, 0.0], [0.0, 0.0, 1.0]]) == true
    @test issubordinate([4.0, 4.0, 2.5], [[1.0, 0.0, 0.0], [0.0, 1.0, 0.0], [0.0, 0.0, 1.0]]) == false
end

@testset "reciprocals" begin
    @test reciprocals([[1.0]]) ≈ [[1.0]]*2pi
    @test reciprocals([[1.0, 0.0], [0.0, 1.0]]) ≈ [[1.0, 0.0], [0.0, 1.0]]*2pi
    @test reciprocals([[1.0, 0.0, 0.0], [0.0, 1.0, 0.0], [0.0, 0.0, 1.0]]) ≈ [[1.0, 0.0, 0.0], [0.0, 1.0, 0.0], [0.0, 0.0, 1.0]]*2pi
end

@testset "translate" begin
    @test translate([1.0 2.0; 3.0 4.0], [1.0, 2.0]) == [2.0 3.0; 5.0 6.0]
end

@testset "rotate" begin
    @test rotate([1.0 2.0; 1.0 2.0], pi/4, axis=([1.0, 1.0], (0.0, 0.0))) ≈ [1.0 1.0; 1.0 1.0+√2.0]
    @test rotate(reshape([2.0, 2.0, 2.0], 3, 1), pi/12, axis=([0.0, 1.0, 1.0], (pi/2, 0.0))) ≈ reshape([2.0, 1.0+√2/2, 1.0+√6/2], 3, 1)
    @test rotate(reshape([2.0, 2.0, 2.0], 3, 1), pi/12, axis=([1.0, 0.0, 1.0], (pi/2, pi/2))) ≈ reshape([1.0+√6/2, 2.0, 1.0+√2/2], 3, 1)
    @test rotate(reshape([2.0, 2.0, 2.0], 3, 1), pi/12, axis=([1.0, 1.0, 0.0], (0.0, 0.0))) ≈ reshape([1.0+√2/2, 1.0+√6/2, 2.0], 3, 1)
end

@testset "tile" begin
    @test tile(reshape([0.0, 0.0], 2, 1), [[1.0, 0.0], [0.0, 1.0]], ((0.5, 0.5), (-0.5, -0.5))) == [0.5 -0.5; 0.5 -0.5]
end

@testset "minimumlengths" begin
    @test minimumlengths(reshape([0.0, 0.0], 2, 1), [[1.0, 0.0], [0.0, 1.0]], 7) ≈ [1.0, √2, 2.0, √5, 2*√2, 3.0, √10]
end

@testset "intralinks" begin
    ps, a1, a2 = [0.0 0.5; 0.0 0.5], [1.0, 0.0], [0.0, 1.0]
    neighbors = Dict(1=>1.0)
    links = [(1, 1, 1, SVector(0.0, -1.0)), (1, 1, 1, SVector(-1.0, 0.0)), (1, 2, 2, SVector(0.0, -1.0)), (1, 2, 2, SVector(-1.0, 0.0))]
    @test intralinks(ps, [a1, a2], neighbors) == links
end

@testset "interlinks" begin
    ps1 = [0.0 1.0; 0.0 0.0]
    ps2 = [0.0 1.0; 1.0 1.0]
    neighbors = Dict(1=>1.0)
    links = [(1, 1, 1), (1, 2, 2)]
    @test interlinks(ps1, ps2, neighbors) == links
end

@testset "PID" begin
    @test PID(1) == PID('T', 1) == PID(scope='T', site=1)
    @test PID(scope="tz", site=1)|>string == "PID(\"tz\", 1)"
end

@testset "Point" begin
    point = Point(PID(0, 1), (0.0, 0.0), (0.0, 0.0))
    @test point|>dtype == point|>typeof|>dtype == Float
    @test point == Point(PID(0, 1), [0.0, 0.0], [0.0, 0.0])
    @test isequal(point, Point{2}(PID(0, 1), [0.0, 0.0], [0.0, 0.0]))
    @test point|>length == point|>typeof|>length == 1
    @test point|>eltype == point|>typeof|>eltype == Point{2, PID{Int}, Float}
    @test point|>rank == point|>typeof|>rank == 1
    @test point|>pidtype == point|>typeof|>pidtype == PID{Int}
    @test point|>dimension == point|>typeof|>dimension == 2
    @test point|>kind == point|>typeof|>kind == 0
    @test point|>string == "Point(PID(0, 1), [0.0, 0.0], [0.0, 0.0])"
    @test point[1] == point
    @test point|>collect == [point]

    @test isintracell(Point(PID(0, 1), (0.0, 0.0), (0.0, 0.0))) == true
    @test isintracell(Point(PID(0, 1), (0.0, 0.0), (1.0, 0.0))) == false
end

@testset "Bond" begin
    bond = Bond(1, Point(PID(1, 1), (0.0, 0.0), (0.0, 0.0)), Point(PID(1, 2), (0.0, 1.0), (0.0, 1.0)))
    @test bond|>deepcopy == bond
    @test isequal(bond|>deepcopy, bond)
    @test bond|>string == "Bond(1, Point(PID(1, 1), [0.0, 0.0], [0.0, 0.0]), Point(PID(1, 2), [0.0, 1.0], [0.0, 1.0]))"
    @test bond|>reverse == Bond(1, Point(PID(1, 2), (0.0, 1.0), (0.0, 1.0)), Point(PID(1, 1), (0.0, 0.0), (0.0, 0.0)))
    @test bond|>length == bond|>typeof|>length == 2
    @test bond|>eltype == bond|>typeof|>eltype == Point{2, PID{Int}, Float}
    @test bond|>rank == bond|>typeof|>rank == 2
    @test bond|>pidtype == bond|>typeof|>pidtype == PID{Int}
    @test bond|>dimension == bond|>typeof|>dimension == 2
    @test bond|>kind == 1
    @test bond|>rcoord == [0.0, 1.0]
    @test bond|>icoord == [0.0, 1.0]
    @test bond|>isintracell == false
    @test bond[1] == bond.epoint && bond[2] == bond.spoint
    @test bond|>collect == [Point(PID(1, 2), (0.0, 1.0), (0.0, 1.0)), Point(PID(1, 1), (0.0, 0.0), (0.0, 0.0))]
end

@testset "Lattice" begin
    lattice = Lattice("Tuanzi", [Point(PID(1, 1), (0.5, 0.5), (0.0, 0.0))], vectors=[[1.0, 0.0], [0.0, 1.0]], neighbors=1)
    @test lattice|>typeof|>contentnames == (:name, :pids, :rcoords, :icoords, :vectors, :reciprocals, :neighbors)
    @test lattice|>deepcopy == lattice
    @test isequal(lattice|>deepcopy, lattice)
    @test lattice|>string == "Lattice(Tuanzi)\n  with 1 point:\n    Point(PID(1, 1), [0.5, 0.5], [0.0, 0.0])\n  with 2 translation vectors:\n    [1.0, 0.0]\n    [0.0, 1.0]\n  with 1 order of nearest neighbors:\n    1 => 1.0\n"
    @test lattice|>length == 1
    @test lattice|>dimension == lattice|>typeof|>dimension == 2
    @test lattice|>keytype == lattice|>typeof|>keytype == PID{Int}
    @test lattice|>valtype == lattice|>typeof|>valtype == Point{2, PID{Int}, Float}
    @test lattice|>dtype == lattice|>typeof|>dtype == Float
    @test lattice[LatticeIndex{'R'}(1)] == lattice[LatticeIndex{'R'}(PID(1, 1))] == SVector(0.5, 0.5)
    @test lattice[LatticeIndex{'I'}(1)] == lattice[LatticeIndex{'I'}(PID(1, 1))] == SVector(0.0, 0.0)
    @test lattice[LatticeIndex{'P'}(1)] == lattice[LatticeIndex{'P'}(PID(1, 1))] == Point(PID(1, 1), (0.5, 0.5), (0.0, 0.0))
    @test lattice|>nneighbor == 1

    zerothbs = [Point(PID(1, 1), (0.5, 0.5), (0.0, 0.0))]
    insidebs = []
    acrossbs = [Bond(1, Point(PID(1, 1), [0.5, 0.5], [0.0, 0.0]), Point(PID(1, 1), [0.5, -0.5], [0.0, -1.0])),
                Bond(1, Point(PID(1, 1), [0.5, 0.5], [0.0, 0.0]), Point(PID(1, 1), [-0.5, 0.5], [-1.0, 0.0]))
    ]
    @test bonds(lattice, zerothbonds) == zerothbs
    @test bonds(lattice, insidebonds) == insidebs
    @test bonds(lattice, acrossbonds) == acrossbs
    @test setdiff(bonds(lattice), [zerothbs; insidebs; acrossbs])|>length == 0

    lattice = Lattice("SuperTuanzi",
                [Lattice("Tuanzi1", [Point(PID(1, 1), (0.0, 0.0))], neighbors=0),
                Lattice("Tuanzi2", [Point(PID(2, 1), (0.5, 0.5))], neighbors=0),
                ],
                vectors=[[1.0, 0.0], [0.0, 1.0]],
                neighbors=2,
                )
    @test setdiff(bonds(lattice), bonds(lattice, zerothbonds, insidebonds, acrossbonds))|>length == 0
end

@testset "SuperLattice" begin
    lattice = SuperLattice("SuperTuanzi",
                [Lattice("TuanziSys", [Point(PID(1, 1), (0.0, 0.0)), Point(PID(1, 2), (0.5, 0.5))], neighbors=Dict(1=>√2/2)),
                Lattice("TuanziEnv", [Point(PID(2, 1), (-0.05, -0.05)), Point(PID(2, 2), (0.55, 0.55))], neighbors=Dict{Int, Float}()),
                ],
                vectors=[[1.0, 0.0], [0.0, 1.0]],
                neighbors=Dict(1=>√2/2, -1=>√2/20)
                )
    @test lattice|>typeof|>contentnames == (:sublattices, :name, :pids, :rcoords, :icoords, :vectors, :reciprocals, :neighbors)
    @test latticetype(lattice) == latticetype(typeof(lattice)) == Lattice{2, PID{Int}, Float}
    @test setdiff(bonds(lattice), bonds(lattice, zerothbonds, insidebonds, acrossbonds))|>length == 0
    @test setdiff(bonds(lattice, insidebonds), bonds(lattice, intrabonds, interbonds))|>length == 0
end

@testset "Cylinder" begin
    cylinder = Cylinder{PID{String}}("Tuanzi", [0.0 0.0; 0.0 1.0], SVector(1.0, 0.0), vector=[0.0, 2.0], neighbors=1)
    @test cylinder|>typeof|>contentnames == (:block, :translation, :name, :pids, :rcoords, :icoords, :vectors, :reciprocals, :neighbors)
    insert!(cylinder, "A", "B")
    insert!(cylinder, "C3", cut=2, scopes=["C1", "C1", "C2", "C2"])
    @test cylinder.pids == [PID("C1", 1), PID("C3", 2), PID("C3", 1), PID("C1", 2), PID("C2", 1), PID("C2", 2)]
    @test cylinder.rcoords == [-1.0 -1.0 0.0 0.0 1.0 1.0; 0.0 1.0 0.0 1.0 0.0 1.0]
    @test cylinder.icoords == [0.0 0.0 0.0 0.0 0.0 0.0; 0.0 0.0 0.0 0.0 0.0 0.0]
    lattice = cylinder("C1", "C2", "C3")
    @test lattice.pids == [PID("C1", 1), PID("C1", 2), PID("C2", 1), PID("C2", 2), PID("C3", 1), PID("C3", 2)]
    @test lattice.rcoords == [-1.0 -1.0 0.0 0.0 1.0 1.0; 0.0 1.0 0.0 1.0 0.0 1.0]
    @test lattice.icoords == [0.0 0.0 0.0 0.0 0.0 0.0; 0.0 0.0 0.0 0.0 0.0 0.0]
end

@testset "LatticeBonds" begin
    @test eltype(typeof(allbonds)) === AbstractBond
    @test eltype(typeof(zerothbonds)) === Point
    @test eltype(typeof(insidebonds)) === eltype(typeof(acrossbonds)) === Bond
    @test eltype(typeof(intrabonds)) === eltype(typeof(interbonds)) === Bond
    @test eltype(AbstractLattice{2, PID{Int}, Float}, allbonds) == eltype(AbstractLattice{2, PID{Int}, Float}, Val(allbonds)) == AbstractBond{2, PID{Int}, Float}
    @test eltype(Lattice{3, PID{Char}, Float}, zerothbonds) == eltype(Lattice{3, PID{Char}, Float}, Val(zerothbonds)) == Point{3, PID{Char}, Float}
    @test eltype(Lattice{1, PID{String}, Float}, insidebonds) == eltype(Lattice{1, PID{String}, Float}, Val(insidebonds)) == Bond{1, PID{String}, Float}

    @test expand(Lattice, Val((allbonds,))) == Lattice|>latticebondsstructure|>leaves|>Tuple
    @test expand(Lattice, Val((zerothbonds,))) == (zerothbonds,)
    @test expand(Lattice, Val((insidebonds,))) == (insidebonds,)
    @test expand(Lattice, Val((acrossbonds,))) == (acrossbonds,)

    @test expand(SuperLattice, Val((allbonds,))) == SuperLattice|>latticebondsstructure|>leaves|>Tuple
    @test expand(SuperLattice, Val((zerothbonds,))) == (zerothbonds,)
    @test expand(SuperLattice, Val((insidebonds,))) == (intrabonds, interbonds)
    @test expand(SuperLattice, Val((acrossbonds,))) == (acrossbonds,)
    @test expand(SuperLattice, Val((intrabonds,))) == (intrabonds,)
    @test expand(SuperLattice, Val((interbonds,))) == (interbonds,)
end

@testset "Bonds" begin
    lattice = Lattice("Tuanzi", [Point(PID(1, 1), (0.5, 0.5), (0.0, 0.0))], vectors=[[1.0, 0.0], [0.0, 1.0]], neighbors=1)
    bs = Bonds(lattice)
    @test bs|>eltype == bs|>typeof|>eltype == AbstractBond{2, PID{Int}, Float}
    @test bs==Bonds(lattice, allbonds) && isequal(bs, Bonds(lattice, allbonds))
    @test (bs|>length == 3) && (bs|>size == (3,))
    @test summary(bs) == "3-element Bonds"
    @test bs|>bondtypes == bs|>typeof|>bondtypes == (zerothbonds, insidebonds, acrossbonds)
    @test bs|>latticetype == bs|>typeof|>latticetype == Lattice{2, PID{Int}, Float}
    @test bs|>rank == bs|>typeof|>rank == 3
    @test bs|>collect == bonds(lattice)
    @test (bs[1] == bs.bonds[1][1]) && (bs[2] == bs.bonds[3][1]) && (bs[3] == bs.bonds[3][2])
    @test filter(allbonds, bs, :include|>Val) == filter((zerothbonds, insidebonds, acrossbonds), bs, :include|>Val) == bs
    @test filter(zerothbonds, bs, :include|>Val) == filter((insidebonds, acrossbonds), bs, :exclude|>Val) == Bonds(lattice, zerothbonds)
    @test filter(insidebonds, bs, :include|>Val) == filter((zerothbonds, acrossbonds), bs, :exclude|>Val) == Bonds(lattice, insidebonds)
    @test filter(acrossbonds, bs, :include|>Val) == filter((zerothbonds, insidebonds), bs, :exclude|>Val) == Bonds(lattice, acrossbonds)
    emptybs = Bonds{(zerothbonds, insidebonds, acrossbonds), Lattice{2, PID{Int}}}((Point{2, PID{Int}, Float}[], Bond{2, PID{Int}, Float}[], Bond{2, PID{Int}, Float}[]))
    @test filter(bond->(dimension(bond) == 3), bs) == empty!(deepcopy(bs)) == empty(bs) == emptybs
    @test reset!(emptybs, lattice) == bs
end

@testset "BrillouinZone" begin
    @test contentnames(BrillouinZone) == (:reciprocals, :contents)
    bz = BrillouinZone(reciprocals([[1.0, 0.0], [0.0, 1.0]]),
            AbelianNumbers('C', [Momentum2D{10}(j-1, i-1) for (i, j) in Iterators.flatten((Iterators.product(1:10, 1:10),))], collect(0:100), :indptr)
            )
    @test getcontent(bz, Val(:contents)) == (bz.momenta,)
    @test dtype(bz) == dtype(typeof(bz)) == Float
end