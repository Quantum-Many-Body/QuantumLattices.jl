using OffsetArrays: OffsetArray
using LinearAlgebra: cross, det, dot, norm
using Plots: plot, savefig, plot!
using QuantumLattices.Spatials
using QuantumLattices: decompose, dimension, expand, rank, shape
using QuantumLattices.QuantumOperators: scalartype
using QuantumLattices.Toolkit: Float, Segment, contentnames
using Random: seed!
using StaticArrays: SVector

@testset "distance" begin
    @test distance([0.0, 0.0], [1.0, 1.0]) ≈ sqrt(2.0)
    @test distance(OffsetArray([0.0, 0.0], -2:-1), OffsetArray([1.0, 1.0], -2:-1)) ≈ sqrt(2.0)
end

@testset "azimuth" begin
    @test azimuth([+1.0]) ≈ 0.0
    @test azimuth([-1.0]) ≈ pi
    @test azimuth([+1.0, +1.0]) ≈ 1//4*pi
    @test azimuth([-1.0, -1.0]) ≈ 5//4*pi
    @test azimuth([+1.0, +1.0, -2.0]) ≈ 1//4*pi
    @test azimuth([+1.0, -1.0, +2.0]) ≈ 7//4*pi

    @test azimuth(OffsetArray([+1.0], -1:-1)) ≈ 0.0
    @test azimuth(OffsetArray([-1.0], -1:-1)) ≈ pi
    @test azimuth(OffsetArray([+1.0, +1.0], -2:-1)) ≈ 1//4*pi
    @test azimuth(OffsetArray([-1.0, -1.0], -2:-1)) ≈ 5//4*pi
    @test azimuth(OffsetArray([+1.0, +1.0, -2.0], -3:-1)) ≈ 1//4*pi
    @test azimuth(OffsetArray([+1.0, -1.0, +2.0], -3:-1)) ≈ 7//4*pi
end

@testset "azimuthd" begin
    @test azimuthd([+1.0]) ≈ 0.0
    @test azimuthd([-1.0]) ≈ 180.0
    @test azimuthd([+1.0, +1.0]) ≈ 45.0
    @test azimuthd([-1.0, -1.0]) ≈ 225.0
    @test azimuthd([+1.0, +1.0, -2.0]) ≈ 45.0
    @test azimuthd([+1.0, -1.0, +2.0]) ≈ 315.0

    @test azimuthd(OffsetArray([+1.0], -1:-1)) ≈ 0.0
    @test azimuthd(OffsetArray([-1.0], -1:-1)) ≈ 180.0
    @test azimuthd(OffsetArray([+1.0, +1.0], -2:-1)) ≈ 45.0
    @test azimuthd(OffsetArray([-1.0, -1.0], -2:-1)) ≈ 225.0
    @test azimuthd(OffsetArray([+1.0, +1.0, -2.0], -3:-1)) ≈ 45.0
    @test azimuthd(OffsetArray([+1.0, -1.0, +2.0], -3:-1)) ≈ 315.0
end

@testset "polar" begin
    @test polar([1.0, 0.0, 1.0]) ≈ 1//4*pi
    @test polar(OffsetArray([1.0, 0.0, 1.0], -3:-1)) ≈ 1//4*pi
end

@testset "polard" begin
    @test polard([1.0, 0.0, 1.0]) ≈ 45.0
    @test polard(OffsetArray([1.0, 0.0, 1.0], -3:-1)) ≈ 45.0
end

@testset "direction" begin
    @test direction('x') == SVector(1, 0, 0)
    @test direction('y') == SVector(0, 1, 0)
    @test direction('z') == SVector(0, 0, 1)

    @test direction(30, :degree) ≈ direction(pi/6, :radian) ≈ direction([√3, 1]) ≈ SVector(√3/2, 1/2)
    @test direction((90, 30), :degree) ≈ direction((pi/2, pi/6), :radian) ≈ direction([√3, 1, 0]) ≈ SVector(√3/2, 1/2, 0)

    @test direction(OffsetArray([√3, 1], -2:-1)) ≈ OffsetArray([√3/2, 1/2], -2:-1)
end

@testset "volume" begin
    m = rand(3, 3)

    @test volume([[m[1, 1]]]) == volume([m[1, 1]]) == m[1, 1]
    @test volume(m[1:2, 1]) ≈ norm(m[1:2, 1])
    @test volume(m[1:3, 1]) ≈ norm(m[1:3, 1])
    
    @test volume([[m[1, 1]], [m[2, 1]]]) == volume([m[1, 1]], [m[2, 1]]) == 0
    @test volume(m[1:2, 1], m[1:2, 2]) ≈ abs(det(m[1:2, 1:2]))
    @test volume(m[1:3, 1], m[1:3, 2]) ≈ norm(cross(m[1:3, 1], m[1:3, 2]))

    @test volume([m[1, 1]], [m[2, 2]], [m[3, 3]]) == volume(m[1:2, 1], m[1:2, 2], m[1:2, 3]) == 0
    @test volume([m[1:3, 1], m[1:3, 2], m[1:3, 3]]) ≈ volume(m[1:3, 1], m[1:3, 2], m[1:3, 3]) ≈ abs(det(m))
end

@testset "decompose" begin
    a, c = randn(3), randn()
    @test decompose(c*a, a) ≈ c
    @test decompose(c*a, [a]) ≈ [c] 

    a₁, a₂ = randn(2), randn(2)
    c₁, c₂ = randn(2)
    @test all(decompose(c₁*a₁+c₂*a₂, a₁, a₂) .≈ (c₁, c₂))
    @test decompose(c₁*a₁+c₂*a₂, [a₁, a₂]) ≈ [c₁, c₂]

    a₁, a₂ = randn(3), randn(3)
    c₁, c₂ = randn(2)
    @test all(decompose(c₁*a₁+c₂*a₂, a₁, a₂) .≈ (c₁, c₂))
    @test decompose(c₁*a₁+c₂*a₂, [a₁, a₂]) ≈ [c₁, c₂]

    a₁, a₂, a₃ = randn(3), randn(3), randn(3)
    c₁, c₂, c₃ = randn(3)
    @test all(decompose(c₁*a₁+c₂*a₂+c₃*a₃, a₁, a₂, a₃) .≈ (c₁, c₂, c₃))
    @test decompose(c₁*a₁+c₂*a₂+c₃*a₃, [a₁, a₂, a₃]) ≈ [c₁, c₂, c₃]

    σ⁰, σ¹, σ², σ³ = [1 0; 0 1], [0 1; 1 0], [0 -1im; 1im 0], [1 0; 0 -1]
    m = σ⁰ + 2σ¹ + 3σ² + 4σ³
    @test decompose(m, (σ⁰, σ¹, σ², σ³)) == (1, 2, 3, 4)
    @test decompose(m, [σ⁰, σ¹, σ², σ³]) == [1, 2, 3, 4]
end

@testset "isintratriangle" begin
    seed!()
    p₁, p₂, p₃ = rand(2), rand(2), rand(2)
    for (i, p) in enumerate((p₁, p₂, p₃))
        @test isintratriangle(p, p₁, p₂, p₃, vertexes=Tuple(j == i ? true : rand(Bool) for j = 1:3), edges=(rand(Bool), rand(Bool), rand(Bool))) == true
        @test isintratriangle(p, p₁, p₂, p₃, vertexes=Tuple(j == i ? false : rand(Bool) for j = 1:3), edges=(rand(Bool), rand(Bool), rand(Bool))) == false
    end
    for (i, p) in enumerate(((p₁+p₂)/2, (p₂+p₃)/2, (p₃+p₁)/2))
        @test isintratriangle(p, p₁, p₂, p₃, vertexes=(rand(Bool), rand(Bool), rand(Bool)), edges=Tuple(j == i ? true : rand(Bool) for j = 1:3)) == true
        @test isintratriangle(p, p₁, p₂, p₃, vertexes=(rand(Bool), rand(Bool), rand(Bool)), edges=Tuple(j == i ? false : rand(Bool) for j = 1:3)) == false
    end
    @test isintratriangle((p₁+p₂+p₃)/3, p₁, p₂, p₃, vertexes=(rand(Bool), rand(Bool), rand(Bool)), edges=(rand(Bool), rand(Bool), rand(Bool))) == true
    @test isintratriangle(p₁+p₂-2p₃, p₁, p₂, p₃, vertexes=(rand(Bool), rand(Bool), rand(Bool)), edges=(rand(Bool), rand(Bool), rand(Bool))) == false
end

@testset "isonline" begin
    seed!()
    p₁, p₂ = rand(2), rand(2)
    @test isonline(p₁, p₁, p₂, ends = (true, rand(Bool))) == true
    @test isonline(p₁, p₁, p₂, ends = (false, rand(Bool))) == false
    @test isonline(p₂, p₁, p₂, ends = (rand(Bool), true)) == true
    @test isonline(p₂, p₁, p₂, ends = (rand(Bool), false)) == false
    @test isonline((p₁+p₂)/2, p₁, p₂, ends = (rand(Bool), rand(Bool))) == true
    @test isonline(p₁+p₂, p₁, p₂, ends = (rand(Bool), rand(Bool))) == false
end

@testset "isparallel" begin
    v = rand(2)
    @test isparallel([0.0, 0.0], v) == 1
    @test isparallel(v, +rand()*v) == +1
    @test isparallel(v, -rand()*v) == -1
    @test isparallel(v, rand(2)) == 0
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
    m = rand(3, 3)
    @test [dot(a, b) for b in reciprocals([[m[1, 1]]]), a in [[m[1, 1]]]] ≈ [2pi;;]
    @test [dot(a, b) for b in reciprocals([m[1:2, 1]]), a in [m[1:2, 1]]] ≈ [2pi;;]
    @test [dot(a, b) for b in reciprocals([m[1:3, 1]]), a in [m[1:3, 1]]] ≈ [2pi;;]
    @test [dot(a, b) for b in reciprocals([m[1:2, 1], m[1:2, 2]]), a in [m[1:2, 1], m[1:2, 2]]] ≈ [2pi 0; 0 2pi]
    @test [dot(a, b) for b in reciprocals([m[1:3, 1], m[1:3, 2]]), a in [m[1:3, 1], m[1:3, 2]]] ≈ [2pi 0; 0 2pi]
    @test [dot(a, b) for b in reciprocals([m[1:3, 1], m[1:3, 2], m[1:3, 3]]), a in [m[1:3, 1], m[1:3, 2], m[1:3, 3]]] ≈ [2pi 0 0; 0 2pi 0; 0 0 2pi]
end

@testset "rotate" begin
    @test rotate([1.0, 0.0], pi/2) ≈ [0.0, 1.0]
    @test rotate([1.0 2.0; 1.0 2.0], pi/4, axis=([1.0, 1.0], (0.0, 0.0))) ≈ [1.0 1.0; 1.0 1.0+√2.0]
    @test rotate([2.0; 2.0; 2.0;;], pi/12, axis=([0.0, 1.0, 1.0], (pi/2, 0.0))) ≈ [2.0; 1.0+√2/2; 1.0+√6/2;;]
    @test rotate([2.0; 2.0; 2.0;;], pi/12, axis=([1.0, 0.0, 1.0], (pi/2, pi/2))) ≈ [1.0+√6/2; 2.0; 1.0+√2/2;;]
    @test rotate([2.0; 2.0; 2.0;;], pi/12, axis=([1.0, 1.0, 0.0], (0.0, 0.0))) ≈ [1.0+√2/2; 1.0+√6/2; 2.0;;]

    @test rotate(OffsetArray([1.0, 0.0], -2:-1), pi/2) ≈ OffsetArray([0.0, 1.0], -2:-1)
    @test rotate(OffsetArray([1.0 2.0; 1.0 2.0], -2:-1, -4:-3), pi/4, axis=(OffsetArray([1.0, 1.0], -2:-1), (0.0, 0.0))) ≈ OffsetArray([1.0 1.0; 1.0 1.0+√2.0], -2:-1, -4:-3)
    @test rotate(OffsetArray([2.0; 2.0; 2.0;;], -3:-1, -4:-4), pi/12, axis=(OffsetArray([0.0, 1.0, 1.0], -3:-1), (pi/2, 0.0))) ≈ OffsetArray([2.0; 1.0+√2/2; 1.0+√6/2;;], -3:-1, -4:-4)
    @test rotate(OffsetArray([2.0; 2.0; 2.0;;], -3:-1, -4:-4), pi/12, axis=(OffsetArray([1.0, 0.0, 1.0], -3:-1), (pi/2, pi/2))) ≈ OffsetArray([1.0+√6/2; 2.0; 1.0+√2/2;;], -3:-1, -4:-4)
    @test rotate(OffsetArray([2.0; 2.0; 2.0;;], -3:-1, -4:-4), pi/12, axis=(OffsetArray([1.0, 1.0, 0.0], -3:-1), (0.0, 0.0))) ≈ OffsetArray([1.0+√2/2; 1.0+√6/2; 2.0;;], -3:-1, -4:-4)
end

@testset "translate" begin
    @test translate([1.0, 2.0], [1.0, 2.0]) == [2.0, 4.0]
    @test translate([1.0 2.0; 3.0 4.0], [1.0, 2.0]) == [2.0 3.0; 5.0 6.0]
    @test translate(OffsetArray([1.0, 2.0], -2:-1), OffsetArray([1.0, 2.0], -2:-1)) == OffsetArray([2.0, 4.0], -2:-1)
    @test translate(OffsetArray([1.0 2.0; 3.0 4.0], -2:-1, -4:-3), OffsetArray([1.0, 2.0], -2:-1)) == OffsetArray([2.0 3.0; 5.0 6.0], -2:-1, -4:-3)
end

@testset "tile" begin
    @test tile([0.0, 0.0], [[1.0, 0.0], [0.0, 1.0]], ((0.5, 0.5), (-0.5, -0.5))) == [0.5 -0.5; 0.5 -0.5]
    @test tile(OffsetArray([0.0, 0.0], -2:-1), [OffsetArray([1.0, 0.0], -2:-1), OffsetArray([0.0, 1.0], -2:-1)], ((0.5, 0.5), (-0.5, -0.5))) == OffsetArray([0.5 -0.5; 0.5 -0.5], -2:-1, 1:2)
end

@testset "minimumlengths" begin
    @test minimumlengths([0.0, 0.0], [[1.0, 0.0], [0.0, 1.0]], 7) ≈ [0.0, 1.0, √2, 2.0, √5, 2*√2, 3.0, √10]
    @test minimumlengths(OffsetArray([0.0, 0.0], -2:-1), [[1.0, 0.0], [0.0, 1.0]], 7) ≈ [0.0, 1.0, √2, 2.0, √5, 2*√2, 3.0, √10]
end

@testset "Neighbors" begin
    neighbors = Neighbors(0=>0.0, 1=>1.0, 2=>√2)
    @test neighbors == Neighbors([0.0, 1.0, √2])
    @test max(neighbors) == √2
    @test nneighbor(neighbors) == 2
    @test nneighbor(Neighbors(Symbol(0)=>0.0, Symbol(1)=>1.0, Symbol(2)=>√2)) == 3
end

@testset "interlinks" begin
    ps₁ = [0.0 1.0; 0.0 0.0]
    ps₂ = [0.0 1.0; 1.0 1.0]
    @test interlinks(ps₁, ps₂, Neighbors(1=>1.0)) == [(1, 1, 1), (1, 2, 2)]

    ps₁ = OffsetArray([0.0 1.0; 0.0 0.0], -2:-1, -4:-3)
    ps₂ = OffsetArray([0.0 1.0; 1.0 1.0], -2:-1, -6:-5)
    @test interlinks(ps₁, ps₂, Neighbors(1=>1.0)) == [(1, -4, -6), (1, -3, -5)]
end

@testset "Point" begin
    point = Point(1, (0.0, 0.0), (0.0, 0.0))
    @test point == Point(1, [0.0, 0.0], [0.0, 0.0])
    @test isequal(point, Point(1, SVector(0.0, 0.0), SVector(0.0, 0.0)))
    @test point|>string == "Point(1, [0.0, 0.0], [0.0, 0.0])"
    @test point|>dimension == point|>typeof|>dimension == 2
    @test point|>scalartype == point|>typeof|>scalartype == Float
    @test isintracell(Point(1, (0.0, 0.0), (0.0, 0.0))) == true
    @test isintracell(Point(1, (0.0, 0.0), (1.0, 0.0))) == false
end

@testset "Bond" begin
    bond = Bond(1, Point(2, (0.0, 1.0), (0.0, 1.0)), Point(1, (0.0, 0.0), (0.0, 0.0)))
    @test bond|>scalartype == bond|>typeof|>scalartype == Float
    @test bond|>eltype == bond|>typeof|>eltype == Point{2, Float}
    @test bond|>length == 2
    @test bond[begin]==bond[1] && bond[end]==bond[2]
    @test bond|>string == "Bond(1, Point(2, [0.0, 1.0], [0.0, 1.0]), Point(1, [0.0, 0.0], [0.0, 0.0]))"
    @test bond|>dimension == bond|>typeof|>dimension == 2
    @test bond|>reverse == Bond(1, Point(1, (0.0, 0.0), (0.0, 0.0)), Point(2, (0.0, 1.0), (0.0, 1.0)))
    @test bond|>collect == [Point(2, (0.0, 1.0), (0.0, 1.0)), Point(1, (0.0, 0.0), (0.0, 0.0))]
    @test bond|>rcoordinate == [0.0, -1.0]
    @test bond|>icoordinate == [0.0, -1.0]
    @test bond|>isintracell == false
    @test Bond(:, Point(2, [1.0], [0.0]), Point(1, [0.0], [0.0]))|>string == "Bond(:, Point(2, [1.0], [0.0]), Point(1, [0.0], [0.0]))"
end

@testset "Lattice" begin
    lattice = Lattice([0.5, 0.5]; name=:Tuanzi, vectors=[[1.0, 0.0], [0.0, 1.0]])
    @test lattice == Lattice([0.5; 0.5;;]; name=:Tuanzi, vectors=[[1.0, 0.0], [0.0, 1.0]])
    @test lattice|>typeof|>contentnames == (:name, :coordinates, :vectors)
    @test lattice|>deepcopy == lattice
    @test isequal(lattice|>deepcopy, lattice)
    @test lattice|>string == "Lattice(Tuanzi)\n  with 1 point:\n    [0.5, 0.5]\n  with 2 translation vectors:\n    [1.0, 0.0]\n    [0.0, 1.0]\n"
    @test lattice|>length == 1
    @test lattice|>dimension == lattice|>typeof|>dimension == 2
    @test lattice|>scalartype == lattice|>typeof|>scalartype == Float
    @test lattice[1] == SVector(0.5, 0.5)
    @test reciprocals(lattice) == reciprocals(lattice.vectors)
    @test Neighbors(lattice, 1) == Neighbors([0.0, 1.0])
    @test Neighbors(lattice, 2) == Neighbors([0.0, 1.0, √2])

    @test bonds(lattice, 1) == bonds(lattice, Neighbors([0.0, 1.0])) == [
        Bond(Point(1, (0.5, 0.5), (0.0, 0.0))),
        Bond(1, Point(1, (0.5, -0.5), (0.0, -1.0)), Point(1, (0.5, 0.5), (0.0, 0.0))),
        Bond(1, Point(1, (-0.5, 0.5), (-1.0, 0.0)), Point(1, (0.5, 0.5), (0.0, 0.0)))
    ]

    unit = Lattice((0.0, 0.0); name=:Square, vectors=[[1.0, 0.0], [0.0, 1.0]])
    @test Lattice(unit, (2, 3), ('P', 'O')) == Lattice([0.0, 0.0], [1.0, 0.0], [0.0, 1.0], [1.0, 1.0], [0.0, 2.0], [1.0, 2.0]; name=Symbol("Square[0:1](0:2)"), vectors=[[2.0, 0.0]])
    @test Lattice(unit, (2, 3), (:periodic, :open); mode=:center) == Lattice([0.0, -1.0], [1.0, -1.0], [0.0, 0.0], [1.0, 0.0], [0.0, 1.0], [1.0, 1.0]; name=Symbol("Square[0:1](-1:1)"), vectors=[[2.0, 0.0]])
    @test bonds(Lattice(unit, (1, 2)), :) == [
        [Bond(:, Point(1, [0.0, 0.0], [0.0, 0.0]), Point(1, [0.0, 0.0], [0.0, 0.0]))] [Bond(:, Point(1, [0.0, 0.0], [0.0, 0.0]), Point(2, [0.0, 1.0], [0.0, 0.0]))];
        [Bond(:, Point(2, [0.0, 1.0], [0.0, 0.0]), Point(1, [0.0, 0.0], [0.0, 0.0]))] [Bond(:, Point(2, [0.0, 1.0], [0.0, 0.0]), Point(2, [0.0, 1.0], [0.0, 0.0]))]
    ]
end

@testset "plot" begin
    lattice = Lattice((0.0, 0.0); name=:Tuanzi, vectors=[[1.0, 0.0], [0.0, 1.0]])
    savefig(plot(lattice, 2), "Lattice.png")
end

@testset "BrillouinZone" begin
    recipls = [[2.0, 0.0], [0.0, 3.0]]
    bz = BrillouinZone{(2, 4)}(recipls)
    @test bz==deepcopy(bz) && isequal(bz, deepcopy(bz))
    @test bz≠BrillouinZone(recipls, 4) && !isequal(bz, BrillouinZone(recipls, 4))
    @test scalartype(bz) == scalartype(typeof(bz)) == Float
    @test dimension(bz) == dimension(typeof(bz)) == 2
    @test label(bz) == label(typeof(bz)) == :k
    @test label(bz, 1) == label(typeof(bz), 1) == "k₁"
    @test label(bz, 2) == label(typeof(bz), 2) == "k₂"
    @test reciprocals(bz) == recipls
    @test rank(bz) == rank(typeof(bz)) == 2
    @test fractionals(bz) == [[0.0, 0.0], [0.0, 0.25], [0.0, 0.5], [0.0, 0.75], [0.5, 0.0], [0.5, 0.25], [0.5, 0.5], [0.5, 0.75]]
    @test collect(bz) ≈ [expand(bz, fractional) for fractional in fractionals(bz)]
    @test shape(bz) == (0:1, 0:3)
    @test hash(bz) == hash(((SVector(2.0, 0.0), SVector(0.0, 3.0)), (2, 4)))
    @test periods(bz) == periods(typeof(bz)) == (2, 4)
    @test period(bz, 1)==period(typeof(bz), 1)==2 && period(bz, 2)==period(typeof(bz), 2)==4
    @test collect(bz) == [[0.0, 0.0], [0.0, 0.75], [0.0, 1.5], [0.0, 2.25], [1.0, 0.0], [1.0, 0.75], [1.0, 1.5], [1.0, 2.25]]
    @test [convert(CartesianIndex, momentum, bz) for momentum in collect(bz)] == CartesianIndex.([(0, 0), (0, 1), (0, 2), (0, 3), (1, 0), (1, 1), (1, 2), (1, 3)])
    @test range(bz, 1)==range(0.0, 0.5, 2)
    @test range(bz, 2)==range(0.0, 0.75, 4)
    @test volume(bz) == 6.0
    @test iscontinuous(bz) == iscontinuous(typeof(bz)) == false
    @test isdiscrete(bz) == isdiscrete(typeof(bz)) == true
    savefig(plot(bz; fractional=false), "BrillouinZone.png")
    savefig(plot(bz; fractional=true), "BrillouinZone-fractional.png")

    bz = BrillouinZone(recipls, Inf)
    @test iscontinuous(bz) == iscontinuous(typeof(bz)) == true
    @test isdiscrete(bz) == isdiscrete(typeof(bz)) == false

    bz = BrillouinZone{(4, Inf)}(recipls)
    @test iscontinuous(bz) == iscontinuous(typeof(bz)) == false
    @test isdiscrete(bz) == isdiscrete(typeof(bz)) == false

    recipls = [[1.0, 0.0, 0.0]]
    @test BrillouinZone{:q, (10,)}(recipls) == BrillouinZone{:q}(recipls, 10)

    recipls = [[1.0, 0.0, 0.0], [0.0, 1.0, 0.0]]
    @test BrillouinZone{(10, 10)}(recipls) == BrillouinZone(recipls, 10)

    recipls = [[1.0, 0.0, 0.0], [0.0, 1.0, 0.0], [0.0, 0.0, 1.0]]
    @test BrillouinZone{(10, 10, 10)}(recipls) == BrillouinZone(recipls, 10)
end

@testset "ReciprocalZone" begin
    rz = ReciprocalZone([[1.0]], length=10)
    @test rz == ReciprocalZone([[1.0]], -1//2=>1//2; length=10)
    @test rz == ReciprocalZone([[1.0]], (-1//2=>1//2,); length=10)
    @test rz == ReciprocalZone([[1.0]], [-1//2=>1//2]; length=10)
    @test all(shrink(rz, 1:5) .== ReciprocalZone([[1.0]], -0.5=>-0.1, length=5, ends=(true, true)))
    @test fractionals(rz) ≈ [[-0.5], [-0.4], [-0.3], [-0.2], [-0.1], [0.0], [0.1], [0.2], [0.3], [0.4]]
    @test collect(rz) ≈ [expand(rz, fractional) for fractional in fractionals(rz)]

    rz = ReciprocalZone{:q}([[1.0]], length=10)
    @test rz == ReciprocalZone{:q}([[1.0]], -1//2=>1//2; length=10)
    @test rz == ReciprocalZone{:q}([[1.0]], (-1//2=>1//2,); length=10)
    @test rz == ReciprocalZone{:q}([[1.0]], [-1//2=>1//2]; length=10)

    b₁, b₂, b₃ = [1.0, 0.0, 0.0], [0.0, 1.0, 0.0], [0.0, 0.0, 1.0]
    @test volume(ReciprocalZone([b₁], -1=>1; length=10)) == 2
    @test volume(ReciprocalZone([b₁, b₂], -1=>1, -1=>1; length=10)) == 4
    @test volume(ReciprocalZone([b₁, b₂, b₃], -1=>1, -1=>1, -1=>1; length=10)) == 8

    rz = ReciprocalZone([b₁, b₂, b₃], -2=>2, -1=>1, -3=>3; length=10)
    @test range(rz, 1) ≈ collect(Segment(-2, 2, 10))
    @test range(rz, 2) ≈ collect(Segment(-1, 1, 10))
    @test range(rz, 3) ≈ collect(Segment(-3, 3, 10))

    bz = BrillouinZone{:q}([[1.0, 0.0], [0.0, 1.0]], 8)
    rz = ReciprocalZone(bz)
    @test rz == ReciprocalZone{:q}([[1.0, 0.0], [0.0, 1.0]], 0=>1, 0=>1; length=8)
    @test collect(rz) == collect(bz)
    savefig(plot(rz; fractional=false), "ReciprocalZone.png")
    savefig(plot(rz; fractional=true), "ReciprocalZone-fractional.png")
end

@testset "ReciprocalScatter" begin
    b₁, b₂ = [2.0, 0.0, 0.0], [0.0, 2.0, 0.0]
    coordinates = [[0.0, 0.5], [0.25, 0.25], [0.5, 0.0], [0.25, -0.25], [0.0, -0.5], [-0.25, -0.25], [-0.5, 0.0], [-0.25, 0.25]]
    rs = ReciprocalScatter([b₁, b₂], coordinates)
    for (i, fractional) in enumerate(fractionals(rs))
        @test rs[i] == expand(rs, fractional) == b₁*fractional[1] + b₂*fractional[2]
    end
    savefig(plot(rs; fractional=false), "ReciprocalScatter.png")
    savefig(plot(rs; fractional=true), "ReciprocalScatter-fractional.png")

    rs = ReciprocalScatter{:q}([b₁, b₂], [[0.0, 0.0], [0.0, 0.25], [0.0, 0.5], [0.0, 0.75], [0.5, 0.0], [0.5, 0.25], [0.5, 0.5], [0.5, 0.75]])
    @test rs == ReciprocalScatter(BrillouinZone{:q}([b₁, b₂], (2, 4)))
    @test rs == ReciprocalScatter(ReciprocalZone{:q}([b₁, b₂], 0=>1, 0=>1; length=(2, 4), ends=(true, false)))
end

@testset "ReciprocalPath" begin
    @test contentnames(ReciprocalPath) == (:contents, :labels)

    b₁, b₂ = [1.0, 0.0], [0.0, 1.0]
    s₁ = (0.0, 0.0) => (0.5, 0.0)
    s₂ = (0.5, 0.0) => (0.5, 0.5)
    s₃ = (0.5, 0.5) => (0.0, 0.0)

    rp = ReciprocalPath([b₁, b₂], s₁, s₂, s₃)
    @test rp == ReciprocalPath(rp.contents, rp.labels)
    @test all(map((x, y)->isapprox(x, y; atol=10^-12), cumsum([step(rp, i) for i=1:length(rp)-1]), [distance(rp, i) for i=2:length(rp)]))

    positions, labels = ticks(rp)
    @test all(map((x, y)->isapprox(x, y; atol=10^-12), positions, [0.0, 0.5, 1.0, 1.0+sqrt(2)/2]))
    @test labels == ["(0.0, 0.0)", "(0.5, 0.0)", "(0.5, 0.5)", "(0.0, 0.0)"]
    savefig(plot(rp), "ReciprocalPath-1.png")

    rp = ReciprocalPath([b₁, b₂], s₁, s₃; labels=("Γ"=>"X", "M"=>"Γ"))
    positions, labels = ticks(rp)
    @test all(map((x, y)->isapprox(x, y; atol=10^-12), positions, [0.0, 0.5, (1+sqrt(2))/2]))
    @test labels == ["Γ", "X / M", "Γ"]
    savefig(plot(rp), "ReciprocalPath-2.png")

    rp = ReciprocalPath{:q}([b₁, b₂], s₁, s₂, s₃)
    @test rp == ReciprocalPath{:q}(rp.contents, rp.labels)

    rp = ReciprocalPath([b₁+b₂], line"X₂-X₁", length=10)
    @test rp ≈ ReciprocalPath([b₁+b₂], -1//2, 1//2; labels=("X₂", "X₁"), length=10)

    rp = ReciprocalPath([b₁, b₂], rectangle"Γ-X-M-Γ", length=10)
    @test rp ≈ ReciprocalPath([b₁, b₂], (0, 0), (1//2, 0), (1//2, 1//2), (0, 0); labels=("Γ", "X", "M", "Γ"), length=10)

    rp = ReciprocalPath{:q}([b₁, b₂], hexagon"Γ-K-M-Γ", length=10)
    @test rp ≈ ReciprocalPath{:q}([b₁, b₂], (0, 0), (2//3, 1//3), (1//2, 1//2), (0, 0); labels=("Γ", "K", "M", "Γ"), length=10)
end

@testset "selectpath" begin
    bz = BrillouinZone([[1.0, 0.0], [0.0, 1.0]], 4)

    path, indexes = selectpath(bz, (0.0, 0.0), (0.5, 0.0), (0.5, 0.5), (0.0, 0.0); atol=1e-9, rtol=1e-9)
    @test indexes == [1, 5, 9, 10, 11, 6, 1]
    @test path == bz[indexes]

    path, indexes = selectpath(bz, rectangle"Γ-X-M-Γ")
    @test indexes == [1, 5, 9, 10, 11, 6, 1]
    @test path == bz[indexes]
    @test path.labels == ("Γ"=>"X", "X"=>"M", "M"=>"Γ")

    path, indexes = selectpath(bz, (-1.5, 0.5)=>(-3.0, 0.0))
    @test collect(path) == [[-1.5, 0.5], [-2.25, 0.25], [-3.0, 0.0]]
    @test bz[indexes] == [[0.5, 0.5], [0.75, 0.25], [0.0, 0.0]]

    plt = plot()
    plot!(plt, bz)
    plot!(plt, path)
    plot!(plt, map(index->Tuple(bz[index]), indexes), seriestype=:scatter)
    savefig(plt, "PickPoint.png")
end

@testset "ReciprocalCurve" begin
    rc = ReciprocalCurve([[0.0, 0.0], [0.5, 0.0], [0.5, 0.5], [0.0, 0.0]])
    @test rc == ReciprocalCurve([(0.0, 0.0), (0.5, 0.0), (0.5, 0.5), (0.0, 0.0)])
    savefig(plot(rc), "ReciprocalCurve.png")

    rp = ReciprocalPath([[1.0, 0.0], [0.0, 1.0]], (0.0, 0.0)=>(0.5, 0.0), (0.5, 0.0)=>(0.5, 0.5), (0.5, 0.5)=>(0.0, 0.0); length=1)
    @test rc == ReciprocalCurve(rp)
end

@testset "utilities" begin
    dlmsave("path.dlm", rand(100), rand(100, 2))
    dlmsave("heatmap.dlm", rand(100), rand(200), rand(200, 100, 3))

    bz = BrillouinZone([[2pi, 0], [0, 2pi]], 200)
    surface = zeros(Float64, (200, 200))
    for (i, k) in enumerate(bz)
        surface[i] = -imag(1/(0.1im+2cos(k[1])+2cos(k[2])))
    end
    savefig(plot(bz, surface), "SingleSurface.png")
    dlmsave("SingleSurface-1.dlm", bz, surface; fractional=false)
    dlmsave("SingleSurface-2.dlm", bz, surface; fractional=true)

    surfaces = zeros(size(surface)..., 2)
    surfaces[:, :, 1] = surface
    surfaces[:, :, 2] = surface
    savefig(plot(bz, surfaces), "MultiSurfaces.png")
    dlmsave("MultiSurfaces-1.dlm", bz, surfaces; fractional=false)
    dlmsave("MultiSurfaces-2.dlm", bz, surfaces; fractional=true)

    rz = ReciprocalZone([[2pi, 0], [0, 2pi]], -2=>2, -1=>1, length=(400, 200))
    surface = zeros(Float64, (200, 400))
    for (i, k) in enumerate(rz)
        surface[i] = -imag(1/(0.1im+2cos(k[1])+2cos(k[2])))
    end
    savefig(plot(rz, surface), "SingleExtendedSurface.png")
    dlmsave("SingleExtendedSurface-1.dlm", rz, surface; fractional=false)
    dlmsave("SingleExtendedSurface-2.dlm", rz, surface; fractional=true)

    surfaces = zeros(size(surface)..., 2)
    surfaces[:, :, 1] = surface
    surfaces[:, :, 2] = surface
    savefig(plot(rz, surfaces), "MultiExtendedSurfaces.png")
    dlmsave("MultiExtendedSurfaces-1.dlm", rz, surfaces; fractional=false)
    dlmsave("MultiExtendedSurfaces-2.dlm", rz, surfaces; fractional=true)

    coordinates = zeros(SVector{2, Float64}, 400)
    coordinates[1:100] = Segment(SVector(0.0, 1.0), SVector(1.0, 0.0), 100)
    coordinates[101:200] = Segment(SVector(1.0, 0.0), SVector(0.0, -1.0), 100)
    coordinates[201:300] = Segment(SVector(0.0, -1.0), SVector(-1.0, 0.0), 100)
    coordinates[301:400] = Segment(SVector(-1.0, 0.0), SVector(0.0, 1.0), 100)
    weights = zeros(400, 2)
    weights[1:100, 1] = fill(1.0, 100)
    weights[101:200, 2] = fill(1.0, 100)
    weights[201:300, 1] = fill(1.0, 100)
    weights[301:400, 2] = fill(1.0, 100)
    rs = ReciprocalScatter([[2pi, 0], [0, 2pi]], coordinates)
    savefig(plot(rs, weights; fractional=false), "Surface-1.png")
    savefig(plot(rs, weights; fractional=true), "Surface-2.png")
    dlmsave("Surface-1.dlm", rs, weights; fractional=false)
    dlmsave("Surface-2.dlm", rs, weights; fractional=true)

    path = ReciprocalPath([[2pi, 0], [0, 2pi]], rectangle"Γ-X-M-Γ")
    band = map(k->-2cos(k[1])-2cos(k[2]), path)
    bands = [band -band]
    savefig(plot(path, bands), "MultiBands.png")
    dlmsave("MultiBands-1.dlm", path, bands; distance=false)
    dlmsave("MultiBands-2.dlm", path, bands; distance=true)

    weights = zeros(length(band), 2, 2)
    weights[:, 1, 1] = abs2.(band)
    weights[:, 2, 2] = abs2.(band)
    savefig(plot(path, bands, weights; weightmultiplier=1.0, weightwidth=2.0, weightcolors=(:blue, :red), weightlabels=("↑", "↓")), "MultiBandsWithWeights.png")
    dlmsave("MultiBandsWithWeights-1.dlm", path, bands, weights; distance=false)
    dlmsave("MultiBandsWithWeights-2.dlm", path, bands, weights; distance=true)

    energies = LinRange(-6.0, 6.0, 401)
    spectrum = [-imag(1/(e+0.1im-b)) for e in energies, b in band]
    savefig(plot(path, energies, spectrum), "SingleSpectrum.png")
    dlmsave("SingleSpectrum-1.dlm", path, energies, spectrum; distance=false)
    dlmsave("SingleSpectrum-2.dlm", path, energies, spectrum; distance=true)

    spectra = zeros(size(spectrum)..., 2)
    spectra[:, :, 1] = spectrum
    spectra[:, :, 2] = spectrum
    savefig(plot(path, energies, spectra), "MultiSpectra.png")
    dlmsave("MultiSpectra-1.dlm", path, energies, spectra; distance=false)
    dlmsave("MultiSpectra-2.dlm", path, energies, spectra; distance=true)
end
