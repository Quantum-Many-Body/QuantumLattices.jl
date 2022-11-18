using Plots: plot, savefig, plot!
using QuantumLattices.Spatials
using QuantumLattices: decompose, dimension, dtype, expand
using QuantumLattices.QuantumNumbers: AbelianNumbers, Momenta, Momentum₁, Momentum₂, Momentum₃
using QuantumLattices.Toolkit: Float, contentnames, getcontent, shape
using Random: seed!
using StaticArrays: SVector

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
    @test volume([[1]]) == volume([1]) == volume([1, 0]) == volume([1, 0, 0]) == 1
    @test volume([1], [2]) == 0
    @test volume([[1, 0], [0, 1]]) == volume([1, 0], [0, 1]) == volume([1, 0, 0], [0, 1, 0]) == 1
    @test volume([1], [2], [3]) == volume([1, 0], [0, 1], [1, 1]) == 0
    @test volume([[1, 0, 0], [0, 1, 0], [0, 0, 1]]) == volume([1, 0, 0], [0, 1, 0], [0, 0, 1]) == 1
end

@testset "isparallel" begin
    @test isparallel([0.0, 0.0], [1.0, 1.0]) == 1
    @test isparallel([2.0, 2.0], [1.0, 1.0]) == 1
    @test isparallel([-2.0, -2.0], [1.0, 1.0]) == -1
    @test isparallel([1.0, -1.0], [1.0, 1.0]) == 0
end

@testset "isonline" begin
    seed!()
    p₁, p₂ = [0.0, 0.0], [1.0, 1.0]
    @test isonline([0.0, 0.0], p₁, p₂, ends = (true, rand(Bool))) == true
    @test isonline([0.0, 0.0], p₁, p₂, ends = (false, rand(Bool))) == false
    @test isonline([1.0, 1.0], p₁, p₂, ends = (rand(Bool), true)) == true
    @test isonline([1.0, 1.0], p₁, p₂, ends = (rand(Bool), false)) == false
    @test isonline([0.5, 0.5], p₁, p₂, ends = (rand(Bool), rand(Bool))) == true
    @test isonline([1.1, 1.1], p₁, p₂, ends = (rand(Bool), rand(Bool))) == false
end

@testset "decompose" begin
    a, c = randn(3), randn()
    @test all(decompose(c*a, a) .≈ (c,))
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
end

@testset "isintratriangle" begin
    seed!()
    p₁, p₂, p₃ = [-1.0, 0.0], [0.0, 1.0], [1.0, 0.0]
    for (i, p) in enumerate((p₁, p₂, p₃))
        @test isintratriangle(p, p₁, p₂, p₃, vertexes=Tuple(j == i ? true : rand(Bool) for j = 1:3), edges=(rand(Bool), rand(Bool), rand(Bool))) == true
        @test isintratriangle(p, p₁, p₂, p₃, vertexes=Tuple(j == i ? false : rand(Bool) for j = 1:3), edges=(rand(Bool), rand(Bool), rand(Bool))) == false
    end
    for (i, p) in enumerate(((p₁+p₂)/2, (p₂+p₃)/2, (p₃+p₁)/2))
        @test isintratriangle(p, p₁, p₂, p₃, vertexes=(rand(Bool), rand(Bool), rand(Bool)), edges=Tuple(j == i ? true : rand(Bool) for j = 1:3)) == true
        @test isintratriangle(p, p₁, p₂, p₃, vertexes=(rand(Bool), rand(Bool), rand(Bool)), edges=Tuple(j == i ? false : rand(Bool) for j = 1:3)) == false
    end
    @assert isintratriangle([0.0, 0.5], p₁, p₂, p₃, vertexes=(rand(Bool), rand(Bool), rand(Bool)), edges=(rand(Bool), rand(Bool), rand(Bool))) == true
    @assert isintratriangle([0.0, 1.5], p₁, p₂, p₃, vertexes=(rand(Bool), rand(Bool), rand(Bool)), edges=(rand(Bool), rand(Bool), rand(Bool))) == false
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
    @test rotate([1.0, 0.0], pi/2) ≈ [0.0, 1.0]
    @test rotate([1.0 2.0; 1.0 2.0], pi/4, axis=([1.0, 1.0], (0.0, 0.0))) ≈ [1.0 1.0; 1.0 1.0+√2.0]
    @test rotate(reshape([2.0, 2.0, 2.0], 3, 1), pi/12, axis=([0.0, 1.0, 1.0], (pi/2, 0.0))) ≈ reshape([2.0, 1.0+√2/2, 1.0+√6/2], 3, 1)
    @test rotate(reshape([2.0, 2.0, 2.0], 3, 1), pi/12, axis=([1.0, 0.0, 1.0], (pi/2, pi/2))) ≈ reshape([1.0+√6/2, 2.0, 1.0+√2/2], 3, 1)
    @test rotate(reshape([2.0, 2.0, 2.0], 3, 1), pi/12, axis=([1.0, 1.0, 0.0], (0.0, 0.0))) ≈ reshape([1.0+√2/2, 1.0+√6/2, 2.0], 3, 1)
end

@testset "tile" begin
    @test tile(reshape([0.0, 0.0], 2, 1), [[1.0, 0.0], [0.0, 1.0]], ((0.5, 0.5), (-0.5, -0.5))) == [0.5 -0.5; 0.5 -0.5]
end

@testset "minimumlengths" begin
    @test minimumlengths(reshape([0.0, 0.0], 2, 1), [[1.0, 0.0], [0.0, 1.0]], 7) ≈ [0.0, 1.0, √2, 2.0, √5, 2*√2, 3.0, √10]
end

@testset "Neighbors" begin
    neighbors = Neighbors(0=>0.0, 1=>1.0, 2=>√2)
    @test neighbors == Neighbors([0.0, 1.0, √2])
    @test max(neighbors) == √2
    @test nneighbor(neighbors) == 2
    @test nneighbor(Neighbors(Symbol(0)=>0.0, Symbol(1)=>1.0, Symbol(2)=>√2)) == 3
end

@testset "interlinks" begin
    ps1 = [0.0 1.0; 0.0 0.0]
    ps2 = [0.0 1.0; 1.0 1.0]
    neighbors = Neighbors(1=>1.0)
    links = [(1, 1, 1), (1, 2, 2)]
    @test interlinks(ps1, ps2, neighbors) == links
end

@testset "Point" begin
    point = Point(1, (0.0, 0.0), (0.0, 0.0))
    @test point == Point(1, [0.0, 0.0], [0.0, 0.0])
    @test isequal(point, Point(1, SVector(0.0, 0.0), SVector(0.0, 0.0)))
    @test point|>string == "Point(1, [0.0, 0.0], [0.0, 0.0])"
    @test point|>dimension == point|>typeof|>dimension == 2
    @test point|>dtype == point|>typeof|>dtype == Float
    @test isintracell(Point(1, (0.0, 0.0), (0.0, 0.0))) == true
    @test isintracell(Point(1, (0.0, 0.0), (1.0, 0.0))) == false
end

@testset "Bond" begin
    bond = Bond(1, Point(2, (0.0, 1.0), (0.0, 1.0)), Point(1, (0.0, 0.0), (0.0, 0.0)))
    @test bond|>deepcopy == bond
    @test isequal(bond|>deepcopy, bond)
    @test bond|>string == "Bond(1, Point(2, [0.0, 1.0], [0.0, 1.0]), Point(1, [0.0, 0.0], [0.0, 0.0]))"
    @test bond|>dimension == bond|>typeof|>dimension == 2
    @test bond|>dtype == bond|>typeof|>dtype == Float
    @test bond|>length == 2
    @test bond|>eltype == bond|>typeof|>eltype == Point{2, Float}
    @test bond|>reverse == Bond(1, Point(1, (0.0, 0.0), (0.0, 0.0)), Point(2, (0.0, 1.0), (0.0, 1.0)))
    @test bond|>collect == [Point(2, (0.0, 1.0), (0.0, 1.0)), Point(1, (0.0, 0.0), (0.0, 0.0))]
    @test bond|>rcoordinate == [0.0, -1.0]
    @test bond|>icoordinate == [0.0, -1.0]
    @test bond|>isintracell == false
end

@testset "Lattice" begin
    lattice = Lattice([0.5, 0.5]; name=:Tuanzi, vectors=[[1.0, 0.0], [0.0, 1.0]])
    @test lattice|>typeof|>contentnames == (:name, :coordinates, :vectors)
    @test lattice|>deepcopy == lattice
    @test isequal(lattice|>deepcopy, lattice)
    @test lattice|>string == "Lattice(Tuanzi)\n  with 1 point:\n    [0.5, 0.5]\n  with 2 translation vectors:\n    [1.0, 0.0]\n    [0.0, 1.0]\n"
    @test lattice|>length == 1
    @test lattice|>dimension == lattice|>typeof|>dimension == 2
    @test lattice|>dtype == lattice|>typeof|>dtype == Float
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
    @test Lattice(unit, (2, 3), ('P', 'O'); mode=:center) == Lattice([0.0, -1.0], [1.0, -1.0], [0.0, 0.0], [1.0, 0.0], [0.0, 1.0], [1.0, 1.0]; name=Symbol("Square[0:1](-1:1)"), vectors=[[2.0, 0.0]])
end

@testset "plot" begin
    lattice = Lattice((0.0, 0.0); name=:Tuanzi, vectors=[[1.0, 0.0], [0.0, 1.0]])
    savefig(plot(lattice, 2), "Lattice.png")
end

@testset "Momentum" begin
    @test expand(Momentum₁{10}(1), [[1.0, 0.0, 0.0]]) == [0.1, 0.0, 0.0]
    @test expand(Momentum₂{10, 100}(1, 1), [[1.0, 0.0, 0.0], [0.0, 1.0, 0.0]]) == [0.1, 0.01, 0.0]
    @test expand(Momentum₃{10, 100, 1000}(1, 1, 1), [[1.0, 0.0, 0.0], [0.0, 1.0, 0.0], [0.0, 0.0, 1.0]]) == [0.1, 0.01, 0.001]

    @test Momentum₁{10}([0.1, 0.0, 0.0], [[1.0, 0.0, 0.0]]) == Momentum₁{10}(1)
    @test Momentum₂{10, 100}([0.1, 0.01, 0.0], [[1.0, 0.0, 0.0], [0.0, 1.0, 0.0]]) == Momentum₂{10, 100}(1, 1)
    @test Momentum₃{10, 100, 1000}([0.1, 0.01, 0.001],[[1.0, 0.0, 0.0], [0.0, 1.0, 0.0], [0.0, 0.0, 1.0]]) == Momentum₃{10, 100, 1000}(1, 1, 1)
end

@testset "BrillouinZone" begin
    recipls = [[1.0, 0.0], [0.0, 1.0]]
    bz = BrillouinZone(Momentum₂{2, 4}, recipls)
    @test dtype(bz) == dtype(typeof(bz)) == Float
    @test shape(bz) == (0:3, 0:1)
    @test keys(bz) == Momenta(Momentum₂{2, 3})
    @test collect(bz) == [ [0.0, 0.0], [0.0, 0.25], [0.0, 0.5], [0.0, 0.75], [0.5, 0.0], [0.5, 0.25], [0.5, 0.5], [0.5, 0.75]]
    savefig(plot(bz), "BrillouinZone.png")

    recipls = [[1.0, 0.0, 0.0]]
    @test BrillouinZone{:q}(Momentum₁{10}, recipls) == BrillouinZone{:q}(recipls, 10)

    recipls = [[1.0, 0.0, 0.0], [0.0, 1.0, 0.0]]
    @test BrillouinZone(Momentum₂{10, 10}, recipls) == BrillouinZone(recipls, 10)

    recipls = [[1.0, 0.0, 0.0], [0.0, 1.0, 0.0], [0.0, 0.0, 1.0]]
    @test BrillouinZone(Momentum₃{10, 10, 10}, recipls) == BrillouinZone(recipls, 10)
end

@testset "ReciprocalZone" begin
    rz = ReciprocalZone([[1.0]], length=10)
    @test rz == ReciprocalZone([[1.0]], -1//2=>1//2; length=10)
    @test rz == ReciprocalZone([[1.0]], (-1//2=>1//2,); length=10)
    @test rz == ReciprocalZone([[1.0]], [-1//2=>1//2]; length=10)

    rz = ReciprocalZone{:q}([[1.0]], length=10)
    @test rz == ReciprocalZone{:q}([[1.0]], -1//2=>1//2; length=10)
    @test rz == ReciprocalZone{:q}([[1.0]], (-1//2=>1//2,); length=10)
    @test rz == ReciprocalZone{:q}([[1.0]], [-1//2=>1//2]; length=10)

    b₁, b₂, b₃ = [1.0, 0.0, 0.0], [0.0, 1.0, 0.0], [0.0, 0.0, 1.0]
    @test ReciprocalZone([b₁], -1=>1; length=10).volume == 2
    @test ReciprocalZone([b₁, b₂], -1=>1, -1=>1; length=10).volume == 4
    @test ReciprocalZone([b₁, b₂, b₃], -1=>1, -1=>1, -1=>1; length=10).volume == 8

    rz = ReciprocalZone(BrillouinZone{:q}(Momentum₂{8, 8}, [[1.0, 0.0], [0.0, 1.0]]))
    @test rz == ReciprocalZone{:q}([[1.0, 0.0], [0.0, 1.0]], 0=>1, 0=>1; length=8)
    savefig(plot(rz), "ReciprocalZone.png")
end

@testset "ReciprocalPath" begin
    @test contentnames(ReciprocalPath) == (:contents,)

    b₁, b₂ = [1.0, 0.0], [0.0, 1.0]
    s₁ = (0.0, 0.0)=>(0.5, 0.0)
    s₂ = (0.5, 0.0)=>(0.5, 0.5)
    s₃ = (0.5, 0.5)=>(0.0, 0.0)

    rp = ReciprocalPath([b₁, b₂], s₁, s₂, s₃)
    @test getcontent(rp, :contents) == rp.momenta
    @test rp == ReciprocalPath(rp.momenta)
    @test rp == ReciprocalPath([b₁, b₂], (s₁, s₂, s₃))
    savefig(plot(rp), "ReciprocalPath.png")

    rp = ReciprocalPath{:q}([b₁, b₂], s₁, s₂, s₃)
    @test rp == ReciprocalPath{:q}(rp.momenta)
    @test rp == ReciprocalPath{:q}([b₁, b₂], (s₁, s₂, s₃))
    @test rp == ReciprocalPath{:q}([b₁, b₂], [s₁, s₂, s₃])

    rp = ReciprocalPath([b₁+b₂], line"X₂-X₁", length=10)
    @test rp == ReciprocalPath([b₁+b₂], (-1//2,)=>(1//2,), length=10)

    rp = ReciprocalPath([b₁, b₂], rectangle"Γ-X-M-Γ", length=10)
    @test rp == ReciprocalPath([b₁, b₂], (0//1, 0//1)=>(1//2, 0//1), (1//2, 0//1)=>(1//2, 1//2), (1//2, 1//2)=>(0//1, 0//1), length=10)

    rp = ReciprocalPath{:q}([b₁, b₂], hexagon"Γ-K-M-Γ", length=10)
    @test rp == ReciprocalPath{:q}([b₁, b₂], (0//1, 0//1)=>(2//3, 1//3), (2//3, 1//3)=>(1//2, 1//2), (1//2, 1//2)=>(0//1, 0//1), length=10)
end

@testset "selectpath" begin
    b₁, b₂ = [1.0, 0.0], [0.0, 1.0]
    bz = BrillouinZone([b₁, b₂], 4)
    segments = ((0.0, 0.0)=>(0.5, 0.0), (0.5, 0.0)=>(0.5, 0.5), (0.5, 0.5)=>(0.0, 0.0))

    path, indexes = selectpath(bz, segments; ends=(true, true), atol=1e-9, rtol=1e-9)
    @test indexes == [1, 5, 9, 10, 11, 6, 1]
    @test path == bz[indexes]

    path, indexes = selectpath(bz, segments; ends=(true, false), atol=1e-9, rtol=1e-9)
    @test indexes == [1, 5, 9, 10, 11, 6]
    @test path == bz[indexes]

    path, indexes = selectpath(bz, (-1.5, 0.5)=>(-3.0, 0.0))
    @test collect(path) == [[-1.5, 0.5], [-2.25, 0.25], [-3.0, 0.0]]
    @test bz[indexes] == [[0.5, 0.5], [0.75, 0.25], [0.0, 0.0]]

    plt = plot()
    plot!(plt, bz)
    plot!(plt, path)
    plot!(plt, map(index->Tuple(bz[index]), indexes), seriestype=:scatter)
    savefig(plt, "PickPoint.png")
end
