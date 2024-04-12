using LaTeXStrings: latexstring
using QuantumLattices: ⊗, ⋅, expand, permute, rank
using QuantumLattices.DegreesOfFreedom: wildcard, AbstractCompositeIndex, CompositeIID, CompositeIndex, Constraint, Coupling, Diagonal, Hilbert, Index, IIDSpace, MatrixCoupling, iidtype, isdefinite, statistics, @indexes
using QuantumLattices.QuantumOperators: Operator, Operators, latexname, matrix, script
using QuantumLattices.QuantumSystems
using QuantumLattices.Spatials: Bond, Point, azimuthd, rcoordinate
using QuantumLattices.Toolkit: Permutations, shape
using SparseArrays: SparseMatrixCSC
using StaticArrays: SVector

@testset "FID" begin
    fid = FID{:f}(1, 1//2, 1)
    @test string(fid) == "FID{:f}(1, 1//2, 1)"
    @test statistics(fid) == statistics(typeof(fid)) == :f
    @test hash(fid) == hash((:f, 1, 1//2, 1))
    @test fid' == replace(fid, nambu=2)
    @test isequal(fid'', replace(fid, nambu=1))
    @test isannihilation(fid) && isannihilation(Index(1, fid)) && isannihilation(CompositeIndex(Index(1, fid), [0.0], [0.0]))
    @test !iscreation(fid) && !iscreation(Index(1, fid)) && !iscreation(CompositeIndex(Index(1, fid), [0.0], [0.0]))

    fid = FID{:b}(1, -1//2, 2)
    @test string(fid) == "FID{:b}(1, -1//2, 2)"
    @test statistics(fid) == statistics(typeof(fid)) == :b
    @test hash(fid) == hash((:b, 1, -1//2, 2))
    @test !isannihilation(fid) && !isannihilation(Index(1, fid)) && !isannihilation(CompositeIndex(Index(1, fid), [0.0], [0.0]))
    @test iscreation(fid) && iscreation(Index(1, fid)) && iscreation(CompositeIndex(Index(1, fid), [0.0], [0.0]))

    fid = FID(1, :α, :)
    @test fid == FID{wildcard}(1, :α, :)
    @test string(fid) == "FID(1, α, :)"
    @test statistics(fid) == wildcard
    @test hash(fid) == hash((wildcard, 1, :α, :))
    @test !isannihilation(fid) && !isannihilation(Index(1, fid)) && !isannihilation(CompositeIndex(Index(1, fid), [0.0], [0.0]))
    @test !iscreation(fid) && !iscreation(Index(1, fid)) && !iscreation(CompositeIndex(Index(1, fid), [0.0], [0.0]))

    @test FID{:f}(1, 1//2, 1)≠FID{:b}(1, 1//2, 1)
    @test isequal(FID{:f}(1, 1//2, 1), FID{:f}(1, 1//2, 1))
    @test !isequal(FID{:f}(1, 1//2, 1), FID{:b}(1, 1//2, 1))

    @test isdefinite(FID{wildcard, Int, Rational{Int}, Int})
    @test !isdefinite(FID{:f, Symbol, typeof(:), Int})
    @test iidtype(FID, Int, typeof(:), Int) == FID{wildcard, Int, typeof(:), Int}
    @test iidtype(FID{:f}, typeof(:), Symbol, Symbol) == FID{:f, typeof(:),  Symbol, Symbol}
end

@testset "Fock" begin
    @test eltype(Fock) == (FID{S, Int, Rational{Int}, Int} where S)
    fock = Fock{:b}(1, 2)
    @test shape(fock) == (1:1, 1:2, 1:2)
    @test CartesianIndex(FID{:b}(1, -1//2, 1), fock) == CartesianIndex(1, 1, 1)
    @test FID(CartesianIndex(1, 1, 1), fock) == FID{:b}(1, -1//2, 1)
    @test collect(fock) == [FID{:b}(1, -1//2, 1), FID{:b}(1, 1//2, 1), FID{:b}(1, -1//2, 2), FID{:b}(1, 1//2, 2)]
    @test statistics(fock) == statistics(typeof(fock)) == :b
    @test string(fock) == "Fock{:b}(norbital=1, nspin=2)"

    @test summary(Fock{:b}(1, 0)) == "0-element Fock{:b}"
    @test summary(Fock{:f}(1, 1)) == "2-element Fock{:f}"

    @test match(FID{wildcard}, Fock{:f}) == match(FID{wildcard}, Fock{:b}) == true
    @test match(FID{:f}, Fock{:f}) == match(FID{:b}, Fock{:b}) == true
    @test match(FID{:b}, Fock{:f}) == match(FID{:f}, Fock{:b}) == false
end

@testset "Fock latex" begin
    @test script(Val(:site), Index(1, FID{:f}(2, 1//2, 1))) == "1"
    @test script(Val(:orbital), FID{:f}(2, 1//2, 1)) == script(Val(:orbital), Index(1, FID{:f}(2, 1//2, 1))) == "2"
    @test script(Val(:spin), FID{:f}(2, 3//2, 1)) == script(Val(:spin), Index(1, FID{:f}(2, 3//2, 1))) == "3//2"
    @test script(Val(:spinsym), FID{:f}(2, 1//2, 1)) == script(Val(:spinsym), Index(1, FID{:f}(2, 1//2, 1))) == "↑"
    @test script(Val(:spinsym), FID{:f}(2, -1//2, 1)) == script(Val(:spinsym), Index(1, FID{:f}(2, -1//2, 1))) == "↓"
    @test script(Val(:nambu), FID{:f}(2, 3//2, 1)) == script(Val(:nambu), Index(1, FID{:f}(2, 3//2, 1))) == ""
    @test script(Val(:nambu), FID{:f}(2, 3//2, 2)) == script(Val(:nambu), Index(1, FID{:f}(2, 3//2, 2))) == "\\dagger"

    @test latexname(Index{<:Union{Int, Colon}, <:FID{:f}}) == Symbol("Index{Union{Int, Colon}, FID{:f}}")
    @test latexname(AbstractCompositeIndex{Index{<:Union{Int, Colon}, <:FID{:f}}}) == Symbol("AbstractCompositeIndex{Index{Union{Int, Colon}, FID{:f}}}")
    @test latexname(FID{:f}) == Symbol("FID{:f}")

    @test latexname(Index{<:Union{Int, Colon}, <:FID{:b}}) == Symbol("Index{Union{Int, Colon}, FID{:b}}")
    @test latexname(AbstractCompositeIndex{Index{<:Union{Int, Colon}, <:FID{:b}}}) == Symbol("AbstractCompositeIndex{Index{Union{Int, Colon}, FID{:b}}}")
    @test latexname(FID{:b}) == Symbol("FID{:b}")

    @test latexname(Index{<:Union{Int, Colon}, <:FID{wildcard}}) == Symbol("Index{Union{Int, Colon}, FID}")
    @test latexname(AbstractCompositeIndex{Index{<:Union{Int, Colon}, <:FID{wildcard}}}) == Symbol("AbstractCompositeIndex{Index{Union{Int, Colon}, FID}}")
    @test latexname(FID{wildcard}) == Symbol("FID")
end

@testset "angle" begin
    @test angle(CompositeIndex(Index(1, FID{:f}(1, 1//2, 1)), [0.0, 0.0], [1.0, 2.0]), [[1.0, 0.0], [0.0, 1.0]], [0.1, 0.0]) ≈ 2pi*0.1
    @test angle(CompositeIndex(Index(1, FID{:f}(1, 1//2, 2)), [0.0, 0.0], [1.0, 2.0]), [[1.0, 0.0], [0.0, 1.0]], [0.0, 0.2]) ≈ -2pi*0.4
end

@testset "Fock Operator" begin
    id₁ = CompositeIndex(Index(2, FID{:f}(1, -1//2, 2)), SVector(0.5, 0.0), SVector(0.0, 0.0))
    id₂ = CompositeIndex(Index(2, FID{:f}(1, -1//2, 1)), SVector(0.5, 0.0), SVector(0.0, 0.0))
    id₃ = CompositeIndex(Index(1, FID{:f}(1, 1//2, 2)), SVector(0.0, 0.0), SVector(0.0, 0.0))
    id₄ = CompositeIndex(Index(1, FID{:f}(1, 1//2, 1)), SVector(0.0, 0.0), SVector(0.0, 0.0))

    opt = Operator(1.0, id₁, id₂)
    @test opt|>isnormalordered

    opt = Operator(1.0, id₁, id₂, id₃, id₄)
    @test opt|>isnormalordered == false
    @test latexstring(opt) == "c^{\\dagger}_{2, 1, ↓}c^{}_{2, 1, ↓}c^{\\dagger}_{1, 1, ↑}c^{}_{1, 1, ↑}"

    op₁ = Operator(1.5, id₁, id₂)
    op₂ = Operator(2.0, id₂, id₁)
    @test op₁*op₂ == Operator(0.0, id₁, id₂, id₂, id₁)

    op₁ = Operator(1.5, id₁, id₂)
    op₂ = Operator(2.0, id₁, id₂)
    @test op₁*op₂ == Operator(3.0, id₁, id₂, id₁, id₂)

    @test permute(id₁, id₂) == (Operator(1), Operator(-1, id₂, id₁))
    @test permute(id₂, id₁) == (Operator(1), Operator(-1, id₁, id₂))
    @test permute(id₁, id₄) == (Operator(-1, id₄, id₁),)
    @test permute(id₄, id₁) == (Operator(-1, id₁, id₄),)


    id₁ = CompositeIndex(Index(2, FID{:b}(1, -1//2, 2)), SVector(0.5, 0.0), SVector(0.0, 0.0))
    id₂ = CompositeIndex(Index(2, FID{:b}(1, -1//2, 1)), SVector(0.5, 0.0), SVector(0.0, 0.0))
    id₃ = CompositeIndex(Index(1, FID{:b}(1, 1//2, 2)), SVector(0.0, 0.0), SVector(0.0, 0.0))
    id₄ = CompositeIndex(Index(1, FID{:b}(1, 1//2, 1)), SVector(0.0, 0.0), SVector(0.0, 0.0))

    opt = Operator(1.0, id₁, id₂)
    @test latexstring(opt) == "b^{\\dagger}_{2, 1, ↓}b^{}_{2, 1, ↓}"

    @test permute(id₁, id₂) == (Operator(-1), Operator(1, id₂, id₁))
    @test permute(id₂, id₁) == (Operator(+1), Operator(1, id₁, id₂))
    @test permute(id₁, id₄) == (Operator(1, id₄, id₁),)
    @test permute(id₄, id₁) == (Operator(1, id₁, id₄),)
end

@testset "Fock IIDSpace" begin
    @test shape(IIDSpace(FID(:α, :σ, 2), Fock{:f}(3, 2))) == (1:3, 1:2, 2:2)
    @test shape(IIDSpace(FID(2, -1//2, 1), Fock{:b}(3, 2))) == (2:2, 1:1, 1:1)
end

@testset "Fock Coupling" begin
    @test CompositeIID(Coupling(2.0, (1, 2), FID, (1, 2), :, (1, 2))) == CompositeIID(FID(1, :, 1), FID(2, :, 2))
    @test CompositeIID(Coupling(2.0, (1, 1, 1, 1), FID, (1, 2, 3, 4), :, :)) == CompositeIID(FID(1, :, 2), FID(2, :, 1), FID(3, :, 2), FID(4, :, 1))

    @test collect(MatrixCoupling(:, FID, :, :, :)) == [
        Coupling(Index(:, FID(:, :, :)), Index(:, FID(:, :, :)))
    ]
    @test collect(MatrixCoupling((1, 2), FID{:f}, :, σ"y", σ"z")) == [
        Coupling(+1im, Index(1, FID{:f}(:, -1//2, 1)), Index(2, FID{:f}(:, 1//2, 2))),
        Coupling(-1im, Index(1, FID{:f}(:, 1//2, 1)), Index(2, FID{:f}(:, -1//2, 2))),
        Coupling(-1im, Index(1, FID{:f}(:, -1//2, 2)), Index(2, FID{:f}(:, 1//2, 1))),
        Coupling(+1im, Index(1, FID{:f}(:, 1//2, 2)), Index(2, FID{:f}(:, -1//2, 1)))
    ]
    @test collect(MatrixCoupling((1, 2), FID{:b}, σ"x", :, σ"0")) == [
        Coupling(Index(1, FID{:b}(2, :, 1)), Index(2, FID{:b}(1, :, 2))),
        Coupling(Index(1, FID{:b}(1, :, 1)), Index(2, FID{:b}(2, :, 2))),
        Coupling(Index(1, FID{:b}(2, :, 2)), Index(2, FID{:b}(1, :, 1))),
        Coupling(Index(1, FID{:b}(1, :, 2)), Index(2, FID{:b}(2, :, 1)))
    ]
    @test collect(MatrixCoupling(:, FID{wildcard}, σ"+", σ"-", :)) == [
        Coupling(Index(:, FID(1, -1//2, :)), Index(:, FID(2, 1//2, :)))
    ]

    fc = Coupling(2.0, (1, 2), FID, (1, 2), :, (2, 1))
    bond = Bond(1, Point(1, SVector(0.0), SVector(0.0)), Point(2, SVector(0.5), SVector(0.0)))
    hilbert = Hilbert(site=>Fock{:f}(2, 2) for site=1:2)
    ex = expand(Val(:Hopping), fc, bond, hilbert)
    @test collect(ex) == [
        Operator(2.0, CompositeIndex(Index(1, FID{:f}(1, -1//2, 2)), SVector(0.0), SVector(0.0)), CompositeIndex(Index(2, FID{:f}(2, -1//2, 1)), SVector(0.5), SVector(0.0))),
        Operator(2.0, CompositeIndex(Index(1, FID{:f}(1, +1//2, 2)), SVector(0.0), SVector(0.0)), CompositeIndex(Index(2, FID{:f}(2, +1//2, 1)), SVector(0.5), SVector(0.0)))
    ]

    fc = Coupling(2.0, (1, 1, 1, 1), FID, :, (1//2, 1//2, -1//2, -1//2), (2, 1, 2, 1))
    point = Point(1, SVector(0.0), SVector(0.0))
    hilbert = Hilbert(point.site=>Fock{:b}(2, 2))
    ex = expand(Val(:term), fc, Bond(point), hilbert)
    @test collect(ex) == [
        Operator(2.0,
                CompositeIndex(Index(1, FID{:b}(1, +1//2, 2)), SVector(0.0), SVector(0.0)),
                CompositeIndex(Index(1, FID{:b}(1, +1//2, 1)), SVector(0.0), SVector(0.0)),
                CompositeIndex(Index(1, FID{:b}(1, -1//2, 2)), SVector(0.0), SVector(0.0)),
                CompositeIndex(Index(1, FID{:b}(1, -1//2, 1)), SVector(0.0), SVector(0.0))
                ),
        Operator(2.0,
                CompositeIndex(Index(1, FID{:b}(2, +1//2, 2)), SVector(0.0), SVector(0.0)),
                CompositeIndex(Index(1, FID{:b}(2, +1//2, 1)), SVector(0.0), SVector(0.0)),
                CompositeIndex(Index(1, FID{:b}(2, -1//2, 2)), SVector(0.0), SVector(0.0)),
                CompositeIndex(Index(1, FID{:b}(2, -1//2, 1)), SVector(0.0), SVector(0.0))
                )
    ]

    fc = Coupling(2.0, @indexes(Index(1, FID(α, 1//2, 2)), Index(1, FID(α, -1//2, 2)), Index(1, FID(β, -1//2, 1)), Index(1, FID(β, 1//2, 1)); constraint=α<β))
    point = Point(1, SVector(0.5), SVector(0.0))
    hilbert = Hilbert(point.site=>Fock{:f}(3, 2))
    ex = expand(Val(:term), fc, Bond(point), hilbert)
    @test collect(ex) == [
        Operator(2.0,
                CompositeIndex(Index(1, FID{:f}(1, +1//2, 2)), SVector(0.5), SVector(0.0)),
                CompositeIndex(Index(1, FID{:f}(1, -1//2, 2)), SVector(0.5), SVector(0.0)),
                CompositeIndex(Index(1, FID{:f}(2, -1//2, 1)), SVector(0.5), SVector(0.0)),
                CompositeIndex(Index(1, FID{:f}(2, +1//2, 1)), SVector(0.5), SVector(0.0))
                ),
        Operator(2.0,
                CompositeIndex(Index(1, FID{:f}(1, +1//2, 2)), SVector(0.5), SVector(0.0)),
                CompositeIndex(Index(1, FID{:f}(1, -1//2, 2)), SVector(0.5), SVector(0.0)),
                CompositeIndex(Index(1, FID{:f}(3, -1//2, 1)), SVector(0.5), SVector(0.0)),
                CompositeIndex(Index(1, FID{:f}(3, +1//2, 1)), SVector(0.5), SVector(0.0))
                ),
        Operator(2.0,
                CompositeIndex(Index(1, FID{:f}(2, +1//2, 2)), SVector(0.5), SVector(0.0)),
                CompositeIndex(Index(1, FID{:f}(2, -1//2, 2)), SVector(0.5), SVector(0.0)),
                CompositeIndex(Index(1, FID{:f}(3, -1//2, 1)), SVector(0.5), SVector(0.0)),
                CompositeIndex(Index(1, FID{:f}(3, +1//2, 1)), SVector(0.5), SVector(0.0))
                )
    ]

    fc₁ = Coupling(+1.0, (1, 1), FID, :, (+1//2, +1//2), (2, 1))
    fc₂ = Coupling(-1.0, (1, 1), FID, :, (-1//2, -1//2), (2, 1))
    point = Point(1, SVector(0.0), SVector(0.0))
    hilbert = Hilbert(point.site=>Fock{:f}(2, 2))
    ex = expand(Val(:term), fc₁*fc₂, Bond(point), hilbert)
    @test collect(ex) == [
        Operator(-1.0,
                CompositeIndex(Index(1, FID{:f}(1, +1//2, 2)), SVector(0.0), SVector(0.0)),
                CompositeIndex(Index(1, FID{:f}(1, +1//2, 1)), SVector(0.0), SVector(0.0)),
                CompositeIndex(Index(1, FID{:f}(1, -1//2, 2)), SVector(0.0), SVector(0.0)),
                CompositeIndex(Index(1, FID{:f}(1, -1//2, 1)), SVector(0.0), SVector(0.0))
                ),
        Operator(-1.0,
                CompositeIndex(Index(1, FID{:f}(2, +1//2, 2)), SVector(0.0), SVector(0.0)),
                CompositeIndex(Index(1, FID{:f}(2, +1//2, 1)), SVector(0.0), SVector(0.0)),
                CompositeIndex(Index(1, FID{:f}(1, -1//2, 2)), SVector(0.0), SVector(0.0)),
                CompositeIndex(Index(1, FID{:f}(1, -1//2, 1)), SVector(0.0), SVector(0.0))
                ),
        Operator(-1.0,
                CompositeIndex(Index(1, FID{:f}(1, +1//2, 2)), SVector(0.0), SVector(0.0)),
                CompositeIndex(Index(1, FID{:f}(1, +1//2, 1)), SVector(0.0), SVector(0.0)),
                CompositeIndex(Index(1, FID{:f}(2, -1//2, 2)), SVector(0.0), SVector(0.0)),
                CompositeIndex(Index(1, FID{:f}(2, -1//2, 1)), SVector(0.0), SVector(0.0))
                ),
        Operator(-1.0,
                CompositeIndex(Index(1, FID{:f}(2, +1//2, 2)), SVector(0.0), SVector(0.0)),
                CompositeIndex(Index(1, FID{:f}(2, +1//2, 1)), SVector(0.0), SVector(0.0)),
                CompositeIndex(Index(1, FID{:f}(2, -1//2, 2)), SVector(0.0), SVector(0.0)),
                CompositeIndex(Index(1, FID{:f}(2, -1//2, 1)), SVector(0.0), SVector(0.0))
                )
    ]
end

@testset "σ" begin
    σ"0" == SparseMatrixCSC([1 0; 0 1])
    σ"x" == SparseMatrixCSC([0 1; 1 0])
    σ"y" == SparseMatrixCSC([0 -1im; 1im 0])
    σ"z" == SparseMatrixCSC([1 0; 0 -1])
    σ"+" == SparseMatrixCSC([0 1; 0 0])
    σ"-" == SparseMatrixCSC([0 0; 1 0])
    σ"11" == SparseMatrixCSC([1 0; 0 0])
    σ"22" == SparseMatrixCSC([0 0; 0 1])
end

@testset "L" begin
    L"x" == SparseMatrixCSC([0 0 0; 0 0 1im; 0 -1im 0])
    L"y" == SparseMatrixCSC([0 0 -1im; 0 0 0; 1im 0 0])
    L"z" == SparseMatrixCSC([0 1im 0; -1im 0 0; 0 0 0])
end

@testset "Onsite" begin
    point = Point(1, (0.5, 0.5), (0.0, 0.0))
    bond = Bond(point)
    hilbert = Hilbert(point.site=>Fock{:f}(2, 2))

    term = Onsite(:mu, 1.5, MatrixCoupling(:, FID, σ"z", σ"x", :))
    operators = Operators(
        Operator(-0.75, CompositeIndex(Index(1, FID{:f}(2, +1//2, 2)), [0.5, 0.5], [0.0, 0.0]), CompositeIndex(Index(1, FID{:f}(2, -1//2, 1)), [0.5, 0.5], [0.0, 0.0])),
        Operator(+0.75, CompositeIndex(Index(1, FID{:f}(1, -1//2, 2)), [0.5, 0.5], [0.0, 0.0]), CompositeIndex(Index(1, FID{:f}(1, +1//2, 1)), [0.5, 0.5], [0.0, 0.0])),
        Operator(-0.75, CompositeIndex(Index(1, FID{:f}(2, -1//2, 2)), [0.5, 0.5], [0.0, 0.0]), CompositeIndex(Index(1, FID{:f}(2, +1//2, 1)), [0.5, 0.5], [0.0, 0.0])),
        Operator(+0.75, CompositeIndex(Index(1, FID{:f}(1, +1//2, 2)), [0.5, 0.5], [0.0, 0.0]), CompositeIndex(Index(1, FID{:f}(1, -1//2, 1)), [0.5, 0.5], [0.0, 0.0]))
    )
    @test expand(term, bond, hilbert, half=true) == operators
    @test expand(term, bond, hilbert, half=false) == operators*2

    term = Onsite(:mu, 1.5, MatrixCoupling(:, FID, σ"z", σ"z", :))
    operators = Operators(
        Operator(+0.75, CompositeIndex(Index(1, FID{:f}(2, -1//2, 2)), [0.5, 0.5], [0.0, 0.0]), CompositeIndex(Index(1, FID{:f}(2, -1//2, 1)), [0.5, 0.5], [0.0, 0.0])),
        Operator(+0.75, CompositeIndex(Index(1, FID{:f}(1, +1//2, 2)), [0.5, 0.5], [0.0, 0.0]), CompositeIndex(Index(1, FID{:f}(1, +1//2, 1)), [0.5, 0.5], [0.0, 0.0])),
        Operator(-0.75, CompositeIndex(Index(1, FID{:f}(2, +1//2, 2)), [0.5, 0.5], [0.0, 0.0]), CompositeIndex(Index(1, FID{:f}(2, +1//2, 1)), [0.5, 0.5], [0.0, 0.0])),
        Operator(-0.75, CompositeIndex(Index(1, FID{:f}(1, -1//2, 2)), [0.5, 0.5], [0.0, 0.0]), CompositeIndex(Index(1, FID{:f}(1, -1//2, 1)), [0.5, 0.5], [0.0, 0.0]))
    )
    @test expand(term, bond, hilbert, half=true) == operators
    @test expand(term, bond, hilbert, half=false) == operators*2
end

@testset "Hopping" begin
    bond = Bond(1, Point(2, (0.0, 0.0), (0.0, 0.0)), Point(1, (0.5, 0.5), (0.0, 0.0)))
    hilbert = Hilbert(site=>Fock{:f}(2, 2) for site=1:2)
    term = Hopping(:t, 1.5, 1)
    operators = Operators(
        Operator(1.5, CompositeIndex(Index(2, FID{:f}(2, +1//2, 2)), [0.0, 0.0], [0.0, 0.0]), CompositeIndex(Index(1, FID{:f}(2, +1//2, 1)), [0.5, 0.5], [0.0, 0.0])),
        Operator(1.5, CompositeIndex(Index(2, FID{:f}(2, -1//2, 2)), [0.0, 0.0], [0.0, 0.0]), CompositeIndex(Index(1, FID{:f}(2, -1//2, 1)), [0.5, 0.5], [0.0, 0.0])),
        Operator(1.5, CompositeIndex(Index(2, FID{:f}(1, -1//2, 2)), [0.0, 0.0], [0.0, 0.0]), CompositeIndex(Index(1, FID{:f}(1, -1//2, 1)), [0.5, 0.5], [0.0, 0.0])),
        Operator(1.5, CompositeIndex(Index(2, FID{:f}(1, +1//2, 2)), [0.0, 0.0], [0.0, 0.0]), CompositeIndex(Index(1, FID{:f}(1, +1//2, 1)), [0.5, 0.5], [0.0, 0.0]))
    )
    @test expand(term, bond, hilbert, half=true) == operators
    @test expand(term, bond, hilbert, half=false) == operators+operators'
end

@testset "Pairing" begin
    bond = Bond(1, Point(2, (0.0, 0.0), (0.0, 0.0)), Point(1, (0.5, 0.5), (0.0, 0.0)))
    hilbert = Hilbert(site=>Fock{:f}(1, 1) for site=1:2)
    term = Pairing(:Δ, 1.5, 1, Coupling((1, 2), FID, :, (0, 0), :); amplitude=bond->(bond|>rcoordinate|>azimuthd ≈ 45 ? 1 : -1))
    operators = Operators(
        Operator(+1.5, CompositeIndex(Index(2, FID{:f}(1, 0, 1)), [0.0, 0.0], [0.0, 0.0]), CompositeIndex(Index(1, FID{:f}(1, 0, 1)), [0.5, 0.5], [0.0, 0.0])),
        Operator(-1.5, CompositeIndex(Index(1, FID{:f}(1, 0, 1)), [0.5, 0.5], [0.0, 0.0]), CompositeIndex(Index(2, FID{:f}(1, 0, 1)), [0.0, 0.0], [0.0, 0.0]))
    )
    @test expand(term, bond, hilbert, half=true) == operators
    @test expand(term, bond, hilbert, half=false) == operators+operators'

    point = Point(1, (0.5, 0.5), (0.0, 0.0))
    hilbert = Hilbert(point.site=>Fock{:f}(1, 2))
    term = Pairing(:Δ, 1.5, 0, MatrixCoupling(:, FID, :, [0 -1; 1 0], :))
    operators = Operators(
        Operator(-1.5, CompositeIndex(Index(1, FID{:f}(1, +1//2, 1)), [0.5, 0.5], [0.0, 0.0]), CompositeIndex(Index(1, FID{:f}(1, -1//2, 1)), [0.5, 0.5], [0.0, 0.0])),
        Operator(+1.5, CompositeIndex(Index(1, FID{:f}(1, -1//2, 1)), [0.5, 0.5], [0.0, 0.0]), CompositeIndex(Index(1, FID{:f}(1, +1//2, 1)), [0.5, 0.5], [0.0, 0.0]))
    )
    @test expand(term, Bond(point), hilbert, half=true) == operators
    @test expand(term, Bond(point), hilbert, half=false) == operators+operators'
end

@testset "Hubbard" begin
    point = Point(1, (0.5, 0.5), (0.0, 0.0))
    bond = Bond(point)
    hilbert = Hilbert(point.site=>Fock{:f}(2, 2))
    term = Hubbard(:H, 2.5)
    operators = Operators(
        Operator(1.25,
            CompositeIndex(Index(1, FID{:f}(1, +1//2, 2)), [0.5, 0.5], [0.0, 0.0]),
            CompositeIndex(Index(1, FID{:f}(1, +1//2, 1)), [0.5, 0.5], [0.0, 0.0]),
            CompositeIndex(Index(1, FID{:f}(1, -1//2, 2)), [0.5, 0.5], [0.0, 0.0]),
            CompositeIndex(Index(1, FID{:f}(1, -1//2, 1)), [0.5, 0.5], [0.0, 0.0])
            ),
        Operator(1.25,
            CompositeIndex(Index(1, FID{:f}(2, +1//2, 2)), [0.5, 0.5], [0.0, 0.0]),
            CompositeIndex(Index(1, FID{:f}(2, +1//2, 1)), [0.5, 0.5], [0.0, 0.0]),
            CompositeIndex(Index(1, FID{:f}(2, -1//2, 2)), [0.5, 0.5], [0.0, 0.0]),
            CompositeIndex(Index(1, FID{:f}(2, -1//2, 1)), [0.5, 0.5], [0.0, 0.0])
            )
    )
    @test expand(term, bond, hilbert, half=true) == operators
    @test expand(term, bond, hilbert, half=false) == operators*2
end

@testset "InterOrbitalInterSpin" begin
    point = Point(1, (0.5, 0.5), (0.0, 0.0))
    bond = Bond(point)
    hilbert = Hilbert(point.site=>Fock{:f}(2, 2))
    term = InterOrbitalInterSpin(:H, 2.5)
    operators = Operators(
        Operator(1.25,
            CompositeIndex(Index(1, FID{:f}(1, +1//2, 2)), [0.5, 0.5], [0.0, 0.0]),
            CompositeIndex(Index(1, FID{:f}(1, +1//2, 1)), [0.5, 0.5], [0.0, 0.0]),
            CompositeIndex(Index(1, FID{:f}(2, -1//2, 2)), [0.5, 0.5], [0.0, 0.0]),
            CompositeIndex(Index(1, FID{:f}(2, -1//2, 1)), [0.5, 0.5], [0.0, 0.0])
            ),
        Operator(1.25,
            CompositeIndex(Index(1, FID{:f}(1, -1//2, 2)), [0.5, 0.5], [0.0, 0.0]),
            CompositeIndex(Index(1, FID{:f}(1, -1//2, 1)), [0.5, 0.5], [0.0, 0.0]),
            CompositeIndex(Index(1, FID{:f}(2, +1//2, 2)), [0.5, 0.5], [0.0, 0.0]),
            CompositeIndex(Index(1, FID{:f}(2, +1//2, 1)), [0.5, 0.5], [0.0, 0.0])
            )
    )
    @test expand(term, bond, hilbert, half=true) == operators
    @test expand(term, bond, hilbert, half=false) == operators*2
end

@testset "InterOrbitalIntraSpin" begin
    point = Point(1, (0.5, 0.5), (0.0, 0.0))
    bond = Bond(point)
    hilbert = Hilbert(point.site=>Fock{:f}(2, 2))
    term = InterOrbitalIntraSpin(:H, 2.5)
    operators = Operators(
        Operator(1.25,
            CompositeIndex(Index(1, FID{:f}(1, -1//2, 2)), [0.5, 0.5], [0.0, 0.0]),
            CompositeIndex(Index(1, FID{:f}(1, -1//2, 1)), [0.5, 0.5], [0.0, 0.0]),
            CompositeIndex(Index(1, FID{:f}(2, -1//2, 2)), [0.5, 0.5], [0.0, 0.0]),
            CompositeIndex(Index(1, FID{:f}(2, -1//2, 1)), [0.5, 0.5], [0.0, 0.0])
            ),
        Operator(1.25,
            CompositeIndex(Index(1, FID{:f}(1, 1//2, 2)), [0.5, 0.5], [0.0, 0.0]),
            CompositeIndex(Index(1, FID{:f}(1, 1//2, 1)), [0.5, 0.5], [0.0, 0.0]),
            CompositeIndex(Index(1, FID{:f}(2, 1//2, 2)), [0.5, 0.5], [0.0, 0.0]),
            CompositeIndex(Index(1, FID{:f}(2, 1//2, 1)), [0.5, 0.5], [0.0, 0.0])
            )
    )
    @test expand(term, bond, hilbert, half=true) == operators
    @test expand(term, bond, hilbert, half=false) == operators*2
end

@testset "SpinFlip" begin
    point = Point(1, (0.5, 0.5), (0.0, 0.0))
    bond = Bond(point)
    hilbert = Hilbert(point.site=>Fock{:f}(2, 2))
    term = SpinFlip(:H, 2.5)
    operators = Operators(
        Operator(2.5,
            CompositeIndex(Index(1, FID{:f}(1, +1//2, 2)), [0.5, 0.5], [0.0, 0.0]),
            CompositeIndex(Index(1, FID{:f}(2, -1//2, 2)), [0.5, 0.5], [0.0, 0.0]),
            CompositeIndex(Index(1, FID{:f}(1, -1//2, 1)), [0.5, 0.5], [0.0, 0.0]),
            CompositeIndex(Index(1, FID{:f}(2, +1//2, 1)), [0.5, 0.5], [0.0, 0.0])
            )
    )
    @test expand(term, bond, hilbert, half=true) == operators
    @test expand(term, bond, hilbert, half=false) == operators+operators'
end

@testset "PairHopping" begin
    point = Point(1, (0.5, 0.5), (0.0, 0.0))
    bond = Bond(point)
    hilbert = Hilbert(point.site=>Fock{:f}(2, 2))
    term = PairHopping(:H, 2.5)
    operators = Operators(
        Operator(2.5,
            CompositeIndex(Index(1, FID{:f}(1, +1//2, 2)), [0.5, 0.5], [0.0, 0.0]),
            CompositeIndex(Index(1, FID{:f}(1, -1//2, 2)), [0.5, 0.5], [0.0, 0.0]),
            CompositeIndex(Index(1, FID{:f}(2, -1//2, 1)), [0.5, 0.5], [0.0, 0.0]),
            CompositeIndex(Index(1, FID{:f}(2, +1//2, 1)), [0.5, 0.5], [0.0, 0.0])
            )
    )
    @test expand(term, bond, hilbert, half=true) == operators
    @test expand(term, bond, hilbert, half=false) == operators+operators'
end

@testset "Coulomb" begin
    bond = Bond(1, Point(2, (0.0, 0.0), (0.0, 0.0)), Point(1, (0.5, 0.5), (0.0, 0.0)))
    hilbert = Hilbert(site=>Fock{:f}(1, 2) for site=1:2)

    term = Coulomb(:V, 2.5, 1, MatrixCoupling(:, FID, :, σ"z", :)^2)
    operators = Operators(
        Operator(-1.25,
            CompositeIndex(Index(2, FID{:f}(1, -1//2, 2)), [0.0, 0.0], [0.0, 0.0]),
            CompositeIndex(Index(2, FID{:f}(1, -1//2, 1)), [0.0, 0.0], [0.0, 0.0]),
            CompositeIndex(Index(1, FID{:f}(1, +1//2, 2)), [0.5, 0.5], [0.0, 0.0]),
            CompositeIndex(Index(1, FID{:f}(1, +1//2, 1)), [0.5, 0.5], [0.0, 0.0])
            ),
        Operator(+1.25,
            CompositeIndex(Index(2, FID{:f}(1, -1//2, 2)), [0.0, 0.0], [0.0, 0.0]),
            CompositeIndex(Index(2, FID{:f}(1, -1//2, 1)), [0.0, 0.0], [0.0, 0.0]),
            CompositeIndex(Index(1, FID{:f}(1, -1//2, 2)), [0.5, 0.5], [0.0, 0.0]),
            CompositeIndex(Index(1, FID{:f}(1, -1//2, 1)), [0.5, 0.5], [0.0, 0.0])
            ),
        Operator(-1.25,
            CompositeIndex(Index(2, FID{:f}(1, +1//2, 2)), [0.0, 0.0], [0.0, 0.0]),
            CompositeIndex(Index(2, FID{:f}(1, +1//2, 1)), [0.0, 0.0], [0.0, 0.0]),
            CompositeIndex(Index(1, FID{:f}(1, -1//2, 2)), [0.5, 0.5], [0.0, 0.0]),
            CompositeIndex(Index(1, FID{:f}(1, -1//2, 1)), [0.5, 0.5], [0.0, 0.0])
            ),
        Operator(+1.25,
            CompositeIndex(Index(2, FID{:f}(1, +1//2, 2)), [0.0, 0.0], [0.0, 0.0]),
            CompositeIndex(Index(2, FID{:f}(1, +1//2, 1)), [0.0, 0.0], [0.0, 0.0]),
            CompositeIndex(Index(1, FID{:f}(1, +1//2, 2)), [0.5, 0.5], [0.0, 0.0]),
            CompositeIndex(Index(1, FID{:f}(1, +1//2, 1)), [0.5, 0.5], [0.0, 0.0])
            )
    )
    @test expand(term, bond, hilbert, half=true) == operators
    @test expand(term, bond, hilbert, half=false) == operators*2

    term = Coulomb(:V, 2.5, 1, MatrixCoupling(:, FID, :, σ"x", :)*MatrixCoupling(:, FID, :, σ"z", :))
    operators = Operators(
        Operator(-1.25,
            CompositeIndex(Index(2, FID{:f}(1, +1//2, 2)), [0.0, 0.0], [0.0, 0.0]),
            CompositeIndex(Index(2, FID{:f}(1, -1//2, 1)), [0.0, 0.0], [0.0, 0.0]),
            CompositeIndex(Index(1, FID{:f}(1, -1//2, 2)), [0.5, 0.5], [0.0, 0.0]),
            CompositeIndex(Index(1, FID{:f}(1, -1//2, 1)), [0.5, 0.5], [0.0, 0.0])
            ),
        Operator(+1.25,
            CompositeIndex(Index(2, FID{:f}(1, -1//2, 2)), [0.0, 0.0], [0.0, 0.0]),
            CompositeIndex(Index(2, FID{:f}(1, +1//2, 1)), [0.0, 0.0], [0.0, 0.0]),
            CompositeIndex(Index(1, FID{:f}(1, +1//2, 2)), [0.5, 0.5], [0.0, 0.0]),
            CompositeIndex(Index(1, FID{:f}(1, +1//2, 1)), [0.5, 0.5], [0.0, 0.0])
            ),
        Operator(+1.25,
            CompositeIndex(Index(2, FID{:f}(1, +1//2, 2)), [0.0, 0.0], [0.0, 0.0]),
            CompositeIndex(Index(2, FID{:f}(1, -1//2, 1)), [0.0, 0.0], [0.0, 0.0]),
            CompositeIndex(Index(1, FID{:f}(1, +1//2, 2)), [0.5, 0.5], [0.0, 0.0]),
            CompositeIndex(Index(1, FID{:f}(1, +1//2, 1)), [0.5, 0.5], [0.0, 0.0])
            ),
        Operator(-1.25,
            CompositeIndex(Index(2, FID{:f}(1, -1//2, 2)), [0.0, 0.0], [0.0, 0.0]),
            CompositeIndex(Index(2, FID{:f}(1, +1//2, 1)), [0.0, 0.0], [0.0, 0.0]),
            CompositeIndex(Index(1, FID{:f}(1, -1//2, 2)), [0.5, 0.5], [0.0, 0.0]),
            CompositeIndex(Index(1, FID{:f}(1, -1//2, 1)), [0.5, 0.5], [0.0, 0.0])
            )
    )
    @test expand(term, bond, hilbert, half=true) == operators
    @test expand(term, bond, hilbert, half=false) == operators*2
end

@testset "SID" begin
    @test SID{1//2}('z')' == SID{1//2}('z')
    @test SID{1//2}('x')' == SID{1//2}('x')
    @test SID{1//2}('y')' == SID{1//2}('y')
    @test SID{1//2}('+')' == SID{1//2}('-')
    @test SID{1//2}('-')' == SID{1//2}('+')

    sid = SID{3//2}('x')
    @test string(sid) == "SID{3//2}('x')"
    @test replace(sid, tag='z') == SID{3//2}('z')
    @test hash(sid) == hash((3//2, 'x'))
    @test totalspin(sid) == totalspin(typeof(sid)) == 3//2
    @test totalspin(Index(:, sid)) == totalspin(CompositeIndex(Index(:, sid), [0], [0])) == 3//2
    @test statistics(sid) == statistics(typeof(sid)) == :b

    sid = SID{wildcard}('z')
    @test sid == SID('z')
    @test string(sid) == "SID('z')"

    @test SID{1//2}('z')≠SID{3//2}('z')
    @test isequal(SID{1//2}('z'), SID{1//2}('z'))
    @test !isequal(SID{1//2}('z'), SID{3//2}('z'))

    @test isdefinite(SID{wildcard, Char})
    @test !isdefinite(SID{1//2, Symbol})
    @test !isdefinite(SID{1, typeof(:)})
    @test iidtype(SID, Char) == SID{wildcard, Char}
    @test iidtype(SID{1//2}, Symbol) == SID{1//2, Symbol}
end

@testset "matrix" begin
    @test isapprox(matrix(SID{1//2}('z')), [[-0.5, 0.0] [0.0, 0.5]])
    @test isapprox(matrix(SID{1//2}('x')), [[0.0, 0.5] [0.5, 0.0]])
    @test isapprox(matrix(SID{1//2}('y')), [[0.0, -0.5im] [0.5im, 0.0]])
    @test isapprox(matrix(SID{1//2}('+')), [[0.0, 1.0] [0.0, 0.0]])
    @test isapprox(matrix(SID{1//2}('-')), [[0.0, 0.0] [1.0, 0.0]])

    @test isapprox(matrix(Index(:, SID{1//2}('z'))), [[-0.5, 0.0] [0.0, 0.5]])
    @test isapprox(matrix(CompositeIndex(Index(:, SID{1//2}('z')), [0], [0])), [[-0.5, 0.0] [0.0, 0.5]])

    @test isapprox(matrix(SID{1}('z')), [[-1.0, 0.0, 0.0] [0.0, 0.0, 0.0] [0.0, 0.0, 1.0]])
    @test isapprox(matrix(SID{1}('x')), [[0.0, √2/2, 0.0] [√2/2, 0.0, √2/2] [0.0, √2/2, 0.0]])
    @test isapprox(matrix(SID{1}('y')), [[0.0, -√2im/2, 0.0] [√2im/2, 0.0, -√2im/2] [0.0, √2im/2, 0.0]])
    @test isapprox(matrix(SID{1}('+')), [[0.0, √2, 0.0] [0.0, 0.0, √2] [0.0, 0.0, 0.0]])
    @test isapprox(matrix(SID{1}('-')), [[0.0, 0.0, 0.0] [√2, 0.0, 0.0] [0.0, √2, 0.0]])
end

@testset "Spin" begin
    @test eltype(Spin) == (SID{S, Char} where S)
    spin = Spin{1}()
    @test shape(spin) == (1:5,)
    @test CartesianIndex(SID{1}('z'), spin) == CartesianIndex(3)
    @test SID(CartesianIndex(1), spin) == SID{1}('x')
    @test summary(spin) == "5-element Spin{1}"
    @test string(spin) == "Spin{1}()"
    @test totalspin(spin) == totalspin(typeof(spin)) == 1
    @test collect(spin) == [SID{1}('x'), SID{1}('y'), SID{1}('z'), SID{1}('+'), SID{1}('-')]

    @test match(SID{wildcard}, Spin{1//2}) == true
    @test match(SID{1//2}, Spin{1//2}) == true
    @test match(SID{1//2}, Spin{1}) == match(SID{1}, Spin{1//2}) == false
end

@testset "latex" begin
    index = Index(1, SID{1//2}('z'))
    @test script(Val(:site), index) == "1"
    @test script(Val(:tag), index.iid) == script(Val(:tag), index) == "z"

    @test latexname(Index{<:Union{Int, Colon}, <:SID}) == Symbol("Index{Union{Int, Colon}, SID}")
    @test latexname(AbstractCompositeIndex{<:Index{<:Union{Int, Colon}, <:SID}}) == Symbol("AbstractCompositeIndex{Index{Union{Int, Colon}, SID}}")
    @test latexname(SID) == Symbol("SID")
end

@testset "SpinOperator" begin
    opt = Operator(
        1.0,
        CompositeIndex(Index(1, SID{1//2}('+')), [0.0, 0.0], [0.0, 0.0]),
        CompositeIndex(Index(1, SID{1//2}('-')), [0.0, 0.0], [0.0, 0.0])
    )
    @test opt' == Operator(
        1.0,
        CompositeIndex(Index(1, SID{1//2}('+')), [0.0, 0.0], [0.0, 0.0]),
        CompositeIndex(Index(1, SID{1//2}('-')), [0.0, 0.0], [0.0, 0.0])
    )
    @test latexstring(opt) == "S^{+}_{1}S^{-}_{1}"
end

@testset "permute" begin
    soptrep(opt::Operator) = opt.value * prod([matrix(opt.id[i].index.iid) for i = 1:rank(opt)])
    for S in (1//2, 1, 3//2)
        indexes = [CompositeIndex(Index(1, SID{S}(tag)), [0.0, 0.0], [0.0, 0.0]) for tag in ('x', 'y', 'z', '+', '-')]
        for (id₁, id₂) in Permutations{2}(indexes)
            left = soptrep(Operator(1, id₁, id₂))
            right = sum([soptrep(opt) for opt in permute(id₁, id₂)])
            @test isapprox(left, right)
        end
    end
    id₁ = CompositeIndex(Index(1, SID{1//2}('z')), [0.0, 0.0], [0.0, 0.0])
    id₂ = CompositeIndex(Index(2, SID{1//2}('z')), [0.0, 0.0], [0.0, 0.0])
    @test permute(id₁, id₂) == (Operator(1, id₂, id₁),)
end

@testset "Spin IIDSpace" begin
    @test shape(IIDSpace(SID{1//2}('x'), Spin{1//2}())) == (1:1,)
    @test shape(IIDSpace(SID{1//2}(:μ), Spin{1//2}())) == (1:5,)
end

@testset "Spin Coupling" begin
    @test collect(MatrixCoupling(:, SID, [1 0 0; 0 1 0; 0 0 1])) == [Coupling(:, SID, ('x', 'x')), Coupling(:, SID, ('y', 'y')), Coupling(:, SID, ('z', 'z'))]

    sc = Coupling(2.0, (1, 2), SID, ('+', '-'))
    bond = Bond(1, Point(1, [0.0], [0.0]), Point(2, [0.5], [0.0]))
    hilbert = Hilbert(site=>Spin{1}() for site=1:2)
    ex = expand(Val(:SpinTerm), sc, bond, hilbert)
    @test collect(ex) == [Operator(2.0, CompositeIndex(Index(1, SID{1}('+')), [0.0], [0.0]), CompositeIndex(Index(2, SID{1}('-')), [0.5], [0.0]))]
end

@testset "Heisenberg" begin
    @test Heisenberg"" == SparseMatrixCSC([1 0 0; 0 1 0; 0 0 1])
end

@testset "Ising" begin
    @test Ising"x" == SparseMatrixCSC([1 0 0; 0 0 0; 0 0 0])
    @test Ising"y" == SparseMatrixCSC([0 0 0; 0 1 0; 0 0 0])
    @test Ising"z" == SparseMatrixCSC([0 0 0; 0 0 0; 0 0 1])
end

@testset "Γ" begin
    @test Γ"x" == SparseMatrixCSC([0 0 0; 0 0 1; 0 1 0])
    @test Γ"y" == SparseMatrixCSC([0 0 1; 0 0 0; 1 0 0])
    @test Γ"z" == SparseMatrixCSC([0 1 0; 1 0 0; 0 0 0])
end

@testset "Γ′" begin
    @test Γ′"x" == SparseMatrixCSC([0 1 1; 1 0 0; 1 0 0])
    @test Γ′"y" == SparseMatrixCSC([0 1 0; 1 0 1; 0 1 0])
    @test Γ′"z" == SparseMatrixCSC([0 0 1; 0 0 1; 1 1 0])
end

@testset "DM" begin
    @test DM"x" == SparseMatrixCSC([0 0 0; 0 0 1; 0 -1 0])
    @test DM"y" == SparseMatrixCSC([0 0 -1; 0 0 0; 1 0 0])
    @test DM"z" == SparseMatrixCSC([0 1 0; -1 0 0; 0 0 0])
end

@testset "SpinTerm" begin
    point = Point(1, (0.5, 0.5), (0.0, 0.0))
    hilbert = Hilbert(point.site=>Spin{1//2}())
    term = SpinTerm(:h, 1.5, 0, Coupling(Index(1, SID('z'))))
    operators = Operators(
        Operator(1.5, CompositeIndex(Index(1, SID{1//2}('z')), [0.5, 0.5], [0.0, 0.0])),
    )
    @test expand(term, Bond(point), hilbert) == operators

    bond = Bond(1, Point(2, (0.5, 0.5), (0.0, 0.0)), Point(1, (0.0, 0.0), (0.0, 0.0)))
    hilbert = Hilbert(site=>Spin{1//2}() for site=1:2)
    term = SpinTerm(:J, 1.5, 1, MatrixCoupling(:, SID, Heisenberg""))
    operators = Operators(
        Operator(1.5, CompositeIndex(Index(2, SID{1//2}('x')), [0.5, 0.5], [0.0, 0.0]), CompositeIndex(Index(1, SID{1//2}('x')), [0.0, 0.0], [0.0, 0.0])),
        Operator(1.5, CompositeIndex(Index(2, SID{1//2}('y')), [0.5, 0.5], [0.0, 0.0]), CompositeIndex(Index(1, SID{1//2}('y')), [0.0, 0.0], [0.0, 0.0])),
        Operator(1.5, CompositeIndex(Index(2, SID{1//2}('z')), [0.5, 0.5], [0.0, 0.0]), CompositeIndex(Index(1, SID{1//2}('z')), [0.0, 0.0], [0.0, 0.0])),
    )
    @test expand(term, bond, hilbert) == operators
end

@testset "PID" begin
    @test PID('u', 'x')' == PID('u', 'x')
    @test PID('p', 'y')' == PID('p', 'y')
    @test statistics(PID('p', 'x')) == statistics(PID) == :b
    @test string(PID('p', :)) == "PID('p', :)"
    @test string(PID('u', 'x')) == "PID('u', 'x')"

    @test isdefinite(PID{Char})
    @test !isdefinite(PID{Symbol})
    @test !isdefinite(PID{typeof(:)})
    @test iidtype(PID, Char, Char) == PID{Char}
    @test iidtype(PID, Char, Symbol) == PID{Symbol}
    @test iidtype(PID, Char, typeof(:)) == PID{typeof(:)}
end

@testset "Phonon" begin
    pn = Phonon(3)
    @test shape(pn) == (1:2, 1:3)
    for i = 1:length(pn)
        @test PID(CartesianIndex(pn[i], pn), pn) == pn[i]
    end
    @test collect(pn) == [PID('u', 'x'), PID('p', 'x'), PID('u', 'y'), PID('p', 'y'), PID('u', 'z'), PID('p', 'z')]
end

@testset "latex" begin
    index = Index(1, PID('u', 'x'))
    @test script(Val(:BD), index.iid, latexofphonons) == "u"
    @test script(Val(:BD), index, latexofphonons) == "u"
    @test script(Val(:BD), CompositeIndex(index, [0.0, 0.0], [0.0, 0.0]), latexofphonons) == "u"
    @test script(Val(:site), index) == "1"
    @test script(Val(:direction), index.iid) == "x"
    @test script(Val(:direction), index) == "x"

    index = Index(2, PID('p', 'y'))
    @test script(Val(:BD), index.iid, latexofphonons) == "p"
    @test script(Val(:BD), index, latexofphonons) == "p"
    @test script(Val(:BD), CompositeIndex(index, [0.0, 0.0], [0.0, 0.0]), latexofphonons) == "p"
    @test script(Val(:site), index) == "2"
    @test script(Val(:direction), index.iid) == "y"
    @test script(Val(:direction), index) == "y"

    @test latexname(Index{<:Union{Int, Colon}, <:PID}) == Symbol("Index{Union{Int, Colon}, PID}")
    @test latexname(AbstractCompositeIndex{<:Index{<:Union{Int, Colon}, <:PID}}) == Symbol("AbstractCompositeIndex{Index{Union{Int, Colon}, PID}}")
    @test latexname(PID) == Symbol("PID")
end

@testset "PhononOperator" begin
    opt = Operator(1.0,
        CompositeIndex(Index(1, PID('p', 'x')), [0.0, 0.0], [0.0, 0.0]),
        CompositeIndex(Index(1, PID('p', 'x')), [0.0, 0.0], [0.0, 0.0])
        )
    @test opt' == Operator(1.0,
        CompositeIndex(Index(1, PID('p', 'x')), [0.0, 0.0], [0.0, 0.0]),
        CompositeIndex(Index(1, PID('p', 'x')), [0.0, 0.0], [0.0, 0.0])
        )
    @test latexstring(opt) == "(p^{x}_{1})^2"
end

@testset "permute" begin
    id₁ = CompositeIndex(Index(1, PID('u', 'x')), [0.0, 0.0], [0.0, 0.0])
    id₂ = CompositeIndex(Index(1, PID('p', 'x')), [0.0, 0.0], [0.0, 0.0])
    @test permute(id₁, id₂) == (Operator(+1im), Operator(1, id₂, id₁))
    @test permute(id₂, id₁) == (Operator(-1im), Operator(1, id₁, id₂))

    id₁ = CompositeIndex(Index(1, PID('u', 'x')), [0.0, 0.0], [0.0, 0.0])
    id₂ = CompositeIndex(Index(1, PID('u', 'x')), [0.0, 0.0], [0.0, 0.0])
    @test permute(id₁, id₂) == (Operator(1, id₂, id₁),)

    id₁ = CompositeIndex(Index(1, PID('u', 'x')), [0.0, 0.0], [0.0, 0.0])
    id₂ = CompositeIndex(Index(1, PID('p', 'y')), [0.0, 0.0], [0.0, 0.0])
    @test permute(id₁, id₂) == (Operator(1, id₂, id₁),)
end

@testset "Phonon IIDSpace" begin
    @test shape(IIDSpace(PID('u', :), Phonon(3))) == (1:1, 1:3)
    @test shape(IIDSpace(PID('u', 'x'), Phonon(3))) == (1:1, 1:1)
    @test shape(IIDSpace(PID('u', 'y'), Phonon(3))) == (1:1, 2:2)
    @test shape(IIDSpace(PID('u', 'z'), Phonon(3))) == (1:1, 3:3)

    @test shape(IIDSpace(PID('p', :), Phonon(2))) == (2:2, 1:2)
    @test shape(IIDSpace(PID('p', 'x'), Phonon(3))) == (2:2, 1:1)
    @test shape(IIDSpace(PID('p', 'y'), Phonon(3))) == (2:2, 2:2)
    @test shape(IIDSpace(PID('p', 'z'), Phonon(3))) == (2:2, 3:3)
end

@testset "Phonon Coupling" begin
    @test collect(MatrixCoupling(:, PID, [1 0 1; 0 1 0; 1 0 1])) == [
        Coupling(Index(:, PID('u', 'x')), Index(:, PID('u', 'x'))),
        Coupling(Index(:, PID('u', 'z')), Index(:, PID('u', 'x'))),
        Coupling(Index(:, PID('u', 'y')), Index(:, PID('u', 'y'))),
        Coupling(Index(:, PID('u', 'x')), Index(:, PID('u', 'z'))),
        Coupling(Index(:, PID('u', 'z')), Index(:, PID('u', 'z')))
    ]

    pnc = Coupling(2.0, @indexes(Index(1, PID('p', μ)), Index(1, PID('p', μ))))
    point = Point(1, [0.5, 0.0], [0.0, 0.0])
    hilbert = Hilbert(point.site=>Phonon(2))
    ex = expand(Val(:Kinetic), pnc, Bond(point), hilbert)
    @test collect(ex) == [
        Operator(2.0, CompositeIndex(Index(1, PID('p', 'x')), [0.5, 0.0], [0.0, 0.0]), CompositeIndex(Index(1, PID('p', 'x')), [0.5, 0.0], [0.0, 0.0])),
        Operator(2.0, CompositeIndex(Index(1, PID('p', 'y')), [0.5, 0.0], [0.0, 0.0]), CompositeIndex(Index(1, PID('p', 'y')), [0.5, 0.0], [0.0, 0.0]))
    ]

    pnc = Coupling(Index(1, PID('u', :)), Index(2, PID('u', :)))
    bond = Bond(1, Point(1, [0.0, 0.0], [0.0, 0.0]), Point(2, [0.5, 0.0], [0.0, 0.0]))
    hilbert = Hilbert(site=>Phonon(2) for site=1:2)
    ex = expand(Val(:Hooke), pnc, bond, hilbert)
    @test shape(ex) == (1:2, 1:2, 1:4)
    @test collect(ex) ==[
        Operator(+1.0, CompositeIndex(Index(1, PID('u', 'x')), [0.0, 0.0], [0.0, 0.0]), CompositeIndex(Index(1, PID('u', 'x')), [0.0, 0.0], [0.0, 0.0])),
        Operator(-0.0, CompositeIndex(Index(1, PID('u', 'y')), [0.0, 0.0], [0.0, 0.0]), CompositeIndex(Index(1, PID('u', 'x')), [0.0, 0.0], [0.0, 0.0])),
        Operator(-0.0, CompositeIndex(Index(1, PID('u', 'x')), [0.0, 0.0], [0.0, 0.0]), CompositeIndex(Index(1, PID('u', 'y')), [0.0, 0.0], [0.0, 0.0])),
        Operator(+0.0, CompositeIndex(Index(1, PID('u', 'y')), [0.0, 0.0], [0.0, 0.0]), CompositeIndex(Index(1, PID('u', 'y')), [0.0, 0.0], [0.0, 0.0])),
        Operator(-1.0, CompositeIndex(Index(1, PID('u', 'x')), [0.0, 0.0], [0.0, 0.0]), CompositeIndex(Index(2, PID('u', 'x')), [0.5, 0.0], [0.0, 0.0])),
        Operator(+0.0, CompositeIndex(Index(1, PID('u', 'y')), [0.0, 0.0], [0.0, 0.0]), CompositeIndex(Index(2, PID('u', 'x')), [0.5, 0.0], [0.0, 0.0])),
        Operator(+0.0, CompositeIndex(Index(1, PID('u', 'x')), [0.0, 0.0], [0.0, 0.0]), CompositeIndex(Index(2, PID('u', 'y')), [0.5, 0.0], [0.0, 0.0])),
        Operator(-0.0, CompositeIndex(Index(1, PID('u', 'y')), [0.0, 0.0], [0.0, 0.0]), CompositeIndex(Index(2, PID('u', 'y')), [0.5, 0.0], [0.0, 0.0])),
        Operator(-1.0, CompositeIndex(Index(2, PID('u', 'x')), [0.5, 0.0], [0.0, 0.0]), CompositeIndex(Index(1, PID('u', 'x')), [0.0, 0.0], [0.0, 0.0])),
        Operator(+0.0, CompositeIndex(Index(2, PID('u', 'y')), [0.5, 0.0], [0.0, 0.0]), CompositeIndex(Index(1, PID('u', 'x')), [0.0, 0.0], [0.0, 0.0])),
        Operator(+0.0, CompositeIndex(Index(2, PID('u', 'x')), [0.5, 0.0], [0.0, 0.0]), CompositeIndex(Index(1, PID('u', 'y')), [0.0, 0.0], [0.0, 0.0])),
        Operator(-0.0, CompositeIndex(Index(2, PID('u', 'y')), [0.5, 0.0], [0.0, 0.0]), CompositeIndex(Index(1, PID('u', 'y')), [0.0, 0.0], [0.0, 0.0])),
        Operator(+1.0, CompositeIndex(Index(2, PID('u', 'x')), [0.5, 0.0], [0.0, 0.0]), CompositeIndex(Index(2, PID('u', 'x')), [0.5, 0.0], [0.0, 0.0])),
        Operator(-0.0, CompositeIndex(Index(2, PID('u', 'y')), [0.5, 0.0], [0.0, 0.0]), CompositeIndex(Index(2, PID('u', 'x')), [0.5, 0.0], [0.0, 0.0])),
        Operator(-0.0, CompositeIndex(Index(2, PID('u', 'x')), [0.5, 0.0], [0.0, 0.0]), CompositeIndex(Index(2, PID('u', 'y')), [0.5, 0.0], [0.0, 0.0])),
        Operator(+0.0, CompositeIndex(Index(2, PID('u', 'y')), [0.5, 0.0], [0.0, 0.0]), CompositeIndex(Index(2, PID('u', 'y')), [0.5, 0.0], [0.0, 0.0]))
        ]
end

@testset "Kinetic" begin
    term = Kinetic(:T, 2.0)
    point = Point(1, [0.5, 0.0], [0.0, 0.0])
    hilbert = Hilbert(point.site=>Phonon(2))
    operators = Operators(
        Operator(2.0, CompositeIndex(Index(1, PID('p', 'x')), [0.5, 0.0], [0.0, 0.0]), CompositeIndex(Index(1, PID('p', 'x')), [0.5, 0.0], [0.0, 0.0])),
        Operator(2.0, CompositeIndex(Index(1, PID('p', 'y')), [0.5, 0.0], [0.0, 0.0]), CompositeIndex(Index(1, PID('p', 'y')), [0.5, 0.0], [0.0, 0.0]))
    )
    @test expand(term, Bond(point), hilbert) == operators
end

@testset "Hooke" begin
    term = Hooke(:V, 2.0, 1)

    bond = Bond(1, Point(1, [0.0, 0.0], [0.0, 0.0]), Point(2, [0.5, 0.0], [0.0, 0.0]))
    hilbert = Hilbert(site=>Phonon(2) for site=1:2)
    operators = Operators(
        Operator(+2.0, CompositeIndex(Index(1, PID('u', 'x')), [0.0, 0.0], [0.0, 0.0]), CompositeIndex(Index(1, PID('u', 'x')), [0.0, 0.0], [0.0, 0.0])),
        Operator(-2.0, CompositeIndex(Index(2, PID('u', 'x')), [0.5, 0.0], [0.0, 0.0]), CompositeIndex(Index(1, PID('u', 'x')), [0.0, 0.0], [0.0, 0.0])),
        Operator(-2.0, CompositeIndex(Index(1, PID('u', 'x')), [0.0, 0.0], [0.0, 0.0]), CompositeIndex(Index(2, PID('u', 'x')), [0.5, 0.0], [0.0, 0.0])),
        Operator(+2.0, CompositeIndex(Index(2, PID('u', 'x')), [0.5, 0.0], [0.0, 0.0]), CompositeIndex(Index(2, PID('u', 'x')), [0.5, 0.0], [0.0, 0.0]))
    )
    @test expand(term, bond, hilbert) == operators

    bond = Bond(1, Point(1, [0.0, 0.0], [0.0, 0.0]), Point(2, [0.0, 0.5], [0.0, 0.0]))
    hilbert = Hilbert(site=>Phonon(2) for site=1:2)
    operators = Operators(
        Operator(+2.0, CompositeIndex(Index(2, PID('u', 'y')), [0.0, 0.5], [0.0, 0.0]), CompositeIndex(Index(2, PID('u', 'y')), [0.0, 0.5], [0.0, 0.0])),
        Operator(+2.0, CompositeIndex(Index(1, PID('u', 'y')), [0.0, 0.0], [0.0, 0.0]), CompositeIndex(Index(1, PID('u', 'y')), [0.0, 0.0], [0.0, 0.0])),
        Operator(-2.0, CompositeIndex(Index(2, PID('u', 'y')), [0.0, 0.5], [0.0, 0.0]), CompositeIndex(Index(1, PID('u', 'y')), [0.0, 0.0], [0.0, 0.0])),
        Operator(-2.0, CompositeIndex(Index(1, PID('u', 'y')), [0.0, 0.0], [0.0, 0.0]), CompositeIndex(Index(2, PID('u', 'y')), [0.0, 0.5], [0.0, 0.0]))
    )
    @test expand(term, bond, hilbert) == operators

    bond = Bond(1, Point(1, [0.0, 0.0], [0.0, 0.0]), Point(2, [0.5, 0.5], [0.0, 0.0]))
    hilbert = Hilbert(site=>Phonon(2) for site=1:2)
    operators = Operators(
        Operator(-1.0, CompositeIndex(Index(2, PID('u', 'x')), [0.5, 0.5], [0.0, 0.0]), CompositeIndex(Index(1, PID('u', 'x')), [0.0, 0.0], [0.0, 0.0])),
        Operator(-1.0, CompositeIndex(Index(1, PID('u', 'y')), [0.0, 0.0], [0.0, 0.0]), CompositeIndex(Index(2, PID('u', 'y')), [0.5, 0.5], [0.0, 0.0])),
        Operator(-1.0, CompositeIndex(Index(1, PID('u', 'x')), [0.0, 0.0], [0.0, 0.0]), CompositeIndex(Index(2, PID('u', 'x')), [0.5, 0.5], [0.0, 0.0])),
        Operator(-1.0, CompositeIndex(Index(1, PID('u', 'y')), [0.0, 0.0], [0.0, 0.0]), CompositeIndex(Index(2, PID('u', 'x')), [0.5, 0.5], [0.0, 0.0])),
        Operator(-1.0, CompositeIndex(Index(2, PID('u', 'x')), [0.5, 0.5], [0.0, 0.0]), CompositeIndex(Index(1, PID('u', 'y')), [0.0, 0.0], [0.0, 0.0])),
        Operator(-1.0, CompositeIndex(Index(2, PID('u', 'y')), [0.5, 0.5], [0.0, 0.0]), CompositeIndex(Index(1, PID('u', 'y')), [0.0, 0.0], [0.0, 0.0])),
        Operator(+1.0, CompositeIndex(Index(1, PID('u', 'x')), [0.0, 0.0], [0.0, 0.0]), CompositeIndex(Index(1, PID('u', 'y')), [0.0, 0.0], [0.0, 0.0])),
        Operator(-1.0, CompositeIndex(Index(1, PID('u', 'x')), [0.0, 0.0], [0.0, 0.0]), CompositeIndex(Index(2, PID('u', 'y')), [0.5, 0.5], [0.0, 0.0])),
        Operator(+1.0, CompositeIndex(Index(1, PID('u', 'y')), [0.0, 0.0], [0.0, 0.0]), CompositeIndex(Index(1, PID('u', 'y')), [0.0, 0.0], [0.0, 0.0])),
        Operator(+1.0, CompositeIndex(Index(2, PID('u', 'y')), [0.5, 0.5], [0.0, 0.0]), CompositeIndex(Index(2, PID('u', 'y')), [0.5, 0.5], [0.0, 0.0])),
        Operator(+1.0, CompositeIndex(Index(1, PID('u', 'y')), [0.0, 0.0], [0.0, 0.0]), CompositeIndex(Index(1, PID('u', 'x')), [0.0, 0.0], [0.0, 0.0])),
        Operator(+1.0, CompositeIndex(Index(2, PID('u', 'y')), [0.5, 0.5], [0.0, 0.0]), CompositeIndex(Index(2, PID('u', 'x')), [0.5, 0.5], [0.0, 0.0])),
        Operator(+1.0, CompositeIndex(Index(2, PID('u', 'x')), [0.5, 0.5], [0.0, 0.0]), CompositeIndex(Index(2, PID('u', 'x')), [0.5, 0.5], [0.0, 0.0])),
        Operator(-1.0, CompositeIndex(Index(2, PID('u', 'y')), [0.5, 0.5], [0.0, 0.0]), CompositeIndex(Index(1, PID('u', 'x')), [0.0, 0.0], [0.0, 0.0])),
        Operator(+1.0, CompositeIndex(Index(1, PID('u', 'x')), [0.0, 0.0], [0.0, 0.0]), CompositeIndex(Index(1, PID('u', 'x')), [0.0, 0.0], [0.0, 0.0])),
        Operator(+1.0, CompositeIndex(Index(2, PID('u', 'x')), [0.5, 0.5], [0.0, 0.0]), CompositeIndex(Index(2, PID('u', 'y')), [0.5, 0.5], [0.0, 0.0]))
    )
    @test expand(term, bond, hilbert) ≈ operators
end

@testset "Elastic" begin
    term = Elastic(:V, 2.0, 1, MatrixCoupling(:, PID, [0 1; 1 0]))

    bond = Bond(1, Point(1, [0.0, 0.0], [0.0, 0.0]), Point(2, [0.5, 0.0], [0.0, 0.0]))
    hilbert = Hilbert(site=>Phonon(2) for site=1:2)
    operators = Operators(
        Operator(1.0, CompositeIndex(Index(1, PID('u', 'x')), [0.0, 0.0], [0.0, 0.0]), CompositeIndex(Index(2, PID('u', 'y')), [0.5, 0.0], [0.0, 0.0])),
        Operator(1.0, CompositeIndex(Index(1, PID('u', 'y')), [0.0, 0.0], [0.0, 0.0]), CompositeIndex(Index(2, PID('u', 'x')), [0.5, 0.0], [0.0, 0.0])),
        Operator(1.0, CompositeIndex(Index(2, PID('u', 'x')), [0.5, 0.0], [0.0, 0.0]), CompositeIndex(Index(1, PID('u', 'y')), [0.0, 0.0], [0.0, 0.0])),
        Operator(1.0, CompositeIndex(Index(2, PID('u', 'y')), [0.5, 0.0], [0.0, 0.0]), CompositeIndex(Index(1, PID('u', 'x')), [0.0, 0.0], [0.0, 0.0])),
    )
    @test expand(term, bond, hilbert) == operators
end
