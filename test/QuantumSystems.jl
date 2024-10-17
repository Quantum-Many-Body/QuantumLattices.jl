using LaTeXStrings: latexstring
using QuantumLattices: ⊗, expand, permute, rank
using QuantumLattices.DegreesOfFreedom: ˢᵗ, ⁿᵈ, AbstractIndex, CompositeIndex, InternalIndexProd, CoordinatedIndex, Coupling, Hilbert, Index, ConstrainedInternal, MatrixCoupling, allequalfields, indextype, isdefinite, patternrule, statistics, @pattern
using QuantumLattices.QuantumOperators: Operator, Operators, latexname, matrix, script
using QuantumLattices.QuantumSystems
using QuantumLattices.Spatials: Bond, Lattice, Neighbors, Point, azimuthd, bonds, rcoordinate, icoordinate
using QuantumLattices.Toolkit: Permutations, shape
using SparseArrays: SparseMatrixCSC
using StaticArrays: SVector

@testset "FockIndex" begin
    index = 𝕗(1, 1//2, 1)
    @test FockIndex{:f, Colon, Colon, Colon}(1, 1//2, 1) == index
    @test statistics(index) == statistics(typeof(index)) == :f
    @test isdefinite(index) == isdefinite(typeof(index)) == true
    @test index' == replace(index, nambu=2)
    @test isequal(index'', replace(index, nambu=1))
    @test hash(index) == hash((:f, 1, 1//2, 1))
    @test string(index) == "𝕗(1, 1//2, 1)"
    @test isannihilation(index) && isannihilation(𝕗(1, 1, 1//2, 1)) && isannihilation(𝕗(1, 1, 1//2, 1, [0.0], [0.0]))
    @test !iscreation(index) && !iscreation(𝕗(1, 1, 1//2, 1)) && !iscreation(𝕗(1, 1, 1//2, 1, [0.0], [0.0]))

    index = 𝕓(1, -1//2, 2)
    @test FockIndex{:b, Colon, Colon, Colon}(1, -1//2, 2) == index
    @test statistics(index) == statistics(typeof(index)) == :b
    @test isdefinite(index) == isdefinite(typeof(index)) == true
    @test hash(index) == hash((:b, 1, -1//2, 2))
    @test string(index) == "𝕓(1, -1//2, 2)"
    @test !isannihilation(index) && !isannihilation(𝕓(1, 1, -1//2, 2)) && !isannihilation(𝕓(1, 1, -1//2, 2, [0.0], [0.0]))
    @test iscreation(index) && iscreation(𝕓(1, 1, -1//2, 2)) && iscreation(𝕓(1, 1, -1//2, 2, [0.0], [0.0]))

    index = 𝔽(1, :α, :)
    @test FockIndex{:, Colon, Colon, Colon}(1, :α, :) == index
    @test statistics(index) == statistics(typeof(index)) == Colon()
    @test isdefinite(index) == isdefinite(typeof(index)) == false
    @test index == FockIndex(1, :α, :)
    @test hash(index) == hash((:, 1, :α, :))
    @test string(index) == "𝔽(1, α, :)"
    @test !isannihilation(index) && !isannihilation(𝔽(1, 1, :α, :)) && !isannihilation(𝔽(1, 1, :α, :, [0.0], [0.0]))
    @test !iscreation(index) && !iscreation(𝔽(1, 1, :α, :)) && !iscreation(𝔽(1, 1, :α, :, [0.0], [0.0]))

    @test 𝕗(1, 1//2, 1) ≠ 𝕓(1, 1//2, 1)
    @test isequal(𝕗(1, 1//2, 1), 𝕗(1, 1//2, 1))
    @test !isequal(𝕗(1, 1//2, 1), 𝕓(1, 1//2, 1))

    @test allequalfields(FockIndex) == (:orbital, :spin)
    @test isdefinite(FockIndex{:, Int, Rational{Int}, Int})
    @test !isdefinite(FockIndex{:f, Symbol, typeof(:), Int})
    @test indextype(FockIndex, Int, typeof(:), Int) == FockIndex{:, Int, typeof(:), Int}
    @test indextype(FockIndex{:f}, typeof(:), Symbol, Symbol) == FockIndex{:f, typeof(:), Symbol, Symbol}
    @test indextype(FockIndex{:b}, typeof(:), Symbol, Symbol) == FockIndex{:b, typeof(:), Symbol, Symbol}

    @test AbstractIndex[FockIndex{:f}] == AbstractIndex[Index{<:FockIndex{:f}}] ==  AbstractIndex[CoordinatedIndex{<:Index{<:FockIndex{:f}}}] == 𝕗
    @test AbstractIndex[FockIndex{:b}] == AbstractIndex[Index{<:FockIndex{:b}}] ==  AbstractIndex[CoordinatedIndex{<:Index{<:FockIndex{:b}}}] == 𝕓
    @test AbstractIndex[FockIndex{:}] == AbstractIndex[Index{<:FockIndex{:}}] ==  AbstractIndex[CoordinatedIndex{<:Index{<:FockIndex{:}}}] == 𝔽
    @test AbstractIndex[𝕗] == FockIndex{:f}
    @test AbstractIndex[𝕓] == FockIndex{:b}
    @test AbstractIndex[𝔽] == FockIndex{:}

    patternrule((:, :, :, :), Val(:), FockIndex, Val(:nambu)) == (2, 1, 2, 1)
end

@testset "Fock latex" begin
    @test script(𝕗(1, 2, 1//2, 1), Val(:site)) == script(𝕗(1, 2, 1//2, 1, [0.0], [0.0]), Val(:site)) == "1"
    @test script(𝕗(2, 1//2, 1), Val(:orbital)) == script(𝕗(1, 2, 1//2, 1), Val(:orbital)) == script(𝕗(1, 2, 1//2, 1, [0.0], [0.0]), Val(:orbital)) == "2"
    @test script(𝕗(2, 3//2, 1), Val(:spin)) == script(𝕗(1, 2, 3//2, 1), Val(:spin)) == script(𝕗(1, 2, 3//2, 1, [0.0], [0.0]), Val(:spin)) == "3//2"
    @test script(𝕗(2, 1//2, 1), Val(:spinsym)) == script(𝕗(1, 2, 1//2, 1), Val(:spinsym)) == script(𝕗(1, 2, 1//2, 1, [0.0], [0.0]), Val(:spinsym)) == "↑"
    @test script(𝕗(2, -1//2, 1), Val(:spinsym)) == script(𝕗(1, 2, -1//2, 1), Val(:spinsym)) == script(𝕗(1, 2, -1//2, 1, [0.0], [0.0]), Val(:spinsym)) == "↓"
    @test script(𝕗(2, 3//2, 1), Val(:nambu)) == script(𝕗(1, 2, 3//2, 1), Val(:nambu)) == script(𝕗(1, 2, 3//2, 1, [0.0], [0.0]), Val(:nambu)) == ""
    @test script(𝕗(2, 3//2, 2), Val(:nambu)) == script(𝕗(1, 2, 3//2, 2), Val(:nambu)) == script(𝕗(1, 2, 3//2, 2, [0.0], [0.0]), Val(:nambu)) == "\\dagger"

    @test latexname(FockIndex{:f}) == Symbol("FockIndex{:f}")
    @test latexname(Index{<:FockIndex{:f}}) == Symbol("Index{FockIndex{:f}}")
    @test latexname(CompositeIndex{Index{<:FockIndex{:f}}}) == Symbol("CompositeIndex{Index{FockIndex{:f}}}")

    @test latexname(FockIndex{:b}) == Symbol("FockIndex{:b}")
    @test latexname(Index{<:FockIndex{:b}}) == Symbol("Index{FockIndex{:b}}")
    @test latexname(CompositeIndex{Index{<:FockIndex{:b}}}) == Symbol("CompositeIndex{Index{FockIndex{:b}}}")

    @test latexname(FockIndex{:}) == Symbol("FockIndex")
    @test latexname(Index{<:FockIndex{:}}) == Symbol("Index{FockIndex}")
    @test latexname(CompositeIndex{Index{<:FockIndex{:}}}) == Symbol("CompositeIndex{Index{FockIndex}}")
end

@testset "Fock" begin
    @test eltype(Fock) == (FockIndex{S, Int, Rational{Int}, Int} where S)

    fock = Fock{:b}(1, 2)
    @test shape(fock) == (1:1, 1:2, 1:2)
    @test convert(CartesianIndex, 𝕓(1, -1//2, 1), fock) == CartesianIndex(1, 1, 1)
    @test convert(FockIndex, CartesianIndex(1, 1, 1), fock) == 𝕓(1, -1//2, 1)
    @test collect(fock) == [𝕓(1, -1//2, 1), 𝕓(1, 1//2, 1), 𝕓(1, -1//2, 2), 𝕓(1, 1//2, 2)]
    @test statistics(fock) == statistics(typeof(fock)) == :b
    @test string(fock) == "Fock{:b}(norbital=1, nspin=2)"

    @test summary(Fock{:b}(1, 0)) == "0-element Fock{:b}"
    @test summary(Fock{:f}(1, 1)) == "2-element Fock{:f}"

    @test match(FockIndex{:}, Fock{:f}) == match(FockIndex{:}, Fock{:b}) == true
    @test match(FockIndex{:f}, Fock{:f}) == match(FockIndex{:b}, Fock{:b}) == true
    @test match(FockIndex{:b}, Fock{:f}) == match(FockIndex{:f}, Fock{:b}) == false

    @test shape(Fock{:f}(3, 2), 𝕗(2, 1//2, 1)) == (2:2, 2:2, 1:1)
    @test shape(Fock{:f}(3, 2), 𝕗(1, :, 2)) ==(1:1, 1:2, 2:2) 
    @test shape(Fock{:f}(3, 2), 𝕗(:, -1//2, 1)) == (1:3, 1:1, 1:1)
    @test shape(Fock{:f}(3, 2), 𝕗(:, :, 2)) == (1:3, 1:2, 2:2)
end

@testset "angle" begin
    @test angle(𝕗(1, 1, 1//2, 1, [0.0, 0.0], [1.0, 2.0]), [[1.0, 0.0], [0.0, 1.0]], [0.1, 0.0]) ≈ 2pi*0.1
    @test angle(𝕗(1, 1, 1//2, 2, [0.0, 0.0], [1.0, 2.0]), [[1.0, 0.0], [0.0, 1.0]], [0.0, 0.2]) ≈ -2pi*0.4
end

@testset "Fock Operator" begin
    id₁ = 𝕗(2, 1, -1//2, 2, SVector(0.5, 0.0), SVector(0.0, 0.0))
    id₂ = 𝕗(2, 1, -1//2, 1, SVector(0.5, 0.0), SVector(0.0, 0.0))
    id₃ = 𝕗(1, 1, 1//2, 2, SVector(0.0, 0.0), SVector(0.0, 0.0))
    id₄ = 𝕗(1, 1, 1//2, 1, SVector(0.0, 0.0), SVector(0.0, 0.0))
    opt = Operator(1.0, id₁, id₂)
    @test opt|>isnormalordered
    opt = Operator(1.0, id₁, id₂, id₃, id₄)
    @test opt|>isnormalordered == false
    @test latexstring(opt) == "c^{\\dagger}_{2,\\,1,\\,↓}c^{}_{2,\\,1,\\,↓}c^{\\dagger}_{1,\\,1,\\,↑}c^{}_{1,\\,1,\\,↑}"
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

    id₁ = 𝕓(2, 1, -1//2, 2, SVector(0.5, 0.0), SVector(0.0, 0.0))
    id₂ = 𝕓(2, 1, -1//2, 1, SVector(0.5, 0.0), SVector(0.0, 0.0))
    id₃ = 𝕓(1, 1, 1//2, 2, SVector(0.0, 0.0), SVector(0.0, 0.0))
    id₄ = 𝕓(1, 1, 1//2, 1, SVector(0.0, 0.0), SVector(0.0, 0.0))
    opt = Operator(1.0, id₁, id₂)
    @test latexstring(opt) == "b^{\\dagger}_{2,\\,1,\\,↓}b^{}_{2,\\,1,\\,↓}"
    @test permute(id₁, id₂) == (Operator(-1), Operator(1, id₂, id₁))
    @test permute(id₂, id₁) == (Operator(+1), Operator(1, id₁, id₂))
    @test permute(id₁, id₄) == (Operator(1, id₄, id₁),)
    @test permute(id₄, id₁) == (Operator(1, id₁, id₄),)

    permute(𝕗(1, -1//2, 2), 𝕓(1, -1//2, 2)) == Operator(1, 𝕓(1, -1//2, 2), 𝕗(1, -1//2, 2))
    permute(𝕓(1, -1//2, 2), 𝕗(1, -1//2, 2)) == Operator(1, 𝕗(1, -1//2, 2), 𝕓(1, -1//2, 2))
end

@testset "Fock Coupling" begin
    @test collect(MatrixCoupling(:, FockIndex, :, :, :)) == collect(MatrixCoupling(:, 𝔽, :, :, :)) == [Coupling(𝔽(:, :, :, :), 𝔽(:, :, :, :))]
    @test collect(MatrixCoupling(:, FockIndex{:}, σ"+", σ"-", :)) == [Coupling(𝔽(:, 1, -1//2, :), 𝔽(:, 2, 1//2, :))]
    @test collect(MatrixCoupling((1ˢᵗ, 2ⁿᵈ), FockIndex{:f}, :, σ"y", σ"z")) == collect(MatrixCoupling((1ˢᵗ, 2ⁿᵈ), 𝕗, :, σ"y", σ"z")) == [
        Coupling(+1im, 𝕗(1ˢᵗ, :, -1//2, 1), 𝕗(2ⁿᵈ, :, 1//2, 2)), Coupling(-1im, 𝕗(1ˢᵗ, :, 1//2, 1), 𝕗(2ⁿᵈ, :, -1//2, 2)),
        Coupling(-1im, 𝕗(1ˢᵗ, :, -1//2, 2), 𝕗(2ⁿᵈ, :, 1//2, 1)), Coupling(+1im, 𝕗(1ˢᵗ, :, 1//2, 2), 𝕗(2ⁿᵈ, :, -1//2, 1))
    ]
    @test collect(MatrixCoupling((1ˢᵗ, 2ⁿᵈ), FockIndex{:b}, σ"x", :, σ"0")) == collect(MatrixCoupling((1ˢᵗ, 2ⁿᵈ), 𝕓, σ"x", :, σ"0")) == [
        Coupling(𝕓(1ˢᵗ, 2, :, 1), 𝕓(2ⁿᵈ, 1, :, 2)), Coupling(𝕓(1ˢᵗ, 1, :, 1), 𝕓(2ⁿᵈ, 2, :, 2)),
        Coupling(𝕓(1ˢᵗ, 2, :, 2), 𝕓(2ⁿᵈ, 1, :, 1)), Coupling(𝕓(1ˢᵗ, 1, :, 2), 𝕓(2ⁿᵈ, 2, :, 1))
    ]

    fc = Coupling(2.0, (1ˢᵗ, 2ⁿᵈ), FockIndex, (1, 2), :, (2, 1))
    bond = Bond(1, Point(1, SVector(0.0), SVector(0.0)), Point(2, SVector(0.5), SVector(0.0)))
    hilbert = Hilbert(site=>Fock{:f}(2, 2) for site=1:2)
    ex = expand(fc, Val(:Hopping), bond, hilbert)
    @test collect(ex) == [
        Operator(2.0, 𝕗(1, 1, -1//2, 2, SVector(0.0), SVector(0.0)), 𝕗(2, 2, -1//2, 1, SVector(0.5), SVector(0.0))),
        Operator(2.0, 𝕗(1, 1, +1//2, 2, SVector(0.0), SVector(0.0)), 𝕗(2, 2, +1//2, 1, SVector(0.5), SVector(0.0)))
    ]

    fc = Coupling(2.0, (1ˢᵗ, 1ˢᵗ, 1ˢᵗ, 1ˢᵗ), 𝔽, :, (1//2, 1//2, -1//2, -1//2), (2, 1, 2, 1))
    point = Point(1, SVector(0.0), SVector(0.0))
    hilbert = Hilbert(point.site=>Fock{:b}(2, 2))
    ex = expand(fc, Val(:term), Bond(point), hilbert)
    @test collect(ex) == [
        Operator(2.0, 𝕓(1, 1, +1//2, 2, SVector(0.0), SVector(0.0)), 𝕓(1, 1, +1//2, 1, SVector(0.0), SVector(0.0)), 𝕓(1, 1, -1//2, 2, SVector(0.0), SVector(0.0)), 𝕓(1, 1, -1//2, 1, SVector(0.0), SVector(0.0))),
        Operator(2.0, 𝕓(1, 2, +1//2, 2, SVector(0.0), SVector(0.0)), 𝕓(1, 2, +1//2, 1, SVector(0.0), SVector(0.0)), 𝕓(1, 2, -1//2, 2, SVector(0.0), SVector(0.0)), 𝕓(1, 2, -1//2, 1, SVector(0.0), SVector(0.0)))
    ]

    fc = Coupling(2.0, @pattern(𝔽(:, α, 1//2, 2), 𝔽(:, α, -1//2, 2), 𝔽(:, β, -1//2, 1), 𝔽(:, β, 1//2, 1); constraint=α<β))
    point = Point(1, SVector(0.5), SVector(0.0))
    hilbert = Hilbert(point.site=>Fock{:f}(3, 2))
    ex = expand(fc, Val(:term), Bond(point), hilbert)
    @test collect(ex) == [
        Operator(2.0, 𝕗(1, 1, +1//2, 2, SVector(0.5), SVector(0.0)), 𝕗(1, 1, -1//2, 2, SVector(0.5), SVector(0.0)), 𝕗(1, 2, -1//2, 1, SVector(0.5), SVector(0.0)), 𝕗(1, 2, +1//2, 1, SVector(0.5), SVector(0.0))),
        Operator(2.0, 𝕗(1, 1, +1//2, 2, SVector(0.5), SVector(0.0)), 𝕗(1, 1, -1//2, 2, SVector(0.5), SVector(0.0)), 𝕗(1, 3, -1//2, 1, SVector(0.5), SVector(0.0)), 𝕗(1, 3, +1//2, 1, SVector(0.5), SVector(0.0))),
        Operator(2.0, 𝕗(1, 2, +1//2, 2, SVector(0.5), SVector(0.0)), 𝕗(1, 2, -1//2, 2, SVector(0.5), SVector(0.0)), 𝕗(1, 3, -1//2, 1, SVector(0.5), SVector(0.0)), 𝕗(1, 3, +1//2, 1, SVector(0.5), SVector(0.0)))
    ]

    fc₁ = Coupling(+1.0, :, 𝔽, :, (+1//2, +1//2), (2, 1))
    fc₂ = Coupling(-1.0, :, 𝔽, :, (-1//2, -1//2), (2, 1))
    point = Point(1, SVector(0.0), SVector(0.0))
    hilbert = Hilbert(point.site=>Fock{:f}(2, 2))
    ex = expand(fc₁*fc₂, Val(:term), Bond(point), hilbert)
    @test collect(ex) == [
        Operator(-1.0, 𝕗(1, 1, +1//2, 2, SVector(0.0), SVector(0.0)), 𝕗(1, 1, +1//2, 1, SVector(0.0), SVector(0.0)), 𝕗(1, 1, -1//2, 2, SVector(0.0), SVector(0.0)), 𝕗(1, 1, -1//2, 1, SVector(0.0), SVector(0.0))),
        Operator(-1.0, 𝕗(1, 2, +1//2, 2, SVector(0.0), SVector(0.0)), 𝕗(1, 2, +1//2, 1, SVector(0.0), SVector(0.0)), 𝕗(1, 1, -1//2, 2, SVector(0.0), SVector(0.0)), 𝕗(1, 1, -1//2, 1, SVector(0.0), SVector(0.0))),
        Operator(-1.0, 𝕗(1, 1, +1//2, 2, SVector(0.0), SVector(0.0)), 𝕗(1, 1, +1//2, 1, SVector(0.0), SVector(0.0)), 𝕗(1, 2, -1//2, 2, SVector(0.0), SVector(0.0)), 𝕗(1, 2, -1//2, 1, SVector(0.0), SVector(0.0))),
        Operator(-1.0, 𝕗(1, 2, +1//2, 2, SVector(0.0), SVector(0.0)), 𝕗(1, 2, +1//2, 1, SVector(0.0), SVector(0.0)), 𝕗(1, 2, -1//2, 2, SVector(0.0), SVector(0.0)), 𝕗(1, 2, -1//2, 1, SVector(0.0), SVector(0.0)))
    ]
end

@testset "σ" begin
    @test σ"0" == SparseMatrixCSC([1 0; 0 1])
    @test σ"x" == SparseMatrixCSC([0 1; 1 0])
    @test σ"y" == SparseMatrixCSC([0 -1im; 1im 0])
    @test σ"z" == SparseMatrixCSC([1 0; 0 -1])
    @test σ"+" == SparseMatrixCSC([0 1; 0 0])
    @test σ"-" == SparseMatrixCSC([0 0; 1 0])
    @test σ"11" == SparseMatrixCSC([1 0; 0 0])
    @test σ"22" == SparseMatrixCSC([0 0; 0 1])
end

@testset "L" begin
    @test L"x" == SparseMatrixCSC([0 0 0; 0 0 1im; 0 -1im 0])
    @test L"y" == SparseMatrixCSC([0 0 -1im; 0 0 0; 1im 0 0])
    @test L"z" == SparseMatrixCSC([0 1im 0; -1im 0 0; 0 0 0])
end

@testset "Onsite" begin
    point = Point(1, (0.5, 0.5), (0.0, 0.0))
    bond = Bond(point)
    hilbert = Hilbert(point.site=>Fock{:f}(2, 2))

    term = Onsite(:mu, 1.5, MatrixCoupling(:, FockIndex, σ"z", σ"x", :))
    operators = Operators(
        Operator(-0.75, 𝕗(1, 2, +1//2, 2, [0.5, 0.5], [0.0, 0.0]), 𝕗(1, 2, -1//2, 1, [0.5, 0.5], [0.0, 0.0])),
        Operator(+0.75, 𝕗(1, 1, -1//2, 2, [0.5, 0.5], [0.0, 0.0]), 𝕗(1, 1, +1//2, 1, [0.5, 0.5], [0.0, 0.0])),
        Operator(-0.75, 𝕗(1, 2, -1//2, 2, [0.5, 0.5], [0.0, 0.0]), 𝕗(1, 2, +1//2, 1, [0.5, 0.5], [0.0, 0.0])),
        Operator(+0.75, 𝕗(1, 1, +1//2, 2, [0.5, 0.5], [0.0, 0.0]), 𝕗(1, 1, -1//2, 1, [0.5, 0.5], [0.0, 0.0]))
    )
    @test expand(term, bond, hilbert, half=true) == operators
    @test expand(term, bond, hilbert, half=false) == operators*2

    term = Onsite(:mu, 1.5, MatrixCoupling(:, FockIndex, σ"z", σ"z", :))
    operators = Operators(
        Operator(+0.75, 𝕗(1, 2, -1//2, 2, [0.5, 0.5], [0.0, 0.0]), 𝕗(1, 2, -1//2, 1, [0.5, 0.5], [0.0, 0.0])),
        Operator(+0.75, 𝕗(1, 1, +1//2, 2, [0.5, 0.5], [0.0, 0.0]), 𝕗(1, 1, +1//2, 1, [0.5, 0.5], [0.0, 0.0])),
        Operator(-0.75, 𝕗(1, 2, +1//2, 2, [0.5, 0.5], [0.0, 0.0]), 𝕗(1, 2, +1//2, 1, [0.5, 0.5], [0.0, 0.0])),
        Operator(-0.75, 𝕗(1, 1, -1//2, 2, [0.5, 0.5], [0.0, 0.0]), 𝕗(1, 1, -1//2, 1, [0.5, 0.5], [0.0, 0.0]))
    )
    @test expand(term, bond, hilbert, half=true) == operators
    @test expand(term, bond, hilbert, half=false) == operators*2
end

@testset "Hopping" begin
    bond = Bond(1, Point(2, (0.0, 0.0), (0.0, 0.0)), Point(1, (0.5, 0.5), (0.0, 0.0)))
    hilbert = Hilbert(site=>Fock{:f}(2, 2) for site=1:2)
    term = Hopping(:t, 1.5, 1)
    operators = Operators(
        Operator(1.5, 𝕗(2, 2, +1//2, 2, [0.0, 0.0], [0.0, 0.0]), 𝕗(1, 2, +1//2, 1, [0.5, 0.5], [0.0, 0.0])),
        Operator(1.5, 𝕗(2, 2, -1//2, 2, [0.0, 0.0], [0.0, 0.0]), 𝕗(1, 2, -1//2, 1, [0.5, 0.5], [0.0, 0.0])),
        Operator(1.5, 𝕗(2, 1, -1//2, 2, [0.0, 0.0], [0.0, 0.0]), 𝕗(1, 1, -1//2, 1, [0.5, 0.5], [0.0, 0.0])),
        Operator(1.5, 𝕗(2, 1, +1//2, 2, [0.0, 0.0], [0.0, 0.0]), 𝕗(1, 1, +1//2, 1, [0.5, 0.5], [0.0, 0.0]))
    )
    @test expand(term, bond, hilbert, half=true) == operators
    @test expand(term, bond, hilbert, half=false) == operators+operators'
end

@testset "Pairing" begin
    bond = Bond(1, Point(2, (0.0, 0.0), (0.0, 0.0)), Point(1, (0.5, 0.5), (0.0, 0.0)))
    hilbert = Hilbert(site=>Fock{:f}(1, 1) for site=1:2)
    term = Pairing(:Δ, 1.5, 1, Coupling{2}(:, 𝔽, :, :, :); amplitude=bond->(bond|>rcoordinate|>azimuthd ≈ 45 ? 1 : -1))
    operators = Operators(
        Operator(+1.5, 𝕗(2, 1, 0, 1, [0.0, 0.0], [0.0, 0.0]), 𝕗(1, 1, 0, 1, [0.5, 0.5], [0.0, 0.0])),
        Operator(-1.5, 𝕗(1, 1, 0, 1, [0.5, 0.5], [0.0, 0.0]), 𝕗(2, 1, 0, 1, [0.0, 0.0], [0.0, 0.0]))
    )
    @test expand(term, bond, hilbert, half=true) == operators
    @test expand(term, bond, hilbert, half=false) == operators+operators'

    point = Point(1, (0.5, 0.5), (0.0, 0.0))
    hilbert = Hilbert(point.site=>Fock{:f}(1, 2))
    term = Pairing(:Δ, 1.5, 0, MatrixCoupling(:, 𝔽, :, [0 -1; 1 0], :))
    operators = Operators(
        Operator(-1.5, 𝕗(1, 1, +1//2, 1, [0.5, 0.5], [0.0, 0.0]), 𝕗(1, 1, -1//2, 1, [0.5, 0.5], [0.0, 0.0])),
        Operator(+1.5, 𝕗(1, 1, -1//2, 1, [0.5, 0.5], [0.0, 0.0]), 𝕗(1, 1, +1//2, 1, [0.5, 0.5], [0.0, 0.0]))
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
        Operator(1.25, 𝕗(1, 1, +1//2, 2, [0.5, 0.5], [0.0, 0.0]), 𝕗(1, 1, +1//2, 1, [0.5, 0.5], [0.0, 0.0]), 𝕗(1, 1, -1//2, 2, [0.5, 0.5], [0.0, 0.0]), 𝕗(1, 1, -1//2, 1, [0.5, 0.5], [0.0, 0.0])),
        Operator(1.25, 𝕗(1, 2, +1//2, 2, [0.5, 0.5], [0.0, 0.0]), 𝕗(1, 2, +1//2, 1, [0.5, 0.5], [0.0, 0.0]), 𝕗(1, 2, -1//2, 2, [0.5, 0.5], [0.0, 0.0]), 𝕗(1, 2, -1//2, 1, [0.5, 0.5], [0.0, 0.0]))
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
        Operator(1.25, 𝕗(1, 1, +1//2, 2, [0.5, 0.5], [0.0, 0.0]), 𝕗(1, 1, +1//2, 1, [0.5, 0.5], [0.0, 0.0]), 𝕗(1, 2, -1//2, 2, [0.5, 0.5], [0.0, 0.0]), 𝕗(1, 2, -1//2, 1, [0.5, 0.5], [0.0, 0.0])),
        Operator(1.25, 𝕗(1, 1, -1//2, 2, [0.5, 0.5], [0.0, 0.0]), 𝕗(1, 1, -1//2, 1, [0.5, 0.5], [0.0, 0.0]), 𝕗(1, 2, +1//2, 2, [0.5, 0.5], [0.0, 0.0]), 𝕗(1, 2, +1//2, 1, [0.5, 0.5], [0.0, 0.0]))
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
        Operator(1.25, 𝕗(1, 1, -1//2, 2, [0.5, 0.5], [0.0, 0.0]), 𝕗(1, 1, -1//2, 1, [0.5, 0.5], [0.0, 0.0]), 𝕗(1, 2, -1//2, 2, [0.5, 0.5], [0.0, 0.0]), 𝕗(1, 2, -1//2, 1, [0.5, 0.5], [0.0, 0.0])),
        Operator(1.25, 𝕗(1, 1, 1//2, 2, [0.5, 0.5], [0.0, 0.0]), 𝕗(1, 1, 1//2, 1, [0.5, 0.5], [0.0, 0.0]), 𝕗(1, 2, 1//2, 2, [0.5, 0.5], [0.0, 0.0]), 𝕗(1, 2, 1//2, 1, [0.5, 0.5], [0.0, 0.0]))
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
        Operator(2.5, 𝕗(1, 1, +1//2, 2, [0.5, 0.5], [0.0, 0.0]), 𝕗(1, 2, -1//2, 2, [0.5, 0.5], [0.0, 0.0]), 𝕗(1, 1, -1//2, 1, [0.5, 0.5], [0.0, 0.0]), 𝕗(1, 2, +1//2, 1, [0.5, 0.5], [0.0, 0.0]))
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
        Operator(2.5, 𝕗(1, 1, +1//2, 2, [0.5, 0.5], [0.0, 0.0]), 𝕗(1, 1, -1//2, 2, [0.5, 0.5], [0.0, 0.0]), 𝕗(1, 2, -1//2, 1, [0.5, 0.5], [0.0, 0.0]), 𝕗(1, 2, +1//2, 1, [0.5, 0.5], [0.0, 0.0]))
    )
    @test expand(term, bond, hilbert, half=true) == operators
    @test expand(term, bond, hilbert, half=false) == operators+operators'
end

@testset "Coulomb" begin
    bond = Bond(1, Point(2, (0.0, 0.0), (0.0, 0.0)), Point(1, (0.5, 0.5), (0.0, 0.0)))
    hilbert = Hilbert(site=>Fock{:f}(1, 2) for site=1:2)

    term = Coulomb(:V, 2.5, 1, MatrixCoupling(:, FockIndex, :, σ"z", :)^2)
    operators = Operators(
        Operator(-1.25, 𝕗(2, 1, -1//2, 2, [0.0, 0.0], [0.0, 0.0]), 𝕗(2, 1, -1//2, 1, [0.0, 0.0], [0.0, 0.0]), 𝕗(1, 1, +1//2, 2, [0.5, 0.5], [0.0, 0.0]), 𝕗(1, 1, +1//2, 1, [0.5, 0.5], [0.0, 0.0])),
        Operator(+1.25, 𝕗(2, 1, -1//2, 2, [0.0, 0.0], [0.0, 0.0]), 𝕗(2, 1, -1//2, 1, [0.0, 0.0], [0.0, 0.0]), 𝕗(1, 1, -1//2, 2, [0.5, 0.5], [0.0, 0.0]), 𝕗(1, 1, -1//2, 1, [0.5, 0.5], [0.0, 0.0])),
        Operator(-1.25, 𝕗(2, 1, +1//2, 2, [0.0, 0.0], [0.0, 0.0]), 𝕗(2, 1, +1//2, 1, [0.0, 0.0], [0.0, 0.0]), 𝕗(1, 1, -1//2, 2, [0.5, 0.5], [0.0, 0.0]), 𝕗(1, 1, -1//2, 1, [0.5, 0.5], [0.0, 0.0])),
        Operator(+1.25, 𝕗(2, 1, +1//2, 2, [0.0, 0.0], [0.0, 0.0]), 𝕗(2, 1, +1//2, 1, [0.0, 0.0], [0.0, 0.0]), 𝕗(1, 1, +1//2, 2, [0.5, 0.5], [0.0, 0.0]), 𝕗(1, 1, +1//2, 1, [0.5, 0.5], [0.0, 0.0]))
    )
    @test expand(term, bond, hilbert, half=true) == operators
    @test expand(term, bond, hilbert, half=false) == operators*2

    term = Coulomb(:V, 2.5, 1, MatrixCoupling(:, FockIndex, :, σ"x", :)*MatrixCoupling(:, FockIndex, :, σ"z", :))
    operators = Operators(
        Operator(-1.25, 𝕗(2, 1, +1//2, 2, [0.0, 0.0], [0.0, 0.0]), 𝕗(2, 1, -1//2, 1, [0.0, 0.0], [0.0, 0.0]), 𝕗(1, 1, -1//2, 2, [0.5, 0.5], [0.0, 0.0]), 𝕗(1, 1, -1//2, 1, [0.5, 0.5], [0.0, 0.0])),
        Operator(+1.25, 𝕗(2, 1, -1//2, 2, [0.0, 0.0], [0.0, 0.0]), 𝕗(2, 1, +1//2, 1, [0.0, 0.0], [0.0, 0.0]), 𝕗(1, 1, +1//2, 2, [0.5, 0.5], [0.0, 0.0]), 𝕗(1, 1, +1//2, 1, [0.5, 0.5], [0.0, 0.0])),
        Operator(+1.25, 𝕗(2, 1, +1//2, 2, [0.0, 0.0], [0.0, 0.0]), 𝕗(2, 1, -1//2, 1, [0.0, 0.0], [0.0, 0.0]), 𝕗(1, 1, +1//2, 2, [0.5, 0.5], [0.0, 0.0]), 𝕗(1, 1, +1//2, 1, [0.5, 0.5], [0.0, 0.0])),
        Operator(-1.25, 𝕗(2, 1, -1//2, 2, [0.0, 0.0], [0.0, 0.0]), 𝕗(2, 1, +1//2, 1, [0.0, 0.0], [0.0, 0.0]), 𝕗(1, 1, -1//2, 2, [0.5, 0.5], [0.0, 0.0]), 𝕗(1, 1, -1//2, 1, [0.5, 0.5], [0.0, 0.0]))
    )
    @test expand(term, bond, hilbert, half=true) == operators
    @test expand(term, bond, hilbert, half=false) == operators*2
end

# @testset "SID" begin
#     @test SID{1//2}('z')' == SID{1//2}('z')
#     @test SID{1//2}('x')' == SID{1//2}('x')
#     @test SID{1//2}('y')' == SID{1//2}('y')
#     @test SID{1//2}('+')' == SID{1//2}('-')
#     @test SID{1//2}('-')' == SID{1//2}('+')

#     sid = SID{3//2}('x')
#     @test string(sid) == "SID{3//2}('x')"
#     @test replace(sid, tag='z') == SID{3//2}('z')
#     @test hash(sid) == hash((3//2, 'x'))
#     @test totalspin(sid) == totalspin(typeof(sid)) == 3//2
#     @test totalspin(Index(:, sid)) == totalspin(CoordinatedIndex(Index(:, sid), [0], [0])) == 3//2
#     @test statistics(sid) == statistics(typeof(sid)) == :b

#     sid = SID{:}('z')
#     @test sid == SID('z')
#     @test string(sid) == "SID('z')"

#     @test SID{1//2}('z')≠SID{3//2}('z')
#     @test isequal(SID{1//2}('z'), SID{1//2}('z'))
#     @test !isequal(SID{1//2}('z'), SID{3//2}('z'))

#     @test isdefinite(SID{:, Char})
#     @test !isdefinite(SID{1//2, Symbol})
#     @test !isdefinite(SID{1, typeof(:)})
#     @test indextype(SID, Char) == SID{:, Char}
#     @test indextype(SID{1//2}, Symbol) == SID{1//2, Symbol}
# end

# @testset "matrix" begin
#     @test isapprox(matrix(SID{1//2}('z')), [[-0.5, 0.0] [0.0, 0.5]])
#     @test isapprox(matrix(SID{1//2}('x')), [[0.0, 0.5] [0.5, 0.0]])
#     @test isapprox(matrix(SID{1//2}('y')), [[0.0, -0.5im] [0.5im, 0.0]])
#     @test isapprox(matrix(SID{1//2}('+')), [[0.0, 1.0] [0.0, 0.0]])
#     @test isapprox(matrix(SID{1//2}('-')), [[0.0, 0.0] [1.0, 0.0]])

#     @test isapprox(matrix(Index(:, SID{1//2}('z'))), [[-0.5, 0.0] [0.0, 0.5]])
#     @test isapprox(matrix(CoordinatedIndex(Index(:, SID{1//2}('z')), [0], [0])), [[-0.5, 0.0] [0.0, 0.5]])

#     @test isapprox(matrix(SID{1}('z')), [[-1.0, 0.0, 0.0] [0.0, 0.0, 0.0] [0.0, 0.0, 1.0]])
#     @test isapprox(matrix(SID{1}('x')), [[0.0, √2/2, 0.0] [√2/2, 0.0, √2/2] [0.0, √2/2, 0.0]])
#     @test isapprox(matrix(SID{1}('y')), [[0.0, -√2im/2, 0.0] [√2im/2, 0.0, -√2im/2] [0.0, √2im/2, 0.0]])
#     @test isapprox(matrix(SID{1}('+')), [[0.0, √2, 0.0] [0.0, 0.0, √2] [0.0, 0.0, 0.0]])
#     @test isapprox(matrix(SID{1}('-')), [[0.0, 0.0, 0.0] [√2, 0.0, 0.0] [0.0, √2, 0.0]])
# end

# @testset "Spin" begin
#     @test eltype(Spin) == (SID{S, Char} where S)
#     spin = Spin{1}()
#     @test shape(spin) == (1:5,)
#     @test CartesianIndex(SID{1}('z'), spin) == CartesianIndex(3)
#     @test SID(CartesianIndex(1), spin) == SID{1}('x')
#     @test summary(spin) == "5-element Spin{1}"
#     @test string(spin) == "Spin{1}()"
#     @test totalspin(spin) == totalspin(typeof(spin)) == 1
#     @test collect(spin) == [SID{1}('x'), SID{1}('y'), SID{1}('z'), SID{1}('+'), SID{1}('-')]

#     @test match(SID{:}, Spin{1//2}) == true
#     @test match(SID{1//2}, Spin{1//2}) == true
#     @test match(SID{1//2}, Spin{1}) == match(SID{1}, Spin{1//2}) == false
# end

# @testset "latex" begin
#     index = Index(1, SID{1//2}('z'))
#     @test script(Val(:site), index) == "1"
#     @test script(Val(:tag), index.iid) == script(Val(:tag), index) == "z"

#     @test latexname(Index{<:<:SID}) == Symbol("Index{SID}")
#     @test latexname(CompositeIndex{<:Index{<:<:SID}}) == Symbol("CompositeIndex{Index{SID}}")
#     @test latexname(SID) == Symbol("SID")
# end

# @testset "SpinOperator" begin
#     opt = Operator(
#         1.0,
#         CoordinatedIndex(Index(1, SID{1//2}('+')), [0.0, 0.0], [0.0, 0.0]),
#         CoordinatedIndex(Index(1, SID{1//2}('-')), [0.0, 0.0], [0.0, 0.0])
#     )
#     @test opt' == Operator(
#         1.0,
#         CoordinatedIndex(Index(1, SID{1//2}('+')), [0.0, 0.0], [0.0, 0.0]),
#         CoordinatedIndex(Index(1, SID{1//2}('-')), [0.0, 0.0], [0.0, 0.0])
#     )
#     @test latexstring(opt) == "S^{+}_{1}S^{-}_{1}"
# end

# @testset "permute" begin
#     soptrep(opt::Operator) = opt.value * prod([matrix(opt.id[i].index.iid) for i = 1:rank(opt)])
#     for S in (1//2, 1, 3//2)
#         indexes = [CoordinatedIndex(Index(1, SID{S}(tag)), [0.0, 0.0], [0.0, 0.0]) for tag in ('x', 'y', 'z', '+', '-')]
#         for (id₁, id₂) in Permutations{2}(indexes)
#             left = soptrep(Operator(1, id₁, id₂))
#             right = sum([soptrep(opt) for opt in permute(id₁, id₂)])
#             @test isapprox(left, right)
#         end
#     end
#     id₁ = CoordinatedIndex(Index(1, SID{1//2}('z')), [0.0, 0.0], [0.0, 0.0])
#     id₂ = CoordinatedIndex(Index(2, SID{1//2}('z')), [0.0, 0.0], [0.0, 0.0])
#     @test permute(id₁, id₂) == (Operator(1, id₂, id₁),)
# end

# @testset "Spin ConstrainedInternal" begin
#     @test shape(ConstrainedInternal(SID{1//2}('x'), Spin{1//2}())) == (1:1,)
#     @test shape(ConstrainedInternal(SID{1//2}(:μ), Spin{1//2}())) == (1:5,)
# end

# @testset "Spin Coupling" begin
#     @test collect(MatrixCoupling(:, SID, [1 0 0; 0 1 0; 0 0 1])) == [Coupling(:, SID, ('x', 'x')), Coupling(:, SID, ('y', 'y')), Coupling(:, SID, ('z', 'z'))]

#     sc = Coupling(2.0, (1, 2), SID, ('+', '-'))
#     bond = Bond(1, Point(1, [0.0], [0.0]), Point(2, [0.5], [0.0]))
#     hilbert = Hilbert(site=>Spin{1}() for site=1:2)
#     ex = expand(Val(:SpinTerm), sc, bond, hilbert)
#     @test collect(ex) == [Operator(2.0, CoordinatedIndex(Index(1, SID{1}('+')), [0.0], [0.0]), CoordinatedIndex(Index(2, SID{1}('-')), [0.5], [0.0]))]
# end

# @testset "Heisenberg" begin
#     @test Heisenberg"" == SparseMatrixCSC([1 0 0; 0 1 0; 0 0 1])
# end

# @testset "Ising" begin
#     @test Ising"x" == SparseMatrixCSC([1 0 0; 0 0 0; 0 0 0])
#     @test Ising"y" == SparseMatrixCSC([0 0 0; 0 1 0; 0 0 0])
#     @test Ising"z" == SparseMatrixCSC([0 0 0; 0 0 0; 0 0 1])
# end

# @testset "Γ" begin
#     @test Γ"x" == SparseMatrixCSC([0 0 0; 0 0 1; 0 1 0])
#     @test Γ"y" == SparseMatrixCSC([0 0 1; 0 0 0; 1 0 0])
#     @test Γ"z" == SparseMatrixCSC([0 1 0; 1 0 0; 0 0 0])
# end

# @testset "Γ′" begin
#     @test Γ′"x" == SparseMatrixCSC([0 1 1; 1 0 0; 1 0 0])
#     @test Γ′"y" == SparseMatrixCSC([0 1 0; 1 0 1; 0 1 0])
#     @test Γ′"z" == SparseMatrixCSC([0 0 1; 0 0 1; 1 1 0])
# end

# @testset "DM" begin
#     @test DM"x" == SparseMatrixCSC([0 0 0; 0 0 1; 0 -1 0])
#     @test DM"y" == SparseMatrixCSC([0 0 -1; 0 0 0; 1 0 0])
#     @test DM"z" == SparseMatrixCSC([0 1 0; -1 0 0; 0 0 0])
# end

# @testset "SpinTerm" begin
#     point = Point(1, (0.5, 0.5), (0.0, 0.0))
#     hilbert = Hilbert(point.site=>Spin{1//2}())
#     term = SpinTerm(:h, 1.5, 0, Coupling(Index(1, SID('z'))))
#     operators = Operators(
#         Operator(1.5, CoordinatedIndex(Index(1, SID{1//2}('z')), [0.5, 0.5], [0.0, 0.0])),
#     )
#     @test expand(term, Bond(point), hilbert) == operators

#     bond = Bond(1, Point(2, (0.5, 0.5), (0.0, 0.0)), Point(1, (0.0, 0.0), (0.0, 0.0)))
#     hilbert = Hilbert(site=>Spin{1//2}() for site=1:2)
#     term = SpinTerm(:J, 1.5, 1, MatrixCoupling(:, SID, Heisenberg""))
#     operators = Operators(
#         Operator(1.5, CoordinatedIndex(Index(2, SID{1//2}('x')), [0.5, 0.5], [0.0, 0.0]), CoordinatedIndex(Index(1, SID{1//2}('x')), [0.0, 0.0], [0.0, 0.0])),
#         Operator(1.5, CoordinatedIndex(Index(2, SID{1//2}('y')), [0.5, 0.5], [0.0, 0.0]), CoordinatedIndex(Index(1, SID{1//2}('y')), [0.0, 0.0], [0.0, 0.0])),
#         Operator(1.5, CoordinatedIndex(Index(2, SID{1//2}('z')), [0.5, 0.5], [0.0, 0.0]), CoordinatedIndex(Index(1, SID{1//2}('z')), [0.0, 0.0], [0.0, 0.0])),
#     )
#     @test expand(term, bond, hilbert) == operators
# end

# @testset "Zeeman" begin
#     point = Point(1, (0.5, 0.5), (0.0, 0.0))
#     hilbert = Hilbert(point.site=>Spin{1//2}())
#     term = Zeeman(:h, 1.5, 'x', 2)
#     operators = Operators(Operator(3.0, CoordinatedIndex(Index(1, SID{1//2}('x')), [0.5, 0.5], [0.0, 0.0])))
#     @test expand(term, Bond(point), hilbert) == operators

#     term = Zeeman(:h, 1.5, [1, 1, 1], 2)
#     operators = Operators(
#         Operator(√3, CoordinatedIndex(Index(1, SID{1//2}('x')), [0.5, 0.5], [0.0, 0.0])),
#         Operator(√3, CoordinatedIndex(Index(1, SID{1//2}('y')), [0.5, 0.5], [0.0, 0.0])),
#         Operator(√3, CoordinatedIndex(Index(1, SID{1//2}('z')), [0.5, 0.5], [0.0, 0.0]))
#     )
#     @test expand(term, Bond(point), hilbert) ≈ operators

#     term = Zeeman(:h, 1.5, [1, -1, 2], [1 0 0; 0 2 0; 0 0 3])
#     operators = Operators(
#         Operator(√6/4, CoordinatedIndex(Index(1, SID{1//2}('x')), [0.5, 0.5], [0.0, 0.0])),
#         Operator(-√6/2, CoordinatedIndex(Index(1, SID{1//2}('y')), [0.5, 0.5], [0.0, 0.0])),
#         Operator(3*√6/2, CoordinatedIndex(Index(1, SID{1//2}('z')), [0.5, 0.5], [0.0, 0.0]))
#     )
#     @test expand(term, Bond(point), hilbert) ≈ operators
# end

# @testset "SingleIonAnisotropy" begin
#     point = Point(1, (0.5, 0.5), (0.0, 0.0))
#     hilbert = Hilbert(point.site=>Spin{1//2}())
#     term = SingleIonAnisotropy(:A, 1.5, 'z')
#     operators = Operators(Operator(1.5, CoordinatedIndex(Index(1, SID{1//2}('z')), [0.5, 0.5], [0.0, 0.0]), CoordinatedIndex(Index(1, SID{1//2}('z')), [0.5, 0.5], [0.0, 0.0])))
#     @test expand(term, Bond(point), hilbert) == operators

#     term = SingleIonAnisotropy(:A, 1.5, [1 0 0; 0 2 0; 0 0 3])
#     operators = Operators(
#         Operator(1.5, CoordinatedIndex(Index(1, SID{1//2}('x')), [0.5, 0.5], [0.0, 0.0]), CoordinatedIndex(Index(1, SID{1//2}('x')), [0.5, 0.5], [0.0, 0.0])),
#         Operator(3.0, CoordinatedIndex(Index(1, SID{1//2}('y')), [0.5, 0.5], [0.0, 0.0]), CoordinatedIndex(Index(1, SID{1//2}('y')), [0.5, 0.5], [0.0, 0.0])),
#         Operator(4.5, CoordinatedIndex(Index(1, SID{1//2}('z')), [0.5, 0.5], [0.0, 0.0]), CoordinatedIndex(Index(1, SID{1//2}('z')), [0.5, 0.5], [0.0, 0.0]))
#     )
#     @test expand(term, Bond(point), hilbert) == operators
# end

# @testset "Ising" begin
#     bond = Bond(1, Point(1, (0.0, 0.0), (0.0, 0.0)), Point(2, (0.5, 0.5), (0.0, 0.0)))
#     hilbert = Hilbert(Spin{1//2}(), 2)
#     term = Ising(:J, 1.5, 1, 'x')
#     operators = Operators(Operator(1.5, CoordinatedIndex(Index(1, SID{1//2}('x')), [0.0, 0.0], [0.0, 0.0]), CoordinatedIndex(Index(2, SID{1//2}('x')), [0.5, 0.5], [0.0, 0.0])))
#     @test expand(term, bond, hilbert) == operators

#     term = Ising(:J, 1.5, 1, 'y')
#     operators = Operators(Operator(1.5, CoordinatedIndex(Index(1, SID{1//2}('y')), [0.0, 0.0], [0.0, 0.0]), CoordinatedIndex(Index(2, SID{1//2}('y')), [0.5, 0.5], [0.0, 0.0])))
#     @test expand(term, bond, hilbert) == operators

#     term = Ising(:J, 1.5, 1, 'z')
#     operators = Operators(Operator(1.5, CoordinatedIndex(Index(1, SID{1//2}('z')), [0.0, 0.0], [0.0, 0.0]), CoordinatedIndex(Index(2, SID{1//2}('z')), [0.5, 0.5], [0.0, 0.0])))
#     @test expand(term, bond, hilbert) == operators
# end

# @testset "Heisenberg" begin
#     bond = Bond(1, Point(1, (0.0, 0.0), (0.0, 0.0)), Point(2, (0.5, 0.5), (0.0, 0.0)))
#     hilbert = Hilbert(Spin{1//2}(), 2)
#     term = Heisenberg(:J, 1.5, 1; form=:xyz)
#     operators = Operators(
#         Operator(1.5, CoordinatedIndex(Index(1, SID{1//2}('x')), [0.0, 0.0], [0.0, 0.0]), CoordinatedIndex(Index(2, SID{1//2}('x')), [0.5, 0.5], [0.0, 0.0])),
#         Operator(1.5, CoordinatedIndex(Index(1, SID{1//2}('y')), [0.0, 0.0], [0.0, 0.0]), CoordinatedIndex(Index(2, SID{1//2}('y')), [0.5, 0.5], [0.0, 0.0])),
#         Operator(1.5, CoordinatedIndex(Index(1, SID{1//2}('z')), [0.0, 0.0], [0.0, 0.0]), CoordinatedIndex(Index(2, SID{1//2}('z')), [0.5, 0.5], [0.0, 0.0]))
#     )
#     @test expand(term, bond, hilbert) == operators

#     term = Heisenberg(:J, 1.5, 1; form=Symbol("+-z"))
#     operators = Operators(
#         Operator(0.75, CoordinatedIndex(Index(1, SID{1//2}('+')), [0.0, 0.0], [0.0, 0.0]), CoordinatedIndex(Index(2, SID{1//2}('-')), [0.5, 0.5], [0.0, 0.0])),
#         Operator(0.75, CoordinatedIndex(Index(1, SID{1//2}('-')), [0.0, 0.0], [0.0, 0.0]), CoordinatedIndex(Index(2, SID{1//2}('+')), [0.5, 0.5], [0.0, 0.0])),
#         Operator(1.5, CoordinatedIndex(Index(1, SID{1//2}('z')), [0.0, 0.0], [0.0, 0.0]), CoordinatedIndex(Index(2, SID{1//2}('z')), [0.5, 0.5], [0.0, 0.0]))
#     )
#     @test expand(term, bond, hilbert) == operators
# end

# @testset "Kitaev" begin
#     lattice = Lattice((0.0, 0.0), (0.0, √3/3); vectors=[[1.0, 0.0], [0.5, √3/2]])
#     bond₁, bond₂, bond₃ = bonds(lattice, Neighbors(1=>1/√3))
#     hilbert = Hilbert(Spin{1//2}(), length(lattice))
#     term = Kitaev(:K, 1.5, 1; x=[90], y=[210], z=[330], unit=:degree)

#     operators = Operators(Operator(1.5, CoordinatedIndex(Index(2, SID{1//2}('x')), -rcoordinate(bond₁), -icoordinate(bond₁)), CoordinatedIndex(Index(1, SID{1//2}('x')), [0.0, 0.0], [0.0, 0.0])))
#     @test expand(term, bond₁, hilbert) == operators
#     @test expand(term, reverse(bond₁), hilbert) == operators'

#     operators = Operators(Operator(1.5, CoordinatedIndex(Index(2, SID{1//2}('y')), -rcoordinate(bond₂), -icoordinate(bond₂)), CoordinatedIndex(Index(1, SID{1//2}('y')), [0.0, 0.0], [0.0, 0.0])))
#     @test expand(term, bond₂, hilbert) == operators
#     @test expand(term, reverse(bond₂), hilbert) == operators'

#     operators = Operators(Operator(1.5, CoordinatedIndex(Index(2, SID{1//2}('z')), -rcoordinate(bond₃), -icoordinate(bond₃)), CoordinatedIndex(Index(1, SID{1//2}('z')), [0.0, 0.0], [0.0, 0.0])))
#     @test expand(term, bond₃, hilbert) == operators
#     @test expand(term, reverse(bond₃), hilbert) == operators'
# end

# @testset "Γ" begin
#     lattice = Lattice((0.0, 0.0), (0.0, √3/3); vectors=[[1.0, 0.0], [0.5, √3/2]])
#     bond₁, bond₂, bond₃ = bonds(lattice, Neighbors(1=>1/√3))
#     hilbert = Hilbert(Spin{1//2}(), length(lattice))
#     term = Γ(:Γ, 1.5, 1; x=[90], y=[210], z=[330], unit=:degree)

#     operators = Operators(
#         Operator(1.5, CoordinatedIndex(Index(2, SID{1//2}('y')), -rcoordinate(bond₁), -icoordinate(bond₁)), CoordinatedIndex(Index(1, SID{1//2}('z')), [0.0, 0.0], [0.0, 0.0])),
#         Operator(1.5, CoordinatedIndex(Index(2, SID{1//2}('z')), -rcoordinate(bond₁), -icoordinate(bond₁)), CoordinatedIndex(Index(1, SID{1//2}('y')), [0.0, 0.0], [0.0, 0.0]))
#     )
#     @test expand(term, bond₁, hilbert) == operators
#     @test expand(term, reverse(bond₁), hilbert) == operators'

#     operators = Operators(
#         Operator(1.5, CoordinatedIndex(Index(2, SID{1//2}('z')), -rcoordinate(bond₂), -icoordinate(bond₂)), CoordinatedIndex(Index(1, SID{1//2}('x')), [0.0, 0.0], [0.0, 0.0])),
#         Operator(1.5, CoordinatedIndex(Index(2, SID{1//2}('x')), -rcoordinate(bond₂), -icoordinate(bond₂)), CoordinatedIndex(Index(1, SID{1//2}('z')), [0.0, 0.0], [0.0, 0.0]))
#     )
#     @test expand(term, bond₂, hilbert) == operators
#     @test expand(term, reverse(bond₂), hilbert) == operators'

#     operators = Operators(
#         Operator(1.5, CoordinatedIndex(Index(2, SID{1//2}('y')), -rcoordinate(bond₃), -icoordinate(bond₃)), CoordinatedIndex(Index(1, SID{1//2}('x')), [0.0, 0.0], [0.0, 0.0])),
#         Operator(1.5, CoordinatedIndex(Index(2, SID{1//2}('x')), -rcoordinate(bond₃), -icoordinate(bond₃)), CoordinatedIndex(Index(1, SID{1//2}('y')), [0.0, 0.0], [0.0, 0.0]))
#     )
#     @test expand(term, bond₃, hilbert) == operators
#     @test expand(term, reverse(bond₃), hilbert) == operators'
# end

# @testset "Γ′" begin
#     lattice = Lattice((0.0, 0.0), (0.0, √3/3); vectors=[[1.0, 0.0], [0.5, √3/2]])
#     bond₁, bond₂, bond₃ = bonds(lattice, Neighbors(1=>1/√3))
#     hilbert = Hilbert(Spin{1//2}(), length(lattice))
#     term = Γ′(:Γ′, 1.5, 1; x=[90], y=[210], z=[330], unit=:degree)

#     operators = Operators(
#         Operator(1.5, CoordinatedIndex(Index(2, SID{1//2}('z')), -rcoordinate(bond₁), -icoordinate(bond₁)), CoordinatedIndex(Index(1, SID{1//2}('x')), [0.0, 0.0], [0.0, 0.0])),
#         Operator(1.5, CoordinatedIndex(Index(2, SID{1//2}('x')), -rcoordinate(bond₁), -icoordinate(bond₁)), CoordinatedIndex(Index(1, SID{1//2}('y')), [0.0, 0.0], [0.0, 0.0])),
#         Operator(1.5, CoordinatedIndex(Index(2, SID{1//2}('x')), -rcoordinate(bond₁), -icoordinate(bond₁)), CoordinatedIndex(Index(1, SID{1//2}('z')), [0.0, 0.0], [0.0, 0.0])),
#         Operator(1.5, CoordinatedIndex(Index(2, SID{1//2}('y')), -rcoordinate(bond₁), -icoordinate(bond₁)), CoordinatedIndex(Index(1, SID{1//2}('x')), [0.0, 0.0], [0.0, 0.0]))
#     )
#     @test expand(term, bond₁, hilbert) == operators
#     @test expand(term, reverse(bond₁), hilbert) == operators'

#     operators = Operators(
#         Operator(1.5, CoordinatedIndex(Index(2, SID{1//2}('z')), -rcoordinate(bond₂), -icoordinate(bond₂)), CoordinatedIndex(Index(1, SID{1//2}('y')), [0.0, 0.0], [0.0, 0.0])),
#         Operator(1.5, CoordinatedIndex(Index(2, SID{1//2}('x')), -rcoordinate(bond₂), -icoordinate(bond₂)), CoordinatedIndex(Index(1, SID{1//2}('y')), [0.0, 0.0], [0.0, 0.0])),
#         Operator(1.5, CoordinatedIndex(Index(2, SID{1//2}('y')), -rcoordinate(bond₂), -icoordinate(bond₂)), CoordinatedIndex(Index(1, SID{1//2}('x')), [0.0, 0.0], [0.0, 0.0])),
#         Operator(1.5, CoordinatedIndex(Index(2, SID{1//2}('y')), -rcoordinate(bond₂), -icoordinate(bond₂)), CoordinatedIndex(Index(1, SID{1//2}('z')), [0.0, 0.0], [0.0, 0.0]))
#     )
#     @test expand(term, bond₂, hilbert) == operators
#     @test expand(term, reverse(bond₂), hilbert) == operators'

#     operators = Operators(
#         Operator(1.5, CoordinatedIndex(Index(2, SID{1//2}('z')), -rcoordinate(bond₃), -icoordinate(bond₃)), CoordinatedIndex(Index(1, SID{1//2}('x')), [0.0, 0.0], [0.0, 0.0])),
#         Operator(1.5, CoordinatedIndex(Index(2, SID{1//2}('y')), -rcoordinate(bond₃), -icoordinate(bond₃)), CoordinatedIndex(Index(1, SID{1//2}('z')), [0.0, 0.0], [0.0, 0.0])),
#         Operator(1.5, CoordinatedIndex(Index(2, SID{1//2}('z')), -rcoordinate(bond₃), -icoordinate(bond₃)), CoordinatedIndex(Index(1, SID{1//2}('y')), [0.0, 0.0], [0.0, 0.0])),
#         Operator(1.5, CoordinatedIndex(Index(2, SID{1//2}('x')), -rcoordinate(bond₃), -icoordinate(bond₃)), CoordinatedIndex(Index(1, SID{1//2}('z')), [0.0, 0.0], [0.0, 0.0]))
#     )
#     @test expand(term, bond₃, hilbert) == operators
#     @test expand(term, reverse(bond₃), hilbert) == operators'
# end

# @testset "DM" begin
#     lattice = Lattice((0.0, 0.0), (0.0, √3/3); vectors=[[1.0, 0.0], [0.5, √3/2]])
#     bond₁, bond₂, bond₃ = bonds(lattice, Neighbors(1=>1/√3))
#     hilbert = Hilbert(Spin{1//2}(), length(lattice))
#     term = DM(:DM, 1.5, 1, [90]=>'x', [210]=>'y', [330]=>'z'; unit=:degree)

#     operators = Operators(
#         Operator(-1.5, CoordinatedIndex(Index(2, SID{1//2}('y')), -rcoordinate(bond₁), -icoordinate(bond₁)), CoordinatedIndex(Index(1, SID{1//2}('z')), [0.0, 0.0], [0.0, 0.0])),
#         Operator(1.5, CoordinatedIndex(Index(2, SID{1//2}('z')), -rcoordinate(bond₁), -icoordinate(bond₁)), CoordinatedIndex(Index(1, SID{1//2}('y')), [0.0, 0.0], [0.0, 0.0]))
#     )
#     @test expand(term, bond₁, hilbert) == operators
#     @test expand(term, reverse(bond₁), hilbert) == operators'

#     operators = Operators(
#         Operator(-1.5, CoordinatedIndex(Index(2, SID{1//2}('z')), -rcoordinate(bond₂), -icoordinate(bond₂)), CoordinatedIndex(Index(1, SID{1//2}('x')), [0.0, 0.0], [0.0, 0.0])),
#         Operator(1.5, CoordinatedIndex(Index(2, SID{1//2}('x')), -rcoordinate(bond₂), -icoordinate(bond₂)), CoordinatedIndex(Index(1, SID{1//2}('z')), [0.0, 0.0], [0.0, 0.0]))
#     )
#     @test expand(term, bond₂, hilbert) == operators
#     @test expand(term, reverse(bond₂), hilbert) == operators'

#     operators = Operators(
#         Operator(1.5, CoordinatedIndex(Index(2, SID{1//2}('y')), -rcoordinate(bond₃), -icoordinate(bond₃)), CoordinatedIndex(Index(1, SID{1//2}('x')), [0.0, 0.0], [0.0, 0.0])),
#         Operator(-1.5, CoordinatedIndex(Index(2, SID{1//2}('x')), -rcoordinate(bond₃), -icoordinate(bond₃)), CoordinatedIndex(Index(1, SID{1//2}('y')), [0.0, 0.0], [0.0, 0.0]))
#     )
#     @test expand(term, bond₃, hilbert) == operators
#     @test expand(term, reverse(bond₃), hilbert) == operators'
# end

# @testset "PID" begin
#     @test PID('u', 'x')' == PID('u', 'x')
#     @test PID('p', 'y')' == PID('p', 'y')
#     @test statistics(PID('p', 'x')) == statistics(PID) == :b
#     @test string(PID('p', :)) == "PID('p', :)"
#     @test string(PID('u', 'x')) == "PID('u', 'x')"

#     @test isdefinite(PID{Char})
#     @test !isdefinite(PID{Symbol})
#     @test !isdefinite(PID{typeof(:)})
#     @test indextype(PID, Char, Char) == PID{Char}
#     @test indextype(PID, Char, Symbol) == PID{Symbol}
#     @test indextype(PID, Char, typeof(:)) == PID{typeof(:)}
# end

# @testset "Phonon" begin
#     pn = Phonon(3)
#     @test shape(pn) == (1:2, 1:3)
#     for i = 1:length(pn)
#         @test PID(CartesianIndex(pn[i], pn), pn) == pn[i]
#     end
#     @test collect(pn) == [PID('u', 'x'), PID('p', 'x'), PID('u', 'y'), PID('p', 'y'), PID('u', 'z'), PID('p', 'z')]
# end

# @testset "latex" begin
#     index = Index(1, PID('u', 'x'))
#     @test script(Val(:BD), index.iid, latexofphonons) == "u"
#     @test script(Val(:BD), index, latexofphonons) == "u"
#     @test script(Val(:BD), CoordinatedIndex(index, [0.0, 0.0], [0.0, 0.0]), latexofphonons) == "u"
#     @test script(Val(:site), index) == "1"
#     @test script(Val(:direction), index.iid) == "x"
#     @test script(Val(:direction), index) == "x"

#     index = Index(2, PID('p', 'y'))
#     @test script(Val(:BD), index.iid, latexofphonons) == "p"
#     @test script(Val(:BD), index, latexofphonons) == "p"
#     @test script(Val(:BD), CoordinatedIndex(index, [0.0, 0.0], [0.0, 0.0]), latexofphonons) == "p"
#     @test script(Val(:site), index) == "2"
#     @test script(Val(:direction), index.iid) == "y"
#     @test script(Val(:direction), index) == "y"

#     @test latexname(Index{<:<:PID}) == Symbol("Index{PID}")
#     @test latexname(CompositeIndex{<:Index{<:<:PID}}) == Symbol("CompositeIndex{Index{PID}}")
#     @test latexname(PID) == Symbol("PID")
# end

# @testset "PhononOperator" begin
#     opt = Operator(1.0,
#         CoordinatedIndex(Index(1, PID('p', 'x')), [0.0, 0.0], [0.0, 0.0]),
#         CoordinatedIndex(Index(1, PID('p', 'x')), [0.0, 0.0], [0.0, 0.0])
#         )
#     @test opt' == Operator(1.0,
#         CoordinatedIndex(Index(1, PID('p', 'x')), [0.0, 0.0], [0.0, 0.0]),
#         CoordinatedIndex(Index(1, PID('p', 'x')), [0.0, 0.0], [0.0, 0.0])
#         )
#     @test latexstring(opt) == "(p^{x}_{1})^2"
# end

# @testset "permute" begin
#     id₁ = CoordinatedIndex(Index(1, PID('u', 'x')), [0.0, 0.0], [0.0, 0.0])
#     id₂ = CoordinatedIndex(Index(1, PID('p', 'x')), [0.0, 0.0], [0.0, 0.0])
#     @test permute(id₁, id₂) == (Operator(+1im), Operator(1, id₂, id₁))
#     @test permute(id₂, id₁) == (Operator(-1im), Operator(1, id₁, id₂))

#     id₁ = CoordinatedIndex(Index(1, PID('u', 'x')), [0.0, 0.0], [0.0, 0.0])
#     id₂ = CoordinatedIndex(Index(1, PID('u', 'x')), [0.0, 0.0], [0.0, 0.0])
#     @test permute(id₁, id₂) == (Operator(1, id₂, id₁),)

#     id₁ = CoordinatedIndex(Index(1, PID('u', 'x')), [0.0, 0.0], [0.0, 0.0])
#     id₂ = CoordinatedIndex(Index(1, PID('p', 'y')), [0.0, 0.0], [0.0, 0.0])
#     @test permute(id₁, id₂) == (Operator(1, id₂, id₁),)
# end

# @testset "Phonon ConstrainedInternal" begin
#     @test shape(ConstrainedInternal(PID('u', :), Phonon(3))) == (1:1, 1:3)
#     @test shape(ConstrainedInternal(PID('u', 'x'), Phonon(3))) == (1:1, 1:1)
#     @test shape(ConstrainedInternal(PID('u', 'y'), Phonon(3))) == (1:1, 2:2)
#     @test shape(ConstrainedInternal(PID('u', 'z'), Phonon(3))) == (1:1, 3:3)

#     @test shape(ConstrainedInternal(PID('p', :), Phonon(2))) == (2:2, 1:2)
#     @test shape(ConstrainedInternal(PID('p', 'x'), Phonon(3))) == (2:2, 1:1)
#     @test shape(ConstrainedInternal(PID('p', 'y'), Phonon(3))) == (2:2, 2:2)
#     @test shape(ConstrainedInternal(PID('p', 'z'), Phonon(3))) == (2:2, 3:3)
# end

# @testset "Phonon Coupling" begin
#     @test collect(MatrixCoupling(:, PID, [1 0 1; 0 1 0; 1 0 1])) == [
#         Coupling(Index(:, PID('u', 'x')), Index(:, PID('u', 'x'))),
#         Coupling(Index(:, PID('u', 'z')), Index(:, PID('u', 'x'))),
#         Coupling(Index(:, PID('u', 'y')), Index(:, PID('u', 'y'))),
#         Coupling(Index(:, PID('u', 'x')), Index(:, PID('u', 'z'))),
#         Coupling(Index(:, PID('u', 'z')), Index(:, PID('u', 'z')))
#     ]

#     pnc = Coupling(2.0, @pattern(Index(1, PID('p', μ)), Index(1, PID('p', μ))))
#     point = Point(1, [0.5, 0.0], [0.0, 0.0])
#     hilbert = Hilbert(point.site=>Phonon(2))
#     ex = expand(Val(:Kinetic), pnc, Bond(point), hilbert)
#     @test collect(ex) == [
#         Operator(2.0, CoordinatedIndex(Index(1, PID('p', 'x')), [0.5, 0.0], [0.0, 0.0]), CoordinatedIndex(Index(1, PID('p', 'x')), [0.5, 0.0], [0.0, 0.0])),
#         Operator(2.0, CoordinatedIndex(Index(1, PID('p', 'y')), [0.5, 0.0], [0.0, 0.0]), CoordinatedIndex(Index(1, PID('p', 'y')), [0.5, 0.0], [0.0, 0.0]))
#     ]

#     pnc = Coupling(Index(1, PID('u', :)), Index(2, PID('u', :)))
#     bond = Bond(1, Point(1, [0.0, 0.0], [0.0, 0.0]), Point(2, [0.5, 0.0], [0.0, 0.0]))
#     hilbert = Hilbert(site=>Phonon(2) for site=1:2)
#     ex = expand(Val(:Hooke), pnc, bond, hilbert)
#     @test shape(ex) == (1:2, 1:2, 1:4)
#     @test collect(ex) ==[
#         Operator(+1.0, CoordinatedIndex(Index(1, PID('u', 'x')), [0.0, 0.0], [0.0, 0.0]), CoordinatedIndex(Index(1, PID('u', 'x')), [0.0, 0.0], [0.0, 0.0])),
#         Operator(-0.0, CoordinatedIndex(Index(1, PID('u', 'y')), [0.0, 0.0], [0.0, 0.0]), CoordinatedIndex(Index(1, PID('u', 'x')), [0.0, 0.0], [0.0, 0.0])),
#         Operator(-0.0, CoordinatedIndex(Index(1, PID('u', 'x')), [0.0, 0.0], [0.0, 0.0]), CoordinatedIndex(Index(1, PID('u', 'y')), [0.0, 0.0], [0.0, 0.0])),
#         Operator(+0.0, CoordinatedIndex(Index(1, PID('u', 'y')), [0.0, 0.0], [0.0, 0.0]), CoordinatedIndex(Index(1, PID('u', 'y')), [0.0, 0.0], [0.0, 0.0])),
#         Operator(-1.0, CoordinatedIndex(Index(1, PID('u', 'x')), [0.0, 0.0], [0.0, 0.0]), CoordinatedIndex(Index(2, PID('u', 'x')), [0.5, 0.0], [0.0, 0.0])),
#         Operator(+0.0, CoordinatedIndex(Index(1, PID('u', 'y')), [0.0, 0.0], [0.0, 0.0]), CoordinatedIndex(Index(2, PID('u', 'x')), [0.5, 0.0], [0.0, 0.0])),
#         Operator(+0.0, CoordinatedIndex(Index(1, PID('u', 'x')), [0.0, 0.0], [0.0, 0.0]), CoordinatedIndex(Index(2, PID('u', 'y')), [0.5, 0.0], [0.0, 0.0])),
#         Operator(-0.0, CoordinatedIndex(Index(1, PID('u', 'y')), [0.0, 0.0], [0.0, 0.0]), CoordinatedIndex(Index(2, PID('u', 'y')), [0.5, 0.0], [0.0, 0.0])),
#         Operator(-1.0, CoordinatedIndex(Index(2, PID('u', 'x')), [0.5, 0.0], [0.0, 0.0]), CoordinatedIndex(Index(1, PID('u', 'x')), [0.0, 0.0], [0.0, 0.0])),
#         Operator(+0.0, CoordinatedIndex(Index(2, PID('u', 'y')), [0.5, 0.0], [0.0, 0.0]), CoordinatedIndex(Index(1, PID('u', 'x')), [0.0, 0.0], [0.0, 0.0])),
#         Operator(+0.0, CoordinatedIndex(Index(2, PID('u', 'x')), [0.5, 0.0], [0.0, 0.0]), CoordinatedIndex(Index(1, PID('u', 'y')), [0.0, 0.0], [0.0, 0.0])),
#         Operator(-0.0, CoordinatedIndex(Index(2, PID('u', 'y')), [0.5, 0.0], [0.0, 0.0]), CoordinatedIndex(Index(1, PID('u', 'y')), [0.0, 0.0], [0.0, 0.0])),
#         Operator(+1.0, CoordinatedIndex(Index(2, PID('u', 'x')), [0.5, 0.0], [0.0, 0.0]), CoordinatedIndex(Index(2, PID('u', 'x')), [0.5, 0.0], [0.0, 0.0])),
#         Operator(-0.0, CoordinatedIndex(Index(2, PID('u', 'y')), [0.5, 0.0], [0.0, 0.0]), CoordinatedIndex(Index(2, PID('u', 'x')), [0.5, 0.0], [0.0, 0.0])),
#         Operator(-0.0, CoordinatedIndex(Index(2, PID('u', 'x')), [0.5, 0.0], [0.0, 0.0]), CoordinatedIndex(Index(2, PID('u', 'y')), [0.5, 0.0], [0.0, 0.0])),
#         Operator(+0.0, CoordinatedIndex(Index(2, PID('u', 'y')), [0.5, 0.0], [0.0, 0.0]), CoordinatedIndex(Index(2, PID('u', 'y')), [0.5, 0.0], [0.0, 0.0]))
#         ]
# end

# @testset "Kinetic" begin
#     term = Kinetic(:T, 2.0)
#     point = Point(1, [0.5, 0.0], [0.0, 0.0])
#     hilbert = Hilbert(point.site=>Phonon(2))
#     operators = Operators(
#         Operator(2.0, CoordinatedIndex(Index(1, PID('p', 'x')), [0.5, 0.0], [0.0, 0.0]), CoordinatedIndex(Index(1, PID('p', 'x')), [0.5, 0.0], [0.0, 0.0])),
#         Operator(2.0, CoordinatedIndex(Index(1, PID('p', 'y')), [0.5, 0.0], [0.0, 0.0]), CoordinatedIndex(Index(1, PID('p', 'y')), [0.5, 0.0], [0.0, 0.0]))
#     )
#     @test expand(term, Bond(point), hilbert) == operators
# end

# @testset "Hooke" begin
#     term = Hooke(:V, 2.0, 1)

#     bond = Bond(1, Point(1, [0.0, 0.0], [0.0, 0.0]), Point(2, [0.5, 0.0], [0.0, 0.0]))
#     hilbert = Hilbert(site=>Phonon(2) for site=1:2)
#     operators = Operators(
#         Operator(+2.0, CoordinatedIndex(Index(1, PID('u', 'x')), [0.0, 0.0], [0.0, 0.0]), CoordinatedIndex(Index(1, PID('u', 'x')), [0.0, 0.0], [0.0, 0.0])),
#         Operator(-2.0, CoordinatedIndex(Index(2, PID('u', 'x')), [0.5, 0.0], [0.0, 0.0]), CoordinatedIndex(Index(1, PID('u', 'x')), [0.0, 0.0], [0.0, 0.0])),
#         Operator(-2.0, CoordinatedIndex(Index(1, PID('u', 'x')), [0.0, 0.0], [0.0, 0.0]), CoordinatedIndex(Index(2, PID('u', 'x')), [0.5, 0.0], [0.0, 0.0])),
#         Operator(+2.0, CoordinatedIndex(Index(2, PID('u', 'x')), [0.5, 0.0], [0.0, 0.0]), CoordinatedIndex(Index(2, PID('u', 'x')), [0.5, 0.0], [0.0, 0.0]))
#     )
#     @test expand(term, bond, hilbert) == operators

#     bond = Bond(1, Point(1, [0.0, 0.0], [0.0, 0.0]), Point(2, [0.0, 0.5], [0.0, 0.0]))
#     hilbert = Hilbert(site=>Phonon(2) for site=1:2)
#     operators = Operators(
#         Operator(+2.0, CoordinatedIndex(Index(2, PID('u', 'y')), [0.0, 0.5], [0.0, 0.0]), CoordinatedIndex(Index(2, PID('u', 'y')), [0.0, 0.5], [0.0, 0.0])),
#         Operator(+2.0, CoordinatedIndex(Index(1, PID('u', 'y')), [0.0, 0.0], [0.0, 0.0]), CoordinatedIndex(Index(1, PID('u', 'y')), [0.0, 0.0], [0.0, 0.0])),
#         Operator(-2.0, CoordinatedIndex(Index(2, PID('u', 'y')), [0.0, 0.5], [0.0, 0.0]), CoordinatedIndex(Index(1, PID('u', 'y')), [0.0, 0.0], [0.0, 0.0])),
#         Operator(-2.0, CoordinatedIndex(Index(1, PID('u', 'y')), [0.0, 0.0], [0.0, 0.0]), CoordinatedIndex(Index(2, PID('u', 'y')), [0.0, 0.5], [0.0, 0.0]))
#     )
#     @test expand(term, bond, hilbert) == operators

#     bond = Bond(1, Point(1, [0.0, 0.0], [0.0, 0.0]), Point(2, [0.5, 0.5], [0.0, 0.0]))
#     hilbert = Hilbert(site=>Phonon(2) for site=1:2)
#     operators = Operators(
#         Operator(-1.0, CoordinatedIndex(Index(2, PID('u', 'x')), [0.5, 0.5], [0.0, 0.0]), CoordinatedIndex(Index(1, PID('u', 'x')), [0.0, 0.0], [0.0, 0.0])),
#         Operator(-1.0, CoordinatedIndex(Index(1, PID('u', 'y')), [0.0, 0.0], [0.0, 0.0]), CoordinatedIndex(Index(2, PID('u', 'y')), [0.5, 0.5], [0.0, 0.0])),
#         Operator(-1.0, CoordinatedIndex(Index(1, PID('u', 'x')), [0.0, 0.0], [0.0, 0.0]), CoordinatedIndex(Index(2, PID('u', 'x')), [0.5, 0.5], [0.0, 0.0])),
#         Operator(-1.0, CoordinatedIndex(Index(1, PID('u', 'y')), [0.0, 0.0], [0.0, 0.0]), CoordinatedIndex(Index(2, PID('u', 'x')), [0.5, 0.5], [0.0, 0.0])),
#         Operator(-1.0, CoordinatedIndex(Index(2, PID('u', 'x')), [0.5, 0.5], [0.0, 0.0]), CoordinatedIndex(Index(1, PID('u', 'y')), [0.0, 0.0], [0.0, 0.0])),
#         Operator(-1.0, CoordinatedIndex(Index(2, PID('u', 'y')), [0.5, 0.5], [0.0, 0.0]), CoordinatedIndex(Index(1, PID('u', 'y')), [0.0, 0.0], [0.0, 0.0])),
#         Operator(+1.0, CoordinatedIndex(Index(1, PID('u', 'x')), [0.0, 0.0], [0.0, 0.0]), CoordinatedIndex(Index(1, PID('u', 'y')), [0.0, 0.0], [0.0, 0.0])),
#         Operator(-1.0, CoordinatedIndex(Index(1, PID('u', 'x')), [0.0, 0.0], [0.0, 0.0]), CoordinatedIndex(Index(2, PID('u', 'y')), [0.5, 0.5], [0.0, 0.0])),
#         Operator(+1.0, CoordinatedIndex(Index(1, PID('u', 'y')), [0.0, 0.0], [0.0, 0.0]), CoordinatedIndex(Index(1, PID('u', 'y')), [0.0, 0.0], [0.0, 0.0])),
#         Operator(+1.0, CoordinatedIndex(Index(2, PID('u', 'y')), [0.5, 0.5], [0.0, 0.0]), CoordinatedIndex(Index(2, PID('u', 'y')), [0.5, 0.5], [0.0, 0.0])),
#         Operator(+1.0, CoordinatedIndex(Index(1, PID('u', 'y')), [0.0, 0.0], [0.0, 0.0]), CoordinatedIndex(Index(1, PID('u', 'x')), [0.0, 0.0], [0.0, 0.0])),
#         Operator(+1.0, CoordinatedIndex(Index(2, PID('u', 'y')), [0.5, 0.5], [0.0, 0.0]), CoordinatedIndex(Index(2, PID('u', 'x')), [0.5, 0.5], [0.0, 0.0])),
#         Operator(+1.0, CoordinatedIndex(Index(2, PID('u', 'x')), [0.5, 0.5], [0.0, 0.0]), CoordinatedIndex(Index(2, PID('u', 'x')), [0.5, 0.5], [0.0, 0.0])),
#         Operator(-1.0, CoordinatedIndex(Index(2, PID('u', 'y')), [0.5, 0.5], [0.0, 0.0]), CoordinatedIndex(Index(1, PID('u', 'x')), [0.0, 0.0], [0.0, 0.0])),
#         Operator(+1.0, CoordinatedIndex(Index(1, PID('u', 'x')), [0.0, 0.0], [0.0, 0.0]), CoordinatedIndex(Index(1, PID('u', 'x')), [0.0, 0.0], [0.0, 0.0])),
#         Operator(+1.0, CoordinatedIndex(Index(2, PID('u', 'x')), [0.5, 0.5], [0.0, 0.0]), CoordinatedIndex(Index(2, PID('u', 'y')), [0.5, 0.5], [0.0, 0.0]))
#     )
#     @test expand(term, bond, hilbert) ≈ operators
# end

# @testset "Elastic" begin
#     term = Elastic(:V, 2.0, 1, MatrixCoupling(:, PID, [0 1; 1 0]))

#     bond = Bond(1, Point(1, [0.0, 0.0], [0.0, 0.0]), Point(2, [0.5, 0.0], [0.0, 0.0]))
#     hilbert = Hilbert(site=>Phonon(2) for site=1:2)
#     operators = Operators(
#         Operator(1.0, CoordinatedIndex(Index(1, PID('u', 'x')), [0.0, 0.0], [0.0, 0.0]), CoordinatedIndex(Index(2, PID('u', 'y')), [0.5, 0.0], [0.0, 0.0])),
#         Operator(1.0, CoordinatedIndex(Index(1, PID('u', 'y')), [0.0, 0.0], [0.0, 0.0]), CoordinatedIndex(Index(2, PID('u', 'x')), [0.5, 0.0], [0.0, 0.0])),
#         Operator(1.0, CoordinatedIndex(Index(2, PID('u', 'x')), [0.5, 0.0], [0.0, 0.0]), CoordinatedIndex(Index(1, PID('u', 'y')), [0.0, 0.0], [0.0, 0.0])),
#         Operator(1.0, CoordinatedIndex(Index(2, PID('u', 'y')), [0.5, 0.0], [0.0, 0.0]), CoordinatedIndex(Index(1, PID('u', 'x')), [0.0, 0.0], [0.0, 0.0])),
#     )
#     @test expand(term, bond, hilbert) == operators
# end
