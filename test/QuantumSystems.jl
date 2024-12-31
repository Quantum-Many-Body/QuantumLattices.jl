using LaTeXStrings: latexstring
using QuantumLattices: expand, kind, permute, rank
using QuantumLattices.DegreesOfFreedom: ˢᵗ, ⁿᵈ, CompositeIndex, CoordinatedIndex, Coupling, Hilbert, Index, InternalIndex, MatrixCoupling, allequalfields, internalindextype, isdefinite, patternrule, statistics, @pattern
using QuantumLattices.QuantumOperators: Operator, OperatorIndex, Operators, latexname, matrix, script
using QuantumLattices.QuantumSystems
using QuantumLattices.Spatials: Bond, Lattice, Neighbors, Point, azimuthd, bonds, rcoordinate, icoordinate
using QuantumLattices.Toolkit: DuplicatePermutations, shape
using SparseArrays: SparseMatrixCSC
using StaticArrays: SVector

@testset "FockIndex" begin
    index = 𝕔(1, 1//2, 1)
    @test FockIndex{:f, Colon, Colon, Colon}(1, 1//2, 1) == index
    @test statistics(index) == statistics(typeof(index)) == :f
    @test isdefinite(index) == isdefinite(typeof(index)) == true
    @test index' == replace(index, nambu=2)
    @test isequal(index'', replace(index, nambu=1))
    @test hash(index) == hash((:f, 1, 1//2, 1))
    @test string(index) == "𝕔(1, 1//2, 1)"
    @test isannihilation(index) && isannihilation(𝕔(1, 1, 1//2, 1)) && isannihilation(𝕔(1, 1, 1//2, 1, [0.0], [0.0]))
    @test !iscreation(index) && !iscreation(𝕔(1, 1, 1//2, 1)) && !iscreation(𝕔(1, 1, 1//2, 1, [0.0], [0.0]))

    index = 𝕓(1, -1//2, 2)
    @test FockIndex{:b, Colon, Colon, Colon}(1, -1//2, 2) == index
    @test statistics(index) == statistics(typeof(index)) == :b
    @test isdefinite(index) == isdefinite(typeof(index)) == true
    @test hash(index) == hash((:b, 1, -1//2, 2))
    @test string(index) == "𝕓(1, -1//2, 2)"
    @test !isannihilation(index) && !isannihilation(𝕓(1, 1, -1//2, 2)) && !isannihilation(𝕓(1, 1, -1//2, 2, [0.0], [0.0]))
    @test iscreation(index) && iscreation(𝕓(1, 1, -1//2, 2)) && iscreation(𝕓(1, 1, -1//2, 2, [0.0], [0.0]))

    index = 𝕕(1, :α, :)
    @test FockIndex{:, Colon, Colon, Colon}(1, :α, :) == index
    @test statistics(index) == statistics(typeof(index)) == Colon()
    @test isdefinite(index) == isdefinite(typeof(index)) == false
    @test index == FockIndex(1, :α, :)
    @test hash(index) == hash((:, 1, :α, :))
    @test string(index) == "𝕕(1, α, :)"
    @test !isannihilation(index) && !isannihilation(𝕕(1, 1, :α, :)) && !isannihilation(𝕕(1, 1, :α, :, [0.0], [0.0]))
    @test !iscreation(index) && !iscreation(𝕕(1, 1, :α, :)) && !iscreation(𝕕(1, 1, :α, :, [0.0], [0.0]))

    @test 𝕔(1, 1//2, 1) ≠ 𝕓(1, 1//2, 1)
    @test !isequal(𝕔(1, 1//2, 1), 𝕓(1, 1//2, 1))

    @test statistics(FockIndex) == statistics(Index{<:FockIndex}) == statistics(CoordinatedIndex{<:Index{<:FockIndex}}) == Symbol(":")
    @test allequalfields(FockIndex) == (:orbital, :spin)
    @test isdefinite(FockIndex{:, Int, Rational{Int}, Int})
    @test !isdefinite(FockIndex{:f, Symbol, typeof(:), Int})
    @test internalindextype(FockIndex, Int, typeof(:), Int) == FockIndex{:, Int, typeof(:), Int}
    @test internalindextype(FockIndex{:f}, typeof(:), Symbol, Symbol) == FockIndex{:f, typeof(:), Symbol, Symbol}
    @test internalindextype(FockIndex{:b}, typeof(:), Symbol, Symbol) == FockIndex{:b, typeof(:), Symbol, Symbol}

    @test OperatorIndex[FockIndex{:f}] == OperatorIndex[Index{<:FockIndex{:f}}] == OperatorIndex[CoordinatedIndex{<:Index{<:FockIndex{:f}}}] == 𝕔
    @test OperatorIndex[FockIndex{:b}] == OperatorIndex[Index{<:FockIndex{:b}}] == OperatorIndex[CoordinatedIndex{<:Index{<:FockIndex{:b}}}] == 𝕓
    @test OperatorIndex[FockIndex{:}] == OperatorIndex[Index{<:FockIndex{:}}] == OperatorIndex[CoordinatedIndex{<:Index{<:FockIndex{:}}}] == 𝕕
    @test OperatorIndex[FockIndex] == OperatorIndex[Index{<:FockIndex}] == OperatorIndex[CoordinatedIndex{<:Index{<:FockIndex}}] == 𝕕
    @test OperatorIndex[𝕔] == FockIndex{:f}
    @test OperatorIndex[𝕓] == FockIndex{:b}
    @test OperatorIndex[𝕕] == FockIndex{:}

    patternrule((:, :, :, :), Val(:), FockIndex, Val(:nambu)) == (2, 1, 2, 1)
end

@testset "Fock latex" begin
    @test script(𝕔(1, 2, 1//2, 1), Val(:site)) == script(𝕔(1, 2, 1//2, 1, [0.0], [0.0]), Val(:site)) == "1"
    @test script(𝕔(2, 1//2, 1), Val(:orbital)) == script(𝕔(1, 2, 1//2, 1), Val(:orbital)) == script(𝕔(1, 2, 1//2, 1, [0.0], [0.0]), Val(:orbital)) == "2"
    @test script(𝕔(2, 3//2, 1), Val(:spin)) == script(𝕔(1, 2, 3//2, 1), Val(:spin)) == script(𝕔(1, 2, 3//2, 1, [0.0], [0.0]), Val(:spin)) == "3//2"
    @test script(𝕔(2, 1//2, 1), Val(:spinsym)) == script(𝕔(1, 2, 1//2, 1), Val(:spinsym)) == script(𝕔(1, 2, 1//2, 1, [0.0], [0.0]), Val(:spinsym)) == "↑"
    @test script(𝕔(2, -1//2, 1), Val(:spinsym)) == script(𝕔(1, 2, -1//2, 1), Val(:spinsym)) == script(𝕔(1, 2, -1//2, 1, [0.0], [0.0]), Val(:spinsym)) == "↓"
    @test script(𝕔(2, 3//2, 1), Val(:nambu)) == script(𝕔(1, 2, 3//2, 1), Val(:nambu)) == script(𝕔(1, 2, 3//2, 1, [0.0], [0.0]), Val(:nambu)) == ""
    @test script(𝕔(2, 3//2, 2), Val(:nambu)) == script(𝕔(1, 2, 3//2, 2), Val(:nambu)) == script(𝕔(1, 2, 3//2, 2, [0.0], [0.0]), Val(:nambu)) == "\\dagger"

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

    @test shape(Fock{:f}(3, 2), 𝕔(2, 1//2, 1)) == (2:2, 2:2, 1:1)
    @test shape(Fock{:f}(3, 2), 𝕔(1, :, 2)) ==(1:1, 1:2, 2:2) 
    @test shape(Fock{:f}(3, 2), 𝕔(:, -1//2, 1)) == (1:3, 1:1, 1:1)
    @test shape(Fock{:f}(3, 2), 𝕔(:, :, 2)) == (1:3, 1:2, 2:2)
end

@testset "angle" begin
    @test angle(𝕔(1, 1, 1//2, 1, [0.0, 0.0], [1.0, 2.0]), [[1.0, 0.0], [0.0, 1.0]], [0.1, 0.0]) ≈ 2pi*0.1
    @test angle(𝕔(1, 1, 1//2, 2, [0.0, 0.0], [1.0, 2.0]), [[1.0, 0.0], [0.0, 1.0]], [0.0, 0.2]) ≈ -2pi*0.4
end

@testset "Fock Operator" begin
    id₁ = 𝕔(2, 1, -1//2, 2, SVector(0.5, 0.0), SVector(0.0, 0.0))
    id₂ = 𝕔(2, 1, -1//2, 1, SVector(0.5, 0.0), SVector(0.0, 0.0))
    id₃ = 𝕔(1, 1, 1//2, 2, SVector(0.0, 0.0), SVector(0.0, 0.0))
    id₄ = 𝕔(1, 1, 1//2, 1, SVector(0.0, 0.0), SVector(0.0, 0.0))
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

    permute(𝕔(1, -1//2, 2), 𝕓(1, -1//2, 2)) == Operator(1, 𝕓(1, -1//2, 2), 𝕔(1, -1//2, 2))
    permute(𝕓(1, -1//2, 2), 𝕔(1, -1//2, 2)) == Operator(1, 𝕔(1, -1//2, 2), 𝕓(1, -1//2, 2))
end

@testset "Fock Coupling" begin
    @test collect(MatrixCoupling(:, FockIndex, :, :, :)) == collect(MatrixCoupling(𝕕, :, :, :, :)) == collect(𝕕⁺𝕕(:, :, :, :)) == [Coupling(𝕕(:, :, :, :), 𝕕(:, :, :, :))]
    @test collect(MatrixCoupling(:, FockIndex{:}, σ"+", σ"-", :)) == [Coupling(𝕕(:, 1, -1//2, :), 𝕕(:, 2, 1//2, :))]
    @test collect(MatrixCoupling((1ˢᵗ, 2ⁿᵈ), FockIndex{:f}, :, σ"y", σ"z")) == collect(MatrixCoupling(𝕔, (1ˢᵗ, 2ⁿᵈ), :, σ"y", σ"z")) == collect(𝕔⁺𝕔((1ˢᵗ, 2ⁿᵈ), :, σ"y", σ"z")) == [
        Coupling(+1im, 𝕔(1ˢᵗ, :, -1//2, 1), 𝕔(2ⁿᵈ, :, 1//2, 2)), Coupling(-1im, 𝕔(1ˢᵗ, :, 1//2, 1), 𝕔(2ⁿᵈ, :, -1//2, 2)),
        Coupling(-1im, 𝕔(1ˢᵗ, :, -1//2, 2), 𝕔(2ⁿᵈ, :, 1//2, 1)), Coupling(+1im, 𝕔(1ˢᵗ, :, 1//2, 2), 𝕔(2ⁿᵈ, :, -1//2, 1))
    ]
    @test collect(MatrixCoupling((1ˢᵗ, 2ⁿᵈ), FockIndex{:b}, σ"x", :, σ"0")) == collect(MatrixCoupling(𝕓, (1ˢᵗ, 2ⁿᵈ), σ"x", :, σ"0")) == collect(𝕓⁺𝕓((1ˢᵗ, 2ⁿᵈ), σ"x", :, σ"0")) == [
        Coupling(𝕓(1ˢᵗ, 2, :, 1), 𝕓(2ⁿᵈ, 1, :, 2)), Coupling(𝕓(1ˢᵗ, 1, :, 1), 𝕓(2ⁿᵈ, 2, :, 2)),
        Coupling(𝕓(1ˢᵗ, 2, :, 2), 𝕓(2ⁿᵈ, 1, :, 1)), Coupling(𝕓(1ˢᵗ, 1, :, 2), 𝕓(2ⁿᵈ, 2, :, 1))
    ]

    fc = Coupling(2.0, (1ˢᵗ, 2ⁿᵈ), FockIndex, (1, 2), :, (2, 1))
    bond = Bond(1, Point(1, SVector(0.0), SVector(0.0)), Point(2, SVector(0.5), SVector(0.0)))
    hilbert = Hilbert(site=>Fock{:f}(2, 2) for site=1:2)
    ex = expand(fc, Val(:Hopping), bond, hilbert)
    @test collect(ex) == [
        Operator(2.0, 𝕔(1, 1, -1//2, 2, SVector(0.0), SVector(0.0)), 𝕔(2, 2, -1//2, 1, SVector(0.5), SVector(0.0))),
        Operator(2.0, 𝕔(1, 1, +1//2, 2, SVector(0.0), SVector(0.0)), 𝕔(2, 2, +1//2, 1, SVector(0.5), SVector(0.0)))
    ]

    fc = Coupling(2.0, 𝕕, (1ˢᵗ, 1ˢᵗ, 1ˢᵗ, 1ˢᵗ), :, (1//2, 1//2, -1//2, -1//2), (2, 1, 2, 1))
    point = Point(1, SVector(0.0), SVector(0.0))
    hilbert = Hilbert(point.site=>Fock{:b}(2, 2))
    ex = expand(fc, Val(:term), Bond(point), hilbert)
    @test collect(ex) == [
        Operator(2.0, 𝕓(1, 1, +1//2, 2, SVector(0.0), SVector(0.0)), 𝕓(1, 1, +1//2, 1, SVector(0.0), SVector(0.0)), 𝕓(1, 1, -1//2, 2, SVector(0.0), SVector(0.0)), 𝕓(1, 1, -1//2, 1, SVector(0.0), SVector(0.0))),
        Operator(2.0, 𝕓(1, 2, +1//2, 2, SVector(0.0), SVector(0.0)), 𝕓(1, 2, +1//2, 1, SVector(0.0), SVector(0.0)), 𝕓(1, 2, -1//2, 2, SVector(0.0), SVector(0.0)), 𝕓(1, 2, -1//2, 1, SVector(0.0), SVector(0.0)))
    ]

    fc = Coupling(2.0, @pattern(𝕕(:, α, 1//2, 2), 𝕕(:, α, -1//2, 2), 𝕕(:, β, -1//2, 1), 𝕕(:, β, 1//2, 1); constraint=α<β))
    point = Point(1, SVector(0.5), SVector(0.0))
    hilbert = Hilbert(point.site=>Fock{:f}(3, 2))
    ex = expand(fc, Val(:term), Bond(point), hilbert)
    @test collect(ex) == [
        Operator(2.0, 𝕔(1, 1, +1//2, 2, SVector(0.5), SVector(0.0)), 𝕔(1, 1, -1//2, 2, SVector(0.5), SVector(0.0)), 𝕔(1, 2, -1//2, 1, SVector(0.5), SVector(0.0)), 𝕔(1, 2, +1//2, 1, SVector(0.5), SVector(0.0))),
        Operator(2.0, 𝕔(1, 1, +1//2, 2, SVector(0.5), SVector(0.0)), 𝕔(1, 1, -1//2, 2, SVector(0.5), SVector(0.0)), 𝕔(1, 3, -1//2, 1, SVector(0.5), SVector(0.0)), 𝕔(1, 3, +1//2, 1, SVector(0.5), SVector(0.0))),
        Operator(2.0, 𝕔(1, 2, +1//2, 2, SVector(0.5), SVector(0.0)), 𝕔(1, 2, -1//2, 2, SVector(0.5), SVector(0.0)), 𝕔(1, 3, -1//2, 1, SVector(0.5), SVector(0.0)), 𝕔(1, 3, +1//2, 1, SVector(0.5), SVector(0.0)))
    ]

    fc₁ = Coupling(+1.0, 𝕕, :, :, (+1//2, +1//2), (2, 1))
    fc₂ = Coupling(-1.0, 𝕕, :, :, (-1//2, -1//2), (2, 1))
    point = Point(1, SVector(0.0), SVector(0.0))
    hilbert = Hilbert(point.site=>Fock{:f}(2, 2))
    ex = expand(fc₁*fc₂, Val(:term), Bond(point), hilbert)
    @test collect(ex) == [
        Operator(-1.0, 𝕔(1, 1, +1//2, 2, SVector(0.0), SVector(0.0)), 𝕔(1, 1, +1//2, 1, SVector(0.0), SVector(0.0)), 𝕔(1, 1, -1//2, 2, SVector(0.0), SVector(0.0)), 𝕔(1, 1, -1//2, 1, SVector(0.0), SVector(0.0))),
        Operator(-1.0, 𝕔(1, 2, +1//2, 2, SVector(0.0), SVector(0.0)), 𝕔(1, 2, +1//2, 1, SVector(0.0), SVector(0.0)), 𝕔(1, 1, -1//2, 2, SVector(0.0), SVector(0.0)), 𝕔(1, 1, -1//2, 1, SVector(0.0), SVector(0.0))),
        Operator(-1.0, 𝕔(1, 1, +1//2, 2, SVector(0.0), SVector(0.0)), 𝕔(1, 1, +1//2, 1, SVector(0.0), SVector(0.0)), 𝕔(1, 2, -1//2, 2, SVector(0.0), SVector(0.0)), 𝕔(1, 2, -1//2, 1, SVector(0.0), SVector(0.0))),
        Operator(-1.0, 𝕔(1, 2, +1//2, 2, SVector(0.0), SVector(0.0)), 𝕔(1, 2, +1//2, 1, SVector(0.0), SVector(0.0)), 𝕔(1, 2, -1//2, 2, SVector(0.0), SVector(0.0)), 𝕔(1, 2, -1//2, 1, SVector(0.0), SVector(0.0)))
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

    term = Onsite(:mu, 1.5, 𝕕⁺𝕕(:, σ"z", σ"x", :))
    operators = Operators(
        Operator(-0.75, 𝕔(1, 2, +1//2, 2, [0.5, 0.5], [0.0, 0.0]), 𝕔(1, 2, -1//2, 1, [0.5, 0.5], [0.0, 0.0])),
        Operator(+0.75, 𝕔(1, 1, -1//2, 2, [0.5, 0.5], [0.0, 0.0]), 𝕔(1, 1, +1//2, 1, [0.5, 0.5], [0.0, 0.0])),
        Operator(-0.75, 𝕔(1, 2, -1//2, 2, [0.5, 0.5], [0.0, 0.0]), 𝕔(1, 2, +1//2, 1, [0.5, 0.5], [0.0, 0.0])),
        Operator(+0.75, 𝕔(1, 1, +1//2, 2, [0.5, 0.5], [0.0, 0.0]), 𝕔(1, 1, -1//2, 1, [0.5, 0.5], [0.0, 0.0]))
    )
    @test expand(term, bond, hilbert, half=true) == operators
    @test expand(term, bond, hilbert, half=false) == operators*2

    term = Onsite(:mu, 1.5, 𝕕⁺𝕕(:, σ"z", σ"z", :))
    operators = Operators(
        Operator(+0.75, 𝕔(1, 2, -1//2, 2, [0.5, 0.5], [0.0, 0.0]), 𝕔(1, 2, -1//2, 1, [0.5, 0.5], [0.0, 0.0])),
        Operator(+0.75, 𝕔(1, 1, +1//2, 2, [0.5, 0.5], [0.0, 0.0]), 𝕔(1, 1, +1//2, 1, [0.5, 0.5], [0.0, 0.0])),
        Operator(-0.75, 𝕔(1, 2, +1//2, 2, [0.5, 0.5], [0.0, 0.0]), 𝕔(1, 2, +1//2, 1, [0.5, 0.5], [0.0, 0.0])),
        Operator(-0.75, 𝕔(1, 1, -1//2, 2, [0.5, 0.5], [0.0, 0.0]), 𝕔(1, 1, -1//2, 1, [0.5, 0.5], [0.0, 0.0]))
    )
    @test expand(term, bond, hilbert, half=true) == operators
    @test expand(term, bond, hilbert, half=false) == operators*2
end

@testset "Hopping" begin
    bond = Bond(1, Point(2, (0.0, 0.0), (0.0, 0.0)), Point(1, (0.5, 0.5), (0.0, 0.0)))
    hilbert = Hilbert(site=>Fock{:f}(2, 2) for site=1:2)
    term = Hopping(:t, 1.5, 1)
    operators = Operators(
        Operator(1.5, 𝕔(2, 2, +1//2, 2, [0.0, 0.0], [0.0, 0.0]), 𝕔(1, 2, +1//2, 1, [0.5, 0.5], [0.0, 0.0])),
        Operator(1.5, 𝕔(2, 2, -1//2, 2, [0.0, 0.0], [0.0, 0.0]), 𝕔(1, 2, -1//2, 1, [0.5, 0.5], [0.0, 0.0])),
        Operator(1.5, 𝕔(2, 1, -1//2, 2, [0.0, 0.0], [0.0, 0.0]), 𝕔(1, 1, -1//2, 1, [0.5, 0.5], [0.0, 0.0])),
        Operator(1.5, 𝕔(2, 1, +1//2, 2, [0.0, 0.0], [0.0, 0.0]), 𝕔(1, 1, +1//2, 1, [0.5, 0.5], [0.0, 0.0]))
    )
    @test expand(term, bond, hilbert, half=true) == operators
    @test expand(term, bond, hilbert, half=false) == operators+operators'
end

@testset "Pairing" begin
    bond = Bond(1, Point(2, (0.0, 0.0), (0.0, 0.0)), Point(1, (0.5, 0.5), (0.0, 0.0)))
    hilbert = Hilbert(site=>Fock{:f}(1, 1) for site=1:2)
    term = Pairing(:Δ, 1.5, 1, Coupling{2}(𝕕, :, :, :, :); amplitude=bond->(bond|>rcoordinate|>azimuthd ≈ 45 ? 1 : -1))
    operators = Operators(
        Operator(+1.5, 𝕔(2, 1, 0, 1, [0.0, 0.0], [0.0, 0.0]), 𝕔(1, 1, 0, 1, [0.5, 0.5], [0.0, 0.0])),
        Operator(-1.5, 𝕔(1, 1, 0, 1, [0.5, 0.5], [0.0, 0.0]), 𝕔(2, 1, 0, 1, [0.0, 0.0], [0.0, 0.0]))
    )
    @test expand(term, bond, hilbert, half=true) == operators
    @test expand(term, bond, hilbert, half=false) == operators+operators'

    point = Point(1, (0.5, 0.5), (0.0, 0.0))
    hilbert = Hilbert(point.site=>Fock{:f}(1, 2))
    term = Pairing(:Δ, 1.5, 0, 𝕕⁺𝕕(:, :, [0 -1; 1 0], :))
    operators = Operators(
        Operator(-1.5, 𝕔(1, 1, +1//2, 1, [0.5, 0.5], [0.0, 0.0]), 𝕔(1, 1, -1//2, 1, [0.5, 0.5], [0.0, 0.0])),
        Operator(+1.5, 𝕔(1, 1, -1//2, 1, [0.5, 0.5], [0.0, 0.0]), 𝕔(1, 1, +1//2, 1, [0.5, 0.5], [0.0, 0.0]))
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
        Operator(1.25, 𝕔(1, 1, +1//2, 2, [0.5, 0.5], [0.0, 0.0]), 𝕔(1, 1, +1//2, 1, [0.5, 0.5], [0.0, 0.0]), 𝕔(1, 1, -1//2, 2, [0.5, 0.5], [0.0, 0.0]), 𝕔(1, 1, -1//2, 1, [0.5, 0.5], [0.0, 0.0])),
        Operator(1.25, 𝕔(1, 2, +1//2, 2, [0.5, 0.5], [0.0, 0.0]), 𝕔(1, 2, +1//2, 1, [0.5, 0.5], [0.0, 0.0]), 𝕔(1, 2, -1//2, 2, [0.5, 0.5], [0.0, 0.0]), 𝕔(1, 2, -1//2, 1, [0.5, 0.5], [0.0, 0.0]))
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
        Operator(1.25, 𝕔(1, 1, +1//2, 2, [0.5, 0.5], [0.0, 0.0]), 𝕔(1, 1, +1//2, 1, [0.5, 0.5], [0.0, 0.0]), 𝕔(1, 2, -1//2, 2, [0.5, 0.5], [0.0, 0.0]), 𝕔(1, 2, -1//2, 1, [0.5, 0.5], [0.0, 0.0])),
        Operator(1.25, 𝕔(1, 1, -1//2, 2, [0.5, 0.5], [0.0, 0.0]), 𝕔(1, 1, -1//2, 1, [0.5, 0.5], [0.0, 0.0]), 𝕔(1, 2, +1//2, 2, [0.5, 0.5], [0.0, 0.0]), 𝕔(1, 2, +1//2, 1, [0.5, 0.5], [0.0, 0.0]))
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
        Operator(1.25, 𝕔(1, 1, -1//2, 2, [0.5, 0.5], [0.0, 0.0]), 𝕔(1, 1, -1//2, 1, [0.5, 0.5], [0.0, 0.0]), 𝕔(1, 2, -1//2, 2, [0.5, 0.5], [0.0, 0.0]), 𝕔(1, 2, -1//2, 1, [0.5, 0.5], [0.0, 0.0])),
        Operator(1.25, 𝕔(1, 1, 1//2, 2, [0.5, 0.5], [0.0, 0.0]), 𝕔(1, 1, 1//2, 1, [0.5, 0.5], [0.0, 0.0]), 𝕔(1, 2, 1//2, 2, [0.5, 0.5], [0.0, 0.0]), 𝕔(1, 2, 1//2, 1, [0.5, 0.5], [0.0, 0.0]))
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
        Operator(2.5, 𝕔(1, 1, +1//2, 2, [0.5, 0.5], [0.0, 0.0]), 𝕔(1, 2, -1//2, 2, [0.5, 0.5], [0.0, 0.0]), 𝕔(1, 1, -1//2, 1, [0.5, 0.5], [0.0, 0.0]), 𝕔(1, 2, +1//2, 1, [0.5, 0.5], [0.0, 0.0]))
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
        Operator(2.5, 𝕔(1, 1, +1//2, 2, [0.5, 0.5], [0.0, 0.0]), 𝕔(1, 1, -1//2, 2, [0.5, 0.5], [0.0, 0.0]), 𝕔(1, 2, -1//2, 1, [0.5, 0.5], [0.0, 0.0]), 𝕔(1, 2, +1//2, 1, [0.5, 0.5], [0.0, 0.0]))
    )
    @test expand(term, bond, hilbert, half=true) == operators
    @test expand(term, bond, hilbert, half=false) == operators+operators'
end

@testset "Coulomb" begin
    bond = Bond(1, Point(2, (0.0, 0.0), (0.0, 0.0)), Point(1, (0.5, 0.5), (0.0, 0.0)))
    hilbert = Hilbert(site=>Fock{:f}(1, 2) for site=1:2)

    term = Coulomb(:V, 2.5, 1, 𝕕⁺𝕕(:, :, σ"z", :)^2)
    operators = Operators(
        Operator(-1.25, 𝕔(2, 1, -1//2, 2, [0.0, 0.0], [0.0, 0.0]), 𝕔(2, 1, -1//2, 1, [0.0, 0.0], [0.0, 0.0]), 𝕔(1, 1, +1//2, 2, [0.5, 0.5], [0.0, 0.0]), 𝕔(1, 1, +1//2, 1, [0.5, 0.5], [0.0, 0.0])),
        Operator(+1.25, 𝕔(2, 1, -1//2, 2, [0.0, 0.0], [0.0, 0.0]), 𝕔(2, 1, -1//2, 1, [0.0, 0.0], [0.0, 0.0]), 𝕔(1, 1, -1//2, 2, [0.5, 0.5], [0.0, 0.0]), 𝕔(1, 1, -1//2, 1, [0.5, 0.5], [0.0, 0.0])),
        Operator(-1.25, 𝕔(2, 1, +1//2, 2, [0.0, 0.0], [0.0, 0.0]), 𝕔(2, 1, +1//2, 1, [0.0, 0.0], [0.0, 0.0]), 𝕔(1, 1, -1//2, 2, [0.5, 0.5], [0.0, 0.0]), 𝕔(1, 1, -1//2, 1, [0.5, 0.5], [0.0, 0.0])),
        Operator(+1.25, 𝕔(2, 1, +1//2, 2, [0.0, 0.0], [0.0, 0.0]), 𝕔(2, 1, +1//2, 1, [0.0, 0.0], [0.0, 0.0]), 𝕔(1, 1, +1//2, 2, [0.5, 0.5], [0.0, 0.0]), 𝕔(1, 1, +1//2, 1, [0.5, 0.5], [0.0, 0.0]))
    )
    @test expand(term, bond, hilbert, half=true) == operators
    @test expand(term, bond, hilbert, half=false) == operators*2

    term = Coulomb(:V, 2.5, 1, 𝕕⁺𝕕(:, :, σ"x", :)*𝕕⁺𝕕(:, :, σ"z", :))
    operators = Operators(
        Operator(-1.25, 𝕔(2, 1, +1//2, 2, [0.0, 0.0], [0.0, 0.0]), 𝕔(2, 1, -1//2, 1, [0.0, 0.0], [0.0, 0.0]), 𝕔(1, 1, -1//2, 2, [0.5, 0.5], [0.0, 0.0]), 𝕔(1, 1, -1//2, 1, [0.5, 0.5], [0.0, 0.0])),
        Operator(+1.25, 𝕔(2, 1, -1//2, 2, [0.0, 0.0], [0.0, 0.0]), 𝕔(2, 1, +1//2, 1, [0.0, 0.0], [0.0, 0.0]), 𝕔(1, 1, +1//2, 2, [0.5, 0.5], [0.0, 0.0]), 𝕔(1, 1, +1//2, 1, [0.5, 0.5], [0.0, 0.0])),
        Operator(+1.25, 𝕔(2, 1, +1//2, 2, [0.0, 0.0], [0.0, 0.0]), 𝕔(2, 1, -1//2, 1, [0.0, 0.0], [0.0, 0.0]), 𝕔(1, 1, +1//2, 2, [0.5, 0.5], [0.0, 0.0]), 𝕔(1, 1, +1//2, 1, [0.5, 0.5], [0.0, 0.0])),
        Operator(-1.25, 𝕔(2, 1, -1//2, 2, [0.0, 0.0], [0.0, 0.0]), 𝕔(2, 1, +1//2, 1, [0.0, 0.0], [0.0, 0.0]), 𝕔(1, 1, -1//2, 2, [0.5, 0.5], [0.0, 0.0]), 𝕔(1, 1, -1//2, 1, [0.5, 0.5], [0.0, 0.0]))
    )
    @test expand(term, bond, hilbert, half=true) == operators
    @test expand(term, bond, hilbert, half=false) == operators*2
end

@testset "SpinIndex" begin
    index = 𝕊{3//2}('x')
    @test statistics(index) == statistics(typeof(index)) == :b
    @test isdefinite(index) == isdefinite(typeof(index)) == true
    @test index == SpinIndex{3//2, Colon}('x')
    @test isequal(index, index')
    @test hash(index) == hash((3//2, 'x'))
    @test replace(index, tag='z') == 𝕊{3//2}('z')
    @test string(index) == "𝕊{3//2}('x')"
    @test totalspin(index) == totalspin(typeof(index)) == 3//2
    @test totalspin(𝕊{3//2}(:, 'x')) == totalspin(typeof(𝕊{3//2}(:, 'x'))) ==3//2
    @test totalspin(𝕊{3//2}(:, 'x', [0], [0])) == totalspin(typeof(𝕊{3//2}(:, 'x', [0], [0]))) == 3//2

    index = 𝕊('z')
    @test index == 𝕊{:}('z')
    @test string(index) == "𝕊('z')"
    @test string(𝕊{:}) == "𝕊"
    @test 𝕊(:, 'z') == Index(:, index)
    @test 𝕊(:, 'z', [0.0], [0.0]) == CoordinatedIndex(Index(:, index), [0.0], [0.0])

    @test 𝕊{1//2}('z') ≠ 𝕊{3//2}('z')
    @test !isequal(𝕊{1//2}('z'), 𝕊{3//2}('z'))

    @test allequal((totalspin(SpinIndex), totalspin(Index{<:SpinIndex}), totalspin(CoordinatedIndex{<:Index{<:SpinIndex}}), NaN))
    @test isdefinite(SpinIndex{:, Char})
    @test !isdefinite(SpinIndex{1//2, Symbol})
    @test !isdefinite(SpinIndex{1, Colon})
    @test internalindextype(SpinIndex, Char) == SpinIndex{:, Char}
    @test internalindextype(SpinIndex{1//2}, Symbol) == SpinIndex{1//2, Symbol}

    @test OperatorIndex[SpinIndex] == OperatorIndex[Index{<:SpinIndex}] == OperatorIndex[CoordinatedIndex{<:Index{<:SpinIndex}}] == 𝕊
    @test OperatorIndex[SpinIndex{:}] == OperatorIndex[Index{<:SpinIndex{:}}] == OperatorIndex[CoordinatedIndex{<:Index{<:SpinIndex{:}}}] == 𝕊
    @test OperatorIndex[SpinIndex{1//2}] == OperatorIndex[Index{<:SpinIndex{1//2}}] == OperatorIndex[CoordinatedIndex{<:Index{<:SpinIndex{1//2}}}] == 𝕊{1//2}
    @test OperatorIndex[𝕊] == SpinIndex{:}
    @test OperatorIndex[𝕊{1//2}] == SpinIndex{1//2}
end

@testset "matrix" begin
    @test isapprox(matrix(𝕊{1//2}('z')), [0.5 0.0; 0.0 -0.5])
    @test isapprox(matrix(𝕊{1//2}('x')), [0.0 0.5; 0.5 0.0])
    @test isapprox(matrix(𝕊{1//2}('y')), [0.0 -0.5im; 0.5im 0.0])
    @test isapprox(matrix(𝕊{1//2}('+')), [0.0 1.0; 0.0 0.0])
    @test isapprox(matrix(𝕊{1//2}('-')), [0.0 0.0; 1.0 0.0])

    @test isapprox(matrix(𝕊{1//2}(:, 'z')), [0.5 0.0; 0.0 -0.5])
    @test isapprox(matrix(𝕊{1//2}(:, 'z', [0], [0])), [0.5 0.0; 0.0 -0.5])

    @test isapprox(matrix(𝕊{1}('z')), [1.0 0.0 0.0; 0.0 0.0 0.0; 0.0 0.0 -1.0])
    @test isapprox(matrix(𝕊{1}('x')), [0.0 √2/2 0.0; √2/2 0.0 √2/2; 0.0 √2/2 0.0])
    @test isapprox(matrix(𝕊{1}('y')), [0.0 -√2im/2 0.0; √2im/2 0.0 -√2im/2; 0.0 √2im/2 0.0])
    @test isapprox(matrix(𝕊{1}('+')), [0.0 √2 0.0; 0.0 0.0 √2; 0.0 0.0 0.0])
    @test isapprox(matrix(𝕊{1}('-')), [0.0 0.0 0.0; √2 0.0 0.0; 0.0 √2 0.0])
end

@testset "Spin latex" begin
    @test script(𝕊{1//2}(1, 'z'), Val(:site)) == script(𝕊{1//2}(1, 'z', [0.0], [0.0]), Val(:site)) == "1"
    @test script(𝕊{1//2}('z'), Val(:tag)) == script(𝕊{1//2}(1, 'z'), Val(:tag)) == script(𝕊{1//2}(1, 'z', [0.0], [0.0]), Val(:tag)) == "z"
    @test script(𝕊{1//2}(:), Val(:tag)) == script(𝕊{1//2}(1, :), Val(:tag)) == script(𝕊{1//2}(1, :, [0.0], [0.0]), Val(:tag)) == ":"

    @test latexname(SpinIndex) == Symbol("SpinIndex")
    @test latexname(Index{<:SpinIndex}) == Symbol("Index{SpinIndex}")
    @test latexname(CompositeIndex{<:Index{<:SpinIndex}}) == Symbol("CompositeIndex{Index{SpinIndex}}")
end

@testset "Spin" begin
    @test eltype(Spin) == (SpinIndex{S, Char} where S)
    spin = Spin{1}()
    @test shape(spin) == (1:3,)
    @test convert(CartesianIndex, 𝕊{1}('z'), spin) == CartesianIndex(3)
    @test convert(SpinIndex, CartesianIndex(1), spin) == 𝕊{1}('x')
    @test summary(spin) == "3-element Spin{1}"
    @test string(spin) == "Spin{1}()"
    @test totalspin(spin) == totalspin(typeof(spin)) == 1
    @test collect(spin) == [𝕊{1}('x'), 𝕊{1}('y'), 𝕊{1}('z')]

    @test match(SpinIndex{:}, Spin{1//2}) == true
    @test match(SpinIndex{1//2}, Spin{1//2}) == true
    @test match(SpinIndex{1//2}, Spin{1}) == match(SpinIndex{1}, Spin{1//2}) == false

    @test shape(Spin{1}(), 𝕊{1}('z')) == (3:3,)
end

@testset "Spin operator" begin
    opt = Operator(1.0, 𝕊{1//2}(1, '+', [0.0, 0.0], [0.0, 0.0]), 𝕊{1//2}(1, '-', [0.0, 0.0], [0.0, 0.0]))
    @test opt' == Operator(1.0, 𝕊{1//2}(1, '+', [0.0, 0.0], [0.0, 0.0]), 𝕊{1//2}(1, '-', [0.0, 0.0], [0.0, 0.0]))
    @test latexstring(opt) == "S^{+}_{1}S^{-}_{1}"

    representation(opt::Operator) = opt.value * prod([matrix(opt[i]) for i = 1:rank(opt)])
    for S in (1//2, 1, 3//2)
        indexes = [𝕊{S}(1, tag, [0.0, 0.0], [0.0, 0.0]) for tag in ('x', 'y', 'z', '+', '-')]
        for (id₁, id₂) in DuplicatePermutations{2}(indexes)
            left = representation(Operator(1, id₁, id₂))
            right = sum([representation(opt) for opt in permute(id₁, id₂)])
            @test isapprox(left, right)
        end
    end
    id₁ = 𝕊{1//2}(1, 'z', [0.0, 0.0], [0.0, 0.0])
    id₂ = 𝕊{1//2}(2, 'z', [0.0, 0.0], [0.0, 0.0])
    @test permute(id₁, id₂) == (Operator(1, id₂, id₁),)
end

@testset "Spin Coupling" begin
    @test collect(MatrixCoupling(:, SpinIndex, [1 0 0; 0 1 0; 0 0 1])) == collect(MatrixCoupling(𝕊, :, [1 0 0; 0 1 0; 0 0 1])) == collect(𝕊ᵀ𝕊(:, [1 0 0; 0 1 0; 0 0 1])) == [
        Coupling(𝕊, :, ('x', 'x')), Coupling(𝕊, :, ('y', 'y')), Coupling(𝕊, :, ('z', 'z'))
    ]

    sc = Coupling(2.0, 𝕊, (1ˢᵗ, 2ⁿᵈ), ('+', '-'))
    bond = Bond(1, Point(1, [0.0], [0.0]), Point(2, [0.5], [0.0]))
    hilbert = Hilbert(Spin{1}(), 2)
    ex = expand(sc, Val(:SpinTerm), bond, hilbert)
    @test collect(ex) == [Operator(2.0, 𝕊{1}(1, '+', [0.0], [0.0]), 𝕊{1}(2, '-', [0.5], [0.0]))]
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
    bond = Bond(Point(1, (0.5, 0.5), (0.0, 0.0)))
    hilbert = Hilbert(Spin{1//2}())
    term = SpinTerm(:h, 1.5, 0, Coupling(𝕊(1ˢᵗ, 'z')))
    operators = Operators(Operator(1.5, 𝕊{1//2}(1, 'z', [0.5, 0.5], [0.0, 0.0])))
    @test expand(term, bond, hilbert) == operators

    bond = Bond(1, Point(2, (0.5, 0.5), (0.0, 0.0)), Point(1, (0.0, 0.0), (0.0, 0.0)))
    hilbert = Hilbert(site=>Spin{1//2}() for site=1:2)
    term = SpinTerm(:J, 1.5, 1, 𝕊ᵀ𝕊(:, Heisenberg""))
    operators = Operators(
        Operator(1.5, 𝕊{1//2}(2, 'x', [0.5, 0.5], [0.0, 0.0]), 𝕊{1//2}(1, 'x', [0.0, 0.0], [0.0, 0.0])),
        Operator(1.5, 𝕊{1//2}(2, 'y', [0.5, 0.5], [0.0, 0.0]), 𝕊{1//2}(1, 'y', [0.0, 0.0], [0.0, 0.0])),
        Operator(1.5, 𝕊{1//2}(2, 'z', [0.5, 0.5], [0.0, 0.0]), 𝕊{1//2}(1, 'z', [0.0, 0.0], [0.0, 0.0])),
    )
    @test expand(term, bond, hilbert) == operators
end

@testset "Zeeman" begin
    point = Point(1, (0.5, 0.5), (0.0, 0.0))
    hilbert = Hilbert(point.site=>Spin{1//2}())
    term = Zeeman(:h, 1.5, 'x', 2)
    operators = Operators(Operator(3.0, 𝕊{1//2}(1, 'x', [0.5, 0.5], [0.0, 0.0])))
    @test expand(term, Bond(point), hilbert) == operators

    term = Zeeman(:h, 1.5, [1, 1, 1], 2)
    operators = Operators(
        Operator(√3, 𝕊{1//2}(1, 'x', [0.5, 0.5], [0.0, 0.0])),
        Operator(√3, 𝕊{1//2}(1, 'y', [0.5, 0.5], [0.0, 0.0])),
        Operator(√3, 𝕊{1//2}(1, 'z', [0.5, 0.5], [0.0, 0.0]))
    )
    @test expand(term, Bond(point), hilbert) ≈ operators

    term = Zeeman(:h, 1.5, [1, -1, 2], [1 0 0; 0 2 0; 0 0 3])
    operators = Operators(
        Operator(√6/4, 𝕊{1//2}(1, 'x', [0.5, 0.5], [0.0, 0.0])),
        Operator(-√6/2, 𝕊{1//2}(1, 'y', [0.5, 0.5], [0.0, 0.0])),
        Operator(3*√6/2, 𝕊{1//2}(1, 'z', [0.5, 0.5], [0.0, 0.0]))
    )
    @test expand(term, Bond(point), hilbert) ≈ operators
end

@testset "SingleIonAnisotropy" begin
    point = Point(1, (0.5, 0.5), (0.0, 0.0))
    hilbert = Hilbert(point.site=>Spin{1//2}())
    term = SingleIonAnisotropy(:A, 1.5, 'z')
    operators = Operators(Operator(1.5, 𝕊{1//2}(1, 'z', [0.5, 0.5], [0.0, 0.0]), 𝕊{1//2}(1, 'z', [0.5, 0.5], [0.0, 0.0])))
    @test expand(term, Bond(point), hilbert) == operators

    term = SingleIonAnisotropy(:A, 1.5, [1 0 0; 0 2 0; 0 0 3])
    operators = Operators(
        Operator(1.5, 𝕊{1//2}(1, 'x', [0.5, 0.5], [0.0, 0.0]), 𝕊{1//2}(1, 'x', [0.5, 0.5], [0.0, 0.0])),
        Operator(3.0, 𝕊{1//2}(1, 'y', [0.5, 0.5], [0.0, 0.0]), 𝕊{1//2}(1, 'y', [0.5, 0.5], [0.0, 0.0])),
        Operator(4.5, 𝕊{1//2}(1, 'z', [0.5, 0.5], [0.0, 0.0]), 𝕊{1//2}(1, 'z', [0.5, 0.5], [0.0, 0.0]))
    )
    @test expand(term, Bond(point), hilbert) == operators
end

@testset "Ising" begin
    bond = Bond(1, Point(1, (0.0, 0.0), (0.0, 0.0)), Point(2, (0.5, 0.5), (0.0, 0.0)))
    hilbert = Hilbert(Spin{1//2}(), 2)
    term = Ising(:J, 1.5, 1, 'x')
    operators = Operators(Operator(1.5, 𝕊{1//2}(1, 'x', [0.0, 0.0], [0.0, 0.0]), 𝕊{1//2}(2, 'x', [0.5, 0.5], [0.0, 0.0])))
    @test expand(term, bond, hilbert) == operators

    term = Ising(:J, 1.5, 1, 'y')
    operators = Operators(Operator(1.5, 𝕊{1//2}(1, 'y', [0.0, 0.0], [0.0, 0.0]), 𝕊{1//2}(2, 'y', [0.5, 0.5], [0.0, 0.0])))
    @test expand(term, bond, hilbert) == operators

    term = Ising(:J, 1.5, 1, 'z')
    operators = Operators(Operator(1.5, 𝕊{1//2}(1, 'z', [0.0, 0.0], [0.0, 0.0]), 𝕊{1//2}(2, 'z', [0.5, 0.5], [0.0, 0.0])))
    @test expand(term, bond, hilbert) == operators
end

@testset "Heisenberg" begin
    bond = Bond(1, Point(1, (0.0, 0.0), (0.0, 0.0)), Point(2, (0.5, 0.5), (0.0, 0.0)))
    hilbert = Hilbert(Spin{1//2}(), 2)
    term = Heisenberg(:J, 1.5, 1; form=:xyz)
    operators = Operators(
        Operator(1.5, 𝕊{1//2}(1, 'x', [0.0, 0.0], [0.0, 0.0]), 𝕊{1//2}(2, 'x', [0.5, 0.5], [0.0, 0.0])),
        Operator(1.5, 𝕊{1//2}(1, 'y', [0.0, 0.0], [0.0, 0.0]), 𝕊{1//2}(2, 'y', [0.5, 0.5], [0.0, 0.0])),
        Operator(1.5, 𝕊{1//2}(1, 'z', [0.0, 0.0], [0.0, 0.0]), 𝕊{1//2}(2, 'z', [0.5, 0.5], [0.0, 0.0]))
    )
    @test expand(term, bond, hilbert) == operators

    term = Heisenberg(:J, 1.5, 1; form=Symbol("+-z"))
    operators = Operators(
        Operator(0.75, 𝕊{1//2}(1, '+', [0.0, 0.0], [0.0, 0.0]), 𝕊{1//2}(2, '-', [0.5, 0.5], [0.0, 0.0])),
        Operator(0.75, 𝕊{1//2}(1, '-', [0.0, 0.0], [0.0, 0.0]), 𝕊{1//2}(2, '+', [0.5, 0.5], [0.0, 0.0])),
        Operator(1.5, 𝕊{1//2}(1, 'z', [0.0, 0.0], [0.0, 0.0]), 𝕊{1//2}(2, 'z', [0.5, 0.5], [0.0, 0.0]))
    )
    @test expand(term, bond, hilbert) == operators
end

@testset "Kitaev" begin
    lattice = Lattice((0.0, 0.0), (0.0, √3/3); vectors=[[1.0, 0.0], [0.5, √3/2]])
    bond₁, bond₂, bond₃ = bonds(lattice, Neighbors(1=>1/√3))
    hilbert = Hilbert(Spin{1//2}(), length(lattice))
    term = Kitaev(:K, 1.5, 1; x=[90], y=[210], z=[330], unit=:degree)

    operators = Operators(Operator(1.5, 𝕊{1//2}(2, 'x', -rcoordinate(bond₁), -icoordinate(bond₁)), 𝕊{1//2}(1, 'x', [0.0, 0.0], [0.0, 0.0])))
    @test expand(term, bond₁, hilbert) == operators
    @test expand(term, reverse(bond₁), hilbert) == operators'

    operators = Operators(Operator(1.5, 𝕊{1//2}(2, 'y', -rcoordinate(bond₂), -icoordinate(bond₂)), 𝕊{1//2}(1, 'y', [0.0, 0.0], [0.0, 0.0])))
    @test expand(term, bond₂, hilbert) == operators
    @test expand(term, reverse(bond₂), hilbert) == operators'

    operators = Operators(Operator(1.5, 𝕊{1//2}(2, 'z', -rcoordinate(bond₃), -icoordinate(bond₃)), 𝕊{1//2}(1, 'z', [0.0, 0.0], [0.0, 0.0])))
    @test expand(term, bond₃, hilbert) == operators
    @test expand(term, reverse(bond₃), hilbert) == operators'
end

@testset "Γ" begin
    lattice = Lattice((0.0, 0.0), (0.0, √3/3); vectors=[[1.0, 0.0], [0.5, √3/2]])
    bond₁, bond₂, bond₃ = bonds(lattice, Neighbors(1=>1/√3))
    hilbert = Hilbert(Spin{1//2}(), length(lattice))
    term = Γ(:Γ, 1.5, 1; x=[90], y=[210], z=[330], unit=:degree)

    operators = Operators(
        Operator(1.5, 𝕊{1//2}(2, 'y', -rcoordinate(bond₁), -icoordinate(bond₁)), 𝕊{1//2}(1, 'z', [0.0, 0.0], [0.0, 0.0])),
        Operator(1.5, 𝕊{1//2}(2, 'z', -rcoordinate(bond₁), -icoordinate(bond₁)), 𝕊{1//2}(1, 'y', [0.0, 0.0], [0.0, 0.0]))
    )
    @test expand(term, bond₁, hilbert) == operators
    @test expand(term, reverse(bond₁), hilbert) == operators'

    operators = Operators(
        Operator(1.5, 𝕊{1//2}(2, 'z', -rcoordinate(bond₂), -icoordinate(bond₂)), 𝕊{1//2}(1, 'x', [0.0, 0.0], [0.0, 0.0])),
        Operator(1.5, 𝕊{1//2}(2, 'x', -rcoordinate(bond₂), -icoordinate(bond₂)), 𝕊{1//2}(1, 'z', [0.0, 0.0], [0.0, 0.0]))
    )
    @test expand(term, bond₂, hilbert) == operators
    @test expand(term, reverse(bond₂), hilbert) == operators'

    operators = Operators(
        Operator(1.5, 𝕊{1//2}(2, 'y', -rcoordinate(bond₃), -icoordinate(bond₃)), 𝕊{1//2}(1, 'x', [0.0, 0.0], [0.0, 0.0])),
        Operator(1.5, 𝕊{1//2}(2, 'x', -rcoordinate(bond₃), -icoordinate(bond₃)), 𝕊{1//2}(1, 'y', [0.0, 0.0], [0.0, 0.0]))
    )
    @test expand(term, bond₃, hilbert) == operators
    @test expand(term, reverse(bond₃), hilbert) == operators'
end

@testset "Γ′" begin
    lattice = Lattice((0.0, 0.0), (0.0, √3/3); vectors=[[1.0, 0.0], [0.5, √3/2]])
    bond₁, bond₂, bond₃ = bonds(lattice, Neighbors(1=>1/√3))
    hilbert = Hilbert(Spin{1//2}(), length(lattice))
    term = Γ′(:Γ′, 1.5, 1; x=[90], y=[210], z=[330], unit=:degree)

    operators = Operators(
        Operator(1.5, 𝕊{1//2}(2, 'z', -rcoordinate(bond₁), -icoordinate(bond₁)), 𝕊{1//2}(1, 'x', [0.0, 0.0], [0.0, 0.0])),
        Operator(1.5, 𝕊{1//2}(2, 'x', -rcoordinate(bond₁), -icoordinate(bond₁)), 𝕊{1//2}(1, 'y', [0.0, 0.0], [0.0, 0.0])),
        Operator(1.5, 𝕊{1//2}(2, 'x', -rcoordinate(bond₁), -icoordinate(bond₁)), 𝕊{1//2}(1, 'z', [0.0, 0.0], [0.0, 0.0])),
        Operator(1.5, 𝕊{1//2}(2, 'y', -rcoordinate(bond₁), -icoordinate(bond₁)), 𝕊{1//2}(1, 'x', [0.0, 0.0], [0.0, 0.0]))
    )
    @test expand(term, bond₁, hilbert) == operators
    @test expand(term, reverse(bond₁), hilbert) == operators'

    operators = Operators(
        Operator(1.5, 𝕊{1//2}(2, 'z', -rcoordinate(bond₂), -icoordinate(bond₂)), 𝕊{1//2}(1, 'y', [0.0, 0.0], [0.0, 0.0])),
        Operator(1.5, 𝕊{1//2}(2, 'x', -rcoordinate(bond₂), -icoordinate(bond₂)), 𝕊{1//2}(1, 'y', [0.0, 0.0], [0.0, 0.0])),
        Operator(1.5, 𝕊{1//2}(2, 'y', -rcoordinate(bond₂), -icoordinate(bond₂)), 𝕊{1//2}(1, 'x', [0.0, 0.0], [0.0, 0.0])),
        Operator(1.5, 𝕊{1//2}(2, 'y', -rcoordinate(bond₂), -icoordinate(bond₂)), 𝕊{1//2}(1, 'z', [0.0, 0.0], [0.0, 0.0]))
    )
    @test expand(term, bond₂, hilbert) == operators
    @test expand(term, reverse(bond₂), hilbert) == operators'

    operators = Operators(
        Operator(1.5, 𝕊{1//2}(2, 'z', -rcoordinate(bond₃), -icoordinate(bond₃)), 𝕊{1//2}(1, 'x', [0.0, 0.0], [0.0, 0.0])),
        Operator(1.5, 𝕊{1//2}(2, 'y', -rcoordinate(bond₃), -icoordinate(bond₃)), 𝕊{1//2}(1, 'z', [0.0, 0.0], [0.0, 0.0])),
        Operator(1.5, 𝕊{1//2}(2, 'z', -rcoordinate(bond₃), -icoordinate(bond₃)), 𝕊{1//2}(1, 'y', [0.0, 0.0], [0.0, 0.0])),
        Operator(1.5, 𝕊{1//2}(2, 'x', -rcoordinate(bond₃), -icoordinate(bond₃)), 𝕊{1//2}(1, 'z', [0.0, 0.0], [0.0, 0.0]))
    )
    @test expand(term, bond₃, hilbert) == operators
    @test expand(term, reverse(bond₃), hilbert) == operators'
end

@testset "DM" begin
    lattice = Lattice((0.0, 0.0), (0.0, √3/3); vectors=[[1.0, 0.0], [0.5, √3/2]])
    bond₁, bond₂, bond₃ = bonds(lattice, Neighbors(1=>1/√3))
    hilbert = Hilbert(Spin{1//2}(), length(lattice))
    term = DM(:DM, 1.5, 1, [90]=>'x', [210]=>'y', [330]=>'z'; unit=:degree)

    operators = Operators(
        Operator(-1.5, 𝕊{1//2}(2, 'y', -rcoordinate(bond₁), -icoordinate(bond₁)), 𝕊{1//2}(1, 'z', [0.0, 0.0], [0.0, 0.0])),
        Operator(1.5, 𝕊{1//2}(2, 'z', -rcoordinate(bond₁), -icoordinate(bond₁)), 𝕊{1//2}(1, 'y', [0.0, 0.0], [0.0, 0.0]))
    )
    @test expand(term, bond₁, hilbert) == operators
    @test expand(term, reverse(bond₁), hilbert) == operators'

    operators = Operators(
        Operator(-1.5, 𝕊{1//2}(2, 'z', -rcoordinate(bond₂), -icoordinate(bond₂)), 𝕊{1//2}(1, 'x', [0.0, 0.0], [0.0, 0.0])),
        Operator(1.5, 𝕊{1//2}(2, 'x', -rcoordinate(bond₂), -icoordinate(bond₂)), 𝕊{1//2}(1, 'z', [0.0, 0.0], [0.0, 0.0]))
    )
    @test expand(term, bond₂, hilbert) == operators
    @test expand(term, reverse(bond₂), hilbert) == operators'

    operators = Operators(
        Operator(1.5, 𝕊{1//2}(2, 'y', -rcoordinate(bond₃), -icoordinate(bond₃)), 𝕊{1//2}(1, 'x', [0.0, 0.0], [0.0, 0.0])),
        Operator(-1.5, 𝕊{1//2}(2, 'x', -rcoordinate(bond₃), -icoordinate(bond₃)), 𝕊{1//2}(1, 'y', [0.0, 0.0], [0.0, 0.0]))
    )
    @test expand(term, bond₃, hilbert) == operators
    @test expand(term, reverse(bond₃), hilbert) == operators'
end

@testset "PhononIndex" begin
    index = 𝕦('x')
    @test statistics(index) == statistics(typeof(index)) == :b
    @test isdefinite(index) == isdefinite(typeof(index)) == true
    @test index == PhononIndex{:u, Colon}('x')
    @test isequal(index, index')
    @test hash(index) == hash((:u, 'x'))
    @test replace(index, direction='y') == 𝕦('y')
    @test string(index) == "𝕦('x')"
    @test kind(index) == kind(typeof(index)) == :u
    @test kind(𝕦(1, 'x')) == kind(typeof(𝕦(1, 'x'))) == :u
    @test kind(𝕦(1, 'x', [0.0], [0.0])) == kind(typeof(𝕦(1, 'x', [0.0], [0.0]))) == :u

    index = 𝕡('x')
    @test statistics(index) == statistics(typeof(index)) == :b
    @test isdefinite(index) == isdefinite(typeof(index)) == true
    @test index == PhononIndex{:p, Colon}('x')
    @test isequal(index, index')
    @test hash(index) == hash((:p, 'x'))
    @test replace(index, direction='y') == 𝕡('y')
    @test string(index) == "𝕡('x')"
    @test kind(index) == kind(typeof(index)) == :p
    @test kind(𝕡(1, 'x')) == kind(typeof(𝕡(1, 'x'))) == :p
    @test kind(𝕡(1, 'x', [0.0], [0.0])) == kind(typeof(𝕡(1, 'x', [0.0], [0.0]))) == :p

    @test 𝕦('x') ≠ 𝕡('x')
    @test !isequal(𝕦('x'), 𝕡('x'))

    @test kind(PhononIndex) == kind(Index{<:PhononIndex}) == kind(CoordinatedIndex{<:Index{<:PhononIndex}}) == Symbol(":")
    @test isdefinite(PhononIndex{:u, Char}) == isdefinite(PhononIndex{:p, Char}) == true
    @test isdefinite(PhononIndex{:u, Symbol}) == isdefinite(PhononIndex{:p, Symbol}) == false
    @test isdefinite(PhononIndex{:u, Colon}) == isdefinite(PhononIndex{:p, Colon}) == false
    @test internalindextype(PhononIndex{:u}, Char) == PhononIndex{:u, Char}
    @test internalindextype(PhononIndex{:p}, Symbol) == PhononIndex{:p, Symbol}
    @test internalindextype(PhononIndex{:}, Colon) == PhononIndex{:, Colon}

    @test OperatorIndex[PhononIndex{:u}] == OperatorIndex[Index{<:PhononIndex{:u}}] == OperatorIndex[CoordinatedIndex{<:Index{<:PhononIndex{:u}}}] == 𝕦
    @test OperatorIndex[PhononIndex{:p}] == OperatorIndex[Index{<:PhononIndex{:p}}] == OperatorIndex[CoordinatedIndex{<:Index{<:PhononIndex{:p}}}] == 𝕡
    @test OperatorIndex[𝕦] == PhononIndex{:u}
    @test OperatorIndex[𝕡] == PhononIndex{:p}
end

@testset "Phonon latex" begin
    index = 𝕦(1, 'x', [0.0, 0.0], [0.0, 0.0])
    @test script(index, latexofphonons, Val(:BD)) == "u"
    @test script(Index(index), latexofphonons, Val(:BD)) == "u"
    @test script(InternalIndex(index), latexofphonons, Val(:BD)) == "u"
    @test script(index, Val(:site)) == script(Index(index), Val(:site)) == "1"
    @test script(index, Val(:direction)) == script(Index(index), Val(:direction)) == script(InternalIndex(index), Val(:direction)) == "x"

    index = 𝕡(2, 'y', [0.0, 0.0], [0.0, 0.0])
    @test script(index, latexofphonons, Val(:BD)) == "p"
    @test script(Index(index), latexofphonons, Val(:BD)) == "p"
    @test script(InternalIndex(index), latexofphonons, Val(:BD)) == "p"
    @test script(index, Val(:site)) == script(Index(index), Val(:site)) == "2"
    @test script(index, Val(:direction)) == script(Index(index), Val(:direction)) == script(InternalIndex(index), Val(:direction)) == "y"

    @test latexname(PhononIndex) == Symbol("PhononIndex")
    @test latexname(Index{<:PhononIndex}) == Symbol("Index{PhononIndex}")
    @test latexname(CompositeIndex{<:Index{<:PhononIndex}}) == Symbol("CompositeIndex{Index{PhononIndex}}")
end

@testset "Phonon" begin
    pn = Phonon{:u}(3)
    @test shape(pn) == (1:3,)
    for i in axes(pn, 1)
        @test convert(PhononIndex, convert(CartesianIndex, pn[i], pn), pn) == pn[i]
    end
    @test summary(pn) == "3-element Phonon{:u}"
    @test string(pn) == "Phonon{:u}(ndirection=3)"
    @test kind(pn) == kind(typeof(pn)) == :u
    @test collect(pn) == [𝕦('x'), 𝕦('y'), 𝕦('z')]

    pn = Phonon{:p}(3)
    @test shape(pn) == (1:3,)
    for i in axes(pn, 1)
        @test convert(PhononIndex, convert(CartesianIndex, pn[i], pn), pn) == pn[i]
    end
    @test summary(pn) == "3-element Phonon{:p}"
    @test string(pn) == "Phonon{:p}(ndirection=3)"
    @test kind(pn) == kind(typeof(pn)) == :p
    @test collect(pn) == [𝕡('x'), 𝕡('y'), 𝕡('z')]

    @test Phonon(3) == Phonon{:}(3)
    @test string(Phonon{:}) == "Phonon{:}"

    @test match(PhononIndex{:u}, Phonon{:}) == match(PhononIndex{:p}, Phonon{:}) == true
    @test match(PhononIndex{:u}, Phonon{:u}) == match(PhononIndex{:p}, Phonon{:p}) == true
    @test match(PhononIndex{:u}, Phonon{:p}) == match(PhononIndex{:p}, Phonon{:u}) == false

    @test filter(PhononIndex{:u}, Phonon(3)) == Phonon{:u}(3)
    @test filter(PhononIndex{:p}, Phonon(3)) == Phonon{:p}(3)
    @test filter(PhononIndex{:u}, Phonon{:}) == Phonon{:u}
    @test filter(PhononIndex{:p}, Phonon{:}) == Phonon{:p}

    @test shape(Phonon{:u}(3), 𝕦(:)) == (1:3,)
    @test shape(Phonon{:u}(3), 𝕦('x')) == (1:1,)
    @test shape(Phonon{:p}(3), 𝕡(:)) == (1:3,)
    @test shape(Phonon{:p}(3), 𝕡('y')) == (2:2,)
end

@testset "PhononOperator" begin
    opt = Operator(1.0, 𝕡(1, 'x', [0.0, 0.0], [0.0, 0.0]), 𝕡(1, 'x', [0.0, 0.0], [0.0, 0.0]))
    @test opt' == Operator(1.0, 𝕡(1, 'x', [0.0, 0.0], [0.0, 0.0]), 𝕡(1, 'x', [0.0, 0.0], [0.0, 0.0]))
    @test latexstring(opt) == "(p^{x}_{1})^2"

    id₁ = 𝕦(1, 'x', [0.0, 0.0], [0.0, 0.0])
    id₂ = 𝕡(1, 'x', [0.0, 0.0], [0.0, 0.0])
    @test permute(id₁, id₂) == (Operator(+1im), Operator(1, id₂, id₁))
    @test permute(id₂, id₁) == (Operator(-1im), Operator(1, id₁, id₂))

    id₁ = 𝕦(1, 'x', [0.0, 0.0], [0.0, 0.0])
    id₂ = 𝕦(1, 'x', [0.0, 0.0], [0.0, 0.0])
    @test permute(id₁, id₂) == (Operator(1, id₂, id₁),)

    id₁ = 𝕦(1, 'x', [0.0, 0.0], [0.0, 0.0])
    id₂ = 𝕡(1, 'y', [0.0, 0.0], [0.0, 0.0])
    @test permute(id₁, id₂) == (Operator(1, id₂, id₁),)
end

@testset "Phonon Coupling" begin
    @test collect(MatrixCoupling(:, PhononIndex{:u}, [1 0 1; 0 1 0; 1 0 1])) == collect(MatrixCoupling(𝕦, :, [1 0 1; 0 1 0; 1 0 1])) == collect(𝕦ᵀ𝕦(:, [1 0 1; 0 1 0; 1 0 1])) == [
        Coupling(𝕦(:, 'x'), 𝕦(:, 'x')), Coupling(𝕦(:, 'z'), 𝕦(:, 'x')), Coupling(𝕦(:, 'y'), 𝕦(:, 'y')), Coupling(𝕦(:, 'x'), 𝕦(:, 'z')), Coupling(𝕦(:, 'z'), 𝕦(:, 'z'))
    ]

    pnc = Coupling(2.0, @pattern(𝕡(:, μ), 𝕡(:, μ)))
    bond = Bond(Point(1, [0.5, 0.0], [0.0, 0.0]))
    hilbert = Hilbert(Phonon(2))
    ex = expand(pnc, Val(:Kinetic), bond, hilbert)
    @test collect(ex) == [
        Operator(2.0, 𝕡(1, 'x', [0.5, 0.0], [0.0, 0.0]), 𝕡(1, 'x', [0.5, 0.0], [0.0, 0.0])),
        Operator(2.0, 𝕡(1, 'y', [0.5, 0.0], [0.0, 0.0]), 𝕡(1, 'y', [0.5, 0.0], [0.0, 0.0]))
    ]

    pnc = Coupling(𝕦(:, :), 𝕦(:, :))
    bond = Bond(1, Point(1, [0.0, 0.0], [0.0, 0.0]), Point(2, [0.5, 0.0], [0.0, 0.0]))
    hilbert = Hilbert(site=>Phonon(2) for site=1:2)
    ex = expand(pnc, Val(:Hooke), bond, hilbert)
    @test shape(ex) == (1:2, 1:2, 1:4)
    @test collect(ex) ==[
        Operator(+1.0, 𝕦(1, 'x', [0.0, 0.0], [0.0, 0.0]), 𝕦(1, 'x', [0.0, 0.0], [0.0, 0.0])),
        Operator(-0.0, 𝕦(1, 'y', [0.0, 0.0], [0.0, 0.0]), 𝕦(1, 'x', [0.0, 0.0], [0.0, 0.0])),
        Operator(-0.0, 𝕦(1, 'x', [0.0, 0.0], [0.0, 0.0]), 𝕦(1, 'y', [0.0, 0.0], [0.0, 0.0])),
        Operator(+0.0, 𝕦(1, 'y', [0.0, 0.0], [0.0, 0.0]), 𝕦(1, 'y', [0.0, 0.0], [0.0, 0.0])),
        Operator(-1.0, 𝕦(1, 'x', [0.0, 0.0], [0.0, 0.0]), 𝕦(2, 'x', [0.5, 0.0], [0.0, 0.0])),
        Operator(+0.0, 𝕦(1, 'y', [0.0, 0.0], [0.0, 0.0]), 𝕦(2, 'x', [0.5, 0.0], [0.0, 0.0])),
        Operator(+0.0, 𝕦(1, 'x', [0.0, 0.0], [0.0, 0.0]), 𝕦(2, 'y', [0.5, 0.0], [0.0, 0.0])),
        Operator(-0.0, 𝕦(1, 'y', [0.0, 0.0], [0.0, 0.0]), 𝕦(2, 'y', [0.5, 0.0], [0.0, 0.0])),
        Operator(-1.0, 𝕦(2, 'x', [0.5, 0.0], [0.0, 0.0]), 𝕦(1, 'x', [0.0, 0.0], [0.0, 0.0])),
        Operator(+0.0, 𝕦(2, 'y', [0.5, 0.0], [0.0, 0.0]), 𝕦(1, 'x', [0.0, 0.0], [0.0, 0.0])),
        Operator(+0.0, 𝕦(2, 'x', [0.5, 0.0], [0.0, 0.0]), 𝕦(1, 'y', [0.0, 0.0], [0.0, 0.0])),
        Operator(-0.0, 𝕦(2, 'y', [0.5, 0.0], [0.0, 0.0]), 𝕦(1, 'y', [0.0, 0.0], [0.0, 0.0])),
        Operator(+1.0, 𝕦(2, 'x', [0.5, 0.0], [0.0, 0.0]), 𝕦(2, 'x', [0.5, 0.0], [0.0, 0.0])),
        Operator(-0.0, 𝕦(2, 'y', [0.5, 0.0], [0.0, 0.0]), 𝕦(2, 'x', [0.5, 0.0], [0.0, 0.0])),
        Operator(-0.0, 𝕦(2, 'x', [0.5, 0.0], [0.0, 0.0]), 𝕦(2, 'y', [0.5, 0.0], [0.0, 0.0])),
        Operator(+0.0, 𝕦(2, 'y', [0.5, 0.0], [0.0, 0.0]), 𝕦(2, 'y', [0.5, 0.0], [0.0, 0.0]))
    ]
end

@testset "Kinetic" begin
    term = Kinetic(:T, 2.0)
    point = Point(1, [0.5, 0.0], [0.0, 0.0])
    hilbert = Hilbert(point.site=>Phonon(2))
    operators = Operators(
        Operator(2.0, 𝕡(1, 'x', [0.5, 0.0], [0.0, 0.0]), 𝕡(1, 'x', [0.5, 0.0], [0.0, 0.0])),
        Operator(2.0, 𝕡(1, 'y', [0.5, 0.0], [0.0, 0.0]), 𝕡(1, 'y', [0.5, 0.0], [0.0, 0.0]))
    )
    @test expand(term, Bond(point), hilbert) == operators
end

@testset "Hooke" begin
    term = Hooke(:V, 2.0, 1)

    bond = Bond(1, Point(1, [0.0, 0.0], [0.0, 0.0]), Point(2, [0.5, 0.0], [0.0, 0.0]))
    hilbert = Hilbert(site=>Phonon(2) for site=1:2)
    operators = Operators(
        Operator(+2.0, 𝕦(1, 'x', [0.0, 0.0], [0.0, 0.0]), 𝕦(1, 'x', [0.0, 0.0], [0.0, 0.0])),
        Operator(-2.0, 𝕦(2, 'x', [0.5, 0.0], [0.0, 0.0]), 𝕦(1, 'x', [0.0, 0.0], [0.0, 0.0])),
        Operator(-2.0, 𝕦(1, 'x', [0.0, 0.0], [0.0, 0.0]), 𝕦(2, 'x', [0.5, 0.0], [0.0, 0.0])),
        Operator(+2.0, 𝕦(2, 'x', [0.5, 0.0], [0.0, 0.0]), 𝕦(2, 'x', [0.5, 0.0], [0.0, 0.0]))
    )
    @test expand(term, bond, hilbert) == operators

    bond = Bond(1, Point(1, [0.0, 0.0], [0.0, 0.0]), Point(2, [0.0, 0.5], [0.0, 0.0]))
    hilbert = Hilbert(site=>Phonon(2) for site=1:2)
    operators = Operators(
        Operator(+2.0, 𝕦(2, 'y', [0.0, 0.5], [0.0, 0.0]), 𝕦(2, 'y', [0.0, 0.5], [0.0, 0.0])),
        Operator(+2.0, 𝕦(1, 'y', [0.0, 0.0], [0.0, 0.0]), 𝕦(1, 'y', [0.0, 0.0], [0.0, 0.0])),
        Operator(-2.0, 𝕦(2, 'y', [0.0, 0.5], [0.0, 0.0]), 𝕦(1, 'y', [0.0, 0.0], [0.0, 0.0])),
        Operator(-2.0, 𝕦(1, 'y', [0.0, 0.0], [0.0, 0.0]), 𝕦(2, 'y', [0.0, 0.5], [0.0, 0.0]))
    )
    @test expand(term, bond, hilbert) == operators

    bond = Bond(1, Point(1, [0.0, 0.0], [0.0, 0.0]), Point(2, [0.5, 0.5], [0.0, 0.0]))
    hilbert = Hilbert(site=>Phonon(2) for site=1:2)
    operators = Operators(
        Operator(-1.0, 𝕦(2, 'x', [0.5, 0.5], [0.0, 0.0]), 𝕦(1, 'x', [0.0, 0.0], [0.0, 0.0])),
        Operator(-1.0, 𝕦(1, 'y', [0.0, 0.0], [0.0, 0.0]), 𝕦(2, 'y', [0.5, 0.5], [0.0, 0.0])),
        Operator(-1.0, 𝕦(1, 'x', [0.0, 0.0], [0.0, 0.0]), 𝕦(2, 'x', [0.5, 0.5], [0.0, 0.0])),
        Operator(-1.0, 𝕦(1, 'y', [0.0, 0.0], [0.0, 0.0]), 𝕦(2, 'x', [0.5, 0.5], [0.0, 0.0])),
        Operator(-1.0, 𝕦(2, 'x', [0.5, 0.5], [0.0, 0.0]), 𝕦(1, 'y', [0.0, 0.0], [0.0, 0.0])),
        Operator(-1.0, 𝕦(2, 'y', [0.5, 0.5], [0.0, 0.0]), 𝕦(1, 'y', [0.0, 0.0], [0.0, 0.0])),
        Operator(+1.0, 𝕦(1, 'x', [0.0, 0.0], [0.0, 0.0]), 𝕦(1, 'y', [0.0, 0.0], [0.0, 0.0])),
        Operator(-1.0, 𝕦(1, 'x', [0.0, 0.0], [0.0, 0.0]), 𝕦(2, 'y', [0.5, 0.5], [0.0, 0.0])),
        Operator(+1.0, 𝕦(1, 'y', [0.0, 0.0], [0.0, 0.0]), 𝕦(1, 'y', [0.0, 0.0], [0.0, 0.0])),
        Operator(+1.0, 𝕦(2, 'y', [0.5, 0.5], [0.0, 0.0]), 𝕦(2, 'y', [0.5, 0.5], [0.0, 0.0])),
        Operator(+1.0, 𝕦(1, 'y', [0.0, 0.0], [0.0, 0.0]), 𝕦(1, 'x', [0.0, 0.0], [0.0, 0.0])),
        Operator(+1.0, 𝕦(2, 'y', [0.5, 0.5], [0.0, 0.0]), 𝕦(2, 'x', [0.5, 0.5], [0.0, 0.0])),
        Operator(+1.0, 𝕦(2, 'x', [0.5, 0.5], [0.0, 0.0]), 𝕦(2, 'x', [0.5, 0.5], [0.0, 0.0])),
        Operator(-1.0, 𝕦(2, 'y', [0.5, 0.5], [0.0, 0.0]), 𝕦(1, 'x', [0.0, 0.0], [0.0, 0.0])),
        Operator(+1.0, 𝕦(1, 'x', [0.0, 0.0], [0.0, 0.0]), 𝕦(1, 'x', [0.0, 0.0], [0.0, 0.0])),
        Operator(+1.0, 𝕦(2, 'x', [0.5, 0.5], [0.0, 0.0]), 𝕦(2, 'y', [0.5, 0.5], [0.0, 0.0]))
    )
    @test expand(term, bond, hilbert) ≈ operators
end

@testset "Elastic" begin
    term = Elastic(:V, 2.0, 1, 𝕦ᵀ𝕦(:, [0 1; 1 0]))
    bond = Bond(1, Point(1, [0.0, 0.0], [0.0, 0.0]), Point(2, [0.5, 0.0], [0.0, 0.0]))
    hilbert = Hilbert(site=>Phonon(2) for site=1:2)
    operators = Operators(
        Operator(1.0, 𝕦(1, 'x', [0.0, 0.0], [0.0, 0.0]), 𝕦(2, 'y', [0.5, 0.0], [0.0, 0.0])),
        Operator(1.0, 𝕦(1, 'y', [0.0, 0.0], [0.0, 0.0]), 𝕦(2, 'x', [0.5, 0.0], [0.0, 0.0])),
        Operator(1.0, 𝕦(2, 'x', [0.5, 0.0], [0.0, 0.0]), 𝕦(1, 'y', [0.0, 0.0], [0.0, 0.0])),
        Operator(1.0, 𝕦(2, 'y', [0.5, 0.0], [0.0, 0.0]), 𝕦(1, 'x', [0.0, 0.0], [0.0, 0.0])),
    )
    @test expand(term, bond, hilbert) == operators
end
