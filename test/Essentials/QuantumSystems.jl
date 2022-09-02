using LaTeXStrings: latexstring
using QuantumLattices.Essentials.DegreesOfFreedom: wildcard, AbstractCompositeIndex, CompositeIndex, Coupling, Couplings, Hilbert, Index, IIDSpace, statistics, @couplings, @subscript_str
using QuantumLattices.Essentials.QuantumOperators: ID, Operator, Operators, latexname, matrix, script
using QuantumLattices.Essentials.QuantumSystems
using QuantumLattices.Essentials.QuantumSystems: dm, gamma, heisenbergpmz, heisenbergxyz
using QuantumLattices.Essentials.Spatials: Bond, Point, azimuthd, rcoordinate
using QuantumLattices.Interfaces: ⊗, ⋅, expand, permute, rank
using QuantumLattices.Prerequisites.Combinatorics: Permutations
using QuantumLattices.Prerequisites.VectorSpaces: shape
using StaticArrays: SVector

@testset "FID" begin
    fid = FID{:f}(orbital=1, spin=1)
    @test string(fid) == "FID{:f}(1, 1, 1)"
    @test statistics(fid) == statistics(typeof(fid)) == :f
    @test fid' == replace(fid, nambu=2)
    @test fid'' == replace(fid, nambu=1)

    fid = FID{:b}(1, 1, 0)
    @test string(fid) == "FID{:b}(1, 1, 0)"
    @test statistics(fid) == statistics(typeof(fid)) == :b
    @test fid' == fid

    @test FID{:f}(1, 1, 1)≠FID{:b}(1, 1, 1)
    @test isequal(FID{:f}(1, 1, 1), FID{:f}(1, 1, 1))
    @test !isequal(FID{:f}(1, 1, 1), FID{:b}(1, 1, 1))
end

@testset "Fock" begin
    @test eltype(Fock) == (FID{S, Int, Int, Int} where S)
    fock = Fock{:b}(norbital=1, nspin=2, nnambu=2)
    @test shape(fock) == (1:1, 1:2, 1:2)
    @test CartesianIndex(FID{:b}(1, 1, 1), fock) == CartesianIndex(1, 1, 1)
    @test FID(CartesianIndex(1, 1, 1), fock) == FID{:b}(1, 1, 1)
    @test collect(fock) == [FID{:b}(1, 1, 1), FID{:b}(1, 2, 1), FID{:b}(1, 1, 2), FID{:b}(1, 2, 2)]
    @test statistics(fock) == statistics(typeof(fock)) == :b
    @test string(fock) == "Fock{:b}(norbital=1, nspin=2, nnambu=2)"

    @test summary(Fock{:b}(nspin=0, nnambu=1)) == "0-element Fock{:b}"
    @test summary(Fock{:b}(nspin=1, nnambu=1)) == "1-element Fock{:b}"
    @test summary(Fock{:f}(nspin=2, nnambu=1)) == "2-element Fock{:f}"

    @test match(FID{wildcard}, Fock{:f}) == match(FID{wildcard}, Fock{:b}) == true
    @test match(FID{:f}, Fock{:f}) == match(FID{:b}, Fock{:b}) == true
    @test match(FID{:b}, Fock{:f}) == match(FID{:f}, Fock{:b}) == false
end

@testset "latex" begin
    @test script(Val(:site), Index(1, FID{:f}(2, 1, 1))) == 1
    @test script(Val(:orbital), FID{:f}(2, 1, 1)) == script(Val(:orbital), Index(1, FID{:f}(2, 1, 1))) == 2
    @test script(Val(:spint), FID{:f}(2, 3, 1)) == script(Val(:spint), Index(1, FID{:f}(2, 3, 1))) == 3
    @test script(Val(:spsym), FID{:f}(2, 2, 1)) == script(Val(:spsym), Index(1, FID{:f}(2, 2, 1))) == "↑"
    @test script(Val(:spsym), FID{:f}(2, 1, 1)) == script(Val(:spsym), Index(1, FID{:f}(2, 1, 1))) == "↓"
    @test script(Val(:nambu), FID{:f}(2, 3, 1)) == script(Val(:nambu), Index(1, FID{:f}(2, 3, 1))) == ""
    @test script(Val(:nambu), FID{:f}(2, 3, 2)) == script(Val(:nambu), Index(1, FID{:f}(2, 3, 2))) == "\\dagger"

    @test latexname(Index{<:FID{:f}}) == Symbol("Index{FID{:f}}")
    @test latexname(AbstractCompositeIndex{Index{<:FID{:f}}}) == Symbol("AbstractCompositeIndex{Index{FID{:f}}}")
    @test latexname(FID{:f}) == Symbol("FID{:f}")
    @test latexname(Index{<:FID{:b}}) == Symbol("Index{FID{:b}}")
    @test latexname(AbstractCompositeIndex{Index{<:FID{:b}}}) == Symbol("AbstractCompositeIndex{Index{FID{:b}}}")
    @test latexname(FID{:b}) == Symbol("FID{:b}")
end

@testset "angle" begin
    @test angle(CompositeIndex(Index(1, FID{:f}(1, 1, 1)), [0.0, 0.0], [1.0, 2.0]), [[1.0, 0.0], [0.0, 1.0]], [0.1, 0.0]) ≈ 2pi*0.1
    @test angle(CompositeIndex(Index(1, FID{:f}(1, 1, 2)), [0.0, 0.0], [1.0, 2.0]), [[1.0, 0.0], [0.0, 1.0]], [0.0, 0.2]) ≈ -2pi*0.4
end

@testset "FockOperator" begin
    id₁ = CompositeIndex(Index(2, FID{:f}(1, 1, 2)), SVector(0.5, 0.0), SVector(0.0, 0.0))
    id₂ = CompositeIndex(Index(2, FID{:f}(1, 1, 1)), SVector(0.5, 0.0), SVector(0.0, 0.0))
    id₃ = CompositeIndex(Index(1, FID{:f}(1, 2, 2)), SVector(0.0, 0.0), SVector(0.0, 0.0))
    id₄ = CompositeIndex(Index(1, FID{:f}(1, 2, 1)), SVector(0.0, 0.0), SVector(0.0, 0.0))

    opt = Operator(1.0, id₁, id₂)
    @test opt|>isnormalordered

    opt = Operator(1.0, id₁, id₂, id₃, id₄)
    @test opt|>isnormalordered == false
    @test latexstring(opt) == "c^{\\dagger}_{2, 1, ↓}c^{}_{2, 1, ↓}c^{\\dagger}_{1, 1, ↑}c^{}_{1, 1, ↑}"

    op₁ = Operator(1.5, id₁, id₂)
    op₂ = Operator(2.0, id₂, id₁)
    @test op₁*op₂ == nothing

    op₁ = Operator(1.5, id₁, id₂)
    op₂ = Operator(2.0, id₁, id₂)
    @test op₁*op₂ == Operator(3.0, id₁, id₂, id₁, id₂)

    @test permute(id₁, id₂) == (Operator(1), Operator(-1, id₂, id₁))
    @test permute(id₂, id₁) == (Operator(1), Operator(-1, id₁, id₂))
    @test permute(id₁, id₄) == (Operator(-1, id₄, id₁),)
    @test permute(id₄, id₁) == (Operator(-1, id₁, id₄),)


    id₁ = CompositeIndex(Index(2, FID{:b}(1, 1, 2)), SVector(0.5, 0.0), SVector(0.0, 0.0))
    id₂ = CompositeIndex(Index(2, FID{:b}(1, 1, 1)), SVector(0.5, 0.0), SVector(0.0, 0.0))
    id₃ = CompositeIndex(Index(1, FID{:b}(1, 2, 2)), SVector(0.0, 0.0), SVector(0.0, 0.0))
    id₄ = CompositeIndex(Index(1, FID{:b}(1, 2, 1)), SVector(0.0, 0.0), SVector(0.0, 0.0))

    opt = Operator(1.0, id₁, id₂)
    @test latexstring(opt) == "b^{\\dagger}_{2, 1, ↓}b^{}_{2, 1, ↓}"

    @test permute(id₁, id₂) == (Operator(+1), Operator(1, id₂, id₁))
    @test permute(id₂, id₁) == (Operator(-1), Operator(1, id₁, id₂))
    @test permute(id₁, id₄) == (Operator(1, id₄, id₁),)
    @test permute(id₄, id₁) == (Operator(1, id₁, id₄),)
end

@testset "FockCoupling" begin
    @test string(FockCoupling{2}(1.0)) == "FockCoupling{2}(value=1.0)"
    @test string(FockCoupling{2}(1.0, spins=(1, 2))) == "FockCoupling{2}(value=1.0, spins=[1 2])"
    @test repr(FockCoupling{2}(2.0)) == "2.0 {2}"

    fc₁ = FockCoupling{2}(1.5, spins=subscript"[x 1]")
    fc₂ = FockCoupling{2}(2.0, orbitals=subscript"[x y](x < y)")
    @test repr(fc₁) == "1.5 sp[x 1]"
    @test repr(fc₂) == "2.0 ob[x y](x < y)"
    fc = fc₁ * fc₂
    @test repr(fc) == "3.0 ob[* *]×[x y](x < y) ⊗ sp[x 1]×[* *]"

    fc₁ = FockCoupling{2}(1.5, spins=subscript"[x 1]")
    fc₂ = FockCoupling{2}(2.0, orbitals=subscript"[x y](x < y)")
    fc = fc₁ ⊗ fc₂
    @test repr(fc) == "3.0 ob[x y](x < y) ⊗ sp[x 1]"

    fc₁ = FockCoupling{2}(1.5, spins=(2, 1))
    fc₂ = FockCoupling{2}(2.0, spins=(1, 2))
    fc = fc₁ ⋅ fc₂
    @test repr(fc) == "3.0 sp[2 2]"

    fc = FockCoupling{2}(2.0, orbitals=(1, 2), nambus=(2, 1))
    bond = Bond(1, Point(1, SVector(0.0), SVector(0.0)), Point(2, SVector(0.5), SVector(0.0)))
    hilbert = Hilbert(site=>Fock{:f}(norbital=2, nspin=2, nnambu=2) for site=1:2)
    ex = expand(fc, bond, hilbert, Val(:Hopping))
    @test collect(ex) == [
        Operator(2.0, CompositeIndex(Index(1, FID{:f}(1, 1, 2)), SVector(0.0), SVector(0.0)), CompositeIndex(Index(2, FID{:f}(2, 1, 1)), SVector(0.5), SVector(0.0))),
        Operator(2.0, CompositeIndex(Index(1, FID{:f}(1, 2, 2)), SVector(0.0), SVector(0.0)), CompositeIndex(Index(2, FID{:f}(2, 2, 1)), SVector(0.5), SVector(0.0)))
    ]

    fc = FockCoupling{4}(2.0, spins=(2, 2, 1, 1), nambus=(2, 1, 2, 1))
    point = Point(1, SVector(0.0), SVector(0.0))
    hilbert = Hilbert(point.site=>Fock{:b}(norbital=2, nspin=2, nnambu=2))
    ex = expand(fc, Bond(point), hilbert, Val(:info))
    @test collect(ex) == [
        Operator(2.0, 
                CompositeIndex(Index(1, FID{:b}(1, 2, 2)), SVector(0.0), SVector(0.0)),
                CompositeIndex(Index(1, FID{:b}(1, 2, 1)), SVector(0.0), SVector(0.0)),
                CompositeIndex(Index(1, FID{:b}(1, 1, 2)), SVector(0.0), SVector(0.0)),
                CompositeIndex(Index(1, FID{:b}(1, 1, 1)), SVector(0.0), SVector(0.0))
                ),
        Operator(2.0,
                CompositeIndex(Index(1, FID{:b}(2, 2, 2)), SVector(0.0), SVector(0.0)),
                CompositeIndex(Index(1, FID{:b}(2, 2, 1)), SVector(0.0), SVector(0.0)),
                CompositeIndex(Index(1, FID{:b}(2, 1, 2)), SVector(0.0), SVector(0.0)),
                CompositeIndex(Index(1, FID{:b}(2, 1, 1)), SVector(0.0), SVector(0.0))
                )
    ]

    fc = FockCoupling{4}(2.0, orbitals=subscript"[α α β β](α < β)", spins=(2, 1, 1, 2), nambus=(2, 2, 1, 1))
    point = Point(1, SVector(0.5), SVector(0.0))
    hilbert = Hilbert(point.site=>Fock{:f}(norbital=3, nspin=2, nnambu=2))
    ex = expand(fc, Bond(point), hilbert, Val(:info))
    @test collect(ex) == [
        Operator(2.0,
                CompositeIndex(Index(1, FID{:f}(1, 2, 2)), SVector(0.5), SVector(0.0)),
                CompositeIndex(Index(1, FID{:f}(1, 1, 2)), SVector(0.5), SVector(0.0)),
                CompositeIndex(Index(1, FID{:f}(2, 1, 1)), SVector(0.5), SVector(0.0)),
                CompositeIndex(Index(1, FID{:f}(2, 2, 1)), SVector(0.5), SVector(0.0))
                ),
        Operator(2.0,
                CompositeIndex(Index(1, FID{:f}(1, 2, 2)), SVector(0.5), SVector(0.0)),
                CompositeIndex(Index(1, FID{:f}(1, 1, 2)), SVector(0.5), SVector(0.0)),
                CompositeIndex(Index(1, FID{:f}(3, 1, 1)), SVector(0.5), SVector(0.0)),
                CompositeIndex(Index(1, FID{:f}(3, 2, 1)), SVector(0.5), SVector(0.0))
                ),
        Operator(2.0,
                CompositeIndex(Index(1, FID{:f}(2, 2, 2)), SVector(0.5), SVector(0.0)),
                CompositeIndex(Index(1, FID{:f}(2, 1, 2)), SVector(0.5), SVector(0.0)),
                CompositeIndex(Index(1, FID{:f}(3, 1, 1)), SVector(0.5), SVector(0.0)),
                CompositeIndex(Index(1, FID{:f}(3, 2, 1)), SVector(0.5), SVector(0.0))
                )
    ]

    fc₁ = FockCoupling{2}(+1.0, spins=(2, 2), nambus=(2, 1))
    fc₂ = FockCoupling{2}(-1.0, spins=(1, 1), nambus=(2, 1))
    point = Point(1, SVector(0.0), SVector(0.0))
    hilbert = Hilbert(point.site=>Fock{:f}(norbital=2, nspin=2, nnambu=2))
    ex = expand(fc₁*fc₂, Bond(point), hilbert, Val(:info))
    @test collect(ex) == [
        Operator(-1.0,
                CompositeIndex(Index(1, FID{:f}(1, 2, 2)), SVector(0.0), SVector(0.0)),
                CompositeIndex(Index(1, FID{:f}(1, 2, 1)), SVector(0.0), SVector(0.0)),
                CompositeIndex(Index(1, FID{:f}(1, 1, 2)), SVector(0.0), SVector(0.0)),
                CompositeIndex(Index(1, FID{:f}(1, 1, 1)), SVector(0.0), SVector(0.0))
                ),
        Operator(-1.0,
                CompositeIndex(Index(1, FID{:f}(2, 2, 2)), SVector(0.0), SVector(0.0)),
                CompositeIndex(Index(1, FID{:f}(2, 2, 1)), SVector(0.0), SVector(0.0)),
                CompositeIndex(Index(1, FID{:f}(1, 1, 2)), SVector(0.0), SVector(0.0)),
                CompositeIndex(Index(1, FID{:f}(1, 1, 1)), SVector(0.0), SVector(0.0))
                ),
        Operator(-1.0,
                CompositeIndex(Index(1, FID{:f}(1, 2, 2)), SVector(0.0), SVector(0.0)),
                CompositeIndex(Index(1, FID{:f}(1, 2, 1)), SVector(0.0), SVector(0.0)),
                CompositeIndex(Index(1, FID{:f}(2, 1, 2)), SVector(0.0), SVector(0.0)),
                CompositeIndex(Index(1, FID{:f}(2, 1, 1)), SVector(0.0), SVector(0.0))
                ),
        Operator(-1.0,
                CompositeIndex(Index(1, FID{:f}(2, 2, 2)), SVector(0.0), SVector(0.0)),
                CompositeIndex(Index(1, FID{:f}(2, 2, 1)), SVector(0.0), SVector(0.0)),
                CompositeIndex(Index(1, FID{:f}(2, 1, 2)), SVector(0.0), SVector(0.0)),
                CompositeIndex(Index(1, FID{:f}(2, 1, 1)), SVector(0.0), SVector(0.0))
                )
    ]
end

@testset "σ⁰" begin
    @test σ⁰"sp" == FockCoupling{2}(1, spins=(1, 1)) + FockCoupling{2}(1, spins=(2, 2))
    @test σ⁰"ob" == FockCoupling{2}(1, orbitals=(1, 1)) + FockCoupling{2}(1, orbitals=(2, 2))
    @test σ⁰"ph" == FockCoupling{2}(1, nambus=(1, 2)) + FockCoupling{2}(1, nambus=(2, 1))
end

@testset "σˣ" begin
    @test σˣ"sp" == FockCoupling{2}(1, spins=(1, 2)) + FockCoupling{2}(1, spins=(2, 1))
    @test σˣ"ob" == FockCoupling{2}(1, orbitals=(1, 2)) + FockCoupling{2}(1, orbitals=(2, 1))
    @test σˣ"ph" == FockCoupling{2}(1, nambus=(1, 1)) + FockCoupling{2}(1, nambus=(2, 2))
end

@testset "σʸ" begin
    @test σʸ"sp" == FockCoupling{2}(1im, spins=(1, 2)) + FockCoupling{2}(-1im, spins=(2, 1))
    @test σʸ"ob" == FockCoupling{2}(1im, orbitals=(1, 2)) + FockCoupling{2}(-1im, orbitals=(2, 1))
    @test σʸ"ph" == FockCoupling{2}(1im, nambus=(1, 1)) + FockCoupling{2}(-1im, nambus=(2, 2))
end

@testset "σᶻ" begin
    @test σᶻ"sp" == FockCoupling{2}(-1, spins=(1, 1)) + FockCoupling{2}(1, spins=(2, 2))
    @test σᶻ"ob" == FockCoupling{2}(-1, orbitals=(1, 1)) + FockCoupling{2}(1, orbitals=(2, 2))
    @test σᶻ"ph" == FockCoupling{2}(-1, nambus=(1, 2)) + FockCoupling{2}(1, nambus=(2, 1))
end

@testset "σ⁺" begin
    @test σ⁺"sp" == Couplings(FockCoupling{2}(1, spins=(2, 1)))
    @test σ⁺"ob" == Couplings(FockCoupling{2}(1, orbitals=(2, 1)))
    @test σ⁺"ph" == Couplings(FockCoupling{2}(1, nambus=(2, 2)))
end

@testset "σ⁻" begin
    @test σ⁻"sp" == Couplings(FockCoupling{2}(1, spins=(1, 2)))
    @test σ⁻"ob" == Couplings(FockCoupling{2}(1, orbitals=(1, 2)))
    @test σ⁻"ph" == Couplings(FockCoupling{2}(1, nambus=(1, 1)))
end

@testset "fockcoupling" begin
    fc = fc"1.0 ob[α α β β](α < β) ⊗ sp[σ γ σ γ](σ ≠ γ) ⊗ ph[2 1 2 1]"
    @test repr(fc) == "1.0 ob[α α β β](α < β) ⊗ sp[σ γ σ γ](σ ≠ γ) ⊗ ph[2 1 2 1]"

    fc = fc"1.0 ob[α α β β] ⊗ sp[σ γ σ γ] ⊗ ph[2 1 2 1]"
    @test repr(fc) == "1.0 ob[α α β β] ⊗ sp[σ γ σ γ] ⊗ ph[2 1 2 1]"

    fc = fc"1.0 ob[α α β β](α < β) ⊗ sp[σ γ σ γ] ⊗ ph[2 1 2 1]"
    @test repr(fc) == "1.0 ob[α α β β](α < β) ⊗ sp[σ γ σ γ] ⊗ ph[2 1 2 1]"

    fc = fc"1.0 ob[α α β β](α < β) ⊗ ph[2 1 2 1]"
    @test repr(fc) == "1.0 ob[α α β β](α < β) ⊗ ph[2 1 2 1]"

    fc = fc"1.0 ob[α α β β](α < β) ⊗ sp[2 1 2 1]"
    @test repr(fc) == "1.0 ob[α α β β](α < β) ⊗ sp[2 1 2 1]"

    fc = fc"1.0 ob[1 1 1 1] ⊗ ph[2 1 2 1]"
    @test repr(fc) == "1.0 ob[1 1 1 1] ⊗ ph[2 1 2 1]"

    fc = fc"1.0 ph[2 1 2 1]"
    @test repr(fc) == "1.0 ph[2 1 2 1]"

    fc = fc"1.0im {2}"
    @test repr(fc) == "1.0im {2}"
end

@testset "Onsite" begin
    point = Point(1, (0.5, 0.5), (0.0, 0.0))
    bond = Bond(point)
    hilbert = Hilbert(point.site=>Fock{:f}(norbital=2, nspin=2, nnambu=2))

    term = Onsite(:mu, 1.5, couplings=σˣ"sp"⊗σᶻ"ob", modulate=true)
    operators = Operators(
        Operator(+0.75, CompositeIndex(Index(1, FID{:f}(2, 2, 2)), [0.5, 0.5], [0.0, 0.0]), CompositeIndex(Index(1, FID{:f}(2, 1, 1)), [0.5, 0.5], [0.0, 0.0])),
        Operator(-0.75, CompositeIndex(Index(1, FID{:f}(1, 1, 2)), [0.5, 0.5], [0.0, 0.0]), CompositeIndex(Index(1, FID{:f}(1, 2, 1)), [0.5, 0.5], [0.0, 0.0])),
        Operator(+0.75, CompositeIndex(Index(1, FID{:f}(2, 1, 2)), [0.5, 0.5], [0.0, 0.0]), CompositeIndex(Index(1, FID{:f}(2, 2, 1)), [0.5, 0.5], [0.0, 0.0])),
        Operator(-0.75, CompositeIndex(Index(1, FID{:f}(1, 2, 2)), [0.5, 0.5], [0.0, 0.0]), CompositeIndex(Index(1, FID{:f}(1, 1, 1)), [0.5, 0.5], [0.0, 0.0]))
    )
    @test expand(term, bond, hilbert, half=true) == operators
    @test expand(term, bond, hilbert, half=false) == operators*2

    term = Onsite(:mu, 1.5, couplings=σᶻ"sp"⊗σᶻ"ob", modulate=true)
    operators = Operators(
        Operator(-0.75, CompositeIndex(Index(1, FID{:f}(2, 1, 2)), [0.5, 0.5], [0.0, 0.0]), CompositeIndex(Index(1, FID{:f}(2, 1, 1)), [0.5, 0.5], [0.0, 0.0])),
        Operator(-0.75, CompositeIndex(Index(1, FID{:f}(1, 2, 2)), [0.5, 0.5], [0.0, 0.0]), CompositeIndex(Index(1, FID{:f}(1, 2, 1)), [0.5, 0.5], [0.0, 0.0])),
        Operator(+0.75, CompositeIndex(Index(1, FID{:f}(2, 2, 2)), [0.5, 0.5], [0.0, 0.0]), CompositeIndex(Index(1, FID{:f}(2, 2, 1)), [0.5, 0.5], [0.0, 0.0])),
        Operator(+0.75, CompositeIndex(Index(1, FID{:f}(1, 1, 2)), [0.5, 0.5], [0.0, 0.0]), CompositeIndex(Index(1, FID{:f}(1, 1, 1)), [0.5, 0.5], [0.0, 0.0]))
    )
    @test expand(term, bond, hilbert, half=true) == operators
    @test expand(term, bond, hilbert, half=false) == operators*2
end

@testset "Hopping" begin
    bond = Bond(1, Point(2, (0.0, 0.0), (0.0, 0.0)), Point(1, (0.5, 0.5), (0.0, 0.0)))
    hilbert = Hilbert(site=>Fock{:f}(norbital=2, nspin=2, nnambu=2) for site=1:2)
    term = Hopping(:t, 1.5, 1)
    operators = Operators(
        Operator(1.5, CompositeIndex(Index(2, FID{:f}(2, 2, 2)), [0.0, 0.0], [0.0, 0.0]), CompositeIndex(Index(1, FID{:f}(2, 2, 1)), [0.5, 0.5], [0.0, 0.0])),
        Operator(1.5, CompositeIndex(Index(2, FID{:f}(2, 1, 2)), [0.0, 0.0], [0.0, 0.0]), CompositeIndex(Index(1, FID{:f}(2, 1, 1)), [0.5, 0.5], [0.0, 0.0])),
        Operator(1.5, CompositeIndex(Index(2, FID{:f}(1, 1, 2)), [0.0, 0.0], [0.0, 0.0]), CompositeIndex(Index(1, FID{:f}(1, 1, 1)), [0.5, 0.5], [0.0, 0.0])),
        Operator(1.5, CompositeIndex(Index(2, FID{:f}(1, 2, 2)), [0.0, 0.0], [0.0, 0.0]), CompositeIndex(Index(1, FID{:f}(1, 2, 1)), [0.5, 0.5], [0.0, 0.0]))
    )
    @test expand(term, bond, hilbert, half=true) == operators
    @test expand(term, bond, hilbert, half=false) == operators+operators'
end

@testset "Pairing" begin
    bond = Bond(1, Point(2, (0.0, 0.0), (0.0, 0.0)), Point(1, (0.5, 0.5), (0.0, 0.0)))
    hilbert = Hilbert(site=>Fock{:f}(norbital=1, nspin=2, nnambu=2) for site=1:2)
    term = Pairing(:Δ, 1.5, 1, couplings=@couplings(FockCoupling{2}(spins=(2, 2))), amplitude=bond->(bond|>rcoordinate|>azimuthd ≈ 45 ? 1 : -1))
    operators = Operators(
        Operator(-1.5, CompositeIndex(Index(2, FID{:f}(1, 2, 1)), [0.0, 0.0], [0.0, 0.0]), CompositeIndex(Index(1, FID{:f}(1, 2, 1)), [0.5, 0.5], [0.0, 0.0])),
        Operator(+1.5, CompositeIndex(Index(1, FID{:f}(1, 2, 1)), [0.5, 0.5], [0.0, 0.0]), CompositeIndex(Index(2, FID{:f}(1, 2, 1)), [0.0, 0.0], [0.0, 0.0]))
    )
    @test expand(term, bond, hilbert, half=true) == operators
    @test expand(term, bond, hilbert, half=false) == operators+operators'

    point = Point(1, (0.5, 0.5), (0.0, 0.0))
    hilbert = Hilbert(point.site=>Fock{:f}(norbital=1, nspin=2, nnambu=2))
    term = Pairing(:Δ, 1.5, 0, couplings=FockCoupling{2}(spins=(2, 1))-FockCoupling{2}(spins=(1, 2)))
    operators = Operators(
        Operator(+1.5, CompositeIndex(Index(1, FID{:f}(1, 2, 1)), [0.5, 0.5], [0.0, 0.0]), CompositeIndex(Index(1, FID{:f}(1, 1, 1)), [0.5, 0.5], [0.0, 0.0])),
        Operator(-1.5, CompositeIndex(Index(1, FID{:f}(1, 1, 1)), [0.5, 0.5], [0.0, 0.0]), CompositeIndex(Index(1, FID{:f}(1, 2, 1)), [0.5, 0.5], [0.0, 0.0]))
    )
    @test expand(term, Bond(point), hilbert, half=true) == operators
    @test expand(term, Bond(point), hilbert, half=false) == operators+operators'
end

@testset "Hubbard" begin
    point = Point(1, (0.5, 0.5), (0.0, 0.0))
    bond = Bond(point)
    hilbert = Hilbert(point.site=>Fock{:f}(norbital=2, nspin=2, nnambu=2))
    term = Hubbard(:H, 2.5)
    operators = Operators(
        Operator(1.25,
            CompositeIndex(Index(1, FID{:f}(1, 2, 2)), [0.5, 0.5], [0.0, 0.0]),
            CompositeIndex(Index(1, FID{:f}(1, 2, 1)), [0.5, 0.5], [0.0, 0.0]),
            CompositeIndex(Index(1, FID{:f}(1, 1, 2)), [0.5, 0.5], [0.0, 0.0]),
            CompositeIndex(Index(1, FID{:f}(1, 1, 1)), [0.5, 0.5], [0.0, 0.0])
            ),
        Operator(1.25,
            CompositeIndex(Index(1, FID{:f}(2, 2, 2)), [0.5, 0.5], [0.0, 0.0]),
            CompositeIndex(Index(1, FID{:f}(2, 2, 1)), [0.5, 0.5], [0.0, 0.0]),
            CompositeIndex(Index(1, FID{:f}(2, 1, 2)), [0.5, 0.5], [0.0, 0.0]),
            CompositeIndex(Index(1, FID{:f}(2, 1, 1)), [0.5, 0.5], [0.0, 0.0])
            )
    )
    @test expand(term, bond, hilbert, half=true) == operators
    @test expand(term, bond, hilbert, half=false) == operators*2
end

@testset "InterOrbitalInterSpin" begin
    point = Point(1, (0.5, 0.5), (0.0, 0.0))
    bond = Bond(point)
    hilbert = Hilbert(point.site=>Fock{:f}(norbital=2, nspin=2, nnambu=2))
    term = InterOrbitalInterSpin(:H, 2.5)
    operators = Operators(
        Operator(1.25,
            CompositeIndex(Index(1, FID{:f}(1, 2, 2)), [0.5, 0.5], [0.0, 0.0]),
            CompositeIndex(Index(1, FID{:f}(1, 2, 1)), [0.5, 0.5], [0.0, 0.0]),
            CompositeIndex(Index(1, FID{:f}(2, 1, 2)), [0.5, 0.5], [0.0, 0.0]),
            CompositeIndex(Index(1, FID{:f}(2, 1, 1)), [0.5, 0.5], [0.0, 0.0])
            ),
        Operator(1.25,
            CompositeIndex(Index(1, FID{:f}(1, 1, 2)), [0.5, 0.5], [0.0, 0.0]),
            CompositeIndex(Index(1, FID{:f}(1, 1, 1)), [0.5, 0.5], [0.0, 0.0]),
            CompositeIndex(Index(1, FID{:f}(2, 2, 2)), [0.5, 0.5], [0.0, 0.0]),
            CompositeIndex(Index(1, FID{:f}(2, 2, 1)), [0.5, 0.5], [0.0, 0.0])
            )
    )
    @test expand(term, bond, hilbert, half=true) == operators
    @test expand(term, bond, hilbert, half=false) == operators*2
end

@testset "InterOrbitalIntraSpin" begin
    point = Point(1, (0.5, 0.5), (0.0, 0.0))
    bond = Bond(point)
    hilbert = Hilbert(point.site=>Fock{:f}(norbital=2, nspin=2, nnambu=2))
    term = InterOrbitalIntraSpin(:H, 2.5)
    operators = Operators(
        Operator(1.25,
            CompositeIndex(Index(1, FID{:f}(1, 1, 2)), [0.5, 0.5], [0.0, 0.0]),
            CompositeIndex(Index(1, FID{:f}(1, 1, 1)), [0.5, 0.5], [0.0, 0.0]),
            CompositeIndex(Index(1, FID{:f}(2, 1, 2)), [0.5, 0.5], [0.0, 0.0]),
            CompositeIndex(Index(1, FID{:f}(2, 1, 1)), [0.5, 0.5], [0.0, 0.0])
            ),
        Operator(1.25,
            CompositeIndex(Index(1, FID{:f}(1, 2, 2)), [0.5, 0.5], [0.0, 0.0]),
            CompositeIndex(Index(1, FID{:f}(1, 2, 1)), [0.5, 0.5], [0.0, 0.0]),
            CompositeIndex(Index(1, FID{:f}(2, 2, 2)), [0.5, 0.5], [0.0, 0.0]),
            CompositeIndex(Index(1, FID{:f}(2, 2, 1)), [0.5, 0.5], [0.0, 0.0])
            )
    )
    @test expand(term, bond, hilbert, half=true) == operators
    @test expand(term, bond, hilbert, half=false) == operators*2
end

@testset "SpinFlip" begin
    point = Point(1, (0.5, 0.5), (0.0, 0.0))
    bond = Bond(point)
    hilbert = Hilbert(point.site=>Fock{:f}(norbital=2, nspin=2, nnambu=2))
    term = SpinFlip(:H, 2.5)
    operators = Operators(
        Operator(2.5,
            CompositeIndex(Index(1, FID{:f}(1, 2, 2)), [0.5, 0.5], [0.0, 0.0]),
            CompositeIndex(Index(1, FID{:f}(2, 1, 2)), [0.5, 0.5], [0.0, 0.0]),
            CompositeIndex(Index(1, FID{:f}(1, 1, 1)), [0.5, 0.5], [0.0, 0.0]),
            CompositeIndex(Index(1, FID{:f}(2, 2, 1)), [0.5, 0.5], [0.0, 0.0])
            )
    )
    @test expand(term, bond, hilbert, half=true) == operators
    @test expand(term, bond, hilbert, half=false) == operators+operators'
end

@testset "PairHopping" begin
    point = Point(1, (0.5, 0.5), (0.0, 0.0))
    bond = Bond(point)
    hilbert = Hilbert(point.site=>Fock{:f}(norbital=2, nspin=2, nnambu=2))
    term = PairHopping(:H, 2.5)
    operators = Operators(
        Operator(2.5,
            CompositeIndex(Index(1, FID{:f}(1, 2, 2)), [0.5, 0.5], [0.0, 0.0]),
            CompositeIndex(Index(1, FID{:f}(1, 1, 2)), [0.5, 0.5], [0.0, 0.0]),
            CompositeIndex(Index(1, FID{:f}(2, 1, 1)), [0.5, 0.5], [0.0, 0.0]),
            CompositeIndex(Index(1, FID{:f}(2, 2, 1)), [0.5, 0.5], [0.0, 0.0])
            )
    )
    @test expand(term, bond, hilbert, half=true) == operators
    @test expand(term, bond, hilbert, half=false) == operators+operators'
end

@testset "Coulomb" begin
    bond = Bond(1, Point(2, (0.0, 0.0), (0.0, 0.0)), Point(1, (0.5, 0.5), (0.0, 0.0)))
    hilbert = Hilbert(site=>Fock{:f}(norbital=1, nspin=2, nnambu=2) for site=1:2)

    term = Coulomb(:V, 2.5, 1, couplings=σᶻ"sp"*σᶻ"sp")
    operators = Operators(
        Operator(-1.25,
            CompositeIndex(Index(2, FID{:f}(1, 1, 2)), [0.0, 0.0], [0.0, 0.0]),
            CompositeIndex(Index(2, FID{:f}(1, 1, 1)), [0.0, 0.0], [0.0, 0.0]),
            CompositeIndex(Index(1, FID{:f}(1, 2, 2)), [0.5, 0.5], [0.0, 0.0]),
            CompositeIndex(Index(1, FID{:f}(1, 2, 1)), [0.5, 0.5], [0.0, 0.0])
            ),
        Operator(+1.25,
            CompositeIndex(Index(2, FID{:f}(1, 1, 2)), [0.0, 0.0], [0.0, 0.0]),
            CompositeIndex(Index(2, FID{:f}(1, 1, 1)), [0.0, 0.0], [0.0, 0.0]),
            CompositeIndex(Index(1, FID{:f}(1, 1, 2)), [0.5, 0.5], [0.0, 0.0]),
            CompositeIndex(Index(1, FID{:f}(1, 1, 1)), [0.5, 0.5], [0.0, 0.0])
            ),
        Operator(-1.25,
            CompositeIndex(Index(2, FID{:f}(1, 2, 2)), [0.0, 0.0], [0.0, 0.0]),
            CompositeIndex(Index(2, FID{:f}(1, 2, 1)), [0.0, 0.0], [0.0, 0.0]),
            CompositeIndex(Index(1, FID{:f}(1, 1, 2)), [0.5, 0.5], [0.0, 0.0]),
            CompositeIndex(Index(1, FID{:f}(1, 1, 1)), [0.5, 0.5], [0.0, 0.0])
            ),
        Operator(+1.25,
            CompositeIndex(Index(2, FID{:f}(1, 2, 2)), [0.0, 0.0], [0.0, 0.0]),
            CompositeIndex(Index(2, FID{:f}(1, 2, 1)), [0.0, 0.0], [0.0, 0.0]),
            CompositeIndex(Index(1, FID{:f}(1, 2, 2)), [0.5, 0.5], [0.0, 0.0]),
            CompositeIndex(Index(1, FID{:f}(1, 2, 1)), [0.5, 0.5], [0.0, 0.0])
            )
    )
    @test expand(term, bond, hilbert, half=true) == operators
    @test expand(term, bond, hilbert, half=false) == operators*2

    term = Coulomb(:V, 2.5, 1, couplings=σˣ"sp"*σᶻ"sp")
    operators = Operators(
        Operator(-1.25,
            CompositeIndex(Index(2, FID{:f}(1, 2, 2)), [0.0, 0.0], [0.0, 0.0]),
            CompositeIndex(Index(2, FID{:f}(1, 1, 1)), [0.0, 0.0], [0.0, 0.0]),
            CompositeIndex(Index(1, FID{:f}(1, 1, 2)), [0.5, 0.5], [0.0, 0.0]),
            CompositeIndex(Index(1, FID{:f}(1, 1, 1)), [0.5, 0.5], [0.0, 0.0])
            ),
        Operator(+1.25,
            CompositeIndex(Index(2, FID{:f}(1, 1, 2)), [0.0, 0.0], [0.0, 0.0]),
            CompositeIndex(Index(2, FID{:f}(1, 2, 1)), [0.0, 0.0], [0.0, 0.0]),
            CompositeIndex(Index(1, FID{:f}(1, 2, 2)), [0.5, 0.5], [0.0, 0.0]),
            CompositeIndex(Index(1, FID{:f}(1, 2, 1)), [0.5, 0.5], [0.0, 0.0])
            ),
        Operator(+1.25,
            CompositeIndex(Index(2, FID{:f}(1, 2, 2)), [0.0, 0.0], [0.0, 0.0]),
            CompositeIndex(Index(2, FID{:f}(1, 1, 1)), [0.0, 0.0], [0.0, 0.0]),
            CompositeIndex(Index(1, FID{:f}(1, 2, 2)), [0.5, 0.5], [0.0, 0.0]),
            CompositeIndex(Index(1, FID{:f}(1, 2, 1)), [0.5, 0.5], [0.0, 0.0])
            ),
        Operator(-1.25,
            CompositeIndex(Index(2, FID{:f}(1, 1, 2)), [0.0, 0.0], [0.0, 0.0]),
            CompositeIndex(Index(2, FID{:f}(1, 2, 1)), [0.0, 0.0], [0.0, 0.0]),
            CompositeIndex(Index(1, FID{:f}(1, 1, 2)), [0.5, 0.5], [0.0, 0.0]),
            CompositeIndex(Index(1, FID{:f}(1, 1, 1)), [0.5, 0.5], [0.0, 0.0])
            )
    )
    @test expand(term, bond, hilbert, half=true) == operators
    @test expand(term, bond, hilbert, half=false) == operators*2
end

@testset "SID" begin
    @test SID{1//2}('z', orbital=1)' == SID{1//2}('z', orbital=1)
    @test SID{1//2}('x', orbital=1)' == SID{1//2}('x', orbital=1)
    @test SID{1//2}('y', orbital=1)' == SID{1//2}('y', orbital=1)
    @test SID{1//2}('+', orbital=1)' == SID{1//2}('-', orbital=1)
    @test SID{1//2}('-', orbital=1)' == SID{1//2}('+', orbital=1)

    sid = SID{3//2}('x')
    @test string(sid) == "SID{3//2}(1, 'x')"
    @test replace(sid, tag='z') == SID{3//2}('z')
    @test totalspin(sid) == totalspin(typeof(sid)) == 3//2
    @test statistics(sid) == statistics(typeof(sid)) == :b

    @test SID{1//2}('z')≠SID{3//2}('z')
    @test isequal(SID{1//2}('z'), SID{1//2}('z'))
    @test !isequal(SID{1//2}('z'), SID{3//2}('z'))
end

@testset "matrix" begin
    @test isapprox(matrix(SID{1//2}(1, 'z')), [[-0.5, 0.0] [0.0, 0.5]])
    @test isapprox(matrix(SID{1//2}(1, 'x')), [[0.0, 0.5] [0.5, 0.0]])
    @test isapprox(matrix(SID{1//2}(1, 'y')), [[0.0, -0.5im] [0.5im, 0.0]])
    @test isapprox(matrix(SID{1//2}(1, '+')), [[0.0, 1.0] [0.0, 0.0]])
    @test isapprox(matrix(SID{1//2}(1, '-')), [[0.0, 0.0] [1.0, 0.0]])

    @test isapprox(matrix(SID{1}(1, 'z')), [[-1.0, 0.0, 0.0] [0.0, 0.0, 0.0] [0.0, 0.0, 1.0]])
    @test isapprox(matrix(SID{1}(1, 'x')), [[0.0, √2/2, 0.0] [√2/2, 0.0, √2/2] [0.0, √2/2, 0.0]])
    @test isapprox(matrix(SID{1}(1, 'y')), [[0.0, -√2im/2, 0.0] [√2im/2, 0.0, -√2im/2] [0.0, √2im/2, 0.0]])
    @test isapprox(matrix(SID{1}(1, '+')), [[0.0, √2, 0.0] [0.0, 0.0, √2] [0.0, 0.0, 0.0]])
    @test isapprox(matrix(SID{1}(1, '-')), [[0.0, 0.0, 0.0] [√2, 0.0, 0.0] [0.0, √2, 0.0]])
end

@testset "Spin" begin
    @test eltype(Spin) == (SID{S, Int, Char} where S)
    spin = Spin{1}(norbital=2)
    @test shape(spin) == (1:2, 1:5)
    @test CartesianIndex(SID{1}(1, 'z'), spin) == CartesianIndex(1, 3)
    @test SID(CartesianIndex(1, 1), spin) == SID{1}(1, 'x')
    @test summary(spin) == "10-element Spin{1}"
    @test string(spin) == "Spin{1}(norbital=2)"
    @test totalspin(spin) == totalspin(typeof(spin)) == 1
    @test collect(spin) == [
        SID{1}(1, 'x'), SID{1}(2, 'x'), SID{1}(1, 'y'), SID{1}(2, 'y'), SID{1}(1, 'z'), SID{1}(2, 'z'),
        SID{1}(1, '+'), SID{1}(2, '+'), SID{1}(1, '-'), SID{1}(2, '-')
    ]
    @test shape(IIDSpace(SID{1//2}(:a, 'x'), Spin{1//2}(2))) == (1:2, 1:1)
    @test shape(IIDSpace(SID{1//2}(2, 'y'), Spin{1//2}(2))) == (2:2, 2:2)
    @test shape(IIDSpace(SID{1//2}(:a, wildcard), Spin{1//2}(2))) == (1:2, 1:5)
    @test shape(IIDSpace(SID{1//2}(2, wildcard), Spin{1//2}(2))) == (2:2, 1:5)

    @test match(SID{wildcard}, Spin{1//2}) == true
    @test match(SID{1//2}, Spin{1//2}) == true
    @test match(SID{1//2}, Spin{1}) == match(SID{1}, Spin{1//2}) == false
end

@testset "latex" begin
    index = Index(1, SID{1//2}(2, 'z'))
    @test script(Val(:site), index) == 1
    @test script(Val(:orbital), index.iid) == script(Val(:orbital), index) == 2
    @test script(Val(:tag), index.iid) == script(Val(:tag), index) == 'z'

    @test latexname(Index{<:SID}) == Symbol("Index{SID}")
    @test latexname(AbstractCompositeIndex{<:Index{<:SID}}) == Symbol("AbstractCompositeIndex{Index{SID}}")
    @test latexname(SID) == Symbol("SID")
end

@testset "SpinOperator" begin
    opt = Operator(1.0,
        CompositeIndex(Index(1, SID{1//2}('+')), [0.0, 0.0], [0.0, 0.0]),
        CompositeIndex(Index(1, SID{1//2}('-')), [0.0, 0.0], [0.0, 0.0])
        )
    @test opt' == Operator(1.0,
        CompositeIndex(Index(1, SID{1//2}('+')), [0.0, 0.0], [0.0, 0.0]),
        CompositeIndex(Index(1, SID{1//2}('-')), [0.0, 0.0], [0.0, 0.0])
        )
    @test latexstring(opt) == "S^{+}_{1}S^{-}_{1}"
end

@testset "permute" begin
    soptrep(opt::Operator) = opt.value * prod([matrix(opt.id[i].index.iid) for i = 1:rank(opt)])
    for S in (1//2, 1, 3//2)
        indexes = [CompositeIndex(Index(1, SID{S}(2, tag)), [0.0, 0.0], [0.0, 0.0]) for tag in ('x', 'y', 'z', '+', '-')]
        for (id₁, id₂) in Permutations{2}(indexes)
            left = soptrep(Operator(1, id₁, id₂))
            right = sum([soptrep(opt) for opt in permute(id₁, id₂)])
            @test isapprox(left, right)
        end
    end
    id₁ = CompositeIndex(Index(1, SID{1//2}(2, 'z')), [0.0, 0.0], [0.0, 0.0])
    id₂ = CompositeIndex(Index(2, SID{1//2}(2, 'z')), [0.0, 0.0], [0.0, 0.0])
    @test permute(id₁, id₂) == (Operator(1, id₂, id₁),)
end

@testset "SpinCoupling" begin
    @test string(SpinCoupling(1.0, ('+', '-'))) == "SpinCoupling(value=1.0, tags=S⁺S⁻)"
    @test string(SpinCoupling(1.0, ('z', 'z'))) == "SpinCoupling(value=1.0, tags=SᶻSᶻ)"
    @test string(SpinCoupling(1.0, ('-', '+'), orbitals=(1, 2))) == "SpinCoupling(value=1.0, tags=S⁻S⁺, orbitals=[1 2])"
    @test repr(SpinCoupling(2.0, ('x', 'y'))) == "2.0 SˣSʸ"

    sc₁ = SpinCoupling(1.5, ('+', '-'), orbitals=subscript"[α β](α > β)")
    sc₂ = SpinCoupling(2.0, ('+', '-'), orbitals=subscript"[α β](α < β)")
    @test repr(sc₁) == "1.5 S⁺S⁻ ob[α β](α > β)"
    @test repr(sc₂) == "2.0 S⁺S⁻ ob[α β](α < β)"

    sc = sc₁ * sc₂
    @test repr(sc) == "3.0 S⁺S⁻S⁺S⁻ ob[α β](α > β)×[α β](α < β)"

    sc = SpinCoupling(2.0, ('+', '-'), orbitals=(1, 2))
    bond = Bond(1, Point(1, [0.0], [0.0]), Point(2, [0.5], [0.0]))
    hilbert = Hilbert(site=>Spin{1}(norbital=2) for site=1:2)
    ex = expand(sc, bond, hilbert, Val(:SpinTerm))
    @test collect(ex) == [
        Operator(2.0,
            CompositeIndex(Index(1, SID{1}(1, '+')), [0.0], [0.0]),
            CompositeIndex(Index(2, SID{1}(2, '-')), [0.5], [0.0])
        )]

    sc = SpinCoupling(2.0, ('+', '-', '+', '-'), orbitals=subscript"[α α β β](α < β)")
    point = Point(1, [0.0], [0.0])
    hilbert = Hilbert(point.site=>Spin{1}(norbital=3))
    ex = expand(sc, Bond(point), hilbert, Val(:info))
    @test collect(ex) == [
        Operator(2.0,
                CompositeIndex(Index(1, SID{1}(1, '+')), [0.0], [0.0]),
                CompositeIndex(Index(1, SID{1}(1, '-')), [0.0], [0.0]),
                CompositeIndex(Index(1, SID{1}(2, '+')), [0.0], [0.0]),
                CompositeIndex(Index(1, SID{1}(2, '-')), [0.0], [0.0])
                ),
        Operator(2.0,
                CompositeIndex(Index(1, SID{1}(1, '+')), [0.0], [0.0]),
                CompositeIndex(Index(1, SID{1}(1, '-')), [0.0], [0.0]),
                CompositeIndex(Index(1, SID{1}(3, '+')), [0.0], [0.0]),
                CompositeIndex(Index(1, SID{1}(3, '-')), [0.0], [0.0])
                ),
        Operator(2.0,
                CompositeIndex(Index(1, SID{1}(2, '+')), [0.0], [0.0]),
                CompositeIndex(Index(1, SID{1}(2, '-')), [0.0], [0.0]),
                CompositeIndex(Index(1, SID{1}(3, '+')), [0.0], [0.0]),
                CompositeIndex(Index(1, SID{1}(3, '-')), [0.0], [0.0])
                )
    ]
end

@testset "heisenberg" begin
    @test heisenberg"+-z ob[1 2]" == heisenberg"ob[1 2]" == heisenbergpmz(orbitals=(1, 2)) == Couplings(
        SpinCoupling(1//1, ('z', 'z'), orbitals=(1, 2)),
        SpinCoupling(1//2, ('+', '-'), orbitals=(1, 2)),
        SpinCoupling(1//2, ('-', '+'), orbitals=(1, 2))
    )
    @test heisenberg"" == heisenbergpmz() == Couplings(SpinCoupling(1//1, ('z', 'z')), SpinCoupling(1//2, ('+', '-')), SpinCoupling(1//2, ('-', '+')))
    @test heisenberg"xyz" == heisenbergxyz() == Couplings(SpinCoupling(1, ('x', 'x')), SpinCoupling(1, ('y', 'y')), SpinCoupling(1, ('z', 'z')))
end

@testset "ising" begin
    @test ising"x" == Couplings(SpinCoupling(1, ('x', 'x')))
    @test ising"y" == Couplings(SpinCoupling(1, ('y', 'y')))
    @test ising"z" == Couplings(SpinCoupling(1, ('z', 'z')))

    @test ising"x ob[α β](α < β)" == Couplings(SpinCoupling(1, ('x', 'x'), orbitals=subscript"[α β](α < β)"))
    @test ising"y ob[α β](α < β)" == Couplings(SpinCoupling(1, ('y', 'y'), orbitals=subscript"[α β](α < β)"))
    @test ising"z ob[α β](α < β)" == Couplings(SpinCoupling(1, ('z', 'z'), orbitals=subscript"[α β](α < β)"))
end

@testset "gamma" begin
    @test gamma"x ob[1 1]" == gamma('y', 'z', orbitals=(1, 1))== SpinCoupling(1, ('y', 'z'), orbitals=(1, 1)) + SpinCoupling(1, ('z', 'y'), orbitals=(1, 1))
    @test gamma"y" == gamma('z', 'x') == SpinCoupling(1, ('z', 'x')) + SpinCoupling(1, ('x', 'z'))
    @test gamma"z" == gamma('x', 'y') == SpinCoupling(1, ('x', 'y')) + SpinCoupling(1, ('y', 'x'))
end

@testset "dm" begin
    @test dm"x ob[1 1]" == dm('y', 'z', orbitals=(1, 1)) == SpinCoupling(1, ('y', 'z'), orbitals=(1, 1)) - SpinCoupling(1, ('z', 'y'), orbitals=(1, 1))
    @test dm"y" == dm('z', 'x') == SpinCoupling(1, ('z', 'x')) - SpinCoupling(1, ('x', 'z'))
    @test dm"z" == dm('x', 'y') == SpinCoupling(1, ('x', 'y')) - SpinCoupling(1, ('y', 'x'))
end

@testset "sᵅ" begin
    @test sˣ"ob[1]" == Couplings(SpinCoupling(1, ('x',), orbitals=(1,)))
    @test sʸ"" == Couplings(SpinCoupling(1, ('y',)))
    @test sᶻ"" == Couplings(SpinCoupling(1, ('z',)))
end

@testset "spincoupling" begin
    sc = sc"1.0 S⁺S⁻ ob[α β](α < β)"
    @test repr(sc) == "1.0 S⁺S⁻ ob[α β](α < β)"

    sc = sc"1.0 S⁺S⁻ ob[α β]"
    @test repr(sc) == "1.0 S⁺S⁻ ob[α β]"

    sc = sc"1.0 S⁺S⁻ ob[1 2]"
    @test repr(sc) == "1.0 S⁺S⁻ ob[1 2]"

    sc = sc"1.0 S⁺S⁻"
    @test repr(sc) == "1.0 S⁺S⁻"
end

@testset "SpinTerm" begin
    point = Point(1, (0.5, 0.5), (0.0, 0.0))
    hilbert = Hilbert(point.site=>Spin{1//2}(norbital=2))
    term = SpinTerm(:h, 1.5, 0, sᶻ"")
    operators = Operators(
        Operator(1.5, CompositeIndex(Index(1, SID{1//2}(1, 'z')), [0.5, 0.5], [0.0, 0.0])),
        Operator(1.5, CompositeIndex(Index(1, SID{1//2}(2, 'z')), [0.5, 0.5], [0.0, 0.0]))
    )
    @test expand(term, Bond(point), hilbert) == operators

    bond = Bond(1, Point(2, (0.5, 0.5), (0.0, 0.0)), Point(1, (0.0, 0.0), (0.0, 0.0)))
    hilbert = Hilbert(site=>Spin{1//2}(norbital=2) for site=1:2)
    term = SpinTerm(:J, 1.5, 1, heisenberg"")
    operators = Operators(
        Operator(1.50, CompositeIndex(Index(2, SID{1//2}(2, 'z')), [0.5, 0.5], [0.0, 0.0]), CompositeIndex(Index(1, SID{1//2}(2, 'z')), [0.0, 0.0], [0.0, 0.0])),
        Operator(0.75, CompositeIndex(Index(2, SID{1//2}(2, '-')), [0.5, 0.5], [0.0, 0.0]), CompositeIndex(Index(1, SID{1//2}(2, '+')), [0.0, 0.0], [0.0, 0.0])),
        Operator(0.75, CompositeIndex(Index(2, SID{1//2}(1, '-')), [0.5, 0.5], [0.0, 0.0]), CompositeIndex(Index(1, SID{1//2}(1, '+')), [0.0, 0.0], [0.0, 0.0])),
        Operator(0.75, CompositeIndex(Index(2, SID{1//2}(1, '+')), [0.5, 0.5], [0.0, 0.0]), CompositeIndex(Index(1, SID{1//2}(1, '-')), [0.0, 0.0], [0.0, 0.0])),
        Operator(1.50, CompositeIndex(Index(2, SID{1//2}(1, 'z')), [0.5, 0.5], [0.0, 0.0]), CompositeIndex(Index(1, SID{1//2}(1, 'z')), [0.0, 0.0], [0.0, 0.0])),
        Operator(0.75, CompositeIndex(Index(2, SID{1//2}(2, '+')), [0.5, 0.5], [0.0, 0.0]), CompositeIndex(Index(1, SID{1//2}(2, '-')), [0.0, 0.0], [0.0, 0.0]))
    )
    @test expand(term, bond, hilbert) == operators
end

@testset "PID" begin
    @test PID('u')' == PID('u')
    @test PID('p')' == PID('p')

    pid = PID('p')
    @test statistics(pid) == statistics(pid) == :b
end

@testset "Phonon" begin
    pn = Phonon(3)
    @test shape(pn) == (1:2, 1:3)
    for i = 1:length(pn)
        @test PID(CartesianIndex(pn[i], pn), pn) == pn[i]
    end
    @test collect(pn) == [PID('u', 'x'), PID('p', 'x'), PID('u', 'y'), PID('p', 'y'), PID('u', 'z'), PID('p', 'z')]

    @test shape(IIDSpace(PID('u'), Phonon(3))) == (1:1, 1:3)
    @test shape(IIDSpace(PID('u', 'x'), Phonon(3))) == (1:1, 1:1)
    @test shape(IIDSpace(PID('u', 'y'), Phonon(3))) == (1:1, 2:2)
    @test shape(IIDSpace(PID('u', 'z'), Phonon(3))) == (1:1, 3:3)

    @test shape(IIDSpace(PID('p'), Phonon(2))) == (2:2, 1:2)
    @test shape(IIDSpace(PID('p', 'x'), Phonon(3))) == (2:2, 1:1)
    @test shape(IIDSpace(PID('p', 'y'), Phonon(3))) == (2:2, 2:2)
    @test shape(IIDSpace(PID('p', 'z'), Phonon(3))) == (2:2, 3:3)
end

@testset "latex" begin
    index = Index(1, PID('u', 'x'))
    @test script(Val(:BD), index.iid, latexofphonons) == 'u'
    @test script(Val(:BD), index, latexofphonons) == 'u'
    @test script(Val(:BD), CompositeIndex(index, [0.0, 0.0], [0.0, 0.0]), latexofphonons) == 'u'
    @test script(Val(:site), index) == 1
    @test script(Val(:direction), index.iid) == 'x'
    @test script(Val(:direction), index) == 'x'

    index = Index(2, PID('p', 'y'))
    @test script(Val(:BD), index.iid, latexofphonons) == 'p'
    @test script(Val(:BD), index, latexofphonons) == 'p'
    @test script(Val(:BD), CompositeIndex(index, [0.0, 0.0], [0.0, 0.0]), latexofphonons) == 'p'
    @test script(Val(:site), index) == 2
    @test script(Val(:direction), index.iid) == 'y'
    @test script(Val(:direction), index) == 'y'

    @test latexname(Index{<:PID}) == Symbol("Index{PID}")
    @test latexname(AbstractCompositeIndex{<:Index{<:PID}}) == Symbol("AbstractCompositeIndex{Index{PID}}")
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

@testset "PhononCoupling" begin
    pnc = PhononCoupling(1.0, ('p', 'p'))
    @test string(pnc) == "PhononCoupling(value=1.0, tags=[p p])"
    @test repr(pnc) == "1.0 [p p]"

    pnc = PhononCoupling(1.0, ('p', 'p'), directions=('x', 'x'))
    @test string(pnc) == "PhononCoupling(value=1.0, tags=[p p], directions=[x x])"
    @test repr(pnc) == "1.0 [p p] dr[x x]"

    pnc = PhononCoupling(2.0, ('p', 'p'))
    point = Point(1, [0.5, 0.0], [0.0, 0.0])
    hilbert = Hilbert(point.site=>Phonon(2))
    ex = expand(pnc, Bond(point), hilbert, Val(:PhononKinetic))
    @test collect(ex) == [
        Operator(2.0, CompositeIndex(Index(1, PID('p', 'x')), [0.5, 0.0], [0.0, 0.0]), CompositeIndex(Index(1, PID('p', 'x')), [0.5, 0.0], [0.0, 0.0])),
        Operator(2.0, CompositeIndex(Index(1, PID('p', 'y')), [0.5, 0.0], [0.0, 0.0]), CompositeIndex(Index(1, PID('p', 'y')), [0.5, 0.0], [0.0, 0.0]))
    ]

    pnc = PhononCoupling(1.0, ('u', 'u'))
    bond = Bond(1, Point(1, [0.0, 0.0], [0.0, 0.0]), Point(2, [0.5, 0.0], [0.0, 0.0]))
    hilbert = Hilbert(site=>Phonon(2) for site=1:2)
    ex = expand(pnc, bond, hilbert, Val(:PhononPotential))
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

    @test kinetic"" == Couplings(PhononCoupling(1, ('p', 'p')))
    @test potential"" == Couplings(PhononCoupling(1, ('u', 'u')))
end

@testset "PhononKinetic" begin
    term = PhononKinetic(:T, 2.0)
    point = Point(1, [0.5, 0.0], [0.0, 0.0])
    hilbert = Hilbert(point.site=>Phonon(2))
    operators = Operators(
        Operator(2.0, CompositeIndex(Index(1, PID('p', 'x')), [0.5, 0.0], [0.0, 0.0]), CompositeIndex(Index(1, PID('p', 'x')), [0.5, 0.0], [0.0, 0.0])),
        Operator(2.0, CompositeIndex(Index(1, PID('p', 'y')), [0.5, 0.0], [0.0, 0.0]), CompositeIndex(Index(1, PID('p', 'y')), [0.5, 0.0], [0.0, 0.0]))
    )
    @test expand(term, Bond(point), hilbert) == operators
end

@testset "PhononPotential" begin
    term = PhononPotential(:V, 2.0, 1)

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
