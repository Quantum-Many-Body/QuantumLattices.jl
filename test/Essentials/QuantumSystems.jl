using Test
using StaticArrays: SVector
using QuantumLattices.Essentials.QuantumSystems
using QuantumLattices.Essentials.QuantumAlgebras: ID
using QuantumLattices.Essentials.Spatials: AbstractPID, PID, CPID, Point, Bond, rcoord, azimuthd
using QuantumLattices.Essentials.DegreesOfFreedom: Index, AbstractCompositeOID, OID, Hilbert, Operator, Operators, script, latexname, isHermitian
using QuantumLattices.Essentials.Terms: IIDSpace, Couplings, abbr, @subscript_str, @couplings, wildcard
using QuantumLattices.Interfaces: ⊗, ⋅, expand, permute, rank
using QuantumLattices.Prerequisites.Combinatorics: Permutations
using QuantumLattices.Prerequisites.VectorSpaces: shape, ndimshape

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
end

@testset "Fock" begin
    fock = Fock{:b}(norbital=1, nspin=2, nnambu=2)
    @test shape(fock) == (1:1, 1:2, 1:2)
    @test ndimshape(fock) == ndimshape(typeof(fock)) == 3
    @test CartesianIndex(FID{:b}(1, 1, 1), fock) == CartesianIndex(1, 1, 1)
    @test FID(CartesianIndex(1, 1, 1), fock) == FID{:b}(1, 1, 1)
    @test collect(fock) == [FID{:b}(1, 1, 1), FID{:b}(1, 2, 1), FID{:b}(1, 1, 2), FID{:b}(1, 2, 2)]
    @test statistics(fock) == statistics(typeof(fock)) == :b
    @test string(fock) == "Fock{:b}(norbital=1, nspin=2, nnambu=2)"

    @test summary(Fock{:b}(nspin=0, nnambu=1)) == "0-element Fock{:b}"
    @test summary(Fock{:b}(nspin=1, nnambu=1)) == "1-element Fock{:b}"
    @test summary(Fock{:f}(nspin=2, nnambu=1)) == "2-element Fock{:f}"
end

@testset "latex" begin
    @test script(Val(:site), Index(PID(1), FID{:f}(2, 1, 1))) == 1
    @test script(Val(:orbital), Index(PID(1), FID{:f}(2, 1, 1))) == 2
    @test script(Val(:spint), Index(PID(1), FID{:f}(2, 3, 1))) == 3
    @test script(Val(:spsym), Index(PID(1), FID{:f}(2, 2, 1))) == "↑"
    @test script(Val(:spsym), Index(PID(1), FID{:f}(2, 1, 1))) == "↓"
    @test script(Val(:nambu), Index(PID(1), FID{:f}(2, 3, 1))) == ""
    @test script(Val(:nambu), Index(PID(1), FID{:f}(2, 3, 2))) == "\\dagger"

    @test latexname(Index{<:AbstractPID, <:FID{:f}}) == Symbol("Index{AbstractPID, FID{:f}}")
    @test latexname(AbstractCompositeOID{Index{<:AbstractPID, <:FID{:f}}}) == Symbol("AbstractCompositeOID{Index{AbstractPID, FID{:f}}}")
    @test latexname(Index{<:AbstractPID, <:FID{:b}}) == Symbol("Index{AbstractPID, FID{:b}}")
    @test latexname(AbstractCompositeOID{Index{<:AbstractPID, <:FID{:b}}}) == Symbol("AbstractCompositeOID{Index{AbstractPID, FID{:b}}}")
end

@testset "angle" begin
    @test angle(OID(Index(PID(1), FID{:f}(1, 1, 1)), [0.0, 0.0], [1.0, 2.0]), [[1.0, 0.0], [0.0, 1.0]], [0.1, 0.0]) ≈ 2pi*0.1
    @test angle(OID(Index(CPID(1, 1), FID{:f}(1, 1, 2)), [0.0, 0.0], [1.0, 2.0]), [[1.0, 0.0], [0.0, 1.0]], [0.0, 0.2]) ≈ -2pi*0.4
end

@testset "FockOperator" begin
    id₁ = OID(Index(PID(2), FID{:f}(1, 1, 2)), SVector(0.5, 0.0), SVector(0.0, 0.0))
    id₂ = OID(Index(PID(2), FID{:f}(1, 1, 1)), SVector(0.5, 0.0), SVector(0.0, 0.0))
    id₃ = OID(Index(PID(1), FID{:f}(1, 2, 2)), SVector(0.0, 0.0), SVector(0.0, 0.0))
    id₄ = OID(Index(PID(1), FID{:f}(1, 2, 1)), SVector(0.0, 0.0), SVector(0.0, 0.0))

    opt = Operator(1.0, ID(id₁, id₂))
    @test opt|>isnormalordered

    opt = Operator(1.0, ID(id₁, id₂, id₃, id₄))
    @test opt|>isnormalordered == false
    @test repr(opt) == "c^{\\dagger}_{2, 1, ↓}c^{}_{2, 1, ↓}c^{\\dagger}_{1, 1, ↑}c^{}_{1, 1, ↑}"

    op₁ = Operator(1.5, ID(id₁, id₂))
    op₂ = Operator(2.0, ID(id₂, id₁))
    @test op₁*op₂ == nothing

    op₁ = Operator(1.5, ID(id₁, id₂))
    op₂ = Operator(2.0, ID(id₁, id₂))
    @test op₁*op₂ == Operator(3.0, ID(id₁, id₂, id₁, id₂))

    @test permute(id₁, id₂) == (Operator(1), Operator(-1, ID(id₂, id₁)))
    @test permute(id₂, id₁) == (Operator(1), Operator(-1, ID(id₁, id₂)))
    @test permute(id₁, id₄) == (Operator(-1, ID(id₄, id₁)),)
    @test permute(id₄, id₁) == (Operator(-1, ID(id₁, id₄)),)


    id₁ = OID(Index(PID(2), FID{:b}(1, 1, 2)), SVector(0.5, 0.0), SVector(0.0, 0.0))
    id₂ = OID(Index(PID(2), FID{:b}(1, 1, 1)), SVector(0.5, 0.0), SVector(0.0, 0.0))
    id₃ = OID(Index(PID(1), FID{:b}(1, 2, 2)), SVector(0.0, 0.0), SVector(0.0, 0.0))
    id₄ = OID(Index(PID(1), FID{:b}(1, 2, 1)), SVector(0.0, 0.0), SVector(0.0, 0.0))

    opt = Operator(1.0, ID(id₁, id₂))
    @test repr(opt) == "b^{\\dagger}_{2, 1, ↓}b^{}_{2, 1, ↓}"

    @test permute(id₁, id₂) == (Operator(+1), Operator(1, ID(id₂, id₁)))
    @test permute(id₂, id₁) == (Operator(-1), Operator(1, ID(id₁, id₂)))
    @test permute(id₁, id₄) == (Operator(1, ID(id₄, id₁)),)
    @test permute(id₄, id₁) == (Operator(1, ID(id₁, id₄)),)
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
    bond = Bond(1, Point(CPID(1, 2), SVector(0.5), SVector(0.0)), Point(CPID(1, 1), SVector(0.0), SVector(0.0)))
    hilbert = Hilbert(pid=>Fock{:f}(norbital=2, nspin=2, nnambu=2) for pid in [bond.epoint.pid, bond.spoint.pid])
    ex = expand(fc, bond, hilbert, Val(:Hopping))
    @test collect(ex) == [
        (2.0, ID(OID(Index(CPID(1, 1), FID{:f}(1, 1, 2)), SVector(0.0), SVector(0.0)), OID(Index(CPID(1, 2), FID{:f}(2, 1, 1)), SVector(0.5), SVector(0.0)))),
        (2.0, ID(OID(Index(CPID(1, 1), FID{:f}(1, 2, 2)), SVector(0.0), SVector(0.0)), OID(Index(CPID(1, 2), FID{:f}(2, 2, 1)), SVector(0.5), SVector(0.0))))
    ]

    fc = FockCoupling{4}(2.0, spins=(2, 2, 1, 1), nambus=(2, 1, 2, 1))
    point = Point(PID(1), SVector(0.0), SVector(0.0))
    hilbert = Hilbert(point.pid=>Fock{:b}(norbital=2, nspin=2, nnambu=2))
    ex = expand(fc, point, hilbert, Val(:info))
    @test collect(ex) == [
        (2.0, ID(OID(Index(PID(1), FID{:b}(1, 2, 2)), SVector(0.0), SVector(0.0)),
                 OID(Index(PID(1), FID{:b}(1, 2, 1)), SVector(0.0), SVector(0.0)),
                 OID(Index(PID(1), FID{:b}(1, 1, 2)), SVector(0.0), SVector(0.0)),
                 OID(Index(PID(1), FID{:b}(1, 1, 1)), SVector(0.0), SVector(0.0))
                 )),
        (2.0, ID(OID(Index(PID(1), FID{:b}(2, 2, 2)), SVector(0.0), SVector(0.0)),
                 OID(Index(PID(1), FID{:b}(2, 2, 1)), SVector(0.0), SVector(0.0)),
                 OID(Index(PID(1), FID{:b}(2, 1, 2)), SVector(0.0), SVector(0.0)),
                 OID(Index(PID(1), FID{:b}(2, 1, 1)), SVector(0.0), SVector(0.0))
                 ))
    ]

    fc = FockCoupling{4}(2.0, orbitals=subscript"[α α β β](α < β)", spins=(2, 1, 1, 2), nambus=(2, 2, 1, 1))
    point = Point(PID(1), SVector(0.5), SVector(0.0))
    hilbert = Hilbert(point.pid=>Fock{:f}(norbital=3, nspin=2, nnambu=2))
    ex = expand(fc, point, hilbert, Val(:info))
    @test collect(ex) == [
        (2.0, ID(OID(Index(PID(1), FID{:f}(1, 2, 2)), SVector(0.5), SVector(0.0)),
                 OID(Index(PID(1), FID{:f}(1, 1, 2)), SVector(0.5), SVector(0.0)),
                 OID(Index(PID(1), FID{:f}(2, 1, 1)), SVector(0.5), SVector(0.0)),
                 OID(Index(PID(1), FID{:f}(2, 2, 1)), SVector(0.5), SVector(0.0))
                 )),
        (2.0, ID(OID(Index(PID(1), FID{:f}(1, 2, 2)), SVector(0.5), SVector(0.0)),
                 OID(Index(PID(1), FID{:f}(1, 1, 2)), SVector(0.5), SVector(0.0)),
                 OID(Index(PID(1), FID{:f}(3, 1, 1)), SVector(0.5), SVector(0.0)),
                 OID(Index(PID(1), FID{:f}(3, 2, 1)), SVector(0.5), SVector(0.0))
                 )),
        (2.0, ID(OID(Index(PID(1), FID{:f}(2, 2, 2)), SVector(0.5), SVector(0.0)),
                 OID(Index(PID(1), FID{:f}(2, 1, 2)), SVector(0.5), SVector(0.0)),
                 OID(Index(PID(1), FID{:f}(3, 1, 1)), SVector(0.5), SVector(0.0)),
                 OID(Index(PID(1), FID{:f}(3, 2, 1)), SVector(0.5), SVector(0.0))
                 ))
    ]

    fc₁ = FockCoupling{2}(+1.0, spins=(2, 2), nambus=(2, 1))
    fc₂ = FockCoupling{2}(-1.0, spins=(1, 1), nambus=(2, 1))
    point = Point(CPID(1, 1), SVector(0.0), SVector(0.0))
    hilbert = Hilbert(point.pid=>Fock{:f}(norbital=2, nspin=2, nnambu=2))
    ex = expand(fc₁*fc₂, point, hilbert, Val(:info))
    @test collect(ex) == [
        (-1.0, ID(OID(Index(CPID(1, 1), FID{:f}(1, 2, 2)), SVector(0.0), SVector(0.0)),
                  OID(Index(CPID(1, 1), FID{:f}(1, 2, 1)), SVector(0.0), SVector(0.0)),
                  OID(Index(CPID(1, 1), FID{:f}(1, 1, 2)), SVector(0.0), SVector(0.0)),
                  OID(Index(CPID(1, 1), FID{:f}(1, 1, 1)), SVector(0.0), SVector(0.0))
                  )),
        (-1.0, ID(OID(Index(CPID(1, 1), FID{:f}(2, 2, 2)), SVector(0.0), SVector(0.0)),
                  OID(Index(CPID(1, 1), FID{:f}(2, 2, 1)), SVector(0.0), SVector(0.0)),
                  OID(Index(CPID(1, 1), FID{:f}(1, 1, 2)), SVector(0.0), SVector(0.0)),
                  OID(Index(CPID(1, 1), FID{:f}(1, 1, 1)), SVector(0.0), SVector(0.0))
                  )),
        (-1.0, ID(OID(Index(CPID(1, 1), FID{:f}(1, 2, 2)), SVector(0.0), SVector(0.0)),
                  OID(Index(CPID(1, 1), FID{:f}(1, 2, 1)), SVector(0.0), SVector(0.0)),
                  OID(Index(CPID(1, 1), FID{:f}(2, 1, 2)), SVector(0.0), SVector(0.0)),
                  OID(Index(CPID(1, 1), FID{:f}(2, 1, 1)), SVector(0.0), SVector(0.0))
                  )),
        (-1.0, ID(OID(Index(CPID(1, 1), FID{:f}(2, 2, 2)), SVector(0.0), SVector(0.0)),
                  OID(Index(CPID(1, 1), FID{:f}(2, 2, 1)), SVector(0.0), SVector(0.0)),
                  OID(Index(CPID(1, 1), FID{:f}(2, 1, 2)), SVector(0.0), SVector(0.0)),
                  OID(Index(CPID(1, 1), FID{:f}(2, 1, 1)), SVector(0.0), SVector(0.0))
                  ))
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
    @test abbr(Onsite) == :st
    @test isnothing(isHermitian(Onsite))

    point = Point(PID(1), (0.5, 0.5), (0.0, 0.0))
    hilbert = Hilbert(point.pid=>Fock{:f}(norbital=2, nspin=2, nnambu=2))

    term = Onsite(:mu, 1.5, couplings=σˣ"sp"⊗σᶻ"ob", modulate=true)
    operators = Operators(
        Operator(+1.5, ID(OID(Index(PID(1), FID{:f}(2, 1, 2)), [0.5, 0.5], [0.0, 0.0]), OID(Index(PID(1), FID{:f}(2, 2, 1)), [0.5, 0.5], [0.0, 0.0]))),
        Operator(-1.5, ID(OID(Index(PID(1), FID{:f}(1, 1, 2)), [0.5, 0.5], [0.0, 0.0]), OID(Index(PID(1), FID{:f}(1, 2, 1)), [0.5, 0.5], [0.0, 0.0])))
    )
    @test expand(term, point, hilbert, true) == operators
    @test expand(term, point, hilbert, false) == operators+operators'

    term = Onsite(:mu, 1.5, couplings=σᶻ"sp"⊗σᶻ"ob", modulate=true)
    operators = Operators(
        Operator(-0.75, ID(OID(Index(PID(1), FID{:f}(2, 1, 2)), [0.5, 0.5], [0.0, 0.0]), OID(Index(PID(1), FID{:f}(2, 1, 1)), [0.5, 0.5], [0.0, 0.0]))),
        Operator(-0.75, ID(OID(Index(PID(1), FID{:f}(1, 2, 2)), [0.5, 0.5], [0.0, 0.0]), OID(Index(PID(1), FID{:f}(1, 2, 1)), [0.5, 0.5], [0.0, 0.0]))),
        Operator(+0.75, ID(OID(Index(PID(1), FID{:f}(2, 2, 2)), [0.5, 0.5], [0.0, 0.0]), OID(Index(PID(1), FID{:f}(2, 2, 1)), [0.5, 0.5], [0.0, 0.0]))),
        Operator(+0.75, ID(OID(Index(PID(1), FID{:f}(1, 1, 2)), [0.5, 0.5], [0.0, 0.0]), OID(Index(PID(1), FID{:f}(1, 1, 1)), [0.5, 0.5], [0.0, 0.0])))
    )
    @test expand(term, point, hilbert, true) == operators
    @test expand(term, point, hilbert, false) == operators+operators'
end

@testset "Hopping" begin
    bond = Bond(1, Point(CPID('a', 1), (0.5, 0.5), (0.0, 0.0)), Point(CPID('b', 2), (0.0, 0.0), (0.0, 0.0)))
    hilbert = Hilbert(pid=>Fock{:f}(norbital=2, nspin=2, nnambu=2) for pid in [bond.spoint.pid, bond.epoint.pid])
    term = Hopping(:t, 1.5, 1)
    operators = Operators(
        Operator(1.5, ID(OID(Index(CPID('b', 2), FID{:f}(2, 2, 2)), [0.0, 0.0], [0.0, 0.0]), OID(Index(CPID('a', 1), FID{:f}(2, 2, 1)), [0.5, 0.5], [0.0, 0.0]))),
        Operator(1.5, ID(OID(Index(CPID('b', 2), FID{:f}(2, 1, 2)), [0.0, 0.0], [0.0, 0.0]), OID(Index(CPID('a', 1), FID{:f}(2, 1, 1)), [0.5, 0.5], [0.0, 0.0]))),
        Operator(1.5, ID(OID(Index(CPID('b', 2), FID{:f}(1, 1, 2)), [0.0, 0.0], [0.0, 0.0]), OID(Index(CPID('a', 1), FID{:f}(1, 1, 1)), [0.5, 0.5], [0.0, 0.0]))),
        Operator(1.5, ID(OID(Index(CPID('b', 2), FID{:f}(1, 2, 2)), [0.0, 0.0], [0.0, 0.0]), OID(Index(CPID('a', 1), FID{:f}(1, 2, 1)), [0.5, 0.5], [0.0, 0.0])))
    )
    @test term|>abbr == term|>typeof|>abbr == :hp
    @test term|>isHermitian == term|>typeof|>isHermitian == false
    @test expand(term, bond, hilbert, true) == operators
    @test expand(term, bond, hilbert, false) == operators+operators'
end

@testset "Pairing" begin
    bond = Bond(1, Point(PID(1), (0.5, 0.5), (0.0, 0.0)), Point(PID(2), (0.0, 0.0), (0.0, 0.0)))
    hilbert = Hilbert(pid=>Fock{:f}(norbital=1, nspin=2, nnambu=2) for pid in [bond.spoint.pid, bond.epoint.pid])
    term = Pairing(:Δ, 1.5, 1, couplings=@couplings(FockCoupling{2}(spins=(2, 2))), amplitude=bond->(bond|>rcoord|>azimuthd ≈ 45 ? 1 : -1))
    operators = Operators(
        Operator(-1.5, ID(OID(Index(PID(2), FID{:f}(1, 2, 1)), [0.0, 0.0], [0.0, 0.0]), OID(Index(PID(1), FID{:f}(1, 2, 1)), [0.5, 0.5], [0.0, 0.0]))),
        Operator(+1.5, ID(OID(Index(PID(1), FID{:f}(1, 2, 1)), [0.5, 0.5], [0.0, 0.0]), OID(Index(PID(2), FID{:f}(1, 2, 1)), [0.0, 0.0], [0.0, 0.0])))
    )
    @test term|>abbr == term|>typeof|>abbr == :pr
    @test term|>isHermitian == term|>typeof|>isHermitian == false
    @test expand(term, bond, hilbert, true) == operators
    @test expand(term, bond, hilbert, false) == operators+operators'

    point = Point(CPID('a', 1), (0.5, 0.5), (0.0, 0.0))
    hilbert = Hilbert(point.pid=>Fock{:f}(norbital=1, nspin=2, nnambu=2))
    term = Pairing(:Δ, 1.5, 0, couplings=FockCoupling{2}(spins=(2, 1))-FockCoupling{2}(spins=(1, 2)))
    operators = Operators(
        Operator(+1.5, ID(OID(Index(CPID('a', 1), FID{:f}(1, 2, 1)), [0.5, 0.5], [0.0, 0.0]), OID(Index(CPID('a', 1), FID{:f}(1, 1, 1)), [0.5, 0.5], [0.0, 0.0]))),
        Operator(-1.5, ID(OID(Index(CPID('a', 1), FID{:f}(1, 1, 1)), [0.5, 0.5], [0.0, 0.0]), OID(Index(CPID('a', 1), FID{:f}(1, 2, 1)), [0.5, 0.5], [0.0, 0.0])))
    )
    @test term|>abbr == term|>typeof|>abbr == :pr
    @test expand(term, point, hilbert, true) == operators
    @test expand(term, point, hilbert, false) == operators+operators'
end

@testset "Hubbard" begin
    point = Point(PID(1), (0.5, 0.5), (0.0, 0.0))
    hilbert = Hilbert(point.pid=>Fock{:f}(norbital=2, nspin=2, nnambu=2))
    term = Hubbard(:H, 2.5)
    operators = Operators(
        Operator(1.25, ID(
            OID(Index(PID(1), FID{:f}(1, 2, 2)), [0.5, 0.5], [0.0, 0.0]),
            OID(Index(PID(1), FID{:f}(1, 2, 1)), [0.5, 0.5], [0.0, 0.0]),
            OID(Index(PID(1), FID{:f}(1, 1, 2)), [0.5, 0.5], [0.0, 0.0]),
            OID(Index(PID(1), FID{:f}(1, 1, 1)), [0.5, 0.5], [0.0, 0.0])
            )),
        Operator(1.25, ID(
            OID(Index(PID(1), FID{:f}(2, 2, 2)), [0.5, 0.5], [0.0, 0.0]),
            OID(Index(PID(1), FID{:f}(2, 2, 1)), [0.5, 0.5], [0.0, 0.0]),
            OID(Index(PID(1), FID{:f}(2, 1, 2)), [0.5, 0.5], [0.0, 0.0]),
            OID(Index(PID(1), FID{:f}(2, 1, 1)), [0.5, 0.5], [0.0, 0.0])
            ))
    )
    @test term|>abbr == term|>typeof|>abbr == :hb
    @test term|>isHermitian == term|>typeof|>isHermitian == true
    @test expand(term, point, hilbert, true) == operators
    @test expand(term, point, hilbert, false) == operators*2
end

@testset "InterOrbitalInterSpin" begin
    point = Point(PID(1), (0.5, 0.5), (0.0, 0.0))
    hilbert = Hilbert(point.pid=>Fock{:f}(norbital=2, nspin=2, nnambu=2))
    term = InterOrbitalInterSpin(:H, 2.5)
    operators = Operators(
        Operator(1.25, ID(
            OID(Index(PID(1), FID{:f}(1, 2, 2)), [0.5, 0.5], [0.0, 0.0]),
            OID(Index(PID(1), FID{:f}(1, 2, 1)), [0.5, 0.5], [0.0, 0.0]),
            OID(Index(PID(1), FID{:f}(2, 1, 2)), [0.5, 0.5], [0.0, 0.0]),
            OID(Index(PID(1), FID{:f}(2, 1, 1)), [0.5, 0.5], [0.0, 0.0])
            )),
        Operator(1.25, ID(
            OID(Index(PID(1), FID{:f}(1, 1, 2)), [0.5, 0.5], [0.0, 0.0]),
            OID(Index(PID(1), FID{:f}(1, 1, 1)), [0.5, 0.5], [0.0, 0.0]),
            OID(Index(PID(1), FID{:f}(2, 2, 2)), [0.5, 0.5], [0.0, 0.0]),
            OID(Index(PID(1), FID{:f}(2, 2, 1)), [0.5, 0.5], [0.0, 0.0])
            ))
    )
    @test term|>abbr == term|>typeof|>abbr == :nons
    @test term|>isHermitian == term|>typeof|>isHermitian == true
    @test expand(term, point, hilbert, true) == operators
    @test expand(term, point, hilbert, false) == operators*2
end

@testset "InterOrbitalIntraSpin" begin
    point = Point(PID(1), (0.5, 0.5), (0.0, 0.0))
    hilbert = Hilbert(point.pid=>Fock{:f}(norbital=2, nspin=2, nnambu=2))
    term = InterOrbitalIntraSpin(:H, 2.5)
    operators = Operators(
        Operator(1.25, ID(
            OID(Index(PID(1), FID{:f}(1, 1, 2)), [0.5, 0.5], [0.0, 0.0]),
            OID(Index(PID(1), FID{:f}(1, 1, 1)), [0.5, 0.5], [0.0, 0.0]),
            OID(Index(PID(1), FID{:f}(2, 1, 2)), [0.5, 0.5], [0.0, 0.0]),
            OID(Index(PID(1), FID{:f}(2, 1, 1)), [0.5, 0.5], [0.0, 0.0])
            )),
        Operator(1.25, ID(
            OID(Index(PID(1), FID{:f}(1, 2, 2)), [0.5, 0.5], [0.0, 0.0]),
            OID(Index(PID(1), FID{:f}(1, 2, 1)), [0.5, 0.5], [0.0, 0.0]),
            OID(Index(PID(1), FID{:f}(2, 2, 2)), [0.5, 0.5], [0.0, 0.0]),
            OID(Index(PID(1), FID{:f}(2, 2, 1)), [0.5, 0.5], [0.0, 0.0])
            ))
    )
    @test term|>abbr == term|>typeof|>abbr == :noes
    @test term|>isHermitian == term|>typeof|>isHermitian == true
    @test expand(term, point, hilbert, true) == operators
    @test expand(term, point, hilbert, false) == operators*2
end

@testset "SpinFlip" begin
    point = Point(PID(1), (0.5, 0.5), (0.0, 0.0))
    hilbert = Hilbert(point.pid=>Fock{:f}(norbital=2, nspin=2, nnambu=2))
    term = SpinFlip(:H, 2.5)
    operators = Operators(
        Operator(2.5, ID(
            OID(Index(PID(1), FID{:f}(1, 2, 2)), [0.5, 0.5], [0.0, 0.0]),
            OID(Index(PID(1), FID{:f}(2, 1, 2)), [0.5, 0.5], [0.0, 0.0]),
            OID(Index(PID(1), FID{:f}(1, 1, 1)), [0.5, 0.5], [0.0, 0.0]),
            OID(Index(PID(1), FID{:f}(2, 2, 1)), [0.5, 0.5], [0.0, 0.0])
            ))
    )
    @test term|>abbr == term|>typeof|>abbr == :sf
    @test term|>isHermitian == term|>typeof|>isHermitian == false
    @test expand(term, point, hilbert, true) == operators
    @test expand(term, point, hilbert, false) == operators+operators'
end

@testset "PairHopping" begin
    point = Point(PID(1), (0.5, 0.5), (0.0, 0.0))
    hilbert = Hilbert(point.pid=>Fock{:f}(norbital=2, nspin=2, nnambu=2))
    term = PairHopping(:H, 2.5)
    operators = Operators(
        Operator(2.5, ID(
            OID(Index(PID(1), FID{:f}(1, 2, 2)), [0.5, 0.5], [0.0, 0.0]),
            OID(Index(PID(1), FID{:f}(1, 1, 2)), [0.5, 0.5], [0.0, 0.0]),
            OID(Index(PID(1), FID{:f}(2, 1, 1)), [0.5, 0.5], [0.0, 0.0]),
            OID(Index(PID(1), FID{:f}(2, 2, 1)), [0.5, 0.5], [0.0, 0.0])
            ))
    )
    @test term|>abbr == term|>typeof|>abbr == :ph
    @test term|>isHermitian == term|>typeof|>isHermitian == false
    @test expand(term, point, hilbert, true) == operators
    @test expand(term, point, hilbert, false) == operators+operators'
end

@testset "Coulomb" begin
    bond = Bond(1, Point(PID(1), (0.5, 0.5), (0.0, 0.0)), Point(PID(2), (0.0, 0.0), (0.0, 0.0)))
    hilbert = Hilbert(pid=>Fock{:f}(norbital=1, nspin=2, nnambu=2) for pid in [bond.spoint.pid, bond.epoint.pid])

    term = Coulomb(:V, 2.5, 1, couplings=σᶻ"sp"*σᶻ"sp")
    operators = Operators(
        Operator(-1.25, ID(
            OID(Index(PID(2), FID{:f}(1, 1, 2)), [0.0, 0.0], [0.0, 0.0]),
            OID(Index(PID(2), FID{:f}(1, 1, 1)), [0.0, 0.0], [0.0, 0.0]),
            OID(Index(PID(1), FID{:f}(1, 2, 2)), [0.5, 0.5], [0.0, 0.0]),
            OID(Index(PID(1), FID{:f}(1, 2, 1)), [0.5, 0.5], [0.0, 0.0])
            )),
        Operator(+1.25, ID(
            OID(Index(PID(2), FID{:f}(1, 1, 2)), [0.0, 0.0], [0.0, 0.0]),
            OID(Index(PID(2), FID{:f}(1, 1, 1)), [0.0, 0.0], [0.0, 0.0]),
            OID(Index(PID(1), FID{:f}(1, 1, 2)), [0.5, 0.5], [0.0, 0.0]),
            OID(Index(PID(1), FID{:f}(1, 1, 1)), [0.5, 0.5], [0.0, 0.0])
            )),
        Operator(-1.25, ID(
            OID(Index(PID(2), FID{:f}(1, 2, 2)), [0.0, 0.0], [0.0, 0.0]),
            OID(Index(PID(2), FID{:f}(1, 2, 1)), [0.0, 0.0], [0.0, 0.0]),
            OID(Index(PID(1), FID{:f}(1, 1, 2)), [0.5, 0.5], [0.0, 0.0]),
            OID(Index(PID(1), FID{:f}(1, 1, 1)), [0.5, 0.5], [0.0, 0.0])
            )),
        Operator(+1.25, ID(
            OID(Index(PID(2), FID{:f}(1, 2, 2)), [0.0, 0.0], [0.0, 0.0]),
            OID(Index(PID(2), FID{:f}(1, 2, 1)), [0.0, 0.0], [0.0, 0.0]),
            OID(Index(PID(1), FID{:f}(1, 2, 2)), [0.5, 0.5], [0.0, 0.0]),
            OID(Index(PID(1), FID{:f}(1, 2, 1)), [0.5, 0.5], [0.0, 0.0])
            ))
    )
    @test term|>abbr == term|>typeof|>abbr == :cl
    @test isnothing(term|>isHermitian) && isnothing(term|>typeof|>isHermitian)
    @test expand(term, bond, hilbert, true) == operators
    @test expand(term, bond, hilbert, false) == operators+operators'

    term = Coulomb(:V, 2.5, 1, couplings=σˣ"sp"*σᶻ"sp")
    operators = Operators(
        Operator(-2.5, ID(
            OID(Index(PID(2), FID{:f}(1, 2, 2)), [0.0, 0.0], [0.0, 0.0]),
            OID(Index(PID(2), FID{:f}(1, 1, 1)), [0.0, 0.0], [0.0, 0.0]),
            OID(Index(PID(1), FID{:f}(1, 1, 2)), [0.5, 0.5], [0.0, 0.0]),
            OID(Index(PID(1), FID{:f}(1, 1, 1)), [0.5, 0.5], [0.0, 0.0])
            )),
        Operator(+2.5, ID(
            OID(Index(PID(2), FID{:f}(1, 1, 2)), [0.0, 0.0], [0.0, 0.0]),
            OID(Index(PID(2), FID{:f}(1, 2, 1)), [0.0, 0.0], [0.0, 0.0]),
            OID(Index(PID(1), FID{:f}(1, 2, 2)), [0.5, 0.5], [0.0, 0.0]),
            OID(Index(PID(1), FID{:f}(1, 2, 1)), [0.5, 0.5], [0.0, 0.0])
            )),
        Operator(+2.5, ID(
            OID(Index(PID(2), FID{:f}(1, 2, 2)), [0.0, 0.0], [0.0, 0.0]),
            OID(Index(PID(2), FID{:f}(1, 1, 1)), [0.0, 0.0], [0.0, 0.0]),
            OID(Index(PID(1), FID{:f}(1, 2, 2)), [0.5, 0.5], [0.0, 0.0]),
            OID(Index(PID(1), FID{:f}(1, 2, 1)), [0.5, 0.5], [0.0, 0.0])
            )),
        Operator(-2.5, ID(
            OID(Index(PID(2), FID{:f}(1, 1, 2)), [0.0, 0.0], [0.0, 0.0]),
            OID(Index(PID(2), FID{:f}(1, 2, 1)), [0.0, 0.0], [0.0, 0.0]),
            OID(Index(PID(1), FID{:f}(1, 1, 2)), [0.5, 0.5], [0.0, 0.0]),
            OID(Index(PID(1), FID{:f}(1, 1, 1)), [0.5, 0.5], [0.0, 0.0])
            ))
    )
    @test term|>abbr == term|>typeof|>abbr == :cl
    @test expand(term, bond, hilbert, true) == operators
    @test expand(term, bond, hilbert, false) == operators+operators'
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
end

@testset "Matrix" begin
    @test isapprox(Matrix(SID{1//2}(1, 'z')), [[-0.5, 0.0] [0.0, 0.5]])
    @test isapprox(Matrix(SID{1//2}(1, 'x')), [[0.0, 0.5] [0.5, 0.0]])
    @test isapprox(Matrix(SID{1//2}(1, 'y')), [[0.0, -0.5im] [0.5im, 0.0]])
    @test isapprox(Matrix(SID{1//2}(1, '+')), [[0.0, 1.0] [0.0, 0.0]])
    @test isapprox(Matrix(SID{1//2}(1, '-')), [[0.0, 0.0] [1.0, 0.0]])

    @test isapprox(Matrix(SID{1}(1, 'z')), [[-1.0, 0.0, 0.0] [0.0, 0.0, 0.0] [0.0, 0.0, 1.0]])
    @test isapprox(Matrix(SID{1}(1, 'x')), [[0.0, √2/2, 0.0] [√2/2, 0.0, √2/2] [0.0, √2/2, 0.0]])
    @test isapprox(Matrix(SID{1}(1, 'y')), [[0.0, -√2im/2, 0.0] [√2im/2, 0.0, -√2im/2] [0.0, √2im/2, 0.0]])
    @test isapprox(Matrix(SID{1}(1, '+')), [[0.0, √2, 0.0] [0.0, 0.0, √2] [0.0, 0.0, 0.0]])
    @test isapprox(Matrix(SID{1}(1, '-')), [[0.0, 0.0, 0.0] [√2, 0.0, 0.0] [0.0, √2, 0.0]])
end

@testset "Spin" begin
    spin = Spin{1}(norbital=2)
    @test shape(spin) == (1:2, 1:5)
    @test ndimshape(spin) == ndimshape(typeof(spin)) == 2
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
end

@testset "latex" begin
    index = Index(PID(1), SID{1//2}(2, 'z'))
    @test script(Val(:site), index) == 1
    @test script(Val(:orbital), index) == 2
    @test script(Val(:tag), index) == 'z'

    @test latexname(Index{<:AbstractPID, <:SID}) == Symbol("Index{AbstractPID, SID}")
    @test latexname(AbstractCompositeOID{<:Index{<:AbstractPID, <:SID}}) == Symbol("AbstractCompositeOID{Index{AbstractPID, SID}}")
end

@testset "SpinOperator" begin
    opt = Operator(1.0, ID(
        OID(Index(PID(1), SID{1//2}('+')), [0.0, 0.0], [0.0, 0.0]),
        OID(Index(PID(1), SID{1//2}('-')), [0.0, 0.0], [0.0, 0.0])
        ))
    @test opt' == Operator(1.0, ID(
        OID(Index(PID(1), SID{1//2}('+')), [0.0, 0.0], [0.0, 0.0]),
        OID(Index(PID(1), SID{1//2}('-')), [0.0, 0.0], [0.0, 0.0])
        ))
    @test isHermitian(opt)
    @test repr(opt) == "S^{+}_{1}S^{-}_{1}"
end

@testset "permute" begin
    soptrep(opt::Operator) = opt.value * prod([Matrix(opt.id[i].index.iid) for i = 1:rank(opt)])
    for S in (1//2, 1, 3//2)
        oids = [OID(Index(PID(1), SID{S}(2, tag)), [0.0, 0.0], [0.0, 0.0]) for tag in ('x', 'y', 'z', '+', '-')]
        for (id₁, id₂) in Permutations{2}(oids)
            left = soptrep(Operator(1, ID(id₁, id₂)))
            right = sum([soptrep(opt) for opt in permute(id₁, id₂)])
            @test isapprox(left, right)
        end
    end
    id₁ = OID(Index(PID(1), SID{1//2}(2, 'z')), [0.0, 0.0], [0.0, 0.0])
    id₂ = OID(Index(PID(2), SID{1//2}(2, 'z')), [0.0, 0.0], [0.0, 0.0])
    @test permute(id₁, id₂) == (Operator(1, ID(id₂, id₁)),)
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
    bond = Bond(1, Point(CPID(1, 2), [0.5], [0.0]), Point(CPID(1, 1), [0.0], [0.0]))
    hilbert = Hilbert(pid=>Spin{1}(norbital=2) for pid in [bond.epoint.pid, bond.spoint.pid])
    ex = expand(sc, bond, hilbert, Val(:SpinTerm))
    @test collect(ex) == [(2.0, ID(
        OID(Index(CPID(1, 1), SID{1}(1, '+')), [0.0], [0.0]),
        OID(Index(CPID(1, 2), SID{1}(2, '-')), [0.5], [0.0])
        ))]

    sc = SpinCoupling(2.0, ('+', '-', '+', '-'), orbitals=subscript"[α α β β](α < β)")
    point = Point(PID(1), [0.0], [0.0])
    hilbert = Hilbert(point.pid=>Spin{1}(norbital=3))
    ex = expand(sc, point, hilbert, Val(:info))
    @test collect(ex) == [
        (2.0, ID(OID(Index(PID(1), SID{1}(1, '+')), [0.0], [0.0]),
                 OID(Index(PID(1), SID{1}(1, '-')), [0.0], [0.0]),
                 OID(Index(PID(1), SID{1}(2, '+')), [0.0], [0.0]),
                 OID(Index(PID(1), SID{1}(2, '-')), [0.0], [0.0])
                 )),
        (2.0, ID(OID(Index(PID(1), SID{1}(1, '+')), [0.0], [0.0]),
                 OID(Index(PID(1), SID{1}(1, '-')), [0.0], [0.0]),
                 OID(Index(PID(1), SID{1}(3, '+')), [0.0], [0.0]),
                 OID(Index(PID(1), SID{1}(3, '-')), [0.0], [0.0])
                 )),
        (2.0, ID(OID(Index(PID(1), SID{1}(2, '+')), [0.0], [0.0]),
                 OID(Index(PID(1), SID{1}(2, '-')), [0.0], [0.0]),
                 OID(Index(PID(1), SID{1}(3, '+')), [0.0], [0.0]),
                 OID(Index(PID(1), SID{1}(3, '-')), [0.0], [0.0])
                 ))
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
    point = Point(PID(1), (0.5, 0.5), (0.0, 0.0))
    hilbert = Hilbert(point.pid=>Spin{1//2}(norbital=2))
    term = SpinTerm{1}(:h, 1.5, 0, couplings=sᶻ"")
    operators = Operators(
        Operator(1.5, ID(OID(Index(PID(1), SID{1//2}(1, 'z')), [0.5, 0.5], [0.0, 0.0]))),
        Operator(1.5, ID(OID(Index(PID(1), SID{1//2}(2, 'z')), [0.5, 0.5], [0.0, 0.0])))
    )
    @test term|>abbr == term|>typeof|>abbr == :sp
    @test term|>isHermitian == term|>typeof|>isHermitian == true
    @test expand(term, point, hilbert) == operators

    bond = Bond(1, Point(CPID('a', 1), (0.0, 0.0), (0.0, 0.0)), Point(CPID('b', 1), (0.5, 0.5), (0.0, 0.0)))
    hilbert = Hilbert(pid=>Spin{1//2}(norbital=2) for pid in [bond.spoint.pid, bond.epoint.pid])
    term = SpinTerm{2}(:J, 1.5, 1, couplings=heisenberg"")
    operators = Operators(
        Operator(1.50, ID(OID(Index(CPID('b', 1), SID{1//2}(2, 'z')), [0.5, 0.5], [0.0, 0.0]), OID(Index(CPID('a', 1), SID{1//2}(2, 'z')), [0.0, 0.0], [0.0, 0.0]))),
        Operator(0.75, ID(OID(Index(CPID('b', 1), SID{1//2}(2, '-')), [0.5, 0.5], [0.0, 0.0]), OID(Index(CPID('a', 1), SID{1//2}(2, '+')), [0.0, 0.0], [0.0, 0.0]))),
        Operator(0.75, ID(OID(Index(CPID('b', 1), SID{1//2}(1, '-')), [0.5, 0.5], [0.0, 0.0]), OID(Index(CPID('a', 1), SID{1//2}(1, '+')), [0.0, 0.0], [0.0, 0.0]))),
        Operator(0.75, ID(OID(Index(CPID('b', 1), SID{1//2}(1, '+')), [0.5, 0.5], [0.0, 0.0]), OID(Index(CPID('a', 1), SID{1//2}(1, '-')), [0.0, 0.0], [0.0, 0.0]))),
        Operator(1.50, ID(OID(Index(CPID('b', 1), SID{1//2}(1, 'z')), [0.5, 0.5], [0.0, 0.0]), OID(Index(CPID('a', 1), SID{1//2}(1, 'z')), [0.0, 0.0], [0.0, 0.0]))),
        Operator(0.75, ID(OID(Index(CPID('b', 1), SID{1//2}(2, '+')), [0.5, 0.5], [0.0, 0.0]), OID(Index(CPID('a', 1), SID{1//2}(2, '-')), [0.0, 0.0], [0.0, 0.0])))
    )
    @test expand(term, bond, hilbert) == operators
end

@testset "PNID" begin
    @test PNID('u')' == PNID('u')
    @test PNID('p')' == PNID('p')
end

@testset "Phonon" begin
    pn = Phonon(3)
    @test shape(pn) == (1:2, 1:3)
    @test ndimshape(pn) == ndimshape(typeof(pn)) == 2
    for i = 1:length(pn)
        @test PNID(CartesianIndex(pn[i], pn), pn) == pn[i]
    end
    @test collect(pn) == [PNID('u', 'x'), PNID('p', 'x'), PNID('u', 'y'), PNID('p', 'y'), PNID('u', 'z'), PNID('p', 'z')]

    @test shape(IIDSpace(PNID('u'), Phonon(3))) == (1:1, 1:3)
    @test shape(IIDSpace(PNID('u', 'x'), Phonon(3))) == (1:1, 1:1)
    @test shape(IIDSpace(PNID('u', 'y'), Phonon(3))) == (1:1, 2:2)
    @test shape(IIDSpace(PNID('u', 'z'), Phonon(3))) == (1:1, 3:3)

    @test shape(IIDSpace(PNID('p'), Phonon(2))) == (2:2, 1:2)
    @test shape(IIDSpace(PNID('p', 'x'), Phonon(3))) == (2:2, 1:1)
    @test shape(IIDSpace(PNID('p', 'y'), Phonon(3))) == (2:2, 2:2)
    @test shape(IIDSpace(PNID('p', 'z'), Phonon(3))) == (2:2, 3:3)
end

@testset "latex" begin
    index = Index(PID(1), PNID('u', 'x'))
    @test script(Val(:BD), index, pndefaultlatex) == 'u'
    @test script(Val(:BD), OID(index, [0.0, 0.0], [0.0, 0.0]), pndefaultlatex) == 'u'
    @test script(Val(:site), index) == 1
    @test script(Val(:dir), index) == 'x'

    index = Index(PID(2), PNID('p', 'y'))
    @test script(Val(:BD), index, pndefaultlatex) == 'p'
    @test script(Val(:BD), OID(index, [0.0, 0.0], [0.0, 0.0]), pndefaultlatex) == 'p'
    @test script(Val(:site), index) == 2
    @test script(Val(:dir), index) == 'y'

    latexname(Index{<:AbstractPID, <:PNID}) == Symbol("Index{AbstractPID, PNID}")
    latexname(AbstractCompositeOID{<:Index{<:AbstractPID, <:PNID}}) == Symbol("AbstractCompositeOID{Index{AbstractPID, PNID}}")
end

@testset "PhononOperator" begin
    opt = Operator(1.0, ID(
        OID(Index(PID(1), PNID('p', 'x')), [0.0, 0.0], [0.0, 0.0]),
        OID(Index(PID(1), PNID('p', 'x')), [0.0, 0.0], [0.0, 0.0])
        ))
    @test opt' == Operator(1.0, ID(
        OID(Index(PID(1), PNID('p', 'x')), [0.0, 0.0], [0.0, 0.0]),
        OID(Index(PID(1), PNID('p', 'x')), [0.0, 0.0], [0.0, 0.0])
        ))
    @test isHermitian(opt)
    @test repr(opt) == "(p^{}_{1x})^2"
end

@testset "permute" begin
    id₁ = OID(Index(PID(1), PNID('u', 'x')), [0.0, 0.0], [0.0, 0.0])
    id₂ = OID(Index(PID(1), PNID('p', 'x')), [0.0, 0.0], [0.0, 0.0])
    @test permute(id₁, id₂) == (Operator(+1im), Operator(1, ID(id₂, id₁)))
    @test permute(id₂, id₁) == (Operator(-1im), Operator(1, ID(id₁, id₂)))

    id₁ = OID(Index(PID(1), PNID('u', 'x')), [0.0, 0.0], [0.0, 0.0])
    id₂ = OID(Index(PID(1), PNID('u', 'x')), [0.0, 0.0], [0.0, 0.0])
    @test permute(id₁, id₂) == (Operator(1, ID(id₂, id₁)),)

    id₁ = OID(Index(PID(1), PNID('u', 'x')), [0.0, 0.0], [0.0, 0.0])
    id₂ = OID(Index(PID(1), PNID('p', 'y')), [0.0, 0.0], [0.0, 0.0])
    @test permute(id₁, id₂) == (Operator(1, ID(id₂, id₁)),)
end

@testset "PhononCoupling" begin
    pnc = PhononCoupling(1.0, ('p', 'p'))
    @test string(pnc) == "PhononCoupling(value=1.0, tags=[p p])"
    @test repr(pnc) == "1.0 [p p]"

    pnc = PhononCoupling(1.0, ('p', 'p'), dirs=('x', 'x'))
    @test string(pnc) == "PhononCoupling(value=1.0, tags=[p p], dirs=[x x])"
    @test repr(pnc) == "1.0 [p p] dr[x x]"

    pnc = PhononCoupling(2.0, ('p', 'p'))
    point = Point(PID(1), [0.5, 0.0], [0.0, 0.0])
    hilbert = Hilbert(point.pid=>Phonon(2))
    ex = expand(pnc, point, hilbert, Val(:PhononKinetic))
    @test collect(ex) == [
        (2.0, ID(OID(Index(PID(1), PNID('p', 'x')), [0.5, 0.0], [0.0, 0.0]), OID(Index(PID(1), PNID('p', 'x')), [0.5, 0.0], [0.0, 0.0]))),
        (2.0, ID(OID(Index(PID(1), PNID('p', 'y')), [0.5, 0.0], [0.0, 0.0]), OID(Index(PID(1), PNID('p', 'y')), [0.5, 0.0], [0.0, 0.0])))
    ]

    pnc = PhononCoupling(1.0, ('u', 'u'))
    bond = Bond(1, Point(PID(2), [0.5, 0.0], [0.0, 0.0]), Point(PID(1), [0.0, 0.0], [0.0, 0.0]))
    hilbert = Hilbert(pid=>Phonon(2) for pid in [bond.epoint.pid, bond.spoint.pid])
    ex = expand(pnc, bond, hilbert, Val(:PhononPotential))
    @test shape(ex) == (1:2, 1:3)
    @test collect(ex) ==[
        (+1.0, ID(OID(Index(PID(1), PNID('u', 'x')), [0.0, 0.0], [0.0, 0.0]), OID(Index(PID(1), PNID('u', 'x')), [0.0, 0.0], [0.0, 0.0]))),
        (+0.0, ID(OID(Index(PID(1), PNID('u', 'y')), [0.0, 0.0], [0.0, 0.0]), OID(Index(PID(1), PNID('u', 'y')), [0.0, 0.0], [0.0, 0.0]))),
        (-2.0, ID(OID(Index(PID(1), PNID('u', 'x')), [0.0, 0.0], [0.0, 0.0]), OID(Index(PID(2), PNID('u', 'x')), [0.5, 0.0], [0.0, 0.0]))), 
        (+0.0, ID(OID(Index(PID(1), PNID('u', 'y')), [0.0, 0.0], [0.0, 0.0]), OID(Index(PID(2), PNID('u', 'y')), [0.5, 0.0], [0.0, 0.0]))),
        (+1.0, ID(OID(Index(PID(2), PNID('u', 'x')), [0.5, 0.0], [0.0, 0.0]), OID(Index(PID(2), PNID('u', 'x')), [0.5, 0.0], [0.0, 0.0]))),
        (+0.0, ID(OID(Index(PID(2), PNID('u', 'y')), [0.5, 0.0], [0.0, 0.0]), OID(Index(PID(2), PNID('u', 'y')), [0.5, 0.0], [0.0, 0.0])))
        ]

    @test kinetic"" == Couplings(PhononCoupling(1, ('p', 'p')))
    @test potential"" == Couplings(PhononCoupling(1, ('u', 'u')))
end

@testset "PhononKinetic" begin
    term = PhononKinetic(:T, 2.0)
    @test abbr(term) == abbr(typeof(term)) == :pnk
    @test isHermitian(term) == isHermitian(typeof(term)) == true

    point = Point(PID(1), [0.5, 0.0], [0.0, 0.0])
    hilbert = Hilbert(point.pid=>Phonon(2))
    operators = Operators(
        Operator(2.0, ID(OID(Index(PID(1), PNID('p', 'x')), [0.5, 0.0], [0.0, 0.0]), OID(Index(PID(1), PNID('p', 'x')), [0.5, 0.0], [0.0, 0.0]))),
        Operator(2.0, ID(OID(Index(PID(1), PNID('p', 'y')), [0.5, 0.0], [0.0, 0.0]), OID(Index(PID(1), PNID('p', 'y')), [0.5, 0.0], [0.0, 0.0])))
    )
    @test expand(term, point, hilbert) == operators
end

@testset "PhononPotential" begin
    term = PhononPotential(:V, 2.0, 1)
    @test abbr(term) == abbr(typeof(term)) == :pnp
    @test isHermitian(term) == isHermitian(typeof(term)) == true

    bond = Bond(1, Point(PID(2), [0.5, 0.0], [0.0, 0.0]), Point(PID(1), [0.0, 0.0], [0.0, 0.0]))
    hilbert = Hilbert(pid=>Phonon(2) for pid in [bond.epoint.pid, bond.spoint.pid])
    operators = Operators(
        Operator(+2.0, ID(OID(Index(PID(1), PNID('u', 'x')), [0.0, 0.0], [0.0, 0.0]), OID(Index(PID(1), PNID('u', 'x')), [0.0, 0.0], [0.0, 0.0]))),
        Operator(-4.0, ID(OID(Index(PID(1), PNID('u', 'x')), [0.0, 0.0], [0.0, 0.0]), OID(Index(PID(2), PNID('u', 'x')), [0.5, 0.0], [0.0, 0.0]))), 
        Operator(+2.0, ID(OID(Index(PID(2), PNID('u', 'x')), [0.5, 0.0], [0.0, 0.0]), OID(Index(PID(2), PNID('u', 'x')), [0.5, 0.0], [0.0, 0.0]))),
    )
    @test expand(term, bond, hilbert) == operators

    bond = Bond(1, Point(PID(2), [0.0, 0.5], [0.0, 0.0]), Point(PID(1), [0.0, 0.0], [0.0, 0.0]))
    hilbert = Hilbert(pid=>Phonon(2) for pid in [bond.epoint.pid, bond.spoint.pid])
    operators = Operators(
        Operator(+2.0, ID(OID(Index(PID(1), PNID('u', 'y')), [0.0, 0.0], [0.0, 0.0]), OID(Index(PID(1), PNID('u', 'y')), [0.0, 0.0], [0.0, 0.0]))),
        Operator(-4.0, ID(OID(Index(PID(1), PNID('u', 'y')), [0.0, 0.0], [0.0, 0.0]), OID(Index(PID(2), PNID('u', 'y')), [0.0, 0.5], [0.0, 0.0]))),
        Operator(+2.0, ID(OID(Index(PID(2), PNID('u', 'y')), [0.0, 0.5], [0.0, 0.0]), OID(Index(PID(2), PNID('u', 'y')), [0.0, 0.5], [0.0, 0.0])))
    )
    @test expand(term, bond, hilbert) == operators

    bond = Bond(1, Point(PID(2), [0.5, 0.5], [0.0, 0.0]), Point(PID(1), [0.0, 0.0], [0.0, 0.0]))
    hilbert = Hilbert(pid=>Phonon(2) for pid in [bond.epoint.pid, bond.spoint.pid])
    operators = Operators(
        Operator(+1.0, ID(OID(Index(PID(1), PNID('u', 'x')), [0.0, 0.0], [0.0, 0.0]), OID(Index(PID(1), PNID('u', 'x')), [0.0, 0.0], [0.0, 0.0]))),
        Operator(+1.0, ID(OID(Index(PID(1), PNID('u', 'y')), [0.0, 0.0], [0.0, 0.0]), OID(Index(PID(1), PNID('u', 'y')), [0.0, 0.0], [0.0, 0.0]))),
        Operator(-2.0, ID(OID(Index(PID(1), PNID('u', 'x')), [0.0, 0.0], [0.0, 0.0]), OID(Index(PID(2), PNID('u', 'x')), [0.5, 0.5], [0.0, 0.0]))), 
        Operator(-2.0, ID(OID(Index(PID(1), PNID('u', 'y')), [0.0, 0.0], [0.0, 0.0]), OID(Index(PID(2), PNID('u', 'y')), [0.5, 0.5], [0.0, 0.0]))),
        Operator(+1.0, ID(OID(Index(PID(2), PNID('u', 'x')), [0.5, 0.5], [0.0, 0.0]), OID(Index(PID(2), PNID('u', 'x')), [0.5, 0.5], [0.0, 0.0]))),
        Operator(+1.0, ID(OID(Index(PID(2), PNID('u', 'y')), [0.5, 0.5], [0.0, 0.0]), OID(Index(PID(2), PNID('u', 'y')), [0.5, 0.5], [0.0, 0.0])))
    )
    @test expand(term, bond, hilbert) ≈ operators
end
