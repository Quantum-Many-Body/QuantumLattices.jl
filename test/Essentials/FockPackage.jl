using Test
using StaticArrays: SVector
using QuantumLattices.Essentials.FockPackage
using QuantumLattices.Essentials.Spatials: Bond, Point, PID, rcoord, azimuthd
using QuantumLattices.Essentials.DegreesOfFreedom: Config, OID, Operators, oidtype, script, latexname, isHermitian
using QuantumLattices.Essentials.Terms: Couplings, Subscripts, @subscripts_str, SubID, abbr, otype
using QuantumLattices.Interfaces: ⊗, ⋅, expand, permute, rank
using QuantumLattices.Prerequisites: Float
using QuantumLattices.Prerequisites.Traits: parameternames, isparameterbound, contentnames, getcontent
using QuantumLattices.Mathematics.AlgebraOverFields: ID

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
    @test Dims(fock) == (1, 2, 2)
    @test CartesianIndex(FID{:b}(1, 1, 1), fock) == CartesianIndex(1, 1, 1)
    @test FID(CartesianIndex(1, 1, 1), fock) == FID{:b}(1, 1, 1)
    @test collect(fock) == [FID{:b}(1, 1, 1), FID{:b}(1, 2, 1), FID{:b}(1, 1, 2), FID{:b}(1, 2, 2)]
    @test statistics(fock) == statistics(typeof(fock)) == :b

    @test summary(Fock{:b}(nspin=0, nnambu=1)) == "0-element Fock{:b}"
    @test summary(Fock{:b}(nspin=1, nnambu=1)) == "1-element Fock{:b}"
    @test summary(Fock{:f}(nspin=2, nnambu=1)) == "2-element Fock{:f}"
end

@testset "FIndex" begin
    @test union(PID{Int}, FID{:f}) == FIndex{:f, Int}

    index = FIndex(PID('f', 1), FID{:f}(1, 1, 1))
    @test string(index) == "FIndex{:f}('f', 1, 1, 1, 1)"
    @test index' == replace(index, nambu=2)
    @test statistics(index) == statistics(typeof(index)) == :f
end

@testset "oidtype" begin
    @test oidtype(Fock{:f}, Point{2, PID{Int}}, Val(:info)) == OID{FIndex{:f, Int}, SVector{2, Float}}
    @test oidtype(Fock{:b}, Point{2, PID{Int}}, Val(:info)) == OID{FIndex{:b, Int}, SVector{2, Float}}
end

@testset "latex" begin
    @test script(Val(:site), FIndex{:f}('c', 1, 2, 1, 1)) == 1
    @test script(Val(:orbital), FIndex{:f}('c', 1, 2, 1, 1)) == 2
    @test script(Val(:spinint), FIndex{:f}('c', 1, 2, 3, 1)) == 3
    @test script(Val(:spinsym), FIndex{:f}('c', 1, 2, 2, 1)) == "↑"
    @test script(Val(:spinsym), FIndex{:f}('c', 1, 2, 1, 1)) == "↓"
    @test script(Val(:nambu), FIndex{:f}('c', 1, 2, 3, 1)) == ""
    @test script(Val(:nambu), FIndex{:f}('c', 1, 2, 3, 2)) == "\\dagger"

    @test latexname(FIndex{:f}) == Symbol("FIndex{:f}")
    @test latexname(OID{FIndex{:f}}) == Symbol("OID{FIndex{:f}}")
    @test latexname(FIndex{:b}) == Symbol("FIndex{:b}")
    @test latexname(OID{FIndex{:b}}) == Symbol("OID{FIndex{:b}}")
end

@testset "angle" begin
    @test angle(OID(FIndex{:f}(1, 1, 1, 1, 1), [0.0, 0.0], [1.0, 2.0]), [[1.0, 0.0], [0.0, 1.0]], [0.1, 0.0]) ≈ 2pi*0.1
    @test angle(OID(FIndex{:f}(1, 1, 1, 1, 2), [0.0, 0.0], [1.0, 2.0]), [[1.0, 0.0], [0.0, 1.0]], [0.0, 0.2]) ≈ -2pi*0.4
end

@testset "FockOperator" begin
    id₁ = OID(FIndex{:f}(1, 2, 1, 1, 2), SVector(0.5, 0.0), SVector(0.0, 0.0))
    id₂ = OID(FIndex{:f}(1, 2, 1, 1, 1), SVector(0.5, 0.0), SVector(0.0, 0.0))
    id₃ = OID(FIndex{:f}(1, 1, 1, 2, 2), SVector(0.0, 0.0), SVector(0.0, 0.0))
    id₄ = OID(FIndex{:f}(1, 1, 1, 2, 1), SVector(0.0, 0.0), SVector(0.0, 0.0))

    opt = FOperator(1.0, ID(id₁, id₂))
    @test opt|>isnormalordered
    @test opt|>statistics == opt|>typeof|>statistics == :f

    opt = FOperator(1.0, ID(id₁, id₂, id₃, id₄))
    @test opt|>isnormalordered == false
    @test repr(opt) == "c^{\\dagger}_{2, 1, ↓}c^{}_{2, 1, ↓}c^{\\dagger}_{1, 1, ↑}c^{}_{1, 1, ↑}"

    op₁ = FOperator(1.5, ID(id₁, id₂))
    op₂ = FOperator(2.0, ID(id₂, id₁))
    @test op₁*op₂ == nothing

    op₁ = FOperator(1.5, ID(id₁, id₂))
    op₂ = FOperator(2.0, ID(id₁, id₂))
    @test op₁*op₂ == FOperator(3.0, ID(id₁, id₂, id₁, id₂))

    @test permute(FOperator, id₁, id₂) == (FOperator(1), FOperator(-1, ID(id₂, id₁)))
    @test permute(FOperator, id₂, id₁) == (FOperator(1), FOperator(-1, ID(id₁, id₂)))
    @test permute(FOperator, id₁, id₄) == (FOperator(-1, ID(id₄, id₁)),)
    @test permute(FOperator, id₄, id₁) == (FOperator(-1, ID(id₁, id₄)),)


    id₁ = OID(FIndex{:b}(1, 2, 1, 1, 2), SVector(0.5, 0.0), SVector(0.0, 0.0))
    id₂ = OID(FIndex{:b}(1, 2, 1, 1, 1), SVector(0.5, 0.0), SVector(0.0, 0.0))
    id₃ = OID(FIndex{:b}(1, 1, 1, 2, 2), SVector(0.0, 0.0), SVector(0.0, 0.0))
    id₄ = OID(FIndex{:b}(1, 1, 1, 2, 1), SVector(0.0, 0.0), SVector(0.0, 0.0))

    opt = BOperator(1.0, ID(id₁, id₂))
    @test opt|>statistics == opt|>typeof|>statistics == :b
    @test repr(opt) == "b^{\\dagger}_{2, 1, ↓}b^{}_{2, 1, ↓}"

    @test permute(BOperator, id₁, id₂) == (BOperator(+1), BOperator(1, ID(id₂, id₁)))
    @test permute(BOperator, id₂, id₁) == (BOperator(-1), BOperator(1, ID(id₁, id₂)))
    @test permute(BOperator, id₁, id₄) == (BOperator(1, ID(id₄, id₁)),)
    @test permute(BOperator, id₄, id₁) == (BOperator(1, ID(id₁, id₄)),)
end

@testset "FockCoupling" begin
    @test parameternames(FockCoupling) == (:value, :atoms, :nambus, :orbitals, :spins, :id)
    @test isparameterbound(FockCoupling, Val(:atoms), Tuple{Int, Int}) == false
    @test isparameterbound(FockCoupling, Val(:atoms), Tuple{Any, Any}) == true
    @test isparameterbound(FockCoupling, Val(:nambus), Tuple{Int, Int}) == false
    @test isparameterbound(FockCoupling, Val(:nambus), Tuple{Any, Any}) == true
    @test isparameterbound(FockCoupling, Val(:orbitals), Subscripts) == true
    @test isparameterbound(FockCoupling, Val(:spins), Subscripts) == true
    @test contentnames(FockCoupling) == (:value, :id, :orbitals, :spins)

    obs, sps = subscripts"[1 1]", subscripts"[σ₁ σ₂](σ₁ + σ₂ == 1)"
    fc = FockCoupling(1.0im, (FCID((1, 1), (2, 1)), SubID(obs), SubID(sps)), obs, sps)
    @test getcontent(fc, Val(:id)) == (FCID(fc.atoms, fc.nambus), SubID(fc.orbitals), SubID(fc.spins))
    @test rank(fc) == rank(typeof(fc)) == 2

    @test string(FockCoupling{2}(1.0)) == "FockCoupling{2}(value=1.0)"
    @test string(FockCoupling{2}(1.0, atoms=(1, 1))) == "FockCoupling{2}(value=1.0, atoms=[1 1])"
    @test string(FockCoupling{2}(1.0, atoms=(1, 1), spins=(1, 2))) == "FockCoupling{2}(value=1.0, atoms=[1 1], spins=[1 2])"
    @test repr(FockCoupling{2}(2.0)) == "2.0 {2}"

    fc₁ = FockCoupling{2}(1.5, atoms=(2, 1), spins=subscripts"[x 1]")
    fc₂ = FockCoupling{2}(2.0, atoms=(1, 2), orbitals=subscripts"[x y](x < y)")
    @test repr(fc₁) == "1.5 sl[2 1] ⊗ sp[x 1]"
    @test repr(fc₂) == "2.0 sl[1 2] ⊗ ob[x y](x < y)"
    fc = fc₁ * fc₂
    @test repr(fc) == "3.0 sl[2 1 1 2] ⊗ ob[* *; x y](:, x < y) ⊗ sp[x 1; * *]"

    fc₁ = FockCoupling{2}(1.5, spins=subscripts"[x 1]")
    fc₂ = FockCoupling{2}(2.0, orbitals=subscripts"[x y](x < y)")
    fc = fc₁ ⊗ fc₂
    @test repr(fc) == "3.0 ob[x y](x < y) ⊗ sp[x 1]"

    fc₁ = FockCoupling{2}(1.5, atoms=(2, 1))
    fc₂ = FockCoupling{2}(2.0, atoms=(1, 2))
    fc = fc₁ ⋅ fc₂
    @test repr(fc) == "3.0 sl[2 2]"

    fc = FockCoupling{2}(2.0, atoms=(1, 1))
    point = Point(PID(1, 1), SVector(0.0, 0.0), SVector(0.0, 0.0))
    fock = Fock{:f}(atom=2, norbital=2, nspin=2, nnambu=2)
    @test collect(expand(fc, (point, point), (fock, fock), Val(:info))) == []

    fc = FockCoupling{2}(2.0, atoms=(1, 2), orbitals=(1, 2), nambus=(2, 1))
    p₁, p₂ = Point(PID(1, 1), SVector(0.0), SVector(0.0)), Point(PID(1, 2), SVector(0.5), SVector(0.0))
    f₁, f₂ = Fock{:f}(atom=1, norbital=2, nspin=2, nnambu=2), Fock{:f}(atom=2, norbital=2, nspin=2, nnambu=2)
    ex = expand(fc, (p₁, p₂), (f₁, f₂), Val(:info))
    @test Dims(ex) == (1, 2)
    @test collect(ex) == [
        (2.0, ID(OID(FIndex{:f}(1, 1, 1, 1, 2), SVector(0.0), SVector(0.0)), OID(FIndex{:f}(1, 2, 2, 1, 1), SVector(0.5), SVector(0.0)))),
        (2.0, ID(OID(FIndex{:f}(1, 1, 1, 2, 2), SVector(0.0), SVector(0.0)), OID(FIndex{:f}(1, 2, 2, 2, 1), SVector(0.5), SVector(0.0))))
    ]

    fc = FockCoupling{4}(2.0, spins=(2, 2, 1, 1), nambus=(2, 1, 2, 1))
    point = Point(PID(1, 1), SVector(0.0), SVector(0.0))
    fock = Fock{:b}(atom=1, norbital=2, nspin=2, nnambu=2)
    ex = expand(fc, (point, point, point, point), (fock, fock, fock, fock), Val(:info))
    @test eltype(ex) == Tuple{Float, ID{OID{FIndex{:b, Int}, SVector{1, Float}}, 4}}
    @test Dims(ex) == (2, 1)
    @test collect(ex) == [
        (2.0, ID(OID(FIndex{:b}(1, 1, 1, 2, 2), SVector(0.0), SVector(0.0)),
                 OID(FIndex{:b}(1, 1, 1, 2, 1), SVector(0.0), SVector(0.0)),
                 OID(FIndex{:b}(1, 1, 1, 1, 2), SVector(0.0), SVector(0.0)),
                 OID(FIndex{:b}(1, 1, 1, 1, 1), SVector(0.0), SVector(0.0))
                 )),
        (2.0, ID(OID(FIndex{:b}(1, 1, 2, 2, 2), SVector(0.0), SVector(0.0)),
                 OID(FIndex{:b}(1, 1, 2, 2, 1), SVector(0.0), SVector(0.0)),
                 OID(FIndex{:b}(1, 1, 2, 1, 2), SVector(0.0), SVector(0.0)),
                 OID(FIndex{:b}(1, 1, 2, 1, 1), SVector(0.0), SVector(0.0))
                 ))
    ]

    fc = FockCoupling{4}(2.0, orbitals=subscripts"[α α β β](α < β)", spins=(2, 1, 1, 2), nambus=(2, 2, 1, 1))
    point = Point(PID(1, 1), SVector(0.5), SVector(0.0))
    fock = Fock{:f}(atom=1, norbital=3, nspin=2, nnambu=2)
    ex = expand(fc, (point, point, point, point), (fock, fock, fock, fock), Val(:info))
    @test Dims(ex) == (3, 1)
    @test collect(ex) == [
        (2.0, ID(OID(FIndex{:f}(1, 1, 1, 2, 2), SVector(0.5), SVector(0.0)),
                 OID(FIndex{:f}(1, 1, 1, 1, 2), SVector(0.5), SVector(0.0)),
                 OID(FIndex{:f}(1, 1, 2, 1, 1), SVector(0.5), SVector(0.0)),
                 OID(FIndex{:f}(1, 1, 2, 2, 1), SVector(0.5), SVector(0.0))
                 )),
        (2.0, ID(OID(FIndex{:f}(1, 1, 1, 2, 2), SVector(0.5), SVector(0.0)),
                 OID(FIndex{:f}(1, 1, 1, 1, 2), SVector(0.5), SVector(0.0)),
                 OID(FIndex{:f}(1, 1, 3, 1, 1), SVector(0.5), SVector(0.0)),
                 OID(FIndex{:f}(1, 1, 3, 2, 1), SVector(0.5), SVector(0.0))
                 )),
        (2.0, ID(OID(FIndex{:f}(1, 1, 2, 2, 2), SVector(0.5), SVector(0.0)),
                 OID(FIndex{:f}(1, 1, 2, 1, 2), SVector(0.5), SVector(0.0)),
                 OID(FIndex{:f}(1, 1, 3, 1, 1), SVector(0.5), SVector(0.0)),
                 OID(FIndex{:f}(1, 1, 3, 2, 1), SVector(0.5), SVector(0.0))
                 ))
    ]

    fc₁ = FockCoupling{2}(+1.0, spins=(2, 2), nambus=(2, 1))
    fc₂ = FockCoupling{2}(-1.0, spins=(1, 1), nambus=(2, 1))
    point = Point(PID(1, 1), SVector(0.0), SVector(0.0))
    fock = Fock{:f}(atom=1, norbital=2, nspin=2, nnambu=2)
    ex = expand(fc₁*fc₂, (point, point, point, point), (fock, fock, fock, fock), Val(:info))
    @test Dims(ex) == (4, 1)
    @test collect(ex) == [
        (-1.0, ID(OID(FIndex{:f}(1, 1, 1, 2, 2), SVector(0.0), SVector(0.0)),
                  OID(FIndex{:f}(1, 1, 1, 2, 1), SVector(0.0), SVector(0.0)),
                  OID(FIndex{:f}(1, 1, 1, 1, 2), SVector(0.0), SVector(0.0)),
                  OID(FIndex{:f}(1, 1, 1, 1, 1), SVector(0.0), SVector(0.0))
                  )),
        (-1.0, ID(OID(FIndex{:f}(1, 1, 2, 2, 2), SVector(0.0), SVector(0.0)),
                  OID(FIndex{:f}(1, 1, 2, 2, 1), SVector(0.0), SVector(0.0)),
                  OID(FIndex{:f}(1, 1, 1, 1, 2), SVector(0.0), SVector(0.0)),
                  OID(FIndex{:f}(1, 1, 1, 1, 1), SVector(0.0), SVector(0.0))
                  )),
        (-1.0, ID(OID(FIndex{:f}(1, 1, 1, 2, 2), SVector(0.0), SVector(0.0)),
                  OID(FIndex{:f}(1, 1, 1, 2, 1), SVector(0.0), SVector(0.0)),
                  OID(FIndex{:f}(1, 1, 2, 1, 2), SVector(0.0), SVector(0.0)),
                  OID(FIndex{:f}(1, 1, 2, 1, 1), SVector(0.0), SVector(0.0))
                  )),
        (-1.0, ID(OID(FIndex{:f}(1, 1, 2, 2, 2), SVector(0.0), SVector(0.0)),
                  OID(FIndex{:f}(1, 1, 2, 2, 1), SVector(0.0), SVector(0.0)),
                  OID(FIndex{:f}(1, 1, 2, 1, 2), SVector(0.0), SVector(0.0)),
                  OID(FIndex{:f}(1, 1, 2, 1, 1), SVector(0.0), SVector(0.0))
                  ))
    ]
end

@testset "σ⁰" begin
    @test σ⁰("sp") == FockCoupling{2}(1, spins=(1, 1)) + FockCoupling{2}(1, spins=(2, 2))
    @test σ⁰("ob") == FockCoupling{2}(1, orbitals=(1, 1)) + FockCoupling{2}(1, orbitals=(2, 2))
    @test σ⁰("sl") == FockCoupling{2}(1, atoms=(1, 1)) + FockCoupling{2}(1, atoms=(2, 2))
    @test σ⁰("ph") == FockCoupling{2}(1, nambus=(1, 2)) + FockCoupling{2}(1, nambus=(2, 1))
end

@testset "σˣ" begin
    @test σˣ("sp") == FockCoupling{2}(1, spins=(1, 2)) + FockCoupling{2}(1, spins=(2, 1))
    @test σˣ("ob") == FockCoupling{2}(1, orbitals=(1, 2)) + FockCoupling{2}(1, orbitals=(2, 1))
    @test σˣ("sl") == FockCoupling{2}(1, atoms=(1, 2)) + FockCoupling{2}(1, atoms=(2, 1))
    @test σˣ("ph") == FockCoupling{2}(1, nambus=(1, 1)) + FockCoupling{2}(1, nambus=(2, 2))
end

@testset "σʸ" begin
    @test σʸ("sp") == FockCoupling{2}(1im, spins=(1, 2)) + FockCoupling{2}(-1im, spins=(2, 1))
    @test σʸ("ob") == FockCoupling{2}(1im, orbitals=(1, 2)) + FockCoupling{2}(-1im, orbitals=(2, 1))
    @test σʸ("sl") == FockCoupling{2}(1im, atoms=(1, 2)) + FockCoupling{2}(-1im, atoms=(2, 1))
    @test σʸ("ph") == FockCoupling{2}(1im, nambus=(1, 1)) + FockCoupling{2}(-1im, nambus=(2, 2))
end

@testset "σᶻ" begin
    @test σᶻ("sp") == FockCoupling{2}(-1, spins=(1, 1)) + FockCoupling{2}(1, spins=(2, 2))
    @test σᶻ("ob") == FockCoupling{2}(-1, orbitals=(1, 1)) + FockCoupling{2}(1, orbitals=(2, 2))
    @test σᶻ("sl") == FockCoupling{2}(-1, atoms=(1, 1)) + FockCoupling{2}(1, atoms=(2, 2))
    @test σᶻ("ph") == FockCoupling{2}(-1, nambus=(1, 2)) + FockCoupling{2}(1, nambus=(2, 1))
end

@testset "σ⁺" begin
    @test σ⁺("sp") == Couplings(FockCoupling{2}(1, spins=(2, 1)))
    @test σ⁺("ob") == Couplings(FockCoupling{2}(1, orbitals=(2, 1)))
    @test σ⁺("sl") == Couplings(FockCoupling{2}(1, atoms=(2, 1)))
    @test σ⁺("ph") == Couplings(FockCoupling{2}(1, nambus=(2, 2)))
end

@testset "σ⁻" begin
    @test σ⁻("sp") == Couplings(FockCoupling{2}(1, spins=(1, 2)))
    @test σ⁻("ob") == Couplings(FockCoupling{2}(1, orbitals=(1, 2)))
    @test σ⁻("sl") == Couplings(FockCoupling{2}(1, atoms=(1, 2)))
    @test σ⁻("ph") == Couplings(FockCoupling{2}(1, nambus=(1, 1)))
end

@testset "fockcoupling" begin
    fc = fc"1.0 sl[1 1 1 1] ⊗ ph[2 1 2 1] ⊗ ob[α α β β](α < β) ⊗ sp[σ γ σ γ](σ ≠ γ)"
    @test repr(fc) == "1.0 sl[1 1 1 1] ⊗ ph[2 1 2 1] ⊗ ob[α α β β](α < β) ⊗ sp[σ γ σ γ](σ ≠ γ)"

    fc = fc"1.0 sl[1 1 1 1] ⊗ ph[2 1 2 1] ⊗ ob[α α β β] ⊗ sp[σ γ σ γ]"
    @test repr(fc) == "1.0 sl[1 1 1 1] ⊗ ph[2 1 2 1] ⊗ ob[α α β β] ⊗ sp[σ γ σ γ]"

    fc = fc"1.0 sl[1 1 1 1] ⊗ ph[2 1 2 1] ⊗ ob[α α β β](α < β) ⊗ sp[σ γ σ γ]"
    @test repr(fc) == "1.0 sl[1 1 1 1] ⊗ ph[2 1 2 1] ⊗ ob[α α β β](α < β) ⊗ sp[σ γ σ γ]"

    fc = fc"1.0 sl[1 1 1 1] ⊗ ph[2 1 2 1] ⊗ ob[α α β β](α < β)"
    @test repr(fc) == "1.0 sl[1 1 1 1] ⊗ ph[2 1 2 1] ⊗ ob[α α β β](α < β)"

    fc = fc"1.0 sl[1 1 1 1] ⊗ ob[α α β β](α < β) ⊗ sp[2 1 2 1]"
    @test repr(fc) == "1.0 sl[1 1 1 1] ⊗ ob[α α β β](α < β) ⊗ sp[2 1 2 1]"

    fc = fc"1.0 sl[1 1 1 1] ⊗ ph[2 1 2 1] ⊗ ob[1 1 1 1]"
    @test repr(fc) == "1.0 sl[1 1 1 1] ⊗ ph[2 1 2 1] ⊗ ob[1 1 1 1]"

    fc = fc"1.0 sl[1 1 1 1] ⊗ ph[2 1 2 1]"
    @test repr(fc) == "1.0 sl[1 1 1 1] ⊗ ph[2 1 2 1]"

    fc = fc"1.0im {2}"
    @test repr(fc) == "1.0im {2}"
end

@testset "fockcouplings" begin
    @test σ⁰"sp" == σ⁰("sp")
    @test σ⁰"ob" == σ⁰("ob")
    @test σ⁰"sl" == σ⁰("sl")
    @test σ⁰"ph" == σ⁰("ph")

    @test σˣ"sp" == σˣ("sp")
    @test σˣ"ob" == σˣ("ob")
    @test σˣ"sl" == σˣ("sl")
    @test σˣ"ph" == σˣ("ph")

    @test σʸ"sp" == σʸ("sp")
    @test σʸ"ob" == σʸ("ob")
    @test σʸ"sl" == σʸ("sl")
    @test σʸ"ph" == σʸ("ph")

    @test σᶻ"sp" == σᶻ("sp")
    @test σᶻ"ob" == σᶻ("ob")
    @test σᶻ"sl" == σᶻ("sl")
    @test σᶻ"ph" == σᶻ("ph")

    @test σ⁺"sp" == σ⁺("sp")
    @test σ⁺"ob" == σ⁺("ob")
    @test σ⁺"sl" == σ⁺("sl")
    @test σ⁺"ph" == σ⁺("ph")

    @test σ⁻"sp" == σ⁻("sp")
    @test σ⁻"ob" == σ⁻("ob")
    @test σ⁻"sl" == σ⁻("sl")
    @test σ⁻"ph" == σ⁻("ph")
end

@testset "Onsite" begin
    @test abbr(Onsite) == :st
    @test isnothing(isHermitian(Onsite))

    point = Point(PID('a', 1), (0.5, 0.5), (0.0, 0.0))
    config = Config{Fock{:f}}(pid->Fock{:f}(atom=pid.site%2, norbital=2, nspin=2, nnambu=2), [point.pid])

    term = Onsite(:mu, 1.5, couplings=σˣ("sp")⊗σᶻ("ob"), modulate=true)
    operators = Operators(
        FOperator(+1.5, ID(OID(FIndex{:f}('a', 1, 2, 1, 2), [0.5, 0.5], [0.0, 0.0]), OID(FIndex{:f}('a', 1, 2, 2, 1), [0.5, 0.5], [0.0, 0.0]))),
        FOperator(-1.5, ID(OID(FIndex{:f}('a', 1, 1, 1, 2), [0.5, 0.5], [0.0, 0.0]), OID(FIndex{:f}('a', 1, 1, 2, 1), [0.5, 0.5], [0.0, 0.0])))
    )
    @test expand(term, point, config, true) == operators
    @test expand(term, point, config, false) == operators+operators'

    term = Onsite(:mu, 1.5, couplings=σᶻ("sp")⊗σᶻ("ob"), modulate=true)
    operators = Operators(
        FOperator(-0.75, ID(OID(FIndex{:f}('a', 1, 2, 1, 2), [0.5, 0.5], [0.0, 0.0]), OID(FIndex{:f}('a', 1, 2, 1, 1), [0.5, 0.5], [0.0, 0.0]))),
        FOperator(-0.75, ID(OID(FIndex{:f}('a', 1, 1, 2, 2), [0.5, 0.5], [0.0, 0.0]), OID(FIndex{:f}('a', 1, 1, 2, 1), [0.5, 0.5], [0.0, 0.0]))),
        FOperator(+0.75, ID(OID(FIndex{:f}('a', 1, 2, 2, 2), [0.5, 0.5], [0.0, 0.0]), OID(FIndex{:f}('a', 1, 2, 2, 1), [0.5, 0.5], [0.0, 0.0]))),
        FOperator(+0.75, ID(OID(FIndex{:f}('a', 1, 1, 1, 2), [0.5, 0.5], [0.0, 0.0]), OID(FIndex{:f}('a', 1, 1, 1, 1), [0.5, 0.5], [0.0, 0.0])))
    )
    @test expand(term, point, config, true) == operators
    @test expand(term, point, config, false) == operators+operators'
end

@testset "Hopping" begin
    bond = Bond(1, Point(PID('a', 1), (0.5, 0.5), (0.0, 0.0)), Point(PID('b', 2), (0.0, 0.0), (0.0, 0.0)))
    config = Config{Fock{:f}}(pid->Fock{:f}(atom=pid.site%2, norbital=2, nspin=2, nnambu=2), [bond.spoint.pid, bond.epoint.pid])
    term = Hopping(:t, 1.5, 1)
    operators = Operators(
        FOperator(1.5, ID(OID(FIndex{:f}('b', 2, 2, 2, 2), [0.0, 0.0], [0.0, 0.0]), OID(FIndex{:f}('a', 1, 2, 2, 1), [0.5, 0.5], [0.0, 0.0]))),
        FOperator(1.5, ID(OID(FIndex{:f}('b', 2, 2, 1, 2), [0.0, 0.0], [0.0, 0.0]), OID(FIndex{:f}('a', 1, 2, 1, 1), [0.5, 0.5], [0.0, 0.0]))),
        FOperator(1.5, ID(OID(FIndex{:f}('b', 2, 1, 1, 2), [0.0, 0.0], [0.0, 0.0]), OID(FIndex{:f}('a', 1, 1, 1, 1), [0.5, 0.5], [0.0, 0.0]))),
        FOperator(1.5, ID(OID(FIndex{:f}('b', 2, 1, 2, 2), [0.0, 0.0], [0.0, 0.0]), OID(FIndex{:f}('a', 1, 1, 2, 1), [0.5, 0.5], [0.0, 0.0])))
    )
    @test term|>abbr == :hp
    @test term|>isHermitian == false
    @test expand(term, bond, config, true) == operators
    @test expand(term, bond, config, false) == operators+operators'
end

@testset "Pairing" begin
    bond = Bond(1, Point(PID('a', 1), (0.5, 0.5), (0.0, 0.0)), Point(PID('b', 2), (0.0, 0.0), (0.0, 0.0)))
    config = Config{Fock{:f}}(pid->Fock{:f}(atom=pid.site%2, norbital=1, nspin=2, nnambu=2), [bond.spoint.pid, bond.epoint.pid])
    term = Pairing(:Δ, 1.5, 1, couplings=FockCoupling{2}(spins=(2, 2)), amplitude=bond->(bond|>rcoord|>azimuthd ≈ 45 ? 1 : -1))
    operators = Operators(
        FOperator(-1.5, ID(OID(FIndex{:f}('b', 2, 1, 2, 1), [0.0, 0.0], [0.0, 0.0]), OID(FIndex{:f}('a', 1, 1, 2, 1), [0.5, 0.5], [0.0, 0.0]))),
        FOperator(+1.5, ID(OID(FIndex{:f}('a', 1, 1, 2, 1), [0.5, 0.5], [0.0, 0.0]), OID(FIndex{:f}('b', 2, 1, 2, 1), [0.0, 0.0], [0.0, 0.0])))
    )
    @test term|>abbr == :pr
    @test term|>isHermitian == false
    @test expand(term, bond, config, true) == operators
    @test expand(term, bond, config, false) == operators+operators'

    point = Point(PID('a', 1), (0.5, 0.5), (0.0, 0.0))
    config = Config{Fock{:f}}(pid->Fock{:f}(atom=pid.site%2, norbital=1, nspin=2, nnambu=2), [point.pid])
    term = Pairing(:Δ, 1.5, 0, couplings=FockCoupling{2}(spins=(2, 1))-FockCoupling{2}(spins=(1, 2)))
    operators = Operators(
        FOperator(+1.5, ID(OID(FIndex{:f}('a', 1, 1, 2, 1), [0.5, 0.5], [0.0, 0.0]), OID(FIndex{:f}('a', 1, 1, 1, 1), [0.5, 0.5], [0.0, 0.0]))),
        FOperator(-1.5, ID(OID(FIndex{:f}('a', 1, 1, 1, 1), [0.5, 0.5], [0.0, 0.0]), OID(FIndex{:f}('a', 1, 1, 2, 1), [0.5, 0.5], [0.0, 0.0])))
    )
    @test term|>abbr == :pr
    @test expand(term, point, config, true) == operators
    @test expand(term, point, config, false) == operators+operators'
end

@testset "Hubbard" begin
    point = Point(PID('a', 1), (0.5, 0.5), (0.0, 0.0))
    config = Config{Fock{:f}}(pid->Fock{:f}(atom=pid.site%2, norbital=2, nspin=2, nnambu=2), [point.pid])
    term = Hubbard(:H, 2.5)
    operators = Operators(
        FOperator(1.25, ID(
            OID(FIndex{:f}('a', 1, 1, 2, 2), [0.5, 0.5], [0.0, 0.0]),
            OID(FIndex{:f}('a', 1, 1, 2, 1), [0.5, 0.5], [0.0, 0.0]),
            OID(FIndex{:f}('a', 1, 1, 1, 2), [0.5, 0.5], [0.0, 0.0]),
            OID(FIndex{:f}('a', 1, 1, 1, 1), [0.5, 0.5], [0.0, 0.0])
            )),
        FOperator(1.25, ID(
            OID(FIndex{:f}('a', 1, 2, 2, 2), [0.5, 0.5], [0.0, 0.0]),
            OID(FIndex{:f}('a', 1, 2, 2, 1), [0.5, 0.5], [0.0, 0.0]),
            OID(FIndex{:f}('a', 1, 2, 1, 2), [0.5, 0.5], [0.0, 0.0]),
            OID(FIndex{:f}('a', 1, 2, 1, 1), [0.5, 0.5], [0.0, 0.0])
            ))
    )
    @test term|>abbr == :hb
    @test term|>isHermitian == true
    @test expand(term, point, config, true) == operators
    @test expand(term, point, config, false) == operators*2
end

@testset "InterOrbitalInterSpin" begin
    point = Point(PID('a', 1), (0.5, 0.5), (0.0, 0.0))
    config = Config{Fock{:f}}(pid->Fock{:f}(atom=pid.site%2, norbital=2, nspin=2, nnambu=2), [point.pid])
    term = InterOrbitalInterSpin(:H, 2.5)
    operators = Operators(
        FOperator(1.25, ID(
            OID(FIndex{:f}('a', 1, 1, 2, 2), [0.5, 0.5], [0.0, 0.0]),
            OID(FIndex{:f}('a', 1, 1, 2, 1), [0.5, 0.5], [0.0, 0.0]),
            OID(FIndex{:f}('a', 1, 2, 1, 2), [0.5, 0.5], [0.0, 0.0]),
            OID(FIndex{:f}('a', 1, 2, 1, 1), [0.5, 0.5], [0.0, 0.0])
            )),
        FOperator(1.25, ID(
            OID(FIndex{:f}('a', 1, 1, 1, 2), [0.5, 0.5], [0.0, 0.0]),
            OID(FIndex{:f}('a', 1, 1, 1, 1), [0.5, 0.5], [0.0, 0.0]),
            OID(FIndex{:f}('a', 1, 2, 2, 2), [0.5, 0.5], [0.0, 0.0]),
            OID(FIndex{:f}('a', 1, 2, 2, 1), [0.5, 0.5], [0.0, 0.0])
            ))
    )
    @test term|>abbr == :nons
    @test term|>isHermitian == true
    @test expand(term, point, config, true) == operators
    @test expand(term, point, config, false) == operators*2
end

@testset "InterOrbitalIntraSpin" begin
    point = Point(PID('a', 1), (0.5, 0.5), (0.0, 0.0))
    config = Config{Fock{:f}}(pid->Fock{:f}(atom=pid.site%2, norbital=2, nspin=2, nnambu=2), [point.pid])
    term = InterOrbitalIntraSpin(:H, 2.5)
    operators = Operators(
        FOperator(1.25, ID(
            OID(FIndex{:f}('a', 1, 1, 1, 2), [0.5, 0.5], [0.0, 0.0]),
            OID(FIndex{:f}('a', 1, 1, 1, 1), [0.5, 0.5], [0.0, 0.0]),
            OID(FIndex{:f}('a', 1, 2, 1, 2), [0.5, 0.5], [0.0, 0.0]),
            OID(FIndex{:f}('a', 1, 2, 1, 1), [0.5, 0.5], [0.0, 0.0])
            )),
        FOperator(1.25, ID(
            OID(FIndex{:f}('a', 1, 1, 2, 2), [0.5, 0.5], [0.0, 0.0]),
            OID(FIndex{:f}('a', 1, 1, 2, 1), [0.5, 0.5], [0.0, 0.0]),
            OID(FIndex{:f}('a', 1, 2, 2, 2), [0.5, 0.5], [0.0, 0.0]),
            OID(FIndex{:f}('a', 1, 2, 2, 1), [0.5, 0.5], [0.0, 0.0])
            ))
    )
    @test term|>abbr == :noes
    @test term|>isHermitian == true
    @test expand(term, point, config, true) == operators
    @test expand(term, point, config, false) == operators*2
end

@testset "SpinFlip" begin
    point = Point(PID('a', 1), (0.5, 0.5), (0.0, 0.0))
    config = Config{Fock{:f}}(pid->Fock{:f}(atom=pid.site%2, norbital=2, nspin=2, nnambu=2), [point.pid])
    term = SpinFlip(:H, 2.5)
    operators = Operators(
        FOperator(2.5, ID(
            OID(FIndex{:f}('a', 1, 1, 2, 2), [0.5, 0.5], [0.0, 0.0]),
            OID(FIndex{:f}('a', 1, 2, 1, 2), [0.5, 0.5], [0.0, 0.0]),
            OID(FIndex{:f}('a', 1, 1, 1, 1), [0.5, 0.5], [0.0, 0.0]),
            OID(FIndex{:f}('a', 1, 2, 2, 1), [0.5, 0.5], [0.0, 0.0])
            ))
    )
    @test term|>abbr == :sf
    @test term|>isHermitian == false
    @test expand(term, point, config, true) == operators
    @test expand(term, point, config, false) == operators+operators'
end

@testset "PairHopping" begin
    point = Point(PID('a', 1), (0.5, 0.5), (0.0, 0.0))
    config = Config{Fock{:f}}(pid->Fock{:f}(atom=pid.site%2, norbital=2, nspin=2, nnambu=2), [point.pid])
    term = PairHopping(:H, 2.5)
    operators = Operators(
        FOperator(2.5, ID(
            OID(FIndex{:f}('a', 1, 1, 2, 2), [0.5, 0.5], [0.0, 0.0]),
            OID(FIndex{:f}('a', 1, 1, 1, 2), [0.5, 0.5], [0.0, 0.0]),
            OID(FIndex{:f}('a', 1, 2, 1, 1), [0.5, 0.5], [0.0, 0.0]),
            OID(FIndex{:f}('a', 1, 2, 2, 1), [0.5, 0.5], [0.0, 0.0])
            ))
    )
    @test term|>abbr == :ph
    @test term|>isHermitian == false
    @test expand(term, point, config, true) == operators
    @test expand(term, point, config, false) == operators+operators'
end

@testset "Coulomb" begin
    bond = Bond(1, Point(PID('a', 1), (0.5, 0.5), (0.0, 0.0)), Point(PID('b', 2), (0.0, 0.0), (0.0, 0.0)))
    config = Config{Fock{:f}}(pid->Fock{:f}(atom=pid.site%2, norbital=1, nspin=2, nnambu=2), [bond.spoint.pid, bond.epoint.pid])

    term = Coulomb(:V, 2.5, 1, couplings=σᶻ("sp")*σᶻ("sp"))
    operators = Operators(
        FOperator(-1.25, ID(
            OID(FIndex{:f}('b', 2, 1, 1, 2), [0.0, 0.0], [0.0, 0.0]),
            OID(FIndex{:f}('b', 2, 1, 1, 1), [0.0, 0.0], [0.0, 0.0]),
            OID(FIndex{:f}('a', 1, 1, 2, 2), [0.5, 0.5], [0.0, 0.0]),
            OID(FIndex{:f}('a', 1, 1, 2, 1), [0.5, 0.5], [0.0, 0.0])
            )),
        FOperator(+1.25, ID(
            OID(FIndex{:f}('b', 2, 1, 1, 2), [0.0, 0.0], [0.0, 0.0]),
            OID(FIndex{:f}('b', 2, 1, 1, 1), [0.0, 0.0], [0.0, 0.0]),
            OID(FIndex{:f}('a', 1, 1, 1, 2), [0.5, 0.5], [0.0, 0.0]),
            OID(FIndex{:f}('a', 1, 1, 1, 1), [0.5, 0.5], [0.0, 0.0])
            )),
        FOperator(-1.25, ID(
            OID(FIndex{:f}('b', 2, 1, 2, 2), [0.0, 0.0], [0.0, 0.0]),
            OID(FIndex{:f}('b', 2, 1, 2, 1), [0.0, 0.0], [0.0, 0.0]),
            OID(FIndex{:f}('a', 1, 1, 1, 2), [0.5, 0.5], [0.0, 0.0]),
            OID(FIndex{:f}('a', 1, 1, 1, 1), [0.5, 0.5], [0.0, 0.0])
            )),
        FOperator(+1.25, ID(
            OID(FIndex{:f}('b', 2, 1, 2, 2), [0.0, 0.0], [0.0, 0.0]),
            OID(FIndex{:f}('b', 2, 1, 2, 1), [0.0, 0.0], [0.0, 0.0]),
            OID(FIndex{:f}('a', 1, 1, 2, 2), [0.5, 0.5], [0.0, 0.0]),
            OID(FIndex{:f}('a', 1, 1, 2, 1), [0.5, 0.5], [0.0, 0.0])
            ))
    )
    @test term|>abbr == :cl
    @test isnothing(term|>isHermitian)
    @test expand(term, bond, config, true) == operators
    @test expand(term, bond, config, false) == operators+operators'

    term = Coulomb(:V, 2.5, 1, couplings=σˣ("sp")*σᶻ("sp"))
    operators = Operators(
        FOperator(-2.5, ID(
            OID(FIndex{:f}('b', 2, 1, 2, 2), [0.0, 0.0], [0.0, 0.0]),
            OID(FIndex{:f}('b', 2, 1, 1, 1), [0.0, 0.0], [0.0, 0.0]),
            OID(FIndex{:f}('a', 1, 1, 1, 2), [0.5, 0.5], [0.0, 0.0]),
            OID(FIndex{:f}('a', 1, 1, 1, 1), [0.5, 0.5], [0.0, 0.0])
            )),
        FOperator(+2.5, ID(
            OID(FIndex{:f}('b', 2, 1, 1, 2), [0.0, 0.0], [0.0, 0.0]),
            OID(FIndex{:f}('b', 2, 1, 2, 1), [0.0, 0.0], [0.0, 0.0]),
            OID(FIndex{:f}('a', 1, 1, 2, 2), [0.5, 0.5], [0.0, 0.0]),
            OID(FIndex{:f}('a', 1, 1, 2, 1), [0.5, 0.5], [0.0, 0.0])
            )),
        FOperator(+2.5, ID(
            OID(FIndex{:f}('b', 2, 1, 2, 2), [0.0, 0.0], [0.0, 0.0]),
            OID(FIndex{:f}('b', 2, 1, 1, 1), [0.0, 0.0], [0.0, 0.0]),
            OID(FIndex{:f}('a', 1, 1, 2, 2), [0.5, 0.5], [0.0, 0.0]),
            OID(FIndex{:f}('a', 1, 1, 2, 1), [0.5, 0.5], [0.0, 0.0])
            )),
        FOperator(-2.5, ID(
            OID(FIndex{:f}('b', 2, 1, 1, 2), [0.0, 0.0], [0.0, 0.0]),
            OID(FIndex{:f}('b', 2, 1, 2, 1), [0.0, 0.0], [0.0, 0.0]),
            OID(FIndex{:f}('a', 1, 1, 1, 2), [0.5, 0.5], [0.0, 0.0]),
            OID(FIndex{:f}('a', 1, 1, 1, 1), [0.5, 0.5], [0.0, 0.0])
            ))
    )
    @test term|>abbr == :cl
    @test expand(term, bond, config, true) == operators
    @test expand(term, bond, config, false) == operators+operators'
end
