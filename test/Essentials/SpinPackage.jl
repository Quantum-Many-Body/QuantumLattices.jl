using Test
using StaticArrays: SVector
using QuantumLattices.Essentials.SpinPackage
using QuantumLattices.Essentials.Spatials: PID, Point, Bond
using QuantumLattices.Essentials.DegreesOfFreedom: Table, OID, isHermitian, Config, Operators, oidtype, script, latexname, iid
using QuantumLattices.Essentials.Terms: Couplings, @subscripts_str, abbr, otype
using QuantumLattices.Interfaces: expand, permute, rank
using QuantumLattices.Prerequisites: Float
using QuantumLattices.Mathematics.Combinatorics: Permutations
using QuantumLattices.Mathematics.AlgebraOverFields: ID

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
    spin = Spin{1}(atom=1, norbital=2)
    @test Dims(spin) == (2, 5)
    @test CartesianIndex(SID{1}(1, 'z'), spin) == CartesianIndex(1, 3)
    @test SID(CartesianIndex(1, 1), spin) == SID{1}(1, 'x')
    @test summary(spin) == "10-element Spin{1}"
    @test totalspin(spin) == totalspin(typeof(spin)) == 1
    @test collect(spin) == [
        SID{1}(1, 'x'), SID{1}(2, 'x'), SID{1}(1, 'y'), SID{1}(2, 'y'), SID{1}(1, 'z'), SID{1}(2, 'z'),
        SID{1}(1, '+'), SID{1}(2, '+'), SID{1}(1, '-'), SID{1}(2, '-')
    ]
end

@testset "SIndex" begin
    index = SIndex{1//2}(1, 1, 1, 'z')
    @test string(index) == "SIndex{1//2}(1, 1, 1, 'z')"
    @test replace(index, tag='x') == SIndex{1//2}(1, 1, 1, 'x')
    @test totalspin(index) == totalspin(typeof(index)) == 1//2
    @test union(PID{Char}, SID{1//2}) == SIndex{1//2, Char}
end

@testset "oidtype" begin
    @test oidtype(Spin{1//2}, Point{2, PID{Int}}, Val(:info)) == OID{SIndex{1//2, Int}, SVector{2, Float}}
end

@testset "latex" begin
    index = SIndex{1//2}('S', 1, 2, 'z')
    @test script(Val(:site), index) == 1
    @test script(Val(:orbital), index) == 2
    @test script(Val(:tag), index) == 'z'

    @test latexname(SIndex) == Symbol("SIndex")
    @test latexname(OID{<:SIndex}) == Symbol("OID{SIndex}")
end

@testset "SOperator" begin
    opt = SOperator(1.0, ID(
        OID(SIndex{1//2}('a', 1, 1, '+'), [0.0, 0.0], [0.0, 0.0]),
        OID(SIndex{1//2}('a', 1, 1, '-'), [0.0, 0.0], [0.0, 0.0])
        ))
    @test opt' == SOperator(1.0, ID(
        OID(SIndex{1//2}('a', 1, 1, '+'), [0.0, 0.0], [0.0, 0.0]),
        OID(SIndex{1//2}('a', 1, 1, '-'), [0.0, 0.0], [0.0, 0.0])
        ))
    @test isHermitian(opt)
    @test repr(opt) == "S^{+}_{1, 1}S^{-}_{1, 1}"
end

@testset "permute" begin
    soptrep(opt::SOperator) = opt.value * prod([Matrix(opt.id[i].index|>iid) for i = 1:rank(opt)])
    for S in (1//2, 1, 3//2)
        oids = [OID(SIndex{S}('S', 1, 2, tag), [0.0, 0.0], [0.0, 0.0]) for tag in ('x', 'y', 'z', '+', '-')]
        for (id₁, id₂) in Permutations{2}(oids)
            left = soptrep(SOperator(1, ID(id₁, id₂)))
            right = sum([soptrep(opt) for opt in permute(SOperator, id₁, id₂)])
            @test isapprox(left, right)
        end
    end
    id₁ = OID(SIndex{1//2}('S', 1, 2, 'z'), [0.0, 0.0], [0.0, 0.0])
    id₂ = OID(SIndex{1//2}('S', 2, 2, 'z'), [0.0, 0.0], [0.0, 0.0])
    @test permute(SOperator, id₁, id₂) == (SOperator(1, ID(id₂, id₁)),)
end

@testset "SpinCoupling" begin
    @test string(SpinCoupling(1.0, ('+', '-'))) == "SpinCoupling(value=1.0, tags=S⁺S⁻)"
    @test string(SpinCoupling(1.0, ('z', 'z'), atoms=(1, 1))) == "SpinCoupling(value=1.0, tags=SᶻSᶻ, atoms=[1 1])"
    @test string(SpinCoupling(1.0, ('-', '+'), atoms=(1, 1), orbitals=(1, 2))) == "SpinCoupling(value=1.0, tags=S⁻S⁺, atoms=[1 1], orbitals=[1 2])"
    @test repr(SpinCoupling(2.0, ('x', 'y'))) == "2.0 SˣSʸ"

    sc₁ = SpinCoupling(1.5, ('+', '-'), atoms=(1, 2), orbitals=subscripts"[x y](x > y)")
    sc₂ = SpinCoupling(2.0, ('+', '-'), atoms=(1, 2), orbitals=subscripts"[x y](x < y)")
    @test repr(sc₁) == "1.5 S⁺S⁻ sl[1 2] ⊗ ob[x y](x > y)"
    @test repr(sc₂) == "2.0 S⁺S⁻ sl[1 2] ⊗ ob[x y](x < y)"

    sc = sc₁ * sc₂
    @test repr(sc) == "3.0 S⁺S⁻S⁺S⁻ sl[1 2 1 2] ⊗ ob[x y; x y](x > y, x < y)"

    sc = SpinCoupling(2.0, ('+', '-'), atoms=(1, 1))
    point = Point(PID(1, 1), [0.0, 0.0], [0.0, 0.0])
    spin = Spin{1}(atom=2, norbital=2)
    ex = expand(sc, (point, point), (spin, spin), Val(:info))
    @test collect(ex) == []

    sc = SpinCoupling(2.0, ('+', '-'), atoms=(1, 2), orbitals=(1, 2))
    p₁, p₂ = Point(PID(1, 1), [0.0], [0.0]), Point(PID(1, 2), [0.5], [0.0])
    s₁, s₂ = Spin{1}(atom=1, norbital=2), Spin{1}(atom=2, norbital=2)
    ex = expand(sc, (p₁, p₂), (s₁, s₂), Val(:info))
    @test Dims(ex) == (1,)
    @test eltype(ex) == Tuple{Float, ID{OID{SIndex{1, Int}, SVector{1, Float}} , 2}}
    @test collect(ex) == [(2.0, ID(
        OID(SIndex{1}(1, 1, 1, '+'), [0.0], [0.0]),
        OID(SIndex{1}(1, 2, 2, '-'), [0.5], [0.0])
        ))]

    sc = SpinCoupling(2.0, ('+', '-', '+', '-'), orbitals=subscripts"[α α β β](α < β)")
    point = Point(PID(1, 1), [0.0], [0.0])
    spin = Spin{1}(norbital=3)
    ex = expand(sc, (point, point, point, point), (spin, spin, spin, spin), Val(:info))
    @test Dims(ex) == (3,)
    @test collect(ex) == [
        (2.0, ID(OID(SIndex{1}(1, 1, 1, '+'), [0.0], [0.0]),
                 OID(SIndex{1}(1, 1, 1, '-'), [0.0], [0.0]),
                 OID(SIndex{1}(1, 1, 2, '+'), [0.0], [0.0]),
                 OID(SIndex{1}(1, 1, 2, '-'), [0.0], [0.0])
                 )),
        (2.0, ID(OID(SIndex{1}(1, 1, 1, '+'), [0.0], [0.0]),
                 OID(SIndex{1}(1, 1, 1, '-'), [0.0], [0.0]),
                 OID(SIndex{1}(1, 1, 3, '+'), [0.0], [0.0]),
                 OID(SIndex{1}(1, 1, 3, '-'), [0.0], [0.0])
                 )),
        (2.0, ID(OID(SIndex{1}(1, 1, 2, '+'), [0.0], [0.0]),
                 OID(SIndex{1}(1, 1, 2, '-'), [0.0], [0.0]),
                 OID(SIndex{1}(1, 1, 3, '+'), [0.0], [0.0]),
                 OID(SIndex{1}(1, 1, 3, '-'), [0.0], [0.0])
                 ))
    ]
end

@testset "Heisenberg" begin
    @test Heisenberg(orbitals=(1, 2)) == Couplings(
        SpinCoupling(1//1, ('z', 'z'), orbitals=(1, 2)),
        SpinCoupling(1//2, ('+', '-'), orbitals=(1, 2)),
        SpinCoupling(1//2, ('-', '+'), orbitals=(1, 2))
    )
    @test Heisenberg("xyz") == Couplings(SpinCoupling(1, ('x', 'x')), SpinCoupling(1, ('y', 'y')), SpinCoupling(1, ('z', 'z')))
end

@testset "Ising" begin
    @test Ising('x', atoms=(1, 2)) == Couplings(SpinCoupling(1, ('x', 'x'), atoms=(1, 2)))
    @test Ising('y', atoms=(1, 2)) == Couplings(SpinCoupling(1, ('y', 'y'), atoms=(1, 2)))
    @test Ising('z', atoms=(1, 2)) == Couplings(SpinCoupling(1, ('z', 'z'), atoms=(1, 2)))
end

@testset "Gamma" begin
    @test Gamma('x', orbitals=(1, 1)) == SpinCoupling(1, ('y', 'z'), orbitals=(1, 1)) + SpinCoupling(1, ('z', 'y'), orbitals=(1, 1))
    @test Gamma('y', atoms=(1, 2)) == SpinCoupling(1, ('z', 'x'), atoms=(1, 2)) + SpinCoupling(1, ('x', 'z'), atoms=(1, 2))
    @test Gamma('z') == SpinCoupling(1, ('x', 'y')) + SpinCoupling(1, ('y', 'x'))
end

@testset "DM" begin
    @test DM('x', orbitals=(1, 1)) == SpinCoupling(1, ('y', 'z'), orbitals=(1, 1)) - SpinCoupling(1, ('z', 'y'), orbitals=(1, 1))
    @test DM('y', atoms=(1, 2)) == SpinCoupling(1, ('z', 'x'), atoms=(1, 2)) - SpinCoupling(1, ('x', 'z'), atoms=(1, 2))
    @test DM('z') == SpinCoupling(1, ('x', 'y')) - SpinCoupling(1, ('y', 'x'))
end

@testset "Sᵅ" begin
    @test Sˣ(atom=1, orbital=1) == Couplings(SpinCoupling(1, ('x',), atoms=(1,), orbitals=(1,)))
    @test Sʸ(atom=1) == Couplings(SpinCoupling(1, ('y',), atoms=(1,)))
    @test Sᶻ(orbital=1) == Couplings(SpinCoupling(1, ('z',), orbitals=(1,)))
end

@testset "spincoupling" begin
    sc = sc"1.0 S⁺S⁻ sl[1 1] ⊗ ob[α β](α < β)"
    @test repr(sc) == "1.0 S⁺S⁻ sl[1 1] ⊗ ob[α β](α < β)"

    sc = sc"1.0 S⁺S⁻ sl[1 1] ⊗ ob[α β]"
    @test repr(sc) == "1.0 S⁺S⁻ sl[1 1] ⊗ ob[α β]"

    sc = sc"1.0 S⁺S⁻ sl[1 1] ⊗ ob[1 2]"
    @test repr(sc) == "1.0 S⁺S⁻ sl[1 1] ⊗ ob[1 2]"

    sc = sc"1.0 S⁺S⁻ sl[1 1]"
    @test repr(sc) == "1.0 S⁺S⁻ sl[1 1]"

    sc = sc"1.0 S⁺S⁻"
    @test repr(sc) == "1.0 S⁺S⁻"
end

@testset "spincouplings" begin
    @test heisenberg"xyz sl[1 1] ⊗ ob[α β](α < β)" == Heisenberg("xyz", atoms=(1, 1), orbitals=subscripts"[α β](α < β)")
    @test heisenberg"sl[1 1] ⊗ ob[1 3]" == Heisenberg(atoms=(1, 1), orbitals=(1, 3))
    @test heisenberg"ob[1 3]" == Heisenberg(orbitals=(1, 3))
    @test heisenberg"sl[1 1]" == Heisenberg(atoms=(1, 1))
    @test heisenberg"" == heisenberg"+-z" == Heisenberg()
    @test heisenberg"xyz" == Heisenberg("xyz")

    @test ising"x" == Ising('x') && ising"x sl[1 1] ⊗ ob[1 3]" == Ising('x', atoms=(1, 1), orbitals=(1, 3))
    @test ising"y" == Ising('y') && ising"y sl[1 1] ⊗ ob[1 3]" == Ising('y', atoms=(1, 1), orbitals=(1, 3))
    @test ising"z" == Ising('z') && ising"z sl[1 1] ⊗ ob[1 3]" == Ising('z', atoms=(1, 1), orbitals=(1, 3))

    @test gamma"x" == Gamma('x') && gamma"x sl[1 1] ⊗ ob[1 3]" == Gamma('x', atoms=(1, 1), orbitals=(1, 3))
    @test gamma"y" == Gamma('y') && gamma"y sl[1 1] ⊗ ob[1 3]" == Gamma('y', atoms=(1, 1), orbitals=(1, 3))
    @test gamma"z" == Gamma('z') && gamma"z sl[1 1] ⊗ ob[1 3]" == Gamma('z', atoms=(1, 1), orbitals=(1, 3))

    @test dm"x" == DM('x') && dm"x sl[1 1] ⊗ ob[1 3]" == DM('x', atoms=(1, 1), orbitals=(1, 3))
    @test dm"y" == DM('y') && dm"y sl[1 1] ⊗ ob[1 3]" == DM('y', atoms=(1, 1), orbitals=(1, 3))
    @test dm"z" == DM('z') && dm"z sl[1 1] ⊗ ob[1 3]" == DM('z', atoms=(1, 1), orbitals=(1, 3))

    @test sˣ"" == Sˣ() && sˣ"sl[1]⊗ob[2]" == Sˣ(atom=1, orbital=2)
    @test sʸ"" == Sʸ() && sʸ"sl[1]⊗ob[2]" == Sʸ(atom=1, orbital=2)
    @test sᶻ"" == Sᶻ() && sᶻ"sl[1]⊗ob[2]" == Sᶻ(atom=1, orbital=2)
end

@testset "SpinTerm" begin
    point = Point(PID('a', 1), (0.5, 0.5), (0.0, 0.0))
    config = Config{Spin{1//2}}(pid->Spin{1//2}(atom=pid.site%2, norbital=2), [point.pid])
    term = SpinTerm{1}(:h, 1.5, 0, couplings=Sᶻ())
    operators = Operators(
        SOperator(1.5, ID(OID(SIndex{1//2}('a', 1, 1, 'z'), [0.5, 0.5], [0.0, 0.0]))),
        SOperator(1.5, ID(OID(SIndex{1//2}('a', 1, 2, 'z'), [0.5, 0.5], [0.0, 0.0])))
    )
    @test term|>abbr == :sp
    @test term|>isHermitian == true
    @test expand(term, point, config) == operators

    bond = Bond(1, Point(PID('a', 1), (0.0, 0.0), (0.0, 0.0)), Point(PID('b', 1), (0.5, 0.5), (0.0, 0.0)))
    config = Config{Spin{1//2}}(pid->Spin{1//2}(atom=pid.site%2, norbital=2), [bond.spoint.pid, bond.epoint.pid])
    term = SpinTerm{2}(:J, 1.5, 1, couplings=Heisenberg())
    operators = Operators(
        SOperator(1.50, ID(OID(SIndex{1//2}('b', 1, 2, 'z'), [0.5, 0.5], [0.0, 0.0]), OID(SIndex{1//2}('a', 1, 2, 'z'), [0.0, 0.0], [0.0, 0.0]))),
        SOperator(0.75, ID(OID(SIndex{1//2}('b', 1, 2, '-'), [0.5, 0.5], [0.0, 0.0]), OID(SIndex{1//2}('a', 1, 2, '+'), [0.0, 0.0], [0.0, 0.0]))),
        SOperator(0.75, ID(OID(SIndex{1//2}('b', 1, 1, '-'), [0.5, 0.5], [0.0, 0.0]), OID(SIndex{1//2}('a', 1, 1, '+'), [0.0, 0.0], [0.0, 0.0]))),
        SOperator(0.75, ID(OID(SIndex{1//2}('b', 1, 1, '+'), [0.5, 0.5], [0.0, 0.0]), OID(SIndex{1//2}('a', 1, 1, '-'), [0.0, 0.0], [0.0, 0.0]))),
        SOperator(1.50, ID(OID(SIndex{1//2}('b', 1, 1, 'z'), [0.5, 0.5], [0.0, 0.0]), OID(SIndex{1//2}('a', 1, 1, 'z'), [0.0, 0.0], [0.0, 0.0]))),
        SOperator(0.75, ID(OID(SIndex{1//2}('b', 1, 2, '+'), [0.5, 0.5], [0.0, 0.0]), OID(SIndex{1//2}('a', 1, 2, '-'), [0.0, 0.0], [0.0, 0.0])))
    )
    @test expand(term, bond, config) == operators
end
