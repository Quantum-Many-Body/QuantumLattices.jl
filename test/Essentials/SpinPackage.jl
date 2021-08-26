using Test
using StaticArrays: SVector
using QuantumLattices.Essentials.SpinPackage
using QuantumLattices.Essentials.Spatials: AbstractPID, PID, CPID, Point, Bond
using QuantumLattices.Essentials.DegreesOfFreedom: OID, AbstractCompositeOID, isHermitian, Config, Index, Operator, Operators, script, latexname
using QuantumLattices.Essentials.Terms: Couplings, @subscripts_str, abbr
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
    spin = Spin{1}(norbital=2)
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

@testset "latex" begin
    index = Index(PID(1), SID{1//2}(2, 'z'))
    @test script(Val(:site), index) == 1
    @test script(Val(:orbital), index) == 2
    @test script(Val(:tag), index) == 'z'

    @test latexname(Index{<:AbstractPID, <:SID}) == Symbol("Index{AbstractPID, SID}")
    @test latexname(AbstractCompositeOID{<:Index{<:AbstractPID, <:SID}}) == Symbol("AbstractCompositeOID{Index{AbstractPID, SID}}")
end

@testset "SOperator" begin
    opt = Operator(1.0, ID(
        OID(Index(PID(1), SID{1//2}(1, '+')), [0.0, 0.0], [0.0, 0.0]),
        OID(Index(PID(1), SID{1//2}(1, '-')), [0.0, 0.0], [0.0, 0.0])
        ))
    @test opt' == Operator(1.0, ID(
        OID(Index(PID(1), SID{1//2}(1, '+')), [0.0, 0.0], [0.0, 0.0]),
        OID(Index(PID(1), SID{1//2}(1, '-')), [0.0, 0.0], [0.0, 0.0])
        ))
    @test isHermitian(opt)
    @test repr(opt) == "S^{+}_{1, 1}S^{-}_{1, 1}"
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

    sc₁ = SpinCoupling(1.5, ('+', '-'), orbitals=subscripts"[x y](x > y)")
    sc₂ = SpinCoupling(2.0, ('+', '-'), orbitals=subscripts"[x y](x < y)")
    @test repr(sc₁) == "1.5 S⁺S⁻ ob[x y](x > y)"
    @test repr(sc₂) == "2.0 S⁺S⁻ ob[x y](x < y)"

    sc = sc₁ * sc₂
    @test repr(sc) == "3.0 S⁺S⁻S⁺S⁻ ob[x y; x y](x > y, x < y)"

    sc = SpinCoupling(2.0, ('+', '-'), orbitals=(1, 2))
    bond = Bond(1, Point(CPID(1, 2), [0.5], [0.0]), Point(CPID(1, 1), [0.0], [0.0]))
    config = Config{Spin{1}}(pid->Spin{1}(norbital=2), [bond.epoint.pid, bond.spoint.pid])
    ex = expand(sc, bond, config, Val(:SpinTerm))
    @test Dims(ex) == (1,)
    @test eltype(ex) == Tuple{Float, ID{OID{Index{CPID{Int}, SID{1}}, SVector{1, Float}}, 2}}
    @test collect(ex) == [(2.0, ID(
        OID(Index(CPID(1, 1), SID{1}(1, '+')), [0.0], [0.0]),
        OID(Index(CPID(1, 2), SID{1}(2, '-')), [0.5], [0.0])
        ))]

    sc = SpinCoupling(2.0, ('+', '-', '+', '-'), orbitals=subscripts"[α α β β](α < β)")
    point = Point(PID(1), [0.0], [0.0])
    config = Config{Spin{1}}(pid->Spin{1}(norbital=3), [point.pid])
    ex = expand(sc, point, config, Val(:info))
    @test eltype(ex) == Tuple{Float, ID{OID{Index{PID, SID{1}}, SVector{1, Float}}, 4}}
    @test Dims(ex) == (3,)
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

@testset "Heisenberg" begin
    @test Heisenberg(orbitals=(1, 2)) == Couplings(
        SpinCoupling(1//1, ('z', 'z'), orbitals=(1, 2)),
        SpinCoupling(1//2, ('+', '-'), orbitals=(1, 2)),
        SpinCoupling(1//2, ('-', '+'), orbitals=(1, 2))
    )
    @test Heisenberg("xyz") == Couplings(SpinCoupling(1, ('x', 'x')), SpinCoupling(1, ('y', 'y')), SpinCoupling(1, ('z', 'z')))
end

@testset "Ising" begin
    @test Ising('x', orbitals=(1, 1)) == Couplings(SpinCoupling(1, ('x', 'x'), orbitals=(1, 1)))
    @test Ising('y') == Couplings(SpinCoupling(1, ('y', 'y')))
    @test Ising('z') == Couplings(SpinCoupling(1, ('z', 'z')))
end

@testset "Gamma" begin
    @test Gamma('x', orbitals=(1, 1)) == SpinCoupling(1, ('y', 'z'), orbitals=(1, 1)) + SpinCoupling(1, ('z', 'y'), orbitals=(1, 1))
    @test Gamma('y') == SpinCoupling(1, ('z', 'x')) + SpinCoupling(1, ('x', 'z'))
    @test Gamma('z') == SpinCoupling(1, ('x', 'y')) + SpinCoupling(1, ('y', 'x'))
end

@testset "DM" begin
    @test DM('x', orbitals=(1, 1)) == SpinCoupling(1, ('y', 'z'), orbitals=(1, 1)) - SpinCoupling(1, ('z', 'y'), orbitals=(1, 1))
    @test DM('y') == SpinCoupling(1, ('z', 'x')) - SpinCoupling(1, ('x', 'z'))
    @test DM('z') == SpinCoupling(1, ('x', 'y')) - SpinCoupling(1, ('y', 'x'))
end

@testset "Sᵅ" begin
    @test Sˣ(orbital=1) == Couplings(SpinCoupling(1, ('x',), orbitals=(1,)))
    @test Sʸ() == Couplings(SpinCoupling(1, ('y',)))
    @test Sᶻ() == Couplings(SpinCoupling(1, ('z',)))
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

@testset "spincouplings" begin
    @test heisenberg"xyz ob[α β](α < β)" == Heisenberg("xyz", orbitals=subscripts"[α β](α < β)")
    @test heisenberg"ob[1 3]" == Heisenberg(orbitals=(1, 3))
    @test heisenberg"" == heisenberg"+-z" == Heisenberg()
    @test heisenberg"xyz" == Heisenberg("xyz")

    @test ising"x" == Ising('x') && ising"x ob[1 3]" == Ising('x', orbitals=(1, 3))
    @test ising"y" == Ising('y') && ising"y ob[1 3]" == Ising('y', orbitals=(1, 3))
    @test ising"z" == Ising('z') && ising"z ob[1 3]" == Ising('z', orbitals=(1, 3))

    @test gamma"x" == Gamma('x') && gamma"x ob[1 3]" == Gamma('x', orbitals=(1, 3))
    @test gamma"y" == Gamma('y') && gamma"y ob[1 3]" == Gamma('y', orbitals=(1, 3))
    @test gamma"z" == Gamma('z') && gamma"z ob[1 3]" == Gamma('z', orbitals=(1, 3))

    @test dm"x" == DM('x') && dm"x ob[1 3]" == DM('x', orbitals=(1, 3))
    @test dm"y" == DM('y') && dm"y ob[1 3]" == DM('y', orbitals=(1, 3))
    @test dm"z" == DM('z') && dm"z ob[1 3]" == DM('z', orbitals=(1, 3))

    @test sˣ"" == Sˣ() && sˣ"ob[2]" == Sˣ(orbital=2)
    @test sʸ"" == Sʸ() && sʸ"ob[2]" == Sʸ(orbital=2)
    @test sᶻ"" == Sᶻ() && sᶻ"ob[2]" == Sᶻ(orbital=2)
end

@testset "SpinTerm" begin
    point = Point(PID(1), (0.5, 0.5), (0.0, 0.0))
    config = Config{Spin{1//2}}(pid->Spin{1//2}(norbital=2), [point.pid])
    term = SpinTerm{1}(:h, 1.5, 0, couplings=Sᶻ())
    operators = Operators(
        Operator(1.5, ID(OID(Index(PID(1), SID{1//2}(1, 'z')), [0.5, 0.5], [0.0, 0.0]))),
        Operator(1.5, ID(OID(Index(PID(1), SID{1//2}(2, 'z')), [0.5, 0.5], [0.0, 0.0])))
    )
    @test term|>abbr == :sp
    @test term|>isHermitian == true
    @test expand(term, point, config) == operators

    bond = Bond(1, Point(CPID('a', 1), (0.0, 0.0), (0.0, 0.0)), Point(CPID('b', 1), (0.5, 0.5), (0.0, 0.0)))
    config = Config{Spin{1//2}}(pid->Spin{1//2}(norbital=2), [bond.spoint.pid, bond.epoint.pid])
    term = SpinTerm{2}(:J, 1.5, 1, couplings=Heisenberg())
    operators = Operators(
        Operator(1.50, ID(OID(Index(CPID('b', 1), SID{1//2}(2, 'z')), [0.5, 0.5], [0.0, 0.0]), OID(Index(CPID('a', 1), SID{1//2}(2, 'z')), [0.0, 0.0], [0.0, 0.0]))),
        Operator(0.75, ID(OID(Index(CPID('b', 1), SID{1//2}(2, '-')), [0.5, 0.5], [0.0, 0.0]), OID(Index(CPID('a', 1), SID{1//2}(2, '+')), [0.0, 0.0], [0.0, 0.0]))),
        Operator(0.75, ID(OID(Index(CPID('b', 1), SID{1//2}(1, '-')), [0.5, 0.5], [0.0, 0.0]), OID(Index(CPID('a', 1), SID{1//2}(1, '+')), [0.0, 0.0], [0.0, 0.0]))),
        Operator(0.75, ID(OID(Index(CPID('b', 1), SID{1//2}(1, '+')), [0.5, 0.5], [0.0, 0.0]), OID(Index(CPID('a', 1), SID{1//2}(1, '-')), [0.0, 0.0], [0.0, 0.0]))),
        Operator(1.50, ID(OID(Index(CPID('b', 1), SID{1//2}(1, 'z')), [0.5, 0.5], [0.0, 0.0]), OID(Index(CPID('a', 1), SID{1//2}(1, 'z')), [0.0, 0.0], [0.0, 0.0]))),
        Operator(0.75, ID(OID(Index(CPID('b', 1), SID{1//2}(2, '+')), [0.5, 0.5], [0.0, 0.0]), OID(Index(CPID('a', 1), SID{1//2}(2, '-')), [0.0, 0.0], [0.0, 0.0])))
    )
    @test expand(term, bond, config) == operators
end
