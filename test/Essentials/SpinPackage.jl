using Test
using StaticArrays: SVector
using QuantumLattices.Essentials.SpinPackage
using QuantumLattices.Essentials.SpinPackage: heisenbergxyz, heisenbergpmz, gamma, dm
using QuantumLattices.Essentials.Spatials: AbstractPID, PID, CPID, Point, Bond
using QuantumLattices.Essentials.DegreesOfFreedom: OID, AbstractCompositeOID, isHermitian, Hilbert, Index, Operator, Operators, script, latexname
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
    hilbert = Hilbert(pid=>Spin{1}(norbital=2) for pid in [bond.epoint.pid, bond.spoint.pid])
    ex = expand(sc, bond, hilbert, Val(:SpinTerm))
    @test Dims(ex) == (1,)
    @test eltype(ex) == Tuple{Float, ID{OID{Index{CPID{Int}, SID{1}}, SVector{1, Float}}, 2}}
    @test collect(ex) == [(2.0, ID(
        OID(Index(CPID(1, 1), SID{1}(1, '+')), [0.0], [0.0]),
        OID(Index(CPID(1, 2), SID{1}(2, '-')), [0.5], [0.0])
        ))]

    sc = SpinCoupling(2.0, ('+', '-', '+', '-'), orbitals=subscripts"[α α β β](α < β)")
    point = Point(PID(1), [0.0], [0.0])
    hilbert = Hilbert(point.pid=>Spin{1}(norbital=3))
    ex = expand(sc, point, hilbert, Val(:info))
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

    @test ising"x ob[α β](α < β)" == Couplings(SpinCoupling(1, ('x', 'x'), orbitals=subscripts"[α β](α < β)"))
    @test ising"y ob[α β](α < β)" == Couplings(SpinCoupling(1, ('y', 'y'), orbitals=subscripts"[α β](α < β)"))
    @test ising"z ob[α β](α < β)" == Couplings(SpinCoupling(1, ('z', 'z'), orbitals=subscripts"[α β](α < β)"))
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
    @test term|>abbr == :sp
    @test term|>isHermitian == true
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
