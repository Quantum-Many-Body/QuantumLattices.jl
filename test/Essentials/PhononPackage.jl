using Test
using QuantumLattices.Essentials.PhononPackage
using QuantumLattices.Essentials.Spatials: PID, AbstractPID, Point, Bond
using QuantumLattices.Essentials.DegreesOfFreedom: Index, script, OID, AbstractCompositeOID, latexname, Operator, isHermitian, Hilbert, Operators
using QuantumLattices.Essentials.Terms: IIDSpace, Couplings, abbr
using QuantumLattices.Interfaces: permute, expand
using QuantumLattices.Mathematics.VectorSpaces: shape, ndimshape
using QuantumLattices.Mathematics.AlgebraOverFields: ID

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
