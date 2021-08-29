using Test
using Printf: @sprintf
using StaticArrays: SVector
using QuantumLattices.Essentials.Terms
using QuantumLattices.Essentials.Spatials: AbstractBond, Point, PID, CPID, Bond, Bonds, Lattice, pidtype, acrossbonds, zerothbonds
using QuantumLattices.Essentials.DegreesOfFreedom: Hilbert, IID, Index, Internal, OID, Table, OIDToTuple, Operator, Operators, plain
using QuantumLattices.Interfaces: rank, expand!
using QuantumLattices.Essentials: kind, update!, reset!
using QuantumLattices.Prerequisites: Float, decimaltostr
using QuantumLattices.Prerequisites.Traits: getcontent, contentnames
using QuantumLattices.Prerequisites.CompositeStructures: NamedContainer
using QuantumLattices.Mathematics.AlgebraOverFields: ID, SimpleID, id, idtype
import QuantumLattices.Interfaces: dimension, expand
import QuantumLattices.Essentials.DegreesOfFreedom: isHermitian
import QuantumLattices.Essentials.Terms: couplingcenters, abbr

struct TID <: IID nambu::Int end
Base.adjoint(sl::TID) = TID(3-sl.nambu)

struct TFock <: Internal{TID} end
Base.Dims(vs::TFock) = (2,)
TID(i::CartesianIndex, vs::TFock) = TID(i.I...)
Base.CartesianIndex(tid::TID, vs::TFock) = CartesianIndex(tid.nambu)

struct TCID <: SimpleID
    nambu::Int
end
Base.fieldnames(::Type{<:TCID}) = (:center, :nambu)

struct TCoupling{V<:Number, I<:ID{TCID}} <: Coupling{V, I}
    value::V
    id::I
end
@inline Base.repr(tc::TCoupling) = @sprintf "%s ph(%s)" decimaltostr(tc.value) join(tc.id.nambus, ", ")
@inline couplingcenters(::TCoupling, ::Bond, ::Val) = (1, 2)
function expand(tc::TCoupling, bond::AbstractBond, hilbert::Hilbert, info::Val)
    points = couplingpoints(tc, bond, info)
    @assert rank(tc)==length(points)
    nambus = tc.id.nambus
    pids = NTuple{rank(tc), pidtype(points|>eltype)}(points[i].pid for i = 1:rank(tc))
    rcoords = NTuple{rank(tc), SVector{dimension(points|>eltype), Float}}(points[i].rcoord for i = 1:rank(tc))
    icoords = NTuple{rank(tc), SVector{dimension(points|>eltype), Float}}(points[i].icoord for i = 1:rank(tc))
    indexes = NTuple{rank(tc), Index{pidtype(points|>eltype), TID}}(Index(pids[i], TID(nambus[i])) for i = 1:rank(tc))
    return ((tc.value, ID(OID, indexes, rcoords, icoords)),)
end

abbr(::Type{<:Term{:Mu}}) = :mu
isHermitian(::Type{<:Term{:Mu}}) = true

abbr(::Type{<:Term{:Hp}}) = :hp
isHermitian(::Type{<:Term{:Hp}}) = false

@testset "Subscripts" begin
    sub = Subscripts(4)
    @test rank(sub) == 1
    @test dimension(sub) == 4
    @test string(sub) == "[* * * *]"
    @test sub((2,)) == (2, 2, 2, 2)
    @test isvalid(sub, (2,)) == true

    sub = Subscripts((1, 2, 2, 1))
    @test rank(sub) == 0
    @test dimension(sub) == 4
    @test string(sub) == "[1 2 2 1]"
    @test sub(()) == (1, 2, 2, 1)
    @test isvalid(sub, ()) == true

    sub = subscripts"[x₁ x₂ x₁ x₂]"
    @test rank(sub) == 2
    @test dimension(sub) == 4
    @test string(sub) == "[x₁ x₂ x₁ x₂]"
    @test sub((1, 2)) == (1, 2, 1, 2)
    @test isvalid(sub, (1, 2)) == true

    sub = subscripts"[x₁ 4 4 x₂](x₁ < x₂)"
    @test rank(sub) == rank(typeof(sub)) == 2
    @test dimension(sub) == dimension(typeof(sub)) == 4
    @test string(sub) == "[x₁ 4 4 x₂](x₁ < x₂)"
    @test sub((1, 2)) == (1, 4, 4, 2)
    @test isvalid(sub, (1, 2)) == true
    @test isvalid(sub, (2, 1)) == false

    sub₁ = @subscripts [1 2 2 1]
    sub₂ = @subscripts [x₁ x₂ x₁ x₂]
    sub₃ = @subscripts [x₁ 4 4 x₂](x₁ < x₂)
    subs = sub₁*sub₂*sub₃
    segs = split(subs)
    @test contentnames(typeof(subs)) == (:contents, :cpattern, :mapping, :constrain)
    @test getcontent(subs, :contents) == subs.opattern
    @test subs == subscripts"[1 2 2 1; x₁ x₂ x₁ x₂; x₁ 4 4 x₂](:, :, x₁ < x₂)" == @subscripts [1 2 2 1; x₁ x₂ x₁ x₂; x₁ 4 4 x₂](:, :, x₁ < x₂)
    @test isequal(deepcopy(subs), subs)
    @test segs == (sub₁, sub₂, sub₃)
    @test string(subs) == "[1 2 2 1; x₁ x₂ x₁ x₂; x₁ 4 4 x₂](:, :, x₁ < x₂)"
    @test expand(subs, (1, 2, 3, 4, 1, 2, 1, 2, 2, 5, 5, 3))|>collect == [
        (1, 2, 2, 1, 1, 1, 1, 1, 1, 4, 4, 2), (1, 2, 2, 1, 1, 2, 1, 2, 1, 4, 4, 2),
        (1, 2, 2, 1, 1, 1, 1, 1, 1, 4, 4, 3), (1, 2, 2, 1, 1, 2, 1, 2, 1, 4, 4, 3),
        (1, 2, 2, 1, 1, 1, 1, 1, 2, 4, 4, 3), (1, 2, 2, 1, 1, 2, 1, 2, 2, 4, 4, 3)
    ]

    sub = segs[1]
    @test sub == subscripts"[1 2 2 1]"
    @test string(sub) == "[1 2 2 1]"
    @test rank(sub) == rank(subs, 1)
    @test dimension(sub) == dimension(subs, 1)
    @test sub(()) == (1, 2, 2, 1)
    @test isvalid(sub, ()) == true

    sub = segs[2]
    @test sub == subscripts"[x₁ x₂ x₁ x₂]"
    @test string(sub) == "[x₁ x₂ x₁ x₂]"
    @test rank(sub) == rank(subs, 2)
    @test dimension(sub) == dimension(subs, 2)
    @test sub((1, 2)) == (1, 2, 1, 2)
    @test isvalid(sub, (1, 2)) == true

    sub = segs[3]
    @test sub == subscripts"[x₁ 4 4 x₂](x₁ < x₂)"
    @test string(sub) == "[x₁ 4 4 x₂](x₁ < x₂)"
    @test rank(sub) == rank(subs, 3)
    @test dimension(sub) == dimension(subs, 3)
    @test sub((1, 2)) == (1, 4, 4, 2)
    @test isvalid(sub, (1, 2)) == true
    @test isvalid(sub, (2, 1)) == false

    subid = SubID((1, 2, 2, 1, :x₁, :x₂, :x₁, :x₂, :x₁, 4, 4, :x₂), ("nonconstrain", "nonconstrain", "x₁ < x₂"), (4, 4, 4))
    @test subid == SubID(subs) == id(subs)
end

@testset "couplings" begin
    @test @couplings(TCoupling(1.0, ID(TCID(1), TCID(1)))) == Couplings(TCoupling(1.0, ID(TCID(1), TCID(1))))
    @test @couplings(TCoupling(1.0, ID(TCID(1)))+TCoupling(1.0, ID(TCID(2)))) == TCoupling(1.0, ID(TCID(1)))+TCoupling(1.0, ID(TCID(2)))

    point = Point(CPID(1, 1), (0.0, 0.0), (0.0, 0.0))
    hilbert = Hilbert{TFock}(pid->TFock(), [CPID(1, 1)])
    tc = TCoupling(1.0, ID(TCID(1), TCID(1)))
    @test couplingcenters(tc, point, Val(:Mu)) == (1, 1)
    @test couplingpoints(tc, point, Val(:Mu)) == (point, point)
    @test couplinginternals(tc, point, hilbert, Val(:Mu)) == (TFock(), TFock())
end

@testset "TermFunction" begin
    ta = TermAmplitude()
    @test ta() == 1

    ta = TermAmplitude(x->x+3.0)
    @test ta(1.0) == 4.0

    tcs = TCoupling(1.0, ID(TCID(1), TCID(1))) + TCoupling(2.0, ID(TCID(2), TCID(2)))
    termcouplings = TermCouplings(tcs)
    @test termcouplings == deepcopy(TermCouplings(tcs))
    @test isequal(termcouplings, deepcopy(TermCouplings(tcs)))
    @test termcouplings() == tcs

    tcs1 = TCoupling(1.0, ID(TCID(1), TCID(1))) + TCoupling(1.0, ID(TCID(2), TCID(2)))
    tcs2 = TCoupling(1.0, ID(TCID(2), TCID(1))) + TCoupling(1.0, ID(TCID(1), TCID(2)))
    termcouplings = TermCouplings(i->(tcs1, tcs2)[(i-1)%2+1])
    @test termcouplings(1) == tcs1
    @test termcouplings(2) == tcs2

    @test ismodulatable(TermModulate{Val{true}}) == ismodulatable(TermModulate{<:Function}) == true
    @test ismodulatable(TermModulate{Val{false}}) == false
    termmodulate = TermModulate(:t)
    @test termmodulate(t=1) == 1
    @test termmodulate(mu=1) == nothing
    @test termmodulate() == nothing

    termmodulate = TermModulate(:t, t->t*2.0)
    @test termmodulate(2) == 4
end

@testset "Term" begin
    point = Point(CPID(1, 1), (0.0, 0.0), (0.0, 0.0))
    hilbert = Hilbert{TFock}(pid->TFock(), [CPID(1, 1)])
    term = Term{:Mu, 2}(:mu, 1.5, 0, couplings=TCoupling(1.0, ID(TCID(2), TCID(1))), amplitude=(bond->3), modulate=false)
    @test term|>kind == term|>typeof|>kind == :Mu
    @test term|>id == term|>typeof|>id == :mu
    @test term|>valtype == term|>typeof|>valtype == Float
    @test term|>rank == term|>typeof|>rank == 2
    @test term|>abbr == term|>typeof|>abbr == :mu
    @test term|>ismodulatable == term|>typeof|>ismodulatable == false
    @test term|>isHermitian == term|>typeof|>isHermitian == true
    @test term == deepcopy(term)
    @test isequal(term, deepcopy(term))
    @test string(term) == "Mu{2}(id=mu, value=1.5, bondkind=0, factor=1.0)"
    @test repr(term, point, hilbert) == "mu: 4.5 ph(2, 1)"
    @test +term == term
    @test -term == term*(-1) == replace(term, factor=-term.factor)
    @test 2*term == term*2 == replace(term, factor=term.factor*2)
    @test term/2 == term*(1/2) == replace(term, factor=term.factor/2)

    p1 = Point(PID(1), (0.0, 0.0), (0.0, 0.0))
    p2 = Point(PID(2), (0.0, 0.0), (0.0, 0.0))
    hilbert = Hilbert{TFock}(pid->TFock(), [PID(1), PID(2)])
    tcs1 = Couplings(TCoupling(1.0, ID(TCID(2), TCID(2))))
    tcs2 = Couplings(TCoupling(1.0, ID(TCID(1), TCID(1))))
    term = Term{:Mu, 2}(:mu, 1.5, 0, couplings=bond->(tcs1, tcs2)[bond.pid.site%2+1], amplitude=bond->3, modulate=true)
    @test term|>ismodulatable == term|>typeof|>ismodulatable == true
    @test repr(term, p1, hilbert) == "mu: 4.5 ph(1, 1)"
    @test repr(term, p2, hilbert) == "mu: 4.5 ph(2, 2)"
    @test one(term) == replace(term, value=1.0)
    @test zero(term) == replace(term, value=0.0)
    @test term.modulate(mu=4.0) == 4.0
    @test term.modulate(t=1.0) == nothing
    @test update!(term, mu=4.25) == replace(term, value=4.25)
    @test term.value == 4.25
end

@testset "expand" begin
    point = Point(PID(1), (0.0, 0.0), (0.0, 0.0))
    hilbert = Hilbert{TFock}(pid->TFock(), [PID(1)])
    term = Term{:Mu, 2}(:mu, 1.5, 0, couplings=TCoupling(1.0, ID(TCID(2), TCID(1))), amplitude=bond->3.0, modulate=true)
    operators = Operators(
        Operator(+2.25, ID(
            OID(Index(PID(1), TID(2)), SVector(0.0, 0.0), SVector(0.0, 0.0)),
            OID(Index(PID(1), TID(1)), SVector(0.0, 0.0), SVector(0.0, 0.0))
            ))
        )
    @test expand(term, point, hilbert, true) == operators
    @test expand(term, point, hilbert, false) == operators*2

    bond = Bond(1, Point(PID(2), (1.5, 1.5), (1.0, 1.0)), Point(PID(1), (0.5, 0.5), (0.0, 0.0)))
    hilbert = Hilbert{TFock}(pid->TFock(), [PID(1), PID(2)])
    term = Term{:Hp, 2}(:t, 1.5, 1, couplings=TCoupling(1.0, ID(TCID(2), TCID(1))), amplitude=bond->3.0, modulate=true)
    operators = Operators(
        Operator(4.5, ID(
            OID(Index(PID(1), TID(2)), SVector(0.5, 0.5), SVector(0.0, 0.0)),
            OID(Index(PID(2), TID(1)), SVector(1.5, 1.5), SVector(1.0, 1.0))
            ))
        )
    @test expand(term, bond, hilbert, true) == operators
    @test expand(term, bond, hilbert, false) == operators+operators'

    lattice = Lattice("Tuanzi", [Point(CPID(1, 1), (0.0, 0.0), (0.0, 0.0))], vectors=[[1.0, 0.0]], neighbors=1)
    bonds = Bonds(lattice)
    hilbert = Hilbert{TFock}(pid->TFock(), lattice.pids)
    term = Term{:Mu, 2}(:mu, 1.5, 0, couplings=TCoupling(1.0, ID(TCID(2), TCID(1))), amplitude=bond->3.0, modulate=true)
    operators = Operators(
        Operator(+2.25, ID(
            OID(Index(CPID(1, 1), TID(2)), SVector(0.0, 0.0), SVector(0.0, 0.0)),
            OID(Index(CPID(1, 1), TID(1)), SVector(0.0, 0.0), SVector(0.0, 0.0))
            ))
        )
    @test expand(term, bonds, hilbert, true) == operators
    @test expand(term, bonds, hilbert, false) == operators*2
end

@testset "Parameters" begin
    ps1 = Parameters{(:t1, :t2, :U)}(1.0im, 1.0, 2.0)
    ps2 = Parameters{(:t1, :U)}(1.0im, 2.0)
    ps3 = Parameters{(:t1, :U)}(1.0im, 2.1)
    @test match(ps1, ps2) == true
    @test match(ps1, ps3) == false
end

@testset "Generator" begin
    @test contentnames(AbstractGenerator) == (:terms, :bonds, :hilbert, :half, :table, :boundary, :operators)

    lattice = Lattice("Tuanzi", [Point(PID(1), (0.0, 0.0), (0.0, 0.0)), Point(PID(2), (0.5, 0.0), (0.0, 0.0))], vectors=[[1.0, 0.0]], neighbors=1)
    bonds = Bonds(lattice)
    hilbert = Hilbert{TFock}(pid->TFock(), lattice.pids)
    table = Table(hilbert, OIDToTuple(:scope, :site))
    t = Term{:Hp, 2}(:t, 2.0, 1, couplings=TCoupling(1.0, ID(TCID(2), TCID(1))))
    μ = Term{:Mu, 2}(:μ, 1.0, 0, couplings=TCoupling(1.0, ID(TCID(2), TCID(1))), modulate=true)
    tops1 = expand(t, filter(acrossbonds, bonds, Val(:exclude)), hilbert, true, table=table)
    tops2 = expand(one(t), filter(acrossbonds, bonds, Val(:include)), hilbert, true, table=table)
    μops = expand(one(μ), filter(zerothbonds, bonds, Val(:include)), hilbert, true, table=table)

    optp = Operator{Float, ID{OID{Index{PID, TID}, SVector{2, Float}}, 2}}
    genops = GenOperators(tops1, NamedContainer{(:μ,)}((μops,)), NamedContainer{(:t, :μ)}((tops2, Operators{optp|>idtype, optp}())))
    @test genops == deepcopy(genops) && isequal(genops, deepcopy(genops))
    @test genops == GenOperators((t, μ), bonds, hilbert, true, table=table)
    @test genops|>eltype == genops|>typeof|>eltype == optp
    @test genops|>idtype == genops|>typeof|>idtype == optp|>idtype
    @test expand!(Operators{idtype(optp), optp}(), genops, plain, t=2.0, μ=1.5) == tops1+tops2*2.0+μops*1.5
    @test empty(genops) == empty!(deepcopy(genops))
    @test empty(genops) == GenOperators(empty(μops), NamedContainer{(:μ,)}((empty(μops),)), NamedContainer{(:t, :μ)}((empty(μops), empty(μops))))
    @test reset!(deepcopy(genops), (t, μ), bonds, hilbert, true, table=table) == genops

    gen = Generator((t, μ), bonds, hilbert; half=true, table=table)
    @test gen == deepcopy(gen) && isequal(gen, deepcopy(gen))
    @test Parameters(gen) == Parameters{(:t, :μ)}(2.0, 1.0)
    @test expand!(Operators{idtype(optp), optp}(), gen) == expand(gen) == tops1+tops2*2.0+μops
    @test expand(gen, :t) == tops1+tops2*2.0
    @test expand(gen, :μ) == μops
    @test expand(gen, 1)+expand(gen, 2)+expand(gen, 3)+expand(gen, 4) == expand(gen)
    @test expand(gen, :μ, 1)+expand(gen, :μ, 2) == μops
    @test expand(gen, :t, 3) == tops1
    @test expand(gen, :t, 4) == tops2*2.0
    @test empty!(deepcopy(gen)) == Generator((t, μ), empty(bonds), empty(hilbert), half=true, table=empty(table), boundary=plain) == empty(gen)
    @test reset!(empty(gen), lattice) == gen
    @test update!(gen, μ=1.5)|>expand == tops1+tops2*2.0+μops*1.5
end
