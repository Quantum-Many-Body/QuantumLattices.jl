using Test
using Printf: @sprintf
using StaticArrays: SVector
using QuantumLattices.Essentials.Terms
using QuantumLattices.Essentials.Spatials: AbstractBond, Point, PID, CPID, Bond, Bonds, Lattice, pidtype, acrossbonds, zerothbonds
using QuantumLattices.Essentials.DegreesOfFreedom: SimpleIID, SimpleInternal, CompositeIID, IIDSpace, IIDConstrain, ConstrainID, Subscript, @subscript_str
using QuantumLattices.Essentials.DegreesOfFreedom: Hilbert, Index, Table, OID, OIDToTuple, Operator, Operators, plain
using QuantumLattices.Interfaces: rank, expand!
using QuantumLattices.Essentials: kind, update!, reset!
using QuantumLattices.Prerequisites: Float, decimaltostr
using QuantumLattices.Prerequisites.Traits: parameternames, contentnames, getcontent
using QuantumLattices.Prerequisites.CompositeStructures: NamedContainer
using QuantumLattices.Mathematics.AlgebraOverFields: ID, SimpleID, id, idtype
import QuantumLattices.Interfaces: dimension, expand
import QuantumLattices.Mathematics.VectorSpaces: shape, ndimshape
import QuantumLattices.Essentials.DegreesOfFreedom: isHermitian
import QuantumLattices.Essentials.Terms: couplingcenters, abbr

struct TID{N<:Union{Int, Symbol}} <: SimpleIID
    nambu::N
end
Base.adjoint(sl::TID) = TID(3-sl.nambu)

struct TFock <: SimpleInternal{TID{Int}}
    nnambu::Int
end
ndimshape(::Type{TFock}) = 1
shape(vs::TFock) = (1:vs.nnambu,)
TID(i::CartesianIndex, vs::TFock) = TID(i.I...)
Base.CartesianIndex(tid::TID, vs::TFock) = CartesianIndex(tid.nambu)
@inline shape(iidspace::IIDSpace{TID{Symbol}, TFock}) = (1:iidspace.internal.nnambu,)
@inline shape(iidspace::IIDSpace{TID{Int}, TFock}) = (iidspace.iid.nambu:iidspace.iid.nambu,)

const TCoupling{V, I<:ID{TID}, C<:IIDConstrain, CI<:ConstrainID} = Coupling{V, I, C, CI}
@inline Base.repr(tc::(Coupling{V, <:ID{TID}} where V)) = @sprintf "%s ph(%s)" decimaltostr(tc.value) join(tc.cid.nambus, ", ")
@inline couplingcenters(::(Coupling{V, <:ID{TID}} where V), ::Bond, ::Val) = (1, 2)
@inline TCoupling(value, nambus::Tuple{Vararg{Int}}) = Coupling(value, ID(TID, nambus), IIDConstrain((nambu=Subscript(nambus),)))
@inline TCoupling(value, nambus::Subscript) = Coupling(value, ID(TID, convert(Tuple, nambus)), IIDConstrain((nambu=nambus,)))

abbr(::Type{<:Term{:Mu}}) = :mu
abbr(::Type{<:Term{:Hp}}) = :hp
isHermitian(::Type{<:Term{:Mu}}) = true
isHermitian(::Type{<:Term{:Hp}}) = false

@testset "couplings" begin
    @test parameternames(Coupling) == (:value, :cid, :constrain, :constrainid)
    @test contentnames(Coupling) == (:value, :id, :constrain)

    tc = TCoupling(2.0, (2,))
    @test id(tc) == ID(CompositeIID(tc.cid), ConstrainID(tc.constrain))
    @test rank(tc) == rank(typeof(tc)) == 1
    @test tc == Coupling(2.0, id(tc), tc.constrain)
    @test ID{SimpleIID}(tc) == tc.cid
    @test IIDConstrain(tc) == tc.constrain

    point = Point(CPID(1, 1), (0.0, 0.0), (0.0, 0.0))
    hilbert = Hilbert(point.pid=>TFock(2))
    @test couplingcenters(tc, point, Val(:Mu)) == (1,)
    @test couplingpoints(tc, point, Val(:Mu)) == (point,)
    @test couplinginternals(tc, point, hilbert, Val(:Mu)) == (TFock(2),)

    tc₁ = TCoupling(1.5, (1, 2))
    tc₂ = TCoupling(2.0, subscript"[a b](a < b)")
    @test tc₁*tc₂ == Coupling(3.0, ID(TID(1), TID(2), TID(:a), TID(:b)), IIDConstrain((nambu=Subscript((1, 2)),), (nambu=subscript"[a b](a < b)",)))

    ex = expand(tc₁, point, hilbert, Val(:info))
    @test eltype(ex) == eltype(typeof(ex)) == Tuple{Float64, NTuple{2, OID{Index{CPID{Int}, TID{Int}}, SVector{2, Float64}}}}
    @test collect(ex) == [
        (1.5, ID(OID(Index(CPID(1, 1), TID(1)), SVector(0.0, 0.0), SVector(0.0, 0.0)),
                 OID(Index(CPID(1, 1), TID(2)), SVector(0.0, 0.0), SVector(0.0, 0.0))
                 ))
        ]

    ex = expand(tc₂, point, hilbert, Val(:info))
    @test eltype(ex) == eltype(typeof(ex)) == Tuple{Float64, NTuple{2, OID{Index{CPID{Int}, TID{Int}}, SVector{2, Float64}}}}
    @test collect(ex) == [
        (2.0, ID(OID(Index(CPID(1, 1), TID(1)), SVector(0.0, 0.0), SVector(0.0, 0.0)),
                    OID(Index(CPID(1, 1), TID(2)), SVector(0.0, 0.0), SVector(0.0, 0.0))
                    ))
        ]

    @test @couplings(TCoupling(1.0, (1, 1))) == Couplings(TCoupling(1.0, (1, 1)))
    @test @couplings(TCoupling(1.0, (1,))+TCoupling(1.0, (2,))) == TCoupling(1.0, (1,))+TCoupling(1.0, (2,))
end

@testset "TermFunction" begin
    ta = TermAmplitude()
    @test ta() == 1

    ta = TermAmplitude(x->x+3.0)
    @test ta(1.0) == 4.0

    tcs = TCoupling(1.0, (1, 1)) + TCoupling(2.0, (2, 2))
    termcouplings = TermCouplings(tcs)
    @test termcouplings == deepcopy(TermCouplings(tcs))
    @test isequal(termcouplings, deepcopy(TermCouplings(tcs)))
    @test termcouplings() == tcs

    tcs1 = TCoupling(1.0, (1, 1)) + TCoupling(1.0, (2, 2))
    tcs2 = TCoupling(1.0, (2, 1)) + TCoupling(1.0, (1, 2))
    termcouplings = TermCouplings(i->(tcs1, tcs2)[(i-1)%2+1])
    @test termcouplings(1) == tcs1
    @test termcouplings(2) == tcs2

    @test ismodulatable(TermModulate{Val{true}}) == ismodulatable(TermModulate{<:Function}) == true
    @test ismodulatable(TermModulate{Val{false}}) == false
    termmodulate = TermModulate(:t)
    @test termmodulate(t=1) == 1
    @test isnothing(termmodulate(mu=1))
    @test isnothing(termmodulate())

    termmodulate = TermModulate(:t, t->t*2.0)
    @test termmodulate(2) == 4
end

@testset "Term" begin
    point = Point(CPID(1, 1), (0.0, 0.0), (0.0, 0.0))
    hilbert = Hilbert(point.pid=>TFock(2))
    term = Term{:Mu, 2}(:mu, 1.5, 0, couplings=@couplings(TCoupling(1.0, (2, 1))), amplitude=(bond->3), modulate=false)
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
    hilbert = Hilbert(pid=>TFock(2) for pid in [PID(1), PID(2)])
    tcs1 = Couplings(TCoupling(1.0, (2, 2)))
    tcs2 = Couplings(TCoupling(1.0, (1, 1)))
    term = Term{:Mu, 2}(:mu, 1.5, 0, couplings=bond->(tcs1, tcs2)[bond.pid.site%2+1], amplitude=bond->3, modulate=true)
    @test term|>ismodulatable == term|>typeof|>ismodulatable == true
    @test repr(term, p1, hilbert) == "mu: 4.5 ph(1, 1)"
    @test repr(term, p2, hilbert) == "mu: 4.5 ph(2, 2)"
    @test one(term) == replace(term, value=1.0)
    @test zero(term) == replace(term, value=0.0)
    @test term.modulate(mu=4.0) == 4.0
    @test isnothing(term.modulate(t=1.0))
    @test update!(term, mu=4.25) == replace(term, value=4.25)
    @test term.value == 4.25
end

@testset "expand" begin
    point = Point(PID(1), (0.0, 0.0), (0.0, 0.0))
    hilbert = Hilbert(point.pid=>TFock(2))
    term = Term{:Mu, 2}(:mu, 1.5, 0, couplings=@couplings(TCoupling(1.0, (2, 1))), amplitude=bond->3.0, modulate=true)
    operators = Operators(
        Operator(+2.25, ID(
            OID(Index(PID(1), TID(2)), SVector(0.0, 0.0), SVector(0.0, 0.0)),
            OID(Index(PID(1), TID(1)), SVector(0.0, 0.0), SVector(0.0, 0.0))
            ))
        )
    @test expand(term, point, hilbert, true) == operators
    @test expand(term, point, hilbert, false) == operators*2

    bond = Bond(1, Point(PID(2), (1.5, 1.5), (1.0, 1.0)), Point(PID(1), (0.5, 0.5), (0.0, 0.0)))
    hilbert = Hilbert(pid=>TFock(2) for pid in [PID(1), PID(2)])
    term = Term{:Hp, 2}(:t, 1.5, 1, couplings=@couplings(TCoupling(1.0, (2, 1))), amplitude=bond->3.0, modulate=true)
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
    hilbert = Hilbert(pid=>TFock(2) for pid in lattice.pids)
    term = Term{:Mu, 2}(:mu, 1.5, 0, couplings=@couplings(TCoupling(1.0, (2, 1))), amplitude=bond->3.0, modulate=true)
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
    hilbert = Hilbert{TFock}(pid->TFock(2), lattice.pids)
    table = Table(hilbert, OIDToTuple(:scope, :site))
    t = Term{:Hp, 2}(:t, 2.0, 1, couplings=@couplings(TCoupling(1.0, (2, 1))))
    μ = Term{:Mu, 2}(:μ, 1.0, 0, couplings=@couplings(TCoupling(1.0, (2, 1))), modulate=true)
    tops1 = expand(t, filter(acrossbonds, bonds, Val(:exclude)), hilbert, true, table=table)
    tops2 = expand(one(t), filter(acrossbonds, bonds, Val(:include)), hilbert, true, table=table)
    μops = expand(one(μ), filter(zerothbonds, bonds, Val(:include)), hilbert, true, table=table)

    optp = Operator{Float, ID{OID{Index{PID, TID{Int}}, SVector{2, Float}}, 2}}
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
