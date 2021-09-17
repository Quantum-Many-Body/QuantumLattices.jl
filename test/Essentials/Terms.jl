using Test
using Printf: @sprintf
using StaticArrays: SVector
using QuantumLattices.Essentials.Terms
using QuantumLattices.Essentials.QuantumAlgebras: ID, id, idtype
using QuantumLattices.Essentials.Spatials: Point, PID, CPID, Bond, Bonds, Lattice, acrossbonds, zerothbonds
using QuantumLattices.Essentials.DegreesOfFreedom: SimpleIID, SimpleInternal, CompositeIID, Hilbert, Index, Table, OID, OIDToTuple, Operator, Operators, plain
using QuantumLattices.Essentials: kind, update!, reset!
using QuantumLattices.Interfaces: rank, expand!, expand, ⊗
using QuantumLattices.Prerequisites: Float, decimaltostr
using QuantumLattices.Prerequisites.Traits: parameternames, isparameterbound, contentnames, getcontent
using QuantumLattices.Prerequisites.CompositeStructures: NamedContainer

import QuantumLattices.Essentials.DegreesOfFreedom: isHermitian
import QuantumLattices.Essentials.Terms: couplingcenters, abbr
import QuantumLattices.Prerequisites.VectorSpaces: shape, ndimshape

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
@inline shape(iidspace::IIDSpace{TID{Symbol}, TFock}) = (1:iidspace.internal.nnambu,)
@inline shape(iidspace::IIDSpace{TID{Int}, TFock}) = (iidspace.iid.nambu:iidspace.iid.nambu,)

const TCoupling{V, I<:ID{TID}, C<:Subscripts, CI<:SubscriptsID} = Coupling{V, I, C, CI}
@inline Base.repr(tc::(Coupling{V, <:ID{TID}} where V)) = @sprintf "%s ph(%s)" decimaltostr(tc.value) join(tc.cid.nambus, ", ")
@inline couplingcenters(::(Coupling{V, <:ID{TID}} where V), ::Bond, ::Val) = (1, 2)
@inline TCoupling(value, nambus::Tuple{Vararg{Int}}) = Coupling(value, ID(TID, nambus), Subscripts((nambu=Subscript(nambus),)))
@inline TCoupling(value, nambus::Subscript) = Coupling(value, ID(TID, convert(Tuple, nambus)), Subscripts((nambu=nambus,)))

abbr(::Type{<:Term{:Mu}}) = :mu
abbr(::Type{<:Term{:Hp}}) = :hp
isHermitian(::Type{<:Term{:Mu}}) = true
isHermitian(::Type{<:Term{:Hp}}) = false

@testset "Subscript" begin
    sub = Subscript(4)
    @test contentnames(typeof(sub)) == (:contents, :rep, :constraint)
    @test getcontent(sub, :contents) == sub.pattern
    @test sub==deepcopy(sub) && isequal(sub, deepcopy(sub))
    @test string(sub) == "[* * * *]"
    @test rank(sub) == rank(typeof(sub)) == 4
    @test match(sub, (2, 2, 2, 2))

    sub = Subscript((1, 2, 3, 4))
    @test string(sub) == "[1 2 3 4]"
    @test sub == subscript"[1 2 3 4]" == subscript"[1 2 3 4]; false"
    @test rank(sub) == 4
    @test match(sub, (1, 3, 3, 4))

    sub = Subscript((1, 2, 2, 1), true)
    @test string(sub) == "[1 2 2 1]"
    @test sub == subscript"[1 2 2 1]; true"
    @test rank(sub) == 4
    @test match(sub, (1, 2, 2, 1)) && !match(sub, (1, 1, 2, 1))

    sub = subscript"[x₁ x₂ x₁ x₂]"
    @test string(sub) == "[x₁ x₂ x₁ x₂]"
    @test rank(sub) == 4
    @test match(sub, (1, 2, 1, 2)) && !match(sub, (1, 2, 2, 1))

    sub = subscript"[x₁ 4 4 x₂](x₁ < x₂)"
    @test string(sub) == "[x₁ 4 4 x₂](x₁ < x₂)"
    @test rank(sub) == rank(typeof(sub)) == 4
    @test match(sub, (1, 4, 4, 2)) && !match(sub, (2, 4, 4, 1))
end

@testset "Subscripts" begin
    subscripts = Subscripts((nambu=subscript"[x y](x > y)",), (nambu=subscript"[x x y y](x < y)",))
    @test string(subscripts) == "nambu[x y](x > y) × nambu[x x y y](x < y)"
    @test repr(subscripts, 1:2, :nambu) == "[x y](x > y)×[x x y y](x < y)"
    @test repr(subscripts, 1, :nambu) == "[x y](x > y)"
    @test repr(subscripts, 2, :nambu) == "[x x y y](x < y)"

    @test rank(subscripts) == rank(typeof(subscripts)) == 6
    @test rank(subscripts, 1) == rank(typeof(subscripts), 1) == 2
    @test rank(subscripts, 2) == rank(typeof(subscripts), 2) == 4
    @test match(subscripts, (TID(2), TID(1), TID(2), TID(2), TID(3), TID(3)))
    @test match(subscripts, TID(2)⊗TID(1)⊗TID(2)⊗TID(2)⊗TID(3)⊗TID(3))
    @test !match(subscripts, TID(1)⊗TID(2)⊗TID(2)⊗TID(2)⊗TID(3)⊗TID(3))
    @test !match(subscripts, TID(2)⊗TID(1)⊗TID(3)⊗TID(3)⊗TID(2)⊗TID(2))
    @test subscripts == Subscripts((nambu=subscript"[x y](x > y)",))*Subscripts((nambu=subscript"[x x y y](x < y)",))

    subscriptsid = SubscriptsID(subscripts)
    @test subscriptsid == SubscriptsID((1:2=>("[x y](x > y)",), 3:6=>("[x x y y](x < y)",)))
    @test idtype(subscripts) == idtype(typeof(subscripts)) == SubscriptsID{NTuple{2, Pair{UnitRange{Int}, Tuple{String}}}}
end

@testset "IIDSpace" begin
    tid₁, tid₂, it = TID(2), TID(:σ), TFock(2)
    iidspace = IIDSpace(tid₁⊗tid₂, it⊗it)
    @test eltype(iidspace) == eltype(typeof(iidspace)) == CompositeIID{Tuple{TID{Int}, TID{Int}}}
    @test kind(iidspace) == kind(typeof(iidspace)) == :info
    @test shape(iidspace) == (2:2, 1:2)
    @test ndimshape(iidspace) == ndimshape(typeof(iidspace)) == 2
    for i = 1:length(iidspace)
        iid = iidspace[i]
        @test iidspace[CartesianIndex(iid, iidspace)] == iid
    end
    @test expand((tid₁, tid₂), (it, it)) == iidspace
    @test collect(iidspace) == [TID(2)⊗TID(1), TID(2)⊗TID(2)]
end

@testset "couplings" begin
    @test parameternames(Coupling) == (:value, :cid, :subscripts, :subscriptsid)
    @test contentnames(Coupling) == (:value, :id, :subscripts)
    @test isparameterbound(Coupling, :cid, Tuple{TID{Int}}) == false
    @test isparameterbound(Coupling, :cid, ID{TID{Int}}) == true
    @test isparameterbound(Coupling, :subscripts, Subscripts{Tuple{NamedTuple{(:nambu,), Tuple{Subscript{Tuple{Int}, typeof(noconstrain)}}}}}) == false
    @test isparameterbound(Coupling, :subscripts, Subscripts) == true
    @test isparameterbound(Coupling, :subscriptsid, SubscriptsID{Tuple{}}) == false
    @test isparameterbound(Coupling, :subscriptsid, SubscriptsID) == true

    tc = TCoupling(2.0, (2,))
    @test id(tc) == ID(CompositeIID(tc.cid), SubscriptsID(tc.subscripts))
    @test rank(tc) == rank(typeof(tc)) == 1
    @test tc == Coupling(2.0, id(tc), tc.subscripts)
    @test ID{SimpleIID}(tc) == tc.cid
    @test Subscripts(tc) == tc.subscripts

    point = Point(CPID(1, 1), (0.0, 0.0), (0.0, 0.0))
    hilbert = Hilbert(point.pid=>TFock(2))
    @test couplingcenters(tc, point, Val(:Mu)) == (1,)
    @test couplingpoints(tc, point, Val(:Mu)) == (point,)
    @test couplinginternals(tc, point, hilbert, Val(:Mu)) == (TFock(2),)

    tc₁ = TCoupling(1.5, (1, 2))
    tc₂ = TCoupling(2.0, subscript"[a b](a < b)")
    @test tc₁*tc₂ == Coupling(3.0, ID(TID(1), TID(2), TID(:a), TID(:b)), Subscripts((nambu=Subscript((1, 2)),), (nambu=subscript"[a b](a < b)",)))

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
    @test valtype(termcouplings) == valtype(typeof(termcouplings)) == typeof(tcs)
    @test termcouplings() == tcs

    fx = i -> i%2==1 ? (TCoupling(1.0, (1, 1))+TCoupling(1.0, (2, 2))) : (TCoupling(1.0, (2, 1))+TCoupling(1.0, (1, 2)))
    termcouplings = TermCouplings(fx)
    @test termcouplings == TermCouplings{typejoin(typeof(fx(1)), typeof(fx(2)))}(fx)
    @test valtype(termcouplings) == valtype(typeof(termcouplings)) == typejoin(typeof(fx(1)), typeof(fx(2)))
    @test termcouplings(1) == fx(1)
    @test termcouplings(2) == fx(2)

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
    term = Term{:Mu}(:mu, 1.5, 0, couplings=@couplings(TCoupling(1.0, (2, 1))), amplitude=(bond->3), modulate=false)
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
    term = Term{:Mu}(:mu, 1.5, 0,
        couplings=bond->bond.pid.site%2==0 ? Couplings(TCoupling(1.0, (2, 2))) : Couplings(TCoupling(1.0, (1, 1))),
        amplitude=bond->3,
        modulate=true
    )
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
    term = Term{:Mu}(:mu, 1.5, 0, couplings=@couplings(TCoupling(1.0, (2, 1))), amplitude=bond->3.0, modulate=true)
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
    term = Term{:Hp}(:t, 1.5, 1, couplings=@couplings(TCoupling(1.0, (2, 1))), amplitude=bond->3.0, modulate=true)
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
    term = Term{:Mu}(:mu, 1.5, 0, couplings=@couplings(TCoupling(1.0, (2, 1))), amplitude=bond->3.0, modulate=true)
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
    @test contentnames(AbstractGenerator) == (:table, :boundary, :operators)
    @test contentnames(AbstractCompleteGenerator) == (:terms, :bonds, :hilbert, :half, :table, :boundary, :operators)
    @test contentnames(AbstractSimplifiedGenerator) == (:parameters, :table, :boundary, :operators)

    lattice = Lattice("Tuanzi", [Point(PID(1), (0.0, 0.0), (0.0, 0.0)), Point(PID(2), (0.5, 0.0), (0.0, 0.0))], vectors=[[1.0, 0.0]], neighbors=1)
    bonds = Bonds(lattice)
    hilbert = Hilbert{TFock}(pid->TFock(2), lattice.pids)
    table = Table(hilbert, OIDToTuple(:scope, :site))
    t = Term{:Hp}(:t, 2.0, 1, couplings=@couplings(TCoupling(1.0, (2, 1))))
    μ = Term{:Mu}(:μ, 1.0, 0, couplings=@couplings(TCoupling(1.0, (2, 1))), modulate=true)
    tops₁ = expand(t, filter(acrossbonds, bonds, Val(:exclude)), hilbert, true, table=table)
    tops₂ = expand(one(t), filter(acrossbonds, bonds, Val(:include)), hilbert, true, table=table)
    μops = expand(one(μ), filter(zerothbonds, bonds, Val(:include)), hilbert, true, table=table)

    optp = Operator{Float, ID{OID{Index{PID, TID{Int}}, SVector{2, Float}}, 2}}
    genops = GenOperators(tops₁, NamedContainer{(:μ,)}((μops,)), NamedContainer{(:t, :μ)}((tops₂, Operators{optp|>idtype, optp}())))
    @test genops == deepcopy(genops) && isequal(genops, deepcopy(genops))
    @test genops == GenOperators((t, μ), bonds, hilbert, true, table=table)
    @test genops|>eltype == genops|>typeof|>eltype == optp
    @test genops|>idtype == genops|>typeof|>idtype == optp|>idtype
    @test expand!(Operators{idtype(optp), optp}(), genops, plain, t=2.0, μ=1.5) == tops₁+tops₂*2.0+μops*1.5
    @test empty(genops) == empty!(deepcopy(genops))
    @test empty(genops) == GenOperators(empty(μops), NamedContainer{(:μ,)}((empty(μops),)), NamedContainer{(:t, :μ)}((empty(μops), empty(μops))))
    @test merge!(empty(genops), genops) == merge(genops, genops) == genops
    @test reset!(deepcopy(genops), (t, μ), bonds, hilbert, true, table=table) == genops

    gen = Generator((t, μ), bonds, hilbert; half=true, table=table)
    @test gen == deepcopy(gen) && isequal(gen, deepcopy(gen))
    @test Parameters(gen) == Parameters{(:t, :μ)}(2.0, 1.0)
    @test expand!(Operators{idtype(optp), optp}(), gen) == expand(gen) == tops₁+tops₂*2.0+μops
    @test expand(gen, :t) == tops₁+tops₂*2.0
    @test expand(gen, :μ) == μops
    @test expand(gen, 1)+expand(gen, 2)+expand(gen, 3)+expand(gen, 4) == expand(gen)
    @test expand(gen, :μ, 1)+expand(gen, :μ, 2) == μops
    @test expand(gen, :t, 3) == tops₁
    @test expand(gen, :t, 4) == tops₂*2.0
    @test empty!(deepcopy(gen)) == Generator((t, μ), empty(bonds), empty(hilbert), half=true, table=empty(table), boundary=plain) == empty(gen)
    @test reset!(empty(gen), lattice) == gen
    @test update!(gen, μ=1.5)|>expand == tops₁+tops₂*2.0+μops*1.5

    params = Parameters{(:t, :μ)}(2.0, 1.0)
    gen = SimplifiedGenerator(params, genops, table=table, boundary=plain)
    @test Parameters(gen) == params
    @test expand!(Operators{idtype(optp), optp}(), gen) == expand(gen) == tops₁+tops₂*2.0+μops
    @test empty!(deepcopy(gen)) == SimplifiedGenerator(params, empty(genops), table=empty(table), boundary=plain) == empty(gen)
    @test reset!(empty(gen), genops, table=table) == gen
    @test update!(gen, μ=1.5)|>expand == tops₁+tops₂*2.0+μops*1.5
end
