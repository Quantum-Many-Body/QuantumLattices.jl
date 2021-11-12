using Test
using Printf: @printf
using LinearAlgebra: tr
using StaticArrays: SVector
using QuantumLattices.Essentials.Frameworks
using QuantumLattices.Essentials: update!, reset!, register!
using QuantumLattices.Essentials.Spatials: Point, PID, Bond, Bonds, Lattice, acrossbonds, zerothbonds
using QuantumLattices.Essentials.DegreesOfFreedom: SimpleIID, SimpleInternal, IIDSpace, Coupling, Subscript, Subscripts, SubscriptsID
using QuantumLattices.Essentials.DegreesOfFreedom: Term, Hilbert, Index, Table, OID, OIDToTuple, plain, @couplings
using QuantumLattices.Essentials.QuantumOperators: ID, Operator, Operators, id, idtype, Identity
using QuantumLattices.Prerequisites: Float
using QuantumLattices.Interfaces:  expand!, expand, add!
using QuantumLattices.Prerequisites.Traits: contentnames
using QuantumLattices.Prerequisites.CompositeStructures: NamedContainer

import QuantumLattices.Essentials.Frameworks: Parameters
import QuantumLattices.Essentials: prepare!, run!, update!
import QuantumLattices.Essentials.DegreesOfFreedom: ishermitian, couplingcenters
import QuantumLattices.Prerequisites.VectorSpaces: shape, ndimshape

struct FID{N<:Union{Int, Symbol}} <: SimpleIID
    nambu::N
end
@inline Base.adjoint(sl::FID{Int}) = FID(3-sl.nambu)

struct FFock <: SimpleInternal{FID{Int}}
    nnambu::Int
end
@inline shape(f::FFock) = (1:f.nnambu,)
@inline ndimshape(::Type{FFock}) = 1
@inline FID(i::CartesianIndex, vs::FFock) = FID(i.I...)
@inline CartesianIndex(did::FID{Int}, vs::FFock) = CartesianIndex(did.nambu)
@inline shape(iidspace::IIDSpace{FID{Symbol}, FFock}) = (1:iidspace.internal.nnambu,)
@inline shape(iidspace::IIDSpace{FID{Int}, FFock}) = (iidspace.iid.nambu:iidspace.iid.nambu,)
@inline shape(iidspace::IIDSpace{FID{Symbol}, FFock}) = (1:iidspace.internal.nnambu,)
@inline shape(iidspace::IIDSpace{FID{Int}, FFock}) = (iidspace.iid.nambu:iidspace.iid.nambu,)

const FCoupling{V, I<:ID{FID}, C<:Subscripts, CI<:SubscriptsID} = Coupling{V, I, C, CI}
@inline couplingcenters(::(Coupling{V, <:ID{FID}} where V), ::Bond, ::Val) = (1, 2)
@inline FCoupling(value, nambus::Tuple{Vararg{Int}}) = Coupling(value, ID(FID, nambus), Subscripts((nambu=Subscript(nambus),)))
@inline FCoupling(value, nambus::Subscript) = Coupling(value, ID(FID, convert(Tuple, nambus)), Subscripts((nambu=nambus,)))

@inline ishermitian(::Type{<:Term{:Mu}}) = true
@inline ishermitian(::Type{<:Term{:Hp}}) = false

@testset "Parameters" begin
    ps1 = Parameters{(:t1, :t2, :U)}(1.0im, 1.0, 2.0)
    ps2 = Parameters{(:t1, :U)}(1.0im, 2.0)
    ps3 = Parameters{(:t1, :U)}(1.0im, 2.1)
    @test match(ps1, ps2) == true
    @test match(ps1, ps3) == false
end

@testset "Generator" begin
    @test contentnames(AbstractGenerator) == (:operators, :table, :boundary)
    @test contentnames(Generator) == (:operators, :terms, :bonds, :hilbert, :half, :table, :boundary)
    @test contentnames(SimplifiedGenerator) == (:parameters, :operators, :table, :boundary)

    lattice = Lattice("Tuanzi", [Point(PID(1), (0.0, 0.0), (0.0, 0.0)), Point(PID(2), (0.5, 0.0), (0.0, 0.0))], vectors=[[1.0, 0.0]], neighbors=1)
    bonds = Bonds(lattice)
    hilbert = Hilbert{FFock}(pid->FFock(2), lattice.pids)
    table = Table(hilbert, OIDToTuple(:scope, :site))
    t = Term{:Hp}(:t, 2.0, 1, couplings=@couplings(FCoupling(1.0, (2, 1))))
    μ = Term{:Mu}(:μ, 1.0, 0, couplings=@couplings(FCoupling(1.0, (2, 1))), modulate=true)
    tops₁ = expand(t, filter(acrossbonds, bonds, Val(:exclude)), hilbert, half=true, table=table)
    tops₂ = expand(one(t), filter(acrossbonds, bonds, Val(:include)), hilbert, half=true, table=table)
    μops = expand(one(μ), filter(zerothbonds, bonds, Val(:include)), hilbert, half=true, table=table)
    i = Identity()

    optp = Operator{Float, ID{OID{Index{PID, FID{Int}}, SVector{2, Float}}, 2}}
    entry = Entry(tops₁, NamedContainer{(:μ,)}((μops,)), NamedContainer{(:t, :μ)}((tops₂, Operators{optp}())))
    @test entry == deepcopy(entry) && isequal(entry, deepcopy(entry))
    @test entry == Entry((t, μ), bonds, hilbert, half=true, table=table)
    @test entry|>valtype == entry|>typeof|>valtype == Operators{optp, idtype(optp)}
    @test expand!(Operators{optp}(), entry, plain, t=2.0, μ=1.5) == tops₁+tops₂*2.0+μops*1.5
    @test empty(entry) == empty!(deepcopy(entry))
    @test empty(entry) == Entry(empty(μops), NamedContainer{(:μ,)}((empty(μops),)), NamedContainer{(:t, :μ)}((empty(μops), empty(μops))))
    @test merge!(empty(entry), entry) == merge(entry, entry) == entry
    @test reset!(deepcopy(entry), (t, μ), bonds, hilbert, half=true, table=table) == entry
    @test i(entry) == entry

    cgen = Generator((t, μ), bonds, hilbert; half=true, table=table, boundary=plain)
    @test cgen == deepcopy(cgen) && isequal(cgen, deepcopy(cgen))
    @test cgen|>eltype == cgen|>typeof|>eltype == optp
    @test cgen|>valtype == cgen|>typeof|>valtype == Operators{optp, idtype(optp)}
    @test collect(cgen) == collect(expand(cgen))
    @test Parameters(cgen) == Parameters{(:t, :μ)}(2.0, 1.0)
    @test expand!(Operators{optp}(), cgen) == expand(cgen) == tops₁+tops₂*2.0+μops
    @test expand(cgen, :t) == tops₁+tops₂*2.0
    @test expand(cgen, :μ) == μops
    @test expand(cgen, 1)+expand(cgen, 2)+expand(cgen, 3)+expand(cgen, 4) == expand(cgen)
    @test expand(cgen, :μ, 1)+expand(cgen, :μ, 2) == μops
    @test expand(cgen, :t, 3) == tops₁
    @test expand(cgen, :t, 4) == tops₂*2.0
    @test empty!(deepcopy(cgen)) == Generator((t, μ), empty(bonds), empty(hilbert), half=true, table=empty(table), boundary=plain) == empty(cgen)
    @test reset!(empty(cgen), lattice) == cgen
    @test update!(cgen, μ=1.5)|>expand == tops₁+tops₂*2.0+μops*1.5

    params = Parameters{(:t, :μ)}(2.0, 1.0)
    sgen = SimplifiedGenerator(params, entry, table=table, boundary=plain)
    @test Parameters(sgen) == params
    @test expand!(Operators{optp}(), sgen) == expand(sgen) == tops₁+tops₂*2.0+μops
    @test empty!(deepcopy(sgen)) == SimplifiedGenerator(params, empty(entry), table=empty(table), boundary=plain) == empty(sgen)
    @test reset!(empty(sgen), entry, table=table) == sgen
    @test update!(sgen, μ=1.5)|>expand == tops₁+tops₂*2.0+μops*1.5
    @test i(cgen) == sgen
end

mutable struct VCA <: Engine
    t::Float
    U::Float
    dim::Int
end
function update!(vca::VCA; kwargs...)
    vca.t = get(kwargs, :t, vca.t)
    vca.U = get(kwargs, :U, vca.U)
    return vca
end
@inline Parameters(vca::VCA) = Parameters{(:t, :U)}(vca.t, vca.U)

mutable struct GF <: Action
    count::Int
end
@inline prepare!(gf::GF, vca::VCA) = Matrix{Float}(undef, vca.dim, vca.dim)
@inline run!(alg::Algorithm{VCA}, assign::Assignment{GF}) = (assign.action.count += 1; assign.data[:, :] .= alg.engine.t+alg.engine.U)

mutable struct DOS <: Action
    mu::Float
end
@inline update!(eb::DOS; kwargs...) = (eb.mu = get(kwargs, :mu, eb.mu); eb)
@inline prepare!(dos::DOS, vca::VCA) = 0.0
@inline function run!(alg::Algorithm{VCA}, assign::Assignment{DOS})
    rundependences!(alg, assign)
    gf = alg.assignments[first(assign.dependences)]
    assign.data = tr(gf.data)
end
@inline dosmap(params::Parameters) = Parameters{(:t, :U, :mu)}(params.t, params.U, -params.U/2)

@testset "Framework" begin
    vca = VCA(1.0, 8.0, 4)
    @test vca == deepcopy(vca)
    @test isequal(vca, deepcopy(vca))
    @test repr(vca) == "VCA"
    @test string(vca) == "VCA"

    gf = GF(0)
    @test gf == deepcopy(gf)
    @test isequal(gf, deepcopy(gf))
    @test update!(gf) == gf

    vca = Algorithm("test", vca)
    @test vca == deepcopy(vca)
    @test isequal(vca, deepcopy(vca))
    @test repr(vca, ≠(:U)) == "test(VCA)_1.0"
    @test repr(vca) == "test(VCA)_1.0_8.0"

    add!(vca, :GF, GF(0))
    rex = r"Action DOS\(DOS\)\: time consumed [0-9]*\.[0-9]*(e[+-][0-9]*)*s."
    @test_logs (:info, rex) register!(vca, :DOS, DOS(-3.5), parameters=(U=7.0,), map=dosmap, dependences=(:GF,))

    dos = vca.assignments[:DOS]
    @test dos == deepcopy(dos)
    @test isequal(dos, deepcopy(dos))
    @test valtype(dos) == valtype(typeof(dos)) == Float
    @test dos.data == 32.0
    @test dos.action.mu == -3.5
    @test nameof(vca, dos) == "test(VCA)_1.0_7.0_DOS"

    update!(dos, U=6.0)
    run!(vca, :DOS, false)
    @test dos.parameters == (t=1.0, U=6.0)
    @test dos.data == 28.0
    @test dos.action.mu == -3.0
end
