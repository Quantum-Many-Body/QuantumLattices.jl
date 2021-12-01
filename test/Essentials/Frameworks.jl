using Test
using Printf: @printf
using LinearAlgebra: tr, dot
using StaticArrays: SVector
using QuantumLattices.Essentials.Frameworks
using QuantumLattices.Essentials: reset!
using QuantumLattices.Essentials.Spatials: Point, AbstractPID, PID, Bond, Bonds, Lattice, acrossbonds, zerothbonds, decompose
using QuantumLattices.Essentials.DegreesOfFreedom: SimpleIID, SimpleInternal, IIDSpace, Coupling, Subscript, Subscripts, SubscriptsID
using QuantumLattices.Essentials.DegreesOfFreedom: Term, Hilbert, Index, Table, OID, OIDToTuple, @couplings, TwistedBoundaryCondition
using QuantumLattices.Essentials.QuantumOperators: ID, Operator, Operators, id, idtype, Identity
using QuantumLattices.Prerequisites: Float
using QuantumLattices.Interfaces:  expand!, expand, add!
using QuantumLattices.Prerequisites.Traits: contentnames

import QuantumLattices.Essentials.DegreesOfFreedom: Parameters
import QuantumLattices.Essentials: prepare!, run!, update!
import QuantumLattices.Essentials.DegreesOfFreedom: ishermitian, couplingcenters
import QuantumLattices.Prerequisites.VectorSpaces: shape, ndimshape

struct FID{N<:Union{Int, Symbol}} <: SimpleIID
    nambu::N
end
@inline Base.adjoint(sl::FID{Int}) = FID(3-sl.nambu)
function Base.angle(id::OID{<:Index{<:AbstractPID, FID{Int}}}, vectors::AbstractVector{<:AbstractVector{Float}}, values::AbstractVector{Float})
    phase=  (length(vectors) == 1) ? 2pi*dot(decompose(id.icoord, vectors[1]), values) :
            (length(vectors) == 2) ? 2pi*dot(decompose(id.icoord, vectors[1], vectors[2]), values) :
            (length(vectors) == 3) ? 2pi*dot(decompose(id.icoord, vectors[1], vectors[2], vectors[3]), values) :
            error("angle error: not supported number of input basis vectors.")
    (id.index.iid.nambu == 1) ? phase : -phase
end

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

@testset "Generator" begin
    @test contentnames(AbstractGenerator) == (:operators, :table)
    @test contentnames(Generator) == (:operators, :terms, :bonds, :hilbert, :half, :table)
    @test contentnames(SimplifiedGenerator) == (:operators, :table, :sourceid)

    lattice = Lattice(:Tuanzi, [Point(PID(1), [0.0]), Point(PID(2), [0.5])], vectors=[[1.0]], neighbors=1)
    bonds = Bonds(lattice)
    hilbert = Hilbert{FFock}(pid->FFock(2), lattice.pids)
    table = Table(hilbert, OIDToTuple(:scope, :site))
    boundary = TwistedBoundaryCondition{(:θ,)}([0.1], lattice.vectors)

    t = Term{:Hp}(:t, 2.0, 1, couplings=@couplings(FCoupling(1.0, (2, 1))))
    μ = Term{:Mu}(:μ, 1.0, 0, couplings=@couplings(FCoupling(1.0, (2, 1))), modulate=true)
    tops₁ = expand(t, filter(acrossbonds, bonds, Val(:exclude)), hilbert, half=true, table=table)
    tops₂ = boundary(expand(one(t), filter(acrossbonds, bonds, Val(:include)), hilbert, half=true, table=table))
    μops = expand(one(μ), filter(zerothbonds, bonds, Val(:include)), hilbert, half=true, table=table)
    i = Identity()

    optp = Operator{Complex{Float}, ID{OID{Index{PID, FID{Int}}, SVector{1, Float}}, 2}}
    entry = Entry(tops₁, (μ=μops,), (t=tops₂, μ=Operators{optp}()), (t=2.0, μ=1.0), boundary)
    @test entry == Entry((t, μ), bonds, hilbert; half=true, table=table, boundary=boundary)
    @test isequal(entry, i(entry))
    @test Parameters(entry) == (t=2.0, μ=1.0, θ=0.1)
    @test entry|>valtype == entry|>typeof|>valtype == Operators{optp, idtype(optp)}
    @test expand(entry) == expand!(Operators{optp}(), entry) ≈ tops₁+tops₂*2.0+μops

    another = Entry(empty(tops₁), (μ=empty(μops),), (t=empty(tops₂), μ=Operators{optp}()), (t=2.0, μ=1.0), boundary)
    @test empty(entry) == empty!(deepcopy(entry)) == another
    @test merge!(deepcopy(another), entry) == merge(another, entry) == entry

    nb = update!(deepcopy(boundary); θ=0.5)
    another = Entry(tops₁, (μ=μops,), (t=nb(tops₂, origin=[0.1]), μ=Operators{optp}()), (t=2.0, μ=1.5), nb)
    dest = deepcopy(entry)
    update!(dest; μ=1.5, θ=0.5)
    @test dest == another
    @test expand(dest) ≈ tops₁+another.boundary(tops₂, origin=[0.1])*2.0+μops*1.5
    dest = deepcopy(entry)
    update!(dest, i, another)
    @test dest == another
    @test expand(dest) ≈ tops₁+another.boundary(tops₂, origin=[0.1])*2.0+μops*1.5

    @test reset!(deepcopy(another), (t, μ), bonds, hilbert, half=true, table=table, boundary=boundary) == entry
    @test reset!(deepcopy(another), i, entry) == entry

    cgen = Generator((t, μ), bonds, hilbert; half=true, table=table, boundary=boundary)
    @test cgen == deepcopy(cgen) && isequal(cgen, deepcopy(cgen))
    @test cgen|>eltype == cgen|>typeof|>eltype == optp
    @test cgen|>valtype == cgen|>typeof|>valtype == Operators{optp, idtype(optp)}
    @test collect(cgen) == collect(expand(cgen))
    @test Parameters(cgen) == (t=2.0, μ=1.0, θ=0.1)
    @test expand!(Operators{optp}(), cgen) == expand(cgen) ≈ tops₁+tops₂*2.0+μops
    @test expand(cgen, :t) ≈ tops₁+tops₂*2.0
    @test expand(cgen, :μ) ≈ μops
    @test expand(cgen, 1)+expand(cgen, 2)+expand(cgen, 3)+expand(cgen, 4) ≈ expand(cgen)
    @test expand(cgen, :μ, 1)+expand(cgen, :μ, 2) ≈ μops
    @test expand(cgen, :t, 3) ≈ tops₁
    @test expand(cgen, :t, 4) ≈ tops₂*2.0
    @test empty!(deepcopy(cgen)) == Generator((t, μ), empty(bonds), empty(hilbert), half=true, table=empty(table), boundary=boundary) == empty(cgen)
    @test reset!(empty(cgen), lattice) == cgen
    @test update!(cgen, μ=1.5)|>expand ≈ tops₁+tops₂*2.0+μops*1.5

    sgen = i(cgen; table=table)
    @test sgen == SimplifiedGenerator(cgen.operators, table=table, sourceid=objectid((i, cgen)))
    @test Parameters(sgen) == (t=2.0, μ=1.5, θ=0.1)
    @test expand!(Operators{optp}(), sgen) == expand(sgen) ≈ tops₁+tops₂*2.0+μops*1.5
    @test empty!(deepcopy(sgen)) == SimplifiedGenerator(empty(cgen.operators), table=empty(table), sourceid=objectid((i, cgen))) == empty(sgen)
    @test update!(sgen, μ=3.5)|>expand ≈ tops₁+tops₂*2.0+μops*3.5
    @test update!(sgen, i, cgen)|>expand ≈ tops₁+tops₂*2.0+μops*1.5
    @test reset!(empty(sgen), cgen.operators; table=table) == sgen
    @test reset!(empty(sgen), i, cgen; table=table) == sgen
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

    vca = Algorithm(:test, vca)
    @test vca == deepcopy(vca)
    @test isequal(vca, deepcopy(vca))
    @test repr(vca, ≠(:U)) == "test(VCA)_1.0"
    @test repr(vca) == "test(VCA)_1.0_8.0"

    add!(vca, :GF, GF(0))
    rex = r"Action DOS\(DOS\)\: time consumed [0-9]*\.[0-9]*(e[+-][0-9]*)*s."
    @test_logs (:info, rex) vca(:DOS, DOS(-3.5), parameters=(U=7.0,), map=dosmap, dependences=(:GF,))

    dos = vca.assignments[:DOS]
    @test dos == deepcopy(dos)
    @test isequal(dos, deepcopy(dos))
    @test valtype(dos) == valtype(typeof(dos)) == Float
    @test dos.data == 32.0
    @test dos.action.mu == -3.5
    @test nameof(vca, dos) == "test(VCA)_1.0_7.0_DOS"

    update!(dos, U=6.0)
    vca(:DOS, info=false)
    @test dos.parameters == (t=1.0, U=6.0)
    @test dos.data == 28.0
    @test dos.action.mu == -3.0
end
