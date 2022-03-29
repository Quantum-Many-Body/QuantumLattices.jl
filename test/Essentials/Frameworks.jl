using Test
using Printf: @printf
using LinearAlgebra: tr, dot
using StaticArrays: SVector
using QuantumLattices.Essentials.Frameworks
using QuantumLattices.Essentials: update, reset!
using QuantumLattices.Essentials.Spatials: Point, AbstractPID, PID, Bond, Bonds, Lattice, acrossbonds, zerothbonds, decompose
using QuantumLattices.Essentials.DegreesOfFreedom: SimpleIID, SimpleInternal, IIDSpace, Coupling, Subscript, Subscripts
using QuantumLattices.Essentials.DegreesOfFreedom: Term, Hilbert, Index, Table, OID, OIDToTuple, @couplings
using QuantumLattices.Essentials.QuantumOperators: ID, Operator, Operators, id, idtype, Identity
using QuantumLattices.Prerequisites: Float
using QuantumLattices.Interfaces:  expand!, expand, add!
using QuantumLattices.Prerequisites.Traits: contentnames, reparameter

import QuantumLattices.Essentials.Frameworks: Parameters
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

const FCoupling{V, I<:ID{FID}, C<:Subscripts} = Coupling{V, I, C}
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
    @test update(ps1; ps3...) == Parameters{(:t1, :t2, :U)}(1.0im, 1.0, 2.1)
end

@testset "Boundary" begin
    op = Operator(4.5,
        OID(Index(PID(1), FID(2)), SVector(0.5, 0.5), SVector(0.0, 0.0)),
        OID(Index(PID(2), FID(1)), SVector(1.5, 1.5), SVector(1.0, 1.0))
        )
    bound = Boundary{(:θ₁, :θ₂)}([0.1, 0.2], [[1.0, 0.0], [0.0, 1.0]])
    M = reparameter(typeof(op), :value, Complex{Float64})
    @test valtype(typeof(bound), typeof(op)) == M
    @test keys(bound) == keys(typeof(bound)) == (:θ₁, :θ₂)
    @test Parameters(bound) == (θ₁=0.1, θ₂=0.2)

    another = Boundary{(:θ₁, :θ₂)}([0.0, 0.0], [[2.0, 0.0], [0.0, 2.0]])
    @test merge!(deepcopy(bound), another) == another

    @test bound(op) ≈ replace(op, 4.5*exp(2im*pi*0.3))
    @test bound(op, origin=[0.05, 0.15]) ≈ replace(op, 4.5*exp(2im*pi*0.1))
    update!(bound, θ₁=0.3)
    @test bound(op) ≈ replace(op, 4.5*exp(2im*pi*0.5))
    @test bound(op, origin=[0.1, 0.1]) ≈ replace(op, 4.5*exp(2im*pi*0.3))

    ops = Operators{M}(op)
    @test valtype(typeof(bound), typeof(ops)) == typeof(ops)
    @test bound(ops) ≈ Operators(replace(op, 4.5*exp(2im*pi*0.5)))
    @test map!(bound, ops) ≈ ops ≈ Operators(replace(op, 4.5*exp(2im*pi*0.5)))

    @test valtype(typeof(plain), typeof(op)) == typeof(op)
    @test valtype(typeof(plain), typeof(ops)) == typeof(ops)
    @test plain(op) == op
    @test plain(ops) == ops
    @test update!(plain) == plain
end

@testset "Formulation" begin
    A(t, μ, Δ; k) = [
          2t*cos(k[1])+2t*cos(k[2])+μ   2im*Δ*sin(k[1])+2Δ*sin(k[2]);
        -2im*Δ*sin(k[1])+2Δ*sin(k[2])   -2t*cos(k[1])-2t*cos(k[2])-μ
    ]
    f = Formulation(A, (t=1.0, μ=0.0, Δ=0.1))
    @test Parameters(f) == (t=1.0, μ=0.0, Δ=0.1)
    @test f(k=[0.0, 0.0]) ≈ [4 0; 0 -4]

    update!(f; μ=0.3)
    @test f(k=[pi/2, pi/2]) ≈ [0.3 0.2+0.2im; 0.2-0.2im -0.3]
end

@testset "Generator" begin
    @test contentnames(CompositeGenerator) == (:operators, :table)
    @test contentnames(Generator) == (:operators, :terms, :bonds, :hilbert, :half, :table)
    @test contentnames(Image) == (:operators, :transformation, :table, :sourceid)

    lattice = Lattice(:Tuanzi, [Point(PID(1), [0.0]), Point(PID(2), [0.5])], vectors=[[1.0]], neighbors=1)
    bonds = Bonds(lattice)
    hilbert = Hilbert{FFock}(pid->FFock(2), lattice.pids)
    table = Table(hilbert, OIDToTuple(:scope, :site))
    boundary = Boundary{(:θ,)}([0.1], lattice.vectors)

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
    @test entry|>eltype == entry|>typeof|>eltype == optp
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
    @test sgen == Image(cgen.operators, i, table, objectid(cgen))
    @test Parameters(sgen) == (t=2.0, μ=1.5, θ=0.1)
    @test expand!(Operators{optp}(), sgen) == expand(sgen) ≈ tops₁+tops₂*2.0+μops*1.5
    @test empty!(deepcopy(sgen)) == Image(empty(cgen.operators), i, empty(table), objectid(cgen)) == empty(sgen)
    @test update!(sgen, μ=3.5)|>expand ≈ tops₁+tops₂*2.0+μops*3.5
    @test update!(sgen, cgen)|>expand ≈ tops₁+tops₂*2.0+μops*1.5
    @test reset!(empty(sgen), cgen; table=table) == sgen
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
