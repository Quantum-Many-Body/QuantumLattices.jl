using LinearAlgebra: dot, tr
using QuantumLattices.Essentials: reset!, update
using QuantumLattices.Essentials.DegreesOfFreedom: Boundary, CompositeIndex, Coupling, Hilbert, IIDSpace, Index, OperatorUnitToTuple, SimpleIID, SimpleInternal, Subscript, Subscripts, Table, Term, @couplings
using QuantumLattices.Essentials.Frameworks
using QuantumLattices.Essentials.QuantumOperators: ID, Identity, Operator, Operators, id, idtype
using QuantumLattices.Essentials.Spatials: Lattice, Point, bonds, decompose, isintracell
using QuantumLattices.Interfaces: add!, expand, expand!
using QuantumLattices.Prerequisites: Float
using QuantumLattices.Prerequisites.Traits: contentnames, reparameter
using StaticArrays: SVector

import QuantumLattices.Essentials.Frameworks: Parameters
import QuantumLattices.Essentials: prepare!, run!, update!
import QuantumLattices.Prerequisites.VectorSpaces: shape

struct FID{N<:Union{Int, Symbol}} <: SimpleIID
    nambu::N
end
@inline Base.adjoint(sl::FID{Int}) = FID(3-sl.nambu)
function Base.angle(id::CompositeIndex{<:Index{FID{Int}}}, vectors::AbstractVector{<:AbstractVector{Float}}, values::AbstractVector{Float})
    phase=  (length(vectors) == 1) ? 2pi*dot(decompose(id.icoordinate, vectors[1]), values) :
            (length(vectors) == 2) ? 2pi*dot(decompose(id.icoordinate, vectors[1], vectors[2]), values) :
            (length(vectors) == 3) ? 2pi*dot(decompose(id.icoordinate, vectors[1], vectors[2], vectors[3]), values) :
            error("angle error: not supported number of input basis vectors.")
    (id.index.iid.nambu == 1) ? phase : -phase
end

struct FFock <: SimpleInternal{FID{Int}}
    nnambu::Int
end
@inline shape(f::FFock) = (1:f.nnambu,)
@inline FID(i::CartesianIndex, vs::FFock) = FID(i.I...)
@inline CartesianIndex(did::FID{Int}, vs::FFock) = CartesianIndex(did.nambu)
@inline shape(iidspace::IIDSpace{FID{Symbol}, FFock}) = (1:iidspace.internal.nnambu,)
@inline shape(iidspace::IIDSpace{FID{Int}, FFock}) = (iidspace.iid.nambu:iidspace.iid.nambu,)
@inline shape(iidspace::IIDSpace{FID{Symbol}, FFock}) = (1:iidspace.internal.nnambu,)
@inline shape(iidspace::IIDSpace{FID{Int}, FFock}) = (iidspace.iid.nambu:iidspace.iid.nambu,)

const FCoupling{V, I<:ID{FID}, C<:Subscripts} = Coupling{V, I, C}
@inline FCoupling(value, nambus::Tuple{Vararg{Int}}) = Coupling(value, ID(FID, nambus), Subscripts((nambu=Subscript(nambus),)))
@inline FCoupling(value, nambus::Subscript) = Coupling(value, ID(FID, convert(Tuple, nambus)), Subscripts((nambu=nambus,)))

@testset "Parameters" begin
    ps1 = Parameters{(:t1, :t2, :U)}(1.0im, 1.0, 2.0)
    ps2 = Parameters{(:t1, :U)}(1.0im, 2.0)
    ps3 = Parameters{(:t1, :U)}(1.0im, 2.1)
    @test match(ps1, ps2) == true
    @test match(ps1, ps3) == false
    @test update(ps1; ps3...) == Parameters{(:t1, :t2, :U)}(1.0im, 1.0, 2.1)

    bound = Boundary{(:θ₁, :θ₂)}([0.1, 0.2], [[1.0, 0.0], [0.0, 1.0]])
    @test Parameters(bound) == (θ₁=0.1, θ₂=0.2)
end

@testset "AnalyticalExpression" begin
    A(t, μ, Δ; k) = [
          2t*cos(k[1])+2t*cos(k[2])+μ   2im*Δ*sin(k[1])+2Δ*sin(k[2]);
        -2im*Δ*sin(k[1])+2Δ*sin(k[2])   -2t*cos(k[1])-2t*cos(k[2])-μ
    ]
    f = AnalyticalExpression(A, (t=1.0, μ=0.0, Δ=0.1))
    @test Parameters(f) == (t=1.0, μ=0.0, Δ=0.1)
    @test f(k=[0.0, 0.0]) ≈ [4 0; 0 -4]

    update!(f; μ=0.3)
    @test f(k=[pi/2, pi/2]) ≈ [0.3 0.2+0.2im; 0.2-0.2im -0.3]
end

@testset "RepresentationGenerator" begin
    @test contentnames(CompositeGenerator) == (:operators, :table)
    @test contentnames(OperatorGenerator) == (:operators, :terms, :bonds, :hilbert, :half, :table)
    @test contentnames(Image) == (:operators, :transformation, :table, :sourceid)

    lattice = Lattice([0.0], [0.5]; vectors=[[1.0]])
    bs = bonds(lattice, 1)
    hilbert = Hilbert(site->FFock(2), 1:length(lattice))
    table = Table(hilbert, OperatorUnitToTuple(:site))
    boundary = Boundary{(:θ,)}([0.1], lattice.vectors)

    t = Term{:Hp}(:t, 2.0, 1, couplings=@couplings(FCoupling(1.0, (2, 1))), ishermitian=false)
    μ = Term{:Mu}(:μ, 1.0, 0, couplings=@couplings(FCoupling(1.0, (2, 1))), ishermitian=true, modulate=true)

    tops₁ = expand(t, filter(bond->isintracell(bond), bs), hilbert; half=true)
    tops₂ = boundary(expand(one(t), filter(bond->!isintracell(bond), bs), hilbert; half=true))
    μops = expand(one(μ), filter(bond->length(bond)==1, bs), hilbert; half=true)
    i = Identity()

    optp = Operator{Complex{Float}, ID{CompositeIndex{Index{FID{Int}}, SVector{1, Float}}, 2}}
    entry = Entry(tops₁, (μ=μops,), (t=tops₂, μ=Operators{optp}()), (t=2.0, μ=1.0), boundary)
    @test entry == Entry((t, μ), bs, hilbert; half=true, boundary=boundary)
    @test isequal(entry, i(entry))
    @test Parameters(entry) == (t=2.0, μ=1.0, θ=0.1)
    @test entry+entry == entry*2 == 2*entry == Entry(tops₁*2, (μ=μops,), (t=tops₂, μ=Operators{optp}()), (t=4.0, μ=2.0), boundary)
    @test entry|>valtype == entry|>typeof|>valtype == Operators{optp, idtype(optp)}
    @test entry|>eltype == entry|>typeof|>eltype == optp
    @test expand(entry) == expand!(Operators{optp}(), entry) ≈ tops₁+tops₂*2.0+μops

    entry₁ = Entry(tops₁, NamedTuple(), (t=tops₂,), (t=2.0,), boundary)
    entry₂ = Entry(zero(tops₁), (μ=μops,), (μ=Operators{optp}(),), (μ=1.0,), boundary)
    @test entry₁+entry₂ == Entry(tops₁, (μ=μops,), (t=tops₂, μ=Operators{optp}()), (t=2.0, μ=1.0), boundary)

    μops₁ = expand(one(μ), bs[1], hilbert; half=true)
    μops₂ = expand(one(μ), bs[2], hilbert; half=true)
    entry₁ = Entry(zero(tops₁), (μ=μops₁,), (μ=Operators{optp}(),), (μ=1.0,), boundary)
    entry₂ = Entry(zero(tops₁), (μ=μops₂,), (μ=Operators{optp}(),), (μ=1.0,), boundary)
    @test entry₁*2+entry₂*2 == Entry(zero(tops₁), (μ=μops,), (μ=Operators{optp}(),), (μ=2.0,), boundary)
    @test entry₁*2+entry₂*3 == Entry(zero(tops₁), (μ=μops₁*2+μops₂*3,), (μ=Operators{optp}(),), (μ=1.0,), boundary)

    another = Entry(empty(tops₁), (μ=empty(μops),), (t=empty(tops₂), μ=Operators{optp}()), (t=2.0, μ=1.0), boundary)
    @test empty(entry) == empty!(deepcopy(entry)) == another

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

    @test reset!(deepcopy(another), (t, μ), bs, hilbert; half=true, boundary=boundary) == entry
    @test reset!(deepcopy(another), i, entry) == entry

    cgen = OperatorGenerator((t, μ), bs, hilbert; half=true, boundary=boundary, table=table)
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
    @test empty!(deepcopy(cgen)) == OperatorGenerator((t, μ), empty(bs), empty(hilbert); half=true, boundary=boundary, table=empty(table)) == empty(cgen)
    @test reset!(empty(cgen), lattice) == cgen
    @test update!(cgen, μ=1.5)|>expand ≈ tops₁+tops₂*2.0+μops*1.5

    sgen = i(cgen; table=table)
    @test sgen == Image(cgen.operators, i, table, objectid(cgen))
    @test Parameters(sgen) == (t=2.0, μ=1.5, θ=0.1)
    @test expand!(Operators{optp}(), sgen) == expand(sgen) ≈ tops₁+tops₂*2.0+μops*1.5
    @test empty!(deepcopy(sgen)) == Image(empty(cgen.operators), i, empty(table), objectid(cgen)) == empty(sgen)
    @test update!(sgen, μ=3.5)|>expand ≈ tops₁+tops₂*2.0+μops*3.5
    @test update!(sgen, cgen)|>expand ≈ tops₁+tops₂*2.0+μops*1.5
    @test reset!(empty(sgen), i, cgen; table=table) == sgen
end

mutable struct VCA <: Frontend
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
@inline run!(alg::Algorithm{VCA}, assign::Assignment{GF}) = (assign.action.count += 1; assign.data[:, :] .= alg.frontend.t+alg.frontend.U)

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
    @test Parameters(vca) == vca.parameters
    @test repr(vca, ≠(:U)) == "test(VCA)_1.0"
    @test repr(vca) == "test(VCA)_1.0_8.0"

    add!(vca, :GF, GF(0))
    rex = r"Action DOS\(DOS\)\: time consumed [0-9]*\.[0-9]*(e[+-][0-9]*)*s."
    @test_logs (:info, rex) vca(:DOS, DOS(-3.5), parameters=(U=7.0,), map=dosmap, dependences=(:GF,))

    dos = vca.assignments[:DOS]
    @test dos == deepcopy(dos)
    @test isequal(dos, deepcopy(dos))
    @test Parameters(dos) == dos.parameters
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
