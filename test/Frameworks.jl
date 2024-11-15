using LinearAlgebra: dot, tr
using QuantumLattices: add!, dimension, dtype, expand, expand!, matrix, reset!, update
using QuantumLattices.DegreesOfFreedom: plain, Boundary, CoordinatedIndex, Coupling, Hilbert, Index, OperatorUnitToTuple, SimpleInternalIndex, SimpleInternal, Term
using QuantumLattices.Frameworks
using QuantumLattices.QuantumOperators: ID, LinearFunction, Operator, Operators, id, idtype
using QuantumLattices.Spatials: Lattice, Point, bonds, decompose, isintracell
using QuantumLattices.Toolkit: Float, reparameter
using StaticArrays: SVector

import QuantumLattices: update!
import QuantumLattices.Frameworks: Parameters, initialize, run!
import QuantumLattices.Toolkit: shape

struct FID{N<:Union{Int, Symbol}} <: SimpleInternalIndex
    nambu::N
end
@inline Base.adjoint(sl::FID{Int}) = FID(3-sl.nambu)
function Base.angle(id::CoordinatedIndex{<:Index{FID{Int}}}, vectors::AbstractVector{<:AbstractVector{Float}}, values::AbstractVector{Float})
    phase=  (length(vectors) == 1) ? 2pi*dot(decompose(id.icoordinate, vectors[1]), values) :
            (length(vectors) == 2) ? 2pi*dot(decompose(id.icoordinate, vectors[1], vectors[2]), values) :
            (length(vectors) == 3) ? 2pi*dot(decompose(id.icoordinate, vectors[1], vectors[2], vectors[3]), values) :
            error("angle error: not supported number of input basis vectors.")
    (id.index.internal.nambu == 1) ? phase : -phase
end

struct FFock <: SimpleInternal{FID{Int}}
    nnambu::Int
end
@inline shape(f::FFock) = (1:f.nnambu,)
@inline Base.convert(::Type{<:FID}, i::CartesianIndex, ::FFock) = FID(i.I...)
@inline Base.convert(::Type{<:CartesianIndex}, fid::FID{Int}, ::FFock) = CartesianIndex(fid.nambu)
@inline shape(::FFock, index::FID{Int}) = (index.nambu:index.nambu,)
@inline shape(internal::FFock, ::FID{Symbol}) = (1:internal.nnambu,)

@testset "Parameters" begin
    ps1 = Parameters{(:t₁, :t₂, :U)}(1.0im, 1.0, 2.0)
    ps2 = Parameters{(:t₁, :U)}(1.0im, 2.0)
    ps3 = Parameters{(:t₁, :U)}(1.0im, 2.1)
    @test match(ps1, ps2) == true
    @test match(ps1, ps3) == false
    @test update(ps1; ps3...) == Parameters{(:t₁, :t₂, :U)}(1.0im, 1.0, 2.1)

    params = Parameters{(:t₁, :t₂)}(1.11111111, 2.2222222222)
    @test repr(params; context=:ndecimal=>2) == "(t₁ = 1.11, t₂ = 2.22)"

    bound = Boundary{(:θ₁, :θ₂)}([0.1, 0.2], [[1.0, 0.0], [0.0, 1.0]])
    @test Parameters(bound) == (θ₁=0.1, θ₂=0.2)

    ops = Operator(1, FID(1)) + Operator(2, FID(2))
    @test Parameters(ops) == NamedTuple()
end

@testset "Formula" begin
    A(t, μ, Δ; k=SVector(0, 0)) = [
        2t*cos(k[1]) + 2t*cos(k[2]) + μ   2im*Δ*sin(k[1]) + 2Δ*sin(k[2]);
        -2im*Δ*sin(k[1]) + 2Δ*sin(k[2])   -2t*cos(k[1]) - 2t*cos(k[2]) - μ
    ]::Matrix{ComplexF64}
    f = Formula(A, (t=1.0, μ=0.0, Δ=0.1))
    @test valtype(f) == Matrix{ComplexF64}
    @test eltype(f) == dtype(f) == ComplexF64
    @test Parameters(f) == (t=1.0, μ=0.0, Δ=0.1)
    @test f(; k=[0.0, 0.0]) ≈ [4 0; 0 -4]

    update!(f; μ=0.3)
    @test f(; k=[pi/2, pi/2]) ≈ [0.3 0.2+0.2im; 0.2-0.2im -0.3]
end

@testset "CategorizedGenerator twist" begin
    lattice = Lattice([0.0], [0.5]; vectors=[[1.0]])
    bs = bonds(lattice, 1)
    hilbert = Hilbert(site=>FFock(2) for site=1:length(lattice))
    boundary = Boundary{(:θ,)}([0.1], lattice.vectors)
    t = Term{:Hp}(:t, 2.0, 1, Coupling(1.0, :, FID, (2, 1)), false; ismodulatable=false)
    μ = Term{:Mu}(:μ, 1.0, 0, Coupling(1.0, :, FID, (2, 1)), true)
    i = LinearFunction(identity)

    optp = Operator{Complex{Float}, ID{CoordinatedIndex{Index{FID{Int}, Int}, SVector{1, Float}}, 2}}
    tops₁ = expand(t, filter(bond->isintracell(bond), bs), hilbert; half=true)
    tops₂ = boundary(expand(one(t), filter(bond->!isintracell(bond), bs), hilbert; half=true))
    μops = expand(one(μ), filter(bond->length(bond)==1, bs), hilbert; half=true)
    μops₁ = expand(one(μ), bs[1], hilbert; half=true)
    μops₂ = expand(one(μ), bs[2], hilbert; half=true)

    cat = CategorizedGenerator(tops₁, (t=Operators{optp}(), μ=μops), (t=tops₂, μ=Operators{optp}()), (t=2.0, μ=1.0), boundary, eager)
    @test cat == deepcopy(cat)
    @test isequal(cat, i(cat))
    @test cat|>valtype == cat|>typeof|>valtype == Operators{optp, idtype(optp)}
    @test cat|>eltype == cat|>typeof|>eltype == optp
    @test cat|>dtype == cat|>typeof|>dtype == Complex{Float}
    @test Parameters(cat) == (t=2.0, μ=1.0, θ=0.1)
    @test !isempty(cat) && isempty(empty(cat))
    @test empty(cat) == empty!(deepcopy(cat)) == CategorizedGenerator(Operators{optp}(), (t=Operators{optp}(), μ=Operators{optp}()), (t=Operators{optp}(), μ=Operators{optp}()), (t=2.0, μ=1.0), boundary, eager)
    @test cat * 2 == 2 * cat == cat + cat == CategorizedGenerator(tops₁*2, (t=Operators{optp}(), μ=μops), (t=tops₂, μ=Operators{optp}()), (t=4.0, μ=2.0), boundary, eager)
    @test expand(cat) == expand!(Operators{optp}(), cat) ≈ tops₁ + tops₂*2.0 + μops
    @test collect(cat) == collect(tops₁ + μops + tops₂*2.0)

    cat₁ = CategorizedGenerator(tops₁, (t=Operators{optp}(),), (t=tops₂,), (t=2.0,), boundary, eager)
    cat₂ = CategorizedGenerator(Operators{optp}(), (μ=μops,), (μ=Operators{optp}(),), (μ=1.0,), boundary, eager)
    @test cat₁+cat₂ == CategorizedGenerator(tops₁, (t=Operators{optp}(), μ=μops), (t=tops₂, μ=Operators{optp}()), (t=2.0, μ=1.0), boundary, eager)

    cat₁ = CategorizedGenerator(Operators{optp}(), (μ=μops₁,), (μ=Operators{optp}(),), (μ=1.0,), boundary, eager)
    cat₂ = CategorizedGenerator(Operators{optp}(), (μ=μops₂,), (μ=Operators{optp}(),), (μ=1.0,), boundary, eager)
    @test cat₁*2 + cat₂*2 == CategorizedGenerator(Operators{optp}(), (μ=μops,), (μ=Operators{optp}(),), (μ=2.0,), boundary, eager)
    @test cat₁*2 + cat₂*3 == CategorizedGenerator(Operators{optp}(), (μ=μops₁*2+μops₂*3,), (μ=Operators{optp}(),), (μ=1.0,), boundary, eager)

    nb = update!(deepcopy(boundary); θ=0.5)
    new = CategorizedGenerator(tops₁, (t=Operators{optp}(), μ=μops), (t=nb(tops₂, origin=[0.1]), μ=Operators{optp}()), (t=2.0, μ=1.5), nb, eager)
    dest = update!(deepcopy(cat); μ=1.5, θ=0.5)
    @test dest == new
    @test expand(dest) ≈ tops₁ + new.boundary(tops₂, origin=[0.1])*2.0 + μops*1.5
    dest = update!(deepcopy(cat), i, new)
    @test dest == new
    @test expand(dest) ≈ tops₁ + new.boundary(tops₂, origin=[0.1])*2.0 + μops*1.5
    @test reset!(deepcopy(new), i, cat) == cat
end

@testset "CategorizedGenerator plain" begin
    lattice = Lattice([0.0], [0.5]; vectors=[[1.0]])
    bs = bonds(lattice, 1)
    hilbert = Hilbert(site=>FFock(2) for site=1:length(lattice))
    t = Term{:Hp}(:t, 2.0, 1, Coupling(1.0, :, FID, (2, 1)), false; ismodulatable=false)
    μ = Term{:Mu}(:μ, 1.0, 0, Coupling(1.0, :, FID, (2, 1)), true)
    i = LinearFunction(identity)
    optp = Operator{Complex{Float}, ID{CoordinatedIndex{Index{FID{Int}, Int}, SVector{1, Float}}, 2}}
    tops = expand(t, bs, hilbert; half=true)
    μops = expand(one(μ), bs, hilbert; half=true)

    cat = CategorizedGenerator(tops, (t=Operators{optp}(), μ=μops), (t=Operators{optp}(), μ=Operators{optp}()), (t=2.0, μ=1.0), plain, eager)
    @test cat == deepcopy(cat)
    @test isequal(cat, i(cat))
    @test Parameters(cat) == (t=2.0, μ=1.0)
    @test cat*2 == 2*cat == cat+cat == CategorizedGenerator(tops*2, (t=Operators{optp}(), μ=μops), (t=Operators{optp}(), μ=Operators{optp}()), (t=4.0, μ=2.0), plain, eager)
    @test expand(cat) == expand!(Operators{optp}(), cat) ≈ tops+μops

    cat₁ = CategorizedGenerator(tops, (t=Operators{optp}(),), (t=Operators{optp}(),), (t=2.0,), plain, eager)
    cat₂ = CategorizedGenerator(Operators{optp}(), (μ=μops,), (μ=Operators{optp}(),), (μ=1.0,), plain, eager)
    @test cat₁+cat₂ == CategorizedGenerator(tops, (t=Operators{optp}(), μ=μops), (t=Operators{optp}(), μ=Operators{optp}()), (t=2.0, μ=1.0), plain, eager)

    new = CategorizedGenerator(Operators{optp}(), (t=Operators{optp}(), μ=Operators{optp}()), (t=Operators{optp}(), μ=Operators{optp}()), (t=2.0, μ=1.0), plain, eager)
    @test empty(cat) == empty!(deepcopy(cat)) == new
    dest = update!(deepcopy(cat); μ=1.5)
    @test expand(dest) ≈ tops + μops*1.5
    dest = update!(deepcopy(cat), i, dest)
    @test expand(dest) ≈ tops + μops*1.5
    @test reset!(deepcopy(new), i, cat) == cat
end

@testset "OperatorGenerator twist" begin
    lattice = Lattice([0.0], [0.5]; vectors=[[1.0]])
    bs = bonds(lattice, 1)
    hilbert = Hilbert(site=>FFock(2) for site=1:length(lattice))
    boundary = Boundary{(:θ,)}([0.1], lattice.vectors)
    t = Term{:Hp}(:t, 2.0, 1, Coupling(1.0, :, FID, (2, 1)), false; ismodulatable=false)
    μ = Term{:Mu}(:μ, 1.0, 0, Coupling(1.0, :, FID, (2, 1)), true)
    optp = Operator{Complex{Float}, ID{CoordinatedIndex{Index{FID{Int}, Int}, SVector{1, Float}}, 2}}
    tops₁ = expand(t, filter(bond->isintracell(bond), bs), hilbert; half=true)
    tops₂ = boundary(expand(one(t), filter(bond->!isintracell(bond), bs), hilbert; half=true))
    μops = expand(one(μ), filter(bond->length(bond)==1, bs), hilbert; half=true)
    cat = CategorizedGenerator(tops₁, (t=Operators{optp}(), μ=μops), (t=tops₂, μ=Operators{optp}()), (t=2.0, μ=1.0), boundary, eager)

    cgen = OperatorGenerator(cat, (t, μ), bs, hilbert, true)
    @test cgen == OperatorGenerator((t, μ), bs, hilbert, boundary; half=true)
    @test isequal(cgen, OperatorGenerator((t, μ), bs, hilbert, boundary; half=true))
    @test cgen|>eltype == cgen|>typeof|>eltype == optp
    @test cgen|>valtype == cgen|>typeof|>valtype == Operators{optp, idtype(optp)}
    @test collect(cgen) == collect(expand(cgen))
    @test Parameters(cgen) == (t=2.0, μ=1.0, θ=0.1)
    @test expand!(Operators{optp}(), cgen) == expand(cgen) ≈ tops₁ + tops₂*2.0 + μops
    @test expand(cgen, :t) ≈ tops₁ + tops₂*2.0
    @test expand(cgen, :μ) ≈ μops
    @test expand(cgen, 1) + expand(cgen, 2) + expand(cgen, 3) + expand(cgen, 4) ≈ expand(cgen)
    @test expand(cgen, :μ, 1) + expand(cgen, :μ, 2) ≈ μops
    @test expand(cgen, :t, 3) ≈ tops₁
    @test expand(cgen, :t, 4) ≈ tops₂*2.0
    @test empty!(deepcopy(cgen)) == OperatorGenerator((t, μ), empty(bs), empty(hilbert), boundary; half=true) == empty(cgen)
    @test !isempty(cgen) && isempty(empty(cgen)) 
    @test reset!(empty(cgen), lattice, hilbert) == cgen
    @test update!(cgen, μ=1.5)|>expand ≈ tops₁ + tops₂*2.0 + μops*1.5
    @test LinearFunction(identity)(cgen) == cgen.operators
end

@testset "OperatorGenerator plain" begin
    lattice = Lattice([0.0], [0.5]; vectors=[[1.0]])
    bs = bonds(lattice, 1)
    hilbert = Hilbert(site=>FFock(2) for site=1:length(lattice))
    t = Term{:Hp}(:t, 2.0, 1, Coupling(1.0, :, FID, (2, 1)), false; ismodulatable=false)
    μ = Term{:Mu}(:μ, 1.0, 0, Coupling(1.0, :, FID, (2, 1)), true)
    optp = Operator{Complex{Float}, ID{CoordinatedIndex{Index{FID{Int}, Int}, SVector{1, Float}}, 2}}
    tops = expand(t, bs, hilbert; half=true)
    μops = expand(one(μ), bs, hilbert; half=true)
    cat = CategorizedGenerator(tops, (t=Operators{optp}(), μ=μops), (t=Operators{optp}(), μ=Operators{optp}()), (t=2.0, μ=1.0), plain, eager)

    cgen = OperatorGenerator(cat, (t, μ), bs, hilbert, true)
    @test cgen == OperatorGenerator((t, μ), bs, hilbert, plain; half=true)
    @test expand(cgen) ≈ tops + μops
    @test expand(cgen, :t) ≈ tops
    @test expand(cgen, :μ) ≈ μops
    @test expand(cgen, 1) + expand(cgen, 2) + expand(cgen, 3) + expand(cgen, 4) ≈ expand(cgen)
    @test expand(cgen, :μ, 1) + expand(cgen, :μ, 2) ≈ μops
    @test expand(cgen, :t, 3) + expand(cgen, :t, 4) ≈ tops
    @test reset!(empty(cgen), lattice, hilbert) == cgen
    @test update!(cgen, μ=1.5)|>expand ≈ tops + μops*1.5
    @test LinearFunction(identity)(cgen) == cgen.operators
end

@testset "Hamiltonian" begin
    A(t, μ, Δ; k=SVector(0, 0)) = [
        2t*cos(k[1]) + 2t*cos(k[2]) + μ   2im*Δ*sin(k[1]) + 2Δ*sin(k[2]);
        -2im*Δ*sin(k[1]) + 2Δ*sin(k[2])   -2t*cos(k[1]) - 2t*cos(k[2]) - μ
    ]::Matrix{ComplexF64}
    h = Hamiltonian(A, (t=1.0, μ=0.0, Δ=0.1))
    @test h == Hamiltonian(Formula(A, (t=1, μ=0, Δ=0.1)))
    @test isequal(h, Hamiltonian(Formula(A, (t=1, μ=0, Δ=0.1))))
    @test h.representation == Formula(A, (t=1.0, μ=0.0, Δ=0.1))
    @test valtype(h) == Matrix{ComplexF64}
    @test eltype(h) == dtype(h) == ComplexF64
    @test Parameters(h) == (t=1.0, μ=0.0, Δ=0.1)
    @test update!(h; μ=0.3) == Hamiltonian(Formula(A, (t=1.0, μ=0.3, Δ=0.1)))
    @test matrix(h; k=[pi/2, pi/2]) ≈ [0.3 0.2+0.2im; 0.2-0.2im -0.3]
    @test dimension(h) == 2

    lattice = Lattice([0.0], [0.5]; vectors=[[1.0]])
    hilbert = Hilbert(site=>FFock(2) for site=1:length(lattice))
    t = Term{:Hp}(:t, 2.0, 1, Coupling(1.0, :, FID, (2, 1)), false; ismodulatable=false)
    μ = Term{:Mu}(:μ, 1.0, 0, Coupling(1.0, :, FID, (2, 1)), true)
    cat = OperatorGenerator((t, μ), bonds(lattice, 1), hilbert, plain, eager; half=true)
    h = Hamiltonian(cat)
    @test valtype(h) == valtype(typeof(h)) == valtype(cat)
    @test eltype(h) == eltype(typeof(h)) == eltype(cat)
    @test dtype(h) == dtype(typeof(h)) == dtype(cat)
    @test Parameters(h) == Parameters(cat)
    @test update!(h; μ=0.2) == Hamiltonian(cat)

    ops = expand(cat)
    h = Hamiltonian(ops)
    @test valtype(h) == valtype(typeof(h)) == typeof(ops)
    @test eltype(h) == eltype(typeof(h)) == eltype(ops)
    @test dtype(h) == dtype(typeof(h)) == dtype(ops)
    @test Parameters(h) == NamedTuple()
    @test update!(h; μ=0.8) == Hamiltonian(ops)
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
@inline initialize(gf::GF, vca::VCA) = Matrix{Float}(undef, vca.dim, vca.dim)
@inline run!(alg::Algorithm{VCA}, assign::Assignment{GF}) = (assign.action.count += 1; assign.data[:, :] .= alg.frontend.t+alg.frontend.U)

mutable struct DOS <: Action
    μ::Float
end
@inline update!(eb::DOS; kwargs...) = (eb.μ = get(kwargs, :μ, eb.μ); eb)
@inline initialize(dos::DOS, vca::VCA) = 0.0
@inline function run!(alg::Algorithm{VCA}, assign::Assignment{DOS})
    prepare!(alg, assign)
    gf = alg.assignments[first(assign.dependencies)]
    assign.data = tr(gf.data)
end
@inline dosmap(params::Parameters) = Parameters{(:t, :U, :μ)}(params.t, params.U, -params.U/2)

@testset "Framework" begin
    vca = VCA(1.0, 8.0, 4)
    @test vca == deepcopy(vca)
    @test isequal(vca, deepcopy(vca))
    @test string(vca) == "VCA"

    gf = GF(0)
    @test gf == deepcopy(gf)
    @test isequal(gf, deepcopy(gf))
    @test update!(gf) == gf

    vca = Algorithm(:test, vca)
    @test vca == deepcopy(vca)
    @test isequal(vca, deepcopy(vca))
    @test Parameters(vca) == vca.parameters
    @test repr(vca; context=:select=>≠(:U)) == "test(VCA)-t(1.0)"
    @test repr(vca) == "test(VCA)-t(1.0)U(8.0)"

    add!(vca, :GF, GF(0))
    rex = r"Action DOS\(DOS\)\: time consumed [0-9]*\.[0-9]*(e[+-][0-9]*)*s."
    @test_logs (:info, rex) vca(:DOS, DOS(-3.5), parameters=(U=7.0,), map=dosmap, dependencies=(:GF,))

    dos = vca.assignments[:DOS]
    @test dos == deepcopy(dos)
    @test isequal(dos, deepcopy(dos))
    @test Parameters(dos) == dos.parameters
    @test valtype(dos) == valtype(typeof(dos)) == Float
    @test dos.data == 32.0
    @test dos.action.μ == -3.5
    @test nameof(vca, dos) == "test(VCA)-t(1.0)U(7.0)-DOS"

    update!(dos, U=6.0)
    vca(:DOS, info=false)
    @test dos.parameters == (t=1.0, U=6.0)
    @test dos.data == 28.0
    @test dos.action.μ == -3.0

    save("path.dat", 1:100, rand(100, 2))
    save("heatmap.dat", 1:100, 1:50, rand(50, 100))
end
