using LinearAlgebra: dot, eigen
using Plots: plot, savefig
using QuantumLattices: expand, expand!, reset!, update
using QuantumLattices.DegreesOfFreedom: plain, Boundary, CoordinatedIndex, Coupling, Hilbert, Index, SimpleInternalIndex, SimpleInternal, Term
using QuantumLattices.Frameworks
using QuantumLattices.Frameworks: seriestype
using QuantumLattices.QuantumOperators: ID, LinearFunction, Operator, Operators, idtype, scalartype
using QuantumLattices.Spatials: BrillouinZone, Lattice, bonds, decompose, isintracell, save
using StaticArrays: SVector, SMatrix, @SMatrix

import QuantumLattices: dimension, update!
import QuantumLattices.Frameworks: Parameters, run!
import QuantumLattices.Toolkit: shape

struct FID{N<:Union{Int, Symbol}} <: SimpleInternalIndex
    nambu::N
end
@inline Base.adjoint(sl::FID{Int}) = FID(3-sl.nambu)
function Base.angle(id::CoordinatedIndex{<:Index{FID{Int}}}, vectors::AbstractVector{<:AbstractVector{Float64}}, values::AbstractVector{Float64})
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
    @test Parameters(ops) == Parameters()
end

@testset "Formula" begin
    A(t, μ, Δ, k=SVector(0, 0)) = @SMatrix([
        2t*cos(k[1]) + 2t*cos(k[2]) + μ   2im*Δ*sin(k[1]) + 2Δ*sin(k[2]);
        -2im*Δ*sin(k[1]) + 2Δ*sin(k[2])   -2t*cos(k[1]) - 2t*cos(k[2]) - μ
    ])
    f = Formula(A, (t=1.0, μ=0.0, Δ=0.1))
    @test valtype(f) == valtype(typeof(f)) == SMatrix{2, 2, ComplexF64, 4}
    @test scalartype(f) == scalartype(typeof(f)) == ComplexF64
    @test Parameters(f) == (t=1.0, μ=0.0, Δ=0.1)
    @test f([0.0, 0.0]) ≈ [4 0; 0 -4]

    update!(f; μ=0.3)
    @test f([pi/2, pi/2]) ≈ [0.3 0.2+0.2im; 0.2-0.2im -0.3]
end

@testset "CategorizedGenerator twist" begin
    lattice = Lattice([0.0], [0.5]; vectors=[[1.0]])
    bs = bonds(lattice, 1)
    hilbert = Hilbert(site=>FFock(2) for site=1:length(lattice))
    boundary = Boundary{(:θ,)}([0.1], lattice.vectors)
    t = Term{:Hp}(:t, 2.0, 1, Coupling(1.0, :, FID, (2, 1)), false; ismodulatable=false)
    μ = Term{:Mu}(:μ, 1.0, 0, Coupling(1.0, :, FID, (2, 1)), true)
    i = LinearFunction(identity)

    optp = Operator{ComplexF64, ID{CoordinatedIndex{Index{FID{Int}, Int}, SVector{1, Float64}}, 2}}
    tops₁ = expand(t, filter(bond->isintracell(bond), bs), hilbert; half=true)
    tops₂ = boundary(expand(one(t), filter(bond->!isintracell(bond), bs), hilbert; half=true))
    μops = expand(one(μ), filter(bond->length(bond)==1, bs), hilbert; half=true)
    μops₁ = expand(one(μ), bs[1], hilbert; half=true)
    μops₂ = expand(one(μ), bs[2], hilbert; half=true)

    cat = CategorizedGenerator(tops₁, (t=Operators{optp}(), μ=μops), (t=tops₂, μ=Operators{optp}()), (t=2.0, μ=1.0), boundary, eager)
    @test cat == deepcopy(cat)
    @test isequal(cat, i(cat))
    @test cat == Generator(tops₁, (t=Operators{optp}(), μ=μops), (t=tops₂, μ=Operators{optp}()), (t=2.0, μ=1.0), boundary, eager)
    @test cat|>valtype == cat|>typeof|>valtype == Operators{optp, idtype(optp)}
    @test cat|>eltype == cat|>typeof|>eltype == optp
    @test cat|>scalartype == cat|>typeof|>scalartype == ComplexF64
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
    optp = Operator{ComplexF64, ID{CoordinatedIndex{Index{FID{Int}, Int}, SVector{1, Float64}}, 2}}
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
    optp = Operator{ComplexF64, ID{CoordinatedIndex{Index{FID{Int}, Int}, SVector{1, Float64}}, 2}}
    tops₁ = expand(t, filter(bond->isintracell(bond), bs), hilbert; half=true)
    tops₂ = boundary(expand(one(t), filter(bond->!isintracell(bond), bs), hilbert; half=true))
    μops = expand(one(μ), filter(bond->length(bond)==1, bs), hilbert; half=true)
    cat = CategorizedGenerator(tops₁, (t=Operators{optp}(), μ=μops), (t=tops₂, μ=Operators{optp}()), (t=2.0, μ=1.0), boundary, eager)

    cgen = OperatorGenerator(cat, bs, hilbert, (t, μ), true)
    @test cgen == OperatorGenerator(bs, hilbert, (t, μ), boundary; half=true)
    @test isequal(cgen, OperatorGenerator(bs, hilbert, (t, μ), boundary; half=true))
    @test cgen == Generator(cat, bs, hilbert, (t, μ), true) == Generator(bs, hilbert, (t, μ), boundary; half=true)
    @test cgen|>valtype == cgen|>typeof|>valtype == Operators{optp, idtype(optp)}
    @test cgen|>eltype == cgen|>typeof|>eltype == optp
    @test cgen|>scalartype == cgen|>typeof|>scalartype == ComplexF64
    @test collect(cgen) == collect(expand(cgen))
    @test Parameters(cgen) == (t=2.0, μ=1.0, θ=0.1)
    @test expand!(Operators{optp}(), cgen) == expand(cgen) ≈ tops₁ + tops₂*2.0 + μops
    @test expand(cgen, :t) ≈ tops₁ + tops₂*2.0
    @test expand(cgen, :μ) ≈ μops
    @test expand(cgen, 1) + expand(cgen, 2) + expand(cgen, 3) + expand(cgen, 4) ≈ expand(cgen)
    @test expand(cgen, :μ, 1) + expand(cgen, :μ, 2) ≈ μops
    @test expand(cgen, :t, 3) ≈ tops₁
    @test expand(cgen, :t, 4) ≈ tops₂*2.0
    @test empty!(deepcopy(cgen)) == OperatorGenerator(empty(bs), empty(hilbert), (t, μ), boundary; half=true) == empty(cgen)
    @test !isempty(cgen) && isempty(empty(cgen)) 
    @test reset!(empty(cgen), bs, hilbert; vectors=lattice.vectors) == cgen
    @test update!(deepcopy(cgen), μ=1.5)|>expand ≈ tops₁ + tops₂*2.0 + μops*1.5
    @test LinearFunction(identity)(cgen) == cgen.operators
    @test reset!(empty(cat), LinearFunction(identity), cgen) == cat
    @test update!(deepcopy(cat), LinearFunction(identity), update!(deepcopy(cgen), μ=1.5))|>expand ≈ tops₁ + tops₂*2.0 + μops*1.5
end

@testset "OperatorGenerator plain" begin
    lattice = Lattice([0.0], [0.5]; vectors=[[1.0]])
    bs = bonds(lattice, 1)
    hilbert = Hilbert(site=>FFock(2) for site=1:length(lattice))
    t = Term{:Hp}(:t, 2.0, 1, Coupling(1.0, :, FID, (2, 1)), false; ismodulatable=false)
    μ = Term{:Mu}(:μ, 1.0, 0, Coupling(1.0, :, FID, (2, 1)), true)
    optp = Operator{ComplexF64, ID{CoordinatedIndex{Index{FID{Int}, Int}, SVector{1, Float64}}, 2}}
    tops = expand(t, bs, hilbert; half=true)
    μops = expand(one(μ), bs, hilbert; half=true)
    cat = CategorizedGenerator(tops, (t=Operators{optp}(), μ=μops), (t=Operators{optp}(), μ=Operators{optp}()), (t=2.0, μ=1.0), plain, eager)

    cgen = OperatorGenerator(cat, bs, hilbert, (t, μ), true)
    @test cgen == OperatorGenerator(bs, hilbert, (t, μ), plain; half=true)
    @test expand(cgen) ≈ tops + μops
    @test expand(cgen, :t) ≈ tops
    @test expand(cgen, :μ) ≈ μops
    @test expand(cgen, 1) + expand(cgen, 2) + expand(cgen, 3) + expand(cgen, 4) ≈ expand(cgen)
    @test expand(cgen, :μ, 1) + expand(cgen, :μ, 2) ≈ μops
    @test expand(cgen, :t, 3) + expand(cgen, :t, 4) ≈ tops
    @test reset!(empty(cgen), bs, hilbert; vectors=lattice.vectors) == cgen
    @test update!(deepcopy(cgen), μ=1.5)|>expand ≈ tops + μops*1.5
    @test LinearFunction(identity)(cgen) == cgen.operators
    @test reset!(empty(cat), LinearFunction(identity), cgen) == cat
    @test update!(deepcopy(cat), LinearFunction(identity), update!(deepcopy(cgen), μ=1.5))|>expand ≈ tops + μops*1.5
end

struct TBA{F<:Formula} <: Frontend
    formula::F
end
@inline Base.show(io::IO, ::TBA) = print(io, "TBA")
@inline Parameters(tba::TBA) = tba.formula.parameters
@inline update!(tba::TBA; parameters...) = (update!(tba.formula; parameters...); tba)

struct EigenSystem{B<:BrillouinZone} <: Action
    brillouinzone::B
end
struct EigenSystemData <: Data
    values::Vector{Vector{Float64}}
    vectors::Vector{Matrix{ComplexF64}}
end
@inline Data(eigensystem::EigenSystem, tba::TBA) = EigenSystemData(Vector{Float64}[], Matrix{ComplexF64}[])
function run!(tba::Algorithm{<:TBA}, eigensystem::Assignment{<:EigenSystem})
    data = eigensystem.data
    empty!(data.values)
    empty!(data.vectors)
    for k in eigensystem.action.brillouinzone
        values, vectors = eigen(tba.frontend.formula(k))
        push!(data.values, values)
        push!(data.vectors, vectors)
    end
end

struct DensityOfStates <:Action
    ne::Int
    σ::Float64
end
struct DensityOfStatesData <: Data
    energies::Vector{Float64}
    values::Vector{Float64}
end
@inline Data(dos::DensityOfStates, ::TBA) = DensityOfStatesData(zeros(dos.ne), zeros(dos.ne))
function run!(tba::Algorithm{<:TBA}, dos::Assignment{<:DensityOfStates})
    @assert isa(dos.dependencies, Tuple{Assignment{<:EigenSystem}}) "run! error: wrong dependencies."
    eigensystem = tba(dos.dependencies[1])
    emin = mapreduce(minimum, min, eigensystem.data.values)
    emax = mapreduce(maximum, max, eigensystem.data.values)
    for (i, ω) in enumerate(range(emin, emax, dos.action.ne))
        dos.data.energies[i] = ω
        dos.data.values[i] = 0.0
        for energies in eigensystem.data.values
            for e in energies
                dos.data.values[i] += exp(-(ω-e)^2/2dos.action.σ^2)
            end
        end
        dos.data.values[i] /= √(2pi)*dos.action.σ
    end
end

@testset "Framework" begin
    A(t, μ, k=SVector(0.0, 0.0); kwargs...) = SMatrix{1, 1}(2t*cos(k[1])+2t*cos(k[2])+μ)
    tba = TBA(Formula(A, (t=1.0, μ=0.5)))
    @test tba==deepcopy(tba) && isequal(tba, deepcopy(tba))

    eigensystem = EigenSystem(BrillouinZone([[2pi, 0], [0, 2pi]], 100))
    @test eigensystem==deepcopy(eigensystem) && isequal(eigensystem, deepcopy(eigensystem))
    @test update!(eigensystem; Parameters(tba)...) == eigensystem
    @test options(typeof(eigensystem)) == Dict{Symbol, String}()
    @test checkoptions(typeof(eigensystem))

    eigensystemdata = Data(eigensystem, tba)
    @test eigensystemdata==deepcopy(eigensystemdata) && isequal(eigensystemdata, deepcopy(eigensystemdata))
    @test Tuple(eigensystemdata) == (Vector{Float64}[], Matrix{ComplexF64}[])

    params(parameters::Parameters) = (t=parameters.t, μ=parameters.U/2)
    eigensystem = Assignment(:eigensystem, eigensystem, (t=1.0, U=2.0), params, (), eigensystemdata, false)
    @test eigensystem==deepcopy(eigensystem) && isequal(eigensystem, deepcopy(eigensystem))
    @test Parameters(eigensystem) == (t=1.0, U=2.0)
    @test valtype(eigensystem) == valtype(typeof(eigensystem)) == EigenSystemData
    update!(eigensystem; U=1.0)
    @test Parameters(eigensystem) == (t=1.0, U=1.0)

    tba = Algorithm(:Square, tba, (t=1.0, U=2.0), params)
    @test tba==deepcopy(tba) && isequal(tba, deepcopy(tba))
    update!(tba; U=1.0)
    @test Parameters(tba) == (t=1.0, U=1.0)
    @test string(tba) == "Square\n  frontend:\n    TBA\n  parameters:\n    t: 1.0\n    U: 1.0"

    dos = tba(:DOS, DensityOfStates(101, 0.1), (t=1.0, U=4.0), params, (eigensystem,))
    @test sum(dos.data.values)*(maximum(dos.data.energies)-minimum(dos.data.energies))/(dos.action.ne-1)/length(eigensystem.action.brillouinzone) ≈ 0.9964676726997486
    summary(dos)
    savefig(plot(dos), "DensityOfStatesU4.png")

    update!(dos; U=8.0)
    savefig(plot(tba(dos)), "DensityOfStatesU8.png")

    @test isnothing(seriestype())
    @test seriestype(dos.data) == seriestype(dos.data.energies, dos.data.values) == :path
    @test seriestype(zeros(0), zeros(0), zeros(0, 0)) == :heatmap

#     save("path.dat", 1:100, rand(100, 2))
#     save("heatmap.dat", 1:100, 1:50, rand(50, 100))
end
