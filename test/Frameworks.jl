using LinearAlgebra: dot, eigen
using Plots: plot, savefig
using QuantumLattices: expand, expand!, reset!, str, update
using QuantumLattices.DegreesOfFreedom: plain, Boundary, CoordinatedIndex, Coupling, Hilbert, Index, SimpleInternalIndex, SimpleInternal, Term
using QuantumLattices.Frameworks
using QuantumLattices.Frameworks: seriestype
using QuantumLattices.QuantumNumbers: periods
using QuantumLattices.QuantumOperators: ID, LinearFunction, Operator, Operators, idtype, scalartype
using QuantumLattices.Spatials: BrillouinZone, Lattice, bonds, decompose, isintracell
using StaticArrays: SVector, SMatrix, @SMatrix

import QuantumLattices: update!
import QuantumLattices.Frameworks: Parameters, options, run!
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
    @test str(params; ndecimal=2) == "t₁(1.11)t₂(2.22)"
    @test str(params; ndecimal=2, select=isequal(:t₁)) == "t₁(1.11)"
    @test str(params; ndecimal=2, select=isequal(:t₁), front="SC") == "SC-t₁(1.11)"
    @test str(params; ndecimal=2, select=isequal(:t₁), rear="[0.0, 0.0]") == "t₁(1.11)-[0.0, 0.0]"
    @test str(params; ndecimal=2, select=isequal(:t₁), front="SC", rear="[0.0, 0.0]") == "SC-t₁(1.11)-[0.0, 0.0]"

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
    @test f == Formula{SMatrix{2, 2, ComplexF64, 4}}(A, (t=1.0, μ=0.0, Δ=0.1))
    @test isequal(f, Formula{SMatrix{2, 2, ComplexF64, 4}}(A, (t=1.0, μ=0.0, Δ=0.1)))
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
Base.show(io::IO, eigensystem::EigenSystem) = print(io, "EigenSystem(", join(periods(keytype(eigensystem.brillouinzone)), "×"), ")")
struct EigenSystemData <: Data
    values::Vector{Vector{Float64}}
    vectors::Vector{Matrix{ComplexF64}}
end
@inline options(::Type{<:Assignment{<:EigenSystem}}) = (
    showinfo = "show the information",
)
function run!(tba::Algorithm{<:TBA}, eigensystem::Assignment{<:EigenSystem}; options...)
    get(options, :showinfo, false) && @info string(eigensystem)
    data = EigenSystemData(Vector{Float64}[], Matrix{ComplexF64}[])
    for k in eigensystem.action.brillouinzone
        values, vectors = eigen(tba.frontend.formula(k))
        push!(data.values, values)
        push!(data.vectors, vectors)
    end
    return data
end

struct DensityOfStates <:Action end
mutable struct DensityOfStatesData <: Data
    energies::Vector{Float64}
    values::Vector{Float64}
end
@inline options(::Type{<:Assignment{<:DensityOfStates}}) = (
    emin = "lower bound of the energy range",
    emax = "upper bound of the energy range",
    ne = "number of sample points in the energy range",
    σ = "broadening factor"
)
function run!(tba::Algorithm{<:TBA}, dos::Assignment{<:DensityOfStates}; emin=nothing, emax=nothing, ne::Int=101, σ=0.1)
    @assert isa(dos.dependencies, Tuple{Assignment{<:EigenSystem}}) "run! error: wrong dependencies."
    eigensystem = first(dos.dependencies)
    isnothing(emin) && (emin = mapreduce(minimum, min, eigensystem.data.values))
    isnothing(emax) && (emax = mapreduce(maximum, max, eigensystem.data.values))
    data = DensityOfStatesData(range(emin, emax, ne), zeros(ne))
    for (i, ω) in enumerate(data.energies)
        data.values[i] = 0.0
        for energies in eigensystem.data.values
            for e in energies
                data.values[i] += exp(-(ω-e)^2/2σ^2)
            end
        end
        data.values[i] /= √(2pi)*σ
    end
    return data
end

A(t, μ, k=SVector(0.0, 0.0); kwargs...) = SMatrix{1, 1}(2t*cos(k[1])+2t*cos(k[2])+μ)

@testset "Frontend & Action & Data" begin
    tba = TBA(Formula(A, (t=1.0, μ=0.5)))
    @test tba==deepcopy(tba) && isequal(tba, deepcopy(tba))

    eigensystem = EigenSystem(BrillouinZone([[2pi, 0], [0, 2pi]], 100))
    @test eigensystem==deepcopy(eigensystem) && isequal(eigensystem, deepcopy(eigensystem))
    @test update!(eigensystem; Parameters(tba)...) == eigensystem

    eigensystemdata = EigenSystemData(Vector{Float64}[], Matrix{ComplexF64}[])
    @test eigensystemdata==deepcopy(eigensystemdata) && isequal(eigensystemdata, deepcopy(eigensystemdata))
    @test Tuple(eigensystemdata) == (Vector{Float64}[], Matrix{ComplexF64}[])
end

params(parameters::Parameters) = (t=parameters.t, μ=parameters.U/2)

@testset "Assignment & Algorithm with map" begin
    tba = Algorithm(:Square, TBA(Formula(A, (t=1.0, μ=0.5))), (t=1.0, U=2.0), params)
    @test tba==deepcopy(tba) && isequal(tba, deepcopy(tba))
    update!(tba; U=1.0)
    @test Parameters(tba) == (t=1.0, U=1.0)
    @test string(tba) == "Square-TBA"
    @test dirname(tba) == "."
    @test basename(tba) == "Square-TBA.qld"
    @test basename(tba; prefix="Prefix") == "Prefix-Square-TBA.qld"
    @test basename(tba; suffix="Suffix") == "Square-TBA-Suffix.qld"
    @test basename(tba; extension="dat") == "Square-TBA.dat"
    @test pathof(tba) == joinpath(dirname(tba), basename(tba))
    @test str(tba) == "Square-TBA-t(1.0)U(1.0)"
    io = IOBuffer()
    show(io, MIME"text/plain"(), tba)
    @test String(take!(io)) == "Square\n  frontend:\n    TBA\n  parameters:\n    t: 1.0\n    U: 1.0"

    @test options(Assignment) == NamedTuple()

    eigensystem = tba(:eigensystem, EigenSystem(BrillouinZone([[2pi, 0], [0, 2pi]], 100)); delay=true)
    @test Parameters(eigensystem) == (t=1.0, U=1.0)
    @test valtype(eigensystem) == valtype(typeof(eigensystem)) == EigenSystemData
    update!(eigensystem; U=2.0)
    @test Parameters(eigensystem) == (t=1.0, U=2.0)
    @test string(eigensystem) == "eigensystem"
    io = IOBuffer()
    show(io, MIME"text/plain"(), eigensystem)
    @test String(take!(io)) == "eigensystem\n  action:\n    EigenSystem(100×100)\n  parameters:\n    t: 1.0\n    U: 2.0"
    @test options(typeof(eigensystem)) == (showinfo="show the information",)
    @test optionsinfo(typeof(eigensystem)) == "Assignment{<:EigenSystem} options:\n  (1) `:showinfo`: show the information.\n"

    dos = tba(:DOS, DensityOfStates(), (t=1.0, U=4.0), (eigensystem,))
    @test dos==deepcopy(dos) && isequal(dos, deepcopy(dos))
    @test options(typeof(dos)) == (emin="lower bound of the energy range", emax="upper bound of the energy range", ne ="number of sample points in the energy range", σ="broadening factor")
    @test optionsinfo(typeof(dos)) == "Assignment{<:DensityOfStates} options:\n  (1) `:emin`: lower bound of the energy range;\n  (2) `:emax`: upper bound of the energy range;\n  (3) `:ne`: number of sample points in the energy range;\n  (4) `:σ`: broadening factor.\n\n  Dependency 1) Assignment{<:EigenSystem} options:\n    (1) `:showinfo`: show the information.\n"
    @test hasoption(typeof(dos), :emin) && hasoption(typeof(dos), :emax) && hasoption(typeof(dos), :ne) && hasoption(typeof(dos), :σ) && hasoption(typeof(dos), :showinfo) && !hasoption(typeof(dos), :hello)
    @test sum(dos.data.values)*(maximum(dos.data.energies)-minimum(dos.data.energies))/(length(dos.data.energies)-1)/length(eigensystem.action.brillouinzone) ≈ 0.9964676726997486
    savefig(plot(dos), "$(string(dos)).png")
    update!(dos; U=8.0)
    tba(dos)
    savefig(plot(tba(dos)), "$(string(dos)).png")
    update!(dos; U=0.0)
    tba(dos; emin=-5.0, emax=5.0)
    savefig(plot(tba(dos)), "$(string(dos)).png")
    summary(tba)
end

@testset "Assignment & Algorithm without map" begin
    tba = Algorithm(:Square, TBA(Formula(A, (t=1.0, μ=2.0))))
    qldsave(tba)
    loaded = qldload(pathof(tba), str(Parameters(tba)))
    @test loaded == qldload(pathof(tba))[str(Parameters(tba))]
    @test all(isequal(loaded), qldload(pathof(tba), str(Parameters(tba)), str(Parameters(tba))))

    eigensystem = loaded(:eigensystem, EigenSystem(BrillouinZone([[2pi, 0], [0, 2pi]], 100)); delay=true)
    dos = loaded(:DOS, DensityOfStates(), (eigensystem,))
    savefig(plot(loaded(dos)), "$(string(dos)).png")

    @test isnothing(seriestype())
    @test seriestype(dos.data) == seriestype(dos.data.energies, dos.data.values) == :path
    @test seriestype(zeros(0), zeros(0), zeros(0, 0)) == :heatmap
end
