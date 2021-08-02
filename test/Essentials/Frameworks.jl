using Test
using Printf: @printf
using LinearAlgebra: tr
using QuantumLattices.Essentials: register!
using QuantumLattices.Essentials.Frameworks
using QuantumLattices.Essentials.Terms: Parameters
using QuantumLattices.Interfaces: id
using QuantumLattices.Prerequisites: Float

import QuantumLattices.Essentials: prepare!, run!, update!
import QuantumLattices.Essentials.Frameworks: dependences

struct FEngine <: Engine end
Base.repr(::FEngine) = "FEngine"
Base.show(io::IO, ::FEngine) = @printf io "FEngine"
struct FApp <: App end

@testset "Assignment" begin
    @test FApp() == FApp()
    @test isequal(FApp(), FApp())
    @test update!(FApp()) == FApp()

    assign = Assignment(:FApp, FApp(), (t=1.0, U=8.0), dependences=(:FApp₁, :FApp₂))
    @test deepcopy(assign) == assign
    @test isequal(deepcopy(assign), assign)
    @test assign|>valtype == assign|>typeof|>valtype == Any
    @test assign|>id == assign|>typeof|>id == :FApp
    update!(assign, t=2.0)
    @test assign.parameters == (t=2.0, U=8.0)
end

@testset "Algorithm" begin
    @test FEngine() == FEngine()
    @test isequal(FEngine(), FEngine())
    @test update!(FEngine()) == FEngine()

    alg = Algorithm("Alg", FEngine(), parameters=(t=1.0, U=8.0))
    @test string(alg) == repr(alg) == "Alg_FEngine_1.0_8.0"
    @test repr(alg, (:U,)) == "Alg_FEngine_1.0"
    @test_logs (:info, r"App FApp₁\(FApp\)\: time consumed [0-9]*\.[0-9]*s.") register!(alg, :FApp₁, FApp(), parameters=(U=5.0,))
    @test_logs (:info, r"App FApp₂\(FApp\)\: time consumed [0-9]*\.[0-9]*s.") register!(alg, :FApp₂, FApp(), parameters=(U=6.0,), dependences=(:FApp₁,))
    @test get(alg, Val(:FApp₁)) == get(alg, :FApp₁) == Assignment(:FApp₁, FApp(), (t=1.0, U=5.0), virgin=false)
    @test get(alg, Val(:FApp₂)) == get(alg, :FApp₂) == Assignment(:FApp₂, FApp(), (t=1.0, U=6.0), dependences=(:FApp₁,), virgin=false)
    @test isnothing(prepare!(alg, get(alg, :FApp₁)))
    @test isnothing(run!(alg, get(alg, :FApp₁)))
    @test dependences(alg, get(alg, :FApp₁)) == ()
    @test dependences(alg, get(alg, :FApp₂)) == (:FApp₁,)
    @test dependences(alg, get(alg, :FApp₂), (:FApp₁,)) == ()
end

mutable struct VCA <: Engine
    t::Float
    U::Float
end
function update!(vca::VCA; kwargs...)
    vca.t = get(kwargs, :t, vca.t)
    vca.U = get(kwargs, :U, vca.U)
    return vca
end
Parameters(vca::VCA) = Parameters{(:t, :U)}(vca.t, vca.U)
gf(alg::Algorithm{VCA}) = get(alg, Val(:_VCAGF_)).data

mutable struct GF <: App
    dim::Int
    count::Int
end

mutable struct DOS <: App
    mu::Float
end
update!(eb::DOS; kwargs...) = (eb.mu = get(kwargs, :mu, eb.mu); eb)

dependences(alg::Algorithm{VCA}, assign::Assignment{GF}, ::Tuple{}=()) = assign.dependences
dependences(alg::Algorithm{VCA}, assign::Assignment, ::Tuple{}=()) = (:_VCAGF_, assign.dependences...)
prepare!(alg::Algorithm{VCA}, assign::Assignment{GF}) = assign.virgin && (assign.data = Matrix{valtype(assign)|>eltype}(undef, assign.app.dim, assign.app.dim))
run!(alg::Algorithm{VCA}, assign::Assignment{GF}) = (assign.app.count += 1; assign.data[:, :] .= alg.engine.t+alg.engine.U)

function run!(alg::Algorithm{VCA}, assign::Assignment{DOS})
    rundependences!(alg, assign)
    assign.data = tr(gf(alg))
end

@testset "Framework" begin
    engine = VCA(1.0, 8.0)
    gf = Assignment(:_VCAGF_, GF(4, 0), Parameters(engine), data=Matrix{Complex{Float}}(undef, 0, 0))
    vca = Algorithm("Test", engine, assignments=(gf,))
    @test_logs (:info, r"App DOS\(DOS\)\: time consumed [0-9]*\.[0-9]*s.") register!(vca, :DOS, DOS(-3.5), parameters=(U=7.0,), map=params::Parameters->Parameters{(:t, :U, :mu)}(params.t, params.U, -params.U/2))
    dos = get(vca, :DOS)
    @test dos.data == 32.0
    @test dos.app.mu == -3.5
    update!(dos, U=6.0)
    run!(vca, :DOS, false)
    @test dos.data == 28.0
    @test dos.app.mu == -3.0
end
