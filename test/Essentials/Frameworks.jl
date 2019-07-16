using Test
using Printf: @printf
using LinearAlgebra: tr
using QuantumLattices.Essentials.Frameworks
using QuantumLattices.Essentials.Terms: Parameters
using QuantumLattices.Interfaces: id,register!
using QuantumLattices.Prerequisites: Float
import QuantumLattices.Interfaces: prepare!,run!,update!
import QuantumLattices.Essentials.Frameworks: dependences

struct FEngine <: Engine end
Base.repr(::FEngine)="FEngine"
Base.show(io::IO,::FEngine)=@printf io "FEngine"
struct FApp <: App end

@testset "Assignment" begin
    @test FApp()==FApp()
    @test isequal(FApp(),FApp())
    @test update!(FApp())==FApp()

    assign=Assignment(:FApp,FApp(),(t=1.0,U=8.0),dependences=(:FApp1,:FApp2))
    @test deepcopy(assign)==assign
    @test isequal(deepcopy(assign),assign)
    @test assign|>valtype==assign|>typeof|>valtype==Any
    @test assign|>id==assign|>typeof|>id==:FApp
    update!(assign,t=2.0)
    @test assign.parameters==(t=2.0,U=8.0)
end

@testset "Algorithm" begin
    @test FEngine()==FEngine()
    @test isequal(FEngine(),FEngine())
    @test update!(FEngine())==FEngine()

    alg=Algorithm("Alg",FEngine(),parameters=(t=1.0,U=8.0))
    @test string(alg)==repr(alg)=="Alg_FEngine_1.0_8.0"
    @test repr(alg,(:U,))=="Alg_FEngine_1.0"
    register!(alg,:FApp1,FApp(),parameters=(U=5.0,))
    register!(alg,:FApp2,FApp(),parameters=(U=6.0,),dependences=(:FApp1,))
    @test get(alg,Val(:FApp1))==get(alg,:FApp1)==Assignment(:FApp1,FApp(),(t=1.0,U=5.0),virgin=false)
    @test get(alg,Val(:FApp2))==get(alg,:FApp2)==Assignment(:FApp2,FApp(),(t=1.0,U=6.0),dependences=(:FApp1,),virgin=false)
    @test prepare!(alg,get(alg,:FApp1))===nothing
    @test run!(alg,get(alg,:FApp1))===nothing
    @test dependences(alg,get(alg,:FApp1))==()
    @test dependences(alg,get(alg,:FApp2))==(:FApp1,)
    @test dependences(alg,get(alg,:FApp2),(:FApp1,))==()
end

mutable struct VCA <: Engine
    t::Float
    U::Float
end
function update!(vca::VCA;kwargs...)
    vca.t=get(kwargs,:t,vca.t)
    vca.U=get(kwargs,:U,vca.U)
    return vca
end
Parameters(vca::VCA)=Parameters{(:t,:U)}(vca.t,vca.U)
gf(alg::Algorithm{VCA})=get(alg,Val(:_VCAGF_)).data

mutable struct GF <: App
    dim::Int
    count::Int
end

mutable struct DOS <: App
    mu::Float
end
update!(eb::DOS;kwargs...)=(eb.mu=get(kwargs,:mu,eb.mu);eb)

dependences(alg::Algorithm{VCA},assign::Assignment{GF},::Tuple{}=())=assign.dependences
dependences(alg::Algorithm{VCA},assign::Assignment,::Tuple{}=())=(:_VCAGF_,assign.dependences...)
prepare!(alg::Algorithm{VCA},assign::Assignment{GF})=assign.virgin && (assign.data=Matrix{valtype(assign)|>eltype}(undef,assign.app.dim,assign.app.dim))
run!(alg::Algorithm{VCA},assign::Assignment{GF})=(assign.app.count+=1; assign.data[:,:].=alg.engine.t+alg.engine.U)

function run!(alg::Algorithm{VCA},assign::Assignment{DOS})
    rundependences!(alg,assign)
    assign.data=tr(gf(alg))
end

@testset "Framework" begin
    engine=VCA(1.0,8.0)
    gf=Assignment(:_VCAGF_,GF(4,0),Parameters(engine),data=Matrix{Complex{Float}}(undef,0,0))
    vca=Algorithm("Test",engine,assignments=(gf,))
    register!(vca,:DOS,DOS(-3.5),parameters=(U=7.0,),map=params::Parameters->Parameters{(:t,:U,:mu)}(params.t,params.U,-params.U/2))
    dos=get(vca,:DOS)
    @test dos.data==32.0
    @test dos.app.mu==-3.5
    update!(dos,U=6.0)
    run!(vca,:DOS,false)
    @test dos.data==28.0
    @test dos.app.mu==-3.0
end
