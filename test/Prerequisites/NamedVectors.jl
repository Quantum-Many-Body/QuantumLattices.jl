using Test
using QuantumLattices.Prerequisites.NamedVectors
using QuantumLattices.Prerequisites: Float

@namedvector mutable struct NHNV
    scope::String
    site::Int
end

@testset "NHNV" begin
    @test NHNV|>fieldnames==(:scope,:site)
    @test NHNV|>length==2

    pid=NHNV("A",0)
    @test pid==convert(NHNV,("A",0))
    @test pid|>string=="NHNV(\"A\",0)"
    @test pid|>length==2
    @test pid|>Iterators.reverse|>collect==[0,"A"]
    @test pid[1]==pid.scope=="A"
    @test pid[2]==pid.site==0
    @test replace(pid,scope="B")==NHNV("B",0)
    @test replace(pid,site=1)==NHNV("A",1)
    @test replace(pid,scope="B",site=1)==NHNV("B",1)
    @test (pid[1]="B";pid[1]=="B")
    @test (pid.site=2;pid.site==2)

    @test NHNV("A",2.0)<NHNV("B",0.0)
    @test NHNV("A",2.0)<NHNV("A",3.0)
    @test isless(NHNV("A",2.0),NHNV("B",0.0))
    @test isless(NHNV("A",2.0),NHNV("A",3.0))

    dict=Dict(NHNV("A",0)=>1,NHNV("A",1)=>2)
    @test dict[NHNV("A",0)]==1
    @test dict[NHNV("A",1)]==2
end

@homonamedvector "FHNV" (:scope,:site) Float64 mutable=true
@homonamedvector "RHNV" (:scope,:site) (<:Real) mutable=true

@testset "FHNV" begin
    @test FHNV|>fieldnames==(:scope,:site)
    @test FHNV|>length==2
    @test FHNV|>eltype==Float64
    @test FHNV|>zero==FHNV(0.0,0.0)
    @test isequal(FHNV|>zero,FHNV(0.0,0.0))

    pid=FHNV(1.0,0.0)
    @test pid==convert(FHNV,(1.0,0.0))
    @test pid|>string=="FHNV(1.0,0.0)"
    @test pid|>length==2
    @test pid|>eltype==Float64
    @test pid|>zero==FHNV(0.0,0.0)
    @test pid[1]==pid.scope==1.0
    @test pid[2]==pid.site==0.0
    @test replace(pid,scope=2.0)==FHNV(2.0,0.0)
    @test replace(pid,site=1.0)==FHNV(1.0,1.0)
    @test replace(pid,scope=2.0,site=1.0)==FHNV(2.0,1.0)
    @test (pid[1]=2.0;pid[1]==2.0)
    @test (pid.site=3.0;pid.site==3.0)

    @test FHNV(1.0,2.0)<FHNV(2.0,0.0)
    @test FHNV(1.0,2.0)<FHNV(1.0,3.0)
    @test isless(FHNV(1.0,2.0),FHNV(2.0,0.0))
    @test isless(FHNV(1.0,2.0),FHNV(1.0,3.0))

    dict=Dict(FHNV(0.0,0.0)=>1,FHNV(0.0,1.0)=>2)
    @test dict[FHNV(0.0,0.0)]==1
    @test dict[FHNV(0.0,1.0)]==2
end

@testset "RHNV" begin
    @test RHNV|>fieldnames==(:scope,:site)
    @test RHNV|>length==2
    @test RHNV{Int}|>eltype==Int
    @test RHNV{Int}|>zero==RHNV(0,0)
    @test RHNV{Float64}|>eltype==Float64
    @test RHNV{Float64}|>zero==RHNV(0.0,0.0)

    pid=RHNV(1,0)
    @test pid==convert(RHNV{Int},(1,0))
    @test pid|>string=="RHNV(1,0)"
    @test pid|>length==2
    @test pid|>eltype==Int
    @test pid|>zero==RHNV(0,0)
    @test pid[1]==pid.scope==1
    @test pid[2]==pid.site==0
    @test replace(pid,scope=2)==FHNV(2,0)
    @test replace(pid,site=1)==FHNV(1,1)
    @test replace(pid,scope=2,site=1)==FHNV(2,1)
    @test (pid[1]=2;pid[1]==2)
    @test (pid.site=3;pid.site==3)

    @test RHNV(1.0,2.0)<RHNV(2.0,0.0)
    @test RHNV(1.0,2.0)<RHNV(1.0,3.0)

    dict=Dict(RHNV(0,0)=>1,RHNV(0,1)=>2)
    @test dict[RHNV(0,0)]==1
    @test dict[RHNV(0,1)]==2
end
