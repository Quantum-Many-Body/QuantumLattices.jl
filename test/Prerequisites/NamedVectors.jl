using Test
using QuantumLattices.Prerequisites.NamedVectors
using QuantumLattices.Prerequisites: Float

mutable struct NHNV <: NamedVector
    scope::String
    site::Int
end

@testset "NHNV" begin
    @test fieldnames(NHNV) == (:scope, :site)
    @test length(NHNV) == 2

    pid = NHNV("A", 0)
    @test pid == convert(NHNV, ("A", 0))
    @test string(pid) == "NHNV(\"A\", 0)"
    @test length(pid) == 2
    @test collect(Iterators.reverse(pid)) == [0, "A"]
    @test pid[1] == pid.scope == "A"
    @test pid[2] == pid.site == 0
    @test replace(pid, scope="B") == NHNV("B", 0)
    @test replace(pid, site=1) == NHNV("A", 1)
    @test replace(pid, scope="B", site=1) == NHNV("B", 1)
    @test (pid[1] = "B"; pid[1] == "B")
    @test (pid.site = 2; pid.site == 2)

    @test NHNV("A", 2.0) < NHNV("B", 0.0)
    @test NHNV("A", 2.0) < NHNV("A", 3.0)
    @test isless(NHNV("A", 2.0), NHNV("B", 0.0))
    @test isless(NHNV("A", 2.0), NHNV("A", 3.0))

    dict = Dict(NHNV("A", 0)=>1, NHNV("A", 1)=>2)
    @test dict[NHNV("A", 0)] == 1
    @test dict[NHNV("A", 1)] == 2
end

mutable struct FHNV <: HomoNamedVector{Float64}
    scope::Float64
    site::Float64
end

mutable struct RHNV{T<:Real} <: HomoNamedVector{T}
    scope::T
    site::T
end

@testset "FHNV" begin
    @test fieldnames(FHNV) == (:scope, :site)
    @test length(FHNV) == 2
    @test eltype(FHNV) == Float64
    @test zero(FHNV) == FHNV(0.0, 0.0)
    @test isequal(zero(FHNV), FHNV(0.0, 0.0))

    pid = FHNV(1.0, 0.0)
    @test pid == convert(FHNV, (1.0, 0.0))
    @test string(pid) == "FHNV(1.0, 0.0)"
    @test length(pid) == 2
    @test eltype(pid) == Float64
    @test zero(pid) == FHNV(0.0, 0.0)
    @test pid[1] == pid.scope == 1.0
    @test pid[2] == pid.site == 0.0
    @test replace(pid, scope=2.0) == FHNV(2.0, 0.0)
    @test replace(pid, site=1.0) == FHNV(1.0, 1.0)
    @test replace(pid, scope=2.0, site=1.0) == FHNV(2.0, 1.0)
    @test (pid[1] = 2.0; pid[1] == 2.0)
    @test (pid.site = 3.0; pid.site == 3.0)

    @test FHNV(1.0, 2.0) < FHNV(2.0, 0.0)
    @test FHNV(1.0, 2.0) < FHNV(1.0, 3.0)
    @test isless(FHNV(1.0, 2.0), FHNV(2.0, 0.0))
    @test isless(FHNV(1.0, 2.0), FHNV(1.0, 3.0))

    dict = Dict(FHNV(0.0, 0.0)=>1, FHNV(0.0, 1.0)=>2)
    @test dict[FHNV(0.0, 0.0)] == 1
    @test dict[FHNV(0.0, 1.0)] == 2
end

@testset "RHNV" begin
    @test fieldnames(RHNV) == (:scope, :site)
    @test length(RHNV) == 2
    @test eltype(RHNV{Int}) == Int
    @test zero(RHNV{Int}) == RHNV(0, 0)
    @test eltype(RHNV{Float64}) == Float64
    @test zero(RHNV{Float64}) == RHNV(0.0, 0.0)

    pid = RHNV(1, 0)
    @test pid == convert(RHNV{Int}, (1, 0))
    @test string(pid) == "RHNV(1, 0)"
    @test length(pid) == 2
    @test eltype(pid) == Int
    @test zero(pid) == RHNV(0, 0)
    @test pid[1] == pid.scope == 1
    @test pid[2] == pid.site == 0
    @test replace(pid, scope=2) == FHNV(2, 0)
    @test replace(pid, site=1) == FHNV(1, 1)
    @test replace(pid, scope=2, site=1) == FHNV(2, 1)
    @test (pid[1] = 2; pid[1] == 2)
    @test (pid.site = 3; pid.site == 3)

    @test RHNV(1.0, 2.0) < RHNV(2.0, 0.0)
    @test RHNV(1.0, 2.0) < RHNV(1.0, 3.0)

    dict = Dict(RHNV(0, 0)=>1, RHNV(0, 1)=>2)
    @test dict[RHNV(0, 0)] == 1
    @test dict[RHNV(0, 1)] == 2
end
