using QuantumLattices.Prerequisites: Float
using QuantumLattices.Prerequisites.NamedVectors

mutable struct NonHomo <: NamedVector
    scope::String
    site::Int
end

@testset "NonHomo" begin
    @test length(NonHomo) == 2

    pid = NonHomo("A", 0)
    @test pid == convert(NonHomo, ("A", 0))
    @test string(pid) == "NonHomo(\"A\", 0)"
    @test length(pid) == 2
    @test collect(Iterators.reverse(pid)) == [0, "A"]
    @test pid[1] == pid.scope == "A"
    @test pid[2] == pid.site == 0
    @test replace(pid, scope="B") == NonHomo("B", 0)
    @test replace(pid, site=1) == NonHomo("A", 1)
    @test replace(pid, scope="B", site=1) == NonHomo("B", 1)
    @test (pid[1] = "B"; pid[1] == "B")
    @test (pid.site = 2; pid.site == 2)

    @test NonHomo("A", 2.0) < NonHomo("B", 0.0)
    @test NonHomo("A", 2.0) < NonHomo("A", 3.0)
    @test isless(NonHomo("A", 2.0), NonHomo("B", 0.0))
    @test isless(NonHomo("A", 2.0), NonHomo("A", 3.0))

    dict = Dict(NonHomo("A", 0)=>1, NonHomo("A", 1)=>2)
    @test dict[NonHomo("A", 0)] == 1
    @test dict[NonHomo("A", 1)] == 2
end

mutable struct Homo{T<:Real} <: HomoNamedVector{T}
    scope::T
    site::T
end

@testset "Homo" begin
    @test length(Homo) == 2
    @test eltype(Homo{Int}) == Int
    @test zero(Homo{Int}) == Homo(0, 0)
    @test eltype(Homo{Float64}) == Float64
    @test zero(Homo{Float64}) == Homo(0.0, 0.0)

    pid = Homo(1, 0)
    @test pid == convert(Homo{Int}, (1, 0))
    @test string(pid) == "Homo(1, 0)"
    @test length(pid) == 2
    @test eltype(pid) == Int
    @test zero(pid) == Homo(0, 0)
    @test pid[1] == pid.scope == 1
    @test pid[2] == pid.site == 0
    @test replace(pid, scope=2) == Homo(2, 0)
    @test replace(pid, site=1) == Homo(1, 1)
    @test replace(pid, scope=2, site=1) == Homo(2, 1)
    @test (pid[1] = 2; pid[1] == 2)
    @test (pid.site = 3; pid.site == 3)

    @test Homo(1.0, 2.0) < Homo(2.0, 0.0)
    @test Homo(1.0, 2.0) < Homo(1.0, 3.0)

    dict = Dict(Homo(0, 0)=>1, Homo(0, 1)=>2)
    @test dict[Homo(0, 0)] == 1
    @test dict[Homo(0, 1)] == 2
end
