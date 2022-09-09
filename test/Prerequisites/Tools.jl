using QuantumLattices.Prerequisites

@testset "allequal" begin
    @test allequal((2, 2, 2, 2))
    @test !allequal((2, 2, 1, 2))
end

@testset "decimaltostr" begin
    @test decimaltostr(:a) == ":a"
    @test decimaltostr(1) == "1"
    @test decimaltostr(10^6) == "1000000"
    @test decimaltostr(1//7) == "1//7"
    @test decimaltostr(1.0) == "1.0"
    @test decimaltostr(10.0^6) == "1.0e+06"
    @test decimaltostr(10^-5) == "1.0e-05"
    @test decimaltostr(1/7) == "0.14286"
    @test decimaltostr(1/7, 8) == "0.14285714"
    @test decimaltostr(0im) == "0.0"
    @test decimaltostr(1//7+0im) == "1//7"
    @test decimaltostr(1.0im) == "1.0im"
    @test decimaltostr(-1.0im) == "-1.0im"
    @test decimaltostr(0.1+0.12im) == "0.1+0.12im"
    @test decimaltostr(0.1-0.12im) == "0.1-0.12im"
end

@testset "ordinal" begin
    @test ordinal(1) == "1st"
    @test ordinal(2) == "2nd"
    @test ordinal(3) == "3rd"
    @test ordinal(4) == "4th"
    @test ordinal(5) == "5th"
end

@testset "delta" begin
    @test delta(1, 2) == 0
    @test delta(1, 1) == 1
end

@testset "concatenate" begin
    @test concatenate((1, 2), (:a, :b)) == (1, 2, :a, :b)
end

@testset "searchsortedfirst" begin
    idx = CartesianIndices((2, 2, 2))
    for (i, index) in enumerate(idx)
        @test searchsortedfirst(idx, index) == i
    end
    @test searchsortedfirst(idx, CartesianIndex(0, 0, 0)) == 1
    @test searchsortedfirst(idx, CartesianIndex(3, 3, 3)) == 9
end