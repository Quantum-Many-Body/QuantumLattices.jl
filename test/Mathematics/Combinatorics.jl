using Test
using QuantumLattices.Mathematics.Combinatorics

@testset "Combinations" begin
    @test Combinations{0}(1:3)|>collect==[()]
    @test Combinations{1}(1:3)|>collect==[(1,),(2,),(3,)]
    @test Combinations{2}(1:3)|>collect==[(1,2),(1,3),(2,3)]
    @test Combinations{3}(1:3)|>collect==[(1,2,3)]
end

@testset "DulCombinations" begin
    @test DulCombinations{0}(1:3)|>collect==[()]
    @test DulCombinations{1}(1:3)|>collect==[(1,),(2,),(3,)]
    @test DulCombinations{2}(1:3)|>collect==[(1,1),(1,2),(1,3),(2,2),(2,3),(3,3)]
    @test DulCombinations{3}(1:3)|>collect==[(1,1,1),(1,1,2),(1,1,3),(1,2,2),(1,2,3),(1,3,3),(2,2,2),(2,2,3),(2,3,3),(3,3,3)]
end

@testset "Permutations" begin
    @test Permutations{0}(1:3)|>collect==[()]
    @test Permutations{1}(1:3)|>collect==[(1,),(2,),(3,)]
    @test Permutations{2}(1:3)|>collect==[(1,2),(1,3),(2,1),(2,3),(3,1),(3,2)]
    @test Permutations{3}(1:3)|>collect==[(1,2,3),(1,3,2),(2,1,3),(2,3,1),(3,1,2),(3,2,1)]
end

@testset "DulPermutations" begin
    @test DulPermutations{0}(1:3)|>collect==[()]
    @test DulPermutations{1}(1:3)|>collect==[(1,),(2,),(3,)]
    @test DulPermutations{2}(1:3)|>collect==[(1,1),(1,2),(1,3),(2,1),(2,2),(2,3),(3,1),(3,2),(3,3)]
    @test DulPermutations{3}(1:3)|>collect==[(i,j,k) for i=1:3 for j=1:3 for k=1:3]
end
