using Test
using QuantumLattices.Mathematics.Combinatorics

@testset "Combinations" begin
    @test Combinations{0}("abc")|>collect == [()]
    @test Combinations{1}("abc")|>collect == [('a',), ('b',), ('c',)]
    @test Combinations{2}("abc")|>collect == [('a', 'b'), ('a', 'c'), ('b', 'c')]
    @test Combinations{3}("abc")|>collect == [('a', 'b', 'c')]
end

@testset "DulCombinations" begin
    @test DulCombinations{0}("abc")|>collect == [()]
    @test DulCombinations{1}("abc")|>collect == [('a',), ('b',), ('c',)]
    @test DulCombinations{2}("abc")|>collect == [('a', 'a'), ('a', 'b'), ('a', 'c'), ('b', 'b'), ('b', 'c'), ('c', 'c')]
    @test DulCombinations{3}("abc")|>collect == [('a', 'a', 'a'), ('a', 'a', 'b'), ('a', 'a', 'c'), ('a', 'b', 'b'), ('a', 'b', 'c'), ('a', 'c', 'c'), ('b', 'b', 'b'), ('b', 'b', 'c'), ('b', 'c', 'c'), ('c', 'c', 'c')]
end

@testset "Permutations" begin
    @test Permutations{0}("abc")|>collect == [()]
    @test Permutations{1}("abc")|>collect == [('a',), ('b',), ('c',)]
    @test Permutations{2}("abc")|>collect == [('a', 'b'), ('a', 'c'), ('b', 'a'), ('b', 'c'), ('c', 'a'), ('c', 'b')]
    @test Permutations{3}("abc")|>collect == [('a', 'b', 'c'), ('a', 'c', 'b'), ('b', 'a', 'c'), ('b', 'c', 'a'), ('c', 'a', 'b'), ('c', 'b', 'a')]
end

@testset "DulPermutations" begin
    @test DulPermutations{0}("abc")|>collect == [()]
    @test DulPermutations{1}("abc")|>collect == [('a',), ('b',), ('c',)]
    @test DulPermutations{2}("abc")|>collect == [('a', 'a'), ('a', 'b'), ('a', 'c'), ('b', 'a'), ('b', 'b'), ('b', 'c'), ('c', 'a'), ('c', 'b'), ('c', 'c')]
    @test DulPermutations{3}("abc")|>collect == [(i, j, k) for i∈"abc" for j∈"abc" for k∈"abc"]
end
