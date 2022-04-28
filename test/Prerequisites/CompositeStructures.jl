using Test
using QuantumLattices.Prerequisites.CompositeStructures

import QuantumLattices.Prerequisites.Traits: contentnames, getcontent

struct CT{S, N, T} <: CompositeNTuple{N, T}
    info::S
    contents::NTuple{N, T}
end
contentnames(::Type{<:CT}) = (:info, :contents)

@testset "CompositeTuple" begin
    @test contentnames(CompositeTuple) == (:contents,)

    t = CT("Info", (1, 2, 3, 4))
    @test contentnames(typeof(t)) == (:info, :contents)
    @test length(t) == 4
    @test length(typeof(t)) == 4
    @test eltype(t) == Int
    @test eltype(typeof(t)) == Int
    @test t == deepcopy(t)
    @test isequal(t, deepcopy(t))
    @test hash(t) == hash(("Info", (1, 2, 3, 4)))
    @test t[1] == 1
    @test t[end] == 4
    @test t[1:3] == CT("Info", (1, 2, 3))
    @test collect(t) == [1, 2, 3, 4]
    @test collect(t|>Iterators.reverse) == [4, 3, 2, 1]
    @test keys(t) == keys((1, 2, 3, 4))
    @test values(t) == values((1, 2, 3, 4))
    @test pairs(t)|>collect == pairs((1, 2, 3, 4))|>collect
    @test reverse(t) == CT("Info", (4, 3, 2, 1))
    @test Tuple(t) == (1, 2, 3, 4)
end

struct CV{S, T} <: CompositeVector{T}
    info::S
    contents::Vector{T}
end
contentnames(::Type{<:CV}) = (:info, :contents)

@testset "CompositeVector" begin
    @test contentnames(CompositeVector) == (:contents,)

    v = CV("Info", [1, 3, 2, 4])
    @test contentnames(typeof(v)) == (:info, :contents)
    @test size(v) == (4,)
    @test size(v, 1) == 4
    @test length(v) == 4
    @test v == deepcopy(v)
    @test isequal(v, deepcopy(v))
    @test v[end] == 4
    @test v[2] == 3
    @test v[CartesianIndex(2)] == 3
    @test v[[1, 3, 4]] == CV("Info", [1, 2, 4])
    @test v[1:3] == CV("Info", [1, 3, 2])
    @test (v[1] = 10; v[2:3] = 11:12; v == CV("Info", [10, 11, 12, 4]))
    @test push!(v, 1) == CV("Info", [10, 11, 12, 4, 1])
    @test push!(v, 2, 3) == CV("Info", [10, 11, 12, 4, 1, 2, 3])
    @test pushfirst!(v, 20, 21) == CV("Info", [20, 21, 10, 11, 12, 4, 1, 2, 3])
    @test insert!(v, 2, 30) == CV("Info", [20, 30, 21, 10, 11, 12, 4, 1, 2, 3])
    @test append!(v, [0, -1]) == CV("Info", [20, 30, 21, 10, 11, 12, 4, 1, 2, 3, 0, -1])
    @test prepend!(v, [8, 9]) == CV("Info", [8, 9, 20, 30, 21, 10, 11, 12, 4, 1, 2, 3, 0, -1])
    @test (splice!(v, 2) == 9) && (v == CV("Info", [8, 20, 30, 21, 10, 11, 12, 4, 1, 2, 3, 0, -1]))
    @test (splice!(v, 1, 9) == 8) && (v == CV("Info", [9, 20, 30, 21, 10, 11, 12, 4, 1, 2, 3, 0, -1]))
    @test (splice!(v, 1:3) == CV("Info", [9, 20, 30])) && (v == CV("Info", [21, 10, 11, 12, 4, 1, 2, 3, 0, -1]))
    @test (splice!(v, 1:3, [8, 7, 6]) == CV("Info", [21, 10, 11])) && (v == CV("Info", [8, 7, 6, 12, 4, 1, 2, 3, 0, -1]))
    @test deleteat!(v, 4) == CV("Info", [8, 7, 6, 4, 1, 2, 3, 0, -1])
    @test deleteat!(v, [1, 2]) == CV("Info", [6, 4, 1, 2, 3, 0, -1])
    @test deleteat!(v, 5:6) == CV("Info", [6, 4, 1, 2, -1])
    @test (pop!(v) == -1) && (v == CV("Info", [6, 4, 1, 2]))
    @test (popfirst!(v) == 6) && (v == CV("Info", [4, 1, 2]))
    @test (empty!(v) == CV("Info", Int[])) && (v == CV("Info", Int[]))

    v = CV("Info", [1, 3, 2, 2, 4])
    @test empty(v) == CV("Info", Int[])
    @test collect(v) == [1, 3, 2, 2, 4]
    @test keys(v) == keys([1, 3, 2, 2, 4])
    @test values(v) == values([1, 3, 2, 2, 4])
    @test pairs(v) == pairs([1, 3, 2, 2, 4])
    @test reverse(v) == CV("Info", [4, 2, 2, 3, 1])
    @test (sort(v) == CV("Info", [1, 2, 2, 3, 4])) && (v == CV("Info", [1, 3, 2, 2, 4]))
    @test (sort!(v) == CV("Info", [1, 2, 2, 3, 4])) && (v == CV("Info", [1, 2, 2, 3, 4]))
    @test (filter(<=(3), v) == CV("Info", [1, 2, 2, 3])) && (v == CV("Info", [1, 2, 2, 3, 4]))
    @test (filter!(<=(3), v) == CV("Info", [1, 2, 2, 3])) && (v == CV("Info", [1, 2, 2, 3]))
    @test findfirst(isequal(2), v) == 2
    @test findlast(isequal(2), v) == 3
    @test findall(isequal(2), v) == [2, 3]
end

struct CD{S, P, I} <: CompositeDict{P, I}
    info::S
    newcontents::Dict{P, I}
end
contentnames(::Type{<:CD}) = (:info, :contents)
getcontent(m::CD, ::Val{:contents}) = getfield(m, :newcontents)

@testset "CompositeDict" begin
    @test contentnames(CompositeDict) == (:contents,)

    d = CD("Info", Dict("a"=>1, "b"=>2))
    @test contentnames(typeof(d)) == (:info, :contents)
    @test d == deepcopy(d)
    @test isequal(d, deepcopy(d))
    @test isempty(d) == false
    @test length(d) == 2
    @test (haskey(d, "a") == true) && (haskey(d, "d") == false)
    @test (Pair("a", 1) ∈ d) && (Pair("d", 4) ∉ d)
    @test get(d, "a", 2) == 1
    @test get(d, "d", 4) == 4
    @test get(()->4, d, "d") == 4
    @test (get!(d, "d", 4) == 4) && (d == CD("Info", Dict("a"=>1, "b"=>2, "d"=>4)))
    @test (get!(()->5, d, "d") == 4) && (d == CD("Info", Dict("a"=>1, "b"=>2, "d"=>4)))
    @test getkey(d, "e", "e") == "e"
    @test d["d"] == 4
    @test (push!(d, Pair("e", 4)) == CD("Info", Dict("a"=>1, "b"=>2, "d"=>4, "e"=>4))) && (d == CD("Info", Dict("a"=>1, "b"=>2, "d"=>4, "e"=>4)))
    @test (d["d"] = 4; d["e"] = 5; d == CD("Info", Dict("a"=>1, "b"=>2, "d"=>4, "e"=>5)))
    @test (pop!(d) == Pair("e", 5)) && (d == CD("Info", Dict("a"=>1, "b"=>2, "d"=>4)))
    @test (pop!(d, "a") == 1) && (d == CD("Info", Dict("b"=>2, "d"=>4)))
    @test (pop!(d, "a", 1) == 1) && (d == CD("Info", Dict("b"=>2, "d"=>4)))
    @test (delete!(d, "b") == CD("Info", Dict("d"=>4))) && (d == CD("Info", Dict("d"=>4)))
    @test (empty!(d) == CD("Info", Dict{String, Int}())) && (d == CD("Info", Dict{String, Int}()))

    d = CD("Info", Dict("a"=>1, "b"=>2, "c"=>3, "d"=>4))
    @test merge(CD("Info", Dict("a"=>1, "b"=>2)), CD("Info", Dict("c"=>3, "d"=>4))) == d
    @test merge(+, CD("Info", Dict("a"=>1, "b"=>2, "c"=>1)), CD("Info", Dict("c"=>2, "d"=>4))) == d
    @test empty(d) == CD("Info", Dict{String, Int}())
    @test collect(d) == collect(Dict("a"=>1, "b"=>2, "c"=>3, "d"=>4))
    @test collect(keys(d)) == collect(keys(Dict("a"=>1, "b"=>2, "c"=>3, "d"=>4)))
    @test collect(values(d)) == collect(values(Dict("a"=>1, "b"=>2, "c"=>3, "d"=>4)))
    @test collect(pairs(d)) == collect(pairs(Dict("a"=>1, "b"=>2, "c"=>3, "d"=>4)))
    @test (filter(p->p.second<=3, d) == CD("Info", Dict("a"=>1, "b"=>2, "c"=>3))) && (d == CD("Info", Dict("a"=>1, "b"=>2, "c"=>3, "d"=>4)))
    @test (filter!(p->p.second<=3, d) == CD("Info", Dict("a"=>1, "b"=>2, "c"=>3))) && (d == CD("Info", Dict("a"=>1, "b"=>2, "c"=>3)))
end

@testset "NamedContainer" begin
    @test NamedContainer{(:a, :b)}((1, 'h')) == (a=1, b='h')
    @test NamedContainer{()}(()) == NamedTuple()
end
