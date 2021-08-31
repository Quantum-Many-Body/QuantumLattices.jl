using Test
using QuantumLattices.Mathematics.VectorSpaces
using QuantumLattices.Interfaces: dimension, rank

import QuantumLattices.Prerequisites.Traits: contentnames, getcontent

struct SimpleVectorSpace{B, N} <: EnumerativeVectorSpace{B}
    sorted::Bool
    table::NTuple{N, B}
    SimpleVectorSpace(sorted::Bool, table::NTuple{N, B}) where {B, N} = new{B, N}(sorted, table)
end
@inline contentnames(::Type{<:SimpleVectorSpace}) = (:sorted, :table)
@inline Base.issorted(vs::SimpleVectorSpace) = vs.sorted

@testset "EnumerativeVectorSpace" begin
    @test contentnames(EnumerativeVectorSpace) == (:table,)

    id0, id4 = (1, 0), (1, 4)
    id1, id2, id3 = (1, 1), (1, 2), (1, 3)
    vs = SimpleVectorSpace(true, (id1, id2, id3))
    @test vs == deepcopy(vs)
    @test isequal(vs, deepcopy(vs))
    @test vs|>size == (3,)
    @test vs|>dimension == 3
    @test vs|>ndims == vs|>typeof|>ndims == 1
    @test vs|>collect == [id1, id2, id3]
    @test vs[1]==id1 && vs[2]==id2 && vs[3]==id3
    @test searchsortedfirst(vs, id0)==1 && searchsortedfirst(vs, id4)==4
    @test searchsortedfirst(vs, id1)==1 && searchsortedfirst(vs, id2)==2 && searchsortedfirst(vs, id3)==3
    @test isnothing(findfirst(id0, vs)) && isnothing(findfirst(id4, vs))
    @test findfirst(id1, vs)==1 && findfirst(id2, vs)==2 && findfirst(id3, vs)==3
    @test (id0∉vs) && (id4 ∉ vs)
    @test (id1∈vs) && (id2∈vs) && (id3 ∈ vs)

    vs = SimpleVectorSpace(false, (id1, id2, id3))
    @test searchsortedfirst(vs, id0)==4 && searchsortedfirst(vs, id4)==4
    @test searchsortedfirst(vs, id1)==1 && searchsortedfirst(vs, id2)==2 && searchsortedfirst(vs, id3)==3
    @test isnothing(findfirst(id0, vs)) && isnothing(findfirst(id4, vs))
    @test findfirst(id1, vs)==1 && findfirst(id2, vs)==2 && findfirst(id3, vs)==3
end

struct SimpleIndices{N} <: CartesianVectorSpace{CartesianIndex{N}}
    dims::NTuple{N, Int}
    SimpleIndices(dims::NTuple{N, Int}) where N = new{N}(dims)
end
SimpleIndices(dims::Int...) = SimpleIndices(dims)
@inline Base.Dims(vs::SimpleIndices) = vs.dims
@inline Base.CartesianIndex(basis::CartesianIndex{N}, ::SimpleIndices{N}) where N = basis

@testset "CartesianVectorSpace" begin
    foi = SimpleIndices(2, 2, 2)
    @test dimension(foi) == 8
    @test issorted(foi) == true
    @test foi|>collect == CartesianIndex.([(1, 1, 1), (2, 1, 1), (1, 2, 1), (2, 2, 1), (1, 1, 2), (2, 1, 2), (1, 2, 2), (2, 2, 2)])
    for (i, index) in enumerate(CartesianIndices((2,2,2)))
        @test foi[i] == foi[index] == index
        @test findfirst(index, foi) == i
        @test searchsortedfirst(foi, index) == i
        @test index∈foi
    end
    i1, i2 = CartesianIndex(0, 0, 0), CartesianIndex(3, 3, 3)
    @test isnothing(findfirst(i1, foi)) && isnothing(findfirst(i2, foi))
    @test searchsortedfirst(foi, i1) == 1 && searchsortedfirst(foi, i2) == 9
    @test i1∉foi && i2∉foi
end

struct VSZNamedVectorSpace{NS, BS<:Tuple, VS<:Tuple{Vararg{Vector}}} <: NamedVectorSpace{:⊕, NS, BS, VS}
    tables::VS
end
@inline getcontent(m::VSZNamedVectorSpace, ::Val{:contents}) = getfield(m, :tables)
@generated function VSZNamedVectorSpace{NS}(contents::Vector...) where NS
    @assert (length(NS) == length(contents)) && isa(NS, Tuple{Vararg{Symbol}})
    BS = Expr(:curly, :Tuple, [contents[i]|>eltype for i = 1:length(NS)]...)
    return quote
        @assert mapreduce(length, ==, contents)
        VSZNamedVectorSpace{NS, $BS, typeof(contents)}(contents)
    end
end

struct VSPNamedVectorSpace{NS, BS<:Tuple, VS<:Tuple{Vararg{Vector}}} <: NamedVectorSpace{:⊗, NS, BS, VS}
    tables::VS
end
@inline getcontent(m::VSPNamedVectorSpace, ::Val{:contents}) = getfield(m, :tables)
@generated function VSPNamedVectorSpace{NS}(contents::Vector...) where NS
    @assert (length(NS) == length(contents)) && isa(NS, Tuple{Vararg{Symbol}})
    BS = Expr(:curly, :Tuple, [contents[i]|>eltype for i = 1:length(NS)]...)
    return :(VSPNamedVectorSpace{NS, $BS, typeof(contents)}(contents))
end

@testset "NamedVectorSpace" begin
    nvs = VSZNamedVectorSpace{(:t, :U)}([1, 2], [8.0, 9.0])
    @test nvs|>keys == nvs|>typeof|>keys == (:t, :U)
    @test nvs|>values == ([1, 2], [8.0, 9.0])
    @test nvs|>pairs|>collect == [:t=>[1, 2], :U=>[8.0, 9.0]]
    @test eltype(nvs, 1) == eltype(nvs|>typeof, 1) == Int
    @test eltype(nvs, 2) == eltype(nvs|>typeof, 2) == Float64
    @test nvs|>typeof|>rank == 1
    @test Dims(nvs) == (2,)
    elements = [(t = 1, U = 8.0), (t = 2, U = 9.0)]
    for i = 1:dimension(nvs)
        @test NamedTuple(CartesianIndex(elements[i], nvs), nvs) == elements[i]
    end
    @test nvs|>collect == elements

    nvs = VSPNamedVectorSpace{(:t, :U)}([1.0, 2.0], [8.0, 9.0])
    @test nvs|>typeof|>rank == 2
    @test Dims(nvs) == (2, 2)
    elements = [(t = 1.0, U = 8.0), (t = 2.0, U = 8.0), (t = 1.0, U = 9.0), (t = 2.0, U = 9.0)]
    for i = 1:dimension(nvs)
        @test NamedTuple(CartesianIndex(elements[i], nvs), nvs) == elements[i]
    end
    @test nvs|>collect == elements
end
