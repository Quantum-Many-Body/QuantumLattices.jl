using QuantumLattices.Interfaces: ⊕, ⊗, dimension, rank
using QuantumLattices.Prerequisites.VectorSpaces

import QuantumLattices.Prerequisites.Traits: contentnames, getcontent
import QuantumLattices.Prerequisites.VectorSpaces: VectorSpaceStyle, shape

struct SimpleVectorSpace{B, N} <: VectorSpace{B}
    sorted::Bool
    table::NTuple{N, B}
    SimpleVectorSpace(sorted::Bool, table::NTuple{N, B}) where {B, N} = new{B, N}(sorted, table)
end
@inline VectorSpaceStyle(::Type{<:SimpleVectorSpace}) = VectorSpaceEnumerative()
@inline contentnames(::Type{<:SimpleVectorSpace}) = (:sorted, :table)
@inline Base.issorted(vs::SimpleVectorSpace) = vs.sorted

@testset "VectorSpaceEnumerative" begin
    id0, id4 = (1, 0), (1, 4)
    id1, id2, id3 = (1, 1), (1, 2), (1, 3)
    vs = SimpleVectorSpace(true, (id1, id2, id3))
    @test vs == deepcopy(vs)
    @test isequal(vs, deepcopy(vs))
    @test vs|>size == (3,)
    @test vs|>dimension == 3
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

struct SimpleIndices{N} <: VectorSpace{CartesianIndex{N}}
    shape::NTuple{N, UnitRange{Int}}
    SimpleIndices(shape::NTuple{N, UnitRange{Int}}) where N = new{N}(shape)
end
@inline SimpleIndices(shape::UnitRange{Int}...) = SimpleIndices(shape)
@inline VectorSpaceStyle(::Type{<:SimpleIndices}) = VectorSpaceCartesian()
@inline shape(vs::SimpleIndices) = vs.shape
@inline Base.CartesianIndex(basis::CartesianIndex{N}, ::SimpleIndices{N}) where N = basis

@testset "VectorSpaceCartesian" begin
    foi = SimpleIndices(2:3, 2:3, 2:3)
    @test dimension(foi) == 8
    @test issorted(foi) == true
    @test foi|>collect == CartesianIndex.([(2, 2, 2), (3, 2, 2), (2, 3, 2), (3, 3, 2), (2, 2, 3), (3, 2, 3), (2, 3, 3), (3, 3, 3)])
    for (i, index) in enumerate(CartesianIndices((2:3, 2:3, 2:3)))
        @test foi[i] == foi[index] == index
        @test findfirst(index, foi) == i
        @test searchsortedfirst(foi, index) == i
        @test index∈foi
    end
    i1, i2 = CartesianIndex(1, 1, 1), CartesianIndex(4, 4, 4)
    @test isnothing(findfirst(i1, foi)) && isnothing(findfirst(i2, foi))
    @test searchsortedfirst(foi, i1) == 1 && searchsortedfirst(foi, i2) == 9
    @test i1∉foi && i2∉foi
end

struct DirectSummedVectorSpace{T<:Tuple, B} <: VectorSpace{B}
    contents::T
    DirectSummedVectorSpace(contents::Tuple) = new{typeof(contents), mapreduce(eltype, typejoin, contents)}(contents)
end
@inline VectorSpaceStyle(::Type{<:DirectSummedVectorSpace}) = VectorSpaceDirectSummed()
@inline DirectSummedVectorSpace(contents...) = DirectSummedVectorSpace(contents)

@testset "VectorSpaceDirectSummed" begin
    a = SimpleIndices(1:3)
    b = SimpleIndices(3:4, 7:7)
    c = DirectSummedVectorSpace(a, b)
    @test dimension(c) == 5
    @test c[1]==CartesianIndex(1) && c[2]==CartesianIndex(2) && c[3]==CartesianIndex(3) && c[4]==CartesianIndex(3, 7) && c[5]==CartesianIndex(4, 7)
end

@testset "NamedVectorSpace" begin
    t = ParameterSpace{:t}(1:2)
    @test contentnames(typeof(t)) == (:content,)
    @test dimension(t) == 2
    @test t≠ParameterSpace{:μ}(1:2) && t==ParameterSpace{:t}(1:2)
    @test !isequal(t, ParameterSpace{:μ}(1:2)) && isequal(t, ParameterSpace{:t}(1:2))
    @test keys(t) == keys(typeof(t)) == (:t,)
    @test t[1]==1 && t[2]==2
    @test t[1]∈t && t[2]∈t
    ps = pairs(t)
    @test size(ps) == (2,)
    @test eltype(ps) == eltype(typeof(ps)) == NamedTuple{(:t,), Tuple{Int}}
    @test ps[1]==(t=1,) && ps[2]==(t=2,)

    U = ParameterSpace{:U}([8.0, 9.0])
    zps = ZippedNamedVectorSpace(t, U)
    @test VectorSpaceStyle(zps) == VectorSpaceZipped()
    @test eltype(zps) == eltype(typeof(zps)) == Tuple{Int, Float64}
    @test keys(zps) == keys(typeof(zps)) == (:t, :U)
    @test dimension(zps) == 2
    @test zps[1]==(1, 8.0) && zps[2]==(2, 9.0)
    ps = pairs(zps)
    @test size(ps) == (2,)
    @test eltype(ps) == eltype(typeof(ps)) == NamedTuple{(:t, :U), Tuple{Int, Float64}}
    @test ps[1]==(t=1, U=8.0) && ps[2]==(t=2, U=9.0)
    pps = DirectProductedNamedVectorSpace(t, U)
    @test VectorSpaceStyle(pps) == VectorSpaceDirectProducted()
    @test eltype(pps) == eltype(typeof(pps)) == Tuple{Int, Float64}
    @test keys(pps) == keys(typeof(pps)) == (:t, :U)
    @test dimension(pps) == 4
    @test pps[1]==(1, 8.0) && pps[2]==(2, 8.0) && pps[3]==(1, 9.0) && pps[4]==(2, 9.0)

    μ = ParameterSpace{:μ}([11, 12])
    @test t⊕U == zps
    @test μ⊕zps == ZippedNamedVectorSpace(μ, t, U)
    @test zps⊕μ == ZippedNamedVectorSpace(t, U, μ)
    @test zps⊕ZippedNamedVectorSpace(μ) == ZippedNamedVectorSpace(t, U, μ)
    @test t⊗U == pps
    @test μ⊗pps == DirectProductedNamedVectorSpace(μ, t, U)
    @test pps⊗μ == DirectProductedNamedVectorSpace(t, U, μ)
    @test pps⊗DirectProductedNamedVectorSpace(μ) == DirectProductedNamedVectorSpace(t, U, μ)
end
