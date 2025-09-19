using DataStructures: OrderedDict
using OffsetArrays: OffsetArray
using QuantumLattices: id, shape, str, value
using QuantumLattices.Toolkit
using StaticArrays: SVector

import QuantumLattices.Toolkit: VectorSpaceStyle, contentnames, contenttype, dissolve, getcontent, isparameterbound, parameternames

@testset "str" begin
    @test str((1, 2)) == "(1, 2)"
    @test str(1) == "1"
    @test str(10^6) == "1000000"
    @test str(1//7) == "1//7"
    @test str(7//1) == "7"
    @test str(1.0) == "1.0"
    @test str(10.0^6) == "1.0e+06"
    @test str(10^-5) == "1.0e-05"
    @test str(1/7) == "0.14286"
    @test str(1/7; ndecimal=8) == "0.14285714"
    @test str(0im) == "0.0"
    @test str(1//7+0im) == "1//7"
    @test str(1.0im) == "1.0im"
    @test str(-1.0im) == "-1.0im"
    @test str(0.1+0.12im) == "0.1+0.12im"
    @test str(0.1-0.12im) == "0.1-0.12im"
    @test str(:a) =="a"
    @test str(:) == ":"
    @test str('a') == "'a'"
end

@testset "superscript" begin
    @test superscript(1234567890) == "¹²³⁴⁵⁶⁷⁸⁹⁰"
end

@testset "subscript" begin
    @test subscript(1234567890) == "₁₂₃₄₅₆₇₈₉₀"
end

@testset "delta" begin
    @test delta(1, 2) == 0
    @test delta(1, 1) == 1
end

@testset "concatenate" begin
    @test concatenate((1, 2), (:a, :b)) == (1, 2, :a, :b)
end

@testset "OrderedDict" begin
    od = OrderedDict(2=>'a', 1=>'b', 3=>'c')
    @test id(od, 1)==2 && id(od, 2)==1 && id(od, 3)==3
    @test value(od, 1)=='a' && value(od, 2)=='b' && value(od, 3)=='c'
end

@testset "searchsortedfirst" begin
    idx = CartesianIndices((2, 2, 2))
    for (i, index) in enumerate(idx)
        @test searchsortedfirst(idx, index) == i
    end
    @test searchsortedfirst(idx, CartesianIndex(0, 0, 0)) == 1
    @test searchsortedfirst(idx, CartesianIndex(3, 3, 3)) == 9
end

@testset "DirectSummedIndices" begin
    indexes₁ = DirectSummedIndices((-3:-1, 1:3, 1:4))
    indexes₂ = collect(indexes₁)
    for i in eachindex(indexes₁, indexes₂)
        @test indexes₁[i] == indexes₂[i]
    end
end

@testset "DirectProductedIndices" begin
    forward = DirectProductedIndices{:forward}((-2:-1, 2:3))
    @test length(forward) == 4
    @test collect(forward) == [CartesianIndex(-2, 2), CartesianIndex(-1, 2), CartesianIndex(-2, 3), CartesianIndex(-1, 3)]
    for i in eachindex(forward)
        @test forward[i] in forward
        @test searchsortedfirst(forward, forward[i]) == i
    end
    @test CartesianIndex(1, 1) ∉ forward

    backward = DirectProductedIndices{:backward}((-2:-1, 2:3))
    @test length(backward) == 4
    @test collect(backward) == [CartesianIndex(-2, 2), CartesianIndex(-2, 3), CartesianIndex(-1, 2), CartesianIndex(-1, 3)]
    for i in eachindex(backward)
        @test backward[i] in backward
        @test searchsortedfirst(backward, backward[i]) == i
    end
    @test CartesianIndex(1, 1) ∉ backward

    @test DirectProductedIndices{:backward}((2, 2:3)) == DirectProductedIndices{:backward}((1:2, 2:3))
end

@testset "Segment" begin
    segment = Segment(SVector(0.0, 0.0), SVector(1.0, 1.0), 10)
    @test segment == Segment((0.0, 0.0), (1.0, 1.0), 10)
    @test isequal(segment,  Segment((0.0, 0.0), (1.0, 1.0), 10))
    @test size(segment) == (10,)
    @test string(segment) == "[p₁, p₂) with p₁=[0.0, 0.0] and p₂=[1.0, 1.0]"
    @test segment[1:4] == Segment(segment[1], segment[4], 4, (true, true))
    @test segment[6:-1:2] == Segment(segment[6], segment[2], 5, (true, true))

    segment = Segment(1, 5, 5, ends=(true, true))
    @test string(segment) == "[1.0, 5.0]"
    @test segment[1]==1 && segment[end]==5
    for (i, seg) in enumerate(segment)
        @test seg ≈ segment[i]
    end

    segment = Segment(1, 5, 4, ends=(true, false))
    @test string(segment) == "[1.0, 5.0)"
    @test segment[1]==1 && segment[end]==4
    for (i, seg) in enumerate(segment)
        @test seg ≈ segment[i]
    end

    segment = Segment(1, 5, 4, ends=(false, true))
    @test string(segment) == "(1.0, 5.0]"
    @test segment[1]==2 && segment[end]==5
    for (i, seg) in enumerate(segment)
        @test seg ≈ segment[i]
    end

    segment = Segment(1, 5, 3, ends=(false, false))
    @test string(segment) == "(1.0, 5.0)"
    @test segment[1]==2 && segment[end]==4
    for (i, seg) in enumerate(segment)
        @test seg ≈ segment[i]
    end
end

@testset "Combinations" begin
    @test Combinations{0}("abc")|>collect == [()]
    @test Combinations{1}("abc")|>collect == [('a',), ('b',), ('c',)]
    @test Combinations{2}("abc")|>collect == [('a', 'b'), ('a', 'c'), ('b', 'c')]
    @test Combinations{3}("abc")|>collect == [('a', 'b', 'c')]

    @test Combinations{0}(OffsetArray(['a', 'b' , 'c'], -3:-1))|>collect == [()]
    @test Combinations{1}(OffsetArray(['a', 'b' , 'c'], -3:-1))|>collect == [('a',), ('b',), ('c',)]
    @test Combinations{2}(OffsetArray(['a', 'b' , 'c'], -3:-1))|>collect == [('a', 'b'), ('a', 'c'), ('b', 'c')]
    @test Combinations{3}(OffsetArray(['a', 'b' , 'c'], -3:-1))|>collect == [('a', 'b', 'c')]
end

@testset "DuplicateCombinations" begin
    @test DuplicateCombinations{0}("abc")|>collect == [()]
    @test DuplicateCombinations{1}("abc")|>collect == [('a',), ('b',), ('c',)]
    @test DuplicateCombinations{2}("abc")|>collect == [('a', 'a'), ('a', 'b'), ('a', 'c'), ('b', 'b'), ('b', 'c'), ('c', 'c')]
    @test DuplicateCombinations{3}("abc")|>collect == [('a', 'a', 'a'), ('a', 'a', 'b'), ('a', 'a', 'c'), ('a', 'b', 'b'), ('a', 'b', 'c'), ('a', 'c', 'c'), ('b', 'b', 'b'), ('b', 'b', 'c'), ('b', 'c', 'c'), ('c', 'c', 'c')]

    @test DuplicateCombinations{0}(OffsetArray(['a', 'b' , 'c'], -3:-1))|>collect == [()]
    @test DuplicateCombinations{1}(OffsetArray(['a', 'b' , 'c'], -3:-1))|>collect == [('a',), ('b',), ('c',)]
    @test DuplicateCombinations{2}(OffsetArray(['a', 'b' , 'c'], -3:-1))|>collect == [('a', 'a'), ('a', 'b'), ('a', 'c'), ('b', 'b'), ('b', 'c'), ('c', 'c')]
    @test DuplicateCombinations{3}(OffsetArray(['a', 'b' , 'c'], -3:-1))|>collect == [('a', 'a', 'a'), ('a', 'a', 'b'), ('a', 'a', 'c'), ('a', 'b', 'b'), ('a', 'b', 'c'), ('a', 'c', 'c'), ('b', 'b', 'b'), ('b', 'b', 'c'), ('b', 'c', 'c'), ('c', 'c', 'c')]
end

@testset "Permutations" begin
    @test Permutations{0}("abc")|>collect == [()]
    @test Permutations{1}("abc")|>collect == [('a',), ('b',), ('c',)]
    @test Permutations{2}("abc")|>collect == [('a', 'b'), ('a', 'c'), ('b', 'a'), ('b', 'c'), ('c', 'a'), ('c', 'b')]
    @test Permutations{3}("abc")|>collect == [('a', 'b', 'c'), ('a', 'c', 'b'), ('b', 'a', 'c'), ('b', 'c', 'a'), ('c', 'a', 'b'), ('c', 'b', 'a')]

    @test Permutations{0}(OffsetArray(['a', 'b' , 'c'], -3:-1))|>collect == [()]
    @test Permutations{1}(OffsetArray(['a', 'b' , 'c'], -3:-1))|>collect == [('a',), ('b',), ('c',)]
    @test Permutations{2}(OffsetArray(['a', 'b' , 'c'], -3:-1))|>collect == [('a', 'b'), ('a', 'c'), ('b', 'a'), ('b', 'c'), ('c', 'a'), ('c', 'b')]
    @test Permutations{3}(OffsetArray(['a', 'b' , 'c'], -3:-1))|>collect == [('a', 'b', 'c'), ('a', 'c', 'b'), ('b', 'a', 'c'), ('b', 'c', 'a'), ('c', 'a', 'b'), ('c', 'b', 'a')]
end

@testset "DuplicatePermutations" begin
    @test DuplicatePermutations{0}("abc")|>collect == [()]
    @test DuplicatePermutations{1}("abc")|>collect == [('a',), ('b',), ('c',)]
    @test DuplicatePermutations{2}("abc")|>collect == [('a', 'a'), ('a', 'b'), ('a', 'c'), ('b', 'a'), ('b', 'b'), ('b', 'c'), ('c', 'a'), ('c', 'b'), ('c', 'c')]
    @test DuplicatePermutations{3}("abc")|>collect == [(i, j, k) for i∈"abc" for j∈"abc" for k∈"abc"]

    @test DuplicatePermutations{0}(OffsetArray(['a', 'b' , 'c'], -3:-1))|>collect == [()]
    @test DuplicatePermutations{1}(OffsetArray(['a', 'b' , 'c'], -3:-1))|>collect == [('a',), ('b',), ('c',)]
    @test DuplicatePermutations{2}(OffsetArray(['a', 'b' , 'c'], -3:-1))|>collect == [('a', 'a'), ('a', 'b'), ('a', 'c'), ('b', 'a'), ('b', 'b'), ('b', 'c'), ('c', 'a'), ('c', 'b'), ('c', 'c')]
    @test DuplicatePermutations{3}(OffsetArray(['a', 'b' , 'c'], -3:-1))|>collect == [(i, j, k) for i∈"abc" for j∈"abc" for k∈"abc"]
end

abstract type FT{T} end
@inline parameternames(::Type{<:FT}) = (:content,)
@inline isparameterbound(::Type{<:FT}, ::Val{:content}, ::Type{D}) where D = !isconcretetype(D)
@inline contentnames(::Type{<:FT}) = (:content,)
@inline dissolve(m::FT, ::Val{:content}, f::Function, args...; kwargs...) = f(getcontent(m, :content), args...; kwargs...)

struct SameNameField <: FT{Vector{Int}}
    info::String
    content::Vector{Int}
end
@inline parameternames(::Type{SameNameField}) = ()
@inline contentnames(::Type{SameNameField}) = (:info, :content)

struct DiffNameField <: FT{Vector{Int}}
    info::String
    table::Vector{Int}
end
@inline parameternames(::Type{DiffNameField}) = ()
@inline contentnames(::Type{DiffNameField}) = (:info, :content)
@inline getcontent(m::DiffNameField, ::Val{:content}) = getfield(m, :table)
@inline contenttype(DF::Type{DiffNameField}, ::Val{:content}) = fieldtype(DF, :table)

@testset "SubTypeTree" begin
    subtypetree = SubTypeTree(Integer)
    @test eltype(subtypetree) == eltype(typeof(subtypetree)) == Type
    @test collect(subtypetree) == [Integer, Bool, Signed, BigInt, Int128, Int16, Int32, Int64, Int8, Unsigned, UInt128, UInt16, UInt32, UInt64, UInt8]
    @test string(subtypetree) == "Integer\n├─ Bool\n├─ Signed\n│  ├─ BigInt\n│  ├─ Int128\n│  ├─ Int16\n│  ├─ Int32\n│  ├─ Int64\n│  └─ Int8\n└─ Unsigned\n   ├─ UInt128\n   ├─ UInt16\n   ├─ UInt32\n   ├─ UInt64\n   └─ UInt8\n"
end

@testset "commontype" begin
    @test commontype(Int) == Int
    @test commontype(Union{Int, Float64}) == Float64
    @test commontype(Union{Int, Float64, Complex{Int}}) == Complex{Float64}

    fx = x::Int->(x==1 ? 1 : x==2 ? 2.0 : 3.0im)
    @test commontype(fx, Tuple{Int}) == Complex{Float64}
end

@testset "abstracttypehelper" begin
    @test DataType(Int) == Int
    T = Vector{<:Integer}
    @test DataType(T) == T.body

    @test supertype(Vector, :Array) == Vector
    @test supertype(Vector, :AbstractArray) == AbstractVector

    @test parametercount(Vector) == 2
    @test parametertype(Vector{<:Real}, 1) == Real
    @test parametertype(Vector{<:Real}, 2) == 1
    @test parametertypes(Vector{<:Real}) == Tuple{Real, 1}
    @test reparameter(Vector, 1, Real, true) == Vector{<:Real}
    @test reparameter(Vector, 1, Real, false) ==  Vector{Real}
    @test reparameter(Vector{Int}, 2, 3, false) == Array{Int, 3}

    @test rawtype(Int) == Int
    @test rawtype(Vector{Int}) == Array
    @test rawtype(Vector) == Array
    @test rawtype(Vector{<:Real}) == Array
    @test fulltype(Array, Tuple{Int, 1}, (false, false)) == Array{Int, 1}
    @test fulltype(Array, Tuple{Real, 2}, (true, false)) == Array{<:Real, 2}

    @test promoteparameters(NamedTuple{(:a, :b), Tuple{Int, Float64}}, NamedTuple{(:b, :c), Tuple{Complex{Float64}, Int}}) == NamedTuple{(:a, :b, :c), Tuple{Int, Complex{Float64}, Int}}

    @test parametercount(FT) == 1
    @test parametername(FT, 1) == :content
    @test parameterorder(FT, :content) == 1
    @test parametertype(FT{Vector}, 1) == parametertype(FT{<:Vector}, :content) == Vector
    @test parametertype(FT{Vector{Int}}, 1) == parametertype(FT{Vector{Int}}, :content) == Vector{Int}
    @test parameterpair(FT{Vector}, :content) == parameterpair(FT{Vector}, 1) == Pair{:content, Vector}
    @test isparameterbound(FT, :content, Vector) == true
    @test isparameterbound(FT, :content, Vector{Int}) == false
    @test isparameterbound(FT, :content, 1) == false
    @test hasparameter(FT, :content) == true

    @test parameternames(FT) == (:content,)
    @test parametertypes(FT{Vector}) == Tuple{Vector}
    @test parameterpairs(FT{Vector}) == NamedTuple{(:content,), Tuple{Vector}}
    @test isparameterbounds(FT, NamedTuple{(:content,), Tuple{Vector}}) == (true,)
    @test isparameterbounds(FT, NamedTuple{(:content,), Tuple{Vector{Int}}}) == (false,)
    @test reparameter(FT{Vector{Int}}, :content, Vector) == FT{<:Vector}
    @test reparameter(FT{Vector{Int}}, :content, Vector{Float64}) == FT{Vector{Float64}}

    @test fulltype(FT, NamedTuple{(:content,), Tuple{Vector}}) == FT{<:Vector}
    @test fulltype(FT, NamedTuple{(:content,), Tuple{Vector{Int}}}) == FT{Vector{Int}}

    @test contentnames(FT) == (:content,)
    @test contentname(FT, 1) == :content
    @test contentorder(FT, :content) == 1
    @test contentcount(FT) == 1
    @test hascontent(FT, :content) == true

    @test contentnames(SameNameField) == (:info, :content)
    @test hasparameter(SameNameField, :content) == false
    @test contenttype(SameNameField, :info) == contenttype(SameNameField, Val(:info)) == String
    @test contenttype(SameNameField, :content) == contenttype(SameNameField, Val(:content)) == Vector{Int}
    @test contenttypes(SameNameField) == Tuple{String, Vector{Int}}

    a = SameNameField("abcd", [1, 2, 3, 4])
    @test getcontent(a, :content) == getcontent(a, 2) == getcontent(a, Val(:content)) == [1, 2, 3, 4]
    @test dissolve(a, getindex, 1) == ("abcd", 1)

    @test contentnames(DiffNameField) == (:info, :content)
    @test hasparameter(DiffNameField, :content) == false
    @test contenttype(DiffNameField, :info) == contenttype(DiffNameField, Val(:info)) == String
    @test contenttype(DiffNameField, :content) == contenttype(DiffNameField, Val(:content)) == Vector{Int}
    @test contenttypes(DiffNameField) == Tuple{String, Vector{Int}}

    b = DiffNameField("abcd", [1, 2, 3, 4])
    @test getcontent(b, :content) == getcontent(b, 2) == getcontent(b, Val(:content)) == [1, 2, 3, 4]
    @test dissolve(b, getindex, 2) == ("abcd", 2)
end

struct EFO{F₁, F₂, F₃}
    f₁::F₁
    f₂::F₂
    f₃::F₃
end
@inline Base.:(==)(fc₁::EFO, fc₂::EFO) = ==(efficientoperations, fc₁, fc₂)
@inline Base.isequal(fc₁::EFO, fc₂::EFO) = isequal(efficientoperations, fc₁, fc₂)
@inline Base.:<(fc₁::EFO, fc₂::EFO) = <(efficientoperations, fc₁, fc₂)
@inline Base.isless(fc₁::EFO, fc₂::EFO) = isless(efficientoperations, fc₁, fc₂)
@inline Base.isapprox(fc₁::EFO, fc₂::EFO; atol::Real=10^-5, rtol::Real=10^-5) = isapprox(efficientoperations, (:f₁, :f₂), fc₁, fc₂; atol=atol, rtol=rtol)

@testset "efficientoperations" begin
    @test ==(efficientoperations, (), ())
    @test isequal(efficientoperations, (), ())
    @test ==(efficientoperations, (1, 2), (1, 2, 3)) == false
    @test isequal(efficientoperations, (1, 2), (1, 2, 3)) == false
    @test isapprox(efficientoperations, (1.0, 2.0), (1, 2)) == true
    @test isapprox(efficientoperations, (), (), ()) == true
    @test isapprox(efficientoperations, (), (), (1,)) == false

    fc₁, fc₂ = EFO(1.0, 2, 3), EFO(1, 2, 3)
    @test fc₁ == fc₂
    @test isequal(fc₁, fc₂)
    fc₁, fc₂ = EFO(1.0, 2, 3), EFO(1, 2, 4.0)
    @test fc₁ < fc₂
    @test isless(fc₁, fc₂)

    fc₁, fc₂ = EFO(1.0+10^-6, 2, 3), EFO(1, 2, 3)
    @test fc₁ ≈ fc₂
    fc₁, fc₂ = EFO(1.0, 2-10^-6, 3), EFO(1, 2, 3)
    @test fc₁ ≈ fc₂
    fc₁, fc₂ = EFO(1.0, 2, 3-10^-6), EFO(1, 2, 3)
    @test fc₁ ≉  fc₂
end

struct CV{S, T} <: CompositeVector{T}
    info::S
    contents::Vector{T}
end
@inline contentnames(::Type{<:CV}) = (:info, :contents)

@testset "CompositeVector" begin
    @test contentnames(CompositeVector) == (:contents,)

    v = CV("Info", [1, 3, 2, 4])
    @test contentnames(typeof(v)) == (:info, :contents)
    @test axes(v) == (Base.OneTo(4),)
    @test size(v) == (4,)
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
@inline contentnames(::Type{<:CD}) = (:info, :contents)
@inline getcontent(m::CD, ::Val{:contents}) = getfield(m, :newcontents)

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

struct SimpleVectorSpace{B, N} <: VectorSpace{B}
    contents::NTuple{N, B}
end
@inline VectorSpaceStyle(::Type{<:SimpleVectorSpace}) = VectorSpaceEnumerative()
@inline SimpleVectorSpace(contents...) = SimpleVectorSpace(contents)

@testset "VectorSpaceEnumerative" begin
    id₀, id₄ = (1, 0), (1, 4)
    id₁, id₂, id₃ = (1, 1), (1, 2), (1, 3)
    vs = SimpleVectorSpace(id₁, id₂, id₃)
    @test vs==deepcopy(vs) && isequal(vs, deepcopy(vs))
    @test axes(vs) == (Base.OneTo(3),)
    @test vs|>size == (3,)
    @test vs|>length == 3
    @test IndexStyle(vs) == IndexLinear()
    @test vs|>collect == [id₁, id₂, id₃]
    @test vs[1]==id₁ && vs[2]==id₂ && vs[3]==id₃
    @test vs[1:3] == [id₁, id₂, id₃]
    @test searchsortedfirst(vs, id₀)==1 && searchsortedfirst(vs, id₄)==4
    @test searchsortedfirst(vs, id₁)==1 && searchsortedfirst(vs, id₂)==2 && searchsortedfirst(vs, id₃)==3
    @test (id₀∉vs) && (id₄∉vs)
    @test (id₁∈vs) && (id₂∈vs) && (id₃∈vs)
end

@testset "VectorSpaceDirectSummed" begin
    a = SimpleVectorSpace(1, 2, 3)
    b = SimpleVectorSpace((3, 7), (4, 7))
    c = DirectSummedVectorSpace(a, b)
    @test length(c) == 5
    @test c[1]==1 && c[2]==2 && c[3]==3 && c[4]==(3, 7) && c[5]==(4, 7)
    for (i, index) in enumerate(DirectSummedIndices((eachindex(a), eachindex(b))))
        @test Int(index, c) == Int(CartesianIndex(index[1], index[2]), c) == i
    end
end

@testset "VectorSpaceDirectProducted" begin
    a = SimpleVectorSpace(1, 2)
    b = SimpleVectorSpace(4, 5)

    forward = DirectProductedVectorSpace{:forward}(a, b)
    @test shape(forward) == (Base.OneTo(2), Base.OneTo(2))
    @test length(forward) == 4
    @test collect(forward) == [(1, 4), (2, 4), (1, 5), (2, 5)]
    for i in eachindex(forward)
        @test forward[i] ∈ forward
        @test searchsortedfirst(forward, forward[i]) == i
    end
    for (i, index) in enumerate(DirectProductedIndices{:forward}(shape(forward)))
        @test Int(index, forward) == i
    end

    backward = DirectProductedVectorSpace{:backward}(a, b)
    @test shape(backward) == (Base.OneTo(2), Base.OneTo(2))
    @test length(backward) == 4
    @test collect(backward) == [(1, 4), (1, 5), (2, 4), (2, 5)]
    for i in eachindex(backward)
        @test backward[i] ∈ backward
        @test searchsortedfirst(backward, backward[i]) == i
    end
    for (i, index) in enumerate(DirectProductedIndices{:backward}(shape(backward)))
        @test Int(index, backward) == i
    end
end

@testset "VectorSpaceZipped" begin
    a = SimpleVectorSpace(1, 2, 3)
    b = SimpleVectorSpace(4, 5, 6)
    c = ZippedVectorSpace(a, b)
    @test length(c) == 3
    @test collect(c) == [(1, 4), (2, 5), (3, 6)]
    for i = 1:length(c)
        @test Int(CartesianIndex(i, i), c) == i
    end
end
