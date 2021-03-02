using Test
using QuantumLattices.Prerequisites.TypeTraits

import QuantumLattices.Prerequisites.TypeTraits: parameternames, isparameterbound

@testset "reflection" begin
    @test DataType(Int) == Int
    T = Vector{<:Integer}
    @test DataType(T) == T.body

    @test supertype(Vector, :Array) == Vector
    @test supertype(Vector, :AbstractArray) == AbstractVector

    @test parametercount(Vector) == 2
    @test parametertype(Vector{<:Real}, 1) == Real
    @test parametertype(Vector{<:Real}, 2) == 1
    @test parametertypes(Vector{<:Real}) == Tuple{Real, 1}

    @test replace(Vector, 1, Real, true) == Vector{<:Real}
    @test replace(Vector, 1, Real, false) ==  Vector{Real}
    @test replace(Vector{Int}, 2, 3, false) == Array{Int, 3}

    @test rawtype(Int) == Int
    @test rawtype(Vector{Int}) == Array
    @test rawtype(Vector) == Array
    @test rawtype(Vector{<:Real}) == Array

    @test fulltype(Array, Tuple{Int, 1}, (false, false)) == Array{Int, 1}
    @test fulltype(Array, Tuple{Real, 2}, (true, false)) == Array{<:Real, 2}
end

abstract type FT{T} end
Base.fieldnames(::Type{Field}, ::Type{<:FT}) = (:content,)
parameternames(::Type{<:FT}) = (:content,)
isparameterbound(::Type{<:FT}, ::Parameter{:content}, D) = !isconcretetype(D)
Base.map(f::Function, m::FT, a::Field{:content}, args::Tuple, kwargs::NamedTuple) = f(getproperty(m, a), args...; kwargs...)

struct SameNameField <: FT{Vector{Int}}
    info::String
    content::Vector{Int}
end
parameternames(::Type{SameNameField}) = ()

struct DiffNameField <: FT{Vector{Int}}
    info::String
    table::Vector{Int}
end
Base.fieldname(::Type{DiffNameField}, ::Field{:content}) = :table
parameternames(::Type{DiffNameField}) = ()

@testset "AbstractTypeHelper" begin
    @test fieldnames(Field, FT) == (:content,)
    @test fieldname(FT, Field(:content)) == :content

    @test promote_type(Parameter, NamedTuple{(:a, :b), Tuple{Int, Float64}}, NamedTuple{(:b, :c), Tuple{Complex{Float64}, Int}}) == NamedTuple{(:a, :b, :c), Tuple{Int, Complex{Float64}, Int}}

    @test parametername(FT, 1) == :content
    @test parametertype(FT{Vector}, 1) == parametertype(FT{Vector}, Parameter(:content)) == Vector
    @test parametertype(FT{Vector{Int}}, 1) == parametertype(FT{Vector{Int}}, Parameter(:content)) == Vector{Int}
    @test parameterpair(FT{Vector}, Parameter(:content)) == Pair{:content, Vector}
    @test isparameterbound(FT, Parameter(:content), Vector) == true
    @test isparameterbound(FT, Parameter(:content), Vector{Int}) == false
    @test hasparameter(FT, Parameter(:content)) == true
    @test parameternames(FT) == (:content,)
    @test parametertypes(FT{Vector}) == Tuple{Vector}
    @test parameterpairs(FT{Vector}) == NamedTuple{(:content,), Tuple{Vector}}
    @test isparameterbounds(FT, NamedTuple{(:content,), Tuple{Vector}}) == (true,)
    @test isparameterbounds(FT, NamedTuple{(:content,), Tuple{Vector{Int}}}) == (false,)
    @test findfirst(Parameter(:content), FT) == 1
    @test replace(FT{Vector{Int}}, Parameter(:content), Vector) == FT{<:Vector}
    @test replace(FT{Vector{Int}}, Parameter(:content), Vector{Float64}) == FT{Vector{Float64}}
    @test fulltype(FT, NamedTuple{(:content,), Tuple{Vector}}) == FT{<:Vector}
    @test fulltype(FT, NamedTuple{(:content,), Tuple{Vector{Int}}}) == FT{Vector{Int}}

    @test fieldnames(Field, SameNameField) == (:content,)
    @test fieldname(SameNameField, Field(:content)) == :content
    @test hasparameter(SameNameField, Parameter(:content)) == false

    a = SameNameField("abcd", [1, 2, 3, 4])
    @test getproperty(a, Field(:content)) == [1, 2, 3, 4]
    @test Tuple(Field, a, getindex, ((1,))) == ("abcd", 1)

    @test fieldnames(Field, DiffNameField) == (:content,)
    @test fieldname(DiffNameField, Field(:content)) == :table
    @test hasparameter(DiffNameField, Parameter(:content)) == false

    b = DiffNameField("abcd", [1, 2, 3, 4])
    @test getproperty(b, Field(:content)) == [1, 2, 3, 4]
    @test Tuple(Field, b, getindex, ((2,))) == ("abcd", 2)
end

struct EFO{F1, F2, F3}
    f1::F1
    f2::F2
    f3::F3
end
Base.:(==)(fc1::EFO, fc2::EFO) = ==(efficientoperations, fc1, fc2)
Base.isequal(fc1::EFO, fc2::EFO) = isequal(efficientoperations, fc1, fc2)
Base.:<(fc1::EFO, fc2::EFO) = <(efficientoperations, fc1, fc2)
Base.isless(fc1::EFO, fc2::EFO) = isless(efficientoperations, fc1, fc2)
Base.isapprox(fc1::EFO, fc2::EFO; atol::Real=10^-5, rtol::Real=10^-5) = isapprox(efficientoperations, Val((:f1, :f2)), fc1, fc2; atol=atol, rtol=rtol)
Base.replace(fc::EFO; kwargs...) = replace(efficientoperations, fc; kwargs...)

@testset "efficientoperations" begin
    @test ==(efficientoperations, (), ())
    @test isequal(efficientoperations, (), ())
    @test ==(efficientoperations, (1, 2), (1, 2, 3)) == false
    @test isequal(efficientoperations, (1, 2), (1, 2, 3)) == false

    fc1, fc2 = EFO(1.0, 2, 3), EFO(1, 2, 3)
    @test fc1 == fc2
    @test isequal(fc1, fc2)
    fc1, fc2 = EFO(1.0, 2, 3), EFO(1, 2, 4.0)
    @test fc1 < fc2
    @test isless(fc1, fc2)

    fc1, fc2 = EFO(1.0+10^-6, 2, 3), EFO(1, 2, 3)
    @test fc1 ≈ fc2
    fc1, fc2 = EFO(1.0, 2-10^-6, 3), EFO(1, 2, 3)
    @test fc1 ≈ fc2
    fc1, fc2 = EFO(1.0, 2, 3-10^-6), EFO(1, 2, 3)
    @test fc1 ≉  fc2

    @test replace(EFO(1, 2, 3), f1='c') == EFO('c', 2, 3)
end

@testset "forder/corder" begin
    dims = (2, 2, 2)
    finds = [(1, 1, 1), (2, 1, 1), (1, 2, 1), (2, 2, 1), (1, 1, 2), (2, 1, 2), (1, 2, 2), (2, 2, 2)]
    cinds = [(1, 1, 1), (1, 1, 2), (1, 2, 1), (1, 2, 2), (2, 1, 1), (2, 1, 2), (2, 2, 1), (2, 2, 2)]
    @test collect(indtosub(dims, i, forder) for i = 1:prod(dims)) == finds
    @test collect(indtosub(dims, i, corder) for i = 1:prod(dims)) == cinds
    @test collect(subtoind(dims, inds, forder) for inds in finds) == collect(1:prod(dims))
    @test collect(subtoind(dims, inds, corder) for inds in cinds) == collect(1:prod(dims))
end
