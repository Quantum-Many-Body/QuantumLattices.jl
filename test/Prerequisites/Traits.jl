using Test
using QuantumLattices.Prerequisites.Traits

import QuantumLattices.Prerequisites.Traits: parameternames, isparameterbound, contentnames, getcontent, dissolve

abstract type FT{T} end
parameternames(::Type{<:FT}) = (:content,)
isparameterbound(::Type{<:FT}, ::Val{:content}, D) = !isconcretetype(D)
contentnames(::Type{<:FT}) = (:content,)
dissolve(m::FT, ::Val{:content}, f::Function, args::Tuple, kwargs::NamedTuple) = f(getcontent(m, :content), args...; kwargs...)

struct SameNameField <: FT{Vector{Int}}
    info::String
    content::Vector{Int}
end
parameternames(::Type{SameNameField}) = ()
contentnames(::Type{SameNameField}) = (:info, :content)

struct DiffNameField <: FT{Vector{Int}}
    info::String
    table::Vector{Int}
end
parameternames(::Type{DiffNameField}) = ()
contentnames(::Type{DiffNameField}) = (:info, :content)
getcontent(m::DiffNameField, ::Val{:content}) = getfield(m, :table)

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
    @test isparameterbound(Vector, 1, Real) == false
    @test isparameterbound(Vector, 2, 1) == false
    @test isparameterbounds(Vector{<:Real}, Tuple{Real, 1}) == (false, false)
    @test reparameter(Vector, 1, Real, true) == Vector{<:Real}
    @test reparameter(Vector, 1, Real, false) ==  Vector{Real}
    @test reparameter(Vector{Int}, 2, 3, false) == Array{Int, 3}

    @test contentnames(Vector) == ()
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
    @test isparameterbound(FT, :nocotent, Vector) == false
    @test isparameterbound(FT, :content, Vector) == true
    @test isparameterbound(FT, :content, Vector{Int}) == false
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

    a = SameNameField("abcd", [1, 2, 3, 4])
    @test getcontent(a, :content) == getcontent(a, 2) == getcontent(a, Val(:content)) == [1, 2, 3, 4]
    @test dissolve(a, getindex, (1,)) == ("abcd", 1)

    @test contentnames(DiffNameField) == (:info, :content)
    @test hasparameter(DiffNameField, :content) == false

    b = DiffNameField("abcd", [1, 2, 3, 4])
    @test getcontent(b, :content) == getcontent(b, 2) == getcontent(b, Val(:content)) == [1, 2, 3, 4]
    @test dissolve(b, getindex, (2,)) == ("abcd", 2)
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
Base.isapprox(fc1::EFO, fc2::EFO; atol::Real=10^-5, rtol::Real=10^-5) = isapprox(efficientoperations, (:f1, :f2), fc1, fc2; atol=atol, rtol=rtol)
Base.replace(fc::EFO; kwargs...) = replace(efficientoperations, fc; kwargs...)

@testset "efficientoperations" begin
    @test ==(efficientoperations, (), ())
    @test isequal(efficientoperations, (), ())
    @test ==(efficientoperations, (1, 2), (1, 2, 3)) == false
    @test isequal(efficientoperations, (1, 2), (1, 2, 3)) == false
    @test isapprox(efficientoperations, (1.0, 2.0), (1, 2)) == true
    @test isapprox(efficientoperations, (), (), ()) == true
    @test isapprox(efficientoperations, (), (), (1,)) == false

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
