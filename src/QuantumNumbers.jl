module QuantumNumbers

using Base: @propagate_inbounds
using Base.Iterators: product
using DataStructures: OrderedDict
using HalfIntegers: HalfInt
using LinearAlgebra: norm
using Printf: @printf
using QuantumLattices: id, value
using Random: seed!
using ..Toolkit: VectorSpace, VectorSpaceCartesian, VectorSpaceDirectProducted, VectorSpaceDirectSummed, VectorSpaceStyle, efficientoperations, subscript

import ..QuantumLattices: ⊕, ⊗, ⊠, decompose, dimension, rank
import ..Toolkit: shape

export Abelian, AbelianQuantumNumber, AbelianGradedSpace, AbelianGradedSpaceProd, AbelianGradedSpaceSum, CompositeAbelianQuantumNumber, Graded, Momenta, RepresentationSpace, SimpleAbelianQuantumNumber
export Momentum, Momentum₁, Momentum₂, Momentum₃, ℕ, 𝕊ᶻ, 𝕌₁, ℤ, ℤ₂, ℤ₃, ℤ₄, findindex, period, periods, regularize, regularize!

"""
    AbelianQuantumNumber

Abstract type of Abelian quantum numbers.

An Abelian quantum number is the label of a irreducible representation of an Abelian group acted on a quantum representation space. 
"""
abstract type AbelianQuantumNumber end
@inline Base.:(==)(qn₁::QN, qn₂::QN) where {QN<:AbelianQuantumNumber} = ==(efficientoperations, qn₁, qn₂)
@inline Base.isequal(qn₁::QN, qn₂::QN) where {QN<:AbelianQuantumNumber} = isequal(efficientoperations, qn₁, qn₂)
@inline Base.:<(qn₁::QN, qn₂::QN) where {QN<:AbelianQuantumNumber} = <(efficientoperations, qn₁, qn₂)
@inline Base.isless(qn₁::QN, qn₂::QN) where {QN<:AbelianQuantumNumber} = isless(efficientoperations, qn₁, qn₂)
@inline Base.hash(qn::AbelianQuantumNumber, h::UInt) = hash(values(qn), h)
@inline Base.show(io::IO, ::Type{T}) where {T<:AbelianQuantumNumber} = @printf io "%s" nameof(T)
@inline Base.iterate(qn::AbelianQuantumNumber) = (qn, nothing)
@inline Base.iterate(::AbelianQuantumNumber, ::Nothing) = nothing

"""
    const Abelian = AbelianQuantumNumber

Type alias for `AbelianQuantumNumber`.
"""
const Abelian = AbelianQuantumNumber

"""
    getindex(::Type{Abelian}, ::Type{T}) where {T<:AbelianQuantumNumber} -> Type{T}

Overloaded `[]` for `Abelian`, i.e., the support of syntax `Abelian[T]` where `T<:AbelianQuantumNumber`, which is helpful for the construction of tensor producted Abelian quantum numbers.
"""
@inline Base.getindex(::Type{Abelian}, ::Type{T}) where {T<:AbelianQuantumNumber} = T

"""
    zero(qn::AbelianQuantumNumber) -> typeof(qn)
    zero(::Type{QN}) where {QN<:AbelianQuantumNumber} -> QN

Get the zero Abelian quantum number.
"""
@inline Base.zero(qn::AbelianQuantumNumber) = zero(typeof(qn))

"""
    ⊗(qns::AbelianQuantumNumber...) -> eltype(qns)

Get the direct product of some `AbelianQuantumNumber`s.
"""
@inline ⊗(qns::AbelianQuantumNumber...) = +(qns...)

"""
    periods(qn::AbelianQuantumNumber) -> Tuple{Vararg{Number}}
    periods(::Type{QN}) where {QN<:AbelianQuantumNumber} -> Tuple{Vararg{Number}}

The periods of Abelian quantum numbers.
"""
@inline periods(qn::AbelianQuantumNumber) = periods(typeof(qn))

"""
    inv(qn::AbelianQuantumNumber, bool::Bool=true) -> typeof(qn)

Get the inverse of an Abelian quantum number `qn` if `bool` is true. Otherwise, return `qn` itself.
"""
@inline Base.inv(qn::AbelianQuantumNumber, bool::Bool=true) = bool ? -qn : qn

"""
    +(qn::AbelianQuantumNumber) -> typeof(qn)

Overloaded `+` operator for `AbelianQuantumNumber`.
"""
@inline Base.:+(qn::AbelianQuantumNumber) = qn

"""
    SimpleAbelianQuantumNumber <: AbelianQuantumNumber

Abstract type of simple Abelian quantum numbers. That is, it contains only one label.
"""
abstract type SimpleAbelianQuantumNumber <: AbelianQuantumNumber end
@inline Base.values(qn::SimpleAbelianQuantumNumber) = getfield(qn, 1)
@inline Base.show(io::IO, qn::SimpleAbelianQuantumNumber) = @printf io "%s(%s)" typeof(qn) values(qn)
@inline Base.zero(::Type{QN}) where {QN<:SimpleAbelianQuantumNumber} = QN(zero(fieldtype(QN, 1)))
@inline periods(::Type{QN}) where {QN<:SimpleAbelianQuantumNumber} = (period(QN),)

"""
    +(qn₁::QN, qn₂::QN, qns::QN...) where {QN<:SimpleAbelianQuantumNumber} -> QN

Overloaded `+` operator for `SimpleAbelianQuantumNumber`.

!!! note
    To ensure type stability, two `SimpleAbelianQuantumNumber` can be added together if and only if they are of the same type.
"""
@inline Base.:+(qn₁::QN, qn₂::QN, qns::QN...) where {QN<:SimpleAbelianQuantumNumber} = QN(sum(map(values, (qn₁, qn₂, qns...))))

"""
    -(qn::SimpleAbelianQuantumNumber) -> typeof(qn)
    -(qn₁::QN, qn₂::QN) where {QN<:SimpleAbelianQuantumNumber} -> QN

Overloaded `-` operator for `SimpleAbelianQuantumNumber`.

!!! note
    To ensure type stability, a `SimpleAbelianQuantumNumber` can be subtracted by another `SimpleAbelianQuantumNumber` if and only if they are of the same type.
"""
@inline Base.:-(qn::SimpleAbelianQuantumNumber) = typeof(qn)(-values(qn))
@inline Base.:-(qn₁::QN, qn₂::QN) where {QN<:SimpleAbelianQuantumNumber} = QN(values(qn₁)-values(qn₂))

"""
    period(qn::SimpleAbelianQuantumNumber) -> Number
    period(::Type{QN}) where {QN<:SimpleAbelianQuantumNumber} -> Number

The period of a simple Abelian quantum number.
"""
@inline period(qn::SimpleAbelianQuantumNumber) = period(typeof(qn))

"""
    𝕌₁ <: SimpleAbelianQuantumNumber

Abstract type of 𝕌₁ quantum numbers.
"""
abstract type 𝕌₁ <: SimpleAbelianQuantumNumber end
@inline period(::Type{<:𝕌₁}) = Inf

"""
    𝕊ᶻ <: 𝕌₁

Concrete Abelian quantum number of the z-component of a spin.
"""
struct 𝕊ᶻ <: 𝕌₁
    charge::HalfInt
end

"""
    ℕ <: <: 𝕌₁

Concrete Abelian quantum number of the particle number.
"""
struct ℕ <: 𝕌₁
    charge::Int
end

"""
    ℤ{N} <: SimpleAbelianQuantumNumber

ℤₙ quantum numbers.
"""
struct ℤ{N} <: SimpleAbelianQuantumNumber
    charge::Int
    function ℤ{N}(charge::Integer) where N
        @assert N>0 "ℤ error: non-positive period ($N)."
        new{N}(mod(charge, N))
    end
end
@inline period(::Type{ℤ{N}}) where N = N
@inline Base.show(io::IO, ::Type{ℤ{N}}) where N = @printf io "ℤ%s" N<5 ? subscript(N) : string("{", N, "}")

"""
    const ℤ₂ = ℤ{2}
    const ℤ₃ = ℤ{3}
    const ℤ₄ = ℤ{4}

Alias for ℤ₂/ℤ₃/ℤ₄ quantum numbers.
"""
const ℤ₂ = ℤ{2}
const ℤ₃ = ℤ{3}
const ℤ₄ = ℤ{4}

"""
    CompositeAbelianQuantumNumber{T<:Tuple{Vararg{SimpleAbelianQuantumNumber}}} <: AbelianQuantumNumber

Deligne tensor product of simple Abelian quantum numbers.
"""
struct CompositeAbelianQuantumNumber{T<:Tuple{Vararg{SimpleAbelianQuantumNumber}}} <: AbelianQuantumNumber
    contents::T
end
@inline Base.values(qn::CompositeAbelianQuantumNumber) = map(values, qn.contents)
@inline Base.show(io::IO, qn::CompositeAbelianQuantumNumber) = @printf io "Abelian[%s]%s" join(fieldtypes(fieldtype(typeof(qn), :contents)), " ⊠ ") values(qn)
@inline Base.zero(::Type{CompositeAbelianQuantumNumber{T}}) where {T<:Tuple{Vararg{SimpleAbelianQuantumNumber}}} = CompositeAbelianQuantumNumber(map(zero,  fieldtypes(T)))
@inline Base.length(qn::CompositeAbelianQuantumNumber) = length(typeof(qn))
@inline Base.length(::Type{<:CompositeAbelianQuantumNumber{T}}) where {T<:Tuple{Vararg{SimpleAbelianQuantumNumber}}} = fieldcount(T)
@inline Base.firstindex(::CompositeAbelianQuantumNumber) = 1
@inline Base.lastindex(qn::CompositeAbelianQuantumNumber) = length(qn)
@inline Base.getindex(qn::CompositeAbelianQuantumNumber, i::Integer) = qn.contents[i]
@inline periods(::Type{CompositeAbelianQuantumNumber{T}}) where {T<:Tuple{Vararg{SimpleAbelianQuantumNumber}}} = map(period, fieldtypes(T))
@inline CompositeAbelianQuantumNumber(contents::SimpleAbelianQuantumNumber...) = CompositeAbelianQuantumNumber(contents)

"""
    CompositeAbelianQuantumNumber{T}(vs::Vararg{Number, N}) where {N, T<:NTuple{N, SimpleAbelianQuantumNumber}}
    CompositeAbelianQuantumNumber{T}(vs::NTuple{N, Number}) where {N, T<:NTuple{N, SimpleAbelianQuantumNumber}}

Construct a Deligne tensor product of simple Abelian quantum numbers by their values.
"""
@inline CompositeAbelianQuantumNumber{T}(vs::Vararg{Number, N}) where {N, T<:NTuple{N, SimpleAbelianQuantumNumber}} = CompositeAbelianQuantumNumber{T}(vs)
@inline CompositeAbelianQuantumNumber{T}(vs::NTuple{N, Number}) where {N, T<:NTuple{N, SimpleAbelianQuantumNumber}} = CompositeAbelianQuantumNumber(map((T, v)->T(v), fieldtypes(T), vs))

"""
    period(qn::CompositeAbelianQuantumNumber, i::Integer) -> Number
    period(::Type{CompositeAbelianQuantumNumber{T}}, i::Integer) where {T<:Tuple{Vararg{SimpleAbelianQuantumNumber}}} -> Number

Get the period of the ith simple Abelian number contained in a Deligne tensor product.
"""
@inline period(qn::CompositeAbelianQuantumNumber, i::Integer) = period(typeof(qn), i)
@inline period(::Type{CompositeAbelianQuantumNumber{T}}, i::Integer) where {T<:Tuple{Vararg{SimpleAbelianQuantumNumber}}} = period(fieldtype(T, i))

"""
    +(qn₁::QN, qn₂::QN, qns::QN...) where {QN<:CompositeAbelianQuantumNumber} -> QN

Overloaded `+` operator for `CompositeAbelianQuantumNumber`.

!!! note
    To ensure type stability, two `CompositeAbelianQuantumNumber` can be added together if and only if they are of the same type.
"""
@inline Base.:+(qn₁::QN, qn₂::QN, qns::QN...) where {QN<:CompositeAbelianQuantumNumber} = QN(mapreduce(values, .+, (qn₁, qn₂, qns...)))

"""
    -(qn::CompositeAbelianQuantumNumber) -> typeof(qn)
    -(qn₁::QN, qn₂::QN) where {QN<:CompositeAbelianQuantumNumber} -> QN

Overloaded `-` operator for `CompositeAbelianQuantumNumber`.

!!! note
    To ensure type stability, a `CompositeAbelianQuantumNumber` can be subtracted by another `CompositeAbelianQuantumNumber` if and only if they are of the same type.
"""
@inline Base.:-(qn::CompositeAbelianQuantumNumber) = CompositeAbelianQuantumNumber(map(-, qn.contents))
@inline Base.:-(qn₁::QN, qn₂::QN) where {QN<:CompositeAbelianQuantumNumber} = QN(map((i₁, i₂)->i₁-i₂, qn₁.contents, qn₂.contents))

"""
    ⊠(qns::SimpleAbelianQuantumNumber...) -> CompositeAbelianQuantumNumber
    ⊠(qn₁::SimpleAbelianQuantumNumber, qn₂::CompositeAbelianQuantumNumber) -> CompositeAbelianQuantumNumber
    ⊠(qn₁::CompositeAbelianQuantumNumber, qn₂::SimpleAbelianQuantumNumber) -> CompositeAbelianQuantumNumber
    ⊠(qn₁::CompositeAbelianQuantumNumber, qn₂::CompositeAbelianQuantumNumber) -> CompositeAbelianQuantumNumber

Deligne tensor product of Abelian quantum numbers.
"""
@inline ⊠(qns::SimpleAbelianQuantumNumber...) = CompositeAbelianQuantumNumber(qns...)
@inline ⊠(qn₁::SimpleAbelianQuantumNumber, qn₂::CompositeAbelianQuantumNumber) = CompositeAbelianQuantumNumber(qn₁, qn₂.contents...)
@inline ⊠(qn₁::CompositeAbelianQuantumNumber, qn₂::SimpleAbelianQuantumNumber) = CompositeAbelianQuantumNumber(qn₁.contents..., qn₂)
@inline ⊠(qn₁::CompositeAbelianQuantumNumber, qn₂::CompositeAbelianQuantumNumber) = CompositeAbelianQuantumNumber(qn₁.contents..., qn₂.contents...)

"""
    ⊠(QNS::Type{<:SimpleAbelianQuantumNumber}...) -> Type{CompositeAbelianQuantumNumber{Tuple{QNS...}}}
    ⊠(::Type{QN}, ::Type{CompositeAbelianQuantumNumber{T}}) where {QN<:SimpleAbelianQuantumNumber, T<:Tuple{Vararg{SimpleAbelianQuantumNumber}}} -> Type{CompositeAbelianQuantumNumber{Tuple{QN, fieldtypes(T)...}}}
    ⊠(::Type{CompositeAbelianQuantumNumber{T}}, ::Type{QN}) where {T<:Tuple{Vararg{SimpleAbelianQuantumNumber}}, QN<:SimpleAbelianQuantumNumber} -> Type{CompositeAbelianQuantumNumber{Tuple{fieldtypes(T)...}, QN}}
    ⊠(::Type{CompositeAbelianQuantumNumber{T₁}}, ::Type{CompositeAbelianQuantumNumber{T₂}}) where {T₁<:Tuple{Vararg{SimpleAbelianQuantumNumber}}, T₂<:Tuple{Vararg{SimpleAbelianQuantumNumber}}} -> Type{CompositeAbelianQuantumNumber{Tuple{fieldtypes(T₁)..., fieldtypes(T₂)...}}}

Deligne tensor product of Abelian quantum numbers.
"""
@inline ⊠(QNS::Type{<:SimpleAbelianQuantumNumber}...) = CompositeAbelianQuantumNumber{Tuple{QNS...}}
@inline ⊠(::Type{QN}, ::Type{CompositeAbelianQuantumNumber{T}}) where {QN<:SimpleAbelianQuantumNumber, T<:Tuple{Vararg{SimpleAbelianQuantumNumber}}} = CompositeAbelianQuantumNumber{Tuple{QN, fieldtypes(T)...}}
@inline ⊠(::Type{CompositeAbelianQuantumNumber{T}}, ::Type{QN}) where {T<:Tuple{Vararg{SimpleAbelianQuantumNumber}}, QN<:SimpleAbelianQuantumNumber} = CompositeAbelianQuantumNumber{Tuple{fieldtypes(T)..., QN}}
@inline ⊠(::Type{CompositeAbelianQuantumNumber{T₁}}, ::Type{CompositeAbelianQuantumNumber{T₂}}) where {T₁<:Tuple{Vararg{SimpleAbelianQuantumNumber}}, T₂<:Tuple{Vararg{SimpleAbelianQuantumNumber}}} = CompositeAbelianQuantumNumber{Tuple{fieldtypes(T₁)..., fieldtypes(T₂)...}}

"""
    const Momentum = CompositeAbelianQuantumNumber{<:Tuple{Vararg{ℤ}}}
    const Momentum₁{N} = CompositeAbelianQuantumNumber{Tuple{ℤ{N}}}
    const Momentum₂{N₁, N₂} = CompositeAbelianQuantumNumber{Tuple{ℤ{N₁}, ℤ{N₂}}}
    const Momentum₃{N₁, N₂, N₃} = CompositeAbelianQuantumNumber{Tuple{ℤ{N₁}, ℤ{N₂}, ℤ{N₃}}}

Type alias for the Abelian quantum numbers of 1d, 2d and 3d momentum.
"""
const Momentum = CompositeAbelianQuantumNumber{<:Tuple{Vararg{ℤ}}}
const Momentum₁{N} = CompositeAbelianQuantumNumber{Tuple{ℤ{N}}}
const Momentum₂{N₁, N₂} = CompositeAbelianQuantumNumber{Tuple{ℤ{N₁}, ℤ{N₂}}}
const Momentum₃{N₁, N₂, N₃} = CompositeAbelianQuantumNumber{Tuple{ℤ{N₁}, ℤ{N₂}, ℤ{N₃}}}
@inline Int(m::Momentum₁) = m[1].charge + 1
@inline Int(m::Momentum₂{N₁, N₂}) where {N₁, N₂} = m[2].charge + m[1].charge*N₂ + 1
@inline Int(m::Momentum₃{N₁, N₂, N₃}) where {N₁, N₂, N₃} = (m[1].charge*N₂+m[2].charge)*N₃ + m[3].charge + 1
@inline Base.show(io::IO, m::Momentum₁{N}) where N = @printf io "Momentum₁{%s}(%s)" N m[1].charge
@inline Base.show(io::IO, m::Momentum₂{N₁, N₂}) where {N₁, N₂} = @printf io "Momentum₂{%s, %s}(%s, %s)" N₁ N₂ m[1].charge m[2].charge
@inline Base.show(io::IO, m::Momentum₃{N₁, N₂, N₃}) where {N₁, N₂, N₃} = @printf io "Momentum₃{%s, %s, %s}(%s, %s, %s)" N₁ N₂ N₃ m[1].charge m[2].charge m[3].charge
@inline Base.show(io::IO, ::Type{Momentum₁{N}}) where N = @printf io "Momentum₁{%s}" N
@inline Base.show(io::IO, ::Type{Momentum₂{N₁, N₂}}) where {N₁, N₂} = @printf io "Momentum₂{%s, %s}" N₁ N₂
@inline Base.show(io::IO, ::Type{Momentum₃{N₁, N₂, N₃}}) where {N₁, N₂, N₃} = @printf io "Momentum₃{%s, %s, %s}" N₁ N₂ N₃

"""
    Momentum₁{N}(k::Integer) where N
    Momentum₂{N}(k₁::Integer, k₂::Integer) where N
    Momentum₂{N₁, N₂}(k₁::Integer, k₂::Integer) where {N₁, N₂}
    Momentum₃{N}(k₁::Integer, k₂::Integer, k₃::Integer) where N
    Momentum₃{N₁, N₂, N₃}(k₁::Integer, k₂::Integer, k₃::Integer) where {N₁, N₂, N₃}

Construct 1d, 2d and 3d momentum.
"""
@inline Momentum₁{N}(k::Integer) where N = CompositeAbelianQuantumNumber(ℤ{N}(k))
@inline Momentum₂{N}(k₁::Integer, k₂::Integer) where N = Momentum₂{N, N}(k₁, k₂)
@inline Momentum₂{N₁, N₂}(k₁::Integer, k₂::Integer) where {N₁, N₂} = CompositeAbelianQuantumNumber(ℤ{N₁}(k₁), ℤ{N₂}(k₂))
@inline Momentum₃{N}(k₁::Integer, k₂::Integer, k₃::Integer) where N =  Momentum₃{N, N, N}(k₁, k₂, k₃)
@inline Momentum₃{N₁, N₂, N₃}(k₁::Integer, k₂::Integer, k₃::Integer) where {N₁, N₂, N₃} = CompositeAbelianQuantumNumber(ℤ{N₁}(k₁), ℤ{N₂}(k₂), ℤ{N₃}(k₃))

"""
    RepresentationSpace{QN<:AbelianQuantumNumber} <: VectorSpace{QN}

Abstract type of quantum representation spaces of Abelian groups.
"""
abstract type RepresentationSpace{QN<:AbelianQuantumNumber} <: VectorSpace{QN} end
@inline Base.show(io::IO, ::MIME"text/plain", rs::RepresentationSpace) = show(io, rs)
@inline dimension(::Type{<:RepresentationSpace}) = Int
@inline Base.range(::Type{<:RepresentationSpace}) = UnitRange{Int}

"""
    dimension(rs::RepresentationSpace, i::Integer) -> Int

Get the degenerate dimension of the ith Abelian quantum number contained in a representation space.
"""
@inline dimension(rs::RepresentationSpace, i::Integer) = dimension(rs, CartesianIndex(i, rs)) 

"""
    range(rs::RepresentationSpace, i::Integer)

Get the slice of the degenerate dimension of the ith Abelian quantum number contained in a representation space.
"""
@inline Base.range(rs::RepresentationSpace, i::Integer) = range(rs, CartesianIndex(i, rs))

"""
    cumsum(rs::RepresentationSpace, i::Union{Integer, CartesianIndex}) -> Int

Get the accumulative degenerate dimension up to the ith Abelian quantum number contained in a representation space.
"""
@inline Base.cumsum(rs::RepresentationSpace, i::Union{Integer, CartesianIndex}) = iszero(i) ? 0 : range(rs, i).stop

"""
    pairs(rs::RepresentationSpace, ::typeof(dimension)) -> RepresentationSpacePairs
    pairs(rs::RepresentationSpace, ::typeof(range)) -> RepresentationSpacePairs

Return an iterator that iterates over the pairs of the Abelian quantum numbers and their corresponding (slices of the) degenerate dimensions contained in a representation space.
"""
@inline Base.pairs(rs::RepresentationSpace, ::typeof(dimension)) = RepresentationSpacePairs{dimension(typeof(rs))}(rs)
@inline Base.pairs(rs::RepresentationSpace, ::typeof(range)) = RepresentationSpacePairs{range(typeof(rs))}(rs)
struct RepresentationSpacePairs{R<:Union{Int, AbstractVector{Int}}, QN<:AbelianQuantumNumber, S<:RepresentationSpace{QN}} <: AbstractVector{Pair{QN, R}}
    space::S
    RepresentationSpacePairs{R}(space::RepresentationSpace) where {R<:Union{Int, AbstractVector{Int}}} = new{R, eltype(space), typeof(space)}(space)
end
@inline Base.size(ps::RepresentationSpacePairs) = (length(ps.space),)
@inline function Base.getindex(ps::RepresentationSpacePairs{Int}, i::Integer)
    index = CartesianIndex(i, ps.space)
    return ps.space[index]=>dimension(ps.space, index)
end
@inline function Base.getindex(ps::RepresentationSpacePairs{<:AbstractVector{Int}}, i::Integer)
    index = CartesianIndex(i, ps.space)
    return ps.space[index]=>range(ps.space, index)
end

"""
    Momenta{P<:Momentum} <: RepresentationSpace{P}

The complete allowed set of momenta.
"""
struct Momenta{P<:Momentum} <: RepresentationSpace{P} end
@inline Momenta(::Type{P}) where {P<:Momentum} = Momenta{P}()
@inline VectorSpaceStyle(::Type{<:Momenta}) = VectorSpaceCartesian()
@inline shape(::Momenta{P}) where {P<:Momentum} = map(period->0:period-1, reverse(periods(P)))
@inline Base.convert(::Type{<:CartesianIndex}, m::P, ::Momenta{P}) where {P<:Momentum} = CartesianIndex(reverse(values(m)))
@inline Base.convert(::Type{P}, index::CartesianIndex, ::Momenta{P}) where {P<:Momentum} = P(reverse(index.I)...)
@inline Base.:(==)(ms₁::Momenta, ms₂::Momenta) = periods(eltype(ms₁))==periods(eltype(ms₂))
@inline Base.isequal(ms₁::Momenta, ms₂::Momenta) = isequal(periods(eltype(ms₁)), periods(eltype(ms₂)))

"""
    regularize!(quantumnumbers::Vector{<:AbelianQuantumNumber}, dimensions::Vector{Int}; check::Bool=false) -> Tuple{typeof(quantumnumbers), typeof(dimensions), Vector{Int}}

In place regularization of the input Abelian quantum numbers and their corresponding degenerate dimensions.

After the regularization, the Abelian quantum numbers will be sorted in the ascending order and duplicates will be merged together.The degenerate dimensions will be processed accordingly. When `check` is `true`, this function also check whether all input degenerate dimensions are positive. The regularized Abelian quantum numbers and degenerate dimensions, as well as the permutation vector that sorts the input Abelian quantum numbers, will be returned. 
"""
function regularize!(quantumnumbers::Vector{<:AbelianQuantumNumber}, dimensions::Vector{Int}; check::Bool=false)
    check && @assert all(>(0), dimensions) "regularize! error: input dimensions must be positive integers."
    perm = sortperm(quantumnumbers)
    permute!(quantumnumbers, perm)
    permute!(dimensions, perm)
    count = 1
    x = first(quantumnumbers)
    for i = 2:length(perm)
        y = quantumnumbers[i]
        d = dimensions[i]
        if y ≠ x
            count += 1
            quantumnumbers[count] = y
            dimensions[count] = d
            x = y
        else
            dimensions[count] += d
        end
    end
    resize!(quantumnumbers, count)
    resize!(dimensions, count)
    return quantumnumbers, dimensions, perm
end

"""
    regularize(quantumnumbers::AbstractVector{<:AbelianQuantumNumber}, dimension::AbstractVector{<:Integer}; check::Bool=false) -> Tuple{Vector{eltype(quantumnumbers)}, Vector{Int}, Vector{Int}}

Regularize of the input Abelian quantum numbers and their corresponding degenerate dimensions.

See [`regularize!`](@ref).
"""
@inline function regularize(quantumnumbers::AbstractVector{<:AbelianQuantumNumber}, dimension::AbstractVector{<:Integer}; check::Bool=false)
    return regularize!(collect(quantumnumbers), collect(Int, dimension); check=check)
end

"""
    AbelianGradedSpace{QN<:AbelianQuantumNumber} <: RepresentationSpace{QN}

A quantum representation space of an Abelian group that has been decomposed into the direct sum of its irreducible representations.
"""
struct AbelianGradedSpace{QN<:AbelianQuantumNumber} <: RepresentationSpace{QN}
    contents::OrderedDict{QN, UnitRange{Int}}
    dual::Bool
end
@inline Base.getindex(::VectorSpaceStyle, gs::AbelianGradedSpace, i::CartesianIndex) = inv(id(gs.contents, i), gs.dual)
function Base.show(io::IO, gs::AbelianGradedSpace)
    @printf io "Graded{%s}(" eltype(gs)
    for (i, qn) in enumerate(keys(gs.contents))
        @printf io "%s=>%s" values(qn) dimension(gs, i)
        i<length(gs) && @printf io "%s" ", "
    end
    @printf io "%s%s" ")" gs.dual ? "'" : ""
end

"""
    const Graded = AbelianGradedSpace

Type alias for `AbelianGradedSpace`.
"""
const Graded = AbelianGradedSpace

"""
    AbelianGradedSpace(quantumnumbers::AbstractVector{<:AbelianQuantumNumber}, dimensions::AbstractVector{<:Integer}, dual::Bool=false; ordercheck::Bool=false, duplicatecheck::Bool=false, degeneracycheck::Bool=false)

Construct an Abelian graded space.

Here:
- `quantumnumbers` specifies the Abelian quantum numbers labeling the irreducible representations of the corresponding Abelian group which must be sorted in the ascending order. Such an ordering should be manually guaranteed by the user. When and only when the keyword argument `ordercheck` is `true`, the constructor will check whether this condition is satisfied and raise an error if it doesn't. Besides, `quantumnumbers` must not contain duplicate Abelian quantum numbers, manually guaranteed by the user as well. This condition can be checked when and only when both `ordercheck==true` and `duplicatecheck==true`. An error will be raised if this check fails.
- `dimensions` specifies the degenerate dimensions of the corresponding Abelian quantum numbers. Apparently, each degenerate dimension must be positive, which should also be manually guaranteed by the user. When and only when the keyword argument `degeneracycheck` is `true`, the constructor will check whether this condition is satisfied and raise an error if it doesn't.
- `dual` specifies whether the graded space is the dual representation of the corresponding Abelian group, which roughly speaking, can be viewed as the direction of the arrow of the Abelian quantum numbers labeling the irreducible representations. We assume `dual==true` corresponds to the in-arrow and `dual==false` corresponds to the out-arrow.
"""
function AbelianGradedSpace(quantumnumbers::AbstractVector{<:AbelianQuantumNumber}, dimensions::AbstractVector{<:Integer}, dual::Bool=false; ordercheck::Bool=false, duplicatecheck::Bool=false, degeneracycheck::Bool=false)
    @assert length(quantumnumbers)==length(dimensions) "AbelianGradedSpace error: mismatched Abelian quantum numbers and accumulative degenerate dimensions."
    if ordercheck
        @assert issorted(quantumnumbers) "AbelianGradedSpace error: Abelian quantum numbers are not sorted."
        duplicatecheck && @assert all(i->quantumnumbers[i]≠quantumnumbers[i+1], firstindex(quantumnumbers):lastindex(quantumnumbers)-1) "AbelianGradedSpace error: duplicate Abelian quantum numbers found."
    end
    degeneracycheck &&  @assert all(>(0), dimensions) "AbelianGradedSpace error: non-positive degenerate dimensions."
    contents = OrderedDict{eltype(quantumnumbers), UnitRange{Int}}()
    start, stop = 1, 0
    @inbounds for (quantumnumber, dimension) in zip(quantumnumbers, dimensions)
        stop += dimension
        contents[quantumnumber] = start:stop
        start += dimension
    end
    return AbelianGradedSpace(contents, dual)
end

"""
    AbelianGradedSpace(pairs; dual::Bool=false)
    AbelianGradedSpace(pairs::Pair...; dual::Bool=false)
    AbelianGradedSpace{QN}(pairs; dual::Bool=false) where {QN<:AbelianQuantumNumber}
    AbelianGradedSpace{QN}(pairs::Pair...; dual::Bool=false) where {QN<:AbelianQuantumNumber}

Construct an Abelian graded space.

In this function, the Abelian quantum numbers will be sorted automatically, therefore, their orders need not be worried. Duplicate and dimension checks of the quantum numbers are also carried out and errors will be raised if either such checks fails.
"""
@inline AbelianGradedSpace(pairs; dual::Bool=false) = AbelianGradedSpace{fieldtype(eltype(pairs), 1)}(pairs; dual=dual)
@inline AbelianGradedSpace(pairs::Pair...; dual::Bool=false) = AbelianGradedSpace(pairs; dual=dual)
@propagate_inbounds function AbelianGradedSpace{QN}(pairs; dual::Bool=false) where {QN<:AbelianQuantumNumber}
    quantumnumbers, dimensions = zeros(QN, length(pairs)), zeros(Int, length(pairs))
    for (i, (qn, dimension)) in enumerate(pairs)
        quantumnumbers[i] = isa(qn, QN) ? qn : QN(qn)
        dimensions[i] = dimension
    end
    quantumnumbers, dimensions = regularize!(quantumnumbers, dimensions; check=true)
    return AbelianGradedSpace(quantumnumbers, dimensions, dual; ordercheck=false, duplicatecheck=false, degeneracycheck=false)
end
@inline AbelianGradedSpace{QN}(pairs::Pair...; dual::Bool=false) where {QN<:AbelianQuantumNumber} = AbelianGradedSpace{QN}(pairs; dual=dual)

"""
    length(gs::AbelianGradedSpace) -> Int

Get the number of inequivalent irreducible representations (i.e., the Abelian quantum numbers) of an Abelian graded space.
"""
@inline Base.length(gs::AbelianGradedSpace) = length(gs.contents)

"""
    getindex(gs::AbelianGradedSpace, indexes::AbstractVector{<:Integer}) -> typeof(gs)
    getindex(gs::AbelianGradedSpace{QN}, quantumnumbers::AbstractVector{QN}) where {QN<:AbelianQuantumNumber} -> AbelianGradedSpace{QN}

Get a subset of an Abelian graded space.
"""
function Base.getindex(gs::AbelianGradedSpace, indexes::AbstractVector{<:Integer})
    issorted(indexes) || (indexes = sort(indexes))
    quantumnumbers, dimensions = zeros(eltype(gs), length(indexes)), zeros(Int, length(indexes))
    for (i, index) in enumerate(indexes)
        quantumnumbers[i] = id(gs.contents, index)
        dimensions[i] = dimension(gs, index)
    end
    return AbelianGradedSpace(quantumnumbers, dimensions, gs.dual; ordercheck=false, duplicatecheck=false, degeneracycheck=false)
end
function Base.getindex(gs::AbelianGradedSpace{QN}, quantumnumbers::AbstractVector{QN}) where {QN<:AbelianQuantumNumber}
    quantumnumbers = [inv(qn, gs.dual) for qn in quantumnumbers]
    sort!(quantumnumbers)
    dimensions = zeros(Int, length(quantumnumbers))
    @inbounds for (i, qn) in enumerate(quantumnumbers)
        dimensions[i] =  dimension(gs, inv(qn, gs.dual))
    end
    return AbelianGradedSpace(quantumnumbers, dimensions, gs.dual; ordercheck=false, duplicatecheck=false, degeneracycheck=false)
end

"""
    in(qn::QN, gs::AbelianGradedSpace{QN}) where {QN<:AbelianQuantumNumber} -> Bool

Check whether an Abelian quantum number is contained in an Abelian graded space.
"""
@inline Base.in(qn::QN, gs::AbelianGradedSpace{QN}) where {QN<:AbelianQuantumNumber} = haskey(gs.contents, inv(qn, gs.dual))

"""
    adjoint(gs::AbelianGradedSpace) -> typeof(gs)

Get the dual of an Abelian graded space.
"""
@inline Base.adjoint(gs::AbelianGradedSpace) = AbelianGradedSpace(gs.contents, !gs.dual)

"""
    dimension(gs::AbelianGradedSpace) -> Int

Get the total dimension of an Abelian graded space.
"""
@inline dimension(gs::AbelianGradedSpace) = last(gs.contents).second.stop

"""
    dimension(gs::AbelianGradedSpace, qn::CartesianIndex) -> Int
    dimension(gs::AbelianGradedSpace{QN}, qn::QN) where {QN<:AbelianQuantumNumber} -> Int

Get the degenerate dimension of an Abelian quantum number contained in an Abelian graded space.
"""
@inline dimension(gs::AbelianGradedSpace, qn::CartesianIndex) = length(value(gs.contents, qn))
@inline dimension(gs::AbelianGradedSpace{QN}, qn::QN) where {QN<:AbelianQuantumNumber} = length(gs.contents[inv(qn, gs.dual)])

"""
    range(gs::AbelianGradedSpace, qn::CartesianIndex) -> UnitRange{Int}
    range(gs::AbelianGradedSpace{QN}, qn::QN) where {QN<:AbelianQuantumNumber} -> UnitRange{Int}

Get the slice of the degenerate dimension of an Abelian quantum number contained in an Abelian graded space.
"""
@inline Base.range(gs::AbelianGradedSpace, qn::CartesianIndex) = value(gs.contents, qn)
@inline Base.range(gs::AbelianGradedSpace{QN}, qn::QN) where {QN<:AbelianQuantumNumber} = gs.contents[inv(qn, gs.dual)]

"""
    cumsum(gs::AbelianGradedSpace{QN}, qn::QN) where {QN<:AbelianQuantumNumber} -> Int

Get the accumulative dimension of an Abelian graded space up to a certain Abelian quantum number contained in an Abelian graded space.
"""
@inline Base.cumsum(gs::AbelianGradedSpace{QN}, qn::QN) where {QN<:AbelianQuantumNumber} = gs.contents[inv(qn, gs.dual)].stop

"""
    findindex(position::Integer, gs::AbelianGradedSpace, guess::Integer) -> Int

Find the index of an Abelian quantum number in an Abelian graded space beginning at `guess` whose position in the complete dimension range is `position`.
"""
@inline findindex(position::Integer, gs::AbelianGradedSpace, guess::Integer) = findnext(i->cumsum(gs, i)>position-1, Base.OneTo(length(gs)), guess)

"""
    CompositeAbelianGradedSpace{N, QN<:AbelianQuantumNumber} <: RepresentationSpace{QN}

Abstract type of composite Abelian graded spaces.
"""
abstract type CompositeAbelianGradedSpace{N, QN<:AbelianQuantumNumber} <: RepresentationSpace{QN} end

"""
    rank(rs::CompositeAbelianGradedSpace) -> Int
    rank(::Type{<:CompositeAbelianGradedSpace{N}}) where N -> Int

Get the number of Abelian graded spaces in the direct sum.
"""
@inline rank(rs::CompositeAbelianGradedSpace) = rank(typeof(rs))
@inline rank(::Type{<:CompositeAbelianGradedSpace{N}}) where N = N

"""
    decompose(rs::CompositeAbelianGradedSpace; expand::Bool=true) -> Tuple{AbelianGradedSpace{eltype(rs)}, Vector{Int}}

Decompose a composite of several Abelian graded spaces to the canonical one.

When `expand` is `false`, the corresponding permutation vector that sorts the Abelian quantum numbers will be returned as well.
When `expand` is `true`, the expanded dimension indexes of the permutation vector that sorts the Abelian quantum numbers will be returned as well.
"""
@propagate_inbounds function decompose(rs::CompositeAbelianGradedSpace; expand::Bool=true)
    len = length(rs)
    quantumnumbers, dimensions = zeros(eltype(rs), len), zeros(Int, len)
    for (i, (qn, dim)) in enumerate(pairs(rs, dimension))
        quantumnumbers[i] = qn
        dimensions[i] = dim
    end
    quantumnumbers, dimensions, perm = regularize!(quantumnumbers, dimensions)
    if expand
        expansion = Int[]
        for index in perm
            append!(expansion, range(rs, index))
        end
        perm = expansion
    end
    return AbelianGradedSpace(quantumnumbers, dimensions, false; ordercheck=false, duplicatecheck=false, degeneracycheck=false), perm
end

"""
    AbelianGradedSpaceSum{N, QN<:AbelianQuantumNumber} <: CompositeAbelianGradedSpace{N, QN}

Direct sum of Abelian graded spaces.
"""
struct AbelianGradedSpaceSum{N, QN<:AbelianQuantumNumber} <: CompositeAbelianGradedSpace{N, QN}
    contents::NTuple{N, AbelianGradedSpace{QN}}
end
@inline AbelianGradedSpaceSum(gses::AbelianGradedSpace...) = AbelianGradedSpaceSum(gses)
@inline VectorSpaceStyle(::Type{<:AbelianGradedSpaceSum}) = VectorSpaceDirectSummed()
function Base.show(io::IO, rs::AbelianGradedSpaceSum)
    for (count, content) in enumerate(rs.contents)
        @printf io "%s" content
        count<rank(rs) && @printf io "%s" " ⊕ "
    end
end

"""
    ⊕(gses::AbelianGradedSpace...) -> AbelianGradedSpaceSum
    ⊕(gs::AbelianGradedSpace, rs::AbelianGradedSpaceSum) -> AbelianGradedSpaceSum
    ⊕(rs::AbelianGradedSpaceSum, gs::AbelianGradedSpace) -> AbelianGradedSpaceSum
    ⊕(rs₁::AbelianGradedSpaceSum, rs₂::AbelianGradedSpaceSum) -> AbelianGradedSpaceSum

Get the direct sum of some Abelian graded spaces.
"""
@inline ⊕(gses::AbelianGradedSpace...) = AbelianGradedSpaceSum(gses...)
@inline ⊕(gs::AbelianGradedSpace, rs::AbelianGradedSpaceSum) = AbelianGradedSpaceSum(gs, rs.contents...)
@inline ⊕(rs::AbelianGradedSpaceSum, gs::AbelianGradedSpace) = AbelianGradedSpaceSum(rs.contents..., gs)
@inline ⊕(rs₁::AbelianGradedSpaceSum, rs₂::AbelianGradedSpaceSum) = AbelianGradedSpaceSum(rs₁.contents..., rs₂.contents...)

"""
    dimension(rs::AbelianGradedSpaceSum) -> Int

Get the total dimension of the direct sum of several Abelian graded spaces.
"""
@inline dimension(rs::AbelianGradedSpaceSum) = sum(dimension, rs.contents)

"""
    dimension(rs::AbelianGradedSpaceSum, i::CartesianIndex) -> Int

Get the degenerate dimension of the ith Abelian quantum number in the direct sum of several Abelian graded spaces.
"""
@inline function dimension(rs::AbelianGradedSpaceSum, i::CartesianIndex)
    m, n = i.I
    return dimension(rs.contents[m], n)
end

"""
    range(rs::AbelianGradedSpaceSum, i::CartesianIndex) -> UnitRange{Int}

Get the slice of the degenerate dimension of the ith Abelian quantum number in the direct sum of several Abelian graded spaces.
"""
@inline function Base.range(rs::AbelianGradedSpaceSum, i::CartesianIndex)
    m, n = i.I
    d = sum(i->dimension(rs.contents[i]), 1:m-1; init=0)
    slice = range(rs.contents[m], n)
    return (slice.start+d):(slice.stop+d)
end

"""
    AbelianGradedSpaceProd{N, QN<:AbelianQuantumNumber} <: CompositeAbelianGradedSpace{N, QN}

Direct product of Abelian graded spaces.
"""
struct AbelianGradedSpaceProd{N, QN<:AbelianQuantumNumber} <: CompositeAbelianGradedSpace{N, QN}
    contents::NTuple{N, AbelianGradedSpace{QN}}
end
@inline AbelianGradedSpaceProd(gses::AbelianGradedSpace...) = AbelianGradedSpaceProd(gses)
@inline VectorSpaceStyle(::Type{<:AbelianGradedSpaceProd}) = VectorSpaceDirectProducted(:backward)
function Base.show(io::IO, rs::AbelianGradedSpaceProd)
    for (count, content) in enumerate(rs.contents)
        @printf io "%s" content
        count<rank(rs) && @printf io "%s" " ⊗ "
    end
end
@inline Base.convert(::Type{QN}, qns::NTuple{N, QN}, rs::AbelianGradedSpaceProd{N, QN}) where {N, QN<:AbelianQuantumNumber} = ⊗(qns...)
@inline Base.range(::Type{<:AbelianGradedSpaceProd{N}}) where N = DirectProductedAbelianGradedSpaceRange{N}

"""
    ⊗(gses::AbelianGradedSpace...) -> AbelianGradedSpaceProd
    ⊗(gs::AbelianGradedSpace, rs::AbelianGradedSpaceProd) -> AbelianGradedSpaceProd
    ⊗(rs::AbelianGradedSpaceProd, gs::AbelianGradedSpace) -> AbelianGradedSpaceProd
    ⊗(rs₁::AbelianGradedSpaceProd, rs₂::AbelianGradedSpaceProd) -> AbelianGradedSpaceProd

Get the direct product of some Abelian graded spaces.
"""
@inline ⊗(gses::AbelianGradedSpace...) = AbelianGradedSpaceProd(gses...)
@inline ⊗(gs::AbelianGradedSpace, rs::AbelianGradedSpaceProd) = AbelianGradedSpaceProd(gs, rs.contents...)
@inline ⊗(rs::AbelianGradedSpaceProd, gs::AbelianGradedSpace) = AbelianGradedSpaceProd(rs.contents..., gs)
@inline ⊗(rs₁::AbelianGradedSpaceProd, rs₂::AbelianGradedSpaceProd) = AbelianGradedSpaceProd(rs₁.contents..., rs₂.contents...)

"""
    dimension(rs::AbelianGradedSpaceProd) -> Int

Get the total dimension of the direct product of several Abelian graded spaces.
"""
@inline dimension(rs::AbelianGradedSpaceProd) = prod(dimension, rs.contents)

"""
    dimension(rs::AbelianGradedSpaceProd, i::CartesianIndex) -> Int

Get the degenerate dimension of the ith Abelian quantum number in the direct product of several Abelian graded spaces.
"""
@inline dimension(rs::AbelianGradedSpaceProd, i::CartesianIndex) = prod(map(dimension, rs.contents, i.I))

"""
    dimension(rs::AbelianGradedSpaceProd{N, QN}, qns::NTuple{N, QN}) where {N, QN<:AbelianQuantumNumber} -> Int

Get the degenerate dimension of the Abelian quantum number fused by `qns` in the direct product of several Abelian graded spaces.
"""
@inline dimension(rs::AbelianGradedSpaceProd{N, QN}, qns::NTuple{N, QN}) where {N, QN<:AbelianQuantumNumber} = prod(map(dimension, rs.contents, qns))

"""
    range(rs::AbelianGradedSpaceProd, i::CartesianIndex) -> AbstractVector{Int}

Get the slice of the degenerate dimension of the ith Abelian quantum number in the direct product of several Abelian graded spaces.
"""
@inline function Base.range(rs::AbelianGradedSpaceProd, i::CartesianIndex)
    contents = reverse(rs.contents)
    linear = LinearIndices(map(dimension, contents))
    cartesian = CartesianIndices(map(range, contents, reverse(i.I)))
    return DirectProductedAbelianGradedSpaceRange(linear, cartesian)
end
struct DirectProductedAbelianGradedSpaceRange{N} <: AbstractVector{Int}
    linear::LinearIndices{N, NTuple{N, Base.OneTo{Int}}}
    cartesian::CartesianIndices{N, NTuple{N, UnitRange{Int}}}
end
@inline Base.size(r::DirectProductedAbelianGradedSpaceRange) = (length(r.cartesian),)
@inline Base.getindex(r::DirectProductedAbelianGradedSpaceRange, i::Integer) = r.linear[r.cartesian[i]]

"""
    range(rs::AbelianGradedSpaceProd{N, QN}, qns::NTuple{N, QN}) where {N, QN<:AbelianQuantumNumber} -> AbstractVector{Int}

Get the slice of the degenerate dimension of the Abelian quantum number fused by `qns` in the direct product of several Abelian graded spaces.
"""
@inline function Base.range(rs::AbelianGradedSpaceProd{N, QN}, qns::NTuple{N, QN}) where {N, QN<:AbelianQuantumNumber}
    contents = reverse(rs.contents)
    linear = LinearIndices(map(dimension, contents))
    cartesian = CartesianIndices(map(range, contents, reverse(qns)))
    return DirectProductedAbelianGradedSpaceRange(linear, cartesian)
end

"""
    merge(rs::AbelianGradedSpaceProd) -> Tuple{AbelianGradedSpace{eltype(rs)}, Dict{eltype(rs), Vector{NTuple{rank(rs), eltype(rs)}}}}

Get the decomposition of the direct product of several Abelian graded spaces and its corresponding fusion processes.

For a set of Abelian graded spaces (gs₁, gs₂, ...), their direct product space can contain several equivalent irreducible representations because for different sets of Abelian quantum numbers (qn₁, qn₂, ...) where qnᵢ∈gsᵢ, the fusion, i.e., `⊗(qn₁, qn₂, ...)` may give the same result `qn`. This function returns the decomposition of the direct product of (gs₁, gs₂, ...) as well as all the fusion processes of each quantum number contained in the decomposition.
"""
function Base.merge(rs::AbelianGradedSpaceProd)
    fusion = Dict{eltype(rs), Vector{NTuple{rank(rs), eltype(rs)}}}()
    record = Dict{eltype(rs), Int}()
    for ps in product(map(content->pairs(content, dimension), reverse(rs.contents))...)
        qns = reverse(map(p->p.first, ps))
        dim = prod(map(p->p.second, ps))
        qn = ⊗(qns...)
        if haskey(fusion, qn)
            push!(fusion[qn], qns)
            record[qn] += dim
        else
            fusion[qn] = [qns]
            record[qn] = dim
        end
    end
    quantumnumbers, dimensions = collect(keys(record)), collect(values(record))
    perm = sortperm(quantumnumbers)
    return AbelianGradedSpace(permute!(quantumnumbers, perm), permute!(dimensions, perm), false; ordercheck=false, duplicatecheck=false, degeneracycheck=false), fusion
end

"""
    split(target::QN, rs::AbelianGradedSpaceProd{N, QN}; nmax::Real=20) where {N, QN<:AbelianQuantumNumber} -> Set{NTuple{N, QN}}

Find a set of splittings of the target Abelian quantum number with respect to the direct product of several Abelian graded spaces.
"""
function Base.split(target::QN, rs::AbelianGradedSpaceProd{N, QN}; nmax::Real=20) where {N, QN<:AbelianQuantumNumber}
    result = Set{NTuple{N, QN}}()
    if isinf(nmax)
        for qns in product(rs.contents...)
            ⊗(qns...)==target && push!(result, qns)
        end
    else
        seed!()
        diff(indexes::NTuple{N, Int}) = norm(values(⊗(map(getindex, rs.contents, indexes)...)-target))
        count = 1
        while true
            indexes = map(content->rand(1:length(content)), rs.contents)
            Δ = diff(indexes)
            @inbounds while Δ > 0
                pos = rand(1:rank(rs))
                index = rand(1:length(rs.contents[pos])-1)
                newindexes = ntuple(i->i==pos ? (index<indexes[pos] ? index : index+1) : indexes[i], Val(N))
                newΔ = diff(newindexes)
                if newΔ<=Δ || exp(Δ-newΔ)>rand()
                    indexes = newindexes
                    Δ = newΔ
                end
            end
            push!(result, map(getindex, rs.contents, indexes))
            count = count + 1
            (length(result)>=nmax || count>nmax*5) && break
        end
    end
    return result
end

end #module
