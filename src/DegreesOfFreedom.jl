module DegreesOfFreedom

using DataStructures: OrderedDict
using Printf: @printf, @sprintf
using SparseArrays: SparseMatrixCSC, nnz
using StaticArrays: SVector
using ..QuantumLattices: add!, decompose, dimension, dtype
using ..QuantumOperators: ID, LinearTransformation, Operator, OperatorPack, Operators, OperatorUnit, QuantumOperator, valuetolatextext
using ..Spatials: Bond, Point
using ..Toolkit: atol, efficientoperations, rtol, CompositeDict, Float, VectorSpace, VectorSpaceCartesian, VectorSpaceDirectProducted, VectorSpaceDirectSummed, VectorSpaceStyle, concatenate, fulltype, parametertype, rawtype, reparameter, tostr

import LaTeXStrings: latexstring
import ..QuantumLattices: ⊕, ⊗, expand, expand!, id, ishermitian, kind, permute, rank, reset!, update!, value
import ..QuantumOperators: idtype, optype, script
import ..Spatials: icoordinate, rcoordinate
import ..Toolkit: contentnames, getcontent, parameternames, shape

export AbstractIndex, AllEqual, Boundary, CompositeInternal, ConstrainedInternal, Internal, InternalIndex, InternalIndexProd, InternalPattern, InternalProd, InternalSum, SimpleInternal, SimpleInternalIndex
export Component, CompositeIndex, CoordinatedIndex, Coupling, Hilbert, Index, MatrixCoupling, MatrixCouplingProd, MatrixCouplingSum, Metric, OperatorUnitToTuple, Ordinal, Pattern, Table, Term, TermAmplitude, TermCoupling, TermFunction
export ˢᵗ, ⁿᵈ, ʳᵈ, ᵗʰ, plain, allequalfields, indextype, isdefinite, partition, patternrule, sitestructure, statistics, @pattern

# AbstractIndex
"""
    AbstractIndex <: OperatorUnit

Abstract type of the index of degrees of freedom.
"""
abstract type AbstractIndex <: OperatorUnit end
@inline Base.getindex(::Type{AbstractIndex}, ::Type{I}) where {I<:AbstractIndex} = I

# InternalIndex and Internal
"""
    InternalIndex <: AbstractIndex

Internal index of the internal degrees of freedom.
"""
abstract type InternalIndex <: AbstractIndex end

"""
    isdefinite(ii::InternalIndex) -> Bool
    isdefinite(::Type{<:InternalIndex}) -> Bool

Judge whether an internal index denotes a definite internal degree of freedom.
"""
@inline isdefinite(ii::InternalIndex) = isdefinite(typeof(ii))

"""
    SimpleInternalIndex <: InternalIndex

Simple internal index, i.e., a complete set of indexes to denote an internal degree of freedom.
"""
abstract type SimpleInternalIndex <: InternalIndex end
@inline isdefinite(::Type{<:SimpleInternalIndex}) = false
@inline Base.show(io::IO, ii::SimpleInternalIndex) = @printf io "%s(%s)" AbstractIndex[typeof(ii)] join(map(tostr, ntuple(i->getfield(ii, i), Val(fieldcount(typeof(ii))))), ", ")

"""
    statistics(ii::SimpleInternalIndex) -> Symbol

Get the statistics of a simple internal index.
"""
@inline statistics(ii::SimpleInternalIndex) = statistics(typeof(ii))

"""
    InternalIndexProd{T<:Tuple{Vararg{SimpleInternalIndex}}} <: InternalIndex

Direct product of several simple internal indexes.
"""
struct InternalIndexProd{T<:Tuple{Vararg{SimpleInternalIndex}}} <: InternalIndex
    contents::T
end
@inline Base.show(io::IO, cii::InternalIndexProd) = @printf io "%s" join((string(cii[i]) for i = 1:rank(cii)), " ⊗ ")
@inline Base.length(cii::InternalIndexProd) = rank(cii)
@inline Base.firstindex(::InternalIndexProd) = 1
@inline Base.lastindex(cii::InternalIndexProd) = length(cii)
@inline Base.getindex(cii::InternalIndexProd, i::Integer) = cii.contents[i]
@inline Base.getindex(cii::InternalIndexProd, slice::AbstractVector{<:Integer}) = InternalIndexProd(cii.contents[slice])
@inline Base.getproperty(cii::InternalIndexProd, name::Symbol) = ciigetproperty(cii, Val(name))
@inline ciigetproperty(cii::InternalIndexProd, ::Val{:contents}) = getfield(cii, :contents)
@inline ciigetproperty(cii::InternalIndexProd, ::Val{name}) where name = getproperty(getfield(cii, :contents), name)
@inline isdefinite(::Type{<:InternalIndexProd{T}}) where {T<:Tuple{Vararg{SimpleInternalIndex}}} = all(map(isdefinite, fieldtypes(T)))

"""
    InternalIndexProd(contents::SimpleInternalIndex...)

Construct the direct product of several simple internal indexes.
"""
@inline InternalIndexProd(contents::SimpleInternalIndex...) = InternalIndexProd(contents)

"""
    rank(cii::InternalIndexProd) -> Int
    rank(::Type{<:InternalIndexProd{T}}) where {T<:Tuple{Vararg{SimpleInternalIndex}}} -> Int

Get the number of simple internal indexes in a direct product.
"""
@inline rank(cii::InternalIndexProd) = rank(typeof(cii))
@inline rank(::Type{<:InternalIndexProd{T}}) where {T<:Tuple{Vararg{SimpleInternalIndex}}} = fieldcount(T)

"""
    ⊗(iis::SimpleInternalIndex...) -> InternalIndexProd
    ⊗(ii::SimpleInternalIndex, cii::InternalIndexProd) -> InternalIndexProd
    ⊗(cii::InternalIndexProd, ii::SimpleInternalIndex) -> InternalIndexProd
    ⊗(cii₁::InternalIndexProd, cii₂::InternalIndexProd) -> InternalIndexProd

Direct product between internal indexes.
"""
@inline ⊗(iis::SimpleInternalIndex...) = InternalIndexProd(iis...)
@inline ⊗(ii::SimpleInternalIndex, cii::InternalIndexProd) = InternalIndexProd(ii, cii.contents...)
@inline ⊗(cii::InternalIndexProd, ii::SimpleInternalIndex) = InternalIndexProd(cii.contents..., ii)
@inline ⊗(cii₁::InternalIndexProd, cii₂::InternalIndexProd) = InternalIndexProd(cii₁.contents..., cii₂.contents...)

"""
    Internal{I<:InternalIndex} <: VectorSpace{I}

Internal space at a single point.
"""
abstract type Internal{I<:InternalIndex} <: VectorSpace{I} end

"""
    SimpleInternal{I<:SimpleInternalIndex} <: Internal{I}

Simple internal space at a single point.
"""
abstract type SimpleInternal{I<:SimpleInternalIndex} <: Internal{I} end
@inline VectorSpaceStyle(::Type{<:SimpleInternal}) = VectorSpaceCartesian()
@inline Base.show(io::IO, i::SimpleInternal) = @printf io "%s(%s)" typeof(i) join(("$name=$(getfield(i, name))" for name in fieldnames(typeof(i))), ", ")

"""
    statistics(i::SimpleInternal) -> Symbol
    statistics(::Type{<:SimpleInternal{I}}) where {I<:SimpleInternalIndex} -> Symbol

Get the statistics of a simple internal space.
"""
@inline statistics(i::SimpleInternal) = statistics(typeof(i))
@inline statistics(::Type{<:SimpleInternal{I}}) where {I<:SimpleInternalIndex} = statistics(I)

"""
    match(ii::SimpleInternalIndex, i::SimpleInternal) -> Bool
    match(::Type{I}, i::SimpleInternal) where {I<:SimpleInternalIndex} -> Bool
    match(ii::SimpleInternalIndex, ::Type{SI}) where {SI<:SimpleInternal} -> Bool
    match(::Type{I}, ::Type{SI}) where {I<:SimpleInternalIndex, SI<:SimpleInternal} -> Bool

Judge whether a simple internal space or the type of a simple internal space matches a simple internal index or the type of a simple internal index.

Here, "match" means that the `eltype` of the simple internal space is consistent with the type of the simple internal index, which usually means that they share the same type name.
"""
@inline Base.match(ii::SimpleInternalIndex, i::SimpleInternal) = match(typeof(ii), typeof(i))
@inline Base.match(::Type{I}, i::SimpleInternal) where {I<:SimpleInternalIndex} = match(I, typeof(i))
@inline Base.match(ii::SimpleInternalIndex, ::Type{SI}) where {SI<:SimpleInternal} = match(typeof(ii), SI)
@inline @generated Base.match(::Type{I}, ::Type{SI}) where {I<:SimpleInternalIndex, SI<:SimpleInternal} = nameof(I)==nameof(eltype(SI))

"""
    filter(ii::SimpleInternalIndex, i::SimpleInternal) -> Union{Nothing, typeof(internal)}
    filter(::Type{I}, i::SimpleInternal) where {I<:SimpleInternalIndex} -> Union{Nothing, typeof(internal)}

Filter a simple internal space with respect to the input simple internal index or the type of a simple internal index.
"""
@inline Base.filter(ii::SimpleInternalIndex, i::SimpleInternal) = filter(typeof(ii), i)
@inline Base.filter(::Type{I}, i::SimpleInternal) where {I<:SimpleInternalIndex} = filterhelper(Val(match(I, typeof(i))), i)
@inline @generated filterhelper(::Val{B}, internal) where B = B ? :(internal) : nothing

"""
    filter(ii::SimpleInternalIndex, ::Type{SI}) where {SI<:SimpleInternal}
    filter(::Type{I}, ::Type{SI}) where {I<:SimpleInternalIndex, SI<:SimpleInternal}

Filter the type of a simple internal space with respect to the input simple internal index or the type of a simple internal index.
"""
@inline Base.filter(ii::SimpleInternalIndex, ::Type{SI}) where {SI<:SimpleInternal} = filter(typeof(ii), SI)
@inline Base.filter(::Type{I}, ::Type{SI}) where {I<:SimpleInternalIndex, SI<:SimpleInternal} = filterhelper(Val(match(I, SI)), SI)

"""
    CompositeInternal{T<:Tuple{Vararg{SimpleInternal}}, I<:InternalIndex} <: Internal{I}

Abstract type of the composition (i.e., direct sum or direct product) of several simple internal spaces.
"""
abstract type CompositeInternal{T<:Tuple{Vararg{SimpleInternal}}, I<:InternalIndex} <: Internal{I} end

"""
    rank(ci::CompositeInternal) -> Int
    rank(::Type{<:CompositeInternal{T}}) where {T<:Tuple{Vararg{SimpleInternal}}} -> Int

Get the number of simple internal spaces in a composite internal space.
"""
@inline rank(ci::CompositeInternal) = rank(typeof(ci))
@inline rank(::Type{<:CompositeInternal{T}}) where {T<:Tuple{Vararg{SimpleInternal}}} = fieldcount(T)

"""
    filter(ii::SimpleInternalIndex, ci::CompositeInternal) -> Union{Nothing, SimpleInternal, CompositeInternal}
    filter(::Type{I}, ci::CompositeInternal) where {I<:SimpleInternalIndex} -> Union{Nothing, SimpleInternal, CompositeInternal}

Filter the composite internal space and select those that match the input simple internal index or the type of a simple internal index.
"""
@inline Base.filter(ii::SimpleInternalIndex, ci::CompositeInternal) = filter(typeof(ii), ci)
@inline Base.filter(::Type{I}, ci::CompositeInternal{T}) where {I<:SimpleInternalIndex, T<:Tuple{Vararg{SimpleInternal}}} = filterhelper(Val(map(SI->match(I, SI), fieldtypes(T))), I, ci)
@inline @generated function filterhelper(::Val{BS}, ::Type{I}, ci::CompositeInternal) where {BS, I<:SimpleInternalIndex}
    exprs = []
    for (i, B) in enumerate(BS)
        B && push!(exprs, :(filter(I, ci.contents[$i])))
    end
    length(exprs)==0 && return
    length(exprs)==1 && return first(exprs)
    return :(rawtype(typeof(ci))($(exprs...)))
end

"""
    filter(ii::SimpleInternalIndex, ::Type{CI}) where {CI<:CompositeInternal}
    filter(::Type{I}, ::Type{CI}) where {I<:SimpleInternalIndex, CI<:CompositeInternal}

Filter the type of a composite internal space and select those that match the input simple internal index or the type of a simple internal index.
"""
@inline Base.filter(ii::SimpleInternalIndex, ::Type{CI}) where {CI<:CompositeInternal} = filter(typeof(ii), CI)
@inline Base.filter(::Type{I}, ::Type{CI}) where {I<:SimpleInternalIndex, T<:Tuple{Vararg{SimpleInternal}}, CI<:CompositeInternal{T}} = filterhelper(Val(map(SI->match(I, SI), fieldtypes(T))), I, CI)
@inline @generated function filterhelper(::Val{BS}, ::Type{I}, ::Type{CI}) where {BS, I<:SimpleInternalIndex, T<:Tuple{Vararg{SimpleInternal}}, CI<:CompositeInternal{T}}
    exprs = []
    for (i, B) in enumerate(BS)
        B && push!(exprs, :(filter(I, fieldtype(T, $i))))
    end
    length(exprs)==0 && return
    length(exprs)==1 && return first(exprs)
    return Expr(:curly, nameof(CI), Expr(:curly, :Tuple, exprs...), :(eltype(rawtype(CI), $(exprs...))))
end

"""
    InternalSum{T<:Tuple{Vararg{SimpleInternal}}, I<:InternalIndex} <: CompositeInternal{T, I}

Direct sum of several single internal spaces.
"""
struct InternalSum{T<:Tuple{Vararg{SimpleInternal}}, I<:InternalIndex} <: CompositeInternal{T, I}
    contents::T
    InternalSum(contents::Tuple{Vararg{SimpleInternal}}) = new{typeof(contents), eltype(InternalSum, fieldtypes(typeof(contents))...)}(contents)
end
@inline VectorSpaceStyle(::Type{<:InternalSum}) = VectorSpaceDirectSummed()
@inline Base.show(io::IO, ci::InternalSum) = @printf io "%s" join((string(ci.contents[i]) for i = 1:rank(ci)), " ⊕ ")
@inline Base.eltype(::Type{InternalSum}, types...) = mapreduce(eltype, typejoin, types; init=Union{})
@inline InternalSum(contents::SimpleInternal...) = InternalSum(contents)

"""
    ⊕(is::SimpleInternal...) -> InternalSum
    ⊕(i::SimpleInternal, ci::InternalSum) -> InternalSum
    ⊕(ci::InternalSum, i::SimpleInternal) -> InternalSum
    ⊕(ci₁::InternalSum, ci₂::InternalSum) -> InternalSum

Direct sum between simple internal spaces and composite internal spaces.
"""
@inline ⊕(is::SimpleInternal...) = InternalSum(is...)
@inline ⊕(i::SimpleInternal, ci::InternalSum) = InternalSum(i, ci.contents...)
@inline ⊕(ci::InternalSum, i::SimpleInternal) = InternalSum(ci.contents..., i)
@inline ⊕(ci₁::InternalSum, ci₂::InternalSum) = InternalSum(ci₁.contents..., ci₂.contents...)

"""
    InternalProd{T<:Tuple{Vararg{SimpleInternal}}, I<:InternalIndex} <: CompositeInternal{T, I}

Direct product of several single internal spaces.
"""
struct InternalProd{T<:Tuple{Vararg{SimpleInternal}}, I<:InternalIndex} <: CompositeInternal{T, I}
    contents::T
    InternalProd(contents::Tuple{Vararg{SimpleInternal}}) = new{typeof(contents), eltype(InternalProd, fieldtypes(typeof(contents))...)}(contents)
end
@inline VectorSpaceStyle(::Type{<:InternalProd}) = VectorSpaceDirectProducted(:forward)
@inline Base.convert(::Type{<:InternalIndexProd{T}}, ii::T, ::InternalProd) where {T<:Tuple{Vararg{SimpleInternalIndex}}} = InternalIndexProd(ii)
@inline Base.show(io::IO, ci::InternalProd) = @printf io "%s" join((string(ci.contents[i]) for i = 1:rank(ci)), " ⊗ ")
@inline Base.eltype(::Type{InternalProd}, types...) = InternalIndexProd{Tuple{map(eltype, types)...}}
@inline InternalProd(contents::SimpleInternal...) = InternalProd(contents)

"""
    ⊗(is::SimpleInternal...) -> InternalProd
    ⊗(i::SimpleInternal, ci::InternalProd) -> InternalProd
    ⊗(ci::InternalProd, i::SimpleInternal) -> InternalProd
    ⊗(ci₁::InternalProd, ci₂::InternalProd) -> InternalProd

Direct product between simple internal spaces and composite internal spaces.
"""
@inline ⊗(is::SimpleInternal...) = InternalProd(is...)
@inline ⊗(i::SimpleInternal, ci::InternalProd) = InternalProd(i, ci.contents...)
@inline ⊗(ci::InternalProd, i::SimpleInternal) = InternalProd(ci.contents..., i)
@inline ⊗(ci₁::InternalProd, ci₂::InternalProd) = InternalProd(ci₁.contents..., ci₂.contents...)

# InternalPattern
"""
    AllEqual{Fields} <: Function

All-equal constraint for the direct product of homogenous simple internal indexes that for every specified field the values of the internal indexes should be all equal.
"""
struct AllEqual{Fields} <: Function
    AllEqual(fields::Tuple{Vararg{Symbol}}) = new{fields}()
end
@inline Base.show(io::IO, ::AllEqual{Fields}) where Fields = @printf io "AllEqual(%s)" join(map(repr, Fields), ", ")

"""
    AllEqual(fields::Tuple{Vararg{Symbol}})
    AllEqual(fields::Symbol...)

Construct an all-equal constraint.
"""
@inline AllEqual(fields::Symbol...) = AllEqual(fields)

"""
    AllEqual(::Type{I}) where {I<:SimpleInternalIndex}

Construct an all-equal constraint based on the type of a simple internal index.
"""
@inline AllEqual(::Type{I}) where {I<:SimpleInternalIndex} = constraint(I, I|>allequalfields|>Val)
@generated function constraint(::Type{I}, ::Val{fields}) where {I<:SimpleInternalIndex, fields}
    result = Symbol[]
    for field in fields
        F = fieldtype(I, field)
        F==Symbol && error("AllEqual error: $I is too complicated for an all-equal constraint to construct since it host the `:$field` field represented by a `Symbol` other than a `Colon`.")
        F==Colon && push!(result, field)
    end
    return AllEqual(result...)
end

"""
    allequalfields(::Type{I}) where {I<:SimpleInternalIndex} -> Tuple{Vararg{Symbol}}

Get the field names that can be subject to all-equal constraint based on the type of a simple internal index.
"""
@inline allequalfields(::Type{I}) where {I<:SimpleInternalIndex} = fieldnames(I)

"""
    (constraint::AllEqual{fields})(index::InternalIndexProd{NTuple{N, I}}) where {fields, N, I<:SimpleInternalIndex} -> Bool

Judge whether the direct product of a set of homogenous simple internal indexes is subject to an all-equal constraint.
"""
@generated function (constraint::AllEqual{fields})(index::InternalIndexProd{NTuple{N, I}}) where {fields, N, I<:SimpleInternalIndex}
    N==0 && return true
    exprs = []
    for field in QuoteNode.(fields)
        push!(exprs, Expr(:call, allequal, Expr(:tuple, [:(getfield(index[$i], $field)) for i=1:N]...)))
    end
    length(exprs)==0 && return true
    return Expr(:call, all, Expr(:tuple, exprs...))
end

"""
    InternalPattern{P, I<:Tuple{Vararg{SimpleInternalIndex}}, N, C<:NTuple{N, Function}} <: QuantumOperator

Internal pattern for the direct product of a set of simple internal indexes with extra constraints.
"""
struct InternalPattern{P, I<:Tuple{Vararg{SimpleInternalIndex}}, N, C<:NTuple{N, Function}} <: QuantumOperator
    index::InternalIndexProd{I}
    constraints::C
    representations::NTuple{N, String}
    function InternalPattern{P}(index::InternalIndexProd{I}, constraints::NTuple{N, Function}, representations::NTuple{N, String}=map(string, constraints))  where {P, I<:Tuple{Vararg{SimpleInternalIndex}}, N}
        @assert isa(P, NTuple{N, Int}) "InternalPattern error: partition ($P) is not $N-tuple of integers."
        @assert sum(P)==rank(index) "InternalPattern error: sum of partition ($(sum(P))) is not equal to the rank ($(rank(index))) of the direct product of simple internal indexes."
        new{P, I, N, typeof(constraints)}(index, constraints, representations)
    end
end
@inline Base.:(==)(pattern₁::InternalPattern{P₁}, pattern₂::InternalPattern{P₂}) where {P₁, P₂} = P₁==P₂ && pattern₁.index==pattern₂.index && pattern₁.representations==pattern₂.representations
@inline Base.:isequal(pattern₁::InternalPattern{P₁}, pattern₂::InternalPattern{P₂}) where {P₁, P₂} = isequal(P₁, P₂) && isequal(pattern₁.index, pattern₂.index) && isequal(pattern₁.representations, pattern₂.representations)
@inline Base.hash(pattern::InternalPattern{P}, h::UInt) where P = hash((P..., pattern.index.contents..., pattern.representations...), h)
function Base.show(io::IO, pattern::InternalPattern)
    len, count = length(pattern.representations), 1
    for i = 1:len
        r = rank(pattern, i)
        start, stop = count, count+r-1
        index = pattern.index[start:stop]
        if i==1
            @printf io "%s" (isdefinite(index) ? "" : "∑")
            (len>1 || !isdefinite(index)) && @printf io "%s" "["
        else
            @printf io " ⊗ %s[" (isdefinite(index) ? "" : "∑")
        end
        @printf io "%s" join(index.contents, " ")
        (len>1 || !isdefinite(index)) && @printf io "%s" "]"
        representation = pattern.representations[i]
        length(representation)==0 || occursin("AllEqual", representation) || @printf io "(%s)" representation
        count = stop+1
    end
end

"""
    InternalPattern(index::InternalIndexProd, constraint::Function, representation::String=string(constraint))
    InternalPattern{P}(index::InternalIndexProd, constraints::NTuple{N, Function}, representations::NTuple{N, String}=map(string, constraints))  where {P, N}

Construct an internal pattern for the direct product of a set of simple internal indexes with 1) only one constraint function, and 2) several constraint functions.
"""
@inline InternalPattern(index::InternalIndexProd, constraint::Function, representation::String=string(constraint)) = InternalPattern{(rank(index),)}(index, (constraint,), (representation,))

"""
    InternalPattern(index::SimpleInternalIndex)
    InternalPattern(index::InternalIndexProd{<:Tuple{Vararg{I}}}) where {I<:SimpleInternalIndex}

Construct an internal pattern whose constraint is an [`AllEqual`](@ref) function for the direct product of a homogeneous set of simple internal indexes.
"""
@inline InternalPattern(index::SimpleInternalIndex) = InternalPattern(InternalIndexProd(index))
@inline InternalPattern(index::InternalIndexProd{<:Tuple{Vararg{I}}}) where {I<:SimpleInternalIndex} = InternalPattern(index, AllEqual(I))

"""
    partition(pattern::InternalPattern) -> Tuple{Vararg{Int}}
    partition(::Type{<:InternalPattern{P}}) where P -> P

Get the partition of the direct product of a set of simple internal indexes on which the extra constraints operate independently.
"""
@inline partition(pattern::InternalPattern) = partition(typeof(pattern))
@inline partition(::Type{<:InternalPattern{P}}) where P = P

"""
    rank(pattern::InternalPattern) -> Int
    rank(::Type{<:InternalPattern{P}}) where P -> Int

Get the rank of the direct product of the simple internal indexes that an internal pattern can apply.
"""
@inline rank(pattern::InternalPattern) = rank(typeof(pattern))
@inline @generated rank(::Type{<:InternalPattern{P}}) where P = sum(P)

"""
    rank(pattern::InternalPattern, i::Integer) -> Int
    rank(::Type{<:InternalPattern{P}}, i::Integer) where P -> Int

Get the rank of the direct product of the simple internal indexes that the ith constraint in an internal pattern can apply.
"""
@inline rank(pattern::InternalPattern, i::Integer) = rank(typeof(pattern), i)
@inline rank(::Type{<:InternalPattern{P}}, i::Integer) where P = P[i]

"""
    match(pattern::InternalPattern, index::InternalIndexProd) -> Bool

Judge whether the direct product of a set of simple internal indexes satisfies an internal pattern.
"""
@generated function Base.match(pattern::InternalPattern{P}, index::InternalIndexProd) where P
    @assert rank(pattern)==rank(index) "match error: mismatched ranks of the direct product of simple internal indexes and internal pattern."
    exprs, count = [], 1
    for (i, r) in enumerate(P)
        start, stop = count, count+r-1
        segment = Expr(:call, :InternalIndexProd, Expr(:tuple, [:(index[$pos]) for pos=start:stop]...))
        push!(exprs, :(pattern.constraints[$i]($segment)::Bool))
        count = stop+1
    end
    return Expr(:call, :all, Expr(:tuple, exprs...))
end

"""
    ⊗(pattern₁::InternalPattern, pattern₂::InternalPattern) -> InternalPattern

Get the combination of two internal patterns.
"""
@inline function ⊗(pattern₁::InternalPattern{P₁}, pattern₂::InternalPattern{P₂}) where {P₁, P₂}
    return InternalPattern{(P₁..., P₂...)}(pattern₁.index⊗pattern₂.index, (pattern₁.constraints..., pattern₂.constraints...), (pattern₁.representations..., pattern₂.representations...))
end

"""
    latexstring(pattern::InternalPattern) -> String

Convert an internal pattern to the latex format.
"""
function latexstring(pattern::InternalPattern)
    result = String[]
    len, count = length(pattern.representations), 1
    for i = 1:len
        r = rank(pattern, i)
        start, stop = count, count+r-1
        index = pattern.index[start:stop]
        representation = pattern.representations[i]
        summation = occursin("AllEqual", representation) ? "" : replace(replace(representation, "&&"=>"\\,\\text{and}\\,"), "||"=>"\\,\\text{or}\\,")
        summation=="" || (summation = join(push!(symbols(index.contents, representation), summation), ",\\,"))
        context = isdefinite(index) ? "" : "\\sum_{$summation}" 
        for content in index.contents
            context = @sprintf "%s%s%s" context (length(context)>0 ? " " : "") latexstring(content)
        end
        push!(result, context)
        count = stop+1
    end
    return join(result, " \\cdot ")
end
function symbols(indexes::Tuple{Vararg{SimpleInternalIndex}}, constraint::String)
    result = String[]
    for index in indexes
        for i = 1:fieldcount(typeof(index))
            attr = getfield(index, i)
            if isa(attr, Symbol)
                value = string(attr)
                occursin(value, constraint) || push!(result, value)
            end
        end
    end
    return sort!(unique!(result))
end

# ConstrainedInternal
"""
    ConstrainedInternal{P<:InternalProd, C<:InternalPattern}

Constrained internal space.
"""
struct ConstrainedInternal{P<:InternalProd, C<:InternalPattern}
    internal::P
    pattern::C
    function ConstrainedInternal(internal::InternalProd, pattern::InternalPattern)
        @assert rank(internal)==rank(pattern) "ConstrainedInternal error: mismatched ranks of the input internal space and internal pattern."
        @assert all(map(match, pattern.index.contents, internal.contents)) "ConstrainedInternal error: mismatched simple internal spaces with respect to simple internal indexes."
        new{typeof(internal), typeof(pattern)}(internal, pattern)
    end
end
@inline Base.eltype(ci::ConstrainedInternal) = eltype(typeof(ci))
@inline Base.eltype(::Type{<:ConstrainedInternal{P}}) where {P<:InternalProd} = eltype(P)
@inline Base.IteratorSize(::Type{<:ConstrainedInternal}) = Base.SizeUnknown()
function Base.iterate(ci::ConstrainedInternal, state=(space=InternalIndexSpace(ci.internal, ci.pattern.index); (space, iterate(space))))
    space, state = state
    while !isnothing(state)
        index, state = state
        state = iterate(space, state)
        match(ci.pattern, index) && return index, (space, state)
    end
    return
end
struct InternalIndexSpace{I<:InternalIndex, E<:InternalIndex, V<:Internal{E}} <: VectorSpace{E}
    internal::V
    index::I
end
@inline VectorSpaceStyle(::Type{<:InternalIndexSpace{<:InternalIndexProd, <:InternalIndexProd, <:CompositeInternal}}) = VectorSpaceDirectProducted(:forward)
@inline function getcontent(iispace::InternalIndexSpace{<:InternalIndexProd, <:InternalIndexProd, <:CompositeInternal}, ::Val{:contents})
    return map((internal, index)->InternalIndexSpace(internal, index), iispace.internal.contents, iispace.index.contents)
end
@inline function Base.convert(::Type{InternalIndexProd{T}}, indexes::T, ::InternalIndexSpace{<:InternalIndexProd, InternalIndexProd{T}, <:CompositeInternal}) where {T<:Tuple{Vararg{SimpleInternalIndex}}}
    return InternalIndexProd(indexes)
end
@inline VectorSpaceStyle(::Type{<:InternalIndexSpace{<:SimpleInternalIndex, <:SimpleInternalIndex, <:SimpleInternal}}) = VectorSpaceCartesian()
@inline shape(iispace::InternalIndexSpace{<:SimpleInternalIndex, <:SimpleInternalIndex, <:SimpleInternal}) = shape(iispace.internal, iispace.index)
@inline function Base.convert(::Type{<:CartesianIndex}, index::E, iispace::InternalIndexSpace{<:SimpleInternalIndex, E, <:SimpleInternal}) where {E<:SimpleInternalIndex}
    return convert(CartesianIndex, index, iispace.internal)
end
@inline function Base.convert(::Type{E}, index::CartesianIndex, iispace::InternalIndexSpace{<:SimpleInternalIndex, E, <:SimpleInternal}) where {E<:SimpleInternalIndex}
    return convert(E, index, iispace.internal)
end

"""
    ConstrainedInternal(internal::SimpleInternal, pattern::InternalPattern)
    ConstrainedInternal(internal::InternalProd, pattern::InternalPattern)

Construct a constrained internal space.
"""
@inline ConstrainedInternal(internal::SimpleInternal, pattern::InternalPattern) = ConstrainedInternal(InternalProd(internal), pattern)

"""
    shape(internal::SimpleInternal, index::SimpleInternalIndex) -> OrdinalRange{Int, Int}

The shape of a simple internal space when a labeled simple internal index are considered.

A constrained internal space need this function to generate all the internal indexes that match the internal pattern, which gets a default implementation, i.e.,
```julia
shape(internal::SimpleInternal, index::SimpleInternalIndex) = shape(internal)
```
It can be overloaded to restrict the shape of a simple internal space based on the input simple internal index to significantly improve efficiency, but this is not necessary.
"""
@inline shape(internal::SimpleInternal, ::SimpleInternalIndex) = shape(internal)

# Index
"""
    Ordinal

Ordinal of an `Int`.
"""
struct Ordinal
    n::Int
end
@inline Base.show(io::IO, nth::Ordinal) = @printf io "%s" nth.n==1 ? "1ˢᵗ" : nth.n==2 ? "2ⁿᵈ" : nth.n==3 ? "3ʳᵈ" : "$(nth.n)ᵗʰ"
@inline Base.:*(i::Int, nth::Ordinal) = Ordinal(i*nth.n)
@inline Base.getindex(v, i::Ordinal) = v[i.n]

"""
    const ˢᵗ = ⁿᵈ = ʳᵈ = ᵗʰ = Ordinal(1)

Constant ordinals.
"""
const ˢᵗ = ⁿᵈ = ʳᵈ = ᵗʰ = Ordinal(1)

"""
    Index(site::Union{Int, Ordinal, Colon}, internal::SimpleInternalIndex)

Index of a degree of freedom, which consist of the spatial part (i.e., the site index) and the internal part (i.e., the internal index).
"""
struct Index{I<:SimpleInternalIndex, S<:Union{Int, Ordinal, Colon}} <: AbstractIndex
    site::S
    internal::I
end
@inline parameternames(::Type{<:Index}) = (:internal, :site)
function Base.show(io::IO, index::Index)
    internal = String[]
    for i = 1:fieldcount(indextype(index))
        value = getfield(index.internal, i)
        push!(internal, tostr(value))
    end
    internal = join(internal, ", ")
    @printf io "%s(%s%s%s)" AbstractIndex[typeof(index)] tostr(index.site) (length(internal)>0 ? ", " : "") internal
end

"""
    isdefinite(index::Index) -> Bool
    isdefinite(::Type{<:Index{I}}) where {I<:SimpleInternalIndex} -> Bool

Determine whether an index denotes a definite degree of freedom.
"""
@inline isdefinite(index::Index) = isdefinite(typeof(index))
@inline isdefinite(::Type{<:Index{I}}) where {I<:SimpleInternalIndex} = isdefinite(I)

"""
    isdefinite(indexes::Tuple{Vararg{Index}}) -> Bool
    isdefinite(::Type{T}) where {T<:Tuple{Vararg{Index}}} -> Bool

Determine whether a tuple of indexes denotes a definite degree of freedom.
"""
@inline isdefinite(indexes::Tuple{Vararg{Index}}) = isdefinite(typeof(indexes))
@inline isdefinite(::Type{T}) where {T<:Tuple{Vararg{Index}}} = all(map(isdefinite, fieldtypes(T)))

"""
    indextype(index::Index)
    indextype(::Type{I}) where {I<:Index}

Get the type of the internal part of an index.
"""
@inline indextype(index::Index) = indextype(typeof(index))
@inline indextype(::Type{I}) where {I<:Index} = parametertype(I, :internal)

"""
    statistics(index::Index) -> Symbol
    statistics(::Type{<:Index{I}}) where {I<:SimpleInternalIndex} -> Symbol

Get the statistics of an index.
"""
@inline statistics(index::Index) = statistics(typeof(index))
@inline statistics(::Type{<:Index{I}}) where {I<:SimpleInternalIndex} = statistics(I)

"""
    adjoint(index::Index) -> typeof(index)

Get the adjoint of an index.
"""
@inline Base.adjoint(index::Index) = rawtype(typeof(index))(index.site, adjoint(index.internal))

"""
    script(index::Index, ::Val{:site}; kwargs...) -> String
    script(index::Index, attr::Val; kwargs...) -> String

Get the required script of an index.
"""
@inline script(index::Index, ::Val{:site}; kwargs...) = tostr(index.site)
@inline script(index::Index, attr::Val; kwargs...) = script(index.internal, attr; kwargs...)

"""
    permute(index₁::Index, index₂::Index) -> Tuple{Vararg{Operator}}

Get the permutation of two indexes.
"""
function permute(index₁::Index, index₂::Index)
    if index₁.site ≠ index₂.site
        return (Operator(statistics(index₁)==:f && statistics(index₂)==:f ? -1 : 1, index₂, index₁),)
    else
        return map(op->Operator(value(op), map(internal->Index(index₁.site, internal), id(op))), permute(index₁.internal, index₂.internal))
    end
end

"""
    indextype(I::Type{<:SimpleInternal})

Get the compatible type of the index based on the type of an internal space.
"""
@inline indextype(I::Type{<:SimpleInternal}) = fulltype(Index, NamedTuple{(:internal, :site), Tuple{eltype(I), Int}})

# CompositeIndex and CoordinatedIndex
"""
    CompositeIndex{I<:Index} <: AbstractIndex

Abstract type of a composite index.
"""
abstract type CompositeIndex{I<:Index} <: AbstractIndex end
@inline contentnames(::Type{<:CompositeIndex}) = (:index,)
@inline parameternames(::Type{<:CompositeIndex}) = (:index,)

"""
    indextype(::CompositeIndex)
    indextype(::Type{<:CompositeIndex})

Get the index type of a composite index.
"""
@inline indextype(index::CompositeIndex) = indextype(typeof(index))
@inline @generated indextype(::Type{I}) where {I<:CompositeIndex} = parametertype(supertype(I, :CompositeIndex), :index)

"""
    statistics(index::CompositeIndex) -> Symbol
    statistics(::Type{<:CompositeIndex{I}}) where {I<:Index} -> Symbol

Get the statistics of a composite index.
"""
@inline statistics(index::CompositeIndex) = statistics(typeof(index))
@inline statistics(::Type{<:CompositeIndex{I}}) where {I<:Index} = statistics(I)

"""
    script(index::CompositeIndex, ::Val{attr}; kwargs...) where attr

Get the `attr` script of a composite index.
"""
@inline script(index::CompositeIndex, ::Val{attr}; kwargs...) where attr = script(getcontent(index, :index), Val(attr); kwargs...)

"""
    CoordinatedIndex{I<:Index, V<:SVector} <: CompositeIndex{I}

Coordinated index, i.e., index with coordinates in the unitcell.
"""
struct CoordinatedIndex{I<:Index, V<:SVector} <: CompositeIndex{I}
    index::I
    rcoordinate::V
    icoordinate::V
    CoordinatedIndex(index::Index, rcoordinate::V, icoordinate::V) where {V<:SVector} = new{typeof(index), V}(index, compositeindexcoordinate(rcoordinate), compositeindexcoordinate(icoordinate))
end
@inline contentnames(::Type{<:CoordinatedIndex}) = (:index, :rcoordinate, :icoordinate)
@inline parameternames(::Type{<:CoordinatedIndex}) = (:index, :coordination)
@inline Base.hash(index::CoordinatedIndex, h::UInt) = hash((index.index, Tuple(index.rcoordinate)), h)
@inline Base.propertynames(::ID{CoordinatedIndex}) = (:indexes, :rcoordinates, :icoordinates)
function Base.show(io::IO, index::CoordinatedIndex)
    internal = String[]
    for i = 1:fieldcount(typeof(index.index.internal))
        value = getfield(index.index.internal, i)
        push!(internal, tostr(value))
    end
    internal = join(internal, ", ")
    @printf io "%s(%s%s%s, %s, %s)" AbstractIndex[typeof(index)] tostr(index.index.site) (length(internal)>0 ? ", " : "") internal index.rcoordinate index.icoordinate
end
@inline compositeindexcoordinate(vector::SVector) = vector
@inline compositeindexcoordinate(vector::SVector{N, Float}) where N = SVector(ntuple(i->vector[i]===-0.0 ? 0.0 : vector[i], Val(N)))

"""
    CoordinatedIndex(index::Index, rcoordinate, icoordinate)
    CoordinatedIndex(index::Index; rcoordinate, icoordinate)

Construct a coordinated index.
"""
@inline CoordinatedIndex(index::Index, rcoordinate, icoordinate) = CoordinatedIndex(index, SVector{length(rcoordinate)}(rcoordinate), SVector{length(icoordinate)}(icoordinate))
@inline CoordinatedIndex(index::Index; rcoordinate, icoordinate) = CoordinatedIndex(index, rcoordinate, icoordinate)

"""
    adjoint(index::CoordinatedIndex) -> typeof(index)

Get the adjoint of a coordinated index.
"""
@inline Base.adjoint(index::CoordinatedIndex) = CoordinatedIndex(index.index', index.rcoordinate, index.icoordinate)

"""
    script(index::CoordinatedIndex, ::Val{:rcoordinate}; kwargs...) -> String
    script(index::CoordinatedIndex, ::Val{:icoordinate}; kwargs...) -> String

Get the rcoordinate/icoordinate script of a coordinated index.
"""
@inline script(index::CoordinatedIndex, ::Val{:rcoordinate}; kwargs...) = @sprintf "[%s]" join(valuetolatextext.(index.rcoordinate), ", ")
@inline script(index::CoordinatedIndex, ::Val{:icoordinate}; kwargs...) = @sprintf "[%s]" join(valuetolatextext.(index.icoordinate), ", ")

"""
    script(index::CoordinatedIndex, ::Val{:integercoordinate}; vectors, kwargs...)

Get the integer coordinate script of a coordinated index.
"""
function script(index::CoordinatedIndex, ::Val{:integercoordinate}; vectors, kwargs...)
    rcoeff = decompose(index.icoordinate, vectors...)
    icoeff = Int.(round.(rcoeff))
    @assert isapprox(efficientoperations, rcoeff, icoeff) "script error: mismatched icoordinate of the input coordinated index and vectors."
    return @sprintf "[%s]" join(icoeff, ", ")
end

"""
    permute(index₁::CoordinatedIndex, index₂::CoordinatedIndex) -> Tuple{Vararg{Operator}}

Get the permutation of two coordinated indexes.
"""
function permute(index₁::CoordinatedIndex, index₂::CoordinatedIndex)
    if index₁.rcoordinate ≠ index₂.rcoordinate || index₁.icoordinate ≠ index₂.icoordinate
        return (Operator(statistics(index₁)==:f && statistics(index₂)==:f ? -1 : 1, index₂, index₁),)
    else
        rcoordinate, icoordinate = index₁.rcoordinate, index₁.icoordinate
        return map(op->Operator(value(op), map(index->CoordinatedIndex(index, rcoordinate, icoordinate), id(op))), permute(index₁.index, index₂.index))
    end
end

"""
    indextype(I::Type{<:SimpleInternal}, P::Type{<:Point})

Get the compatible type of the coordinated index based on the type of an internal space and the type of a point.
"""
@inline indextype(I::Type{<:SimpleInternal}, P::Type{<:Point}) = fulltype(CoordinatedIndex, NamedTuple{(:index, :coordination), Tuple{indextype(I), SVector{dimension(P), dtype(P)}}})

"""
    rcoordinate(opt::Operator{<:Number, <:ID{CoordinatedIndex}}) -> SVector

Get the whole rcoordinate of an operator.
"""
@inline function rcoordinate(opt::Operator{<:Number, <:ID{CoordinatedIndex}})
    rank(opt)==1 && return id(opt)[1].rcoordinate
    rank(opt)==2 && return id(opt)[2].rcoordinate-id(opt)[1].rcoordinate
    error("rcoordinate error: not supported rank($(rank(opt))) of $(nameof(opt)).")
end

"""
    icoordinate(opt::Operator{<:Number, <:ID{CoordinatedIndex}}) -> SVector

Get the whole icoordinate of an operator.
"""
@inline function icoordinate(opt::Operator{<:Number, <:ID{CoordinatedIndex}})
    rank(opt)==1 && return id(opt)[1].icoordinate
    rank(opt)==2 && return id(opt)[2].icoordinate-id(opt)[1].icoordinate
    error("icoordinate error: not supported rank($(rank(opt))) of $(nameof(opt)).")
end

# Hilbert
"""
    Hilbert{I<:Internal} <: CompositeDict{Int, I}

Hilbert space at a lattice.
"""
struct Hilbert{I<:Internal} <: CompositeDict{Int, I}
    contents::OrderedDict{Int, I}
    Hilbert(contents::OrderedDict{Int, <:Internal}) = new{valtype(contents)}(contents)
end

"""
    Hilbert(ps::Pair...)
    Hilbert(kv)

Construct a Hilbert space the same way as an `OrderedDict`.
"""
@inline Hilbert(ps::Pair...) = Hilbert(ps)
@inline Hilbert(kv) = Hilbert(OrderedDict(kv))

"""
    Hilbert(internals::Internal...)
    Hilbert(internals::Tuple{Vararg{Internal}})
    Hilbert(internals::AbstractVector{<:Internal})

Construct a Hilbert space with the given internal spaces.
"""
@inline Hilbert(internals::Internal...) = Hilbert(internals)
@inline Hilbert(internals::Tuple{Vararg{Internal}}) = Hilbert(i=>internal for (i, internal) in enumerate(internals))
@inline Hilbert(internals::AbstractVector{<:Internal}) = Hilbert(i=>internal for (i, internal) in enumerate(internals))

"""
    Hilbert(internal::Internal, num::Integer)

Construct a Hilbert space with all same internal spaces.
"""
@inline Hilbert(internal::Internal, num::Integer) = Hilbert(i=>internal for i in 1:num)

# Pattern
"""
    Pattern{S<:Tuple{Vararg{Union{Ordinal, Colon}}}, C<:InternalPattern} <: QuantumOperator

Coupling pattern.
"""
struct Pattern{S<:Tuple{Vararg{Union{Ordinal, Colon}}}, C<:InternalPattern} <: QuantumOperator
    sites::S
    internal::C
end
@inline parameternames(::Type{<:Pattern}) = (:sites, :internal)
@inline Base.hash(pattern::Pattern, h::UInt) = hash((pattern.sites, pattern.internal), h)
function Base.show(io::IO, pattern::Pattern)
    len, count = length(pattern.internal.representations), 1
    for i = 1:len
        r = rank(pattern.internal, i)
        start, stop = count, count+r-1
        indexes = map(Index, pattern.sites[start:stop], pattern.internal.index[start:stop].contents)
        if i==1
            @printf io "%s" (isdefinite(indexes) ? "" : "∑")
            (len>1 || !isdefinite(indexes)) && @printf io "%s" "["
        else
            @printf io " ⊗ %s[" (isdefinite(indexes) ? "" : "∑")
        end
        @printf io "%s" join(indexes, " ")
        (len>1 || !isdefinite(indexes)) && @printf io "%s" "]"
        representation = pattern.internal.representations[i]
        length(representation)==0 || occursin("AllEqual", representation) || @printf io "(%s)" representation
        count = stop+1
    end
end

"""
    Pattern(indexes::Index...)
    Pattern(indexes::Tuple{Vararg{Index}})

Construct a coupling pattern from a set of indexes.
"""
@inline Pattern(indexes::Index...) = Pattern(indexes)
@inline Pattern(indexes::Tuple{Vararg{Index}}) = Pattern(indexes.sites, InternalPattern(InternalIndexProd(indexes.internals)))

"""
    rank(pattern::Pattern) -> Int
    rank(::Type{<:Pattern{S}}) where {S<:Tuple{Vararg{Union{Ordinal, Colon}}} -> Int

Get the rank of a coupling pattern.
"""
@inline rank(pattern::Pattern) = rank(typeof(pattern))
@inline rank(::Type{<:Pattern{S}}) where {S<:Tuple{Vararg{Union{Ordinal, Colon}}}} = fieldcount(S)

"""
    ⊗(pattern₁::Pattern, pattern₂::Pattern) -> Pattern

Get the combination of two coupling patterns.
"""
@inline ⊗(pattern₁::Pattern, pattern₂::Pattern) = Pattern((pattern₁.sites..., pattern₂.sites...), pattern₁.internal⊗pattern₂.internal)

"""
    latexstring(pattern::Pattern) -> String

Convert a coupling pattern to the latex format.
"""
function latexstring(pattern::Pattern)
    result = String[]
    len, count = length(pattern.internal.representations), 1
    for i = 1:len
        r = rank(pattern.internal, i)
        start, stop = count, count+r-1
        indexes = map(Index, pattern.sites[start:stop], pattern.internal.index[start:stop].contents)
        representation = pattern.internal.representations[i]
        summation = occursin("AllEqual", representation) ? "" : replace(replace(representation, "&&"=>"\\,\\text{and}\\,"), "||"=>"\\,\\text{or}\\,")
        summation=="" || (summation = join(push!(symbols(indexes.internals, representation), summation), ",\\,"))
        context = isdefinite(indexes) ? "" : "\\sum_{$summation}"
        for index in indexes
            context = @sprintf "%s%s%s" context (length(context)>0 ? " " : "") latexstring(index)
        end
        push!(result, context)
        count = stop+1
    end
    return join(result, " \\cdot ")
end

"""
    @pattern index₁ index₂ ...
    @pattern(index₁, index₂, ...; constraint=...)

Construct a coupling pattern according to the pattern of the input indexes and an optional constraint.
"""
macro pattern(exprs...)
    if exprs[1].head==:parameters
        @assert(
            length(exprs[1].args)==1 && exprs[1].args[1].head==:kw && exprs[1].args[1].args[1]==:constraint,
            "@pattern error: constraint must be specified by the keyword argument `constraint`."
        )
        constraint = exprs[1].args[1].args[2]
        patterns = exprs[2:end]
    else
        constraint = true
        patterns = exprs
    end
    indexes = []
    blocks = []
    for (i, pattern) in enumerate(patterns)
        if pattern.head==:call && length(pattern.args)==3 && isa(pattern.args[3], Expr) && pattern.args[3].head==:call && pattern.args[3].args[1]≠://
            candidates = pattern.args[3].args[2:end]
            flag = 1
        elseif pattern.head==:call && length(pattern.args)>=3
            candidates = pattern.args[3:end]
            flag = 2
        else
            error("@pattern error: wrong pattern.")
        end
        attrs = []
        for (j, attr) in enumerate(candidates)
            isa(attr, Expr) && !(attr.head==:call && length(attr.args)==3 && attr.args[1]==:// && isa(attr.args[2], Int) && isa(attr.args[3], Int)) && (attr = Symbol(attr))
            if isa(attr, Symbol)
                @assert attr≠Symbol(":") "@pattern error: Colon(:) cannot be used for internal indexes in @pattern."
                @assert Base.isidentifier(attr) "@pattern error: wrong pattern."
                if !isdefined(__module__, attr)
                    push!(blocks, quote
                        local $attr
                        if @isdefined($attr)
                            ($attr)!=getfield(index[$i], $j) && return false
                        else
                            $attr = getfield(index[$i], $j)
                        end
                    end)
                    attr = QuoteNode(attr)
                end
            else
                push!(blocks, quote
                    (getfield(index[$i], $j)!=$attr) && return false
                end)
            end
            push!(attrs, attr)
        end
        if flag == 1
            index = Expr(:call, pattern.args[1], pattern.args[2], Expr(:call, pattern.args[3].args[1], attrs...))
        elseif flag == 2
            index = Expr(:call, pattern.args[1], pattern.args[2], attrs...)
        end
        push!(indexes, index)
    end
    N = length(indexes)
    indexes = Expr(:tuple, indexes...)
    blocks = Expr(:block, vcat([block.args for block in blocks]...)...)
    iname, fname = gensym("indexes"), gensym("constraint")
    representation = constraint==true ? "" : string(constraint)
    return quote
        function ($fname)(index::InternalIndexProd{<:NTuple{$N, SimpleInternalIndex}})
            $blocks
            return $constraint
        end
        local $iname = $(esc(indexes))
        Pattern($iname.sites, InternalPattern(InternalIndexProd($iname.internals), $fname, $representation))
    end
end

# patternrule
"""
    patternrule(value, ::Val{Name}, args...; kwargs...) where Name

By overriding this function, a named rule for the value of some attributes in a pattern can be defined so that the value can be transformed into the desired one.

In most cases, such overridden functions define the default rules for the attributes in a pattern when they take on the default value `:`. 
"""
function patternrule end

"""
    patternrule(value, ::Val, args...; kwargs...)-> typeof(value)

Default pattern rule unless otherwise specified.
"""
@inline patternrule(value, ::Val, args...; kwargs...) = value

"""
    patternrule(index::InternalIndexProd{T}, ::Val{Name}) where {N, T<:NTuple{N, SimpleInternalIndex}, Name} -> InternalIndexProd

Define the default rule for the internal index in an internal pattern.
"""
@inline @generated function patternrule(index::InternalIndexProd{T}, ::Val{Name}) where {N, T<:NTuple{N, SimpleInternalIndex}, Name}
    allequal(map(nameof, fieldtypes(T))) || return :(index)
    exprs = []
    for name in QuoteNode.(fieldnames(eltype(T)))
        vs = Expr(:tuple, [:(getfield(index[$i], $name)) for i = 1:N]...)
        push!(exprs, :(patternrule($vs, Val($(QuoteNode(Name))), eltype(T), Val($name))))
    end
    return :(InternalIndexProd(map(apply, fieldtypes(T), $(exprs...))))
end
@inline apply(::Type{F}, args...) where F = F(args...)

"""
    patternrule(sites::NTuple{N, Colon}, ::Val, bondlength::Integer) where N -> NTuple{N, Ordinal}

Define the default rule for the sites of a set of indexes in a coupling pattern.
"""
function patternrule(::NTuple{N, Colon}, ::Val, bondlength::Integer) where N
    bondlength==1 && return ntuple(i->1ˢᵗ, Val(N))
    bondlength==2 && begin
        N==2 && return (1ˢᵗ, 2ⁿᵈ)
        N==4 && return (1ˢᵗ, 1ˢᵗ, 2ⁿᵈ, 2ⁿᵈ)
        error("patternrule error: not implemented.")
    end
    error("patternrule error: not implemented.")
end

# Coupling
"""
    Coupling{V, P<:Pattern} <: OperatorPack{V, P}

Coupling among internal degrees of freedom at different or same lattice points.
"""
struct Coupling{V, P<:Pattern} <: OperatorPack{V, P}
    value::V
    pattern::P
end
@inline getcontent(coupling::Coupling, ::Val{:id}) = coupling.pattern
@inline Base.length(coupling::Coupling) = length(typeof(coupling))
@inline Base.length(::Type{<:Coupling}) = 1
@inline Base.eltype(coupling::Coupling) = eltype(typeof(coupling))
@inline Base.eltype(C::Type{<:Coupling}) = C
@inline Base.iterate(coupling::Coupling) = (coupling, nothing)
@inline Base.iterate(::Coupling, ::Nothing) = nothing
@inline Base.show(io::IO, coupling::Coupling) = @printf io "%s%s" (coupling.value≈1 ? "" : coupling.value≈-1 ? "- " : string(tostr(coupling.value), " ")) coupling.pattern

"""
    Coupling(indexes::Index...)
    Coupling(value, indexes::Index...)
    Coupling(value, indexes::Tuple{Vararg{Index}})

Construct a coupling with the input indexes as the pattern.
"""
@inline Coupling(indexes::Index...) = Coupling(1, indexes...)
@inline Coupling(value, indexes::Index...) = Coupling(value, indexes)
@inline Coupling(value, indexes::Tuple{Vararg{Index}}) = Coupling(value, Pattern(indexes))

"""
    Coupling(sites::Union{NTuple{N, Ordinal}, Colon}, ::Type{I}, fields::Union{NTuple{N}, Colon}...) where {N, I<:SimpleInternalIndex}
    Coupling(value, sites::Union{NTuple{N, Ordinal}, Colon}, ::Type{I}, fields::Union{NTuple{N}, Colon}...) where {N, I<:SimpleInternalIndex}
    Coupling{N}(sites::Union{NTuple{N, Ordinal}, Colon}, ::Type{I}, fields::Union{NTuple{N}, Colon}...) where {N, I<:SimpleInternalIndex}
    Coupling{N}(value, sites::Union{NTuple{N, Ordinal}, Colon}, ::Type{I}, fields::Union{NTuple{N}, Colon}...) where {N, I<:SimpleInternalIndex}

    Coupling(sites::Union{NTuple{N, Ordinal}, Colon}, f::Function, fields::Union{NTuple{N}, Colon}...) where N
    Coupling(value, sites::Union{NTuple{N, Ordinal}, Colon}, f::Function, fields::Union{NTuple{N}, Colon}...) where N
    Coupling{N}(sites::Union{NTuple{N, Ordinal}, Colon}, f::Function, fields::Union{NTuple{N}, Colon}...) where N
    Coupling{N}(value, sites::Union{NTuple{N, Ordinal}, Colon}, f::Function, fields::Union{NTuple{N}, Colon}...) where N

Construct a `Coupling` with the input sites and the fields of a kind of simple internal index.
"""
@inline Coupling(sites::Union{NTuple{N, Ordinal}, Colon}, f::Function, fields::Union{NTuple{N}, Colon}...) where N = Coupling{N}(sites, f, fields...)
@inline Coupling(value, sites::Union{NTuple{N, Ordinal}, Colon}, f::Function, fields::Union{NTuple{N}, Colon}...) where N = Coupling{N}(value, sites, f, fields...)
@inline Coupling{N}(sites::Union{NTuple{N, Ordinal}, Colon}, f::Function, fields::Union{NTuple{N}, Colon}...) where N = Coupling{N}(1, sites, f, fields...)
@inline Coupling{N}(value, sites::Union{NTuple{N, Ordinal}, Colon}, f::Function, fields::Union{NTuple{N}, Colon}...) where N = Coupling{N}(value, sites, AbstractIndex[f], fields...)
@inline Coupling(sites::Union{NTuple{N, Ordinal}, Colon}, ::Type{I}, fields::Union{NTuple{N}, Colon}...) where {N, I<:SimpleInternalIndex} = Coupling{N}(sites, I, fields...)
@inline Coupling(value, sites::Union{NTuple{N, Ordinal}, Colon}, ::Type{I}, fields::Union{NTuple{N}, Colon}...) where {N, I<:SimpleInternalIndex} = Coupling{N}(value, sites, I, fields...)
@inline Coupling{N}(sites::Union{NTuple{N, Ordinal}, Colon}, ::Type{I}, fields::Union{NTuple{N}, Colon}...) where {N, I<:SimpleInternalIndex} = Coupling{N}(1, sites, I, fields...)
@inline function Coupling{N}(value, sites::Union{NTuple{N, Ordinal}, Colon}, ::Type{I}, fields::Union{NTuple{N}, Colon}...) where {N, I<:SimpleInternalIndex}
    return Coupling(value, map(Index, default(sites, Val(N)), map(I, map(field->default(field, Val(N)), fields)...)))
end
@inline default(fields, ::Val) = fields
@inline default(::Colon, N::Val) = ntuple(i->:, N)

"""
    rank(coupling::Coupling) -> Int
    rank(::Type{M}) where {M<:Coupling} -> Int

Get the rank of a coupling.
"""
@inline rank(coupling::Coupling) = rank(typeof(coupling))
@inline rank(::Type{<:Coupling{V, P} where V}) where {P<:Pattern} = rank(P)

"""
    *(cp₁::Coupling, cp₂::Coupling) -> Coupling
    ⊗(cp₁::Coupling, cp₂::Coupling) -> Coupling

Get the multiplication between two coupling.
"""
@inline ⊗(cp₁::Coupling, cp₂::Coupling) = cp₁ * cp₂
@inline Base.:*(cp₁::Coupling, cp₂::Coupling) = Coupling(cp₁.value*cp₂.value, cp₁.pattern⊗cp₂.pattern)

"""
    latexstring(coupling::Coupling) -> String

Convert a coupling to the latex format.
"""
@inline latexstring(coupling::Coupling) = @sprintf "%s%s" (coupling.value≈1 ? "" : coupling.value≈-1 ? "- " : string(tostr(coupling.value), "\\,")) latexstring(coupling.pattern)

"""
    expand(coupling::Coupling, ::Val{Rule}, bond::Bond, hilbert::Hilbert) where Rule

Expand a coupling with the given bond and Hilbert space under a given named pattern rule.
"""
function expand(coupling::Coupling, ::Val{Rule}, bond::Bond, hilbert::Hilbert) where Rule
    sites = patternrule(coupling.pattern.sites, Val(Rule), length(bond))
    index = patternrule(coupling.pattern.internal.index, Val(Rule))
    points = ntuple(i->bond[sites[i]], Val(rank(coupling)))
    internal = InternalProd(map((index, internal)->filter(index, internal), index.contents, ntuple(i->hilbert[points[i].site], Val(rank(coupling)))))
    pattern = InternalPattern{partition(coupling.pattern.internal)}(index, coupling.pattern.internal.constraints, coupling.pattern.internal.representations)
    return CExpand(coupling.value, points, ConstrainedInternal(internal, pattern))
end
struct CExpand{V, N, SV<:SVector, S<:ConstrainedInternal}
    value::V
    sites::NTuple{N, Int}
    rcoordinates::NTuple{N, SV}
    icoordinates::NTuple{N, SV}
    internal::S
end
function CExpand(value, points::NTuple{N, Point}, internal::ConstrainedInternal) where N
    sites = ntuple(i->points[i].site, Val(N))
    rcoordinates = ntuple(i->points[i].rcoordinate, Val(N))
    icoordinates = ntuple(i->points[i].icoordinate, Val(N))
    return CExpand(value, sites, rcoordinates, icoordinates, internal)
end
@inline Base.eltype(ex::CExpand) = eltype(typeof(ex))
@inline @generated function Base.eltype(::Type{<:CExpand{V, N, SV, S}}) where {V, N, SV<:SVector, S<:ConstrainedInternal}
    return Operator{V, Tuple{map(I->CoordinatedIndex{Index{I, Int}, SV}, fieldtypes(fieldtype(eltype(S), :contents)))...}}
end
@inline Base.IteratorSize(::Type{<:CExpand}) = Base.SizeUnknown()
function Base.iterate(ex::CExpand, state=iterate(ex.internal))
    while !isnothing(state)
        index, state = state
        state = iterate(ex.internal, state)
        return Operator(ex.value, map(CoordinatedIndex, map(Index, ex.sites, index.contents), ex.rcoordinates, ex.icoordinates)), state
    end
    return
end

# MatrixCoupling
"""
    Component{T₁, T₂} <: VectorSpace{Tuple{T₁, T₁, T₂}}

A component of a matrix coupling, i.e., a matrix acting on a separated internal space.
"""
struct Component{T₁, T₂, V<:AbstractVector{T₁}} <: VectorSpace{Tuple{T₁, T₁, T₂}}
    left::V
    right::V
    matrix::SparseMatrixCSC{T₂, Int}
    function Component(left::V, right::V, matrix::AbstractMatrix{T}) where {V<:AbstractVector, T}
        @assert length(left)==size(matrix)[1] && length(right)==size(matrix)[2] "Component error: mismatched inputs."
        new{eltype(V), T, V}(left, right, convert(SparseMatrixCSC{T, Int}, matrix))
    end
end
@inline Base.length(component::Component) = nnz(component.matrix)
@inline function Base.getindex(component::Component, i::Integer)
    row = component.matrix.rowval[i]
    col = searchsortedlast(component.matrix.colptr, i)
    return (component.left[row], component.right[col], component.matrix.nzval[i])
end
@inline Base.promote_rule(::Type{Component{T, T₁, V}}, ::Type{Component{T, T₂, V}}) where {T, T₁, T₂, V<:AbstractVector{T}} = Component{T, promote_type(T₁, T₂), V}
@inline function Base.convert(::Type{Component{T, T₁, V}}, component::Component{T, T₂, V}) where {T, T₁, T₂, V<:AbstractVector{T}}
    return Component(component.left, component.right, convert(SparseMatrixCSC{T₁, Int}, component.matrix))
end

"""
    MatrixCoupling{I}(sites::Union{NTuple{2, Ordinal}, NTuple{2, Colon}}, contents::Tuple{Vararg{Component}}) where {I<:SimpleInternalIndex}
    MatrixCoupling(sites::Union{NTuple{2, Ordinal}, Colon}, ::Type{I}, contents::Component...) where {I<:SimpleInternalIndex}

Matrix coupling, i.e., a set of couplings whose coefficients are specified by matrices acting on separated internal spaces.
"""
struct MatrixCoupling{I<:SimpleInternalIndex, S<:Union{Ordinal, Colon}, C<:Tuple{Vararg{Component}}} <: VectorSpace{Coupling}
    sites::Tuple{S, S}
    contents::C
    function MatrixCoupling{I}(sites::Union{NTuple{2, Ordinal}, NTuple{2, Colon}}, contents::Tuple{Vararg{Component}}) where {I<:SimpleInternalIndex}
        @assert fieldcount(I)==length(contents) "MatrixCoupling error: mismatched type of internal index ($nameof(I)) and components (len=$length(contents))."
        new{I, eltype(sites), typeof(contents)}(sites, contents)
    end
end
@inline MatrixCoupling(sites::Union{NTuple{2, Ordinal}, Colon}, ::Type{I}, contents::Component...) where {I<:SimpleInternalIndex} = MatrixCoupling{I}(default(sites, Val(2)), contents)
@inline @generated function Base.eltype(::Type{MC}) where {MC<:MatrixCoupling}
    types = fieldtypes(fieldtype(MC, :contents))
    V = Expr(:call, :promote_type, [:(parametertype($C, 2)) for C in types]...)
    S = parametertype(MC, 2)
    I = Expr(:call, :indextype, parametertype(MC, 1), [:(parametertype($C, 1)) for C in types]...)
    C = :(InternalPattern{(2,), Tuple{$I, $I}, 1, Tuple{typeof(AllEqual($I))}})
    return :(Coupling{$V, Pattern{Tuple{$S, $S}, $C}})
end
@inline VectorSpaceStyle(::Type{<:MatrixCoupling}) = VectorSpaceDirectProducted(:forward)
function Base.convert(::Type{<:Coupling}, contents::Tuple, mc::MatrixCoupling{I}) where {I<:SimpleInternalIndex}
    value = mapreduce(content->content[3], *, contents, init=1)
    index₁ = Index(mc.sites[1], I(map(content->content[1], contents)...))
    index₂ = Index(mc.sites[2], I(map(content->content[2], contents)...))
    return Coupling(value, index₁, index₂)
end
@inline function Base.promote_rule(::Type{MatrixCoupling{I, S, C₁}}, ::Type{MatrixCoupling{I, S, C₂}}) where {I<:SimpleInternalIndex, S<:Union{Ordinal, Colon}, C₁<:Tuple{Vararg{Component}}, C₂<:Tuple{Vararg{Component}}}
    return MatrixCoupling{I, S, _promote_type_(C₁, C₂)}
end
@inline Base.convert(::Type{MatrixCoupling{I, S, C}}, mc::MatrixCoupling{I, S}) where {I<:SimpleInternalIndex, S<:Union{Ordinal, Colon}, C<:Tuple{Vararg{Component}}} = MatrixCoupling{I}(mc.sites, convert(C, mc.contents))
@inline @generated function _promote_type_(::Type{T₁}, ::Type{T₂}) where {N, T₁<:NTuple{N, Any}, T₂<:NTuple{N, Any}}
    return Expr(:curly, Tuple, [:(promote_type(fieldtype(T₁, $i), fieldtype(T₂, $i))) for i=1:N]...)
end

"""
    MatrixCouplingProd{V<:Number, C<:Tuple{Vararg{MatrixCoupling}}} <: VectorSpace{Coupling}

Product of matrix couplings together with an overall coefficient.
"""
struct MatrixCouplingProd{V<:Number, C<:Tuple{Vararg{MatrixCoupling}}} <: VectorSpace{Coupling}
    value::V
    contents::C
end
@inline MatrixCouplingProd(contents::MatrixCoupling...) = MatrixCouplingProd(1, contents)
@inline MatrixCouplingProd(value::Number, contents::MatrixCoupling...) = MatrixCouplingProd(value, contents)
@inline Base.eltype(::Type{MCP}) where {MCP<:MatrixCouplingProd} = _eltype_(fieldtype(MCP, :value), _eltypes_(MCP))
@inline @generated _eltypes_(::Type{MCP}) where {MCP<:MatrixCouplingProd} = Expr(:curly, :Tuple, [:(eltype($C)) for C in fieldtypes(fieldtype(MCP, :contents))]...)
@inline @generated function _eltype_(::Type{V}, ::Type{TS}) where {V<:Number, TS<:Tuple}
    types = fieldtypes(TS)
    MV = promote_type(V, map(valtype, types)...)
    SS = Tuple{concatenate(map(C->fieldtypes(parametertype(idtype(C), :sites)), types)...)...}
    IS = Tuple{concatenate(map(C->fieldtypes(parametertype(parametertype(idtype(C), :internal), 2)), types)...)...}
    FS = Tuple{concatenate(map(C->fieldtypes(parametertype(parametertype(idtype(C), :internal), 4)), types)...)...}
    N = length(types)
    RS = ntuple(i->2, Val(N))
    P = InternalPattern{RS, IS, N, FS}
    return Coupling{MV, Pattern{SS, P}}
end
@inline VectorSpaceStyle(::Type{<:MatrixCouplingProd}) = VectorSpaceDirectProducted(:forward)
@inline Base.convert(::Type{<:Coupling}, contents::Tuple{Vararg{Coupling}}, mcp::MatrixCouplingProd) = prod(contents; init=mcp.value)
@inline function Base.promote_rule(::Type{MatrixCouplingProd{V₁, C₁}}, ::Type{MatrixCouplingProd{V₂, C₂}}) where {V₁<:Number, C₁<:Tuple{Vararg{MatrixCoupling}}, V₂<:Number, C₂<:Tuple{Vararg{MatrixCoupling}}}
    return MatrixCouplingProd{promote_type(V₁, V₂), _promote_type_(C₁, C₂)}
end
@inline Base.convert(::Type{MatrixCouplingProd{V, C}}, mcp::MatrixCouplingProd) where {V<:Number, C<:Tuple{Vararg{MatrixCoupling}}} = MatrixCouplingProd(convert(V, mcp.value), convert(C, mcp.contents))

"""
    MatrixCouplingSum{C<:MatrixCouplingProd, N} <: VectorSpace{Coupling}

Sum of the products of matrix couplings.
"""
struct MatrixCouplingSum{C<:MatrixCouplingProd, N} <: VectorSpace{Coupling}
    contents::NTuple{N, C}
end
@inline MatrixCouplingSum(contents::Union{MatrixCoupling, MatrixCouplingProd}...) = MatrixCouplingSum(promote(map(m->1*m, contents)...))
@inline Base.eltype(::Type{<:MatrixCouplingSum{C}}) where {C<:MatrixCouplingProd} = eltype(C)
@inline VectorSpaceStyle(::Type{<:MatrixCouplingSum}) = VectorSpaceDirectSummed()

"""
    *(mc₁::MatrixCoupling, mc₂::MatrixCoupling) -> MatrixCouplingProd
    *(mc::MatrixCoupling, mcp::MatrixCouplingProd) -> MatrixCouplingProd
    *(mcp::MatrixCouplingProd, mc::MatrixCoupling) -> MatrixCouplingProd
    *(mcp₁::MatrixCouplingProd, mcp₂::MatrixCouplingProd) -> MatrixCouplingProd
    *(mc::MatrixCoupling, factor::Number) -> MatrixCouplingProd
    *(factor::Number, mc::MatrixCoupling) -> MatrixCouplingProd
    *(factor::Number, mcp::MatrixCouplingProd) -> MatrixCouplingProd
    *(mcp::MatrixCouplingProd, factor::Number) -> MatrixCouplingProd
    *(mcs::MatrixCouplingSum, element::Union{Number, MatrixCoupling, MatrixCouplingProd}) -> MatrixCouplingSum
    *(element::Union{Number, MatrixCoupling, MatrixCouplingProd}, mcs::MatrixCouplingSum) -> MatrixCouplingSum
    *(mcs₁::MatrixCouplingSum, mcs₂::MatrixCouplingSum) -> MatrixCouplingSum

Product between `MatrixCoupling`s and `MatrixCouplingProd`s.
"""
@inline Base.:*(mc₁::MatrixCoupling, mc₂::MatrixCoupling) = MatrixCouplingProd(mc₁, mc₂)
@inline Base.:*(mc::MatrixCoupling, mcp::MatrixCouplingProd) = MatrixCouplingProd(mcp.value, mc, mcp.contents...)
@inline Base.:*(mcp::MatrixCouplingProd, mc::MatrixCoupling) = MatrixCouplingProd(mcp.value, mcp.contents..., mc)
@inline Base.:*(mcp₁::MatrixCouplingProd, mcp₂::MatrixCouplingProd) = MatrixCouplingProd(mcp₁.value*mcp₂.value, mcp₁.contents..., mcp₂.contents...)
@inline Base.:*(factor::Number, mc::MatrixCoupling) = MatrixCouplingProd(factor, mc)
@inline Base.:*(mc::MatrixCoupling, factor::Number) = MatrixCouplingProd(factor, mc)
@inline Base.:*(factor::Number, mcp::MatrixCouplingProd) = MatrixCouplingProd(factor*mcp.value, mcp.contents)
@inline Base.:*(mcp::MatrixCouplingProd, factor::Number) = MatrixCouplingProd(factor*mcp.value, mcp.contents)
@inline Base.:*(mcs::MatrixCouplingSum, factor::Number) = MatrixCouplingSum(map(m->m*factor, mcs.contents))
@inline Base.:*(factor::Number, mcs::MatrixCouplingSum) = MatrixCouplingSum(map(m->m*factor, mcs.contents))
@inline Base.:*(mcs::MatrixCouplingSum, element::Union{MatrixCoupling, MatrixCouplingProd}) = MatrixCouplingSum(map(m->m*element, mcs.contents))
@inline Base.:*(element::Union{MatrixCoupling, MatrixCouplingProd}, mcs::MatrixCouplingSum) = MatrixCouplingSum(map(m->element*m, mcs.contents))
@inline Base.:*(mcs₁::MatrixCouplingSum, mcs₂::MatrixCouplingSum) = MatrixCouplingSum(concatenate(map(m₁->map(m₂->m₁*m₂, mcs₂.contents), mcs₁.contents)...))

"""
    +(mc₁::Union{MatrixCoupling, MatrixCouplingProd}, mc₂::Union{MatrixCoupling, MatrixCouplingProd}) -> MatrixCouplingSum
    +(mc::Union{MatrixCoupling, MatrixCouplingProd}, mcs::MatrixCouplingSum) -> MatrixCouplingSum
    +(mcs::MatrixCouplingSum, mc::Union{MatrixCoupling, MatrixCouplingProd}) -> MatrixCouplingSum
    +(mcs₁::MatrixCouplingSum, mcs₂::MatrixCouplingSum) -> MatrixCouplingSum

Addition between `MatrixCoupling`s and `MatrixCouplingProd`s.
"""
@inline Base.:+(mc₁::Union{MatrixCoupling, MatrixCouplingProd}, mc₂::Union{MatrixCoupling, MatrixCouplingProd}) = MatrixCouplingSum(mc₁, mc₂)
@inline Base.:+(mc::Union{MatrixCoupling, MatrixCouplingProd}, mcs::MatrixCouplingSum) = MatrixCouplingSum(mc, mcs.contents...)
@inline Base.:+(mcs::MatrixCouplingSum, mc::Union{MatrixCoupling, MatrixCouplingProd}) = MatrixCouplingSum(mcs.contents..., mc)
@inline Base.:+(mcs₁::MatrixCouplingSum, mcs₂::MatrixCouplingSum) = MatrixCouplingSum(mcs₁.contents..., mcs₂.contents...)

"""
    ^(mc::Union{MatrixCoupling, MatrixCouplingProd, MatrixCouplingSum}, n::Integer) -> Union{MatrixCoupling, MatrixCouplingProd, MatrixCouplingSum}

Get the nth power of a `MatrixCoupling`/`MatrixCouplingProd`/`MatrixCouplingSum`.
"""
@inline Base.:^(mc::Union{MatrixCoupling, MatrixCouplingProd, MatrixCouplingSum}, n::Integer) = prod(ntuple(i->mc, Val(n)); init=1)

"""
    /(mcp::MatrixCouplingProd, factor::Number) -> MatrixCouplingProd
    /(mcs::MatrixCouplingSum, factor::Number) -> MatrixCouplingSum
    /(mc::MatrixCoupling, factor::Number) -> MatrixCouplingProd
    //(mcp::MatrixCouplingProd, factor::Number) -> MatrixCouplingProd
    //(mcs::MatrixCouplingSum, factor::Number) -> MatrixCouplingSum
    //(mc::MatrixCoupling, factor::Number) -> MatrixCouplingProd
    -(mc::MatrixCoupling) -> MatrixCouplingProd
    -(mcp::MatrixCouplingProd) -> MatrixCouplingProd
    -(mcs::MatrixCouplingSum) -> MatrixCouplingSum
    -(mc₁::Union{MatrixCoupling, MatrixCouplingProd}, mc₂::Union{MatrixCoupling, MatrixCouplingProd}) -> MatrixCouplingSum
    -(mc::Union{MatrixCoupling, MatrixCouplingProd}, mcs::MatrixCouplingSum) -> MatrixCouplingSum
    -(mcs::MatrixCouplingSum, mc::Union{MatrixCoupling, MatrixCouplingProd}) -> MatrixCouplingSum
    -(mcs₁::MatrixCouplingSum, mcs₂::MatrixCouplingSum) -> MatrixCouplingSum

Define right-division, minus and subtraction operator for a `MatrixCoupling`/`MatrixCouplingProd`/`MatrixCouplingSum`.
"""
@inline Base.:/(mcp::MatrixCouplingProd, factor::Number) = MatrixCouplingProd(mcp.value/factor, mcp.contents)
@inline Base.:/(mcs::MatrixCouplingSum, factor::Number) = MatrixCouplingSum(map(m->m/factor, mcs.contents))
@inline Base.:/(mc::MatrixCoupling, factor::Number) = MatrixCouplingProd(1/factor, mc)
@inline Base.://(mcp::MatrixCouplingProd, factor::Number) = MatrixCouplingProd(mcp.value//factor, mcp.contents)
@inline Base.://(mcs::MatrixCouplingSum, factor::Number) = MatrixCouplingSum(map(m->m//factor, mcs.contents))
@inline Base.://(mc::MatrixCoupling, factor::Number) = MatrixCouplingProd(1//factor, mc)
@inline Base.:-(mc::MatrixCoupling) = MatrixCouplingProd(-1, mc)
@inline Base.:-(mcp::MatrixCouplingProd) = MatrixCouplingProd(-1*mcp.value, mcp.contents)
@inline Base.:-(mcs::MatrixCouplingSum) = MatrixCouplingSum(map(m->-m, mcs.contents))
@inline Base.:-(mc₁::Union{MatrixCoupling, MatrixCouplingProd}, mc₂::Union{MatrixCoupling, MatrixCouplingProd}) = MatrixCouplingSum(mc₁, -mc₂)
@inline Base.:-(mc::Union{MatrixCoupling, MatrixCouplingProd}, mcs::MatrixCouplingSum) = MatrixCouplingSum(mc, map(m->-m, mcs.contents)...)
@inline Base.:-(mcs::MatrixCouplingSum, mc::Union{MatrixCoupling, MatrixCouplingProd}) = MatrixCouplingSum(mcs.contents..., -mc)
@inline Base.:-(mcs₁::MatrixCouplingSum, mcs₂::MatrixCouplingSum) = MatrixCouplingSum(mcs₁.contents..., map(m->-m, mcs₂.contents)...)

# Term
"""
    TermFunction{F} <: Function

Abstract type for concrete term functions.
"""
abstract type TermFunction{F} <: Function end
@inline Base.:(==)(tf₁::TermFunction{F₁}, tf₂::TermFunction{F₂}) where {F₁, F₂} = F₁==F₂ && ==(efficientoperations, tf₁, tf₂)
@inline Base.isequal(tf₁::TermFunction{F₁}, tf₂::TermFunction{F₂}) where {F₁, F₂} = isequal(F₁, F₂) && isequal(efficientoperations, tf₁, tf₂)

"""
    TermAmplitude{F} <: TermFunction{F}

The function for the amplitude of a term.
"""
struct TermAmplitude{F} <: TermFunction{F}
    TermAmplitude(amplitude::Union{Function, Nothing}) = new{amplitude}()
end
@inline (::TermAmplitude{nothing})(::Bond) = 1
@inline (::TermAmplitude{F})(bond::Bond) where F = F(bond)
@inline Base.valtype(tf::TermFunction, bond::Bond) = valtype(typeof(tf), typeof(bond))
@inline Base.valtype(::Type{TermAmplitude{nothing}}, ::Type{<:Bond}) = Int
@inline Base.valtype(::Type{TermAmplitude{F}}, ::Type{B}) where {F, B<:Bond} = Core.Compiler.return_type(F, Tuple{B})

"""
    TermCoupling{C<:Coupling, F} <: TermFunction{F}

The function for the coupling of a term.
"""
struct TermCoupling{C<:Coupling, F} <: TermFunction{F}
    coupling::F
    TermCoupling(coupling) = new{eltype(coupling), typeof(coupling)}(coupling)
    TermCoupling(coupling::Function) = new{eltype(Core.Compiler.return_type(coupling, Tuple{Bond})), typeof(coupling)}(coupling)
end
@inline Base.valtype(termcoupling::TermCoupling) = valtype(typeof(termcoupling))
@inline Base.valtype(::Type{<:TermCoupling{C}}) where {C<:Coupling} = C
@inline (termcoupling::TermCoupling)(::Bond) = termcoupling.coupling
@inline (termcoupling::TermCoupling{<:Coupling, <:Function})(bond::Bond) = termcoupling.coupling(bond)

"""
    Term{K, I, V, B, C<:TermCoupling, A<:TermAmplitude}

Term of a quantum lattice system.
"""
mutable struct Term{K, I, V, B, C<:TermCoupling, A<:TermAmplitude}
    value::V
    const bondkind::B
    const coupling::C
    const amplitude::A
    const ishermitian::Bool
    const ismodulatable::Bool
    const factor::V
    function Term{K, I}(value, bondkind, coupling::TermCoupling, amplitude::TermAmplitude, ishermitian::Bool, ismodulatable::Bool, factor) where {K, I}
        @assert isa(K, Symbol) "Term error: kind must be a Symbol."
        @assert isa(I, Symbol) "Term error: id must be a Symbol."
        @assert value==value' "Term error: only real values are allowed. Complex values always have the positive directions, therefore, they must be specified by the amplitude function."
        new{K, I, typeof(value), typeof(bondkind), typeof(coupling), typeof(amplitude)}(value, bondkind, coupling, amplitude, ishermitian, ismodulatable, factor)
    end
end
@inline Base.:(==)(term₁::Term, term₂::Term) = ==(efficientoperations, term₁, term₂)
@inline Base.isequal(term₁::Term, term₂::Term) = isequal(efficientoperations, term₁, term₂)

"""
    Term{K}(id::Symbol, value, bondkind, coupling, ishermitian::Bool; amplitude::Union{Function, Nothing}=nothing, ismodulatable::Bool=true) where K

Construct a term.
"""
@inline function Term{K}(id::Symbol, value, bondkind, coupling, ishermitian::Bool; amplitude::Union{Function, Nothing}=nothing, ismodulatable::Bool=true) where K
    return Term{K, id}(value, bondkind, TermCoupling(coupling), TermAmplitude(amplitude), ishermitian, ismodulatable, 1)
end

"""
    kind(term::Term) -> Symbol
    kind(::Type{<:Term) -> Symbol

Get the kind of a term.
"""
@inline kind(term::Term) = kind(typeof(term))
@inline kind(::Type{<:Term{K}}) where K = K

"""
    id(term::Term) -> Symbol
    id(::Type{<:Term) -> Symbol

Get the id of a term.
"""
@inline id(term::Term) = id(typeof(term))
@inline id(::Type{<:Term{K, I} where K}) where I = I

"""
    valtype(term::Term)
    valtype(::Type{<:Term)

Get the value type of a term.
"""
@inline Base.valtype(term::Term) = valtype(typeof(term))
@inline Base.valtype(::Type{<:Term{K, I, V} where {K, I}}) where V = V

"""
    rank(term::Term) -> Int
    rank(::Type{<:Term) -> Int

Get the rank of a term.
"""
@inline rank(term::Term) = rank(typeof(term))
@inline rank(::Type{<:Term{K, I, V, B, C} where {K, I, V, B}}) where {C<:TermCoupling} = rank(valtype(C))

"""
    string(term::Term, bond::Bond, hilbert::Hilbert) -> String

Get the string representation of a term on a bond with a given Hilbert space.
"""
function Base.string(term::Term, bond::Bond, hilbert::Hilbert)
    cache = String[]
    if term.bondkind == bond.kind
        value = term.value * term.amplitude(bond) * term.factor
        if !isapprox(value, 0.0, atol=atol, rtol=rtol)
            for coupling in term.coupling(bond)
                if !isnothing(iterate(expand(coupling, Val(kind(term)), bond, hilbert)))
                    representation = string(value * coupling)
                    term.ishermitian || (representation = string(representation, " + h.c."))
                    push!(cache, @sprintf "%s" representation)
                end
            end
        end
    end
    return join(cache, "\n")
end

"""
    replace(term::Term; kwargs...) -> Term

Replace some attributes of a term with key word arguments.
"""
@inline @generated function Base.replace(term::Term; kwargs...)
    exprs = [:(get(kwargs, $name, getfield(term, $name))) for name in QuoteNode.(term|>fieldnames)]
    return :(Term{kind(term), id(term)}($(exprs...)))
end

"""
    one(term::Term) -> Term

Get a unit term.
"""
@inline Base.one(term::Term) = replace(term, value=one(term.value))

"""
    zero(term::Term) -> Term

Get a zero term.
"""
@inline Base.zero(term::Term) = replace(term, value=zero(term.value))

"""
    update!(term::Term, args...; kwargs...) -> Term

Update the value of a term if it is ismodulatable.
"""
function update!(term::Term, args...; kwargs...)
    @assert term.ismodulatable "update! error: not modulatable term."
    term.value = get(kwargs, id(term), term.value)
    return term
end

"""
    optype(::Type{T}, ::Type{H}, ::Type{B}) where {T<:Term, H<:Hilbert, B<:Bond}

Get the compatible `Operator` type from the type of a term, a Hilbert space and a bond.
"""
@inline function optype(::Type{T}, ::Type{H}, ::Type{B}) where {T<:Term, H<:Hilbert, B<:Bond}
    C = valtype(fieldtype(T, :coupling))
    @assert C<:Coupling "optype error: not supported."
    V, V′, V′′ = valtype(T), valtype(C), valtype(fieldtype(T, :amplitude), B)
    isconcretetype(V′) && (V = promote_type(V, V′))
    isconcretetype(V′′) && (V = promote_type(V, V′′))
    indextypes = ntuple(i->indextype(filter(fieldtype(parametertype(parametertype(idtype(C), :internal), 2), i), valtype(H)), eltype(B)), Val(rank(C)))
    return fulltype(Operator, NamedTuple{(:value, :id), Tuple{V, Tuple{indextypes...}}})
end

"""
    expand!(operators::Operators, term::Term, bond::Bond, hilbert::Hilbert; half::Bool=false) -> Operators
    expand!(operators::Operators, term::Term, bonds, hilbert::Hilbert; half::Bool=false) -> Operators

Expand the operators of a term on a bond/set-of-bonds with a given Hilbert space.

The `half` parameter determines the behavior of generating operators, which falls into the following two categories
* `true`: "Hermitian half" of the generated operators
* `false`: "Hermitian whole" of the generated operators
"""
function expand!(operators::Operators, term::Term, bond::Bond, hilbert::Hilbert; half::Bool=false)
    if term.bondkind == bond.kind
        value = term.value * term.amplitude(bond) * term.factor
        if !isapprox(value, 0.0, atol=atol, rtol=rtol)
            for coupling in term.coupling(bond)
                for opt in expand(coupling, Val(kind(term)), bond, hilbert)
                    isapprox(opt.value, 0.0, atol=atol, rtol=rtol) && continue
                    if half
                        add!(operators, opt*valtype(eltype(operators))(value/(term.ishermitian ? 2 : 1)))
                    else
                        opt = opt*valtype(eltype(operators))(value)
                        add!(operators, opt)
                        term.ishermitian || add!(operators, opt')
                    end
                end
            end
        end
    end
    return operators
end
@inline function expand!(operators::Operators, term::Term, bonds, hilbert::Hilbert; half::Bool=false)
    for bond in bonds
        expand!(operators, term, bond, hilbert; half=half)
    end
    return operators
end

"""
    expand(term::Term, bond::Bond, hilbert::Hilbert; half::Bool=false) -> Operators
    expand(term::Term, bonds, hilbert::Hilbert; half::Bool=false) -> Operators

Expand the operators of a term on a bond/set-of-bonds with a given Hilbert space.
"""
@inline function expand(term::Term, bond::Bond, hilbert::Hilbert; half::Bool=false)
    M = optype(term|>typeof, hilbert|>typeof, bond|>typeof)
    expand!(Operators{M}(), term, bond, hilbert; half=half)
end
@inline function expand(term::Term, bonds, hilbert::Hilbert; half::Bool=false)
    M = optype(term|>typeof, hilbert|>typeof, bonds|>eltype)
    expand!(Operators{M}(), term, bonds, hilbert; half=half)
end

# Metric and Table
"""
    Metric <: Function

The rules for measuring an operator unit so that different operator units can be compared.

As a function, every instance should accept only one positional argument, i.e. the operator unit to be measured.
"""
abstract type Metric <: Function end
@inline Base.:(==)(m₁::T, m₂::T) where {T<:Metric} = ==(efficientoperations, m₁, m₂)
@inline Base.isequal(m₁::T, m₂::T) where {T<:Metric} = isequal(efficientoperations, m₁, m₂)
@inline (M::Type{<:Metric})(::Type{I}) where {I<:CompositeIndex} = M(indextype(I))
@inline (M::Type{<:Metric})(::Type{H}) where {H<:Hilbert} = M(Index{H|>valtype|>eltype, Int})
@inline (metric::Metric)(index::CompositeIndex) = metric(getcontent(index, :index))
@inline Base.valtype(::Type{M}, ::Type{I}) where {M<:Metric, I<:CompositeIndex} = valtype(M, indextype(I))

"""
    OperatorUnitToTuple{Fields} <: Metric

A rule that converts an operator unit to a tuple by iterating over a set of selected fields in a specific order.
"""
struct OperatorUnitToTuple{Fields} <: Metric
    OperatorUnitToTuple(fields::Tuple{Vararg{Symbol}}) = new{fields}()
end
@inline OperatorUnitToTuple(fields::Symbol...) = OperatorUnitToTuple(fields)

"""
    keys(::OperatorUnitToTuple{Fields}) where Fields -> Fields
    keys(::Type{<:OperatorUnitToTuple{Fields}}) where Fields -> Fields

Get the names of the selected fields.
"""
@inline Base.keys(::OperatorUnitToTuple{Fields}) where Fields = Fields
@inline Base.keys(::Type{<:OperatorUnitToTuple{Fields}}) where Fields = Fields

"""
    OperatorUnitToTuple(::Type{I}) where {I<:Index}

Construct the metric rule from the information of the `Index` type.
"""
@inline OperatorUnitToTuple(::Type{I}) where {I<:Index} = OperatorUnitToTuple(:site, (fieldnames(indextype(I)))...)

"""
    valtype(::Type{<:OperatorUnitToTuple}, ::Type{<:Index})

Get the valtype of applying an `OperatorUnitToTuple` rule to an `Index`.
"""
@inline @generated function Base.valtype(::Type{M}, ::Type{I}) where {M<:OperatorUnitToTuple, I<:Index}
    types = []
    for field in keys(M)
        if field==:site
            push!(types, Int)
        elseif hasfield(indextype(I), field)
            push!(types, fieldtype(indextype(I), field))
        end
    end
    return  Expr(:curly, :Tuple, types...)
end

"""
    (operatorunittotuple::OperatorUnitToTuple)(index::Index) -> Tuple

Convert an index to a tuple.
"""
@inline @generated function (operatorunittotuple::OperatorUnitToTuple)(index::Index)
    exprs = []
    for name in keys(operatorunittotuple)
        field = QuoteNode(name)
        if name==:site
            push!(exprs, :(index.site))
        elseif hasfield(indextype(index), name)
            push!(exprs, :(getfield(index.internal, $field)))
        end
    end
    return Expr(:tuple, exprs...)
end

"""
    Table{I, B<:Metric} <: CompositeDict{I, Int}

The table of operator unit v.s. sequence pairs.
"""
struct Table{I, B<:Metric} <: CompositeDict{I, Int}
    by::B
    contents::OrderedDict{I, Int}
end
@inline contentnames(::Type{<:Table}) = (:by, :contents)
@inline Table{I}(by::Metric) where {I<:OperatorUnit} = Table(by, OrderedDict{valtype(typeof(by), I), Int}())
@inline vec2dict(vs::AbstractVector) = OrderedDict{eltype(vs), Int}(v=>i for (i, v) in enumerate(vs))

"""
    getindex(table::Table, operatorunit::OperatorUnit) -> Int

Inquiry the sequence of an operator unit.
"""
@inline Base.getindex(table::Table, operatorunit::OperatorUnit) = table[table.by(operatorunit)]

"""
    haskey(table::Table, operatorunit::OperatorUnit) -> Bool
    haskey(table::Table, operatorunits::ID{OperatorUnit}) -> Tuple{Vararg{Bool}}

Judge whether a single operator unit or a set of operator units have been assigned with sequences in table.
"""
@inline Base.haskey(table::Table, operatorunit::OperatorUnit) = haskey(table, table.by(operatorunit))
@inline Base.haskey(table::Table, operatorunits::ID{OperatorUnit}) = map(operatorunit->haskey(table, operatorunit), operatorunits)

"""
    Table(operatorunits::AbstractVector{<:OperatorUnit}, by::Metric=OperatorUnitToTuple(eltype(operatorunits)))

Convert a set of operator units to the corresponding table of operator unit vs. sequence pairs.

The input operator units are measured by the input `by` function with the duplicates removed. The resulting unique values are sorted, which determines the sequence of the input `operatorunits`. Note that two operator units have the same sequence if their converted values are equal to each other.
"""
@inline Table(operatorunits::AbstractVector{<:OperatorUnit}, by::Metric=OperatorUnitToTuple(eltype(operatorunits))) = Table(by, [by(operatorunit) for operatorunit in operatorunits]|>unique!|>sort!|>vec2dict)

"""
    Table(hilbert::Hilbert, by::Metric=OperatorUnitToTuple(typeof(hilbert))) -> Table

Get the index-sequence table of a Hilbert space.
"""
function Table(hilbert::Hilbert, by::Metric=OperatorUnitToTuple(typeof(hilbert)))
    result = Index{hilbert|>valtype|>eltype, Int}[]
    for (site, internal) in hilbert
        for index in internal
            push!(result, (result|>eltype)(site, index))
        end
    end
    return Table(result, by)
end

"""
    union(tables::Table...) -> Table

Unite several operator unit vs. sequence tables.
"""
function Base.union(tables::Table...)
    @assert allequal(map(table->table.by, tables)) "union error: all input tables should have the same `by` attribute."
    indices = (tables|>eltype|>keytype)[]
    for table in tables
        for index in keys(table)
            push!(indices, index)
        end
    end
    return Table(tables[1].by, indices|>unique!|>sort!|>vec2dict)
end

"""
    reset!(table::Table, operatorunits::AbstractVector{<:OperatorUnit}) -> Table

Reset a table by a new set of operatorunits.
"""
function reset!(table::Table, operatorunits::AbstractVector{<:OperatorUnit})
    empty!(table)
    for (i, id) in enumerate([table.by(operatorunit) for operatorunit in operatorunits]|>unique!|>sort!)
        table[id] = i
    end
    return table
end

"""
    reset!(table::Table, hilbert::Hilbert) -> Table

Reset a table by a Hilbert space.
"""
function reset!(table::Table, hilbert::Hilbert)
    indices = Index{hilbert|>valtype|>eltype, Int}[]
    for (site, internal) in pairs(hilbert)
        for internalindex in internal
            push!(indices, (indices|>eltype)(site, internalindex))
        end
    end
    return reset!(table, indices)
end

"""
    findall(select::Function, hilbert::Hilbert, table::Table) -> Vector{Int}

Find all the sequences of indexes contained in a Hilbert space according to a table and a select function.
"""
function Base.findall(select::Function, hilbert::Hilbert, table::Table)
    result = Int[]
    for (site, internal) in pairs(hilbert)
        for internalindex in internal
            index = Index(site, internalindex)
            if select(index)
                push!(result, table[index])
            end
        end
    end
    return sort!(unique!(result))
end

# Boundary
"""
    Boundary{Names}(values::AbstractVector{<:Number}, vectors::AbstractVector{<:AbstractVector{<:Number}}) where Names

Boundary twist of operators.
"""
struct Boundary{Names, D<:Number, V<:AbstractVector} <: LinearTransformation
    values::Vector{D}
    vectors::Vector{V}
    function Boundary{Names}(values::AbstractVector{<:Number}, vectors::AbstractVector{<:AbstractVector{<:Number}}) where Names
        @assert length(Names)==length(values)==length(vectors) "Boundary error: mismatched names, values and vectors."
        datatype = promote_type(eltype(values), Float)
        new{Names, datatype, eltype(vectors)}(convert(Vector{datatype}, values), vectors)
    end
end
@inline Base.:(==)(bound₁::Boundary, bound₂::Boundary) = keys(bound₁)==keys(bound₂) && ==(efficientoperations, bound₁, bound₂)
@inline Base.isequal(bound₁::Boundary, bound₂::Boundary) = isequal(keys(bound₁), keys(bound₂)) && isequal(efficientoperations, bound₁, bound₂)
@inline Base.valtype(::Type{<:Boundary}, M::Type{<:Operator}) = reparameter(M, :value, promote_type(Complex{Int}, valtype(M)))
@inline Base.valtype(B::Type{<:Boundary}, MS::Type{<:Operators}) = (M = valtype(B, eltype(MS)); Operators{M, idtype(M)})

"""
    keys(bound::Boundary) -> Tuple{Vararg{Symbol}}
    keys(::Type{<:Boundary{Names}}) where Names -> Names

Get the names of the boundary parameters.
"""
@inline Base.keys(bound::Boundary) = keys(typeof(bound))
@inline Base.keys(::Type{<:Boundary{Names}}) where Names = Names

"""
    (bound::Boundary)(operator::Operator; origin::Union{AbstractVector, Nothing}=nothing) -> Operator

Get the boundary twisted operator.
"""
@inline function (bound::Boundary)(operator::Operator; origin::Union{AbstractVector, Nothing}=nothing)
    values = isnothing(origin) ? bound.values : bound.values-origin
    return replace(operator, operator.value*exp(1im*mapreduce(u->angle(u, bound.vectors, values), +, id(operator))))
end

"""
    update!(bound::Boundary; parameters...) -> Boundary

Update the values of the boundary twisted phase.
"""
@inline @generated function update!(bound::Boundary; parameters...)
    exprs = []
    for (i, name) in enumerate(QuoteNode.(keys(bound)))
        push!(exprs, :(bound.values[$i] = get(parameters, $name, bound.values[$i])))
    end
    return Expr(:block, exprs..., :(return bound))
end

"""
    merge!(bound::Boundary, another::Boundary) -> typeof(bound)

Merge the values and vectors of the twisted boundary condition from another one.
"""
@inline function Base.merge!(bound::Boundary, another::Boundary)
    @assert keys(bound)==keys(another) "merge! error: mismatched names of boundary parameters."
    bound.values .= another.values
    bound.vectors .= another.vectors
    return bound
end

"""
    replace(bound::Boundary; values=bound.values, vectors=bound.vectors) -> Boundary

Replace the values or vectors of a twisted boundary condition and get the new one.

!!! note
    The plain boundary condition keeps plain even when replaced with new values or new vectors.
"""
@inline Base.replace(bound::Boundary; values=bound.values, vectors=bound.vectors) = Boundary{keys(bound)}(values, vectors)

"""
    plain

Plain boundary condition without any twist.
"""
const plain = Boundary{()}(Float[], SVector{0, Float}[])
@inline Base.valtype(::Type{typeof(plain)}, M::Type{<:Operator}) = M
@inline Base.valtype(::Type{typeof(plain)}, M::Type{<:Operators}) = M
@inline (::typeof(plain))(operator::Operator; kwargs...) = operator
@inline Base.replace(::(typeof(plain)); kwargs...) = plain

end #module
