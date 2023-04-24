module QuantumNumbers

using Base.Iterators: Reverse, flatten, product
using DataStructures: OrderedDict
using LinearAlgebra: norm
using Printf: @printf, @sprintf
using Random: MersenneTwister, seed!, shuffle!
using ..Toolkit: Combinations, HomoNamedVector, VectorSpace, VectorSpaceCartesian, VectorSpaceEnumerative, VectorSpaceStyle

import ..QuantumLattices: ⊕, ⊗, decompose, dimension, expand, permute
import ..Toolkit: shape

export AbelianNumber, AbelianNumbers, regularize, regularize!, periods, @abeliannumber
export Momenta, Momentum, Momentum₁, Momentum₂, Momentum₃, ParticleNumber, SpinfulParticle, SpinZ, particlenumbers, spinfulparticles, spinzs

"""
    positives(inputs::NTuple{N, Any}) where N -> NTuple{N, Int}

Return a tuple of all positive signs.
"""
@inline positives(inputs::NTuple{N, Any}) where N = ntuple(i->1, Val(N))

"""
    AbelianNumber{T<:Real} <: HomoNamedVector{T}

Abstract type for all concrete quantum numbers for a single basis.
"""
abstract type AbelianNumber{T<:Real} <: HomoNamedVector{T} end

"""
    periods(qn::AbelianNumber)

The periods of the components of a concrete `AbelianNumber`.
"""
@inline periods(qn::AbelianNumber) = periods(typeof(qn))

"""
    dimension(::Type{<:AbelianNumber}) -> Int
    dimension(::AbelianNumber) -> Int

The dimension of the Hilbert space an `AbelianNumber` represents. Apparently, this is always 1.
"""
@inline dimension(::Type{<:AbelianNumber}) = 1
@inline dimension(::AbelianNumber) = 1

"""
    +(qn::AbelianNumber) -> typeof(qn)
    +(qn::QN, qns::QN...) where {QN<:AbelianNumber} -> QN

Overloaded `+` operator for `AbelianNumber`.

!!! note
    To ensure type stability, two `AbelianNumber` can be added together if and only if they are of the same type.
"""
@inline Base.:+(qn::AbelianNumber) = qn
@inline Base.:+(qn::QN, qns::QN...) where {QN<:AbelianNumber} = map(+, qn, qns...)

"""
    -(qn::AbelianNumber) -> typeof(qn)
    -(qn₁::QN, qn₂::QN) where {QN<:AbelianNumber} -> QN

Overloaded `-` operator for `AbelianNumber`.

!!! note
    To ensure type stability, an `AbelianNumber` can be subtracted by another `AbelianNumber` if and only if they are of the same type.
"""
@inline Base.:-(qn::AbelianNumber) = map(-, qn)
@inline Base.:-(qn₁::QN, qn₂::QN) where {QN<:AbelianNumber} = map(-, qn₁, qn₂)

"""
    *(factor::Integer, qn::AbelianNumber) -> typeof(qn)
    *(qn::AbelianNumber, factor::Integer) -> typeof(qn)

Overloaded `*` operator for the multiplication between an integer and an `AbelianNumber`.
"""
@inline Base.:*(factor::Integer, qn::AbelianNumber) = qn * factor
@inline Base.:*(qn::AbelianNumber, factor::Integer) = typeof(qn)(map(num->num*factor, convert(Tuple, qn))...)

"""
    ^(qn::AbelianNumber, factor::Integer) -> typeof(qn)

Overloaded `^` operator for `AbelianNumber`.

This operation translates into the direct product of `factor` copies of `qn`.
"""
@inline Base.:^(qn::AbelianNumber, factor::Integer) = kron(NTuple{factor, typeof(qn)}(qn for i = 1:factor)...)

"""
    ⊗(qns::AbelianNumber...) -> eltype(qns)

Get the direct product of some `AbelianNumber`s.
"""
@inline ⊗(qns::AbelianNumber...) = kron(qns...)

"""
    kron(qns::Vararg{AbelianNumber}; signs=positives(qns))-> eltype(qns)

Get the direct product of some `AbelianNumber`s.
!!! note
    Physically, the direct product of a couple of `AbelianNumber`s are defined through the direct product of the bases of the Hilbert spaces they represent. Apparently, the result is still an `AbelianNumber` whose dimension is 1. At the same time, each component of the result is obtained by a summation of the corresponding components of the inputs with the correct signs. This is a direct consequence of the Abelian nature of our quantum numbers.
"""
@inline Base.kron(qns::Vararg{AbelianNumber}; signs=positives(qns)) = sum(sign==1 ? qn : -qn for (sign, qn) in zip(signs, qns))

"""
    regularize!(::Type{QN}, array::AbstractVector{<:Real}) where {QN<:AbelianNumber} -> typeof(array)
    regularize!(::Type{QN}, array::AbstractMatrix{<:Real}) where {QN<:AbelianNumber} -> typeof(array)

Regularize the elements of an array in place so that it can represent quantum numbers.
"""
function regularize!(::Type{QN}, array::AbstractVector{<:Real}) where {QN<:AbelianNumber}
    @assert array|>length == QN|>length "regularize! error: inconsistent shape of input array and $QN."
    for (i, period) in enumerate(QN|>periods)
        period==Inf || @inbounds begin
            remainder = array[i] % period
            array[i] = remainder<0 ? remainder+period : remainder
        end
    end
    return array
end
function regularize!(::Type{QN}, array::AbstractMatrix{<:Real}) where {QN<:AbelianNumber}
    @assert size(array, 1) == QN|>length "regularize! error: inconsistent shape of input array and $QN."
    for (i, period) in enumerate(QN|>periods)
        period==Inf || @inbounds for j = 1:size(array, 2)
            remainder = array[i, j] % period
            array[i, j] = remainder<0 ? remainder+period : remainder
        end
    end
    return array
end

"""
    regularize(::Type{QN}, array::Union{AbstractVector{<:Real}, AbstractMatrix{<:Real}}) where {QN<:AbelianNumber} -> typeof(array)

Regularize the elements of an array and return a copy that can represent quantum numbers.
"""
@inline regularize(::Type{QN}, array::Union{AbstractVector{<:Real}, AbstractMatrix{<:Real}}) where {QN<:AbelianNumber} = regularize!(QN, deepcopy(array))

"""
    @abeliannumber typename T fields periods

Construct a concrete `AbelianNumber` with the type name being `typename`, fieldtype specified by `T`, fieldnames specified by `fields`, and periods specified by `periods`.
"""
macro abeliannumber(typename, T, fields, periods)
    typename = Symbol(typename)
    fields = tuple(eval(fields)...)
    periods = tuple(eval(periods)...)
    arguments = ntuple(i->Symbol(:v, i), length(fields))
    @assert length(fields) == length(periods) "@abeliannumber error: number of fields($(length(fields))) and periods($(length(periods))) not equal."
    @assert all(isa(name, Symbol) for name in fields) "@abeliannumber error: all field names should be Symbol."
    @assert all(periods .> 0) "@abeliannumber error: all field periods should be positive."
    if all(periods .== Inf)
        title = Expr(:call, :($(esc(typename))), arguments...)
        body = Expr(:call, :new, arguments...)
    else
        title = Expr(:call, :($(esc(typename))), arguments..., Expr(:kw, :(regularize::Bool), :true))
        body = Expr(:call, :new, (p==Inf ? arg : :(regularize ? (r=$arg%$p; r<0 ? r+$p : r) : $arg) for (p, arg) in zip(periods, arguments))...)
    end
    newtype = Expr(:struct, false, :($(esc(typename))<:AbelianNumber{$(esc(T))}), Expr(:block, (:($field::$(esc(T))) for field in fields)..., Expr(:(=), title, body)))
    functions = Expr(:block, :(periods(::Type{<:$(esc(typename))}) = $periods))
    return Expr(:block, :(global periods), :(Base.@__doc__($newtype)), functions)
end

"""
    AbelianNumbers{QN<:AbelianNumber} <: VectorSpace{QN}

The whole quantum numbers of the total bases of a Hilbert space.
"""
struct AbelianNumbers{QN<:AbelianNumber} <: VectorSpace{QN}
    form::Char
    contents::Vector{QN}
    indptr::Vector{Int}
    function AbelianNumbers(form::Char, contents::Vector{<:AbelianNumber}, indptr::Vector{Int})
        @assert uppercase(form)∈('G', 'U', 'C') "AbelianNumbers error: 'form'($form) is not 'G', 'U' or 'C'."
        @assert length(contents)+1==length(indptr) "AbelianNumbers error: mismatch shapes of contents and indptr ($(length(contents))+1!=$length(indptr))."
        return new{contents|>eltype}(form|>uppercase, contents, indptr)
    end
end
@inline VectorSpaceStyle(::Type{<:AbelianNumbers}) = VectorSpaceEnumerative()
@inline dimension(qns::AbelianNumbers) = @inbounds(qns.indptr[end])
@inline Base.length(qns::AbelianNumbers) = length(qns.contents)
@inline Base.issorted(qns::AbelianNumbers) = qns.form=='C'
@inline Base.string(qns::AbelianNumbers) = @sprintf "QNS(%s, %s)" qns|>length qns|>dimension
@inline Base.show(io::IO, qns::AbelianNumbers) = @printf io "QNS(%s)" join((@sprintf("%s=>%s", qn, slice) for (qn, slice) in pairs(qns, :indptr)), ", ")

"""
    AbelianNumbers(form::Char, contents::Vector{<:AbelianNumber}, indptr::Vector{Int}, choice::Symbol)
    AbelianNumbers(form::Char, contents::Vector{<:AbelianNumber}, counts::Vector{Int}, ::Val{:counts})
    AbelianNumbers(form::Char, contents::Vector{<:AbelianNumber}, indptr::Vector{Int}, ::Val{:indptr})

Construct an `AbelianNumbers` from a vector of concrete quantum numbers and an vector containing their counts or indptr.
"""
@inline AbelianNumbers(form::Char, contents::Vector{<:AbelianNumber}, indptr::Vector{Int}, choice::Symbol) = AbelianNumbers(form, contents, indptr, choice|>Val)
@inline AbelianNumbers(form::Char, contents::Vector{<:AbelianNumber}, counts::Vector{Int}, ::Val{:counts}) = AbelianNumbers(form|>uppercase, contents, [0, cumsum(counts)...])
@inline AbelianNumbers(form::Char, contents::Vector{<:AbelianNumber}, indptr::Vector{Int}, ::Val{:indptr}) = AbelianNumbers(form|>uppercase, contents, indptr)

"""
    AbelianNumbers(form::Char, contents::Vector{<:AbelianNumber})

Construct an `AbelianNumbers` with a set of quantum numbers whose counts are all one.
"""
@inline function AbelianNumbers(form::Char, contents::Vector{<:AbelianNumber})
    return AbelianNumbers(form, contents, collect(0:length(contents)), Val(:indptr))
end

"""
    AbelianNumbers(qn::AbelianNumber, count::Int=1)

Construct an `AbelianNumbers` with one unique quantum number which occurs `count` times.
"""
@inline AbelianNumbers(qn::AbelianNumber, count::Int=1) = AbelianNumbers('C', [qn], [0, count], :indptr)

"""
    AbelianNumbers(od::OrderedDict{<:AbelianNumber, <:Integer})

Construct an `AbelianNumbers` from an ordered dict containing concrete quantum numbers and their counts.
"""
function AbelianNumbers(od::OrderedDict{<:AbelianNumber, <:Integer})
    contents = Vector{od|>keytype}(undef, length(od))
    indptr = zeros(Int, length(od)+1)
    for (i, (qn, count)) in enumerate(od)
        @inbounds contents[i] = qn
        @inbounds indptr[i+1] = indptr[i] + count
    end
    return AbelianNumbers('U', contents, indptr, :indptr)
end

"""
    AbelianNumbers(od::OrderedDict{<:AbelianNumber, <:UnitRange{<:Integer}})

Construct an `AbelianNumbers` from an ordered dict containing concrete quantum numbers and their slices.
"""
function AbelianNumbers(od::OrderedDict{<:AbelianNumber, <:UnitRange{<:Integer}})
    contents = Vector{od|>keytype}(undef, length(od))
    indptr = zeros(Int, length(od)+1)
    for (i, (qn, slice)) in enumerate(od)
        @inbounds contents[i] = qn
        @inbounds @assert indptr[i]+1==slice.start "AbelianNumbers error: slice not consistent."
        @inbounds indptr[i+1] = slice.stop
    end
    return AbelianNumbers('U', contents, indptr, :indptr)
end

"""
    getindex(qns::AbelianNumbers, slice::UnitRange{Int}) -> AbelianNumbers
    getindex(qns::AbelianNumbers, indices::Vector{Int}) -> AbelianNumbers

Overloaded `[]` operator.
!!! note
    For an `AbelianNumbers`, all the `getindex` functions act on its `contents`, i.e. its compressed data, but not on its expansion, i.e. the uncompressed data. This definition is consistent with the `length` of an `AbelianNumbers`.
"""
function Base.getindex(qns::AbelianNumbers, slice::UnitRange{Int})
    contents = qns.contents[slice]
    indptr = qns.indptr[slice.start:slice.stop+1]
    indptr .= indptr .- qns.indptr[slice.start]
    return AbelianNumbers(qns.form, contents, indptr, :indptr)
end
function Base.getindex(qns::AbelianNumbers, indices::Vector{Int})
    contents = Vector{qns|>eltype}(undef, length(indices))
    indptr = zeros(Int, length(indices)+1)
    for (i, index) in enumerate(indices)
        contents[i] = qns.contents[index]
        indptr[i+1] = indptr[i] + qns.indptr[index+1] - qns.indptr[index]
    end
    return AbelianNumbers('G', contents, indptr, :indptr)
end

"""
    keys(qns::AbelianNumbers) -> Vector{qns|>eltype}

Iterate over the concrete `AbelianNumber`s contained in an `AbelianNumbers`.
"""
@inline Base.keys(qns::AbelianNumbers) = qns.contents

"""
    values(qns::AbelianNumbers, choice::Symbol)
    values(qns::AbelianNumbers, ::Val{:indptr})
    values(qns::AbelianNumbers, ::Val{:counts})

Iterate over the slices/counts of the `AbelianNumbers`.
"""
@inline Base.values(qns::AbelianNumbers, choice::Symbol) = values(qns, choice|>Val)
@inline @views Base.values(qns::AbelianNumbers, ::Val{:indptr}) = ((start+1):stop for (start, stop) in zip(qns.indptr[1:end-1], qns.indptr[2:end]))
@inline @views Base.values(qns::AbelianNumbers, ::Val{:counts}) = (stop-start for (start, stop) in zip(qns.indptr[1:end-1], qns.indptr[2:end]))

"""
    pairs(qns::AbelianNumbers, choice::Symbol)
    pairs(qns::AbelianNumbers, ::Val{:indptr})
    pairs(qns::AbelianNumbers, ::Val{:counts})

Iterate over the `AbelianNumber=>slice` or `AbelianNumber=>count` pairs.
"""
@inline Base.pairs(qns::AbelianNumbers, choice::Symbol) = pairs(qns, choice|>Val)
@inline Base.pairs(qns::AbelianNumbers, choice::Union{Val{:indptr}, Val{:counts}}) = Base.Generator(=>, keys(qns), values(qns, choice))

"""
    OrderedDict(qns::AbelianNumbers, choice::Symbol) -> OrderedDict
    OrderedDict(qns::AbelianNumbers, ::Val{:indptr}) -> OrderedDict{qns|>eltype, UnitRange{Int}}
    OrderedDict(qns::AbelianNumbers, ::Val{:counts}) -> OrderedDict{qns|>eltype, Int}

Convert an `AbelianNumbers` to an ordered dict.
"""
@inline OrderedDict(qns::AbelianNumbers, choice::Symbol) = OrderedDict(qns, choice|>Val)
function OrderedDict(qns::AbelianNumbers, ::Val{:indptr})
    @assert qns.form≠'G' "OrderedDict error: input `AbelianNumbers` cannot be `G` formed."
    result = OrderedDict{qns|>eltype, UnitRange{Int}}()
    @inbounds for i = 1:length(qns)
        result[qns.contents[i]] = qns.indptr[i]+1:qns.indptr[i+1]
    end
    return result
end
function OrderedDict(qns::AbelianNumbers, ::Val{:counts})
    @assert qns.form≠'G' "OrderedDict error: input `AbelianNumbers` cannot be `G` formed."
    result = OrderedDict{qns|>eltype, Int}()
    @inbounds for i = 1:length(qns)
        result[qns.contents[i]] = qns.indptr[i+1] - qns.indptr[i]
    end
    return result
end

"""
    sort(qns::AbelianNumbers) -> Tuple{AbelianNumbers, Vector{Int}}

Sort the quantum numbers of an `AbelianNumbers`, return the sorted `AbelianNumbers` and the permutation array that sorts the expansion of the original `AbelianNumbers`.
"""
function Base.sort(qns::AbelianNumbers)
    ctpts = sortperm(qns.contents, alg=Base.Sort.QuickSort)
    masks = Vector{Bool}(undef, length(qns.contents))
    masks[1] = true
    nonduplicate = 1
    @inbounds for i = 2:length(ctpts)
        masks[i] = qns.contents[ctpts[i]]≠qns.contents[ctpts[i-1]]
        masks[i] && (nonduplicate += 1)
    end
    contents = Vector{qns|>eltype}(undef, nonduplicate)
    indptr = zeros(Int, nonduplicate+1)
    permutation = Vector{Int}(undef, dimension(qns))
    qncount, ptcount = 0, 0
    @inbounds for (mask, index) in zip(masks, ctpts)
        mask && (qncount += 1; contents[qncount] = qns.contents[index]; indptr[qncount+1] = indptr[qncount])
        indptr[qncount+1] += qns.indptr[index+1] - qns.indptr[index]
        for (i, p) in enumerate(qns.indptr[index]+1:qns.indptr[index+1])
            permutation[ptcount+i] = p
        end
        ptcount += qns.indptr[index+1] - qns.indptr[index]
    end
    return AbelianNumbers('C', contents, indptr, :indptr), permutation
end

"""
    findall(target::Union{QN, Tuple{Vararg{QN}}}, qns::AbelianNumbers{QN}, choice::Symbol) where {QN<:AbelianNumber} -> Vector{Int}
    findall(target::Union{QN, Tuple{Vararg{QN}}}, qns::AbelianNumbers{QN}, ::Val{:compression}) where {QN<:AbelianNumber} -> Vector{Int})
    findall(target::Union{QN, Tuple{Vararg{QN}}}, qns::AbelianNumbers{QN}, ::Val{:expansion}) where {QN<:AbelianNumber} -> Vector{Int}

Find all the indices of the target quantum numbers in the contents or the expansion of an `AbelianNumbers`.
"""
@inline Base.findall(target::Union{QN, Tuple{Vararg{QN}}}, qns::AbelianNumbers{QN}, choice::Symbol) where {QN<:AbelianNumber} = findall(target, qns, choice|>Val)
@inline Base.findall(target::QN, qns::AbelianNumbers{QN}, ::Val{:compression}) where {QN<:AbelianNumber} = findall((target,), qns, :compression|>Val)
@inline Base.findall(target::QN, qns::AbelianNumbers{QN}, ::Val{:expansion}) where {QN<:AbelianNumber} = findall((target,), qns, :expansion|>Val)
function Base.findall(targets::Tuple{Vararg{QN}}, qns::AbelianNumbers{QN}, ::Val{:compression}) where {QN<:AbelianNumber}
    result = Int[]
    if issorted(qns)
        for qn in targets
            range = searchsorted(qns.contents, qn)
            range.start<=range.stop && push!(result, range.start)
        end
    else
        for (i, qn) in enumerate(qns)
            qn∈targets && push!(result, i)
        end
    end
    return result
end
function Base.findall(targets::Tuple{Vararg{QN}}, qns::AbelianNumbers{QN}, ::Val{:expansion}) where {QN<:AbelianNumber}
    result = Int[]
    @inbounds for index in findall(targets, qns, :compression)
        append!(result, qns.indptr[index]+1:qns.indptr[index+1])
    end
    return result
end

"""
    filter(target::QN, qns::AbelianNumbers{QN}) where {QN<:AbelianNumber} -> AbelianNumbers{QN}
    filter(targets::Tuple{Vararg{QN}}, qns::AbelianNumbers{QN}) where {QN<:AbelianNumber} -> AbelianNumbers{QN}

Find a subset of an `AbelianNumbers` by picking out the target quantum numbers.
"""
@inline Base.filter(target::QN, qns::AbelianNumbers{QN}) where {QN<:AbelianNumber} = qns[findall(target, qns, :compression)]
@inline Base.filter(targets::Tuple{Vararg{QN}}, qns::AbelianNumbers{QN}) where {QN<:AbelianNumber} = qns[findall(targets, qns, :compression)]

"""
    permute(qns::AbelianNumbers, permutation::Vector{Int}, choice::Symbol) -> AbelianNumbers
    permute(qns::AbelianNumbers, permutation::Vector{Int}, ::Val{:compression}) -> AbelianNumbers
    permute(qns::AbelianNumbers, permutation::Vector{Int}, ::Val{:expansion}) -> AbelianNumbers

Reorder the quantum numbers contained in an `AbelianNumbers` with a permutation and return the new one.

For `:compression` case, the permutation is for the compressed contents of the original `AbelianNumbers` while for `:expansion` case, the permutation is for the expanded contents of the original `AbelianNumbers`.
"""
@inline permute(qns::AbelianNumbers, permutation::Vector{Int}, choice::Symbol) = permute(qns, permutation, choice|>Val)
@inline permute(qns::AbelianNumbers, permutation::Vector{Int}, ::Val{:compression}) = qns[permutation]
function permute(qns::AbelianNumbers, permutation::Vector{Int}, ::Val{:expansion})
    contents = Vector{qns|>eltype}(undef, length(permutation))
    indptr = zeros(Int, length(permutation)+1)
    expansion = expand(qns, :indices)
    @inbounds for (i, p) in enumerate(permutation)
        contents[i] = qns.contents[expansion[p]]
        indptr[i+1] = indptr[i] + 1
    end
    return AbelianNumbers('G', contents, indptr, :indptr)
end

"""
    expand(qns::AbelianNumbers, choice::Symbol)
    expand(qns::AbelianNumbers, ::Val{:contents}) -> Vector{eltype(qns)}
    expand(qns::AbelianNumbers, ::Val{:indices}) -> Vector{Int}

Expand the contents or indices of an `AbelianNumbers` to the uncompressed form.
"""
@inline expand(qns::AbelianNumbers, choice::Symbol) = expand(qns, choice|>Val)
function expand(qns::AbelianNumbers, ::Val{:contents})
    result = Vector{qns|>eltype}(undef, dimension(qns))
    @inbounds for i = 1:length(qns)
        for j = qns.indptr[i]+1:qns.indptr[i+1]
            result[j] = qns.contents[i]
        end
    end
    return result
end
function expand(qns::AbelianNumbers, ::Val{:indices})
    result = Vector{Int}(undef, dimension(qns))
    @inbounds for i = 1:length(qns)
        for j = qns.indptr[i]+1:qns.indptr[i+1]
            result[j] = i
        end
    end
    return result
end

"""
    +(qns::AbelianNumbers) -> AbelianNumbers
    +(qn::QN, qns::AbelianNumbers{QN}) where {QN<:AbelianNumber} -> AbelianNumbers{QN}
    +(qns::AbelianNumbers{QN}, qn::QN) where {QN<:AbelianNumber} -> AbelianNumbers{QN}

Overloaded `+` operator for `AbelianNumber`/`AbelianNumbers`.
!!! note
    1. The addition between an `AbelianNumbers` and an `AbelianNumber` is just a global shift of the contents of the `AbelianNumbers` by the `AbelianNumber`, therefore, the result is an `AbelianNumbers`.
    2. `+` cannot be used between two `AbelianNumbers` because the result is ambiguous. Instead, use `⊕` for direct sum and `⊗` for direct product.
    3. To ensure type stability, an `AbelianNumber` and an `AbelianNumbers` can be added together if and only if the former's type is the same with the latter's eltype.
"""
@inline Base.:+(qns::AbelianNumbers) = qns
@inline Base.:+(qn::QN, qns::AbelianNumbers{QN}) where {QN<:AbelianNumber} = qns + qn
@inline Base.:+(qns::AbelianNumbers{QN}, qn::QN) where {QN<:AbelianNumber} = AbelianNumbers(qns.form=='G' ? 'G' : 'U', [iqn+qn for iqn in qns], qns.indptr, :indptr)

"""
    -(qns::AbelianNumbers) -> AbelianNumbers
    -(qn::QN, qns::AbelianNumbers{QN}) where {QN<:AbelianNumber} -> AbelianNumbers{QN}
    -(qns::AbelianNumbers{QN}, qn::QN) where {QN<:AbelianNumber} -> AbelianNumbers{QN}

Overloaded `-` operator for `AbelianNumber`/`AbelianNumbers`.
!!! note
    1. The subtraction between an `AbelianNumbers` and an `AbelianNumber` is just a global shift of the contents of the `AbelianNumbers` by the `AbelianNumber`, therefore, the result is an `AbelianNumbers`.
    2. `-` cannot be used between two `AbelianNumbers` because the result is ambiguous. Instead, use `⊕` with signs for direct sum and `⊗` with signs for direct product.
    3. To ensure type stability, an `AbelianNumber` can be subtracted by an `AbelianNumbers` or vice versa if and only if the former's type is the same with the latter's eltype.
"""
@inline Base.:-(qns::AbelianNumbers) = AbelianNumbers(qns.form=='G' ? 'G' : 'U', -qns.contents, qns.indptr, :indptr)
@inline Base.:-(qn::QN, qns::AbelianNumbers{QN}) where {QN<:AbelianNumber} = AbelianNumbers(qns.form=='G' ? 'G' : 'U', [qn-iqn for iqn in qns], qns.indptr, :indptr)
@inline Base.:-(qns::AbelianNumbers{QN}, qn::QN) where {QN<:AbelianNumber} = AbelianNumbers(qns.form=='G' ? 'G' : 'U', [iqn-qn for iqn in qns], qns.indptr, :indptr)

"""
    *(factor::Integer, qns::AbelianNumbers) -> AbelianNumbers
    *(qns::AbelianNumbers, factor::Integer) -> AbelianNumbers

Overloaded `*` operator for the multiplication between an integer and an `AbelianNumbers`.
"""
@inline Base.:*(factor::Integer, qns::AbelianNumbers) = qns * factor
@inline Base.:*(qns::AbelianNumbers, factor::Integer) = AbelianNumbers(qns.form=='G' ? 'G' : 'U', [qn*factor for qn in qns.contents], qns.indptr, :indptr)

"""
    ^(qns::AbelianNumbers, factor::Integer) -> AbelianNumbers

Overloaded `^` operator for `AbelianNumbers`.

This operation translates into the direct product of `factor` copies of `qns`.
"""
@inline Base.:^(qns::AbelianNumbers, factor::Integer) = kron(NTuple{factor, qns|>typeof}(qns for i = 1:factor)...)

"""
    ⊕(qns::AbelianNumber...) -> AbelianNumbers{qns|>eltype}
    ⊕(qnses::AbelianNumbers...) -> qnses|>eltype

Get the direct sum of some `AbelianNumber`s or `AbelianNumbers`es.
"""
@inline ⊕(qns::AbelianNumber...) = union(qns...)
@inline ⊕(qnses::AbelianNumbers...) = union(qnses...)

"""
    ⊗(qnses::AbelianNumbers...) -> eltype(qnses)

Get the direct product of some `AbelianNumbers`es.
"""
@inline ⊗(qnses::AbelianNumbers...) = kron(qnses...)

"""
    union(qns::AbelianNumber...; signs=positives(qns)) -> AbelianNumbers
    union(qnses::AbelianNumbers{QN}...; signs=positives(qnses)) where {QN<:AbelianNumber} -> AbelianNumbers{QN}

Get the direct sum of some `AbelianNumber`s or `AbelianNumbers`es.
!!! note
    1. Physically, the direct sum of a couple of `AbelianNumber`s or `AbelianNumbers`es is defined by the direct sum of the bases of the Hilbert spaces they represent. Therefore, the dimension of the result equals the summation of those of the inputs. As a consequence, even for `AbelianNumber`s, the result will be an `AbelianNumbers` because the dimension of the result is greater than 1.
    2. Signs of `AbelianNumber`s or `AbelianNumbers`es can be provided when getting their direct sums.
"""
@inline function Base.union(qns::AbelianNumber...; signs=positives(qns))
    contents = [sign==1 ? qn : -qn for (sign, qn) in zip(signs, qns)]
    indptr = collect(0:length(qns))
    return AbelianNumbers('G', contents, indptr, :indptr)
end
function Base.union(qnses::AbelianNumbers{QN}...; signs=positives(qnses)) where {QN<:AbelianNumber}
    lengths = NTuple{fieldcount(typeof(qnses)), Int}(length(qns) for qns in qnses)
    contents = Vector{QN}(undef, sum(lengths))
    indptr = zeros(Int, sum(lengths)+1)
    count = 0
    @inbounds for (sign, qns, length) in zip(signs, qnses, lengths)
        for i = 1:length
            contents[count+i] = (sign == 1) ? qns.contents[i] : -qns.contents[i]
            indptr[count+i+1] = indptr[count+i] + qns.indptr[i+1] - qns.indptr[i]
        end
        count += length
    end
    return AbelianNumbers('G', contents, indptr, :indptr)
end

"""
    kron(qnses::AbelianNumbers{QN}...; signs=positives(qnses)) where {QN<:AbelianNumber} -> AbelianNumbers{QN}

Get the direct product of some `AbelianNumbers`es.
!!! note
    Physically, the direct product of a couple of `AbelianNumbers`es are defined by the direct product of the bases of the Hilbert spaces they represent. Therefore, the dimension of the result equals the product of those of the inputs. Meanwhile, each quantum number in the contents of the result is obtained by a summation of the corresponding quantum numbers of the inputs with the correct signs. This is a direct consequence of the Abelian nature of our quantum numbers.
"""
function Base.kron(qnses::AbelianNumbers{QN}...; signs=positives(qnses)) where {QN<:AbelianNumber}
    N = fieldcount(typeof(qnses))
    lengths = NTuple{N, Int}(i<N ? dimension(qns) : length(qns) for (i, qns) in enumerate(qnses))
    contents = Vector{QN}(undef, prod(lengths))
    indptr = zeros(Int, prod(lengths)+1)
    cache = Vector{QN}(undef, N)
    @inbounds expansions = NTuple{N-1, Vector{Int}}(expand(qnses[i], :indices) for i = 1:N-1)
    @inbounds for (i, indices) in enumerate(product(NTuple{N, UnitRange{Int64}}(1:length for length in reverse(lengths))...))
        for (j, (sign, qns, index)) in enumerate(zip(signs, qnses, reverse(indices)))
            pos = j<N ? expansions[j][index] : index
            cache[j] = sign==1 ? qns.contents[pos] : -qns.contents[pos]
        end
        indptr[i+1] = indptr[i] + qnses[end].indptr[indices[1]+1] - qnses[end].indptr[indices[1]]
        contents[i] = sum(cache)
    end
    return AbelianNumbers('G', contents, indptr, :indptr)
end

"""
    prod(qnses::AbelianNumbers{QN}...; signs=positives(qnses)) where {QN<:AbelianNumber} -> AbelianNumbers{QN}, Dict{QN, Dict{NTuple{length(qnses), QN}, UnitRange{Int}}}

Unitary Kronecker product of several `AbelianNumbers`es. The product result as well as the records of the product will be returned.
!!! note
    1. All input `AbelianNumbers` must be 'U' formed or 'C' formed.
    2. Since duplicate quantum number are not allowed in 'U' formed and 'C' formed `AbelianNumbers`es, in general, there exists a merge process of duplicate quantum numbers in the result. Therefore, records are needed to keep track of this process, which will be returned along with the product result. The records are stored in a `Dict{QN, Dict{NTuple{NTuple{length(qnses), QN}, UnitRange{Int}}}` typed dict, in which, for each nonduplicate quantum number `qn` in the result, there exist a record `Dict((qn₁, qn₂, ...)=>start:stop, ...)` telling what quantum numbers `(qn₁, qn₂, ...)` a merged duplicate `qn` comes from and what slice `start:stop` this merged duplicate corresponds in the result.
"""
function Base.prod(qnses::AbelianNumbers{QN}...; signs=positives(qnses)) where {QN<:AbelianNumber}
    @assert all((qns.form == 'U') || (qns.form == 'C') for qns in qnses) "prod error: all input qnses should be 'U' formed or 'C' formed."
    N = fieldcount(typeof(qnses))
    lengths = NTuple{N, Int}(length(qns) for qns in qnses)
    cache = Vector{QN}(undef, N)
    container = OrderedDict{QN, Int}()
    records = Dict{QN, Dict{NTuple{N, QN}, UnitRange{Int}}}()
    for indices in product(NTuple{N, UnitRange{Int64}}(1:length for length in reverse(lengths))...)
        qn, count = QN|>zero, 1
        @inbounds for (j, (sign, qns, index)) in enumerate(zip(signs, qnses, reverse(indices)))
            cache[j] = qns.contents[index]
            qn = sign==1 ? (qn+cache[j]) : (qn-cache[j])
            count = count * (qns.indptr[index+1]-qns.indptr[index])
        end
        container[qn] = get(container, qn, 0) + count
        !haskey(records, qn) && (records[qn] = Dict{NTuple{N, QN}, UnitRange{Int}}())
        records[qn][NTuple{N, QN}(cache)] = (container[qn]-count+1):container[qn]
    end
    contents = QN[qn for qn in keys(container)]
    sort!(contents)
    indptr = zeros(Int, length(container)+1)
    @inbounds for (i, qn) in enumerate(contents)
        indptr[i+1] = indptr[i] + container[qn]
    end
    return AbelianNumbers('C', contents, indptr, :indptr), records
end

"""
    decompose(target::QN, qnses::AbelianNumbers{QN}...; signs=positives(qnses), method=:montecarlo, nmax=20) where {QN<:AbelianNumber} -> Vector{NTuple{length(qnses), Int}}

Find a couple of decompositions of `target` with respect to `qnses`.
!!! note
    A tuple of integers `(i₁, i₂, ...)` is called a decomposition of a given `target` with respect to the given `qnses` if and only if they satisfy the "decomposition rule":
    ```math
    \\sum_\\text{j} \\text{signs}[\\text{j}]\\times\\text{qnses}[\\text{j}][\\text{i}_{\\text{j}}]==\\text{target}
    ```
    This equation is in fact a set of restricted [linear Diophantine equations](https://en.wikipedia.org/wiki/Diophantine_equation#Linear_Diophantine_equations). Indeed, our quantum numbers are always discrete Abelian ones and all instances of a concrete `AbelianNumber` forms a [module](https://en.wikipedia.org/wiki/Module_(mathematics)) over the [ring](https://en.wikipedia.org/wiki/Ring_(mathematics)) of integers. Therefore, each quantum number can be represented as a integral multiple of the unit element of the Abelian module, which results in the final reduction of the above equation to a set of linear Diophantine equations. Then finding a decomposition is equivalent to find a solution of the reduced linear Diophantine equations, with the restriction that the quantum numbers constructed from the solution should be in the corresponding `qnses`. Here we provide two methods to find such decompositions, one is by brute force, and the other is by Monte Carlo simulations.
"""
@inline decompose(target::QN, qnses::AbelianNumbers{QN}...; signs=positives(qnses), method=:montecarlo, nmax=20) where {QN<:AbelianNumber} = _decompose(Val(method), target, qnses...; signs=signs, nmax=nmax)
function _decompose(::Val{:bruteforce}, target::QN, qnses::AbelianNumbers{QN}...; signs, nmax) where {QN<:AbelianNumber}
    N = fieldcount(typeof(qnses))
    result = Set{NTuple{N, Int}}()
    cache = Vector{Int}(undef, N)
    dimensions = NTuple{N, Int}(dimension(qns) for qns in reverse(qnses))
    indices = findall(target, kron(qnses...; signs=signs), :expansion)
    nmax<length(indices) && (shuffle!(MersenneTwister(), indices); indices = @views indices[1:nmax])
    for index in indices
        @inbounds for (i, dimension) in enumerate(dimensions)
            cache[end+1-i] = (index-1)%dimension + 1
            index = (index-1)÷dimension + 1
        end
        push!(result, NTuple{N, Int}(cache))
    end
    return collect(NTuple{N, Int}, result)
end
function _decompose(::Val{:montecarlo}, target::QN, qnses::AbelianNumbers{QN}...; signs, nmax) where {QN<:AbelianNumber}
    seed!()
    N = fieldcount(typeof(qnses))
    result = Set{NTuple{N, Int}}()
    expansions = Vector{QN}[expand(qns, :contents) for qns in qnses]
    dimensions = NTuple{N, Int}(dimension(qns) for qns in qnses)
    diff = indices->norm(sum(sign==1 ? expansion[index] : -expansion[index] for (sign, expansion, index) in zip(signs, expansions, indices))-target)
    count = 1
    while true
        oldindices, newindices = Int[rand(1:dimension) for dimension in dimensions], ones(Int, N)
        olddiff, newdiff = diff(oldindices), diff(newindices)
        @inbounds while newdiff > 0
            pos = rand(1:N)
            index = rand(1:dimensions[pos]-1)
            newindices[pos] = index<newindices[pos] ? index : index+1
            newdiff = diff(newindices)
            if newdiff<=olddiff || exp(olddiff-newdiff)>rand()
                oldindices[:] = newindices
                olddiff = newdiff
            end
            newindices[:] = oldindices
        end
        count = count + 1
        push!(result, NTuple{N, Int}(newindices))
        (length(result)>=nmax || count>nmax*5) && break
    end
    return collect(NTuple{N, Int}, result)
end

"""
    SpinZ(Sz::Real)

The concrete `AbelianNumber` of a quantum system with spin z-component `Sz` conserved.
"""
@abeliannumber "SpinZ" Float64 (:Sz,) (Inf,)

"""
    ParticleNumber(N::Real)

The concrete `AbelianNumber` of a quantum system with particle number `N` conserved.
"""
@abeliannumber "ParticleNumber" Int (:N,) (Inf,)

"""
    SpinfulParticle(N::Real, Sz::Real)

The concrete `AbelianNumber` of a quantum system with both particle number `N` and spin z-component `Sz` conserved.
"""
@abeliannumber "SpinfulParticle" Float64 (:N, :Sz) (Inf, Inf)

"""
    spinzs(S::Real) -> AbelianNumbers{SpinZ}

Construct the `AbelianNumbers` of the Hilbert space of a single spin `S`.
"""
@inline spinzs(S::Real) = AbelianNumbers('C', [SpinZ(sz) for sz = -S:S], collect(0:Int(2*S+1)), :indptr)

"""
    particlenumbers(N::Real) -> AbelianNumbers{ParticleNumber}

Construct the `AbelianNumbers` of the Hilbert space of a single-particle state with at most `N` identical particles.
"""
@inline particlenumbers(N::Real) = AbelianNumbers('C', [ParticleNumber(np) for np = 0:N], collect(0:Int(N)+1), :indptr)

"""
    spinfulparticles(S::Real) -> AbelianNumbers{SpinfulParticle}

Construct the `AbelianNumbers` of the Hilbert space of a single site with internal degrees of freedom that can be ascribed to a spin `S`.
"""
function spinfulparticles(S::Real)
    contents = [SpinfulParticle(0.0, 0.0)]
    for n = 1:Int(2*S+1)
        for szs in Combinations{n}(-S:S)
            push!(contents, SpinfulParticle(n, sum(szs)))
        end
    end
    return sort(AbelianNumbers('G', contents, collect(0:length(contents)), :indptr))[1]
end

"""
    Momentum <: AbelianNumber{Int}

Abstract type for momentum.
"""
abstract type Momentum <: AbelianNumber{Int} end

"""
    Momentum₁{N}(k::Integer) where N

One dimensional momentum.
"""
struct Momentum₁{N} <: Momentum
    k::Int
    function Momentum₁{N}(k::Integer) where N
        @assert N>0 "Momentum₁ error: wrong period ($N)."
        remainder = k % N
        remainder < 0 && (remainder = remainder + N)
        new{N}(remainder)
    end
end
@inline periods(::Type{<:Momentum₁{N}}) where N = (N,)
@inline Int(m::Momentum₁) = m.k + 1

"""
    Momentum₂{N}(k₁::Integer, k₂::Integer) where N
    Momentum₂{N₁, N₂}(k₁::Integer, k₂::Integer) where {N₁, N₂}

Two dimensional momentum.
"""
struct Momentum₂{N₁, N₂} <: Momentum
    k₁::Int
    k₂::Int
    Momentum₂{N}(k₁::Integer, k₂::Integer) where N =  Momentum₂{N, N}(k₁, k₂)
    function Momentum₂{N₁, N₂}(k₁::Integer, k₂::Integer) where {N₁, N₂}
        @assert N₁>0 "Momentum₂ error: wrong 1st period ($N₁)."
        @assert N₂>0 "Momentum₂ error: wrong 2nd period ($N₂)."
        r₁ = k₁ % N₁
        r₂ = k₂ % N₂
        r₁ < 0 && (r₁ = r₁ + N₁)
        r₂ < 0 && (r₂ = r₂ + N₂)
        new{N₁, N₂}(r₁, r₂)
    end
end
@inline periods(::Type{<:Momentum₂{N₁, N₂}}) where {N₁, N₂} = (N₁, N₂)
@inline Int(m::Momentum₂{N₁, N₂}) where {N₁, N₂} = m.k₂ + m.k₁*N₂ + 1

"""
    Momentum₃{N}(k₁::Integer, k₂::Integer, k₃::Integer) where N
    Momentum₃{N₁, N₂, N₃}(k₁::Integer, k₂::Integer, k₃::Integer) where {N₁, N₂, N₃}

Three dimensional momentum.
"""
struct Momentum₃{N₁, N₂, N₃} <: Momentum
    k₁::Int
    k₂::Int
    k₃::Int
    Momentum₃{N}(k₁::Integer, k₂::Integer, k₃::Integer) where N =  Momentum₃{N, N, N}(k₁, k₂, k₃)
    function Momentum₃{N₁, N₂, N₃}(k₁::Integer, k₂::Integer, k₃::Integer) where {N₁, N₂, N₃}
        @assert N₁>0 "Momentum₃ error: wrong 1st period ($N₁)."
        @assert N₂>0 "Momentum₃ error: wrong 2nd period ($N₂)."
        @assert N₃>0 "Momentum₃ error: wrong 3rd period ($N₃)."
        r₁ = k₁ % N₁
        r₂ = k₂ % N₂
        r₃ = k₃ % N₃
        r₁ < 0 && (r₁ = r₁ + N₁)
        r₂ < 0 && (r₂ = r₂ + N₂)
        r₃ < 0 && (r₃ = r₃ + N₃)
        new{N₁, N₂, N₃}(r₁, r₂, r₃)
    end
end
@inline periods(::Type{<:Momentum₃{N₁, N₂, N₃}}) where {N₁, N₂, N₃} = (N₁, N₂, N₃)
@inline Int(m::Momentum₃{N₁, N₂, N₃}) where {N₁, N₂, N₃} = (m.k₁*N₂+m.k₂)*N₃ + m.k₃ + 1

"""
    Momenta{P<:Momentum} <: VectorSpace{P}

The allowed set of momenta.
"""
struct Momenta{P<:Momentum} <: VectorSpace{P} end
@inline Momenta(::Type{P}) where {P<:Momentum} = Momenta{P}()
@inline VectorSpaceStyle(::Type{<:Momenta}) = VectorSpaceCartesian()
@inline shape(::Momenta{P}) where {P<:Momentum} = map(period->0:period-1, reverse(periods(P)))
@inline Base.CartesianIndex(m::P, ::Momenta{P}) where {P<:Momentum} = CartesianIndex(Int.(reverse(values(m))))
@inline (::Type{<:Momentum})(index::CartesianIndex, ::Momenta{P}) where {P<:Momentum} = P(reverse(index.I)...)
@inline Base.:(==)(ms₁::Momenta, ms₂::Momenta) = periods(eltype(ms₁))==periods(eltype(ms₂))
@inline Base.isequal(ms₁::Momenta, ms₂::Momenta) = isequal( periods(eltype(ms₁)), periods(eltype(ms₂)))

end #module
