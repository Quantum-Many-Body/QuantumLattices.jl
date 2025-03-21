module Toolkit

using Base: @propagate_inbounds
using DataStructures: OrderedDict
using Format: FormatSpec, pyfmt
using InteractiveUtils: subtypes
using Printf: @printf
using StaticArrays: SVector

import QuantumLattices: id, shape, value

# Utilities
export atol, rtol, Float
export DirectProductedIndices, DirectSummedIndices, Segment, concatenate, delta, subscript, superscript, tostr

# Combinatorics
export Combinatorics, Combinations, DuplicateCombinations, DuplicatePermutations, Permutations

# Traits
export SubTypeTree, commontype, dissolve, fulltype, rawtype
export hasparameter, isparameterbound, isparameterbounds, parametercount, parametername, parameternames, parameterorder, parameterpair, parameterpairs, parametertype, parametertypes, promoteparameters, reparameter 
export contentcount, contentname, contentnames, contentorder, contenttype, contenttypes, getcontent, hascontent
export efficientoperations

# Composite structures
export CompositeDict, CompositeVector

# Vector spaces
export VectorSpace, VectorSpaceDirectProducted, VectorSpaceDirectSummed, VectorSpaceEnumerative, VectorSpaceGeneral, VectorSpaceStyle, VectorSpaceZipped
export DirectProductedVectorSpace, DirectSummedVectorSpace, ZippedVectorSpace

# Utilities
"Absolute tolerance for float numbers."
const atol = 5 * 10^-14

"Relative tolerance for float numbers."
const rtol = √atol

"Default float type."
const Float = Float64

"""
    tostr(number, ::Integer=5) -> String
    tostr(number::Integer, n::Integer=5) -> String
    tostr(number::Rational, n::Integer=5) -> String
    tostr(number::AbstractFloat, n::Integer=5) -> String
    tostr(number::Complex, n::Integer=5) -> String

Convert a number to a string with at most `n` decimal places.
"""
@inline tostr(number, ::Integer=5) = repr(number)
@inline tostr(number::Integer, ::Integer=5) = string(number)
@inline tostr(number::Rational, ::Integer=5) = number.den==1 ? repr(number.num) : repr(number)
function tostr(number::AbstractFloat, n::Integer=5)
    if number == 0.0
        result = "0.0"
    elseif 10^-5 < abs(number) < 10^6
        result = rstrip(pyfmt(FormatSpec(".$(n)f"), number), '0')
        (result[end] == '.') && (result = result * '0')
    else
        result = pyfmt(FormatSpec(".$(n)e"), number)
        epos = findfirst(isequal('e'), result)
        temp = rstrip(result[1:epos-1], '0')
        result = (temp[end] == '.') ? (temp * "0" * result[epos:end]) : (temp * result[epos:end])
    end
    return result
end
function tostr(number::Complex, n::Integer=5)
    sreal = (real(number) == 0) ? "0" : tostr(real(number), n)
    simag = (imag(number) == 0) ? "0" : tostr(imag(number), n)
    result = ""
    (sreal == "0") || (result = result * sreal)
    (simag == "0") || (result = ((simag[1] == '-') ? (result * simag) : (length(result) == 0) ? simag : (result * "+" * simag)) * "im")
    (length(result) == 0) && (result = "0.0")
    return result
end

"""
    tostr(value::Symbol) -> String
    tostr(value::Colon) -> String
    tostr(value::Char) -> String

Convert a single value to string.
"""
@inline tostr(value::Symbol) = string(value)
@inline tostr(::Colon) = ":"
@inline tostr(value::Char) = repr(value)

"""
    superscript(i::Integer) -> String

Convert an integer to the superscript.
"""
@inline superscript(i::Integer) = join(d==1 ? '¹' : d==2 ? '²' : d==3 ? '³' : '⁰'+d for d in Iterators.reverse(digits(i)))

"""
    subscript(i::Integer) -> String

Convert an integer to the subscript.
"""
@inline subscript(i::Integer) = join('₀'+d for d in Iterators.reverse(digits(i)))

"""
    delta(i, j) -> Int

Kronecker delta function.
"""
@inline delta(i, j) = (i == j) ? 1 : 0

"""
    concatenate(ts::Tuple...) -> Tuple

Concatenate tuples.
"""
@inline @generated concatenate(ts::Tuple...) = Expr(:tuple, map(i->:(ts[$i]...), 1:length(ts))...)

"""
    id(od::OrderedDict, i) -> keytype(od)

Get the ith key of an `OrderedDict`.
"""
@inline id(od::OrderedDict, i) = od.keys[i]

"""
    value(od::OrderedDict, i) -> valtype(od)

Get the ith value of an `OrderedDict`.
"""
@inline value(od::OrderedDict, i) = od.vals[i]

"""
    searchsortedfirst(table, basis; by=identity, lt=isless, rev=false) -> Int

Use the binary search method to find the position of a basis in a sorted table so that the order is preserved if the basis in inserted in that position.
"""
function Base.searchsortedfirst(table, basis; by=identity, lt=isless, rev=false)
    lo, hi = firstindex(table)-1, lastindex(table)+1
    @inbounds while lo < hi-1
        m = (lo+hi) >>> 1
        if lt(by(table[m]), by(basis)) ≠ rev
            lo = m 
        else
            hi = m
        end
    end
    return hi
end

"""
    DirectSummedIndices{T<:Tuple{Vararg{OrdinalRange{Int, Int}}}} <: AbstractVector{CartesianIndex{3}}

Direct sum of indexes.
"""
struct DirectSummedIndices{T<:Tuple{Vararg{OrdinalRange{Int, Int}}}} <: AbstractVector{CartesianIndex{3}}
    indices::T
end
@inline Base.size(indexes::DirectSummedIndices) = (sum(map(length, indexes.indices)),)
@inline Base.iterate(indexes::DirectSummedIndices) = CartesianIndex(1, first(first(indexes.indices)), 0), (1, first(first(indexes.indices)), 0)
function Base.iterate(indexes::DirectSummedIndices, state::NTuple{3, Int})
    m, n, len = state[1], state[2]+1, state[3]
    if n>last(indexes.indices[m])
        m += 1
        m>length(indexes.indices) && return
        n = first(indexes.indices[m])
        len += length(indexes.indices[m-1])
    end
    return CartesianIndex(m, n, len), (m, n, len)
end
function Base.getindex(indexes::DirectSummedIndices, i::Integer)
    count = 0
    @inbounds for (m, index) in enumerate(indexes.indices)
        len = length(index)
        if i<=len
            return CartesianIndex(m, index[i], count)
        else
            i -= len
            count += len
        end
    end
    error("getindex error: out of range.")
end

"""
    DirectSummedIndices(indexes::Tuple{Vararg{Union{<:Integer, OrdinalRange{<:Integer, <:Integer}}}})

Construct a `DirectSummedIndices`.
"""
@inline DirectSummedIndices(indexes::Tuple{Vararg{Union{<:Integer, OrdinalRange{<:Integer, <:Integer}}}}) = DirectSummedIndices(map(index->convert2ind(index), indexes))
@inline convert2ind(index::Integer) = Base.OneTo(index)
@inline convert2ind(index::AbstractUnitRange{<:Integer}) = first(index):last(index)
@inline convert2ind(index::OrdinalRange{<:Integer, <:Integer}) = first(index):step(index):last(index)

"""
    DirectProductedIndices{Order, N, T<:NTuple{N, OrdinalRange{Int, Int}}} <: AbstractVector{CartesianIndex{N}}

Direct product of indexes.

Here, the type parameter `Order` must be either `:forward` or `:backward`:
1) `:forward`: the direct product iterates over the first sub-index first like a Julia array;
2) `:backward`: the direct product iterates over the last sub-index first like a C/C++ array.
"""
struct DirectProductedIndices{Order, N, T<:NTuple{N, OrdinalRange{Int, Int}}} <: AbstractVector{CartesianIndex{N}}
    indices::T
    function DirectProductedIndices{Order}(indices::NTuple{N, OrdinalRange{Int, Int}}) where {Order, N}
        @assert Order∈(:forward, :backward) "DirectProductedIndices error: `Order` must be either `:forward` or `:backward`."
        new{Order, N, typeof(indices)}(indices)
    end
end
@inline Base.size(indexes::DirectProductedIndices) = (mapreduce(length, *, indexes.indices),)
@inline Base.getindex(indexes::DirectProductedIndices{:forward}, i::Integer) = CartesianIndices(indexes.indices)[i]
@inline Base.getindex(indexes::DirectProductedIndices{:backward}, i::Integer) = CartesianIndex(reverse(CartesianIndices(reverse(indexes.indices))[i].I)...)
@inline function Base.searchsortedfirst(indexes::DirectProductedIndices{:forward, N}, index::CartesianIndex{N}) where N
    return LinearIndices(indexes.indices)[index-first(indexes)+oneunit(typeof(index))]
end
@inline function Base.searchsortedfirst(indexes::DirectProductedIndices{:backward, N}, index::CartesianIndex{N}) where N
    return LinearIndices(reverse(indexes.indices))[CartesianIndex(reverse((index-first(indexes)+oneunit(typeof(index))).I)...)]
end
@inline Base.in(index::CartesianIndex, indexes::DirectProductedIndices) = index in CartesianIndices(indexes.indices)

"""
    DirectProductedIndices{Order}(indexes::Tuple{Vararg{Union{<:Integer, OrdinalRange{<:Integer, <:Integer}}}}) where Order

Construct a `DirectProductedIndices`.
"""
@inline function DirectProductedIndices{Order}(indexes::Tuple{Vararg{Union{<:Integer, OrdinalRange{<:Integer, <:Integer}}}}) where Order
    return DirectProductedIndices{Order}(map(index->convert2ind(index), indexes))
end

"""
    Segment{S} <: AbstractVector{S}

A segment.
"""
struct Segment{S} <: AbstractVector{S}
    start::S
    stop::S
    length::Int
    ends::Tuple{Bool, Bool}
end
@inline Base.:(==)(s₁::Segment, s₂::Segment) = ==(efficientoperations, s₁, s₂)
@inline Base.isequal(s₁::Segment, s₂::Segment) = isequal(efficientoperations, s₁, s₂)
@inline Base.size(segment::Segment) = (segment.length,)
function Base.getindex(segment::Segment, i::Integer)
    length = segment.length + count(isequal(false), segment.ends) - 1
    step = convert(eltype(segment), (segment.stop-segment.start)/length)
    start = segment.ends[1] ? segment.start : segment.start+step
    return start+(i-1)*step
end
@inline Base.getindex(segment::Segment, range::OrdinalRange{<:Integer}) = Segment(segment[first(range)], segment[last(range)], length(range), (true, true))
function Base.iterate(segment::Segment)
    segment.length==0 && return
    length = segment.length + count(isequal(false), segment.ends) - 1
    step = convert(eltype(segment), (segment.stop-segment.start)/length)
    start = segment.ends[1] ? segment.start : segment.start+step
    return start, (1, start, step)
end
function Base.iterate(segment::Segment, state)
    i, middle, step = state
    i==segment.length && return
    middle = middle + step
    return middle, (i+1, middle, step)
end
function Base.show(io::IO, segment::Segment{<:Number})
    left, right = (segment.ends[1] ? "[" : "("), (segment.ends[2] ? "]" : ")")
    @printf io "%s%s, %s%s" left segment.start segment.stop right
end
function Base.show(io::IO, segment::Segment)
    left, right = (segment.ends[1] ? "[" : "("), (segment.ends[2] ? "]" : ")")
    @printf io "%sp₁, p₂%s with p₁=%s and p₂=%s" left right segment.start segment.stop
end

"""
    Segment(start::Number, stop::Number, length::Integer; ends::Tuple{Bool, Bool}=(true, false))
    Segment(start::AbstractVector{<:Number}, stop::AbstractVector{<:Number}, length::Integer; ends::Tuple{Bool, Bool}=(true, false))
    Segment(start::NTuple{N, Number}, stop::NTuple{N, Number}, length::Integer; ends::Tuple{Bool, Bool}=(true, false)) where N

Construct a segment.
"""
function Segment(start::Number, stop::Number, length::Integer; ends::Tuple{Bool, Bool}=(true, false))
    @assert length>=0 "Segment error: length must be non-negative."
    dtype = promote_type(typeof(start), typeof(stop), Float)
    return Segment(convert(dtype, start), convert(dtype, stop), length, ends)
end
function Segment(start::AbstractVector{<:Number}, stop::AbstractVector{<:Number}, length::Integer; ends::Tuple{Bool, Bool}=(true, false))
    @assert length>=0 "Segment error: length must be non-negative."
    @assert Base.length(start)==Base.length(stop) "Segment error: start and stop should have equal length."
    dtype = SVector{Base.length(start), promote_type(eltype(start), eltype(stop), Float)}
    return Segment(convert(dtype, start), convert(dtype, stop), length, ends)
end
function Segment(start::NTuple{N, Number}, stop::NTuple{N, Number}, length::Integer; ends::Tuple{Bool, Bool}=(true, false)) where N
    @assert length>=0 "Segment error: length must be non-negative."
    dtype = SVector{N, promote_type(reduce(promote_type, fieldtypes(typeof(start))), reduce(promote_type, fieldtypes(typeof(start))), Float)}
    return Segment(convert(dtype, start), convert(dtype, stop), length, ends)
end

# Combinatorics
"""
    Combinatorics{M, C}

Abstract combinatorial algorithms.
"""
abstract type Combinatorics{M, C} end
@inline Base.eltype(::Type{<:Combinatorics{M, C}}) where {M, C} = NTuple{M, eltype(C)}

"""
    Combinations{M}(contents::C) where {M, C}

Combinations of `M` elements from `contents`. Duplicates are not allowed.
"""
struct Combinations{M, C} <: Combinatorics{M, C}
    contents::C
    firstindex::Int
    lastindex::Int
    length::Int
    function Combinations{M}(contents::C) where {M, C}
        fi, li, len = firstindex(contents), lastindex(contents), length(contents)
        @assert M>=0 && len==li-fi+1 "Combinations error: not supported."
        new{M, C}(contents, fi, li, len)
    end
end
@inline Base.length(com::Combinations{M}) where M = binomial(com.length, M)
@propagate_inbounds function Base.iterate(com::Combinations{M}) where M
    M>com.length && return
    M==0 && return (), [com.lastindex+2]
    return ntuple(i->com.contents[com.firstindex+i-1], Val(M)), nextmstate!(collect(com.firstindex:com.firstindex+M-1), com.lastindex, M)
end 
@propagate_inbounds function Base.iterate(com::Combinations{M}, state) where M
    state[1]>com.lastindex-M+1 && return
    return ntuple(i->com.contents[state[i]], Val(M)), nextmstate!(state, com.lastindex, M)
end
@propagate_inbounds function nextmstate!(state::Vector{Int}, lastindex::Int, M::Int)
    for i = M:-1:1
        state[i] += 1
        state[i]>lastindex-(M-i) && continue
        for j = i+1:M
            state[j] = state[j-1] + 1
        end
        break
    end
    return state
end

"""
    DuplicateCombinations{M}(contents::C) where {M, C}

Combinations of `M` elements from `contents`. Duplicates are allowed.
"""
struct DuplicateCombinations{M, C} <: Combinatorics{M, C}
    contents::C
    firstindex::Int
    lastindex::Int
    length::Int
    function DuplicateCombinations{M}(contents::C) where {M, C}
        fi, li, len = firstindex(contents), lastindex(contents), length(contents)
        @assert M>=0 && li-fi+1==len "DuplicateCombinations error: not supported."
        new{M, C}(contents, fi, li, len)
    end
end
@inline Base.length(com::DuplicateCombinations{M}) where M = binomial(com.length+M-1, com.length-1)
@propagate_inbounds function Base.iterate(com::DuplicateCombinations{M}) where M
    M==0 && return (), [com.length+1]
    return ntuple(i->com.contents[com.firstindex], Val(M)), nextdmstate!(fill(com.firstindex, M), com.lastindex, M)
end
@propagate_inbounds function Base.iterate(com::DuplicateCombinations{M}, state) where M
    state[1]>com.lastindex && return
    return ntuple(i->com.contents[state[i]], Val(M)), nextdmstate!(state, com.lastindex, M)
end
@propagate_inbounds function nextdmstate!(state::Vector{Int}, lastindex::Int, M::Int)
    for i = M:-1:1
        state[i] += 1
        state[i]>lastindex && continue
        for j = i+1:M
            state[j] = state[j-1]
        end
        break
    end
    return state
end

"""
    Permutations{M}(contents::C) where {M, C}

Permutations of `M` elements from `contents`. Duplicates are not allowed.
"""
struct Permutations{M, C} <: Combinatorics{M, C}
    contents::C
    firstindex::Int
    lastindex::Int
    length::Int
    function Permutations{M}(contents::C) where {M, C}
        fi, li, len = firstindex(contents), lastindex(contents), length(contents)
        @assert M>=0 && li-fi+1==len "Permutations error: not supported."
        new{M, C}(contents, fi, li, len)
    end
end
@inline Base.length(perm::Permutations{M}) where M = M<=perm.length ? prod((perm.length-M+1):perm.length) : 0
@propagate_inbounds function Base.iterate(perm::Permutations{M}) where M
    perm.length<M && return
    M==0 && return (), [perm.lastindex+1]
    state = collect(perm.firstindex:perm.lastindex)
    return ntuple(i->perm.contents[state[i]], Val(M)), nextpstate!(state, perm.length, M)
end
@propagate_inbounds function Base.iterate(perm::Permutations{M}, state) where M
    perm.lastindex<state[1] && return
    return ntuple(i->perm.contents[state[i]], Val(M)), nextpstate!(state, perm.length, M)
end
@propagate_inbounds function nextpstate!(state::Vector{Int}, length::Int, M::Int)
    if M < length
        j = M + 1
        while j<=length && state[M]>=state[j]
            j += 1
        end
    end
    if M<length && j<=length
        state[M], state[j] = state[j], state[M]
    else
        M<length && reverse!(state, M+1)
        i = M - 1
        while i>=1 && state[i]>=state[i+1]
            i -= 1
        end
        if i > 0
            j = length
            while j>i && state[i]>=state[j]
                j -= 1
            end
            state[i], state[j] = state[j], state[i]
            reverse!(state, i+1)
        else
            state[1] = typemax(eltype(state))
        end
    end
    return state
end

"""
    DuplicatePermutations{M}(contents::C) where {M, C}

Permutations of `M` elements from `contents`. Duplicates are allowed.
"""
struct DuplicatePermutations{M, C} <: Combinatorics{M, C}
    contents::C
    firstindex::Int
    lastindex::Int
    length::Int
    function DuplicatePermutations{M}(contents::C) where {M, C}
        fi, li, len = firstindex(contents), lastindex(contents), length(contents)
        @assert M>=0 && li-fi+1==len "DuplicatePermutations error: not supported."
        new{M, C}(contents, fi, li, len)
    end
end
@inline Base.length(perm::DuplicatePermutations{M}) where M = perm.length^M
@propagate_inbounds function Base.iterate(perm::DuplicatePermutations{M}) where M
    indices = CartesianIndices(ntuple(i->perm.firstindex:perm.lastindex, Val(M)))
    current = iterate(indices)
    isnothing(current) && return
    index, state = reverse(current[1].I), current[2]
    return ntuple(i->perm.contents[index[i]], Val(M)), (indices, state)
end
@propagate_inbounds function Base.iterate(perm::DuplicatePermutations{M}, state) where M
    current = iterate(state[1], state[2])
    isnothing(current) && return
    indices, index, state = state[1], reverse(current[1].I), current[2]
    return ntuple(i->perm.contents[index[i]], Val(M)), (indices, state)
end

# Traits
## with type themselves
"""
    SubTypeTree(root::Type)

Construct the complete depth-first subtype tree of a root `Type`.
"""
struct SubTypeTree
    root::Type
end
@inline Base.eltype(subtypetree::SubTypeTree) = eltype(typeof(subtypetree))
@inline Base.eltype(::Type{SubTypeTree}) = Type
@inline Base.IteratorSize(::Type{<:SubTypeTree}) = Base.SizeUnknown()
@inline function Base.iterate(subtypetree::SubTypeTree, state=Type[subtypetree.root])
    length(state)==0 && return
    current = popfirst!(state)
    prepend!(state, subtypes(current))
    return current, state
end
function Base.show(io::IO, subtypetree::SubTypeTree)
    println(io, subtypetree.root)
    depth = get(io, :depth, 0)
    prefix = get(io, :prefix, "")
    middle, terminator, skip, dash = "├", "└", "│", "─"
    children = subtypes(subtypetree.root)
    for (i, child) in enumerate(children)
        print(io, prefix)
        extra = if i == length(children)
            print(io, terminator)
            " " ^ (textwidth(skip)+textwidth(dash)+1)
        else
            print(io, middle)
            skip * " "^(textwidth(dash)+1)
        end
        print(io, dash, " ")
        show(IOContext(io, :depth=>depth+1, :prefix=>prefix*extra), SubTypeTree(child))
    end
end

"""
    commontype(f::Function, types, ::Type{T}=Any) where T

Find the common return type of a function.
"""
function commontype(f::Function, types, ::Type{T}=Any) where T
    result = Union{}
    for rtype in Base.return_types(f, types)
        result = promote_type(result, commontype(rtype))
    end
    @assert result<:T "commontype error: wrong common type."
    return result
end
@inline commontype(t) = t
@inline commontype(t::Union) = promote_type(commontype(t.a), commontype(t.b))

"""
    DataType(T::DataType) -> DataType
    DataType(T::UnionAll) -> DataType

Get the DataType.
"""
@inline Base.DataType(T::DataType) = T
@inline Base.DataType(T::UnionAll) = DataType(T.body)

"""
    supertype(T, termination::Symbol) -> DataType

Get the supertype of `T` till termination.
"""
@inline Base.supertype(::Type{T}, termination::Symbol) where T = _supertype(T, Val(termination))
@inline @generated function _supertype(::Type{T}, ::Val{termination}) where {T, termination}
    result = T
    while nameof(result)≠termination && nameof(result)≠:Any
        result = supertype(result)
    end
    nameof(result)≠termination && error("_supertype error: termination is not the name of a valid supertype of $(nameof(T)).")
    return result
end

@inline _rawtype(T::DataType) = T.name.wrapper
@inline _rawtype(T::UnionAll) = _rawtype(T.body)
"""
    rawtype(::Type{T}) where T -> DataType/UnionAll

Get the "raw part" of a type. That is, the type without all its type parameters.
"""
@inline @generated rawtype(::Type{T}) where T = _rawtype(T)

"""
    map(f::Function, ::Type{T}) where {T<:Tuple} -> Type{Tuple}

Apply the `f` function on the fieldtypes of a type of tuple `T`.
"""
@inline @generated Base.map(f::Function, ::Type{T}) where {T<:Tuple} = Expr(:curly, Tuple, [:(f($type)) for type in fieldtypes(T)]...)

## with type parameters
"""
    parametercount(::Type{T}) where T -> Int

For a type `T`, get the number of its type parameters.
"""
@inline @generated parametercount(::Type{T}) where T = length(DataType(T).parameters)

"""
    parametername(::Type{T}, i::Integer) where T -> Symbol

For a type `T`, get the name of its ith type parameter.
"""
@inline parametername(::Type{T}, i::Integer) where T = _parametername(parameternames(T)|>Val, Val(i))
@inline @generated _parametername(::Val{names}, ::Val{i}) where {names, i} = QuoteNode(names[i])

"""
    parameterorder(::Type{T}, name::Symbol) where T -> Int

For a type `T`, get the order of one of its type parameters.
"""
@inline parameterorder(::Type{T}, name::Symbol) where T = _order(Val(name), parameternames(T)|>Val)
@inline @generated _order(::Val{name}, ::Val{names}) where {name, names} = findfirst(isequal(name), names)::Int

"""
    parametertype(::Type{T}, name::Symbol) where T
    parametertype(::Type{T}, i::Integer) where T

For a type `T`, get the type of one of its type parameters.
"""
@inline parametertype(::Type{T}, name::Symbol) where T = parametertype(T, parameterorder(T, name))
@inline parametertype(::Type{T}, i::Integer) where T = _parametertype(T, Val(i))
@inline @generated function _parametertype(::Type{T}, ::Val{i}) where {T, i}
    result = DataType(T).parameters[i]
    return isa(result, TypeVar) ? result.ub : isa(result, Symbol) ? QuoteNode(result) : result
end

"""
    parameterpair(::Type{T}, name::Symbol) where T
    parameterpair(::Type{T}, i::Integer) where T

For type `T`, get the name-type pair of one of its type parameters.

The result is stored in the type parameters of a `Pair`.
"""
@inline parameterpair(::Type{T}, name::Symbol) where T = Pair{name, parametertype(T, name)}
@inline parameterpair(::Type{T}, i::Integer) where T = Pair{parametername(T, i), parametertype(T, i)}

"""
    isparameterbound(::Type{T}, i::Integer, D) where T -> Bool
    isparameterbound(::Type{T}, name::Symbol, D) where T -> Bool

For a type `T`, judge whether a type `D` should be considered as the upper bound of one of its type parameters.

!!! note
    The default implementations of this function is
    ```julia
    isparameterbound(::Type{}, ::Val{}, ::Type{D}) where D = !isconcretetype(D)
    isparameterbound(::Type{}, ::Val{}, ::Any) = false
    ```
"""
@inline isparameterbound(::Type{T}, i::Integer, D) where T = isparameterbound(T, Val(i), D)
@inline isparameterbound(::Type{T}, name::Symbol, D) where T = isparameterbound(T, Val(name), D)
@inline isparameterbound(::Type{}, ::Val{}, ::Type{D}) where D = !isconcretetype(D)
@inline isparameterbound(::Type{}, ::Val{}, ::Any) = false

"""
    hasparameter(::Type{T}, name::Symbol) where T -> Bool

For type `T`, judge whether it has a type parameter specified by `name`.
"""
@inline hasparameter(::Type{T}, name::Symbol) where T = _hasparameter(Val(name), parameternames(T)|>Val)
@inline _hasparameter(::Val{name}, ::Val{names}) where {name, names} = name∈names

"""
    parameternames(::Type{T}) where T -> Tuple{Vararg{Symbol}}

For a type `T`, get the names of all its type parameters.
"""
@inline parameternames(::Type{T}) where T = error("parameternames error: not defined for $(nameof(T)).")

"""
    parametertypes(::Type{T}) where T

For a type `T`, get the types of all its type parameters.

The returned types are stored in the type parameters of a `Tuple`.
"""
@inline parametertypes(::Type{T}) where T = _parametertypes(T, parametercount(T)|>Val)
@inline @generated function _parametertypes(::Type{T}, ::Val{C}) where {T, C}
    exprs = []
    for i = 1:C
        push!(exprs, :(parametertype(T, $i)))
    end
    return Expr(:curly, :Tuple, exprs...)
end

"""
    parameterpairs(::Type{T}) where T

For a type `T`, get the name-type pairs of all its type parameters.

The return types are stored in the type parameters of a `NamedTuple`.
"""
@inline parameterpairs(::Type{T}) where T = NamedTuple{parameternames(T), parametertypes(T)}

"""
    isparameterbounds(::Type{T}, ::Type{PS}) where {T, PS<:Tuple} -> Tuple{Vararg{Bool}}
    isparameterbounds(::Type{T}, ::Type{PS}) where {T, PS<:NamedTuple} -> Tuple{Vararg{Bool}}

For a type `T`, judge whether the types specified by `PS` should be considered as the upper bounds of its corresponding type parameters.
"""
@inline isparameterbounds(::Type{T}, ::Type{PS}) where {T, PS<:Tuple} = _isparameterbounds(T, PS, parametercount(T)|>Val)
@inline isparameterbounds(::Type{T}, ::Type{PS}) where {T, PS<:NamedTuple} = _isparameterbounds(T, PS, parameternames(T)|>Val)
@inline @generated function _isparameterbounds(::Type{T}, ::Type{PS}, ::Val{C}) where {T, PS<:Tuple, C}
    exprs = []
    for i = 1:C
        N = Val(i)
        P = fieldtype(PS, i)
        push!(exprs, :(isparameterbound(T, $N, $P)))
    end
    return Expr(:tuple, exprs...)
end
@inline @generated function _isparameterbounds(::Type{T}, ::Type{PS}, ::Val{names}) where {T, PS<:NamedTuple, names}
    exprs = []
    for i = 1:length(names)
        N = Val(names[i])
        P = fieldtype(PS, names[i])
        push!(exprs, :(isparameterbound(T, $N, $P)))
    end
    return Expr(:tuple, exprs...)
end

"""
    reparameter(::Type{T}, i::Integer, P, ub::Bool=isparameterbound(T, i, P)) where T
    reparameter(::Type{T}, name::Symbol, P, ub::Bool=isparameterbound(T, name, P)) where T

For a type `T`, replace the type of its ith type parameter with `P`. Here, `ub` determines whether `P` should be considered as the upper bound. 
"""
@inline reparameter(::Type{T}, i::Integer, P, ub::Bool=isparameterbound(T, i, P)) where T = _reparameter(T, Val(i), Tuple{P}, Val(ub))
@inline reparameter(::Type{T}, name::Symbol, P, ub::Bool=isparameterbound(T, name, P)) where T = _reparameter(T, parameterorder(T, name)|>Val, Tuple{P}, Val(ub))
@inline @generated function _reparameter(::Type{T}, ::Val{i}, ::Type{P}, ::Val{ub}) where {T, i, P, ub}
    params = collect(DataType(T).parameters)
    params[i] = ub ? TypeVar(gensym(), fieldtype(P, 1)) : fieldtype(P, 1)
    V = Core.apply_type(rawtype(T), params...)
    for k = 1:length(params)
        isa(params[k], TypeVar) && (V = UnionAll(params[k], V))
    end
    return V
end

"""
    promoteparameters(::Type{T₁}, ::Type{T₂}) where {T₁<:NamedTuple, T₂<:NamedTuple}

Promote the types specified by two named tuples with the same names accordingly.

The result is stored in the type parameters of a `NamedTuple`.
"""
@inline @generated function promoteparameters(::Type{T₁}, ::Type{T₂}) where {T₁<:NamedTuple, T₂<:NamedTuple}
    exprs = []
    names = Tuple(union(fieldnames(T₁), fieldnames(T₂)))
    for name in names
        hasfield(T₁, name) && (F₁ = fieldtype(T₁, name))
        hasfield(T₂, name) && (F₂ = fieldtype(T₂, name))
        if hasfield(T₁, name) && hasfield(T₂, name)
            push!(exprs, :(promote_type($F₁, $F₂)))
        else
            push!(exprs, hasfield(T₁, name) ? F₁ : F₂)
        end
    end
    return Expr(:curly, :NamedTuple, names, Expr(:curly, :Tuple, exprs...))
end

"""
    fulltype(::Type{T}, ::Type{PS}, ubs::Tuple{Vararg{Bool}}=isparameterbounds(T, PS)) where {T, PS<:Tuple}
    fulltype(::Type{T}, ::Type{PS}, ubs::Tuple{Vararg{Bool}}=isparameterbounds(T, PS)) where {T, PS<:NamedTuple}

Get the full type of type `T` with the type parameters replaced by those of `PS`.

Here, `ubs` determines whether the new type parameter should be considered as the upper bound accordingly.
"""
@inline function fulltype(::Type{T}, ::Type{PS}, ubs::Tuple{Vararg{Bool}}=isparameterbounds(T, PS)) where {T, PS<:Tuple}
    @assert parametercount(T) == fieldcount(PS) == length(ubs) "fulltype error: length-mismatched input parameters."
    return _fulltype(T, PS, Val(ubs))
end
@inline function fulltype(::Type{T}, ::Type{PS}, ubs::Tuple{Vararg{Bool}}=isparameterbounds(T, PS)) where {T, PS<:NamedTuple}
    _fulltype(T, PS, parameternames(T)|>Val, Val(ubs))
end
@inline @generated function _fulltype(::Type{T}, ::Type{PS}, ::Val{ubs}) where {T, PS, ubs}
    VS = [ubs[i] ? TypeVar(gensym(), fieldtype(PS, i)) : fieldtype(PS, i) for i=1:length(ubs)]
    V = Core.apply_type(rawtype(T), VS...)
    for i = 1:length(VS)
        ubs[i] && (V = UnionAll(VS[i], V))
    end
    return V
end
@inline @generated function _fulltype(::Type{T}, ::Type{PS}, ::Val{names}, ::Val{ubs}) where {T, PS<:NamedTuple, names, ubs}
    VS = [(ubs[i] ? TypeVar(gensym(), fieldtype(PS, names[i])) : fieldtype(PS, names[i])) for i=1:length(names)]
    V = Core.apply_type(rawtype(T), VS...)
    for i = 1:length(VS)
        ubs[i] && (V = UnionAll(VS[i], V))
    end
    return V
end

## with type contents
"""
    contentcount(::Type{T}) where T -> Int

For a type `T`, get the number of its predefined contents.
"""
@inline contentcount(::Type{T}) where T = length(contentnames(T))

"""
    contentname(::Type{T}, i::Integer) where T -> Symbol

For a type `T`, get the name of its ith predefined content.
"""
@inline contentname(::Type{T}, i::Integer) where T = contentnames(T)[i]

"""
    contentorder(::Type{T}, name::Symbol) where T -> Int

For a type `T`, get the position order of a predefined content by the name.
"""
@inline contentorder(::Type{T}, name::Symbol) where T = _order(Val(name), contentnames(T)|>Val)

"""
    contenttype(::Type{T}, name::Symbol) where T
    contenttype(::Type{T}, ::Val{name}) where {T, name}

For a type `T`, get the type of a predefined content by the name.
"""
@inline contenttype(::Type{T}, name::Symbol) where T = contenttype(T, Val(name))
@inline contenttype(::Type{T}, ::Val{name}) where {T, name} = fieldtype(T, name)

"""
    hascontent(::Type{T}, name::Symbol) where T -> Bool

For a type `T`, judge whether it has a predefined content specified by `name`.
"""
@inline hascontent(::Type{T}, name::Symbol) where T = _hascontent(name|>Val, contentnames(T)|>Val)
@inline @generated _hascontent(::Val{name}, ::Val{names}) where {name, names} = name∈names

"""
    getcontent(m, i::Integer)
    getcontent(m, name::Symbol)
    getcontent(m, ::Val{name}) where name

Get the value of the predefined content of `m`. 
"""
@inline getcontent(m, i::Integer) = getcontent(m, contentname(typeof(m), i))
@inline getcontent(m, name::Symbol) = getcontent(m, Val(name))
@inline getcontent(m, ::Val{name}) where name = getfield(m, name)

"""
    contentnames(::Type{T}) where T -> Tuple{Vararg{Symbol}}

For a type `T`, define the names of its predefined contents.
"""
@inline @generated contentnames(::Type{T}) where T = fieldnames(T)

"""
    contenttypes(::Type{T}) where T

For a type `T`, get the types of its predefined contents.
"""
@inline contenttypes(::Type{T}) where T = _contenttypes(T, contentnames(T)|>Val)
@inline @generated function _contenttypes(::Type{T}, ::Val{names}) where {T, names}
    exprs = []
    for name in QuoteNode.(names)
        push!(exprs, :(contenttype(T, $name)))
    end
    return Expr(:curly, :Tuple, exprs...)
end

"""
    dissolve(m, f::Function=identity, args...; kwargs...) -> Tuple

Convert `m` to a tuple by the function `f` applied elementally to its contents with the extra positional arguments (`args`) and keyword arguments (`kwargs`). 

To each content of `m`, the underlying interface of the `dissolve` function when `f` is applied is as follows:
```julia
dissolve(m, Val(name), f, args...; kwargs...)
```
Here, `name` is the name of the corresponding content of `m`.

Basically, the rule of how `f` operates on each field of `m` can be overridden by redefining the above `dissolve` function.

!!! note
    The default `dissolve` function ignores the operation of function `f` and just return the content value of `m`.
"""
@inline dissolve(m, f::Function=identity, args...; kwargs...) = dissolvehelper(m, Val(contentnames(typeof(m))), f, args...; kwargs...)
@inline @generated function dissolvehelper(m, ::Val{names}, f::Function, args...; kwargs...) where names
    exprs = []
    for name in names
        name = Val(name)
        push!(exprs, :(dissolve(m, $name, f, args...; kwargs...)))
    end
    return Expr(:tuple, exprs...)
end

"""
    dissolve(m, ::Val{name}, f::Function, args...; kwargs...) where name

Dissolve the content specified by `name` of `m` by the function `f` applied with the extra positional arguments (`args`) and keyword arguments (`kwargs`).
"""
@inline dissolve(m, ::Val{name}, ::Function, args...; kwargs...) where name = getcontent(m, Val(name))

## with type operations
struct EfficientOperations end
"""
    efficientoperations

Indicate that the efficient operations, i.e. "=="/"isequal", "<"/"isless" or "replace", will be used.
"""
const efficientoperations = EfficientOperations()

"""
    ==(::EfficientOperations, o₁, o₂) -> Bool

Compare two objects and judge whether they are equivalent to each other.
"""
@inline @generated function Base.:(==)(::EfficientOperations, o₁, o₂)
    fcount = fieldcount(o₁)
    if fcount == fieldcount(o₂)
        if fcount == 0
            return true
        else
            expr = :(getfield(o₁, 1) == getfield(o₂, 1))
            for i = 2:fcount
                expr = Expr(:&&, expr, :(getfield(o₁, $i) == getfield(o₂, $i)))
            end
            return expr
        end
    else
        return false
    end
end

"""
    isequal(::EfficientOperations, o₁, o₂) -> Bool

Compare two objects and judge whether they are equivalent to each other.
"""
@inline @generated function Base.isequal(::EfficientOperations, o₁, o₂)
    fcount = fieldcount(o₁)
    if fcount == fieldcount(o₂)
        if fcount == 0
            return true
        else
            expr = :(isequal(getfield(o₁, 1), getfield(o₂, 1)))
            for i = 2:fcount
                expr = Expr(:&&, expr, :(isequal(getfield(o₁, $i), getfield(o₂, $i))))
            end
            return expr
        end
    else
        return false
    end
end

"""
    <(::EfficientOperations, o₁, o₂) -> Bool

Compare two objects and judge whether the first is less than the second.
"""
@inline @generated function Base.:<(::EfficientOperations, o₁, o₂)
    n₁, n₂ = fieldcount(o₁), fieldcount(o₂)
    N = min(n₁, n₂)
    expr = (n₁ < n₂) ? Expr(:if, :(getfield(o₁, $N) == getfield(o₂, $N)), true, false) : false
    expr = Expr(:if, :(getfield(o₁, $N) < getfield(o₂, $N)), true, expr)
    for i in range(N-1, stop=1, step=-1)
        expr = Expr(:if, :(getfield(o₁, $i) > getfield(o₂, $i)), false, expr)
        expr = Expr(:if, :(getfield(o₁, $i) < getfield(o₂, $i)), true, expr)
    end
    return expr
end

"""
    isless(::EfficientOperations, o₁, o₂) -> Bool

Compare two objects and judge whether the first is less than the second.
"""
@inline @generated function Base.isless(::EfficientOperations, o₁, o₂)
    n₁, n₂ = fieldcount(o₁), fieldcount(o₂)
    N = min(n₁, n₂)
    expr = n₁ < n₂ ? Expr(:if, :(getfield(o₁, $N) == getfield(o₂, $N)), true, false) : false
    expr = Expr(:if, :(isless(getfield(o₁, $N), getfield(o₂, $N))), true, expr)
    for i in range(N-1, stop=1, step=-1)
        expr = Expr(:if, :(isless(getfield(o₂, $i), getfield(o₁, $i))), false, expr)
        expr = Expr(:if, :(isless(getfield(o₁, $i), getfield(o₂, $i))), true, expr)
    end
    return expr
end

"""
    isapprox(::EfficientOperations, o₁, o₂; atol=atol, rtol=rtol) -> Bool
    isapprox(::EfficientOperations, fields::Union{Nothing, Union{Integer, Symbol}, Tuple{Vararg{Union{Integer, Symbol}}}}, o₁, o₂; atol=atol, rtol=rtol) -> Bool
    isapprox(::EfficientOperations, ::Val{fields}, o₁, o₂; atol=atol, rtol=rtol) where fields -> Bool

Compare two objects and judge whether they are inexactly equivalent to each other.
"""
@inline Base.isapprox(::EfficientOperations, o₁, o₂; atol=atol, rtol=rtol) = isapprox(efficientoperations, nothing, o₁, o₂, atol=atol, rtol=rtol)
@inline function Base.isapprox(::EfficientOperations, fields::Union{Nothing, Union{Integer, Symbol}, Tuple{Vararg{Union{Integer, Symbol}}}}, o₁, o₂; atol=atol, rtol=rtol)
    isapprox(efficientoperations, fields|>Val, o₁, o₂; atol=atol, rtol=rtol)
end
@inline @generated function Base.isapprox(::EfficientOperations, ::Val{fields}, o₁, o₂; atol=atol, rtol=rtol) where fields
    (fieldcount(o₁)≠fieldcount(o₂) || fieldnames(o₁)≠fieldnames(o₂)) && return false
    isnothing(fields) && (fields = ntuple(i->i, fieldcount(o₁)))
    isa(fields, Union{Integer, Symbol}) && (fields = (fields,))
    fields = Set(isa(field, Symbol) ? findfirst(isequal(field), fieldnames(o₁)) : field for field in fields)
    exprs = []
    for i = 1:fieldcount(o₁)
        if i∈fields
            push!(exprs, :(isapprox(getfield(o₁, $i), getfield(o₂, $i); atol=atol, rtol=rtol)::Bool))
        else
            push!(exprs, :(isequal(getfield(o₁, $i), getfield(o₂, $i))::Bool))
        end
    end
    return Expr(:&&, exprs...)
end

# Composite structures
"""
    CompositeVector{T}

A composite vector can be considered as a vector that is implemented by including a concrete subtype of `AbstractVector` as its data attribute.
"""
abstract type CompositeVector{T} <: AbstractVector{T} end
@inline contentnames(::Type{<:CompositeVector}) = (:contents,)
@inline dissolve(cv::CompositeVector, ::Val{:contents}, f::Function, args...; kwargs...) = f(getcontent(cv, :contents), args...; kwargs...)
@inline Base.axes(cv::CompositeVector) = axes(getcontent(cv, :contents))
@inline Base.size(cv::CompositeVector) = size(getcontent(cv, :contents))
@inline Base.:(==)(cv1::CompositeVector, cv2::CompositeVector) = ==(efficientoperations, cv1, cv2)
@inline Base.isequal(cv1::CompositeVector, cv2::CompositeVector) = isequal(efficientoperations, cv1, cv2)
@inline Base.getindex(cv::CompositeVector, i::Union{<:Integer, CartesianIndex}) = getcontent(cv, :contents)[i]
@inline Base.getindex(cv::CompositeVector, inds) = rawtype(typeof(cv))(dissolve(cv, getindex, inds)...)
@inline Base.setindex!(cv::CompositeVector, value, inds) = (getcontent(cv, :contents)[inds] = value)
@inline Base.push!(cv::CompositeVector, values...) = (push!(getcontent(cv, :contents), values...); cv)
@inline Base.pushfirst!(cv::CompositeVector, values...) = (pushfirst!(getcontent(cv, :contents), values...); cv)
@inline Base.insert!(cv::CompositeVector, index::Integer, value) = (insert!(getcontent(cv, :contents), index, value); cv)
@inline Base.append!(cv::CompositeVector, values) = (append!(getcontent(cv, :contents), values); cv)
@inline Base.prepend!(cv::CompositeVector, values) = (prepend!(getcontent(cv, :contents), values); cv)
@inline Base.splice!(cv::CompositeVector, index::Integer, replacement=Base._default_splice) = splice!(getcontent(cv, :contents), index, replacement)
@inline Base.splice!(cv::CompositeVector, range::UnitRange{<:Integer}, replacement=Base._default_splice) = rawtype(typeof(cv))(dissolve(cv, splice!, range, replacement)...)
@inline Base.deleteat!(cv::CompositeVector, indices) = (deleteat!(getcontent(cv, :contents), indices); cv)
@inline Base.pop!(cv::CompositeVector) = pop!(getcontent(cv, :contents))
@inline Base.popfirst!(cv::CompositeVector) = popfirst!(getcontent(cv, :contents))
@inline Base.empty!(cv::CompositeVector) = (empty!(getcontent(cv, :contents)); cv)
@inline Base.empty(cv::CompositeVector) = rawtype(typeof(cv))(dissolve(cv, empty)...)
@inline Base.reverse(cv::CompositeVector) =  rawtype(typeof(cv))(dissolve(cv, reverse)...)
@inline Base.similar(cv::CompositeVector, dtype::Type=eltype(cv), dims::Tuple{Vararg{Int}}=size(cv)) = rawtype(typeof(cv))(dissolve(cv, similar, dtype, dims)...)

"""
    CompositeDict{K, V}

A composite dict can be considered as a dict that is implemented by including a concrete subtype of `AbstractDict` as its data attribute.
"""
abstract type CompositeDict{K, V} <: AbstractDict{K, V} end
@inline contentnames(::Type{<:CompositeDict}) = (:contents,)
@inline dissolve(cd::CompositeDict, ::Val{:contents}, f::Function, args...; kwargs...) = f(getcontent(cd, :contents), args...; kwargs...)
@inline Base.isempty(cd::CompositeDict) = isempty(getcontent(cd, :contents))
@inline Base.length(cd::CompositeDict) = length(getcontent(cd, :contents))
@inline Base.haskey(cd::CompositeDict, key) = haskey(getcontent(cd, :contents), key)
@inline Base.in(p::Pair, cd::CompositeDict, valcmp=(==)) = in(p, getcontent(cd, :contents), valcmp)
@inline Base.:(==)(cd1::CompositeDict, cd2::CompositeDict) = ==(efficientoperations, cd1, cd2)
@inline Base.isequal(cd1::CompositeDict, cd2::CompositeDict) = isequal(efficientoperations, cd1, cd2)
@inline Base.get(cd::CompositeDict, key, default) = get(getcontent(cd, :contents), key, default)
@inline Base.get(f::Base.Callable, cd::CompositeDict, key) = get(f, getcontent(cd, :contents), key)
@inline Base.getkey(cd::CompositeDict, key, default) = getkey(getcontent(cd, :contents), key, default)
@inline Base.getindex(cd::CompositeDict{K, V}, index::K) where {K, V} = getcontent(cd, :contents)[index]
@inline Base.push!(cd::CompositeDict, ps::Pair...) = (push!(getcontent(cd, :contents), ps...); cd)
@inline Base.get!(cd::CompositeDict, key, default) = get!(getcontent(cd, :contents), key, default)
@inline Base.get!(default::Union{Function, Type}, cd::CompositeDict, key) = get!(default, getcontent(cd, :contents), key)
@inline Base.setindex!(cd::CompositeDict{K, V}, value::V, index::K) where {K, V} = (getcontent(cd, :contents)[index] = value)
@inline Base.pop!(cd::CompositeDict) = pop!(getcontent(cd, :contents))
@inline Base.pop!(cd::CompositeDict, key) = pop!(getcontent(cd, :contents), key)
@inline Base.pop!(cd::CompositeDict, key, default) = pop!(getcontent(cd, :contents), key, default)
@inline Base.delete!(cd::CompositeDict, key) = (delete!(getcontent(cd, :contents), key); cd)
@inline Base.empty!(cd::CompositeDict) = (empty!(getcontent(cd, :contents)); cd)
@inline Base.merge(cd::CD, others::CD...) where CD <: CompositeDict = merge!(empty(cd), cd, others...)
@inline Base.merge(combine::Function, cd::CD, others::CD...) where {CD<:CompositeDict} = merge!(combine, empty(cd), cd, others...)
@inline Base.empty(cd::CompositeDict) = rawtype(typeof(cd))(dissolve(cd, empty)...)
@inline Base.iterate(cd::CompositeDict) = iterate(getcontent(cd, :contents))
@inline Base.iterate(cd::CompositeDict, state) = iterate(getcontent(cd, :contents), state)
@inline Base.keys(cd::CompositeDict) = keys(getcontent(cd, :contents))
@inline Base.values(cd::CompositeDict) = values(getcontent(cd, :contents))
@inline Base.pairs(cd::CompositeDict) = pairs(getcontent(cd, :contents))

# Vector spaces
"""
    VectorSpace{B} <: AbstractVector{B}

Abstract vector space.
"""
abstract type VectorSpace{B} <: AbstractVector{B} end
@inline Base.:(==)(vs₁::VS, vs₂::VS) where {VS<:VectorSpace} = ==(efficientoperations, vs₁, vs₂)
@inline Base.:(==)(::VS₁, ::VS₂) where {VS₁<:VectorSpace, VS₂<:VectorSpace} = false
@inline Base.isequal(vs₁::VS, vs₂::VS) where {VS<:VectorSpace} = isequal(efficientoperations, vs₁, vs₂)
@inline Base.isequal(::VS₁, ::VS₂) where {VS₁<:VectorSpace, VS₂<:VectorSpace} = false
@inline Base.axes(vs::VectorSpace) = (Base.OneTo(length(vs)),)
@inline Base.size(vs::VectorSpace) = (length(vs),)
@inline Base.IndexStyle(::Type{<:VectorSpace}) = IndexLinear()
@inline Base.getindex(vs::VectorSpace, i::Integer) = getindex(vs, CartesianIndex(i, vs))

"""
    VectorSpaceStyle

Style of a concrete type of vector space.
"""
abstract type VectorSpaceStyle end
@inline VectorSpaceStyle(vs::VectorSpace) = VectorSpaceStyle(typeof(vs))
@inline Base.length(vs::VectorSpace) = length(VectorSpaceStyle(vs), vs)
@inline Base.CartesianIndex(i::Integer, vs::VectorSpace) = CartesianIndex(VectorSpaceStyle(vs), i, vs)
@inline Base.getindex(vs::VectorSpace, i::CartesianIndex) = getindex(VectorSpaceStyle(vs), vs, i)
@inline Base.searchsortedfirst(vs::VectorSpace{B}, basis::B) where B = searchsortedfirst(VectorSpaceStyle(vs), vs, basis)
@inline Base.in(basis::B, vs::VectorSpace{B}) where B = in(VectorSpaceStyle(vs), basis, vs)

"""
    VectorSpaceGeneral <: VectorSpaceStyle

Default vector space style.
"""
struct VectorSpaceGeneral <: VectorSpaceStyle end
@inline VectorSpaceStyle(::Type{<:VectorSpace}) = VectorSpaceGeneral()
@inline Base.CartesianIndex(::VectorSpaceStyle, i::Integer, ::VectorSpace) = CartesianIndex(i)
@inline Base.searchsortedfirst(::VectorSpaceStyle, vs::VectorSpace{B}, basis::B) where B = invoke(searchsortedfirst, Tuple{Any, Any}, vs, basis)
@inline Base.in(::VectorSpaceStyle, basis::B, vs::VectorSpace{B}) where B = invoke(in, Tuple{Any, Any}, basis, vs)

"""
    VectorSpaceEnumerative <: VectorSpaceStyle

Enumerative vector space style, which indicates that the vector space has a predefined content named `contents` that contains all its bases.
"""
struct VectorSpaceEnumerative <: VectorSpaceStyle end
@inline Base.length(::VectorSpaceEnumerative, vs::VectorSpace) = length(getcontent(vs, :contents))
@inline function Base.getindex(::VectorSpaceEnumerative, vs::VectorSpace, i::CartesianIndex)
    contents = getcontent(vs, :contents)
    return contents[first(axes(contents))[i]]
end

"""
    VectorSpaceDirectSummed <: VectorSpaceStyle

Vector space style which indicates that a vector space is the direct sum of its sub-components.
"""
struct VectorSpaceDirectSummed <: VectorSpaceStyle end
@inline Base.length(::VectorSpaceDirectSummed, vs::VectorSpace) = mapreduce(length, +, getcontent(vs, :contents))
@inline Base.CartesianIndex(::VectorSpaceDirectSummed, i::Integer, vs::VectorSpace) = DirectSummedIndices(map(eachindex, getcontent(vs, :contents)))[i]
@inline function Base.getindex(::VectorSpaceDirectSummed, vs::VectorSpace, i::CartesianIndex{3})
    contents = getcontent(vs, :contents)
    m, n = i.I
    return contents[m][n]
end

"""
    VectorSpaceDirectProducted{Order} <: VectorSpaceStyle

Vector space style which indicates that a vector space is the direct product of its sub-components.

The type parameter `Order` must be either `:forward` or `:backward`:
1) `:forward`: the direct product iterates over the first sub-component first like a Julia array;
2) `:backward`: the direct product iterates over the last sub-component first like a C/C++ array.
"""
struct VectorSpaceDirectProducted{Order} <: VectorSpaceStyle end
@inline function VectorSpaceDirectProducted(order::Symbol)
    @assert order∈(:forward, :backward) "VectorSpaceDirectProducted error: not supported order."
    return VectorSpaceDirectProducted{order}()
end
@inline shape(vs::VectorSpace) = concatenate(map(axes, getcontent(vs, :contents))...)
@inline Base.length(::VectorSpaceDirectProducted, vs::VectorSpace) = mapreduce(length, *, shape(vs))
@inline Base.CartesianIndex(::VectorSpaceDirectProducted{Order}, i::Integer, vs::VectorSpace) where Order = DirectProductedIndices{Order}(shape(vs))[i]
@inline Base.getindex(::VectorSpaceDirectProducted, vs::VectorSpace, index::CartesianIndex) = convert(eltype(vs), index, vs)
function Base.searchsortedfirst(::VectorSpaceDirectProducted{Order}, vs::VectorSpace{B}, basis::B) where {Order, B}
    index = convert(CartesianIndex, basis, vs)
    indexes = DirectProductedIndices{Order}(shape(vs))
    return searchsortedfirst(indexes, index)
end
@inline Base.in(::VectorSpaceDirectProducted{Order}, basis::B, vs::VectorSpace{B}) where {Order, B} = convert(CartesianIndex, basis, vs) in DirectProductedIndices{Order}(shape(vs))

"""
    VectorSpaceZipped <: VectorSpaceStyle

Vector space style which indicates that a vector space is the zip of its sub-components.
"""
struct VectorSpaceZipped <: VectorSpaceStyle end
@inline Base.length(::VectorSpaceZipped, vs::VectorSpace) = mapreduce(length, min, getcontent(vs, :contents))
@inline function Base.CartesianIndex(::VectorSpaceZipped, i::Integer, vs::VectorSpace)
    contents = getcontent(vs, :contents)
    index = map(content->first(axes(content))[i], contents)
    return CartesianIndex(index)
end
@inline Base.getindex(::VectorSpaceZipped, vs::VectorSpace, index::CartesianIndex) = convert(eltype(vs), index, vs)

"""
    DirectSummedVectorSpace{B, T<:Tuple} <: VectorSpace{B}

A simple direct summed vector space.
"""
struct DirectSummedVectorSpace{B, T<:Tuple} <: VectorSpace{B}
    contents::T
    DirectSummedVectorSpace(contents::Tuple) = new{mapreduce(eltype, typejoin, contents), typeof(contents)}(contents)
end
@inline VectorSpaceStyle(::Type{<:DirectSummedVectorSpace}) = VectorSpaceDirectSummed()

"""
    DirectSummedVectorSpace(contents...)
    DirectSummedVectorSpace(contents::Tuple)

Construct a simple direct summed vector space.
"""
@inline DirectSummedVectorSpace(contents...) = DirectSummedVectorSpace(contents)

"""
    DirectProductedVectorSpace{Order, B, T<:Tuple} <: VectorSpace{B}

A simple direct producted vector space.
"""
struct DirectProductedVectorSpace{Order, B, T<:Tuple} <: VectorSpace{B}
    contents::T
    DirectProductedVectorSpace{Order, B}(contents::Tuple) where {Order, B} = new{Order, B, typeof(contents)}(contents)
end
@inline VectorSpaceStyle(::Type{<:DirectProductedVectorSpace{Order}}) where Order = VectorSpaceDirectProducted(Order)
@inline Base.convert(::Type{B}, index::CartesianIndex, vs::DirectProductedVectorSpace{Order, B}) where {B, Order} = B(map(getindex, vs.contents, index.I))
@inline Base.convert(::Type{<:CartesianIndex}, basis::B, vs::DirectProductedVectorSpace{Order, B}) where {B<:Tuple, Order} = CartesianIndex(map(searchsortedfirst, vs.contents, basis)...)

"""
    DirectProductedVectorSpace{Order}(contents...) where Order
    DirectProductedVectorSpace{Order}(contents::Tuple) where Order
    DirectProductedVectorSpace{Order, B}(contents...) where {Order, B}
    DirectProductedVectorSpace{Order, B}(contents::Tuple) where {Order, B}

Construct a simple direct producted vector space.
"""
@inline DirectProductedVectorSpace{Order}(contents...) where Order = DirectProductedVectorSpace{Order}(contents)
@inline DirectProductedVectorSpace{Order}(contents::Tuple) where Order = DirectProductedVectorSpace{Order, Tuple{map(eltype, fieldtypes(typeof(contents)))...}}(contents...)
@inline DirectProductedVectorSpace{Order, B}(contents...) where {Order, B} = DirectProductedVectorSpace{Order, B}(contents)

"""
    ZippedVectorSpace{B, T<:Tuple} <: VectorSpace{B}

A simple zipped vector space.
"""
struct ZippedVectorSpace{B, T<:Tuple} <: VectorSpace{B}
    contents::T
    ZippedVectorSpace{B}(contents::Tuple) where B = new{B, typeof(contents)}(contents)
end
@inline VectorSpaceStyle(::Type{<:ZippedVectorSpace}) = VectorSpaceZipped()
@inline Base.convert(::Type{B}, index::CartesianIndex, vs::ZippedVectorSpace{B}) where B = B(map(getindex, vs.contents, index.I))

"""
    ZippedVectorSpace(contents...)
    ZippedVectorSpace(contents::Tuple)
    ZippedVectorSpace{B}(contents...) where B
    ZippedVectorSpace{B}(contents::Tuple) where B

Construct a simple zipped vector space.
"""
@inline ZippedVectorSpace(contents...) = ZippedVectorSpace(contents)
@inline ZippedVectorSpace(contents::Tuple) = ZippedVectorSpace{Tuple{map(eltype, fieldtypes(typeof(contents)))...}}(contents...)
@inline ZippedVectorSpace{B}(contents...) where B = ZippedVectorSpace{B}(contents)

end # module
