module Toolkit

using Formatting: FormatSpec, fmt
using Printf: @printf
using StaticArrays: SVector

import QuantumLattices: ⊕, ⊗

# Utilities
export atol, rtol, Float
export Segment, concatenate, decimaltostr, delta, ordinal

# Combinatorics
export Combinatorics, Combinations, DuplicateCombinations, DuplicatePermutations, Permutations

# Traits
export commontype, dissolve, fulltype, rawtype
export hasparameter, isparameterbound, isparameterbounds, parametercount, parametername, parameternames, parameterorder, parameterpair, parameterpairs, parametertype, parametertypes, promoteparameters, reparameter 
export contentcount, contentname, contentnames, contentorder, contenttype, contenttypes, getcontent, hascontent
export efficientoperations

# Composite structures
export CompositeDict, CompositeNTuple, CompositeTuple, CompositeVector, NamedContainer

# Named vectors
export HomoNamedVector, NamedVector

# Vector spaces
export CompositeNamedVectorSpace, DirectProductedNamedVectorSpace, NamedVectorSpace, ParameterSpace, SimpleNamedVectorSpace, VectorSpace, ZippedNamedVectorSpace
export VectorSpaceCartesian, VectorSpaceDirectProducted, VectorSpaceDirectSummed, VectorSpaceEnumerative, VectorSpaceStyle, VectorSpaceZipped
export shape

# Simple trees
export simpletreedepth, simpletreewidth, AbstractSimpleTree, SimpleTree, SimpleTreeCore
export addnode!, ancestor, children, deletenode!, descendants, isleaf, leaves, level, move!, root, siblings, subtree

# Utilities
"Absolute tolerance for float numbers."
const atol = 5 * 10^-14

"Relative tolerance for float numbers."
const rtol = √atol

"Default float type."
const Float = Float64

"""
    decimaltostr(number, ::Int=5)
    decimaltostr(number::Integer, n::Int=5)
    decimaltostr(number::Rational, n::Int=5)
    decimaltostr(number::AbstractFloat, n::Int=5)
    decimaltostr(number::Complex, n::Int=5)

Convert a number to a string with at most `n` decimal places.
"""
@inline decimaltostr(number, ::Int=5) = repr(number)
@inline decimaltostr(number::Integer, ::Int=5) = string(number)
@inline decimaltostr(number::Rational, ::Int=5) = string(number)
function decimaltostr(number::AbstractFloat, n::Int=5)
    if number == 0.0
        result = "0.0"
    elseif 10^-5 < abs(number) < 10^6
        result = rstrip(fmt(FormatSpec(".$(n)f"), number), '0')
        (result[end] == '.') && (result = result * '0')
    else
        result = fmt(FormatSpec(".$(n)e"), number)
        epos = findfirst(isequal('e'), result)
        temp = rstrip(result[1:epos-1], '0')
        result = (temp[end] == '.') ? (temp * "0" * result[epos:end]) : (temp * result[epos:end])
    end
    result
end
function decimaltostr(number::Complex, n::Int=5)
    sreal = (real(number) == 0) ? "0" : decimaltostr(real(number), n)
    simag = (imag(number) == 0) ? "0" : decimaltostr(imag(number), n)
    result = ""
    (sreal == "0") || (result = result * sreal)
    (simag == "0") || (result = ((simag[1] == '-') ? (result * simag) : (length(result) == 0) ? simag : (result * "+" * simag)) * "im")
    (length(result) == 0) && (result = "0.0")
    return result
end

"""
    ordinal(number::Integer)

Convert a positive number to its corresponding ordinal.
"""
@inline function ordinal(number::Integer)
    @assert number > 0 "ordinal error: input number must be positive."
    (number == 1) ? "1st" : (number == 2) ? "2nd" : (number == 3) ? "3rd" : "$(number)th"
end

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
    searchsortedfirst(table, basis; compare=<) -> Int

Use the binary search method to find the position of a basis in a sorted table so that the order is preserved if the basis in inserted in that position.
"""
function Base.searchsortedfirst(table, basis; compare=<)
    lo, hi = 0, length(table)+1
    @inbounds while lo < hi-1
        m = (lo+hi) >>> 1
        if compare(table[m], basis)
            lo = m
        else
            hi = m
        end
    end
    return hi
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

Combinations of M elements from contents. Duplicates are not allowed.
"""
struct Combinations{M, C} <: Combinatorics{M, C}
    contents::C
    N::Int
    Combinations{M}(contents::C) where {M, C} = new{M, C}(contents, length(contents))
end
@inline Base.length(c::Combinations{M}) where M = binomial(c.N, M)
Base.iterate(c::Combinations{M}) where M = (M > c.N) ? nothing : (M == 0) ? ((), [c.N+2]) : (ntuple(i->c.contents[i], Val(M)), nextmstate!(collect(1:M), c.N, M))
Base.iterate(c::Combinations{M}, state) where M = (state[1] > c.N-M+1) ? nothing : (ntuple(i->c.contents[state[i]], Val(M)), nextmstate!(state, c.N, M))
function nextmstate!(state::Vector{Int}, N::Int, M::Int)
    for i = M:-1:1
        state[i] += 1
        (state[i] > N-(M-i)) && continue
        for j = i+1:M
            state[j] = state[j-1] + 1
        end
        break
    end
    state
end

"""
    DuplicateCombinations{M}(contents::C) where {M, C}

Combinations of M elements from contents. Duplicates are allowed.
"""
struct DuplicateCombinations{M, C} <: Combinatorics{M, C}
    contents::C
    N::Int
    DuplicateCombinations{M}(contents::C) where {M, C} = new{M, C}(contents, length(contents))
end
@inline Base.length(c::DuplicateCombinations{M}) where M = binomial(c.N+M-1, c.N-1)
function Base.iterate(c::DuplicateCombinations{M}) where M
    (M == 0) ? ((), [c.N+1]) : (ntuple(i->c.contents[1], Val(M)), nextdmstate!(collect(ntuple(i->1, M)), c.N, M))
end
Base.iterate(c::DuplicateCombinations{M}, state) where M = (state[1] > c.N) ? nothing : (ntuple(i->c.contents[state[i]], Val(M)), nextdmstate!(state, c.N, M))
function nextdmstate!(state::Vector{Int}, N::Int, M::Int)
    for i = M:-1:1
        state[i] += 1
        (state[i] > N) && continue
        for j = i+1:M
            state[j] = state[j-1]
        end
        break
    end
    state
end

"""
    Permutations{M}(contents::C) where {M, C}

Permutations of M elements from contents. Duplicates are not allowed.
"""
struct Permutations{M, C} <: Combinatorics{M, C}
    contents::C
    N::Int
    Permutations{M}(contents::C) where {M, C} = new{M, C}(contents, length(contents))
end
@inline Base.length(p::Permutations{M}) where M = (0 <= M <= p.N) ? prod((p.N-M+1):p.N) : 0
function Base.iterate(p::Permutations{M}) where M
    ((p.N == 0) && (M > 0) || (0 < p.N < M)) ? nothing : (state = collect(1:p.N); (ntuple(i->p.contents[state[i]], Val(M)), nextpstate!(state, p.N, M)))
end
function Base.iterate(p::Permutations{M}, state) where M
    ((p.N == 0) && (M > 0) || (0 < p.N < max(state[1], M))) ? nothing : (ntuple(i->p.contents[state[i]], Val(M)), nextpstate!(state, p.N, M))
end
function nextpstate!(state::Vector{Int}, N::Int, M::Int)
    (M <= 0) && return [N+1]
    (M < N) && (j = M + 1; while (j <= N) && (state[M] >= state[j]); j += 1; end)
    if (M < N) && (j <= N)
        state[M], state[j] = state[j], state[M]
    else
        (M < N) && reverse!(state, M+1)
        i = M - 1
        while (i >= 1) && (state[i] >= state[i+1]); i -= 1; end
        if i > 0
            j = N
            while (j > i) && (state[i] >= state[j]); j -= 1; end
            state[i], state[j] = state[j], state[i]
            reverse!(state, i+1)
        else
            state[1] = N + 1
        end
    end
    return state
end

"""
    DuplicatePermutations{M}(contents::C) where {M, C}

Permutations of M elements from contents. Duplicates are allowed.
"""
struct DuplicatePermutations{M, C} <: Combinatorics{M, C}
    contents::C
    N::Int
    DuplicatePermutations{M}(contents::C) where {M, C} = new{M, C}(contents, length(contents))
end
@inline Base.length(p::DuplicatePermutations{M}) where M = p.N^M
function Base.iterate(p::DuplicatePermutations{M}) where M
    indices = CartesianIndices(ntuple(i->p.N, M|>Val))
    index = iterate(indices)
    isnothing(index) && return nothing
    return ntuple(i->p.contents[reverse(index[1].I)[i]], M|>Val), (indices, index[2])
end
function Base.iterate(p::DuplicatePermutations{M}, state) where M
    index = iterate(state[1], state[2])
    isnothing(index) && return nothing
    return ntuple(i->p.contents[reverse(index[1].I)[i]], M|>Val), (state[1], index[2])
end

# Traits
## with type themselves
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
    while nameof(result) != termination && nameof(result) != :Any
        result = supertype(result)
    end
    (nameof(result) != termination) && error("_supertype error: termination is not the name of a valid supertype of $(nameof(T)).")
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
    return isa(result, TypeVar) ? result.ub : result
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
    isparameterbound(::Type{}, ::Val{}, ::Any) -> Bool

For a type `T`, judge whether a type `D` should be considered as the upper bound of one of its type parameters.
"""
@inline isparameterbound(::Type{T}, i::Integer, D) where T = isparameterbound(T, Val(i), D)
@inline isparameterbound(::Type{T}, name::Symbol, D) where T = isparameterbound(T, Val(name), D)
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
    promoteparameters(::Type{T1}, ::Type{T2}) where {T1<:NamedTuple, T2<:NamedTuple}

Promote the types specified by two named tuples with the same names accordingly.

The result is stored in the type parameters of a `NamedTuple`.
"""
@inline @generated function promoteparameters(::Type{T1}, ::Type{T2}) where {T1<:NamedTuple, T2<:NamedTuple}
    exprs = []
    names = Tuple(union(fieldnames(T1), fieldnames(T2)))
    for (i, name) in enumerate(names)
        hasfield(T1, name) && (F1 = fieldtype(T1, name))
        hasfield(T2, name) && (F2 = fieldtype(T2, name))
        if hasfield(T1, name) && hasfield(T2, name)
            push!(exprs, :(promote_type($F1, $F2)))
        else
            push!(exprs, hasfield(T1, name) ? F1 : F2)
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
    dissolve(m, f::Function=identity, args::Tuple=(), kwargs::NamedTuple=NamedTuple()) -> Tuple

Convert `m` to a tuple by the function `f` applied elementally to its contents with the extra positional arguments (`args`) and keyword arguments (`kwargs`). 

The underlying called interface is the `dissolve` function when `f` is applied to each content of `m`:
```julia
dissolve(m, Val(name), f, args, kwargs)
```
Here, `name` is the name of a content of `m`.

Basically, the rule of how `f` operates on each field of `m` can be overridden by redefining the above `dissolve` function.
!!!note
   The default `dissolve` function ignores the operation of function `f` and just return the content value of `m`.
"""
@inline function dissolve(m, f::Function=identity, args::Tuple=(), kwargs::NamedTuple=NamedTuple())
    dissolvehelper(m, f, args, kwargs, m|>typeof|>contentnames|>Val)
end
@inline @generated function dissolvehelper(m, f::Function, args::Tuple, kwargs::NamedTuple, ::Val{names}) where names
    exprs = []
    for name in names
        name = Val(name)
        push!(exprs, :(dissolve(m, $name, f, args, kwargs)))
    end
    return Expr(:tuple, exprs...)
end

"""
    dissolve(m, ::Val{name}, f::Function, args::Tuple, kwargs::NamedTuple) where name

Dissolve the content specified by `name` of `m` by the function `f` applied with the extra positional arguments (`args`) and keyword arguments (`kwargs`).
"""
@inline dissolve(m, ::Val{name}, f::Function, args::Tuple, kwargs::NamedTuple) where name = getcontent(m, name|>Val)

## with type operations
struct EfficientOperations end
"""
    efficientoperations

Indicate that the efficient operations, i.e. "=="/"isequal", "<"/"isless" or "replace", will be used.
"""
const efficientoperations = EfficientOperations()

"""
    ==(::EfficientOperations, o1, o2) -> Bool

Compare two objects and judge whether they are equivalent to each other.
"""
@inline @generated function Base.:(==)(::EfficientOperations, o1, o2)
    fcount = fieldcount(o1)
    if fcount == fieldcount(o2)
        if fcount == 0
            return true
        else
            expr = :(getfield(o1, 1) == getfield(o2, 1))
            for i = 2:fcount
                expr = Expr(:&&, expr, :(getfield(o1, $i) == getfield(o2, $i)))
            end
            return expr
        end
    else
        return false
    end
end

"""
    isequal(::EfficientOperations, o1, o2) -> Bool

Compare two objects and judge whether they are equivalent to each other.
"""
@inline @generated function Base.isequal(::EfficientOperations, o1, o2)
    fcount = fieldcount(o1)
    if fcount == fieldcount(o2)
        if fcount == 0
            return true
        else
            expr = :(isequal(getfield(o1, 1), getfield(o2, 1)))
            for i = 2:fcount
                expr = Expr(:&&, expr, :(isequal(getfield(o1, $i), getfield(o2, $i))))
            end
            return expr
        end
    else
        return false
    end
end

"""
    <(::EfficientOperations, o1, o2) -> Bool

Compare two objects and judge whether the first is less than the second.
"""
@inline @generated function Base.:<(::EfficientOperations, o1, o2)
    n1, n2 = fieldcount(o1), fieldcount(o2)
    N = min(n1, n2)
    expr = (n1 < n2) ? Expr(:if, :(getfield(o1, $N) == getfield(o2, $N)), true, false) : false
    expr = Expr(:if, :(getfield(o1, $N) < getfield(o2, $N)), true, expr)
    for i in range(N-1, stop=1, step=-1)
        expr = Expr(:if, :(getfield(o1, $i) > getfield(o2, $i)), false, expr)
        expr = Expr(:if, :(getfield(o1, $i) < getfield(o2, $i)), true, expr)
    end
    return expr
end

"""
    isless(::EfficientOperations, o1, o2) -> Bool

Compare two objects and judge whether the first is less than the second.
"""
@inline @generated function Base.isless(::EfficientOperations, o1, o2)
    n1, n2 = fieldcount(o1), fieldcount(o2)
    N = min(n1, n2)
    expr = n1 < n2 ? Expr(:if, :(getfield(o1, $N) == getfield(o2, $N)), true, false) : false
    expr = Expr(:if, :(isless(getfield(o1, $N), getfield(o2, $N))), true, expr)
    for i in range(N-1, stop=1, step=-1)
        expr = Expr(:if, :(isless(getfield(o2, $i), getfield(o1, $i))), false, expr)
        expr = Expr(:if, :(isless(getfield(o1, $i), getfield(o2, $i))), true, expr)
    end
    return expr
end

"""
    isapprox(::EfficientOperations, o1, o2; atol=atol, rtol=rtol) -> Bool
    isapprox(::EfficientOperations, fields::Union{Union{Nothing, Integer, Symbol}, Tuple{Vararg{Union{Integer, Symbol}}}}, o1, o2; atol=atol, rtol=rtol) -> Bool
    isapprox(::EfficientOperations, ::Val{fields}, o1, o2; atol=atol, rtol=rtol) where fields -> Bool

Compare two objects and judge whether they are inexactly equivalent to each other.
"""
@inline Base.isapprox(::EfficientOperations, o1, o2; atol=atol, rtol=rtol) = isapprox(efficientoperations, nothing, o1, o2, atol=atol, rtol=rtol)
@inline function Base.isapprox(::EfficientOperations, fields::Union{Nothing, Union{Integer, Symbol}, Tuple{Vararg{Union{Integer, Symbol}}}}, o1, o2; atol=atol, rtol=rtol)
    isapprox(efficientoperations, fields|>Val, o1, o2; atol=atol, rtol=rtol)
end
@inline @generated function Base.isapprox(::EfficientOperations, ::Val{fields}, o1, o2; atol=atol, rtol=rtol) where fields
    (fieldcount(o1)≠fieldcount(o2) || fieldnames(o1)≠fieldnames(o2)) && return false
    isnothing(fields) && (fields = ntuple(i->i, fieldcount(o1)))
    isa(fields, Union{Integer, Symbol}) && (fields = (fields,))
    fields = Set(isa(field, Symbol) ? findfirst(isequal(field), fieldnames(o1)) : field for field in fields)
    exprs = []
    for i = 1:fieldcount(o1)
        if i∈fields
            push!(exprs, :(isapprox(getfield(o1, $i), getfield(o2, $i); atol=atol, rtol=rtol)::Bool))
        else
            push!(exprs, :(isequal(getfield(o1, $i), getfield(o2, $i))::Bool))
        end
    end
    return Expr(:&&, exprs...)
end

"""
    replace(::EfficientOperations, o; kwargs...) -> typeof(o)

Return a copy of the input object with some of the field values replaced by the keyword arguments.
"""
@inline @generated function Base.replace(::EfficientOperations, o; kwargs...)
    exprs = [:(get(kwargs, $name, getfield(o, $name))) for name in QuoteNode.(fieldnames(o))]
    return :(rawtype(typeof(o))($(exprs...)))
end

# Composite structures
"""
    CompositeTuple{T<:Tuple}

A composite tuple can be considered as a tuple that is implemented by including an ordinary `Tuple` as its data attribute.
"""
abstract type CompositeTuple{T<:Tuple} end
@inline contentnames(::Type{<:CompositeTuple}) = (:contents,)
@inline dissolve(ct::CompositeTuple, ::Val{:contents}, f::Function, args::Tuple, kwargs::NamedTuple) = f(getcontent(ct, :contents), args...; kwargs...)

@inline Base.length(::CompositeTuple{T}) where T = fieldcount(T)
@inline Base.length(::Type{<:CompositeTuple{T}}) where T = fieldcount(T)
@inline Base.eltype(::CompositeTuple{T}) where T = eltype(T)
@inline Base.eltype(::Type{<:CompositeTuple{T}}) where T = eltype(T)
@inline Base.hash(ct::CompositeTuple, h::UInt) = hash(dissolve(ct), h)
@inline Base.:(==)(ct1::CompositeTuple, ct2::CompositeTuple) = ==(efficientoperations, ct1, ct2)
@inline Base.isequal(ct1::CompositeTuple, ct2::CompositeTuple) = isequal(efficientoperations, ct1, ct2)
@inline Base.getindex(ct::CompositeTuple, i::Union{<:Integer, CartesianIndex}) = getcontent(ct, :contents)[i]
@inline Base.getindex(ct::CompositeTuple, inds) = rawtype(typeof(ct))(dissolve(ct, getindex, (inds,))...)
@inline Base.lastindex(ct::CompositeTuple) = lastindex(getcontent(ct, :contents))
@inline Base.iterate(ct::CompositeTuple) = iterate(getcontent(ct, :contents))
@inline Base.iterate(ct::CompositeTuple, state) = iterate(getcontent(ct, :contents), state)
@inline Base.iterate(rv::Iterators.Reverse{<:CompositeTuple}, state=length(rv.itr)) = (state < 1) ? nothing : (rv.itr[state], state-1)
@inline Base.keys(ct::CompositeTuple) = keys(getcontent(ct, :contents))
@inline Base.values(ct::CompositeTuple) = values(getcontent(ct, :contents))
@inline Base.pairs(ct::CompositeTuple) = pairs(getcontent(ct, :contents))
@inline Base.reverse(ct::CompositeTuple) = rawtype(typeof(ct))(dissolve(ct, reverse)...)
@inline Base.Tuple(ct::CompositeTuple) = getcontent(ct, :contents)

"""
    CompositeNTuple{N, T}

A composite ntuple can be considered as a ntuple that is implemented by including an ordinary `NTuple` as its data attribute.

Alias for `CompositeTuple{NTuple{N, T}}`.
"""
const CompositeNTuple{N, T} = CompositeTuple{NTuple{N, T}}

"""
    CompositeVector{T}

A composite vector can be considered as a vector that is implemented by including a concrete subtype of `AbstractVector` as its data attribute.
"""
abstract type CompositeVector{T} <: AbstractVector{T} end
@inline contentnames(::Type{<:CompositeVector}) = (:contents,)
@inline dissolve(cv::CompositeVector, ::Val{:contents}, f::Function, args::Tuple, kwargs::NamedTuple) = f(getcontent(cv, :contents), args...; kwargs...)

@inline Base.size(cv::CompositeVector) = size(getcontent(cv, :contents))
@inline Base.size(cv::CompositeVector, i) = size(getcontent(cv, :contents), i)
@inline Base.length(cv::CompositeVector) = length(getcontent(cv, :contents))
@inline Base.:(==)(cv1::CompositeVector, cv2::CompositeVector) = ==(efficientoperations, cv1, cv2)
@inline Base.isequal(cv1::CompositeVector, cv2::CompositeVector) = isequal(efficientoperations, cv1, cv2)
@inline Base.getindex(cv::CompositeVector, i::Union{<:Integer, CartesianIndex}) = getcontent(cv, :contents)[i]
@inline Base.getindex(cv::CompositeVector, inds) = rawtype(typeof(cv))(dissolve(cv, getindex, (inds,))...)
@inline Base.lastindex(cv::CompositeVector) = lastindex(getcontent(cv, :contents))
@inline Base.setindex!(cv::CompositeVector, value, inds) = (getcontent(cv, :contents)[inds] = value)
@inline Base.push!(cv::CompositeVector, values...) = (push!(getcontent(cv, :contents), values...); cv)
@inline Base.pushfirst!(cv::CompositeVector, values...) = (pushfirst!(getcontent(cv, :contents), values...); cv)
@inline Base.insert!(cv::CompositeVector, index::Integer, value) = (insert!(getcontent(cv, :contents), index, value); cv)
@inline Base.append!(cv::CompositeVector, values) = (append!(getcontent(cv, :contents), values); cv)
@inline Base.prepend!(cv::CompositeVector, values) = (prepend!(getcontent(cv, :contents), values); cv)
@inline Base.splice!(cv::CompositeVector, index::Integer, replacement=Base._default_splice) = splice!(getcontent(cv, :contents), index, replacement)
@inline Base.splice!(cv::CompositeVector, range::UnitRange{<:Integer}, replacement=Base._default_splice) = rawtype(typeof(cv))(dissolve(cv, splice!, (range, replacement))...)
@inline Base.deleteat!(cv::CompositeVector, indices) = (deleteat!(getcontent(cv, :contents), indices); cv)
@inline Base.pop!(cv::CompositeVector) = pop!(getcontent(cv, :contents))
@inline Base.popfirst!(cv::CompositeVector) = popfirst!(getcontent(cv, :contents))
@inline Base.empty!(cv::CompositeVector) = (empty!(getcontent(cv, :contents)); cv)
@inline Base.empty(cv::CompositeVector) = rawtype(typeof(cv))(dissolve(cv, empty)...)
@inline Base.reverse(cv::CompositeVector) =  rawtype(typeof(cv))(dissolve(cv, reverse)...)
@inline Base.similar(cv::CompositeVector, dtype::Type=eltype(cv), dims::Tuple{Vararg{Int}}=size(cv)) = rawtype(typeof(cv))(dissolve(cv, similar, (dtype, dims))...)
@inline Base.iterate(cv::CompositeVector, state=1) = iterate(getcontent(cv, :contents), state)
@inline Base.keys(cv::CompositeVector) = keys(getcontent(cv, :contents))
@inline Base.values(cv::CompositeVector) = values(getcontent(cv, :contents))
@inline Base.pairs(cv::CompositeVector) = pairs(getcontent(cv, :contents))

"""
    CompositeDict{K, V}

A composite dict can be considered as a dict that is implemented by including a concrete subtype of `AbstractDict` as its data attribute.
"""
abstract type CompositeDict{K, V} <: AbstractDict{K, V} end
@inline contentnames(::Type{<:CompositeDict}) = (:contents,)
@inline dissolve(cd::CompositeDict, ::Val{:contents}, f::Function, args::Tuple, kwargs::NamedTuple) = f(getcontent(cd, :contents), args...; kwargs...)

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

"""
    NamedContainer{T, Names} = NamedTuple{Names, <:Tuple{Vararg{T}}}

NamedContainer is just a wrapper of Julia NamedTuple, but not a composite type.
"""
const NamedContainer{T, Names} = NamedTuple{Names, <:Tuple{Vararg{T}}}
@inline NamedContainer{Names}(contents::Tuple) where Names = NamedTuple{Names}(contents)

# Named vectors
"""
    NamedVector

Abstract type for all named vectors.
"""
abstract type NamedVector end
@inline Base.:(==)(nv₁::NamedVector, nv₂::NamedVector) = keys(nv₁)==keys(nv₂) && values(nv₁)==values(nv₂)
@inline Base.isequal(nv₁::NamedVector, nv₂::NamedVector) = isequal(keys(nv₁), keys(nv₂)) && isequal(values(nv₁), values(nv₂))
@inline Base.getindex(nv::NamedVector, index::Int) = getfield(nv, index)
@inline Base.setindex!(nv::NamedVector, value, index::Int) = setfield!(nv, index, value)
@inline Base.:<(nv₁::NamedVector, nv₂::NamedVector) = <(efficientoperations, nv₁, nv₂)
@inline Base.isless(nv₁::NamedVector, nv₂::NamedVector) = isless(efficientoperations, nv₁, nv₂)
@inline Base.hash(nv::NamedVector, h::UInt) = hash(values(nv), h)
@inline Base.length(::Type{NV}) where {NV<:NamedVector} = fieldcount(NV)
@inline Base.length(nv::NamedVector) = length(typeof(nv))
@inline Base.iterate(nv::NamedVector, state=1) = (state > length(nv)) ? nothing : (getfield(nv, state), state+1)
@inline Base.iterate(rv::Iterators.Reverse{<:NamedVector}, state=length(rv.itr)) = (state < 1) ? nothing : (getfield(rv.itr, state), state-1)
@inline Base.keys(nv::NamedVector) = fieldnames(typeof(nv))
@inline Base.values(nv::NamedVector) = convert(Tuple, nv)
@inline Base.pairs(nv::NamedVector) = Base.Generator(=>, keys(nv), values(nv))
Base.show(io::IO, nv::NamedVector) = @printf io "%s(%s)" nameof(typeof(nv)) join(repr.(values(nv)), ", ")

"""
    convert(::Type{Tuple}, nv::NamedVector) -> Tuple
    convert(::Type{NV}, nv::Tuple) where {NV<:NamedVector} -> NV

Convert a named vector to tuple and vice versa.
"""
@generated Base.convert(::Type{Tuple}, nv::NamedVector) = Expr(:tuple, (:(getfield(nv, $i)) for i = 1:fieldcount(nv))...)
function Base.convert(::Type{NV}, nv::Tuple) where NV<:NamedVector
    @assert fieldcount(NV) == length(nv) "convert error: mismatched length between $NV($(fieldcount(NV))) and input tuple($(length(nv)))."
    return NV(nv...)
end

"""
    zero(::Type{NV}) where NV<:NamedVector -> NV
    zero(nv::NamedVector) -> typeof(nv)

Get a concrete `NamedVector` with all values being zero.
"""
@inline @generated function Base.zero(::Type{NV}) where NV<:NamedVector
    zeros = (zero(fieldtype(NV, i)) for i = 1:fieldcount(NV))
    return :(NV($(zeros...)))
end
@inline Base.zero(nv::NamedVector) = zero(typeof(nv))

"""
    replace(nv::NamedVector; kwargs...) -> typeof(nv)

Return a copy of a concrete `NamedVector` with some of the field values replaced by the keyword arguments.
"""
@inline Base.replace(nv::NamedVector; kwargs...) = replace(efficientoperations, nv; kwargs...)

"""
    map(f, nvs::NV...) where NV<:NamedVector -> NV

Apply function `f` element-wisely on the input named vectors.
"""
@generated function Base.map(f, nvs::NV...) where NV<:NamedVector
    exprs = Vector{Expr}(undef, fieldcount(NV))
    for i = 1:length(exprs)
        tmp = [:(nvs[$j][$i]) for j = 1:length(nvs)]
        exprs[i] = :(f($(tmp...)))
    end
    return :(($NV)($(exprs...)))
end

"""
    HomoNamedVector{T}

Abstract type for all homogeneous named vectors.
"""
abstract type HomoNamedVector{T} <: NamedVector end
@inline Base.eltype(::Type{<:HomoNamedVector{T}}) where T = T
@inline Base.eltype(nv::HomoNamedVector) = eltype(typeof(nv))

# Vector spaces
"""
    VectorSpace{B} <: AbstractVector{B}

Abstract vector space.
"""
abstract type VectorSpace{B} <: AbstractVector{B} end
@inline Base.:(==)(vs₁::VectorSpace, vs₂::VectorSpace) = ==(efficientoperations, vs₁, vs₂)
@inline Base.isequal(vs₁::VectorSpace, vs₂::VectorSpace) = isequal(efficientoperations, vs₁, vs₂)
@inline Base.size(vs::VectorSpace) = (length(vs),)

"""
    VectorSpaceStyle

The style of a concrete type of vector space.
"""
abstract type VectorSpaceStyle end
@inline VectorSpaceStyle(vs::VectorSpace) = VectorSpaceStyle(typeof(vs))

@inline Base.length(vs::VectorSpace) = length(VectorSpaceStyle(vs), vs)
@inline Base.getindex(vs::VectorSpace, i) = getindex(VectorSpaceStyle(vs), vs, i)
@inline Base.getindex(vs::VectorSpace, i::CartesianIndex{1}) = vs[i[1]]
@inline Base.getindex(style::VectorSpaceStyle, vs::VectorSpace, indexes) = map(index->getindex(style, vs, index), indexes)
@inline Base.issorted(vs::VectorSpace) = issorted(VectorSpaceStyle(vs), vs)
@inline Base.findfirst(basis::B, vs::VectorSpace{B}) where B = findfirst(VectorSpaceStyle(vs), basis, vs)
@inline Base.searchsortedfirst(vs::VectorSpace{B}, basis::B) where B = searchsortedfirst(VectorSpaceStyle(vs), vs, basis)
@inline Base.in(basis::B, vs::VectorSpace{B}) where B = in(VectorSpaceStyle(vs), basis, vs)

@inline Base.issorted(::VectorSpaceStyle, vs::VectorSpace) = false
@inline function Base.findfirst(::VectorSpaceStyle, basis::B, vs::VectorSpace{B}) where B
    if issorted(vs)
        index = invoke(searchsortedfirst, Tuple{Any, Any}, vs, basis)
        return (0<index<=length(vs) && vs[index]==basis) ? index : nothing
    end
    for i = 1:length(vs)
        vs[i]==basis && return i
    end
end
@inline function Base.searchsortedfirst(::VectorSpaceStyle, vs::VectorSpace{B}, basis::B) where B
    issorted(vs) && return invoke(searchsortedfirst, Tuple{Any, Any}, vs, basis)
    for i = 1:length(vs) 
        vs[i]==basis && return i
    end
    return length(vs)+1
end
@inline Base.in(::VectorSpaceStyle, basis::B, vs::VectorSpace{B}) where B = isa(findfirst(basis, vs), Integer)

"""
    VectorSpaceEnumerative <: VectorSpaceStyle

Enumerative vector space style, which indicates that the vector space has a predefined content named `contents` that contains all its bases.
"""
struct VectorSpaceEnumerative <: VectorSpaceStyle end
@inline Base.length(::VectorSpaceEnumerative, vs::VectorSpace) = length(getcontent(vs, :contents))
@inline Base.getindex(::VectorSpaceEnumerative, vs::VectorSpace, i::Integer) = getcontent(vs, :contents)[i]
@inline Base.in(::VectorSpaceEnumerative, vs::VectorSpace{B}, basis::B) where B = in(basis, getcontent(vs, :contents))

"""
    VectorSpaceCartesian <: VectorSpaceStyle

Cartesian vector space style, which indicates that every basis in it could be represented by a Cartesian index.
"""
struct VectorSpaceCartesian <: VectorSpaceStyle end
function shape end
@inline Base.length(::VectorSpaceCartesian, vs::VectorSpace) = mapreduce(length, *, shape(vs))
@inline Base.getindex(::VectorSpaceCartesian, vs::VectorSpace, i::Integer) = getindex(VectorSpaceCartesian(), vs, CartesianIndices(shape(vs))[i])
@inline Base.getindex(::VectorSpaceCartesian, vs::VectorSpace, i::CartesianIndex) = rawtype(eltype(vs))(i, vs)
@inline Base.issorted(::VectorSpaceCartesian, vs::VectorSpace) = true
@inline function Base.findfirst(::VectorSpaceCartesian, basis::B, vs::VectorSpace{B}) where B
    index, idx = CartesianIndex(basis, vs), CartesianIndices(shape(vs))
    index∉idx && return nothing
    return LinearIndices(shape(vs))[index-first(idx)+oneunit(typeof(index))]
end
@inline function Base.searchsortedfirst(::VectorSpaceCartesian, vs::VectorSpace{B}, basis::B) where B
    index, idx = CartesianIndex(basis, vs), CartesianIndices(shape(vs))
    index∈idx && return LinearIndices(shape(vs))[index-first(idx)+oneunit(typeof(index))]
    return searchsortedfirst(idx, index)
end
@inline Base.in(::VectorSpaceCartesian, basis::B, vs::VectorSpace{B}) where B = CartesianIndex(basis, vs)∈CartesianIndices(shape(vs))

"""
    VectorSpaceDirectSummed <: VectorSpaceStyle

Vector space style which indicates that a vector space is the direct sum of its sub-components.
"""
struct VectorSpaceDirectSummed <: VectorSpaceStyle end
@inline Base.length(::VectorSpaceDirectSummed, vs::VectorSpace) = mapreduce(length, +, getcontent(vs, :contents))
@inline function Base.getindex(::VectorSpaceDirectSummed, vs::VectorSpace, i::Integer)
    contents = getcontent(vs, :contents)
    dimsum = cumsum(map(length, contents))
    m = searchsortedfirst(dimsum, i)
    n = m>1 ? (i-dimsum[m-1]) : i
    return contents[m][n]
end

"""
    VectorSpaceDirectProducted <: VectorSpaceStyle

Vector space style which indicates that a vector space is the direct product of its sub-components.
"""
struct VectorSpaceDirectProducted <: VectorSpaceStyle end
@inline Base.length(::VectorSpaceDirectProducted, vs::VectorSpace) = mapreduce(length, *, getcontent(vs, :contents))
@inline function Base.getindex(::VectorSpaceDirectProducted, vs::VectorSpace, i::Integer)
    contents = getcontent(vs, :contents)
    index = CartesianIndices(map(length, contents))[i]
    return rawtype(eltype(vs))(map(getindex, contents, Tuple(index)), vs)
end

"""
    VectorSpaceZipped <: VectorSpaceStyle

Vector space style which indicates that a vector space is the zip of its sub-components.
"""
struct VectorSpaceZipped <: VectorSpaceStyle end
@inline Base.length(::VectorSpaceZipped, vs::VectorSpace) = mapreduce(length, min, getcontent(vs, :contents))
@inline function Base.getindex(::VectorSpaceZipped, vs::VectorSpace, i::Integer)
    return rawtype(eltype(vs))(map(content->getindex(content, i), getcontent(vs, :contents)), vs)
end

"""
    NamedVectorSpace{B} <: VectorSpace{B}

Abstract named vector space.
"""
abstract type NamedVectorSpace{B} <: VectorSpace{B} end
@inline Base.names(vs::NamedVectorSpace) = names(typeof(vs))
@inline Base.:(==)(vs₁::NamedVectorSpace, vs₂::NamedVectorSpace) = names(vs₁)==names(vs₂) && invoke(==, Tuple{VectorSpace, VectorSpace}, vs₁, vs₂)
@inline Base.isequal(vs₁::NamedVectorSpace, vs₂::NamedVectorSpace) = isequal(names(vs₁), names(vs₂)) && invoke(isequal, Tuple{VectorSpace, VectorSpace}, vs₁, vs₂)
@inline Base.pairs(vs::NamedVectorSpace) = NamedVectorSpacePairIteration(vs)
struct NamedVectorSpacePairIteration{V<:NamedVectorSpace, B<:NamedTuple} <: AbstractVector{B}
    namedvectorspace::V
    NamedVectorSpacePairIteration(vs::NamedVectorSpace) = new{typeof(vs), pairtype(typeof(vs))}(vs)
end
@inline Base.size(ps::NamedVectorSpacePairIteration) = size(ps.namedvectorspace)

"""
    SimpleNamedVectorSpace{N, B} <: NamedVectorSpace{B}

Abstract simple named vector space.
"""
abstract type SimpleNamedVectorSpace{N, B} <: NamedVectorSpace{B} end
@inline Base.names(::Type{<:SimpleNamedVectorSpace{N}}) where N = (N,)
@inline pairtype(::Type{V}) where {V<:SimpleNamedVectorSpace} = NamedTuple{names(V), Tuple{eltype(V)}}
@inline Base.getindex(ps::NamedVectorSpacePairIteration{V}, i::Integer) where {V<:SimpleNamedVectorSpace} = NamedTuple{names(V)}((ps.namedvectorspace[i],))

"""
    ParameterSpace{N, T, B} <: SimpleNamedVectorSpace{N, B}

Parameter space.
"""
struct ParameterSpace{N, T, B} <: SimpleNamedVectorSpace{N, B}
    contents::T
    function ParameterSpace{N}(contents) where N
        @assert isa(N, Symbol) "ParameterSpace error: $(N) is not a Symbol."
        new{N, typeof(contents), eltype(contents)}(contents)
    end
end
@inline VectorSpaceStyle(::Type{<:ParameterSpace}) = VectorSpaceEnumerative()

"""
    CompositeNamedVectorSpace{T<:Tuple{Vararg{SimpleNamedVectorSpace}}, B<:Tuple} <: NamedVectorSpace{B}

Abstract composite named vector space.
"""
abstract type CompositeNamedVectorSpace{T<:Tuple{Vararg{SimpleNamedVectorSpace}}, B<:Tuple} <: NamedVectorSpace{B} end
@inline @generated Base.names(::Type{<:CompositeNamedVectorSpace{T}}) where {T<:Tuple{Vararg{SimpleNamedVectorSpace}}} = map(first, map(names, fieldtypes(T)))
@inline pairtype(::Type{V}) where {V<:CompositeNamedVectorSpace} = NamedTuple{names(V), eltype(V)}
@inline Tuple(basis::Tuple, ::CompositeNamedVectorSpace) = basis
@inline Base.getindex(ps::NamedVectorSpacePairIteration{V}, i::Integer) where {V<:CompositeNamedVectorSpace} = NamedTuple{names(V)}(ps.namedvectorspace[i])

"""
    ZippedNamedVectorSpace{T<:Tuple{Vararg{SimpleNamedVectorSpace}}, B<:Tuple} <: CompositeNamedVectorSpace{T, B}

Zipped named vector space.
"""
struct ZippedNamedVectorSpace{T<:Tuple{Vararg{SimpleNamedVectorSpace}}, B<:Tuple} <: CompositeNamedVectorSpace{T, B}
    contents::T
    ZippedNamedVectorSpace(contents::Tuple{Vararg{SimpleNamedVectorSpace}}) = new{typeof(contents), map(eltype, typeof(contents))}(contents)
end
@inline ZippedNamedVectorSpace(contents::SimpleNamedVectorSpace...) = ZippedNamedVectorSpace(contents)
@inline VectorSpaceStyle(::Type{<:ZippedNamedVectorSpace}) = VectorSpaceZipped()

"""
    DirectProductedNamedVectorSpace{T<:Tuple{Vararg{SimpleNamedVectorSpace}}, B<:Tuple} <: CompositeNamedVectorSpace{T, B}

Direct producted named vector space.
"""
struct DirectProductedNamedVectorSpace{T<:Tuple{Vararg{SimpleNamedVectorSpace}}, B<:Tuple} <: CompositeNamedVectorSpace{T, B}
    contents::T
    DirectProductedNamedVectorSpace(contents::Tuple{Vararg{SimpleNamedVectorSpace}}) = new{typeof(contents), map(eltype, typeof(contents))}(contents)
end
@inline DirectProductedNamedVectorSpace(contents::SimpleNamedVectorSpace...) = DirectProductedNamedVectorSpace(contents)
@inline VectorSpaceStyle(::Type{<:DirectProductedNamedVectorSpace}) = VectorSpaceDirectProducted()

"""
    ⊕(vs₁::SimpleNamedVectorSpace, vs₂::SimpleNamedVectorSpace) -> ZipNamedVectorSpace
    ⊕(vs₁::SimpleNamedVectorSpace, vs₂::ZippedNamedVectorSpace) -> ZipNamedVectorSpace
    ⊕(vs₁::ZippedNamedVectorSpace, vs₂::SimpleNamedVectorSpace) -> ZipNamedVectorSpace
    ⊕(vs₁::ZippedNamedVectorSpace, vs₂::ZippedNamedVectorSpace) -> ZipNamedVectorSpace

The zip of named vector spaces.
"""
@inline ⊕(vs₁::SimpleNamedVectorSpace, vs₂::SimpleNamedVectorSpace) = ZippedNamedVectorSpace(vs₁, vs₂)
@inline ⊕(vs₁::SimpleNamedVectorSpace, vs₂::ZippedNamedVectorSpace) = ZippedNamedVectorSpace(vs₁, vs₂.contents...)
@inline ⊕(vs₁::ZippedNamedVectorSpace, vs₂::SimpleNamedVectorSpace) = ZippedNamedVectorSpace(vs₁.contents..., vs₂)
@inline ⊕(vs₁::ZippedNamedVectorSpace, vs₂::ZippedNamedVectorSpace) = ZippedNamedVectorSpace(vs₁.contents..., vs₂.contents...)

"""
    ⊗(vs₁::SimpleNamedVectorSpace, vs₂::SimpleNamedVectorSpace) -> DirectProductedNamedVectorSpace
    ⊗(vs₁::SimpleNamedVectorSpace, vs₂::DirectProductedNamedVectorSpace) -> DirectProductedNamedVectorSpace
    ⊗(vs₁::DirectProductedNamedVectorSpace, vs₂::SimpleNamedVectorSpace) -> DirectProductedNamedVectorSpace
    ⊗(vs₁::DirectProductedNamedVectorSpace, vs₂::DirectProductedNamedVectorSpace) -> DirectProductedNamedVectorSpace

The direct product of named vector spaces.
"""
@inline ⊗(vs₁::SimpleNamedVectorSpace, vs₂::SimpleNamedVectorSpace) = DirectProductedNamedVectorSpace(vs₁, vs₂)
@inline ⊗(vs₁::SimpleNamedVectorSpace, vs₂::DirectProductedNamedVectorSpace) = DirectProductedNamedVectorSpace(vs₁, vs₂.contents...)
@inline ⊗(vs₁::DirectProductedNamedVectorSpace, vs₂::SimpleNamedVectorSpace) = DirectProductedNamedVectorSpace(vs₁.contents..., vs₂)
@inline ⊗(vs₁::DirectProductedNamedVectorSpace, vs₂::DirectProductedNamedVectorSpace) = DirectProductedNamedVectorSpace(vs₁.contents..., vs₂.contents...)

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
    length = segment.length+count(isequal(false), segment.ends)-1
    step = convert(eltype(segment), (segment.stop-segment.start)/length)
    start = segment.ends[1] ? segment.start : segment.start+step
    return start+(i-1)*step
end
@inline Base.getindex(segment::Segment, range::OrdinalRange{<:Integer}) = Segment(segment[first(range)], segment[last(range)], length(range), (true, true))
function Base.iterate(segment::Segment)
    segment.length==0 && return
    length = segment.length+count(isequal(false), segment.ends)-1
    step = convert(eltype(segment), (segment.stop-segment.start)/length)
    start = segment.ends[1] ? segment.start : segment.start+step
    return start, (1, start, step)
end
function Base.iterate(segment::Segment, state)
    i, middle, step = state
    i==segment.length && return
    middle = middle+step
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
    Segment(start::AbstractVector, stop::AbstractVector, length::Integer; ends::Tuple{Bool, Bool}=(true, false))
    Segment(start::NTuple{N, Number}, stop::NTuple{N, Number}, length::Integer; ends::Tuple{Bool, Bool}=(true, false)) where N

Construct a segment.
"""
function Segment(start::Number, stop::Number, length::Integer; ends::Tuple{Bool, Bool}=(true, false))
    @assert length>=0 "Segment error: length must be non-negative."
    dtype = promote_type(typeof(start), typeof(stop), Float64)
    return Segment(convert(dtype, start), convert(dtype, stop), length, ends)
end
function Segment(start::AbstractVector, stop::AbstractVector, length::Integer; ends::Tuple{Bool, Bool}=(true, false))
    @assert length>=0 "Segment error: length must be non-negative."
    @assert Base.length(start)==Base.length(stop) "Segment error: start and stop should have equal length."
    dtype = SVector{Base.length(start), promote_type(eltype(start), eltype(stop), Float64)}
    return Segment(convert(dtype, start), convert(dtype, stop), length, ends)
end
function Segment(start::NTuple{N, Number}, stop::NTuple{N, Number}, length::Integer; ends::Tuple{Bool, Bool}=(true, false)) where N
    @assert length>=0 "Segment error: length must be non-negative."
    dtype = SVector{N, promote_type(eltype(start), eltype(stop), Float64)}
    return Segment(convert(dtype, start), convert(dtype, stop), length, ends)
end

# Simple trees
"""
    SimpleTreeCore()

The core of a simple tree.
"""
mutable struct SimpleTreeCore{N, D}
    root::Union{N, Nothing}
    const contents::Dict{N, D}
    const parent::Dict{N, N}
    const children::Dict{N, Vector{N}}
    SimpleTreeCore{N, D}() where {N, D} = new{N, D}(nothing, Dict{N, D}(), Dict{N, N}(), Dict{N, Vector{N}}())
end
@inline Base.:(==)(tc1::TC, tc2::TC) where {TC<:SimpleTreeCore} = ==(efficientoperations, tc1, tc2)
@inline Base.isequal(tc1::TC, tc2::TC) where {TC<:SimpleTreeCore} = isequal(efficientoperations, tc1, tc2)

abstract type SimpleTreeIteration end
struct SimpleTreeDepth <: SimpleTreeIteration end
struct SimpleTreeWidth <: SimpleTreeIteration end

"""
    simpletreedepth

Indicate that the iteration over a tree is depth-first.
"""
const simpletreedepth = SimpleTreeDepth()

"""
    simpletreewidth

Indicate that the iteration over a tree is width-first.
"""
const simpletreewidth = SimpleTreeWidth()

"""
    AbstractSimpleTree{Node, Data}

Abstract type for all concrete trees.
"""
abstract type AbstractSimpleTree{N, D} end
@inline contentnames(::Type{<:AbstractSimpleTree}) = (:TREECORE,)
@inline Base.:(==)(t₁::T, t₂::T) where {T<:AbstractSimpleTree} = ==(efficientoperations, t₁, t₂)
@inline Base.isequal(t₁::T, t₂::T) where {T<:AbstractSimpleTree} = isequal(efficientoperations, t₁, t₂)

"""
    SimpleTreeCore(tree::AbstractSimpleTree) -> SimpleTreeCore

Get the core of a simple tree.
"""
@inline SimpleTreeCore(tree::AbstractSimpleTree) = getcontent(tree, :TREECORE)

"""
    keytype(tree::AbstractSimpleTree)
    keytype(::Type{T}) where {T<:AbstractSimpleTree}

Get a tree's node type.
"""
@inline Base.keytype(tree::AbstractSimpleTree) = keytype(typeof(tree))
@inline @generated Base.keytype(::Type{T}) where {T<:AbstractSimpleTree} = parametertype(supertype(T, :AbstractSimpleTree), 1)

"""
    valtype(tree::AbstractSimpleTree)
    valtype(::Type{T}) where {T<:AbstractSimpleTree}

Get a tree's data type.
"""
@inline Base.valtype(tree::AbstractSimpleTree) = valtype(typeof(tree))
@inline @generated Base.valtype(::Type{T}) where {T<:AbstractSimpleTree} = parametertype(supertype(T, :AbstractSimpleTree), 2)

"""
    eltype(tree::AbstractSimpleTree)
    eltype(::Type{T}) where {T<:AbstractSimpleTree}

Get the eltype of a tree.
"""
@inline Base.eltype(tree::AbstractSimpleTree) = eltype(typeof(tree))
@inline Base.eltype(::Type{T}) where {T<:AbstractSimpleTree} = Pair{keytype(T), valtype(T)}

"""
    root(tree::AbstractSimpleTree) -> Union{keytype(tree), Nothing}

Get a tree's root node.
"""
@inline root(tree::AbstractSimpleTree) = SimpleTreeCore(tree).root

"""
    haskey(tree::AbstractSimpleTree{N}, node::N) where N -> Bool

Check whether a node is in a tree.
"""
@inline Base.haskey(tree::AbstractSimpleTree{N}, node::N) where N = haskey(SimpleTreeCore(tree).contents, node)

"""
    length(tree::AbstractSimpleTree) -> Int

Get the number of a tree's nodes.
"""
@inline Base.length(tree::AbstractSimpleTree) = length(SimpleTreeCore(tree).contents)

"""
    parent(tree::AbstractSimpleTree{N}, node::N, superparent::Union{N, Nothing}=nothing) where N -> Union{N, Nothing}

Get the parent of a tree's node. When `node` is the tree's root, return `superparent`.
"""
@inline Base.parent(tree::AbstractSimpleTree{N}, node::N, superparent::Union{N, Nothing}=nothing) where {N} = (node == root(tree) ? superparent : SimpleTreeCore(tree).parent[node])

"""
    children(tree::AbstractSimpleTree) -> Vector{keytype(tree)}
    children(tree::AbstractSimpleTree, ::Nothing) -> Vector{keytype(tree)}
    children(tree::AbstractSimpleTree{N}, node::N) where N -> Vector{N}

Get the children of a tree's node.
"""
@inline children(tree::AbstractSimpleTree) = children(tree, nothing)
@inline children(tree::AbstractSimpleTree, ::Nothing) = (root(tree) === nothing) ? error("children error: empty tree!") : [root(tree)]
@inline children(tree::AbstractSimpleTree{N}, node::N) where N = SimpleTreeCore(tree).children[node]

"""
    addnode!(tree::AbstractSimpleTree{N}, node::N) where N} -> typeof(tree)
    addnode!(tree::AbstractSimpleTree{N}, ::Nothing, node::N) where N -> typeof(tree)
    addnode!(tree::AbstractSimpleTree{N}, parent::N, node::N) where N -> typeof(tree)

Update the structure of a tree by adding a node. When the parent is `nothing`, the input tree must be empty and the input node becomes the tree's root.
"""
@inline addnode!(tree::AbstractSimpleTree{N}, node::N) where {N} = addnode!(tree, nothing, node)
function addnode!(tree::AbstractSimpleTree{N}, ::Nothing, node::N) where N
    @assert root(tree) === nothing "addnode! error: not empty tree."
    SimpleTreeCore(tree).root = node
    SimpleTreeCore(tree).children[node] = N[]
    return tree
end
function addnode!(tree::AbstractSimpleTree{N}, parent::N, node::N) where N
    @assert haskey(tree, parent) "addnode! error: parent($parent) not in tree."
    @assert !haskey(tree, node) "addnode! error: node($node) already in tree."
    push!(SimpleTreeCore(tree).children[parent], node)
    SimpleTreeCore(tree).parent[node] = parent
    SimpleTreeCore(tree).children[node] = N[]
    return tree
end

"""
    deletenode!(tree::AbstractSimpleTree{N}, node::N) where N -> typeof(tree)

Update the structure of a tree by deleting a node.
"""
function deletenode!(tree::AbstractSimpleTree{N}, node::N) where N
    @assert haskey(tree, node) "deletenode! error: node($node) not in tree."
    if node == root(tree)
        SimpleTreeCore(tree).root = nothing
    else
        pnode = parent(tree, node)
        haskey(SimpleTreeCore(tree).children, pnode) && filter!(!=(node), SimpleTreeCore(tree).children[pnode])
    end
    pop!(SimpleTreeCore(tree).contents, node)
    pop!(SimpleTreeCore(tree).parent, node, nothing)
    pop!(SimpleTreeCore(tree).children, node)
    return tree
end

"""
    getindex(tree::AbstractSimpleTree{N}, node::N) where N -> N

Get the data of a tree's node.
"""
@inline Base.getindex(tree::AbstractSimpleTree{N}, node::N) where {N} = SimpleTreeCore(tree).contents[node]

"""
    setindex!(tree::AbstractSimpleTree{N, D}, data::D, node::N) where {N, D}

Set the data of a tree's node.
"""
@inline Base.setindex!(tree::AbstractSimpleTree{N, D}, data::D, node::N) where {N, D} = (SimpleTreeCore(tree).contents[node] = data)

"""
    empty(tree::AbstractSimpleTree)

Construct an empty tree of the same type with the input one.
"""
@inline Base.empty(tree::AbstractSimpleTree) = rawtype(typeof(tree))(dissolve(tree, empty)...)
@inline dissolve(tree::AbstractSimpleTree{N, D}, ::Val{:TREECORE}, ::typeof(empty), args::Tuple, kwargs::NamedTuple) where {N, D} = SimpleTreeCore{N, D}()

"""
    keys(tree::AbstractSimpleTree{N}, ::SimpleTreeDepth, node::Union{N, Nothing}=root(tree)) where N
    keys(tree::AbstractSimpleTree{N}, ::SimpleTreeWidth, node::Union{N, Nothing}=root(tree)) where N

Iterate over a tree's nodes starting from a certain `node` by depth first search or width first search.
"""
function Base.keys(tree::AbstractSimpleTree{N}, ti::SimpleTreeIteration, node::Union{N, Nothing}=root(tree)) where N
    TreeKeys{typeof(ti), N, valtype(tree), typeof(tree)}(tree, node)
end
struct TreeKeys{P<:SimpleTreeIteration, N, D, T<:AbstractSimpleTree{N, D}}
    tree::T
    node::Union{N, Nothing}
end
Base.eltype(::Type{<:TreeKeys{<:SimpleTreeIteration, N, D, <:AbstractSimpleTree{N, D}} where D}) where {N} = N
Base.IteratorSize(::Type{<:TreeKeys}) = Base.SizeUnknown()
Base.iterate(tk::TreeKeys) = (root(tk.tree) === nothing) ? nothing : (tk.node, copy(children(tk.tree, tk.node)))
function Base.iterate(tk::TreeKeys{SimpleTreeDepth, N}, state::Vector{N}) where N
    (length(state) == 0) ? nothing : (node = popfirst!(state); prepend!(state, children(tk.tree, node)); (node, state))
end
function Base.iterate(tk::TreeKeys{SimpleTreeWidth, N}, state::Vector{N}) where N
    (length(state) == 0) ? nothing : (node = popfirst!(state); append!(state, children(tk.tree, node)); (node, state))
end

"""
    values(tree::AbstractSimpleTree{N}, ::SimpleTreeDepth, node::Union{N, Nothing}=root(tree)) where N
    values(tree::AbstractSimpleTree{N}, ::SimpleTreeWidth, node::Union{N, Nothing}=root(tree)) where N

Iterate over a tree's data starting from a certain `node` by depth first search or width first search.
"""
@inline Base.values(tree::AbstractSimpleTree{N}, ti::SimpleTreeIteration, node::Union{N, Nothing}=root(tree)) where {N} = (tree[key] for key in keys(tree, ti, node))

"""
    pairs(tree::AbstractSimpleTree{N}, ::SimpleTreeDepth, node::Union{N, Nothing}=root(tree)) where N
    pairs(tree::AbstractSimpleTree{N}, ::SimpleTreeWidth, node::Union{N, Nothing}=root(tree)) where N

Iterate over a tree's (node, data) pairs starting from a certain `node` by depth first search or width first search.
"""
@inline Base.pairs(tree::AbstractSimpleTree{N}, ti::SimpleTreeIteration, node::Union{N, Nothing}=root(tree)) where {N} = ((key, tree[key]) for key in keys(tree, ti, node))

"""
    isleaf(tree::AbstractSimpleTree{N}, node::N) where N -> Bool

Judge whether a tree's node is a leaf (a node without children) or not.
"""
@inline isleaf(tree::AbstractSimpleTree{N}, node::N) where {N} = length(children(tree, node)) == 0

"""
    level(tree::AbstractSimpleTree{N}, node::N) where N -> Int

Get the level of tree's node.
"""
@inline function level(tree::AbstractSimpleTree{N}, node::N) where N
    result = 1
    while node != root(tree)
        result += 1
        node = parent(tree, node)
    end
    return result
end

"""
    ancestor(tree::AbstractSimpleTree{N}, node::N, generation::Int=1) where N -> N

Get the ancestor of a tree's node of the n-th generation.
"""
@inline function ancestor(tree::AbstractSimpleTree{N}, node::N, generation::Int=1) where N
    @assert generation >= 0 "ancestor error: generation($generation) must be non-negative."
    result = node
    for i = 1:generation
        result = parent(tree, result)
    end
    result
end

"""
    descendants(tree::AbstractSimpleTree{N}, node::N, generation::Int=1) where N -> Vector{N}

Get the descendants of a tree's node of the nth generation.
"""
function descendants(tree::AbstractSimpleTree{N}, node::N, generation::Int=1) where N
    @assert generation >= 0 "descendants error: generation($generation) must be non-negative."
    result = N[node]
    for i = 1:generation
        result = vcat((children(tree, node) for node in result)...)
    end
    result
end

"""
    siblings(tree::AbstractSimpleTree{N}, node::N) where N -> Vector{N}

Get the siblings (other nodes sharing the same parent) of a tree's node.
"""
@inline siblings(tree::AbstractSimpleTree{N}, node::N) where N = filter(!=(node), children(tree, parent(tree, node)))

"""
    leaves(tree::AbstractSimpleTree) -> Vector{keytype(tree)}

Get a tree's leaves.
"""
@inline leaves(tree::AbstractSimpleTree) = keytype(tree)[node for node in keys(tree, simpletreedepth) if isleaf(tree, node)]

"""
    push!(tree::AbstractSimpleTree{N, D}, node::N, data::D) where {N, D} -> typeof(tree)
    push!(tree::AbstractSimpleTree{N, D}, parent::Union{N, Nothing}, node::N, data::D) where {N, D} -> typeof(tree)

Push a new node to a tree. When `parent` is `nothing`, this function set the root node of an empty tree.
"""
Base.push!(tree::AbstractSimpleTree{N, D}, node::N, data::D) where {N, D} = push!(tree, nothing, node, data)
function Base.push!(tree::AbstractSimpleTree{N, D}, parent::Union{N, Nothing}, node::N, data::D) where {N, D}
    addnode!(tree, parent, node)
    tree[node] = data
    return tree
end

"""
    append!(tree::AbstractSimpleTree{N, D}, subtree::AbstractSimpleTree{N, D}) where {N, D} -> typeof(tree)
    append!(tree::AbstractSimpleTree{N, D}, node::Union{N, Nothing}, subtree::AbstractSimpleTree{N, D}) where {N, D} -> typeof(tree)

Append a subtree to a tree.
"""
Base.append!(tree::AbstractSimpleTree{N, D}, subtree::AbstractSimpleTree{N, D}) where {N, D} = append!(tree, nothing, subtree)
function Base.append!(tree::AbstractSimpleTree{N, D}, node::Union{N, Nothing}, subtree::AbstractSimpleTree{N, D}) where {N, D}
    for (key, value) in pairs(subtree, simpletreewidth)
        push!(tree, parent(subtree, key, node), key, value)
    end
    return tree
end

"""
    delete!(tree::AbstractSimpleTree{N}, node::N) where N -> typeof(tree)

Delete a node and all its descendants from a tree.
"""
function Base.delete!(tree::AbstractSimpleTree{N}, node::N) where N
    for key in collect(N, keys(tree, simpletreedepth, node))
        deletenode!(tree, key)
    end
    return tree
end

"""
    empty!(tree::AbstractSimpleTree) -> typeof(tree)

Empty a tree.
"""
@inline Base.empty!(tree::AbstractSimpleTree) = delete!(tree, root(tree))

"""
    subtree(tree::AbstractSimpleTree{N}, node::N) where N -> typeof(tree)

Get a subtree whose root is `node`.
"""
function subtree(tree::AbstractSimpleTree{N}, node::N) where N
    result = empty(tree)
    for (i, (key, value)) in enumerate(pairs(tree, simpletreedepth, node))
        push!(result, (i == 1) ? nothing : parent(tree, key), key, value)
    end
    return result
end

"""
    move!(tree::AbstractSimpleTree{N}, node::N, parent::N) where N -> typeof(tree)

Move a subtree to a new position.
"""
function move!(tree::AbstractSimpleTree{N}, node::N, parent::N) where N
    sub = subtree(tree, node)
    delete!(tree, node)
    append!(tree, parent, sub)
    return tree
end

"""
    SimpleTree{N, D} <: AbstractSimpleTree{N, D}

The minimum tree structure that implements all the default tree methods.
"""
struct SimpleTree{N, D} <: AbstractSimpleTree{N, D}
    TREECORE::SimpleTreeCore{N, D}
end

"""
    SimpleTree{N, D}() where {N, D}

Construct an empty simple tree.
"""
SimpleTree{N, D}() where {N, D} = SimpleTree(SimpleTreeCore{N, D}())

end # module
