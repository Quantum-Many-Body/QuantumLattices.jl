module QuantumOperators

using DataStructures: OrderedDict
using Latexify: unicode2latex
using Printf: @printf, @sprintf
using ..Toolkit: atol, efficientoperations, rtol, contentorder, decimaltostr, fulltype, getcontent, parameterpairs, parametertype, promoteparameters, rawtype, reparameter

import LaTeXStrings: latexstring
import LinearAlgebra: dot
import ..QuantumLattices: ⊗, add!, div!, dtype, id, ishermitian, mul!, permute, rank, sub!, value
import ..Toolkit: contentnames, dissolve, isparameterbound, parameternames, subscript, superscript

export ID, LaTeX, Operator, OperatorPack, OperatorProd, Operators, OperatorSet, OperatorSum, OperatorUnit, QuantumOperator, Representation
export LinearFunction, LinearTransformation, Matrixization, Permutation, RankFilter, TabledUnitSubstitution, Transformation, UnitSubstitution
export idtype, ishermitian, isscalartype, latexname, latexformat, matrix, optype, script, sequence

# Generic representation of (a part of) a quantum lattice system
"""
    Representation

The abstract type of a representation of (a part of) a quantum lattice system.
"""
abstract type Representation end
@inline Base.:(==)(r₁::Representation, r₂::Representation) = ==(efficientoperations, r₁, r₂)
@inline Base.isequal(r₁::Representation, r₂::Representation) = isequal(efficientoperations, r₁, r₂)

# Generic quantum operator
"""
    QuantumOperator <: Representation

The abstract type of any quantum operator.
"""
abstract type QuantumOperator <: Representation end
@inline Base.zero(m::QuantumOperator) = zero(typeof(m))
@inline dot(m₁::QuantumOperator, m₂::QuantumOperator) = conj(m₁) * m₂
@inline dot(m::QuantumOperator, c::Number) = conj(m) * c
@inline dot(c::Number, m::QuantumOperator) = conj(c) * m
@inline ⊗(m::QuantumOperator, c::Number) = m * c
@inline ⊗(c::Number, m::QuantumOperator) = c * m

"""
    replace(m::QuantumOperator; kwargs...) -> typeof(m)

Return a copy of a concrete `QuantumOperator` with some of the field values replaced by the keyword arguments.
"""
@inline Base.replace(m::QuantumOperator; kwargs...) = replace(efficientoperations, m; kwargs...)

# Operator unit
"""
    OperatorUnit <: QuantumOperator

An operator unit is the irreducible symbolic unit to completely represent a quantum operator.

It plays the role of the symbols as in usual computer algebras while it can host internal structures, which is convenient for quantum operators in representative of the internal degrees of freedom.
"""
abstract type OperatorUnit <: QuantumOperator end
@inline Base.show(io::IO, u::OperatorUnit) = @printf io "%s(%s)" nameof(typeof(u)) join(map(repr, ntuple(i->getfield(u, i), Val(fieldcount(typeof(u))))), ", ")
@inline @generated Base.hash(u::OperatorUnit, h::UInt) = Expr(:call, :hash, Expr(:tuple, [:(getfield(u, $i)) for i=1:fieldcount(u)]...), :h)

# ID of a composite quantum operator
"""
    ID{U<:OperatorUnit, N}

The id of a composite quantum operator, which is an ordered set of operator units.

Type alias for `NTuple{N, U} where {U<:OperatorUnit}`.
"""
const ID{U<:OperatorUnit, N} = NTuple{N, U}
@inline Base.promote_rule(::Type{Tuple{}}, I::Type{<:Tuple{OperatorUnit, Vararg{OperatorUnit}}}) = ID{I|>eltype}
@inline Base.promote_rule(I::Type{<:Tuple{OperatorUnit, Vararg{OperatorUnit}}}, ::Type{Tuple{}}) = ID{I|>eltype}

"""
    ID(id::OperatorUnit...)
    ID(u::OperatorUnit, id::ID{OperatorUnit})
    ID(id::ID{OperatorUnit}, u::OperatorUnit)
    ID(id₁::ID{OperatorUnit}, id₂::ID{OperatorUnit})

Get the id from operator units/ids.
"""
@inline ID(id::OperatorUnit...) = id
@inline ID(u::OperatorUnit, id::ID{OperatorUnit}) = ID(u, id...)
@inline ID(id::ID{OperatorUnit}, u::OperatorUnit) = ID(id..., u)
@inline ID(id₁::ID{OperatorUnit}, id₂::ID{OperatorUnit}) = ID(id₁..., id₂...)

"""
    ID(::Type{U}, attrs::Vararg{NTuple{N}, M}) where {U<:OperatorUnit, N, M}

Get the composite id from the components of singular ids.
"""
@inline @generated function ID(::Type{U}, attrs::Vararg{NTuple{N, Any}, M}) where {U<:OperatorUnit, N, M}
    exprs = []
    for i = 1:N
        args = [:(attrs[$j][$i]) for j = 1:M]
        push!(exprs, :(U($(args...))))
    end
    return :(ID($(exprs...)))
end

"""
    propertynames(::Type{I}) where I<:ID{OperatorUnit} -> Tuple{Vararg{Symbol}}

Get the property names of a composite id.
"""
@inline @generated function Base.propertynames(I::ID{OperatorUnit})
    exprs = [QuoteNode(Symbol(name, 's')) for name in I|>eltype|>fieldnames]
    return Expr(:tuple, exprs...)
end

"""
    getproperty(id::ID{OperatorUnit}, name::Symbol)

Get the property of a composite id.
"""
@inline Base.getproperty(id::ID{OperatorUnit}, name::Symbol) = idgetproperty(id, Val(name), Val(id|>propertynames))
@inline @generated function idgetproperty(id::ID{OperatorUnit}, ::Val{name}, ::Val{names}) where {name, names}
    index = findfirst(isequal(name), names)::Int
    exprs = [:(getfield(id[$i], $index)) for i = 1:fieldcount(id)]
    return Expr(:tuple, exprs...)
end

"""
    rank(id::ID{OperatorUnit}) -> Int
    rank(::Type{<:ID{OperatorUnit, N}}) where N -> Int

Get the rank of an id.
"""
@inline rank(id::ID{OperatorUnit}) = rank(typeof(id))
@inline rank(::Type{<:ID{OperatorUnit, N}}) where N = N

"""
    adjoint(id::ID{OperatorUnit}) -> ID

Get the adjoint of an id.
"""
@inline Base.adjoint(id::ID{OperatorUnit}) = map(adjoint, reverse(id))

"""
    ishermitian(id::ID{OperatorUnit}) -> Bool

Judge whether an id is Hermitian.
"""
function ishermitian(id::ID{OperatorUnit})
    for i = 1:((rank(id)+1)÷2)
        id[i]'==id[rank(id)+1-i] || return false
    end
    return true
end

# Operator pack
"""
    OperatorPack{V, I} <: QuantumOperator

The entity that represent the pack of a number and an id of a quantum operator.

Basically, a concrete subtype should contain two predefined contents:
- `value::V`: the coefficient of the pack
- `id::I`: the id of the pack
"""
abstract type OperatorPack{V, I} <: QuantumOperator end
@inline contentnames(::Type{<:OperatorPack}) = (:value, :id)
@inline parameternames(::Type{<:OperatorPack}) = (:value, :id)
@inline isparameterbound(::Type{<:OperatorPack}, ::Val{:value}, ::Type{V}) where V = false
@inline isscalartype(::Type{<:OperatorPack}) = false
@inline function Base.promote_rule(::Type{M₁}, ::Type{M₂}) where {M₁<:OperatorPack, M₂<:OperatorPack}
    M₁<:M₂ && return M₂
    M₂<:M₁ && return M₁
    s₁, s₂ = isscalartype(M₁), isscalartype(M₂)
    M = s₂ ? rawtype(M₁) : s₁ ? rawtype(M₂) : typejoin(rawtype(M₁), rawtype(M₂))
    return fulltype(M, promoteparameters(parameterpairs(M₁), parameterpairs(M₂)))
end

"""
    valtype(m::OperatorPack)
    valtype(::Type{T}) where {T<:OperatorPack}

Get the type of the value of an `OperatorPack`.
"""
@inline Base.valtype(m::OperatorPack) = valtype(typeof(m))
@inline @generated Base.valtype(::Type{T}) where {T<:OperatorPack} = parametertype(supertype(T, :OperatorPack), :value)

"""
    idtype(m::OperatorPack)
    idtype(::Type{T}) where {T<:OperatorPack}

The type of the id of an `OperatorPack`.
"""
@inline idtype(m::OperatorPack) = idtype(typeof(m))
@inline @generated idtype(::Type{T}) where {T<:OperatorPack} = parametertype(supertype(T, :OperatorPack), :id)

"""
    dtype(m::OperatorPack)
    dtype(::Type{T}) where {T<:OperatorPack}

The data type of the coefficient of an `OperatorPack`.
"""
@inline dtype(m::OperatorPack) = dtype(typeof(m))
@inline dtype(::Type{T}) where {T<:OperatorPack} = valtype(T)

"""
    value(m::OperatorPack) -> valtype(m)

Get the value of an `OperatorPack`.
"""
@inline value(m::OperatorPack) = getcontent(m, :value)

"""
    id(m::OperatorPack) -> idtype(m)

Get the id of an `OperatorPack`.
"""
@inline id(m::OperatorPack) = getcontent(m, :id)

"""
    replace(m::OperatorPack, v) -> OperatorPack

Replace the value of an `OperatorPack`.
"""
@inline Base.replace(m::OperatorPack, v) = rawtype(typeof(m))(dissolve(m, replace, v)...)
@inline dissolve(::OperatorPack, ::Val{:value}, ::typeof(replace), v) = v

"""
    isapprox(m₁::OperatorPack, m₂::OperatorPack; atol::Real=atol, rtol::Real=rtol) -> Bool

Compare two `OperatorPack`s and judge whether they are approximate to each other.
"""
@inline function Base.isapprox(m₁::OperatorPack, m₂::OperatorPack; atol::Real=atol, rtol::Real=rtol)
    return isapprox(efficientoperations, Val(contentorder(typeof(m₁), :value)), dissolve(m₁), dissolve(m₂); atol=atol, rtol=rtol)::Bool
end

"""
    convert(::Type{M}, m::OperatorPack) where {M<:OperatorPack}

Convert a quantum operator from one type to another.
"""
@inline function Base.convert(::Type{M}, m::OperatorPack) where {M<:OperatorPack}
    typeof(m)<:M && return m
    @assert convertible(M, m) "convert error: $(typeof(m)) cannot be converted to $M."
    return replace(m, convert(valtype(M), value(m)))
end
function convertible(::Type{M}, m::OperatorPack) where {M<:OperatorPack}
    !(rawtype(typeof(m)) <: rawtype(M)) && return false
    S = supertype(typeof(m), nameof(M))
    for (i, name) in enumerate(parameternames(M))
        (name != :value) && !(parametertype(S, i) <: parametertype(M, i)) && return false
    end
    return true
end

# Operator prod
"""
    OperatorProd{V, I<:Tuple} <: OperatorPack{V, I}

A special kind of `OperatorPack`, where the relation between the coefficient and each component of the opeator id can be viewed as product.
"""
abstract type OperatorProd{V, I<:Tuple} <: OperatorPack{V, I} end
@inline Base.eltype(m::OperatorProd) = eltype(typeof(m))
@inline Base.eltype(::Type{M}) where {M<:OperatorProd} = eltype(idtype(M))
@inline Base.iterate(m::OperatorProd) = iterate(id(m))
@inline Base.iterate(m::OperatorProd, state) = iterate(id(m), state)
@inline Base.firstindex(m::OperatorProd) = 1
@inline Base.lastindex(m::OperatorProd) = rank(m)
@inline Base.promote_rule(::Type{M}, ::Type{N}) where {M<:OperatorProd, N<:Number} = reparameter(M, :value, promote_type(valtype(M), N))
@inline isscalartype(::Type{<:OperatorProd{V, Tuple{}} where V}) = true

"""
    rank(m::OperatorProd) -> Int
    rank(::Type{M}) where {M<:OperatorProd} -> Int

Get the rank of an `OperatorProd`.
"""
@inline rank(m::OperatorProd) = rank(typeof(m))
@inline rank(::Type{M}) where {M<:OperatorProd} = fieldcount(idtype(M))

"""
    length(m::OperatorProd) -> Int

Get the length of an `OperatorProd`.
"""
@inline Base.length(m::OperatorProd) = rank(m)

"""
    getindex(m::OperatorProd, i::Integer) -> OperatorUnit
    getindex(m::OperatorProd, slice) -> OperatorProd

Overloaded `[]`.
"""
@inline Base.getindex(m::OperatorProd, i::Integer) = id(m)[i]
@inline Base.getindex(m::OperatorProd, slice) = rawtype(typeof(m))(dissolve(m, getindex, slice)...)
@inline dissolve(m::OperatorProd, ::Val{:value}, ::typeof(getindex), ::Any) = one(valtype(m))
@inline dissolve(m::OperatorProd, ::Val{:id}, ::typeof(getindex), slice) = id(m)[slice]

"""
    split(m::OperatorProd) -> Tuple{valtype(m), Vararg{Any}}

Split an `OperatorProd` into the coefficient and a sequence of the components of its id.
"""
@inline Base.split(m::OperatorProd) = (value(m), id(m)...)

"""
    one(::Type{M}) where {M<:OperatorProd}
    one(m::OperatorProd)

Get the identity quantum operator.
"""
@inline function Base.one(::Type{M}) where {M<:OperatorProd}
    rtype = rawtype(M)
    vtype = isconcretetype(valtype(M)) ? valtype(M) : Int
    @assert fieldnames(rtype) == (:value, :id) "one error: not supported type($(nameof(rtype)))."
    return rtype(one(vtype), ())
end
@inline Base.one(m::OperatorProd) = rawtype(typeof(m))(dissolve(m, one)...)
@inline dissolve(m::OperatorProd, ::Val{:value}, ::typeof(one)) = one(valtype(m))
@inline dissolve(::OperatorProd, ::Val{:id}, ::typeof(one)) = ()

"""
    convert(::Type{M}, m::Number) where {M<:OperatorProd}

Convert a number to a quantum operator.
"""
@inline function Base.convert(::Type{M}, m::Number) where {M<:OperatorProd}
    @assert isa(one(M), M) "convert error: not convertible."
    return one(M)*convert(dtype(M), m)
end

"""
    sequence(m::OperatorProd, table) -> NTuple{rank(m), Int}

Get the sequence of the id of a quantum operator according to a table.
"""
@inline sequence(m::OperatorProd, table) = map(u->table[u], id(m))

# Operator
"""
    Operator{V<:Number, I<:ID{OperatorUnit}} <: OperatorProd{V, I}

Operator.
"""
struct Operator{V<:Number, I<:ID{OperatorUnit}} <: OperatorProd{V, I}
    value::V
    id::I
end
@inline Operator(value::Number, id::OperatorUnit...) = Operator(value, id)
function Base.show(io::IO, m::Operator)
    @printf io "%s(%s%s%s)" nameof(typeof(m)) decimaltostr(value(m)) (rank(m)>0 ? ", " : "") join(id(m), ", ")
end

"""
    adjoint(m::Operator) -> Operator

Get the adjoint of an operator.
"""
@inline Base.adjoint(m::Operator) = rawtype(typeof(m))(value(m)', id(m)')

"""
    ishermitian(m::Operator) -> Bool

Judge whether an operator is Hermitian.
"""
@inline ishermitian(m::Operator) = isreal(value(m)) && ishermitian(id(m))

"""
    convert(::Type{M}, u::OperatorUnit) where {M<:Operator}

Convert an operator unit to an operator.
"""
@inline function Base.convert(::Type{M}, u::OperatorUnit) where {M<:Operator}
    @assert Tuple{typeof(u)} <: idtype(M) "convert error: not convertible."
    return Operator(one(valtype(M)), ID(u))
end

# Operator set
"""
    OperatorSet{M<:OperatorPack} <: QuantumOperator

The set of `OperatorPack`s.

1) The relation between two `OperatorPack`s in an `OperatorSet` can be viewed as addition.
2) But in general, only iteration over `OperatorPack`s and length are supported.
3) To use arithmetic operations, please refer to its subtype, `OperatorSum`.
"""
abstract type OperatorSet{M<:OperatorPack} <: QuantumOperator end
@inline Base.eltype(ms::OperatorSet) = eltype(typeof(ms))
@inline Base.eltype(::Type{<:OperatorSet{M}}) where {M<:OperatorPack} = M
function Base.show(io::IO, ms::OperatorSet)
    @printf io "%s with %s %s\n" summary(ms) length(ms) nameof(eltype(ms))
    for m in ms
        @printf io "  %s\n" m
    end
end

# Operator sum
"""
    OperatorSum{M<:OperatorPack, I} <: OperatorSet{M}

The sum of `OperatorPack`s.

Similar items are automatically merged with the aid of the id system.
"""
struct OperatorSum{M<:OperatorPack, I} <: OperatorSet{M}
    contents::OrderedDict{I, M}
    OperatorSum(contents::OrderedDict{I, M}) where {I, M<:OperatorPack} = new{M, I}(contents)
end
@inline Base.iterate(ms::OperatorSum) = iterate(values(ms.contents))
@inline Base.iterate(ms::OperatorSum, state) = iterate(values(ms.contents), state)
@inline Base.firstindex(::OperatorSum) = 1
@inline Base.lastindex(ms::OperatorSum) = length(ms)
@inline function Base.promote_rule(::Type{MS₁}, ::Type{MS₂}) where {MS₁<:OperatorSum, MS₂<:OperatorSum}
    M = promote_type(eltype(MS₁), eltype(MS₂))
    return OperatorSum{M, idtype(M)}
end
@inline function Base.convert(::Type{MS}, ms::OperatorSum) where {MS<:OperatorSum}
    @assert eltype(MS)>:eltype(ms) "convert error: cannot convert an object of $(typeof(ms)) to an object of $(MS)."
    return add!(zero(MS), ms)
end

"""
    OperatorSum(ms)
    OperatorSum(ms::QuantumOperator...)
    OperatorSum{M}(ms) where {M<:OperatorPack}
    OperatorSum{M}(ms::QuantumOperator...) where {M<:OperatorPack}

Get the sum of `OperatorPack`s.
"""
@inline OperatorSum(ms::QuantumOperator...) = OperatorSum{eltype(ms)}(ms)
@inline OperatorSum(ms) = OperatorSum{eltype(ms)}(ms)
@inline OperatorSum{M}(ms::QuantumOperator...) where {M<:OperatorPack} = OperatorSum{M}(ms)
function OperatorSum{M}(ms) where {M<:OperatorPack}
    result = OperatorSum(OrderedDict{idtype(M), M}())
    for m in ms
        value(m)≠0 && add!(result, m)
    end
    return result
end

"""
    length(ms::OperatorSum) -> Int

Get the number of `OperatorPack`s contained in an `OperatorSum`.
"""
@inline Base.length(ms::OperatorSum) = length(ms.contents)

"""
    haskey(ms::OperatorSum, id) -> Bool

Judge whether an `OperatorSum` contains an `OperatorPack` with the given `id`. 
"""
@inline Base.haskey(ms::OperatorSum, id) = haskey(ms.contents, id)

"""
    getindex(ms::OperatorSum, index::Integer) -> eltype(ms)
    getindex(ms::OperatorSum, indexes::AbstractVector{<:Integer}) -> typeof(ms)
    getindex(ms::OperatorSum, ::Colon) -> typeof(ms)

Overloaded `[]`.
"""
@inline Base.getindex(ms::OperatorSum, index::Integer) = value(ms.contents, index)
@inline Base.getindex(ms::OperatorSum, indexes::AbstractVector{<:Integer}) = OperatorSum{eltype(ms)}(ms[index] for index in indexes)
@inline Base.getindex(ms::OperatorSum, ::Colon) = OperatorSum(copy(ms.contents))

"""
    empty(ms::OperatorSum) -> typeof(ms)
    empty!(ms::OperatorSum) -> typeof(ms)

Get an empty copy or empty an `OperatorSum`.
"""
@inline Base.empty(ms::OperatorSum) = OperatorSum(empty(ms.contents))
@inline Base.empty!(ms::OperatorSum) = (empty!(ms.contents); ms)

"""
    isapprox(ms₁::OperatorSum, ms₂::OperatorSum; atol::Real=atol, rtol::Real=rtol) -> Bool

Compare two `OperatorSum`s and judge whether they are approximate to each other.
"""
@inline function Base.isapprox(ms₁::OperatorSum, ms₂::OperatorSum; atol::Real=atol, rtol::Real=rtol)
    for m in ms₁
        isapprox(value(m), 0; atol=atol, rtol=rtol) && continue
        k = id(m)
        haskey(ms₂, k) || return false
        isapprox(m, ms₂.contents[k]; atol=atol, rtol=rtol) || return false
    end
    for m in ms₂
        isapprox(value(m), 0; atol=atol, rtol=rtol) && continue
        k = id(m)
        haskey(ms₁, k) || return false
        isapprox(m, ms₁.contents[k]; atol=atol, rtol=rtol) || return false
    end
    return true
end

"""
    add!(ms::OperatorSum) -> typeof(ms)
    add!(ms::OperatorSum, m::Union{Number, OperatorUnit, OperatorPack}) -> typeof(ms)
    add!(ms::OperatorSum, mms::OperatorSum) -> typeof(ms)

Get the in-place addition of quantum operators.
"""
@inline add!(ms::OperatorSum) = ms
@inline function add!(ms::OperatorSum, m::Union{Number, OperatorUnit, OperatorPack})
    isa(m, Number) && m==0 && return ms
    m = convert(eltype(ms), m)
    old = get(ms.contents, id(m), nothing)
    new = isnothing(old) ? m : replace(old, value(old)+value(m))
    value(new)==0 ? delete!(ms.contents, id(m)) : (ms.contents[id(m)] = new)
    return ms
end
@inline function add!(ms::OperatorSum, mms::OperatorSum)
    for m in mms
        add!(ms, m)
    end
    return ms
end

"""
    sub!(ms::OperatorSum) -> typeof(ms)
    sub!(ms::OperatorSum, m::Union{Number, OperatorUnit, OperatorPack}) -> typeof(ms)
    sub!(ms::OperatorSum, mms::OperatorSum) -> typeof(ms)

Get the in-place subtraction of quantum operators.
"""
@inline sub!(ms::OperatorSum) = ms
@inline function sub!(ms::OperatorSum, m::Union{Number, OperatorUnit, OperatorPack})
    isa(m, Number) && abs(m)==0 && return ms
    m = convert(eltype(ms), m)
    old = get(ms.contents, id(m), nothing)
    new = isnothing(old) ? -m : replace(old, value(old)-value(m))
    value(new)==0 ? delete!(ms.contents, id(m)) : (ms.contents[id(m)] = new)
    return ms
end
@inline function sub!(ms::OperatorSum, mms::OperatorSum)
    for m in mms
        sub!(ms, m)
    end
    return ms
end

"""
    mul!(ms::OperatorSum, factor::Number) -> OperatorSum

Get the in-place multiplication of an `OperatorSum` with a number.
"""
function mul!(ms::OperatorSum, factor::Number)
    abs(factor)==0 && return empty!(ms)
    for m in ms
        ms.contents[id(m)] = replace(m, value(m)*factor)
    end
    return ms
end

"""
    div!(ms::OperatorSum, factor::Number) -> OperatorSum

Get the in-place division of an `OperatorSum` with a number.
"""
function div!(ms::OperatorSum, factor::Number)
    isinf(factor) && return empty!(ms)
    for m in ms
        ms.contents[id(m)] = replace(m, value(m)/factor)
    end
    return ms
end

"""
    map!(f::Function, ms::OperatorSum; kwargs...) -> typeof(ms)

In place map of an `OperatorSum` by the function `f` elementally.
"""
function Base.map!(f::Function, ms::OperatorSum; kwargs...)
    for m in ms
        ms.contents[id(m)] = f(m; kwargs...)
    end
    return ms
end

# Operators
"""
    Operators{O<:Operator, I<:ID{OperatorUnit}}

A set of operators.

Type alias for `OperatorSum{O<:Operator, I<:ID{OperatorUnit}}`.
"""
const Operators{O<:Operator, I<:ID{OperatorUnit}} = OperatorSum{O, I}
@inline Base.summary(io::IO, ::Operators) = @printf io "Operators"

"""
    Operators(opts::Operator...)
    Operators{M}(opts::Operator...)

Get a set of operators.
"""
@inline Operators(opts::Operator...) = OperatorSum(opts)
@inline Operators{M}(opts::Operator...) where {M<:Operator} = OperatorSum{M}(opts)

"""
    adjoint(opts::Operators) -> Operators

Get the adjoint of a set of operators.
"""
function Base.adjoint(opts::Operators)
    result = zero(opts)
    for opt in opts
        add!(result, opt')
    end
    return result
end

"""
    ishermitian(opts::Operators) -> Bool

Judge whether a set of operators as a whole is Hermitian.
"""
@inline ishermitian(opts::Operators) = opts == opts'

# Arithmetic of quantum operators
"""
    optype(m::QuantumOperator)
    optype(::Type{<:QuantumOperator})

Get the corresponding `OperatorPack` type of a generic quantum operator.
"""
@inline optype(m::QuantumOperator) = optype(typeof(m))
@inline optype(::Type{M}) where {M<:OperatorUnit} = fulltype(Operator, NamedTuple{(:value, :id), Tuple{Int, Tuple{M}}})
@inline optype(::Type{M}) where {M<:OperatorPack} = M
@inline optype(::Type{M}) where {M<:OperatorSet} = eltype(M)

"""
    zero(::Type{M}) where {M<:QuantumOperator} -> OperatorSum

Get the zero sum.
"""
@inline Base.zero(::Type{M}) where {M<:OperatorUnit} = zero(optype(M))
@inline Base.zero(::Type{M}) where {M<:OperatorPack} = OperatorSum{M}()
@inline Base.zero(::Type{M}) where {M<:OperatorSum} = OperatorSum{eltype(M)}()

"""
    conj(m::OperatorUnit) -> OperatorUnit
    conj(m::OperatorPack) -> OperatorPack
    conj(m::OperatorSum) -> OperatorSum

Get the conjugation.
"""
@inline Base.conj(m::OperatorUnit) = m
@inline Base.conj(m::OperatorPack) = replace(m, conj(value(m)))
function Base.conj(ms::OperatorSum)
    result = zero(ms)
    for m in ms
        add!(result, conj(m))
    end
    return result
end

"""
    +(m::QuantumOperator) -> typeof(m)
    +(m₁::QuantumOperator, m₂::QuantumOperator) -> OperatorSum
    +(factor::Number, m::QuantumOperator) -> OperatorSum
    +(m::QuantumOperator, factor::Number) -> OperatorSum

Overloaded `+` between quantum operators.
"""
@inline Base.:+(m::QuantumOperator) = m
@inline function Base.:+(m₁::QuantumOperator, m₂::QuantumOperator)
    M = promote_type(optype(m₁), optype(m₂))
    result = OperatorSum{M}()
    add!(result, m₁)
    add!(result, m₂)
    return result
end
@inline Base.:+(factor::Number, m::QuantumOperator) = one(optype(m))*factor + m
@inline Base.:+(m::QuantumOperator, factor::Number) = m + one(optype(m))*factor

"""
    -(m::QuantumOperator) -> QuantumOperator
    -(m₁::QuantumOperator, m₂::QuantumOperator) -> OperatorSum
    -(factor::Number, m::QuantumOperator) -> OperatorSum
    -(m::QuantumOperator, factor::Number) -> OperatorSum

Overloaded `-` between quantum operators.
"""
@inline Base.:-(m::QuantumOperator) = m*(-1)
@inline function Base.:-(m₁::QuantumOperator, m₂::QuantumOperator)
    M = promote_type(optype(m₁), optype(m₂))
    result = OperatorSum{M}()
    add!(result, m₁)
    sub!(result, m₂)
    return result
end
@inline Base.:-(factor::Number, m::QuantumOperator) = one(optype(m))*factor - m
@inline Base.:-(m::QuantumOperator, factor::Number) = m - one(optype(m))*factor

"""
    *(factor::Number, m::OperatorUnit) -> Operator
    *(m::OperatorUnit, factor::Number) -> Operator
    *(m₁::OperatorUnit, m₂::OperatorUnit) -> Operator
    *(factor::Number, m::OperatorPack) -> OperatorPack
    *(m::OperatorPack, factor::Number) -> OperatorPack
    *(m₁::OperatorPack, m₂::OperatorUnit) -> OperatorPack
    *(m₁::OperatorUnit, m₁::OperatorPack) -> OperatorPack
    *(m₁::OperatorPack, m₂::OperatorPack) -> OperatorPack
    *(factor::Number, ms::OperatorSum) -> OperatorSum
    *(ms::OperatorSum, factor::Number) -> OperatorSum
    *(m::OperatorUnit, ms::OperatorSum) -> OperatorSum
    *(ms::OperatorSum, m::OperatorUnit) -> OperatorSum
    *(m::OperatorPack, ms::OperatorSum) -> OperatorSum
    *(ms::OperatorSum, m::OperatorPack) -> OperatorSum
    *(ms₁::OperatorSum, ms₂::OperatorSum) -> OperatorSum

Overloaded `*` between quantum operators or a quantum operator and a number.
"""
@inline Base.:*(factor::Number, m::OperatorUnit) = Operator(factor, m)
@inline Base.:*(m::OperatorUnit, factor::Number) = Operator(factor, m)
@inline Base.:*(m₁::OperatorUnit, m₂::OperatorUnit) = Operator(1, m₁, m₂)
@inline Base.:*(factor::Number, m::OperatorPack) = replace(m, factor*value(m))
@inline Base.:*(m::OperatorPack, factor::Number) = replace(m, value(m)*factor)
@inline Base.:*(m₁::OperatorPack, m₂::OperatorUnit) = m₁*Operator(1, m₂)
@inline Base.:*(m₁::OperatorUnit, m₂::OperatorPack) = Operator(1, m₁)*m₂
@inline function Base.:*(m₁::OperatorPack, m₂::OperatorPack)
    M₁, M₂ = typeof(m₁), typeof(m₂)
    @assert nameof(M₁)==nameof(M₂) && contentnames(M₁)==(:value, :id)==contentnames(M₂) "\"*\" error: not implemented between $(nameof(M₁)) and $(nameof(M₂))."
    return rawtype(M₁)(value(m₁)*value(m₂), ID(id(m₁), id(m₂)))
end
@inline Base.:*(factor::Number, ms::OperatorSum) = ms * factor
function Base.:*(ms::OperatorSum, factor::Number)
    result = OperatorSum{promote_type(eltype(ms), typeof(factor))}()
    abs(factor)==0 && return result
    for m in ms
        add!(result, m*factor)
    end
    return result
end
@inline Base.:*(m::Union{OperatorPack, OperatorUnit}, ms::OperatorSum) = OperatorSum(collect(m*mm for mm in ms))
@inline Base.:*(ms::OperatorSum, m::Union{OperatorPack, OperatorUnit}) = OperatorSum(collect(mm*m for mm in ms))
@inline Base.:*(ms₁::OperatorSum, ms₂::OperatorSum) = OperatorSum(collect(m₁*m₂ for m₁ in ms₁ for m₂ in ms₂))

"""
    /(m::QuantumOperator, factor::Number) -> QuantumOperator

Overloaded `/` between a quantum operator and a number.
"""
@inline Base.:/(m::QuantumOperator, factor::Number) = m * (one(dtype(optype(m)))/factor)

"""
    //(m::QuantumOperator, factor::Number) -> QuantumOperator

Overloaded `//` between a quantum operator and a number.
"""
@inline Base.://(m::QuantumOperator, factor::Number) = m * (1//factor)

"""
    ^(m::QuantumOperator, n::Integer) -> QuantumOperator

Overloaded `^` between a quantum operator and an integer.
"""
@inline Base.:^(m::QuantumOperator, n::Integer) = (@assert n>0 "^ error: non-positive integers are not allowed."; prod(ntuple(i->m, Val(n)), init=1))

# LaTeX format of quantum operators
"""
    LaTeX{SP, SB}(body, spdelimiter::String=", ", sbdelimiter::String=", "; options...) where {SP, SB}

LaTeX string representation of quantum operators.
"""
struct LaTeX{SP, SB, B<:Function, O}
    body::B
    spdelimiter::String
    sbdelimiter::String
    options::O
    function LaTeX{SP, SB}(body::Function, spdelimiter::String=", ", sbdelimiter::String=", "; options...) where {SP, SB}
        @assert isa(SP, Tuple{Vararg{Symbol}}) && isa(SB, Tuple{Vararg{Symbol}}) "LaTeX error: SP and SB must be tuple of symbols."
        new{SP, SB, typeof(body), typeof(options)}(body, spdelimiter, sbdelimiter, options)
    end
end
@inline superscript(::Type{<:LaTeX{SP}}) where SP = SP
@inline subscript(::Type{<:LaTeX{SP, SB} where SP}) where SB = SB
@inline LaTeX{SP, SB}(body, spdelimiter::String=", ", sbdelimiter::String=", "; options...) where {SP, SB} = LaTeX{SP, SB}(::OperatorUnit->body, spdelimiter, sbdelimiter; options...)

"""
    latexname(T::Type{<:OperatorUnit}) -> Symbol

Get the name of a type of `OperatorUnit` in the latex format lookups.
"""
@inline latexname(T::Type{<:OperatorUnit}) = nameof(T)

const latexformats = Dict{Symbol, LaTeX}()
"""
    latexformat(T::Type{<:OperatorUnit}) -> LaTeX
    latexformat(T::Type{<:OperatorUnit}, l::LaTeX) -> LaTeX

Get/Set the LaTeX format for a subtype of `OperatorUnit`.
"""
@inline latexformat(T::Type{<:OperatorUnit}) = latexformats[latexname(T)]
@inline latexformat(T::Type{<:OperatorUnit}, l::LaTeX) = latexformats[latexname(T)] = l

"""
    script(u::OperatorUnit, l::LaTeX, ::Val{:BD}) -> Any
    script(u::OperatorUnit, l::LaTeX, ::Val{:SP}) -> Tuple
    script(u::OperatorUnit, l::LaTeX, ::Val{:SB}) -> Tuple

Get the body/superscript/subscript of the LaTeX string representation of an operator unit.
"""
@inline script(u::OperatorUnit, l::LaTeX, ::Val{:BD}) = l.body(u)
@inline @generated script(u::OperatorUnit, l::LaTeX, ::Val{:SP}) = Expr(:tuple, [:(script(u, Val($sup); l.options...)) for sup in QuoteNode.(l|>superscript)]...)
@inline @generated script(u::OperatorUnit, l::LaTeX, ::Val{:SB}) = Expr(:tuple, [:(script(u, Val($sub); l.options...)) for sub in QuoteNode.(l|>subscript)]...)

"""
    latexstring(u::OperatorUnit) -> String

LaTeX string representation of an operator unit.
"""
@inline function latexstring(u::OperatorUnit)
    l = latexformat(typeof(u))
    body = script(u, l, Val(:BD))
    superscript = join((str for str in script(u, l, Val(:SP)) if length(str)>0), l.spdelimiter)
    subscript = join((str for str in script(u, l, Val(:SB)) if length(str)>0), l.sbdelimiter)
    return @sprintf "%s^{%s}_{%s}" body superscript subscript
end

"""
    latexstring(opt::OperatorProd) -> String

Get the string representation of an operator in the LaTeX format.
"""
function latexstring(opt::OperatorProd)
    rank(opt)==0 && return replace(valuetolatextext(value(opt)), " "=>"")
    poses = Int[]
    push!(poses, 1)
    for i = 2:rank(opt)
        opt[i]≠opt[i-1] && push!(poses, i)
    end
    push!(poses, rank(opt)+1)
    result = valuetostr(value(opt))
    for i = 1:length(poses)-1
        order = poses[i+1] - poses[i]
        if order == 1
            result = @sprintf "%s%s" result latexstring(opt[poses[i]])
        else
            result = @sprintf "%s(%s)^%s" result latexstring(opt[poses[i]]) order
        end
    end
    return result
end
function valuetostr(v)
    v==+1 && return ""
    v==-1 && return "-"
    result = valuetolatextext(v)
    if occursin(" ", result) || (isa(v, Complex) && real(v)≠0 && imag(v)≠0)
        bra, ket = occursin("(", result) ? ("[", "]") : ("(", ")")
        result = @sprintf "%s%s%s" bra replace(result, " "=>"") ket
    end
    return result
end
@inline valuetolatextext(value::Union{Real, Complex}) = decimaltostr(value)
function valuetolatextext(value)
    io = IOBuffer()
    if showable(MIME"text/latex"(), value)
        show(IOContext(io, :limit=>false), MIME"text/latex"(), value)
    else
        show(IOContext(io, :limit=>false), MIME"text/plain"(), value)
    end
    return replace(replace(replace(replace(String(take!(io)), "\\begin{equation*}" => ""), "\\end{equation*}" => ""), "\$" => ""), "\n" => "")
end

"""
    latexstring(opts::OperatorSet) -> String

Get the string representation of a set of operators in the LaTeX format.
"""
function latexstring(opts::OperatorSet)
    result = String[]
    for (i, opt) in enumerate(opts)
        rep = latexstring(opt)
        i>1 && rep[1]≠'-' && push!(result, "+")
        push!(result, rep)
    end
    return join(result, "")
end

"""
    show(io::IO, ::MIME"text/latex", m::QuantumOperator)

Show a quantum operator.
"""
@inline Base.show(io::IO, ::MIME"text/latex", m::QuantumOperator) = show(io, MIME"text/latex"(), latexstring(unicode2latex(latexstring(m))))

# Transformations
"""
    Transformation <: Function

Abstract transformation on quantum operators.
"""
abstract type Transformation <: Function end
@inline Base.:(==)(transformation₁::Transformation, transformation₂::Transformation) = ==(efficientoperations, transformation₁, transformation₂)
@inline Base.isequal(transformation₁::Transformation, transformation₂::Transformation) = isequal(efficientoperations, transformation₁, transformation₂)

"""
    LinearTransformation <: Transformation

Abstract linear transformation on quantum operators.
"""
abstract type LinearTransformation <: Transformation end
@inline Base.valtype(transformation::LinearTransformation, m::QuantumOperator) = valtype(typeof(transformation), typeof(m))
@inline Base.zero(transformation::LinearTransformation, m::QuantumOperator) = zero(valtype(transformation, m))

"""
    (transformation::LinearTransformation)(ms::OperatorSet; kwargs...) -> OperatorSet

Get the linear transformed quantum operators.
"""
function (transformation::LinearTransformation)(ms::OperatorSet; kwargs...)
    result = zero(transformation, ms)
    for m in ms
        add!(result, transformation, m; kwargs...)
    end
    return result
end

"""
    add!(destination, transformation::LinearTransformation, op::OperatorPack; kwargs...) -> typeof(destination)
    add!(destination, transformation::LinearTransformation, op::OperatorSet; kwargs...) -> typeof(destination)

Add the result of the linear transformation on a quantum operator to the destination.
"""
@inline function add!(destination, transformation::LinearTransformation, op::OperatorPack; kwargs...)
    return add!(destination, transformation(op; kwargs...))
end
@inline function add!(destination, transformation::LinearTransformation, op::OperatorSet; kwargs...)
    for m in op
        add!(destination, transformation, m; kwargs...)
    end
    return destination
end

"""
    LinearFunction{F<:Function} <: LinearTransformation

Wrapper a function to be a linear transformation.
"""
struct LinearFunction{F<:Function} <: LinearTransformation
    f::F
end
@inline (f::LinearFunction)(op::OperatorPack; kwargs...) = f.f(op; kwargs...)
@inline (f::LinearFunction)(op::OperatorSet; kwargs...) = f.f(op; kwargs...)

"""
    Matrixization <: LinearTransformation

The matrixization transformation.
"""
abstract type Matrixization <: LinearTransformation end

"""
    matrix

Generic matrix representation.
"""
function matrix end

"""
    Permutation{T} <: LinearTransformation

The permutation transformation.
"""
struct Permutation{T} <: LinearTransformation
    table::T
end
@inline function Base.valtype(::Type{<:Permutation}, M::Type{<:OperatorProd})
    M = reparameter(M, :id, ID{eltype(M)})
    return OperatorSum{M, idtype(M)}
end
@inline Base.valtype(P::Type{<:Permutation}, M::Type{<:OperatorSum}) = valtype(P, eltype(M))

"""
    (permutation::Permutation)(m::OperatorProd; rev::Bool=false, kwargs...) -> OperatorSum

Permute the operator units of an `OperatorProd` to the descending order according to the table contained in `permutation`.

!!! note
    To use this function, the user must implement a method of `permute`, which computes the result of the permutation of two operator units:
    ```julia
    permute(u₁::OperatorUnit, u₂::OperatorUnit) -> Union{OperatorProd, OperatorSum}
    ```
    Here, `u₁` and `u₂` are two arbitrary operator units contained in `id(m)`.
"""
@inline function (permutation::Permutation)(m::OperatorProd; rev::Bool=false, kwargs...)
    result = zero(permutation, m)
    cache = eltype(result)[m]
    while length(cache) > 0
        current = pop!(cache)
        pos = operatorprodcommuteposition(sequence(current, permutation.table); rev=rev)
        if isnothing(pos)
            add!(result, current)
        else
            left = current[1:pos-1] * value(m)
            right = current[pos+2:end]
            for middle in permute(id(current)[pos], id(current)[pos+1])
                temp = left * middle * right
                isnothing(temp) || push!(cache, temp)
            end
        end
    end
    return result
end
function operatorprodcommuteposition(sequences; rev::Bool)
    pos = 1
    while pos < length(sequences)
        ((rev && sequences[pos]>sequences[pos+1]) || (!rev && sequences[pos]<sequences[pos+1])) && return pos
        pos += 1
    end
    return nothing
end

"""
    UnitSubstitution{U<:OperatorUnit, S<:OperatorSum} <: LinearTransformation

The "unit substitution" transformation, which substitutes each `OperatorUnit` in the old quantum operators to a new expression represented by an `OperatorSum`.
"""
abstract type UnitSubstitution{U<:OperatorUnit, S<:OperatorSum} <: LinearTransformation end
@inline function Base.valtype(::Type{<:UnitSubstitution{U, S}}, M::Type{<:OperatorProd}) where {U<:OperatorUnit, S<:OperatorSum}
    @assert U<:eltype(eltype(M)) "valtype error: mismatched unit transformation."
    V = fulltype(eltype(S), NamedTuple{(:value, :id), Tuple{promote_type(valtype(eltype(S)), valtype(M)), ID{eltype(eltype(S))}}})
    return OperatorSum{V, idtype(V)}
end
@inline Base.valtype(P::Type{<:UnitSubstitution}, M::Type{<:OperatorSum}) = valtype(P, eltype(M))

"""
    (unitsubstitution::UnitSubstitution)(m::OperatorProd; kwargs...) -> OperatorSum

Substitute every `OperatorUnit` in an `OperatorProd` with a new `OperatorSum`.
"""
function (unitsubstitution::UnitSubstitution)(m::OperatorProd; kwargs...)
    return prod(ntuple(i->unitsubstitution(m[i]; kwargs...), Val(rank(m))), init=value(m))
end

"""
    TabledUnitSubstitution{U<:OperatorUnit, S<:OperatorSum, T<:AbstractDict{U, S}} <: UnitSubstitution{U, S}

A concrete "unit substitution" transformation, which stores every substitution of the old `OperatorUnit`s in its table as a dictionary.
"""
struct TabledUnitSubstitution{U<:OperatorUnit, S<:OperatorSum, T<:AbstractDict{U, S}} <: UnitSubstitution{U, S}
    table::T
    function TabledUnitSubstitution(table::AbstractDict{<:OperatorUnit, <:OperatorSum})
        new{keytype(table), valtype(table), typeof(table)}(table)
    end
end
(unitsubstitution::TabledUnitSubstitution)(m::OperatorUnit; kwargs...) = unitsubstitution.table[m]

"""
    RankFilter{R} <: LinearTransformation

Rank filter, which filters out the `OperatorProd` with a given rank `R`.
"""
struct RankFilter{R} <: LinearTransformation
    function RankFilter(rank::Int)
        @assert rank>=0 "RankFilter error: the wanted rank must be non-negative."
        new{rank}()
    end
end
@inline Base.valtype(::Type{RankFilter{R}}, M::Type{<:OperatorProd}) where R = reparameter(M, :id, ID{eltype(M), R})
@inline function Base.valtype(R::Type{<:RankFilter}, M::Type{<:OperatorSum})
    V = valtype(R, eltype(M))
    return OperatorSum{V, idtype(V)}
end
@inline rank(rf::RankFilter) = rank(typeof(rf))
@inline rank(::Type{RankFilter{R}}) where R = R
@inline @generated (rf::RankFilter)(m::OperatorProd; kwargs...) = rank(m)==rank(rf) ? :(m) : 0

end #module
