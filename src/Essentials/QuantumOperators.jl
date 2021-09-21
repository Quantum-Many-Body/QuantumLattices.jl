module QuantumOperators

using Printf: @printf, @sprintf
using ...Prerequisites: atol, rtol
using ...Prerequisites.NamedVectors: NamedVector
using ...Prerequisites.Traits: efficientoperations, getcontent, contentorder
using ...Prerequisites.Traits: rawtype, fulltype, parametertype, parameterpairs, reparameter, promoteparameters
using ...Prerequisites.CompositeStructures: CompositeDict

import ...Interfaces: id, value, rank, add!, sub!, mul!, div!, ⊗, ⋅, permute
import ...Prerequisites.Traits: contentnames, dissolve, isparameterbound, parameternames

export SingularID, ID, QuantumOperator, OperatorProd, OperatorSum, Scalar, idtype, sequence
export Transformation, Identity, Numericalization

"""
    SingularID <: NamedVector

A singular id is the building block of the id system of quantum operators.
"""
abstract type SingularID <: NamedVector end

"""
    ID{I<:SingularID, N}

The (composite) id system of quantum operators.

Type alias for `NTuple{N, I} where {N, I<:SingularID}`.
"""
const ID{I<:SingularID, N} = NTuple{N, I}
@inline Base.show(io::IO, cid::Tuple{SingularID, Vararg{SingularID}}) = @printf io "ID(%s)" join(cid, ", ")

"""
    ID(ids::SingularID...)
    ID(ids::NTuple{N, SingularID}) where N

Get the composite id from singular ids.
"""
@inline ID(ids::SingularID...) = ids
@inline ID(ids::NTuple{N, SingularID}) where {N} = ids

"""
    ID(::Type{SID}, attrs::Vararg{NTuple{N}, M}) where {SID<:SingularID, N, M}

Get the composite id from the components of singular ids.
"""
@inline @generated function ID(::Type{SID}, attrs::Vararg{NTuple{N, Any}, M}) where {SID<:SingularID, N, M}
    exprs = []
    for i = 1:N
        args = [:(attrs[$j][$i]) for j = 1:M]
        push!(exprs, :(SID($(args...))))
    end
    return :(ID($(exprs...)))
end

"""
    propertynames(::Type{I}) where I<:ID{SingularID} -> Tuple{Vararg{Symbol}}

Get the property names of a composite id.
"""
@inline @generated function Base.propertynames(I::ID{SingularID})
    exprs = [QuoteNode(Symbol(name, 's')) for name in I|>eltype|>fieldnames]
    return Expr(:tuple, exprs...)
end

"""
    getproperty(cid::ID{SingularID}, name::Symbol)

Get the property of a composite id.
"""
@inline Base.getproperty(cid::ID{SingularID}, name::Symbol) = idgetproperty(cid, Val(name), Val(cid|>propertynames))
@inline @generated function idgetproperty(cid::ID{SingularID}, ::Val{name}, ::Val{names}) where {name, names}
    index = findfirst(isequal(name), names)::Int
    exprs = [:(getfield(cid[$i], $index)) for i = 1:fieldcount(cid)]
    return Expr(:tuple, exprs...)
end

"""
    promote_type(::Type{Tuple{}}, I::Type{<:ID{SingularID, N}}) where N
    promote_type(I::Type{<:ID{SingularID, N}}, ::Type{Tuple{}}) where N

Define the promote rule for ID types.
"""
@inline Base.promote_type(::Type{Tuple{}}, I::Type{<:Tuple{SingularID, Vararg{SingularID}}}) = ID{I|>eltype}
@inline Base.promote_type(I::Type{<:Tuple{SingularID, Vararg{SingularID}}}, ::Type{Tuple{}}) = ID{I|>eltype}

"""
    isless(::Type{<:SingularID}, cid₁::ID{SingularID}, cid₂::ID{SingularID}) -> Bool

Compare two ids and judge whether the first is less than the second.

The comparison rule are as follows:
1. ids with smaller ranks are always less than those with higher ranks;
2. if two ids are of the same rank, the comparison goes just like that between tuples.
"""
@inline function Base.isless(::Type{<:SingularID}, cid₁::ID{SingularID}, cid₂::ID{SingularID})
    r₁, r₂ = cid₁|>rank, cid₂|>rank
    (r₁ < r₂) ? true : (r₁ > r₂) ? false : isless(cid₁, cid₂)
end

"""
    rank(id::ID{SingularID}) -> Int
    rank(::Type{<:ID{SingularID}}) -> Any
    rank(::Type{<:ID{SingularID, N}}) where N -> Int

Get the rank of a composite id.
"""
@inline rank(id::ID{SingularID}) = rank(typeof(id))
@inline rank(::Type{<:ID{SingularID}}) = Any
@inline rank(::Type{<:ID{SingularID, N}}) where N = N

"""
    *(sid₁::SingularID, sid₂::SingularID) -> ID{SingularID}
    *(sid::SingularID, cid::ID{SingularID}) -> ID{SingularID}
    *(cid::ID{SingularID}, sid::SingularID) -> ID{SingularID}
    *(cid₁::ID{SingularID}, cid₂::ID{SingularID}) -> ID{SingularID}

Get the product of the id system.
"""
@inline Base.:*(sid₁::SingularID, sid₂::SingularID) = ID(sid₁, sid₂)
@inline Base.:*(sid::SingularID, cid::ID{SingularID}) = ID(sid, cid...)
@inline Base.:*(cid::ID{SingularID}, sid::SingularID) = ID(cid..., sid)
@inline Base.:*(cid₁::ID{SingularID}, cid₂::ID{SingularID}) = ID(cid₁..., cid₂...)

"""
    QuantumOperator{I, V} <:CompositeDict{I, V}

The abstract type of any quantum operator.
"""
abstract type QuantumOperator{I, V} <:CompositeDict{I, V} end
@inline Base.:(==)(m₁::QuantumOperator, m₂::QuantumOperator) = ==(efficientoperations, m₁, m₂)
@inline Base.isequal(m₁::QuantumOperator, m₂::QuantumOperator) = isequal(efficientoperations, m₁, m₂)
@inline parameternames(::Type{<:QuantumOperator}) = (:id, :value)
@inline isparameterbound(::Type{<:QuantumOperator}, ::Val{:id}, ::Type{I}) where {I<:ID{SingularID}} = !isconcretetype(I)
@inline isparameterbound(::Type{<:QuantumOperator}, ::Val{:value}, ::Type{V}) where V = !isconcretetype(V) 

"""
    idtype(m::QuantumOperator)
    idtype(::Type{T}) where {T<:QuantumOperator}

The type of the id of a `QuantumOperator`.
"""
@inline idtype(m::QuantumOperator) = idtype(typeof(m))
@inline @generated idtype(::Type{T}) where {T<:QuantumOperator} = parametertype(supertype(T, :QuantumOperator), :id)

"""
    valtype(m::QuantumOperator)
    valtype(::Type{T}) where {T<:QuantumOperator}

Get the type of the value of a `QuantumOperator`.
"""
@inline Base.valtype(m::QuantumOperator) = valtype(typeof(m))
@inline @generated Base.valtype(::Type{T}) where {T<:QuantumOperator} = parametertype(supertype(T, :QuantumOperator), :value)

"""
    rank(m::QuantumOperator) -> Int
    rank(::Type{M}) where M<:QuantumOperator -> Int

Get the rank of a `QuantumOperator`.
"""
@inline rank(m::QuantumOperator) = rank(typeof(m))
@inline rank(::Type{M}) where {M<:QuantumOperator} = rank(idtype(M))

"""
    OperatorProd{V, I<:ID{SingularID}} <: QuantumOperator{I, V}

The entity that represent the product of several quantum operators.

Basically, a concrete subtype should contain two predefined contents:
- `value::V`: the coefficient of the product
- `id::I`: the total id of the product
"""
abstract type OperatorProd{V, I<:ID{SingularID}} <: QuantumOperator{I, V} end
@inline contentnames(::Type{<:OperatorProd}) = (:value, :id)
@inline parameternames(::Type{<:OperatorProd}) = (:value, :id)
@inline isparameterbound(::Type{<:OperatorProd}, ::Val{:value}, ::Type{V}) where V = false

@inline newvalue(m::OperatorProd, v::Number) = v
"""
    replace(m::OperatorProd, v::Number) -> OperatorProd

Replace the value of an `OperatorProd`.
"""
@inline Base.replace(m::OperatorProd, v::Number) = rawtype(typeof(m))(dissolve(m, newvalue, (v,))...)
@inline dissolve(m::OperatorProd, ::Val{:value}, ::typeof(newvalue), args::Tuple, kwargs::NamedTuple) = newvalue(m, args...; kwargs...)

"""
    value(m::OperatorProd) -> valtype(m)

Get the value of an `OperatorProd`.
"""
@inline value(m::OperatorProd) = getcontent(m, :value)

"""
    id(m::OperatorProd) -> idtype(m)

Get the id of an `OperatorProd`.
"""
@inline id(m::OperatorProd) = getcontent(m, :id)

"""
    Scalar{V}

Scalar quantum operator.
"""
const Scalar{V} = OperatorProd{V, Tuple{}}

"""
    isapprox(m₁::OperatorProd, m₂::OperatorProd; atol::Real=atol, rtol::Real=rtol) -> Bool

Compare two `OperatorProd`s and judge whether they are approximate to each other.
"""
@inline function Base.isapprox(m₁::OperatorProd, m₂::OperatorProd; atol::Real=atol, rtol::Real=rtol)
    isapprox(efficientoperations, contentorder(typeof(m₁), :value)|>Val, dissolve(m₁), dissolve(m₂); atol=atol, rtol=rtol)::Bool
end

"""
    replace(m::OperatorProd; kwargs...) -> typeof(m)

Return a copy of a concrete `OperatorProd` with some of the field values replaced by the keyword arguments.
"""
@inline Base.replace(m::OperatorProd; kwargs...) = replace(efficientoperations, m; kwargs...)

"""
    promote_rule(::Type{M₁}, ::Type{M₂}) where {M₁<:OperatorProd, M₂<:OperatorProd}

Define the promote rule for `OperatorProd` types.
"""
@inline function Base.promote_rule(::Type{M₁}, ::Type{M₂}) where {M₁<:OperatorProd, M₂<:OperatorProd}
    M₁<:M₂ && return M₂
    M₂<:M₁ && return M₁
    r₁, r₂ = rank(M₁), rank(M₂)
    M = r₂==0 ? rawtype(M₁) : r₁==0 ? rawtype(M₂) : typejoin(rawtype(M₁), rawtype(M₂))
    return fulltype(M, promoteparameters(parameterpairs(M₁), parameterpairs(M₂)))
end

"""
    promote_type(::Type{M}, ::Type{V}, ::Val{:*}) where {M<:OperatorProd, V<:Number}
    promote_type(::Type{V}, ::Type{M}, ::Val{:*}) where {M<:OperatorProd, V<:Number}

Define the promote rule for the multiplication between an `OperatorProd` and a scalar.
"""
@inline Base.promote_type(::Type{M}, ::Type{V}, ::Val{:*}) where {M<:OperatorProd, V<:Number} = promote_type(V, M, Val(:*))
@inline function Base.promote_type(::Type{V}, ::Type{M}, ::Val{:*}) where {M<:OperatorProd, V<:Number}
    return reparameter(M, :value, promote_type(valtype(M), V))
end

"""
    getindex(m::OperatorProd, i) -> OperatorProd

Overloaded `[]`.
"""
@inline Base.getindex(m::OperatorProd, i) = rawtype(typeof(m))(dissolve(m, getindex, (i,))...)
@inline dissolve(m::OperatorProd, ::Val{:value}, ::typeof(getindex), ::Tuple{Any}, ::NamedTuple) = one(valtype(m))
@inline dissolve(m::OperatorProd, ::Val{:id}, ::typeof(getindex), i::Tuple{Any}, ::NamedTuple) = ID(id(m)[first(i)])

"""
    length(m::OperatorProd) -> Int

Get the length of an `OperatorProd`.
"""
@inline Base.length(m::OperatorProd) = rank(m)
@inline Base.firstindex(m::OperatorProd) = 1
@inline Base.lastindex(m::OperatorProd) = rank(m)

"""
    one(::Type{M}) where M<:OperatorProd
    one(m::OperatorProd)

Get the identity quantum operator.
"""
@inline function Base.one(::Type{M}) where M<:OperatorProd
    rtype = rawtype(M)
    vtype = isconcretetype(valtype(M)) ? valtype(M) : Int
    @assert fieldnames(rtype) == (:value, :id) "one error: not supported type($(nameof(rtype)))."
    return rtype(one(vtype), ID())
end
@inline Base.one(m::OperatorProd) = rawtype(typeof(m))(dissolve(m, one)...)
@inline dissolve(m::OperatorProd, ::Val{:value}, ::typeof(one), ::Tuple, ::NamedTuple) = one(valtype(m))
@inline dissolve(::OperatorProd, ::Val{:id}, ::typeof(one), ::Tuple, ::NamedTuple) = ID()

"""
    convert(::Type{M}, m::Scalar) where {M<:Scalar}
    convert(::Type{M}, m::Number) where {M<:Scalar}
    convert(::Type{M}, m::OperatorProd) where {M<:OperatorProd}

1) Convert a scalar quantum operator from one type to another;
2) Convert a scalar to a scalar quantum operator;
3) Convert a quantum operator from one type to another.
"""
@inline Base.convert(::Type{M}, m::Scalar) where {M<:Scalar} = typeof(m)<:M ? m : one(M)*value(m)
@inline Base.convert(::Type{M}, m::Number) where {M<:Scalar} = one(M)*m
@inline function Base.convert(::Type{M}, m::OperatorProd) where {M<:OperatorProd}
    (typeof(m) <: M) && return m
    @assert convertible(M, m) "convert error: $(typeof(m)) cannot be converted to $M."
    return replace(m, convert(valtype(M), value(m)))
end
function convertible(::Type{M}, m::OperatorProd) where {M<:OperatorProd}
    !(rawtype(typeof(m)) <: rawtype(M)) && return false
    S = supertype(typeof(m), nameof(M))
    for (i, name) in enumerate(parameternames(M))
        (name != :value) && !(parametertype(S, i) <: parametertype(M, i)) && return false
    end
    return true
end

"""
    split(m::OperatorProd) -> Tuple{Any, Vararg{OperatorProd}}

Split an `OperatorProd` into the coefficient and a sequence of rank-1 `OperatorProd`s.
"""
@inline @generated function Base.split(m::OperatorProd)
    exprs = [:(value(m))]
    for i = 1:rank(m)
        push!(exprs, :(m[$i]))
    end
    return Expr(:tuple, exprs...)
end

"""
    OperatorSum{I<:ID{SingularID}, M<:OperatorProd} <: QuantumOperator{I, M}

The sum of `OperatorProd`s.

Similar items are automatically merged with the aid of the id system.
"""
struct OperatorSum{I<:ID{SingularID}, M<:OperatorProd} <: QuantumOperator{I, M}
    contents::Dict{I, M}
    OperatorSum(contents::Dict{<:ID{SingularID}, <:OperatorProd}) = new{keytype(contents), valtype(contents)}(contents)
end
@inline function Base.promote_rule(::Type{MS₁}, ::Type{MS₂}) where {MS₁<:OperatorSum, MS₂<:OperatorSum}
    I = promote_type(idtype(MS₁), idtype(MS₂))
    M = promote_type(valtype(MS₁), valtype(MS₂))
    return OperatorSum{I, M}
end
function Base.show(io::IO, ms::OperatorSum)
    @printf io "OperatorSum with %s %s:\n" length(ms) nameof(valtype(ms))
    for m in values(ms)
        @printf io "  %s\n" m
    end
end
@inline Base.zero(ms::OperatorSum) = zero(typeof(ms))
@inline Base.zero(::Type{MS}) where {MS<:OperatorSum} = MS()

"""
    OperatorSum(ms::OperatorProd...)
    OperatorSum{I, M}(ms::OperatorProd...) where {I<:ID{SingularID}, M<:OperatorProd}

Get the sum of `OperatorProd`s with similar items merged.
"""
@inline OperatorSum(ms::OperatorSum) = OperatorSum{idtype(ms), valtype(ms)}(ms)
@inline OperatorSum(ms::OperatorProd...) = OperatorSum{idtype(eltype(ms)), eltype(ms)}(ms...)
@inline  OperatorSum{I, M}(ms::OperatorSum) where {I<:ID{SingularID}, M<:OperatorProd} = OperatorSum(Dict{I, M}(ms.contents))
function OperatorSum{I, M}(ms::OperatorProd...) where {I<:ID{SingularID}, M<:OperatorProd}
    contents = Dict{I, M}()
    for m in ms
        abs(value(m))==0 || (contents[id(m)] = m)
    end
    return OperatorSum(contents)
end

"""
    repr(ms::OperatorSum) -> String

Get the repr representation of a sum of `OperatorProd`s.
"""
function Base.repr(ms::OperatorSum)
    cache = [@sprintf("OperatorSum with %s %s:", length(ms), nameof(valtype(ms)))]
    for m in ms|>values
        push!(cache, @sprintf("  %s", repr(m)))
    end
    return join(cache, "\n")
end 

"""
    isapprox(ms₁::OperatorSum, ms₂::OperatorSum; atol::Real=atol, rtol::Real=rtol) -> Bool

Compare two `OperatorSum`s and judge whether they are approximate to each other.
"""
@inline function Base.isapprox(ms₁::OperatorSum, ms₂::OperatorSum; atol::Real=atol, rtol::Real=rtol)
    for (k, m) in ms₁
        isapprox(value(m), 0, atol=atol, rtol=rtol) && continue
        haskey(ms₂, k) || return false
        isapprox(m, ms₂[k]) || return false
    end
    for (k, m) in ms₂
        isapprox(value(m), 0,  atol=atol, rtol=rtol) && continue
        haskey(ms₁, k) || return false
        isapprox(m, ms₁[k]) || return false
    end
    return true
end

"""
    add!(ms::OperatorSum) -> typeof(ms)
    add!(ms::OperatorSum, m::Number) -> typeof(ms)
    add!(ms::OperatorSum, m::OperatorProd) -> typeof(ms)
    add!(ms::OperatorSum, mms::OperatorSum) -> typeof(ms)

Get the in-place addition of quantum operators.
"""
@inline add!(ms::OperatorSum) = ms
@inline add!(ms::OperatorSum, m::Number) = add!(ms, one(valtype(ms))*m)
@inline function add!(ms::OperatorSum, m::OperatorProd)
    m = convert(ms|>valtype, m)
    old = get(ms, id(m), nothing)
    new = isnothing(old) ? m : replace(old, value(old)+value(m))
    abs(value(new))==0.0 ? delete!(ms, id(m)) : (ms[id(m)] = new)
    return ms
end
@inline function add!(ms::OperatorSum, mms::OperatorSum)
    for m in mms|>values
        add!(ms, m)
    end
    return ms
end

"""
    sub!(ms::OperatorSum) -> typeof(ms)
    sub!(ms::OperatorSum, m::Number) -> typeof(ms)
    sub!(ms::OperatorSum, m::OperatorProd) -> typeof(ms)
    sub!(ms::OperatorSum, mms::OperatorSum) -> typeof(ms)

Get the in-place subtraction of quantum operators.
"""
@inline sub!(ms::OperatorSum) = ms
@inline sub!(ms::OperatorSum, m::Number) = add!(ms, one(valtype(ms))*(-m))
@inline function sub!(ms::OperatorSum, m::OperatorProd)
    m = convert(ms|>valtype, m)
    old = get(ms, id(m), nothing)
    new = isnothing(old) ? -m : replace(old, value(old)-value(m))
    abs(value(new))==0.0 ? delete!(ms, id(m)) : (ms[id(m)] = new)
    return ms
end
@inline function sub!(ms::OperatorSum, mms::OperatorSum)
    for m in mms|>values
        sub!(ms, m)
    end
    return ms
end

"""
    mul!(ms::OperatorSum, factor::Scalar) -> OperatorSum
    mul!(ms::OperatorSum, factor::Number) -> OperatorSum

Get the in-place multiplication of an `OperatorSum` with a scalar.
"""
@inline mul!(ms::OperatorSum, factor::Scalar) = mul!(ms, value(factor))
function mul!(ms::OperatorSum, factor::Number)
    @assert isa(one(ms|>valtype|>valtype)*factor, ms|>valtype|>valtype) "mul! error: mismatched type, $(ms|>valtype) and $(factor|>typeof)."
    for m in values(ms)
        ms[id(m)] = replace(m, value(m)*factor)
    end
    return ms
end

"""
    div!(ms::OperatorSum, factor::Scalar) -> OperatorSum
    div!(ms::OperatorSum, factor::Number) -> OperatorSum

Get the in-place division of an `OperatorSum` with a scalar.
"""
@inline div!(ms::OperatorSum, factor::Scalar) = div!(ms, value(factor))
function div!(ms::OperatorSum, factor::Number)
    @assert isa(one(ms|>valtype|>valtype)/factor, ms|>valtype|>valtype) "div! error: mismatched type, $(ms|>valtype) and $(factor|>typeof)."
    for m in values(ms)
        ms[id(m)] = replace(m, value(m)/factor)
    end
    return ms
end

"""
    +(m::OperatorProd) -> typeof(m)
    +(m::OperatorProd, factor::Number) -> OperatorSum
    +(factor::Number, m::OperatorProd) -> OperatorSum
    +(m₁::OperatorProd, m₂::OperatorProd) -> OperatorSum
    +(ms::OperatorSum) -> typeof(ms)
    +(ms::OperatorSum, factor::Number) -> OperatorSum
    +(factor::Number, ms::OperatorSum) -> OperatorSum
    +(ms::OperatorSum, m::OperatorProd) -> OperatorSum
    +(m::OperatorProd, ms::OperatorSum) -> OperatorSum
    +(ms₁::OperatorSum, ms₂::OperatorSum) -> OperatorSum

Overloaded `+` between quantum operators.
"""
@inline Base.:+(m::OperatorProd) = m
@inline Base.:+(ms::OperatorSum) = ms
@inline Base.:+(factor::Number, m::OperatorProd) = m + one(m)*factor
@inline Base.:+(m::OperatorProd, factor::Number) = m + one(m)*factor
@inline Base.:+(factor::Number, m::Scalar) = replace(m, value(m)+factor)
@inline Base.:+(m::Scalar, factor::Number) = replace(m, value(m)+factor)
@inline Base.:+(m₁::Scalar, m₂::Scalar) = replace(m₁, value(m₁)+value(m₂))
@inline Base.:+(factor::Number, ms::OperatorSum) = ms + one(valtype(ms))*factor
@inline Base.:+(ms::OperatorSum, factor::Number) = ms + one(valtype(ms))*factor
@inline function Base.:+(m₁::OperatorProd, m₂::OperatorProd)
    I = promote_type(m₁|>idtype, m₂|>idtype)
    M = promote_type(typeof(m₁), typeof(m₂))
    return OperatorSum{I, M}(m₁, m₂)
end
@inline Base.:+(m::OperatorProd, ms::OperatorSum) = ms + m
@inline function Base.:+(ms::OperatorSum, m::OperatorProd)
    I = promote_type(ms|>idtype, m|>idtype)
    M = promote_type(valtype(ms), typeof(m))
    return add!(OperatorSum{I, M}(ms), m)
end
@inline function Base.:+(ms₁::OperatorSum, ms₂::OperatorSum)
    I = promote_type(ms₁|>idtype, ms₂|>idtype)
    M = promote_type(valtype(ms₁), valtype(ms₂))
    return add!(OperatorSum{I, M}(ms₁), ms₂)
end

"""
    *(factor::Number, m::OperatorProd) -> OperatorProd
    *(m::OperatorProd, factor::Number) -> OperatorProd
    *(m₁::OperatorProd, m₂::OperatorProd) -> OperatorProd
    *(factor::Number, ms::OperatorSum) -> OperatorSum
    *(ms::OperatorSum, factor::Number) -> OperatorSum
    *(m::OperatorProd, ms::OperatorSum) -> OperatorSum
    *(ms::OperatorSum, m::OperatorProd) -> OperatorSum
    *(ms₁::OperatorSum, ms₂::OperatorSum) -> OperatorSum

Overloaded `*` between quantum operators or a quantum operator and a scalar.
"""
@inline Base.:*(factor::Scalar, m::OperatorProd) = m * value(factor)
@inline Base.:*(m::OperatorProd, factor::Scalar) = m * value(factor)
@inline Base.:*(m₁::Scalar, m₂::Scalar) = replace(m₁, value(m₁)*value(m₂))
@inline Base.:*(factor::Number, m::OperatorProd) = m * factor
@inline Base.:*(m::OperatorProd, factor::Number) = replace(m, factor*value(m))
@inline Base.:*(factor::Scalar, ms::OperatorSum) = ms * value(factor)
@inline Base.:*(ms::OperatorSum, factor::Scalar) = ms * value(factor)
@inline Base.:*(factor::Number, ms::OperatorSum) = ms * factor
@inline function Base.:*(ms::OperatorSum, factor::Number)
    abs(factor)==0 && return empty(ms)
    result = OperatorSum{idtype(ms), promote_type(valtype(ms), typeof(factor), Val(:*))}()
    for (id, m) in ms
        result[id] = m * factor
    end
    return result
end
@inline Base.:*(m::OperatorProd, ms::OperatorSum) = OperatorSum((m * mm for mm in ms|>values)...)
@inline Base.:*(ms::OperatorSum, m::OperatorProd) = OperatorSum((mm * m for mm in ms|>values)...)
@inline Base.:*(ms₁::OperatorSum, ms₂::OperatorSum) = OperatorSum((m₁ * m₂ for m₁ in ms₁|>values for m₂ in ms₂|>values)...)
@inline function Base.:*(m₁::OperatorProd, m₂::OperatorProd)
    @assert((m₁|>typeof|>nameof == m₂|>typeof|>nameof) && (m₁|>typeof|>fieldcount == m₂|>typeof|>fieldcount == 2),
            "\"*\" error: not implemented between $(m₁|>typeof|>nameof) and $(m₂|>typeof|>nameof)."
            )
    return rawtype(typeof(m₁))(value(m₁)*value(m₂), id(m₁)*id(m₂))
end

"""
    -(m::OperatorProd) -> typeof(m)
    -(m::OperatorProd, factor::Number) -> OperatorSum
    -(factor::Number, m::OperatorProd) -> OperatorSum
    -(m₁::OperatorProd, m₂::OperatorProd) -> OperatorSum
    -(ms::OperatorSum) -> typeof(ms)
    -(ms::OperatorSum, factor::Number) -> OperatorSum
    -(factor::Number, ms::OperatorSum) -> OperatorSum
    -(m::OperatorProd, ms::OperatorSum) -> OperatorSum
    -(ms::OperatorSum, m::OperatorProd) -> OperatorSum
    -(ms₁::OperatorSum, ms₂::OperatorSum) -> OperatorSum

Overloaded `-` between quantum operators.
"""
@inline Base.:-(m::OperatorProd) = m * (-1)
@inline Base.:-(ms::OperatorSum) = ms * (-1)
@inline Base.:-(factor::Number, m::OperatorProd) = one(m)*factor - m
@inline Base.:-(m::OperatorProd, factor::Number) = m + one(m)*(-factor)
@inline Base.:-(factor::Number, m::Scalar) = replace(m, factor-value(m))
@inline Base.:-(m::Scalar, factor::Number) = replace(m, value(m)-factor)
@inline Base.:-(m₁::Scalar, m₂::Scalar) = replace(m₁, value(m₁)-value(m₂))
@inline Base.:-(factor::Number, ms::OperatorSum) = one(valtype(ms))*factor - ms
@inline Base.:-(ms::OperatorSum, factor::Number) = ms + one(valtype(ms))*(-factor)
@inline function Base.:-(m₁::OperatorProd, m₂::OperatorProd)
    I = promote_type(m₁|>idtype, m₂|>idtype)
    M = promote_type(typeof(m₁), typeof(m₂))
    return OperatorSum{I, M}(m₁, -m₂)
end
@inline function Base.:-(m::OperatorProd, ms::OperatorSum)
    I = promote_type(m|>idtype, ms|>idtype)
    M = promote_type(typeof(m), valtype(ms))
    return sub!(OperatorSum{I, M}(m), ms)
end
@inline function Base.:-(ms::OperatorSum, m::OperatorProd)
    I = promote_type(ms|>idtype, m|>idtype)
    M = promote_type(valtype(ms), typeof(m))
    return sub!(OperatorSum{I, M}(ms), m)
end
@inline function Base.:-(ms₁::OperatorSum, ms₂::OperatorSum)
    I = promote_type(ms₁|>idtype, ms₂|>idtype)
    M = promote_type(valtype(ms₁), valtype(ms₂))
    return sub!(OperatorSum{I, M}(ms₁), ms₂)
end

"""
    /(m::OperatorProd, factor::Number) -> OperatorProd
    /(m::OperatorProd, factor::Scalar) -> OperatorProd
    /(ms::OperatorSum, factor::Number) -> OperatorSum
    /(ms::OperatorSum, factor::Scalar) -> OperatorSum

Overloaded `/` between a quantum operator and a scalar.
"""
Base.:/(m::OperatorProd, factor::Number) = m * (one(valtype(m))/factor)
Base.:/(m::OperatorProd, factor::Scalar) = m * (one(valtype(m))/value(factor))
Base.:/(ms::OperatorSum, factor::Number) = ms * (one(valtype(valtype(ms)))/factor)
Base.:/(ms::OperatorSum, factor::Scalar) = ms * (one(valtype(valtype(ms)))/value(factor))

"""
    //(m::OperatorProd, factor::Integer) -> OperatorProd
    //(m::OperatorProd, factor::Scalar) -> OperatorProd
    //(ms::OperatorSum, factor::Integer) ->  OperatorSum
    //(ms::OperatorSum, factor::Scalar) -> OperatorSum

Overloaded `//` between a quantum operator and a scalar.
"""
Base.://(m::OperatorProd, factor::Integer) = m * (1//factor)
Base.://(m::OperatorProd, factor::Scalar) = m * (1//value(factor))
Base.://(ms::OperatorSum, factor::Number) = ms * (1//factor)
Base.://(ms::OperatorSum, factor::Scalar) = ms * (1//value(factor))

"""
    ^(m::OperatorProd, n::Integer) -> OperatorProd
    ^(ms::OperatorSum, n::Integer) -> OperatorSum

Overloaded `^` between a quantum operator and an integer.
"""
Base.:^(m::OperatorProd, n::Integer) = (@assert n>0 "^ error: non-positive integers are not allowed."; prod(ntuple(i->m, Val(n))))
Base.:^(ms::OperatorSum, n::Integer) = (@assert n>0 "^ error: non-positive integers are not allowed."; prod(ntuple(i->ms, Val(n))))

"""
    ⊗(m::OperatorProd, ms::OperatorSum) -> OperatorSum
    ⊗(ms::OperatorSum, m::OperatorProd) -> OperatorSum
    ⊗(ms₁::OperatorSum, ms₂::OperatorSum) -> OperatorSum

Overloaded `⊗` between quantum operators.
"""
⊗(m::OperatorProd, ms::OperatorSum) = OperatorSum((m ⊗ mm for mm in ms|>values)...)
⊗(ms::OperatorSum, m::OperatorProd) = OperatorSum((mm ⊗ m for mm in ms|>values)...)
⊗(ms₁::OperatorSum, ms₂::OperatorSum) = OperatorSum((m₁ ⊗ m₂ for m₁ in ms₁|>values for m₂ in ms₂|>values)...)

"""
    ⋅(m::OperatorProd, ms::OperatorSum) -> OperatorSum
    ⋅(ms::OperatorSum, m::OperatorProd) -> OperatorSum
    ⋅(ms₁::OperatorSum, ms₂::OperatorSum) -> OperatorSum

Overloaded `⋅` between quantum operators.
"""
⋅(m::OperatorProd, ms::OperatorSum) = OperatorSum((m ⋅ mm for mm in ms|>values)...)
⋅(ms::OperatorSum, m::OperatorProd) = OperatorSum((mm ⋅ m for mm in ms|>values)...)
⋅(ms₁::OperatorSum, ms₂::OperatorSum) = OperatorSum((m₁ ⋅ m₂ for m₁ in ms₁|>values for m₂ in ms₂|>values)...)

"""
    sequence(m::OperatorProd, table) -> NTuple{rank(m), Int}

Get the sequence of the id of a quantum operator according to a table.
"""
@generated sequence(m::OperatorProd, table) = Expr(:tuple, [:(get(table, id(m)[$i], nothing)) for i = 1:rank(m)]...)

"""
    replace(m::OperatorProd, pairs::Pair{<:SingularID, <:Union{OperatorProd, OperatorSum}}...) -> OperatorProd/OperatorSum
    replace(ms::OperatorSum, pairs::Pair{<:SingularID, <:Union{OperatorProd, OperatorSum}}...) -> OperatorSum

Replace the rank-1 components of a quantum operator with new quantum operators.
"""
function Base.replace(m::OperatorProd, pairs::Pair{<:SingularID, <:Union{OperatorProd, OperatorSum}}...)
    replacedids = NTuple{length(pairs), idtype(m)|>eltype}(pair.first for pair in pairs)
    ms = split(m)
    result = ms[1]
    for i = 1:rank(m)
        index = findfirst(isequal(id(m)[i]), replacedids)
        result = result * (isa(index, Int) ? pairs[index].second : ms[i+1])
    end
    return result
end
function Base.replace(ms::OperatorSum, pairs::Pair{<:SingularID, <:Union{OperatorProd, OperatorSum}}...)
    result = operatorsumtype(pairs...)()
    for m in values(ms)
        add!(result, replace(m, pairs...))
    end
    return result
end
@generated function operatorsumtype(ms::Vararg{Pair{<:SingularID, <:Union{OperatorProd, OperatorSum}}, N}) where N
    M = Union{}
    for i = 1:N
        E = fieldtype(ms[i], 2)
        M = promote_type(E<:OperatorSum ? valtype(E) : E, M)
    end
    M = reparameter(M, :id, ID{eltype(idtype(M))})
    return OperatorSum{idtype(M), M}
end

"""
    permute(id₁::SingularID, id₂::SingularID) -> Tuple{Vararg{OperatorProd}}

Permutation rule of two singular ids.
"""
permute(::SingularID, ::SingularID) = error("permute error: not implemented for $(nameof(T)).")

"""
    permute!(result::OperatorSum, m::OperatorProd, table) -> OperatorSum
    permute!(result::OperatorSum, ms::OperatorSum, table) -> OperatorSum

Permute the singular ids of a quantum operator to the descending order according to a table, and store the permuted quantum operators in result.

!!! note
    To use this function, the user must implement a method of `permute`, which computes the result of the permutation of two ids and takes the following interface:
    ```julia
    permute(id₁::SingularID, id₂::SingularID) -> Union{OperatorProd, OperatorSum}
    ```
    Here, `id₁` and `id₂` are two arbitrary singular ids contained in `id(m)`.
"""
function Base.permute!(result::OperatorSum, m::OperatorProd, table)
    cache = valtype(result)[m]
    while length(cache) > 0
        current = pop!(cache)
        pos = operatorprodcommuteposition(sequence(current, table))
        if isa(pos, Nothing)
            add!(result, current)
        else
            left = current[1:pos-1] * value(m)
            right = current[pos+2:end]
            for middle in permute(id(current)[pos], id(current)[pos+1])
                temp = left * middle * right
                (temp === nothing) || push!(cache, temp)
            end
        end
    end
    return result
end
function Base.permute!(result::OperatorSum, ms::OperatorSum, table)
    for m in values(ms)
        permute!(result, m, table)
    end
    return result
end
function operatorprodcommuteposition(sequences)
    pos = 1
    while pos < length(sequences)
        (sequences[pos] < sequences[pos+1]) && return pos
        pos += 1
    end
    return nothing
end

"""
    permute(m::OperatorProd, table) -> OperatorSum
    permute(ms::OperatorSum, table) -> OperatorSum

Permute the singular ids of a quantum operator to the descending order according to a table.
"""
@inline function permute(m::OperatorProd, table)
    M = reparameter(typeof(m), :id, ID{m|>idtype|>eltype})
    permute!(OperatorSum{idtype(M), M}(), m, table)
end
@inline function permute(ms::OperatorSum, table)
    M = reparameter(valtype(ms), :id, ID{ms|>idtype|>eltype})
    permute!(OperatorSum{idtype(M), M}(), ms, table)
end

"""
    Transformation <: Function

Abstract transformation that could transform the quantum operators.
"""
abstract type Transformation <: Function end
@inline Base.:(==)(transformation₁::Transformation, transformation₂::Transformation) = ==(efficientoperations, transformation₁, transformation₂)
@inline Base.isequal(transformation₁::Transformation, transformation₂::Transformation) = isequal(efficientoperations, transformation₁, transformation₂)
@inline Base.valtype(transformation::Transformation, m::OperatorProd) = valtype(typeof(transformation), typeof(m))
@inline Base.valtype(transformation::Transformation, ms::OperatorSum) = valtype(typeof(transformation), typeof(ms))
@inline Base.valtype(T::Type{<:Transformation}, MS::Type{<:OperatorSum}) = OperatorSum{idtype(valtype(T, valtype(MS))), valtype(T, valtype(MS))}

"""
    (transformation::Transformation)(ms::OperatorSum) -> OperatorSum

Get the transformed quantum operators.
"""
function (transformation::Transformation)(ms::OperatorSum)
    result = valtype(transformation, ms)()
    for m in values(ms)
        add!(result, transformation(m))
    end
    return result
end

"""
    Identity <: Transformation

The identity transformation.
"""
struct Identity <: Transformation end
@inline Base.valtype(::Type{Identity}, M::Type{<:OperatorProd}) = M
@inline (i::Identity)(m::OperatorProd) = m

"""
    Numericalization{T<:Number} <: Transformation

The numericalization transformation, which converts the value of a quantum operator to a number of type `T`.
"""
struct Numericalization{T<:Number} <: Transformation end
@inline Base.valtype(num::Numericalization) = valtype(typeof(num))
@inline Base.valtype(::Type{<:Numericalization{T}}) where {T<:Number} = T
@inline Base.valtype(T::Type{<:Numericalization}, M::Type{<:OperatorProd}) = reparameter(M, :value, valtype(T))
@inline (n::Numericalization)(m::OperatorProd) = convert(valtype(n, m), m)

end #module
