module AlgebraOverFields

using Printf: @printf, @sprintf
using ...Prerequisites: atol, rtol
using ...Prerequisites.NamedVectors: NamedVector
using ...Prerequisites.TypeTraits: efficientoperations, rawtype, fulltype, Field, Parameter, parametertype, parameterpairs

import ...Interfaces: rank, scalar, id, add!, sub!, mul!, div!, ⊗, ⋅, sequence, permute
import ...Prerequisites.TypeTraits: isparameterbound, parameternames

export SimpleID, ID
export Element, Scalar, Elements
export idtype

"""
    SimpleID <: NamedVector

A simple id is the building block of the id system of an algebra over a field.
"""
abstract type SimpleID <: NamedVector end

"""
    ID{I<:SimpleID, N}

The (composite) id system of an algebra over a field.

Type alias for `NTuple{N, I} where {N, I<:SimpleID}`.
"""
const ID{I<:SimpleID, N} = NTuple{N, I}

"""
    ID(ids::SimpleID...)
    ID(ids::NTuple{N, SimpleID}) where N

Get the composite id from simple ids.
"""
ID(ids::SimpleID...) = ids
ID(ids::NTuple{N, SimpleID}) where {N} = ids

"""
    ID(::Type{SID}, attrs::Vararg{NTuple{N}, M}) where {SID<:SimpleID, N, M}

Get the composite id from the components of simple ids.
"""
@generated function ID(::Type{SID}, attrs::Vararg{NTuple{N, Any}, M}) where {SID<:SimpleID, N, M}
    exprs = []
    for i = 1:N
        args = [:(attrs[$j][$i]) for j = 1:M]
        push!(exprs, :(SID($(args...))))
    end
    return :(ID($(exprs...)))
end

"""
    propertynames(::Type{I}) where I<:ID{SimpleID} -> Tuple{Vararg{Symbol}}

Get the property names of a composite id.
"""
@generated function Base.propertynames(I::ID{SimpleID})
    exprs = [QuoteNode(Symbol(name, 's')) for name in I|>eltype|>fieldnames]
    return Expr(:tuple, exprs...)
end

"""
    getproperty(cid::ID{SimpleID}, name::Symbol)

Get the property of a composite id.
"""
Base.getproperty(cid::ID{SimpleID}, name::Symbol) = idgetproperty(cid, Val(name), Val(cid|>propertynames))
@generated function idgetproperty(cid::ID{SimpleID}, ::Val{name}, ::Val{names}) where {name, names}
    index = findfirst(isequal(name), names)::Int
    exprs = [:(getfield(cid[$i], $index)) for i = 1:fieldcount(cid)]
    return Expr(:tuple, exprs...)
end

"""
    show(io::IO, cid::Tuple{SimpleID, Vararg{SimpleID}})

Show a composite id.
"""
Base.show(io::IO, cid::Tuple{SimpleID, Vararg{SimpleID}}) = @printf io "ID(%s)" join(cid, ", ")

"""
    promote_type(::Type{Tuple{}}, I::Type{<:ID{SimpleID, N}}) where N
    promote_type(I::Type{<:ID{SimpleID, N}}, ::Type{Tuple{}}) where N

Define the promote rule for ID types.
"""
Base.promote_type(::Type{Tuple{}}, I::Type{<:Tuple{SimpleID, Vararg{SimpleID}}}) = ID{I|>eltype}
Base.promote_type(I::Type{<:Tuple{SimpleID, Vararg{SimpleID}}}, ::Type{Tuple{}}) = ID{I|>eltype}

"""
    isless(::Type{<:SimpleID}, cid1::ID{SimpleID}, cid2::ID{SimpleID}) -> Bool

Compare two ids and judge whether the first is less than the second.

The comparison rule are as follows:
1. ids with smaller ranks are always less than those with higher ranks;
2. if two ids are of the same rank, the comparison goes just like that between tuples.
"""
function Base.isless(::Type{<:SimpleID}, cid1::ID{SimpleID}, cid2::ID{SimpleID})
    r1, r2 = cid1|>rank, cid2|>rank
    (r1 < r2) ? true : (r1 > r2) ? false : isless(cid1, cid2)
end

"""
    rank(id::ID{SimpleID}) -> Int
    rank(::Type{<:ID{SimpleID}}) -> Any
    rank(::Type{<:ID{SimpleID, N}}) where N -> Int

Get the rank of a composite id.
"""
rank(id::ID{SimpleID}) = id |> typeof |> rank
rank(::Type{<:ID{SimpleID}}) = Any
rank(::Type{<:ID{SimpleID, N}}) where {N} = N

"""
    *(sid1::SimpleID, sid2::SimpleID) -> ID{SimpleID}
    *(sid::SimpleID, cid::ID{SimpleID}) -> ID{SimpleID}
    *(cid::ID{SimpleID}, sid::SimpleID) -> ID{SimpleID}
    *(cid1::ID{SimpleID}, cid2::ID{SimpleID}) -> ID{SimpleID}

Get the product of the id system.
"""
Base.:*(sid1::SimpleID, sid2::SimpleID) = ID(sid1, sid2)
Base.:*(sid::SimpleID, cid::ID{SimpleID}) = ID(sid, cid...)
Base.:*(cid::ID{SimpleID}, sid::SimpleID) = ID(cid..., sid)
Base.:*(cid1::ID{SimpleID}, cid2::ID{SimpleID}) = ID(cid1..., cid2...)

"""
    Element{V, I<:ID{SimpleID}}

An element of an algebra over a field.

Basically, a concrete subtype should contain two attributes:
- `value::V`: the coefficient of the element
- `id::I`: the id of the element
"""
abstract type Element{V, I<:ID{SimpleID}} end
Base.fieldnames(::Type{Field}, ::Type{<:Element}) = (:value, :id)
parameternames(::Type{<:Element}) = (:value, :id)
isparameterbound(::Type{<:Element}, ::Parameter{:value}, ::Type{V}) where V = false
isparameterbound(::Type{<:Element}, ::Parameter{:id}, ::Type{I}) where I<:ID{SimpleID} = !isconcretetype(I) 

"""
    replace(m::Element, v::Number) -> Element

Replace the value of an element.
"""
Base.replace(m::Element, v::Number) = _replace(m, v, fieldname(typeof(m), Field(:value))|>Val)
@generated _replace(m::Element, v::Number, ::Val{name}) where name = Expr(:call, :replace, :m, Expr(:kw, name, :v))

"""
    scalar(m::Element) -> valtype(m)

Get the value of an element.
"""
scalar(m::Element) = getproperty(m, Field(:value))

"""
    id(m::Element) -> idtype(m)

Get the id of an element.
"""
id(m::Element) = getproperty(m, Field(:id))

"""
    Scalar{V}

Scalar element.
"""
const Scalar{V} = Element{V, Tuple{}}

"""
    valtype(m::Element)
    valtype(::Type{T}) where {T<:Element}

Get the type of the value of an element.

The result is also the type of the field over which the algebra is defined.
"""
Base.valtype(m::Element) = m |> typeof |> valtype
@generated Base.valtype(::Type{T}) where {T<:Element} = parametertype(supertype(T, :Element), 1)

"""
    idtype(m::Element)
    idtype(::Type{T}) where {T<:Element}

The type of the id of an element.
"""
idtype(m::Element) = m |> typeof |> idtype
@generated idtype(::Type{T}) where {T<:Element} = parametertype(supertype(T, :Element), 2)

"""
    rank(m::Element) -> Int
    rank(::Type{M}) where M<:Element -> Int

Get the rank of an element.
"""
rank(m::Element) = m |> typeof |> rank
rank(::Type{M}) where {M<:Element} = M |> idtype |> rank

"""
    ==(m1::Element, m2::Element) -> Bool

Compare two elements and judge whether they are equal to each other.
"""
Base.:(==)(m1::Element, m2::Element) = ==(efficientoperations, m1, m2)

"""
    isequal(m1::Element, m2::Element) -> Bool

Compare two elements and judge whether they are equal to each other.
"""
Base.isequal(m1::Element, m2::Element) = isequal(efficientoperations, m1, m2)

"""
    isapprox(m1::Element, m2::Element; atol::Real=atol, rtol::Real=rtol) -> Bool

Compare two elements and judge whether they are inexactly equivalent to each other.
"""
Base.isapprox(m1::Element, m2::Element; atol::Real=atol, rtol::Real=rtol) = isapprox(efficientoperations, (fieldname(typeof(m1), Field(:value)),)|>Val, m1, m2; atol=atol, rtol=rtol)

"""
    replace(m::Element; kwargs...) -> typeof(m)

Return a copy of a concrete `Element` with some of the field values replaced by the keyword arguments.
"""
Base.replace(m::Element; kwargs...) = replace(efficientoperations, m; kwargs...)

"""
    promote_rule(::Type{M1}, ::Type{M2}) where {M1<:Element, M2<:Element}

Define the promote rule for Element types.
"""
function Base.promote_rule(::Type{M1}, ::Type{M2}) where {M1<:Element, M2<:Element}
    (M1 <: M2) && return M2
    (M2 <: M1) && return M1
    r1, r2 = M1|>rank, M2|>rank
    M = (r2 == 0) ? rawtype(M1) : (r1 == 0) ? rawtype(M2) : typejoin(rawtype(M1), rawtype(M2))
    return fulltype(M, promote_type(Parameter, parameterpairs(M1), parameterpairs(M2)))
end

"""
    promote_type(::Type{M}, ::Type{V}, ::Val{:*}) where {M<:Element, V<:Number}
    promote_type(::Type{V}, ::Type{M}, ::Val{:*}) where {M<:Element, V<:Number}

Define the promote rule for the multiplication between an Element and a scalar.
"""
Base.promote_type(::Type{M}, ::Type{V}, ::Val{:*}) where {M<:Element, V<:Number} = promote_type(V, M, Val(:*))
function Base.promote_type(::Type{V}, ::Type{M}, ::Val{:*}) where {M<:Element, V<:Number}
    return replace(M, Parameter(:value), promote_type(valtype(M), V))
end

"""
    getindex(m::Element, i) -> Element

Overloaded `[]` operator.
"""
Base.getindex(m::Element, i) = rawtype(typeof(m))(Tuple(Field, m, getindex, (i,))...)
Base.map(::typeof(getindex), m::Element, ::Field{:value}, ::Tuple{Any}, ::NamedTuple) = one(valtype(m))
Base.map(::typeof(getindex), m::Element, ::Field{:id}, i::Tuple{Any}, ::NamedTuple) = ID(id(m)[i[1]])

"""
    length(m::Element) -> Int

Get the length of an element.
"""
Base.length(m::Element) = rank(m)
Base.firstindex(m::Element) = 1
Base.lastindex(m::Element) = rank(m)

"""
    one(::Type{M}) where M<:Element
    one(m::Element)

Get the identity operator.
"""
function Base.one(::Type{M}) where M<:Element
    rtype = rawtype(M)
    vtype = isconcretetype(valtype(M)) ? valtype(M) : Int
    @assert fieldnames(rtype) == (:value, :id) "one error: not supproted type($(nameof(rtype)))."
    return rtype(one(vtype), ID())
end
Base.one(m::Element) = rawtype(typeof(m))(Tuple(Field, m, one)...)
Base.map(::typeof(one), m::Element, ::Field{:value}, ::Tuple, ::NamedTuple) = one(valtype(m))
Base.map(::typeof(one), ::Element, ::Field{:id}, ::Tuple, ::NamedTuple) = ID()

"""
    convert(::Type{M}, m::Scalar) where M<:Scalar
    convert(::Type{M}, m::Number) where M<:Scalar
    convert(::Type{M}, m::Element) where M<:Element

1) Convert a scalar element from one type to another;
2) Convert a scalar to a scalar element;
3) Convert an element from one type to another.
"""
Base.convert(::Type{M}, m::Scalar) where M<:Scalar = (typeof(m) <: M) ? m : one(M)*scalar(m)
Base.convert(::Type{M}, m::Number) where M<:Scalar = one(M)*m
function Base.convert(::Type{M}, m::Element) where M<:Element
    (typeof(m) <: M) && return m
    @assert convertable(M, m) "convert error: $(nameof(typeof(m))) cannot be converted to $(nameof(M))."
    return replace(m, convert(valtype(M), scalar(m)))
end
function convertable(::Type{M}, m::Element) where M<:Element
    !(rawtype(typeof(m)) <: rawtype(M)) && return false
    S = supertype(typeof(m), nameof(M))
    for (i, name) in enumerate(parameternames(M))
        (name != :value) && !(parametertype(S, i) <: parametertype(M, i)) && return false
    end
    return true
end

"""
    split(m::Element) -> Tuple{Any, Vararg{Element}}

Split an element into the coefficient and a sequence of rank-1 elements.
"""
@generated function Base.split(m::Element)
    exprs = [:(scalar(m))]
    for i = 1:rank(m)
        push!(exprs, :(m[$i]))
    end
    return Expr(:tuple, exprs...)
end

"""
    Elements{I<:ID{SimpleID}, M<:Element} <: AbstractDict{I, M}

An set of elements of an algebra over a field.

Type alias for `Dict{I<:ID{SimpleID}, M<:Element}`.
Similar iterms are automatically merged thanks to the id system.
"""
const Elements{I<:ID{SimpleID}, M<:Element} = Dict{I, M}

"""
    Elements(ms)
    Elements(ms::Pair{I, M}...) where {I<:ID{SimpleID}, M<:Element}
    Elements(ms::Element...)

Get the set of elements with similar items merged.
"""
Elements(ms) = Base.dict_with_eltype((K, V)->Dict{K, V}, ms, eltype(ms))
function Elements(ms::Pair{I, M}...) where {I<:ID{SimpleID}, M<:Element}
    result = Elements{I, M}()
    for (id, m) in ms
        result[id] = m
    end
    return result
end
function Elements(ms::Element...)
    result = Elements{ms|>eltype|>idtype, ms|>eltype}()
    for m in ms add!(result, m) end
    return result
end

"""
    show(io::IO, ms::Elements)

Show a set of elements.
"""
function Base.show(io::IO, ms::Elements)
    @printf io "Elements with %s entries:\n" length(ms)
    for m in values(ms)
        @printf io "  %s\n" m
    end
end

"""
    repr(ms::Elements) -> String

Get the repr representation of a set of elements.
"""
function Base.repr(ms::Elements)
    cache = [@sprintf("Elements with %s entries:", length(ms))]
    for m in ms|>values
        push!(cache, @sprintf("  %s", repr(m)))
    end
    return join(cache, "\n")
end

"""
    zero(ms::Elements) -> Nothing
    zero(::Type{<:Elements}) -> Nothing

Get a zero set of elements.

A zero set of elements is defined to be the one with no elements.
"""
Base.zero(ms::Elements) = ms |> typeof |> zero
Base.zero(::Type{<:Elements}) = nothing

"""
    ==(ms::Elements, ::Nothing) -> Bool
    ==(::Nothing, ms::Elements) -> Bool

Judge whether a set of elements is identically empty.
"""
Base.:(==)(ms::Elements, ::Nothing) = length(ms) == 0
Base.:(==)(::Nothing, ms::Elements) = length(ms) == 0

"""
    isequal(ms::Elements, ::Nothing) -> Bool
    isequal(::Nothing, ms::Elements) -> Bool

Judge whether a set of elements is identically empty.
"""
Base.isequal(ms::Elements, ::Nothing) = length(ms) == 0
Base.isequal(::Nothing, ms::Elements) = length(ms) == 0

"""
    add!(ms::Elements) -> typeof(ms)
    add!(ms::Elements, ::Nothing) -> typeof(ms)
    add!(ms::Elements, m::Number) -> typeof(ms)
    add!(ms::Elements, m::Element) -> typeof(ms)
    add!(ms::Elements, mms::Elements) -> typeof(ms)

Get the inplace addition of elements to a set.
"""
add!(ms::Elements) = ms
add!(ms::Elements, ::Nothing) = ms
add!(ms::Elements, m::Number) = add!(ms, one(valtype(ms))*m)
function add!(ms::Elements, m::Element)
    m = convert(ms|>valtype, m)
    old = get(ms, id(m), nothing)
    new = (old === nothing) ? m : replace(old, scalar(old)+scalar(m))
    (abs(scalar(new)) == 0.0) ? delete!(ms, id(m)) : (ms[id(m)] = new)
    return ms
end
add!(ms::Elements, mms::Elements) = (for m in mms|>values add!(ms, m) end; ms)

"""
    sub!(ms::Elements) -> typeof(ms)
    sub!(ms::Elements, ::Nothing) -> typeof(ms)
    add!(ms::Elements, m::Number) -> typeof(ms)
    sub!(ms::Elements, m::Element) -> typeof(ms)
    sub!(ms::Elements, mms::Elements) -> typeof(ms)

Get the inplace subtraction of elements from a set.
"""
sub!(ms::Elements) = ms
sub!(ms::Elements, ::Nothing) = ms
sub!(ms::Elements, m::Number) = add!(ms, one(valtype(ms))*(-m))
function sub!(ms::Elements, m::Element)
    m = convert(ms|>valtype, m)
    old = get(ms, id(m), nothing)
    new = (old === nothing) ? -m : replace(old, scalar(old)-scalar(m))
    (abs(scalar(new)) == 0.0) ? delete!(ms, id(m)) : (ms[id(m)] = new)
    return ms
end
sub!(ms::Elements, mms::Elements) = (for m in mms|>values sub!(ms, m) end; ms)

"""
    mul!(ms::Elements, factor::Scalar) -> Elements
    mul!(ms::Elements, factor::Number) -> Elements

Get the inplace multiplication of elements with a scalar.
"""
mul!(ms::Elements, factor::Scalar) = mul!(ms, scalar(factor))
function mul!(ms::Elements, factor::Number)
    @assert isa(one(ms|>valtype|>valtype)*factor, ms|>valtype|>valtype) "mul! error: dismatched type, $(ms|>valtype) and $(factor|>typeof)."
    for m in values(ms)
        ms[id(m)] = replace(m, scalar(m)*factor)
    end
    return ms
end

"""
    div!(ms::Elements, factor::Scalar) -> Elements
    div!(ms::Elements, factor::Number) -> Elements

Get the inplace division of element with a scalar.
"""
div!(ms::Elements, factor::Scalar) = div!(ms, scalar(factor))
function div!(ms::Elements, factor::Number)
    @assert isa(one(ms|>valtype|>valtype)/factor, ms|>valtype|>valtype) "div! error: dismatched type, $(ms|>valtype) and $(factor|>typeof)."
    for m in values(ms)
        ms[id(m)] = replace(m, scalar(m)/factor)
    end
    return ms
end

"""
    +(m::Element) -> typeof(m)
    +(m::Element, ::Nothing) -> typeof(m)
    +(::Nothing, m::Element) -> typeof(m)
    +(m::Element, factor::Number) -> Elements
    +(factor::Number, m::Element) -> Elements
    +(m1::Element, m2::Element) -> Elements
    +(ms::Elements) -> typeof(ms)
    +(ms::Elements, ::Nothing) -> typeof(ms)
    +(::Nothing, ms::Elements) -> typeof(ms)
    +(ms::Elements, factor::Number) -> Elements
    +(factor::Number, ms::Elements) -> Elements
    +(ms::Elements, m::Element) -> Elements
    +(m::Element, ms::Elements) -> Elements
    +(ms1::Elements, ms2::Elements) -> Elements

Overloaded `+` operator between elements of an algebra over a field.
"""
Base.:+(m::Element) = m
Base.:+(ms::Elements) = ms
Base.:+(m::Element, ::Nothing) = m
Base.:+(::Nothing, m::Element) = m
Base.:+(ms::Elements, ::Nothing) = ms
Base.:+(::Nothing, ms::Elements) = ms
Base.:+(factor::Number, m::Element) = m + one(m)*factor
Base.:+(m::Element, factor::Number) = m + one(m)*factor
Base.:+(factor::Number, m::Scalar) = replace(m, scalar(m)+factor)
Base.:+(m::Scalar, factor::Number) = replace(m, scalar(m)+factor)
Base.:+(m1::Scalar, m2::Scalar) = replace(m1, scalar(m1)+scalar(m2))
Base.:+(factor::Number, ms::Elements) = ms + one(valtype(ms))*factor
Base.:+(ms::Elements, factor::Number) = ms + one(valtype(ms))*factor
function Base.:+(m1::Element, m2::Element)
    I = promote_type(m1|>idtype, m2|>idtype)
    M = promote_type(typeof(m1), typeof(m2))
    return add!(Elements{I, M}(id(m1)=>m1), m2)
end
Base.:+(m::Element, ms::Elements) = ms + m
function Base.:+(ms::Elements, m::Element)
    I = promote_type(ms|>keytype, m|>idtype)
    M = promote_type(valtype(ms), typeof(m))
    return add!(Elements{I, M}(ms), m)
end
function Base.:+(ms1::Elements, ms2::Elements)
    I = promote_type(ms1|>keytype, ms2|>keytype)
    M = promote_type(valtype(ms1), valtype(ms2))
    return add!(Elements{I, M}(ms1), ms2)
end

"""
    *(m::Element, ::Nothing) -> Nothing
    *(::Nothing, m::Element) -> Nothing
    *(factor::Number, m::Element) -> Element
    *(m::Element, factor::Number) -> Element
    *(m1::Element, m2::Element) -> Element
    *(ms::Elements, ::Nothing) -> Nothing
    *(::Nothing, ms::Elements) -> Nothing
    *(factor::Number, ms::Elements) -> Elements
    *(ms::Elements, factor::Number) -> Elements
    *(m::Element, ms::Elements) -> Elements
    *(ms::Elements, m::Element) -> Elements
    *(ms1::Elements, ms2::Elements) -> Elements

Overloaded `*` operator for element-scalar multiplications and element-element multiplications of an algebra over a field.
"""
Base.:*(m::Element, ::Nothing) = nothing
Base.:*(::Nothing, m::Element) = nothing
Base.:*(ms::Elements, ::Nothing) = nothing
Base.:*(::Nothing, ms::Elements) = nothing
Base.:*(factor::Scalar, m::Element) = m * scalar(factor)
Base.:*(m::Element, factor::Scalar) = m * scalar(factor)
Base.:*(m1::Scalar, m2::Scalar) = replace(m1, scalar(m1)*scalar(m2))
Base.:*(factor::Number, m::Element) = m * factor
Base.:*(m::Element, factor::Number) = replace(m, factor*scalar(m))
Base.:*(factor::Scalar, ms::Elements) = ms * scalar(factor)
Base.:*(ms::Elements, factor::Scalar) = ms * scalar(factor)
Base.:*(factor::Number, ms::Elements) = ms * factor
function Base.:*(ms::Elements, factor::Number)
    (abs(factor) == 0) && return zero(Elements)
    result = Elements{keytype(ms), promote_type(valtype(ms), typeof(factor), Val(:*))}()
    for (id, m) in ms
        result[id] = m * factor
    end
    return result
end
Base.:*(m::Element, ms::Elements) = Elements((m * mm for mm in ms|>values)...)
Base.:*(ms::Elements, m::Element) = Elements((mm * m for mm in ms|>values)...)
Base.:*(ms1::Elements, ms2::Elements) = Elements((m1 * m2 for m1 in ms1|>values for m2 in ms2|>values)...)
function Base.:*(m1::Element, m2::Element)
    @assert((m1|>typeof|>nameof == m2|>typeof|>nameof) && (m1|>typeof|>fieldcount == m2|>typeof|>fieldcount == 2),
            "\"*\" error: not implemented between $(m1|>typeof|>nameof) and $(m2|>typeof|>nameof)."
            )
    rawtype(typeof(m1))(scalar(m1)*scalar(m2), id(m1)*id(m2))
end

"""
    -(m::Element) -> typeof(m)
    -(m::Element, ::Nothing) -> typeof(m)
    -(::Nothing, m::Element) -> typeof(m)
    -(m::Element, factor::Number) -> Elements
    -(factor::Number, m::Element) -> Elements
    -(m1::Element, m2::Element) -> Elements
    -(ms::Elements) -> typeof(ms)
    -(ms::Elements, ::Nothing) -> typeof(ms)
    -(::Nothing, ms::Elements) -> typeof(ms)
    -(ms::Elements, factor::Number) -> Elements
    -(factor::Number, ms::Elements) -> Elements
    -(m::Element, ms::Elements) -> Elements
    -(ms::Elements, m::Element) -> Elements
    -(ms1::Elements, ms2::Elements) -> Elements

Overloaded `-` operator between elements of an algebra over a field.
"""
Base.:-(m::Element) = m * (-1)
Base.:-(ms::Elements) = ms * (-1)
Base.:-(m::Element, ::Nothing) = m
Base.:-(::Nothing, m::Element) = -m
Base.:-(ms::Elements, ::Nothing) = ms
Base.:-(::Nothing, ms::Elements) = -ms
Base.:-(factor::Number, m::Element) = one(m)*factor - m
Base.:-(m::Element, factor::Number) = m + one(m)*(-factor)
Base.:-(factor::Number, m::Scalar) = replace(m, factor-scalar(m))
Base.:-(m::Scalar, factor::Number) = replace(m, scalar(m)-factor)
Base.:-(m1::Scalar, m2::Scalar) = replace(m1, scalar(m1)-scalar(m2))
Base.:-(factor::Number, ms::Elements) = one(valtype(ms))*factor - ms
Base.:-(ms::Elements, factor::Number) = ms + one(valtype(ms))*(-factor)
function Base.:-(m1::Element, m2::Element)
    I = promote_type(m1|>idtype, m2|>idtype)
    M = promote_type(typeof(m1), typeof(m2))
    return sub!(Elements{I, M}(id(m1)=>m1), m2)
end
function Base.:-(m::Element, ms::Elements)
    I = promote_type(m|>idtype, ms|>keytype)
    M = promote_type(typeof(m), valtype(ms))
    return sub!(Elements{I, M}(id(m)=>m), ms)
end
function Base.:-(ms::Elements, m::Element)
    I = promote_type(ms|>keytype, m|>idtype)
    M = promote_type(valtype(ms), typeof(m))
    return sub!(Elements{M|>idtype, M}(ms), m)
end
function Base.:-(ms1::Elements, ms2::Elements)
    I = promote_type(ms1|>keytype, ms2|>keytype)
    M = promote_type(valtype(ms1), valtype(ms2))
    return sub!(Elements{I, M}(ms1), ms2)
end

"""
    /(m::Element, factor::Number) -> Element
    /(m::Element, factor::Scalar) -> Element
    /(ms::Elements, factor::Number) -> Elements
    /(ms::Elements, factor::Scalar) -> Elements

Overloaded `/` operator for element-scalar division of an algebra over a field.
"""
Base.:/(m::Element, factor::Number) = m * (1/factor)
Base.:/(m::Element, factor::Scalar) = m * (1/scalar(factor))
Base.:/(ms::Elements, factor::Number) = ms * (1/factor)
Base.:/(ms::Elements, factor::Scalar) = ms * (1/scalar(factor))

"""
    //(m::Element, factor::Integer) -> Element
    //(m::Element, factor::Scalar) -> Element
    //(ms::Elements, factor::Integer) ->  Elements
    //(ms::Elements, factor::Scalar) -> Elements

Overloaded `//` operator for element-scalar division of an algebra over a field.
"""
Base.://(m::Element, factor::Integer) = m * (1//factor)
Base.://(m::Element, factor::Scalar) = m * (1//scalar(factor))
Base.://(ms::Elements, factor::Number) = ms * (1//factor)
Base.://(ms::Elements, factor::Scalar) = ms * (1//scalar(factor))

"""
    ^(m::Element, n::Integer) -> Element
    ^(ms::Elements, n::Integer) -> Elements

Overloaded `^` operator for element-integer power of an algebra over a field.
"""
Base.:^(m::Element, n::Integer) = (@assert n > 0 "^ error: non-positive integers are not allowed."; prod(ntuple(i->m, Val(n))))
Base.:^(ms::Elements, n::Integer) = (@assert n > 0 "^ error: non-positive integers are not allowed."; prod(ntuple(i->ms, Val(n))))

"""
    ⊗(m::Element, ms::Elements) -> Elements
    ⊗(ms::Elements, m::Element) -> Elements
    ⊗(ms1::Elements, ms2::Elements) -> Elements

Overloaded `⊗` operator for element-element multiplications of an algebra over a field.
"""
⊗(m::Element, ms::Elements) = Elements((m ⊗ mm for mm in ms|>values)...)
⊗(ms::Elements, m::Element) = Elements((mm ⊗ m for mm in ms|>values)...)
⊗(ms1::Elements, ms2::Elements) = Elements((m1 ⊗ m2 for m1 in ms1|>values for m2 in ms2|>values)...)

"""
    ⋅(m::Element, ms::Elements) -> Elements
    ⋅(ms::Elements, m::Element) -> Elements
    ⋅(ms1::Elements, ms2::Elements) -> Elements

Overloaded `⋅` operator for element-element multiplications of an algebra over a field.
"""
⋅(m::Element, ms::Elements) = Elements((m ⋅ mm for mm in ms|>values)...)
⋅(ms::Elements, m::Element) = Elements((mm ⋅ m for mm in ms|>values)...)
⋅(ms1::Elements, ms2::Elements) = Elements((m1 ⋅ m2 for m1 in ms1|>values for m2 in ms2|>values)...)

"""
    sequence(m::Element, table=nothing) -> NTuple{rank(m), Int}

Get the sequence of the ids of an element according to a table.
"""
@generated sequence(m::Element, table) = Expr(:tuple, [:(get(table, id(m)[$i], nothing)) for i = 1:rank(m)]...)

"""
    replace(m::Element, pairs::Pair{<:SimpleID, <:Union{Element, Elements}}...) -> Element/Elements
    replace(ms::Elements, pairs::Pair{<:SimpleID, <:Union{Element, Elements}}...) -> Elements

Replace the rank-1 components of an element with new element/elements.
"""
function Base.replace(m::Element, pairs::Pair{<:SimpleID, <:Union{Element, Elements}}...)
    replacedids = NTuple{length(pairs), idtype(m)|>eltype}(pair.first for pair in pairs)
    ms = split(m)
    result = ms[1]
    for i = 1:rank(m)
        index = findfirst(isequal(id(m)[i]), replacedids)
        result = result * (isa(index, Int) ? pairs[index].second : ms[i+1])
    end
    return result
end
function Base.replace(ms::Elements, pairs::Pair{<:SimpleID, <:Union{Element, Elements}}...)
    result = elementstype(pairs...)()
    for m in values(ms)
        add!(result, replace(m, pairs...))
    end
    return result
end
@generated function elementstype(ms::Vararg{Pair{<:SimpleID, <:Union{Element, Elements}}, N}) where N
    M = Union{}
    for i = 1:N
        E = fieldtype(ms[i], 2)
        M = promote_type(E<:Elements ? valtype(E) : E, M)
    end
    M = replace(M, Parameter(:id), ID{eltype(idtype(M))})
    return Elements{idtype(M), M}
end

"""
    permute(::Type{T}, id1::SimpleID, id2::SimpleID, table::Any) where T<:Element -> Tuple{Vararg{Element}}

Permutation rule of two ids of type `T` according to an extra table.
"""
permute(::Type{T}, ::SimpleID, ::SimpleID, ::Any) where T<:Element = error("permute error: not implemented for $(nameof(T)).")

"""
    permute!(result::Elements, m::Element, table=nothing) -> Elements
    permute!(result::Elements, ms::Elements, table=nothing) -> Elements

Permute the ids of an-element/a-set-of-elements to the descending order according to a table, and store the permuted elements in result.

!!! note
    To use this function, the user must implement a method of `permute`, which computes the result of the permutation of a rank-2 element and takes the following interface:
    ```julia
    permute(::Type{M}, id1::SimpleID, id2::SimpleID, table=nothing) -> Union{M, Elements}
    ```
    Here, `M` is the type of the input element `m` in `permute!`, `id1` and `id2` are two arbitary simple ids contained in `id(m)`.
"""
function Base.permute!(result::Elements, m::Element, table=nothing)
    cache = valtype(result)[m]
    while length(cache) > 0
        current = pop!(cache)
        pos = elementcommuteposition(sequence(current, table))
        if isa(pos, Nothing)
            add!(result, current)
        else
            left = current[1:pos-1] * scalar(m)
            right = current[pos+2:end]
            for middle in permute(typeof(current), id(current)[pos], id(current)[pos+1], table)
                temp = left * middle * right
                (temp === nothing) || push!(cache, temp)
            end
        end
    end
    return result
end
function Base.permute!(result::Elements, ms::Elements, table=nothing)
    for m in values(ms)
        permute!(result, m, table)
    end
    return result
end
function elementcommuteposition(seqs)
    pos = 1
    while pos < length(seqs)
        (seqs[pos] < seqs[pos+1]) && return pos
        pos += 1
    end
    return nothing
end

"""
    permute(m::Element, table=nothing) -> Elements
    permute(ms::Elements, table=nothing) -> Elements

Permute the ids of an-element/a-set-of-elements to the descending order according to a table.
"""
function permute(m::Element, table=nothing)
    M = replace(typeof(m), Parameter(:id), ID{m|>idtype|>eltype})
    permute!(Elements{idtype(M), M}(), m, table)
end
function permute(ms::Elements, table=nothing)
    M = replace(valtype(ms), Parameter(:id), ID{ms|>keytype|>eltype})
    permute!(Elements{idtype(M), M}(), ms, table)
end

end #module
