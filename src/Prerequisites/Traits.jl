module Traits

using ..Prerequisites: atol, rtol

export parametercount, parametername, parameterorder, parametertype, parameterpair, isparameterbound, hasparameter
export parameternames, parametertypes, parameterpairs, isparameterbounds, reparameter, promoteparameters, rawtype, fulltype
export contentcount, contentname, contentorder, hascontent, getcontent, contentnames, dissolve
export efficientoperations

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
    @assert parametercount(T) == fieldcount(PS) == length(ubs) "fulltype error: length-dismatched input parameters."
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
    dissolve(m, f::Function=identity, args::Tuple=(), kwargs::NamedTuple=NamedTuple()) -> Tuple

Convert `m` to a tuple by the function `f` applied elementally to its contents with the extra positional arguments (`args`) and keyword arguments (`kwargs`). 

The underlying called interface is the `dissolve` function when `f` is applied to each content of `m`:
```julia
dissolve(m, Val(name), f, args, kwargs)
```
Here, `name` is the name of a content of `m`.

Basically, the rule of how `f` operates on each field of `m` can be overriden by redefining the above `dissolve` function.
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

Disolve the content specified by `name` of `m` by the function `f` applied with the extra positional arguments (`args`) and keyword arguments (`kwargs`).
"""
@inline dissolve(m, ::Val{name}, f::Function, args::Tuple, kwargs::NamedTuple) where name = getcontent(m, name|>Val)

struct EfficientOperations end
"""
    efficientoperations

Indicate that the efficient operations, i.e. "=="/"isequal", "<"/"isless" or "replace", will be used.
"""
const efficientoperations = EfficientOperations()

"""
    ==(::EfficientOperations, o1, o2) -> Bool

Compare two objects and judge whether they are eqaul to each other.
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

Compare two objects and judge whether they are eqaul to each other.
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
    isapprox(::EfficientOperations, fields::Union{Union{Integer, Symbol}, Tuple{Vararg{Union{Integer, Symbol}}}}, o1, o2; atol=atol, rtol=rtol) -> Bool
    isapprox(::EfficientOperations, ::Val{fields}, o1, o2; atol=atol, rtol=rtol) where fields -> Bool

Compare two objects and judge whether they are inexactly equivalent to each other.
"""
@inline function Base.isapprox(::EfficientOperations, fields::Union{Union{Integer, Symbol}, Tuple{Vararg{Union{Integer, Symbol}}}}, o1, o2; atol=atol, rtol=rtol)
    isapprox(efficientoperations, fields|>Val, o1, o2; atol=atol, rtol=rtol)
end
@inline @generated function Base.isapprox(::EfficientOperations, ::Val{fields}, o1, o2; atol=atol, rtol=rtol) where fields
    (fieldcount(o1)≠fieldcount(o2) || fieldnames(o1)≠fieldnames(o2)) && return false
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

end # module
