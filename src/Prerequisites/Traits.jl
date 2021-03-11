module Traits

using ..Prerequisites: atol, rtol

export parametercount, parametername, parameterorder, parametertype, parameterpair, isparameterbound, hasparameter
export parameternames, parametertypes, parameterpairs, isparameterbounds, reparameter, promoteparameters, rawtype, fulltype
export contentcount, contentname, hascontent, fieldnameofcontent, getcontent, contentnames, disolve
export efficientoperations
export MemoryOrder, FOrder, COrder, forder, corder, indtosub, subtoind

"""
    DataType(T::DataType) -> DataType
    DataType(T::UnionAll) -> DataType

Get the DataType.
"""
Base.DataType(T::DataType) = T
Base.DataType(T::UnionAll) = DataType(T.body)

"""
    supertype(T, termination::Symbol) -> DataType

Get the supertype of `T` till termination.
"""
Base.supertype(::Type{T}, termination::Symbol) where T = _supertype(T, Val(termination))
@generated function _supertype(::Type{T}, ::Val{termination}) where {T, termination}
    result = T
    while nameof(result) != termination && nameof(result) != :Any
        result = supertype(result)
    end
    (nameof(result) != termination) && error("_supertype error: termination is not the name of a valid supertype of $(nameof(T)).")
    return result
end

"""
    parametercount(::Type{T}) where T -> Int

For a type `T`, get the number of its type parameters.
"""
@generated parametercount(::Type{T}) where T = length(DataType(T).parameters)

"""
    parametername(::Type{T}, i::Integer) where T -> Symbol

For a type `T`, get the name of its ith type parameter.
"""
parametername(::Type{T}, i::Integer) where T = _parametername(parameternames(T)|>Val, Val(i))
@generated _parametername(::Val{names}, ::Val{i}) where {names, i} = QuoteNode(names[i])

"""
    parameterorder(::Type{T}, name::Symbol) where T -> Int

For a type `T`, get the order of one of its type parameters.
"""
parameterorder(::Type{T}, name::Symbol) where T = _parameterorder(Val(name), parameternames(T)|>Val)
@generated _parameterorder(::Val{name}, ::Val{names}) where {name, names} = findfirst(isequal(name), names)::Int

"""
    parametertype(::Type{T}, name::Symbol) where T
    parametertype(::Type{T}, i::Integer) where T

For a type `T`, get the type of one of its type parameters.
"""
parametertype(::Type{T}, name::Symbol) where T = parametertype(T, parameterorder(T, name))
parametertype(::Type{T}, i::Integer) where T = _parametertype(T, Val(i))
@generated function _parametertype(::Type{T}, ::Val{i}) where {T, i}
    result = DataType(T).parameters[i]
    return isa(result, TypeVar) ? result.ub : result
end

"""
    parameterpair(::Type{T}, name::Symbol) where T
    parameterpair(::Type{T}, i::Integer) where T

For type `T`, get the name-type pair of one of its type parameters.

The result is stored in the type parameters of a `Pair`.
"""
parameterpair(::Type{T}, name::Symbol) where T = Pair{name, parametertype(T, name)}
parameterpair(::Type{T}, i::Integer) where T = Pair{parametername(T, i), parametertype(T, i)}

"""
    isparameterbound(::Type{T}, name::Symbol, D) where T -> Bool
    isparameterbound(::Type{T}, i::Integer, D) where T -> Bool
    isparameterbound(::Type{T}, ::Val{name}, D) where {T, name} -> Bool

For a type `T`, judge whether a type `D` should be considered as the upper bound of one of its type parameters.
"""
isparameterbound(::Type{}, ::Integer, ::Any) = false
isparameterbound(::Type{T}, name::Symbol, D) where T = isparameterbound(T, Val(name), D)
isparameterbound(::Type{}, ::Val{}, ::Any) = false

"""
    hasparameter(::Type{T}, name::Symbol) where T -> Bool

For type `T`, judge whether it has a type parameter specified by `name`.
"""
hasparameter(::Type{T}, name::Symbol) where T = _hasparameter(Val(name), parameternames(T)|>Val)
_hasparameter(::Val{name}, ::Val{names}) where {name, names} = name∈names

"""
    parameternames(::Type{T}) where T -> Tuple{Vararg{Symbol}}

For a type `T`, get the names of all its type parameters.
"""
parameternames(::Type{T}) where T = error("parameternames error: not defined for $(nameof(T)).")

"""
    parametertypes(::Type{T}) where T

For a type `T`, get the types of all its type parameters.

The returned types are stored in the type parameters of a `Tuple`.
"""
parametertypes(::Type{T}) where T = _parametertypes(T, parametercount(T)|>Val)
@generated function _parametertypes(::Type{T}, ::Val{C}) where {T, C}
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
parameterpairs(::Type{T}) where T = NamedTuple{parameternames(T), parametertypes(T)}

"""
    isparameterbounds(::Type{T}, ::Type{PS}) where {T, PS<:Tuple}
    isparameterbounds(::Type{T}, ::Type{PS}) where {T, PS<:NamedTuple}

For a type `T`, judge whether the types specified by `PS` should be considered as the upper bounds of its corresponding type parameters.
"""
isparameterbounds(::Type{T}, ::Type{PS}) where {T, PS<:Tuple} = _isparameterbounds(T, PS, parametercount(T)|>Val)
isparameterbounds(::Type{T}, ::Type{PS}) where {T, PS<:NamedTuple} = _isparameterbounds(T, PS, parameternames(T)|>Val)
@generated function _isparameterbounds(::Type{T}, ::Type{PS}, ::Val{C}) where {T, PS<:Tuple, C}
    exprs = []
    for i = 1:C
        P = fieldtype(PS, i)
        push!(exprs, :(isparameterbound(T, $i, $P)))
    end
    return Expr(:tuple, exprs...)
end
@generated function _isparameterbounds(::Type{T}, ::Type{PS}, ::Val{names}) where {T, PS<:NamedTuple, names}
    exprs = []
    for i = 1:length(names)
        N = QuoteNode(names[i])
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
reparameter(::Type{T}, i::Integer, P, ub::Bool=isparameterbound(T, i, P)) where T = _replace(T, Val(i), Tuple{P}, Val(ub))
reparameter(::Type{T}, name::Symbol, P, ub::Bool=isparameterbound(T, name, P)) where T = _replace(T, parameterorder(T, name)|>Val, Tuple{P}, Val(ub))
@generated function _replace(::Type{T}, ::Val{i}, ::Type{P}, ::Val{ub}) where {T, i, P, ub}
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
@generated function promoteparameters(::Type{T1}, ::Type{T2}) where {T1<:NamedTuple, T2<:NamedTuple}
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

_rawtype(T::DataType) = T.name.wrapper
_rawtype(T::UnionAll) = _rawtype(T.body)
"""
    rawtype(::Type{T}) where T -> DataType/UnionAll

Get the "raw part" of a type. That is, the type without all its type parameters.
"""
@generated rawtype(::Type{T}) where T = _rawtype(T)

"""
    fulltype(::Type{T}, ::Type{PS}, ubs::Tuple{Vararg{Bool}}=isparameterbounds(T, PS)) where {T, PS<:Tuple}
    fulltype(::Type{T}, ::Type{PS}, ubs::Tuple{Vararg{Bool}}=isparameterbounds(T, PS)) where {T, PS<:NamedTuple}

Get the full type of type `T` with the type parameters replaced by those of `PS`.

Here, `ubs` determines whether the new type parameter should be considered as the upper bound accordingly.
"""
function fulltype(::Type{T}, ::Type{PS}, ubs::Tuple{Vararg{Bool}}=isparameterbounds(T, PS)) where {T, PS<:Tuple}
    @assert parametercount(T) == fieldcount(PS) == length(ubs) "fulltype error: length-dismatched input parameters."
    return _fulltype(T, PS, Val(ubs))
end
function fulltype(::Type{T}, ::Type{PS}, ubs::Tuple{Vararg{Bool}}=isparameterbounds(T, PS)) where {T, PS<:NamedTuple}
    _fulltype(T, PS, parameternames(T)|>Val, Val(ubs))
end
@generated function _fulltype(::Type{T}, ::Type{PS}, ::Val{ubs}) where {T, PS, ubs}
    VS = [ubs[i] ? TypeVar(gensym(), fieldtype(PS, i)) : fieldtype(PS, i) for i=1:length(ubs)]
    V = Core.apply_type(T, VS...)
    for i = 1:length(VS)
        ubs[i] && (V = UnionAll(VS[i], V))
    end
    return V
end
@generated function _fulltype(::Type{T}, ::Type{PS}, ::Val{names}, ::Val{ubs}) where {T, PS<:NamedTuple, names, ubs}
    VS = [(ubs[i] ? TypeVar(gensym(), fieldtype(PS, names[i])) : fieldtype(PS, names[i])) for i=1:length(names)]
    V = Core.apply_type(T, VS...)
    for i = 1:length(VS)
        ubs[i] && (V = UnionAll(VS[i], V))
    end
    return V
end

"""
    contentcount(::Type{T}) where T -> Int

For a type `T`, get the number of its predefined contents.
"""
contentcount(::Type{T}) where T = length(contentnames(T))

"""
    contentname(::Type{T}, i::Integer) where T -> Symbol

For a type `T`, get the name of its ith predefined content.
"""
contentname(::Type{T}, i::Integer) where T = contentnames(T)[i]

"""
    hascontent(::Type{T}, name::Symbol) where T -> Bool

For a type `T`, judge whether it has a predefined content specified by `name`.
"""
hascontent(::Type{T}, name::Symbol) where T = _hascontent(name|>Val, contentnames(T)|>Val)
@generated _hascontent(::Val{name}, ::Val{names}) where {name, names} = name∈names

"""
    fieldnameofcontent(::Type{T}, i::Integer) where T -> Symbol
    fieldnameofcontent(::Type{T}, name::Symbol) where T -> Symbol
    fieldnameofcontent(::Type{T}, ::Val{name}) where {T, name} -> Symbol

For a type `T`, get the field name of one of its predefined contents.
"""
fieldnameofcontent(::Type{T}, i::Integer) where T = fieldnameofcontent(T, contentname(T, i))
fieldnameofcontent(::Type{T}, name::Symbol) where T = fieldnameofcontent(T, Val(name))
fieldnameofcontent(::Type{}, ::Val{name}) where name = name

"""
    getcontent(m, i::Integer)
    getcontent(m, name::Symbol)

Get the value of the predefined content of `m`. 
"""
getcontent(m, i::Integer) = getcontent(m, contentname(typeof(m), i))
getcontent(m, name::Symbol) = getfield(m, fieldnameofcontent(typeof(m), name))

"""
    contentnames(::Type{T}) where T -> Tuple{Vararg{Symbol}}

For a type `T`, define the names of its predefined contents.
"""
contentnames(::Type{}) = ()

"""
    disolve(m, f::Function=identity, args::Tuple=(), kwargs::NamedTuple=NamedTuple()) -> Tuple

Convert `m` to a tuple by the function `f` applied elementally to its fields with the extra positional arguments (`args`) and keyword arguments (`kwargs`). 

The underlying called interface is the `map` function when `f` is applied to each field of `m`:
```julia
disolve(m, Val(name), f, args, kwargs)
```
Here, `name` is the name of a field of `m`:
1. When the field corresponds to a predefined one, `name` should be the predefined name but not the actual one;
2. When the field does not correspond to a predefined one, `name` is the actual one.

Basically, the rule of how `f` operates on each field of `m` can be overriden by redefining the above `disolve` function.
!!!note
   The default `disolve` function ignores the operation of function `f` and just return the field value of `m`.
"""
function disolve(m, f::Function=identity, args::Tuple=(), kwargs::NamedTuple=NamedTuple())
    tuplehelper(m, f, args, kwargs, fieldnamelookups(typeof(m), contentnames(typeof(m))|>Val))
end
@generated function fieldnamelookups(::Type{T}, ::Val{names}) where {T, names}
    exprs = []
    for name in QuoteNode.(names)
        push!(exprs, :(fieldnameofcontent(T, $name)))
    end
    return Expr(:curly, :NamedTuple, Expr(:tuple, exprs...), Expr(:curly, :Tuple, QuoteNode.(names)...))
end
@generated function tuplehelper(m, f::Function, args::Tuple, kwargs::NamedTuple, ::Type{LP}) where LP<:NamedTuple
    exprs = []
    for name in fieldnames(m)
        name = Val(hasfield(LP, name) ? fieldtype(LP, name) : name)
        push!(exprs, :(disolve(m, $name, f, args, kwargs)))
    end
    return Expr(:tuple, exprs...)
end

"""
    disolve(m, ::Val{name}, f::Function, args::Tuple, kwargs::NamedTuple) where name

Disolve the field specified by `name` of `m` by the function `f` applied with the extra positional arguments (`args`) and keyword arguments (`kwargs`).
"""
function disolve(m, ::Val{name}, f::Function, args::Tuple, kwargs::NamedTuple) where name
    _disolve(m, Val(name), f, args, kwargs, Val(hascontent(typeof(m), name)))
end
_disolve(m, ::Val{name}, f::Function, args::Tuple, kwargs::NamedTuple, ::Val{true}) where name = f(getcontent(m, name), args...; kwargs...)
_disolve(m, ::Val{name}, ::Function, ::Tuple, ::NamedTuple, ::Val{false}) where name = getfield(m, name)

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
@generated function Base.:(==)(::EfficientOperations, o1, o2)
    fcount = fieldcount(o1)
    if fcount == fieldcount(o2)
        if fcount == 0
            return :(true)
        else
            expr = :(getfield(o1, 1) == getfield(o2, 1))
            for i = 2:fcount
                expr = Expr(:&&, expr, :(getfield(o1, $i) == getfield(o2, $i)))
            end
            return expr
        end
    else
        return :(false)
    end
end

"""
    isequal(::EfficientOperations, o1, o2) -> Bool

Compare two objects and judge whether they are eqaul to each other.
"""
@generated function Base.isequal(::EfficientOperations, o1, o2)
    fcount = fieldcount(o1)
    if fcount == fieldcount(o2)
        if fcount == 0
            return :(true)
        else
            expr = :(isequal(getfield(o1, 1), getfield(o2, 1)))
            for i = 2:fcount
                expr = Expr(:&&, expr, :(isequal(getfield(o1, $i), getfield(o2, $i))))
            end
            return expr
        end
    else
        return :(false)
    end
end

"""
    <(::EfficientOperations, o1, o2) -> Bool

Compare two objects and judge whether the first is less than the second.
"""
@generated function Base.:<(::EfficientOperations, o1, o2)
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
@generated function Base.isless(::EfficientOperations, o1, o2)
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
    isapprox(::EfficientOperations, ::Val{Names}, o1, o2; atol::Real=atol, rtol::Real=rtol) where Names -> Bool

Compare two objects and judge whether they are inexactly equivalent to each other.
"""
@generated function Base.isapprox(::EfficientOperations, ::Val{Names}, o1, o2; atol::Real=atol, rtol::Real=rtol) where Names
    fcount = fieldcount(o1)
    if fcount == fieldcount(o2)
        if fcount == 0
            return :(true)
        else
            exprs = [:(isapprox(getfield(o1, $name), getfield(o2, $name); atol=atol, rtol=rtol)) for name in QuoteNode.(Names)]
            for name in fieldnames(o1)
                (name ∉ Names) && (name = QuoteNode(name); push!(exprs, :(isequal(getfield(o1, $name), getfield(o2, $name)))))
            end
            return Expr(:&&, exprs...)
        end
    else
        return :(false)
    end
end

"""
    replace(::EfficientOperations, o; kwargs...) -> typeof(o)

Return a copy of the input object with some of the field values replaced by the keyword arguments.
"""
@generated function Base.replace(::EfficientOperations, o; kwargs...)
    exprs = [:(get(kwargs, $name, getfield(o, $name))) for name in QuoteNode.(fieldnames(o))]
    return :(rawtype(typeof(o))($(exprs...)))
end

abstract type MemoryOrder end
struct FOrder <: MemoryOrder end
struct COrder <: MemoryOrder end
"""
    forder

Indicate that the convertion between Cartesian index and linear index is using the Fortran order.
"""
const forder = FOrder()
"""
    corder

Indicate that the convertion between Cartesian index and linear index is using the C/C++ order.
"""
const corder = COrder()

@generated tail(ts::NTuple{N}) where N = Expr(:tuple, [:(ts[$i]) for i = 2:N]...)
@generated head(ts::NTuple{N}) where N = Expr(:tuple, [:(ts[$i]) for i = 1:N-1]...)

"""
    indtosub(dims::Tuple, ind::Int, order::FOrder) -> Tuple
    indtosub(dims::Tuple, ind::Int, order::COrder) -> Tuple

Convert an linear index to Cartesian index. Fortran-order or C-order can be assigned.
"""
function indtosub(dims::Tuple, ind::Int, order::FOrder)
    (length(dims) == 0) && return ()
    return ((ind-1)%dims[1]+1, indtosub(tail(dims), (ind-1)÷dims[1]+1, order)...)
end
function indtosub(dims::Tuple, ind::Int, order::COrder)
    (length(dims) == 0) && return ()
    return (indtosub(head(dims), (ind-1)÷dims[end]+1, order)..., (ind-1)%dims[end]+1)
end

"""
    subtoind(dims::NTuple{N, Int}, inds::NTuple{N, Int}, order::FOrder) where N -> Int
    subtoind(dims::NTuple{N, Int}, inds::NTuple{N, Int}, order::COrder) where N -> Int

Convert an Cartesian index to linear index. Fortran-order or C-order can be assigned.
"""
function subtoind(dims::NTuple{N, Int}, inds::NTuple{N, Int}, order::FOrder) where N
    length(dims) == 0 && return 1
    return (subtoind(tail(dims), tail(inds), order)-1)*dims[1]+inds[1]
end
function subtoind(dims::NTuple{N, Int}, inds::NTuple{N, Int}, order::COrder) where N
    length(dims) == 0 && return 1
    return (subtoind(head(dims), head(inds), order)-1)*dims[end]+inds[end]
end

end # module
