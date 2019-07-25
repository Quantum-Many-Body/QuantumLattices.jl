```@meta
CurrentModule=QuantumLattices.Prerequisites.Factories
```

```@setup factories
push!(LOAD_PATH,"../../../../src/")
using QuantumLattices.Prerequisites.Factories
```

# Factories

The aim of `Factories` is to provide tools to hack into Julia codes without knowing the details of their abstract syntax trees and regularize the mechanism to "escape" variables in `Expr` expressions, so that users can manipulate the existing codes, modify them and generate new ones in macros. In particular, a factory in this module means the representation of certain blocks of Julia codes by a usual Julia struct. This representation is much easier to comprehend than the canonical `Expr` representation and makes it far more convenient to define macros. In general, we propose the following requirements that any factory must satisfy:
* *DECOMPOSITION* - An `Expr` expression can be decomposed into its corresponding factory by the factory's constructor.
* *COMPOSITION* - A factory can compose its corresponding `Expr` expression by calling itself.
* *ESCAPE* - A variable should be or not be escaped in the composed `Expr` expression by a factory depends on predefined escape mechanisms.
These three requirements also define the basic interfaces to interact with factories. In practice, we combine the second and third in a single interface, i.e. by passing an instance of certain concrete [`EscapeMechanism`](@ref) as the only argument of calling a factory, the needed `Expr` expression with variables correctly escaped can be obtained.

## Escape mechanisms

We adopt Julia structs to denote escape mechanisms so that we can utilize Julia's multi-dispatch to implement different mechanisms whereas keeping the same interface.

### EscapeMechanism

[`EscapeMechanism`](@ref) is the abstract type for all concrete escape mechanisms.

### Escaped

[`Escaped`](@ref) has only one attribute:
* `names::NTuple{N,Symbol} where N`: the names of variables to be escaped
Apparently, a variable should be escaped if its name is in the `names` of an `Escaped`.
This mechanism suits a factory whose variables should be unescaped by default.

### UnEscaped

[`UnEscaped`](@ref) also has only on attribute:
* `names::NTuple{N,Symbol} where N`: the names of variables not to be escaped
Obviously, on the contrary to [`Escaped`](@ref), a variable should be escaped if its name is not in the `names` of an `UnEscaped`.
This mechanism suits a factory whose variables should be escaped by default.

### MixEscaped

[`MixEscaped`](@ref) has two attributes:
* `escaped::Escaped`: the escaped part of the mixed mechanism
* `unescaped::UnEscaped`: the UnEscaped part of the mixed mechanism
This mechanism suits complex factories that parts of it suit the "escaped" mechanism while others suit the "unescaped" mechanism.

### RawExpr

[`RawExpr`](@ref) has no attributes and it means "raw expression without any variable escaped". This mechanism is used for the print of all factories by default.

## Concrete factories

Out of practical purposes, we implement 7 kinds of factories, i.e. **[`Inference`](@ref)**, **[`Argument`](@ref)**, **[`Parameter`](@ref)**, **[`Field`](@ref)**, **[`Block`](@ref)**, **[`FunctionFactory`](@ref)** and **[`TypeFactory`](@ref)**, which represent **a type inference**, **a function argument**, **a method or type parameter**, **a struct field**, **a `begin ... end` block**, **a function itself** and **a struct itself**, respectively. Some of the basic methods making the above three requirements fulfilled with these types are based on the powerful functions defined in [`MacroTools`](https://github.com/MikeInnes/MacroTools.jl).

We want to give a remark that although the types and functions provided in this module helps a lot for the definition of macros, macros should not be abused. On the one hand, some macros may change the language specifications, which makes it hard to understand the codes, and even splits the community; on the one hand, macros usually increases the precompiling/jit time, which means enormous uses of macros in a module may lead to an extremely long load time. Besides, due to the limited ability of the author, the codes in this module are not optimal, which adds to the jit overhead. Any promotion that keeps the interfaces unchanged is welcomed.

### Inference

An [`Inference`](@ref) has 3 attributes:
* `head::Union{Symbol, Nothing}`: the head of the type inference, which must be one of `(nothing,:(<:),:curly)`
* `name::Union{Symbol, Nothing}`: the name of the type inference
* `params::Union{Inference, Vector{Inference},Nothing}`: the parameters of the type inference

All valid expressions representing type inferences can be passed to the constructor.

A type variable represented by a Symbol:
```@example factories
Inference(:T)
```

A range of types represented by a UnionAll:
```@example factories
Inference(:(<:Number))
```

A type with formal type parameters:
```@example factories
Inference(:(Vector{T}))
```

A definite parametric type:
```@example factories
Inference(:(Vector{Tuple{String,Int}}))
```

A type with type parameters being UnionAlls:
```@example factories
Inference(:(Type{<:Number}))
```

On the other hand, you can use the macro [`@inference`](@ref) to construct an `Inference` directly from a type inference:
```@example factories
@inference Vector{Tuple{String,Int}}
```
!!! note
    1. `Inference` is a recursive struct, i.e. it recursively decomposes a type inference until the final type inference is just a `Symbol`.
    2. When the input expression is a `Symbol`, the `head` and `params` attributes of the resulting `Inference` is `nothing`. Otherwise, its `head` is the same with that of the input expression, and the `args` of the input expression will be further decomposed, whose result will be stored in `params`.
    3. When the head of the input expression is `:(<:)`, the `params` is an `Inference` whereas when the head of the input expression is `:curly`, the `params` is a `Vector{Inference}`.

[`Inference`](@ref) uses the [`UnEscaped`](@ref) mechanism to escape variables, i.e. you have to specify the names that are not to be escaped.

Empty `UnEscaped()` escapes everything:
```@example factories
Inference(:(Vector{T}))(UnEscaped())
```

Unescape formal type parameters:
```@example factories
Inference(:(Vector{T}))(UnEscaped(:T))
```

Unescape both type names and type parameters:
```@example factories
Inference(:(Vector{T}))(UnEscaped(:Vector,:T))
```

### Argument

An [`Argument`](@ref) has 4 attributes:
* `name::Union{Symbol,Nothing}`: the name of the argument
* `type::Inference`: the type inference of the argument
* `slurp::Bool`: whether the argument should be expanded by `...`
* `default::Any`: the default value of the argument, `nothing` for those with no default values

All valid expressions representing the arguments of functions can be passed to the constructor.

A formal argument represented by a Symbol:
```@example factories
Argument(:arg)
```

A formal argument with its type specified:
```@example factories
Argument(:(arg::ArgType))
```

A slurped formal argument:
```@example factories
Argument(:(arg::ArgType...))
```

A formal argument with a default value
```@example factories
Argument(:(arg::ArgType=default))
```

Or you can use the macro [`@argument`](@ref) for a direct construction from an argument declaration:
```@example factories
@argument arg::ArgType=default
```
The construction from such expressions is based on the `MacroTools.splitarg` function.

[`Argument`](@ref) uses the [`MixEscaped`](@ref) mechanism to escape variables, with the [`UnEscaped`](@ref) mechanism for `type` and [`Escaped`](@ref) mechanism for `default`, e.g.
```@example factories
Argument(:(arg::Real=zero(Int)))(MixEscaped(UnEscaped(),Escaped(:zero,:Int)))
```
It can be seen that the name of an argument will never be escaped, which is obvious since the name of a function argument is always local. By the way, the composition of an [`Argument`](@ref) expression is based on the `MacroTools.combinearg` function.

### Parameter

A [`Parameter`](@ref) has 3 attributes:
* `head::Union{Symbol,Nothing}`: the head of the parameter, must be either `:(<:)` or `:(::)` if not `nothing`
* `name::Union{Symbol,Nothing}`: the name of the parameter
* `type::Union{Inference,Nothing}`: the type inference of the parameter

All expressions that represent type parameters or method parameters are allowed to be passed to the constructor.

A formal type/function parameter represented by a Symbol:
```@example factories
Parameter(:T)
```

A formal type/function parameter represented by a UnionAll with the name omitted:
```@example factories
Parameter(:(<:Number))
```

A formal type/function parameter represented by a UnionAll with the name specified:
```@example factories
Parameter(:(T<:Number))
```

The macro [`@parameter`](@ref) completes the construction directly from a parameter declaration:
```@example factories
@parameter T<:Number
```
!!! note
    1. We use `nothing` to denote a missing `name` or `type.`
    2. Two subtle situations of type/method parameters, e.g. `MyType{T}` and `MyType{Int}`, should be distinguished by `Parameter(:T)` and `Parameter(:(::Int))`.

[`Parameter`](@ref) uses the [`UnEscaped`](@ref) mechanism to escape variables, too, e.g.
```@example factories
Parameter(:(N<:Vector{T}))(UnEscaped(:T))
```
As is similar to [`Argument`](@ref), the `name` of a method/type parameter will never be escaped because of its local scope.

### Field

A [`Field`](@ref) has 2 attributes:
* `name::Symbol`: the name of the field
* `type::Inference`: the type inference of the field

Legal expressions can be used to construct a `Field` instance by its constructor.

A field name specified by a Symbol:
```@example factories
Field(:field)
```

A field with both the name and type specified:
```@example factories
Field(:(field::FieldType))
```

A field with the type being parametric:
```@example factories
Field(:(field::ParametricType{T}))
```

The macro [`@field`](@ref) is also provided to help the construction directly from a field declaration:
```@example factories
@field field::FieldType
```
The construction from these expressions is based on the `MacroTools.splitarg` function.

[`Field`](@ref) uses the [`UnEscaped`](@ref) mechanism to escape variables as well, e.g.
```@example factories
Field(:(field::Dict{N,D}))(UnEscaped(:N,:D))
```
The name of a struct will never be escaped either because it is a local variable tightly binding to a struct. It is noted that the composition of field expressions is based on the `MacroTools.combinefield` function.

### Block

A [`Block`](@ref) has only one attribute:
* `body::Vector{Any}`: the body of the `begin ... end` block

Any expression can be passed to the constructor of `Block`.

A single expression:
```@example factories
Block(:(x=1))
```

A semicolon expression:
```@example factories
Block(:(x=1;y=2))
```

A `begin ... end` block:
```@example factories
Block(:(begin x=1 end))
```

A `quote` block:
```@example factories
Block(quote
        x=1
        y=2
    end)
```

Or you can construct a `Block` instance directly from any code by the macro [`@block`](@ref):
```@example factories
@block x=1 y=2
```
The body of a `block` can also be extended by the [`push!`](@ref) function or the [`@push!`](@ref) macro.
!!! note
    1. The body of a `Block` is somewhat "flattened", i.e. it contains no `begin ... end` blocks. During the initialization, any such input block will be unblocked and added to the body part by part. So is the [`push!`](@ref) and [`@push!`](@ref) procedures.
    2. All `LineNumberNode`s generated by the input codes will also be included in the block's body. However, you can use [`rmlines!`](@ref) or [`@rmlines!`](@ref) to remove them from the body of an existing `Block`, or use [`rmlines`](@ref) or [`@rmlines`](@ref) to get a copy with them removed in the body.

[`Block`](@ref) uses the [`Escaped`](@ref) mechanism to escape variables. This is because variables in a block are often local ones and should not be escaped. Therefore, only those defined in other modules should be noted and escaped, which usually constitute the minority. For example,
```@example factories
Block(:(x=1;y=2;z=Int[1,2,3]))(Escaped(:Int))
```

### FunctionFactory

A [`FunctionFactory`](@ref) has 7 attributes:
* `name::Union{Symbol,Expr}`: the name of the function
* `params::Vector{Inference}`: the method parameters of the function
* `args::Vector{Argument}`: the positional arguments of the function
* `kwargs::Vector{Argument}`: the keyword arguments of the function
* `rtype::Inference`: the return type of the function
* `whereparams::Vector{Parameter}`: the method parameters specified by the `where` keyword
* `body::Block`: the body of the function

All expressions that represent functions are allowed to be passed to the constructor.

A function with no arguments:
```@example factories
FunctionFactory(:(f()=nothing))
```

A function with arguments:
```@example factories
FunctionFactory(:(f(x)=x))
```

A function with arguments that have default values:
```@example factories
FunctionFactory(:(f(x::Int, y::Int; choice::Function=sum)=choice(x, y)))
```

A function with arguments whose type are specified by `where` keyword:
```@example factories
FunctionFactory(:(f(x::T, y::T; choice::Function=sum) where T<:Number=choice(x, y)))
```

A function with return type:
```@example factories
FunctionFactory(:((f(x::T, y::T; choice::Function=sum)::T) where T<:Number=choice(x, y)))
```

A multiline function:
```@example factories
FunctionFactory(:(
    function (f(x::T, y::T; choice::Function=sum)::T) where T<:Number
        choice(x, y)
    end
))
```

A function in a `quote` block:
```@example factories
FunctionFactory(
    quote
        function (f(x::T, y::T; choice::Function=sum)::T) where T<:Number
            choice(x, y)
        end
    end
)
```

Similarly, an instance can also be constructed from the macro [`@functionfactory`](@ref):
```@example factories
@functionfactory (f(x::T, y::T; choice::Function=sum)::T) where T<:Number=choice(x, y)
```
The construction from such expressions are based on the `MacroTools.splitdef` function.
!!! note
    1. Since Julia 0.7, the form `MyType{D}(data::D) where D` only appears in struct constructors, therefore, the attribute `:params` of a function factory is nonempty only when this factory aims to represent a struct constructor.
    2. Usually, the name of a function factory is a `Symbol`. However, if the factory aims to extend some methods of a function defined in another module, e.g., `Base.eltype`, the name will be an `Expr`.

[`FunctionFactory`](@ref) adopts the [`MixEscaped`](@ref) mechanism to escape variables, with [`UnEscaped`](@ref) for `params`, `args`, `kwargs`, `rtype` and `whereparams` while [`Escaped`](@ref) for `name` and `body`. It is worth to emphasize that the name of a function factory belongs to the `Escaped` part. Therefore, when it is an `Expr`, it will never be escaped because an `Expr` cannot be a element of a `NTuple{N,Symbol} where N`. See following examples.

Escape the function name:
```@example factories
FunctionFactory(:(
    (f(x::T, y::T; choice::Function=sum)::T) where T<:Number=max(x, y, choice(x, y))
    ))(MixEscaped(UnEscaped(:T),Escaped(:f, :max,)))
```

Do not escape the function name:
```@example factories
FunctionFactory(:(
    (f(x::T, y::T; choice::Function=sum)::T) where T<:Number=max(x, y, choice(x, y))
    ))(MixEscaped(UnEscaped(:T),Escaped(:max)))
```

The compositions of function expressions are based on the `MacroTools.combinedef` function.

Other features include:
* Positional arguments can be added by [`addargs!`](@ref) or [`@addargs!`](@ref)
* Keyword arguments can be added by [`addkwargs!`](@ref) or [`@addkwargs!`](@ref)
* Where parameters can be added by [`addwhereparams!`](@ref) or [`@addwhereparams!`](@ref)
* Body can be extended by [`extendbody!`](@ref) or [`@extendbody!`](@ref)

### TypeFactory

A [`TypeFactory`](@ref) has 6 attributes:
* `name::Symbol`: the name of the struct
* `mutable::Bool`: whether or not the struct is mutable
* `params::Vector{Parameter}`: the type parameters of the struct
* `supertype::Inference`: the supertype of the struct
* `fields::Vector{Field}`: the fields of the struct
* `constructors::Vector{FunctionFactory}`: the inner constructors of the struct

Any expression representing valid struct definitions can be passed to the constructor.

A struct with no fields:
```@example factories
TypeFactory(:(struct StructName end))
```

A parametric struct with no fields:
```@example factories
TypeFactory(:(struct StructName{T} end))
```

An empty struct with a supertype:
```@example factories
TypeFactory(:(struct Child{T} <: Parent{T} end))
```

A fielded struct with a supertype
```@example factories
TypeFactory(:(
    struct Child{T<:Number} <: Parent{T}
        field1::T
        field2::T
    end
))
```

A struct in a `quote` block:
```@example factories
TypeFactory(
    quote
        struct Child{T<:Number} <: Parent{T}
            field1::T
            field2::T
        end
    end
)
```

Also, the macro [`@typefactory`](@ref) supports the construction directly from a type definition:
```@example factories
@typefactory struct Child{T<:Number} <: Parent{T}
                field1::T
                field2::T
                Child(field1::T,field2::T=zero(T)) where T=new{T}(field1,field2)
            end
```
The construction from these expressions is based on the `MacroTools.splitstructdef` function.

[`TypeFactory`](@ref) also uses the [`MixEscaped`](@ref) mechanism to escape variables, with the [`UnEscaped`](@ref) part for `params`, `supertype` and `fields`, the [`Escaped`](@ref) part for `name`, and both for `constructors`. For example,
```@example factories
@typefactory(struct Child{T<:Number} <: Parent{T}
    field::T
    Child(field::T) where T=new{T}(field)
end)(MixEscaped(UnEscaped(:T),Escaped(:Child)))
```
The composition of a type expression is based on the `MacroTools.combinestructdef` function.

Other features include:
* Fields can be added by [`addfields!`](@ref) or [`@addfields!`](@ref)
* Type parameters can be added by [`addparams!`](@ref) or [`@addparams!`](@ref)
* Inner constructors can be added by [`addconstructors!`](@ref) or [`@addconstructors!`](@ref)

## Manual

```@autodocs
Modules = [Factories]
Order =  [:module, :constant, :type, :macro, :function]
```
