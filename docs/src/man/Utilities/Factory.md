```@meta
CurrentModule=Hamiltonian.Utilities.Factory
```

```@setup factory
push!(LOAD_PATH,"../../../../src/")
using Hamiltonian.Utilities.Factory
```

# Factory

The aim of `Factory` is to provide tools to hack into Julia codes without knowing the details of their abstract syntax trees and regularize the mechanism to "escape" variables in `Expr` expressions, so that users can manipulate the existing codes, modify them and generate new ones in macros. In particular, `Factory` in this module means the representation of certain blocks of Julia codes by a usual Julia struct. This representation is much easier to comprehend than the canonical `Expr` representation and makes it far more convenient to define macros. In general, we propose two basic requirements that any factory must satisfy:
* **FEATURE-1**: A concrete factory can be constructed from a legal `Expr` expression that represents the block of codes it aims to represent;
* **FEATURE-2**: The canonical `Expr` expression of a block of codes that a concrete factory represents can be obtained by "calling" the factory itself.
Also, we propose a simple mechanism to implement the escape of variables in factories:
* **FEATURE-3**: During the construction from `Expr` expressions, a keyword argument `unescaped`, which is a tuple of `Symbol`s, can be passed to the constructor, and all variables whose names are not in this keyword argument will be escaped.
Different from the former two requirements, this mechanism is just an optional feature of concrete factories. Meanwhile, all of these features also defines the basic interfaces to interact with factories.

Out of practical purposes, we implemente 7 kinds of factories, i.e. *[`Inference`](@ref)*, *[`Argument`](@ref)*, *[`Parameter`](@ref)*, *[`Field`](@ref)*, *[`Block`](@ref)*, *[`FunctionFactory`](@ref)* and *[`TypeFactory`](@ref)*, which represent *a type inference*, *a function argument*, *a method or type parameter*, *a struct field*, *a `begin ... end` block*, *a function itself* and *a struct itself*, respectively. Some of the basic methods making the above requirements fulfilled with these types are based on the powerful functions defined in [`MacroTools`](https://github.com/MikeInnes/MacroTools.jl).

## Inference

An [`Inference`](@ref) has 4 attributes:
* `head::Union{Symbol,Nothing}`: the head of the type inference, which must be one of `(nothing,:(::),:(<:),:curly)`
* `name::Union{Symbol,Nothing}`: the name of the type inference
* `params::Union{Inference,Vector{Inference},Nothing}`: the parameters of the type inference
* `escape::Bool`: whether the name of the type inference should be escaped when calling the factory
All valid expressions representing type inferences can be passed to the constructor:
```@repl factory
Inference(:T,unescaped=(:T,))
Inference(:(::T),unescaped=(:T,))
Inference(:(<:Number))
Inference(:(Vector{T}),unescaped=(:T,))
Inference(:(Vector{Tuple{String,Int}}))
Inference(:(Type{<:Number}))
```
On the other hand, you can use the macro [`@inference`](@ref) to construct an `Inference` directly from a type inference:
```@repl factory
@inference Vector{Tuple{String,Int}}
```
!!! note
    1. `Inference` is a recursive struct, i.e. it recursively decomposes a type inference until the final type inference is just a `Symbol`.
    2. When the input expression representing a type inference is a `Symbol`, the `head` and `params` attributes of the resulting `Inference` is `nothing`. Otherwise, its `head` is the same with that of the input expression, and the `args` of the input expression will be further decomposed, whose result will be stored in `params`.
    3. When the head of the input expression is `:(::)` or `:(<:)`, the `params` is an `Inference` whereas when the head of the input expression is `:curly`, the `params` is a `Vector{Inference}`.

## Argument

An [`Argument`](@ref) has 4 attributes:
* `name::Union{Symbol,Nothing}`: the name of the argument
* `type::Inference`: the type inference of the argument
* `slurp::Bool`: whether the argument should be expanded by `...`
* `default::Any`: the default value of the argument, `nothing` for those with no default values
All valid expressions representing the arguments of functions can be passed to the constructor:
```@repl factory
Argument(:arg)
Argument(:(arg::ArgType),unescaped=(:ArgType,))
Argument(:(arg::ArgType...),unescaped=(:ArgType,))
Argument(:(arg::ArgType=default),unescaped=(:ArgType,:default))
```
Or you can use the macro [`@argument`](@ref) for construction directly from an argument declaration:
```@repl factory
@argument arg::ArgType=default (:ArgType,:default)
```
The construction from such expressions is based on the the `MacroTools.splitarg` function. On the other hand, calling an instance of `Argument` will get the corresponding `Expr` expression, e.g.,
```julia-repl
julia> Argument(:(arg::ArgType=default),unescaped=(:ArgType,:default))()
:(arg::ArgType=default)
```
This feature is based on the `MacroTools.combinearg` function.

## Parameter

A [`Parameter`](@ref) has 2 attributes:
* `name::Union{Symbol,Nothing}`: the name of the parameter
* `type::Union{Inference,Nothing}`: the type inference of the parameter
All expressions that represent a type parameter or a method parameter are allowed to be passed to the constructor:
```@repl factory
Parameter(:T)
Parameter(:(<:Number))
Parameter(:(T<:Number))
Parameter(:(::Int))
```
The macro [`@parameter`](@ref) completes the construction directly from a parameter declaration:
```@repl factory
@parameter T<:Number
```
!!! note
    1. We use `nothing` to denote a missing `name` or `type.`
    2. Two subtle situations of type/method parameters, e.g. `MyType{T}` and `MyType{Int}`, are distinguished by `Parameter(:T)` and `Parameter(:(::Int))`. For the first case, `Parameter(:T).name==:T` and `Parameter(:T).type==:nothing` while for the second case, `Parameter(:(::Int)).name==:nothing` and `Parameter(:(::Int)).name==:Int`. Moreover, the callings of them return different forms, e.g. `Parameter(:T)()==:T` and `Parameter(:(::Int))==:(<:$(esc(Int)))`.

## Field

A [`Field`](@ref) has 2 attributes:
* `name::Symbol`: the name of the field
* `type::Inference`: the type inference of the field
Legal expressions can be used to construct a `Field` instance by its constructor:
```@repl factory
Field(:field)
Field(:(field::FieldType),unescaped=(:FieldType,))
Field(:(field::ParametricType{T}),unescaped=(:ParametricType,:T))
```
The macro [`@field`](@ref) is also provided to help the construction directly from a field declaration:
```@repl factory
@field field::FieldType (:FieldType,)
```
The construction from these expressions is based on the `MacroTools.splitarg` function and the convertion to these expressions is based on the `MacroTools.combinefield` function.

## Block

A [`Block`](@ref) has only one attribute:
* `body::Vector{Any}`: the body of the `begin ... end` block
Any expression can be passed to the constructor of `Block`:
```@repl factory
Block(:(x=1))
Block(:(x=1;y=2))
Block(:(begin x=1 end))
Block(quote
        x=1
        y=2
    end)
```
Or you can construct a `Block` instance directly from any code by the macro [`@block`](@ref):
```@repl factory
@block x=1 y=2
```
The body of a `block` can also be extended by the [`push!`](@ref) function or the [`@push!`](@ref) macro.
!!! note
    1. The body of a `Block` is somewhat "flattened", i.e. it contains no `begin ... end` blocks. During the initialization, any such input block will be unblocked and added to the body part by part. So is the [`push!`](@ref) and [`@push!`](@ref) processes.
    2. All `LineNumberNode`s generated by the input codes will also be included in the block's body. However, you can use [`rmlines!`](@ref) or [`@rmlines!`](@ref) to remove them from the body of an existing `Block`, or use [`rmlines`](@ref) or [`@rmlines`](@ref) to get a copy with them removed in the body.
    3. `Block` is the only factory that **DOES NOT** cope with the "escape" mechanism we propose, because usually variables in a block are local ones and should not be escaped. However, the function [`escape`](@ref) is provided to escape the variables in a `Expr` expression with a reverse mechanism, where a tuple of variables to be escaped should be transmitted.

## FunctionFactory

A [`FunctionFactory`](@ref) has 6 attributes:
* `name::Symbol`: the name of the function
* `args::Vector{Argument}`: the positional arguments of the function
* `kwargs::Vector{Argument}`: the keyword arguments of the function
* `rtype::Inference`: the return type of the function
* `params::Vector{Parameter}`: the method parameters specified by the `where` keyword
* `body::Block`: the body of the function
All expressions that represent functions are allowed to be passed to the constructor:
```@repl factory
FunctionFactory(:(f()=nothing))
FunctionFactory(:(f(x)=x))
FunctionFactory(:(f(x::Int,y::Int;choice::Function=sum)=choice(x,y)))
FunctionFactory(:(f(x::T,y::T;choice::Function=sum) where T<:Number=choice(x,y)),unescaped=(:T,))
FunctionFactory(:((f(x::T,y::T;choice::Function=sum)::T) where T<:Number=choice(x,y)),unescaped=(:T,))
FunctionFactory(:(
    function (f(x::T,y::T;choice::Function=sum)::T) where T<:Number
        choice(x,y)
    end
    ),
    unescaped=(:T,)
)
FunctionFactory(
    quote
        function (f(x::T,y::T;choice::Function=sum)::T) where T<:Number
            choice(x,y)
        end
    end,
    unescaped=(:T,)
)
```
Similarly, an instance can also be constructed from the macro [`@functionfactory`](@ref):
```@repl factory
@functionfactory (f(x::T,y::T;choice::Function=sum)::T) where T<:Number=choice(x,y) (:T,)
```
The construction from and the convertion to such expressions are based on the `MacroTools.splitdef` and `MacroTools.combinedef` functions, respectively.
!!! note
    Because the form `f{T}(x::T,y::T;choice::Function=sum)` has no longer been supported since Julia 0.7, the entry `:params` in the returned dict by `MacroTools.splitarg` is always missing. Therefore, we abandon its corresponding field in `FunctionFactory` but use the attribute `:params` to denote the `:whereparams` entry.

Other features include:
* Positional arguments can be added by [`addargs!`](@ref) or [`@addargs!`](@ref)
* Keyword arguments can be added by [`addkwargs!`](@ref) or [`@addkwargs!`](@ref)
* Method parameters can be added by [`addparams!`](@ref) or [`@addparams!`](@ref)
* Body can be extended by [`extendbody!`](@ref) or [`@extendbody!`](@ref)

## TypeFactory

A [`TypeFactory`](@ref) has 6 attributes:
* `name::Symbol`: the name of the struct
* `mutable::Bool`: whether or not the struct is mutable
* `params::Vector{Parameter}`: the type parameters of the struct
* `supertype::Inference`: the supertype of the struct
* `fields::Vector{Field}`: the fields of the struct
* `constructors::Vector{FunctionFactory}`: the inner constructors of the struct
Any expression representing valid struct definitions can be passed to the constructor:
```@repl factory
TypeFactory(:(struct StructName end))
TypeFactory(:(struct StructName{T} end),unescaped=(:T,))
TypeFactory(:(struct Child{T} <: Parent{T} end),unescaped=(:T,))
TypeFactory(:(
    struct Child{T<:Number} <: Parent{T}
        field1::T
        field2::T
    end
    ),
    unescaped=(:T,)
)
TypeFactory(
    quote
        struct Child{T<:Number} <: Parent{T}
            field1::T
            field2::T
        end
    end,
    unescaped=(:T,)
)
```
Also, the macro [`@typefactory`](@ref) supports the construction directly from a type definition:
```@repl factory
@typefactory struct Child{T<:Number} <: Parent{T} field1::T; field2::T; Child(field1::T,field2::T=zero(T)) where T=new{T}(field1,field2) end (:T,)
```
The construction from these expressions is based on the `MacroTools.splitstructdef` function. Meanwhile, the convertion to the corresponding expression from a `TypeFactory` is based on the `MacroTools.combinestructdef` function.

Other features include:
* Fields can be added by [`addfields!`](@ref) or [`@addfields!`](@ref)
* Type parameters can be added by [`addparams!`](@ref) or [`@addparams!`](@ref)
* Inner constructors can be added by [`addconstructors!`](@ref) or [`@addconstructors!`](@ref)

## More about escape

As you may have noticed, during the construction of factoies, we often provide a keyword argument `unescaped` to tell what variables should not be escaped. The criterion to determine such variables is quite direct and simple: if a variable is local to current module, or to be defined, or to be declared, it should not be escaped, while if it has been defined or declared in other modules, it should be escaped. Therefore, the name of a function argument, or a type/method parameter, or a struct field, should not be escaped because the first is a local variable and the latter two are to-be-declared ones. In fact, by design, these names will **NEVER** be escaped even when they are not provided in the `unescaped` keyword argument. Another common situation where unescaped variables exist occurs when methods or structs have parameters. In this case, the type inferences specified by these method/type parameters should not be escaped because they are also local ones. Note that this situation is not automatically handled by the above factories, because the method/type parameters can be modified after the factoies are constructed, which means, you have to offer the `unescaped` tuple by hand. To help decorate existing methods or structs, we define an inquiry function [`paramnames`](@ref) to obtain the method/type parameter names. By the way, the factory [`Block`](@ref) does not support variable escape by design, because it is usually used as the body of a function, which is always a local environment. This is in sharp contrast to the method/type declaration, where type inferences are usually involved with variables that are defined or declared in other modules. This also explains the philosophy behind our design that for such factories, variables that are not escaped should be assigned while for the function [`escape`](@ref), variables to be escaped should be assigned.

## Manual

```@autodocs
Modules=[Factory]
Order=  [:module,:constant,:type,:macro,:function]
```
