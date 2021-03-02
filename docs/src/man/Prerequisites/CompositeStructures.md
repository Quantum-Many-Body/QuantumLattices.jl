```@meta
CurrentModule = QuantumLattices.Prerequisites.CompositeStructures
```

```@setup compositestructrues
push!(LOAD_PATH, "../../../../src/")
using QuantumLattices.Prerequisites.CompositeStructures
```

# Composite structures

In principle, Julia is not an object-oriented programming language. For example, only abstract types can be inherited so that subtype cannot inherit fields from their parents. Therefore, Julia prefers composition over inheritance. However, to make a new concrete type behaves much alike another one, tedious repetitions of redefining the generic interfaces are usually not avoidable, especially for the basic types in Julia base. In this module, we implement three such composited types, [`CompositeTuple`](@ref), [`CompositeVector`](@ref) and [`CompositeDict`](@ref), for the sake of future usages. Besides, [`NamedContainer`](@ref), as a wrapper of Julia `NamedTuple`, is also provided here so that the construction of a Julia `NamedTuple` can be more flexible over the standard `(name=value, ... )` syntax.

## CompositeTuple and CompositeNTuple

A composite tuple (ntuple) can be considered as a tuple (ntuple) that is implemented by including an ordinary [`Tuple`](https://docs.julialang.org/en/v1/base/base/#Core.Tuple)([`NTuple`](https://docs.julialang.org/en/v1/manual/types/#Vararg-Tuple-Types)) as its data attribute.

To take full advantages of the Julia base, the following interfaces are defined:
* inquiry of info: `length`, `eltype`, `hash`
* comparison between objects: `==`, `isequal`
* obtainment of old elements: `getindex`
* iteration: `iterate`, `keys`, `values`, `pairs`
* construction of new objects: `reverse`

Note that arithmetic operations and logical operations excluding `==` and `isequal` are not supported. Besides, a composite tuple is **not** a tuple since Julia has no abstract tuples.

## CompositeVector

A composite vector can be considered as a vector that is implemented by including a concrete subtype of [`AbstractVector`](https://docs.julialang.org/en/v1/base/arrays/#Base.AbstractVector) as its data attribute, and it itself is a subtype of [`AbstractVector`](https://docs.julialang.org/en/v1/base/arrays/#Base.AbstractVector).

To take full advantages of the Julia base, the following interfaces are redefined:
* inquiry of info: `size`, `length`
* comparison between objects: `==`, `isequal`
* obtainment of old elements: `getindex`
* operation of old elements: `setindex!`
* addition of new elements: `push!`, `pushfirst!`, `insert!`, `append!`, `prepend!`
* removal of old elements: `splice!`, `deleteat!`, `pop!`, `popfirst!`, `empty!`
* construction of new objects: `empty`, `reverse`, `similar`
* iteration: `iterate`, `keys`, `values`, `pairs`

Note that arithmetic operations and logical operations excluding `==` and `isequal` are not supported.

## CompositeDict

A composite dict can be considered as a dict that is implemented by including a concrete subtype of [`AbstractDict`](https://docs.julialang.org/en/v1/base/collections/#Base.AbstractDict) as its data attribute and it itself is a subtype of [`AbstractDict`](https://docs.julialang.org/en/v1/base/collections/#Base.AbstractDict).

To take full advantages of the Julia base, the following interfaces are redefined:
* inquiry of info: `isempty`, `length`, `haskey`, `in`
* comparison between objects: `==`, `isequal`
* obtainment of old elements: `get`, `getkey`, `getindex`
* operation and addition of elements: `push!`, `get!`, `setindex!`
* removal of old elements: `pop!`, `delete!`, `empty!`
* construction of new objects: `merge`, `empty`
* iteration: `iterate`, `keys`, `values`, `pairs`

## NamedContainer

[`NamedContainer`](@ref) is just a wrapper (type alias) of Julia NamedTuple, but not a composite type.

Julia NamedTuple is useful to keep type stability of codes when we deal with inhomogeneous immutable dict-like objects, but its default constructor is not so convenient because the names and contents must be assigned pair by pair in a pair of parentheses explicitly. Therefore, we define a type alias of NamedTuple under the name of [`NamedContainer`](@ref), so that we can construct a NamedTuple by the usual-formed constructor [`NamedContainer`](@ref), e.g.
```@example compositestructrues
NamedContainer{(:a, :b)}((1, 2))
```

## Manual

```@autodocs
Modules = [CompositeStructures]
Order = [:module, :constant, :type, :macro, :function]
```
