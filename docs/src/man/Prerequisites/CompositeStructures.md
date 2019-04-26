```@meta
CurrentModule=QuantumLattices.Prerequisites.CompositeStructures
```

# Composite structures

In principle, Julia is not an object-oriented programming language. For example, only abstract types can be inherited so that subtype cannot inherit fields from their parents. Therefore, Julia prefers composition over inheritance. However, to make a new concrete type behaves much alike another one, tedious reputitions of redifining the generic interfaces are usually not avoidable, especially for the basic types in Julia base. In this module, we implement three such composited types, [`CompositeTuple`](@ref), [`CompositeVector`](@ref) and [`CompositeDict`](@ref), for the sake of future usages.

## CompositeTuple and CompositeNTuple

A composite tuple (ntuple) is a tuple (ntuple) that is implemented by including an ordinary `Tuple` (`NTuple`) as one of its attributes with the name `:contents`.

To take full advantages of the Julia base, the following interfaces are defined:
* inquiry of info: `length`, `eltype`, `hash`
* comparison between objects: `==`, `isequal`
* obtainment of old elements: `getindex`
* iteration: `iterate`, `keys`, `values`, `pairs`
* construction of new objects: `reverse`

Composite ntuples are suited for the situations where other attributes are not affected by the modification of the elements. Note that arithmatic operations and logical operations excluding `==` and `isequal` are not supported. Besides, a composite ntuple is **not** a tuple since Julia has no abstract tuples.

## CompositeVector

A composite vector is a vector that is implemented by including an ordinary `Vector` as one of its attributes with the name `:contents`.

To take full advantages of the Julia base, the following interfaces are redefined:
* inquiry of info: `size`, `length`
* comparison between objects: `==`, `isequal`
* obtainment of old elements: `getindex`
* modification of old elements: `setindex!`
* addition of new elements: `push!`, `pushfirst!`, `insert!`, `append!`, `prepend!`
* removal of old elements: `splice!`, `deleteat!`, `pop!`, `popfirst!`, `empty!`
* construction of new objects: `empty`, `reverse`, `similar`
* iteration: `iterate`, `keys`, `values`, `pairs`

Composite vectors are suited for the situations where other attributes are not affected by the modification of the elements. Note that arithmatic operations and logical operations excluding `==` and `isequal` are not supported.

## CompositeDict

A composite dict is a dict that is implemented by including an ordinary `Dict` as one of its attributes with the name `:contents`.

To take full advantages of the Julia base, the following interfaces are redefined:
* inquiry of info: `isempty`, `length`, `haskey`, `in`
* comparison between objects: `==`, `isequal`
* obtainment of old elements: `get`, `getkey`, `getindex`
* modification and addition of elements: `push!`, `get!`, `setindex!`
* removal of old elements: `pop!`, `delete!`, `empty!`
* construction of new objects: `merge`, `empty`
* iteration: `iterate`, `keys`, `values`, `pairs`

As is similar to composite vectors, composite dicts are suited for the situations where other attributes are not affected by the modification of the elements.

## NamedContainer

[`NamedContainer`](@ref) is just a wrapper (type alias) of Julia NamedTuple, but not a composite type.

Julia NamedTuple is useful to keep type stability of codes when we deal with unhomogenous immutable dict-like objects, but its default constructor is not so convenient becase the names and contents must be assigned pair by pair in a pair of parentheses explicitly. Therefore, we define a type alias of NamedTuple under the name of [`NamedContainer`](@ref), so that we can construct a NamedTuple by the usual-formed constructor [`NamedContainer`](@ref).

## Manul

```@autodocs
Modules=[CompositeStructures]
Order=  [:module,:constant,:type,:macro,:function]
```
