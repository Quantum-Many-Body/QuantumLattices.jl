```@meta
CurrentModule=QuantumLattices.Prerequisites.TypeTraits
```

# Type traits

This module defines generic type traits that are useful to the package.

Julia does not support multi-inheritance, which is sometimes not convenient. A way around this is to use traits, i.e. by utilizing the dispatch on certain (singleton) types known as traits to simulate multi-inheritance. Although this method cannot aviod small repetitive codes, it suits methods well that are complicated and lengthy.

## EfficientOperations

`EfficientOperations` defines efficient operations such as `==/isequal`, `</isless`, `isapprox`, `replace`, etc, that ensure type stability.

Type stability is the key of Julia to improve the code efficiency. However, it cannot be ensured in some unexpected cases, especially where an iterator is involved. For example, the following codes appears type unstable:
```julia
function Base.:(==)(o1::AbstractType,o2::AbstractType)
    n1,n2=o1|>typeof|>fieldcount,o2|>typeof|>fieldcount
    n1==n2 ? all(getfield(o1,i)==getfield(o2,i) for i=1:n1) : false
end
```
Methods like above are common when we design abstract types, but they are not type stable. To get rid of it, the generated function trick can be used:
```julia
@generated function Base.:(==)(o1::AbstractType,o2::AbstractType)
    n1,n2=o1|>typeof|>fieldcount,o2|>typeof|>fieldcount
    if n1==n2
        expr=:(getfield(o1,1)==getfield(o2,1))
        for i=2:fcount
            expr=Expr(:&&,expr,:(getfield(o1,$i)==getfield(o2,$i)))
        end
        return expr
    else
        return :(false)
    end
end
```
Then type stability can be ensured. We use this trick to implement the methods such as `==/isequal`, `</isless`, `isapprox`, `replace`, etc, with the trait `efficientoperations::EfficientOperations`. Other types can resort to these methods by passing [`efficientoperations`](@ref) as the first argument.

## MemoryOrder

`MemoryOrder` provides the convertions, [`subtoind`](@ref) and [`indtosub`](@ref), between a Cartesian index represented by a tuple and a linear index represented by an integer. C/C++ order or Fortran order can be specified, though the constant instances [`corder`](@ref) or [`forder`](@ref) of singleton types `COrder` and `FOrder`, which are both subtypes of the abstract type `MemoryOrder`.

## Manual

```@autodocs
Modules=[TypeTraits]
Order=  [:module,:constant,:type,:macro,:function]
```
