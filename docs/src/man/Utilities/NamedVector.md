```@meta
CurrentModule=Hamiltonian.Utilities.NamedVector
```

# Named vector

A named vector is similiar to a named tuple, which associate each of its values with a name. Although the names of a named vector cannot be changed, the values can be modified if needed. In contrast to the predefined [`NamedTuple`](https://docs.julialang.org/en/v1/base/base/#Core.NamedTuple) in Julia, which employs the names as type parameters, we just implement a named vector as a composite struct equipped with the `getindex` and `setindex!` functions, with the fieldnames being its names. This simple implementation makes it possible to define your own concrete named vector with any of your preferred type names, and ensures that all instances of a certain concrete named vector share the same names. Therefore, if you are familiar with Python, you will find that our named vector is more qualified to be the counterpart of the [`namedtuple`](https://docs.python.org/3.7/library/collections.html#collections.namedtuple) in Python than the default Julia implementation. Specifically, we also define a macro [`@namedvector`](@ref) as the type factory to help users to define their own concrete named vectors. Last but not least important, it is also worth noted that **a named vector is not a vector**, as is similar to that a named tuple is not a tuple in Julia. This results from our basic expectation that a named vector should be more like a tuple other than a vector so that not all operations valid to vectors are also valid to named vectors.

## AbstractNamedVector
[`AbstractNamedVector`](@ref) defines the abstract type for all concrete named vectors.

Main features include:
* Values can be accessed or modified either by the `.` operator or by the `[]` operator.
* Comparisons, such as `≡`, `≢`, `==`, `≠`, `>`, `<`, `≥`, `≤` are supported. Therefore a vector of named vectors can be sorted by the default `sort` function.
* Hash is supported by `hash`. Therefore, a named vector can be used as the key of a dict or set.
* Iteration over its fieldnames is supported by `keys`, over its values is supported by `values`, over its field-value pairs is supported by `pairs`. A reverse iteration is also supported.

To subtype it, please note:
1. A concrete type can be either mutable or immutable as you need, but all its fields should be of the same type. A recommended template for the subtype is
   ```julia
   [mutable] struct YourNamedVector{T} <: AbstractNamedVector{T}
       filedname1::T
       filedname2::T
       ...
   end
   ```
2. It is recommended to overload the `Base.fieldnames` function for concrete subtypes to ensure type stability and improve efficiency, which though is not a necessity. A template for such an overloading is
   ```julia
   Base.fieldnames(Type{YourNamedVector})=(:fieldname1,:fieldname2,...)
   ```
3. For all concrete subtypes, if inner constructors are defined, the one which has the same interface with the default one must be implemented. Otherwise, some functionalities will not work.
4. Arithmetic operations, such as `+`, `-`, `*`, `/`, `%`, `÷`, etc. are **NOT** supported. However, an efficient [`map`](@ref NamedVector.map) function  is implemented, which can help users do the overloadings of these operations.

## Manual
```@autodocs
Modules=[NamedVector]
Order=  [:module,:constant,:type,:macro,:function]
```
