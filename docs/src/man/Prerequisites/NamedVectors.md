```@meta
CurrentModule=QuantumLattices.Prerequisites.NamedVectors
```

```@setup namedvectors
push!(LOAD_PATH,"../../../../src/")
using QuantumLattices.Prerequisites.NamedVectors
```

# Named vectors

A named vector is similar to a named tuple, which associate each of its values with a name. Although the names of a named vector cannot be changed, the values can be modified if needed. In contrast to the predefined [`NamedTuple`](https://docs.julialang.org/en/v1/base/base/#Core.NamedTuple) in Julia, which employs the names as type parameters, we just implement a named vector as a composite struct equipped with the `getindex` and `setindex!` functions, with the fieldnames being its names. This simple implementation makes it possible to define your own concrete named vector with any of your preferred type names, and ensures that all instances of a certain concrete named vector share the same names. Therefore, if you are familiar with Python, you will find that our named vector is more qualified to be the counterpart of the [`namedtuple`](https://docs.python.org/3.7/library/collections.html#collections.namedtuple) in Python than the default Julia implementation. Last but not least important, it is also worth noted that **a named vector is not a vector**, as is similar to that a named tuple is not a tuple in Julia. This results from our basic expectation that a named vector should be more like a tuple other than a vector so that not all operations valid to vectors are also valid to named vectors.

## NamedVector

[`NamedVector`](@ref) defines the abstract type for all concrete named vectors.

Main features include:
* Values can be accessed or modified either by the `.` operator or by the `[]` operator.
* Comparisons, such as `≡`, `≢`, `==`, `≠`, `>`, `<`, `≥`, `≤` are supported. Therefore a vector of named vectors can be sorted by the default `sort` function.
* Hash is supported by `hash`. Therefore, a named vector can be used as the key of a dict or set.
* Iteration over its fieldnames is supported by `keys`, over its values is supported by `values`, over its field-value pairs is supported by `pairs`.
  A reverse iteration is also supported.

To subtype it, please note:
1. A concrete type can be either mutable or immutable as you need, which is different from tuples.
2. The fields of a concrete type can be of the same type or not. For the former, we denote the named vector as "homogeneous" while for the latter as "inhomogeneous".
   For homogeneous ones, we define a sub abstract type, [`HomoNamedVector`](@ref) for further optimization of the default methods. See [HomoNamedVector](@ref) below.
3. It is recommended to overload the `Base.fieldnames` function for concrete subtypes to ensure type stability and improve efficiency, which though is not a necessity.
   A template for such an overloading is
   ```julia
   Base.fieldnames(Type{<:YourNamedVector})=(:fieldname1,:fieldname2,...)
   ```
4. For all concrete subtypes, if inner constructors are defined, the one which has the same interface with the default one must be implemented.
   Otherwise, some functionalities will not work.
5. Arithmetic operations, such as `+`, `-`, `*`, `/`, `%`, `÷`, etc. are **not** supported.
   However, the function [`map`](@ref) is implemented, which can help users do the overloadings of these operations.

We define a macro [`@namedvector`](@ref) as the type factory to decorate a "raw" struct to be a subtype of `NamedVector`. Here, "raw" means the struct to be decorated has no explicit supertype other than `Any`, neither inner constructors as well. For example,
```@example namedvectors
# decorate an immutable "raw" struct
@namedvector struct InHomoIMNV
    scope::String
    site::Int
end

# decorate a mutable "raw" struct
@namedvector mutable struct InHomoMNV
    scope::String
    site::Int
end
```
This macro encapsulate the overloading of `Base.fieldnames`, and you have no need to do this by hand any more.

## HomoNamedVector

[`HomoNamedVector`](@ref) is the subtype of [`NamedVector`](@ref) that of all its fields share the same type. Compared to `NamedVector`, one more default method is implemented with [`HomoNamedVector`](@ref), i.e. `eltype`, which returns the type of its fields. This function ensures the type stability of all the methods that involves an iteration of the field values of a named vector. Therefore, homogeneous named vector are usually more efficient than inhomogeneous ones. Use homogeneous ones as much as possible unless the code efficiency does not matter.

To subtype [`HomoNamedVector`](@ref), all the suggestions mentioned in the previous subsection for [`NamedVector`](@ref) also applies. A recommended template for a subtype is
```julia
[mutable] struct YourNamedVector{T} <: HomoNamedVector{T}
    filed1::T
    filed2::T
    ...
end
```
We also provide a macro [`@homonamedvector`](@ref) to help the definition of concrete homogeneous named vector, where you only need specify the type name, field names, data type and optionally whether the subtype is mutable. For example,
```@example namedvectors
# construct an immutable homogeneous named vector without type parameters
@homonamedvector HomoNVWithoutParameter (:scope,:site) Int mutable=false

# construct a mutable homogeneous named vector with type parameters
@homonamedvector HomoNVWithParameter (:scope,:site) (<:Real) mutable=true
```
This macro also integrates the `Base.fieldnames` function, thus its overloading by hand is on longer needed.

## Manual

```@autodocs
Modules = [NamedVectors]
Order =  [:module, :constant, :type, :macro, :function]
```
