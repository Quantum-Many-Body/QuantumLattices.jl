```@meta
CurrentModule = QuantumLattices.Essentials.QuantumAlgebras
```

```@setup quantumalgebras
push!(LOAD_PATH, "../../../../src/")
using QuantumLattices.Essentials.QuantumAlgebras
```

# Quantum algebras

An algebra over a field is a vector space over that field, in which a bilinear operator (often called the "multiplication") between vectors is defined.

With the help of the structure constants of the algebra, the result of the bilinear operation between any arbitrary two vectors can be expressed by a sum of individual ones. Therefore, in principle, an algebra can be represented by the complete basis set of its corresponding vector space and a rank-3 tensor encapsulating its structure constants. It is noted that the "bilinear operation" is not restricted to the usual multiplication only. For example, it is the commutator, which is a composition of the usual multiplication and subtraction (for any A and B, the commutator [A, B] is defined as [A, B]≝AB-BA) that serves as the bilinear operator for Lie algebras. In this module, for scalars in the field and elements in the algebra, we only provide the interfaces of the scalar multiplication (including the scalar division) between a scalar and an element, the addition (including the subtraction) and the usual multiplication between two elements. Other complicated operations should be composed from these basic ones.

## SimpleID and ID

[`SimpleID`](@ref) is the building block of the id system of an algebra over a field, while [`ID`](@ref) defines the specific identifier of an element in that algebra.

Generally, the usual multiplication between two elements of an algebra is not commutable, and the rank of the multiplication is just the add-up before the simplification with the help of the algebra structure. We thus use a simple id to denote a single basis of the corresponding vector space, and an id to denote the identifier of an element. With the help of the direct product (`⊗`) of two ids, an over complete id system designed for the whole algebra is constructed. This id system is redundant because it does not reflects the structure constants of the algebra, which reduces independent basis elements. Extra mechanisms should be provided to kill this redundancy, which goes beyond the current module. Users should define them themselves.

## Element and Elements

[`Element`](@ref) defines a single element of an algebra while [`Elements`](@ref) defines an expression composed of several elements from an algebra.

An [`Element`](@ref) must have two predefined contents:
- `value::Number`: the coefficient of the element
- `id::ID`: the id of the element

Arithmetic operations (`+`, `-`, `*`, `/`) between a scalar, an [`Element`](@ref) or an [`Elements`](@ref) is defined. See Manual for details.

## Manual

```@autodocs
Modules = [QuantumAlgebras]
Order = [:module, :constant, :type, :macro, :function]
```
