```@meta
CurrentModule = QuantumLattices.Essentials.QuantumOperators
```

```@setup quantumalgebras
push!(LOAD_PATH, "../../../../src/")
using QuantumLattices.Essentials.QuantumOperators
```

# Quantum operators

Quantum operators form an algebra over a field, which are vector spaces with a bilinear operation (often called the "multiplication") between vectors defined.

With the help of the structure constants of the algebra, the result of the bilinear operation between any arbitrary two vectors can be expressed by a sum of individual ones. Therefore, in principle, an algebra can be represented by the complete basis set of its corresponding vector space and a rank-3 tensor encapsulating its structure constants. It is noted that the "bilinear operation" is not restricted to the usual multiplication. For example, it is the commutator, which is a composition of the usual multiplication and subtraction (for any A and B, the commutator [A, B] is defined as [A, B]≝AB-BA) that serves as the bilinear operator for Lie algebras.

In general, there are three basic operations on quantum operators, i.e. the scalar multiplication between a scalar and a quantum operator, the usual addition and the usual multiplication between quantum operators. Other complicated operations can be composed from these basic ones. These basic operations are implemented in this module.

## SingularID and ID

[`SingularID`](@ref) is the building block of the id system of quantum operators, which specifies the basis of the vector space of the corresponding algebra. On the other hand, [`ID`](@ref) is the direct product (`⊗`) of several [`SingularID`](@ref)s, which could specify the product of the basis quantum operators. This defines an over complete id system for the quantum operators because it does not reflects the structure constants of the corresponding algebra, which would reduce the number of independent basis quantum operators. Extra mechanisms should be provided to kill this redundancy, which goes beyond the current module. Users should define them themselves.

## OperatorProd and OperatorSum

[`OperatorProd`](@ref) defines the product operator as an entity of basis quantum operators while [`OperatorSum`](@ref) defines the summation as an entity of [`OperatorProd`](@ref)s. Both of them are subtypes of [`QuantumOperator`](@ref), which is the abstract type for all quantum operators.

An [`OperatorProd`](@ref) must have two predefined contents:
- `value::Number`: the coefficient of the quantum operator
- `id::ID`: the id of the quantum operator

Arithmetic operations (`+`, `-`, `*`, `/`) between a scalar, an [`OperatorProd`](@ref) or an [`OperatorSum`](@ref) is defined. See Manual for details.

## Manual

```@autodocs
Modules = [QuantumOperators]
Order = [:module, :constant, :type, :macro, :function]
```
