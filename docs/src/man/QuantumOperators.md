```@meta
CurrentModule = QuantumLattices.QuantumOperators
```

# Quantum Operators

*Module for the operator algebra of quantum lattice systems.*

The `QuantumOperators` module implements the concrete operator algebra, the types that represent individual operators, their sums, products, and collections. These types form the algebraic backbone described in tutorial chapter 5.

**Key types:** [`Operator`](@ref), [`OperatorSum`](@ref), [`OperatorProd`](@ref), [`OperatorSet`](@ref), [`Operators`](@ref).

**Key operations:** Addition, scalar multiplication, and operator multiplication are supported between all operator types. Linear transformations ([`LinearTransformation`](@ref), [`Permutation`](@ref), [`UnitSubstitution`](@ref), [`Matrixization`](@ref)) act systematically on the algebra.

Quantum operators form an algebra over a field, i.e., a vector space equipped with a bilinear operation (multiplication) defined between vectors. There are three basic operations: scalar multiplication between a scalar and a quantum operator, the usual addition, and the usual multiplication between quantum operators. More complicated operations can be composed from these basic ones.

## OperatorIndex

[`OperatorIndex`](@ref) is the building block of quantum operators, which specifies the basis of the vector space of the corresponding algebra.

## OperatorProd and OperatorSum

[`OperatorProd`](@ref) defines the product operator as an entity of basis quantum operators while [`OperatorSum`](@ref) defines the summation as an entity of [`OperatorProd`](@ref)s. Both of them are subtypes of [`QuantumOperator`](@ref), which is the abstract type for all quantum operators.

An [`OperatorProd`](@ref) must have two predefined contents:
- `value::Number`: the coefficient of the quantum operator
- `id::ID`: the id of the quantum operator

Arithmetic operations (`+`, `-`, `*`, `/`) between a scalar, an [`OperatorProd`](@ref) or an [`OperatorSum`](@ref) are defined. See Manual for details.

In addition to [`Operator`](@ref) and [`Operators`](@ref) (covered in [Chapter 4](@ref TutorialOperators)), this module provides specialized container types: [`OperatorSum`](@ref) (explicit sum), [`OperatorProd`](@ref) (product), and [`OperatorSet`](@ref) (unordered collection). It also defines the linear transformation framework: [`LinearTransformation`](@ref), [`Permutation`](@ref), [`UnitSubstitution`](@ref), [`TabledUnitSubstitution`](@ref), and [`Matrixization`](@ref).

## Manual

```@autodocs
Modules = [QuantumOperators]
Order = [:module, :constant, :type, :macro, :function]
```
