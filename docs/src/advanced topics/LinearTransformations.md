```@meta
CurrentModule = QuantumLattices
DocTestFilters = [r"im +[-\+]0\.0[-\+]"]
DocTestSetup = quote
    push!(LOAD_PATH, "../../../src/")
    using QuantumLattices
end
```

# [Linear Transformations](@id LinearTransformations)

[`LinearTransformation`](@ref) was introduced in [Chapter 4](@ref TutorialOperators) as the abstract framework for linear maps on the operator algebra, with [`Permutation`](@ref) as a concrete example. Here we discuss advanced usage: how [`LinearTransformation`](@ref) interacts with [`OperatorGenerator`](@ref) and its internal [`CategorizedGenerator`](@ref).

## CategorizedGenerator

[`CategorizedGenerator`](@ref) is the internal representation used by [`OperatorGenerator`](@ref). It groups the operators of a quantum lattice system into three categories:
- **Constant operators** (`constops`): independent of tunable parameters
- **Alterable operators** (`alterops`): depend on tunable parameters
- **Boundary operators** (`boundops`): associated with boundary twists

When a [`LinearTransformation`](@ref) is applied to an [`OperatorGenerator`](@ref), the [`CategorizedGenerator`](@ref) structure is preserved: each category is transformed independently. This design keeps the internal structure clean and enables efficient parameter updates: only the alterable part needs recomputation when parameters change.

## UnitSubstitution and TabledUnitSubstitution

[`UnitSubstitution`](@ref) and [`TabledUnitSubstitution`](@ref) are [`LinearTransformation`](@ref) subtypes that replace operator indices with concrete values (or other indices). This is the mechanism behind term expansion (Chapter 5): the placeholder orbitals, spins, and sites in coupling patterns are substituted with the concrete values determined by the Hilbert space and bond geometry.

## Matrixization

[`Matrixization`](@ref) (abstract `<: LinearTransformation`) converts operators to matrix representations in a specified basis. The [`matrix`](@ref) function uses a [`Table`](@ref) ([Section 3.5](@ref TutorialTableAndMetric)) to determine the basis ordering.
