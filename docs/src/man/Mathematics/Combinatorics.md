```@meta
CurrentModule=QuantumLattices.Mathematics.Combinatorics
```

# Combinatorics

This module implements the combinations and permutations of an indexable object, with duplicate elements allowed or not. Compared to another Julia package [Combinatorics](https://github.com/JuliaMath/Combinatorics.jl), the iterators return tuples instead of vectors, which greatly decreases the momory allocation times and improves the code efficiency.

## AbstractCombinatorics

[`AbstractCombinatorics{M,C}`](@ref) is the abstract type of all combinatoric algorithms. It has two type parameters:
* `M`: the number of elements to be taken
* `C`: the type of the collection of candidate elements
To avoid momery allocation, the iteration of a concrete combinatoric algorithm returns a tuple, whose length is `M` and eltype is `eltype(C)`.

### Combinations and DulCombinations

[`Combinations{M,C}`](@ref) and [`DulCombinations`](@ref) generate all the combinations of `M` elements from an indexable collection whose type is `C`, with the differences being that the former forbids duplicate elements in the combinations while the latter allows.

### Permutations and DulPermutations

[`Permutations{M,C}`](@ref) and [`DulPermutations`](@ref) generate all the permutations of `M` elements from an indexable collection whose type is `C`, with the differences being that the former forbids duplicate elements in the permutations while the latter allows.

## Manul

```@autodocs
Modules=[Combinatorics]
Order=  [:module,:constant,:type,:macro,:function]
```
