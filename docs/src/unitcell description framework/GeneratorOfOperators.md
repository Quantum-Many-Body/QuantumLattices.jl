```@meta
CurrentModule = QuantumLattices
DocTestFilters = [r"im +[-\+]0\.0[-\+]"]
DocTestSetup = quote
    push!(LOAD_PATH, "../../../src/")
    using QuantumLattices
end
```

```@setup gen
using QuantumLattices
using SymPy: Sym, symbols
```

# Generator of operators

In previous pages we have introduced the essential elements in detail of a quantum lattice system in the unitcell description framework. In this page, we will discuss how such elements can be incorporated to get the operator representation of the lattice Hamiltonian.

## Operator-formed lattice Hamiltonians

**The essential elements to obtain an operator-formed lattice Hamiltonian are 1) terms, 2) bonds and 3) Hilbert space**.

Accordingly, [`OperatorGenerator`](@ref) is the type to incorporate all these elements:
```julia
OperatorGenerator(terms::Tuple{Vararg{Term}}, bonds::Vector{<:Bond}, hilbert::Hilbert)
```
Then, the operator representation of the lattice Hamiltonian can be obtained by the [`expand`](@ref) function:
```julia
expand(gen::OperatorGenerator)
```
which is based on the expansion of terms introduced in the last section of the previous page [Couplings among different degrees of freedom](@ref).

Now, let's go back to the example proposed in the page of [Introduction](@ref UnitcellDescriptionIntroduction):
```@repl gen
lattice = Lattice([zero(Sym)], [one(Sym)]);
hilbert = Hilbert(site=>Fock{:f}(1, 2) for site=1:length(lattice));
t = Hopping(:t, symbols("t", real=true), 1);
U = Hubbard(:U, symbols("U", real=true));
operators = expand(OperatorGenerator((t, U), bonds(lattice, 1), hilbert))
```
Note to run the above codes, `SymPy.Sym` and `SymPy.symbols` should be imported first.

The LaTeX formatted outputs will be discussed in the page of [LaTeX formatted outputs](@ref).

Here, we have two comments:
* In the construction of [`OperatorGenerator`](@ref), a vector of [`Bond`](@ref)s other than a [`Lattice`](@ref) is required. Indeed, it is more flexible to control the generated operators by passing a vector of [`Bond`](@ref)s. By restricting the range of the input bonds, a certain subset of the lattice Hamiltonian can be easily obtained, which is sometimes useful in some quantum many-body algorithms.
* The Hilbert space is essential because in general not all information of the local generators of a concrete quantum lattice system is contained in the terms. Remind: 1) the statistics of particles or the total spin of local spins is usually omitted in the coupling pattern of a term, and 2) the total number of local orbital/spin degrees of freedom, or the dimension of lattice vibrations, cannot be determined solely by the information provided by the coupling pattern of a term. Therefore, in order to obtain the lattice Hamiltonian, such incomplete information must be supplemented by the Hilbert space.

## Parameters tuning of lattice Hamiltonians

It is customary to tune the parameters of a lattice Hamiltonian. This can be achieved by the [`update!`](@ref) function exported by this package:
```julia
update!(gen::OperatorGenerator; termid₁=termvalue₁, termid₁=termvalue₂, ...)
```
Here, `termidᵢ` (i = 1, 2, ...) is the id of a term in the lattice Hamiltonian, and `termvalueᵢ` is the new overall coefficient of this term.

The parameters of the terms in an [`OperatorGenerator`](@ref) can be requested by the type [`Parameters`](@ref) (a type alias of `NamedTuple`):
```julia
Parameters(gen::OperatorGenerator) -> Parameters
```

Let's see an example.
```jldoctest
julia> lattice = Lattice([0.0], [1.0]);

julia> hilbert = Hilbert(site=>Fock{:f}(1, 2) for site=1:length(lattice));

julia> t = Hopping(:t, 1.0, 1);

julia> gen = OperatorGenerator((t,), bonds(lattice, 1), hilbert);

julia> Parameters(gen)
(t = 1.0,)

julia> expand(gen)
Operators with 4 Operator
  Operator(1.0, 𝕗(2, 1, -1//2, 2, [1.0], [0.0]), 𝕗(1, 1, -1//2, 1, [0.0], [0.0]))
  Operator(1.0, 𝕗(1, 1, -1//2, 2, [0.0], [0.0]), 𝕗(2, 1, -1//2, 1, [1.0], [0.0]))
  Operator(1.0, 𝕗(2, 1, 1//2, 2, [1.0], [0.0]), 𝕗(1, 1, 1//2, 1, [0.0], [0.0]))
  Operator(1.0, 𝕗(1, 1, 1//2, 2, [0.0], [0.0]), 𝕗(2, 1, 1//2, 1, [1.0], [0.0]))

julia> update!(gen; t=2.0);

julia> Parameters(gen)
(t = 2.0,)

julia> expand(gen)
Operators with 4 Operator
  Operator(2.0, 𝕗(2, 1, -1//2, 2, [1.0], [0.0]), 𝕗(1, 1, -1//2, 1, [0.0], [0.0]))
  Operator(2.0, 𝕗(1, 1, -1//2, 2, [0.0], [0.0]), 𝕗(2, 1, -1//2, 1, [1.0], [0.0]))
  Operator(2.0, 𝕗(2, 1, 1//2, 2, [1.0], [0.0]), 𝕗(1, 1, 1//2, 1, [0.0], [0.0]))
  Operator(2.0, 𝕗(1, 1, 1//2, 2, [0.0], [0.0]), 𝕗(2, 1, 1//2, 1, [1.0], [0.0]))
```
