```@meta
CurrentModule=Hamiltonian.Utilities.GoodQuantumNumber
```

# Good quantum numbers

Good qunatum numbers can be considered as the conserved labels for the bases of a Hilbert space when a quantum system hosts some symmetries. Here we only implement Abelian good quantum numbers because non-Abelian ones are far more complicated yet much less used. In practice, good quantum numbers can be integers or half integers, therefore, we use real numbers to denote them in this module for simplicity. Independent quantum numbers, such as the particle number and the spin z-component, can coexist at the same time. We use type `QuantumNumber` to represent the complete set of independent ones for a single basis of a Hilbert space, and type `QuantumNumbers` to represent the whole quantum numbers for the total bases.

## QuantumNumber

The abstract type for the complete set of independent good quantum numbers for a single basis.

Main features include:
* function `fieldnames`: get the names of the quantum numbers
* function `periods`: get the periods of the quantum numbers
* arithmetic operations: `+`,`-`,`*`,`⊕`
* hashable: concrete instances can be used as keys for a dict or a set
* iterable: concrete instances are iterable over their values

In particular, `QuantumNumber <: AbstractNamedVector{Float64}`, all features supported by `AbstractNamedVector` are also available for `QuantumNumber`. See also [`AbstractNamedVector`](@ref).

For convenience, **4** kinds of good quantum numbers are predefined in this module, i.e.
* `SQN`: for spin z-component reserved systems
* `PQN`: for particle number reserved systems
* `SPQN`: for both particle number and spin-z component reserved systems
* `Z2QN`: for systems with a ``Z_2`` conservation quantum number

Users who want to define their own ``Z_N``-like quantum numbers must handle the periodicities in the construction function, otherwise, wrong results will be get when arithmetic operations, such as `+` or `-`, are involved. It is highly recommended to use the macro `@quantumnumber` to define your own concrete `QuantumNumber`s.

## QuantumNumbers

The whole quantum numbers for the total bases.

To achieve high efficiency:
* The quantum numbers are stored in a compressed form similiar to that of a CSC/CSR sparse matrix.
* The contents of a `QuantumNumbers` are not an array of some kind of concrete `QuantumNumber`s, but an 2d array of floats with the columns being the values of those concrete `QuantumNumber`s.

Main features include:
* function `eltype`: get the concrete type of the quantum numbers it contains
* arithmetic operations: `+`,`-`,`*`,`^`,`⊗`,`⊕`
* iteration: concrete `QuantumNumber`s it contains will be constructed and returned

## Manual
```@autodocs
Modules=[GoodQuantumNumber]
Order=  [:module,:constant,:type,:macro,:function]
```
