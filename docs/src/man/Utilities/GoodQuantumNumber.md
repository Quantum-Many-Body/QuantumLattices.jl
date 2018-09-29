```@meta
CurrentModule=Hamiltonian.Utilities.GoodQuantumNumber
```

# Good quantum numbers

Good qunatum numbers can be considered as the conserved labels for the bases of a Hilbert space when a quantum system hosts some symmetries. Here we only implement Abelian good quantum numbers because non-Abelian ones are far more complicated yet much less used. In practice, good quantum numbers are usually integers or half integers. Whereas we use real numbers to represent half integers in this module for simplicity. Independent quantum numbers, such as the particle number and the spin z-component, can coexist at the same time. We use type `QuantumNumber` to represent the complete set of independent ones for a single basis of a Hilbert space, and type `QuantumNumbers` to represent the whole quantum numbers for the total bases.

## QuantumNumber

The base type for the complete set of independent good quantum numbers for a single basis. Main features include:
* function `names`: get the names of the quantum numbers
* function `periods`: get the periods of the quantum numbers
* arithmetic operations: `+`,`-`,`*`
* hashable: concrete instances can be used as keys for a dict

## QuantumNumbers

The whole quantum numbers for the total bases.

To achieve high efficiency:
* The quantum numbers are stored in a compressed form similiar to that of a CSC/CSR sparse matrix.
* The contents of a `QuantumNumbers` are not an array of some kind of concrete `QuantumNumber`s, but an 2d array of integers or floats with the columns being the values of those concrete `QuantumNumber`s.

Main features include:
* function `qntype`: get the concrete type of the quantum numbers it contains
* arithmetic operations: `+`,`-`,`*`,`^`,`⊗`,`⊕`
* iteration:

```@autodocs
Modules=[GoodQuantumNumber]
Order=  [:module,:constant,:type,:macro,:function]
```
