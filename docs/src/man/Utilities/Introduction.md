```@meta
CurrentModule=Hamiltonian.Utilities
```

# Introduction

This module contains the utilities of the Hamiltonian package.

The constants, types, macros, functions or submodules defined in this module will **NOT** be exported by the package. Instead, they serve as the prerequisites and fundamentals.
The range of the utilities' contents are quite wide, but basically, they fall into four categories:
* Global constants and miscellaneous tiny useful functions;
* Generic functions that are extended by other parts of the package;
* Basic data structures as supplements to the Julia.Base and other common packages;
* Basic mathematical algorithms or structures that are needed by or specific to this package.
The first category is contained in the main body of this module, while the others come in separate submodules.

## Constants and functions

All the following constants and functions in this section are defined in the main body and are exported by this module.

```@docs
atol
rtol
Float
decimaltostr
ordinal
delta
```

## Generic interfaces

For details, please refer to the page [Interfaces](@ref):
```@contents
Pages=[
    "Interfaces.md",
    ]
Depth=2
```

## Basic structures

Here lists the table of contents of the basic data structures that are supplements to the Julia.Base and other common packages:
```@contents
Pages=[
    "TypeTraits.md",
    "Factories.md",
    "CompositeStructures.md",
    "Trees.md",
    "NamedVectors.md",
    ]
Depth=2
```

## Mathematical supplements

The following lists the submodules containing the mathematical algorithms or structures that are needed by or specific to this package:
```@contents
Pages=[
    "Combinatorics.md",
    "AlgebraOverFields.md",
    "QuantumNumbers.md",
    ]
Depth=2
```
