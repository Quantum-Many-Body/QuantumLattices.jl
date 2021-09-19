```@meta
CurrentModule=QuantumLattices
```

# QuantumLattices

*Julia package for the construction of quantum lattice systems.*

Welcome to [QuantumLattices](https://github.com/Quantum-Many-Body/QuantumLattices.jl). Here we provide a general framework to construct the **second-quantized operator-formed Hamiltonian of any quantum lattice system**, with the inputs as simple as its description by the natural language. Combined with [SymPy](https://github.com/JuliaPy/SymPy.jl), this operator-formed Hamiltonian supports **complete symbolic computations**, making it a convenient prerequisite of quantum many-body algorithms, such as **TBA**(tight-bind approximation), **SCMF**(self-consistent mean field theory), **ED**(exact diagonalization), **CPT/VCA**(cluster perturbation theory and variational cluster approach ), **DMRG**(density matrix renormalization group), etc. Generic interfaces are defined to provide a unified access to these algorithms which would be implemented in separate packages. Only minor modifications need be made when users alter from one algorithm to another.

## Introduction

The core of the package is the construction of the **operator representation of a lattice Hamiltonian**. This is based on the following mathematical observations that the operators in a lattice Hamiltonian:
* **act on local Hilbert spaces**, and
* **form an algebra over the complex field**.

The first observation is the starting point of our [**unitcell description framework**](https://quantum-many-body.github.io/QuantumLattices.jl/dev/tutorials/UnitcellDescription/Introduction/) and the second is the mathematical foundation of our **symbolic computing system**.

It is noted that our implementation of the symbolic computation only involves
* the mathematical operations between a scalar and an operator, and
* the mathematical operations between two operators.

The symbolic operations between two scalars are **not** implemented because:
* in condensed matter physics, for many cases, only the numerical values of operators are important because the analytical expressions can be too complicated to analyze or they may even not exist;
* our construction process of the operators and their mathematical operations are **compatible with the [SymPy](https://github.com/JuliaPy/SymPy.jl) package**, therefore, a fully symbolic computation can be achieved by a simple combination of both.

Another major aim of this package is to provide unified interfaces to access all quantum many-body algorithms. Much of the job can be done by the construction of the operator-formed Hamiltonian, which serves as a common input for different algorithms. The remaining stuff concerns mainly with project management, such as result recording, data caching, parameter updating, code logging, dependency managing, etc. Utilities are provided to handle these tasks.

## Package Features

* **Unitcell Description Framework**: by telling the information of the quantum lattice system within a unitcell, the construction of the symbolic representation of the Hamiltonian is just as simple as describing the system in a usual research paper.
* **Complete Symbolic Computation**: with only this package, symbolic computation between operators is realized whereas the coefficient of any operator remains numeric; by integrating it with [SymPy](https://github.com/JuliaPy/SymPy.jl), complete symbolic computation can be achieved and no modifications need be made on the methods in this package.
* **Generic Many-Body Algorithmic Interfaces**: quantum many-body algorithms can be initialized in quite similar ways with only minor modifications needed. Moreover, automatic project management is realized, including that of result recording, data caching, parameter updating, code logging, dependency managing, etc.

## Supported Systems

Four common kinds of systems in condensed matter physics are supported:
* **canonical fermionic systems**
* **canonical/hard-core bosonic systems**
* **SU(2) spin systems**
* **Phononic systems**

Furthermore, other systems can be supported easily by extending the generic "protocols" provided in this package.

## Supported Algorithms

Concrete algorithms are implemented in separate packages (still in progress):
* **[TBA](https://github.com/Quantum-Many-Body/TightBindingApproximation.jl)**: tight-binding approximation for fermionic/bosonic systems;
* **SCMF**: self-consistent mean field theory for fermionic systems;
* **[ED](https://github.com/Quantum-Many-Body/ExactDiagonalization.jl)**: exact diagonalization for fermionic/hard-core-bosonic/spin systems;
* **CPT/VCA**: cluster perturbation theory and variational cluster approach for fermionic systems;
* **DMRG**: density matrix renormalization group for fermionic/hard-core-bosonic/spin systems;
* **[SWT](https://github.com/Quantum-Many-Body/SpinWaveTheory.jl)**: spin wave theory for local spin systems.

## Getting Started
[Tutorials: unitcell description](https://quantum-many-body.github.io/QuantumLattices.jl/dev/tutorials/UnitcellDescription/)

## Python counterpart
[HamiltonianPy](https://github.com/waltergu/HamiltonianPy): in fact, the authors of this Julia package worked on the python package at first and only turned to Julia later.
