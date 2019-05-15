# QuantumLattices.jl

[![Build Status](https://api.travis-ci.org/Quantum-Many-Body/QuantumLattices.jl.svg?branch=master)](https://travis-ci.org/Quantum-Many-Body/QuantumLattices.jl)
[![codecov](https://codecov.io/gh/Quantum-Many-Body/QuantumLattices.jl/branch/master/graph/badge.svg)](https://codecov.io/gh/Quantum-Many-Body/QuantumLattices.jl)
[![][docs-latest-img]][docs-latest-url]
[![][docs-stable-img]][docs-stable-url]
[![996.icu](https://img.shields.io/badge/link-996.icu-red.svg)](https://996.icu)
[![LICENSE](https://img.shields.io/badge/License-Apache%202.0-blue.svg)](https://opensource.org/licenses/Apache-2.0)
[![LICENSE](https://img.shields.io/badge/license-Anti%20996-blue.svg)](https://github.com/996icu/996.ICU/blob/master/LICENSE)

*Julia package for the construction of quantum lattice systems.*

Welcome to [QuantumLattices](https://github.com/Quantum-Many-Body/QuantumLattices.jl). Here we provide a general framework to construct the **second quantized operator formed Hamiltonian of any quantum lattice system**, with the inputs as simple as its description by natural languages. Combined with [SymPy](https://github.com/JuliaPy/SymPy.jl), this operator formed Hamiltonian supports **complete symbolic computations**, making it a convenient prerequisite of quantum many-body algorithms, such as **TBA**(tight-bind approximation), **SCMF**(self-consistent mean field theory), **ED**(exact diagonalizaiton), **CPT/VCA**(cluster perturbation theory and variational cluster approach ), **DMRG**(density matrix renormalization group), etc. Generic interfaces are defined to give a unified access to these algorithms although their real implementations come in seperate packages. Only minor modifications need be made when users alter from one algorithm to another.

## Installation

In Julia **v1.1+**, please type `]` in the REPL to use the package mode, then type this command:

```julia
pkg> add QuantumLattices
```

## Introduction

The core of the package is the construction of the **operator representations of lattice Hamiltonians**. This is based on the following mathematical observations that the operators in a lattice Hamiltonian:
* **act on local Hilbert spaces**, and
* **form algebras over the complex field**.

The first observation is the starting point of our [**unitcell description framework**](https://quantum-many-body.github.io/QuantumLattices.jl/dev/tutorials/UnitcellDescription/) and the second is the mathmatical foundation of our **symbolic computing system**.

It is noted that our implementation of the symbolic computation only involves
* the mathematical operations between a scalar and an operator, and
* the mathematical operations between two operators.

The symbolic operations between two scalars are **not** implemented becase:
* in condensed matter physics, for many cases, only the numerical values of operators are important because the analytical expressions can be too complicated to analyze or they may even not exist;
* our construction process of the operators and their mathematical operations are **completely compatible with the [SymPy](https://github.com/JuliaPy/SymPy.jl) package**, therefore, a fully symbolic computation can be acheived by a simple combination of both.

Another major aim of this package is to provide unified interfaces to access all quantum-many algorithms. Much of the job can be done by the construction of the operator-formed Hamiltonians, which serves as a common input for different algorithms. The remaining stuff concerns mainly with project management, such as result recording, data caching, parameter updating, code logging, dependency managing, etc. Utilities are provided to handle these tasks.

## Package Features

* **Unitcell Description Framework**: by telling the information of the quantum lattice system within a unitcell, the construction of the symbolic representation of the Hamiltonian is just as simple as describing the system in a usual research paper.
* **Complete symbolic computation**: with only this package, symbolic computation between operators is realized whereas the coeffcient of any operator remains numeric; by integrating it with [SymPy](https://github.com/JuliaPy/SymPy.jl), complete symbolic computation can be acheived and no modifications need be made on the methods in this package.
* **Generic Engine-App Interfaces**: by regarding the relation between algorithms and tasks as that between engines and apps, automatic project management is realized, including that of result recording, data caching, parameter updating, code logging, dependency managing, etc, moreover, all algorithms are initialized in quite similiar ways with only minor modifications needed.

## Supported Systems

Three common kinds of systems in condensed matter physics are perfectly supported:
* **canonical fermionic systems**
* **canonical/hard-core bosonic systems**
* **SU(2) spin systems**

Furthermore, other systems can be supported easily by extending the generic "protocols" provided in this package.

## Supported Algorithms

Concrete algorithms are implemented in seperate packages (still in progess):
* **TBA**: tight-binding approximation for fermionic/bosonic systems;
* **SCMF**: self-consistent mean field theory for fermionic systems;
* **[ED](https://github.com/Quantum-Many-Body/ExactDiagonalization.jl)**: exact diagonalizaiton for fermionic/hard-core-bosonic/spin systems;
* **CPT/VCA**: cluster perturbation theory and variational cluster approach for fermionic systems;
* **DMRG**: density matrix renormalization group for fermionic/hard-core-bosonic/spin systems;
* **LSWT**: linear spin wave theory for local spin systems;
* **FBFMSW**: spin wave theory for flatband ferromagnets.

## Getting Started
[Tutorials: unitcell description](https://quantum-many-body.github.io/QuantumLattices.jl/dev/tutorials/UnitcellDescription/)

## Documentation
- [**LATEST**][docs-latest-url] &mdash; **documentation of the latest version.**
- [**STABLE**][docs-stable-url] &mdash; **documentation of the stable version.**

## Note

Due to the fast development of this package, releases with different minor version numbers are **not** guaranteed to be compatible with previous ones **before** the release of v1.0.0. Comments are welcomed in the github issues.

## Contact
waltergu1989@gmail.com

## Python counterpart
[HamiltonianPy](https://github.com/waltergu/HamiltonianPy): in fact, the authors of this Julia package worked on the python package at first and only turned to Julia later.

[docs-latest-img]: https://img.shields.io/badge/docs-latest-blue.svg
[docs-latest-url]: https://quantum-many-body.github.io/QuantumLattices.jl/latest/
[docs-stable-img]: https://img.shields.io/badge/docs-stable-blue.svg
[docs-stable-url]: https://quantum-many-body.github.io/QuantumLattices.jl/stable/
