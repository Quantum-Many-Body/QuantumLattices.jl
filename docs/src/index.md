# QuantumLattices.jl

[![CI](https://github.com/Quantum-Many-Body/QuantumLattices.jl/actions/workflows/CI.yml/badge.svg)](https://github.com/Quantum-Many-Body/QuantumLattices.jl/actions/workflows/CI.yml)
[![codecov](https://codecov.io/gh/Quantum-Many-Body/QuantumLattices.jl/branch/master/graph/badge.svg)](https://codecov.io/gh/Quantum-Many-Body/QuantumLattices.jl)
[![](https://img.shields.io/badge/docs-latest-blue.svg)](https://quantum-many-body.github.io/QuantumLattices.jl/latest/)
[![](https://img.shields.io/badge/docs-stable-blue.svg)](https://quantum-many-body.github.io/QuantumLattices.jl/stable/)
[![996.icu](https://img.shields.io/badge/link-996.icu-red.svg)](https://996.icu)
[![LICENSE](https://img.shields.io/badge/License-Apache%202.0-blue.svg)](https://opensource.org/licenses/Apache-2.0)
[![LICENSE](https://img.shields.io/badge/license-Anti%20996-blue.svg)](https://github.com/996icu/996.ICU/blob/master/LICENSE)
[![Code Style: Blue](https://img.shields.io/badge/code%20style-blue-4495d1.svg)](https://github.com/invenia/BlueStyle)
[![ColPrac: Contributor's Guide on Collaborative Practices for Community Packages](https://img.shields.io/badge/ColPrac-Contributor's%20Guide-blueviolet)](https://github.com/SciML/ColPrac)

*Julia package for the construction of quantum lattice systems.*

Welcome to **[QuantumLattices](https://github.com/Quantum-Many-Body/QuantumLattices.jl)**. Here we provide a general framework to construct the **second-quantized operator-formed Hamiltonian of any quantum lattice system**, with the inputs as simple as its description by the natural language. This operator-formed Hamiltonian supports **complete symbolic computations** when combined with **[SymPy](https://github.com/JuliaPy/SymPy.jl)**, and can serve as a convenient **frontend of quantum many-body algorithms**, such as **[TBA (tight-bind approximation)](https://github.com/Quantum-Many-Body/TightBindingApproximation.jl)**, **[LSWT (linear spin wave theory)](https://github.com/Quantum-Many-Body/SpinWaveTheory.jl)**, **SCMF (self-consistent mean field theory)**, **[ED (exact diagonalization)](https://github.com/Quantum-Many-Body/ExactDiagonalization.jl)**, **[CPT/VCA (cluster perturbation theory / variational cluster approach)](https://github.com/Quantum-Many-Body/QuantumClusterTheories.jl)**, **DMRG (density matrix renormalization group)**, etc. Generic interfaces are defined to provide a unified access to these algorithms with automatic project management.

## Installation

In Julia **v1.6+**, please type `]` in the REPL to use the package mode, then type this command:

```julia
pkg> add QuantumLattices
```

## Package Features

The mathematical foundations of our package is that the operators in a lattice Hamiltonian:
* **act on local Hilbert spaces**, and
* **form an algebra over the complex field**.

Based on this, the package has the following features:
* **Unitcell Description Framework**: the Hamiltonian can be constructed based on the unitcell of a lattice with the information of the local algebra acting on the local Hilbert space living on each point and the terms that couples different degrees of freedom on the same or different points. Such information can be input into the program as simple as describing the quantum system in a usual research paper.

* **Complete Symbolic Computation**: with only this package, symbolic computation between operators is realized while the coefficient of any operator remains numeric; by integrating it with [SymPy](https://github.com/JuliaPy/SymPy.jl), complete symbolic computation can be achieved and no modifications need be made on the methods in this package.

* **Generic Frontend of Many-Body Algorithms**: quantum many-body algorithms can be initialized in quite similar ways with only minor modifications needed. Moreover, automatic project management is realized, including that of result recording, data caching, parameter updating, information logging, dependency managing, etc.

## Supported Systems

Four common categories of quantum lattice systems in condensed matter physics are supported:
* **canonical fermionic systems**
* **canonical/hard-core bosonic systems**
* **SU(2) spin systems**
* **Phononic systems**

Furthermore, other systems can be supported easily by extending the generic protocols provided in this package.

## Supported Algorithms

Concrete algorithms could be considered as the "backend" of quantum lattice systems. They are developed in separate packages (still in progress):
* **[TBA](https://github.com/Quantum-Many-Body/TightBindingApproximation.jl)**: tight-binding approximation for fermionic/bosonic systems;
* **SCMF**: self-consistent mean field theory for fermionic systems;
* **[ED](https://github.com/Quantum-Many-Body/ExactDiagonalization.jl)**: exact diagonalization for fermionic/hard-core-bosonic/spin systems;
* **[CPT/VCA](https://github.com/Quantum-Many-Body/QuantumClusterTheories.jl)**: cluster perturbation theory and variational cluster approach for fermionic systems;
* **DMRG**: density matrix renormalization group for fermionic/hard-core-bosonic/spin systems;
* **[LSWT](https://github.com/Quantum-Many-Body/SpinWaveTheory.jl)**: linear spin wave theory for magnetically ordered local spin systems.

## Getting Started
* [Tutorials: unitcell description](@ref UnitcellDescriptionIntroduction)
* [Tutorials: advanced usage](@ref AdvancedUsageIntroduction)

## Note

Due to the fast development of this package, releases with different minor version numbers are **not** guaranteed to be compatible with previous ones **before** the release of v1.0.0. Comments are welcomed in the GitHub issues.

## Contact
waltergu1989@gmail.com

## Python counterpart
[HamiltonianPy](https://github.com/waltergu/HamiltonianPy): in fact, the authors of this Julia package worked on the python package at first and only turned to Julia later.