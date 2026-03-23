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

Welcome to **[QuantumLattices](https://github.com/Quantum-Many-Body/QuantumLattices.jl)**. Here we provide a general framework to construct the **operator-formed Hamiltonian of any quantum lattice system**, with the inputs as simple as its description by the natural language. This operator-formed Hamiltonian supports **complete symbolic computations** when combined with **[SymPy](https://github.com/JuliaPy/SymPy.jl)**, and can serve as a convenient **frontend of quantum many-body algorithms**, such as **[TBA (tight-bind approximation)](https://github.com/Quantum-Many-Body/TightBindingApproximation.jl)**, **[LSWT (linear spin wave theory)](https://github.com/Quantum-Many-Body/SpinWaveTheory.jl)**, **[SCMF (self-consistent mean field theory)](https://github.com/Quantum-Many-Body/MeanFieldTheory.jl)**, **[ED (exact diagonalization)](https://github.com/Quantum-Many-Body/ExactDiagonalization.jl)**, **[CPT/VCA (cluster perturbation theory / variational cluster approach)](https://github.com/Quantum-Many-Body/QuantumClusterTheories.jl)**, **[DMRG (density matrix renormalization group)](https://github.com/ZongYongyue/DynamicalCorrelators.jl)**, **[RPA (random phase approximation)](https://github.com/Quantum-Many-Body/RandomPhaseApproximation.jl)**, etc. Generic interfaces are defined to provide a unified access to these algorithms with automatic project management.

## Installation

In Julia **v1.8+**, please type `]` in the REPL to use the package mode, then type this command:

```julia
pkg> add QuantumLattices
```

## Quick Start

Build your first quantum lattice system in just a few lines:

```julia
using QuantumLattices
using SymPy: Sym, symbols

# 1. Define the lattice (1D chain with 2 sites)
lattice = Lattice([zero(Sym)], [one(Sym)])

# 2. Define the internal degrees of freedom (spin-1/2 fermions)
hilbert = Hilbert(site => Fock{:f}(1, 2) for site in eachindex(lattice))

# 3. Define the Hamiltonian terms
t = Hopping(:t, symbols("t", real=true), 1) # nearest-neighbor hopping
U = Hubbard(:U, symbols("U", real=true))    # Hubbard interaction

# 4. Generate the Hamiltonian operators
operators = expand(OperatorGenerator(bonds(lattice, 1), hilbert, (t, U)))
```

This generates all the operators in the Hamiltonian:
```math
H = t c^\dagger_{1\uparrow} c_{2\uparrow} + t c^\dagger_{2\uparrow} c_{1\uparrow} + t c^\dagger_{1\downarrow} c_{2\downarrow} + t c^\dagger_{2\downarrow} c_{1\downarrow} + U n_{1\uparrow} n_{1\downarrow} + U n_{2\uparrow} n_{2\downarrow}
```

## Package Features

The mathematical foundations of our package are that the operators in a lattice Hamiltonian:
* **act on local Hilbert spaces**, and
* **form an algebra over the complex field**.

Based on this, the package provides the following features:
* **Unitcell Description Framework**: The Hamiltonian can be constructed based on the unitcell of a lattice, using information about the local algebra acting on the local Hilbert space at each point and the terms that couple different degrees of freedom at the same or different points. This information can be provided to the program in the same way as describing the quantum system in a research paper.

* **Complete Symbolic Computation**: With this package alone, symbolic computation between operators is supported while keeping the coefficients of any operator numeric. By integrating it with [SymPy](https://github.com/JuliaPy/SymPy.jl), complete symbolic computation can be achieved without requiring any modifications to the methods in this package.

* **Generic Frontend of Many-Body Algorithms**: Using the operator-formed Hamiltonian as a foundation, quantum many-body algorithms can be initialized in a consistent manner with minimal modifications. Moreover, automatic project management is provided, including result recording, data caching, parameter updating, information logging, dependency management, and more.

## Supported Systems

Four common categories of quantum lattice systems in condensed matter physics are supported:
* **Canonical complex fermionic systems**
* **Canonical complex and hard-core bosonic systems**
* **SU(2) spin systems**
* **Phononic systems**

Furthermore, other systems can be easily supported by extending the generic protocols provided in this package.

## Supported Algorithms

Concrete algorithms can be considered as the "backend" of quantum lattice systems. They are developed in separate packages:
* **[TBA](https://github.com/Quantum-Many-Body/TightBindingApproximation.jl)**: Tight-binding approximation for complex-fermionic/complex-bosonic/phononic systems.
* **[LSWT](https://github.com/Quantum-Many-Body/SpinWaveTheory.jl)**: Linear spin wave theory for magnetically ordered local-spin systems.
* **[SCMF](https://github.com/Quantum-Many-Body/MeanFieldTheory.jl)**: Self-consistent mean field theory for complex fermionic systems.
* **[ED](https://github.com/Quantum-Many-Body/ExactDiagonalization.jl)**: Exact diagonalization for complex-fermionic/hard-core-bosonic/local-spin systems.
* **[CPT/VCA](https://github.com/Quantum-Many-Body/QuantumClusterTheories.jl)**: Cluster perturbation theory and variational cluster approach for complex fermionic and local spin systems.
* **[DMRG](https://github.com/ZongYongyue/DynamicalCorrelators.jl)**: Density matrix renormalization group for complex-fermionic/hard-core-bosonic/local-spin systems based on [TensorKit](https://github.com/Jutho/TensorKit.jl) and [MPSKit](https://github.com/QuantumKitHub/MPSKit.jl).
* **[RPA](https://github.com/Quantum-Many-Body/RandomPhaseApproximation.jl)**: Random phase approximation for complex fermionic systems.

## Getting Started

* [Tutorial: Unitcell Description](@ref UnitcellDescriptionIntroduction)
* [Tutorial: Advanced Topics](@ref AdvancedTopicsIntroduction)

## Note

Due to the rapid development of this package, releases with different minor version numbers are **not** guaranteed to be compatible with previous ones **before** the release of v1.0.0. Comments are welcome in the GitHub issues.

## Contact

* Email: waltergu1989@gmail.com

## Python Counterpart

[HamiltonianPy](https://github.com/waltergu/HamiltonianPy): The authors of this Julia package initially worked on a Python package before transitioning to Julia.