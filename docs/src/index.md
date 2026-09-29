# QuantumLattices.jl

[![CI](https://github.com/Quantum-Many-Body/QuantumLattices.jl/actions/workflows/CI.yml/badge.svg)](https://github.com/Quantum-Many-Body/QuantumLattices.jl/actions/workflows/CI.yml)
[![codecov](https://codecov.io/gh/Quantum-Many-Body/QuantumLattices.jl/branch/master/graph/badge.svg)](https://codecov.io/gh/Quantum-Many-Body/QuantumLattices.jl)
[![](https://img.shields.io/badge/docs-latest-blue.svg)](https://quantum-many-body.github.io/QuantumLattices.jl/latest/)
[![](https://img.shields.io/badge/docs-stable-blue.svg)](https://quantum-many-body.github.io/QuantumLattices.jl/stable/)
[![LICENSE](https://img.shields.io/badge/License-Apache%202.0-blue.svg)](https://opensource.org/licenses/Apache-2.0)
[![LICENSE](https://img.shields.io/badge/license-Anti%20996-blue.svg)](https://github.com/996icu/996.ICU/blob/master/LICENSE)
[![Code Style: Blue](https://img.shields.io/badge/code%20style-blue-4495d1.svg)](https://github.com/invenia/BlueStyle)
[![ColPrac: Contributor's Guide on Collaborative Practices for Community Packages](https://img.shields.io/badge/ColPrac-Contributor's%20Guide-blueviolet)](https://github.com/SciML/ColPrac)

*Julia package for the construction of quantum lattice systems.*

Welcome to [QuantumLattices](https://github.com/Quantum-Many-Body/QuantumLattices.jl). This package provides a general framework to construct arbitrary quantum lattice systems and derive their Hamiltonians from inputs as simple as a natural-language description. The standard definition of a lattice model consists of three components: the lattice, the internal degrees of freedom, and the coupling terms. Then the operator-based Hamiltonian can be generated automatically and fed directly into any supported quantum many-body algorithm, including [TBA (tight-binding approximation)](https://github.com/Quantum-Many-Body/TightBindingApproximation.jl), [ED (exact diagonalization)](https://github.com/Quantum-Many-Body/ExactDiagonalization.jl), [CPT/VCA (cluster perturbation theory / variational cluster approach)](https://github.com/Quantum-Many-Body/QuantumClusterTheories.jl), [DMRG (density matrix renormalization group)](https://github.com/ZongYongyue/DynamicalCorrelators.jl), [LSWT (linear spin wave theory)](https://github.com/Quantum-Many-Body/SpinWaveTheory.jl), [SCMF (self-consistent mean field theory)](https://github.com/Quantum-Many-Body/MeanFieldTheory.jl), [RPA (random phase approximation)](https://github.com/Quantum-Many-Body/RandomPhaseApproximation.jl), and more. When combined with [SymPy](https://github.com/JuliaPy/SymPy.jl), the entire Hamiltonian, including all the coefficients and algebraic operations, can be kept fully symbolic and manipulated analytically, with the numerical evaluation postponed until explicitly requested. None of this requires any modification to the package methods. Beyond this standard pipeline, the Hamiltonian can also be supplied directly as a Julia function or prebuilt operators, enabling seamless interfacing with external inputs such as [Wannier90](https://github.com/wannier-developers/wannier90) Hamiltonians. Automatic project management, including result recording, data caching, parameter updating, and dependency tracking, is provided uniformly regardless of which representation is used.

## Installation

In Julia **v1.10+**, type `]` in the REPL to enter the package mode, then:

```julia
pkg> add QuantumLattices
```

## Quick Start

Build your first quantum lattice system in just a few lines:

```@example quickstart
using QuantumLattices
using SymPy: Sym, symbols    # install with: pkg> add SymPy

# 1. Define the lattice (1D chain with 2 sites)
lattice = Lattice([zero(Sym)], [one(Sym)])

# 2. Define the internal degrees of freedom (spin-1/2 fermions)
hilbert = Hilbert(site => Fock{:f}(1, 2) for site in eachindex(lattice))

# 3. Define the Hamiltonian terms
t = Hopping(:t, symbols("t", real=true), 1)     # nearest-neighbor hopping
U = Hubbard(:U, symbols("U", real=true))        # Hubbard interaction

# 4. Build the model and display the Hamiltonian
model = LatticeModel(lattice, hilbert, (t, U))
expand(model)
```

## Package Features

The package is built on three integrated features:

1. **Unitcell Description Framework**: In the standard workflow, a quantum lattice system is fully specified by its lattice, its internal degrees of freedom, and its coupling terms: the same three ingredients shown in the [Quick Start](#Quick-Start) above. These are bundled into a [`LatticeModel`](@ref), which automatically generates the complete Hamiltonian, much as one would write it in a research paper.

2. **Symbolic Operator Algebra**: Operators form an algebra over the complex field. Sums, products, scalar multiplication, and (anti-)commutation relations are all supported natively. Combined with [SymPy](https://github.com/JuliaPy/SymPy.jl), the entire Hamiltonian, including all the coefficients and algebraic operations, remains symbolic until numeric evaluation is explicitly requested. No modifications to the package methods are required.

3. **Algorithm Interface**: [`LatticeModel`](@ref) provides a uniform interface to all supported quantum many-body algorithms. Automatic project management, including result recording, data caching, parameter updating, and dependency tracking, is built in directly. Beyond the standard unitcell-based pipeline, [`LatticeModel`](@ref) can also wrap a Julia function or prebuilt operators, enabling interfaces with external inputs such as [Wannier90](https://github.com/wannier-developers/wannier90) Hamiltonians.

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
* **[ED](https://github.com/Quantum-Many-Body/ExactDiagonalization.jl)**: Exact diagonalization for complex-fermionic/hard-core-bosonic/local-spin systems.
* **[CPT/VCA](https://github.com/Quantum-Many-Body/QuantumClusterTheories.jl)**: Cluster perturbation theory and variational cluster approach for complex fermionic and local spin systems.
* **[DMRG](https://github.com/ZongYongyue/DynamicalCorrelators.jl)**: Density matrix renormalization group for complex-fermionic/hard-core-bosonic/local-spin systems based on [TensorKit](https://github.com/Jutho/TensorKit.jl) and [MPSKit](https://github.com/QuantumKitHub/MPSKit.jl).
* **[LSWT](https://github.com/Quantum-Many-Body/SpinWaveTheory.jl)**: Linear spin wave theory for magnetically ordered local-spin systems.
* **[SCMF](https://github.com/Quantum-Many-Body/MeanFieldTheory.jl)**: Self-consistent mean field theory for complex fermionic systems.
* **[RPA](https://github.com/Quantum-Many-Body/RandomPhaseApproximation.jl)**: Random phase approximation for complex fermionic systems.

## Getting Started

* **Tutorials**: A step-by-step guide from first principles to building models and connecting to algorithms. Start with the [Introduction](@ref TutorialIntroduction).
* **Advanced Topics**: Specialized guides on advanced topics, including [hybrid systems](@ref HybridSystems), [boundary conditions](@ref BoundaryConditions), and [linear transformations](@ref LinearTransformations).
* **Manual**: Module-by-module reference with concept introductions, useful for looking up specific types and functions.

## Note

Due to the rapid development of this package, releases with different minor version numbers are **not** guaranteed to be compatible with previous ones **before** the release of v1.0.0. Comments are welcome in the GitHub issues.

## Contact

* Email: waltergu1989@gmail.com

## Python Counterpart

[HamiltonianPy](https://github.com/waltergu/HamiltonianPy): The authors of this Julia package initially worked on a Python package before transitioning to Julia.
