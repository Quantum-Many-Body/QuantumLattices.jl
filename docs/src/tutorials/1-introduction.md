```@meta
CurrentModule = QuantumLattices
```

# [1. Introduction](@id TutorialIntroduction)

A large class of problems in quantum many-body physics involves particles that interact on a lattice. The Hamiltonian of such a system, which encodes its equilibrium and dynamical properties, is built from the operators that form an algebra acting on the local Hilbert space at each lattice point. Writing down this Hamiltonian correctly and efficiently is the first step of any theoretical investigation, yet it is surprisingly easy to get wrong by hand when the system is beyond the simplest cases.

[QuantumLattices](https://github.com/Quantum-Many-Body/QuantumLattices.jl) automates this process. In the standard workflow, a quantum lattice system is completely defined based on its unitcell, which essentially requires three types of information:

1. The **spatial information**, such as the coordinates of the points contained in the unitcell of the lattice;
2. The **internal degrees of freedom**, such as the operator algebra acting on the local Hilbert space at each point;
3. The **couplings among different degrees of freedom**, such as the terms present in the Hamiltonian.

These three ingredients are combined into a [`LatticeModel`](@ref), a unifying abstraction that automatically generates the Hamiltonian and bundles it together with the model parameters into a single object. From there, the model can be passed directly to any supported quantum many-body algorithm.

Beyond this standard pipeline, which is covered in detail throughout this tutorial, [`LatticeModel`](@ref) can also wrap a Julia function or precomputed operators, making it possible to interface with external inputs such as [Wannier90](https://github.com/wannier-developers/wannier90) Hamiltonians. This flexibility, together with the symbolic algebra powered by [SymPy](https://github.com/JuliaPy/SymPy.jl), makes [QuantumLattices](https://github.com/Quantum-Many-Body/QuantumLattices.jl) a truly generic platform for quantum many-body computations.

## 1.1 Standard workflow: unitcell description framework

Let's start with a simple example, *"the single orbital electronic Hubbard model with only nearest neighbor hopping on a one dimensional lattice with only two sites"*. Here, *"one dimensional lattice with only two sites"* describes the spatial information, *"single orbital electronic"* defines the local Hilbert space and thus the local operator algebra, and *"Hubbard model with only nearest neighbor hopping"* expresses the terms present in the Hamiltonian. From this description, we can derive that the Hamiltonian of the system is

```math
H=tc^\dagger_{1\uparrow}c_{2\uparrow}+tc^\dagger_{2\uparrow}c_{1\uparrow}+tc^\dagger_{1\downarrow}c_{2\downarrow}+tc^\dagger_{2\downarrow}c_{1\downarrow}+Uc^\dagger_{1\uparrow}c_{1\uparrow}c^\dagger_{1\downarrow}c_{1\downarrow}+Uc^\dagger_{2\uparrow}c_{2\uparrow}c^\dagger_{2\downarrow}c_{2\downarrow}
```

where ``t`` is the hopping amplitude, ``U`` is the Hubbard interaction strength, and the electronic creation/annihilation operator ``c^\dagger_{i\sigma}/c_{i\sigma}`` carries a site index ``i`` (``i=1, 2``) and a spin index ``\sigma`` (``\sigma=\uparrow, \downarrow``). The **unitcell description framework** follows exactly this train of thought. The same system can be constructed with the following code:

```@example intro
using QuantumLattices
using SymPy: Sym, symbols

# Define the unitcell
lattice = Lattice([zero(Sym)], [one(Sym)])

# Define the internal degrees of freedom, i.e., the single-orbital spin-1/2 fermionic algebra
hilbert = Hilbert(site => Fock{:f}(1, 2) for site in eachindex(lattice))

# Define the terms
t = Hopping(:t, symbols("t", real=true), 1)
U = Hubbard(:U, symbols("U", real=true))

# Build the model and display the Hamiltonian
model = LatticeModel(lattice, hilbert, (t, U))
expand(model)
```

Note that on this page, the Hamiltonian is rendered in LaTeX because the documentation is generated with Documenter; in a plain REPL, `expand(model)` prints the same object as plain text. Here, in the subscript of the electronic annihilation/creation operator, an extra orbital index is also displayed. Let's walk through each component:

- **`Lattice([zero(Sym)], [one(Sym)])`** creates a 1D lattice that only contains 2 sites at positions 0 and 1. The coordinates are given as SymPy symbols only so that the whole construction stays symbolic; plain numbers such as `Lattice([0.0], [1.0])` work equally well, even together with symbolic coefficients. Note that no translation vectors are given here: a `Lattice` without translation vectors describes a finite cluster, for which the "unitcell" coincides with the whole system. Supplying the `vectors` keyword, as will be shown in [Chapter 2](@ref TutorialLattice), turns it into the unitcell of an infinite lattice.
- **`Hilbert(site => Fock{:f}(1, 2) for site in eachindex(lattice))`** assigns a single-orbital (1), spin-1/2 (2) fermionic (`:f`) Fock space to each site.
- **`Hopping(:t, symbols("t", real=true), 1)`** defines a nearest-neighbor hopping term. The symbol `:t` serves as its identifier, and the coefficient is the [SymPy](https://github.com/JuliaPy/SymPy.jl) symbol `t`.
- **`Hubbard(:U, symbols("U", real=true))`** defines an on-site Hubbard interaction with symbolic amplitude `U`.
- **`LatticeModel(lattice, hilbert, (t, U))`** combines the three ingredients. The terms themselves determine which bonds are needed; the constructor handles this automatically.
- **`expand(model)`** generates the resulting Hamiltonian operators.

Detailed explanations will be covered in the following chapters.

## 1.2 How This Tutorial Is Organized

The tutorial follows the natural structure of a quantum lattice system. Each chapter builds on the previous ones:

| Chapter | Content |
|---------|---------|
| 2 | Spatial structure: lattices, points, bonds, and reciprocal space |
| 3 | Internal degrees of freedom: algebra of Fock/Spin/Phonon systems, and index ordering |
| 4 | Operator algebra: Operator, Operators, LaTeX output, and linear transformations |
| 5 | Couplings: Terms, coupling patterns, and automatic operator expansion |
| 6 | LatticeModel in depth: representations, parameters, and worked examples |
| 7 | The algorithm interface, user guide: connecting models to solvers, project management |
| 8 | The algorithm interface, developer guide: writing an algorithm package |

By the end, you will be able to define an arbitrary quantum lattice model as a [`LatticeModel`](@ref) and connect it to numerical algorithms.

## 1.3 Comparison with Existing Tools

Several excellent packages exist for working with quantum lattice systems. Here we highlight the distinctive features of [QuantumLattices](https://github.com/Quantum-Many-Body/QuantumLattices.jl):

- [Python/QuSpin](https://quspin.github.io/QuSpin/): QuSpin constructs Hamiltonians from operator strings and site-wise operators, which is expressive but requires manual enumeration of all terms. [QuantumLattices](https://github.com/Quantum-Many-Body/QuantumLattices.jl) specifies terms by bond kind and coupling patterns, and is fundamentally symbolic rather than matrix-based.
- [Julia/ITensor](https://github.com/ITensor/ITensors.jl): ITensor focuses on tensor-network representations; it provides building blocks for MPS/MPO algorithms. [QuantumLattices](https://github.com/Quantum-Many-Body/QuantumLattices.jl) operates at a higher level: you define the model, obtain a [`LatticeModel`](@ref), and that model can be passed to ITensor-based solvers through the algorithm interface.
- [Python/HamiltonianPy](https://github.com/waltergu/HamiltonianPy): This is the Python predecessor of [QuantumLattices](https://github.com/Quantum-Many-Body/QuantumLattices.jl). The Julia version benefits from Julia's type system and multiple dispatch, resulting in a cleaner API and competitive performance.
