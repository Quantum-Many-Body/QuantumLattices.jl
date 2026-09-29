```@meta
CurrentModule = QuantumLattices.QuantumSystems
```

# Quantum Systems

*Module for pre-built quantum lattice systems and term libraries.*

The `QuantumSystems` module provides concrete implementations of the common categories of quantum lattice systems: fermionic, bosonic, spin, and phononic. It also defines the standard term types (Hopping, Hubbard, Heisenberg, etc.) and their coupling patterns.

**Key types for algebras:** [`Fock`](@ref), [`FockIndex`](@ref), [`Spin`](@ref), [`SpinIndex`](@ref), [`Phonon`](@ref), [`PhononIndex`](@ref).

**Key types for terms:** [`Hopping`](@ref), [`Hubbard`](@ref), [`Heisenberg`](@ref), [`Ising`](@ref), [`Kitaev`](@ref), [`Zeeman`](@ref), [`Pairing`](@ref), [`Coulomb`](@ref), [`Onsite`](@ref), and more.

**Convenience functions:** `𝕔`, `𝕔⁺`, `𝕒`, `𝕒⁺`, `𝕕`, `𝕕⁺` for fermionic/bosonic indices; `𝕊` for spin indices; `σˣ`, `σʸ`, `σᶻ` for Pauli matrices.

```@autodocs
Modules = [QuantumSystems]
Order = [:module, :constant, :type, :macro, :function]
```
