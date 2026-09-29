```@meta
CurrentModule = QuantumLattices.DegreesOfFreedom
```

# Degrees of Freedom

*Module for internal degrees of freedom of quantum lattice systems.*

The `DegreesOfFreedom` module provides the types for specifying the internal structure of a quantum lattice system: the local algebras acting on each site, their generators, and the hierarchical index system that labels them. These concepts are explained in tutorial chapters 4 and 6.

**Key types:** [`Hilbert`](@ref), [`Internal`](@ref), `Fock`, `Spin`, `Phonon`, [`Index`](@ref), [`CoordinatedIndex`](@ref), [`Coupling`](@ref), [`Term`](@ref).

```@autodocs
Modules = [DegreesOfFreedom]
Order = [:module, :constant, :type, :macro, :function]
```
