```@meta
CurrentModule = QuantumLattices.Spatials
```

# Spatials

*Module for spatial structures of quantum lattice systems.*

The `Spatials` module provides the fundamental types for describing the geometry of a quantum lattice system: lattices, points, bonds, and reciprocal space. It implements the mathematical concepts introduced in the Tutorials.

**Key types:** [`Lattice`](@ref), [`Point`](@ref), [`Bond`](@ref), [`Neighbors`](@ref), [`BrillouinZone`](@ref), [`ReciprocalZone`](@ref), [`ReciprocalPath`](@ref).

**Key functions:** [`bonds`](@ref), [`isparallel`](@ref), [`reciprocals`](@ref), [`tile`](@ref), [`minimumlengths`](@ref).

```@autodocs
Modules = [Spatials]
Order = [:module, :constant, :type, :macro, :function]
```
