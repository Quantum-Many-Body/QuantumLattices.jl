```@meta
CurrentModule = QuantumLattices.Frameworks
```

# Frameworks

*Module for high-level abstractions: models, algorithms, and project management.*

The `Frameworks` module provides the top-level abstractions that tie everything together. Its centerpiece is [`LatticeModel`](@ref), the unifying container that wraps a lattice, Hilbert space, Hamiltonian, and parameters into a single object and serves as the universal input to quantum many-body algorithms.

**Model layer:**
- **[`LatticeModel`](@ref)**: the complete description of a quantum lattice system. Acts as a constructor factory dispatching to [`Formula`](@ref), [`StaticGenerator`](@ref), or [`OperatorGenerator`](@ref).
- **[`Formula`](@ref)**: function-based Hamiltonian representation.
- **[`Generator`](@ref)** (abstract): operator-based Hamiltonian representation.
- **[`StaticGenerator`](@ref)**: wraps a fixed [`OperatorSet`](@ref).
- **[`ParametricGenerator`](@ref)** (abstract): supports parameter updates.
- **[`CategorizedGenerator`](@ref)**: groups operators into constant, alterable, and boundary categories.
- **[`OperatorGenerator`](@ref)**: generates operators from bonds, Hilbert space, and terms.
- **[`Embedding`](@ref)**: replicates unitcell operators onto larger lattices.

**Algorithm interface:**
- **[`Frontend`](@ref)**, **[`Algorithm`](@ref)**, **[`Action`](@ref)**, **[`Data`](@ref)**, **[`Assignment`](@ref)**: the standard protocol for connecting models to solvers.

**Utilities:**
- **[`Parameters`](@ref)**, **[`Boundary`](@ref)**: parameter and boundary condition management.
- **[`ParametricGenerator`](@ref)**, **[`StaticGenerator`](@ref)**: parameter sweep utilities.

**Configuration management:** [`config`](@ref), [`stamp`](@ref), [`qlsave`](@ref)/[`qlload`](@ref), [`qlcsave`](@ref)/[`qldclean`](@ref) provide automatic caching, persistence, and reproducibility.

See tutorial chapters 6 and 7 for detailed narratives.

```@autodocs
Modules = [Frameworks]
Order = [:module, :constant, :type, :macro, :function]
```
