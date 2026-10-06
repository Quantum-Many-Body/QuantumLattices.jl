```@meta
CurrentModule = QuantumLattices.Frameworks
```

# Frameworks

*Module for high-level abstractions: models, algorithms, and project management.*

The `Frameworks` module provides the top-level abstractions that tie everything together. Its root is [`FrameworkElement`](@ref), the common supertype of everything the module manages: it carries the unified protocol of parameters ([`Parameters`](@ref)/[`update!`](@ref)), naming ([`str`](@ref)/`basename`/`pathof`), persistence ([`config`](@ref)/[`stamp`](@ref)/[`qlsave`](@ref)) and display, plus the optional `valtype` machinery from which `scalartype`/`eltype` are derived. Its three branches are [`LatticeModel`](@ref), the description of a quantum lattice system; [`Algorithm`](@ref), the execution context; and [`Assignment`](@ref), a single computation task assigned to an algorithm.

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
- **frontend** (a role, not a type): a [`LatticeModel`](@ref) subtype defined by an algorithm package, holding the model prepared in the form that algorithm needs.
- **[`Algorithm`](@ref)**: a frontend together with the algorithm-level parameters, the `map` and the timer. It is a transparent proxy of its frontend for pure query interfaces, which are mirrored onto it by [`@delegate`](@ref).
- **task** (a role, not a type): a plain struct specifying what to compute; [`dependencytypes`](@ref) declares the expected shape of its dependencies, and [`datatype`](@ref) infers the type of its result from the return type of [`run!`](@ref).
- **[`Data`](@ref)**: the result of a computation, with a tuple view used for output and plotting.
- **[`Assignment`](@ref)**: a named computation task paired with parameters, dependencies and a result slot; created by calling an [`Algorithm`](@ref) and executed by [`run!`](@ref).

**Utilities:**
- **[`Parameters`](@ref)**, **[`Boundary`](@ref)**: parameter and boundary condition management.
- **[`ParametricGenerator`](@ref)**, **[`StaticGenerator`](@ref)**: parameter sweep utilities.

**Configuration management:** [`config`](@ref), [`stamp`](@ref), [`qlsave`](@ref)/[`qlload`](@ref), [`qlcsave`](@ref)/[`qldclean`](@ref) provide automatic caching, persistence, and reproducibility.

See tutorial chapters 6, 7 and 8 for detailed narratives.

```@autodocs
Modules = [Frameworks]
Order = [:module, :constant, :type, :macro, :function]
```
