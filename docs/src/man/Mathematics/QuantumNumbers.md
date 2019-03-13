```@meta
CurrentModule=QuantumLattices.Mathematics.QuantumNumbers
```

# Quantum numbers

Qunatum numbers can be considered as the conserved labels for the bases of a Hilbert space when a quantum system hosts some symmetries. Here we only implement Abelian quantum numbers because non-Abelian ones are far more complicated yet much less used. In practice, quantum numbers can be integers or half integers, therefore, we use real numbers to denote them in this module for simplicity. Independent quantum numbers, such as the particle number and the spin z-component, can coexist at the same time. We use type [`AbelianNumber`](@ref) to represent the complete set of independent ones for a single basis of a Hilbert space, and type [`AbelianNumbers`](@ref) to represent the whole quantum numbers for the total bases.

## AbelianNumber

The abstract type for the complete set of independent quantum numbers for a single basis.

Main features include:
* function `fieldnames`: get the names of the quantum numbers
* function `periods`: get the periods of the quantum numbers
* arithmetic operations: `+`, `-`, `*`, `^`, `⊕`, `⊗`
* hashable: concrete instances can be used as keys for a dict or a set
* iterable: concrete instances are iterable over their values
* comparable: two concrete instances can be compared

In particular, `AbelianNumber <: HomoNamedVector{Float}`, all features supported by `HomoNamedVector` are also available for `HomoNamedVector`. See also [HomoNamedVector](@ref).

For convenience, **4** kinds of quantum numbers are predefined in this module, i.e.
* [`SQN`](@ref): for spin z-component reserved systems
* [`PQN`](@ref): for particle number reserved systems
* [`SPQN`](@ref): for both particle number and spin-z component reserved systems
* [`Z2QN`](@ref): for systems with a ``Z_2`` conservation quantum number

Users who want to define their own ``Z_N``-like quantum numbers must handle the periodicities in the construction function, otherwise, wrong results will be get when arithmetic operations, such as `+` or `-`, are involved. It is highly recommended to use the macro [`@quantumnumber`](@ref) to define your own concrete `AbelianNumber`s.

## AbelianNumbers

The whole quantum numbers for the total bases.

By design, a [`AbelianNumbers{QN}`](@ref) has one type parameter:
* `QN<:AbelianNumber`: the type of the quantum numbers contained in it
And 3 attributes:
* `form::Char`: Its form, whose value must be one of the followings
  - `'G'`: the general form, which has no restriction for its `contents`
  - `'U'`: the unitary form, which requires no duplicates in its `contents`
  - `'C'`: the canonical form, which requires both no duplicates and accending-order in its `contents`
  Usually, G-formed and U-formed `AbelianNumbers`es can be transformed to the corresponding C-formed ones by the [`sort`](@ref) function.
* `contents::Vector{QN}`: The quantum numbers contained in it. To achieve high efficiency, it is required to be an homogenous array of a certain kind of concrete `AbelianNumber`.
* `indptr::Vector{Int}`: The indptr of the quantum numbers contained in it, which is similar to the `colptr` attribute of a [CSC sparse matrix](https://docs.julialang.org/en/v1/stdlib/SparseArrays/#man-csc-1) and records the compression info of its `contents`.

Main features include:
* function `eltype`: get the concrete type of the quantum numbers it contains
* index access: get the contents directly by the `getindex` function
* arithmetic operations: `+`, `-`, `*`, `^`, `⊗`, `⊕`
* iterable: various iteration supports, including functions such as `iterate`, `keys`, `values` and `pairs`
For a complete summation of its features, please refer to the [manual](@ref qnmanual).

For convenience, **5** functions are predefined to generate the `AbelianNumbers` of common physical systems, i.e.
* [`SQNS`](@ref): a signle spin
* [`PQNS`](@ref): a single-particle state with at most `N` identical particles
* [`SzPQNS`](@ref): a single-paritcle state with at most one particle whose spin-z component is `Sz`
* [`SPQNS`](@ref): a single site with internal degrees of freedom that can be ascribed to a spin
* [`Z2QNS`](@ref): any ``Z_2`` Hilbert space

## [Manual](@id qnmanual)

```@autodocs
Modules=[QuantumNumbers]
Order=  [:module,:constant,:type,:macro,:function]
```
