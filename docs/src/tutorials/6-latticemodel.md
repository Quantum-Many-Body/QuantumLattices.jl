```@meta
CurrentModule = QuantumLattices
DocTestSetup = quote
    using QuantumLattices
end
```

```@setup latticemodel
using QuantumLattices
using SymPy: Sym, symbols
```

# [6. LatticeModel: the Unifying Abstraction](@id TutorialLatticeModel)

We have now covered the entire pipeline: lattices ([Chapter 2](@ref TutorialLattice)), internal degrees of freedom ([Chapter 3](@ref TutorialDOF)), operators ([Chapter 4](@ref TutorialOperators)), and coupling terms ([Chapter 5](@ref TutorialCouplings)). The final piece is the container that ties everything together into a single object: [`LatticeModel`](@ref).

## 6.1 What a LatticeModel Is

[`LatticeModel`](@ref) is an abstract supertype that aims at generating the Hamiltonian of a quantum lattice system. It is the object that algorithms accept as their input. Its concrete representations differ in what they store. The standard workflow of this tutorial uses the first of them:

* [`OperatorGenerator`](@ref): the standard, operator-based representation. It stores the bonds, the Hilbert space, the coupling terms and the parameters, and generates the operator Hamiltonian on demand.

The standard workflow is not always convenient, however: a Hamiltonian may come from an external code, such as a [Wannier90](https://github.com/wannier-developers/wannier90)-exported tight-binding model, which can be efficiently converted into operator form, but would be inefficient to reconstruct from a lattice, a Hilbert space and coupling terms; or it may be known analytically, so that regenerating it term by term would be wasteful. For these cases, two further representations exist:

* [`StaticGenerator`](@ref): wraps a fixed set of prebuilt operators, e.g. a Hamiltonian exported by an external code.
* [`Formula`](@ref): wraps a Julia function of the parameters, e.g. an analytically known Hamiltonian.

Despite the differences in what they store, all three representations support a common set of generic interfaces, which are grouped below by their purposes. These groups make a model interchangeable between algorithms, so that one model can be passed to any algorithm and one algorithm can accept any model:

* **Hamiltonian access**: [`expand`](@ref) expands the Hamiltonian into operators, and the collection interface (`length`, `iterate`, `isempty`) inspects the result (a [`Formula`](@ref), whose value may be of any type rather than operators, supports neither).
* **Parameter management**: [`Parameters`](@ref) returns the model parameters, and [`update!`](@ref) changes the tunable ones.
* **Value types**: [`valtype`](@ref), [`eltype`](@ref) and [`scalartype`](@ref) report the types of the Hamiltonian, of its elements, and of its coefficients, the static information an algorithm needs before it runs.

This uniformity of the interface is precisely the benefit of unifying the three representations under [`LatticeModel`](@ref). Which representation is created is decided by the arguments passed to the [`LatticeModel`](@ref) constructor, as the following sections show.

## 6.2 Standard Construction

### 6.2.1 From Lattice, Hilbert Space, and Terms

The standard way to build a model is from a lattice, a Hilbert space, and the coupling terms:

```julia
LatticeModel(lattice, hilbert, terms) -> OperatorGenerator
```

The bonds on which the terms act are inferred automatically from the terms themselves, via [`nneighbor`](@ref): each [`Term`](@ref) ([Chapter 5](@ref TutorialCouplings)) declares a bond kind, and the constructor collects the bonds up to the highest order required. Let us build the two-site Hubbard model symbolically, and obtain the expanded Hamiltonian by [`expand`](@ref):

```@example latticemodel
lattice = Lattice([zero(Sym)], [one(Sym)])
hilbert = Hilbert(site => Fock{:f}(1, 2) for site in eachindex(lattice))
t = Hopping(:t, symbols("t", real=true), 1)
U = Hubbard(:U, symbols("U", real=true))
model = LatticeModel(lattice, hilbert, (t, U))
expand(model)
```

Being an [`OperatorGenerator`](@ref), the model supports the Hamiltonian-access interface of Section 6.1: the collection interface (`length`, `iterate`, `isempty`) inspects the expanded operators, e.g. `length(model)` counts them:

```@example latticemodel
length(model) |> println
```

### 6.2.2 Restricting the Bonds

For fine-grained control over which bonds are included, the bonds can be passed explicitly:

```julia
LatticeModel(bonds, hilbert, terms) -> OperatorGenerator
```

Instead of the whole lattice, one may pass only a selected set of bonds, thereby obtaining a subset of the lattice Hamiltonian. For instance, on a two-site chain, passing only the 1st-neighbor bonds yields a model that contains hopping operators alone, and passing only the one-point bonds yields one that contains Hubbard operators alone:

```@example latticemodel
lattice = Lattice([zero(Sym)], [one(Sym)])
hilbert = Hilbert(site => Fock{:f}(1, 2) for site in eachindex(lattice))
t = Hopping(:t, symbols("t", real=true), 1)
U = Hubbard(:U, symbols("U", real=true))

bonds_hop = bonds(lattice, Neighbors(1=>1))
model_hop = LatticeModel(bonds_hop, hilbert, (t, U))
expand(model_hop)
```

```@example latticemodel
bonds_hub = bonds(lattice, Neighbors(0=>0))
model_hub = LatticeModel(bonds_hub, hilbert, (t, U))
expand(model_hub)
```

The first model expands into 4 hopping operators and the second into 2 Hubbard operators, although both models accept the hopping and Hubbard terms simultaneously. Restricting the bonds this way is sometimes useful in quantum many-body algorithms.

### 6.2.3 Parameter Management and Value Types

This section demonstrates two more interface groups of Section 6.1, namely parameter management and value types, with an [`OperatorGenerator`](@ref).

#### Parameter management

[`Parameters`](@ref) gives the parameter name-value pairs of a model, and [`update!`](@ref) changes them in place while keeping the structure of the expansion intact:

```@example latticemodel
lattice = Lattice([zero(Sym)], [one(Sym)])
hilbert = Hilbert(site => Fock{:f}(1, 2) for site in eachindex(lattice))
t = Hopping(:t, symbols("t", real=true), 1)
U = Hubbard(:U, symbols("U", real=true))
model = LatticeModel(lattice, hilbert, (t, U))
Parameters(model) |> println
```

```@example latticemodel
update!(model; t = 2, U = 4)
expand(model)
```

!!! note
    The [`update!`](@ref) function acts in place, so that the coefficients of the terms stored in an [`OperatorGenerator`](@ref) are updated. Since a model stores the very terms you pass it, which are not copied at construction by default, this also updates those terms. In the example above, this means that the `t` and `U` defined outside `model` are changed as well:

    ```@example latticemodel
    (t.value, U.value)
    ```

    To avoid this, you can explicitly `deepcopy` the terms before passing them to the construction function:

    ```@example latticemodel
    t = Hopping(:t, symbols("t", real=true), 1)
    U = Hubbard(:U, symbols("U", real=true))
    model = LatticeModel(lattice, hilbert, deepcopy((t, U)))
    update!(model; t = 2, U = 4)
    (t.value, U.value)
    ```

    In general, a model never copies its terms; if terms are reused or modified afterwards, pass copies.

#### Value types

[`valtype`](@ref), [`eltype`](@ref) and [`scalartype`](@ref) report the type of the Hamiltonian represented by the model, the type of its elements, and the type of the coefficients:

```@example latticemodel
lattice = Lattice([0.0], [1.0])
hilbert = Hilbert(site => Fock{:f}(1, 2) for site in eachindex(lattice))
t = Hopping(:t, -1.0, 1)
U = Hubbard(:U, 8.0)
model = LatticeModel(lattice, hilbert, (t, U))

(valtype(model) <: Operators, eltype(model) <: Operator, scalartype(model))
```

## 6.3 Pre-Built Construction

Another way to build a model is from a fixed set of operators:

```julia
LatticeModel(operators::Operators) -> StaticGenerator
```

This path is meant for integrating external Hamiltonians: any Hamiltonian produced by another code, for instance a tight-binding Hamiltonian exported by [Wannier90](https://github.com/wannier-developers/wannier90), can be converted into the [`Operators`](@ref)-like objects of this package, wrapped into a [`LatticeModel`](@ref), and then passed to any algorithm written against the common interface of Section 6.1. Note that neither a lattice nor a Hilbert space is needed here: the model stores the operators as given. To see the actual integration with Wannier90 and the related technical details, see the [Interface With Wannier90](https://quantum-many-body.github.io/TightBindingApproximation.jl/dev/wannier90/) chapter of the [TightBindingApproximation.jl](https://github.com/Quantum-Many-Body/TightBindingApproximation.jl) package and its manual.

As a simple illustration, let us write down by hand the operators of a two-site hopping Hamiltonian, which contains the hopping in the two directions:

```@example latticemodel
t = symbols("t", real=true)
operators = Operators(
    Operator(t, 𝕔⁺(2, 1, 0, [one(Sym)/2], [zero(Sym)]), 𝕔(1, 1, 0, [zero(Sym)], [zero(Sym)])),
    Operator(t, 𝕔⁺(1, 1, 0, [zero(Sym)], [zero(Sym)]), 𝕔(2, 1, 0, [one(Sym)/2], [zero(Sym)])),
)
model = LatticeModel(operators)
expand(model)
```

The model is a [`StaticGenerator`](@ref) that shares the common interface of Section 6.1. `length(model)` gives the number of operators it stores, and `Parameters(model)` its parameters:

```@example latticemodel
length(model) |> println
Parameters(model) |> println
```

It is noted that a [`StaticGenerator`](@ref) carries no parameters; this is why it is called static. This holds even when the stored operators contain symbolic coefficients, as in the example above: [`Parameters`](@ref) returns an empty named tuple, and the symbols can never be substituted numerically through the model. If tunable parameters are needed, an [`OperatorGenerator`](@ref) or a [`Formula`](@ref) should be used instead. Correspondingly, [`update!`](@ref) on a [`StaticGenerator`](@ref) is a no-op:

```@example latticemodel
update!(model; t=5)
expand(model)
```

## 6.4 Function-Based Construction

Yet another way is to give the Hamiltonian directly by a type-stable Julia function of the parameters:

```julia
LatticeModel(expression::Function, parameters::Parameters) -> Formula
```

The function receives the parameter values positionally, in the order in which they appear in the parameter named tuple. Here, [`Parameters`](@ref) is an alias of `NamedTuple` whose values are restricted to numbers (SymPy symbols included, since `Sym` is a subtype of `Number`), so a plain named-tuple literal such as `(t=t, μ=μ, Δ=Δ)` below is already a valid instance. For example, the spinless p+ip-wave BdG Hamiltonian on the square lattice, which depends on the parameters ``t``, ``\mu`` and ``\Delta``, can be given directly by an analytic expression:

```@example latticemodel
ham(t, μ, Δ, k=[0, 0]) = [2t*cos(k[1]) + 2t*cos(k[2]) + μ  2im*Δ*sin(k[1]) + 2Δ*sin(k[2]);
                          -2im*Δ*sin(k[1]) + 2Δ*sin(k[2])  -2t*cos(k[1]) - 2t*cos(k[2]) - μ]
t, μ, Δ = symbols("t μ Δ", real=true)
model = LatticeModel(ham, (t=t, μ=μ, Δ=Δ))
nothing # hide
```

Besides the parameter values, the expression may take extra positional arguments, provided that they come after the parameter values, are given default values, and keep the expression type-stable. Such a default is required, because the constructor inspects the expression by applying it to the parameter values alone. In the example above, the momentum `k`, which defaults to the Γ point, is such an argument. The model then can be evaluated by supplying these extra arguments:

```@example latticemodel
k₁, k₂ = symbols("k₁ k₂", real=true)
model([k₁, k₂])
```

As for an [`OperatorGenerator`](@ref), the parameters of a [`Formula`](@ref) are obtained with [`Parameters`](@ref):

```@example latticemodel
Parameters(model)
```

With [`update!`](@ref), the parameters can be updated, and the model then can be re-evaluated with the new values:

```@example latticemodel
update!(model; μ=0)
model([k₁, k₂])
```

Since a [`Formula`](@ref) does not generate a list of operators, it does not support [`expand`](@ref) or the collection interface.

The generic interfaces of Section 6.2.3 apply to a [`Formula`](@ref) as well, with one peculiarity: the value type is the return type of the expression:

```@example latticemodel
valtype(model)
```

## 6.5 Worked Examples

Here we present some common examples. Except the Kitaev model, which must be defined on a honeycomb lattice, all of them are built on the one-dimensional periodic chain with a single point in the unitcell and a translation vector, just for simplicity. The `latexformat` lines in the following examples serve only to keep the typeset Hamiltonians compact on this page ([Section 4.3](@ref TutorialOperators)); they can be skipped when the examples are run in a REPL.

### 6.5.1 Fermionic and Bosonic Systems

#### Hubbard Model on a Chain

```@example latticemodel
lattice = Lattice([zero(Sym)]; vectors=[[one(Sym)]], name=:Chain)
hilbert = Hilbert(site => Fock{:f}(1, 2) for site in eachindex(lattice))
t = Hopping(:t, symbols("t", real=true), 1)
U = Hubbard(:U, symbols("U", real=true))
model = LatticeModel(lattice, hilbert, (t, U))

# set the custom latexformat
# one point per unitcell and one orbital: only the spin and the rcoordinate are needed
fock = CoordinatedIndex{<:Index{<:FockIndex{:f}}}
old = latexformat(fock)
latexformat(fock, LaTeX{(:nambu,), (:spinsym, :rcoordinate)}('c', "", ""))

# show the Hamiltonian
expand(model)
```

```@example latticemodel
# restore the default latexformat
latexformat(fock, old)
nothing # hide
```

The hopping term produces four operators on the single 1st-neighbor bond, one for each spin species in each of the two directions, and the Hubbard term produces one operator on the one-point bond, so that the Hamiltonian contains 5 operators in total.

#### Multi-Orbital Hubbard Model

For systems with multi-orbital degrees of freedom, increase the number of orbitals in the [`Fock`](@ref) constructor:

```@example latticemodel
lattice = Lattice([zero(Sym)]; vectors=[[one(Sym)]], name=:Chain)
hilbert = Hilbert(site => Fock{:f}(2, 2) for site in eachindex(lattice))
t = Hopping(:t, symbols("t", real=true), 1)
U = Hubbard(:U, symbols("U", real=true))
model = LatticeModel(lattice, hilbert, (t, U))

# set the custom latexformat
# one point per unitcell and two orbitals: the orbital is displayed as well
fock = CoordinatedIndex{<:Index{<:FockIndex{:f}}}
old = latexformat(fock)
latexformat(fock, LaTeX{(:nambu,), (:orbital, :spinsym, :rcoordinate)}('c', "", ""))

# show the Hamiltonian
expand(model)
```

```@example latticemodel
# restore the default latexformat
latexformat(fock, old)
nothing # hide
```

Every operator now carries an orbital index in addition to the spin one, so that the hopping term produces eight operators on the 1st-neighbor bond, and the Hubbard term two operators on the one-point bond, one for each orbital. The Hamiltonian contains 10 operators in total.

#### Bose-Hubbard Model

The same construction applies to bosonic systems, which differ only in the statistics of the [`Fock`](@ref) space and, since bosons usually carry no spin, in the number of spin species. The interaction between bosons is a four-operator coupling, which [`Onsite`](@ref) accepts as an explicit [`Coupling`](@ref):

```@example latticemodel
lattice = Lattice([zero(Sym)]; vectors=[[one(Sym)]], name=:Chain)
hilbert = Hilbert(site => Fock{:b}(1, 1) for site in eachindex(lattice))
t = Hopping(:t, symbols("t", real=true), 1)
U = Onsite(:U, symbols("U", real=true), Coupling(𝕒⁺(:, :, :), 𝕒⁺(:, :, :), 𝕒(:, :, :), 𝕒(:, :, :)))
model = LatticeModel(lattice, hilbert, (t, U))

# set the custom latexformat
# one point per unitcell and one spin species: only the rcoordinate is needed
boson = CoordinatedIndex{<:Index{<:FockIndex{:b}}}
old = latexformat(boson)
latexformat(boson, LaTeX{(:nambu,), (:rcoordinate,)}('b', "", ""))

# show the Hamiltonian
expand(model)
```

```@example latticemodel
# restore the default latexformat
latexformat(boson, old)
nothing # hide
```

The hopping term produces two operators, and the interaction term produces the single operator ``U\sum_i b^\dagger_i b^\dagger_i b_i b_i`` on the one-point bond, so that the Hamiltonian contains 3 operators in total.

### 6.5.2 Spin Systems

#### Heisenberg Model on a Chain

```@example latticemodel
lattice = Lattice([zero(Sym)]; vectors=[[one(Sym)]], name=:Chain)
hilbert = Hilbert(site => Spin{1//2}() for site in eachindex(lattice))
J = Heisenberg(:J, symbols("J", real=true), 1)
model = LatticeModel(lattice, hilbert, J)

# set the custom latexformat
# one point per unitcell: only the rcoordinate is needed
spin = CoordinatedIndex{<:Index{<:SpinIndex}}
old = latexformat(spin)
latexformat(spin, LaTeX{(:tag,), (:rcoordinate,)}('S', "", ""))

# show the Hamiltonian
expand(model)
```

The Heisenberg term expands on the 1st-neighbor bond into three operators, giving 3 in total. Further terms are combined in the same way. Adding, for instance, a magnetic field along ``z`` gives 4 operators:

```@example latticemodel
h = Zeeman(:h, symbols("h", real=true), 'z')
model = LatticeModel(lattice, hilbert, (J, h))
expand(model)
```

```@example latticemodel
# restore the default latexformat
latexformat(spin, old)
nothing # hide
```

#### Kitaev Model on the Honeycomb Lattice

The Kitaev model is defined on the honeycomb lattice. Its unitcell contains two points and its 1st-neighbor bonds point along three different directions. Since the coupling of each bond depends on its direction, the term takes the keyword arguments `x`, `y`, `z` that assign the three kinds of bonds by their azimuthal angles:

```@example latticemodel
lattice = Lattice(
    [zero(Sym), zero(Sym)], [zero(Sym), sqrt(Sym(3))/3];
    vectors=[[one(Sym), zero(Sym)], [one(Sym)/2, sqrt(Sym(3))/2]],
    name=:Honeycomb
)
hilbert = Hilbert(Spin{1//2}(), length(lattice))
K = Kitaev(:K, symbols("K", real=true), 1; x=[90], y=[210], z=[330])
model = LatticeModel(lattice, hilbert, K)

# set the custom latexformat
# two points per unitcell: the site is needed, and the default subscript delimiter is kept
spin = CoordinatedIndex{<:Index{<:SpinIndex}}
old = latexformat(spin)
latexformat(spin, LaTeX{(:tag,), (:site, :rcoordinate)}('S'))

# show the Hamiltonian
expand(model)
```

```@example latticemodel
# restore the default latexformat
latexformat(spin, old)
nothing # hide
```

Each of the three 1st-neighbor bonds carries a single Ising-type operator, so that the Hamiltonian contains 3 operators.

### 6.5.3 Phononic Systems

Phononic systems are built in the same way, from a phononic Hilbert space and phononic terms such as the kinetic and Hooke terms:

```@example latticemodel
lattice = Lattice([zero(Sym)]; vectors=[[one(Sym)]], name=:Chain)
hilbert = Hilbert(site => Phonon(1) for site in eachindex(lattice))
T = Kinetic(:T, 1)
V = Hooke(:k, symbols("k", real=true), 1)
model = LatticeModel(lattice, hilbert, (T, V))

# set the custom latexformat
# one point per unitcell and one vibration direction: only the rcoordinate is needed
phonon = CoordinatedIndex{<:Index{<:PhononIndex}}
old = latexformat(phonon)
latexformat(phonon, LaTeX{(), (:rcoordinate,)}(old.body, "", ""))

# show the Hamiltonian
expand(model)
```

```@example latticemodel
# restore the default latexformat
latexformat(phonon, old)
nothing # hide
```

The kinetic term acts on the one-point bonds and the Hooke term on the 1st-neighbor bonds; for a chain with a single vibration direction, they produce 5 operators in total.

## Summary

[`LatticeModel`](@ref) is the highest-level abstraction of the package:

- It bundles the spatial structure, the Hilbert space, the Hamiltonian, and the parameters into a single object.
- Its Hamiltonian may be represented by an [`OperatorGenerator`](@ref), a [`StaticGenerator`](@ref), or a [`Formula`](@ref), and algorithms work with any of them uniformly.
- Parameters are managed uniformly through [`Parameters`](@ref) and [`update!`](@ref).
- All representations share the generic interfaces of Section 6.1.

To learn how a [`LatticeModel`](@ref) connects to numerical algorithms, continue to [Chapter 7](@ref TutorialAlgorithmInterface).
