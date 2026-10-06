```@meta
CurrentModule = QuantumLattices
DocTestFilters = [r"im +[-\+]0\.0[-\+]", r"QuantumLattices(?:\.[A-Za-z_]\w*)*\."]
DocTestSetup = quote
    using QuantumLattices
    using StaticArrays: SVector
end
```

# [3. Internal Degrees of Freedom](@id TutorialDOF)

With the spatial structure defined, we now turn to the internal degrees of freedom.

## 3.1 Hierarchy of the Internal Degrees of Freedom

In general, a lattice Hamiltonian can be expressed in terms of the generators of an [algebra](https://en.wikipedia.org/wiki/Algebra_over_a_field) acting on the system's Hilbert space. For complex fermionic (bosonic) systems, the Hilbert space is the [Fock space](https://en.wikipedia.org/wiki/Fock_space), and the Hamiltonian is built from the generators of the fermionic (bosonic) algebra, i.e., the creation and annihilation operators $\{c^\dagger_\alpha, c_\alpha\}$ $\left(\{b^\dagger_\alpha, b_\alpha\}\right)$. For local spin-1/2 systems, the Hilbert space is $\otimes_\alpha\{\lvert\uparrow\rangle, \lvert\downarrow\rangle\}_\alpha$, and the Hamiltonian is built from the generators of the SU(2) spin algebra, i.e., the spin operators $\{S^x_\alpha, S^y_\alpha, S^z_\alpha\}$ or $\{S^+_\alpha, S^-_\alpha, S^z_\alpha\}$. In both cases, the subscript $\alpha$ represents a complete set of internal indices. Identifying the algebra and its generators is therefore the central task in constructing the operator representation of a lattice Hamiltonian.

The global Hilbert space of a lattice system can usually decompose into the direct product of local internal spaces living on individual lattice points, consequently, the global algebra decomposes similarly into local algebras. To incorporate the translation symmetry of the lattice, an intermediate representation is needed for the internal degrees of freedom that are equivalent under translation within the origin unitcell. Thus, from the microscopic to the macroscopic, the internal degrees of freedom are organized into a three-level hierarchy: local, unitcell, and global.

* **Local (individual-point) level**: the algebra is represented by [`Internal`](@ref), and a generator of that algebra by [`InternalIndex`](@ref). Both are abstract types; their concrete subtypes represent the specific algebras and generators of different quantum lattice systems.
* **Unitcell level**: the algebra is represented by [`Hilbert`](@ref), which assigns a concrete local algebra to each point within the origin unitcell. Accordingly, [`Index`](@ref) combines a site index with an [`InternalIndex`](@ref) to specify a translation-equivalent generator within the origin unitcell.
* **Global (whole-lattice) level**: a representation of the algebra is not needed, but the generators must be specified outside the origin unitcell when a bond crosses unitcell boundaries. The type [`CoordinatedIndex`](@ref) combines an [`Index`](@ref) with the coordinates $\mathbf{R}_{im}$ and $\mathbf{R}_m$ of the underlying point to represent a generator at this level.

It is noted that [`InternalIndex`](@ref), [`Index`](@ref) and [`CoordinatedIndex`](@ref) are all subtypes of the same abstract type, [`OperatorIndex`](@ref), the common root of the whole hierarchy of algebra generators. This design enables operators and their algebraic operations to be defined uniformly over all three levels in the next chapter.

The above discussions can be summarized by the following table, which also displays how the spatial part of a quantum lattice system is represented:

|            | Local (Individual-Point) Level                     | Unitcell Level    | Global (Whole-Lattice) Level |
|:----------:|:--------------------------------------------------:|:-----------------:|:----------------------------:|
| Spatial    | [`Point`](@ref)                                    | [`Lattice`](@ref) |                              |
| Algebra    | [`Internal`](@ref) and its concrete subtypes       | [`Hilbert`](@ref) |                              |
| Generator  | [`InternalIndex`](@ref)  and its concrete subtypes | [`Index`](@ref)   | [`CoordinatedIndex`](@ref)   |

## 3.2 Fermionic and Bosonic Systems

### 3.2.1 Local Level: Fock and FockIndex

Roughly speaking, these systems share similar internal structures of local Hilbert spaces termed as the [Fock space](https://en.wikipedia.org/wiki/Fock_space) where the generators of local algebras are the annihilation and creation operators. Besides the nambu index to distinguish whether it is an annihilation one or a creation one, such a generator usually adopts an orbital index and a spin index. Thus, the type [`FockIndex`](@ref)`<:`[`InternalIndex`](@ref), which specifies a certain local generator of a local Fock algebra, has the following attributes:
* `orbital::Int`: the orbital index
* `spin::Rational{Int}`: the spin index, which must be a half integer or integer
* `nambu::Int`: the nambu index, which must be 1 (annihilation) or 2 (creation).
Correspondingly, the type [`Fock`](@ref)`<:`[`Internal`](@ref), which specifies the local algebra acting on the local [Fock space](https://en.wikipedia.org/wiki/Fock_space), has the following attributes:
* `norbital::Int`: the number of allowed orbital indices
* `nspin::Int`: the number of allowed spin indices

To distinguish whether the system is a fermionic one or a bosonic one, [`FockIndex`](@ref) and [`Fock`](@ref) take a symbol `:f`(for fermionic) or `:b`(for bosonic) to be their first type parameters.

Now let's see some examples.

A [`FockIndex`](@ref) instance can be initialized by giving all its three attributes:
```jldoctest FFF
julia> FockIndex{:f}(2, 1//2, 1)
𝕔(2, 1//2)

julia> FockIndex{:f}(2, 1//2, 2)
𝕔⁺(2, 1//2)

julia> FockIndex{:b}(2, 0, 1)
𝕒(2, 0)

julia> FockIndex{:b}(2, 0, 2)
𝕒⁺(2, 0)
```

Here, `𝕔` (\bbc<tab>), `𝕔⁺` (\bbc<tab>\^+<tab>), `𝕒` (\bba<tab>) and `𝕒⁺` (\bba<tab>\^+<tab>) are functions that are convenient to construct and display instances of `FockIndex{:f}` and `FockIndex{:b}`, respectively.
```jldoctest FFF
julia> 𝕔(2, 1//2) isa FockIndex{:f}
true

julia> 𝕔⁺(2, 1//2) isa FockIndex{:f}
true

julia> 𝕒(2, 0) isa FockIndex{:b}
true

julia> 𝕒⁺(2, 0) isa FockIndex{:b}
true
```

The adjoint of a [`FockIndex`](@ref) instance is also defined:
```jldoctest FFF
julia> 𝕔(3, 3//2)'
𝕔⁺(3, 3//2)

julia> 𝕒⁺(3, 1)'
𝕒(3, 1)
```
Apparently, this operation is nothing but the "Hermitian conjugate".

A [`Fock`](@ref) instance can be initialized by giving all its attributes:
```jldoctest FFF
julia> Fock{:f}(1, 2)
4-element Fock{:f}:
 𝕔(1, -1//2)
 𝕔(1, 1//2)
 𝕔⁺(1, -1//2)
 𝕔⁺(1, 1//2)

julia> Fock{:b}(1, 1)
2-element Fock{:b}:
 𝕒(1, 0)
 𝕒⁺(1, 0)
```
As can be seen, a [`Fock`](@ref) instance behaves like a vector (because the parent type [`Internal`](@ref) is a subtype of `AbstractVector`), and its iteration just generates all the allowed [`FockIndex`](@ref) instances on its associated spatial point:
```jldoctest FFF
julia> fck = Fock{:f}(2, 1);

julia> fck |> typeof |> eltype
FockIndex{:f, Int64, Rational{Int64}}

julia> fck |> length
4

julia> [fck[1], fck[2], fck[3], fck[4]]
4-element Vector{FockIndex{:f, Int64, Rational{Int64}}}:
 𝕔(1, 0)
 𝕔(2, 0)
 𝕔⁺(1, 0)
 𝕔⁺(2, 0)

julia> fck |> collect
4-element Vector{FockIndex{:f, Int64, Rational{Int64}}}:
 𝕔(1, 0)
 𝕔(2, 0)
 𝕔⁺(1, 0)
 𝕔⁺(2, 0)
```
This is isomorphic to the mathematical fact that a local algebra is a vector space of the local generators.

The statistics of a [`FockIndex`](@ref)/[`Fock`](@ref) instance can be obtained by the [`statistics`](@ref) function:
```jldoctest FFF
julia> 𝕔(1, 0) |> statistics
:f

julia> 𝕒(1, 0) |> statistics
:b

julia> Fock{:f}(2, 2) |> statistics
:f

julia> Fock{:b}(2, 2) |> statistics
:b
```

### 3.2.2 Unitcell Level: Hilbert and Index

To specify the Fock algebra at the unitcell level, [`Hilbert`](@ref) associates each point within the origin unitcell with an instance of [`Fock`](@ref):
```jldoctest FFF
julia> Hilbert(1=>Fock{:f}(1, 2), 2=>Fock{:f}(1, 2))
Hilbert{Fock{:f}} with 2 entries:
  1 => Fock{:f}(norbital=1, nspin=2)
  2 => Fock{:f}(norbital=1, nspin=2)

julia> Hilbert(site=>Fock{:f}(2, 2) for site=1:2)
Hilbert{Fock{:f}} with 2 entries:
  1 => Fock{:f}(norbital=2, nspin=2)
  2 => Fock{:f}(norbital=2, nspin=2)

julia> Hilbert(Fock{:f}(2, 2), 2)
Hilbert{Fock{:f}} with 2 entries:
  1 => Fock{:f}(norbital=2, nspin=2)
  2 => Fock{:f}(norbital=2, nspin=2)

julia> Hilbert([Fock{:f}(2, 2), Fock{:f}(2, 2)])
Hilbert{Fock{:f}} with 2 entries:
  1 => Fock{:f}(norbital=2, nspin=2)
  2 => Fock{:f}(norbital=2, nspin=2)
```

In general, at different sites, the local Fock algebra could be different:
```jldoctest FFF
julia> Hilbert(site=>Fock{:f}(iseven(site) ? 2 : 1, 1) for site=1:2)
Hilbert{Fock{:f}} with 2 entries:
  1 => Fock{:f}(norbital=1, nspin=1)
  2 => Fock{:f}(norbital=2, nspin=1)

julia> Hilbert(1=>Fock{:f}(1, 2), 2=>Fock{:b}(1, 2))
Hilbert{Fock} with 2 entries:
  1 => Fock{:f}(norbital=1, nspin=2)
  2 => Fock{:b}(norbital=1, nspin=2)
```

[`Hilbert`](@ref) itself is a subtype of `AbstractDict`, the iteration over the keys gives the sites, and the iteration over the values gives the local algebras:
```jldoctest FFF
julia> hilbert = Hilbert(site=>Fock{:f}(iseven(site) ? 2 : 1, 1) for site=1:2);

julia> collect(keys(hilbert))
2-element Vector{Int64}:
 1
 2

julia> collect(values(hilbert))
2-element Vector{Fock{:f}}:
 Fock{:f}(norbital=1, nspin=1)
 Fock{:f}(norbital=2, nspin=1)

julia> collect(hilbert)
2-element Vector{Pair{Int64, Fock{:f}}}:
 1 => Fock{:f}(norbital=1, nspin=1)
 2 => Fock{:f}(norbital=2, nspin=1)

julia> [hilbert[1], hilbert[2]]
2-element Vector{Fock{:f}}:
 Fock{:f}(norbital=1, nspin=1)
 Fock{:f}(norbital=2, nspin=1)
```

To specify a translation-equivalent generator of the Fock algebra within the unitcell, [`Index`](@ref) just combines a `site::Int` attribute and an `internal::FockIndex` attribute:
```jldoctest FFF
julia> index = Index(1, FockIndex{:f}(1, -1//2, 2))
𝕔⁺(1, 1, -1//2)

julia> index.site
1

julia> index.internal
𝕔⁺(1, -1//2)
```

Here, the functions `𝕔`, `𝕔⁺`, `𝕒` and `𝕒⁺` can also construct and display instances of `Index{<:FockIndex{:f}}` and `Index{<:FockIndex{:b}}`, respectively.
```jldoctest FFF
julia> 𝕔(1, 1, -1//2) isa Index{<:FockIndex{:f}}
true

julia> 𝕔⁺(1, 1, -1//2) isa Index{<:FockIndex{:f}}
true

julia> 𝕒(1, 1, -1//2) isa Index{<:FockIndex{:b}}
true

julia> 𝕒⁺(1, 1, -1//2) isa Index{<:FockIndex{:b}}
true
```

The Hermitian conjugate and statistics of an [`Index`](@ref) instance is also defined:
```jldoctest FFF
julia> 𝕔⁺(1, 1, -1//2)'
𝕔(1, 1, -1//2)

julia> 𝕒(1, 1, -1//2)'
𝕒⁺(1, 1, -1//2)

julia> 𝕔(1, 1, -1//2) |> statistics
:f

julia> 𝕒(1, 1, -1//2) |> statistics
:b
```

### 3.2.3 Global Level: CoordinatedIndex

Since the local algebra of a quantum lattice system can be defined point by point, the global algebra can be completely compressed into the origin unitcell. However, generators outside the origin unitcell cannot be avoided because we have to use them to compose the Hamiltonian on the bonds that go across the unitcell boundaries. This situation is similar to the case of [`Lattice`](@ref) and [`Point`](@ref). Therefore, we take a similar solution for the generators to that adopted for the [`Point`](@ref), i.e., we include the $\mathbf{R}_{im}$ coordinate (by the `rcoordinate` attribute) and the $\mathbf{R}_m$ coordinate (by the `icoordinate` attribute) of the underlying point together with the `index::Index` attribute in the [`CoordinatedIndex`](@ref) type to represent a generator that could be inside or outside the origin unitcell:
```jldoctest FFF
julia> index = CoordinatedIndex(Index(1, FockIndex{:f}(1, 0, 2)), [0.5, 0.0], [0.0, 0.0])
𝕔⁺(1, 1, 0, [0.5, 0.0], [0.0, 0.0])

julia> index.index
𝕔⁺(1, 1, 0)

julia> index.rcoordinate
2-element SVector{2, Float64} with indices SOneTo(2):
 0.5
 0.0

julia> index.icoordinate
2-element SVector{2, Float64} with indices SOneTo(2):
 0.0
 0.0

julia> index'
𝕔(1, 1, 0, [0.5, 0.0], [0.0, 0.0])

julia> index |> statistics
:f
```

Here, as can be expected, the functions `𝕔`, `𝕔⁺`, `𝕒` and `𝕒⁺` can construct and display instances of `CoordinatedIndex{<:Index{<:FockIndex{:f}}}` and `CoordinatedIndex{<:Index{<:FockIndex{:b}}}`, respectively, as well.
```jldoctest FFF
julia> 𝕔(1, 1, 0, [0.5, 0.0], [0.0, 0.0]) isa CoordinatedIndex{<:Index{<:FockIndex{:f}}}
true

julia> 𝕔⁺(1, 1, 0, [0.5, 0.0], [0.0, 0.0]) isa CoordinatedIndex{<:Index{<:FockIndex{:f}}}
true

julia> 𝕒(1, 1, 0, [0.5, 0.0], [0.0, 0.0]) isa CoordinatedIndex{<:Index{<:FockIndex{:b}}}
true

julia> 𝕒⁺(1, 1, 0, [0.5, 0.0], [0.0, 0.0]) isa CoordinatedIndex{<:Index{<:FockIndex{:b}}}
true
```

## 3.3 Spin Systems

### 3.3.1 Local Level: Spin and SpinIndex

[`Spin`](@ref)`<:`[`Internal`](@ref) and [`SpinIndex`](@ref)`<:`[`InternalIndex`](@ref) are designed to deal with SU(2) spin systems at the local level.

Although spin systems are essentially bosonic, the commonly-used local Hilbert space is distinct from that of a usual bosonic system: it is the space spanned by the eigenstates of a local $S^z$ operator rather than a [Fock space](https://en.wikipedia.org/wiki/Fock_space). At the same time, a spin Hamiltonian is usually expressed by local spin operators, such as $S^x$, $S^y$, $S^z$, $S^+$ and $S^-$, instead of creation and annihilation operators. Therefore, it is convenient to define another set of concrete subtypes for spin systems.

To specify which one of the five $\{S^x, S^y, S^z, S^+, S^-\}$ a local spin operator is, the type [`SpinIndex`](@ref) has the following attribute:
* `tag::Char`: the tag, which must be `'x'`, `'y'`, `'z'`, `'+'` or `'-'`.
Correspondingly, the type [`Spin`](@ref), which defines the local SU(2) spin algebra, does not need any attribute.

For [`SpinIndex`](@ref) and [`Spin`](@ref), it is also necessary to know what the total spin is, which is taken as their first type parameters and should be a half-integer or an integer.

Now let's see examples.

A [`SpinIndex`](@ref) instance can be initialized as follows:
```jldoctest SSS
julia> SpinIndex{3//2}('x')
𝕊{3//2}('x')

julia> SpinIndex{1//2}('z')
𝕊{1//2}('z')

julia> SpinIndex{1}('+')
𝕊{1}('+')
```

Here, the type `𝕊` (\bbS<tab>) plays a similar role in spin systems as `𝕔`/`𝕔⁺` and `𝕒`/`𝕒⁺` in Fock systems.
```jldoctest SSS
julia> 𝕊{3//2}('x') isa SpinIndex{3//2}
true
```

The "Hermitian conjugate" of a [`SpinIndex`](@ref) instance can be obtained by the adjoint operation:
```jldoctest SSS
julia> 𝕊{3//2}('x')'
𝕊{3//2}('x')

julia> 𝕊{3//2}('y')'
𝕊{3//2}('y')

julia> 𝕊{3//2}('z')'
𝕊{3//2}('z')

julia> 𝕊{3//2}('+')'
𝕊{3//2}('-')

julia> 𝕊{3//2}('-')'
𝕊{3//2}('+')
```

The local spin space is determined by the total spin. The standard matrix representation of a [`SpinIndex`](@ref) instance on this local spin space can be obtained by the [`matrix`](@ref) function exported by this package:
```jldoctest SSS
julia> 𝕊{1//2}('x') |> matrix
2×2 Matrix{ComplexF64}:
 0.0+0.0im  0.5+0.0im
 0.5+0.0im  0.0+0.0im

julia> 𝕊{1//2}('y') |> matrix
2×2 Matrix{ComplexF64}:
  0.0-0.0im  0.0-0.5im
 -0.0+0.5im  0.0-0.0im

julia> 𝕊{1//2}('z') |> matrix
2×2 Matrix{ComplexF64}:
  0.5+0.0im   0.0+0.0im
 -0.0+0.0im  -0.5+0.0im

julia> 𝕊{1//2}('+') |> matrix
2×2 Matrix{ComplexF64}:
 0.0+0.0im  1.0+0.0im
 0.0+0.0im  0.0+0.0im

julia> 𝕊{1//2}('-') |> matrix
2×2 Matrix{ComplexF64}:
 0.0+0.0im  0.0+0.0im
 1.0+0.0im  0.0+0.0im
```

A [`Spin`](@ref) instance can be initialized as follows:
```jldoctest SSS
julia> Spin{1}()
3-element Spin{1}:
 𝕊{1}('x')
 𝕊{1}('y')
 𝕊{1}('z')

julia> Spin{1//2}()
3-element Spin{1//2}:
 𝕊{1//2}('x')
 𝕊{1//2}('y')
 𝕊{1//2}('z')
```

Similar to [`Fock`](@ref), a [`Spin`](@ref) instance behaves like a vector whose iteration generates the [`SpinIndex`](@ref) instances on its associated spatial point:
```jldoctest SSS
julia> sp = Spin{1}();

julia> sp |> typeof |> eltype
SpinIndex{1, Char}

julia> sp |> length
3

julia> [sp[1], sp[2], sp[3]]
3-element Vector{SpinIndex{1, Char}}:
 𝕊{1}('x')
 𝕊{1}('y')
 𝕊{1}('z')

julia> sp |> collect
3-element Vector{SpinIndex{1, Char}}:
 𝕊{1}('x')
 𝕊{1}('y')
 𝕊{1}('z')
```
It is noted that a [`Spin`](@ref) instance only generates [`SpinIndex`](@ref) instances limited to those of $S^x$, $S^y$, $S^z$ but not those of $S^+$ and $S^-$ because the former three have already formed a complete set of the generators of the local SU(2) spin algebra. This doesn't matter if you also want to construct your spin Hamiltonian by use of $S^+$ and $S^-$. For more details, see [Chapter 5](@ref TutorialCouplings).

The total spin of a [`SpinIndex`](@ref)/[`Spin`](@ref) instance can be obtained by the [`totalspin`](@ref) function:
```jldoctest SSS
julia> 𝕊{1//2}('x') |> totalspin
1//2

julia> Spin{1}() |> totalspin
1
```

### 3.3.2 Unitcell and Global Levels

At the unitcell and global levels to construct the SU(2) spin algebra and spin generators, it is completely the same to that of the Fock algebra and Fock generators as long as we replace [`Fock`](@ref) and [`FockIndex`](@ref) with [`Spin`](@ref) and [`SpinIndex`](@ref), respectively:
```jldoctest SSS
julia> Hilbert(1=>Spin{1//2}(), 2=>Spin{1}())
Hilbert{Spin} with 2 entries:
  1 => Spin{1//2}()
  2 => Spin{1}()

julia> index = Index(1, SpinIndex{1//2}('+'))
𝕊{1//2}(1, '+')

julia> index |> totalspin
1//2

julia> 𝕊{1//2}(1, '+') isa Index{<:SpinIndex{1//2}}
true

julia> index = CoordinatedIndex(Index(1, SpinIndex{1//2}('-')), [0.5, 0.5], [1.0, 1.0])
𝕊{1//2}(1, '-', [0.5, 0.5], [1.0, 1.0])

julia> index |> totalspin
1//2

julia> 𝕊{1//2}(1, '-', [0.5, 0.5], [1.0, 1.0]) isa CoordinatedIndex{<:Index{<:SpinIndex{1//2}}}
true
```

## 3.4 Phononic Systems

### 3.4.1 Local Level: Phonon and PhononIndex

Phononic systems are also bosonic systems. However, the canonical creation and annihilation operators of phonons depend on the eigenvalues and eigenvectors of the dynamical matrix, making them difficult to define locally at each point. Instead, we resort to the displacement ($\mathbf{u}$) and momentum ($\mathbf{p}$) operators of lattice vibrations as the generators, which can be easily defined locally. The type [`PhononIndex`](@ref)`<:`[`InternalIndex`](@ref) can specify such a local generator, which has the following attributes:
* `direction::Char`: the direction, which must be one of `'x'`, `'y'` and `'z'`, to indicate which spatial directional component of the generator it is
Correspondingly, the type [`Phonon`](@ref)`<:`[`Internal`](@ref), which defines the local $\{\mathbf{u}, \mathbf{p}\}$ algebra of the lattice vibrations, has the following attributes:
* `ndirection::Int`: the spatial dimension of the lattice vibrations, which must be 1, 2, or 3.

For [`PhononIndex`](@ref) and [`Phonon`](@ref), it is also necessary to distinguish whether it is for the displacement ($\mathbf{u}$) or for the momentum ($\mathbf{p}$). Their first type parameters are designed to solve this problem, with `:u` and `:p` denoting $\mathbf{u}$ and $\mathbf{p}$, respectively.

Now let's see examples:
```jldoctest PPP
julia> PhononIndex{:u}('x')
𝕦('x')

julia> PhononIndex{:p}('x')
𝕡('x')

julia> # one-dimensional lattice vibration only has the x component

julia> Phonon{:u}(1)
1-element Phonon{:u}:
 𝕦('x')

julia> Phonon{:p}(1)
1-element Phonon{:p}:
 𝕡('x')

julia> # two-dimensional lattice vibration only has the x and y components

julia> Phonon{:u}(2) 
2-element Phonon{:u}:
 𝕦('x')
 𝕦('y')

julia> Phonon{:p}(2)
2-element Phonon{:p}:
 𝕡('x')
 𝕡('y')

julia> # three-dimensional lattice vibration has the x, y and z components

julia> Phonon{:u}(3)
3-element Phonon{:u}:
 𝕦('x')
 𝕦('y')
 𝕦('z')

julia> Phonon{:p}(3)
3-element Phonon{:p}:
 𝕡('x')
 𝕡('y')
 𝕡('z')
```

As is usual, we define functions `𝕦` (\bbu<tab>) and `𝕡` (\bbp<tab>) to construct and display instances of `PhononIndex{:u}` and `PhononIndex{:p}` for convenience, respectively.
```jldoctest PPP
julia> 𝕦('x') isa PhononIndex{:u}
true

julia> 𝕡('x') isa PhononIndex{:p}
true
```

The kind (`:u` for displacement and `:p` for momentum) of a [`PhononIndex`](@ref)/[`Phonon`](@ref) instance can be obtained by the [`kind`](@ref) function:
```jldoctest PPP
julia> 𝕦('x') |> kind
:u

julia> 𝕡('x') |> kind
:p
```

```jldoctest PPP
julia> Phonon{:u}(3) |> kind
:u

julia> Phonon{:p}(3) |> kind
:p
```

### 3.4.2 Unitcell and Global Levels

At the unitcell and global levels, lattice-vibration algebras and generators are the same as in the previous cases by replacing [`Fock`](@ref) and [`FockIndex`](@ref) with [`Phonon`](@ref) and [`PhononIndex`](@ref):
```jldoctest PPP
julia> Hilbert(site=>Phonon(2) for site=1:3)
Hilbert{Phonon{:}} with 3 entries:
  1 => Phonon(ndirection=2)
  2 => Phonon(ndirection=2)
  3 => Phonon(ndirection=2)

julia> index = Index(1, PhononIndex{:u}('x'))
𝕦(1, 'x')

julia> index |> kind
:u

julia> 𝕦(1, 'x') isa Index{<:PhononIndex{:u}}
true

julia> index = Index(1, PhononIndex{:p}('x'))
𝕡(1, 'x')

julia> index |> kind
:p

julia> 𝕡(1, 'x') isa Index{<:PhononIndex{:p}}
true

julia> index = CoordinatedIndex(Index(1, PhononIndex{:u}('x')), [0.5, 0.5], [1.0, 1.0])
𝕦(1, 'x', [0.5, 0.5], [1.0, 1.0])

julia> index |> kind
:u

julia> 𝕦(1, 'x', [0.5, 0.5], [1.0, 1.0]) isa CoordinatedIndex{<:Index{<:PhononIndex{:u}}}
true

julia> index = CoordinatedIndex(Index(1, PhononIndex{:p}('x')), [0.5, 0.5], [1.0, 1.0])
𝕡(1, 'x', [0.5, 0.5], [1.0, 1.0])

julia> index |> kind
:p

julia> 𝕡(1, 'x', [0.5, 0.5], [1.0, 1.0]) isa CoordinatedIndex{<:Index{<:PhononIndex{:p}}}
true
```
It is noted that `Phonon{:}` is a special kind of [`Phonon`](@ref), which are used to specify the Hilbert space of lattice vibrations. In this way, both the displacement ($\mathbf{u}$) and the momentum ($\mathbf{p}$) degrees of freedom can be incorporated.

## [3.5 Table and Metric: Index Orders](@id TutorialTableAndMetric)

To construct matrix representations of operators, we need an ordering of their indices, which maps each [`Index`](@ref)/[`CoordinatedIndex`](@ref) to an integer sequence number. This is the role of [`Table`](@ref) and [`Metric`](@ref).

### 3.5.1 Metric and OperatorIndexToTuple

[`Metric`](@ref) (abstract `<: Function`) is a rule that converts an operator index into a value that can be compared and sorted. The concrete subtype provided by the package is [`OperatorIndexToTuple`](@ref), which converts an [`Index`](@ref)/[`CoordinatedIndex`](@ref) to a tuple element-by-element, as specified by the type parameter `Fields`. Each field can be either a `Symbol` (an attribute name of the internal index) or a `Function` (a trait function).

For fermionic/bosonic operators, common fields include:
- `:site` -- the site index
- `:orbital` -- the orbital index
- `:spin` -- the spin index
- `:nambu` -- the nambu index
- [`statistics`](@ref) -- the function that obtains the statistics of a fermionic/bosonic operator

For spin operators, common fields include:
- `:tag` -- the spin tag
- [`totalspin`](@ref) -- the function that obtains the total spin of a spin operator

For phononic operators, common fields include:
- `:direction` -- the direction
- [`kind`](@ref) -- the function that obtains the kind of a phononic operator

Let's see examples of constructing a [`Metric`](@ref) and applying it to an [`Index`](@ref)/[`CoordinatedIndex`](@ref):
```jldoctest table-metric
julia> metric = OperatorIndexToTuple(statistics, :nambu, :spin, :orbital, :site);

julia> keys(metric)
(statistics, :nambu, :spin, :orbital, :site)

julia> 𝕔⁺(1, 1, -1//2) |> metric
(:f, 2, -1//2, 1, 1)

julia> 𝕔⁺(1, 1, -1//2, [0.0], [0.0]) |> metric
(:f, 2, -1//2, 1, 1)
```

```jldoctest table-metric
julia> metric = OperatorIndexToTuple(totalspin, :site, :tag);

julia> keys(metric)
(totalspin, :site, :tag)

julia> 𝕊{1//2}(1, 'x') |> metric
(1//2, 1, 'x')

julia> 𝕊{1//2}(1, 'x', [0.0], [0.0]) |> metric
(1//2, 1, 'x')
```

```jldoctest table-metric
julia> metric = OperatorIndexToTuple(kind, :direction);

julia> keys(metric)
(kind, :direction)

julia> 𝕦(1, 'x') |> metric
(:u, 'x')

julia> 𝕦(1, 'x', [0.0], [0.0]) |> metric
(:u, 'x')

julia> 𝕡(1, 'y') |> metric
(:p, 'y')

julia> 𝕡(1, 'y', [0.0], [0.0]) |> metric
(:p, 'y')
```

Different `Fields` combinations produce different ordering conventions, which can be chosen to match the conventions of specific algorithms or physical setups.

### 3.5.2 Table

[`Table{I, B<:Metric}`](@ref) is an `OrderedDict` that maps each [`Index`](@ref)/[`CoordinatedIndex`](@ref) to an integer sequence number. The construction proceeds as:
1. Convert each [`Index`](@ref)/[`CoordinatedIndex`](@ref) via the [`Metric`](@ref) to obtain comparable values
2. Take the unique values and sort them
3. Assign sequence numbers based on the sorted order

A [`Table`](@ref) can be built from different sources.

**From a [`Hilbert`](@ref):**
```jldoctest table-metric
julia> hilbert = Hilbert(site=>Fock{:f}(1, 2) for site=1:2);

julia> metric = OperatorIndexToTuple(:spin, :site, :orbital);

julia> table = Table(hilbert, metric);

julia> length(table)
4

julia> table[𝕔(1, 1, 1//2)]
3

julia> table[𝕔⁺(1, 1, 1//2)]
3

julia> table[𝕔(2, 1, -1//2, [0.0], [0.0])]
2

julia> table[𝕔⁺(2, 1, -1//2, [0.0], [0.0])]
2
```

Note that with the above metric, the creation and annihilation indexes of the same mode share one sequence number: `𝕔(1, 1, 1//2)` and `𝕔⁺(1, 1, 1//2)` both map to `3`, because `OperatorIndexToTuple(:spin, :site, :orbital)` contains no `:nambu` field. This is usually the desired behavior when the sequence numbers label the basis of a Hamiltonian matrix, since both operators act on the same local Fock state. When the nambu index must be distinguished, include `:nambu` in the metric, as in the following example.

The lookup of a [`Table`](@ref) compares the metric values rather than the object identities. This is why a table built from a [`Hilbert`](@ref), whose keys are unitcell-level [`Index`](@ref)es, can be queried directly with a global-level [`CoordinatedIndex`](@ref), as in the example above, since the coordinates are simply invisible to the metric.

**From an `AbstractVector` of [`Index`](@ref)es:**
```jldoctest table-metric
julia> indexes = [𝕔(1, 1, 1//2), 𝕔(1, 1, -1//2), 𝕔⁺(1, 1, 1//2), 𝕔⁺(1, 1, -1//2)];

julia> metric = OperatorIndexToTuple(:nambu, :spin, :site, :orbital);

julia> table = Table(indexes, metric);

julia> length(table)
4

julia> table[𝕔(1, 1, 1//2)]
2

julia> table[𝕔⁺(1, 1, 1//2)]
4

julia> table[𝕔(1, 1, -1//2, [0.0], [0.0])]
1

julia> table[𝕔⁺(1, 1, -1//2, [0.0], [0.0])]
3
```

Two [`Table`](@ref) instances sharing the same [`Metric`](@ref) can union together:
```jldoctest table-metric
julia> metric = OperatorIndexToTuple(:nambu, :spin, :site, :orbital);

julia> table₁ = Table([𝕔(1, 1, 1//2), 𝕔(1, 1, -1//2)], metric);

julia> table₂ = Table([𝕔⁺(1, 1, 1//2), 𝕔⁺(1, 1, -1//2)], metric);

julia> merged = union(table₁, table₂);

julia> length(merged)
4
```

## Summary

The internal degrees of freedom are organized by the three-level hierarchy. At each level, the relevant types are summarized as follows:

| Level    | Fermionic/Bosonic                     | Spin                                  | Phononic                               | Purpose                                                |
|----------|---------------------------------------|---------------------------------------|----------------------------------------|--------------------------------------------------------|
| Local    | [`Fock`](@ref)/[`FockIndex`](@ref)    | [`Spin`](@ref)/[`SpinIndex`](@ref)    | [`Phonon`](@ref)/[`PhononIndex`](@ref) | Algebra and its generators at one point                |
| Unitcell | [`Hilbert`](@ref) and [`Index`](@ref) | [`Hilbert`](@ref) and [`Index`](@ref) | [`Hilbert`](@ref) and [`Index`](@ref)  | Translation-equivalent algebra and generators          |
| Global   | [`CoordinatedIndex`](@ref)            | [`CoordinatedIndex`](@ref)            | [`CoordinatedIndex`](@ref)             | Generators at arbitrary points of the infinite lattice |

The [`Hilbert`](@ref), [`Index`](@ref) and [`CoordinatedIndex`](@ref) types are system-agnostic: they apply to fermions/bosons, spins and phonons alike, with the concrete local types as their parameters. To build matrix representations or fix an ordering, the [`Metric`](@ref) rule converts an [`Index`](@ref)/[`CoordinatedIndex`](@ref) instance into a sortable value, and the [`Table`](@ref) assigns a sequence number to each [`Index`](@ref)/[`CoordinatedIndex`](@ref) instance.

With both the spatial structure and the internal degrees of freedom in hand, we are ready to construct operators out of these generators, which is the subject of [Chapter 4](@ref TutorialOperators).
