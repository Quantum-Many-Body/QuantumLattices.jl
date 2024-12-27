```@meta
CurrentModule = QuantumLattices
DocTestFilters = [r"im +[-\+]0\.0[-\+]"]
DocTestSetup = quote
    push!(LOAD_PATH, "../../../src/")
    using QuantumLattices
end
```

# Internal degrees of freedom

Now let's move to the second step, the internal degrees of freedom.

## Hierarchy of the internal degrees of freedom

In general, a lattice Hamiltonian can be expressed by the generators of the [algebra](https://en.wikipedia.org/wiki/Algebra_over_a_field) that acts on the Hilbert space of the system. For example for the complex fermionic (bosonic) system, the Hilbert space is the [Fock space](https://en.wikipedia.org/wiki/Fock_space), and the lattice Hamiltonian can be expressed by the generators of the fermionic (bosonic) algebra, i.e., the the creation and annihilation operators $\{c^\dagger_\alpha, c_\alpha\}$ $\left(\{b^\dagger_\alpha, b_\alpha\}\right)$. For another example for the local spin-1/2 system, the Hilbert space is the ``\otimes_\alpha\{\lvert\uparrow\rangle, \lvert\downarrow\rangle\}_\alpha`` space, and the lattice Hamiltonian can be expressed by the generators of the SU(2) spin algebra, i.e., the spin operators $\{S^x_\alpha, S^y_\alpha, S^z_\alpha\}$ or $\{S^+_\alpha, S^-_\alpha, S^z_\alpha\}$. In both examples, the subscript $\alpha$ denotes a complete set of indexes of the internal degrees of freedom of the quantum system. Therefore, the determination of the algebra acting on the system's Hilbert space and its corresponding generators lies at the center of the constructions of the operator representations of lattice Hamiltonians.

The global Hilbert space of a lattice system can be decomposed into the direct product of the local internal spaces "living" on individual points, leading to a similar decomposition of the global algebra into local ones. To incorporate with the unitcell construction of the lattice, an extra intermediate representation of the translation-equivalent internal degrees of freedom within the origin unitcell is also needed. Thus, from the microscopic to the macroscopic, we arrive at a three level hierarchy, namely the local-unitcell-global hierarchy, of the internal degrees of freedom.

At the local or the individual-point level, the local algebra is represented by the type [`Internal`](@ref), and a local generator of the local algebra is represented by the type [`InternalIndex`](@ref). Both types are abstract types with their concrete subtypes to represent concrete local algebras and concrete local generators of different quantum lattice systems with different internal structures, respectively.

At the unitcell level, the algebra of the system is represented by the type [`Hilbert`](@ref), which defines the concrete local algebras point by point within the origin unitcell. Accordingly, the type [`Index`](@ref), which combines a site index and an instance of [`InternalIndex`](@ref), could specify a translation-equivalent generator within the origin unitcell.

At the global or the whole-lattice level, we do not actually need a representation of the algebra of the system, but really do for the generators because we have to specify them outside the origin unitcell when the bond goes across the unitcell boundaries. The type [`CoordinatedIndex`](@ref), which combines an instance of [`Index`](@ref) and the coordinates $\mathbf{R}$ and $\mathbf{R}_i$ of the underlying point, represents a generator at such a level.

The above discussions can be summarized by the following table, which also displays how the spatial part of a quantum lattice system is represented:

|           | local (individual-point) level                    | unitcell level    | global (whole-lattice) level |
|:---------:|:-------------------------------------------------:|:-----------------:|:----------------------------:|
| spatial   | [`Point`](@ref)                                   | [`Lattice`](@ref) |                              |
| algebra   | [`Internal`](@ref) and its concrete subtypes      | [`Hilbert`](@ref) |                              |
| generator | [`InternalIndex`](@ref) and its concrete subtypes | [`Index`](@ref)   | [`CoordinatedIndex`](@ref)   |

## Quantum lattice systems with different internal structures

In this section, we will explain in detail for the common categories of quantum lattice systems implemented in this package about how their algebras and generators are organized according to the above three level hierarchy.

### Canonical complex fermionic, canonical complex bosonic and hard-core bosonic systems

#### Local level: Fock and FockIndex

Roughly speaking, these systems share similar internal structures of local Hilbert spaces termed as the [Fock space](https://en.wikipedia.org/wiki/Fock_space) where the generators of local algebras are the annihilation and creation operators. Besides the nambu index to distinguish whether it is an annihilation one or a creation one, such a generator usually adopts an orbital index and a spin index. Thus, the type [`FockIndex`](@ref)`<:`[`InternalIndex`](@ref), which specifies a certain local generator of a local Fock algebra, has the following attributes:
* `orbital::Int`: the orbital index
* `spin::Rational{Int}`: the spin index, which must be a half integer
* `nambu::Int`: the nambu index, which must be 1(annihilation) or 2(creation).
Correspondingly, the type [`Fock`](@ref)`<:`[`Internal`](@ref), which specifies the local algebra acting on the local [Fock space](https://en.wikipedia.org/wiki/Fock_space), has the following attributes:
* `norbital::Int`: the number of allowed orbital indices
* `nspin::Int`: the number of allowed spin indices

To distinguish whether the system is a fermionic one or a bosonic one, [`FockIndex`](@ref) and [`Fock`](@ref) take a symbol `:f`(for fermionic) or `:b`(for bosonic) to be their first type parameters.

Now let's see some examples.

An [`FockIndex`](@ref) instance can be initialized by giving all its three attributes:
```jldoctest FFF
julia> FockIndex{:f}(2, 1//2, 1)
𝕔(2, 1//2, 1)

julia> FockIndex{:b}(2, 0, 1)
𝕓(2, 0, 1)
```

Here, `𝕔` (\bbc<tab>) and `𝕓` (\bbb<tab>) are two functions that are convenient to construct and display instances of `FockIndex{:f}` and `FockIndex{:b}`, respectively.
```jldoctest FFF
julia> 𝕔(2, 1//2, 1) isa FockIndex{:f}
true

julia> 𝕓(2, 0, 1) isa FockIndex{:b}
true
```

The adjoint of an [`FockIndex`](@ref) instance is also defined:
```jldoctest FFF
julia> 𝕔(3, 3//2, 1)'
𝕔(3, 3//2, 2)

julia> 𝕓(3, 3//2, 2)'
𝕓(3, 3//2, 1)
```
Apparently, this operation is nothing but the "Hermitian conjugate".

A [`Fock`](@ref) instance can be initialized by giving all its attributes:
```jldoctest FFF
julia> Fock{:f}(1, 2)
4-element Fock{:f}:
 𝕔(1, -1//2, 1)
 𝕔(1, 1//2, 1)
 𝕔(1, -1//2, 2)
 𝕔(1, 1//2, 2)

julia> Fock{:b}(1, 1)
2-element Fock{:b}:
 𝕓(1, 0, 1)
 𝕓(1, 0, 2)
```
As can be seen, a [`Fock`](@ref) instance behaves like a vector (because the parent type [`Internal`](@ref) is a subtype of `AbstractVector`), and its iteration just generates all the allowed [`FockIndex`](@ref) instances on its associated spatial point:
```jldoctest FFF
julia> fck = Fock{:f}(2, 1);

julia> fck |> typeof |> eltype
FockIndex{:f, Int64, Rational{Int64}, Int64}

julia> fck |> length
4

julia> [fck[1], fck[2], fck[3], fck[4]]
4-element Vector{FockIndex{:f, Int64, Rational{Int64}, Int64}}:
 𝕔(1, 0, 1)
 𝕔(2, 0, 1)
 𝕔(1, 0, 2)
 𝕔(2, 0, 2)

julia> fck |> collect
4-element Vector{FockIndex{:f, Int64, Rational{Int64}, Int64}}:
 𝕔(1, 0, 1)
 𝕔(2, 0, 1)
 𝕔(1, 0, 2)
 𝕔(2, 0, 2)
```
This is isomorphic to the mathematical fact that a local algebra is a vector space of the local generators.

#### Unitcell level: Hilbert and Index

To specify the Fock algebra at the unitcell level, [`Hilbert`](@ref) associate each point within the origin unitcell with an instance of [`Fock`](@ref):
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
𝕔(1, 1, -1//2, 2)

julia> index.site
1

julia> index.internal
𝕔(1, -1//2, 2)
```

Here, the functions `𝕔` and `𝕓` can also construct and display instances of `Index{<:FockIndex{:f}}` and `Index{<:FockIndex{:b}}`, respectively.
```jldoctest FFF
julia> 𝕔(1, 1, -1//2, 2) isa Index{<:FockIndex{:f}}
true

julia> 𝕓(1, 1, -1//2, 2) isa Index{<:FockIndex{:b}}
true
```

The Hermitian conjugate of an [`Index`](@ref) is also defined:
```jldoctest FFF
julia> 𝕔(1, 1, -1//2, 2)'
𝕔(1, 1, -1//2, 1)

julia> 𝕓(1, 1, -1//2, 1)'
𝕓(1, 1, -1//2, 2)
```

#### Global level: CoordinatedIndex

Since the local algebra of a quantum lattice system can be defined point by point, the global algebra can be completely compressed into the origin unitcell. However, the generator outside the origin unitcell cannot be avoided because we have to use them to compose the Hamiltonian on the bonds that goes across the unitcell boundaries. This situation is similar to the case of [`Lattice`](@ref) and [`Point`](@ref). Therefore, we take a similar solution for the generators to that is adopted for the [`Point`](@ref), i.e., we include the $\mathbf{R}$ coordinate (by the `rcoordinate` attribute) and the $\mathbf{R}_i$ coordinate (by the `icoordinate` attribute) of the underlying point together with the `index::Index` attribute in the [`CoordinatedIndex`](@ref) type to represent a generator that could be inside or outside the origin unitcell:
```jldoctest FFF
julia> index = CoordinatedIndex(Index(1, FockIndex{:f}(1, 0, 2)), [0.5, 0.0], [0.0, 0.0])
𝕔(1, 1, 0, 2, [0.5, 0.0], [0.0, 0.0])

julia> index.index
𝕔(1, 1, 0, 2)

julia> index.rcoordinate
2-element StaticArraysCore.SVector{2, Float64} with indices SOneTo(2):
 0.5
 0.0

julia> index.icoordinate
2-element StaticArraysCore.SVector{2, Float64} with indices SOneTo(2):
 0.0
 0.0

julia> index' # the Hermitian conjugate of a CoordinatedIndex is also defined
𝕔(1, 1, 0, 1, [0.5, 0.0], [0.0, 0.0])
```

Here, as can be expected, the functions `𝕔` and `𝕓` can construct and display instances of `CoordinatedIndex{<:Index{<:FockIndex{:f}}}` and `CoordinatedIndex{<:Index{<:FockIndex{:b}}}`, respectively, as well.
```jldoctest FFF
julia> 𝕔(1, 1, 0, 2, [0.5, 0.0], [0.0, 0.0]) isa CoordinatedIndex{<:Index{<:FockIndex{:f}}}
true

julia> 𝕓(1, 1, 0, 2, [0.5, 0.0], [0.0, 0.0]) isa CoordinatedIndex{<:Index{<:FockIndex{:b}}}
true
```

### SU(2) spin systems

#### Local level: Spin and SpinIndex

[`Spin`](@ref)`<:`[`Internal`](@ref) and [`SpinIndex`](@ref)`<:`[`InternalIndex`](@ref) are designed to deal with SU(2) spin systems at the local level.

Although spin systems are essentially bosonic, the commonly-used local Hilbert space is distinct from that of an usual bosonic system: it is the space spanned by the eigenstates of a local $S^z$ operator rather than a [Fock space](https://en.wikipedia.org/wiki/Fock_space). At the same time, a spin Hamiltonian is usually expressed by local spin operators, such as $S^x$, $S^y$, $S^z$, $S^+$ and $S^-$, instead of creation and annihilation operators. Therefore, it is convenient to define another set of concrete subtypes for spin systems.

To specify which one of the five $\{S^x, S^y, S^z, S^+, S^-\}$ a local spin operator is, the type [`SpinIndex`](@ref) has the following attribute:
* `tag::Char`: the tag, which must be `'x'`, `'y'`, `'z'`, `'+'` or `'-'`.
Correspondingly, the type [`Spin`](@ref), which defines the local SU(2) spin algebra, does not need any attribute.

For [`SpinIndex`](@ref) and [`Spin`](@ref), it is also necessary to know what the total spin is, which is taken as their first type parameters and should be a half-integer or an integer.

Now let's see examples.

An [`SpinIndex`](@ref) instance can be initialized as follows
```jldoctest SSS
julia> SpinIndex{3//2}('x')
𝕊{3//2}('x')

julia> SpinIndex{1//2}('z')
𝕊{1//2}('z')

julia> SpinIndex{1}('+')
𝕊{1}('+')
```

Here, the type `𝕊` (\bbS<tab>) plays a similar role in spin systems as `𝕔` and `𝕓` in Fock system.
```jldoctest SSS
julia> 𝕊{3//2}('x') isa SpinIndex{3//2}
true
```

The "Hermitian conjugate" of an [`SpinIndex`](@ref) instance can be obtained by the adjoint operation:
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

The local spin space is determined by the total spin. The standard matrix representation of an [`SpinIndex`](@ref) instance on this local spin space can be obtained by the [`matrix`](@ref) function exported by this package:
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
5-element Spin{1}:
 𝕊{1}('x')
 𝕊{1}('y')
 𝕊{1}('z')
 𝕊{1}('+')
 𝕊{1}('-')

julia> Spin{1//2}()
5-element Spin{1//2}:
 𝕊{1//2}('x')
 𝕊{1//2}('y')
 𝕊{1//2}('z')
 𝕊{1//2}('+')
 𝕊{1//2}('-')
```

Similar to [`Fock`](@ref), a [`Spin`](@ref) instance behaves like a vector whose iteration generates all the allowed [`SpinIndex`](@ref) instances on its associated spatial point:
```jldoctest SSS
julia> sp = Spin{1}();

julia> sp |> typeof |> eltype
SpinIndex{1, Char}

julia> sp |> length
5

julia> [sp[1], sp[2], sp[3], sp[4], sp[5]]
5-element Vector{SpinIndex{1, Char}}:
 𝕊{1}('x')
 𝕊{1}('y')
 𝕊{1}('z')
 𝕊{1}('+')
 𝕊{1}('-')

julia> sp |> collect
5-element Vector{SpinIndex{1, Char}}:
 𝕊{1}('x')
 𝕊{1}('y')
 𝕊{1}('z')
 𝕊{1}('+')
 𝕊{1}('-')
```
It is noted that a [`Spin`](@ref) instance generates [`SpinIndex`](@ref) instances not only limited to those corresponding to $S^x$, $S^y$, $S^z$, but also those to $S^+$ and $S^-$ although the former three already forms a complete set of the generators of the local SU(2) spin algebra. This overcomplete feature is for the convenience to the construction of spin Hamiltonians.

#### Unitcell and global levels: Hilbert, Index and CoordinatedIndex

At the unitcell and global levels to construct the SU(2) spin algebra and spin generators, it is completely the same to that of the Fock algebra and Fock generators as long as we replace [`Fock`](@ref) and [`FockIndex`](@ref) with [`Spin`](@ref) and [`SpinIndex`](@ref), respectively:
```jldoctest SSS
julia> Hilbert(1=>Spin{1//2}(), 2=>Spin{1}())
Hilbert{Spin} with 2 entries:
  1 => Spin{1//2}()
  2 => Spin{1}()

julia> Index(1, SpinIndex{1//2}('+'))
𝕊{1//2}(1, '+')

julia> 𝕊{1//2}(1, '+') isa Index{<:SpinIndex{1//2}}
true

julia> CoordinatedIndex(Index(1, SpinIndex{1//2}('-')), [0.5, 0.5], [1.0, 1.0])
𝕊{1//2}(1, '-', [0.5, 0.5], [1.0, 1.0])

julia> 𝕊{1//2}(1, '-', [0.5, 0.5], [1.0, 1.0]) isa CoordinatedIndex{<:Index{<:SpinIndex{1//2}}}
true
```

### Phononic systems

#### Local level: Phonon and PhononIndex

Phononic systems are also bosonic systems. However, the canonical creation and annihilation operators of phonons depends on the eigenvalues and eigenvectors of the dynamical matrix, making them difficult to be defined locally at each point. Instead, we resort to the displacement ($\mathbf{u}$) and momentum ($\mathbf{p}$) operators of lattice vibrations as the generators, which can be easily defined locally. The type [`PhononIndex`](@ref)`<:`[`InternalIndex`](@ref) could specify such a local generator, which has the following attributes:
* `direction::Char`: the direction, which must be one of `'x'`, `'y'` and `'z'`, to indicate which spatial directional component of the generator it is
Correspondingly, the type [`Phonon`](@ref)`<:`[`Internal`](@ref), which defines the local $\{\mathbf{u}, \mathbf{p}\}$ algebra of the lattice vibrations, has the following attributes:
* `ndirection::Int`: the spatial dimension of the lattice vibrations, which must be 1, 2, or 3.

For [`PhononIndex`](@ref) and [`Phonon`](@ref), it is also necessary to distinguish whether it is for the displacement ($\mathbf{u}$) or for the momentum ($\mathbf{p}$). Their first type parameters are designed to solve this problem, with `:u` and `:b` denoting $\mathbf{u}$ and $\mathbf{p}$, respectively.

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

#### Unitcell and global levels: Hilbert, Index and CoordinatedIndex

At the unitcell and global levels, lattice-vibration algebras and generators are the same to previous situations by [`Phonon`](@ref) and [`PhononIndex`](@ref) replaced with in the corresponding types:
```jldoctest PPP
julia> Hilbert(site=>Phonon(2) for site=1:3)
Hilbert{Phonon{:}} with 3 entries:
  1 => Phonon(ndirection=2)
  2 => Phonon(ndirection=2)
  3 => Phonon(ndirection=2)

julia> Index(1, PhononIndex{:u}('x'))
𝕦(1, 'x')

julia> 𝕦(1, 'x') isa Index{<:PhononIndex{:u}}
true

julia> Index(1, PhononIndex{:p}('x'))
𝕡(1, 'x')

julia> 𝕡(1, 'x') isa Index{<:PhononIndex{:p}}
true

julia> CoordinatedIndex(Index(1, PhononIndex{:u}('x')), [0.5, 0.5], [1.0, 1.0])
𝕦(1, 'x', [0.5, 0.5], [1.0, 1.0])

julia> 𝕦(1, 'x', [0.5, 0.5], [1.0, 1.0]) isa CoordinatedIndex{<:Index{<:PhononIndex{:u}}}
true

julia> CoordinatedIndex(Index(1, PhononIndex{:p}('x')), [0.5, 0.5], [1.0, 1.0])
𝕡(1, 'x', [0.5, 0.5], [1.0, 1.0])

julia> 𝕡(1, 'x', [0.5, 0.5], [1.0, 1.0]) isa CoordinatedIndex{<:Index{<:PhononIndex{:p}}}
true
```
It is noted that `Phonon{:}` is a special kind of [`Phonon`](@ref), which are used to specify the Hilbert space of lattice vibrations. In this way, both the displacement ($\mathbf{u}$) and the momentum ($\mathbf{p}$) degrees of freedom can be incorporated.

## Operator and Operators

Now we arrive at the core types of this package, the [`Operator`](@ref) and [`Operators`](@ref). They are defined to deal with the mathematical operations, i.e., the `+`/`-`/`*` operations between two elements of the algebra acting on the Hilbert space, and the scalar multiplication between an element of the algebra and a number. Specifically, an [`Operator`](@ref) represents a **product** of several generators of the algebra specified at any of the three levels along with a coefficient, and an [`Operators`](@ref) represents the **sum** of several instances of [`Operator`](@ref).

[`Operator`](@ref) can be initialized by two ways:
```jldoctest OO
julia> Operator(2, 𝕔(1, -1//2, 2), 𝕔(1, -1//2, 1), 𝕊{1//2}('z'))
Operator(2, 𝕔(1, -1//2, 2), 𝕔(1, -1//2, 1), 𝕊{1//2}('z'))

julia> 2 * 𝕔(1, 1, -1//2, 2) * 𝕊{1//2}(2, 'z')
Operator(2, 𝕔(1, 1, -1//2, 2), 𝕊{1//2}(2, 'z'))
```
It is noted that the number of the generators can be any natural number.

Although generators at different levels can be producted to make an [`Operator`](@ref), it is not recommended to do so because the logic will be muddled:
```jldoctest OO
julia> Operator(2, 𝕔(1, 0, 2), 𝕔(2, 1, 0, 1, [0.0], [0.0])) # never do this !!!
Operator(2, 𝕔(1, 0, 2), 𝕔(2, 1, 0, 1, [0.0], [0.0]))
```

[`Operator`](@ref) can be iterated and indexed by integers, which will give the corresponding generators in the product:
```jldoctest OO
julia> op = Operator(2, 𝕔(1, 1//2, 2), 𝕔(1, 1//2, 1));

julia> length(op)
2

julia> [op[1], op[2]]
2-element Vector{FockIndex{:f, Int64, Rational{Int64}, Int64}}:
 𝕔(1, 1//2, 2)
 𝕔(1, 1//2, 1)

julia> collect(op)
2-element Vector{FockIndex{:f, Int64, Rational{Int64}, Int64}}:
 𝕔(1, 1//2, 2)
 𝕔(1, 1//2, 1)
```

To get the coefficient of an [`Operator`](@ref) or all its individual generators as a whole, use the [`value`](@ref) and [`id`](@ref) function exported by this package, respectively:
```jldoctest OO
julia> op = Operator(2, 𝕔(1, 0, 2), 𝕔(1, 0, 1));

julia> value(op)
2

julia> id(op)
(𝕔(1, 0, 2), 𝕔(1, 0, 1))
```

The product between two [`Operator`](@ref)s, or the scalar multiplication between a number and an [`Operator`](@ref) is also an [`Operator`](@ref):
```jldoctest OO
julia> Operator(2, 𝕔(1, 1//2, 2)) * Operator(3, 𝕔(1, 1//2, 1))
Operator(6, 𝕔(1, 1//2, 2), 𝕔(1, 1//2, 1))

julia> 3 * Operator(2, 𝕔(1, 1//2, 2))
Operator(6, 𝕔(1, 1//2, 2))

julia> Operator(2, 𝕔(1, 1//2, 2)) * 3
Operator(6, 𝕔(1, 1//2, 2))
```

The Hermitian conjugate of an [`Operator`](@ref) can be obtained by the adjoint operator:
```jldoctest OO
julia> op = Operator(6, 𝕔(2, 1//2, 2), 𝕔(1, 1//2, 1));

julia> op'
Operator(6, 𝕔(1, 1//2, 2), 𝕔(2, 1//2, 1))
```

There also exists a special [`Operator`](@ref), which only has the coefficient:
```jldoctest OO
julia> Operator(2)
Operator(2)
```

[`Operators`](@ref) can be initialized by two ways:
```jldoctest OO
julia> Operators(Operator(2, 𝕔(1, 1//2, 1)), Operator(3, 𝕔(1, 1//2, 2)))
Operators with 2 Operator
  Operator(2, 𝕔(1, 1//2, 1))
  Operator(3, 𝕔(1, 1//2, 2))

julia> Operator(2, 𝕔(1, 1//2, 1)) - Operator(3, 𝕓(1, 1//2, 2))
Operators with 2 Operator
  Operator(2, 𝕔(1, 1//2, 1))
  Operator(-3, 𝕓(1, 1//2, 2))
```

Similar items are automatically merged during the construction of [`Operators`](@ref):
```jldoctest OO
julia> Operators(Operator(2, 𝕔(1, 1//2, 1)), Operator(3, 𝕔(1, 1//2, 1)))
Operators with 1 Operator
  Operator(5, 𝕔(1, 1//2, 1))

julia> Operator(2, 𝕔(1, 1//2, 1)) + Operator(3, 𝕔(1, 1//2, 1))
Operators with 1 Operator
  Operator(5, 𝕔(1, 1//2, 1))
```

The multiplication between two [`Operators`](@ref)es, or between an [`Operators`](@ref) and an [`Operator`](@ref), or between a number and an [`Operators`](@ref) are defined:
```jldoctest OO
julia> ops = Operator(2, 𝕔(1, 1//2, 1)) + Operator(3, 𝕔(1, 1//2, 2));

julia> op = Operator(2, 𝕔(2, 1//2, 1));

julia> ops * op
Operators with 2 Operator
  Operator(4, 𝕔(1, 1//2, 1), 𝕔(2, 1//2, 1))
  Operator(6, 𝕔(1, 1//2, 2), 𝕔(2, 1//2, 1))

julia> op * ops
Operators with 2 Operator
  Operator(4, 𝕔(2, 1//2, 1), 𝕔(1, 1//2, 1))
  Operator(6, 𝕔(2, 1//2, 1), 𝕔(1, 1//2, 2))

julia> another = Operator(2, 𝕔(1, 1//2, 1)) + Operator(3, 𝕔(1, 1//2, 2));

julia> ops * another 
Operators with 2 Operator
  Operator(6, 𝕔(1, 1//2, 1), 𝕔(1, 1//2, 2))
  Operator(6, 𝕔(1, 1//2, 2), 𝕔(1, 1//2, 1))

julia> 2 * ops
Operators with 2 Operator
  Operator(4, 𝕔(1, 1//2, 1))
  Operator(6, 𝕔(1, 1//2, 2))

julia> ops * 2
Operators with 2 Operator
  Operator(4, 𝕔(1, 1//2, 1))
  Operator(6, 𝕔(1, 1//2, 2))
```
It is noted that in the result, the distributive law automatically applies. Besides, the fermion operator relation $c^2=c\dagger^2=0$ is also used.

As is usual, the Hermitian conjugate of an [`Operators`](@ref) can be obtained by the adjoint operator:
```jldoctest
julia> op₁ = Operator(6, 𝕔(1, 1//2, 2), 𝕔(2, 1//2, 1));

julia> op₂ = Operator(4, 𝕔(1, 1//2, 1), 𝕔(2, 1//2, 1));

julia> ops = op₁ + op₂;

julia> ops'
Operators with 2 Operator
  Operator(6, 𝕔(2, 1//2, 2), 𝕔(1, 1//2, 1))
  Operator(4, 𝕔(2, 1//2, 2), 𝕔(1, 1//2, 2))
```

[`Operators`](@ref) can be iterated and indexed:
```jldoctest OO
julia> ops = Operator(2, 𝕔(1, 1//2, 1)) + Operator(3, 𝕔(1, 1//2, 2));

julia> collect(ops)
2-element Vector{Operator{Int64, Tuple{FockIndex{:f, Int64, Rational{Int64}, Int64}}}}:
 Operator(2, 𝕔(1, 1//2, 1))
 Operator(3, 𝕔(1, 1//2, 2))

julia> ops[1]
Operator(2, 𝕔(1, 1//2, 1))

julia> ops[2]
Operator(3, 𝕔(1, 1//2, 2))
```
The index order of an [`Operators`](@ref) is the insertion order of the operators it contains.
