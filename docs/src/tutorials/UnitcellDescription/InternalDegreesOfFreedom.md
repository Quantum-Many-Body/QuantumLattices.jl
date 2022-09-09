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

In general, a lattice Hamiltonian can be expressed by the generators of the [algebra](https://en.wikipedia.org/wiki/Algebra_over_a_field) that acts on the Hilbert space of the system. For example for the complex fermionic (bosonic) system, the Hilbert space is the [Fock space](https://en.wikipedia.org/wiki/Fock_space), and the lattice Hamiltonian can be expressed by the generators of the fermionic (bosonic) algebra, i.e., the the creation and annihilation operators $\{c^\dagger_\alpha, c_\alpha\}$ $\left(\{b^\dagger_\alpha, b_\alpha\}\right)$. For another example for the local spin-1/2 system, the Hilbert space is the ``\otimes_\alpha\{\lvert\uparrow\rangle, \lvert\downarrow\rangle\}_\alpha`` space, and the lattice Hamiltonian can be expressed by the generators of the SU(2) spin algebra, i.e., the spin operators $\{S^x_\alpha, S^y_\alpha, S^z_\alpha\}$ or $\{S^+_\alpha, S^-_\alpha, S^z_\alpha\}$. In both examples, the subscript $\alpha$ denotes a complete set of indexes of the internal degrees of freedom of the quantum system. Therefore, the determination of the algebra acting on the system's Hilbert space and its corresponding generators lies at the center of the construction of lattice Hamiltonians.

The global Hilbert space of a lattice system can be decomposed into the direct product of the local internal spaces "living" on individual points, leading to a similar decomposition of the global algebra into local ones. To incorporate with the unitcell construction of the lattice, an extra intermediate representation of the translation-equivalent internal degrees of freedom within the origin unitcell is also needed. Thus, from the microscopic to the macroscopic, we arrive at a three level hierarchy, namely the local-unitcell-global hierarchy, of the internal degrees of freedom.

At the local or the individual-point level, the local algebra is represented by the type [`Internal`](@ref), and a local generator of the local algebra is represented by the type [`IID`](@ref). Both types are abstract types with their concrete subtypes to represent concrete local algebras and concrete local generators of different quantum lattice systems with different internal structures, respectively.

At the unitcell level, the algebra of the system is represented by the type [`Hilbert`](@ref), which defines the concrete local algebras point by point within the origin unitcell. Accordingly, the type [`Index`](@ref), which combines a site index and an instance of [`IID`](@ref), could specify a translation-equivalent generator within the origin unitcell.

At the global or the whole-lattice level, we do not actually need a representation of the algebra of the system, but really do for the generators because we have to specify them outside the origin unitcell when the bond goes across the unitcell boundaries. The type [`CompositeIndex`](@ref), which combines an instance of [`Index`](@ref) and the coordinates $\mathbf{R}$ and $\mathbf{R}_i$ of the underlying point, represents a generator at such a level.

The above discussions can be summarized by the following table, which also displays how the spatial part of a quantum lattice system is represented:

|           | local (individual-point) level               | unitcell level    | global (whole-lattice) level |
|:---------:|:--------------------------------------------:|:-----------------:|:----------------------------:|
| spatial   | [`Point`](@ref)                              | [`Lattice`](@ref) |                              |
| algebra   | [`Internal`](@ref) and its concrete subtypes | [`Hilbert`](@ref) |                              |
| generator | [`IID`](@ref) and its concrete subtypes      | [`Index`](@ref)   | [`CompositeIndex`](@ref)     |

## Quantum lattice systems with different internal structures

In this section, we will explain in detail for three common categories of quantum lattice systems implemented in this package about how their algebras and generators are organized according to the above three level hierarchy.

### Canonical complex fermionic, canonical bosonic and hard-core bosonic systems

#### Local level: Fock and FID

Roughly speaking, these systems share similar internal structures of local Hilbert spaces termed as the [Fock space](https://en.wikipedia.org/wiki/Fock_space) where the generators of local algebras are the annihilation and creation operators. Besides the nambu index to distinguish whether it is an annihilation one or a creation one, such a generator usually adopts an orbital index and a spin index. Thus, the type [`FID`](@ref)`<:`[`IID`](@ref), which specifies a certain local generator of a local Fock algebra, has the following attributes:
* `orbital::Int`: the orbital index
* `spin::Int`: the spin index
* `nambu::Int`: the nambu index, which must be 1(annihilation) or 2(creation).
Correspondingly, the type [`Fock`](@ref)`<:`[`Internal`](@ref), which specifies the local algebra acting on the local [Fock space](https://en.wikipedia.org/wiki/Fock_space), has the following attributes:
* `norbital::Int`: the number of allowed orbital indices
* `nspin::Int`: the number of allowed spin indices

To distinguish whether the system is a fermionic one or a bosonic one, [`FID`](@ref) and [`Fock`](@ref) take a symbol `:f`(for fermionic) or `:b`(for bosonic) to be their first type parameters.

Now let's see some examples.

An [`FID`](@ref) instance can be initialized by giving all its three attributes:
```jldoctest FFF
julia> FID{:f}(2, 2, 1)
FID{:f}(2, 2, 1)

julia> FID{:b}(2, 2, 1)
FID{:b}(2, 2, 1)
```

The adjoint of an [`FID`](@ref) instance is also defined:
```jldoctest FFF
julia> FID{:f}(3, 4, 1)'
FID{:f}(3, 4, 2)

julia> FID{:b}(3, 4, 2)'
FID{:b}(3, 4, 1)
```
Apparently, this operation is nothing but the "Hermitian conjugate".

A [`Fock`](@ref) instance can be initialized by giving all its attributes or by key word arguments to specify those beyond the default values:
```jldoctest FFF
julia> Fock{:f}(1, 2)
4-element Fock{:f}:
 FID{:f}(1, 1, 1)
 FID{:f}(1, 2, 1)
 FID{:f}(1, 1, 2)
 FID{:f}(1, 2, 2)
```
As can be seen, a [`Fock`](@ref) instance behaves like a vector (because the parent type [`Internal`](@ref) is a subtype of `AbstractVector`), and its iteration just generates all the allowed [`FID`](@ref) instances on its associated spatial point:
```jldoctest FFF
julia> fck = Fock{:f}(2, 1);

julia> fck |> typeof |> eltype
FID{:f, Int64, Int64, Int64}

julia> fck |> length
4

julia> [fck[1], fck[2], fck[3], fck[4]]
4-element Vector{FID{:f, Int64, Int64, Int64}}:
 FID{:f}(1, 1, 1)
 FID{:f}(2, 1, 1)
 FID{:f}(1, 1, 2)
 FID{:f}(2, 1, 2)

julia> fck |> collect
4-element Vector{FID{:f, Int64, Int64, Int64}}:
 FID{:f}(1, 1, 1)
 FID{:f}(2, 1, 1)
 FID{:f}(1, 1, 2)
 FID{:f}(2, 1, 2)
```
This is isomorphic to the mathematical fact that a local algebra is a vector space of the local generators.

#### Unitcell level: Hilbert{<:Fock} and Index{<:FID}

To specify the Fock algebra at the unitcell level, [`Hilbert`](@ref) associate each point within the origin unitcell with an instance of [`Fock`](@ref):
```jldoctest FFF
julia> Hilbert(1=>Fock{:f}(1, 2), 2=>Fock{:f}(1, 2))
Hilbert{Fock{:f}} with 2 entries:
  2 => Fock{:f}(norbital=1, nspin=2)
  1 => Fock{:f}(norbital=1, nspin=2)

julia> Hilbert(site=>Fock{:f}(2, 2) for site=1:2)
Hilbert{Fock{:f}} with 2 entries:
  2 => Fock{:f}(norbital=2, nspin=2)
  1 => Fock{:f}(norbital=2, nspin=2)

julia> Hilbert(site->Fock{:f}(1, 1), 1:2)
Hilbert{Fock{:f}} with 2 entries:
  2 => Fock{:f}(norbital=1, nspin=1)
  1 => Fock{:f}(norbital=1, nspin=1)
```

In general, at different sites, the local Fock algebra could be different:
```jldoctest FFF
julia> Hilbert(site=>Fock{:f}(iseven(site) ? 2 : 1, 1) for site=1:2)
Hilbert{Fock{:f}} with 2 entries:
  2 => Fock{:f}(norbital=2, nspin=1)
  1 => Fock{:f}(norbital=1, nspin=1)

julia> Hilbert(1=>Fock{:f}(1, 2), 2=>Fock{:b}(1, 2))
Hilbert{Fock} with 2 entries:
  2 => Fock{:b}(norbital=1, nspin=2)
  1 => Fock{:f}(norbital=1, nspin=2)
```

[`Hilbert`](@ref) itself is a subtype of `AbstractDict`, the iteration over the keys gives the sites, and the iteration over the values gives the local algebras:
```jldoctest FFF
julia> hilbert = Hilbert(site=>Fock{:f}(iseven(site) ? 2 : 1, 1) for site=1:2);

julia> collect(keys(hilbert))
2-element Vector{Int64}:
 2
 1

julia> collect(values(hilbert))
2-element Vector{Fock{:f}}:
 Fock{:f}(norbital=2, nspin=1)
 Fock{:f}(norbital=1, nspin=1)

julia> collect(hilbert)
2-element Vector{Pair{Int64, Fock{:f}}}:
 2 => Fock{:f}(norbital=2, nspin=1)
 1 => Fock{:f}(norbital=1, nspin=1)

julia> [hilbert[1], hilbert[2]]
2-element Vector{Fock{:f}}:
 Fock{:f}(norbital=1, nspin=1)
 Fock{:f}(norbital=2, nspin=1)
```

To specify a translation-equivalent generator of the Fock algebra within the unitcell, [`Index`](@ref) just combines a `:site::Int` attribute and an `:iid::FID` attribute:
```jldoctest FFF
julia> index = Index(1, FID{:f}(1, 1, 2))
Index(1, FID{:f}(1, 1, 2))

julia> index.site
1

julia> index.iid
FID{:f}(1, 1, 2)
```

The Hermitian conjugate of an [`Index`](@ref) is also defined:
```jldoctest FFF
julia> Index(1, FID{:f}(1, 1, 2))'
Index(1, FID{:f}(1, 1, 1))
```

#### Global level: CompositeIndex{<:Index{<:FID}}

Since the local algebra of a quantum lattice system can be defined point by point, the global algebra can be completely compressed into the origin unitcell. However, the generator outside the origin unitcell cannot be avoided because we have to use them to compose the Hamiltonian on the bonds that goes across the unitcell boundaries. This situation is the same to the case of [`Lattice`](@ref) and [`Point`](@ref). Therefore, we take a similar solution for the generators to that is used for the [`Point`](@ref), i.e., we include the $\mathbf{R}$ coordinate (the `:rcoordinate` attribute) and the $\mathbf{R}_i$ coordinate (the `:icoordinate` attribute) of the underlying point together with the `:index::Index` attribute in the [`CompositeIndex`](@ref) type to represent a generator that could be inside or outside the origin unitcell:
```jldoctest FFF
julia> index = CompositeIndex(Index(1, FID{:f}(1, 1, 2)), [0.5, 0.0], [0.0, 0.0])
CompositeIndex(Index(1, FID{:f}(1, 1, 2)), [0.5, 0.0], [0.0, 0.0])

julia> index.index
Index(1, FID{:f}(1, 1, 2))

julia> index.rcoordinate
2-element StaticArrays.SVector{2, Float64} with indices SOneTo(2):
 0.5
 0.0

julia> index.icoordinate
2-element StaticArrays.SVector{2, Float64} with indices SOneTo(2):
 0.0
 0.0

julia> index' # the Hermitian conjugate of a CompositeIndex is also defined
CompositeIndex(Index(1, FID{:f}(1, 1, 1)), [0.5, 0.0], [0.0, 0.0])
```

### SU(2) spin systems

#### Local level: Spin and SID

[`Spin`](@ref)`<:`[`Internal`](@ref) and [`SID`](@ref)`<:`[`IID`](@ref) are designed to deal with SU(2) spin systems at the local level.

Although spin systems are essentially bosonic, the commonly-used local Hilbert space is distinct from that of a canonical bosonic system: it is the space spanned by the eigenstates of a local $S^z$ operator rather than a [Fock space](https://en.wikipedia.org/wiki/Fock_space). At the same time, a spin Hamiltonian is usually expressed by local spin operators, such as $S^x$, $S^y$, $S^z$, $S^+$ and $S^-$, instead of creation and annihilation operators. Therefore, it is convenient to define another set of concrete subtypes for spin systems.

To specify which one of the five $\{S^x, S^y, S^z, S^+, S^-\}$ a local spin operator is, the type [`SID`](@ref) has the following attribute:
* `tag::Char`: the tag, which must be `'x'`, `'y'`, `'z'`, `'+'` or `'-'`.
Correspondingly, the type [`Spin`](@ref), which defines the local SU(2) spin algebra, has no attributes.

For [`SID`](@ref) and [`Spin`](@ref), it is also necessary to know what the total spin is, which is taken as their first type parameters and should be a half-integer or an integer.

Now let's see examples.

An [`SID`](@ref) instance can be initialized as follows
```jldoctest SSS
julia> SID{3//2}('x')
SID{3//2}('x')

julia> SID{1//2}('z')
SID{1//2}('z')

julia> SID{1}('+')
SID{1}('+')
```

The "Hermitian conjugate" of an [`SID`](@ref) instance can be obtained by the adjoint operation:
```jldoctest SSS
julia> SID{3//2}('x')'
SID{3//2}('x')

julia> SID{3//2}('y')'
SID{3//2}('y')

julia> SID{3//2}('z')'
SID{3//2}('z')

julia> SID{3//2}('+')'
SID{3//2}('-')

julia> SID{3//2}('-')'
SID{3//2}('+')
```

The local spin space is determined by the total spin. The standard matrix representation of an [`SID`](@ref) instance on this local spin space can be obtained by the [`matrix`](@ref) function exported by this package:
```jldoctest SSS
julia> SID{1//2}('x') |> matrix
2×2 Matrix{ComplexF64}:
 0.0+0.0im  0.5+0.0im
 0.5+0.0im  0.0+0.0im

julia> SID{1//2}('y') |> matrix
2×2 Matrix{ComplexF64}:
 0.0-0.0im  -0.0+0.5im
 0.0-0.5im   0.0-0.0im

julia> SID{1//2}('z') |> matrix
2×2 Matrix{ComplexF64}:
 -0.5+0.0im  -0.0+0.0im
  0.0+0.0im   0.5+0.0im

julia> SID{1//2}('+') |> matrix
2×2 Matrix{ComplexF64}:
 0.0+0.0im  0.0+0.0im
 1.0+0.0im  0.0+0.0im

julia> SID{1//2}('-') |> matrix
2×2 Matrix{ComplexF64}:
 0.0+0.0im  1.0+0.0im
 0.0+0.0im  0.0+0.0im
```

A [`Spin`](@ref) instance can be initialized as follows:
```jldoctest SSS
julia> Spin{1}()
5-element Spin{1}:
 SID{1}('x')
 SID{1}('y')
 SID{1}('z')
 SID{1}('+')
 SID{1}('-')

julia> Spin{1//2}()
5-element Spin{1//2}:
 SID{1//2}('x')
 SID{1//2}('y')
 SID{1//2}('z')
 SID{1//2}('+')
 SID{1//2}('-')
```

Similar to [`Fock`](@ref), a [`Spin`](@ref) instance behaves like a vector whose iteration generates all the allowed [`SID`](@ref) instances on its associated spatial point:
```jldoctest SSS
julia> sp = Spin{1}();

julia> sp |> typeof |> eltype
SID{1, Char}

julia> sp |> length
5

julia> [sp[1], sp[2], sp[3], sp[4], sp[5]]
5-element Vector{SID{1, Char}}:
 SID{1}('x')
 SID{1}('y')
 SID{1}('z')
 SID{1}('+')
 SID{1}('-')

julia> sp |> collect
5-element Vector{SID{1, Char}}:
 SID{1}('x')
 SID{1}('y')
 SID{1}('z')
 SID{1}('+')
 SID{1}('-')
```
It is noted that a [`Spin`](@ref) instance generates [`SID`](@ref) instances not only limited to those corresponding to $S^x$, $S^y$, $S^z$, but also those to $S^+$ and $S^-$ although the former three already forms a complete set of the generators of the local SU(2) spin algebra. This overcomplete feature is for the convenience to the construction of spin Hamiltonians.

#### Unitcell and global levels: Hilbert{<:Spin}, Index{<:SID} and CompositeIndex{<:Index{<:SID}}

At the unitcell and global levels to construct the SU(2) spin algebra and spin generators, it is completely the same to that of the Fock algebra and Fock generators as long as we replace [`Fock`](@ref) and [`FID`](@ref) with [`Spin`](@ref) and [`SID`](@ref), respectively:
```jldoctest SSS
julia> Hilbert(1=>Spin{1//2}(), 2=>Spin{1}())
Hilbert{Spin} with 2 entries:
  2 => Spin{1}()
  1 => Spin{1//2}()

julia> Index(1, SID{1//2}('+'))
Index(1, SID{1//2}('+'))

julia> CompositeIndex(Index(1, SID{1//2}('-')), [0.5, 0.5], [1.0, 1.0])
CompositeIndex(Index(1, SID{1//2}('-')), [0.5, 0.5], [1.0, 1.0])
```

### Phononic systems

#### Local level: Phonon and PID

Phononic systems are also bosonic systems. However, the canonical creation and annihilation operators of phonons depends on the eigenvalues and eigenvectors of the dynamical matrix, making them difficult to be defined locally at each point. Instead, we resort to the momentum $\mathbf{p}$ and displacement $\mathbf{u}$ operators of lattice vibrations as the generators, which can be easily defined locally. The type [`PID`](@ref)`<:`[`IID`](@ref) could specify such a local generator, which has the following attributes:
* `:tag::Char`: the tag, which must be either `'p'` or `'u'`, to specify whether it is the momentum or the displacement operator, respectively
* `:direction::Char`: the direction, which must be one of `'x'`, `'y'` and `'z'`, to indicate which spatial directional component of the generator is
Correspondingly, the type [`Phonon`](@ref)`<:`[`Internal`](@ref), which defines the local $\{\mathbf{u}, \mathbf{p}\}$ algebra of the lattice vibrations, has the following attributes:
* `ndirection::Int`: the spatial dimension of the lattice vibrations, which must be 1, 2, or 3

Now let's see examples:
```jldoctest PPP
julia> PID('u', 'x')
PID('u', 'x')

julia> PID('p', 'x')
PID('p', 'x')

julia> Phonon(1) # one-dimensional lattice vibration only has the x component
2-element Phonon:
 PID('u', 'x')
 PID('p', 'x')

julia> Phonon(2) # two-dimensional lattice vibration only has the x and y components
4-element Phonon:
 PID('u', 'x')
 PID('p', 'x')
 PID('u', 'y')
 PID('p', 'y')

julia> Phonon(3) # three-dimensional lattice vibration has the x, y and z components
6-element Phonon:
 PID('u', 'x')
 PID('p', 'x')
 PID('u', 'y')
 PID('p', 'y')
 PID('u', 'z')
 PID('p', 'z')
```

#### Unitcell and global levels: Hilbert{<:Phonon}, Index{<:PID} and CompositeIndex{<:Index{<:PID}}

At the unitcell and global levels, lattice-vibration algebras and generators are the same to previous situations with [`Phonon`](@ref) and [`PID`](@ref) replaced with:
```jldoctest PPP
julia> Hilbert(site=>Phonon(2) for site=1:3)
Hilbert{Phonon} with 3 entries:
  2 => Phonon(ndirection=2)
  3 => Phonon(ndirection=2)
  1 => Phonon(ndirection=2)

julia> Index(1, PID('u', 'x'))
Index(1, PID('u', 'x'))

julia> CompositeIndex(Index(1, PID('u', 'x')), [0.5, 0.5], [1.0, 1.0])
CompositeIndex(Index(1, PID('u', 'x')), [0.5, 0.5], [1.0, 1.0])
```

## Operator and Operators

Now we arrive at the core types of this package, the [`Operator`](@ref) and [`Operators`](@ref). They are defined to deal with the mathematical operations, i.e., the `+`/`-`/`*` operations between two elements of the algebra, and the scalar multiplication between a element of the algebra and a number. Specifically, an [`Operator`](@ref) represents a **product** of several generators specified at any of the three levels along with a coefficient, and an [`Operators`](@ref) represents the **sum** of several instances of [`Operator`](@ref).

[`Operator`](@ref) can be initialized by two ways:
```jldoctest OO
julia> Operator(2, FID{:f}(1, 1, 2), FID{:f}(1, 1, 1), SID{1//2}('z'))
Operator(2, FID{:f}(1, 1, 2), FID{:f}(1, 1, 1), SID{1//2}('z'))

julia> 2 * Index(1, FID{:f}(1, 1, 2)) * Index(2, SID{1//2}('z'))
Operator(2, Index(1, FID{:f}(1, 1, 2)), Index(2, SID{1//2}('z')))
```
It is noted that the number of the generators can be any natural number.

Although generators at different levels can be producted to make a generator, it is not recommended to do so because the logic will be muddled:
```jldoctest OO
julia> Operator(2, FID{:f}(1, 1, 2), CompositeIndex(Index(2, FID{:f}(1, 1, 1)), [0.0, 0.0], [0.0, 0.0])) # never do this !!!
Operator(2, FID{:f}(1, 1, 2), CompositeIndex(Index(2, FID{:f}(1, 1, 1)), [0.0, 0.0], [0.0, 0.0]))
```

[`Operator`](@ref) can be iterated and indexed by an integer, which will give the corresponding generator in the product:
```jldoctest OO
julia> op = Operator(2, FID{:f}(1, 1, 2), FID{:f}(1, 1, 1));

julia> length(op)
2

julia> [op[1], op[2]]
2-element Vector{FID{:f, Int64, Int64, Int64}}:
 FID{:f}(1, 1, 2)
 FID{:f}(1, 1, 1)

julia> collect(op)
2-element Vector{FID{:f, Int64, Int64, Int64}}:
 FID{:f}(1, 1, 2)
 FID{:f}(1, 1, 1)
```

To get the coefficient of an [`Operator`](@ref) or all its individual generators as a whole, use the [`value`](@ref) and [`id`](@ref) function exported by this package, respectively:
```jldoctest OO
julia> op = Operator(2, FID{:f}(1, 1, 2), FID{:f}(1, 1, 1));

julia> value(op)
2

julia> id(op)
(FID{:f}(1, 1, 2), FID{:f}(1, 1, 1))
```

The product between two [`Operator`](@ref)s, or the scalar multiplication between a number and an [`Operator`](@ref) is also an [`Operator`](@ref):
```jldoctest OO
julia> Operator(2, FID{:f}(1, 1, 2)) * Operator(3, FID{:f}(1, 1, 1))
Operator(6, FID{:f}(1, 1, 2), FID{:f}(1, 1, 1))

julia> 3 * Operator(2, FID{:f}(1, 1, 2))
Operator(6, FID{:f}(1, 1, 2))

julia> Operator(2, FID{:f}(1, 1, 2)) * 3
Operator(6, FID{:f}(1, 1, 2))
```

The Hermitian conjugate of an [`Operator`](@ref) can be obtained by the adjoint operator:
```jldoctest OO
julia> op = Operator(6, FID{:f}(2, 1, 2), FID{:f}(1, 1, 1));

julia> op'
Operator(6, FID{:f}(1, 1, 2), FID{:f}(2, 1, 1))
```

There also exists a special [`Operator`](@ref), which only has the coefficient:
```jldoctest OO
julia> Operator(2)
Operator(2)
```

[`Operators`](@ref) can be initialized by two ways:
```jldoctest OO
julia> Operators(Operator(2, FID{:f}(1, 1, 1)), Operator(3, FID{:f}(1, 1, 2)))
Operators with 2 Operator
  Operator(3, FID{:f}(1, 1, 2))
  Operator(2, FID{:f}(1, 1, 1))

julia> Operator(2, FID{:f}(1, 1, 1)) + Operator(3, FID{:f}(1, 1, 2)) - Operator(3, FID{:b}(1, 1, 2))
Operators with 3 Operator
  Operator(3, FID{:f}(1, 1, 2))
  Operator(2, FID{:f}(1, 1, 1))
  Operator(-3, FID{:b}(1, 1, 2))
```

Similar items are automatically merged during the construction of [`Operators`](@ref):
```jldoctest OO
julia> Operators(Operator(2, FID{:f}(1, 1, 1)), Operator(3, FID{:f}(1, 1, 1)))
Operators with 1 Operator
  Operator(5, FID{:f}(1, 1, 1))

julia> Operator(2, FID{:f}(1, 1, 1)) + Operator(3, FID{:f}(1, 1, 1))
Operators with 1 Operator
  Operator(5, FID{:f}(1, 1, 1))
```

The multiplication between two [`Operators`](@ref)es, or between an [`Operators`](@ref) and an [`Operator`](@ref), or between a number and a [`Operators`](@ref) are defined:
```jldoctest OO
julia> ops = Operator(2, FID{:f}(1, 1, 1)) + Operator(3, FID{:f}(1, 1, 2));

julia> op = Operator(2, FID{:f}(2, 1, 1));

julia> ops * op
Operators with 2 Operator
  Operator(6, FID{:f}(1, 1, 2), FID{:f}(2, 1, 1))
  Operator(4, FID{:f}(1, 1, 1), FID{:f}(2, 1, 1))

julia> op *  ops
Operators with 2 Operator
  Operator(6, FID{:f}(2, 1, 1), FID{:f}(1, 1, 2))
  Operator(4, FID{:f}(2, 1, 1), FID{:f}(1, 1, 1))

julia> another = Operator(2, FID{:f}(1, 1, 1)) + Operator(3, FID{:f}(1, 1, 2));

julia> ops * another 
Operators with 4 Operator
  Operator(6, FID{:f}(1, 1, 1), FID{:f}(1, 1, 2))
  Operator(4, FID{:f}(1, 1, 1), FID{:f}(1, 1, 1))
  Operator(9, FID{:f}(1, 1, 2), FID{:f}(1, 1, 2))
  Operator(6, FID{:f}(1, 1, 2), FID{:f}(1, 1, 1))

julia> 2 * ops
Operators with 2 Operator
  Operator(6, FID{:f}(1, 1, 2))
  Operator(4, FID{:f}(1, 1, 1))

julia> ops * 2
Operators with 2 Operator
  Operator(6, FID{:f}(1, 1, 2))
  Operator(4, FID{:f}(1, 1, 1))
```
It is noted that in the result, the distributive law automatically applies.

As is usual, the Hermitian conjugate of an [`Operators`](@ref) can be obtained by the adjoint operator:
```jldoctest
julia> ops = Operator(6, FID{:f}(1, 1, 2), FID{:f}(2, 1, 1)) + Operator(4, FID{:f}(1, 1, 1), FID{:f}(2, 1, 1));

julia> ops'
Operators with 2 Operator
  Operator(6, FID{:f}(2, 1, 2), FID{:f}(1, 1, 1))
  Operator(4, FID{:f}(2, 1, 2), FID{:f}(1, 1, 2))
```

[`Operators`](@ref) can be iterated, but cannot be indexed:
```jldoctest OO
julia> ops = Operator(2, FID{:f}(1, 1, 1)) + Operator(3, FID{:f}(1, 1, 2));

julia> collect(ops)
2-element Vector{Operator{Int64, Tuple{FID{:f, Int64, Int64, Int64}}}}:
 Operator(3, FID{:f}(1, 1, 2))
 Operator(2, FID{:f}(1, 1, 1))

julia> ops[1]
ERROR: MethodError: no method matching getindex(::Operators{Operator{Int64, Tuple{FID{:f, Int64, Int64, Int64}}}, Tuple{FID{:f, Int64, Int64, Int64}}}, ::Int64)
[...]
```
The different behaviors between [`Operator`](@ref) and [`Operators`](@ref) when they are indexed result from their underlying implementations: [`Operator`](@ref) is something like a tuple while [`Operators`](@ref) is something like a dictionary. Do not ask why not to implement them based on [expression trees](https://en.wikipedia.org/wiki/Binary_expression_tree). If you have to ask, then the answer will be the authors of the package don't know how. So if you are not satisfied with this implementation, feel free to post a pull request.