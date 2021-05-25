```@meta
CurrentModule = QuantumLattices
DocTestFilters = [r"[ +-]0\.0[ +-]"]
DocTestSetup = quote
    push!(LOAD_PATH, "../../../src/")
    using QuantumLattices
end
```

# Internal degrees of freedom

Now let's move to the second step. A main feature of a quantum lattice system is that one can decompose the whole Hilbert space into local ones. Correspondingly, the Hamiltonian of the system can always be expressed by products/sums of **local generators of the algebras acting on these local Hilbert spaces**. Thus, a parallel to the spatial info of the unitcell of the lattice system is the internal degrees of freedom that "lives" in each spatial point.

## Internal, IID and Index

Different kinds of quantum lattice systems can have different species of local Hilbert spaces. For example, the local Hilbert space of a complex fermionic/bosonic system is the [Fock space](https://en.wikipedia.org/wiki/Fock_space) whereas that of a spin-1/2 system is the two-dimensional ``\{\lvert\uparrow\rangle, \lvert\downarrow\rangle\}`` space. The abstract type [`Internal`](@ref) is an abstraction of local Hilbert spaces. More precisely, it is an abstraction of **local algebras acting on these local Hilbert spaces**. To specify a generator of such a local algebra, two sets of tags are needed: one identifies which the local algebra is and the other represents which the internal degree of freedom is. The former is already encoded by a [`PID`](@ref) object, and the latter will be stored in an object of a concrete subtype of the abstract type [`IID`](@ref). These two sets of tags are combined to be the [`Index`](@ref) type, which contains all the tags of a local generator. Since different [`IID`](@ref)s can have different tags, [`Index`](@ref) is also an abstract type. For every concrete [`IID`](@ref) subtype, there corresponds a concrete [`Index`](@ref) subtype.

In this package, we implement two groups of concrete subtypes to handle with the following two sets of systems, respectively:
* Canonical fermionic, canonical bosonic and hard-core bosonic systems,
* SU(2) spin systems.

### Fock, FID and FIndex

This group of concrete subtypes are designed to deal with canonical fermionic, canonical bosonic and hard-core bosonic systems.

Roughly speaking, these systems share similar internal structures of local Hilbert spaces termed as the [Fock space](https://en.wikipedia.org/wiki/Fock_space) where the generators of local algebras are the annihilation and creation operators. Besides the spatial index, such an operator usually adopts an orbital index and a spin index. It also needs a nambu index to distinguish whether it is an annihilation one or a creation one. Thus, the type [`FID`](@ref)`<:IID`, which specifies a certain internal degree of freedom of a local [Fock space](https://en.wikipedia.org/wiki/Fock_space), has the following attributes:
* `orbital::Int`: the orbital index
* `spin::Int`: the spin index
* `nambu::Int`: the nambu index, which must be 0, 1(annihilation) or 2(creation).
Correspondingly, the type [`Fock`](@ref)`<:Internal`, which specifies the whole internal structure of a local [Fock space](https://en.wikipedia.org/wiki/Fock_space), has the following attributes:
* `atom::Int`: the atom index associated with a local Hilbert space
* `norbital::Int`: the number of allowed orbital indices
* `nspin::Int`: the number of allowed spin indices
* `nnambu::Int`: the number of allowed nambu indices, which must be either 1 or 2.
It is noted that we also associate an atom index with each [`Fock`](@ref) instance, which proves to be helpful in future usages.

One more remark. The `:nambu` attribute of an [`FID`](@ref) instance can be 1 or 2, which means it represents an annihilation operator or a creation operator, respectively. This corresponds to a usual complex fermionic/bosonic system. The `:nambu` attribute can also be 0. In this case, it corresponds to a real fermionic/bosonic system where annihilation and creation operators are identical to each other, e.g. a Majorana fermionic system. Accordingly, The `:nnambu` attribute of a [`Fock`](@ref) instance can be either 2 or 1. Being 2, it allows usual complex annihilation/creation operators, while being 1 it only allows real fermionic/bosonic operators.

The type [`FIndex`](@ref)`<:Index` gathers all the tags in [`PID`](@ref) and [`FID`](@ref), which apparently has the following attributes:
* `scope::Any`: the scope of a point
* `site::Int`: the site index of a point
* `orbital::Int`: the orbital index
* `spin::Int`: the spin index
* `nambu::Int`: the nambu index, which must be 0, 1(annihilation) or 2(creation).
It is the complete set of tags to specify a generator (the annihilation/creation operator) of local algebras of such systems.

To distinguish whether the system is a fermionic one or a bosonic one, [`FID`](@ref), [`Fock`](@ref) and [`FIndex`](@ref) take a symbol `:f`(for fermionic) or `:b`(for bosonic) to be their first type parameters.

Now let's see some examples.

An [`FID`](@ref) instance can be initialized by giving all its three attributes:
```jldoctest FFF
julia> FID{:f}(2, 2, 1)
FID{:f}(2, 2, 1)

julia> FID{:b}(2, 2, 1)
FID{:b}(2, 2, 1)
```
Or you can only specify part of the three by key word arguments and the others will get default values:
```jldoctest FFF
julia> FID{:f}() # default values: orbital=1, spin=1, nambu=1
FID{:f}(1, 1, 1)

julia> FID{:b}(spin=2, nambu=0)
FID{:b}(1, 2, 0)
```
The adjoint of an [`FID`](@ref) instance is also defined:
```jldoctest FFF
julia> FID{:f}(orbital=3, spin=4, nambu=1)'
FID{:f}(3, 4, 2)

julia> FID{:b}(orbital=3, spin=4, nambu=2)'
FID{:b}(3, 4, 1)

julia> FID{:b}(orbital=3, spin=4, nambu=0)'
FID{:b}(3, 4, 0)
```
Apparently, this operation is nothing but the "Hermitian conjugate".

An [`FIndex`](@ref) instance can be initialized as usual by giving all its attributes:
```jldoctest FFF
julia> FIndex{:f}('a', 1, 3, 4, 0)
FIndex{:f}('a', 1, 3, 4, 0)
```
Or it can be initialized by a [`PID`](@ref) instance and an [`FID`](@ref) instance:
```jldoctest FFF
julia> FIndex(PID('a', 1), FID{:b}(3, 4, 0))
FIndex{:b}('a', 1, 3, 4, 0)
```
Note at this time, the statistics of the system is omitted because it can be deduced from the input [`FID`](@ref). Conversely, the corresponding [`PID`](@ref) instance and [`FID`](@ref) instance can be extracted from an [`FIndex`](@ref) instance by the [`pid`](@ref) and [`iid`](@ref) method, respectively:
```jldoctest FFF
julia> fidx = FIndex(PID('a', 1), FID{:b}(3, 4, 0));

julia> fidx |> pid
PID('a', 1)

julia> fidx |> iid
FID{:b}(3, 4, 0)
```
Similar to [`FID`](@ref), the adjoint of an [`FIndex`](@ref) instance is also defined:
```jldoctest FFF
julia> FIndex(PID('a', 1), FID{:f}(3, 4, 1))'
FIndex{:f}('a', 1, 3, 4, 2)

julia> FIndex(PID('a', 1), FID{:f}(3, 4, 2))'
FIndex{:f}('a', 1, 3, 4, 1)

julia> FIndex(PID('a', 1), FID{:b}(3, 4, 0))'
FIndex{:b}('a', 1, 3, 4, 0)
```

A [`Fock`](@ref) instance can be initialized by giving all its attributes or by key word arguments to specify those beyond the default values:
```jldoctest FFF
julia> Fock{:f}(1, 1, 2, 2)
4-element Fock{:f}:
 FID{:f}(1, 1, 1)
 FID{:f}(1, 2, 1)
 FID{:f}(1, 1, 2)
 FID{:f}(1, 2, 2)

julia> Fock{:b}() # default values: atom=1, norbital=1, nspin=2, nnambu=2
4-element Fock{:b}:
 FID{:b}(1, 1, 1)
 FID{:b}(1, 2, 1)
 FID{:b}(1, 1, 2)
 FID{:b}(1, 2, 2)

julia> Fock{:b}(atom=2, norbital=2, nspin=1, nnambu=1)
2-element Fock{:b}:
 FID{:b}(1, 1, 0)
 FID{:b}(2, 1, 0)
```
As can be seen, a [`Fock`](@ref) instance behaves like a vector (because the parent type [`Internal`](@ref) is a subtype of `AbstractVector`), and its iteration just generates all the allowed [`FID`](@ref) instances on its associated spatial point:
```jldoctest FFF
julia> fck = Fock{:f}(atom=1, norbital=2, nspin=1, nnambu=2);

julia> fck |> typeof |> eltype
FID{:f}

julia> fck |> length
4

julia> [fck[1], fck[2], fck[3], fck[4]]
4-element Vector{FID{:f}}:
 FID{:f}(1, 1, 1)
 FID{:f}(2, 1, 1)
 FID{:f}(1, 1, 2)
 FID{:f}(2, 1, 2)

julia> fck |> collect
4-element Vector{FID{:f}}:
 FID{:f}(1, 1, 1)
 FID{:f}(2, 1, 1)
 FID{:f}(1, 1, 2)
 FID{:f}(2, 1, 2)
```

### Spin, SID and SIndex

```@setup SSS
push!(LOAD_PATH, "../../../src/")
using QuantumLattices
```

This group of concrete subtypes are designed to deal with SU(2) spin systems.

Although spin systems are essentially bosonic, the commonly-used local Hilbert space is distinct from that of a canonical bosonic system: it is the space spanned by the eigenstates of a local ``S^z`` operator rather than a "[Fock space](https://en.wikipedia.org/wiki/Fock_space)". At the same time, a spin Hamiltonian is usually expressed by local spin operators, such as ``S^x``, ``S^y``, ``S^z``, ``S^+`` and ``S^-``, instead of creation and annihilation operators. Therefore, it is convenient to define another set of concrete subtypes for spin systems.

Apart from the spatial indices, a local spin operator can have an orbital index. Besides, it also needs an index to indicate which spin operator it is. Thus, the type [`SID`](@ref), which specifies the internal degree of freedom of a local spin operator, has the following attributes:
* `orbital::Int`: the orbital index
* `tag::Char`: the tag, which must be `'x'`, `'y'`, `'z'`, `'+'` or `'-'`.
Correspondingly, the type [`Spin`](@ref), which defines the whole internal structure of a local spin space, has the following attributes:
* `atom::Int`: the atom index associated with a local spin space
* `norbital::Int`: the number of allowed orbital indices
Similarly, the type [`SIndex`](@ref), which encapsulates all indices needed to specify a local spin operator, is just a combination of the indices in [`PID`](@ref) and [`SID`](@ref). Its attributes are as follows:
* `scope::Any`: the scope of a point
* `site::Int`: the site index of a point
* `orbital::Int`: the orbital index
* `tag::Char`: the tag, which must be `'x'`, `'y'`, `'z'`, `'+'` or `'-'`.

For [`SID`](@ref), [`Spin`](@ref) and [`SIndex`](@ref), it is also necessary to know what the total spin is, which is taken as their first type parameters and should be a half-integer or an integer.

Now let's see examples.

An [`SID`](@ref) instance can be initialized as follows
```jldoctest SSS
julia> SID{3//2}(1, 'x')
SID{3//2}(1, 'x')

julia> SID{1//2}('z') # default values: orbital=1
SID{1//2}(1, 'z')

julia> SID{1}('+', orbital=2)
SID{1}(2, '+')
```
The "Hermitian conjugate" of an [`SID`](@ref) instance can be obtained by the adjoint operation:
```jldoctest SSS
julia> SID{3//2}(1, 'x')'
SID{3//2}(1, 'x')

julia> SID{3//2}(1, 'y')'
SID{3//2}(1, 'y')

julia> SID{3//2}(1,'z')'
SID{3//2}(1, 'z')

julia> SID{3//2}(1, '+')'
SID{3//2}(1, '-')

julia> SID{3//2}(1, '-')'
SID{3//2}(1, '+')
```
Regardless of the local orbital space, the local spin space is determined by the total spin. The standard matrix representation of an [`SID`](@ref) instance on this local spin space can be obtained by the `Matrix` function:
```jldoctest SSS
julia> SID{1//2}(1, 'x') |> Matrix
2×2 Matrix{ComplexF64}:
 0.0+0.0im  0.5+0.0im
 0.5+0.0im  0.0+0.0im

julia> SID{1//2}(1, 'y') |> Matrix
2×2 Matrix{ComplexF64}:
 0.0-0.0im  -0.0+0.5im
 0.0-0.5im   0.0-0.0im

julia> SID{1//2}(1, 'z') |> Matrix
2×2 Matrix{ComplexF64}:
 -0.5+0.0im  0.0+0.0im
  0.0+0.0im  0.5+0.0im

julia> SID{1//2}(1, '+') |> Matrix
2×2 Matrix{ComplexF64}:
 0.0+0.0im  0.0+0.0im
 1.0+0.0im  0.0+0.0im

julia> SID{1//2}(1, '-') |> Matrix
2×2 Matrix{ComplexF64}:
 0.0+0.0im  1.0+0.0im
 0.0+0.0im  0.0+0.0im
```

The usage of [`SIndex`](@ref) is completely parallel to that of [`FIndex`](@ref):
```jldoctest SSS
julia> sidx = SIndex{1//2}('a', 1, 3, '-')
SIndex{1//2}('a', 1, 3, '-')

julia> sidx == SIndex(PID('a', 1), SID{1//2}(3, '-'))
true

julia> sidx |> pid
PID('a', 1)

julia> sidx |> iid
SID{1//2}(3, '-')

julia> sidx'
SIndex{1//2}('a', 1, 3, '+')
```

A [`Spin`](@ref) instance can be initialized as follows:
```jldoctest SSS
julia> Spin{1}(1, 2)
10-element Spin{1}:
 SID{1}(1, 'x')
 SID{1}(2, 'x')
 SID{1}(1, 'y')
 SID{1}(2, 'y')
 SID{1}(1, 'z')
 SID{1}(2, 'z')
 SID{1}(1, '+')
 SID{1}(2, '+')
 SID{1}(1, '-')
 SID{1}(2, '-')

julia> Spin{1//2}() # default values: atom=1, norbital=1
5-element Spin{1//2}:
 SID{1//2}(1, 'x')
 SID{1//2}(1, 'y')
 SID{1//2}(1, 'z')
 SID{1//2}(1, '+')
 SID{1//2}(1, '-')

julia> Spin{1}(atom=2, norbital=1)
5-element Spin{1}:
 SID{1}(1, 'x')
 SID{1}(1, 'y')
 SID{1}(1, 'z')
 SID{1}(1, '+')
 SID{1}(1, '-')
```
Similar to [`Fock`](@ref), a [`Spin`](@ref) instance behaves like a vector whose iteration generates all the allowed [`SID`](@ref) instances on its associated spatial point:
```jldoctest SSS
julia> sp = Spin{1}(atom=1, norbital=1);

julia> sp |> typeof |> eltype
SID{1}

julia> sp |> length
5

julia> [sp[1], sp[2], sp[3], sp[4], sp[5]]
5-element Vector{SID{1}}:
 SID{1}(1, 'x')
 SID{1}(1, 'y')
 SID{1}(1, 'z')
 SID{1}(1, '+')
 SID{1}(1, '-')

julia> sp |> collect
5-element Vector{SID{1}}:
 SID{1}(1, 'x')
 SID{1}(1, 'y')
 SID{1}(1, 'z')
 SID{1}(1, '+')
 SID{1}(1, '-')
```
It is noted that a [`Spin`](@ref) instance generates [`SID`](@ref) instances corresponding to not only ``S^x``, ``S^y`` and ``S^z``,  but also ``S^+`` and ``S^-`` although the former three already forms a complete set of generators of local spin algebras. This overcomplete feature is for the convenience to the construction of spin Hamiltonians.

## Config

## Operator and Operators

### OID

### FOperator and BOperator

### latex-formatted output
