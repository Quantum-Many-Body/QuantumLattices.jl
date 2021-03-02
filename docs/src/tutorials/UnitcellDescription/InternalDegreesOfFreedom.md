```@meta
CurrentModule = QuantumLattices
```

# Internal degrees of freedom

Now let's move to the second step. A main feature of a quantum lattice system is that one can decompose the whole Hilbert space into local ones. Correspondingly, the Hamiltonian of the system can always be expressed by products/sums of local generators of the algebras acting on these local Hilbert spaces. Thus, a parallel to the spatial info of the unitcell of the lattice system is the internal degrees of freedom that "lives" in each spatial point.

## Internal, IID and Index

Different kinds of quantum lattice systems can have different species of local Hilbert spaces. For example, the local Hilbert space of a complex fermionic/bosonic system is the [Fock space](https://en.wikipedia.org/wiki/Fock_space) whereas that of a spin-1/2 system is the two-dimensional ``\{\lvert\uparrow\rangle, \lvert\downarrow\rangle\}`` space. The abstract type [`Internal`](@ref) is an abstraction of local Hilbert spaces. More precisely, it is an abstraction of local algebras acting on these local Hilbert spaces. To specify a generator of such a local algebra, two sets of tags are needed: one identifies which the local algebra is and the other represents which the internal degree of freedom is. The former is already encoded by a [`PID`](@ref) object, and the latter will be stored in an object of a concrete subtype of the abstract type [`IID`](@ref). These two sets of tags are combined to be the [`Index`](@ref) type. Since different [`IID`](@ref)s can have different tags, [`Index`](@ref) is also an abstract type. For every concrete [`IID`](@ref) subtype, there corresponds a concrete [`Index`](@ref) subtype. In short, a concrete [`Internal`](@ref) object generates all concrete [`IID`](@ref) objects, which combined with the spatial tags provided by a [`PID`](@ref) object to be concrete [`Index`](@ref) objects, define the complete indices needed to specify the generators of the local algebras to construct the Hamiltonian.

In this package, we implement two groups of concrete subtypes to handle with the following two sets of systems, respectively:
* Canonical fermionic, canonical bosonic and hard-core bosonic systems,
* SU(2) spin systems.

### Fock, FID and FIndex

```@setup FFF
push!(LOAD_PATH, "../../../src/")
using QuantumLattices
```

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
It is the complete set of tags to specify a generator of local algebras of such systems.

Now let's see some examples.

An [`FID`](@ref) instance can be initialized by giving all its three attributes:
```@repl FFF
FID(2, 2, 1)
```
Or you can only specify part of the three by key word arguments and the others will get default values:
```@repl FFF
FID() # default values: orbital=1, spin=1, nambu=1
FID(spin=2, nambu=0)
```
The adjoint of an [`FID`](@ref) instance is also defined:
```@repl FFF
FID(orbital=3, spin=4, nambu=1)'
FID(orbital=3, spin=4, nambu=2)'
FID(orbital=3, spin=4, nambu=0)'
```
Apparently, this operation is nothing but the "Hermitian conjugate".

An [`FIndex`](@ref) instance can be initialized as usual by giving all its attributes:
```@repl FFF
FIndex('a', 1, 3, 4, 0)
```
Or it can be initialized by a [`PID`](@ref) instance and an [`FID`](@ref) instance:
```@repl FFF
FIndex(PID('a', 1), FID(3, 4, 0))
```
Conversely, the corresponding [`PID`](@ref) instance and [`FID`](@ref) instance can be extracted from an [`FIndex`](@ref) instance by the [`pid`](@ref) and [`iid`](@ref) method, respectively:
```@repl FFF
fidx = FIndex(PID('a', 1), FID(3, 4, 0));
fidx |> pid
fidx |> iid
```
Similar to [`FID`](@ref), the adjoint of an [`FIndex`](@ref) instance is also defined:
```@repl FFF
FIndex(PID('a', 1), FID(3, 4, 1))'
FIndex(PID('a', 1), FID(3, 4, 2))'
FIndex(PID('a', 1), FID(3, 4, 0))'
```

A [`Fock`](@ref) instance can be initialized by giving all its attributes or by key word arguments to specify those beyond the default values:
```@repl FFF
Fock(1, 1, 2, 2)
Fock() # default values: atom=1, norbital=1, nspin=2, nnambu=2
Fock(atom=2, norbital=2, nspin=1, nnambu=1)
```
As can be seen, a [`Fock`](@ref) instance behaves like a vector (because the parent type [`Internal`](@ref) is a subtype of `AbstractVector`), and its iteration just generates all the allowed [`FID`](@ref) instances on its associated spatial point:
```@repl FFF
fck = Fock(atom=1, norbital=2, nspin=1, nnambu=2);
fck |> typeof |> eltype
fck |> length
[fck[1], fck[2], fck[3], fck[4]]
fck |> collect
```

### Spin, SID and SIndex

```@setup SSS
push!(LOAD_PATH, "../../../src/")
using QuantumLattices
```

This group of concrete subtypes are designed to deal with SU(2) spin systems.

Although spin systems are essentially bosonic, the local Hilbert space often used is distinct from that of a canonical bosonic system: it is the space spanned by the eigenstates of a local ``S^z`` operator rather than a "[Fock space](https://en.wikipedia.org/wiki/Fock_space)". At the same time, a spin Hamiltonian is usually expressed by local spin operators, such as ``S^x``, ``S^y``, ``S^z``, ``S^+`` and ``S^-``, instead of creation and annihilation operators. Therefore, it is convenient to define another set of concrete subtypes for spin systems.

Apart from the spatial indices, a local spin operator can have an orbital index. It is also necessary to know what the total spin is. Besides, it needs an index to indicate which spin operator it is. Thus, the type [`SID`](@ref), which specifies the internal degree of freedom of a local spin operator, has the following attributes:
* `orbital::Int`: the orbital index
* `spin::Float64`: the total spin
* `tag::Char`: the tag, which must be `'x'`, `'y'`, `'z'`, `'+'` or `'-'`.
Correspondingly, the type [`Spin`](@ref), which defines the whole internal structure of a local spin space, has the following attributes:
* `atom::Int`: the atom index associated with a local spin space
* `norbital::Int`: the number of allowed orbital indices
* `spin::Float64`: the total spin
Similarly, the type [`SIndex`](@ref), which encapsulates all indices needed to specify a local spin operator, is just a combination of the indices in [`PID`](@ref) and [`SID`](@ref). Its attributes are as follows:
* `scope::Any`: the scope of a point
* `site::Int`: the site index of a point
* `orbital::Int`: the orbital index
* `spin::Float64`: the total spin
* `tag::Char`: the tag, which must be `'x'`, `'y'`, `'z'`, `'+'` or `'-'`.

Now let's see examples.

An [`SID`](@ref) instance can be initialized as follows
```@repl SSS
SID(1, 1.5, 'x')
SID() # default values: orbital=1, spin=0.5, tag='z'
SID(orbital=2, spin=1.0, tag='+')
```
The "Hermitian conjugate" of an [`SID`](@ref) instance can be obtained by the adjoint operation:
```@repl SSS
SID(1, 1.5, 'x')'
SID(1, 1.5, 'y')'
SID(1, 1.5, 'z')'
SID(1, 1.5, '+')'
SID(1, 1.5, '-')'
```
Regardless of the local orbital space, the local spin space is determined by the total spin. The standard matrix representation of an [`SID`](@ref) instance on this local spin space can be obtained by the `matrix` function:
```@repl SSS
SID(1, 0.5, 'x')' |> matrix
SID(1, 0.5, 'y')' |> matrix
SID(1, 0.5, 'z')' |> matrix
SID(1, 0.5, '+')' |> matrix
SID(1, 0.5, '-')' |> matrix
```

The usage of [`SIndex`](@ref) is completely parallel to that of [`FIndex`](@ref):
```@repl SSS
sidx = SIndex('a', 1, 3, 0.5, '-')
sidx == SIndex(PID('a', 1), SID(3, 0.5, '-'))
sidx |> pid
sidx |> iid
sidx'
```

A [`Spin`](@ref) instance can be initialized as follows:
```@repl SSS
Spin(1, 2, 1.0)
Spin() # default values: atom=1, orbital=1, spin=0.5
Spin(atom=2, norbital=1, spin=1.0)
```
Similar to [`Fock`](@ref), a [`Spin`](@ref) instance behaves like a vector whose iteration generates all the allowed [`SID`](@ref) instances on its associated spatial point:
```@repl FFF
sp = Spin(atom=1, norbital=1, spin=1.0);
sp |> typeof |> eltype
sp |> length
[sp[1], sp[2], sp[3], sp[4], sp[5]]
sp |> collect
```
It is noted that a [`Spin`](@ref) instance generates [`SID`](@ref) instances corresponding to not only ``S^x``, ``S^y`` and ``S^z``,  but also ``S^+`` and ``S^-`` although the former three already forms a complete set of generators of local spin algebras. This overcomplete feature is for the convenience to the construction of spin Hamiltonians.

## Config and Table

## Operator and Operators

### OID

### FOperator and BOperator

### latex-formatted output
