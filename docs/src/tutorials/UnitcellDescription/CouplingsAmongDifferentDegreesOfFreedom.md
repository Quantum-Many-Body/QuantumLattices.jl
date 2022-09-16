```@meta
CurrentModule = QuantumLattices
DocTestFilters = [r"im +[-\+]0\.0[-\+]"]
DocTestSetup = quote
    push!(LOAD_PATH, "../../../src/")
    using QuantumLattices
    using SymPy: symbols
end
```

# Couplings among different degrees of freedom

Now we arrive at the final step toward the complete description of a quantum lattice system, i.e., the terms that represent the couplings among different degrees of freedom.

## Ingredients of terms in Hamiltonians

In this package, the type [`Term`](@ref) is the representation of a term in lattice Hamiltonians.

As is well-known, different quantum lattice models have different terms. For example, the [Hubbard model](https://en.wikipedia.org/wiki/Hubbard_model) consist of a hopping term and a Hubbard term while the [transverse-field Ising model](https://en.wikipedia.org/wiki/Transverse-field_Ising_model) contains an Ising term and a transverse-field term. Despite the rich diversity of the terms in quantum lattice models, they host common ingredients as well:

1) **Overall coefficient**: every term has an overall coefficient, e.g., the hopping amplitude for the hopping term, the Hubbard interaction strength for the Hubbard term, etc.
2) **Kind of bonds to be summed over**: as the natural result of the lattice symmetry, every term contains a summation over some kind of generic bonds, e.g., the hopping term sums over a certain order of nearest-neighbor bonds, the Hubbard term sums over all individual points (namely the 1-point bonds), etc.
3) **Coupling pattern of local generators in the summation body**: apart from the site indexes, every term contains in the summation body over bonds a coupling pattern that can be represented by a certain combination of local generators, e.g., in the summation body with the site indexes dropped, the hopping term can be represented by $c^†c$, the Hubbard term can be represented by $c^†_↑ c_↑ c^†_↓c_↓$, the [XXX Heisenberg term](https://en.wikipedia.org/wiki/Quantum_Heisenberg_model) can be represented by $S^xS^x+S^yS^y+S^zS^z$ or $\frac{1}{2}S^+S^-+\frac{1}{2}S^-S^++S^zS^z$, etc.
4) **Site structure merging with local generators**: for every kind of term, the site indexes in the summation body over bonds has a definite structure, e.g., in the hopping term the $c^†$ operator is always associated with the first site of a 2-point bond while the $c$ operator is always associated with the second, in the Hubbard term all the four operators are always associated with the same site, etc.
5) **Hermiticity**: to guarantee the Hamiltonian to be Hermitian, the Hermitian conjugate of non-Hermitian terms must be added, e.g., the Hermitian conjugate of the hopping term must be added in the expression of the lattice Hamiltonian while that of the Hubbard term need not.
6) **Amplitude dependency on bonds** (optional): the amplitude of a term can be dependent on the generic bonds, e.g., the staggered local chemical potential depends on the sublattice index of a point, the $p+ip$ pairing potential depends on the direction of a bond, etc.

Such common ingredients determine the underlying organization of [`Term`](@ref). In fact, almost all of them manifest themselves in the basic construction function of [`Term`](@ref) shown as follows:
```julia
Term{termkind}(
    id::Symbol, coefficient, bondkind, coupling, ishermitian::Bool;
    amplitude::Union{Function, Nothing}=nothing
) where termkind
```
where `termkind` must be a `Symbol`, `coupling` specifies the coupling pattern of the term which can accept an instance of [`Coupling`](@ref), or an iterator of [`Coupling`](@ref)s, or a function that returns a [`Coupling`](@ref) or an iterator of [`Coupling`](@ref)s, and the keyword argument `amplitude` specifies the amplitude dependency on bonds of the term if it is not `nothing`. Here, the new type [`Coupling`](@ref) is the building block of the coupling pattern, which will be discussed in detail later. It is also noted that an extra id is also assigned with each term which can be used for fast lookup for later convenience.

### Site structures

The only ingredient that is not included in the above construction function of a [`Term`](@ref) is the 4th in the list, i.e., the site structure merging with local generators.

First of all, we want to comment on why the coupling pattern of a term does not include the site indexes but instead the site structure is needed:
* On the one hand, the lattice symmetry of the Hamiltonian guarantees the isomorphism of the coupling patterns of the local generators on symmetry-equivalent bonds with possibly different site indexes. Therefore, it is reasonable to exclude the site indexes to obtain a symmetry-equivalent coupling pattern.
* On the other hand, the dropped site indexes can be easily supplemented back as long as the site structure merging with local generators is provided, furthermore, such site structure is also invariant for all symmetry-equivalent bonds, making it more suitable for a term that sums over the symmetry-equivalent bonds.

Apparently, the site structure depends only on the type but not on the concrete instances of [`Term`](@ref). Therefore, it is desirable to define it separately outside the initialization of a term, which is implemented by the [`sitestructure`](@ref) function in this package:
```julia
sitestructure(::Val{termkind}, ::Val{termrank}, bondlength::Int) where {termkind, termrank} -> NTuple{termrank, Int}
```
Here, `termkind` is the kind of the term used for dispatch, `termrank` is the number of local generators in the coupling pattern, `bondlength` is the number of points contained in the bonds to be summed over of the term, and the returning tuple contains the sequences of the points that correspond to the local generators when the site indexes are supplemented back.

The methods aiming at the common used kinds of terms in condensed matter physics have been predefined by default in this package. You need not care about [`sitestructure`](@ref) unless an exception is encountered. Then you will have to fix the exception by simply define another method of [`sitestructure`](@ref) specialized for your own kind of [`Term`](@ref). For illustration, we provide an example of [`sitestructure`](@ref) as follows:
```julia
function sitestructure(::Val{:Hopping}, ::Val{2}, bondlength::Integer)
    @assert(
        bondlength==2,
        "sitestructure error: the bonds summed over for hopping terms must contain two points."
    )
    return (1, 2)
end
```

### Coupling patterns

Before the further discussion of [`Term`](@ref), we at first turn to the coupling pattern, which lies at the center of the construction of a [`Term`](@ref).

#### Coupling: building block of coupling patterns

[`Coupling`](@ref), as a combination of a set of local generators together with a coefficient, is the building block of the coupling pattern.

Let's see a typical example, which represents the coupling pattern $c^†c$ of an usual hopping term:
```jldoctest HM
julia> Coupling(1, FID(:, :, 2), FID(:, :, 1))
∑[FID(:, :, 2) FID(:, :, 1)]
```
Here the coefficient of the [`Coupling`](@ref) is $1$, which can be omitted in the construction function:
```jldoctest HM
julia> Coupling(FID(:, :, 2), FID(:, :, 1))
∑[FID(:, :, 2) FID(:, :, 1)]
```
It is noted that the [`FID`](@ref) instances used here are a bit different from those introduced in the previous section of [Internal degrees of freedom](@ref):
1) They don't have the type parameter `:f` or `:b` to specify the statistics of the generator, which means they apply to both fermionic and bosonic quantum lattice system.
2) Some of their fields, such as the `:orbital` and `:spin` attributes, are initialized by the `:` operator rather than an integer, which means such indexes are omitted and will be handled by the default mechanism that **omitted indexes are summed diagonally**, i.e., $c^†c≡\sum_{α, σ}c^†_{α, σ}c_{α, σ}$ (note there are no site indexes!). This rule is in fact a tradition in the literature of condensed matter physics.

The implicit summation in the construction of a [`Coupling`](@ref) is made explicit in its string representation by the "∑" symbol, as can be seen in above examples.

Of course, it also supports usual [`FID`](@ref)s to initialize more specific coupling patterns:
```jldoctest HM
julia> Coupling(-1, FID{:f}(1, 1, 2), FID{:f}(1, 1, 1))
- FID{:f}(1, 1, 2) FID{:f}(1, 1, 1)
```
which corresponds to $-c^†_{1, ↓}c_{1, ↓}$. In this coupling pattern, all indexes are definite, therefore, there is no summation symbol "∑" in the string representation.

The local generators can be of different types, which corresponds to a hybrid quantum lattice system that couples different categories of internal degrees of freedom:
```jldoctest HM
julia> Coupling(FID{:f}(1, 1, 2), FID{:f}(1, 1, 1), SID{1//2}('z'))
FID{:f}(1, 1, 2) FID{:f}(1, 1, 1) SID{1//2}('z')
```
Here, local spins are coupled to itinerant fermions. For more discussions on hybrid systems, please refer to the section of [Hybrid systems](@ref).

When all local generators are of the same type, a [`Coupling`](@ref) can be initialized in a different simpler way:
```jldoctest HM
julia> Coupling{FID}(:, :, (2, 2, 1, 1))
∑[FID(:, :, 2) FID(:, :, 2) FID(:, :, 1) FID(:, :, 1)]

julia> Coupling{FID}(2, :, :, (2, 2, 1, 1))
2 ∑[FID(:, :, 2) FID(:, :, 2) FID(:, :, 1) FID(:, :, 1)]

julia> Coupling{SID}(('z', 'z'))
SID('z') SID('z')

julia> Coupling{SID}(-1, ('z', 'z'))
- SID('z') SID('z')

julia> Coupling{PID}(('p', 'p'), :)
∑[PID('p', :) PID('p', :)]

julia> Coupling{PID}(1//2, ('p', 'p'), :)
1//2 ∑[PID('p', :) PID('p', :)]
```
Here, in the 3rd and 4th examples, the total spin of an [`SID`](@ref) is omitted, meaning it could represent any allowable value of total spins, as is similar to [`FID`](@ref) when the statistics is omitted.

The coefficient and the local generators of a [`Coupling`](@ref) are stored in the `:value` and `:iids` attributes, respectively:
```jldoctest
julia> coupling = Coupling(1//2, SID('+'), SID('-'));

julia> coupling.value
1//2

julia> coupling.iids
(SID('+'), SID('-'))
```

A [`Coupling`](@ref) can be multiplied with a number:
```jldoctest
julia> coupling = Coupling(1//2, SID('+'), SID('-'));

julia> coupling * 3
3//2 SID('+') SID('-')

julia> 3 * coupling
3//2 SID('+') SID('-')
```

Two [`Coupling`](@ref)s can be multiplied together:
```jldoctest
julia> cp₁ = Coupling{FID}(2, (1, 1), :, (2, 1));

julia> cp₂ = Coupling{FID}(2, (2, 2), :, (2, 1));

julia> cp₁ * cp₂
4 ∑[FID(1, :, 2) FID(1, :, 1)] ⋅ ∑[FID(2, :, 2) FID(2, :, 1)]
```

It is noted that due to the implicit summation in the coupling pattern, the above product is not equal to `Coupling{FID}(4, (1, 1, 2, 2), :, (2, 1, 2, 1))`:
```jldoctest
julia> cp₁, cp₂ = Coupling{FID}(2, (1, 1), :, (2, 1)), Coupling{FID}(2, (2, 2), :, (2, 1));

julia> cp = Coupling{FID}(4, (1, 1, 2, 2), :, (2, 1, 2, 1));

julia> cp == cp₁ * cp₂
false
```

#### Coupling pattern with constraints


#### MatrixCoupling: coupling pattern with a matrix acting on the internal (sub)space


#### Bond dependent coupling patterns


### Amplitude dependency on bonds


## Specialized terms

For each certain kind of terms, some of the input parameters of the basic construction function are in fact fixed or have default values, e.g., the hopping term is always non-Hermitian while the Hubbard term is always Hermitian. Therefore, for each common kind of terms in condensed matter physics, it is more convenient to define the specialized construction function. All such predefinitions can be referenced in the section of [Advanced terms](@ref). Here, we only introduce [`Hopping`](@ref) and [`Hubbard`](@ref) specialized for the hopping terms and Hubbard terms, respectively, for illustration:
```julia
Hopping(id::Symbol, value, bondkind; coupling=Coupling{FID}(:, :, (2, 1)))
Hubbard(id::Symbol, value)
```

It can be seen that the terms constructed by the generic [`Term`](@ref) function or the specialized [`Hopping`](@ref) or [`Hubbard`](@ref) functions are the same, but the specialized functions are much more convenient:
```jldoctest
julia> t₁ =  Term{:Hopping}(:t, 2.0, 1, Coupling{FID}(:, :, (2, 1)), false);

julia> t₂ = Hopping(:t, 2.0, 1);

julia> t₁ == t₂
true

julia> U₁ = Term{:Hubbard}(:U, 1.5, 0, Coupling{FID}(:, (2, 2, 1, 1), (2, 1, 2, 1)), true);

julia> U₂ = Hubbard(:U, 1.5);

julia> U₁ == U₂
true
```

It is highly recommended to use the specialized functions to construct a term unless such a specialization does not exist.

## Expand terms to obtain operators

To obtain the operators of a [`Term`](@ref), the [`expand`](@ref) function exported by this package can be used as follows:
```julia
expand(term::Term, bond::Bond, hilbert::Hilbert) -> Operators
expand(term::Term, bonds::Vector{<:Bond}, hilbert::Hilbert) -> Operators
```

Let's see some examples.

For the hopping term:
```jldoctest
julia> t = Hopping(:t, 2.0, 1);

julia> bond = Bond(1, Point(1, [0.0], [0.0]), Point(2, [0.5], [0.0]));

julia> hilbert = Hilbert(1=>Fock{:f}(1, 2), 2=>Fock{:f}(1, 2));

julia> expand(t, bond, hilbert)
Operators with 4 Operator
  Operator(2.0, CompositeIndex(Index(2, FID{:f}(1, 2, 2)), [0.5], [0.0]), CompositeIndex(Index(1, FID{:f}(1, 2, 1)), [0.0], [0.0]))
  Operator(2.0, CompositeIndex(Index(2, FID{:f}(1, 1, 2)), [0.5], [0.0]), CompositeIndex(Index(1, FID{:f}(1, 1, 1)), [0.0], [0.0]))
  Operator(2.0, CompositeIndex(Index(1, FID{:f}(1, 2, 2)), [0.0], [0.0]), CompositeIndex(Index(2, FID{:f}(1, 2, 1)), [0.5], [0.0]))
  Operator(2.0, CompositeIndex(Index(1, FID{:f}(1, 1, 2)), [0.0], [0.0]), CompositeIndex(Index(2, FID{:f}(1, 1, 1)), [0.5], [0.0]))
```

When a bond and a term do not match each other, the [`expand`](@ref) function will return an empty [`Operators`](@ref):
```jldoctest
julia> t = Hopping(:t, 1.0, 1);

julia> bond = Bond(2, Point(1, [0.0], [0.0]), Point(1, [1.0], [1.0]));

julia> hilbert = Hilbert(1=>Fock{:f}(1, 2), 2=>Fock{:f}(1, 2));

julia> expand(t, bond, hilbert)
Operators with 0 Operator
```

In the [`expand`](@ref) function, a `Vector` of [`Bond`](@ref)s can also be provided to get all the operators expanded on such bonds:
```jldoctest
julia> t = Hopping(:t, 1.0, 1);

julia> lattice = Lattice([0.0], [0.5]; vectors=[[1.0]]);

julia> hilbert = Hilbert(1=>Fock{:f}(1, 2), 2=>Fock{:f}(1, 2));

julia> expand(t, bonds(lattice, 1), hilbert)
Operators with 8 Operator
  Operator(1.0, CompositeIndex(Index(1, FID{:f}(1, 1, 2)), [0.0], [0.0]), CompositeIndex(Index(2, FID{:f}(1, 1, 1)), [-0.5], [-1.0]))
  Operator(1.0, CompositeIndex(Index(2, FID{:f}(1, 2, 2)), [0.5], [0.0]), CompositeIndex(Index(1, FID{:f}(1, 2, 1)), [0.0], [0.0]))
  Operator(1.0, CompositeIndex(Index(2, FID{:f}(1, 1, 2)), [0.5], [0.0]), CompositeIndex(Index(1, FID{:f}(1, 1, 1)), [0.0], [0.0]))
  Operator(1.0, CompositeIndex(Index(1, FID{:f}(1, 2, 2)), [0.0], [0.0]), CompositeIndex(Index(2, FID{:f}(1, 2, 1)), [-0.5], [-1.0]))
  Operator(1.0, CompositeIndex(Index(1, FID{:f}(1, 2, 2)), [0.0], [0.0]), CompositeIndex(Index(2, FID{:f}(1, 2, 1)), [0.5], [0.0]))
  Operator(1.0, CompositeIndex(Index(2, FID{:f}(1, 2, 2)), [-0.5], [-1.0]), CompositeIndex(Index(1, FID{:f}(1, 2, 1)), [0.0], [0.0]))
  Operator(1.0, CompositeIndex(Index(2, FID{:f}(1, 1, 2)), [-0.5], [-1.0]), CompositeIndex(Index(1, FID{:f}(1, 1, 1)), [0.0], [0.0]))
  Operator(1.0, CompositeIndex(Index(1, FID{:f}(1, 1, 2)), [0.0], [0.0]), CompositeIndex(Index(2, FID{:f}(1, 1, 1)), [0.5], [0.0]))
```
