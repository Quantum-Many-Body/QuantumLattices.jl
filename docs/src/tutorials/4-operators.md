```@meta
CurrentModule = QuantumLattices
DocTestSetup = quote
    using QuantumLattices
end
```

```@setup operators
using QuantumLattices
using SymPy: symbols
```

# [4. Operator Algebra](@id TutorialOperators)

The previous two chapters provided the two ingredients out of which the Hamiltonian of a quantum lattice system is composed: the spatial structure of the lattice ([Chapter 2](@ref TutorialLattice)) and the generators acting on the local Hilbert space at each lattice point ([Chapter 3](@ref TutorialDOF)). This chapter turns to the operator representation of the Hamiltonian, which is built from these generators.

A Hamiltonian is a sum of operators, and each operator is a scalar coefficient multiplying a product of generators,

```math
t\,\hat{g}_{\mu_1}\hat{g}_{\mu_2}\cdots\hat{g}_{\mu_k},
```

where ``t`` is the scalar coefficient and each ``\hat{g}_{\mu_j}`` is a generator specified by the labels ``\mu_j``. Following the local-unitcell-global hierarchy of Section 3.1, these labels include its site together with its internal indices, supplemented by the coordinates when the generator lies outside the origin unitcell. In this package every generator is stored, together with its labels, as a single object of type [`OperatorIndex`](@ref), so that each component of the product above is such a whole generator rather than a bare set of labels. The Hamiltonian is then the sum of all such operators,

```math
H = \sum_j t_j\,\hat{g}_{\mu^{(j)}_1}\hat{g}_{\mu^{(j)}_2}\cdots\hat{g}_{\mu^{(j)}_{k_j}}.
```

In this way a single operator and the whole Hamiltonian are realized by two types: each operator of the form above corresponds to an [`Operator`](@ref), and their sum, i.e., the Hamiltonian itself, to an [`Operators`](@ref). Together they form an algebra over ``\mathbb{C}``: operators can be added, multiplied, and scaled by complex numbers, yielding other operators or their sums. This chapter introduces these concrete types, and their algebraic operations and transformations.

!!! note
    In this chapter the word "operator" denotes a single summand of the Hamiltonian, i.e., a scalar coefficient times a product of generators. This is different from the coupling *terms* to be introduced in [Chapter 5](@ref TutorialCouplings): a coupling term is a physical specification that expands into many operators.

## 4.1 Operator and Operators

[`Operator`](@ref) represents a product of generators with a coefficient. It can be initialized in two ways:
```jldoctest Op
julia> Operator(2, 𝕔⁺(1, -1//2), 𝕔(1, -1//2), 𝕊{1//2}('z'))
Operator(2, 𝕔⁺(1, -1//2), 𝕔(1, -1//2), 𝕊{1//2}('z'))

julia> 2 * 𝕔⁺(1, 1, -1//2) * 𝕊{1//2}(2, 'z')
Operator(2, 𝕔⁺(1, 1, -1//2), 𝕊{1//2}(2, 'z'))
```
Note that the number of generators can be any natural number.

Although generators at different levels can be multiplied to form an [`Operator`](@ref), it is not recommended to do so because the logic will become confusing:
```jldoctest Op
julia> Operator(2, 𝕔⁺(1, 0), 𝕔(2, 1, 0, [0.0], [0.0])) # never do this !!!
Operator(2, 𝕔⁺(1, 0), 𝕔(2, 1, 0, [0.0], [0.0]))
```

[`Operator`](@ref) can be iterated and indexed by integers, which will give the corresponding generators in the product:
```jldoctest Op
julia> op = Operator(2, 𝕔⁺(1, 1//2), 𝕔(1, 1//2));

julia> length(op)
2

julia> [op[1], op[2]]
2-element Vector{FockIndex{:f, Int64, Rational{Int64}}}:
 𝕔⁺(1, 1//2)
 𝕔(1, 1//2)

julia> collect(op)
2-element Vector{FockIndex{:f, Int64, Rational{Int64}}}:
 𝕔⁺(1, 1//2)
 𝕔(1, 1//2)
```

To get the coefficient of an [`Operator`](@ref) or all its individual generators as a whole, use the [`value`](@ref) and [`id`](@ref) function, respectively:
```jldoctest Op
julia> op = Operator(2, 𝕔⁺(1, 0), 𝕔(1, 0));

julia> value(op)
2

julia> id(op)
(𝕔⁺(1, 0), 𝕔(1, 0))
```

The product between two [`Operator`](@ref)s, or the scalar multiplication between a number and an [`Operator`](@ref) is also an [`Operator`](@ref):
```jldoctest Op
julia> Operator(2, 𝕔⁺(1, 1//2)) * Operator(3, 𝕔(1, 1//2))
Operator(6, 𝕔⁺(1, 1//2), 𝕔(1, 1//2))

julia> 3 * Operator(2, 𝕔⁺(1, 1//2))
Operator(6, 𝕔⁺(1, 1//2))

julia> Operator(2, 𝕔⁺(1, 1//2)) * 3
Operator(6, 𝕔⁺(1, 1//2))
```

The Hermitian conjugate of an [`Operator`](@ref) can be obtained by the adjoint operator:
```jldoctest Op
julia> op = Operator(6, 𝕔⁺(2, 1//2), 𝕔(1, 1//2));

julia> op'
Operator(6, 𝕔⁺(1, 1//2), 𝕔(2, 1//2))
```

There also exists a special [`Operator`](@ref), which only has the coefficient:
```jldoctest Op
julia> Operator(2)
Operator(2)
```

[`Operators`](@ref) represents a sum of [`Operator`](@ref) instances. It can be initialized in two ways:
```jldoctest Op
julia> Operators(Operator(2, 𝕔(1, 1//2)), Operator(3, 𝕔⁺(1, 1//2)))
Operators with 2 Operator
  Operator(2, 𝕔(1, 1//2))
  Operator(3, 𝕔⁺(1, 1//2))

julia> Operator(2, 𝕔(1, 1//2)) - Operator(3, 𝕒⁺(1, 1//2))
Operators with 2 Operator
  Operator(2, 𝕔(1, 1//2))
  Operator(-3, 𝕒⁺(1, 1//2))
```

Similar items are automatically merged during the construction of [`Operators`](@ref):
```jldoctest Op
julia> Operators(Operator(2, 𝕔(1, 1//2)), Operator(3, 𝕔(1, 1//2)))
Operators with 1 Operator
  Operator(5, 𝕔(1, 1//2))

julia> Operator(2, 𝕔(1, 1//2)) + Operator(3, 𝕔(1, 1//2))
Operators with 1 Operator
  Operator(5, 𝕔(1, 1//2))
```

The multiplication between two [`Operators`](@ref)es, or between an [`Operators`](@ref) and an [`Operator`](@ref), or between a number and an [`Operators`](@ref) is defined:
```jldoctest Op
julia> ops = Operator(2, 𝕔(1, 1//2)) + Operator(3, 𝕔⁺(1, 1//2));

julia> op = Operator(2, 𝕔(2, 1//2));

julia> ops * op
Operators with 2 Operator
  Operator(4, 𝕔(1, 1//2), 𝕔(2, 1//2))
  Operator(6, 𝕔⁺(1, 1//2), 𝕔(2, 1//2))

julia> op * ops
Operators with 2 Operator
  Operator(4, 𝕔(2, 1//2), 𝕔(1, 1//2))
  Operator(6, 𝕔(2, 1//2), 𝕔⁺(1, 1//2))

julia> another = Operator(2, 𝕔(1, 1//2)) + Operator(3, 𝕔⁺(1, 1//2));

julia> ops * another
Operators with 2 Operator
  Operator(6, 𝕔(1, 1//2), 𝕔⁺(1, 1//2))
  Operator(6, 𝕔⁺(1, 1//2), 𝕔(1, 1//2))

julia> 2 * ops
Operators with 2 Operator
  Operator(4, 𝕔(1, 1//2))
  Operator(6, 𝕔⁺(1, 1//2))

julia> ops * 2
Operators with 2 Operator
  Operator(4, 𝕔(1, 1//2))
  Operator(6, 𝕔⁺(1, 1//2))
```
Note that in the result, the distributive law automatically applies. Besides, the fermion operator relation ``c^2=(c^\dagger)^2=0`` is also used.

As is usual, the Hermitian conjugate of an [`Operators`](@ref) can be obtained by the adjoint operator:
```jldoctest OpHerm
julia> op₁ = Operator(6, 𝕔⁺(1, 1//2), 𝕔(2, 1//2));

julia> op₂ = Operator(4, 𝕔(1, 1//2), 𝕔(2, 1//2));

julia> ops = op₁ + op₂;

julia> ops'
Operators with 2 Operator
  Operator(6, 𝕔⁺(2, 1//2), 𝕔(1, 1//2))
  Operator(4, 𝕔⁺(2, 1//2), 𝕔⁺(1, 1//2))
```

[`Operators`](@ref) can be iterated and indexed:
```jldoctest Op
julia> ops = Operator(2, 𝕔(1, 1//2)) + Operator(3, 𝕔⁺(1, 1//2));

julia> collect(ops)
2-element Vector{Operator{Int64, Tuple{FockIndex{:f, Int64, Rational{Int64}}}}}:
 Operator(2, 𝕔(1, 1//2))
 Operator(3, 𝕔⁺(1, 1//2))

julia> ops[1]
Operator(2, 𝕔(1, 1//2))

julia> ops[2]
Operator(3, 𝕔⁺(1, 1//2))
```
The index order of an [`Operators`](@ref) is the insertion order of the operators it contains.

## 4.2 Symbolic Computation

When combined with [SymPy](https://github.com/JuliaPy/SymPy.jl), the operator algebra supports fully symbolic manipulation. The coefficient of any operator can be a SymPy symbol, and the algebraic operations of Section 4.1 preserve this symbolic nature:

```@example operators
# symbolic coefficients
t, U = symbols("t U", real=true)

# hopping and Hubbard operators
op_hop = Operator(t, 𝕔⁺(1, 1, 1//2), 𝕔(2, 1, 1//2))
op_hub = Operator(U, 𝕔⁺(1, 1, 1//2), 𝕔(1, 1, 1//2), 𝕔⁺(1, 1, -1//2), 𝕔(1, 1, -1//2))

# their sum keeps the coefficients symbolic
op_hop + op_hub
```

Here, the last result is displayed in LaTeX (see Section 4.3): the coefficients `t` and `U` remain as symbols throughout the addition. Only when you explicitly substitute numeric values does an operator become numeric. This separation of symbolic manipulation from numeric evaluation is a core design principle.

## 4.3 LaTeX Output

Operators are displayed in LaTeX format by default in Jupyter notebooks and in Documenter-generated documentation such as this one (but not in the plain REPL), as has been seen above. This is extremely useful for verifying generated Hamiltonians against hand-written expressions. The typesetting is recursive: to render an [`Operators`](@ref), each of its constituent [`Operator`](@ref)s is rendered in turn and the results are joined by ``+``/``-``; to render an [`Operator`](@ref), its coefficient is rendered first, followed by each of its generators in order. What remains to be specified is the LaTeX form of a single generator, and this is customizable: the form of each generator type is set manually through a [`LaTeX`](@ref) format stored per type and accessed through [`latexformat`](@ref). Such a format consists of a body, the literal symbol of the generator, together with two tuples of symbols that single out the attributes entering its typeset form:

```julia
LaTeX{SP, SB}(body, spdelimiter=",\\,", sbdelimiter=",\\,")
```

where `body` is the symbol of the operator (any literal works, e.g. the string `"c"` or the character `'c'` for a fermionic operator), the attributes named in `SP` form the **super**script, and those named in `SB` form the **sub**script, with the entries of each tuple joined by its own delimiter. For instance, the default format of the fermionic `Index{<:FockIndex{:f}}` type is

```julia
LaTeX{(:nambu,), (:site, :orbital, :spinsym)}("c")
```

That is, the `nambu` attribute is superscripted, e.g., a creation operator acquires a ``\dagger``, while the `site`, the `orbital`, and the spin symbol `spinsym` (``↓`` for `-1//2` and ``↑`` for `1//2`) are subscripted. The rendering of a concrete operator under this default format is:

```@example operators
# the default format typesets the spin as a symbol, i.e. ↓ for -1//2 and ↑ for 1//2
op = Operator(2, 𝕔⁺(1, 1, -1//2), 𝕔(2, 1, 1//2))
```

To customize the appearance of a given generator type, register a new [`LaTeX`](@ref) format with the two-argument [`latexformat`](@ref). The body, the superscript, and the subscript are set independently, so one can change the symbol of the generator, drop some of its labels, or alter how a particular label is typeset. For example, suppose the fermionic operator should be written with the body `"f"` instead of `"c"`, should not show its orbital, and should spell the spin out numerically rather than as a symbol; this is achieved by

```@example operators
# redefine the format of the Index{<:FockIndex{:f}} type
old = latexformat(Index{<:FockIndex{:f}})
latexformat(Index{<:FockIndex{:f}}, LaTeX{(:nambu,), (:site, :spin)}("f"))
op
```

The effect is threefold: the operator is now typeset with an `f`, the orbital label no longer appears in the subscript, and the spin is shown as a number (`-1//2` instead of ``\downarrow`` and `1//2` instead of ``\uparrow``). Since the registration is global, the format is restored afterwards so that the rest of the documentation is unaffected:

```@example operators
latexformat(Index{<:FockIndex{:f}}, old)
nothing # hide
```

Two further attributes can enter the typeset form, `rcoordinate` and `icoordinate`, but they are special in that they exist only for [`CoordinatedIndex`](@ref), the generators that carry the coordinates of their underlying point ([Section 3.2](@ref TutorialDOF)). Both are typeset as their bracketed value, and displaying them is what makes the generators on translationally equivalent bonds distinguishable. Adding them to the subscript of the coordinated fermionic index gives, for instance:

```@example operators
# rcoordinate and icoordinate exist only for CoordinatedIndex
oldcoord = latexformat(CoordinatedIndex{<:Index{<:FockIndex{:f}}})
latexformat(
    CoordinatedIndex{<:Index{<:FockIndex{:f}}},
    LaTeX{(:nambu,), (:site, :orbital, :spinsym, :rcoordinate, :icoordinate)}("c")
)
Operator(2, 𝕔⁺(1, 1, -1//2, [0.5, 0.0], [0.0, 0.0]), 𝕔(2, 1, 1//2, [1.5, 0.0], [1.0, 0.0]))
```

Naming them in the format of the [`Index`](@ref) type itself, by contrast, has no effect, because an [`Index`](@ref) carries no coordinates and the two entries are simply dropped from the subscript:

```@example operators
# naming them for the non-coordinated type has no effect
oldindex = latexformat(Index{<:FockIndex{:f}})
latexformat(
    Index{<:FockIndex{:f}},
    LaTeX{(:nambu,), (:site, :orbital, :spinsym, :rcoordinate, :icoordinate)}("c")
)
Operator(2, 𝕔⁺(1, 1, -1//2), 𝕔(2, 1, 1//2))
```

```@example operators
latexformat(Index{<:FockIndex{:f}}, oldindex)
latexformat(CoordinatedIndex{<:Index{<:FockIndex{:f}}}, oldcoord)
nothing # hide
```

## 4.4 Linear Transformation

A key feature of the operator algebra is that operators can be systematically transformed. The abstract type [`LinearTransformation`](@ref) represents any linear map on the operator algebra:

```math
T(t₁O₁ + t₂O₂) = t₁T(O₁) + t₂T(O₂)
```

All subtypes implement the callable interface `(transformation::LinearTransformation)(operator::Operator)`, making them composable and uniformly applicable to [`Operators`](@ref), and other containers. For a simple illustration, one can wrap an arbitrary function as a [`LinearTransformation`](@ref) via [`LinearFunction`](@ref). Linearity means that applying `f` to each operator of a sum yields the same result as applying it to the sum directly, so acting on an [`Operators`](@ref) simply transforms each constituent [`Operator`](@ref):

```jldoctest LT
julia> f = LinearFunction(op -> 2 * op);

julia> op = Operator(3, 𝕔(1, 1, 1//2));

julia> f(op)
Operator(6, 𝕔(1, 1, 1//2))

julia> ops = Operator(3, 𝕔(1, 1, 1//2)) + Operator(5, 𝕔⁺(1, 1, 1//2));

julia> f(ops)
Operators with 2 Operator
  Operator(6, 𝕔(1, 1, 1//2))
  Operator(10, 𝕔⁺(1, 1, 1//2))
```

## 4.5 Permutation

[`Permutation`](@ref) is a concrete [`LinearTransformation`](@ref) that reorders the generators in an operator product according to a reference [`Table`](@ref) (introduced in [Section 3.5](@ref TutorialTableAndMetric)).

### 4.5.1 Commutation Relations

When two adjacent generators are swapped, the result depends on the algebraic relations of the system, which differ from one category of system to another. In each case below, ``\alpha`` stands for the labels carried by a generator, and the function [`permute`](@ref) returns the outcome of exchanging the two generators.

#### **Fermions**

The generators of a fermionic algebra anticommute,

```math
\{c_\alpha, c_\beta\} = 0, \qquad \{c^\dagger_\alpha, c^\dagger_\beta\} = 0, \qquad \{c_\alpha, c^\dagger_\beta\} = \delta_{\alpha\beta}.
```

Whenever the two generators carry *different* labels, the swap merely introduces a factor of ``-1``, regardless of whether they are two annihilations, two creations, or one of each. The only non-trivial case is an annihilation and a creation with the *same* label, where ``c_\alpha c^\dagger_\alpha = 1 - c^\dagger_\alpha c_\alpha``: the swap then produces the constant operator ``1`` alongside the reordered ``-c^\dagger_\alpha c_\alpha``. Two identical annihilations or creations would instead vanish, since ``c_\alpha^2 = (c^\dagger_\alpha)^2 = 0``:

```jldoctest Permutation
julia> permute(𝕔(1, 1//2), 𝕔(1, -1//2))
(Operator(-1, 𝕔(1, -1//2), 𝕔(1, 1//2)),)

julia> permute(𝕔(1, 1//2), 𝕔⁺(1, -1//2))
(Operator(-1, 𝕔⁺(1, -1//2), 𝕔(1, 1//2)),)

julia> permute(𝕔(1, 1//2), 𝕔⁺(1, 1//2))
(Operator(1), Operator(-1, 𝕔⁺(1, 1//2), 𝕔(1, 1//2)))

julia> permute(𝕔⁺(1, 1//2), 𝕔⁺(1, -1//2))
(Operator(-1, 𝕔⁺(1, -1//2), 𝕔⁺(1, 1//2)),)
```

#### **Bosons**

The generators of a bosonic algebra commute,

```math
[a_\alpha, a_\beta] = 0, \qquad [a^\dagger_\alpha, a^\dagger_\beta] = 0, \qquad [a_\alpha, a^\dagger_\beta] = \delta_{\alpha\beta}.
```

Generators with *different* labels therefore swap without any sign change. As for fermions, the sole non-trivial case is an annihilation and a creation with the *same* label, where ``a_\alpha a^\dagger_\alpha = 1 + a^\dagger_\alpha a_\alpha``: the swap produces the constant operator ``1`` together with the reordered ``a^\dagger_\alpha a_\alpha``:

```jldoctest Permutation
julia> permute(𝕒(1, 0), 𝕒(1, 1))
(Operator(1, 𝕒(1, 1), 𝕒(1, 0)),)

julia> permute(𝕒(1, 0), 𝕒⁺(1, 1))
(Operator(1, 𝕒⁺(1, 1), 𝕒(1, 0)),)

julia> permute(𝕒(1, 0), 𝕒⁺(1, 0))
(Operator(1), Operator(1, 𝕒⁺(1, 0), 𝕒(1, 0)))

julia> permute(𝕒⁺(1, 0), 𝕒⁺(1, 1))
(Operator(1, 𝕒⁺(1, 1), 𝕒⁺(1, 0)),)
```

#### **Spins**

On the same lattice site the spin operators obey the SU(2) algebra,

```math
[S^\alpha, S^\beta] = i\varepsilon^{\alpha\beta\gamma} S^\gamma,
```

so that swapping two *different* components of the same site produces the third one, with the sign fixed by the orientation of the pair, ``S^\alpha S^\beta = S^\beta S^\alpha + i\varepsilon^{\alpha\beta\gamma} S^\gamma``. Two identical components commute, and so do components at different sites; in both cases the swap merely reorders:

```jldoctest Permutation
julia> permute(𝕊{1//2}('x'), 𝕊{1//2}('y'))
(Operator(1im, 𝕊{1//2}('z')), Operator(1, 𝕊{1//2}('y'), 𝕊{1//2}('x')))

julia> permute(𝕊{1//2}('y'), 𝕊{1//2}('x'))
(Operator(-1im, 𝕊{1//2}('z')), Operator(1, 𝕊{1//2}('x'), 𝕊{1//2}('y')))

julia> permute(𝕊{1//2}('x'), 𝕊{1//2}('x'))
(Operator(1, 𝕊{1//2}('x'), 𝕊{1//2}('x')),)

julia> permute(𝕊{1//2}(1, 'x'), 𝕊{1//2}(2, 'y'))
(Operator(1, 𝕊{1//2}(2, 'y'), 𝕊{1//2}(1, 'x')),)
```

#### **Phonons**

The displacement and momentum operators obey

```math
[u^\mu, p^\nu] = i\delta^{\mu\nu}.
```

The outcome depends on the two generators being the same kind or not, and on their directions. When they are a ``u`` and a ``p`` of the *same* direction, swapping them produces the constant operator ``\pm i`` together with the reordered one; all other pairs (two generators of the same kind, or generators of different directions) commute, and the swap merely reorders them:

```jldoctest Permutation
julia> permute(𝕦('x'), 𝕡('x'))
(Operator(1im), Operator(1, 𝕡('x'), 𝕦('x')))

julia> permute(𝕡('x'), 𝕦('x'))
(Operator(-1im), Operator(1, 𝕦('x'), 𝕡('x')))

julia> permute(𝕦('x'), 𝕦('y'))
(Operator(1, 𝕦('y'), 𝕦('x')),)

julia> permute(𝕦('x'), 𝕡('y'))
(Operator(1, 𝕡('y'), 𝕦('x')),)
```

### 4.5.2 Permutation in Action

`Permutation(table)` holds a [`Table`](@ref) as the target ordering. When called on an operator, it performs a bubble-sort-like scan: find the first adjacent pair that is out of order (according to the Table sequences), call [`permute`](@ref) on that pair, and continue on each resulting operator until all are ordered. Applied to an [`Operators`](@ref), the permutation acts on each constituent [`Operator`](@ref) independently. By default the product is arranged into the *descending* order of the Table sequences; the keyword `rev` requests the ascending order instead. To keep the prescribed order transparent, we build the [`Table`](@ref) explicitly from the generators that occur, together with a chosen [`Metric`](@ref), and the [`Table`](@ref) must contain every generator appearing in the operators to be permuted.

#### **Fermions**

Suppose we order the fermionic generators by the fields `(site, orbital, spin, nambu)`. All generators below share the same orbital ``1`` and spin ``1//2``, so they are distinguished by site and nambu alone. Consider a sum of two operators: the first one is an annihilation followed by a creation of the *same* label (same site, orbital, and spin), of sequences ``1`` and ``2``; the second one is a hopping-like product ``c_1 c^\dagger_2`` of sequences ``1`` and ``3``:

```jldoctest Permutation
julia> metric = OperatorIndexToTuple(:site, :orbital, :spin, :nambu);

julia> indexes = [𝕔(1, 1, 1//2), 𝕔⁺(1, 1, 1//2), 𝕔⁺(2, 1, 1//2)];

julia> table = Table(indexes, metric);

julia> table[𝕔(1, 1, 1//2)]
1

julia> table[𝕔⁺(1, 1, 1//2)]
2

julia> table[𝕔⁺(2, 1, 1//2)]
3

julia> P = Permutation(table);

julia> ops = Operator(2, 𝕔(1, 1, 1//2), 𝕔⁺(1, 1, 1//2)) + Operator(3, 𝕔(1, 1, 1//2), 𝕔⁺(2, 1, 1//2));

julia> P(ops)
Operators with 3 Operator
  Operator(-2, 𝕔⁺(1, 1, 1//2), 𝕔(1, 1, 1//2))
  Operator(2)
  Operator(-3, 𝕔⁺(2, 1, 1//2), 𝕔(1, 1, 1//2))
```

Each operator is permuted independently: reordering the same-label pair ``c_1 c^\dagger_1`` into ``c^\dagger_1 c_1`` uses ``c c^\dagger = 1 - c^\dagger c``, yielding the reordered ``-2\,c^\dagger_1 c_1`` together with the constant ``2``; the hopping-like operator ``c_1 c^\dagger_2`` is merely exchanged into ``-c^\dagger_2 c_1``.

#### **Bosons**

Suppose we order the bosonic generators by the fields `(site, orbital, spin, nambu)`, with all spins set to ``0``. The first one is again a same-label pair, an annihilation followed by a creation of the same label (sequences ``1`` and ``2``), while the second one is the hopping-like ``a_1 a^\dagger_2`` of sequences ``1`` and ``3``:

```jldoctest Permutation
julia> metric = OperatorIndexToTuple(:site, :orbital, :spin, :nambu);

julia> indexes = [𝕒(1, 1, 0), 𝕒⁺(1, 1, 0), 𝕒⁺(2, 1, 0)];

julia> table = Table(indexes, metric);

julia> table[𝕒(1, 1, 0)]
1

julia> table[𝕒⁺(1, 1, 0)]
2

julia> table[𝕒⁺(2, 1, 0)]
3

julia> P = Permutation(table);

julia> ops = Operator(2, 𝕒(1, 1, 0), 𝕒⁺(1, 1, 0)) + Operator(3, 𝕒(1, 1, 0), 𝕒⁺(2, 1, 0));

julia> P(ops)
Operators with 3 Operator
  Operator(2, 𝕒⁺(1, 1, 0), 𝕒(1, 1, 0))
  Operator(2)
  Operator(3, 𝕒⁺(2, 1, 0), 𝕒(1, 1, 0))
```

Both operators are permuted independently. The same-label pair is reordered through ``a a^\dagger = 1 + a^\dagger a``, giving the reordered ``2\,a^\dagger_1 a_1`` together with the constant ``2``; the hopping-like operator ``a_1 a^\dagger_2`` is exchanged into ``a^\dagger_2 a_1`` without any sign change, as befits bosonic statistics.

#### **Spins**

Suppose we order the spin generators by `(site, tag)`. Let's take a sum whose first operator is the same-site product ``S^x S^y``, of sequences ``1`` and ``2``, together with an ``S^x`` at another site. Reordering ``S^x S^y`` into ``S^y S^x`` now uses the SU(2) relation ``S^x S^y = S^y S^x + \mathrm{i} S^z``, so the result acquires an extra ``\mathrm{i} S^z`` operator. A complex coefficient is used for the operators so that this complex operator can be accumulated:

```jldoctest Permutation
julia> metric = OperatorIndexToTuple(:site, :tag);

julia> indexes = [𝕊{1//2}(1, 'x'), 𝕊{1//2}(1, 'y'), 𝕊{1//2}(1, 'z'), 𝕊{1//2}(2, 'x')];

julia> table = Table(indexes, metric);

julia> table[𝕊{1//2}(1, 'x')]
1

julia> table[𝕊{1//2}(1, 'y')]
2

julia> P = Permutation(table);

julia> ops = Complex(2)*𝕊{1//2}(1, 'x')*𝕊{1//2}(1, 'y') + Complex(1)*𝕊{1//2}(2, 'x');

julia> P(ops)
Operators with 3 Operator
  Operator(2, 𝕊{1//2}(1, 'y'), 𝕊{1//2}(1, 'x'))
  Operator(2im, 𝕊{1//2}(1, 'z'))
  Operator(1, 𝕊{1//2}(2, 'x'))
```

#### **Phonons**

Suppose we order the phononic generators by `(site, direction)` and, to distinguish a ``u`` from a ``p`` of the same site and direction, include the trait function [`kind`](@ref) in the metric. A ``p\,u`` product of the same site and direction then lies in the ascending order and must be exchanged into ``u\,p``, which through ``[u, p] = \mathrm{i}`` produces a constant operator ``-\mathrm{i}``:

```jldoctest Permutation
julia> metric = OperatorIndexToTuple(kind, :site, :direction);

julia> indexes = [𝕡(1, 'x'), 𝕦(1, 'x'), 𝕦(2, 'x')];

julia> table = Table(indexes, metric);

julia> table[𝕡(1, 'x')]
1

julia> table[𝕦(1, 'x')]
2

julia> P = Permutation(table);

julia> ops = Operator(Complex(3), 𝕡(1, 'x'), 𝕦(1, 'x')) + Operator(Complex(1), 𝕦(2, 'x'));

julia> P(ops)
Operators with 3 Operator
  Operator(3, 𝕦(1, 'x'), 𝕡(1, 'x'))
  Operator(-3im)
  Operator(1, 𝕦(2, 'x'))
```

## Summary

| Type | Represents |
|------|-----------|
| [`Operator`](@ref) | Scalar × product of local generators |
| [`Operators`](@ref) | Sum of `Operator` instances |
| [`LinearTransformation`](@ref) | Linear map on the operator algebra |
| [`Permutation`](@ref) | Reordering generators according to a reference Table |

Now that we understand the operator algebra, the next step is to learn how these operators are generated automatically from physical coupling terms ([Chapter 5](@ref TutorialCouplings)).
