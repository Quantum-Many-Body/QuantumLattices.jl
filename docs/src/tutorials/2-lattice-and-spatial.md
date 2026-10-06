```@meta
CurrentModule = QuantumLattices
DocTestSetup = quote
    using QuantumLattices
    using StaticArrays: SVector
end
```

# [2. Lattice and Spatial Structure](@id TutorialLattice)

In the standard workflow to define a quantum lattice system, we must first specify where the sites are and how they are connected. This chapter introduces the spatial building blocks: the lattice, its points, and the bonds between them. Besides, the reciprocal-space types essential for momentum-space analysis are also covered.

## 2.1 Constructing a Lattice

### 2.1.1 Spatial information of a unitcell

A lattice in condensed matter physics is a set of points in $\mathbb{R}^d$ closed under translations by a discrete set of vectors. Translation symmetry introduces an equivalence relation: two points related by multiples of the translation vectors are physically equivalent. This observation leads to the **unitcell construction**: it suffices to specify a finite set of points $\{\mathbf{r}_1, \ldots, \mathbf{r}_n\}$ within the origin unitcell, together with a set of translation vectors $\{\mathbf{a}_1, \ldots, \mathbf{a}_d\}$.

[`Lattice`](@ref) is the concrete type representing the unitcell of an infinite lattice. It has three attributes:
* `name::Symbol`: the name of the lattice
* `coordinates::Matrix{<:Number}`: the coordinates of the points within the origin unitcell
* `vectors::SVector{N, <:SVector}`: the translation vectors of the lattice
Here, [`SVector`](https://juliaarrays.github.io/StaticArrays.jl/stable/api/#SVector) is an immutable vector defined and exported by the [StaticArrays](https://github.com/JuliaArrays/StaticArrays.jl) package.

[`Lattice`](@ref) can be constructed by providing the coordinates, with optional keyword arguments to specify its name and translation vectors:
```jldoctest unitcell
julia> Lattice([0.0])
Lattice(lattice)
  with 1 point:
    [0.0]

julia> Lattice((0.0, 0.0), (0.5, 0.5); vectors=[[1.0, 0.0], [0.0, 1.0]], name=:Square)
Lattice(Square)
  with 2 points:
    [0.0, 0.0]
    [0.5, 0.5]
  with 2 translation vectors:
    [1.0, 0.0]
    [0.0, 1.0]

julia> Lattice(
           (0.0, 0.0, 0.0);
           name=:Cube,
           vectors=[[1.0, 0.0, 0.0], [0.0, 1.0, 0.0], [0.0, 0.0, 1.0]]
       )
Lattice(Cube)
  with 1 point:
    [0.0, 0.0, 0.0]
  with 3 translation vectors:
    [1.0, 0.0, 0.0]
    [0.0, 1.0, 0.0]
    [0.0, 0.0, 1.0]
```
The coordinates can be specified using vectors or tuples. When the `vectors` keyword is omitted, as in the first example above, the lattice is a finite cluster: it possesses no translation symmetry, and the "unitcell" coincides with the whole system.

Iteration over a lattice yields the coordinates of the points in it:
```jldoctest unitcell
julia> lattice = Lattice((0.0, 0.0), (0.5, 0.5); vectors=[[1.0, 0.0], [0.0, 1.0]]);

julia> length(lattice)
2

julia> [lattice[1], lattice[2]]
2-element Vector{SVector{2, Float64}}:
 [0.0, 0.0]
 [0.5, 0.5]

julia> collect(lattice)
2-element Vector{SVector{2, Float64}}:
 [0.0, 0.0]
 [0.5, 0.5]
```

The reciprocal translation vectors of the dual lattice can be obtained by [`reciprocals`](@ref):
```jldoctest unitcell
julia> lattice = Lattice((0.0, 0.0); vectors=[[1.0, 0.0], [0.0, 1.0]]);

julia> reciprocals(lattice)
2-element SVector{2, SVector{2, Float64}} with indices SOneTo(2):
 [6.283185307179586, -0.0]
 [-0.0, 6.283185307179586]
```

### 2.1.2 Building Larger Lattices

Often we need a supercell, a lattice obtained by translating a unitcell several times along each direction. This is done by passing the base lattice and the number of repetitions (or explicit ranges) to the `Lattice` constructor:
```jldoctest lattice-ch3
julia> lattice = Lattice((0.0, 0.0); vectors=[[1.0, 0.0], [0.0, 1.0]], name=:Square);

julia> Lattice(lattice, (2, 2))  # 2x2 open boundaries
Lattice(Square(0:1)(0:1))
  with 4 points:
    [0.0, 0.0]
    [1.0, 0.0]
    [0.0, 1.0]
    [1.0, 1.0]

julia> Lattice(lattice, (1:2, 1:2), ('P', 'P'))  # 2x2 periodic
Lattice(Square[1:2][1:2])
  with 4 points:
    [1.0, 1.0]
    [2.0, 1.0]
    [1.0, 2.0]
    [2.0, 2.0]
  with 2 translation vectors:
    [2.0, 0.0]
    [0.0, 2.0]
```
When the third argument `boundaries` is provided, `'P'` (or `:periodic`) retains the translation vector in that direction (scaled by the number of repetitions), while `'O'` (or `:open`) drops it. In the automatically generated name of the new lattice, each direction is shown as `[range]` when it is periodic and as `(range)` when it is open, as in `Square[1:2][1:2]` and `Square(0:1)(0:1)` above.

## 2.2 Points: Specifying Positions Beyond the Unitcell

With translation symmetry, all points of a lattice are equivalent to those within the origin unitcell. However, things become complicated when bonds are requested. Bonds between different unitcells cannot be compressed into a single unitcell. Therefore, even in the unitcell construction framework, it is necessary to specify a point outside the origin unitcell, which requires extra information beyond a single coordinate if we wish to simultaneously keep track of which point it is equivalent to within the origin unitcell.

In literature, it is customary to express the coordinate $\mathbf{R}_{i\mathbf{m}}$ of a point in the infinite lattice as

```math
\mathbf{R}_{i\mathbf{m}} = \mathbf{r}_i + \mathbf{R}_\mathbf{m}, \qquad \mathbf{R}_{\mathbf{m}}\equiv\sum_{\alpha=1}^{d} m_\alpha \mathbf{a}_\alpha, \; m_\alpha \in \mathbb{Z}
```

Here, $\mathbf{R}_\mathbf{m}$ is the integral coordinate of the unitcell the point belongs to and $\mathbf{r}_i$ is the relative displacement of the point in the unitcell, where $i \in \{1, \ldots, n\}$ identifies the site index and $\mathbf{m} = (m_1, \ldots, m_d)$ identifies the unitcell. Any two of these three coordinates are sufficient to get the full information. In this package, we choose $\mathbf{R}_{i\mathbf{m}}$ and $\mathbf{R}_\mathbf{m}$ as the complete set for an individual lattice point. Additionally, we also include the `site` index for fast lookup of the equivalence of a point within the origin unitcell, although it is redundant in theory. Thus, the [`Point`](@ref) defined in this package has three attributes as follows:
* `site::Int`: the site index of a point that specifies the equivalent point within the origin unitcell
* `rcoordinate::`[`SVector`](https://github.com/JuliaArrays/StaticArrays.jl): the **r**eal **coordinate** of the point ($\mathbf{R}_{i\mathbf{m}}$)
* `icoordinate::`[`SVector`](https://github.com/JuliaArrays/StaticArrays.jl): the **i**ntegral **coordinate** of the unitcell the point belongs to ($\mathbf{R}_\mathbf{m}$)

When constructing a [`Point`](@ref), `rcoordinate` and `icoordinate` can accept tuples or standard vectors as inputs, such as:
```jldoctest unitcell
julia> Point(1, [0.0], [0.0])
Point(1, [0.0], [0.0])

julia> Point(1, (1.5, 0.0), (1.0, 0.0))
Point(1, [1.5, 0.0], [1.0, 0.0])
```
`icoordinate` can be omitted, then it will be initialized by a zero [`SVector`](https://github.com/JuliaArrays/StaticArrays.jl):
```jldoctest unitcell
julia> Point(1, [0.0, 0.5])
Point(1, [0.0, 0.5], [0.0, 0.0])
```

## 2.3 Request for the Bonds of a Lattice

### 2.3.1 Generic bonds

A bond in the narrow sense consists of two points. However, in quantum lattice systems, it is common to refer to generic bonds with only one or more than two points. Additionally, it is convenient to associate a bond with kind information, such as the order of the nearest neighbors of the bond. Thus, the [`Bond`](@ref) is defined as follows:
* `kind`: the kind information of a generic bond
* `points::AbstractVector{<:Point}`: the points a generic bond contains
```jldoctest unitcell
julia> Bond(Point(1, [0.0, 0.0], [0.0, 0.0])) # 1-point bond
Bond(0, Point(1, [0.0, 0.0], [0.0, 0.0]))

julia> Bond(2, Point(1, [0.0, 0.0], [0.0, 0.0]), Point(1, [1.0, 1.0], [1.0, 1.0])) # 2-point bond
Bond(2, Point(1, [0.0, 0.0], [0.0, 0.0]), Point(1, [1.0, 1.0], [1.0, 1.0]))

julia> Bond(:plaquette, Point(1, [0.0, 0.0]), Point(2, [1.0, 0.0]), Point(3, [1.0, 1.0]), Point(4, [0.0, 1.0])) # generic bond with 4 points
Bond(:plaquette, Point(1, [0.0, 0.0], [0.0, 0.0]), Point(2, [1.0, 0.0], [0.0, 0.0]), Point(3, [1.0, 1.0], [0.0, 0.0]), Point(4, [0.0, 1.0], [0.0, 0.0]))
```
Note that the `kind` attribute of a bond with only one point is set to 0.

Iteration over a bond will yield the points it contains:
```jldoctest unitcell
julia> bond = Bond(2, Point(1, [0.0, 0.0], [0.0, 0.0]), Point(2, [1.0, 0.0], [0.0, 0.0]));

julia> length(bond)
2

julia> [bond[1], bond[2]]
2-element Vector{Point{2, Float64}}:
 Point(1, [0.0, 0.0], [0.0, 0.0])
 Point(2, [1.0, 0.0], [0.0, 0.0])

julia> collect(bond)
2-element Vector{Point{2, Float64}}:
 Point(1, [0.0, 0.0], [0.0, 0.0])
 Point(2, [1.0, 0.0], [0.0, 0.0])
```

The coordinate of a bond as a whole is also defined for those that only contain one or two points. The coordinate of a 1-point bond is defined to be the corresponding coordinate of this point, and the coordinate of a 2-point bond is defined to be the corresponding coordinate of the second point minus that of the first:
```jldoctest unitcell
julia> bond = Bond(Point(1, [2.0], [1.0]));

julia> rcoordinate(bond)
1-element SVector{1, Float64} with indices SOneTo(1):
 2.0

julia> icoordinate(bond)
1-element SVector{1, Float64} with indices SOneTo(1):
 1.0

julia> bond = Bond(1, Point(1, [1.0, 1.0], [1.0, 1.0]), Point(2, [0.5, 0.5], [0.0, 0.0]));

julia> rcoordinate(bond)
2-element SVector{2, Float64} with indices SOneTo(2):
 -0.5
 -0.5

julia> icoordinate(bond)
2-element SVector{2, Float64} with indices SOneTo(2):
 -1.0
 -1.0
```

### 2.3.2 Generation of 1-point and 2-point bonds

In this package, we provide the function [`bonds`](@ref) to get the 1-point and 2-point bonds of a lattice:
```julia
bonds(lattice::Lattice, nneighbor::Int) -> Vector{<:Bond}
bonds(lattice::Lattice, neighbors::Neighbors) -> Vector{<:Bond}
```
This function is based on the `KDTree` type provided by the [`NearestNeighbors.jl`](https://github.com/KristofferC/NearestNeighbors.jl) package. In the first method, all bonds up to the `nneighbor`th nearest neighbors are returned, including the 1-point bonds:
```jldoctest unitcell
julia> lattice = Lattice([0.0, 0.0]; vectors=[[1.0, 0.0], [0.0, 1.0]]);

julia> bonds(lattice, 2)
5-element Vector{Bond{Int64, Point{2, Float64}, Vector{Point{2, Float64}}}}:
 Bond(0, Point(1, [0.0, 0.0], [0.0, 0.0]))
 Bond(2, Point(1, [-1.0, -1.0], [-1.0, -1.0]), Point(1, [0.0, 0.0], [0.0, 0.0]))
 Bond(1, Point(1, [0.0, -1.0], [0.0, -1.0]), Point(1, [0.0, 0.0], [0.0, 0.0]))
 Bond(2, Point(1, [1.0, -1.0], [1.0, -1.0]), Point(1, [0.0, 0.0], [0.0, 0.0]))
 Bond(1, Point(1, [-1.0, 0.0], [-1.0, 0.0]), Point(1, [0.0, 0.0], [0.0, 0.0]))
```
Note that for a lattice with translation symmetry, only one representative bond is returned for each class of bonds related by lattice translations, and every two-point representative is oriented so that its **second** point lies in the origin unitcell (with zero `icoordinate`). Every other bond is a translated copy of one of these representatives, possibly with the two points exchanged, so listing them all would be redundant. The returned bonds are not sorted by neighbor order: the single-point bonds come first, and the two-point bonds follow in the order in which they are found by the underlying `KDTree` search.

However, this method is not very efficient, as `KDTree` only searches for bonds with lengths less than a given value, and it does not know the bond lengths for each order of nearest neighbors. This information must be computed first. Therefore, in the second method, [`bonds`](@ref) can accept a new type, [`Neighbors`](@ref), as its second positional parameter to improve efficiency, as it can tell the program the bond length information a priori:
```jldoctest unitcell
julia> lattice = Lattice([0.0, 0.0]; vectors=[[1.0, 0.0], [0.0, 1.0]]);

julia> bonds(lattice, Neighbors(0=>0.0, 1=>1.0, 2=>√2))
5-element Vector{Bond{Int64, Point{2, Float64}, Vector{Point{2, Float64}}}}:
 Bond(0, Point(1, [0.0, 0.0], [0.0, 0.0]))
 Bond(2, Point(1, [-1.0, -1.0], [-1.0, -1.0]), Point(1, [0.0, 0.0], [0.0, 0.0]))
 Bond(1, Point(1, [0.0, -1.0], [0.0, -1.0]), Point(1, [0.0, 0.0], [0.0, 0.0]))
 Bond(2, Point(1, [1.0, -1.0], [1.0, -1.0]), Point(1, [0.0, 0.0], [0.0, 0.0]))
 Bond(1, Point(1, [-1.0, 0.0], [-1.0, 0.0]), Point(1, [0.0, 0.0], [0.0, 0.0]))
```
Meanwhile, an instance of [`Neighbors`](@ref) can also serve as a filter for the generated bonds, selecting those bonds with the given bond lengths:
```jldoctest unitcell
julia> lattice = Lattice([0.0, 0.0]; vectors=[[1.0, 0.0], [0.0, 1.0]]);

julia> bonds(lattice, Neighbors(2=>√2))
2-element Vector{Bond{Int64, Point{2, Float64}, Vector{Point{2, Float64}}}}:
 Bond(2, Point(1, [-1.0, -1.0], [-1.0, -1.0]), Point(1, [0.0, 0.0], [0.0, 0.0]))
 Bond(2, Point(1, [1.0, -1.0], [1.0, -1.0]), Point(1, [0.0, 0.0], [0.0, 0.0]))
```

To obtain generic bonds containing more points, users are encouraged to implement their own `bonds` methods. Pull requests are welcome.

## 2.4 Reciprocal Space

The reciprocal space dual to a lattice is essential for momentum-space analysis. This package provides three common types for different use cases, all of which are subtypes of [`AbstractVector`](https://docs.julialang.org/en/v1/base/arrays/#Base.AbstractVector): the discrete Brillouin zone ([`BrillouinZone`](@ref)), the rectangular reciprocal zone ([`ReciprocalZone`](@ref)), and the high-symmetry reciprocal path ([`ReciprocalPath`](@ref)). All of them support [`length`](https://docs.julialang.org/en/v1/base/arrays/#Base.length-Tuple{AbstractArray}), integer indexing, [`collect`](https://docs.julialang.org/en/v1/base/collections/#Base.collect-Tuple{Any}), and iteration. Iteration and indexing yield the coordinates of the ``\mathbf{k}``-points. These types will be put to use in [Chapter 7](@ref TutorialAlgorithmInterface), where a [`BrillouinZone`](@ref) supplies the ``\mathbf{k}``-point grid for a band-structure calculation, and a [`ReciprocalPath`](@ref) defines the high-symmetry path along which bands are plotted.

### 2.4.1 BrillouinZone

[`BrillouinZone`](@ref) represents the first Brillouin zone of a lattice, parameterized by a grid of discrete ``\mathbf{k}``-points. It can be constructed from a lattice or directly from reciprocal vectors:

```jldoctest reciprocal
julia> lattice = Lattice((0.0, 0.0); vectors=[[1.0, 0.0], [0.0, 1.0]]);

julia> bz = BrillouinZone(lattice, 2)
4-element BrillouinZone{:k, (2, 2), 2, SVector{2, Float64}}:
 [0.0, 0.0]
 [0.0, 3.141592653589793]
 [3.141592653589793, 0.0]
 [3.141592653589793, 3.141592653589793]

julia> bz == BrillouinZone(reciprocals(lattice), 2)
true

julia> length(bz)
4

julia> collect(bz) == [bz[1], bz[2], bz[3], bz[4]]
true
```

The periods along different directions can be different:
```jldoctest reciprocal
julia> lattice = Lattice((0.0, 0.0); vectors=[[1.0, 0.0], [0.0, 1.0]]);

julia> BrillouinZone(lattice, (2, 4))
8-element BrillouinZone{:k, (2, 4), 2, SVector{2, Float64}}:
 [0.0, 0.0]
 [0.0, 1.5707963267948966]
 [0.0, 3.141592653589793]
 [0.0, 4.71238898038469]
 [3.141592653589793, 0.0]
 [3.141592653589793, 1.5707963267948966]
 [3.141592653589793, 3.141592653589793]
 [3.141592653589793, 4.71238898038469]
```

Query methods include [`periods`](@ref), [`dimension`](@ref), and [`scalartype`](@ref):

```jldoctest reciprocal
julia> lattice = Lattice((0.0, 0.0); vectors=[[1.0, 0.0], [0.0, 1.0]]);

julia> bz = BrillouinZone(lattice, (2, 4));

julia> periods(bz)
(2, 4)

julia> dimension(bz)
2

julia> scalartype(bz)
Float64
```

### 2.4.2 ReciprocalZone

[`ReciprocalZone`](@ref) provides a rectangular zone defined by per-axis fractional bounds. It can be constructed from reciprocal vectors with explicit bounds and grid resolution:

```jldoctest reciprocal
julia> recipls = [[1.0, 0.0], [0.0, 1.0]];

julia> rz = ReciprocalZone(recipls, -1//2=>1//2, -1//2=>1//2; length=2)
4-element ReciprocalZone{:k, 2, SVector{2, Float64}, Float64}:
 [-0.5, -0.5]
 [-0.5, 0.0]
 [0.0, -0.5]
 [0.0, 0.0]

julia> length(rz)
4

julia> dimension(rz)
2

julia> scalartype(rz)
Float64
```

Note that the bounds are given in fractional coordinates: the actual ``\mathbf{k}``-points are the fractional values multiplied by the corresponding reciprocal vectors, which in this example happen to be the unit vectors.

Like [`BrillouinZone`](@ref), a [`ReciprocalZone`](@ref) can also be constructed from a lattice, and the lengths along different directions can be different:

```jldoctest reciprocal
julia> lattice = Lattice((0.0, 0.0); vectors=[[1.0, 0.0], [0.0, 1.0]]);

julia> rz = ReciprocalZone(lattice, -1//2=>1//2, -1//2=>1//2; length=(2, 4))
8-element ReciprocalZone{:k, 2, SVector{2, Float64}, Float64}:
 [-3.141592653589793, -3.141592653589793]
 [-3.141592653589793, -1.5707963267948966]
 [-3.141592653589793, 0.0]
 [-3.141592653589793, 1.5707963267948966]
 [0.0, -3.141592653589793]
 [0.0, -1.5707963267948966]
 [0.0, 0.0]
 [0.0, 1.5707963267948966]
```

### 2.4.3 ReciprocalPath

[`ReciprocalPath`](@ref) defines a piecewise-linear path through high-symmetry points, useful for computing and displaying band structures:

```jldoctest reciprocal
julia> recipls = [[1.0, 0.0], [0.0, 1.0]];

julia> rp = ReciprocalPath(
           recipls, (0, 0)=>(1//2, 0), (1//2, 0)=>(1//2, 1//2), (1//2, 1//2)=>(0, 0);
           length=2
       )
7-element ReciprocalPath{:k, SVector{2, Float64}, 3, Tuple{Rational{Int64}, Rational{Int64}}}:
 [0.0, 0.0]
 [0.25, 0.0]
 [0.5, 0.0]
 [0.5, 0.25]
 [0.5, 0.5]
 [0.25, 0.25]
 [0.0, 0.0]

julia> length(rp)
7

julia> dimension(rp)
2
```

A path can also be constructed directly from a lattice:
```jldoctest reciprocal
julia> lattice = Lattice((0.0, 0.0); vectors=[[1.0, 0.0], [0.0, 1.0]]);

julia> rp = ReciprocalPath(
           lattice, (0, 0)=>(1//2, 0), (1//2, 0)=>(1//2, 1//2), (1//2, 1//2)=>(0, 0);
           length=2
       )
7-element ReciprocalPath{:k, SVector{2, Float64}, 3, Tuple{Rational{Int64}, Rational{Int64}}}:
 [0.0, 0.0]
 [1.5707963267948966, 0.0]
 [3.141592653589793, 0.0]
 [3.141592653589793, 1.5707963267948966]
 [3.141592653589793, 3.141592653589793]
 [1.5707963267948966, 1.5707963267948966]
 [0.0, 0.0]
```

### 2.4.4 High-symmetry points: the string macros

For the common one-dimensional, rectangular, and hexagonal Brillouin zones, the standard high-symmetry points can be referred to by name through the exported string macros [`@line_str`](@ref), [`@rectangle_str`](@ref) and [`@hexagon_str`](@ref), with the point names separated by `-` (spaces are ignored), e.g.:

```julia
ReciprocalPath(lattice, line"Γ-X-...")
ReciprocalPath(lattice, rectangle"Γ-X-M-...")
ReciprocalPath(lattice, hexagon"Γ-K-M-...")
```

Here, the supported names and their positions as fractional coordinates in the reciprocal basis are:

* **`line"P₁-P₂-..."`**

  | Names | `Γ`/`Γ₁` | `X`/`X₁` | `Γ₂` | `X₂` |
  |-------|-----------|-----------|------|------|
  | Position | 0 | 1/2 | 1 | -1/2 |

* **`rectangle"P₁-P₂-..."`**

  | Names | `Γ` | `X`/`X₁`  | `X₂` | `Y`/`Y₁`  | `Y₂` | `M`/`M₁`  | `M₂` | `M₃` | `M₄` |
  |-------|-----|-----------|------|-----------|------|-----------|------|------|------|
  | Position | (0, 0) | (1/2, 0) | (-1/2, 0) | (0, 1/2) | (0, -1/2) | (1/2, 1/2) | (-1/2, 1/2) | (-1/2, -1/2) | (1/2, -1/2) |

* **`hexagon"P₁-P₂-..."`**: the angle between the two reciprocal vectors is 120° by default; append `, 60°` (e.g., `hexagon"Γ-K-M-Γ, 60°"`) when it is 60°.

  | Names | `Γ` | `K`/`K₁` | `K₂` | `K₃` | `K₄` | `K₅` | `K₆` | `M`/`M₁` | `M₂` | `M₃` | `M₄` | `M₅` | `M₆` |
  |-------|-----|-----------|------|------|------|------|------|-----------|------|------|------|------|------|
  | 120°  | (0, 0) | (2/3, 1/3) | (1/3, 2/3) | (1/3, -1/3) | (-2/3, -1/3) | (-1/3, -2/3) | (-1/3, 1/3) | (1/2, 1/2) | (1/2, 0) | (0, -1/2) | (-1/2, -1/2) | (-1/2, 0) | (0, 1/2) |
  | 60°   | (0, 0) | (1/3, 1/3) | (2/3, -1/3) | (1/3, -2/3) | (-1/3, -1/3) | (-2/3, 1/3) | (-1/3, 2/3) | (0, 1/2) | (1/2, 0) | (1/2, -1/2) | (0, -1/2) | (-1/2, 0) | (-1/2, 1/2) |

The subscripted names enumerate the symmetry-equivalent copies of a point, so that a path can visit a specific one: for example, `line"Γ₁-Γ₂"` sweeps a full period of the reciprocal lattice instead of staying at Γ, and `hexagon"Γ-K₁-K₂-Γ"` visits along the edges of the small equilateral triangle contained in the hexagonal Brillouin zone.

## Summary

| Type | Purpose |
|------|---------|
| [`Lattice`](@ref) | Unitcell description with translation vectors |
| [`Point`](@ref) | Labeled point with real and integral coordinates |
| [`Bond`](@ref) | Ordered set of points, classified by kind |
| [`Neighbors`](@ref) | Map from neighbor order to bond length |
| [`BrillouinZone`](@ref) | Discrete Brillouin zone grid |
| [`ReciprocalZone`](@ref) | Rectangular reciprocal zone grid |
| [`ReciprocalPath`](@ref) | High-symmetry path in reciprocal space |

With the spatial structure in place, we now turn to [Chapter 3](@ref TutorialDOF), where we define what lives at each lattice point.
