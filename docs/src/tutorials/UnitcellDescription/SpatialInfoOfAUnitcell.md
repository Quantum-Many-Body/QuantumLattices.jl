```@meta
CurrentModule = QuantumLattices
DocTestSetup = quote
    push!(LOAD_PATH, "../../../src/")
    using QuantumLattices
end
```

# Spatial information of a unitcell

The first step toward the complete description of a quantum lattice system is the understanding of the spatial information of a unitcell.

## Construction of a lattice

In general, a lattice has translation symmetry. This symmetry introduces an equivalence relation for the points in a lattice when they can be translated into each other by multiple times of the translation vectors. This observation sets the mathematical foundation of the unitcell construction. Therefore, it is enough for a lattice to restrict all points within the origin unitcell together with the translation vectors.

[`Lattice`](@ref) is the simplest structure to encode all the spatial information within the origin unitcell. Apparently, it must contain all the coordinates of the points in the origin unitcell and the translation vectors of the lattice. Other stuff appears to be useful as well, such as the name of the lattice and the reciprocals dual to the translation vectors. Therefore, in this package, [`Lattice`](@ref) has four attributes:
* `name::Symbol`: the name of the lattice
* `coordinates::Matrix{<:Number}`: the coordinates of the points within the origin unitcell
* `vectors::Vector{<:StaticArrays.SVector}`: the translation vectors of the lattice
* `reciprocals::Vector{<:StaticArrays.SVector}`: the reciprocals dual to the translation vectors

[`Lattice`](@ref) can be constructed by offering the coordinates, with optional keyword arguments to specify its name and translation vectors:
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

julia> Lattice((0.0, 0.0, 0.0); name=:Cube, vectors=[[1.0, 0.0, 0.0], [0.0, 1.0, 0.0], [0.0, 0.0, 1.0]])
Lattice(Cube)
  with 1 point:
    [0.0, 0.0, 0.0]
  with 3 translation vectors:
    [1.0, 0.0, 0.0]
    [0.0, 1.0, 0.0]
    [0.0, 0.0, 1.0]
```
The coordinates could be specified by vectors or tuples. It is noted that the `:reciprocals` attribute need not be assigned because it can be deduced from the input `:vectors`.

Iteration over a lattice will get the coordinates of the points in it:
```jldoctest unitcell
julia> lattice = Lattice((0.0, 0.0), (0.5, 0.5); vectors=[[1.0, 0.0], [0.0, 1.0]]);

julia> length(lattice)
2

julia> [lattice[1], lattice[2]]
2-element Vector{StaticArrays.SVector{2, Float64}}:
 [0.0, 0.0]
 [0.5, 0.5]

julia> collect(lattice)
2-element Vector{StaticArrays.SVector{2, Float64}}:
 [0.0, 0.0]
 [0.5, 0.5]
```

## Request for the bonds of a lattice

Before the introduction of how to obtain the bonds of a lattice, let's discuss more about the unitcell construction to clarify the logic behind the definitions of the [`Point`](@ref) type and the [`Bond`](@ref) type in this package.

### Point

With the translation symmetry, all points of a lattice are equivalent to those within the origin unitcell. However, it becomes complicated when the bonds are requested. The bonds inter different unitcells cannot be compressed into a single unitcell. Therefore, even in the unitcell construction framework, it turns out to be unavoidable to specify a point outside the origin unitcell, which requires extra information beyond a single coordinate if we want to remember which point it is equivalent to within the origin unitcell at the same time. In fact, it is customary in literature to express the coordinate $\mathbf{R}$ of a point in a lattice as $\mathbf{R}=\mathbf{R}_i+\mathbf{r}$, where $\mathbf{R}_i$ is the integral coordinate of the unitcell the point belongs to and $\mathbf{r}$ is the relative displacement of the point in the unitcell. Apparently, any two of these three coordinates are complete to get the full information. In this package, we choose $\mathbf{R}$ and $\mathbf{R}_i$ as the complete set for a individual lattice point. Besides, we also associate a `:site` index with a point for the fast lookup for its equivalence within the origin unitcell although it is redundant in theory. Thus, the [`Point`](@ref) defined in this package has three attributes as follows:
* `site::Int`: the site index of a point that specifies the equivalent point within the origin unitcell
* `rcoordinate::`[`StaticArrays.SVector`](https://github.com/JuliaArrays/StaticArrays.jl): the **r**eal **coordinate** of the point ($\mathbf{R}$)
* `icoordinate::`[`StaticArrays.SVector`](https://github.com/JuliaArrays/StaticArrays.jl): the **i**ntegral **coordinate** of the unitcell the point belongs to ($\mathbf{R}_i$)

At the construction of a [`Point`](@ref), `:rcoordinate` and `:icoordinate` can accept tuples or usual vectors as inputs, such as
```jldoctest unitcell
julia> Point(1, [0.0], [0.0])
Point(1, [0.0], [0.0])

julia> Point(1, (1.5, 0.0), (1.0, 0.0))
Point(1, [1.5, 0.0], [1.0, 0.0])
```
`:icoordinate` can be omitted, then it will be initialized by a zero [`StaticArrays.SVector`](https://github.com/JuliaArrays/StaticArrays.jl):
```jldoctest unitcell
julia> Point(1, [0.0, 0.5])
Point(1, [0.0, 0.5], [0.0, 0.0])
```

### Bond

A bond in the narrow sense consist of two points. However, in quantum lattice systems, it is common to refer to generic bonds with only one or more than two points. In addition, it is also convenient to associate a bond with a kind information, such as the order of the nearest neighbors of the bond. Thus, the [`Bond`](@ref) is defined as follows:
* `kind`: the kind information of a generic bond
* `points::Vector{<:Point}`: the points a generic bond contains
```jldoctest unitcell
julia> Bond(Point(1, [0.0, 0.0], [0.0, 0.0])) # 1-point bond
Bond(0, Point(1, [0.0, 0.0], [0.0, 0.0]))

julia> Bond(2, Point(1, [0.0, 0.0], [0.0, 0.0]), Point(1, [1.0, 1.0], [1.0, 1.0])) # 2-point bond
Bond(2, Point(1, [0.0, 0.0], [0.0, 0.0]), Point(1, [1.0, 1.0], [1.0, 1.0]))

julia> Bond(:plaquette, Point(1, [0.0, 0.0]), Point(2, [1.0, 0.0]), Point(3, [1.0, 1.0]), Point(4, [0.0, 1.0])) # generic bond with 4 points
Bond(:plaquette, Point(1, [0.0, 0.0], [0.0, 0.0]), Point(2, [1.0, 0.0], [0.0, 0.0]), Point(3, [1.0, 1.0], [0.0, 0.0]), Point(4, [0.0, 1.0], [0.0, 0.0]))
```
It is noted that the `:kind` attribute of a bond with only one point is set to be 0.

Iteration over a bond will get the points it contains:
```jldoctest unitcell
julia> bond = Bond(:plaquette, Point(1, [0.0, 0.0]), Point(2, [1.0, 0.0]), Point(3, [1.0, 1.0]), Point(4, [0.0, 1.0]));

julia> length(bond)
4

julia> [bond[1], bond[2], bond[3], bond[4]]
4-element Vector{Point{2, Float64}}:
 Point(1, [0.0, 0.0], [0.0, 0.0])
 Point(2, [1.0, 0.0], [0.0, 0.0])
 Point(3, [1.0, 1.0], [0.0, 0.0])
 Point(4, [0.0, 1.0], [0.0, 0.0])

julia> collect(bond)
4-element Vector{Point{2, Float64}}:
 Point(1, [0.0, 0.0], [0.0, 0.0])
 Point(2, [1.0, 0.0], [0.0, 0.0])
 Point(3, [1.0, 1.0], [0.0, 0.0])
 Point(4, [0.0, 1.0], [0.0, 0.0])
```

The coordinate of a bond as a whole is also defined for those that only contains one or two points. The coordinate of a 1-point bond is defined to be the corresponding coordinate of this point, and the coordinate of a 2-point bond is defined to be the corresponding coordinate of the second point minus that of the first:
```jldoctest unitcell
julia> bond1p = Bond(Point(1, [2.0], [1.0]));

julia> rcoordinate(bond1p)
1-element StaticArrays.SVector{1, Float64} with indices SOneTo(1):
 2.0

julia> icoordinate(bond1p)
1-element StaticArrays.SVector{1, Float64} with indices SOneTo(1):
 1.0

julia> bond2p = Bond(1, Point(1, [1.0, 1.0], [1.0, 1.0]), Point(2, [0.5, 0.5], [0.0, 0.0]));

julia> rcoordinate(bond2p)
2-element StaticArrays.SVector{2, Float64} with indices SOneTo(2):
 -0.5
 -0.5

julia> icoordinate(bond2p)
2-element StaticArrays.SVector{2, Float64} with indices SOneTo(2):
 -1.0
 -1.0
```

### Generation of 1-point and 2-point bonds of a lattice

In this package, we provide the function [`bonds`](@ref) to get the 1-point and 2-point bonds of a lattice:
```julia
bonds(lattice::Lattice, nneighbor::Int) -> Vector{<:Bond}
bonds(lattice::Lattice, neighbors::Neighbors) -> Vector{<:Bond}
```
which is based on the `KDTree` type provided by the [`NearestNeighbors.jl`](https://github.com/KristofferC/NearestNeighbors.jl) package. In the first method, all the bonds up to the `nneighbor`th nearest neighbors are returned, including the 1-point bonds:
```jldoctest unitcell
julia> lattice = Lattice([0.0, 0.0]; vectors=[[1.0, 0.0], [0.0, 1.0]]);

julia> bonds(lattice, 2)
5-element Vector{Bond{Int64, Point{2, Float64}}}:
 Bond(0, Point(1, [0.0, 0.0], [0.0, 0.0]))
 Bond(2, Point(1, [-1.0, -1.0], [-1.0, -1.0]), Point(1, [0.0, 0.0], [0.0, 0.0]))
 Bond(1, Point(1, [0.0, -1.0], [0.0, -1.0]), Point(1, [0.0, 0.0], [0.0, 0.0]))
 Bond(2, Point(1, [1.0, -1.0], [1.0, -1.0]), Point(1, [0.0, 0.0], [0.0, 0.0]))
 Bond(1, Point(1, [-1.0, 0.0], [-1.0, 0.0]), Point(1, [0.0, 0.0], [0.0, 0.0]))
```
However, this method is not so efficient, as `KDTree` only searches the bonds with the lengths less than a value, and it does not know the bond lengths for each order of nearest neighbors. Such information must be computed as first. Therefore, in the second method, [`bonds`](@ref) can accept a new type, the [`Neighbors`](@ref), as its second positional parameter to improve the efficiency, which could tell the program the information of the bond lengths in priori:
```jldoctest unitcell
julia> lattice = Lattice([0.0, 0.0]; vectors=[[1.0, 0.0], [0.0, 1.0]]);

julia> bonds(lattice, Neighbors(0=>0.0, 1=>1.0, 2=>√2))
5-element Vector{Bond{Int64, Point{2, Float64}}}:
 Bond(0, Point(1, [0.0, 0.0], [0.0, 0.0]))
 Bond(2, Point(1, [-1.0, -1.0], [-1.0, -1.0]), Point(1, [0.0, 0.0], [0.0, 0.0]))
 Bond(1, Point(1, [0.0, -1.0], [0.0, -1.0]), Point(1, [0.0, 0.0], [0.0, 0.0]))
 Bond(2, Point(1, [1.0, -1.0], [1.0, -1.0]), Point(1, [0.0, 0.0], [0.0, 0.0]))
 Bond(1, Point(1, [-1.0, 0.0], [-1.0, 0.0]), Point(1, [0.0, 0.0], [0.0, 0.0]))
```
Meanwhile, an instance of [`Neighbors`](@ref) could also serve as a filter of the generated bonds, which select those bonds with the given bond lengths:
```jldoctest unitcell
julia> lattice = Lattice([0.0, 0.0]; vectors=[[1.0, 0.0], [0.0, 1.0]]);

julia> bonds(lattice, Neighbors(2=>√2))
2-element Vector{Bond{Int64, Point{2, Float64}}}:
 Bond(2, Point(1, [-1.0, -1.0], [-1.0, -1.0]), Point(1, [0.0, 0.0], [0.0, 0.0]))
 Bond(2, Point(1, [1.0, -1.0], [1.0, -1.0]), Point(1, [0.0, 0.0], [0.0, 0.0]))
```
 
 To obtain generic bonds containing more points, user are encouraged to implement their own `bonds` methods. Pull requests are welcomed.
