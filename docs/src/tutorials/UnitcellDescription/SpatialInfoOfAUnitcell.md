```@meta
CurrentModule = QuantumLattices
DocTestSetup = quote
    push!(LOAD_PATH, "../../../src/")
    using QuantumLattices
end
```

# Spatial info of a unitcell

The first step toward the complete description of a quantum lattice system is the understanding of the spatial info of a unitcell.

## Point

The basic data structure encoding the spatial info of a unitcell is [`Point`](@ref).

Theoretically, the only information that is needed to determine a point in a lattice is its coordinates in the real space. Now that coordinates sometimes are complicated real numbers and are not convenient for lookup, it is desirable to attach to each point with a sensible id. Then you may agree that the most appropriate data structure representing a point should contain two parts, the id part and the coordinate part. But the story does not end up here. This is because we want to compress the whole spatial info of a lattice into its unitcells. For each lattice, there exists a freedom to choose its unitcells. Sometimes we even need enlarged unitcells. Therefore, something must be adopted to keep the info of which unitcell a point belongs to in the lattice. This information is useful even when we only keep the data of points within a single unitcell because at this time we usually have to obtain those bonds across the unitcell boundaries that must contain a point in other unitcells. Now we arrive at the final structure, just as the [`Point`](@ref) defined in this package, which has three attributes:
* `pid::`[`PID`](@ref): the id of a point
* `rcoord::`[`StaticArrays.SVector`](https://github.com/JuliaArrays/StaticArrays.jl): the coordinates of the point in the real space
* `icoord::`[`StaticArrays.SVector`](https://github.com/JuliaArrays/StaticArrays.jl): the coordinates of the unitcell the point belongs to in the real space

Here [`PID`](@ref) contains two attributes:
* `scope::Any`: the scope of a point
* `site::Int`: the site index of a point

The `:site` attribute is necessary and easy to understand for a point id. Yet sometimes it is more convenient if we can assign extra information to a point id, e.g., a priori knowledge of the groupings of lattice points. Therefore, we provide another attribute, `:scope`, to act as the supplement to the `:site` attribute, which can be anything you want.

Let's see some examples.

You can specify both the `:scope` attribute and the `:site` attribute during the initialization of a [`PID`](@ref):
```jldoctest unitcell
julia> PID("WhateverYouWant", 1)
PID("WhateverYouWant", 1)
```
Or, you can omit the `:scope` attribute:
```jldoctest unitcell
julia> PID(1)
PID('T', 1)
```
Then the `:scope` attribute get a default value `'T'`, which is short for the nick name of my wife.

At the construction of a [`Point`](@ref), `:rcoord` and `:icoord` can accept tuples or usual vectors as inputs, such as
```jldoctest unitcell
julia> Point(PID(1), (0.0,), (0.0,))
Point(PID('T', 1), [0.0], [0.0])

julia> Point(PID(1), [0.0], [0.0])
Point(PID('T', 1), [0.0], [0.0])
```
If the `:icoord` is omitted, it will be initialized by a zero [`StaticArrays.SVector`](https://github.com/JuliaArrays/StaticArrays.jl):
```jldoctest unitcell
julia> Point(PID(1), [0.0])
Point(PID('T', 1), [0.0], [0.0])
```

## Lattice

[`Lattice`](@ref) is the simplest structure to encode all the spatial info within a unitcell. Apparently, it must contain all the points of a unitcell. Besides, a unitcell can assume either open or periodic boundary for every spatial dimension, thus a [`Lattice`](@ref) should also contain the translation vectors. Other stuff also appears to be useful, such as the name, the reciprocals dual to the translation vectors, and the bond length of each order of nearest neighbors. Therefore, in this package, [`Lattice`](@ref) gets seven attributes:
* `name::String`: the name of the lattice
* `pids::Vector{<:PID}`: the pids of the lattice
* `rcoords::Matrix{Float64}`: the rcoords of the lattice
* `icoords::Matrix{Float64}`: the icoords of the lattice
* `vectors::Vector{<:StaticArrays.SVector}`: the translation vectors of the lattice
* `reciprocals::Vector{<:StaticArrays.SVector}`: the reciprocals of the lattice
* `neighbors::Dict{Int, Float64}`: the order-distance map of the nearest neighbors of the lattice
Here, the `:pids`, `:rcoords` and `:icoords` attributes decompose the points in a lattice, which makes it convenient for global operations on the lattice.

Points can be used directly to construct a lattice, whereas `:vectors` and `:neighbors` can be assigned by keyword arguments:
```jldoctest unitcell
julia> Lattice("L2P", [Point(PID(1), [0.0]), Point(PID(2), [1.0])],
           vectors=[[2.0]],
           neighbors=Dict(1=>1.0, 2=>2.0)
           )
Lattice(L2P)
  with 2 points:
    Point(PID('T', 1), [0.0], [0.0])
    Point(PID('T', 2), [1.0], [0.0])
  with 1 translation vector:
    [2.0]
  with 2 orders of nearest neighbors:
    2 => 2.0
    1 => 1.0
```

The `:neighbors` keyword argument can also be a natural number, which sets the highest order of nearest neighbors, and the order-distance map of nearest neighbors can be computed automatically by the construction function:
```jldoctest unitcell
julia> Lattice("L2P", [Point(PID(1), [0.0]), Point(PID(2), [1.0])],
           vectors=[[2.0]],
           neighbors=2
           )
Lattice(L2P)
  with 2 points:
    Point(PID('T', 1), [0.0], [0.0])
    Point(PID('T', 2), [1.0], [0.0])
  with 1 translation vector:
    [2.0]
  with 2 orders of nearest neighbors:
    2 => 2.0
    1 => 1.0
```

It is noted that the `:vectors` and `:neighbors` attributes can also be omitted at the initialization, then `:vectors` will be set to be empty and `:neighbors` to be 1 upon the call of the construction function:
```jldoctest unitcell
julia> Lattice("L2P", [Point(PID(1), [0.0]), Point(PID(2), [1.0])])
Lattice(L2P)
  with 2 points:
    Point(PID('T', 1), [0.0], [0.0])
    Point(PID('T', 2), [1.0], [0.0])
  with 1 order of nearest neighbors:
    1 => 1.0
```

In all cases, the `:reciprocals` attributes need not be assigned because it can be deduced from the input `:vectors`.

## Bonds

One of the most important functions of a lattice is to inquiry the bonds it contains.

A usual bond contains two points, the start point and the end point. This structure is implemented as [`Bond`](@ref), which has three attributes:
* `neighbor::Int`: the nearest neighbor order of the bond
* `spoint::Point`: the start point of the bond
* `epoint::Point`: the end point of the bond
The `:neighbor` attribute provides the a priori info of the nearest neighbor order of a bond, which proves to be quite advantageous in future uses.

There are other types of generalized bonds. In fact, a single point can also be viewed as a kind of bond, namely, the zeroth order nearest neighbor bond. We can also have more complex generalized bonds, such as a plaquette (the minimum four-site square) in the square lattice. All these generalized bonds gather under the abstract type, [`AbstractBond`](@ref), and the generation from a lattice of such generalized bonds can be managed by the type [`Bonds`](@ref). In this package, we only implement two types of concrete generalized bonds, i.e. [`Point`](@ref) and [`Bond`](@ref). Users interested in other types can define them themselves by extending our protocols. In this way, the management of the generation of these user extended bonds can be utilized by [`Bonds`](@ref) without extra modifications. See [`Bonds`](@ref) for more details.

Now let's see a simple example:
```jldoctest unitcell
julia> lattice = Lattice("L2P", [Point(PID(1), [0.0]), Point(PID(2), [1.0])],
           vectors=[[2.0]],
           neighbors=2
           )
Lattice(L2P)
  with 2 points:
    Point(PID('T', 1), [0.0], [0.0])
    Point(PID('T', 2), [1.0], [0.0])
  with 1 translation vector:
    [2.0]
  with 2 orders of nearest neighbors:
    2 => 2.0
    1 => 1.0

julia> Bonds(lattice)
6-element Bonds:
 Point(PID('T', 1), [0.0], [0.0])
 Point(PID('T', 2), [1.0], [0.0])
 Bond(1, Point(PID('T', 1), [0.0], [0.0]), Point(PID('T', 2), [1.0], [0.0]))
 Bond(2, Point(PID('T', 1), [0.0], [0.0]), Point(PID('T', 1), [-2.0], [-2.0]))
 Bond(1, Point(PID('T', 1), [0.0], [0.0]), Point(PID('T', 2), [-1.0], [-2.0]))
 Bond(2, Point(PID('T', 2), [1.0], [0.0]), Point(PID('T', 2), [-1.0], [-2.0]))
```
By default, `Bonds(lattice::Lattice)` generates all the generalized bonds with orders of nearest neighbors specified by the attribute `:neighbors` of the input lattice, including the individual points and the bonds across the periodic boundaries. Note that the bonds whose lengths are not present in the `:neighbors` attribute of the input lattice won't be included in the result, even when their lengths are shorter:
```jldoctest unitcell
julia> lattice = Lattice("L2P", [Point(PID(1), [0.0]), Point(PID(2), [1.0])],
           vectors=[[2.0]],
           neighbors=Dict(2=>2.0)
           )
Lattice(L2P)
  with 2 points:
    Point(PID('T', 1), [0.0], [0.0])
    Point(PID('T', 2), [1.0], [0.0])
  with 1 translation vector:
    [2.0]
  with 1 order of nearest neighbors:
    2 => 2.0

julia> Bonds(lattice)
4-element Bonds:
 Point(PID('T', 1), [0.0], [0.0])
 Point(PID('T', 2), [1.0], [0.0])
 Bond(2, Point(PID('T', 1), [0.0], [0.0]), Point(PID('T', 1), [-2.0], [-2.0]))
 Bond(2, Point(PID('T', 2), [1.0], [0.0]), Point(PID('T', 2), [-1.0], [-2.0]))
```
In other words, the `:neighbors` attribute can be viewed as a filter of the generated bonds (but this filter only affects the [`Bond`](@ref) typed but not the [`Point`](@ref) typed generalized bonds). When the input lattice has no translation vectors, the generated bonds will only contain the individual points and the intra-unitcell bonds, just as expected:
```jldoctest unitcell
julia> lattice = Lattice("L2P", [Point(PID(1), [0.0]), Point(PID(2), [1.0])])
Lattice(L2P)
  with 2 points:
    Point(PID('T', 1), [0.0], [0.0])
    Point(PID('T', 2), [1.0], [0.0])
  with 1 order of nearest neighbors:
    1 => 1.0

julia> Bonds(lattice)
3-element Bonds:
 Point(PID('T', 1), [0.0], [0.0])
 Point(PID('T', 2), [1.0], [0.0])
 Bond(1, Point(PID('T', 1), [0.0], [0.0]), Point(PID('T', 2), [1.0], [0.0]))
```
