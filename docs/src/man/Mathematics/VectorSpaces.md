```@meta
CurrentModule=QuantumLattices.Mathematics.VectorSpaces
```

# Vector spaces

A [vector space](https://en.wikipedia.org/wiki/Vector_space) is a linear space, in which the addition of vectors and multiplication of a vector by a scalar are defined.

Vector spaces are frequently encountered in physics, e.g. the Hilbert space in quantum mechanics. In this submodule, we only implement those with finite dimensions. We want to remark that in our implementation, a vector space is a subtype of an abstract vector, therefore, the bases always possess a order, which means, two vector spaces are not considered to be equal to each other even if their corresponding actual mathmatical spaces are the same but the the orders of the bases are different.

## VectorSpace

[`VectorSpace{B}`](@ref) is the abstaction of a vector space, which has only one type parameter:
* `B<:Any`: the type of the bases of the vector space

Basically, a subtype should implement the following 3 methods:
1) ```julia
   dimension(vs::VectorSpace) -> Int
   ```
   Get the dimension of a vector space
2) ```julia
   Base.getindex(vs::VectorSpace{B},i::Int)  where B -> B
   ```
   Get the ith basis of a vector space
3) ```julia
   Base.searchsortedfirst(vs::VectorSpace{B},basis::B) where B -> Int
   ```
   Search the index of a basis in a vector space

However, we provide several interfaces, including type traits and methods to deal with common situations:
1) A vector space whose bases are stored in a table under the attribute name `:table` can be ascribed to the [`HasTable`](@ref) trait and the [`TableSorted`](@ref) trait.
   Specifically, the first trait must be implemented as
   ```julia
   HasTable(::Type{SubType})=HasTable(true)
   ```
   While, if the table is unsorted, the second trait should be implemented as
   ```julia
   TableSorted(::Type{SubType})=TableSorted(false)
   ```
   and if the table is sorted, the second trait should be implemented as
   ```julia
   TableSorted(::Type{SubType})=TableSorted(true)
   ```
2) A vector space whose bases may be represented by a multiindex (Cartesian index) can be ascribed to the traits [`IsMultiIndexable`](@ref) and [`MultiIndexOrderStyle`](@ref).
   Specifically, the first trait must be implemented as
   ```julia
   IsMultiIndexable(::Type{SubType})=IsMultiIndexable(true)
   ```
   While, if the order style of the multiindex is C/C++ like, the second trait shoule be implemented as
   ```julia
   MultiIndexOrderStyle(::Type{SubType})=MultiIndexOrderStyle('C')
   ```
   and if the order style is Fortran like, the second trait shoule be implemented as
   ```julia
   MultiIndexOrderStyle(::Type{SubType})=MultiIndexOrderStyle('F')
   ```
   Furthermore, it should implement the following methods
   * ```julia
     rank(::Type{SubType}) -> Int
     ```
     Get the rank of a multiindexable vector space.
   * ```julia
     dims(vs::SubType) -> NTuple{vs|>typeof|>rank,Int}
     ```
     Get the dimensions along each axes of a multiindexable vector space.
   * ```julia
     inds(basis,vs::SubType) ->  NTuple{vs|>typeof|>rank,Int}
     ```
     Get the Cartesian index representation of a basis in a multiindexable vector space.
   * ```julia
     eltype(SubType).name.wrapper(index::NTuple{N,Int},vs::SubType)
     ```
     Construct a basis from a tuple of integers and a multiindexable vector space.
   Note that a multiindexable vector space can also have a sorted or unsorted table. But then the trait [`MultiIndexOrderStyle`](@ref) takes no effects and the sequences of its bases will be completely determined by its attribute `:table`.
If the type taits and methods are defined properly as stated above, the `dimension`, `getindex` and `searchsortedfirst` methods get default implementations. No need to worry about them any more.

Other features include
* comparison: `==` and `isequal`
* iteration: `iterate`
* inquiry: `size`, `findfirst` and `in`

## SimpleVectorSpace

[`SimpleVectorSpace{S,B,N}`](@ref) is the simplest vector space, whose bases are stored in the attribute `:table::NTuple{N,B}` as an ntuple.

The `:table` attribute can be sorted or unsorted, which is determined by the type parameter `S`, with `'T'` for sorted and `'F'` for unsorted.

## OrderedIndices

[`OrderedIndices{N}`](@ref) defines the simplest abstract class of multiindexable vector spaces, whose bases are just tuples of integers.

This class of vector spaces must have the following attribute:
`dims::NTuple{N,Int}`: the dimesnions of the Cartesian indices along all axes

### SimpleIndices

[`SimpleIndices{M,N}`](@ref) is the simple ordered Cartesian indices.

It is worth noting that
1) It can be C/C++ ordered or Fortran ordered depending on the first type parameter `M`, with `'C'` for the former and `'F'` the latter.
2) For its bases (Cartesian indices), there is no restriction except that they should be in the proper range defined by its `dims`.

### TabledIndices

[`TabledIndices{S,N}`](@ref) defines the tabled ordered Cartesian indices.

Compared to [`SimpleIndices`](@ref), the bases of this kind of vector spaces are stored in the attribute `:table`, which must be a vector of tuple of integers. The `:table` attribute can be sorted or unsorted, which is determined by the type parameter `S`, with `'T'` for sorted and `'F'` for unsorted. This type suits the situations when the Cartesian indices are restricted by extra conditions except that propoesed by the attribute `:dims`.

## NamedVectorSpace

[`NamedVectorSpace{M,NS,BS,VS}`](@ref) defines a multiindexable vector space, each of whose indexable dimensions is associated with a name.

It has four type parameters:
* `M`: mode of the named vector space. It specifies how the indexable dimensions are combined to form the bases of the named vector space, and must take one of the following values:
  - `:zip`: elements from each indexable dimensions are zipped together to form the bases,
  - `:product`: elements from each indexable dimensions are direct producted together to form the bases.
For the `:zip` mode, all the indexable dimensions should have the same number of elements, and the number of formed bases is equal to this number; for the `:product` mode, there are no restriction on the numbers of the indexable dimensions, and the number of the final bases is equal to their product.
* `NS::Tuple{Vararg{Symbol}}`: the names of the indexable dimensions
* `BS<:Tuple`: the eltypes of the indexable dimensions
* `VS<:Tuple{Vararg{AbstractVector}}`: the contents of the indexable dimensions

The concrete types must have the following attribute:
* `:contents::VS`: storage of the contents of the indexable dimensions

By default, a named vector space uses C order for the indexable dimensions when the mode is `:product`. You can change it to F order for your own subtypes by defining the [`MultiIndexOrderStyle`](@ref) trait.

## Manul

```@autodocs
Modules=[VectorSpaces]
Order=  [:module,:constant,:type,:macro,:function]
```
