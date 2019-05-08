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

## GradedVectorSpace

[`GradedVectorSpace{G,B,V,T}`](@ref) defines the abstract type of graded vector spaces, which are vector spaces that have the extra structure of a grading, which is a decomposition of the vector space into a direct sum of vector subspaces.

It has 4 type parameters
* `G`: the type of the grades
* `B`: the eltype of the subspaces
* `V<:VectorSpace`: the type of the subspaces
* `T<:GradedTables{G,V}`: the type of the subspaces' contents

Concrete subtypes must have the following attribute:
* `:tables::T`: the contents of the subspaces, which must be a [`GradedTables`](@ref).

Specifically, the `dimension`, `getindex` and `searchsortedfirst` methods are overloaded in support of various purposes. For details, please refer to the manual.

## Manul

```@autodocs
Modules=[VectorSpaces]
Order=  [:module,:constant,:type,:macro,:function]
```
