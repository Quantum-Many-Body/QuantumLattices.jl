```@meta
CurrentModule = QuantumLattices.Mathematics.VectorSpaces
```

# Vector spaces

A [vector space](https://en.wikipedia.org/wiki/Vector_space) is a linear space, in which the addition of vectors and multiplication of a vector by a scalar are defined.

Vector spaces are frequently encountered in physics, e.g. the Hilbert space in quantum mechanics. In this submodule, we only implement those with finite dimensions. We want to remark that in our implementation, a vector space is a subtype of an abstract vector, therefore, the bases always possess a order, which means, two vector spaces are not considered to be equal to each other even if their corresponding actual mathematical spaces are the same but the orders of the bases are different.

## VectorSpace

[`VectorSpace{B}`](@ref) is the abstraction of a vector space, which has only one type parameter:
* `B<:Any`: the type of the bases of the vector space

Basically, a subtype should implement the following 3 methods:
1) ```julia
   dimension(vs::VectorSpace) -> Int
   ```
   Get the dimension of a vector space
2) ```julia
   Base.getindex(vs::VectorSpace{B}, i::Int)  where B -> B
   ```
   Get the ith basis of a vector space
3) ```julia
   Base.searchsortedfirst(vs::VectorSpace{B}, basis::B) where B -> Int
   ```
   Search the index of a basis in a vector space

Other features include
* comparison: `==` and `isequal`
* iteration: `iterate`
* inquiry: `size`, `findfirst` and `in`

## EnumerativeVectorSpace

[`EnumerativeVectorSpace`](@ref) is the simplest vector space, whose bases are stored in the predefined content `:table`.

## CartesianVectorSpace

[`CartesianVectorSpace`](@ref) defines the abstract class of multiindexable vector spaces, whose bases can be accessed by a Cartesian index.

## NamedVectorSpace

[`NamedVectorSpace{M, NS, BS, VS}`](@ref) defines a multiindexable vector space, each of whose indexable dimensions is associated with a name.

It has four type parameters:
* `M`: mode of the named vector space. It specifies how the indexable dimensions are combined to form the bases of the named vector space, and must take one of the following values:
  - `:⊕`: elements from each indexable dimensions are zipped together to form the bases,
  - `:⊗`: elements from each indexable dimensions are direct producted together to form the bases.
For the `:⊕` mode, all the indexable dimensions should have the same number of elements, and the number of formed bases is equal to this number; for the `:⊗` mode, there are no restriction on the numbers of the indexable dimensions, and the number of the final bases is equal to their product.
* `NS::Tuple{Vararg{Symbol}}`: the names of the indexable dimensions
* `BS<:Tuple`: the eltypes of the indexable dimensions
* `VS<:Tuple{Vararg{AbstractVector}}`: the contents of the indexable dimensions

The concrete types must have the following predefined content:
* `:contents::VS`: storage of the contents of the indexable dimensions

## Manual

```@autodocs
Modules = [VectorSpaces]
Order = [:module, :constant, :type, :macro, :function]
```
