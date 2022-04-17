```@meta
CurrentModule = QuantumLattices.Prerequisites.VectorSpaces
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

## Manual

```@autodocs
Modules = [VectorSpaces]
Order = [:module, :constant, :type, :macro, :function]
```
