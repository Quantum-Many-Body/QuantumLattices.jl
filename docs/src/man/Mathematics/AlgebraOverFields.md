```@meta
CurrentModule=QuantumLattices.Mathematics.AlgebraOverFields
```

```@setup algebraoverfields
push!(LOAD_PATH,"../../../../src/")
using QuantumLattices.Mathematics.AlgebraOverFields
```

# Algebra over fields

An algebra over a field is a vector space over that field, in which a bilinear operator (often called the "multiplication") between vectors is defined.

With the help of the structure constants of the algebra, the result of the bilinear operation between any arbitary two vectors can be expressed by a sum of individual ones. Therefore, in principle, an algebra can be represented by the complete basis set of its corresponding vector space and a rank-3 tensor encapsulating its structure constants. It is noted that the "bilinear operation" is not restriced to the usual multiplication only. For example, it is the commutator, which is a composition of the usual multiplication and subtraction (for any A and B, the commutator [A,B] is defined as [A,B]≝AB-BA) that serves as the bilinear operator for Lie algebras. In this module, for scalars in the field and elements in the algebra, we only provide the interfaces of the scalar multiplication (including the scalar division) bwteen a sclar and an element, the addition (including the subtraction) and the usual multiplication between two elements. Other complicated operations should be composed from these basic ones.

## SimpleID and ID

[`SimpleID`](@ref) is the building block of the id system of an algebra over a field, while [`ID`](@ref) defines the specific identifier of an element in that algebra.

Generally, the usual multiplication between two elements of an algebra is not commutable, and the rank of the multiplication is just the add-up before the simplication with the help of the algebra structure. We thus use a simple id to denote a single basis of the corresponding vector space, and an id to denote the identifier of an element. With the help of the direct product (`⊗`) of two ids, an over complete id system designed for the whole algebra is constructed. This id system is redundant because it does not reflects the structure constants of the algebra, which reduces independent basis elements. Extra mechanisms should be provided to kill this redundancy, which goes beyond the current module. Users should define them themselves.

## IdSpace

[`IdSpace`](@ref) defines the complete vector space that corresponds to an algebra.

An id space uses a set of simple id bases to generated all ranked ids, with the former represented by an instance of [`SimpleVectorSpace`](@ref) and stored in the attribute `:sids`, while the latter by an instance of [`GradedTables`](@ref) in the attribute `:tables`. For the sick of the orderings among different ids, an id will be represented by a multiindex, which are determined by the corresponding sequences of its simple ids in the basis set. The id system can be complete or over complete, depending on whether the attribute `:tables` contain multiindices that actually represent the same basis of the algebra. Besides, the methods
```julia
Base.getindex(idspace::IdSpace,i::Int) -> ID
```
and
```julia
Base.findfirst(id::ID,idspace::IdSpace) -> Int
```
are provided to get the ith id and to find the index of an id of a idspace, respectively.

## Element and Elements

[`Element`](@ref) defines a single element of an algebra while [`Elements`](@ref) defines an exprssion composed of several elements from an algebra.

The first and second attributes of an [`Element`](@ref) must be
- `value::Number`: the coefficient of the element
- `id::ID`: the id of the element

Arithmetic operations (`+`, `-`, `*`, `/`) bwteen a scalar, an [`Element`](@ref) or an [`Elements`](@ref) is defined. See Manual for details.

## Manual

```@autodocs
Modules=[AlgebraOverFields]
Order=  [:module,:constant,:type,:macro,:function]
```
