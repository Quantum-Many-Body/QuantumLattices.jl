```@meta
CurrentModule = QuantumLattices.Toolkit
DocTestFilters = [r"var\".*\"", r"generic function with [0-9]* method", r".*s \(.*\% GC\)", r"evals/sample:.*"]
DocTestSetup = quote
    push!(LOAD_PATH, "../../../../src/")
    using QuantumLattices.Toolkit
    import QuantumLattices.Toolkit: isparameterbound, parameternames, contentnames, getcontent
end
```

```@setup toolkit
push!(LOAD_PATH, "../../../src/")
using QuantumLattices.Toolkit
```

# Toolkit

*This module contains the toolkit of the package.*

The constants, types, macros, functions defined in this module will **not** be exported by the package. Instead, they serve as the prerequisites. The range of the contents are quite wide, but basically, they fall into two categories:
* Utilities, such as global constants and miscellaneous tiny useful functions;
* Basic data structures as supplements to the `Julia.Base` and other common packages.

## Utilities

```@docs
atol
rtol
Float
concatenate
tostr
subscript
superscript
delta
id
value
DirectSummedIndices
Segment
```

## Combinatorics

The combinations and permutations of an indexable object are implemented, with duplicate elements allowed or not. Compared to another Julia package [Combinatorics](https://github.com/JuliaMath/Combinatorics.jl), the iterators return tuples instead of vectors, which could greatly decrease the memory allocation times and improves the code efficiency.

[`Combinatorics{M, C}`](@ref) is the abstract type of all combinatorial algorithms. It has two type parameters:
* `M`: the number of elements to be taken
* `C`: the type of the collection of candidate elements
To avoid memory allocation, the iteration of a concrete combinatorial algorithm returns a tuple, whose length is `M` and eltype is `eltype(C)`.

### Combinations and DuplicateCombinations

[`Combinations{M, C}`](@ref) and [`DuplicateCombinations{M, C}`](@ref) generate all the combinations of `M` elements from an indexable collection whose type is `C`, with the differences being that the former forbids duplicate elements in the combinations while the latter allows.

All combinations of 2 integers taken from 1 to 3 without duplicate:
```@example toolkit
Combinations{2}(1:3) |> collect
```

All combinations of 2 integers taken from 1 to 3 with duplicate allowed:
```@example toolkit
DuplicateCombinations{2}(1:3) |> collect
```

### Permutations and DuplicatePermutations

[`Permutations{M, C}`](@ref) and [`DuplicatePermutations{M, C}`](@ref) generate all the permutations of `M` elements from an indexable collection whose type is `C`, with the differences being that the former forbids duplicate elements in the permutations while the latter allows.

All permutations of 2 integers taken from 1 to 3 without duplicate:
```@example toolkit
Permutations{2}(1:3) |> collect
```

All permutations of 2 integers taken from 1 to 3 with duplicate allowed:
```@example toolkit
DuplicatePermutations{2}(1:3) |> collect
```

### Manual

```@docs
Combinatorics
Combinations
DuplicateCombinations
Permutations
DuplicatePermutations
```

## Traits

Trait functions and trait types that are useful to the package are defined.

Generally speaking, traits in Julia could fall into two categories according to their usages, the first may be term as "type helpers" and the second are usually called "Holy traits" named after [Tim Holy](https://github.com/timholy). Type helpers aim at the inquiry, alteration and computation of the compile-time information of types, while Holy traits can be applied as an alternative to multi-inheritance by use of the Julia multidispatch feature.

### Type helpers

Type helpers are important for the generic programming in Julia, especially in the design of generic interfaces and abstract types.

Let's see a simple situation, i.e. the elemental addition of two vectors of numbers. The numbers can assume different types and the type of the result depends on both of them, for example, the result between two vectors of integers is a vector of integers while that between a vector of integers and a vector of floats is a vector of floats. Of course, one can explicitly define every elemental addition function between any two different types of vectors of numbers, like this:
```julia
# wrong design pattern

function elementaladdition(v₁::Vector{Int64}, v₂::Vector{Int64})
    result = Int[]
    ...
end
function elementaladdition(v₁::Vector{Int64}, v₂::Vector{Float64})
    result = Float64[]
    ...
end
...
...
```
Writing down all such methods is already a heavy repetition. What's worse, you will quickly find that a lot more functions, such as the elemental subtraction, elemental multiplication and elemental division, are waiting for you to implement. This is a total disaster.

The correct strategy is to define the promotion rule of any two types of numbers and use it to define the type of the result:
```julia
# correct design pattern

promotion(::Type{Int64}, ::Type{Int64}) = Int64
promotion(::Type{Int64}, ::Type{Float64}) = Float64
...
...

function elementaladdition(v₁::Vector{T₁}, v₂::Vector{T₂}) where {T₁<:Number, T₂<:Number}
    result = promotion(T₁, T₂)[]
    ...
end
function elementalsubtraction(v₁::Vector{T₁}, v₂::Vector{T₂}) where {T₁<:Number, T₂<:Number}
    result = promotion(T₁, T₂)[]
    ...
end
...
...
```
The promotion rule applies equally to all the arithmetic operations on numbers. Therefore, tedious code repetition could be avoided with it. In fact, similar promotion rules have already been defined in Julia base, and the default implementations of arithmetic operations in Julia are indeed based on them (see [`Base.promote_rule`](https://docs.julialang.org/en/v1/base/base/#Base.promote_type) and [`Base.promote_type`](https://docs.julialang.org/en/v1/base/base/#Base.promote_rule)). When new user-defined numeric types are introduced, the only things you need to do is to add new promotion rules and implement a few basic arithmetic functions for these new types. Then quite a lot of generic codes could apply to them without any modification.

#### Type helpers with type parameters

The input and output types of a promotion rule are known at compile time, thus, the promotion rule is a trait function aiming at the computation of compile-time information of types. Trait functions dealing with the inquiries of compile-time information of types are also widely used in Julia, such as the [`eltype`](https://docs.julialang.org/en/v1/base/collections/#Base.eltype) function of [`Vector`](https://docs.julialang.org/en/v1/base/arrays/#Base.Vector):
```jldoctest traits
julia> eltype(Vector{String})
String
```

For a user-defined parametric type, it is also useful to provide an inquiry function to access to the type parameters:
```jldoctest traits
julia> struct Hi{T<:Number}
           content::T
       end
       contenttype(::Type{Hi{T}}) where T<:Number = T
       contenttype(Hi{Int64})
Int64
```

However, the above defined function `contenttype` could not apply to a [`UnionAll`](https://docs.julialang.org/en/v1/manual/types/#UnionAll-Types) type, such as `Hi{<:Real}`:
```jldoctest traits
julia> contenttype(Hi{<:Real})
ERROR: MethodError: no method matching contenttype(::Type{Hi{<:Real}})
[...]
```

In fact in Julia base, all such inquiry functions, e.g., the `eltype` function, work poor for the `UnionAll` types:
```jldoctest traits
julia> eltype(Vector{<:Real})
Any
```
In concept, `eltype(Vector{<:Real}` should return `Real` instead of `Any` as every element in `Vector{<:Real}` is indeed a real number. Similarly, we expect that `contenttype(Hi{<:Real})` should also give us `Real`. Unfortunately, functions defined in the similar form like this could never achieve such goals. Julia base doesn't provide generic functions to access or change the information of the parameters of a type. In this module, we try to fill this gap with a set of generic trait functions.

##### Access or change the type parameters by their position orders

The most direct information of the parameters of a type is their position orders. We provide [`parametertype`](@ref) to access to them by such information:
```jldoctest traits
julia> parametertype(Hi{<:Real}, 1)
Real

julia> parametertype(Vector{<:Real}, 1)
Real

julia> parametertype(Vector{<:Real}, 2)
1
```

You can use [`parametercount`](@ref) to inquire the total number of the parameters of a type:
```jldoctest traits
julia> parametercount(Hi)
1

julia> parametercount(Vector)
2
```
It is noted that `Vector` has 2 type parameters because it is just a type alias for `Array{T, 1} where T`.

To change the parameters of a type, [`reparameter`](@ref) can be used:
```jldoctest traits
julia> reparameter(Hi{Int64}, 1, Real)
Hi{<:Real}

julia> reparameter(Vector{Int64}, 1, Real)
Vector{<:Real} (alias for Array{<:Real, 1})

julia> reparameter(Vector{<:Real}, 2, 3)
Array{<:Real, 3}

julia> reparameter(Hi{Int64}, 1, Real, false)
Hi{Real}

julia> reparameter(Hi{Int64}, 1, Real, true)
Hi{<:Real}
```
We want to remark that by providing the fourth positional argument with the `true` value, a `UnionAll` type could be generated. When the fourth positional argument is omitted, it is actually determined by another trait function, i.e., [`isparameterbound`](@ref). This function judges whether an input type should be considered as the upper bound of the new parameter of a type. By default, it is defined to be
```julia
isparameterbound(::Type{}, ::Val{}, ::Type{D}) where D = !isconcretetype(D)
isparameterbound(::Type{}, ::Val{}, ::Any) = false
```
This function can be overloaded to change the behavior for a certain type:
```jldoctest traits
julia> isparameterbound(::Type{<:Vector}, ::Val{1}, ::Type{D}) where D = false;

julia> reparameter(Vector, 1, Real)
Vector{Real} (alias for Array{Real, 1})
```
The second positional argument of [`isparameterbound`](@ref) must be of type `Val` because in principle you should be able to assign different behaviors for different parameters of a type separately. If it is of type `Integer`, a single overloading would change the behaviors for all.

Besides, you can inquire all the parameters of a type by [`parametertypes`](@ref):
```jldoctest traits
julia> parametertypes(Hi{<:Real})
Tuple{Real}

julia> parametertypes(Vector{Int64})
Tuple{Int64, 1}
```
The obtained type parameters are stored as those of a `Tuple`.

At the same time, you can change all the parameters of a type by [`fulltype`](@ref):
```jldoctest traits
julia> fulltype(Hi{Int64}, Tuple{Real})
Hi{<:Real}

julia> fulltype(Hi{Int64}, Tuple{Real}, (false,))
Hi{Real}

julia> fulltype(Vector{Int64}, Tuple{Real, 2})
Matrix{Real} (alias for Array{Real, 2})

julia> fulltype(Vector{Int64}, Tuple{Real, 2}, (true, false))
Matrix{<:Real} (alias for Array{<:Real, 2})
```
Like [`reparameter`](@ref), the last positional argument of [`fulltype`](@ref) could determine whether the corresponding types specified by the type parameters of the input `Tuple` should be considered as the upper bounds of the new parameters of a type. When this argument is omitted, it is determined by another trait function [`isparameterbounds`](@ref), which successively calls the [`isparameterbound`](@ref) function to determine the behaviors for all the parameters of a type as the literal indicates.

##### Associate type parameters with names

Sometimes, it is more convenient to associate names with the parameters of a type, and then access or change them by their names. This can be done by overloading the [`parameternames`](@ref) trait function for a certain type:
```jldoctest traits
julia> parameternames(::Type{<:Hi}) = (:content,)
parameternames (generic function with 3 methods)
```

Now, you can inquire the name of a type parameter by [`parametername`](@ref) with the given position order or vice versa by [`parameterorder`](@ref):
```jldoctest traits
julia> parametername(Hi, 1)
:content

julia> parameterorder(Hi, :content)
1
```

You can also inquire whether a type has a parameter with the given name by [`hasparameter`](@ref):
```jldoctest traits
julia> hasparameter(Hi, :content)
true

julia> hasparameter(Hi, :others)
false
```

And [`parametertype`](@ref) and [`reparameter`](@ref) can be applied by the name of a type parameter instead of its position order:
```jldoctest traits
julia> parametertype(Hi{<:Real}, :content)
Real

julia> reparameter(Hi{Int}, :content, Real)
Hi{<:Real}

julia> reparameter(Hi{Int}, :content, Real, false)
Hi{Real}
```

To change the [`reparameter`](@ref) behavior when its last positional argument is omitted, you should overload the [`isparameterbound`](@ref) function accordingly, e.g.:
```jldoctest traits
julia> isparameterbound(::Type{<:Hi}, ::Val{:content}, ::Type{D}) where D = false;

julia> reparameter(Hi{Int}, :content, Real)
Hi{Real}
```
!!! note
    Accessing or altering a parameter of a type by its name is independent from that by its position order. Thus, even the following method
    ```julia
    isparameterbound(::Type{<:Hi}, ::Val{1}, D)
    ```
    has been overloaded, it doesn't affect the result of the function call like
    ```julia
    reparameter(Hi{Int}, :content, Real)
    ```
    Rather, it only affect the result of the function call like
    ```julia
    reparameter(Hi{Int}, 1, Real)
    ```
    To change the default behavior of the former function call, you must overload the following method manually as well
    ```julia
    isparameterbound(::Type{<:Hi}, ::Val{:content}, D)
    ```

A new trait function [`parameterpair`](@ref) is provided to inquire the name-type pair of a parameter of a type:
```jldoctest traits
julia> parameterpair(Hi{<:Real}, 1)
Pair{:content, Real}

julia> parameterpair(Hi{<:Real}, :content)
Pair{:content, Real}
```

And a new trait function [`parameterpairs`](@ref) can be used to inquire all the name-type pairs of the parameters of a type:
```jldoctest traits
julia> parameterpairs(Hi{<:Real})
@NamedTuple{content::Real}
```

The parameters of a type can be altered all at once by giving the name-type pairs to [`fulltype`](@ref):
```jldoctest traits
julia> fulltype(Hi{Int}, NamedTuple{(:content,), Tuple{Real}}, (false,))
Hi{Real}

julia> fulltype(Hi{Int}, NamedTuple{(:content,), Tuple{Real}}, (true,))
Hi{<:Real}

julia> fulltype(Hi{Int}, NamedTuple{(:content,), Tuple{Real}})
Hi{Real}
```
Here, the last positional argument can be omitted whose default value would be determined by the [`isparameterbounds`](@ref) function which successively calls the [`isparameterbound`](@ref) function on each of the named parameter. Note that similar to the situation of the [`reparameter`](@ref) function in this subsubsection, the [`isparameterbound`](@ref) function called here is also the version that takes the parameter name as the input rather than that of the position order.

#### Type helpers with predefined contents

Julia abstract types don't have any field or attribute. They are only tags on the type tree. However, we may expect sometimes an abstract type to possess some kind of predefined content so that the design of some methods could be highly simplified. For example, we may need an abstract type that describes a composite vector. Apparently, it should have a field that is the vector contained in it. Of course, we can appoint a fixed field name with it and force every concrete subtype must contain such a field. In such a design pattern, the name of this field in every concrete subtype must be kept unchanged, which may be annoying when it conflicts with that of another field. What's worse, a predefined content of an abstract type is not always limited to a certain field. Maybe we need more than one fields to construct such a content. The just mentioned design pattern cannot deal with such situations.

Here, we provide a set of trait functions to help the design of abstract types with predefined contents. We take the case of composite vector for illustration, and the generalization to other situations is straightforward. First, the trait function [`contentnames`](@ref) should be overloaded to define the names of the predefined contents:
```jldoctest traits
julia> abstract type CompositeVector{T} end
       contentnames(::Type{<:CompositeVector}) = (:content,);
```

Then you can inquire the total number of predefined contents by [`contentcount`](@ref), inquire the name of a predefined content with its position order by [`contentname`](@ref), and judge whether a type has a predefined content with a given name by [`hascontent`](@ref):
```jldoctest traits
julia> contentcount(CompositeVector)
1

julia> contentname(CompositeVector, 1)
:content

julia> hascontent(CompositeVector, :content)
true

julia> hascontent(CompositeVector, :value)
false
```

The key is the interface [`getcontent`](@ref), which defines how to get the value of the predefined content. For the simple case when the predefined content just corresponds to a field, and also the field name of the predefined content coincides with the content name, the overloading of [`getcontent`](@ref) can be omitted, e.g.:
```jldoctest traits
julia> struct AnotherCompositeVector{T} <: CompositeVector{T}
           content::Vector{T}
       end;

julia> v = AnotherCompositeVector([1, 2, 3])
       getcontent(v, :content)
3-element Vector{Int64}:
 1
 2
 3
```

For the cases when a predefined contents does not share the same name with a certain field, or even it is not limited to only one certain field, you must implement your own [`getcontent`](@ref) manually. Let's see two typical examples:
```jldoctest traits
julia> struct DifferentFieldName{T} <: CompositeVector{T}
           data::Vector{T}
       end
       getcontent(v::DifferentFieldName, ::Val{:content}) = v.data;

julia> struct BeyondSingleField{T} <: CompositeVector{T}
           firsthalf::Vector{T}
           secondhalf::Vector{T}
       end
       getcontent(v::BeyondSingleField, ::Val{:content}) = [v.firsthalf; v.secondhalf];

julia> v = DifferentFieldName([1, 2, 3])
       getcontent(v, :content)
3-element Vector{Int64}:
 1
 2
 3

julia> v = BeyondSingleField([1, 2, 3], [4, 5, 6])
       getcontent(v, :content)
6-element Vector{Int64}:
 1
 2
 3
 4
 5
 6
```
Note that for the method overloading of [`getcontent`](@ref), the second argument is of type `Val{:content}`. This is convenient because in principle an abstract type could have more than only one predefined content, thus, the behaviors of the [`getcontent`](@ref) function could be defined separately for different predefined contents in this way. In fact, the function call `getcontent(m, contentname)` is just an alias for `getcontent(m, contentname|>Val)`.

### Holy traits

As an emergent feature of Julia, basically speaking, a Holy trait is a Julia type that could direct the generic function of a user-defined type to a certain implementation based on the Julia multi-dispatch mechanism. For different user-defined types, they could be assigned with different Holy traits, leading to different implementations of the same generic interface. Since the information of Holy traits are known at compile time, such design pattern doesn't affect the runtime efficiency as long as type stability is ensured.

#### Alternative of multi-inheritance

Maybe the most common application of Holy traits is to serve as the alternative of multi-inheritance. Let's see a simple scenario. You have defined an abstract type. It is natural to demand that for every concrete subtype of it, a pair of instances could be compared and judge whether they are equivalent to each other by the value. Unfortunately, for a new user-defined type, the default `==` function in Julia actually judges whether they are the same object, but not equal to each other by the value. Therefore, you need to define your own `==` function for this abstract type. However, you may need define a lot of abstract types when you are developing a Julia package. It is annoying if such simple functions must be written for each of them. In other languages like Python, this could be solved with the help of multi-inheritance. But Julia does not support multi-inheritance. The common way is to use Holy traits. For example, the above issue could be solved like this:
```@example traits
struct Equivalence end
const equivalence = Equivalence()
function Base.:(==)(::Equivalence, o₁, o₂)
    n₁, n₂ = fieldcount(typeof(o₁)), fieldcount(typeof(o₂))
    n₁≠n₂ && return false
    for i = 1:n₁
        getfield(o₁, i)≠getfield(o₂, i) && return false
    end
    return true
end

abstract type TypeWithEquivalence end
Base.:(==)(o₁::TypeWithEquivalence, o₂::TypeWithEquivalence) = ==(equivalence, o₁, o₂);

struct ConcreteTypeWithEquivalence{F₁, F₂} <: TypeWithEquivalence
    f₁::F₁
    f₂::F₂
end;

a₁ = ConcreteTypeWithEquivalence(("a", "b", "c"), [1, 2, 3])
a₂ = ConcreteTypeWithEquivalence(("a", "b", "c"), [1.0, 2.0, 3.0])
a₁ == a₂
```
Here, the type `Equivalence` is the Holy trait that helps the abstract type `TypeWithEquivalence` to implement the `==` function, which applies equally to any other types.

#### Type stability and the generated function trick

However, the story does not end up here. If you are concerned about the code efficiency, you may find that the above implementation is not type stable:
```@example traits
using BenchmarkTools
@benchmark $a₁ == $a₂
```
The memory allocation occurs when the `==` function tries to compare the values of `getfield(o₁, i)` and `getfield(o₂, i)` because in principle the types of these values depend on the runtime value of the variable `i`. To ensure type stability, the generated function trick can be utilized:
```@example traits
struct EfficientEquivalence end
const efficientequivalence = EfficientEquivalence()
@generated function Base.:(==)(::EfficientEquivalence, o₁, o₂)
    n₁, n₂ = fieldcount(o₁), fieldcount(o₂)
    n₁≠n₂ && return false
    expr = :(getfield(o₁, 1) == getfield(o₂, 1))
    for i = 2:n₁
        expr = Expr(:&&, expr, :(getfield(o₁, $i) == getfield(o₂, $i)))
    end
    return expr
end

abstract type TypeWithEfficientEquivalence end
function Base.:(==)(o₁::TypeWithEfficientEquivalence, o₂::TypeWithEfficientEquivalence)
    return ==(efficientequivalence, o₁, o₂)
end

struct ConcreteTypeWithEfficientEquivalence{F₁, F₂} <: TypeWithEfficientEquivalence
    f₁::F₁
    f₂::F₂
end

a₁ = ConcreteTypeWithEfficientEquivalence(("a", "b", "c"), [1, 2, 3])
a₂ = ConcreteTypeWithEfficientEquivalence(("a", "b", "c"), [1.0, 2.0, 3.0])
a₁ == a₂
```
```@example traits
@benchmark $a₁ == $a₂
```
At runtime of the generated `==` function, it compares the values of `getfield(o₁, 1)` and `getfield(o₂, 1)`, `getfield(o₁, 2)` and `getfield(o₂, 2)`, etc., whose types are known at compile time. Therefore, type stability could be ensured.

#### EfficientOperations

`EfficientOperations` is a Holy trait defined in this module that packs several common operations, such as `==/isequal`, `</isless`, `isapprox` and `replace`, to help other (abstract) types to implement such functions by passing [`efficientoperations`](@ref) as the first argument, just as illustrated above. See the manual for more detailed information.

### Manual

For traits with types themselves:
```@docs
commontype
fulltype
rawtype
DataType
supertype
```

For traits with type parameters:
```@docs
hasparameter
isparameterbound
isparameterbounds
parametercount
parametername
parameternames
parameterorder
parameterpair
parameterpairs
parametertype
parametertypes
promoteparameters
reparameter
```

For traits with type contents:
```@docs
contentcount
contentname
contentnames
contentorder
contenttype
contenttypes
getcontent
hascontent
dissolve
```

For traits with type operations:
```@docs
efficientoperations
```

## Composite structures

In principle, Julia is not an object-oriented programming language. For example, only abstract types can be inherited so that subtype cannot inherit fields from their parents. Therefore, Julia prefers composition over inheritance. However, to make a new concrete type behaves much alike another one, tedious repetitions of redefining the generic interfaces are usually not avoidable, especially for the basic types in Julia base. In this module, we implement to such composited types, [`CompositeVector`](@ref) and [`CompositeDict`](@ref), for the sake of future usages.

### CompositeVector

A composite vector can be considered as a vector that is implemented by including a concrete subtype of [`AbstractVector`](https://docs.julialang.org/en/v1/base/arrays/#Base.AbstractVector) as its data attribute, and it itself is a subtype of [`AbstractVector`](https://docs.julialang.org/en/v1/base/arrays/#Base.AbstractVector).

To take full advantages of the Julia base, the following interfaces are redefined:
* inquiry of info: `size`, `length`
* comparison between objects: `==`, `isequal`
* obtainment of old elements: `getindex`
* operation of old elements: `setindex!`
* addition of new elements: `push!`, `pushfirst!`, `insert!`, `append!`, `prepend!`
* removal of old elements: `splice!`, `deleteat!`, `pop!`, `popfirst!`, `empty!`
* construction of new objects: `empty`, `reverse`, `similar`
* iteration: `iterate`, `keys`, `values`, `pairs`

Note that arithmetic operations and logical operations excluding `==` and `isequal` are not supported.

### CompositeDict

A composite dict can be considered as a dict that is implemented by including a concrete subtype of [`AbstractDict`](https://docs.julialang.org/en/v1/base/collections/#Base.AbstractDict) as its data attribute and it itself is a subtype of [`AbstractDict`](https://docs.julialang.org/en/v1/base/collections/#Base.AbstractDict).

To take full advantages of the Julia base, the following interfaces are redefined:
* inquiry of info: `isempty`, `length`, `haskey`, `in`
* comparison between objects: `==`, `isequal`
* obtainment of old elements: `get`, `getkey`, `getindex`
* operation and addition of elements: `push!`, `get!`, `setindex!`
* removal of old elements: `pop!`, `delete!`, `empty!`
* construction of new objects: `merge`, `empty`
* iteration: `iterate`, `keys`, `values`, `pairs`

### Manual

```@docs
CompositeDict
CompositeVector
```

## Vector spaces

A [vector space](https://en.wikipedia.org/wiki/Vector_space) is a linear space, in which the addition of vectors and multiplication of a vector by a scalar are defined.

Vector spaces are frequently encountered in physics, e.g. the Hilbert space in quantum mechanics. In this submodule, we only implement those with finite dimensions. We want to remark that in our implementation, a vector space is a subtype of an abstract vector, therefore, the bases always possess a order, which means, two vector spaces are not considered to be equal to each other even if their corresponding actual mathematical spaces are the same but the orders of the bases are different.

### VectorSpace

[`VectorSpace{B}`](@ref) is the abstraction of a vector space, which has only one type parameter:
* `B<:Any`: the type of the bases of the vector space

Basically, a subtype should implement the following 2 methods:
1) ```julia
   Base.length(vs::VectorSpace) -> Int
   ```
   Get the dimension of a vector space
2) ```julia
   Base.getindex(vs::VectorSpace{B}, i::Int)  where B -> B
   ```
   Get the ith basis of a vector space

Other features include
* comparison: `==` and `isequal`
* iteration: `iterate`
* inquiry: `size` and `in`
* search: `searchsortedfirst`

### Manual

Predefined types of vector spaces:
```@docs
VectorSpace
NamedVectorSpace
SimpleNamedVectorSpace
ParameterSpace
NamedVectorSpaceProd
NamedVectorSpaceZip
```

Predefined types of vector space style:
```@docs
VectorSpaceStyle
VectorSpaceGeneral
VectorSpaceEnumerative
VectorSpaceCartesian
VectorSpaceDirectProducted
VectorSpaceDirectSummed
VectorSpaceZipped
```
