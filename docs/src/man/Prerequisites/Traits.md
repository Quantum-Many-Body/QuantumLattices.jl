```@meta
CurrentModule = QuantumLattices.Prerequisites.Traits
DocTestFilters = [r"var\".*\"", r"generic function with [0-9]* method", r".*s \(.*\% GC\)", r"evals/sample:.*"]
DocTestSetup = quote
    push!(LOAD_PATH, "../../../../src/")
    using QuantumLattices.Prerequisites.Traits
    import QuantumLattices.Prerequisites.Traits: isparameterbound, parameternames, contentnames, getcontent
end
```

# Traits

This module defines some trait functions and trait types that are useful to the package.

Generally speaking, traits in Julia could fall into two categories according to their usages, the first may be term as "type helpers" and the second are usually called "Holy traits" named after [Tim Holy](https://github.com/timholy). Type helpers aim at the inquiry, alteration and computation of the compile-time information of types, while Holy traits can be applied as an alternative to multi-inheritance by use of the Julia multidispatch feature.

## Type helpers

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

### Type helpers with type parameters

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
ERROR: MethodError: no method matching contenttype(::Type{Hi{var"#s1"} where var"#s1"<:Real})
[...]
```

In fact in Julia base, all such inquiry functions, e.g., the `eltype` function, work poor for the `UnionAll` types:
```jldoctest traits
julia> eltype(Vector{<:Real})
Any
```
In concept, `eltype(Vector{<:Real}` should return `Real` instead of `Any` as every element in `Vector{<:Real}` is indeed a real number. Similarly, we expect that `contenttype(Hi{<:Real})` should also give us `Real`. Unfortunately, functions defined in the similar form like this could never achieve such goals. Julia base doesn't provide generic functions to access or change the information of the parameters of a type. In this module, we try to fill this gap with a set of generic trait functions.

#### Access or change the type parameters by their position orders

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
Hi{Real}

julia> reparameter(Vector{Int64}, 1, Real)
Vector{Real} (alias for Array{Real, 1})

julia> reparameter(Vector{<:Real}, 2, 3)
Array{var"#s1",3} where var"#s1"<:Real

julia> reparameter(Hi{Int64}, 1, Real, false)
Hi{Real}

julia> reparameter(Hi{Int64}, 1, Real, true)
Hi{var"#s2"} where var"#s2"<:Real
```
We want to remark that by providing the fourth positional argument with the `true` value, a `UnionAll` type could be generated. When the fourth positional argument is omitted, it is actually determined by another trait function, i.e., [`isparameterbound`](@ref). This function judges whether an input type should be considered as the upper bound of the new parameter of a type. By default, it is always defined to be `false`. This function can be overloaded to change the behavior for a certain type:
```jldoctest traits
julia> isparameterbound(::Type{<:Hi}, ::Val{1}, D) = !isconcretetype(D);

julia> reparameter(Hi, 1, Real)
Hi{var"#s3"} where var"#s3"<:Real
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
Hi{var"#s4"} where var"s#4"<:Real

julia> fulltype(Hi{Int64}, Tuple{Real}, (false,))
Hi{Real}

julia> fulltype(Vector{Int64}, Tuple{Real, 2})
Matrix{Real} (alias for Array{Real, 2})

julia> fulltype(Vector{Int64}, Tuple{Real, 2}, (true, false))
Matrix{var"#s5"} where var"#s5"<:Real (alias for Array{var"#s5", 2} where var"#s5"<:Real)
```
Like [`reparameter`](@ref), the last positional argument of [`fulltype`](@ref) could determine whether the corresponding types specified by the type parameters of the input `Tuple` should be considered as the upper bounds of the new parameters of a type. When this argument is omitted, it is determined by another trait function [`isparameterbounds`](@ref), which successively calls the [`isparameterbound`](@ref) function to determine the behaviors for all the parameters of a type as the literal indicates.

#### Associate type parameters with names

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
Hi{Real}

julia> reparameter(Hi{Int}, :content, Real, true)
Hi{var"#s6"} where var"#s6"<:Real
```

To change the [`reparameter`](@ref) behavior when its last positional argument is omitted, you should overload the [`isparameterbound`](@ref) function accordingly, e.g.:
```jldoctest traits
julia> isparameterbound(::Type{<:Hi}, ::Val{:content}, D) = !isconcretetype(D);

julia> reparameter(Hi{Int}, :content, Real)
Hi{var"#s7"} where var"#s7"<:Real
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
NamedTuple{(:content,), Tuple{Real}}
```

The parameters of a type can be altered all at once by giving the name-type pairs to [`fulltype`](@ref):
```jldoctest traits
julia> fulltype(Hi{Int}, NamedTuple{(:content,), Tuple{Real}}, (false,))
Hi{Real}

julia> fulltype(Hi{Int}, NamedTuple{(:content,), Tuple{Real}}, (true,))
Hi{var"#s8"} where var"#s8"<:Real

julia> fulltype(Hi{Int}, NamedTuple{(:content,),Tuple{Real}})
Hi{var"#s9"} where var"#s9"<:Real
```
Here, the last positional argument can be omitted whose default value would be determined by the [`isparameterbounds`](@ref) function which successively calls the [`isparameterbound`](@ref) function on each of the named parameter. Note that similar to the situation of the [`reparameter`](@ref) function in this subsection, the [`isparameterbound`](@ref) function called here is also the version that takes the parameter name as the input rather than that of the position order.

### Type helpers with predefined contents

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

## Holy traits

As an emergent feature of Julia, basically speaking, a Holy trait is a Julia type that could direct the generic function of a user-defined type to a certain implementation based on the Julia multi-dispatch mechanism. For different user-defined types, they could be assigned with different Holy traits, leading to different implementations of the same generic interface. Since the information of Holy traits are known at compile time, such design pattern doesn't affect the runtime efficiency as long as type stability is ensured.

### Alternative of multi-inheritance

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

### Type stability and the generated function trick

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

### EfficientOperations

`EfficientOperations` is a Holy trait defined in this module that packs several common operations, such as `==/isequal`, `</isless`, `isapprox` and `replace`, to help other (abstract) types to implement such functions by passing [`efficientoperations`](@ref) as the first argument, just as illustrated above. See the manual for more detailed information.

## Manual

```@autodocs
Modules = [Traits]
Order = [:module, :constant, :type, :macro, :function]
```

```@meta
DocTestSetup = nothing
```
