```@meta
CurrentModule = QuantumLattices.Prerequisites.Traits
DocTestFilters = [r"var\".*\"", r"generic function with [0-9]* method"]
DocTestSetup = quote
    push!(LOAD_PATH, "../../../../src/")
    using QuantumLattices.Prerequisites.Traits
    import QuantumLattices.Prerequisites.Traits: isparameterbound, parameternames, contentnames, fieldnameofcontent
end
```

# Traits

This module defines some trait functions and trait types that are useful to the package.

In general, traits in Julia fall into two categories according to their usages, the first I term as "type helpers" and the second called "Holy traits" named after [Tim Holy](https://github.com/timholy). Type helpers are aimed at the inquiry, alteration and computation of the compile-time information of types, while Holy traits can be applied as an alternative to multi-inheritance to enhance polymorphism.

## Type helpers

Type helpers are important for the generic programming in Julia, especially in the design of generic interfaces and abstract types.

Let's see a simple situation, i.e. the elemental addition of two vectors of numbers. The numbers can assume different types and the type of the result depends on both of them, for example, the result between two vectors of integers is a vector of integers while that between a vector of integers and a vector of floats is a vector of floats. Of course, one can explicitly define every elemental addition function between any two different types of vectors of numbers, like this:
```julia
# wrong design pattern

function elementaladdition(v1::Vector{Int64}, v2::Vector{Int64})
    result = Int[]
    ...
end
function elementaladdition(v1::Vector{Int64}, v2::Vector{Float64})
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

function elementaladdition(v1::Vector{T1}, v2::Vector{T2}) where {T1<:Number, T2<:Number}
    result = promotion(T1, T2)[]
    ...
end
function elementalsubtraction(v1::Vector{T1}, v2::Vector{T2}) where {T1<:Number, T2<:Number}
    result = promotion(T1, T2)[]
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
Array{Real,1}

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
The second positional argument must be of type `Val` because in principle you should be able to assign different behaviors for different parameters of a type separately. If it is of type `Integer`, a single overloading would change the behaviors for all.

Besides, you can inquire all the parameters of a type by [`parametertypes`](@ref):
```jldoctest traits
julia> parametertypes(Hi{<:Real})
Tuple{Real}

julia> parametertypes(Vector{Int64})
Tuple{Int64,1}
```
The obtained type parameters are stored in those of a `Tuple`.

At the same time, you can change all the parameters of a type by [`fulltype`](@ref):
```jldoctest traits
julia> fulltype(Hi{Int64}, Tuple{Real})
Hi{var"#s4"} where var"s#4"<:Real

julia> fulltype(Hi{Int64}, Tuple{Real}, (false,))
Hi{Real}

julia> fulltype(Vector{Int64},Tuple{Real, 2})
Array{Real,2}

julia> fulltype(Vector{Int64}, Tuple{Real, 2}, (true, false))
Array{var"#s5",2} where var"#s5"<:Real
```
Like [`reparameter`](@ref), the last positional argument of [`fulltype`](@ref) could determine whether the corresponding types specified by the type parameters of the input `Tuple` should be considered as the upper bounds of the new parameters of a type. When this argument is omitted, it is determined by another trait function [`isparameterbounds`](@ref), which successively calls the [`isparameterbound`](@ref) function to determine the behaviors for all the parameters of a type as the literal indicates.

#### Associate type parameters with names

Sometimes, it is more convenient to associate names with the parameters of a type, and then access or change them by their names. This can be done by overloading the [`parameternames`](@ref) trait function for a certain type:
```jldoctest traits
julia> parameternames(::Type{<:Hi}) = (:content,)
parameternames (generic function with 3 methods)
```

Now, you can inquire the parameter name by [`parametername`](@ref) with the given position order or vice versa by [`parameterorder`](@ref):
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

And [`parametertype`](@ref) and [`reparameter`](@ref) can be applied by the parameter name instead of the position order:
```jldoctest traits
julia> parametertype(Hi{<:Real}, :content)
Real

julia> reparameter(Hi{Int}, :content, Real)
Hi{Real}

julia> reparameter(Hi{Int}, :content, Real, true)
Hi{var"#s6"} where var"#s6"<:Real
```

To change the [`reparameter`](@ref) behavior when its last positional argument is omitted, you should overload the [`isparameterbound`](@ref) function, e.g.:
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
    To change the default behavior of the former function call, you must overload the following method
    ```julia
    isparameterbound(::Type{<:Hi}, ::Val{:content}, D)
    ```

A new trait function [`parameterpair`](@ref) is provided to inquire the name-type pair of a parameter of a type:
```jldoctest traits
julia> parameterpair(Hi{<:Real}, 1)
Pair{:content,Real}

julia> parameterpair(Hi{<:Real}, :content)
Pair{:content,Real}
```

And a new trait function [`parameterpairs`](@ref) can be used to inquire all the name-type pairs of the parameters of a type:
```jldoctest traits
julia> parameterpairs(Hi{<:Real})
NamedTuple{(:content,),Tuple{Real}}
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

Julia abstract types don't have any field or attribute. They are only tags on the type tree. However, we may expect sometimes an abstract type to possess some kind of predefined content so that the design of some methods could be highly simplified. For example, we may want an abstract type that describes a composite vector. Apparently, it should have a field that is the vector contained in it. Of course, we can appoint a fixed field name with it and force every concrete type must contain such a field. In such a design pattern, the name of this field in every concrete type must be kept unchanged, which may be annoying when it conflicts with that of another field. On the other hand, a predefined content of an abstract type is not always limited to a certain field. Maybe we need more than one fields to construct such a content. The just mentioned design pattern cannot deal with such situations.

Here, we provide a set of trait functions to help the design of abstract types with predefined contents. We take the case of composite vector for illustration, and the generalization to other situations is straightforward. First, the trait function [`contentnames`](@ref) should be overloaded to define the names of the predefined contents:
```jldoctest traits
julia> abstract type CompositeVector{T} end
       contentnames(::Type{<:CompositeVector}) = (:content,);
```

Then you can inquire the total number of predefined contents by [`contentcount`](@ref), inquire the name of a predefined content with its position order by [`contentname`](@ref), and judge whether a type has a predefined content with a name by [`hascontent`](@ref):
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

The key is the interface [`getcontent`](@ref), which defines how to get the value of the predefined content. For the simple case when the predefined content just corresponds to a field, you can use the trait function [`fieldnameofcontent`](@ref) to specify the field name of the content, then [`getcontent`](@ref) gets a  default implementation, and you need not overload it explicitly, e.g.:
```jldoctest traits
julia> struct ConcreteCompositeVector{T} <: CompositeVector{T}
           vector::Vector{T}
       end
       fieldnameofcontent(::Type{<:ConcreteCompositeVector}, ::Val{:content}) = :vector;

julia> v = ConcreteCompositeVector([1, 2, 3])
       getcontent(v, :content)
3-element Array{Int64,1}:
 1
 2
 3
```

If the field name of the predefined content coincides with the content name, the overloading of [`fieldnameofcontent`](@ref) can also be omitted, e.g.:
```jldoctest traits
julia> struct AnotherCompositeVector{T} <: CompositeVector{T}
           content::Vector{T}
       end;

julia> v = AnotherCompositeVector([1, 2, 3])
       getcontent(v, :content)
3-element Array{Int64,1}:
 1
 2
 3
```

For the case that a predefined content is not limited to a certain field, you must implement your own [`getcontent`](@ref) manually.

## Holy traits

As an emergent feature of Julia, basically speaking, a Holy trait is a Julia type that could direct the generic function of a user-defined type to a certain implementation based on the Julia multi-dispatch mechanism. For different user-defined types, they could be assigned with different Holy traits, leading to different implementations of the same generic interface. Since the information of Holy traits are known at compile time, such design pattern doesn't affect the runtime efficiency as long as type stability is ensured.

### EfficientOperations

`EfficientOperations` defines efficient operations such as `==/isequal`, `</isless`, `isapprox`, `replace`, etc., that ensure type stability.

Type stability is the key of Julia to improve the code efficiency. However, it cannot be ensured in some unexpected cases, especially where an iterator is involved. For example, the following codes appears type unstable:
```julia
function Base.:(==)(o1::AbstractType, o2::AbstractType)
    n1, n2 = fieldcount(typeof(o1)), fieldcount(typeof(o2))
    (n1 == n2) ? all(getfield(o1, i) == getfield(o2, i) for i = 1:n1) : false
end
```
Methods like above are common when we design abstract types, but they are not type stable. To get rid of it, the generated function trick can be used:
```julia
@generated function Base.:(==)(o1::AbstractType, o2::AbstractType)
    n1, n2 = fieldcount(o1), fieldcount(o2)
    if n1 == n2
        expr=:(getfield(o1, 1) == getfield(o2, 1))
        for i = 2:n1
            expr = Expr(:&&, expr, :(getfield(o1, $i) == getfield(o2, $i)))
        end
        return expr
    else
        return :(false)
    end
end
```
Then type stability can be ensured. We use this trick to implement the methods such as `==/isequal`, `</isless`, `isapprox`, `replace`, etc., with the trait `efficientoperations::EfficientOperations`. Other types can resort to these methods by passing [`efficientoperations`](@ref) as the first argument.

### MemoryOrder

`MemoryOrder` provides the conversions, [`subtoind`](@ref) and [`indtosub`](@ref), between a Cartesian index represented by a tuple and a linear index represented by an integer. C/C++ order or Fortran order can be specified by use of the constants [`corder`](@ref) or [`forder`](@ref), which are instances of singleton types `COrder` and `FOrder` that are both subtypes of the abstract type `MemoryOrder`.

## Manual

```@autodocs
Modules = [Traits]
Order = [:module, :constant, :type, :macro, :function]
```

```@meta
DocTestSetup = nothing
```