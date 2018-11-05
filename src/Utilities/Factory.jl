module Factory

using MacroTools: block
using MacroTools: splitarg,combinefield,combinearg
using MacroTools: splitstructdef,combinestructdef
using MacroTools: splitdef,combinedef
using Printf: @printf

export Block,@block,@push!
export Argument,@argument
export FunctionFactory,@functionfactory,addargs!,@addargs!,addkwargs!,@addkwargs!,extendbody!,@extendbody!
export Field,@field
export TypeFactory,@typefactory,addfields!,@addfields!,addconstructors!,@addconstructors!

"Factory expression types."
const FExpr=Union{Symbol,Expr}

"""
    AbstractFactory

Abstract type for all concrete factories.
"""
abstract type AbstractFactory end

"""
    ==(f1::F,f2::F) where F<:AbstractFactory

Overloaded `==` operator.
"""
Base.:(==)(f1::F,f2::F) where F<:AbstractFactory=all(getfield(f1,name)==getfield(f2,name) for name in F|>fieldnames)

"""
    show(io::IO,f::AbstractFactory)

Show a concrete `AbstractFactory`.
"""
Base.show(io::IO,f::AbstractFactory)=@printf io "%s(\n%s\n)" f|>typeof|>Base.typename join(("$name: $(getfield(f,name))" for name in f|>typeof|>fieldnames),'\n')

"""
    replace(f::AbstractFactory;kwargs...)

Return a copy of a concrete `AbstractFactory` with some of the field values replaced by the keyword arguments.
"""
Base.replace(f::AbstractFactory;kwargs...)=(f|>typeof)((get(kwargs,key,getfield(f,name)) for name in f|>typeof|>fieldnames)...)

"""
    Block(parts::FExpr...)

The struct to describe a `begin ... end` block.
"""
struct Block <: AbstractFactory
    body::Vector{Any}
    Block(parts::FExpr...)=new(vcat((block(part).args for part in parts)...))
end

"""
    @block parts::FExpr...

Construct a `Block` directly from a `begin ... end` block definition.
"""
macro block(parts::FExpr...) :(Block($parts...)) end

"""
    (b::Block)()

Convert a `Block` to the `Expr` representation of the `begin ... end` block it describes.
"""
(b::Block)()=Expr(:block,b.body...)

"""
    push!(b::Block,parts::FExpr...)
    push!(b::Block,parts::Block...)

Push other parts into the body of a block.
"""
Base.push!(b::Block,parts::FExpr...)=push!(b,Block.(parts)...)
function Base.push!(b::Block,parts::Block...)
    for part in parts
        append!(b.body,part.body)
    end
end

"""
    @push! b parts::FExpr...

Push other parts into the body of a block.
"""
macro push!(b,parts::FExpr...) :(push!($(esc(b)),$parts...)) end

"""
    Argument(name::Symbol,type::FExpr,slurp::Bool,default::Any)
    Argument(name::Symbol;type::FExpr=:Any,slurp::Bool=false,default::Any=nothing)
    Argument(expr::Expr)

The struct to describe a argument of a `function`.
"""
struct Argument <: AbstractFactory
    name::Symbol
    type::FExpr
    slurp::Bool
    default::Any
    Argument(name::Symbol,type::FExpr,slurp::Bool,default::Any)=new(name,type,slurp,default)
    Argument(name::Symbol;type::FExpr=:Any,slurp::Bool=false,default::Any=nothing)=new(name,type,slurp,default)
    Argument(expr::Expr)=(cache=splitarg(expr);new(cache[1],cache[2],cache[3],cache[4]))
end

"""
    @argument expr::FExpr

Construct an `Argument` directly from an argument statement.
"""
macro argument(expr::FExpr)
    expr=[expr]
    :(Argument($expr...))
end

"""
    (a::Argument)()

Convert a `Argument` to the `Expr` representation of the argument it describes.
"""
(a::Argument)()=combinearg((getfield(a,name) for name in a|>typeof|>fieldnames)...)

"""
    FunctionFactory(name::Symbol,args::Vector,kwargs::Vector,rtype::FExpr,params::Vector,body::FExpr)
    FunctionFactory(name::Symbol;args::Vector=Any[],kwargs::Vector=Any[],rtype::FExpr=:Any,params::Vector=Any[],body::FExpr)
    FunctionFactory(expr::Expr)

The struct to describe a `function`.
"""
struct FunctionFactory <: AbstractFactory
    name::Symbol
    args::Vector{FExpr}
    kwargs::Vector{FExpr}
    rtype::FExpr
    params::Vector{FExpr}
    body::Expr
    FunctionFactory(name::Symbol,args::Vector,kwargs::Vector,rtype::FExpr,params::Vector,body::FExpr)=new(name,args,kwargs,rtype,params,block(body))
    function FunctionFactory(name::Symbol;args::Vector=Any[],kwargs::Vector=Any[],rtype::FExpr=:Any,params::Vector=Any[],body::FExpr)
        new(name,args,kwargs,rtype,params,block(body))
    end
    function FunctionFactory(expr::Expr)
        dict=splitdef(expr)
        FunctionFactory(dict[:name],dict[:args],dict[:kwargs],get(dict,:rtype,:Any),[v for v in dict[:whereparams]],dict[:body])
    end
end

"""
    @functionfactory expr::Expr

Construct a `FunctionFactory` directly from a function definition.
"""
macro functionfactory(expr::Expr)
    expr=[expr]
    :(FunctionFactory($expr...))
end

"""
    (ff::FunctionFactory)()

Convert a `FunctionFactory` to the `Expr` representation of the `function` it describes.
"""
(ff::FunctionFactory)()=combinedef(Dict(name=>getfield(ff,name) for name in ff|>typeof|>fieldnames))

"""
    addargs!(ff::FunctionFactory,args::FExpr...)
    addargs!(ff::FunctionFactory,args::Argument...)

Add a couple of positional arguments to a function factory.
"""
addargs!(ff::FunctionFactory,args::FExpr...)=addargs!(ff,Argument.(args)...)
addargs!(ff::FunctionFactory,args::Argument...)=push!(ff.args,(arg() for arg in args)...)

"""
    @addargs! ff args::FExpr...

Add a couple of positional arguments to a function factory.
"""
macro addargs!(ff,args::FExpr...) :(addargs!($(esc(ff)),$args...)) end

"""
    addkwargs!(ff::FunctionFactory,kwargs::FExpr...)
    addkwargs!(ff::FunctionFactory,kwargs::Argument...)

Add a couple of keyword arguments to a function factory.
"""
addkwargs!(ff::FunctionFactory,kwargs::FExpr...)=addkwargs!(ff,Argument.(kwargs)...)
addkwargs!(ff::FunctionFactory,kwargs::Argument...)=push!(ff.kwargs,(kwarg() for kwarg in kwargs)...)

"""
    @addkwargs! ff kwargs::FExpr...

Add a couple of keyword arguments to a function factory.
"""
macro addkwargs!(ff,kwargs::FExpr...) :(addkwargs!($(esc(ff)),$kwargs...)) end

"""
    extendbody!(ff::FunctionFactory,parts::FExpr...)
    extendbody!(ff::FunctionFactory,parts::Block...)

Extend the body of a function factory.
"""
extendbody!(ff::FunctionFactory,parts::FExpr...)=extendbody!(ff,Block.(parts)...)
function extendbody!(ff::FunctionFactory,parts::Block...)
    for part in parts
        append!(ff.body.args,part.body)
    end
end

"""
    @extendbody! ff parts::FExpr...

Extend the body of a function factory.
"""
macro extendbody!(ff,parts::FExpr...) :(extendbody!($(esc(ff)),$parts...)) end

"""
    Field(name::Symbol,type::FExpr)
    Field(name::Symbol;type::FExpr=:Any)
    Field(expr::Expr)

The struct to describe a field of a `struct`.
"""
struct Field <: AbstractFactory
    name::Symbol
    type::FExpr
    Field(name::Symbol,type::FExpr)=new(name,type)
    Field(name::Symbol;type::FExpr=:Any)=new(name,type)
    Field(expr::Expr)=(cache=splitarg(expr);new(cache[1],cache[2]))
end

"""
    @field expr::FExpr

Construct a `Field` directly from a field statement.
"""
macro field(expr::FExpr)
    expr=[expr]
    :(Field($expr...))
end

"""
    (f::Field)()

Convert a `Field` to the `Expr` representation of the field it describes.
"""
(f::Field)()=combinefield((f.name,f.type))

"""
    TypeFactory(name::Symbol,mutable::Bool,params::Vector,supertype::FExpr,fields::Vector,constructors::Vector)
    TypeFactory(name::Symbol;mutable::Bool=false,params::Vector=Any[],supertype::FExpr=:Any,fields::Vector=Any[],constructors::Vector=Any[])
    TypeFactory(expr::Expr)

The struct to describe a `struct`.
"""
struct TypeFactory <: AbstractFactory
    name::Symbol
    mutable::Bool
    params::Vector{FExpr}
    supertype::FExpr
    fields::Vector{Tuple{Symbol,FExpr}}
    constructors::Vector{Expr}
    TypeFactory(name::Symbol,mutable::Bool,params::Vector,supertype::FExpr,fields::Vector,constructors::Vector)=new(name,mutable,params,supertype,fields,constructors)
    function TypeFactory(name::Symbol;mutable::Bool=false,params::Vector=Any[],supertype::FExpr=:Any,fields::Vector=Any[],constructors::Vector=Any[])
        new(mutable,name,params,supertype,fields,constructors)
    end
    function TypeFactory(expr::Expr)
        dict=splitstructdef(expr)
        new(dict[:name],dict[:mutable],dict[:params],dict[:supertype],dict[:fields],dict[:constructors])
    end
end

"""
    @typefactory expr::Expr

Construct a `TypeFactory` directly from a type definition.
"""
macro typefactory(expr::Expr)
    expr=[expr]
    :(TypeFactory($expr...))
end

"""
    (tf::TypeFactory)()

Convert a `TypeFactory` to the `Expr` representation of the `struct` it describes.
"""
(tf::TypeFactory)()=combinestructdef(Dict(name=>getfield(tf,name) for name in tf|>typeof|>fieldnames))

"""
    addfields!(rf::TypeFactory,fields::FExpr...)
    addfields!(rf::TypeFactory,fields::Field...)

Add a couple of fields to a type factory.
"""
addfields!(rf::TypeFactory,fields::FExpr...)=addfields!(rf,Field.(fields)...)
addfields!(rf::TypeFactory,fields::Field...)=push!(rf.fields,((field.name,field.type) for field in fields)...)

"""
    @addfields! rf fields::FExpr...

Add a couple of fields to a type factory.
"""
macro addfields!(rf,fields::FExpr...) :(addfields!($(esc(rf)),$fields...)) end

"""
    addconstructors!(rf::TypeFactory,constructors::Expr...)
    addconstructors!(rf::TypeFactory,constructors::FunctionFactory...)

Add a couple of constructors to a type factory.
"""
addconstructors!(rf::TypeFactory,constructors::Expr...)=addconstructors!(rf,FunctionFactory.(constructors)...)
addconstructors!(rf::TypeFactory,constructors::FunctionFactory...)=push!(rf.constructors,(constructor() for constructor in constructors)...)

"""
    @addconstructors! rf constructors::Expr...

Add a couple of constructors to a type factory.
"""
macro addconstructors!(rf,constructors::Expr...) :(addconstructors!($(esc(rf)),$constructors...)) end

end #module
