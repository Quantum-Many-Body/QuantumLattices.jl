module Factory

using MacroTools: block
using MacroTools: splitarg,combinefield,combinearg
using MacroTools: splitstructdef,combinestructdef
using MacroTools: splitdef,combinedef
using Printf: @printf,@sprintf
import MacroTools: rmlines

export Argument,@argument
export Parameter,@parameter
export Field,@field
export Block,@block,@push!,rmlines!,@rmlines!,rmlines,@rmlines
export FunctionFactory,@functionfactory,addargs!,@addargs!,addkwargs!,@addkwargs!,extendbody!,@extendbody!
export TypeFactory,@typefactory,addfields!,@addfields!,addconstructors!,@addconstructors!
export addparams!,@addparams!

@generated function maxfieldnamelength(::T) where T
    result=1
    for name in T|>fieldnames
        result=max(result,name|>string|>length)
    end
    return :($result)
end

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
function Base.show(io::IO,f::AbstractFactory)
    @printf io "%s(\n" f|>typeof|>nameof
    for name in f|>typeof|>fieldnames
        nstr=" "^2*string(name)*":"*" "^((f|>maxfieldnamelength)+1-(name|>string|>length))
        value=getfield(f,name)
        if isa(value,Nothing)
            vstr="nothing"
        elseif isa(value,AbstractFactory)
            vstr=string(value())
        elseif isa(value,Vector{<:AbstractFactory})
            eltypename=value|>eltype|>nameof
            strings=Tuple(string(af()) for af in value)
            if (length(strings)>0 ? sum(map(length,strings)) : 0)<90-(f|>maxfieldnamelength)
                vstr=@sprintf "%s[%s]" eltypename join(strings,", ")
            else
                vstr=replace((@sprintf "%s[\n%s\n]" eltypename join(strings,"\n")),"\n"=>"\n  ")
            end
        else
            vstr=string(value)
        end
        '\n' in vstr && (vstr=replace(vstr,"\n"=>"\n"*" "^((f|>maxfieldnamelength)+4)))
        @printf(io,"%s%s\n",nstr,vstr)
    end
    @printf io ")"
end

"""
    replace(f::AbstractFactory;kwargs...)

Return a copy of a concrete `AbstractFactory` with some of the field values replaced by the keyword arguments.
"""
Base.replace(f::AbstractFactory;kwargs...)=(f|>typeof)((get(kwargs,key,getfield(f,name)) for name in f|>typeof|>fieldnames)...)

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
macro argument(expr::FExpr) expr=[expr];:(Argument($expr...)) end

"""
    (a::Argument)()

Convert an `Argument` to the `Expr` representation of the argument it describes.
"""
(a::Argument)()=combinearg((getfield(a,name) for name in a|>typeof|>fieldnames)...)

"""
    Parameter(name::Symbol,type::FExpr)
    Parameter(name::Symbol;type::FExpr=:Any)
    Parameter(expr::FExpr)

The struct to describe a parameter of a `function` or a `type`.
"""
struct Parameter <: AbstractFactory
    name::Symbol
    type::FExpr
    Parameter(name::Symbol,type::FExpr)=new(name,type)
    Parameter(name::Symbol;type::FExpr=:nothing)=new(name,type)
    function Parameter(expr::FExpr)
        if isa(expr,Symbol)
            name,type=expr,:nothing
        else
            @assert expr.head==:(<:) "Parameter error: wrong formed input expression."
            (name,type)=length(expr.args)==2 ? expr.args : (:nothing,expr.args[1])
        end
        new(name,type)
    end
end

"""
    @parameter expr::FExpr

Construct a `Parameter` directly from an parameter statement.
"""
macro parameter(expr::FExpr) expr=[expr];:(Parameter($expr...)) end

"""
    (p::Parameter)()

Convert a `Parameter` to the `Expr` representation of the parameter it describes.
"""
(p::Parameter)()=p.name==:nothing ? :(<:$(p.type)) : p.type==:nothing ? :($(p.name)) : :($(p.name)<:$(p.type))

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
macro field(expr::FExpr) expr=[expr];:(Field($expr...)) end

"""
    (f::Field)()

Convert a `Field` to the `Expr` representation of the field it describes.
"""
(f::Field)()=combinefield((f.name,f.type))

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
    rmlines!(b::Block)

Remove line number nodes in the body of a block.
"""
rmlines!(b::Block)=filter!(part->!isa(part,LineNumberNode),b.body)

"""
    @rmlines! b::Expr

Remove line number nodes in the body of a block.
"""
macro rmlines!(b::Expr) :(rmlines!($(esc(b)))) end

"""
    rmlines(b::Block)

Return a copy of a block with the line number nodes removed.
"""
rmlines(b::Block)=Block(filter(part->!isa(part,LineNumberNode),b.body)...)

"""
    @rmlines b::Expr

Return a copy of a block with the line number nodes removed.
"""
macro rmlines(b::Expr) :(rmlines($(esc(b)))) end

"""
    FunctionFactory(name::Symbol,args::Vector{Argument},kwargs::Vector{Argument},rtype::FExpr,params::Vector{Parameter},body::Block)
    FunctionFactory(    name::Symbol;
                        args::Vector{Argument}=Argument[],
                        kwargs::Vector{Argument}=Argument[],
                        rtype::FExpr=:Any,
                        params::Vector{Parameter}=Parameter[],
                        body::Block=Block()
                        )
    FunctionFactory(expr::Expr)

The struct to describe a `function`.
"""
struct FunctionFactory <: AbstractFactory
    name::Symbol
    args::Vector{Argument}
    kwargs::Vector{Argument}
    rtype::FExpr
    params::Vector{Parameter}
    body::Block
    function FunctionFactory(name::Symbol,args::Vector{Argument},kwargs::Vector{Argument},rtype::FExpr,params::Vector{Parameter},body::Block)
        new(name,args,kwargs,rtype,params,body)
    end
    function FunctionFactory(   name::Symbol;
                                args::Vector{Argument}=Argument[],
                                kwargs::Vector{Argument}=Argument[],
                                rtype::FExpr=:Any,
                                params::Vector{Parameter}=Parameter[],
                                body::Block=Block()
                                )
        new(name,args,kwargs,rtype,params,body)
    end
    function FunctionFactory(expr::Expr)
        dict=splitdef(expr)
        new(    dict[:name],
                Argument.(dict[:args]),
                Argument.(dict[:kwargs]),
                get(dict,:rtype,:Any),
                Parameter[Parameter(param) for param in dict[:whereparams]],
                Block(dict[:body])
                )
    end
end

"""
    @functionfactory expr::Expr

Construct a `FunctionFactory` directly from a function definition.
"""
macro functionfactory(expr::Expr) expr=[expr];:(FunctionFactory($expr...)) end

"""
    (ff::FunctionFactory)()

Convert a `FunctionFactory` to the `Expr` representation of the `function` it describes.
"""
function (ff::FunctionFactory)()
    combinedef(Dict(
        :name           =>      ff.name,
        :args           =>      [arg() for arg in ff.args],
        :kwargs         =>      [kwarg() for kwarg in ff.kwargs],
        :rtype          =>      ff.rtype,
        :whereparams    =>      [param() for param in ff.params],
        :body           =>      ff.body()
    ))
end

"""
    addargs!(ff::FunctionFactory,args::FExpr...)
    addargs!(ff::FunctionFactory,args::Argument...)

Add a couple of positional arguments to a function factory.
"""
addargs!(ff::FunctionFactory,args::FExpr...)=addargs!(ff,Argument.(args)...)
addargs!(ff::FunctionFactory,args::Argument...)=push!(ff.args,args...)

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
addkwargs!(ff::FunctionFactory,kwargs::Argument...)=push!(ff.kwargs,kwargs...)

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
extendbody!(ff::FunctionFactory,parts::Block...)=push!(ff.body,parts...)

"""
    @extendbody! ff parts::FExpr...

Extend the body of a function factory.
"""
macro extendbody!(ff,parts::FExpr...) :(extendbody!($(esc(ff)),$parts...)) end

"""
    TypeFactory(name::Symbol,mutable::Bool,params::Vector{Parameter},supertype::FExpr,fields::Vector{Field},constructors::Vector{FunctionFactory})
    TypeFactory(    name::Symbol;
                    mutable::Bool=false,
                    params::Vector{Parameter}=Parameter[],
                    supertype::FExpr=:Any,
                    fields::Vector{Field}=Field[],
                    constructors::Vector{FunctionFactory}=FunctionFactory[]
                    )
    TypeFactory(expr::Expr)

The struct to describe a `struct`.
"""
struct TypeFactory <: AbstractFactory
    name::Symbol
    mutable::Bool
    params::Vector{Parameter}
    supertype::FExpr
    fields::Vector{Field}
    constructors::Vector{FunctionFactory}
    function TypeFactory(name::Symbol,mutable::Bool,params::Vector{Parameter},supertype::FExpr,fields::Vector{Field},constructors::Vector{FunctionFactory})
        new(name,mutable,params,supertype,fields,constructors)
    end
    function TypeFactory(   name::Symbol;
                            mutable::Bool=false,
                            params::Vector{Parameter}=Parameter[],
                            supertype::FExpr=:Any,
                            fields::Vector{Field}=Field[],
                            constructors::Vector{FunctionFactory}=FunctionFactory[]
                            )
        new(mutable,name,params,supertype,fields,constructors)
    end
    function TypeFactory(expr::Expr)
        dict=splitstructdef(expr)
        new(    dict[:name],
                dict[:mutable],
                Parameter.(dict[:params]),
                dict[:supertype],
                Field[Field(field...) for field in dict[:fields]],
                FunctionFactory.(dict[:constructors])
        )
    end
end

"""
    @typefactory expr::Expr

Construct a `TypeFactory` directly from a type definition.
"""
macro typefactory(expr::Expr) expr=[expr];:(TypeFactory($expr...)) end

"""
    (tf::TypeFactory)()

Convert a `TypeFactory` to the `Expr` representation of the `struct` it describes.
"""
function (tf::TypeFactory)()
    combinestructdef(Dict(
        :name           =>      tf.name,
        :mutable        =>      tf.mutable,
        :params         =>      [param() for param in tf.params],
        :supertype      =>      tf.supertype,
        :fields         =>      [(field.name,field.type) for field in tf.fields],
        :constructors   =>      [constructor() for constructor in tf.constructors]
        ))
end

"""
    addfields!(tf::TypeFactory,fields::FExpr...)
    addfields!(tf::TypeFactory,fields::Field...)

Add a couple of fields to a type factory.
"""
addfields!(tf::TypeFactory,fields::FExpr...)=addfields!(tf,Field.(fields)...)
addfields!(tf::TypeFactory,fields::Field...)=push!(tf.fields,fields...)

"""
    @addfields! tf fields::FExpr...

Add a couple of fields to a type factory.
"""
macro addfields!(tf,fields::FExpr...) :(addfields!($(esc(tf)),$fields...)) end

"""
    addconstructors!(tf::TypeFactory,constructors::Expr...)
    addconstructors!(tf::TypeFactory,constructors::FunctionFactory...)

Add a couple of constructors to a type factory.
"""
addconstructors!(tf::TypeFactory,constructors::Expr...)=addconstructors!(tf,FunctionFactory.(constructors)...)
addconstructors!(tf::TypeFactory,constructors::FunctionFactory...)=push!(tf.constructors,constructors...)

"""
    @addconstructors! tf constructors::Expr...

Add a couple of constructors to a type factory.
"""
macro addconstructors!(tf,constructors::Expr...) :(addconstructors!($(esc(tf)),$constructors...)) end

"""
    addparams!(f::Union{FunctionFactory,TypeFactory},params::FExpr...)
    addparams!(f::Union{FunctionFactory,TypeFactory},params::Parameter...)

Add a couple of method parameters to a function factory or a type factory.
"""
addparams!(f::Union{FunctionFactory,TypeFactory},params::FExpr...)=addparams!(f,Parameter.(params)...)
addparams!(f::Union{FunctionFactory,TypeFactory},params::Parameter...)=push!(f.params,params...)

"""
    @addparams! f params::FExpr...

Add a couple of method parameters to a function factory or a type factory.
"""
macro addparams!(ff,params::FExpr...) :(addparams!($(esc(ff)),$params...)) end

end #module
