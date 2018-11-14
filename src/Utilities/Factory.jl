module Factory

using MacroTools: block
using MacroTools: splitarg,combinefield,combinearg
using MacroTools: splitstructdef,combinestructdef
using MacroTools: splitdef,combinedef
using Printf: @printf,@sprintf
import MacroTools: rmlines

export FExpr
export symbols,escape
export Inference,@inference
export Argument,@argument
export Parameter,@parameter
export Field,@field
export Block,@block,@push!,rmlines!,@rmlines!,rmlines,@rmlines
export FunctionFactory,@functionfactory,addargs!,@addargs!,addkwargs!,@addkwargs!,extendbody!,@extendbody!
export TypeFactory,@typefactory,addfields!,@addfields!,addconstructors!,@addconstructors!
export addparams!,@addparams!
export paramnames

@generated function maxfieldnamelength(::T) where T
    result=1
    for name in T|>fieldnames
        result=max(result,name|>string|>length)
    end
    return :($result)
end

"Factory expression types, which is defined as `Union{Symbol,Expr}`."
const FExpr=Union{Symbol,Expr}

"""
    symbols(expr)
"""
function symbols(expr)
    result=[]
    if isa(expr,Symbol)
        push!(result,expr)
    elseif isa(expr,Expr)
        for arg in expr.args
            append!(result,symbols(arg))
        end
    end
    return result
end

"""
    escape(expr,::NTuple{N,Symbol}=()) where N
    escape(expr::Symbol,escaped::NTuple{N,Symbol}=()) where N
    escape(expr::Expr,escaped::NTuple{N,Symbol}=()) where N

Escape the symbols sepecified by `escaped` in the input expression.
"""
escape(expr,::NTuple{N,Symbol}=()) where N=expr
escape(expr::Symbol,escaped::NTuple{N,Symbol}=()) where N=expr ∈ escaped ? :($(esc(expr))) : expr
escape(expr::Expr,escaped::NTuple{N,Symbol}=()) where N=(expr.head==:escape ? expr : Expr(expr.head,(escape(arg,escaped) for arg in expr.args)...))

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
    SubTypeofAbstractFactory.(exprs::Vector{FExpr};unescaped::NTuple{N,Symbol}=()) where N
    SubTypeofAbstractFactory.(exprs::NTuple{N,FExpr};unescaped::NTuple{M,Symbol}=()) where {N,M}

Broadcast construction of concrete subtypes of `AbstractFactory`.
"""
Base.Broadcast.broadcasted(f::Type{<:AbstractFactory},exprs::Vector{FExpr};unescaped::NTuple{N,Symbol}=()) where N=[f(expr,unescaped=unescaped) for expr in exprs]
Base.Broadcast.broadcasted(f::Type{<:AbstractFactory},exprs::NTuple{N,FExpr};unescaped::NTuple{M,Symbol}=()) where {N,M}=NTuple{N,f}(f(expr,unescaped=unescaped) for expr in exprs)

"""
    Inference(head::Union{Symbol,Nothing},name::Union{Symbol,Nothing},params::Union{Inference,Vector{Inference},Nothing},escape::Bool)
    Inference(;
            head::Union{Symbol,Nothing}=nothing,
            name::Union{Symbol,Nothing}=nothing,
            params::Union{Inference,Vector{Inference},Nothing}=nothing,
            escape::Bool=true
            )
    Inference(expr::FExpr;unescaped::NTuple{N,Symbol}=()) where N

The struct to describe a type inference.
"""
mutable struct Inference <: AbstractFactory
    head::Union{Symbol,Nothing}
    name::Union{Symbol,Nothing}
    params::Union{Inference,Vector{Inference},Nothing}
    escape::Bool
    function Inference(head::Union{Symbol,Nothing},name::Union{Symbol,Nothing},params::Union{Inference,Vector{Inference},Nothing},escape::Bool)
        @assert (head===nothing)==(params===nothing) "Inference error: `head` and `params` must be or not be `nothing` at the same time."
        new(head,name,params,escape)
    end
end
Inference(;
        head::Union{Symbol,Nothing}=nothing,
        name::Union{Symbol,Nothing}=nothing,
        params::Union{Inference,Vector{Inference},Nothing}=nothing,
        escape::Bool=true
        )=Inference(head,name,params,escape)
function Inference(expr::FExpr;unescaped::NTuple{N,Symbol}=()) where N
    if isa(expr,Symbol)
        return Inference(nothing,expr,nothing,expr ∉ unescaped)
    elseif expr.head==:(<:) || expr.head==:(::)
        @assert length(expr.args)==1 "Inference error: wrong input expr."
        return Inference(expr.head,nothing,Inference(expr.args[1],unescaped=unescaped),false)
    elseif expr.head==:curly
        return Inference(expr.head,expr.args[1],Inference.(expr.args[2:end],unescaped=unescaped),expr.args[1] ∉ unescaped)
    else
        error("Inference error: wrong input expr.")
    end
end

"""
    @inference expr::FExpr unescaped::FExpr=:()

Construct an `Inference` directly from a type inference.
"""
macro inference(expr::FExpr,unescaped::FExpr=:()) expr=[expr];:(Inference($expr...;unescaped=$(esc(unescaped)))) end

"""
    (i::Inference)()

Convert a `Inference` to the `Expr` representation of the type inference it describes.
"""
function (i::Inference)()
    if i.head===nothing
        return i.escape ? :($(esc(i.name))) : i.name
    elseif i.head==:(::) || i.head==:(<:)
        return Expr(i.head,i.params())
    elseif i.head==:curly
        return Expr(i.head,i.escape ? :($(esc(i.name))) : i.name,[param() for param in i.params]...)
    end
end

"""
    Argument(name::Union{Symbol,Nothing},type::Inference,slurp::Bool,default::Any)
    Argument(;name::Union{Symbol,Nothing}=nothing,type::Inference=Inference(:Any),slurp::Bool=false,default::Any=nothing)
    Argument(expr::FExpr;unescaped::NTuple{N,Symbol}=()) where N)

The struct to describe a argument of a `function`.
"""
mutable struct Argument <: AbstractFactory
    name::Union{Symbol,Nothing}
    type::Inference
    slurp::Bool
    default::Any
end
Argument(;name::Union{Symbol,Nothing}=nothing,type::Inference=Inference(:Any),slurp::Bool=false,default::Any=nothing)=Argument(name,type,slurp,default)
function Argument(expr::FExpr;unescaped::NTuple{N,Symbol}=()) where N
    cache=splitarg(expr)
    Argument(cache[1],Inference(cache[2],unescaped=unescaped),cache[3],escape(cache[4],tuple(setdiff(symbols(cache[4]),unescaped)...)))
end

"""
    @argument expr::FExpr unescaped::FExpr=:()

Construct an `Argument` directly from an argument statement.
"""
macro argument(expr::FExpr,unescaped::FExpr=:()) expr=[expr];:(Argument($expr...;unescaped=$(esc(unescaped)))) end

"""
    (a::Argument)()

Convert an `Argument` to the `Expr` representation of the argument it describes.
"""
(a::Argument)()=combinearg(a.name,a.type(),a.slurp,a.default)

"""
    Parameter(name::Union{Symbol,Nothing},type::Union{Inference,Nothing})
    Parameter(;name::Union{Symbol,Nothing}=nothing,type::Union{Inference,Nothing}=nothing)
    Parameter(expr::FExpr;unescaped::NTuple{N,Symbol}=()) where N

The struct to describe a parameter of a `function` or a `type`.
"""
mutable struct Parameter <: AbstractFactory
    name::Union{Symbol,Nothing}
    type::Union{Inference,Nothing}
end
Parameter(;name::Union{Symbol,Nothing}=nothing,type::Union{Inference,Nothing}=nothing)=Parameter(name,type)
function Parameter(expr::FExpr;unescaped::NTuple{N,Symbol}=()) where N
    if isa(expr,Symbol)
        name,type=expr,nothing
    elseif expr.head==:(<:)
        (name,type)=length(expr.args)==2 ? expr.args : (nothing,expr.args[1])
        type=Inference(type,unescaped=unescaped)
    elseif expr.head==:(::)
        @assert length(expr.args)==1 "Parameter error: wrong input expr."
        name,type=nothing,Inference(expr.args[1],unescaped=unescaped)
    else
        error("Parameter error: wrong input expr.")
    end
    Parameter(name,type)
end

"""
    @parameter expr::FExpr unescaped::FExpr=:()

Construct a `Parameter` directly from an parameter statement.
"""
macro parameter(expr::FExpr,unescaped::FExpr=:()) expr=[expr];:(Parameter($expr...;unescaped=$(esc(unescaped)))) end

"""
    (p::Parameter)()

Convert a `Parameter` to the `Expr` representation of the parameter it describes.
"""
(p::Parameter)()=p.name===nothing ? Expr(:(<:),p.type()) : p.type===nothing ? p.name : Expr(:(<:),p.name,p.type())

"""
    Field(name::Symbol,type::Inference)
    Field(name::Symbol;type::FExpr=Inference(:Any))
    Field(expr::Expr;unescaped::NTuple{N,Symbol}=()) where N

The struct to describe a field of a `struct`.
"""
mutable struct Field <: AbstractFactory
    name::Symbol
    type::Inference
end
Field(name::Symbol,type::FExpr=Inference(:Any))=Field(name,type)
Field(expr::FExpr;unescaped::NTuple{N,Symbol}=()) where N=(cache=splitarg(expr);Field(cache[1],Inference(cache[2],unescaped=unescaped)))

"""
    @field expr::FExpr unescaped::FExpr=:()

Construct a `Field` directly from a field statement.
"""
macro field(expr::FExpr,unescaped::FExpr=:()) expr=[expr];:(Field($expr...;unescaped=$(esc(unescaped)))) end

"""
    (f::Field)()

Convert a `Field` to the `Expr` representation of the field it describes.
"""
(f::Field)()=combinefield((f.name,f.type()))

"""
    convert(::Type{Tuple},f::Field)
    convert(::Type{Tuple{Symbol,FExpr}},f::Field)

Convert a `Field` to tuple.
"""
Base.convert(::Type{Tuple},f::Field)=(f.name,f.type())
Base.convert(::Type{Tuple{Symbol,FExpr}},f::Field)=convert(Tuple,f)

"""
    Block(parts::FExpr...)

The struct to describe a `begin ... end` block.
"""
mutable struct Block <: AbstractFactory
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
Base.push!(::Block)=nothing
Base.push!(b::Block,parts::FExpr...)=push!(b,Block(parts...))
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
    FunctionFactory(name::Symbol,args::Vector{Argument},kwargs::Vector{Argument},rtype::Inference,params::Vector{Parameter},body::Block)
    FunctionFactory(    name::Symbol;
                        args::Vector{Argument}=Argument[],
                        kwargs::Vector{Argument}=Argument[],
                        rtype::Inference=Inference(:Any),
                        params::Vector{Parameter}=Parameter[],
                        body::Block=Block()
                        )
    FunctionFactory(expr::Expr;unescaped::NTuple{N,Symbol}=()) where N

The struct to describe a `function`.
"""
mutable struct FunctionFactory <: AbstractFactory
    name::Symbol
    args::Vector{Argument}
    kwargs::Vector{Argument}
    rtype::Inference
    params::Vector{Parameter}
    body::Block
end
function FunctionFactory(   name::Symbol;
                            args::Vector{Argument}=Argument[],
                            kwargs::Vector{Argument}=Argument[],
                            rtype::Inference=Inference(:Any),
                            params::Vector{Parameter}=Parameter[],
                            body::Block=Block(),
                            )
    FunctionFactory(name,args,kwargs,rtype,params,body)
end
function FunctionFactory(expr::Expr;unescaped::NTuple{N,Symbol}=()) where N
    dict=splitdef(expr)
    FunctionFactory(    dict[:name],
                        Argument.(dict[:args],unescaped=unescaped),
                        Argument.(dict[:kwargs],unescaped=unescaped),
                        Inference(get(dict,:rtype,:Any),unescaped=unescaped),
                        Parameter.(collect(FExpr,dict[:whereparams]),unescaped=unescaped),
                        Block(dict[:body]),
                        )
end

"""
    @functionfactory expr::FExpr unescaped::FExpr=:()

Construct a `FunctionFactory` directly from a function definition.
"""
macro functionfactory(expr::Expr,unescaped::FExpr=:()) expr=[expr];:(FunctionFactory($expr...;unescaped=$(esc(unescaped)))) end

"""
    (ff::FunctionFactory)()

Convert a `FunctionFactory` to the `Expr` representation of the `function` it describes.
"""
function (ff::FunctionFactory)()
    combinedef(Dict(
        :name           =>      ff.name,
        :args           =>      [arg() for arg in ff.args],
        :kwargs         =>      [kwarg() for kwarg in ff.kwargs],
        :rtype          =>      ff.rtype(),
        :whereparams    =>      [param() for param in ff.params],
        :body           =>      ff.body()
    ))
end

"""
    addargs!(ff::FunctionFactory,args::Argument...)
    addargs!(ff::FunctionFactory,args::FExpr...;unescaped::NTuple{N,Symbol}=()) where N

Add a couple of positional arguments to a function factory.
"""
addargs!(::FunctionFactory)=nothing
addargs!(ff::FunctionFactory,args::Argument...)=push!(ff.args,args...)
addargs!(ff::FunctionFactory,args::FExpr...;unescaped::NTuple{N,Symbol}=()) where N=addargs!(ff,Argument.(args,unescaped=unescaped)...)

"""
    @addargs! ff args::FExpr...

Add a couple of positional arguments to a function factory.
"""
macro addargs!(ff,args::FExpr...) :(addargs!($(esc(ff)),$args...)) end

"""
    addkwargs!(ff::FunctionFactory,kwargs::Argument...)
    addkwargs!(ff::FunctionFactory,kwargs::FExpr...;unescaped::NTuple{N,Symbol}=()) where N

Add a couple of keyword arguments to a function factory.
"""
addkwargs!(::FunctionFactory)=nothing
addkwargs!(ff::FunctionFactory,kwargs::Argument...)=push!(ff.kwargs,kwargs...)
addkwargs!(ff::FunctionFactory,kwargs::FExpr...;unescaped::NTuple{N,Symbol}=()) where N=addkwargs!(ff,Argument.(kwargs,unescaped=unescaped)...)

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
extendbody!(::FunctionFactory)=nothing
extendbody!(ff::FunctionFactory,parts::FExpr...)=extendbody!(ff,Block(parts...))
extendbody!(ff::FunctionFactory,parts::Block...)=push!(ff.body,parts...)

"""
    @extendbody! ff parts::FExpr...

Extend the body of a function factory.
"""
macro extendbody!(ff,parts::FExpr...) :(extendbody!($(esc(ff)),$parts...)) end

"""
    TypeFactory(name::Symbol,mutable::Bool,params::Vector{Parameter},supertype::Inference,fields::Vector{Field},constructors::Vector{FunctionFactory})
    TypeFactory(    name::Symbol;
                    mutable::Bool=false,
                    params::Vector{Parameter}=Parameter[],
                    supertype::Inference=Inference(:Any),
                    fields::Vector{Field}=Field[],
                    constructors::Vector{FunctionFactory}=FunctionFactory[],
                    )
    TypeFactory(expr::Expr;unescaped::NTuple{N,Symbol}=()) where N

The struct to describe a `struct`.
"""
mutable struct TypeFactory <: AbstractFactory
    name::Symbol
    mutable::Bool
    params::Vector{Parameter}
    supertype::Inference
    fields::Vector{Field}
    constructors::Vector{FunctionFactory}
end
function TypeFactory(   name::Symbol;
                        mutable::Bool=false,
                        params::Vector{Parameter}=Parameter[],
                        supertype::Inference=Inference(:Any),
                        fields::Vector{Field}=Field[],
                        constructors::Vector{FunctionFactory}=FunctionFactory[],
                        )
    TypeFactory(name,mutable,params,supertype,fields,constructors)
end
function TypeFactory(expr::Expr;unescaped::NTuple{N,Symbol}=()) where N
    dict=splitstructdef(expr)
    TypeFactory(    dict[:name],
                    dict[:mutable],
                    Parameter.(dict[:params],unescaped=unescaped),
                    Inference(dict[:supertype],unescaped=unescaped),
                    [Field(field[1],Inference(field[2],unescaped=unescaped)) for field in dict[:fields]],
                    FunctionFactory.(dict[:constructors],unescaped=unescaped)
    )
end

"""
    @typefactory expr::Expr unescaped::FExpr=:()

Construct a `TypeFactory` directly from a type definition.
"""
macro typefactory(expr::Expr,unescaped::FExpr=:()) expr=[expr];:(TypeFactory($expr...;unescaped=$(esc(unescaped)))) end

"""
    (tf::TypeFactory)()

Convert a `TypeFactory` to the `Expr` representation of the `struct` it describes.
"""
function (tf::TypeFactory)()
    combinestructdef(Dict(
        :name           =>      tf.name,
        :mutable        =>      tf.mutable,
        :params         =>      [param() for param in tf.params],
        :supertype      =>      tf.supertype(),
        :fields         =>      [convert(Tuple,field) for field in tf.fields],
        :constructors   =>      [constructor() for constructor in tf.constructors]
        ))
end

"""
    addfields!(tf::TypeFactory,fields::Field...)
    addfields!(tf::TypeFactory,fields::FExpr...;unescaped::NTuple{N,Symbol}=()) where N

Add a couple of fields to a type factory.
"""
addfields!(::TypeFactory)=nothing
addfields!(tf::TypeFactory,fields::Field...)=push!(tf.fields,fields...)
addfields!(tf::TypeFactory,fields::FExpr...;unescaped::NTuple{N,Symbol}=()) where N=addfields!(tf,Field.(fields,unescaped=unescaped)...)

"""
    @addfields! tf fields::FExpr...

Add a couple of fields to a type factory.
"""
macro addfields!(tf,fields::FExpr...) :(addfields!($(esc(tf)),$fields...)) end

"""
    addconstructors!(tf::TypeFactory,constructors::FunctionFactory...)
    addconstructors!(tf::TypeFactory,constructors::Expr...;unescaped::NTuple{N,Symbol}=()) where N

Add a couple of constructors to a type factory.
"""
addconstructors!(::TypeFactory)=nothing
addconstructors!(tf::TypeFactory,constructors::FunctionFactory...)=push!(tf.constructors,constructors...)
addconstructors!(tf::TypeFactory,constructors::Expr...;unescaped::NTuple{N,Symbol}=()) where N=addconstructors!(tf,FunctionFactory.(constructors,unescaped=unescaped)...)

"""
    @addconstructors! tf constructors::Expr...

Add a couple of constructors to a type factory.
"""
macro addconstructors!(tf,constructors::Expr...) :(addconstructors!($(esc(tf)),$constructors...)) end

"""
    addparams!(f::Union{FunctionFactory,TypeFactory},params::Parameter...)
    addparams!(f::Union{FunctionFactory,TypeFactory},params::FExpr...;unescaped::NTuple{N,Symbol}=()) where N

Add a couple of method parameters to a function factory or a type factory.
"""
addparams!(::Union{FunctionFactory,TypeFactory})=nothing
addparams!(f::Union{FunctionFactory,TypeFactory},params::Parameter...)=push!(f.params,params...)
addparams!(f::Union{FunctionFactory,TypeFactory},params::FExpr...;unescaped::NTuple{N,Symbol}=()) where N=addparams!(f,Parameter.(params,unescaped=unescaped)...)

"""
    @addparams! f params::FExpr...

Add a couple of method parameters to a function factory or a type factory.
"""
macro addparams!(ff,params::FExpr...) :(addparams!($(esc(ff)),$params...)) end

"""
    paramnames(expr::Union{Expr,FunctionFactory,TypeFactory}) -> NTuple{N,Symbol} where N

Get the type/method parameters names.
"""
paramnames(f::Union{FunctionFactory,TypeFactory})=NTuple{f.params|>length,Symbol}(param.name for param in f.params)
paramnames(expr::Expr)=paramnames(expr.head==:struct ? TypeFactory(expr) : FunctionFactory(expr))

end #module
