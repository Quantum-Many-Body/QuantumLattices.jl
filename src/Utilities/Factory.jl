module Factory

using MacroTools: block
using MacroTools: splitarg,combinefield,combinearg
using MacroTools: splitstructdef,combinestructdef
using MacroTools: splitdef,combinedef
using Printf: @printf,@sprintf
import MacroTools: rmlines

export FExpr,RawExpr
export escape
export Inference,@inference
export Argument,@argument
export Parameter,@parameter
export Field,@field
export Block,@block,@push!,rmlines!,@rmlines!,rmlines,@rmlines
export FunctionFactory,@functionfactory,addargs!,@addargs!,addkwargs!,@addkwargs!,addwhereparams!,@addwhereparams!,extendbody!,@extendbody!
export TypeFactory,@typefactory,addfields!,@addfields!,addconstructors!,@addconstructors!
export addparams!,@addparams!

@generated function maxfieldnamelength(::T) where T
    result=1
    for name in T|>fieldnames
        result=max(result,name|>string|>length)
    end
    return :($result)
end

"Factory expression types, which is defined as `Union{Symbol,Expr}`."
const FExpr=Union{Symbol,Expr}

"Whether or not to show raw expressions of factoies."
const RawExpr=Val(true)

"""
    names(expr::Union{Symbol,Expr}) -> Vector{Symbol}

Get all the variable names in an `Expr`.
"""
Base.names(expr)=[]
Base.names(expr::Symbol)=[expr]
function Base.names(expr::Expr)
    result=[]
    for arg in expr.args
        append!(result,names(arg))
    end
    result
end

"""
    escape(expr,::NTuple{N,Symbol}=()) where N -> Any
    escape(expr::Symbol,escaped::NTuple{N,Symbol}=()) where N -> FExpr
    escape(expr::Expr,escaped::NTuple{N,Symbol}=()) where N -> Expr

Escape the variables sepecified by `escaped` in the input expression.
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
    ==(f1::F,f2::F) where F<:AbstractFactory -> Bool

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
            vstr=string(value(RawExpr))
        elseif isa(value,Vector{<:AbstractFactory})
            eltypename=value|>eltype|>nameof
            strings=Tuple(string(af(RawExpr)) for af in value)
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
    replace(f::AbstractFactory;kwargs...) -> typeof(f)

Return a copy of a concrete `AbstractFactory` with some of the field values replaced by the keyword arguments.
"""
Base.replace(f::AbstractFactory;kwargs...)=(f|>typeof)((get(kwargs,key,getfield(f,name)) for name in f|>typeof|>fieldnames)...)

"""
    Inference(head::Union{Symbol,Nothing},name::Union{Symbol,Nothing},params::Union{Inference,Vector{Inference},Nothing})
    Inference(;
            head::Union{Symbol,Nothing}=nothing,
            name::Union{Symbol,Nothing}=nothing,
            params::Union{Inference,Vector{Inference},Nothing}=nothing,
            )
    Inference(expr::FExpr)

The struct to describe a type inference.
"""
mutable struct Inference <: AbstractFactory
    head::Union{Symbol,Nothing}
    name::Union{Symbol,Nothing}
    params::Union{Inference,Vector{Inference},Nothing}
    function Inference(head::Union{Symbol,Nothing},name::Union{Symbol,Nothing},params::Union{Inference,Vector{Inference},Nothing})
        @assert (head===nothing)==(params===nothing) "Inference error: `head` and `params` must be or not be `nothing` at the same time."
        new(head,name,params)
    end
end
Inference(;
        head::Union{Symbol,Nothing}=nothing,
        name::Union{Symbol,Nothing}=nothing,
        params::Union{Inference,Vector{Inference},Nothing}=nothing
        )=Inference(head,name,params)
function Inference(expr::FExpr)
    if isa(expr,Symbol)
        return Inference(nothing,expr,nothing)
    elseif expr.head==:(<:) || expr.head==:(::)
        @assert length(expr.args)==1 "Inference error: wrong input expr."
        return Inference(expr.head,nothing,Inference(expr.args[1]))
    elseif expr.head==:curly
        return Inference(expr.head,expr.args[1],Inference.(expr.args[2:end]))
    else
        error("Inference error: wrong input expr.")
    end
end

"""
    @inference expr::FExpr

Construct an `Inference` directly from a type inference.
"""
macro inference(expr::FExpr) expr=[expr];:(Inference($expr...)) end

"""
    (i::Inference)(::typeof(RawExpr)) -> FExpr
    (i::Inference)(;unescaped::NTuple{N,Symbol}=()) where N -> FExpr

Convert a `Inference` to the `Expr` representation of the type inference it describes.
"""
function (i::Inference)(::typeof(RawExpr))
    if i.head===nothing
        return i.name
    elseif i.head==:(::) || i.head==:(<:)
        return Expr(i.head,i.params(RawExpr))
    else
        return Expr(i.head,i.name,[param(RawExpr) for param in i.params]...)
    end
end
function (i::Inference)(;unescaped::NTuple{N,Symbol}=()) where N
    if i.head===nothing
        return i.name ∉ unescaped ? :($(esc(i.name))) : i.name
    elseif i.head==:(::) || i.head==:(<:)
        return Expr(i.head,i.params(unescaped=unescaped))
    else
        return Expr(i.head,i.name ∉ unescaped ? :($(esc(i.name))) : i.name,[param(unescaped=unescaped) for param in i.params]...)
    end
end

"""
    Argument(name::Union{Symbol,Nothing},type::Inference,slurp::Bool,default::Any)
    Argument(;name::Union{Symbol,Nothing}=nothing,type::Inference=Inference(:Any),slurp::Bool=false,default::Any=nothing)
    Argument(expr::FExpr)

The struct to describe a argument of a `function`.
"""
mutable struct Argument <: AbstractFactory
    name::Union{Symbol,Nothing}
    type::Inference
    slurp::Bool
    default::Any
end
Argument(;name::Union{Symbol,Nothing}=nothing,type::Inference=Inference(:Any),slurp::Bool=false,default::Any=nothing)=Argument(name,type,slurp,default)
Argument(expr::FExpr)=(cache=splitarg(expr);Argument(cache[1],Inference(cache[2]),cache[3],cache[4]))

"""
    @argument expr::FExpr

Construct an `Argument` directly from an argument statement.
"""
macro argument(expr::FExpr) expr=[expr];:(Argument($expr...)) end

"""
    (a::Argument)(::typeof(RawExpr)) -> Expr
    (a::Argument)(;unescaped::NTuple{N,Symbol}=()) where N -> Expr

Convert an `Argument` to the `Expr` representation of the argument it describes.
"""
(a::Argument)(::typeof(RawExpr))=combinearg(a.name,a.type(RawExpr),a.slurp,a.default)
(a::Argument)(;unescaped::NTuple{N,Symbol}=()) where N=combinearg(a.name,a.type(unescaped=unescaped),a.slurp,escape(a.default,tuple(setdiff(names(a.default),unescaped)...)))

"""
    Parameter(name::Union{Symbol,Nothing},type::Union{Inference,Nothing})
    Parameter(;name::Union{Symbol,Nothing}=nothing,type::Union{Inference,Nothing}=nothing)
    Parameter(expr::FExpr)

The struct to describe a parameter of a `function` or a `type`.
"""
mutable struct Parameter <: AbstractFactory
    name::Union{Symbol,Nothing}
    type::Union{Inference,Nothing}
end
Parameter(;name::Union{Symbol,Nothing}=nothing,type::Union{Inference,Nothing}=nothing)=Parameter(name,type)
function Parameter(expr::FExpr)
    if isa(expr,Symbol)
        name,type=expr,nothing
    elseif expr.head==:(<:)
        (name,type)=length(expr.args)==2 ? expr.args : (nothing,expr.args[1])
        type=Inference(type)
    elseif expr.head==:(::)
        @assert length(expr.args)==1 "Parameter error: wrong input expr."
        name,type=nothing,Inference(expr.args[1])
    else
        error("Parameter error: wrong input expr.")
    end
    Parameter(name,type)
end

"""
    @parameter expr::FExpr

Construct a `Parameter` directly from an parameter statement.
"""
macro parameter(expr::FExpr) expr=[expr];:(Parameter($expr...)) end

"""
    (p::Parameter)(::typeof(RawExpr)) -> FExpr
    (p::Parameter)(;unescaped::NTuple{N,Symbol}=()) where N -> FExpr

Convert a `Parameter` to the `Expr` representation of the parameter it describes.
"""
(p::Parameter)(::typeof(RawExpr))=p.name===nothing ? Expr(:(<:),p.type(RawExpr)) : p.type===nothing ? p.name : Expr(:(<:),p.name,p.type(RawExpr))
function (p::Parameter)(;unescaped::NTuple{N,Symbol}=()) where N
    if p.name===nothing
        Expr(:(<:),p.type(unescaped=unescaped))
    elseif p.type===nothing
        p.name
    else
        Expr(:(<:),p.name,p.type(unescaped=unescaped))
    end
end

"""
    Field(name::Symbol,type::Inference)
    Field(;name::Symbol,type::FExpr=Inference(:Any))
    Field(expr::FExpr)

The struct to describe a field of a `struct`.
"""
mutable struct Field <: AbstractFactory
    name::Symbol
    type::Inference
end
Field(;name::Symbol,type::FExpr=Inference(:Any))=Field(name,type)
Field(expr::FExpr)=(cache=splitarg(expr);Field(cache[1],Inference(cache[2])))

"""
    @field expr::FExpr

Construct a `Field` directly from a field statement.
"""
macro field(expr::FExpr) expr=[expr];:(Field($expr...)) end

"""
    (f::Field)(::typeof(RawExpr)) -> Expr
    (f::Field)(;unescaped::NTuple{N,Symbol}=()) where N -> Expr

Convert a `Field` to the `Expr` representation of the field it describes.
"""
(f::Field)(::typeof(RawExpr))=combinefield((f.name,f.type(RawExpr)))
(f::Field)(;unescaped::NTuple{N,Symbol}=()) where N=combinefield((f.name,f.type(unescaped=unescaped)))

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
    (b::Block)(::typeof(RawExpr)) -> Expr
    (b::Block)(;escaped::NTuple{N,Symbol}=()) where N -> Expr

Convert a `Block` to the `Expr` representation of the `begin ... end` block it describes.
"""
(b::Block)(::typeof(RawExpr))=Expr(:block,b.body...)
(b::Block)(;escaped::NTuple{N,Symbol}=()) where N=Expr(:block,(escape(part,escaped) for part in b.body)...)

"""
    push!(b::Block,parts::FExpr...) -> Block
    push!(b::Block,parts::Block...) -> Block

Push other parts into the body of a block.
"""
Base.push!(b::Block)=b
Base.push!(b::Block,parts::FExpr...)=push!(b,Block(parts...))
function Base.push!(b::Block,parts::Block...)
    for part in parts
        append!(b.body,part.body)
    end
    return b
end

"""
    @push! b parts::FExpr...

Push other parts into the body of a block.
"""
macro push!(b,parts::FExpr...) :(push!($(esc(b)),$parts...)) end

"""
    rmlines!(b::Block) -> Block

Remove line number nodes in the body of a block.
"""
rmlines!(b::Block)=(filter!(part->!isa(part,LineNumberNode),b.body);b)

"""
    @rmlines! b::Expr

Remove line number nodes in the body of a block.
"""
macro rmlines!(b::Expr) :(rmlines!($(esc(b)))) end

"""
    rmlines(b::Block) -> Block

Return a copy of a block with the line number nodes removed.
"""
rmlines(b::Block)=Block(filter(part->!isa(part,LineNumberNode),b.body)...)

"""
    @rmlines b::Expr

Return a copy of a block with the line number nodes removed.
"""
macro rmlines(b::Expr) :(rmlines($(esc(b)))) end

"""
    FunctionFactory(name::FExpr,params::Vector{Inference},args::Vector{Argument},kwargs::Vector{Argument},rtype::Inference,whereparams::Vector{Parameter},body::Block)
    FunctionFactory(    ;name::FExpr,
                        params::Vector{Inference}=Inference[],
                        args::Vector{Argument}=Argument[],
                        kwargs::Vector{Argument}=Argument[],
                        rtype::Inference=Inference(:Any),
                        whereparams::Vector{Parameter}=Parameter[],
                        body::Block=Block()
                        )
    FunctionFactory(expr::Expr)

The struct to describe a `function`.
"""
mutable struct FunctionFactory <: AbstractFactory
    name::FExpr
    params::Vector{Inference}
    args::Vector{Argument}
    kwargs::Vector{Argument}
    rtype::Inference
    whereparams::Vector{Parameter}
    body::Block
end
function FunctionFactory(   ;name::FExpr,
                            params::Vector{Inference}=Inference[],
                            args::Vector{Argument}=Argument[],
                            kwargs::Vector{Argument}=Argument[],
                            rtype::Inference=Inference(:Any),
                            whereparams::Vector{Parameter}=Parameter[],
                            body::Block=Block(),
                            )
    FunctionFactory(name,params,args,kwargs,rtype,whereparams,body)
end
function FunctionFactory(expr::Expr)
    dict=splitdef(expr)
    FunctionFactory(    dict[:name],
                        Inference.(get(dict,:params,[])),
                        Argument.(dict[:args]),
                        Argument.(dict[:kwargs]),
                        Inference(get(dict,:rtype,:Any)),
                        Parameter.(collect(FExpr,dict[:whereparams])),
                        Block(dict[:body]),
                        )
end

"""
    @functionfactory expr::FExpr

Construct a `FunctionFactory` directly from a function definition.
"""
macro functionfactory(expr::Expr) expr=[expr];:(FunctionFactory($expr...)) end

"""
    (ff::FunctionFactory)(::typeof(RawExpr)) -> Expr
    (ff::FunctionFactory)(;unescaped::NTuple{N,Symbol}=(),escaped::NTuple{M,Symbol}=()) where {N,M} -> Expr

Convert a `FunctionFactory` to the `Expr` representation of the `function` it describes.
"""
function (ff::FunctionFactory)(::typeof(RawExpr))
    combinedef(Dict(
        :name           =>      ff.name,
        :params         =>      [param(RawExpr) for param in ff.params],
        :args           =>      [arg(RawExpr) for arg in ff.args],
        :kwargs         =>      [kwarg(RawExpr) for kwarg in ff.kwargs],
        :rtype          =>      ff.rtype(RawExpr),
        :whereparams    =>      [whereparam(RawExpr) for whereparam in ff.whereparams],
        :body           =>      ff.body(RawExpr)
    ))
end
function (ff::FunctionFactory)(;unescaped::NTuple{N,Symbol}=(),escaped::NTuple{M,Symbol}=()) where {N,M}
    combinedef(Dict(
        :name           =>      ff.name ∈ escaped ? :($(esc(ff.name))) : ff.name,
        :params         =>      [param(unescaped=unescaped) for param in ff.params],
        :args           =>      [arg(unescaped=unescaped) for arg in ff.args],
        :kwargs         =>      [kwarg(unescaped=unescaped) for kwarg in ff.kwargs],
        :rtype          =>      ff.rtype(unescaped=unescaped),
        :whereparams    =>      [whereparam(unescaped=unescaped) for whereparam in ff.whereparams],
        :body           =>      ff.body(escaped=escaped)
    ))
end

"""
    addargs!(ff::FunctionFactory,args::Argument...) -> FunctionFactory
    addargs!(ff::FunctionFactory,args::FExpr...) -> FunctionFactory

Add a couple of positional arguments to a function factory.
"""
addargs!(ff::FunctionFactory)=ff
addargs!(ff::FunctionFactory,args::Argument...)=(push!(ff.args,args...);ff)
addargs!(ff::FunctionFactory,args::FExpr...)=addargs!(ff,Argument.(args)...)

"""
    @addargs! ff args::FExpr...

Add a couple of positional arguments to a function factory.
"""
macro addargs!(ff,args::FExpr...) :(addargs!($(esc(ff)),$args...)) end

"""
    addkwargs!(ff::FunctionFactory,kwargs::Argument...) -> FunctionFactory
    addkwargs!(ff::FunctionFactory,kwargs::FExpr...) -> FunctionFactory

Add a couple of keyword arguments to a function factory.
"""
addkwargs!(ff::FunctionFactory)=ff
addkwargs!(ff::FunctionFactory,kwargs::Argument...)=(push!(ff.kwargs,kwargs...);ff)
addkwargs!(ff::FunctionFactory,kwargs::FExpr...)=addkwargs!(ff,Argument.(kwargs)...)

"""
    @addkwargs! ff kwargs::FExpr...

Add a couple of keyword arguments to a function factory.
"""
macro addkwargs!(ff,kwargs::FExpr...) :(addkwargs!($(esc(ff)),$kwargs...)) end

"""
    extendbody!(ff::FunctionFactory,parts::FExpr...) -> FunctionFactory
    extendbody!(ff::FunctionFactory,parts::Block...) -> FunctionFactory

Extend the body of a function factory.
"""
extendbody!(ff::FunctionFactory)=ff
extendbody!(ff::FunctionFactory,parts::FExpr...)=extendbody!(ff,Block(parts...))
extendbody!(ff::FunctionFactory,parts::Block...)=(push!(ff.body,parts...);ff)

"""
    @extendbody! ff parts::FExpr...

Extend the body of a function factory.
"""
macro extendbody!(ff,parts::FExpr...) :(extendbody!($(esc(ff)),$parts...)) end

"""
    addwhereparams!(ff::FunctionFactory,whereparams::Parameter...) -> FunctionFactory
    addwhereparams!(ff::FunctionFactory,whereparams::FExpr...) -> FunctionFactory

Add a couple of method where parameters to a function factory or a type factory.
"""
addwhereparams!(ff::FunctionFactory)=ff
addwhereparams!(ff::FunctionFactory,whereparams::Parameter...)=(push!(ff.whereparams,whereparams...);ff)
addwhereparams!(ff::FunctionFactory,whereparams::FExpr...)=addwhereparams!(ff,Parameter.(whereparams)...)

"""
    @addwhereparams! ff whereparams::FExpr...

Add a couple of method parameters to a function factory or a type factory.
"""
macro addwhereparams!(ff,whereparams::FExpr...) :(addwhereparams!($(esc(ff)),$whereparams...)) end

"""
    TypeFactory(name::Symbol,mutable::Bool,params::Vector{Parameter},supertype::Inference,fields::Vector{Field},constructors::Vector{FunctionFactory})
    TypeFactory(    ;name::Symbol,
                    mutable::Bool=false,
                    params::Vector{Parameter}=Parameter[],
                    supertype::Inference=Inference(:Any),
                    fields::Vector{Field}=Field[],
                    constructors::Vector{FunctionFactory}=FunctionFactory[],
                    )
    TypeFactory(expr::Expr)

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
function TypeFactory(   ;name::Symbol,
                        mutable::Bool=false,
                        params::Vector{Parameter}=Parameter[],
                        supertype::Inference=Inference(:Any),
                        fields::Vector{Field}=Field[],
                        constructors::Vector{FunctionFactory}=FunctionFactory[],
                        )
    TypeFactory(name,mutable,params,supertype,fields,constructors)
end
function TypeFactory(expr::Expr)
    dict=splitstructdef(expr)
    TypeFactory(    dict[:name],
                    dict[:mutable],
                    Parameter.(dict[:params]),
                    Inference(dict[:supertype]),
                    [Field(field[1],Inference(field[2])) for field in dict[:fields]],
                    FunctionFactory.(dict[:constructors])
    )
end

"""
    @typefactory expr::Expr

Construct a `TypeFactory` directly from a type definition.
"""
macro typefactory(expr::Expr) expr=[expr];:(TypeFactory($expr...)) end

"""
    (tf::TypeFactory)(::typeof(RawExpr)) -> Expr
    (tf::TypeFactory)(;unescaped::NTuple{N,Symbol}=()) where N -> Expr

Convert a `TypeFactory` to the `Expr` representation of the `struct` it describes.
"""
function (tf::TypeFactory)(::typeof(RawExpr))
    combinestructdef(Dict(
        :name           =>      tf.name,
        :mutable        =>      tf.mutable,
        :params         =>      [param(RawExpr) for param in tf.params],
        :supertype      =>      tf.supertype(RawExpr),
        :fields         =>      [(field.name,field.type(RawExpr)) for field in tf.fields],
        :constructors   =>      [constructor(RawExpr) for constructor in tf.constructors]
        ))
end
function (tf::TypeFactory)(;unescaped::NTuple{N,Symbol}=(),escaped::NTuple{M,Symbol}=()) where {N,M}
    combinestructdef(Dict(
        :name           =>      tf.name ∈ escaped ? :($(esc(tf.name))) : tf.name,
        :mutable        =>      tf.mutable,
        :params         =>      [param(unescaped=unescaped) for param in tf.params],
        :supertype      =>      tf.supertype(unescaped=unescaped),
        :fields         =>      [(field.name,field.type(unescaped=unescaped)) for field in tf.fields],
        :constructors   =>      [constructor(unescaped=unescaped,escaped=escaped) for constructor in tf.constructors]
        ))
end

"""
    addfields!(tf::TypeFactory,fields::Field...) -> TypeFactory
    addfields!(tf::TypeFactory,fields::FExpr...) -> TypeFactory

Add a couple of fields to a type factory.
"""
addfields!(tf::TypeFactory)=tf
addfields!(tf::TypeFactory,fields::Field...)=(push!(tf.fields,fields...);tf)
addfields!(tf::TypeFactory,fields::FExpr...)=addfields!(tf,Field.(fields)...)

"""
    @addfields! tf fields::FExpr...

Add a couple of fields to a type factory.
"""
macro addfields!(tf,fields::FExpr...) :(addfields!($(esc(tf)),$fields...)) end

"""
    addconstructors!(tf::TypeFactory,constructors::FunctionFactory...) -> TypeFactory
    addconstructors!(tf::TypeFactory,constructors::Expr...) -> TypeFactory

Add a couple of constructors to a type factory.
"""
addconstructors!(tf::TypeFactory)=tf
addconstructors!(tf::TypeFactory,constructors::FunctionFactory...)=(push!(tf.constructors,constructors...);tf)
addconstructors!(tf::TypeFactory,constructors::Expr...)=addconstructors!(tf,FunctionFactory.(constructors)...)

"""
    @addconstructors! tf constructors::Expr...

Add a couple of constructors to a type factory.
"""
macro addconstructors!(tf,constructors::Expr...) :(addconstructors!($(esc(tf)),$constructors...)) end

"""
    addparams!(f::FunctionFactory,params::Inference...) ->FunctionFactory
    addparams!(f::FunctionFactory,params::FExpr...) -> FunctionFactory
    addparams!(f::TypeFactory,params::Parameter...) -> TypeFactory
    addparams!(f::TypeFactory,params::FExpr...) -> TypeFactory


Add a couple of parameters to a function factory or a type factory.
"""
addparams!(f::Union{FunctionFactory,TypeFactory})=f
addparams!(f::FunctionFactory,params::Inference...)=(push!(f.params,params...);f)
addparams!(f::FunctionFactory,params::FExpr...)=addparams!(f,Inference.(params)...)
addparams!(f::TypeFactory,params::Parameter...)=(push!(f.params,params...);f)
addparams!(f::TypeFactory,params::FExpr...)=addparams!(f,Parameter.(params)...)

"""
    @addparams! f params::FExpr...

Add a couple of method parameters to a function factory or a type factory.
"""
macro addparams!(f,params::FExpr...) :(addparams!($(esc(f)),$params...)) end

end #module
