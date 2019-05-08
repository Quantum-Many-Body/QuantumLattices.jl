module Factories

using MacroTools: block
using MacroTools: splitarg,combinefield,combinearg
using MacroTools: splitstructdef,combinestructdef
using MacroTools: splitdef,combinedef
using Printf: @printf,@sprintf
using ..TypeTraits: efficientoperations

import MacroTools: rmlines

export FExpr,escape
export Escaped,UnEscaped,MixEscaped
export RawExpr,rawexpr
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

"Abstract escape mechanism."
abstract type EscapeMechanism end

"""
    Escaped(names::Symbol...)

Indicate that symbols of a factory should be escaped if they are in `names`.
"""
struct Escaped{N} <: EscapeMechanism
    names::NTuple{N,Symbol}
    Escaped(names::Symbol...)=new{names|>length}(names)
end

"""
    UnEscaped(names::Symbol...)

IIndicate that symbols of a factory should be escaped if they are not in `names`.
"""
struct UnEscaped{N} <: EscapeMechanism
    names::NTuple{N,Symbol}
    UnEscaped(names::Symbol...)=new{names|>length}(names)
end

"""
    MixEscaped(escaped::Escaped)
    MixEscaped(unescaped::UnEscaped)
    MixEscaped(escaped::Escaped,unescaped::UnEscaped)
    MixEscaped(unescaped::UnEscaped,escaped::Escaped)

Indicate that some parts of a factory use the `Escaped` mechanism while other parts use the `UnEscaped` mechanism.
"""
struct MixEscaped{N,M} <: EscapeMechanism
    escaped::Escaped{N}
    unescaped::UnEscaped{M}
end
MixEscaped()=MixEscaped(Escaped(),UnEscaped())
MixEscaped(escaped::Escaped)=MixEscaped(escaped,UnEscaped())
MixEscaped(unescaped::UnEscaped)=MixEscaped(Escaped(),unescaped)
MixEscaped(unescaped::UnEscaped,escaped::Escaped)=MixEscaped(escaped,unescaped)

"Raw expression without any variable escaped."
struct RawExpr <: EscapeMechanism end
"""
    rawexpr

Indicate that no variable in a factory should be escaped.
"""
const rawexpr=RawExpr()

"""
    escape(expr,::RawExpr) -> Any
    escape(expr,::Escaped) -> Any
    escape(expr,::UnEscaped) -> Any
    escape(expr::Symbol,em::Escaped) -> FExpr
    escape(expr::Expr,em::Escaped) -> Expr
    escape(expr::Symbol,em::UnEscaped) -> FExpr
    escape(expr::Expr,em::UnEscaped) -> Expr

Escape the variables in the input expression.
"""
escape(expr,::RawExpr)=expr
escape(expr,::Escaped)=expr
escape(expr,::UnEscaped)=expr
escape(expr::Symbol,em::Escaped)=expr ∈ em.names ? :($(esc(expr))) : expr
escape(expr::Expr,em::Escaped)=(expr.head==:escape ? expr : Expr(expr.head,(escape(arg,em) for arg in expr.args)...))
escape(expr::Symbol,em::UnEscaped)=expr ∉ em.names ? :($(esc(expr))) : expr
function escape(expr::Expr,em::UnEscaped)
    if expr.head==:escape
        @assert length(expr.args)==1 "escape error: internal error..."
        escape(expr.args[1],em)
    else
        Expr(expr.head,(escape(arg,em) for arg in expr.args)...)
    end
end

"""
    names(expr::Symbol) -> Vector{Symbol}
    names(expr::Expr) -> Vector{Symbol}

Get all the symbols in an expression.
"""
Base.names(expr::Symbol)=[expr]
function Base.names(expr::Expr)
    result=Symbol[]
    for arg in expr.args
        (isa(arg,Symbol)|| isa(arg,Expr)) && append!(result,names(arg))
    end
    return result
end

"""
    AbstractFactory

Abstract type for all concrete factories.
"""
abstract type AbstractFactory end

"""
    ==(f1::F,f2::F) where F<:AbstractFactory -> Bool
    isequal(f1::F,f2::F) where F<:AbstractFactory -> Bool

Overloaded equivalent operator.
"""
Base.:(==)(f1::F,f2::F) where F<:AbstractFactory = ==(efficientoperations,f1,f2)
Base.isequal(f1::F,f2::F) where F<:AbstractFactory=isequal(efficientoperations,f1,f2)

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
            vstr=string(value(rawexpr))
        elseif isa(value,Vector{<:AbstractFactory})
            eltypename=value|>eltype|>nameof
            strings=Tuple(string(af(rawexpr)) for af in value)
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
Base.replace(f::AbstractFactory;kwargs...)=replace(efficientoperations,f;kwargs...)

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
    elseif expr.head==:(<:)
        @assert length(expr.args)==1 "Inference error: wrong input expr."
        return Inference(expr.head,nothing,Inference(expr.args[1]))
    elseif expr.head==:curly
        return Inference(expr.head,expr.args[1],Inference.(@views(expr.args[2:end])))
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
    (i::Inference)(em::RawExpr) -> FExpr
    (i::Inference)(em::UnEscaped) -> FExpr
    (i::Inference)(em::MixEscaped) -> FExpr

Convert a `Inference` to the `Expr` representation of the type inference it describes.
"""
(i::Inference)(em::MixEscaped)=i(em.unescaped)
function (i::Inference)(em::Union{RawExpr,<:UnEscaped})
    if i.head===nothing
        return escape(i.name,em)
    elseif i.head==:(<:)
        return Expr(i.head,i.params(em))
    else
        return Expr(i.head,escape(i.name,em),[param(em) for param in i.params]...)
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
    (a::Argument)(em::RawExpr) -> Expr
    (a::Argument)(em::MixEscaped) -> Expr

Convert an `Argument` to the `Expr` representation of the argument it describes.
"""
(a::Argument)(em::RawExpr)=combinearg(a.name,a.type(em),a.slurp,escape(a.default,em))
(a::Argument)(em::MixEscaped)=combinearg(a.name,a.type(em.unescaped),a.slurp,escape(a.default,em.escaped))

"""
    Parameter(head::Union{Symbol,Nothing},name::Union{Symbol,Nothing},type::Union{Inference,Nothing})
    Parameter(;head::Union{Symbol,Nothing}=nothing,name::Union{Symbol,Nothing}=nothing,type::Union{Inference,Nothing}=nothing)
    Parameter(expr::FExpr)

The struct to describe a parameter of a `function` or a `type`.
"""
mutable struct Parameter <: AbstractFactory
    head::Union{Symbol,Nothing}
    name::Union{Symbol,Nothing}
    type::Union{Inference,Nothing}
end
Parameter(;head::Union{Symbol,Nothing}=nothing,name::Union{Symbol,Nothing}=nothing,type::Union{Inference,Nothing}=nothing)=Parameter(head,name,type)
function Parameter(expr::FExpr)
    if isa(expr,Symbol)
        head,name,type=nothing,expr,nothing
    elseif expr.head==:(<:)
        head=:(<:)
        (name,type)=length(expr.args)==2 ? expr.args : (nothing,expr.args[1])
        type=Inference(type)
    elseif expr.head==:(::)
        @assert length(expr.args)==1 "Parameter error: wrong input expr."
        head,name,type=:(::),nothing,Inference(expr.args[1])
    else
        error("Parameter error: wrong input expr.")
    end
    Parameter(head,name,type)
end

"""
    @parameter expr::FExpr

Construct a `Parameter` directly from an parameter statement.
"""
macro parameter(expr::FExpr) expr=[expr];:(Parameter($expr...)) end

"""
    (p::Parameter)(em::RawExpr) -> FExpr
    (p::Parameter)(em::UnEscaped) -> FExpr
    (p::Parameter)(em::MixEscaped) -> FExpr

Convert a `Parameter` to the `Expr` representation of the parameter it describes.
"""
function (p::Parameter)(em::Union{RawExpr,<:UnEscaped,<:MixEscaped})
    p.name===nothing ? (p.head==:(::) ? p.type(em) : Expr(:(<:),p.type(em))) : p.type===nothing ? p.name : Expr(:(<:),p.name,p.type(em))
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
Field(;name::Symbol,type::Inference=Inference(:Any))=Field(name,type)
Field(expr::FExpr)=(cache=splitarg(expr);Field(cache[1],Inference(cache[2])))

"""
    @field expr::FExpr

Construct a `Field` directly from a field statement.
"""
macro field(expr::FExpr) expr=[expr];:(Field($expr...)) end

"""
    (f::Field)(em::RawExpr) -> Expr
    (f::Field)(em::UnEscaped) -> Expr
    (f::Field)(em::MixEscaped) -> Expr

Convert a `Field` to the `Expr` representation of the field it describes.
"""
(f::Field)(em::Union{RawExpr,<:UnEscaped,<:MixEscaped})=combinefield((f.name,f.type(em)))

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
    (b::Block)(em::RawExpr) -> Expr
    (b::Block)(em::Escaped) -> Expr
    (b::Block)(em::MixEscaped) -> Expr

Convert a `Block` to the `Expr` representation of the `begin ... end` block it describes.
"""
(b::Block)(em::MixEscaped)=b(em.escaped)
(b::Block)(em::Union{RawExpr,<:Escaped})=Expr(:block,(escape(part,em) for part in b.body)...)

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
    (ff::FunctionFactory)(em::RawExpr) -> Expr
    (ff::FunctionFactory)(em::MixEscaped) -> Expr

Convert a `FunctionFactory` to the `Expr` representation of the `function` it describes.
"""
function (ff::FunctionFactory)(em::Union{RawExpr,<:MixEscaped})
    combinedef(Dict(
        :name           =>      escape(ff.name,isa(em,RawExpr) ? em : em.escaped),
        :params         =>      [param(em) for param in ff.params],
        :args           =>      [arg(em) for arg in ff.args],
        :kwargs         =>      [kwarg(em) for kwarg in ff.kwargs],
        :rtype          =>      ff.rtype(em),
        :whereparams    =>      [whereparam(em) for whereparam in ff.whereparams],
        :body           =>      ff.body(em)
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
    (tf::TypeFactory)(em::RawExpr) -> Expr
    (tf::TypeFactory)(em::MixEscaped) -> Expr

Convert a `TypeFactory` to the `Expr` representation of the `struct` it describes.
"""
function (tf::TypeFactory)(em::Union{RawExpr,<:MixEscaped})
    combinestructdef(Dict(
        :name           =>      escape(tf.name,isa(em,RawExpr) ? em : em.escaped),
        :mutable        =>      tf.mutable,
        :params         =>      [param(em) for param in tf.params],
        :supertype      =>      tf.supertype(em),
        :fields         =>      [(field.name,field.type(em)) for field in tf.fields],
        :constructors   =>      [constructor(em) for constructor in tf.constructors]
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
addparams!(f::FunctionFactory)=f
addparams!(f::FunctionFactory,params::Inference...)=(push!(f.params,params...);f)
addparams!(f::FunctionFactory,params::FExpr...)=addparams!(f,Inference.(params)...)
addparams!(f::TypeFactory)=f
addparams!(f::TypeFactory,params::Parameter...)=(push!(f.params,params...);f)
addparams!(f::TypeFactory,params::FExpr...)=addparams!(f,Parameter.(params)...)

"""
    @addparams! f params::FExpr...

Add a couple of method parameters to a function factory or a type factory.
"""
macro addparams!(f,params::FExpr...) :(addparams!($(esc(f)),$params...)) end

end #module
