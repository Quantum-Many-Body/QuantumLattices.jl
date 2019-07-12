module Frameworks

using Printf: @printf
using TimerOutputs: TimerOutputs,TimerOutput,@timeit
using ..Terms: Parameters
using ...Prerequisites: decimaltostr
using ...Prerequisites.CompositeStructures: NamedContainer
using ...Prerequisites.TypeTraits: efficientoperations

import ...Interfaces: id,update!,prepare!,register!,run!,add!

export App,Engine
export Assignment,Algorithm
export dependences,rundependences!

"""
    App

Abstract type for all apps.
"""
abstract type App end

"""
    ==(app1::App,app2::App) -> Bool
    isequal(app1::App,app2::App) -> Bool

Judge whether two apps are equivalent to each other.
"""
Base.:(==)(app1::App,app2::App)= ==(efficientoperations,app1,app2)
Base.isequal(app1::App,app2::App)=isequal(efficientoperations,app1,app2)

"""
    update!(app::App;kwargs...) -> App

Update the status of an app.
"""
update!(app::App;kwargs...)=app

"""
    Engine

Abstract type for all engines.
"""
abstract type Engine end

"""
    ==(engine1::Engine,engine2::Engine) -> Bool
    isequal(engine1::Engine,engine2::Engine) -> Bool

Judge whether two engines are equivalent to each other.
"""
Base.:(==)(engine1::Engine,engine2::Engine)= ==(efficientoperations,engine1,engine2)
Base.isequal(engine1::Engine,engine2::Engine)=isequal(efficientoperations,engine1,engine2)

"""
    update!(engine::Engine;kwargs...) -> Engine

Update the status of an engine.
"""
update!(engine::Engine;kwargs...)=engine

"""
    Assignment( id::Symbol,app::App,parameters::Parameters;
                map::Function=identity,
                dependences::Tuple{Vararg{Symbol}}=(),
                data::Any=nothing,
                savedata::Bool=true,
                virgin::Bool=true,
                kwargs...
                )

An assignment associated with an app.
"""
mutable struct Assignment{A<:App,P<:Parameters,M<:Function,D<:Tuple{Vararg{Symbol}},R<:Any,I}
    app::A
    parameters::P
    map::M
    dependences::D
    data::R
    savedata::Bool
    virgin::Bool
end
function Assignment(id::Symbol,app::App,parameters::Parameters;
                    map::Function=identity,
                    dependences::Tuple{Vararg{Symbol}}=(),
                    data::Any=nothing,
                    savedata::Bool=true,
                    virgin::Bool=true,
                    kwargs...
                    )
    A,P,M,D=app|>typeof,parameters|>typeof,map|>typeof,dependences|>typeof
    R=data===nothing ? Any : data|>typeof
    return Assignment{A,P,M,D,R,id}(app,parameters,map,dependences,data,savedata,virgin)
end

"""
    ==(assign1::Assignment,assign2::Assignment) -> Bool
    isequal(assign1::Assignment,assign2::Assignment) -> Bool

Judge whether two assignments are equivalent to each other.
"""
Base.:(==)(assign1::Assignment,assign2::Assignment)= ==(efficientoperations,assign1,assign2)
Base.isequal(assign1::Assignment,assign2::Assignment)=isequal(efficientoperations,assign1,assign2)

"""
    valtype(assign::Assignment)
    valtype(::Type{<:Assignment{<:App,<:Parameters,<:Function,<:Tuple{Vararg{Symbol}},R}}) where R

The type of the data(result) of an assignment.
"""
Base.valtype(assign::Assignment)=assign|>typeof|>valtype
Base.valtype(::Type{<:Assignment{<:App,<:Parameters,<:Function,<:Tuple{Vararg{Symbol}},R}}) where R=R

"""
    id(assign::Assignment) -> Symbol
    id(::Type{<:Assignment{<:App,<:Parameters,<:Function,<:Tuple{Vararg{Symbol}},<:Any,I}}) where I -> Symbol

The id of an assignment.
"""
id(assign::Assignment)=assign|>typeof|>id
id(::Type{<:Assignment{<:App,<:Parameters,<:Function,<:Tuple{Vararg{Symbol}},<:Any,I}}) where I=I

"""
    update!(assign::Assignment;kwargs...) -> Assignment

Update the parameters of an assignment and the status of its associated app.
"""
@generated function update!(assign::Assignment;kwargs...)
    exprs=[]
    names=fieldnames(fieldtype(assign,:parameters))
    for (i,name) in enumerate(names)
        name=QuoteNode(name)
        push!(exprs,:(get(kwargs,$name,getfield(assign.parameters,$i))))
    end
    return quote
        assign.parameters=Parameters{$names}($(exprs...))
        update!(assign.app;assign.map(assign.parameters)...)
        return assign
    end
end

"""
    Algorithm(  name::String,engine::Engine;
                din::String=".",
                dout::String=".",
                parameters::Union{Parameters,Nothing}=nothing,
                map::Function=identity,
                assignments::Tuple{Vararg{Assignment}}=(),
                kwargs...
                )

An algorithm associated with an engine.
"""
mutable struct Algorithm{E<:Engine,P<:Parameters,M<:Function,S<:NamedContainer{Assignment}}
    name::String
    engine::E
    din::String
    dout::String
    parameters::P
    map::M
    sassignments::S
    dassignments::Dict{Symbol,Assignment}
    timer::TimerOutput
end
function Algorithm( name::String,engine::Engine;
                    din::String=".",
                    dout::String=".",
                    parameters::Union{Parameters,Nothing}=nothing,
                    map::Function=identity,
                    assignments::Tuple{Vararg{Assignment}}=(),
                    kwargs...
                    )
    parameters===nothing && (parameters=Parameters(engine))
    assignments=namedassignments(assignments)
    return Algorithm(name,engine,din,dout,parameters,map,assignments,Dict{Symbol,Assignment}(),TimerOutput())
end
@generated function namedassignments(assignments::Tuple{Vararg{Assignment}})
    names,values=[],[]
    for i=1:fieldcount(assignments)
        push!(names,fieldtype(assignments,1)|>id)
        push!(values,:(assignments[$i]))
    end
    names=NTuple{fieldcount(assignments),Symbol}(names)
    return :(NamedContainer{$names}($(values...)))
end

"""
    update!(alg::Algorithm;kwargs...) -> Algorithm

Update the parameters of an algorithm and its associated engine.
"""
@generated function update!(alg::Algorithm;kwargs...)
    exprs=[]
    names=fieldnames(fieldtype(alg,:parameters))
    for (i,name) in enumerate(names)
        name=QuoteNode(name)
        push!(exprs,:(get(kwargs,$name,getfield(alg.parameters,$i))))
    end
    return quote
        alg.parameters=Parameters{$names}($(exprs...))
        update!(alg.engine;alg.map(alg.parameters)...)
        return alg
    end
end

"""
    repr(alg::Algorithm,mask::Tuple{Vararg{Symbol}}=();ndecimal::Int=10) -> String

Get the repr representation of an algorithm.

Optionally, some parameters of the algorithm can be masked. Besides, the maximum number of decimals of the parameters can also be specified.
"""
function Base.repr(alg::Algorithm,mask::Tuple{Vararg{Symbol}}=();ndecimal::Int=10)
    result=[string(alg.name),repr(alg.engine)]
    for (name,value) in pairs(alg.parameters)
        name∉mask && push!(result,decimaltostr(value,ndecimal))
    end
    return join(result,"_")
end

"""
    show(io::IO,alg::Algorithm)

Show an algorithm.
"""
function Base.show(io::IO,alg::Algorithm)
    @printf io "%s_%s" alg.name alg.engine
    for (name,value) in pairs(alg.parameters)
        @printf io "_%s" decimaltostr(value,10)
    end
end

"""
    get(alg::Algorithm,id::Symbol) -> Assignment
    get(alg::Algorithm,::Val{id}) where id -> Assignment

Find the assignment registered on a algorithm by its id.
"""
Base.get(alg::Algorithm,id::Symbol)=id∈keys(alg.sassignments) ? getfield(alg.sassignments,id) : alg.dassignments[id]
@generated Base.get(alg::Algorithm,::Val{id}) where id=id∈fieldnames(fieldtype(alg,:sassignments)) ? :(getfield(alg.sassignments,id)) : :(alg.dassignments[id])

"""
    summary(alg::Algorithm)

Provide a summary of an algorithm.
"""
function Base.summary(alg::Algorithm)
    @info "Summary of $(alg.name)($(nameof(typeof(alg.engine)))):"
    @info string(alg.timer)
end

"""
    prepare!(alg::Algorithm,assign::Assignment) -> Nothing

Prepare an assignment registered on a algorithm.
"""
prepare!(alg::Algorithm,assign::Assignment)=nothing

"""
    run!(alg::Algorithm,assign::Assignment) -> Nothing

Run an assignment registered on a algorithm.
"""
run!(alg::Algorithm,assign::Assignment)=nothing

"""
    register!(alg::Algorithm,id::Symbol,app::App;kwargs...) -> Algorithm

Add an assignment on a algorithm by providing the contents of the assignment, and run this assignment.
"""
function register!(alg::Algorithm,id::Symbol,app::App;kwargs...)
    add!(alg,id,app;kwargs...)
    run!(alg,id,true)
end

"""
    add!(alg::Algorithm,id::Symbol,app::App;kwargs...) -> Algorithm

Add an assignment on a algorithm by providing the contents of the assignment.

The difference between `add!` and `register!` is that the `add!` function does not run the newly added assignment but the `register!` function does.
"""
function add!(alg::Algorithm,id::Symbol,app::App;kwargs...)
    @assert id ∉ keys(alg.sassignments) "add! error: id($id) conflict."
    alg.dassignments[id]=Assignment(id,app,merge(alg.parameters,get(kwargs,:parameters,Parameters{()}()));kwargs...)
    return alg
end

"""
    run!(alg::Algorithm,id::Symbol,timing::Bool=true) -> Algorithm
    run!(alg::Algorithm,::Val{id},timing::Bool=true) where id -> Algorithm

Run an assignment with the given id registered on an algorithm. Optionally, the run process can be timed by setting the `timing` argument to be `true`.
"""
run!(alg::Algorithm,id::Symbol,timing::Bool=true)=run!(alg,Val(id),timing)
function run!(alg::Algorithm,::Val{id},timing::Bool=true) where id
    assign=get(alg,Val(id))
    if timing
        @timeit alg.timer string(id) algrunassignment!(alg,assign)
        @info "App $id($(nameof(assign.app|>typeof))): time consumed $(TimerOutputs.time(alg.timer[string(id)])/10^9)s."
    else
        algrunassignment!(alg,assign)
    end
    return alg
end
function algrunassignment!(alg::Algorithm,assign::Assignment)
    ismatched=match(assign.parameters,alg.parameters)
    !assign.virgin && ismatched && return
    !ismatched && update!(alg;assign.parameters...)
    prepare!(alg,assign)
    run!(alg,assign)
    assign.virgin=false
end

"""
    dependences(alg::Algorithm,assign::Assignment,::Tuple{}=()) -> Tuple{Vararg{Symbol}}
    dependences(alg::Algorithm,assign::Assignment,mask::Tuple{Vararg{Symbol}}) -> Tuple{Vararg{Symbol}}

Get the dependences of an assignment and return their ids.
"""
dependences(alg::Algorithm,assign::Assignment,::Tuple{}=())=assign.dependences
dependences(alg::Algorithm,assign::Assignment,mask::Tuple{Vararg{Symbol}})=Tuple(filter(x->x∉mask,collect(assign.dependences)))

"""
    rundependences!(alg::Algorithm,assign::Assignment,mask::Tuple{Vararg{Symbol}}=()) -> Algorithm

Run the dependences of an assignment. Optionally, some dependences can be jumped by specifying the `mask` argument.
"""
function rundependences!(alg::Algorithm,assign::Assignment,mask::Tuple{Vararg{Symbol}}=())
    for id in dependences(alg,assign,mask)
        assign=get(alg,id)
        algrundependence!(alg,assign)
    end
    return alg
end
function algrundependence!(alg::Algorithm,assign::Assignment)
    ismatched=match(alg.parameters,assign.parameters)
    !assign.virgin && ismatched && return
    !ismatched && update!(assign;alg.parameters...)
    prepare!(alg,assign)
    run!(alg,assign)
    assign.virgin=false
end

end  #module
