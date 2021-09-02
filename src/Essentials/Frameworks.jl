module Frameworks

using Printf: @printf
using TimerOutputs: TimerOutputs, TimerOutput, @timeit
using ..Terms: Parameters
using ...Prerequisites: decimaltostr
using ...Prerequisites.CompositeStructures: NamedContainer
using ...Prerequisites.Traits: efficientoperations

import ...Interfaces: id, add!
import ...Essentials: update!

export App, Engine, Assignment, Algorithm
export prepare!, register!, run!, dependences, rundependences!

"""
    App

Abstract type for all apps.
"""
abstract type App end
@inline Base.:(==)(app₁::App, app₂::App) = ==(efficientoperations, app₁, app₂)
@inline Base.isequal(app₁::App, app₂::App) = isequal(efficientoperations, app₁, app₂)

"""
    update!(app::App; kwargs...) -> App

Update the status of an app.
"""
@inline update!(app::App; kwargs...) = app

"""
    Engine

Abstract type for all engines.
"""
abstract type Engine end
@inline Base.:(==)(engine₁::Engine, engine₂::Engine) = ==(efficientoperations, engine₁, engine₂)
@inline Base.isequal(engine₁::Engine, engine₂::Engine) = isequal(efficientoperations, engine₁, engine₂)

"""
    update!(engine::Engine; kwargs...) -> Engine

Update the status of an engine.
"""
@inline update!(engine::Engine; kwargs...) = engine

"""
    Assignment{A<:App, P<:Parameters, M<:Function, D<:Tuple{Vararg{Symbol}}, R<:Any, I}

An assignment associated with an app.
"""
mutable struct Assignment{A<:App, P<:Parameters, M<:Function, D<:Tuple{Vararg{Symbol}}, R<:Any, I}
    app::A
    parameters::P
    map::M
    dependences::D
    data::R
    savedata::Bool
    virgin::Bool
end
@inline Base.:(==)(assign₁::Assignment, assign₂::Assignment) = ==(efficientoperations, assign₁, assign₂)
@inline Base.isequal(assign₁::Assignment, assign₂::Assignment) = isequal(efficientoperations, assign₁, assign₂)

"""
    Assignment(id::Symbol, app::App, parameters::Parameters;
        map::Function=identity,
        dependences::Tuple{Vararg{Symbol}}=(),
        data::Any=nothing,
        savedata::Bool=true,
        virgin::Bool=true,
        kwargs...
        )

Construct an assignment.
"""
@inline function Assignment(id::Symbol, app::App, parameters::Parameters;
        map::Function=identity,
        dependences::Tuple{Vararg{Symbol}}=(),
        data::Any=nothing,
        savedata::Bool=true,
        virgin::Bool=true,
        kwargs...)
    A, P, M, D, R = typeof(app), typeof(parameters), typeof(map), typeof(dependences), isnothing(data) ? Any : typeof(data)
    return Assignment{A, P, M, D, R, id}(app, parameters, map, dependences, data, savedata, virgin)
end

"""
    valtype(assign::Assignment)
    valtype(::Type{<:Assignment})

The type of the data(result) of an assignment.
"""
@inline Base.valtype(assign::Assignment) = valtype(typeof(assign))
@inline Base.valtype(::Type{<:Assignment{<:App, <:Parameters, <:Function, <:Tuple{Vararg{Symbol}}, R}}) where R = R

"""
    id(assign::Assignment) -> Symbol
    id(::Type{<:Assignment}) -> Symbol

The id of an assignment.
"""
@inline id(assign::Assignment) = id(typeof(assign))
@inline id(::Type{<:Assignment{<:App, <:Parameters, <:Function, <:Tuple{Vararg{Symbol}}, <:Any, I}}) where I = I

"""
    update!(assign::Assignment; kwargs...) -> Assignment

Update the parameters of an assignment and the status of its associated app.
"""
@generated function update!(assign::Assignment; kwargs...)
    exprs = []
    names = fieldnames(fieldtype(assign, :parameters))
    for (i, name) in enumerate(names)
        name = QuoteNode(name)
        push!(exprs, :(get(kwargs, $name, getfield(assign.parameters, $i))))
    end
    return quote
        assign.parameters = Parameters{$names}($(exprs...))
        update!(assign.app; assign.map(assign.parameters)...)
        return assign
    end
end

"""
    Algorithm{E<:Engine, P<:Parameters, M<:Function, S<:NamedContainer{Assignment}}

An algorithm associated with an engine.
"""
mutable struct Algorithm{E<:Engine, P<:Parameters, M<:Function, S<:NamedContainer{Assignment}}
    name::String
    engine::E
    din::String
    dout::String
    parameters::P
    map::M
    sassignments::S
    dassignments::Dict{Symbol, Assignment}
    timer::TimerOutput
end
@inline Base.:(==)(algorithm₁::Algorithm, algorithm₂::Algorithm) = ==(efficientoperations, algorithm₁, algorithm₂)
@inline Base.isequal(algorithm₁::Algorithm, algorithm₂::Algorithm) = isequal(efficientoperations, algorithm₁, algorithm₂)
function Base.show(io::IO, alg::Algorithm)
    @printf io "%s_%s" alg.name alg.engine
    for (name, value) in pairs(alg.parameters)
        @printf io "_%s" decimaltostr(value, 10)
    end
end

"""
    Algorithm(name::String, engine::Engine;
        din::String=".",
        dout::String=".",
        parameters::Union{Parameters, Nothing}=nothing,
        map::Function=identity,
        assignments::Tuple{Vararg{Assignment}}=(),
        kwargs...
        )

Construct an algorithm.
"""
@inline function Algorithm(name::String, engine::Engine;
        din::String=".",
        dout::String=".",
        parameters::Parameters=Parameters(engine),
        map::Function=identity,
        assignments::Tuple{Vararg{Assignment}}=(),
        kwargs...
        )
    return Algorithm(name, engine, din, dout, parameters, map, namedassignments(assignments), Dict{Symbol, Assignment}(), TimerOutput())
end
@generated function namedassignments(assignments::Tuple{Vararg{Assignment}})
    names, values = [], []
    for i = 1:fieldcount(assignments)
        push!(names, fieldtype(assignments, i)|>id)
        push!(values, :(assignments[$i]))
    end
    names = NTuple{fieldcount(assignments), Symbol}(names)
    values = Expr(:tuple, values...)
    return :(NamedContainer{$names}($values))
end

"""
    update!(alg::Algorithm; kwargs...) -> Algorithm

Update the parameters of an algorithm and its associated engine.
"""
@generated function update!(alg::Algorithm; kwargs...)
    exprs = []
    names = fieldnames(fieldtype(alg, :parameters))
    for (i, name) in enumerate(names)
        name = QuoteNode(name)
        push!(exprs, :(get(kwargs, $name, getfield(alg.parameters, $i))))
    end
    return quote
        alg.parameters = Parameters{$names}($(exprs...))
        update!(alg.engine; alg.map(alg.parameters)...)
        return alg
    end
end

"""
    repr(alg::Algorithm, mask::Tuple{Vararg{Symbol}}=(); ndecimal::Int=10) -> String

Get the repr representation of an algorithm.

Optionally, some parameters of the algorithm can be masked. Besides, the maximum number of decimals of the parameters can also be specified.
"""
function Base.repr(alg::Algorithm, mask::Tuple{Vararg{Symbol}}=(); ndecimal::Int=10)
    result = [string(alg.name), repr(alg.engine)]
    for (name, value) in pairs(alg.parameters)
        name∉mask && push!(result, decimaltostr(value, ndecimal))
    end
    return join(result, "_")
end

"""
    get(alg::Algorithm, id::Symbol) -> Assignment
    get(alg::Algorithm, ::Val{id}) where id -> Assignment

Find the assignment registered on a algorithm by its id.
"""
@inline Base.get(alg::Algorithm, id::Symbol) = get(alg, id|>Val)
@generated function Base.get(alg::Algorithm, ::Val{id}) where id
    return id∈fieldnames(fieldtype(alg, :sassignments)) ? :(getfield(alg.sassignments, id)) : :(alg.dassignments[id])
end

"""
    summary(alg::Algorithm)

Provide a summary of an algorithm.
"""
function Base.summary(alg::Algorithm)
    @info "Summary of $(alg.name)($(nameof(typeof(alg.engine)))):"
    @info string(alg.timer)
end

"""
    prepare!(alg::Algorithm, assign::Assignment) -> Nothing

Prepare an assignment registered on a algorithm.
"""
@inline prepare!(alg::Algorithm, assign::Assignment) = nothing

"""
    run!(alg::Algorithm, assign::Assignment) -> Nothing

Run an assignment registered on a algorithm.
"""
@inline run!(alg::Algorithm, assign::Assignment) = nothing

"""
    register!(alg::Algorithm, id::Symbol, app::App; kwargs...) -> Algorithm

Add an assignment on a algorithm by providing the contents of the assignment, and run this assignment.
"""
@inline function register!(alg::Algorithm, id::Symbol, app::App; kwargs...)
    add!(alg, id, app; kwargs...)
    run!(alg, id, true)
end

"""
    add!(alg::Algorithm, id::Symbol, app::App; kwargs...) -> Algorithm

Add an assignment on a algorithm by providing the contents of the assignment.

The difference between `add!` and `register!` is that the `add!` function does not run the newly added assignment but the `register!` function does.
"""
function add!(alg::Algorithm, id::Symbol, app::App; kwargs...)
    @assert id∉keys(alg.sassignments) "add! error: id($id) conflict."
    alg.dassignments[id] = Assignment(id, app, merge(alg.parameters, get(kwargs, :parameters, Parameters{()}())); kwargs...)
    return alg
end

"""
    run!(alg::Algorithm, id::Symbol, timing::Bool=true) -> Algorithm
    run!(alg::Algorithm, ::Val{id}, timing::Bool=true) where id -> Algorithm

Run an assignment with the given id registered on an algorithm. Optionally, the run process can be timed by setting the `timing` argument to be `true`.
"""
@inline run!(alg::Algorithm, id::Symbol, timing::Bool=true) = run!(alg, Val(id), timing)
function run!(alg::Algorithm, ::Val{id}, timing::Bool=true) where id
    assign = get(alg, Val(id))
    if timing
        @timeit alg.timer string(id) algrunassignment!(alg, assign)
        @info "App $id($(nameof(assign.app|>typeof))): time consumed $(TimerOutputs.time(alg.timer[string(id)]) / 10^9)s."
    else
        algrunassignment!(alg, assign)
    end
    return alg
end
function algrunassignment!(alg::Algorithm, assign::Assignment)
    ismatched = match(assign.parameters, alg.parameters)
    !assign.virgin && ismatched && return
    !ismatched && update!(alg; assign.parameters...)
    prepare!(alg, assign)
    run!(alg, assign)
    assign.virgin = false
end

"""
    dependences(alg::Algorithm, assign::Assignment, ::Tuple{}=()) -> Tuple{Vararg{Symbol}}
    dependences(alg::Algorithm, assign::Assignment, mask::Tuple{Vararg{Symbol}}) -> Tuple{Vararg{Symbol}}

Get the dependences of an assignment and return their ids.
"""
@inline dependences(alg::Algorithm, assign::Assignment, ::Tuple{}=()) = assign.dependences
@inline dependences(alg::Algorithm, assign::Assignment, mask::Tuple{Vararg{Symbol}}) = Tuple(filter(x->x∉mask, collect(assign.dependences)))

"""
    rundependences!(alg::Algorithm, assign::Assignment, mask::Tuple{Vararg{Symbol}}=()) -> Algorithm

Run the dependences of an assignment. Optionally, some dependences can be jumped by specifying the `mask` argument.
"""
function rundependences!(alg::Algorithm, assign::Assignment, mask::Tuple{Vararg{Symbol}}=())
    for id in dependences(alg, assign, mask)
        assign = get(alg, id)
        algrundependence!(alg, assign)
    end
    return alg
end
function algrundependence!(alg::Algorithm, assign::Assignment)
    ismatched = match(alg.parameters, assign.parameters)
    !assign.virgin && ismatched && return
    !ismatched && update!(assign; alg.parameters...)
    prepare!(alg, assign)
    run!(alg, assign)
    assign.virgin = false
end

end  #module
