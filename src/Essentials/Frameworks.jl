module Frameworks

using Printf: @printf, @sprintf
using LaTeXStrings: latexstring
using Serialization: serialize
using TimerOutputs: TimerOutputs, TimerOutput, @timeit
using RecipesBase: RecipesBase, @recipe, @series
using ..QuantumOperators: Operators, Transformation, idtype, optype
using ..Spatials: Bonds, AbstractLattice, acrossbonds, isintracell
using ..DegreesOfFreedom: Hilbert, Table, Boundary, plain, Term, ismodulatable
using ...Prerequisites: atol, rtol, decimaltostr
using ...Prerequisites.Traits: efficientoperations
using ...Prerequisites.CompositeStructures: NamedContainer

import ...Interfaces: id, add!, expand, expand!
import ...Essentials: update, update!, reset!
import ...Prerequisites.Traits: contentnames, getcontent

export Parameters, Entry, Engine, AbstractGenerator, Generator, SimplifiedGenerator, Action, Assignment, Algorithm
export prepare!, run!, save, rundependences!

"""
    Parameters{Names}(values::Number...) where Names

A NamedTuple that contain the key-value pairs.
"""
const Parameters{Names} = NamedTuple{Names, <:Tuple{Vararg{Number}}}
@inline Parameters{Names}(values::Number...) where {Names} = NamedTuple{Names}(values)
@inline @generated function update(params::Parameters; kwargs...)
    names = fieldnames(params)
    values = [:(get(kwargs, $name, getfield(params, $name))) for name in QuoteNode.(names)]
    return :(Parameters{$names}($(values...)))
end

"""
    match(params₁::Parameters, params₂::Parameters; atol=atol, rtol=rtol) -> Bool

Judge whether the second set of parameters matches the first.
"""
function Base.match(params₁::Parameters, params₂::Parameters; atol=atol, rtol=rtol)
    for name in keys(params₂)
        haskey(params₁, name) && !isapprox(getfield(params₁, name), getfield(params₂, name); atol=atol, rtol=rtol) && return false
    end
    return true
end

"""
    Entry{C, A<:NamedTuple, B<:NamedTuple}

An entry of quantum operators related to (part of) a quantum lattice system.
"""
struct Entry{C, A<:NamedTuple, B<:NamedTuple}
    constops::C
    alterops::A
    boundops::B
end
@inline Base.:(==)(entry₁::Entry, entry₂::Entry) = ==(efficientoperations, entry₁, entry₂)
@inline Base.isequal(entry₁::Entry, entry₂::Entry) = isequal(efficientoperations, entry₁, entry₂)

"""
    valtype(entry::Entry)
    valtype(::Type{<:Entry})

Get the valtype of an entry of quantum operators.
"""
@inline Base.valtype(entry::Entry) = valtype(typeof(entry))
@inline @generated function Base.valtype(::Type{<:Entry{C, A, B}}) where {C, A<:NamedTuple, B<:NamedTuple}
    optp = C
    (fieldcount(A) > 0) && (optp = reduce(promote_type, fieldtypes(A), init=optp))
    (fieldcount(B) > 0) && (optp = reduce(promote_type, fieldtypes(B), init=optp))
    return optp
end

"""
    empty!(entry::Entry) -> Entry

Empty an entry of quantum operators.
"""
function Base.empty!(entry::Entry)
    empty!(entry.constops)
    map(empty!, values(entry.alterops))
    map(empty!, values(entry.boundops))
    return entry
end

"""
    empty(entry::Entry) -> Entry

Get an empty entry of quantum operators.
"""
function Base.empty(entry::Entry)
    constops = empty(entry.constops)
    alterops = NamedTuple{keys(entry.alterops)}(map(empty, values(entry.alterops)))
    boundops = NamedTuple{keys(entry.boundops)}(map(empty, values(entry.boundops)))
    return Entry(constops, alterops, boundops)
end

"""
    merge!(entry::Entry, another::Entry) -> Entry

Merge the entry of quantum operators by another one.
"""
function Base.merge!(entry::Entry, another::Entry)
    merge!(entry.constops, another.constops)
    for name in fieldnames(typeof(entry.alterops))
        merge!(getfield(entry.alterops, name), getfield(another.alterops, name))
    end
    for name in fieldnames(typeof(entry.boundops))
        merge!(getfield(entry.boundops, name), getfield(another.boundops, name))
    end
    return entry
end

"""
    merge(entry::Entry, another::Entry) -> Entry

Get a new entry of quantum operators by merging two entries.
"""
@inline Base.merge(entry::Entry, another::Entry) = merge!(empty(entry), another)

"""
    (transformation::Transformation)(entry::Entry; kwargs...) -> Entry

Get the transformed entry of quantum operators.
"""
function (transformation::Transformation)(entry::Entry; kwargs...)
    wrapper(m) = transformation(m; kwargs...)
    constops = wrapper(entry.constops)
    alterops = NamedTuple{keys(entry.alterops)}(map(wrapper, values(entry.alterops)))
    boundops = NamedTuple{keys(entry.boundops)}(map(wrapper, values(entry.boundops)))
    return Entry(constops, alterops, boundops)
end

"""
    Entry(terms::Tuple{Vararg{Term}}, bonds::Bonds, hilbert::Hilbert; half::Bool=false, table::Union{Nothing, Table}=nothing)

Construct an entry of operators based on the input terms, bonds and Hilbert space.
"""
function Entry(terms::Tuple{Vararg{Term}}, bonds::Bonds, hilbert::Hilbert; half::Bool=false, table::Union{Nothing, Table}=nothing)
    constterms, alterterms, choosedterms = termclassifier(terms)
    innerbonds = filter(acrossbonds, bonds, Val(:exclude))
    boundbonds = filter(acrossbonds, bonds, Val(:include))
    constops = Operators{mapreduce(term->optype(typeof(term), typeof(hilbert), eltype(bonds)), promote_type, choosedterms)}()
    map(term->expand!(constops, term, innerbonds, hilbert, half=half, table=table), constterms)
    alterops = NamedTuple{map(id, alterterms)}(map(term->expand(one(term), innerbonds, hilbert, half=half, table=table), alterterms))
    boundops = NamedTuple{map(id, terms)}(map(term->expand(one(term), boundbonds, hilbert, half=half, table=table), terms))
    return Entry(constops, alterops, boundops)
end
@generated function termclassifier(terms::Tuple{Vararg{Term}})
    constterms, alterterms = [], []
    for (i, term) in enumerate(fieldtypes(terms))
        ismodulatable(term) ? push!(alterterms, :(terms[$i])) : push!(constterms, :(terms[$i]))
    end
    constterms, alterterms = Expr(:tuple, constterms...), Expr(:tuple, alterterms...)
    return Expr(:tuple, constterms, alterterms, (length(constterms.args)>0 ? constterms : alterterms))
end

"""
    expand!(operators::Operators, entry::Entry, boundary::Boundary; kwargs...) -> Operators

Expand an entry of operators with the given boundary twist and term coefficients.
"""
@generated function expand!(operators::Operators, entry::Entry, boundary::Boundary; kwargs...)
    exprs = [:(add!(operators, entry.constops))]
    for name in QuoteNode.(fieldnames(fieldtype(entry, :alterops)))
        push!(exprs, :(value = get(kwargs, $name, nothing)))
        push!(exprs, :(for opt in getfield(entry.alterops, $name)
            add!(operators, opt*value)
        end))
    end
    for name in QuoteNode.(fieldnames(fieldtype(entry, :boundops)))
        push!(exprs, :(value = get(kwargs, $name, nothing)))
        push!(exprs, :(for opt in getfield(entry.boundops, $name)
            add!(operators, boundary(opt)*value)
        end))
    end
    push!(exprs, :(return operators))
    return Expr(:block, exprs...)
end

"""
    reset!(entry::Entry{<:Operators}, terms::Tuple{Vararg{Term}}, bonds::Bonds, hilbert::Hilbert; half::Bool=false, table::Union{Nothing, Table}=nothing) -> Entry

Reset an entry of operators by the new terms, bonds and Hilbert space.
"""
function reset!(entry::Entry{<:Operators}, terms::Tuple{Vararg{Term}}, bonds::Bonds, hilbert::Hilbert; half::Bool=false, table::Union{Nothing, Table}=nothing)
    empty!(entry)
    constterms, alterterms, _ = termclassifier(terms)
    innerbonds = filter(acrossbonds, bonds, Val(:exclude))
    boundbonds = filter(acrossbonds, bonds, Val(:include))
    map(term->expand!(entry.constops, term, innerbonds, hilbert, half=half, table=table), constterms)
    map(term->expand!(getfield(entry.alterops, id(term)), one(term), innerbonds, hilbert, half=half, table=table), alterterms)
    map(term->expand!(getfield(entry.boundops, id(term)), one(term), boundbonds, hilbert, half=half, table=table), terms)
    return entry
end

"""
    Engine

Abstract type for all engines.
"""
abstract type Engine end
@inline Base.:(==)(engine₁::Engine, engine₂::Engine) = ==(efficientoperations, engine₁, engine₂)
@inline Base.isequal(engine₁::Engine, engine₂::Engine) = isequal(efficientoperations, engine₁, engine₂)
@inline update!(engine::Engine; kwargs...) = engine
@inline Base.repr(engine::Engine) = String(nameof(typeof(engine)))
@inline Base.show(io::IO, engine::Engine) = @printf io "%s" nameof(typeof(engine))

"""
    AbstractGenerator{E<:Entry, T<:Union{Table, Nothing}, B<:Boundary} <: Engine

Abstract generator for quantum operators.

By protocol, a concrete generator should have the following predefined contents:
* `operators::E`: the entry for the generated quantum operators
* `table::T`: the index-sequence table if it is not nothing
* `boundary::B`: boundary twist for the generated operators
"""
abstract type AbstractGenerator{E<:Entry, T<:Union{Table, Nothing}, B<:Boundary} <: Engine end
@inline contentnames(::Type{<:AbstractGenerator}) = (:operators, :table, :boundary)
@inline Base.valtype(gen::AbstractGenerator) = valtype(typeof(gen))
@inline Base.valtype(::Type{<:AbstractGenerator{E}}) where {E<:Entry} = valtype(E)
@inline Base.eltype(gen::AbstractGenerator) = eltype(typeof(gen))
@inline Base.eltype(::Type{<:AbstractGenerator{E}}) where {E<:Entry} = eltype(valtype(E))
@inline Base.IteratorSize(::Type{<:AbstractGenerator}) = Base.SizeUnknown()
@inline function Base.iterate(gen::AbstractGenerator)
    ops = expand(gen)
    index = iterate(ops)
    isnothing(index) && return nothing
    return index[1], (ops, index[2])
end
@inline function Base.iterate(gen::AbstractGenerator, state)
    index = iterate(state[1], state[2])
    isnothing(index) && return nothing
    return index[1], (state[1], index[2])
end
Base.show(io::IO, ::MIME"text/latex", gen::AbstractGenerator) = show(io, MIME"text/latex"(), latexstring(latexstring(expand(gen))))

"""
    expand!(result, gen::AbstractGenerator) -> typeof(result)

Expand the generator.
"""
@inline expand!(result, gen::AbstractGenerator) = expand!(result, getcontent(gen, :operators), getcontent(gen, :boundary); Parameters(gen)...)

"""
    expand(gen::AbstractGenerator) -> valtype(gen)

Expand the generator.
"""
@inline expand(gen::AbstractGenerator) = expand!(zero(valtype(gen)), gen)

"""
    Generator{E<:Entry, TS<:Tuple{Vararg{Term}}, BS<:Bonds, H<:Hilbert, T<:Union{Table, Nothing}, B<:Boundary} <: AbstractGenerator{E, T, B}

A generator of operators based on terms, bonds, Hilbert space, and boundary twist.
"""
struct Generator{E<:Entry, TS<:Tuple{Vararg{Term}}, BS<:Bonds, H<:Hilbert, T<:Union{Table, Nothing}, B<:Boundary} <: AbstractGenerator{E, T, B}
    operators::E
    terms::TS
    bonds::BS
    hilbert::H
    half::Bool
    table::T
    boundary::B
end
@inline contentnames(::Type{<:Generator}) = (:operators, :terms, :bonds, :hilbert, :half, :table, :boundary)

"""
    Generator(terms::Tuple{Vararg{Term}}, bonds::Bonds, hilbert::Hilbert;
        half::Bool=false, table::Union{Table,Nothing}=nothing, boundary::Boundary=plain
        )

Construct a generator of operators.
"""
@inline function Generator(terms::Tuple{Vararg{Term}}, bonds::Bonds, hilbert::Hilbert;
        half::Bool=false, table::Union{Table,Nothing}=nothing, boundary::Boundary=plain
        )
    return Generator(Entry(terms, bonds, hilbert; half=half, table=table), terms, bonds, hilbert, half, table, boundary)
end

"""
    Parameters(gen::Generator) -> Parameters

Get the parameters of the terms of a generator.
"""
@inline Parameters(gen::Generator) = NamedTuple{map(id, gen.terms)}(map(term->term.value, gen.terms))

"""
    expand(gen::Generator, name::Symbol) -> Operators
    expand(gen::Generator, i::Int) -> Operators
    expand(gen::Generator, name::Symbol, i::Int) -> Operators

Expand the operators of a generator:
1) the operators of a specific term;
2) the operators on a specific bond;
3) the operators of a specific term on a specific bond.
"""
function expand(gen::Generator, name::Symbol)
    term = get(gen.terms, Val(name))
    result = Operators{optype(typeof(term), typeof(gen.hilbert), eltype(gen.bonds))}()
    if !ismodulatable(term)
        expand!(result, term, filter(acrossbonds, gen.bonds, Val(:exclude)), gen.hilbert, half=gen.half, table=gen.table)
    else
        for opt in getfield(gen.operators.alterops, name)
            add!(result, opt*term.value)
        end
    end
    for opt in getfield(gen.operators.boundops, name)
        add!(result, gen.boundary(opt)*term.value)
    end
    return result
end
function expand(gen::Generator, i::Int)
    bond = gen.bonds[i]
    result = zero(valtype(gen))
    map(term->expand!(result, term, bond, gen.hilbert, half=gen.half, table=gen.table), gen.terms)
    isintracell(bond) || for opt in result
        result[id(opt)] = gen.boundary(opt)
    end
    return result
end
function expand(gen::Generator, name::Symbol, i::Int)
    bond = gen.bonds[i]
    result = expand(get(gen.terms, Val(name)), bond, gen.hilbert, half=gen.half, table=gen.table)
    isintracell(bond) || for opt in result
        result.contents[id(opt)] = gen.boundary(opt)
    end
    return result
end
@inline @generated function Base.get(terms::Tuple{Vararg{Term}}, ::Val{Name}) where Name
    i = findfirst(isequal(Name), map(id, fieldtypes(terms)))::Int
    return :(terms[$i])
end

"""
    update!(gen::Generator; kwargs...) -> typeof(gen)

Update the coefficients of the terms in a generator.
"""
function update!(gen::Generator; kwargs...)
    update!(gen.boundary; kwargs...)
    map(term->(ismodulatable(term) ? update!(term; kwargs...) : term), gen.terms)
    return gen
end

"""
    empty!(gen::Generator) -> Generator

Empty the :bonds, :hilbert, :table and :operators of a generator.
"""
function Base.empty!(gen::Generator)
    empty!(gen.bonds)
    empty!(gen.hilbert)
    isnothing(gen.table) || empty!(gen.table)
    empty!(gen.operators)
    return gen
end

"""
    empty(gen::Generator) -> Generator

Get an empty copy of a generator.
"""
@inline function Base.empty(gen::Generator)
    Generator(
        empty(gen.operators),
        gen.terms,
        empty(gen.bonds),
        empty(gen.hilbert),
        gen.half,
        isnothing(gen.table) ? nothing : empty(gen.table),
        gen.boundary,
        )
end

"""
    reset!(gen::Generator, lattice::AbstractLattice) -> Generator

Reset a generator by a new lattice.
"""
function reset!(gen::Generator, lattice::AbstractLattice)
    reset!(gen.bonds, lattice)
    reset!(gen.hilbert, lattice.pids)
    isnothing(gen.table) || reset!(gen.table, gen.hilbert)
    reset!(gen.operators, gen.terms, gen.bonds, gen.hilbert, half=gen.half, table=gen.table)
    return gen
end

"""
    SimplifiedGenerator{P<:Parameters, T<:Union{Table, Nothing}, B<:Boundary, E<:Entry} <: AbstractGenerator{T, B, E}

The simplified generator for the entry of quantum operators.
"""
mutable struct SimplifiedGenerator{P<:Parameters, E<:Entry, T<:Union{Table, Nothing}, B<:Boundary} <: AbstractGenerator{E, T, B}
    parameters::P
    operators::E
    table::T
    boundary::B
end
@inline contentnames(::Type{<:SimplifiedGenerator}) = (:parameters, :operators, :table, :boundary)

"""
    SimplifiedGenerator(parameters::Parameters, operators::Entry; table::Union{Table,Nothing}=nothing, boundary::Boundary=plain)

Construct an simplified generator of operators.
"""
@inline function SimplifiedGenerator(parameters::Parameters, operators::Entry; table::Union{Table,Nothing}=nothing, boundary::Boundary=plain)
    return SimplifiedGenerator(parameters, operators, table, boundary)
end

"""
    Parameters(gen::SimplifiedGenerator) -> Parameters

Get the parameters of the terms of a generator.
"""
@inline Parameters(gen::SimplifiedGenerator) = gen.parameters

"""
    empty!(gen::SimplifiedGenerator) -> SimplifiedGenerator

Empty the operators of an simplified generator.
"""
@inline function Base.empty!(gen::SimplifiedGenerator)
    empty!(gen.operators)
    isnothing(gen.table) || empty!(gen.table)
    return gen
end

"""
    empty(gen::SimplifiedGenerator) -> SimplifiedGenerator

Get an empty copy of an simplified generator.
"""
@inline function Base.empty(gen::SimplifiedGenerator)
    return SimplifiedGenerator(gen.parameters, empty(gen.operators), isnothing(gen.table) ? nothing : empty(gen.table), gen.boundary)
end

"""
    update!(gen::SimplifiedGenerator; kwargs...) -> typeof(gen)

Update the parameters of a generator.
"""
function update!(gen::SimplifiedGenerator; kwargs...)
    gen.parameters = update(gen.parameters; kwargs...)
    update!(gen.boundary; kwargs...)
    return gen
end

"""
    reset!(gen::SimplifiedGenerator, operators::Entry; table::Union{Table, Nothing}=nothing) -> SimplifiedGenerator

Reset a generator by a new lattice.
"""
function reset!(gen::SimplifiedGenerator, operators::Entry; table::Union{Table, Nothing}=nothing)
    isnothing(gen.table) || isnothing(table) || merge!(empty!(gen.table), table)
    merge!(empty!(gen.operators), operators)
    return gen
end

"""
    (transformation::Transformation)(gen::AbstractGenerator;
        table::Union{Table, Nothing}=getcontent(gen, :table),
        boundary::Boundary=getcontent(gen, :boundary),
        kwargs...
        ) -> SimplifiedGenerator

Get the result of a transformation on the generator for the entry of quantum operators.
"""
function (transformation::Transformation)(gen::AbstractGenerator;
        table::Union{Table, Nothing}=getcontent(gen, :table),
        boundary::Boundary=getcontent(gen, :boundary),
        kwargs...
        )
    return SimplifiedGenerator(Parameters(gen), transformation(getcontent(gen, :operators); kwargs...), table=table, boundary=boundary)
end

"""
    Action

Abstract type for all actions.
"""
abstract type Action end
@inline Base.:(==)(action₁::Action, action₂::Action) = ==(efficientoperations, action₁, action₂)
@inline Base.isequal(action₁::Action, action₂::Action) = isequal(efficientoperations, action₁, action₂)
@inline prepare!(action::Action, engine::Engine) = nothing
@inline update!(action::Action; kwargs...) = action

"""
    Assignment{A<:Action, P<:Parameters, M<:Function, D<:Tuple{Vararg{Symbol}}, R} <: Function

An assignment associated with an action.
"""
mutable struct Assignment{A<:Action, P<:Parameters, M<:Function, D<:Tuple{Vararg{Symbol}}, R} <: Function
    id::Symbol
    action::A
    parameters::P
    map::M
    dependences::D
    data::R
    ismatched::Bool
    save::Bool
end
@inline Base.:(==)(assign₁::Assignment, assign₂::Assignment) = ==(efficientoperations, assign₁, assign₂)
@inline Base.isequal(assign₁::Assignment, assign₂::Assignment) = isequal(efficientoperations, assign₁, assign₂)

"""
    valtype(assign::Assignment)
    valtype(::Type{<:Assignment})

The type of the data(result) of an assignment.
"""
@inline Base.valtype(assign::Assignment) = valtype(typeof(assign))
@inline Base.valtype(::Type{<:Assignment{<:Action, <:Parameters, <:Function, <:Tuple{Vararg{Symbol}}, R}}) where R = R

"""
    update!(assign::Assignment; kwargs...) -> Assignment

Update the parameters of an assignment and the status of its associated action.
"""
function update!(assign::Assignment; kwargs...)
    if length(kwargs)>0
        assign.parameters = update(assign.parameters; kwargs...)
        update!(assign.action; assign.map(assign.parameters)...)
    end
    return assign
end

"""
    Algorithm{E<:Engine, P<:Parameters, M<:Function} <: Function

An algorithm associated with an engine.
"""
mutable struct Algorithm{E<:Engine, P<:Parameters, M<:Function} <: Function
    name::Symbol
    engine::E
    din::String
    dout::String
    parameters::P
    map::M
    assignments::Dict{Symbol, Assignment}
    timer::TimerOutput
end
@inline run!(alg::Algorithm, assign::Assignment) = nothing
@inline function Base.:(==)(algorithm₁::Algorithm, algorithm₂::Algorithm)
    return ==(
        (algorithm₁.name, algorithm₁.engine, algorithm₁.din, algorithm₁.dout, algorithm₁.parameters, algorithm₁.map, algorithm₁.assignments),
        (algorithm₂.name, algorithm₂.engine, algorithm₂.din, algorithm₂.dout, algorithm₂.parameters, algorithm₂.map, algorithm₂.assignments),
    )
end
@inline function Base.isequal(algorithm₁::Algorithm, algorithm₂::Algorithm)
    return isequal(
        (algorithm₁.name, algorithm₁.engine, algorithm₁.din, algorithm₁.dout, algorithm₁.parameters, algorithm₁.map, algorithm₁.assignments),
        (algorithm₂.name, algorithm₂.engine, algorithm₂.din, algorithm₂.dout, algorithm₂.parameters, algorithm₂.map, algorithm₂.assignments),
    )
end
function Base.show(io::IO, alg::Algorithm)
    @printf io "%s(%s)" alg.name alg.engine
    for (name, value) in pairs(alg.parameters)
        @printf io "_%s" decimaltostr(value, 10)
    end
end

"""
    Algorithm(name::Symbol, engine::Engine; din::String=".", dout::String=".", parameters::Parameters=Parameters(engine), map::Function=identity)

Construct an algorithm.
"""
@inline function Algorithm(name::Symbol, engine::Engine; din::String=".", dout::String=".", parameters::Parameters=Parameters(engine), map::Function=identity)
    return Algorithm(name, engine, din, dout, parameters, map, Dict{Symbol, Assignment}(), TimerOutput())
end

"""
    update!(alg::Algorithm; kwargs...) -> Algorithm

Update the parameters of an algorithm and its associated engine.
"""
function update!(alg::Algorithm; kwargs...)
    if length(kwargs)>0
        alg.parameters = update(alg.parameters; kwargs...)
        update!(alg.engine; alg.map(alg.parameters)...)
    end
    return alg
end

"""
    repr(alg::Algorithm, f::Function=param->true; ndecimal::Int=10) -> String

Get the repr representation of an algorithm.

Optionally, some parameters of the algorithm can be filtered by specifying the `f` function.
Besides, the maximum number of decimals of the parameters can also be specified.
"""
function Base.repr(alg::Algorithm, f::Function=param->true; ndecimal::Int=10)
    result = [@sprintf "%s(%s)" alg.name alg.engine]
    for (name, value) in pairs(alg.parameters)
        f(name) && push!(result, decimaltostr(value, ndecimal))
    end
    return join(result, "_")
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
    nameof(alg::Algorithm, assign::Assignment) -> String

Get the name of the combination of an algorithm and an assignment.
"""
@inline Base.nameof(alg::Algorithm, assign::Assignment) = @sprintf "%s_%s" repr(alg) assign.id

"""
    add!(alg::Algorithm, id::Symbol, action::Action; parameters::Parameters=Parameters{()}(), kwargs...) -> Tuple{Algorithm, Assignment}

Add an assignment on a algorithm by providing the contents of the assignment without the execution of it.
"""
function add!(alg::Algorithm, id::Symbol, action::Action; parameters::Parameters=Parameters{()}(), kwargs...)
    parameters = merge(alg.parameters, parameters)
    map = get(kwargs, :map, identity)
    dependences = get(kwargs, :dependences, ())
    data = prepare!(action, alg.engine)
    save = get(kwargs, :save, false)
    assign = Assignment(id, action, parameters, map, dependences, data, false, save)
    alg.assignments[id] = assign
    return (alg, assign)
end

"""
    (alg::Algorithm)(assign::Assignment) -> Tuple{Algorithm, Assignment}
    (assign::Assignment)(alg::Algorithm) -> Tuple{Algorithm, Assignment}

Run an assignment based on an algorithm.

The difference between these two methods is that the first uses the parameters of `assign` as the current parameters while the second uses those of `alg`.
"""
function (alg::Algorithm)(assign::Assignment)
    ismatched = match(assign.parameters, alg.parameters)
    if !(assign.ismatched && ismatched)
        ismatched || update!(alg; assign.parameters...)
        run!(alg, assign)
        assign.ismatched = true
    end
    return (alg, assign)
end
function (assign::Assignment)(alg::Algorithm)
    ismatched = match(alg.parameters, assign.parameters)
    if !(assign.ismatched && ismatched)
        ismatched || update!(assign; alg.parameters...)
        run!(alg, assign)
        assign.ismatched = true
    end
    return (alg, assign)
end

"""
    (alg::Algorithm)(id::Symbol, action::Action; info::Bool=true, parameters::Parameters=Parameters{()}(), kwargs...) -> Tuple{Algorithm, Assignment}

Add an assignment on a algorithm by providing the contents of the assignment, and run this assignment.
"""
@inline function (alg::Algorithm)(id::Symbol, action::Action; info::Bool=true, parameters::Parameters=Parameters{()}(), kwargs...)
    add!(alg, id, action; parameters=parameters, kwargs...)
    alg(id, info=info)
end

"""
    (alg::Algorithm)(id::Symbol; info::Bool=true, parameters::Parameters=Parameters{()}()) -> Tuple{Algorithm, Assignment}

Run an assignment specified by its id.

Optionally, the time of the run process can be informed by setting the `info` argument to be `true`.
"""
function (alg::Algorithm)(id::Symbol; info::Bool=true, parameters::Parameters=Parameters{()}())
    assign = alg.assignments[id]
    update!(assign; parameters...)
    @timeit alg.timer string(id) alg(assign)
    info && @info "Action $id($(nameof(assign.action|>typeof))): time consumed $(TimerOutputs.time(alg.timer[string(id)]) / 10^9)s."
    assign.save && save(alg, assign)
    return (alg, assign)
end

"""
    save(alg::Algorithm, assign::Assignment) -> Tuple{Algorithm, Assignment}

Save the data of an assignment registered on an algorithm.
"""
function save(alg::Algorithm, assign::Assignment)
    serialize(@sprintf("%s/%s.dat", alg.dout, nameof(alg, assign)), data)
    return (alg, assign)
end

"""
    rundependences!(alg::Algorithm, assign::Assignment, f::Function=assign->true) -> Tuple{Algorithm, Assignment}

Run the dependences of an assignment.

Optionally, some dependences can be filtered by specifying the `f` function.
"""
function rundependences!(alg::Algorithm, assign::Assignment, f::Function=assign->true)
    for id in assign.dependences
        dependence = alg.assignments[id]
        f(dependence) && dependence(alg)
    end
    return (alg, assign)
end

"""
    @recipe plot(pack::Tuple{Algorithm, Assignment})

Define the recipe for the visualization of an assignment of an algorithm.
"""
@recipe function plot(pack::Tuple{Algorithm, Assignment})
    title --> nameof(pack[1], pack[2])
    titlefontsize --> 10
    legend --> false
    seriestype --> (isa(pack[2].data, Tuple{Any, Any, Any}) ? :heatmap : :path)
    pack[2].data
end

end  # module
