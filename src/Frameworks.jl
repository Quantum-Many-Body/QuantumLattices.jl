module Frameworks

using DelimitedFiles: writedlm
using LaTeXStrings: latexstring
using Printf: @printf, @sprintf
using RecipesBase: RecipesBase, @recipe
using Serialization: serialize
using TimerOutputs: TimerOutput, TimerOutputs, @timeit
using ..DegreesOfFreedom: plain, Boundary, Hilbert, Table, Term, ismodulatable
using ..QuantumOperators: Operators, LinearFunction, LinearTransformation, Transformation, identity, optype
using ..Spatials: AbstractLattice, Bond, Neighbors, bonds!, isintracell
using ..Toolkit: atol, efficientoperations, rtol, decimaltostr

import ..QuantumLattices: add!, expand, expand!, id, reset!, update, update!
import ..Spatials: save
import ..Toolkit: contentnames, getcontent

export Action, Algorithm, AnalyticalExpression, Assignment, CompositeGenerator, Entry, Frontend, Image, OperatorGenerator, Parameters, RepresentationGenerator, initialize, prepare!, run!, save

"""
    Parameters{Names}(values::Number...) where Names

A NamedTuple that contain the key-value pairs.
"""
const Parameters{Names} = NamedTuple{Names, <:Tuple{Vararg{Number}}}
@inline Parameters{Names}(values::Number...) where {Names} = NamedTuple{Names}(values)
@inline @generated function update(params::NamedTuple; parameters...)
    names = fieldnames(params)
    values = Expr(:tuple, [:(get(parameters, $name, getfield(params, $name))) for name in QuoteNode.(names)]...)
    return :(NamedTuple{$names}($values))
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
    Parameters(bound::Boundary)

Get the parameters of the twisted boundary condition.
"""
@inline Parameters(bound::Boundary) = NamedTuple{keys(bound)}(ntuple(i->bound.values[i], Val(fieldcount(typeof(keys(bound))))))

"""
    Frontend

The frontend of algorithms applied to a quantum lattice system.
"""
abstract type Frontend end
@inline Base.:(==)(frontend₁::Frontend, frontend₂::Frontend) = ==(efficientoperations, frontend₁, frontend₂)
@inline Base.isequal(frontend₁::Frontend, frontend₂::Frontend) = isequal(efficientoperations, frontend₁, frontend₂)
@inline Base.repr(frontend::Frontend) = String(nameof(typeof(frontend)))
@inline Base.show(io::IO, frontend::Frontend) = @printf io "%s" nameof(typeof(frontend))
@inline Base.valtype(frontend::Frontend) = valtype(typeof(frontend))
@inline update!(frontend::Frontend; kwargs...) = error("update! error: not implemented for $(nameof(typeof(frontend))).")
@inline Parameters(frontend::Frontend) = error("Parameters error: not implemented for $(nameof(typeof(frontend))).")

"""
    RepresentationGenerator <: Frontend

Representation generator of a quantum lattice system.
"""
abstract type RepresentationGenerator <: Frontend end
@inline Base.eltype(gen::RepresentationGenerator) = eltype(typeof(gen))
@inline Base.IteratorSize(::Type{<:RepresentationGenerator}) = Base.SizeUnknown()
@inline function Base.iterate(gen::RepresentationGenerator)
    ops = expand(gen)
    index = iterate(ops)
    isnothing(index) && return nothing
    return index[1], (ops, index[2])
end
@inline function Base.iterate(gen::RepresentationGenerator, state)
    index = iterate(state[1], state[2])
    isnothing(index) && return nothing
    return index[1], (state[1], index[2])
end
@inline Base.show(io::IO, ::MIME"text/latex", gen::RepresentationGenerator) = show(io, MIME"text/latex"(), latexstring(latexstring(expand(gen))))

"""
    expand(gen::RepresentationGenerator) -> valtype(gen)
    expand!(result, gen::RepresentationGenerator) -> typeof(result)

Expand the generator to get the representation of the quantum lattice system (or some part of it).
"""
@inline expand(gen::RepresentationGenerator) = expand!(zero(valtype(gen)), gen)
@inline expand!(result, gen::RepresentationGenerator) = error("expand! error: not implemented for $(nameof(typeof(gen))).")

"""
    AnalyticalExpression{F<:Function, P<:Parameters} <: RepresentationGenerator

Representation of a quantum lattice system by an analytical expression.
"""
mutable struct AnalyticalExpression{F<:Function, P<:Parameters} <: RepresentationGenerator
    const expression::F
    parameters::P
end
@inline function update!(expression::AnalyticalExpression; parameters...)
    expression.parameters = update(expression.parameters; parameters...)
    update!(expression.expression; parameters...)
    return expression
end
@inline update!(expression::Function; parameters...) = expression
@inline Parameters(expression::AnalyticalExpression) = expression.parameters
@inline (expression::AnalyticalExpression)(; kwargs...) = expression.expression(values(expression.parameters)...; kwargs...)

"""
    Entry{C, A<:NamedTuple, B<:NamedTuple, P<:Parameters, D<:Boundary} <: RepresentationGenerator

The basic representation generator of a quantum lattice system that records the quantum operators or a representation of the quantum operators related to (part of) the system.
"""
mutable struct Entry{C, A<:NamedTuple, B<:NamedTuple, P<:Parameters, D<:Boundary} <: RepresentationGenerator
    const constops::C
    const alterops::A
    const boundops::B
    parameters::P
    const boundary::D
end
@inline Entry(entry::Entry) = entry
@inline Base.eltype(E::Type{<:Entry}) = eltype(valtype(E))

"""
    Entry(terms::Tuple{Vararg{Term}}, bonds::Vector{<:Bond}, hilbert::Hilbert; half::Bool=false, boundary::Boundary=plain)

Construct an entry of quantum operators based on the input terms, bonds, Hilbert space and (twisted) boundary condition.
"""
function Entry(terms::Tuple{Vararg{Term}}, bonds::Vector{<:Bond}, hilbert::Hilbert; half::Bool=false, boundary::Boundary=plain)
    constterms, alterterms, choosedterms = termclassifier(terms)
    innerbonds = filter(bond->isintracell(bond), bonds)
    boundbonds = filter(bond->!isintracell(bond), bonds)
    constops = Operators{mapreduce(term->optype(typeof(term), typeof(hilbert), eltype(bonds)), promote_type, choosedterms)}()
    map(term->expand!(constops, term, innerbonds, hilbert; half=half), constterms)
    alterops = NamedTuple{map(id, alterterms)}(map(term->expand(one(term), innerbonds, hilbert; half=half), alterterms))
    boundops = NamedTuple{map(id, terms)}(map(term->map!(boundary, expand!(Operators{valtype(typeof(boundary), optype(typeof(term), typeof(hilbert), eltype(bonds)))}(), one(term), boundbonds, hilbert, half=half)), terms))
    parameters = NamedTuple{map(id, terms)}(map(term->term.value, terms))
    return Entry(constops, alterops, boundops, parameters, boundary)
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
    *(entry::Entry, factor) -> Entry
    *(factor, entry::Entry) -> Entry

Multiply an entry of quantum operators with a factor.
"""
@inline Base.:*(entry::Entry, factor) = factor * entry
@inline function Base.:*(factor, entry::Entry)
    parameters = NamedTuple{keys(entry.parameters)}(map(value->factor*value, values(entry.parameters)))
    return Entry(factor*entry.constops, entry.alterops, entry.boundops, parameters, entry.boundary)
end

"""
    +(entry₁::Entry, entry₂::Entry) -> Entry

Addition of two entries of quantum operators.
"""
function Base.:+(entry₁::Entry, entry₂::Entry)
    @assert entry₁.boundary==entry₂.boundary "+ error: in order to be added, two entries must share the same boundary condition (including the twist angles at the boundary)."
    constops = entry₁.constops + entry₂.constops
    alls, allshares = totalkeys(entry₁.parameters, entry₂.parameters), sharedkeys(entry₁.parameters, entry₂.parameters)
    allmatches = NamedTuple{keymaps(allshares)}(map(key->opsmatch(entry₁.alterops, entry₂.alterops, key) && opsmatch(entry₁.boundops, entry₂.boundops, key), allshares))
    parameters = NamedTuple{keymaps(alls)}(map(key->combinevalue(entry₁.parameters, entry₂.parameters, allmatches, key), alls))
    alteralls, altershares = totalkeys(entry₁.alterops, entry₂.alterops), sharedkeys(entry₁.alterops, entry₂.alterops)
    boundalls, boundshares = totalkeys(entry₁.boundops, entry₂.boundops), sharedkeys(entry₁.boundops, entry₂.boundops)
    altermatches = NamedTuple{keymaps(altershares)}(map(((::Val{key}) where key)->getfield(allmatches, key), altershares))
    boundmatches = NamedTuple{keymaps(boundshares)}(map(((::Val{key}) where key)->getfield(allmatches, key), boundshares))
    alterops = NamedTuple{keymaps(alteralls)}(map(key->combineops(entry₁.alterops, entry₁.parameters, entry₂.alterops, entry₂.parameters, altermatches, key), alteralls))
    boundops = NamedTuple{keymaps(boundalls)}(map(key->combineops(entry₁.boundops, entry₁.parameters, entry₂.boundops, entry₂.parameters, boundmatches, key), boundalls))
    return Entry(constops, alterops, boundops, parameters, deepcopy(entry₁.boundary))
end
@inline keymaps(keys) = map(((::Val{key}) where key)->key, keys)
@generated totalkeys(content₁::NamedTuple, content₂::NamedTuple) = map(Val, Tuple(unique((fieldnames(content₁)..., fieldnames(content₂)...))))
@generated sharedkeys(content₁::NamedTuple, content₂::NamedTuple) = map(Val, Tuple(intersect(fieldnames(content₁), fieldnames(content₂))))
@inline opsmatch(ops₁::NamedTuple, ops₂::NamedTuple, ::Val{key}) where key = opsmatch(get(ops₁, key, nothing), get(ops₂, key, nothing))
@inline opsmatch(ops₁, ops₂) = ops₁==ops₂ || zero(ops₁)==ops₁ || zero(ops₂)==ops₂
@inline opsmatch(ops, ::Nothing) = true
@inline opsmatch(::Nothing, ops) = true
@inline opsmatch(::Nothing, ::Nothing) = true
@inline combinevalue(params₁::Parameters, params₂::Parameters, matches::NamedTuple, ::Val{key}) where key = combinevalue(get(params₁, key, nothing), get(params₂, key, nothing), get(matches, key, nothing))
@inline combinevalue(value₁, value₂, flag::Bool) = flag ? value₁+value₂ : value₁==value₂ ? value₁ : promote(one(value₁), one(value₂))[1]
@inline combinevalue(value, ::Nothing, ::Nothing) = value
@inline combinevalue(::Nothing, value, ::Nothing) = value
@inline function combineops(ops₁::NamedTuple, params₁::Parameters, ops₂::NamedTuple, params₂::Parameters, matches::NamedTuple, ::Val{key}) where key
    combineops(get(ops₁, key, nothing), get(params₁, key, nothing), get(ops₂, key, nothing), get(params₂, key, nothing), get(matches, key, nothing))
end
@inline combineops(ops₁, value₁, ops₂, value₂, flag::Bool) = flag ? deepcopy(ops₁) : value₁==value₂ ? ops₁+ops₂ : ops₁*value₁+ops₂*value₂
@inline combineops(ops, value, ::Nothing, ::Nothing, ::Nothing) = deepcopy(ops)
@inline combineops(::Nothing, ::Nothing, ops, value, ::Nothing) = deepcopy(ops)

"""
    (transformation::LinearTransformation)(entry::Entry; kwargs...) -> Entry

Apply a linear transformation to an entry of (representations of) quantum operators.
"""
function (transformation::LinearTransformation)(entry::Entry; kwargs...)
    wrapper(m) = transformation(m; kwargs...)
    constops = wrapper(entry.constops)
    alterops = NamedTuple{keys(entry.alterops)}(map(wrapper, values(entry.alterops)))
    boundops = NamedTuple{keys(entry.boundops)}(map(wrapper, values(entry.boundops)))
    return Entry(constops, alterops, boundops, entry.parameters, deepcopy(entry.boundary))
end

"""
    valtype(::Type{<:Entry})

Get the valtype of an entry of (representations of) quantum operators.
"""
@inline @generated function Base.valtype(::Type{<:Entry{C, A, B}}) where {C, A<:NamedTuple, B<:NamedTuple}
    optp = C
    (fieldcount(A) > 0) && (optp = reduce(promote_type, fieldtypes(A), init=optp))
    (fieldcount(B) > 0) && (optp = reduce(promote_type, fieldtypes(B), init=optp))
    return optp
end

"""
    expand!(result, entry::Entry) -> typeof(result)

Expand an entry to get the (representation of) quantum operators related to a quantum lattice system.
"""
@generated function expand!(result, entry::Entry)
    exprs = [:(add!(result, entry.constops))]
    for name in QuoteNode.(fieldnames(fieldtype(entry, :alterops)))
        push!(exprs, :(add!(result, get(entry.parameters, $name, 1)*getfield(entry.alterops, $name))))
    end
    for name in QuoteNode.(fieldnames(fieldtype(entry, :boundops)))
        push!(exprs, :(add!(result, get(entry.parameters, $name, 1)*getfield(entry.boundops, $name))))
    end
    push!(exprs, :(return result))
    return Expr(:block, exprs...)
end

"""
    update!(entry::Entry{<:Operators}; parameters...) -> Entry

Update the parameters (including the boundary parameters) of an entry of quantum operators.

!!! Note
    The coefficients of `boundops` are also updated due to the change of the boundary parameters.
"""
function update!(entry::Entry{<:Operators}; parameters...)
    entry.parameters = update(entry.parameters; parameters...)
    if !match(Parameters(entry.boundary), NamedTuple{keys(parameters)}(values(parameters)))
        old = copy(entry.boundary.values)
        update!(entry.boundary; parameters...)
        map(ops->map!(LinearFunction(op->entry.boundary(op, origin=old)), ops), values(entry.boundops))
    end
    return entry
end

"""
    update!(entry::Entry, transformation::LinearTransformation, source::Entry{<:Operators}; kwargs...) -> Entry

Update the parameters (including the boundary parameters) of an entry based on its source entry of quantum operators and the corresponding linear transformation.

!!! Note
    The coefficients of `boundops` are also updated due to the change of the boundary parameters.
"""
function update!(entry::Entry, transformation::LinearTransformation, source::Entry{<:Operators}; kwargs...)
    entry.parameters = update(entry.parameters; source.parameters...)
    if !match(Parameters(entry.boundary), Parameters(source.boundary))
        update!(entry.boundary; Parameters(source.boundary)...)
        map((dest, ops)->map!(LinearFunction(op->transformation(op; kwargs...)), empty!(dest), ops), values(entry.boundops), values(source.boundops))
    end
    return entry
end

"""
    empty(entry::Entry) -> Entry
    empty!(entry::Entry) -> Entry

Get an empty copy of an entry or empty an entry of (representations of) quantum operators.
"""
@inline function Base.empty(entry::Entry)
    constops = empty(entry.constops)
    alterops = NamedTuple{keys(entry.alterops)}(map(empty, values(entry.alterops)))
    boundops = NamedTuple{keys(entry.boundops)}(map(empty, values(entry.boundops)))
    return Entry(constops, alterops, boundops, entry.parameters, deepcopy(entry.boundary))
end
@inline function Base.empty!(entry::Entry)
    empty!(entry.constops)
    map(empty!, values(entry.alterops))
    map(empty!, values(entry.boundops))
    return entry
end

"""
    reset!(entry::Entry{<:Operators}, terms::Tuple{Vararg{Term}}, bonds::Vector{<:Bond}, hilbert::Hilbert; half::Bool=false, boundary::Boundary=entry.boundary) -> Entry

Reset an entry of quantum operators by the new terms, bonds, Hilbert space and (twisted) boundary condition.
"""
function reset!(entry::Entry{<:Operators}, terms::Tuple{Vararg{Term}}, bonds::Vector{<:Bond}, hilbert::Hilbert; half::Bool=false, boundary::Boundary=entry.boundary)
    empty!(entry)
    constterms, alterterms, _ = termclassifier(terms)
    innerbonds = filter(bond->isintracell(bond), bonds)
    boundbonds = filter(bond->!isintracell(bond), bonds)
    map(term->expand!(entry.constops, term, innerbonds, hilbert; half=half), constterms)
    map(term->expand!(getfield(entry.alterops, id(term)), one(term), innerbonds, hilbert; half=half), alterterms)
    map(term->map!(boundary, expand!(getfield(entry.boundops, id(term)), one(term), boundbonds, hilbert; half=half)), terms)
    entry.parameters = NamedTuple{map(id, terms)}(map(term->term.value, terms))
    merge!(entry.boundary, boundary)
    return entry
end

"""
    reset!(entry::Entry, transformation::LinearTransformation, source::Entry{<:Operators}; kwargs...)

Reset an entry by its source entry of quantum operators and the corresponding linear transformation.
"""
function reset!(entry::Entry, transformation::LinearTransformation, source::Entry{<:Operators}; kwargs...)
    wrapper = LinearFunction(op->transformation(op; kwargs...))
    map!(wrapper, empty!(entry.constops), source.constops)
    map((dest, ops)->map!(wrapper, empty!(dest), ops), values(entry.alterops), values(source.alterops))
    map((dest, ops)->map!(wrapper, empty!(dest), ops), values(entry.boundops), values(source.boundops))
    entry.parameters = update(entry.parameters; source.parameters...)
    merge!(entry.boundary, source.boundary)
    return entry
end

"""
    Parameters(entry::Entry)

Get the complete set of parameters of an entry of (representations of) quantum operators.
"""
@inline Parameters(entry::Entry) = merge(entry.parameters, Parameters(entry.boundary))

"""
    CompositeGenerator{E<:Entry, T<:Union{Table, Nothing}} <: RepresentationGenerator

Abstract type for a composite representation generator of a quantum lattice system.

By protocol, it must have the following predefined contents:
* `operators::E`: the entry for the generated (representations of) quantum operators
* `table::T`: the index-sequence table if it is not nothing
"""
abstract type CompositeGenerator{E<:Entry, T<:Union{Table, Nothing}} <: RepresentationGenerator end
@inline contentnames(::Type{<:CompositeGenerator}) = (:operators, :table)
@inline Base.valtype(::Type{<:CompositeGenerator{E}}) where {E<:Entry} = valtype(E)
@inline Base.eltype(::Type{<:CompositeGenerator{E}}) where {E<:Entry} = eltype(E)
@inline expand!(result, gen::CompositeGenerator) = expand!(result, getcontent(gen, :operators))
@inline Parameters(gen::CompositeGenerator) = Parameters(getcontent(gen, :operators))
@inline Entry(gen::CompositeGenerator) = getcontent(gen, :operators)

"""
    OperatorGenerator{E<:Entry{<:Operators}, TS<:Tuple{Vararg{Term}}, B<:Bond, H<:Hilbert, T<:Union{Table, Nothing}} <: CompositeGenerator{E, T}

A generator of operators based on the terms, bonds and Hilbert space of a quantum lattice system.
"""
struct OperatorGenerator{E<:Entry{<:Operators}, TS<:Tuple{Vararg{Term}}, B<:Bond, H<:Hilbert, T<:Union{Table, Nothing}} <: CompositeGenerator{E, T}
    operators::E
    terms::TS
    bonds::Vector{B}
    hilbert::H
    half::Bool
    table::T
end
@inline contentnames(::Type{<:OperatorGenerator}) = (:operators, :terms, :bonds, :hilbert, :half, :table)

"""
    OperatorGenerator(terms::Tuple{Vararg{Term}}, bonds::Vector{<:Bond}, hilbert::Hilbert; half::Bool=false, boundary::Boundary=plain, table::Union{Table,Nothing}=nothing)

Construct a generator of operators.
"""
@inline function OperatorGenerator(terms::Tuple{Vararg{Term}}, bonds::Vector{<:Bond}, hilbert::Hilbert; half::Bool=false, boundary::Boundary=plain, table::Union{Table,Nothing}=nothing)
    return OperatorGenerator(Entry(terms, bonds, hilbert; half=half, boundary=boundary), terms, bonds, hilbert, half, table)
end

"""
    update!(gen::OperatorGenerator; parameters...) -> typeof(gen)

Update the coefficients of the terms in a generator.
"""
@inline function update!(gen::OperatorGenerator; parameters...)
    update!(gen.operators; parameters...)
    map(term->(ismodulatable(term) ? update!(term; parameters...) : term), gen.terms)
    return gen
end

"""
    empty(gen::OperatorGenerator) -> OperatorGenerator
    empty!(gen::OperatorGenerator) -> OperatorGenerator

Get an empty copy of or empty an operator generator.
"""
@inline function Base.empty(gen::OperatorGenerator)
    return OperatorGenerator(empty(gen.operators), gen.terms, empty(gen.bonds), empty(gen.hilbert), gen.half, isnothing(gen.table) ? nothing : empty(gen.table))
end
function Base.empty!(gen::OperatorGenerator)
    empty!(gen.bonds)
    empty!(gen.hilbert)
    isnothing(gen.table) || empty!(gen.table)
    empty!(gen.operators)
    return gen
end

"""
    reset!(gen::OperatorGenerator, lattice::AbstractLattice, hilbert::Hilbert; neighbors=max(map(term->term.bondkind, gen.terms)...)) -> OperatorGenerator

Reset an operator generator by a new lattice and the corresponding new hilbert space.
"""
function reset!(gen::OperatorGenerator, lattice::AbstractLattice, hilbert::Hilbert; neighbors=max(map(term->term.bondkind, gen.terms)...))
    isa(neighbors, Neighbors) || (neighbors = Neighbors(lattice, neighbors))
    bonds!(empty!(gen.bonds), lattice, neighbors)
    merge!(empty!(gen.hilbert), hilbert)
    isnothing(gen.table) || reset!(gen.table, gen.hilbert)
    reset!(gen.operators, gen.terms, gen.bonds, gen.hilbert; half=gen.half, boundary=replace(gen.operators.boundary; vectors=lattice.vectors))
    return gen
end

"""
    expand(gen::OperatorGenerator, name::Symbol) -> Operators
    expand(gen::OperatorGenerator, i::Int) -> Operators
    expand(gen::OperatorGenerator, name::Symbol, i::Int) -> Operators

Expand an operator generator to get:
1) the operators of a specific term;
2) the operators on a specific bond;
3) the operators of a specific term on a specific bond.
"""
function expand(gen::OperatorGenerator, name::Symbol)
    result = zero(valtype(gen))
    term = get(gen.terms, Val(name))
    if !ismodulatable(term)
        expand!(result, term, filter(bond->isintracell(bond), gen.bonds), gen.hilbert; half=gen.half)
    else
        for opt in getfield(gen.operators.alterops, name)
            add!(result, opt*term.value)
        end
    end
    for opt in getfield(gen.operators.boundops, name)
        add!(result, opt*term.value)
    end
    return result
end
function expand(gen::OperatorGenerator, i::Int)
    bond = gen.bonds[i]
    result = zero(valtype(gen))
    map(term->expand!(result, term, bond, gen.hilbert; half=gen.half), gen.terms)
    isintracell(bond) || for opt in result
        result[id(opt)] = gen.operators.boundary(opt)
    end
    return result
end
function expand(gen::OperatorGenerator, name::Symbol, i::Int)
    bond = gen.bonds[i]
    term = get(gen.terms, Val(name))
    result = expand!(zero(valtype(gen)), term, bond, gen.hilbert; half=gen.half)
    isintracell(bond) || for opt in result
        result[id(opt)] = gen.operators.boundary(opt)
    end
    return result
end
@inline @generated function Base.get(terms::Tuple{Vararg{Term}}, ::Val{Name}) where Name
    i = findfirst(isequal(Name), map(id, fieldtypes(terms)))::Int
    return :(terms[$i])
end

"""
    Image{E<:Entry, H<:Transformation, T<:Union{Table, Nothing}} <: CompositeGenerator{E, T}

The image of a transformation applied to a representation of a quantum lattice system.
"""
mutable struct Image{E<:Entry, H<:Transformation, T<:Union{Table, Nothing}} <: CompositeGenerator{E, T}
    const operators::E
    transformation::H
    const table::T
    const sourceid::UInt
end
@inline contentnames(::Type{<:Image}) = (:operators, :transformation, :table, :sourceid)

"""
    (transformation::Transformation)(gen::RepresentationGenerator; table::Union{Table, Nothing}=nothing, kwargs...) -> Image

Get the image of a transformation applied to a representation of a quantum lattice system.
"""
@inline function (transformation::Transformation)(gen::RepresentationGenerator; table::Union{Table, Nothing}=nothing, kwargs...)
    return Image(transformation(Entry(gen); kwargs...), transformation, table, objectid(gen))
end

"""
    empty(gen::Image) -> Image
    empty!(gen::Image) -> Image

Get an empty copy of or empty the image of a transformation applied to a representation.
"""
@inline Base.empty(gen::Image) = Image(empty(gen.operators), gen.transformation, isnothing(gen.table) ? nothing : empty(gen.table), gen.sourceid)
@inline function Base.empty!(gen::Image)
    empty!(gen.operators)
    isnothing(gen.table) || empty!(gen.table)
    return gen
end

"""
    update!(gen::Image; parameters...) -> typeof(gen)

Update the parameters of the image of a transformation applied to a representation.
"""
@inline function update!(gen::Image; parameters...)
    update!(gen.operators; parameters...)
    return gen
end

"""
    update!(gen::Image, source::RepresentationGenerator; kwargs...) -> Image

Update the parameters of the image based on its source representation.
"""
@inline function update!(gen::Image, source::RepresentationGenerator; kwargs...)
    @assert gen.sourceid==objectid(source) "update! error: mismatched image, transformation and source representation."
    update!(gen.operators, gen.transformation, Entry(source); kwargs...)
    return gen
end

"""
    reset!(gen::Image, transformation::Transformation, source::RepresentationGenerator; table::Union{Table, Nothing}=nothing, kwargs...) -> Image

Reset the image of a transformation applied to a representation.
"""
@inline function reset!(gen::Image, transformation::Transformation, source::RepresentationGenerator; table::Union{Table, Nothing}=nothing, kwargs...)
    @assert gen.sourceid==objectid(source) "reset! error: mismatched image, transformation and source representation."
    reset!(gen.operators, transformation, Entry(source); kwargs...)
    gen.transformation = transformation
    isnothing(table) && (table = isa(source, CompositeGenerator) ? getcontent(source, :table) : gen.table)
    isnothing(table) || gen.table===table || merge!(empty!(gen.table), table)
    return gen
end

"""
    Action

Abstract type for all actions.
"""
abstract type Action end
@inline Base.:(==)(action₁::Action, action₂::Action) = ==(efficientoperations, action₁, action₂)
@inline Base.isequal(action₁::Action, action₂::Action) = isequal(efficientoperations, action₁, action₂)
@inline initialize(action::Action, frontend::Frontend) = error("initialize error: not implemented.")
@inline update!(action::Action; parameters...) = action

"""
    Assignment{A<:Action, P<:Parameters, M<:Function, N, D} <: Function

An assignment associated with an action.
"""
mutable struct Assignment{A<:Action, P<:Parameters, M<:Function, N, D} <: Function
    const id::Symbol
    const action::A
    parameters::P
    const map::M
    const dependences::NTuple{N, Symbol}
    data::D
    ismatched::Bool
end
@inline Base.:(==)(assign₁::Assignment, assign₂::Assignment) = ==(efficientoperations, assign₁, assign₂)
@inline Base.isequal(assign₁::Assignment, assign₂::Assignment) = isequal(efficientoperations, assign₁, assign₂)
@inline Parameters(assignment::Assignment) = assignment.parameters

"""
    valtype(assign::Assignment)
    valtype(::Type{<:Assignment})

The type of the data(result) of an assignment.
"""
@inline Base.valtype(assign::Assignment) = valtype(typeof(assign))
@inline Base.valtype(::Type{<:Assignment{<:Action, <:Parameters, <:Function, N, R} where N}) where R = R

"""
    update!(assign::Assignment; parameters...) -> Assignment

Update the parameters of an assignment and the status of its associated action.
"""
function update!(assign::Assignment; parameters...)
    if length(parameters)>0
        assign.parameters = update(assign.parameters; parameters...)
        update!(assign.action; assign.map(assign.parameters)...)
    end
    return assign
end

"""
    Algorithm{F<:Frontend, P<:Parameters, M<:Function} <: Function

An algorithm associated with an frontend.
"""
mutable struct Algorithm{F<:Frontend, P<:Parameters, M<:Function} <: Function
    const name::Symbol
    const frontend::F
    const din::String
    const dout::String
    parameters::P
    const map::M
    const assignments::Dict{Symbol, Assignment}
    const timer::TimerOutput
end
@inline run!(alg::Algorithm, assign::Assignment) = nothing
@inline function Base.:(==)(algorithm₁::Algorithm, algorithm₂::Algorithm)
    return ==(
        (algorithm₁.name, algorithm₁.frontend, algorithm₁.din, algorithm₁.dout, algorithm₁.parameters, algorithm₁.map, algorithm₁.assignments),
        (algorithm₂.name, algorithm₂.frontend, algorithm₂.din, algorithm₂.dout, algorithm₂.parameters, algorithm₂.map, algorithm₂.assignments),
    )
end
@inline function Base.isequal(algorithm₁::Algorithm, algorithm₂::Algorithm)
    return isequal(
        (algorithm₁.name, algorithm₁.frontend, algorithm₁.din, algorithm₁.dout, algorithm₁.parameters, algorithm₁.map, algorithm₁.assignments),
        (algorithm₂.name, algorithm₂.frontend, algorithm₂.din, algorithm₂.dout, algorithm₂.parameters, algorithm₂.map, algorithm₂.assignments),
    )
end
function Base.show(io::IO, alg::Algorithm)
    @printf io "%s(%s)" alg.name alg.frontend
    for (name, value) in pairs(alg.parameters)
        @printf io "_%s" decimaltostr(value, 10)
    end
end
@inline Parameters(algorithm::Algorithm) = algorithm.parameters

"""
    Algorithm(name::Symbol, frontend::Frontend; din::String=".", dout::String=".", parameters::Parameters=Parameters(frontend), map::Function=identity)

Construct an algorithm.
"""
@inline function Algorithm(name::Symbol, frontend::Frontend; din::String=".", dout::String=".", parameters::Parameters=Parameters(frontend), map::Function=identity)
    return Algorithm(name, frontend, din, dout, parameters, map, Dict{Symbol, Assignment}(), TimerOutput())
end

"""
    update!(alg::Algorithm; parameters...) -> Algorithm

Update the parameters of an algorithm and its associated frontend.
"""
function update!(alg::Algorithm; parameters...)
    if length(parameters)>0
        alg.parameters = update(alg.parameters; parameters...)
        update!(alg.frontend; alg.map(alg.parameters)...)
    end
    return alg
end

"""
    repr(alg::Algorithm, f::Function=param->true; ndecimal::Int=10) -> String

Get the repr representation of an algorithm.

Optionally, some parameters of the algorithm can be filtered by specifying the `f` function. Besides, the maximum number of decimals of the parameters can also be specified by the keyword argument `ndecimal`.
"""
function Base.repr(alg::Algorithm, f::Function=param->true; ndecimal::Int=10)
    result = String[]
    for (name, value) in pairs(alg.parameters)
        f(name) && push!(result, @sprintf "%s(%s)" name decimaltostr(value, ndecimal))
    end
    return @sprintf "%s(%s)-%s" alg.name alg.frontend join(result, "")
end

"""
    summary(alg::Algorithm)

Provide a summary of an algorithm.
"""
function Base.summary(alg::Algorithm)
    @info "Summary of $(alg.name)($(nameof(typeof(alg.frontend)))):"
    @info string(alg.timer)
end

"""
    nameof(alg::Algorithm, assign::Assignment) -> String

Get the name of the combination of an algorithm and an assignment.
"""
@inline Base.nameof(alg::Algorithm, assign::Assignment) = @sprintf "%s-%s" repr(alg) assign.id

"""
    add!(alg::Algorithm, id::Symbol, action::Action; parameters::Parameters=Parameters{()}(), map::Function=identity, dependences::Tuple=(), kwargs...) -> Tuple{Algorithm, Assignment}

Add an assignment on an algorithm by providing the contents of the assignment without the execution of it.
"""
function add!(alg::Algorithm, id::Symbol, action::Action; parameters::Parameters=Parameters{()}(), map::Function=identity, dependences::Tuple=(), kwargs...)
    assign = Assignment(id, action, merge(alg.parameters, parameters), map, dependences, initialize(action, alg.frontend), false)
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
    (alg::Algorithm)(id::Symbol, action::Action; info::Bool=true, kwargs...) -> Tuple{Algorithm, Assignment}

Add an assignment on a algorithm by providing the contents of the assignment, and run this assignment.
"""
@inline function (alg::Algorithm)(id::Symbol, action::Action; info::Bool=true, kwargs...)
    add!(alg, id, action; kwargs...)
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
    return (alg, assign)
end

"""
    save(alg::Algorithm, assign::Assignment; delimited=false) -> Tuple{Algorithm, Assignment}

Save the data of an assignment registered on an algorithm.
"""
@inline function save(alg::Algorithm, assign::Assignment; delimited=false)
    filename = @sprintf("%s/%s.dat", alg.dout, nameof(alg, assign))
    delimited ? save(filename, assign.data) : serialize(filename, assign.data)
    return (alg, assign)
end
@inline save(filename::AbstractString, data::Tuple) = save(filename, data...)
function save(filename::AbstractString, x::AbstractVector{<:Number}, y::Union{AbstractVector{<:Number}, AbstractMatrix{<:Number}})
    @assert length(x)==size(y)[1] "save error: mismatched size of x and y."
    open(filename, "w") do f
        writedlm(f, [x y])
    end
end
function save(filename::AbstractString, x::AbstractVector{<:Number}, y::AbstractVector{<:Number}, z::Union{AbstractMatrix{<:Number}, AbstractArray{<:Number, 3}})
    @assert size(z)[1:2]==(length(y), length(x)) "save error: mismatched size of x, y and z."
    open(filename, "w") do f
        new_x = kron(x, ones(length(y)))
        new_y = kron(ones(length(x)), y)
        new_z = reshape(z, length(x)*length(y), :)
        writedlm(f, [new_x new_y new_z])
    end
end

"""
    prepare!(alg::Algorithm, assign::Assignment, f::Function=assign->true) -> Tuple{Algorithm, Assignment}

Run the dependences of an assignment.

Optionally, some dependences can be filtered by specifying the `f` function.
"""
function prepare!(alg::Algorithm, assign::Assignment, f::Function=assign->true)
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
    title --> nameof(pack...)
    titlefontsize --> 10
    attr = seriestype(pack)
    isnothing(attr) || begin
        seriestype --> attr
        attr==:path && begin
            legend --> false
            minorgrid --> true
            xminorticks --> 10
            yminorticks --> 10
        end
    end
    pack[2].data
end
@inline seriestype(_...) = nothing
@inline seriestype(pack::Tuple{Algorithm, Assignment}) = seriestype(pack[2].data...)
@inline seriestype(::AbstractVector{<:Number}, ::Union{AbstractVector{<:Number}, AbstractMatrix{<:Number}}, _...) = :path
@inline seriestype(::AbstractVector{<:Number}, ::AbstractVector{<:Number}, ::Union{AbstractMatrix{<:Number}, AbstractArray{<:Number, 3}}, _...) = :heatmap

end  # module
