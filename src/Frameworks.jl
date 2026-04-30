module Frameworks

using Base: @propagate_inbounds
using Base.Iterators: flatten, repeated
using IndentWrappers: indent
using JLD2: jldopen, loadtodict!
using Latexify: latexify
using Serialization: deserialize, serialize
using TimerOutputs: @timeit, TimerOutput, time
using ..DegreesOfFreedom: Boundary, Hilbert, Term, plain
using ..QuantumLattices: OneOrMore, ZeroAtLeast, ZeroOrMore, value
using ..QuantumOperators: LinearTransformation, OperatorPack, OperatorSet, OperatorSum, Operators, identity, operatortype
using ..Spatials: Bond, isintracell
using ..Toolkit: atol, efficientoperations, parametertype, rtol

import ..QuantumLattices: add!, expand, expand!, id, reset!, str, update, update!
import ..QuantumOperators: scalartype
import ..Spatials: dlmsave

export Action, Algorithm, Assignment, CategorizedGenerator, Data, Eager, ExpansionStyle, Formula, Frontend, Generator, LatticeModel, Lazy, OperatorGenerator, Parameters, ParametricGenerator, StaticGenerator
export checkoptions, datatype, eager, fingerprint, hasoption, lazy, options, optionsinfo, qldload, qldsave, run!, seriestype

"""
    Parameters{Names}(values::Number...) where Names

A NamedTuple that contains the key-value pairs.
"""
const Parameters{Names, T<:ZeroAtLeast{Number}} = NamedTuple{Names, T}
@inline Parameters() = NamedTuple()
@inline Parameters{Names}(values::Number...) where {Names} = NamedTuple{Names}(values)
function Base.show(io::IO, params::Parameters)
    haskey(io, :ndecimal) && (params = NamedTuple{keys(params)}(map(value->round(value; digits=io[:ndecimal]), values(params))))
    invoke(show, Tuple{IO, NamedTuple}, io, params)
end

"""
    update(params::NamedTuple; parameters...) -> Parameters

Update a set of `Parameters` and return the updated one.
"""
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
    str(params::Parameters; ndecimal::Int=10, select::Function=name::Symbol->true, front::String="", rear::String="") -> String

Convert a set of `Parameters` to a string with each number rounded to at most `ndecimal` decimal places. The `select` function can be used to filter which key-value pairs to include.
"""
function str(params::Parameters; ndecimal::Int=10, select::Function=name::Symbol->true, front::String="", rear::String="")
    result = String[]
    for (name, value) in pairs(params)
        if select(name)
            push!(result, string(name, "(", str(value; ndecimal=ndecimal), ")"))
        end
    end
    return string(append(front, "-"), join(result), prepend(rear, "-"))
end
@inline append(s::String, suffix::String) = isempty(s) ? s : string(s, suffix)
@inline prepend(s::String, prefix::String) = isempty(s) ? s : string(prefix, s)

"""
    Parameters(bound::Boundary)

Get the parameters of the twisted boundary condition.
"""
@inline Parameters(bound::Boundary) = NamedTuple{keys(bound)}(ntuple(i->bound.values[i], Val(fieldcount(typeof(keys(bound))))))

"""
    LatticeModel

Abstract supertype for all representations of a quantum lattice system.

Subtypes must implement `valtype`. `Parameters` and `update!` should also be implemented as applicable.
"""
abstract type LatticeModel end
@inline Base.:(==)(model₁::LatticeModel, model₂::LatticeModel) = ==(efficientoperations, model₁, model₂)
@inline Base.isequal(model₁::LatticeModel, model₂::LatticeModel) = isequal(efficientoperations, model₁, model₂)

"""
    valtype(model::LatticeModel)
    valtype(::Type{<:LatticeModel})

Get the valtype of a `LatticeModel`. Subtypes must implement the type-level method.
"""
@inline Base.valtype(model::LatticeModel) = valtype(typeof(model))

"""
    scalartype(model::LatticeModel)
    scalartype(::Type{T}) where {T<:LatticeModel}

Get the scalar type of a `LatticeModel`.
"""
@inline scalartype(model::LatticeModel) = scalartype(typeof(model))
@inline scalartype(::Type{T}) where {T<:LatticeModel} = scalartype(valtype(T))

"""
    eltype(model::LatticeModel)
    eltype(::Type{T}) where {T<:LatticeModel}

Get the eltype of a `LatticeModel`.
"""
@inline Base.eltype(model::LatticeModel) = eltype(typeof(model))
@inline Base.eltype(::Type{T}) where {T<:LatticeModel} = eltype(valtype(T))

"""
    Formula{V, F<:Function, P<:Parameters}

Representation of a quantum lattice system with an explicit analytical formula.
"""
mutable struct Formula{V, F<:Function, P<:Parameters} <: LatticeModel
    const expression::F
    parameters::P
    function Formula(expression::Function, parameters::Parameters)
        V = Core.Compiler.return_type(expression, parametertype(typeof(parameters), 2))
        @assert isconcretetype(V) "Formula error: input expression is not type-stable."
        new{V, typeof(expression), typeof(parameters)}(expression, parameters)
    end
    function Formula{V}(expression::Function, parameters::Parameters) where V
        new{V, typeof(expression), typeof(parameters)}(expression, parameters)
    end
end
@inline Base.valtype(::Type{<:Formula{V}}) where V = V

"""
    update!(formula::Formula; parameters...) -> Formula

Update the parameters of a `Formula` in place and return itself after update.
"""
@inline function update!(formula::Formula; parameters...)
    formula.parameters = update(formula.parameters; parameters...)
    update!(formula.expression; parameters...)
    return formula
end
@inline update!(expression::Function; parameters...) = expression

"""
    Parameters(formula::Formula) -> Parameters

Get the parameters of a `Formula`.
"""
@inline Parameters(formula::Formula) = formula.parameters

"""
    (formula::Formula)(args...; kwargs...) -> valtype(formula)

Get the result of a `Formula`.
"""
@inline @generated function (formula::Formula)(args...; kwargs...)
    exprs = [:(getfield(formula.parameters, $i)) for i = 1:fieldcount(fieldtype(formula, :parameters))]
    return :(formula.expression($(exprs...), args...; kwargs...))
end

"""
    ExpansionStyle

Expansion style of a generator of (representations of) quantum operators. It has two singleton subtypes, [`Eager`](@ref) and [`Lazy`](@ref).
"""
abstract type ExpansionStyle end

"""
    Eager <: ExpansionStyle

Eager expansion style with eager computation so that similar terms are combined in the final result.
"""
struct Eager <: ExpansionStyle end

"""
    const eager = Eager()

Singleton instance of [`Eager`](@ref).
"""
const eager = Eager()

"""
    Lazy <: ExpansionStyle

Lazy expansion style with lazy computation so that similar terms are not combined in the final result.
"""
struct Lazy <: ExpansionStyle end

"""
    const lazy = Lazy()

Singleton instance of [`Lazy`](@ref).
"""
const lazy = Lazy()

"""
    Generator{V} <: LatticeModel

Abstract supertype for lattice models that are represented by generators of operators.

It has three branches: [`StaticGenerator`](@ref), and [`ParametricGenerator`](@ref) (which includes [`CategorizedGenerator`](@ref) and [`OperatorGenerator`](@ref)).

`Generator` also serves as a constructor factory, dispatching to the appropriate concrete subtype based on the arguments.
"""
abstract type Generator{V} <: LatticeModel end
@inline Base.valtype(::Type{<:Generator{V}}) where V = V

"""
    expand(gen::Generator) -> valtype(gen)
    expand(gen::Generator, ::Eager) -> valtype(gen)
    expand(gen::Generator, ::Lazy) -> valtype(gen)

Expand a `Generator`. Defaults to eager expansion which combines similar terms;
pass `lazy` to preserve separate terms.

Subtypes must implement `expand(gen, ::Lazy)`.
"""
@inline expand(gen::Generator) = expand(gen, eager)
@inline expand(gen::Generator, ::Eager) = expand!(zero(valtype(gen)), gen)

"""
    expand!(result, gen::Generator) -> typeof(result)

Expand the generator lazily and add the operators to `result`.
"""
function expand!(result, gen::Generator)
    for op in expand(gen, lazy)
        add!(result, op)
    end
    return result
end

"""
    iterate(gen::Generator)
    iterate(::Generator, state)

Iterate over a `Generator` lazily.
"""
@propagate_inbounds function Base.iterate(gen::Generator)
    ops = expand(gen, lazy)
    index = iterate(ops)
    isnothing(index) && return nothing
    return index[1], (ops, index[2])
end
@propagate_inbounds function Base.iterate(::Generator, state)
    index = iterate(state[1], state[2])
    isnothing(index) && return nothing
    return index[1], (state[1], index[2])
end

"""
    length(gen::Generator) -> Int

Get the number of operators after lazy expansion.
"""
@inline Base.length(gen::Generator) = length(expand(gen, lazy))

"""
    isempty(gen::Generator) -> Bool

Judge whether a `Generator` is empty after lazy expansion.
"""
@inline Base.isempty(gen::Generator) = iszero(length(gen))

"""
    StaticGenerator{M<:OperatorSet} <: Generator{M}

A `LatticeModel` that wraps a static `OperatorSet`.

Unlike `ParametricGenerator`, the operators have no parameter-update mechanism — they are treated as a fixed set.
"""
struct StaticGenerator{M<:OperatorSet} <: Generator{M}
    operators::M
end
@inline expand(gen::StaticGenerator, ::Lazy) = gen.operators

"""
    (transformation::LinearTransformation)(gen::StaticGenerator; kwargs...) -> StaticGenerator

Apply a linear transformation to a static generator of (representations of) quantum operators.
"""
@inline (transformation::LinearTransformation)(gen::StaticGenerator; kwargs...) = StaticGenerator(transformation(gen.operators; kwargs...))

"""
    Parameters(gen::StaticGenerator) -> Parameters

Get the parameters of a static generator, which are always empty.
"""
@inline Parameters(gen::StaticGenerator) = Parameters()

"""
    empty(gen::StaticGenerator) -> StaticGenerator
    empty!(gen::StaticGenerator) -> StaticGenerator

Get an empty copy of a static generator or empty it in place.
"""
@inline Base.empty(gen::StaticGenerator) = StaticGenerator(empty(gen.operators))
@inline function Base.empty!(gen::StaticGenerator)
    empty!(gen.operators)
    return gen
end

"""
    update!(gen::StaticGenerator; parameters...) -> StaticGenerator

Update the parameters of a static generator. Since a static generator has no parameters, this is a no-op that returns the generator unchanged.
"""
@inline update!(gen::StaticGenerator; parameters...) = gen

"""
    update!(gen::StaticGenerator, transformation::LinearTransformation, source::StaticGenerator; kwargs...) -> StaticGenerator

Update the parameters of a static generator from a source. Since a static generator has no parameters, this is a no-op that returns the generator unchanged.
"""
@inline update!(gen::StaticGenerator, transformation::LinearTransformation, source::StaticGenerator; kwargs...) = gen

"""
    reset!(gen::StaticGenerator, transformation::LinearTransformation, source::StaticGenerator; kwargs...) -> StaticGenerator

Reset a static generator from a source by applying the linear transformation to the source operators.
"""
@inline function reset!(gen::StaticGenerator, transformation::LinearTransformation, source::StaticGenerator; kwargs...)
    add!(empty!(gen.operators), transformation, source.operators; kwargs...)
    return gen
end

"""
    ParametricGenerator{V} <: Generator{V}

Abstract supertype for generators that carry parameters and support [`update!`](@ref).

Its concrete subtypes are [`CategorizedGenerator`](@ref) and [`OperatorGenerator`](@ref).
"""
abstract type ParametricGenerator{V} <: Generator{V} end

"""
    CategorizedGenerator{V, C, A<:NamedTuple, B<:NamedTuple, P<:Parameters, D<:Boundary} <: ParametricGenerator{V}

Parameterized generator that groups the (representations of) quantum operators in a quantum lattice system into three categories, i.e., the constant, the alterable, and the boundary.
"""
mutable struct CategorizedGenerator{V, C, A<:NamedTuple, B<:NamedTuple, P<:Parameters, D<:Boundary} <: ParametricGenerator{V}
    const constops::C
    const alterops::A
    const boundops::B
    parameters::P
    const boundary::D
    function CategorizedGenerator(constops, alterops::NamedTuple, boundops::NamedTuple, parameters::Parameters, boundary::Boundary)
        C, A, B = typeof(constops), typeof(alterops), typeof(boundops)
        new{commontype(C, A, B), C, A, B, typeof(parameters), typeof(boundary)}(constops, alterops, boundops, parameters, boundary)
    end
end
@inline @generated function commontype(::Type{C}, ::Type{A}, ::Type{B}) where {C, A<:NamedTuple, B<:NamedTuple}
    exprs = [:(optp = C)]
    fieldcount(A)>0 && append!(exprs, [:(optp = promote_type(optp, $T)) for T in fieldtypes(A)])
    fieldcount(B)>0 && append!(exprs, [:(optp = promote_type(optp, $T)) for T in fieldtypes(B)])
    push!(exprs, :(return optp))
    return Expr(:block, exprs...)
end
function expand(cat::CategorizedGenerator, ::Lazy)
    params = (one(eltype(cat.parameters)), values(cat.parameters, keys(cat.alterops)|>Val)..., values(cat.parameters, keys(cat.alterops)|>Val)...)
    counts = (length(cat.constops), map(length, values(cat.alterops))..., map(length, values(cat.boundops))...)
    ops = (cat.constops, values(cat.alterops)..., values(cat.boundops)...)
    return CategorizedGeneratorExpand{eltype(cat)}(flatten(map((param, count)->repeated(param, count), params, counts)), flatten(ops))
end
@inline @generated function Base.values(parameters::Parameters, ::Val{KS}) where KS
    exprs = [:(getfield(parameters, $name)) for name in QuoteNode.(KS)]
    return Expr(:tuple, exprs...)
end
struct CategorizedGeneratorExpand{M<:OperatorPack, VS, OS} <: OperatorSet{M}
    values::VS
    ops::OS
    CategorizedGeneratorExpand{M}(values, ops) where {M<:OperatorPack} = new{M, typeof(values), typeof(ops)}(values, ops)
end
@inline Base.length(ee::CategorizedGeneratorExpand) = mapreduce(length, +, ee.values.it)
@propagate_inbounds function Base.iterate(ee::CategorizedGeneratorExpand, state=((), ()))
    v = iterate(ee.values, state[1])
    isnothing(v) && return nothing
    op = iterate(ee.ops, state[2])
    return op[1]*v[1], (v[2], op[2])
end

"""
    (transformation::LinearTransformation)(cat::CategorizedGenerator; kwargs...) -> CategorizedGenerator

Apply a linear transformation to a categorized generator of (representations of) quantum operators.
"""
function (transformation::LinearTransformation)(cat::CategorizedGenerator; kwargs...)
    wrapper(m) = transformation(m; kwargs...)
    constops = wrapper(cat.constops)
    alterops = NamedTuple{keys(cat.alterops)}(map(wrapper, values(cat.alterops)))
    boundops = NamedTuple{keys(cat.boundops)}(map(wrapper, values(cat.boundops)))
    return CategorizedGenerator(constops, alterops, boundops, cat.parameters, deepcopy(cat.boundary))
end

"""
    Parameters(cat::CategorizedGenerator)

Get the complete set of parameters of a categorized generator of (representations of) quantum operators.
"""
@inline Parameters(cat::CategorizedGenerator) = merge(cat.parameters, Parameters(cat.boundary))

"""
    empty(cat::CategorizedGenerator) -> CategorizedGenerator
    empty!(cat::CategorizedGenerator) -> CategorizedGenerator

Get an empty copy of a categorized generator or empty a categorized generator of (representations of) quantum operators.
"""
@inline function Base.empty(cat::CategorizedGenerator)
    constops = empty(cat.constops)
    alterops = NamedTuple{keys(cat.alterops)}(map(empty, values(cat.alterops)))
    boundops = NamedTuple{keys(cat.boundops)}(map(empty, values(cat.boundops)))
    return CategorizedGenerator(constops, alterops, boundops, cat.parameters, deepcopy(cat.boundary))
end
@inline function Base.empty!(cat::CategorizedGenerator)
    empty!(cat.constops)
    map(empty!, values(cat.alterops))
    map(empty!, values(cat.boundops))
    return cat
end

"""
    update!(cat::CategorizedGenerator{<:OperatorSum}; parameters...) -> CategorizedGenerator

Update the parameters (including the boundary parameters) of a categorized generator of (representations of) quantum operators.

!!! Note
    The coefficients of `boundops` are also updated due to the change of the boundary parameters.
"""
function update!(cat::CategorizedGenerator{<:OperatorSum}; parameters...)
    cat.parameters = update(cat.parameters; parameters...)
    if !match(Parameters(cat.boundary), NamedTuple{keys(parameters)}(values(parameters)))
        old = copy(cat.boundary.values)
        update!(cat.boundary; parameters...)
        map(ops->map!(op->cat.boundary(op, origin=old), ops), values(cat.boundops))
    end
    return cat
end

"""
    update!(cat::CategorizedGenerator, transformation::LinearTransformation, source::CategorizedGenerator; kwargs...) -> CategorizedGenerator

Update the parameters (including the boundary parameters) of a categorized generator based on its source categorized generator of (representations of) quantum operators and the corresponding linear transformation.

!!! Note
    The coefficients of `boundops` are also updated due to the change of the boundary parameters.
"""
function update!(cat::CategorizedGenerator, transformation::LinearTransformation, source::CategorizedGenerator; kwargs...)
    cat.parameters = update(cat.parameters; source.parameters...)
    if !match(Parameters(cat.boundary), Parameters(source.boundary))
        update!(cat.boundary; Parameters(source.boundary)...)
        map((dest, ops)->add!(empty!(dest), transformation, ops; kwargs...), values(cat.boundops), values(source.boundops))
    end
    return cat
end

"""
    reset!(cat::CategorizedGenerator, transformation::LinearTransformation, source::CategorizedGenerator; kwargs...)

Reset a categorized generator by its source categorized generator of (representations of) quantum operators and the corresponding linear transformation.
"""
function reset!(cat::CategorizedGenerator, transformation::LinearTransformation, source::CategorizedGenerator; kwargs...)
    add!(empty!(cat.constops), transformation, source.constops; kwargs...)
    map((dest, ops)->add!(empty!(dest), transformation, ops; kwargs...), values(cat.alterops), values(source.alterops))
    map((dest, ops)->add!(empty!(dest), transformation, ops; kwargs...), values(cat.boundops), values(source.boundops))
    cat.parameters = update(cat.parameters; source.parameters...)
    merge!(cat.boundary, source.boundary)
    return cat
end

"""
    OperatorGenerator{V<:Operators, CG<:CategorizedGenerator{V}, B<:Bond, H<:Hilbert, TS<:ZeroAtLeast{Term}} <: ParametricGenerator{V}

A generator of operators based on the terms, bonds and Hilbert space of a quantum lattice system.
"""
struct OperatorGenerator{V<:Operators, CG<:CategorizedGenerator{V}, B<:Bond, H<:Hilbert, TS<:ZeroAtLeast{Term}} <: ParametricGenerator{V}
    operators::CG
    bonds::Vector{B}
    hilbert::H
    terms::TS
    half::Bool
    function OperatorGenerator(operators::CategorizedGenerator, bonds::Vector{<:Bond}, hilbert::Hilbert, terms::OneOrMore{Term}, half::Bool)
        terms = OneOrMore(terms)
        new{valtype(operators), typeof(operators), eltype(bonds), typeof(hilbert), typeof(terms)}(operators, bonds, hilbert, terms, half)
    end
end
@inline expand(gen::OperatorGenerator, ::Lazy) = expand(gen.operators, lazy)

"""
    OperatorGenerator(bonds::Vector{<:Bond}, hilbert::Hilbert, terms::OneOrMore{Term}, boundary::Boundary=plain; half::Bool=false)

Construct a generator of quantum operators based on the input bonds, Hilbert space, terms and (twisted) boundary condition.

When the boundary condition is [`plain`](@ref), the boundary operators will be set to be empty for simplicity and efficiency.
"""
function OperatorGenerator(bonds::Vector{<:Bond}, hilbert::Hilbert, terms::OneOrMore{Term}, boundary::Boundary=plain; half::Bool=false)
    emptybonds = eltype(bonds)[]
    innerbonds, boundbonds = if boundary == plain
        bonds, eltype(bonds)[]
    else
        filter(isintracell, bonds), filter((!)∘isintracell, bonds)
    end
    terms = OneOrMore(terms)
    constops = Operators{mapreduce(term->operatortype(eltype(bonds), typeof(hilbert), typeof(term)), promote_type, terms)}()
    map(term->expand!(constops, term, term.ismodulatable ? emptybonds : innerbonds, hilbert; half=half), terms)
    alterops = NamedTuple{map(id, terms)}(expansion(terms, emptybonds, innerbonds, hilbert, scalartype(constops); half=half))
    boundops = NamedTuple{map(id, terms)}(expansion(terms, boundbonds, hilbert, boundary, scalartype(constops); half=half))
    parameters = NamedTuple{map(id, terms)}(map(value, terms))
    return OperatorGenerator(CategorizedGenerator(constops, alterops, boundops, parameters, boundary), bonds, hilbert, terms, half)
end
function expansion(terms::ZeroAtLeast{Term}, emptybonds::Vector{<:Bond}, innerbonds::Vector{<:Bond}, hilbert::Hilbert, ::Type{V}; half) where V
    return map(terms) do term
        expand(replace(term, one(V)), term.ismodulatable ? innerbonds : emptybonds, hilbert; half=half)
    end
end
function expansion(terms::ZeroAtLeast{Term}, bonds::Vector{<:Bond}, hilbert::Hilbert, boundary::Boundary, ::Type{V}; half) where V
    return map(terms) do term
        O = promote_type(valtype(typeof(boundary), operatortype(eltype(bonds), typeof(hilbert), typeof(term))), V)
        map!(boundary, expand!(Operators{O}(), one(term), bonds, hilbert, half=half))
    end
end

"""
    Parameters(gen::OperatorGenerator) -> Parameters

Get the parameters of an `OperatorGenerator`.
"""
@inline Parameters(gen::OperatorGenerator) = Parameters(gen.operators)

"""
    update!(gen::OperatorGenerator; parameters...) -> typeof(gen)

Update the coefficients of the terms in a generator.
"""
@inline function update!(gen::OperatorGenerator; parameters...)
    update!(gen.operators; parameters...)
    map(term->(term.ismodulatable ? update!(term; parameters...) : term), gen.terms)
    return gen
end

"""
    empty(gen::OperatorGenerator) -> OperatorGenerator
    empty!(gen::OperatorGenerator) -> OperatorGenerator

Get an empty copy of or empty an operator generator.
"""
@inline function Base.empty(gen::OperatorGenerator)
    return OperatorGenerator(empty(gen.operators), empty(gen.bonds), empty(gen.hilbert), gen.terms, gen.half)
end
function Base.empty!(gen::OperatorGenerator)
    empty!(gen.operators)
    empty!(gen.bonds)
    empty!(gen.hilbert)
    return gen
end

"""
    reset!(gen::OperatorGenerator, bonds::AbstractVector{<:Bond}, hilbert::Hilbert; vectors::AbstractVector{<:AbstractVector}=gen.operators.boundary.vectors) -> OperatorGenerator

Reset an operator generator by a new lattice and the corresponding hilbert space.
"""
function reset!(gen::OperatorGenerator, bonds::AbstractVector{<:Bond}, hilbert::Hilbert; vectors::AbstractVector{<:AbstractVector}=gen.operators.boundary.vectors)
    append!(empty!(gen.bonds), bonds)
    merge!(empty!(gen.hilbert), hilbert)
    empty!(gen.operators)
    reset!(gen.operators.boundary, vectors)
    emptybonds = eltype(gen.bonds)[]
    innerbonds, boundbonds = if gen.operators.boundary == plain
        gen.bonds, eltype(gen.bonds)[]
    else
        filter(bond->isintracell(bond), gen.bonds), filter(bond->!isintracell(bond), gen.bonds)
    end
    map(term->expand!(gen.operators.constops, term, term.ismodulatable ? emptybonds : innerbonds, gen.hilbert; half=gen.half), gen.terms)
    map(term->expand!(getfield(gen.operators.alterops, id(term)), one(term), term.ismodulatable ? innerbonds : emptybonds, gen.hilbert; half=gen.half), gen.terms)
    map(term->map!(gen.operators.boundary, expand!(getfield(gen.operators.boundops, id(term)), one(term), boundbonds, gen.hilbert; half=gen.half)), gen.terms)
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
    for bond in gen.bonds
        if isintracell(bond)
            expand!(result, term, bond, gen.hilbert; half=gen.half)
        else
            for opt in expand(term, bond, gen.hilbert; half=gen.half)
                add!(result, gen.operators.boundary(opt))
            end
        end
    end
    return result
end
function expand(gen::OperatorGenerator, i::Int)
    bond = gen.bonds[i]
    result = zero(valtype(gen))
    map(term->expand!(result, term, bond, gen.hilbert; half=gen.half), gen.terms)
    isintracell(bond) || map!(gen.operators.boundary, result)
    return result
end
function expand(gen::OperatorGenerator, name::Symbol, i::Int)
    bond = gen.bonds[i]
    term = get(gen.terms, Val(name))
    result = expand!(zero(valtype(gen)), term, bond, gen.hilbert; half=gen.half)
    isintracell(bond) || map!(gen.operators.boundary, result)
    return result
end
@inline @generated function Base.get(terms::ZeroAtLeast{Term}, ::Val{Name}) where Name
    i = findfirst(isequal(Name), map(id, fieldtypes(terms)))::Int
    return :(terms[$i])
end

"""
    (transformation::LinearTransformation)(gen::OperatorGenerator; kwargs...) -> CategorizedGenerator

Get the transformation applied to a generator of quantum operators.
"""
@inline (transformation::LinearTransformation)(gen::OperatorGenerator; kwargs...) = transformation(gen.operators; kwargs...)

"""
    update!(cat::CategorizedGenerator, transformation::LinearTransformation, source::OperatorGenerator; kwargs...) -> CategorizedGenerator

Update the parameters (including the boundary parameters) of a categorized generator based on its source operator generator of (representations of) quantum operators and the corresponding linear transformation.

!!! Note
    The coefficients of `boundops` are also updated due to the change of the boundary parameters.
"""
@inline update!(cat::CategorizedGenerator, transformation::LinearTransformation, source::OperatorGenerator; kwargs...) = update!(cat, transformation, source.operators; kwargs...)

"""
    reset!(cat::CategorizedGenerator, transformation::LinearTransformation, source::OperatorGenerator; kwargs...)

Reset a categorized generator by its source operator generator of (representations of) quantum operators and the corresponding linear transformation.
"""
@inline reset!(cat::CategorizedGenerator, transformation::LinearTransformation, source::OperatorGenerator; kwargs...) = reset!(cat, transformation, source.operators; kwargs...)

"""
    Generator(ops::OperatorSet) -> StaticGenerator
    Generator(constops, alterops::NamedTuple, boundops::NamedTuple, parameters::Parameters, boundary::Boundary) -> CategorizedGenerator
    Generator(ops::CategorizedGenerator{<:Operators}, bonds::Vector{<:Bond}, hilbert::Hilbert, terms::OneOrMore{Term}, half::Bool) -> OperatorGenerator
    Generator(bonds::Vector{<:Bond}, hilbert::Hilbert, terms::OneOrMore{Term}, boundary::Boundary=plain; half::Bool=false) -> OperatorGenerator

Factory constructor for `Generator` subtypes.

Dispatches on the argument types to construct the appropriate concrete representation.
Unlike [`LatticeModel`](@ref), this factory only constructs `Generator` subtypes (not `Formula`).
"""
@inline Generator(ops::OperatorSet) = StaticGenerator(ops)
@inline Generator(constops, alterops::NamedTuple, boundops::NamedTuple, parameters::Parameters, boundary::Boundary) = CategorizedGenerator(constops, alterops, boundops, parameters, boundary)
@inline Generator(ops::CategorizedGenerator{<:Operators}, bonds::Vector{<:Bond}, hilbert::Hilbert, terms::OneOrMore{Term}, half::Bool) = OperatorGenerator(ops, bonds, hilbert, terms, half)
@inline Generator(bonds::Vector{<:Bond}, hilbert::Hilbert, terms::OneOrMore{Term}, boundary::Boundary=plain; half::Bool=false) = OperatorGenerator(bonds, hilbert, terms, boundary; half=half)

"""
    LatticeModel(expression::Function, parameters::Parameters) -> Formula
    LatticeModel(ops::OperatorSet) -> StaticGenerator
    LatticeModel(constops, alterops::NamedTuple, boundops::NamedTuple, parameters::Parameters, boundary::Boundary) -> CategorizedGenerator
    LatticeModel(ops::CategorizedGenerator{<:Operators}, bonds::Vector{<:Bond}, hilbert::Hilbert, terms::OneOrMore{Term}, half::Bool) -> OperatorGenerator
    LatticeModel(bonds::Vector{<:Bond}, hilbert::Hilbert, terms::OneOrMore{Term}, boundary::Boundary=plain; half::Bool=false) -> OperatorGenerator

Unified factory constructor for `LatticeModel` subtypes.

Dispatches on the argument types to construct the appropriate concrete representation.
`Frontend` subtypes and `Algorithm` are not constructed via this factory.
"""
@inline LatticeModel(expression::Function, parameters::Parameters) = Formula(expression, parameters)
@inline LatticeModel(ops::OperatorSet) = StaticGenerator(ops)
@inline LatticeModel(constops, alterops::NamedTuple, boundops::NamedTuple, parameters::Parameters, boundary::Boundary) = CategorizedGenerator(constops, alterops, boundops, parameters, boundary)
@inline LatticeModel(ops::CategorizedGenerator{<:Operators}, bonds::Vector{<:Bond}, hilbert::Hilbert, terms::OneOrMore{Term}, half::Bool) = OperatorGenerator(ops, bonds, hilbert, terms, half)
@inline LatticeModel(bonds::Vector{<:Bond}, hilbert::Hilbert, terms::OneOrMore{Term}, boundary::Boundary=plain; half::Bool=false) = OperatorGenerator(bonds, hilbert, terms, boundary; half=half)

"""
    Frontend <: LatticeModel

Frontend of algorithms applied to a quantum lattice system.
"""
abstract type Frontend <: LatticeModel end

"""
    Action

Abstract type for all actions.
"""
abstract type Action end
@inline Base.:(==)(action₁::Action, action₂::Action) = ==(efficientoperations, action₁, action₂)
@inline Base.isequal(action₁::Action, action₂::Action) = isequal(efficientoperations, action₁, action₂)

"""
    update!(action::Action; parameters...) -> Action

Update the parameters of an action.
"""
@inline update!(action::Action; parameters...) = action

"""
    Data

Abstract type for the data of an action.
"""
abstract type Data end
@inline Base.:(==)(data₁::Data, data₂::Data) = ==(efficientoperations, data₁, data₂)
@inline Base.isequal(data₁::Data, data₂::Data) = isequal(efficientoperations, data₁, data₂)

"""
    Tuple(data::Data)

Convert `Data` to `Tuple`.
"""
@inline @generated function Base.Tuple(data::Data)
    exprs = [:(getfield(data, $i)) for i in 1:fieldcount(data)]
    return Expr(:tuple, exprs...)
end

"""
    Assignment{A<:Action, P<:Parameters, M<:Function, N, D<:Data} <: Function

An assignment associated with an action.
"""
mutable struct Assignment{A<:Action, P<:Parameters, M<:Function, T<:Tuple, D<:Data} <: Function
    const dir::String
    const name::Symbol
    const action::A
    parameters::P
    const map::M
    const dependencies::T
    data::D
    function Assignment(::Type{D}, dir::String, name::Symbol, action::Action, parameters::Parameters, map::Function, dependencies::ZeroAtLeast{Assignment}) where {D<:Data}
        new{typeof(action), typeof(parameters), typeof(map), typeof(dependencies), D}(dir, name, action, parameters, map, dependencies)
    end
end
@inline Base.:(==)(assign₁::Assignment, assign₂::Assignment) = ==(efficientoperations, assign₁, assign₂)
@inline Base.isequal(assign₁::Assignment, assign₂::Assignment) = isequal(efficientoperations, assign₁, assign₂)

"""
    Parameters(assignment::Assignment) -> Parameters

Get the parameters of an assignment.
"""
@inline Parameters(assignment::Assignment) = assignment.parameters

"""
    valtype(assign::Assignment)
    valtype(::Type{<:Assignment})

Type of the data (result) of an assignment.
"""
@inline Base.valtype(assign::Assignment) = valtype(typeof(assign))
@inline Base.valtype(::Type{<:Assignment{<:Action, <:Parameters, <:Function, <:Tuple, D}}) where {D<:Data} = D

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
    show(io::IO, assign::Assignment)

Show an assignment.
"""
@inline Base.show(io::IO, assign::Assignment) = print(io, assign.name)

"""
    show(io::IO, ::MIME"text/plain", assign::Assignment)

Show an assignment.

Optionally, some parameters of the algorithm can be filtered by specifying the `:select` context in `io`. Besides, the maximum number of decimals of the parameters can also be specified by the `:ndecimal` context in `io`.
"""
function Base.show(io::IO, ::MIME"text/plain", assign::Assignment)
    io₁ = indent(io, 2)
    io₂ = indent(io, 4)
    print(io, assign.name)
    print(io₁, '\n', "action:")
    print(io₂, '\n', assign.action)
    print(io₁, '\n', "parameters:")
    select = get(io, :select, param->true)
    ndecimal = get(io, :ndecimal, 10)
    for (name, value) in pairs(assign.parameters)
        if select(name)
            print(io₂, '\n', name, ": ", str(value; ndecimal=ndecimal))
        end
    end
end

"""
    options(::Type{<:Assignment}) -> NamedTuple

Get the options of a certain type of `Assignment`.
"""
@inline options(::Type{<:Assignment}) = NamedTuple()

"""
    optionsinfo(::Type{A}; level::Int=1) where {A<:Assignment} -> String

Get the complete info of the options of a certain type of `Assignment`, including that of its dependencies.
"""
function optionsinfo(::Type{A}; level::Int=1) where {A<:Assignment}
    io = IOBuffer()
    indent = repeat(" ", 2level)
    print(io, "Assignment{<:$(nameof(fieldtype(A, :action)))} options:\n")
    options = Frameworks.options(A)
    for (i, (key, value)) in enumerate(pairs(options))
        print(io, indent, "($i) `:$key`: $value", i<length(options) ? ";" : ".", '\n')
    end
    for (i, D) in enumerate(fieldtypes(fieldtype(A, :dependencies)))
        print(io, '\n', indent, "Dependency $i) ", optionsinfo(D; level=level+1))
    end
    return String(take!(io))
end

"""
    hasoption(::Type{A}, option::Symbol) where {A<:Assignment} -> Bool

Judge whether a certain type of `Assignment` has an option.
"""
function hasoption(::Type{A}, option::Symbol) where {A<:Assignment}
    haskey(options(A), option) && return true
    for D in fieldtypes(fieldtype(A, :dependencies))
        hasoption(D, option) && return true
    end
    return false
end

"""
    checkoptions(::Type{A}; options...) where {A<:Assignment}

Check whether the keyword arguments are legal options of a certain type of `Assignment`.
"""
@inline function checkoptions(::Type{A}; options...) where {A<:Assignment}
    for candidate in keys(options)
        @assert(hasoption(A, candidate), "checkoptions error: improper option(`:$candidate`). See following.\n$(optionsinfo(A))")
    end
end

"""
    Algorithm{F<:Frontend, P<:Parameters, M<:Function} <: LatticeModel

An algorithm associated with a frontend.
"""
mutable struct Algorithm{F<:Frontend, P<:Parameters, M<:Function} <: LatticeModel
    const dir::String
    const name::Symbol
    const frontend::F
    parameters::P
    const map::M
    const timer::TimerOutput
end
@inline function Base.:(==)(algorithm₁::Algorithm, algorithm₂::Algorithm)
    return ==(
        (algorithm₁.dir, algorithm₁.name, algorithm₁.frontend, algorithm₁.parameters, algorithm₁.map),
        (algorithm₂.dir, algorithm₂.name, algorithm₂.frontend, algorithm₂.parameters, algorithm₂.map),
    )
end
@inline function Base.isequal(algorithm₁::Algorithm, algorithm₂::Algorithm)
    return isequal(
        (algorithm₁.dir, algorithm₁.name, algorithm₁.frontend, algorithm₁.parameters, algorithm₁.map),
        (algorithm₂.dir, algorithm₂.name, algorithm₂.frontend, algorithm₂.parameters, algorithm₂.map),
    )
end
@inline Base.valtype(::Type{<:Algorithm{F}}) where {F<:Frontend} = valtype(F)

"""
    Algorithm(name::Symbol, frontend::Frontend, parameters::Parameters=Parameters(frontend), map::Function=identity; dir::String=".", timer::TimerOutput=TimerOutput())

Construct an algorithm.
"""
@inline function Algorithm(name::Symbol, frontend::Frontend, parameters::Parameters=Parameters(frontend), map::Function=identity; dir::String=".", timer::TimerOutput=TimerOutput())
    return Algorithm(dir, name, frontend, parameters, map, timer)
end

"""
    datatype(::Type{A}, ::Type{F}) where {A<:Action, F<:Frontend}

Get the concrete subtype of `Data` according to the types of an `Action` and a `Frontend`.
"""
@inline function datatype(::Type{A}, ::Type{F}) where {A<:Action, F<:Frontend}
    D = Core.Compiler.return_type(run!, Tuple{Algorithm{F}, Assignment{A}})
    @assert isconcretetype(D) && D<:Data "datatype error: failure ($D) of the default method for type $A and type $F."
    return D
end

"""
    Parameters(algorithm::Algorithm)

Get the parameters of an algorithm.
"""
@inline Parameters(algorithm::Algorithm) = algorithm.parameters

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
    show(io::IO, alg::Algorithm)

Show an algorithm.
"""
@inline Base.show(io::IO, alg::Algorithm) = print(io, alg.name, "-", nameof(typeof(alg.frontend)))

"""
    show(io::IO, ::MIME"text/plain", alg::Algorithm)

Show an algorithm.

Optionally, some parameters of the algorithm can be filtered by specifying the `:select` context in `io`. Besides, the maximum number of decimals of the parameters can also be specified by the `:ndecimal` context in `io`.
"""
function Base.show(io::IO, ::MIME"text/plain", alg::Algorithm)
    io₁ = indent(io, 2)
    io₂ = indent(io, 4)
    print(io, alg.name)
    print(io₁, '\n', "frontend:")
    print(io₂, '\n', alg.frontend)
    print(io₁, '\n', "parameters:")
    select = get(io, :select, param->true)
    ndecimal = get(io, :ndecimal, 10)
    for (name, value) in pairs(alg.parameters)
        if select(name)
            print(io₂, '\n', name, ": ", str(value; ndecimal=ndecimal))
        end
    end
end

"""
    id(obj::Union{Assignment, Algorithm})

Get the identifier of an assignment/algorithm.
"""
@inline id(obj::Union{Assignment, Algorithm}) = string(obj)

"""
    summary(alg::Algorithm)

Provide a summary of an algorithm.
"""
function Base.summary(alg::Algorithm)
    @info "Summary of $(alg.name) with $(nameof(typeof(alg.frontend))) frontend:"
    @info string(alg.timer)
end

"""
    (alg::Algorithm)(assign::Assignment, checkoptions::Bool=true; options...) -> Assignment
    (assign::Assignment)(alg::Algorithm, checkoptions::Bool=true; options...) -> Assignment

Run an assignment based on an algorithm.

The difference between these two methods is that the first uses the parameters of `assign` as the current parameters while the second uses those of `alg`.
"""
function (alg::Algorithm)(assign::Assignment, checkoptions::Bool=true; options...)
    @timeit alg.timer id(assign) begin
        checkoptions && Frameworks.checkoptions(typeof(assign); options...)
        ismatched = match(assign.parameters, alg.parameters)
        if !(isdefined(assign, :data) && ismatched)
            ismatched || update!(alg; assign.parameters...)
            map(dependency->dependency(alg, false; options...), assign.dependencies)
            @timeit alg.timer "run!" (assign.data = run!(alg, assign; options...))
        end
    end
    return assign
end
function (assign::Assignment)(alg::Algorithm, checkoptions::Bool=true; options...)
    @timeit alg.timer id(assign) begin
        checkoptions && Frameworks.checkoptions(typeof(assign); options...)
        ismatched = match(alg.parameters, assign.parameters)
        if !(isdefined(assign, :data) && ismatched)
            ismatched || update!(assign; alg.parameters...)
            map(dependency->dependency(alg, false; options...), assign.dependencies)
            @timeit alg.timer "run!" (assign.data = run!(alg, assign; options...))
        end
    end
    return assign
end

"""
    (alg::Algorithm)(name::Symbol, action::Action, dependencies::ZeroOrMore{Assignment}; dir::String=alg.dir, delay::Bool=false, options...) -> Assignment
    (alg::Algorithm)(name::Symbol, action::Action, parameters::Parameters, dependencies::ZeroOrMore{Assignment}; dir::String=alg.dir, delay::Bool=false, options...) -> Assignment
    (alg::Algorithm)(name::Symbol, action::Action, parameters::Parameters=Parameters(), map::Function=identity, dependencies::ZeroOrMore{Assignment}=(); dir::String=alg.dir, delay::Bool=false, options...) -> Assignment

Add an assignment on an algorithm by providing the contents of the assignment, and run this assignment.
"""
@inline function (alg::Algorithm)(name::Symbol, action::Action, dependencies::ZeroOrMore{Assignment}; dir::String=alg.dir, delay::Bool=false, options...)
    return alg(name, action, Parameters(), dependencies; delay=delay, dir=dir, options...)
end
@inline function (alg::Algorithm)(name::Symbol, action::Action, parameters::Parameters, dependencies::ZeroOrMore{Assignment}; dir::String=alg.dir, delay::Bool=false, options...)
    return alg(name, action, parameters, identity, dependencies; delay=delay, dir=dir, options...)
end
@inline function (alg::Algorithm)(name::Symbol, action::Action, parameters::Parameters=Parameters(), map::Function=identity, dependencies::ZeroOrMore{Assignment}=(); dir::String=alg.dir, delay::Bool=false, options...)
    assign = Assignment(datatype(typeof(action), typeof(alg.frontend)), dir, name, action, merge(alg.parameters, parameters), map, ZeroOrMore(dependencies))
    delay || begin
        alg(assign; options...)
        @info "Assignment $name: time consumed $(time(alg.timer[id(assign)])/10^9)s."
    end
    return assign
end

"""
    run!(alg::Algorithm, assign::Assignment; options...)

Run an assignment based on an algorithm.
"""
function run! end

"""
    dirname(obj::Union{Assignment, Algorithm}) -> String

Get the dirname of the data file of an assignment/algorithm.
"""
@inline Base.dirname(obj::Union{Assignment, Algorithm}) = obj.dir

"""
    basename(obj::Union{Assignment, Algorithm}; prefix::String="", suffix::String="", extension::String="qld") -> String

Get the basename of the data file of an assignment/algorithm.
"""
@inline function Base.basename(obj::Union{Assignment, Algorithm}; prefix::String="", suffix::String="", extension::String="qld")
    return string(append(prefix, "-"), id(obj), prepend(suffix, "-"), prepend(extension, "."))
end

"""
    pathof(obj::Union{Assignment, Algorithm}; prefix::String="", suffix::String="", extension::String="qld") -> String

Get the path of the data file of an assignment/algorithm.
"""
@inline function Base.pathof(obj::Union{Assignment, Algorithm}; prefix::String="", suffix::String="", extension::String="qld")
    return joinpath(dirname(obj), basename(obj; prefix=prefix, suffix=suffix, extension=extension))
end

"""
    str(obj::Union{Assignment, Algorithm}; prefix::String="", suffix::String="", ndecimal::Int=10, select::Function=name::Symbol->true, front::String="", rear::String="") -> String

Get the string representation of an assignment/algorithm.
"""
@inline function str(obj::Union{Assignment, Algorithm}; prefix::String="", suffix::String="", ndecimal::Int=10, select::Function=name::Symbol->true, front::String="", rear::String="")
    base = basename(obj; prefix=prefix, suffix=suffix, extension="")
    parameters = str(Parameters(obj); ndecimal=ndecimal, select=select, front=front, rear=rear)
    return string(base, prepend(parameters, "-"))
end

"""
    qldsave(obj::Union{Assignment, Algorithm}, objs::Union{Assignment, Algorithm}...; mode::String="a+", prefix::String="", suffix::String="", ndecimal::Int=10, select::Function=name::Symbol->true, front::String="", rear::String="")

Save a series of assignments/algorithms to qld files.
"""
function qldsave(obj::Union{Assignment, Algorithm}, objs::Union{Assignment, Algorithm}...; mode::String="a+", prefix::String="", suffix::String="", ndecimal::Int=10, select::Function=name::Symbol->true, front::String="", rear::String="")
    @assert mode∈("a+", "w") "qldsave error: mode ($(repr(mode))) is not \"a+\" or \"w\"."
    map((obj, objs...)) do object
        qldsave(pathof(object; prefix=prefix, suffix=suffix), str(Parameters(object); ndecimal=ndecimal, select=select, front=front, rear=rear), object; mode=mode)
    end
end

"""
    qldsave(filename::String, args...; mode::String="a+")

Save arbitrary data as key-value pairs to a qld file.
"""
function qldsave(filename::String, args...; mode::String="a+")
    @assert mode∈("a+", "w") "qldsave error: mode ($(repr(mode))) is not \"a+\" or \"w\"."
    @assert iseven(length(args)) "qldsave error: wrong formed input data."
    jldopen(filename, mode) do file
        io = IOBuffer()
        for i in 1:2:length(args)
            serialize(io, args[i+1])
            file[string(args[i])] = take!(io)
        end
    end
end

"""
    qldload(filename::String) -> Dict{String, Any}
    qldload(filename::String, name::String) -> Any
    qldload(filename::String, name₁::String, name₂::String, names::String...) -> Tuple

Load data from a qld file.
With a single filename argument, returns all entries as a flat `Dict`.
With a key, returns that single entry.
With multiple keys, returns a `Tuple`.
"""
function qldload(filename::String)
    raw = jldopen(filename, "r") do file
        loadtodict!(Dict{String, Any}(), file)
    end
    result = Dict{String, Any}()
    for (key, bytes) in raw
        bytes isa Vector{UInt8} || continue
        result[key] = deserialize(IOBuffer(bytes))
    end
    return result
end
function qldload(filename::String, name::String)
    return jldopen(filename, "r") do file
        haskey(file, name) || error("qldload error: key '$name' not found in '$filename'.")
        data = file[name]
        data isa Vector{UInt8} || error("qldload error: key '$name' is not a data entry.")
        deserialize(IOBuffer(data))
    end
end
function qldload(filename::String, name₁::String, name₂::String, names::String...)
    return jldopen(filename, "r") do file
        map((name₁, name₂, names...)) do name
            haskey(file, name) || error("qldload error: key '$name' not found.")
            data = file[name]
            data isa Vector{UInt8} || error("qldload error: key '$name' is not a data entry.")
            deserialize(IOBuffer(data))
        end
    end
end

"""
    fingerprint(obj; ndecimal::Int=10) -> String

Generate a deterministic string fingerprint for an object, combining the string representation of `id(obj)` with that of `Parameters(obj)` (with numeric values rounded to `ndecimal` digits). When `id(obj)`/`Parameters(obj)` are not defined, the type name and an empty parameter list are used as fallbacks respectively.
"""
function fingerprint(obj; ndecimal::Int=10)
    id_str = try
        string(id(obj))
    catch
        string(nameof(typeof(obj)))
    end
    p = try
        Parameters(obj)
    catch
        Parameters()
    end
    if !isempty(p)
        return str(p; ndecimal=ndecimal, front=id_str)
    else
        return id_str
    end
end

"""
    dlmsave(assignment::Assignment, delim='\t'; prefix::String="", suffix::String="", ndecimal::Int=10, select::Function=name::Symbol->true, front::String="", rear::String="")

Save the data of an assignment to a delimited file.
"""
@inline function dlmsave(assignment::Assignment, delim='\t'; prefix::String="", suffix::String="", ndecimal::Int=10, select::Function=name::Symbol->true, front::String="", rear::String="")
    dlmsave(joinpath(dirname(assignment), string(str(assignment; prefix=prefix, suffix=suffix, ndecimal=ndecimal, select=select, front=front, rear=rear), ".dlm")), Tuple(assignment.data)..., delim)
end

"""
    seriestype(data::Data) -> Union{Symbol, Nothing}

Determine the plot seriestype for data. Returns `:path`, `:heatmap`, or `nothing`.
"""
@inline seriestype(data::Data) = seriestype(Tuple(data)...)
@inline seriestype(_...) = nothing
@inline seriestype(::AbstractVector{<:Number}, ::Union{AbstractVector{<:Number}, AbstractMatrix{<:Number}}, _...) = :path
@inline seriestype(::AbstractVector{<:Number}, ::AbstractVector{<:Number}, ::Union{AbstractMatrix{<:Number}, AbstractArray{<:Number, 3}}, _...) = :heatmap

end  # module
