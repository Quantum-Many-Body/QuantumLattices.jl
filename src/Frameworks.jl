module Frameworks

using Base: @propagate_inbounds
using Base.Iterators: flatten, repeated
using DelimitedFiles: writedlm
using LaTeXStrings: latexstring
using Printf: @printf, @sprintf
using RecipesBase: RecipesBase, @recipe
using Serialization: serialize
using TimerOutputs: TimerOutput, TimerOutputs, @timeit
using ..DegreesOfFreedom: plain, Boundary, Hilbert, Term
using ..QuantumLattices: id, value
using ..QuantumOperators: OperatorPack, Operators, OperatorSet, OperatorSum, LinearTransformation, Transformation, identity, operatortype
using ..Spatials: Bond, isintracell
using ..Toolkit: atol, efficientoperations, rtol, parametertype, tostr

import ..QuantumLattices: add!, expand, expand!, reset!, update, update!
import ..QuantumOperators: scalartype
import ..Spatials: save

export Action, Algorithm, Assignment, CategorizedGenerator, Eager, ExpansionStyle, Formula, Frontend, Generator, Lazy, OneOrMore, OperatorGenerator, Parameters
export checkoptions, eager, lazy, initialize, options, prepare!, run!

"""
    const OneOrMore{A} = Union{A, Tuple{A, Vararg{A}}}

One or more something.
"""
const OneOrMore{A} = Union{A, Tuple{A, Vararg{A}}}

"""
    OneOrMore(x) -> Tuple{typeof(x)}
    OneOrMore(x::Tuple) -> typeof(x)

If `x` is a tuple, return itself; if not, return `(x,)`.
"""
@inline OneOrMore(x) = (x,)
@inline OneOrMore(xs::Tuple) = xs

"""
    Parameters{Names}(values::Number...) where Names

A NamedTuple that contains the key-value pairs.
"""
const Parameters{Names, T<:Tuple{Vararg{Number}}} = NamedTuple{Names, T}
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
    Parameters(bound::Boundary)

Get the parameters of the twisted boundary condition.
"""
@inline Parameters(bound::Boundary) = NamedTuple{keys(bound)}(ntuple(i->bound.values[i], Val(fieldcount(typeof(keys(bound))))))

"""
    Parameters(ops::OperatorSet) -> NamedTuple{(), Tuple{}}

Get the parameters of an `OperatorSet`, which is defined to be an empty `NamedTuple`.
"""
@inline Parameters(ops::OperatorSet) = NamedTuple()

"""
    Formula{V, F<:Function, P<:Parameters}

Representation of a quantum lattice system with an explicit analytical formula.
"""
mutable struct Formula{V, F<:Function, P<:Parameters}
    const expression::F
    parameters::P
    function Formula(expression::Function, parameters::Parameters)
        V = Core.Compiler.return_type(expression, parametertype(typeof(parameters), 2))
        @assert isconcretetype(V) "Formula error: input expression is not type-stable."
        new{V, typeof(expression), typeof(parameters)}(expression, parameters)
    end
end
@inline Base.:(==)(formula₁::Formula, formula₂::Formula) = ==(efficientoperations, formula₁, formula₂)
@inline Base.isequal(formula₁::Formula, formula₂::Formula) = isequal(efficientoperations, formula₁, formula₂)

"""
    valtype(formula::Formula)
    valtype(::Type{<:Formula{V}})

Get the valtype of a `Formula`.
"""
@inline Base.valtype(formula::Formula) = valtype(typeof(formula))
@inline Base.valtype(::Type{<:Formula{V}}) where V = V

"""
    scalartype(::Type{F}) where {F<:Formula}

Get the scalar type of a `Formula`.
"""
@inline scalartype(::Type{F}) where {F<:Formula} = scalartype(valtype(F))

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
    Generator{V}

Generator of (representations of) quantum operators in a quantum lattice system.
"""
abstract type Generator{V} end
@inline Base.:(==)(gen₁::Generator, gen₂::Generator) = ==(efficientoperations, gen₁, gen₂)
@inline Base.isequal(gen₁::Generator, gen₂::Generator) = isequal(efficientoperations, gen₁, gen₂)
@inline Base.show(io::IO, ::MIME"text/latex", gen::Generator) = show(io, MIME"text/latex"(), latexstring(latexstring(expand(gen))))
@inline ExpansionStyle(gen::Generator) = ExpansionStyle(typeof(gen))
@inline Base.IteratorSize(::Type{<:Generator}) = Base.SizeUnknown()

"""
    valtype(gen::Generator)
    valtype(::Type{<:Generator{V}}) where V

Get the valtype of a `Generator`.
"""
@inline Base.valtype(gen::Generator) = valtype(typeof(gen))
@inline Base.valtype(::Type{<:Generator{V}}) where V = V

"""
    eltype(gen::Generator)
    eltype(::Type{T}) where {T<:Generator}

Get the eltype of a `Generator`.
"""
@inline Base.eltype(gen::Generator) = eltype(typeof(gen))
@inline Base.eltype(::Type{T}) where {T<:Generator} = eltype(valtype(T))

"""
    scalartype(::Type{T}) where {T<:Generator}

Get the scalar type of a `Generator`.
"""
@inline scalartype(::Type{T}) where {T<:Generator} = scalartype(valtype(T))

"""
    iterate(gen::Generator)
    iterate(::Generator, state)

Iterate over a `Generator`.
"""
@propagate_inbounds function Base.iterate(gen::Generator)
    ops = expand(gen)
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
    expand(gen::Generator)
    expand(gen::Generator, ::Eager)
    expand(gen::Generator, ::Lazy)

Expand the generator to get the (representations of) quantum operators in a quantum lattice system.
"""
@inline expand(gen::Generator) = expand(gen, ExpansionStyle(gen))
@inline expand(gen::Generator, ::Eager) = expand!(zero(valtype(gen)), gen)
@inline expand(gen::Generator, ::Lazy) = error("expand! error: not implemented for $(nameof(typeof(gen))).")

"""
    expand!(result, gen::Generator) -> typeof(result)

Expand the generator to add the (representations of) quantum operators in a quantum lattice system to `result`.
"""
function expand!(result, gen::Generator)
    for op in expand(gen, lazy)
        add!(result, op)
    end
    return result
end

"""
    CategorizedGenerator{V, C, A<:NamedTuple, B<:NamedTuple, P<:Parameters, D<:Boundary, S<:ExpansionStyle} <: Generator{V}

Categorized generator that groups the (representations of) quantum operators in a quantum lattice system into three categories, i.e., the constant, the alterable, and the boundary.
"""
mutable struct CategorizedGenerator{V, C, A<:NamedTuple, B<:NamedTuple, P<:Parameters, D<:Boundary, S<:ExpansionStyle} <: Generator{V}
    const constops::C
    const alterops::A
    const boundops::B
    parameters::P
    const boundary::D
    function CategorizedGenerator(constops, alterops::NamedTuple, boundops::NamedTuple, parameters::Parameters, boundary::Boundary, style::ExpansionStyle)
        C, A, B = typeof(constops), typeof(alterops), typeof(boundops)
        new{commontype(C, A, B), C, A, B, typeof(parameters), typeof(boundary), typeof(style)}(constops, alterops, boundops, parameters, boundary)
    end
end
@inline ExpansionStyle(::Type{<:CategorizedGenerator{V, C, <:NamedTuple, <:NamedTuple, <:Parameters, <:Boundary, S} where {V, C}}) where {S<:ExpansionStyle} = S()
@inline @generated function commontype(::Type{C}, ::Type{A}, ::Type{B}) where {C, A<:NamedTuple, B<:NamedTuple}
    exprs = [:(optp = C)]
    fieldcount(A)>0 && append!(exprs, [:(optp = promote_type(optp, $T)) for T in fieldtypes(A)])
    fieldcount(B)>0 && append!(exprs, [:(optp = promote_type(optp, $T)) for T in fieldtypes(B)])
    push!(exprs, :(return optp))
    return Expr(:block, exprs...)
end

"""
    Generator(constops, alterops::NamedTuple, boundops::NamedTuple, parameters::Parameters, boundary::Boundary, style::ExpansionStyle) -> CategorizedGenerator

Construct a `CategorizedGenerator`.
"""
@inline function Generator(constops, alterops::NamedTuple, boundops::NamedTuple, parameters::Parameters, boundary::Boundary, style::ExpansionStyle)
    return CategorizedGenerator(constops, alterops, boundops, parameters, boundary, style)
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
    return CategorizedGenerator(constops, alterops, boundops, cat.parameters, deepcopy(cat.boundary), ExpansionStyle(cat))
end

"""
    Parameters(cat::CategorizedGenerator)

Get the complete set of parameters of a categorized generator of (representations of) quantum operators.
"""
@inline Parameters(cat::CategorizedGenerator) = merge(cat.parameters, Parameters(cat.boundary))

"""
    isempty(cat::CategorizedGenerator) -> Bool

Judge whether a categorized generator is empty.
"""
@inline Base.isempty(cat::CategorizedGenerator) = isempty(cat.constops) && all(map(isempty, values(cat.alterops))) && all(map(isempty, values(cat.boundops)))

"""
    empty(cat::CategorizedGenerator) -> CategorizedGenerator
    empty!(cat::CategorizedGenerator) -> CategorizedGenerator

Get an empty copy of a categorized generator or empty a categorized generator of (representations of) quantum operators.
"""
@inline function Base.empty(cat::CategorizedGenerator)
    constops = empty(cat.constops)
    alterops = NamedTuple{keys(cat.alterops)}(map(empty, values(cat.alterops)))
    boundops = NamedTuple{keys(cat.boundops)}(map(empty, values(cat.boundops)))
    return CategorizedGenerator(constops, alterops, boundops, cat.parameters, deepcopy(cat.boundary), ExpansionStyle(cat))
end
@inline function Base.empty!(cat::CategorizedGenerator)
    empty!(cat.constops)
    map(empty!, values(cat.alterops))
    map(empty!, values(cat.boundops))
    return cat
end

"""
    *(cat::CategorizedGenerator, factor) -> CategorizedGenerator
    *(factor, cat::CategorizedGenerator) -> CategorizedGenerator

Multiply a categorized generator of (representations of) quantum operators with a factor.
"""
@inline Base.:*(cat::CategorizedGenerator, factor) = factor * cat
@inline function Base.:*(factor, cat::CategorizedGenerator)
    parameters = NamedTuple{keys(cat.parameters)}(map(value->factor*value, values(cat.parameters)))
    return CategorizedGenerator(factor*cat.constops, cat.alterops, cat.boundops, parameters, cat.boundary, ExpansionStyle(cat))
end

"""
    +(cat₁::CategorizedGenerator, cat₂::CategorizedGenerator) -> CategorizedGenerator

Addition of two categorized generators of (representations of) quantum operators.
"""
function Base.:+(cat₁::CategorizedGenerator, cat₂::CategorizedGenerator)
    @assert cat₁.boundary==cat₂.boundary "+ error: in order to be added, two entries must share the same boundary condition (including the twist angles at the boundary)."
    @assert ExpansionStyle(cat₁)==ExpansionStyle(cat₂) "+ error: in order to be added, two entries must share the same expansion style."
    constops = cat₁.constops + cat₂.constops
    alls, allshares = totalkeys(cat₁.parameters, cat₂.parameters), sharedkeys(cat₁.parameters, cat₂.parameters)
    allmatches = NamedTuple{keymaps(allshares)}(map(key->opsmatch(cat₁.alterops, cat₂.alterops, key) && opsmatch(cat₁.boundops, cat₂.boundops, key), allshares))
    parameters = NamedTuple{keymaps(alls)}(map(key->combinevalue(cat₁.parameters, cat₂.parameters, allmatches, key), alls))
    alteralls, altershares = totalkeys(cat₁.alterops, cat₂.alterops), sharedkeys(cat₁.alterops, cat₂.alterops)
    boundalls, boundshares = totalkeys(cat₁.boundops, cat₂.boundops), sharedkeys(cat₁.boundops, cat₂.boundops)
    altermatches = NamedTuple{keymaps(altershares)}(map(((::Val{key}) where key)->getfield(allmatches, key), altershares))
    boundmatches = NamedTuple{keymaps(boundshares)}(map(((::Val{key}) where key)->getfield(allmatches, key), boundshares))
    alterops = NamedTuple{keymaps(alteralls)}(map(key->combineops(cat₁.alterops, cat₁.parameters, cat₂.alterops, cat₂.parameters, altermatches, key), alteralls))
    boundops = NamedTuple{keymaps(boundalls)}(map(key->combineops(cat₁.boundops, cat₁.parameters, cat₂.boundops, cat₂.parameters, boundmatches, key), boundalls))
    return CategorizedGenerator(constops, alterops, boundops, parameters, deepcopy(cat₁.boundary), ExpansionStyle(cat₁))
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
    expand(cat::CategorizedGenerator, ::Lazy)

Expand a categorized generator to get the (representations of) quantum operators in a quantum lattice system.
"""
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
    update!(cat::CategorizedGenerator, transformation::LinearTransformation, source::CategorizedGenerator{<:Operators}; kwargs...) -> CategorizedGenerator

Update the parameters (including the boundary parameters) of a categorized generator based on its source categorized generator of (representations of) quantum operators and the corresponding linear transformation.

!!! Note
    The coefficients of `boundops` are also updated due to the change of the boundary parameters.
"""
function update!(cat::CategorizedGenerator, transformation::LinearTransformation, source::CategorizedGenerator{<:Operators}; kwargs...)
    cat.parameters = update(cat.parameters; source.parameters...)
    if !match(Parameters(cat.boundary), Parameters(source.boundary))
        update!(cat.boundary; Parameters(source.boundary)...)
        map((dest, ops)->add!(empty!(dest), transformation, ops; kwargs...), values(cat.boundops), values(source.boundops))
    end
    return cat
end

"""
    reset!(cat::CategorizedGenerator, transformation::LinearTransformation, source::CategorizedGenerator{<:Operators}; kwargs...)

Reset a categorized generator by its source categorized generator of (representations of) quantum operators and the corresponding linear transformation.
"""
function reset!(cat::CategorizedGenerator, transformation::LinearTransformation, source::CategorizedGenerator{<:Operators}; kwargs...)
    add!(empty!(cat.constops), transformation, source.constops; kwargs...)
    map((dest, ops)->add!(empty!(dest), transformation, ops; kwargs...), values(cat.alterops), values(source.alterops))
    map((dest, ops)->add!(empty!(dest), transformation, ops; kwargs...), values(cat.boundops), values(source.boundops))
    cat.parameters = update(cat.parameters; source.parameters...)
    merge!(cat.boundary, source.boundary)
    return cat
end

"""
    OperatorGenerator{V<:Operators, CG<:CategorizedGenerator{V}, B<:Bond, H<:Hilbert, TS<:Tuple{Vararg{Term}}} <: Generator{V}

A generator of operators based on the terms, bonds and Hilbert space of a quantum lattice system.
"""
struct OperatorGenerator{V<:Operators, CG<:CategorizedGenerator{V}, B<:Bond, H<:Hilbert, TS<:Tuple{Vararg{Term}}} <: Generator{V}
    operators::CG
    bonds::Vector{B}
    hilbert::H
    terms::TS
    half::Bool
end
@inline ExpansionStyle(::Type{<:OperatorGenerator{<:Operators, CG}}) where {CG<:CategorizedGenerator} = ExpansionStyle(CG)
@inline Base.valtype(::Type{<:OperatorGenerator{V}}) where {V<:Operators} = V
@inline expand(gen::OperatorGenerator, ::Lazy) = expand(gen.operators, lazy)

"""
    OperatorGenerator(bonds::Vector{<:Bond}, hilbert::Hilbert, terms::Tuple{Vararg{Term}}, boundary::Boundary=plain, style::ExpansionStyle=eager; half::Bool=false)

Construct a generator of quantum operators based on the input bonds, Hilbert space, terms and (twisted) boundary condition.

When the boundary condition is [`plain`](@ref), the boundary operators will be set to be empty for simplicity and efficiency.
"""
function OperatorGenerator(bonds::Vector{<:Bond}, hilbert::Hilbert, terms::Tuple{Vararg{Term}}, boundary::Boundary=plain, style::ExpansionStyle=eager; half::Bool=false)
    emptybonds = eltype(bonds)[]
    innerbonds, boundbonds = if boundary == plain
        bonds, eltype(bonds)[]
    else
        filter(isintracell, bonds), filter((!)∘isintracell, bonds)
    end
    constops = Operators{mapreduce(term->operatortype(eltype(bonds), typeof(hilbert), typeof(term)), promote_type, terms)}()
    map(term->expand!(constops, term, term.ismodulatable ? emptybonds : innerbonds, hilbert; half=half), terms)
    alterops = NamedTuple{map(id, terms)}(expansion(terms, emptybonds, innerbonds, hilbert, valtype(eltype(constops)); half=half))
    boundops = NamedTuple{map(id, terms)}(expansion(terms, boundbonds, hilbert, boundary, valtype(eltype(constops)); half=half))
    parameters = NamedTuple{map(id, terms)}(map(value, terms))
    return OperatorGenerator(CategorizedGenerator(constops, alterops, boundops, parameters, boundary, style), bonds, hilbert, terms, half)
end
function expansion(terms::Tuple{Vararg{Term}}, emptybonds::Vector{<:Bond}, innerbonds::Vector{<:Bond}, hilbert::Hilbert, ::Type{V}; half) where V
    return map(terms) do term
        expand(replace(term, one(V)), term.ismodulatable ? innerbonds : emptybonds, hilbert; half=half)
    end
end
function expansion(terms::Tuple{Vararg{Term}}, bonds::Vector{<:Bond}, hilbert::Hilbert, boundary::Boundary, ::Type{V}; half) where V
    return map(terms) do term
        O = promote_type(valtype(typeof(boundary), operatortype(eltype(bonds), typeof(hilbert), typeof(term))), V)
        map!(boundary, expand!(Operators{O}(), one(term), bonds, hilbert, half=half))
    end
end

"""
    Generator(operators::CategorizedGenerator{<:Operators}, bonds::Vector{<:Bond}, hilbert::Hilbert, terms::Tuple{Vararg{Term}}, half::Bool) -> OperatorGenerator
    Generator(bonds::Vector{<:Bond}, hilbert::Hilbert, terms::Tuple{Vararg{Term}}, boundary::Boundary=plain, style::ExpansionStyle=eager; half::Bool=false) -> OperatorGenerator

Construct an `OperatorGenerator`.
"""
@inline function Generator(operators::CategorizedGenerator{<:Operators}, bonds::Vector{<:Bond}, hilbert::Hilbert, terms::Tuple{Vararg{Term}}, half::Bool)
    return OperatorGenerator(operators, bonds, hilbert, terms, half)
end
@inline function Generator(bonds::Vector{<:Bond}, hilbert::Hilbert, terms::Tuple{Vararg{Term}}, boundary::Boundary=plain, style::ExpansionStyle=eager; half::Bool=false)
    return OperatorGenerator(bonds, hilbert, terms, boundary, style; half=half)
end

"""
    Parameters(gen::OperatorGenerator) -> Parameters

Get the parameters of an `OperatorGenerator`.
"""
@inline Parameters(gen::OperatorGenerator) = Parameters(gen.operators)

"""
    isempty(gen::OperatorGenerator) -> Bool

Judge whether an `OperatorGenerator` is empty.
"""
@inline Base.isempty(gen::OperatorGenerator) = isempty(gen.operators)

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

Reset an operator generator by a new lattice and the corresponding new hilbert space.
"""
function reset!(gen::OperatorGenerator, bonds::AbstractVector{<:Bond}, hilbert::Hilbert; vectors::AbstractVector{<:AbstractVector}=gen.operators.boundary.vectors)
    append!(empty!(gen.bonds), bonds)
    merge!(empty!(gen.hilbert), hilbert)
    empty!(gen.operators)
    merge!(gen.operators.boundary, replace(gen.operators.boundary; vectors=vectors))
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
@inline @generated function Base.get(terms::Tuple{Vararg{Term}}, ::Val{Name}) where Name
    i = findfirst(isequal(Name), map(id, fieldtypes(terms)))::Int
    return :(terms[$i])
end

"""
    (transformation::Transformation)(gen::OperatorGenerator; kwargs...) -> CategorizedGenerator

Get the transformation applied to a generator of quantum operators.
"""
@inline (transformation::Transformation)(gen::OperatorGenerator; kwargs...) = transformation(gen.operators; kwargs...)

"""
    Frontend

Frontend of algorithms applied to a quantum lattice system.
"""
abstract type Frontend end
@inline Base.:(==)(frontend₁::Frontend, frontend₂::Frontend) = ==(efficientoperations, frontend₁, frontend₂)
@inline Base.isequal(frontend₁::Frontend, frontend₂::Frontend) = isequal(efficientoperations, frontend₁, frontend₂)
@inline Base.show(io::IO, frontend::Frontend) = @printf io "%s" nameof(typeof(frontend))
@inline update!(frontend::Frontend; kwargs...) = error("update! error: not implemented for $(nameof(typeof(frontend))).")
@inline Parameters(frontend::Frontend) = error("Parameters error: not implemented for $(nameof(typeof(frontend))).")

"""
    Action

Abstract type for all actions.
"""
abstract type Action end
@inline Base.:(==)(action₁::Action, action₂::Action) = ==(efficientoperations, action₁, action₂)
@inline Base.isequal(action₁::Action, action₂::Action) = isequal(efficientoperations, action₁, action₂)
@inline initialize(action::Action, frontend::Frontend) = error("initialize error: not implemented.")
@inline update!(action::Action; parameters...) = action
@inline options(::Type{<:Action}) = Dict{Symbol, String}()
@inline function checkoptions(::Type{A}; kwargs...) where A
    legals = options(A)
    for candidate in keys(kwargs)
        @assert candidate∈keys(legals) "checkoptions error: improper option(`:$candidate`) for `$(nameof(A))`, which should be\n$(join(("$i) `:$key`: $value" for (i, (key, value)) in enumerate(legals)), ";\n"))."
    end
end

"""
    Assignment{A<:Action, P<:Parameters, M<:Function, N, D} <: Function

An assignment associated with an action.
"""
mutable struct Assignment{A<:Action, P<:Parameters, M<:Function, N, D} <: Function
    const id::Symbol
    const action::A
    parameters::P
    const map::M
    const dependencies::NTuple{N, Symbol}
    data::D
    ismatched::Bool
end
@inline Base.:(==)(assign₁::Assignment, assign₂::Assignment) = ==(efficientoperations, assign₁, assign₂)
@inline Base.isequal(assign₁::Assignment, assign₂::Assignment) = isequal(efficientoperations, assign₁, assign₂)
@inline Parameters(assignment::Assignment) = assignment.parameters

"""
    valtype(assign::Assignment)
    valtype(::Type{<:Assignment})

Type of the data (result) of an assignment.
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
@inline Parameters(algorithm::Algorithm) = algorithm.parameters

"""
    Algorithm(name::Symbol, frontend::Frontend; din::String=".", dout::String=".", parameters::Parameters=Parameters(frontend), map::Function=identity)

Construct an algorithm.
"""
@inline function Algorithm(name::Symbol, frontend::Frontend; din::String=".", dout::String=".", parameters::Parameters=Parameters(frontend), map::Function=identity, timer::TimerOutput=TimerOutput())
    return Algorithm(name, frontend, din, dout, parameters, map, Dict{Symbol, Assignment}(), timer)
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
    show(io::IO, alg::Algorithm)

Show an algorithm.

Optionally, some parameters of the algorithm can be filtered by specifying the `:select` context in `io`. Besides, the maximum number of decimals of the parameters can also be specified by the `:ndecimal` context in `io`.
"""
function Base.show(io::IO, alg::Algorithm)
    select = get(io, :select, param->true)
    ndecimal = get(io, :ndecimal, 10)
    @printf io "%s(%s)" alg.name alg.frontend
    flag = false
    for (name, value) in pairs(alg.parameters)
        if select(name)
            flag || @printf io "-"
            flag = true
            @printf io "%s(%s)" name tostr(value, ndecimal)
        end
    end
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
@inline Base.nameof(alg::Algorithm, assign::Assignment) = @sprintf "%s-%s" alg assign.id

"""
    add!(alg::Algorithm, id::Symbol, action::Action; parameters::Parameters=Parameters{()}(), map::Function=identity, dependencies::Tuple=(), kwargs...) -> Tuple{Algorithm, Assignment}

Add an assignment on an algorithm by providing the contents of the assignment without the execution of it.
"""
function add!(alg::Algorithm, id::Symbol, action::Action; parameters::Parameters=Parameters{()}(), map::Function=identity, dependencies::Tuple=(), kwargs...)
    assign = Assignment(id, action, merge(alg.parameters, parameters), map, dependencies, initialize(action, alg.frontend), false)
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

Run the dependencies of an assignment.

Optionally, some dependencies can be filtered by specifying the `f` function.
"""
function prepare!(alg::Algorithm, assign::Assignment, f::Function=assign->true)
    for id in assign.dependencies
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
