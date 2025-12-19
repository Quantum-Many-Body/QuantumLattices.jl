module Frameworks

using Base: @propagate_inbounds
using Base.Iterators: flatten, repeated
using DelimitedFiles: writedlm
using IndentWrappers: indent
using JLD2: jldopen, loadtodict!
using Latexify: latexify
using RecipesBase: RecipesBase, @recipe
using Serialization: deserialize, serialize
using TimerOutputs: TimerOutput, time, @timeit
using ..DegreesOfFreedom: plain, Boundary, Hilbert, Term
using ..QuantumLattices: OneOrMore, ZeroAtLeast, ZeroOrMore, id, value
using ..QuantumOperators: OperatorPack, Operators, OperatorSet, OperatorSum, LinearTransformation, identity, operatortype
using ..Spatials: Bond, isintracell
using ..Toolkit: atol, efficientoperations, rtol, parametertype

import ..QuantumLattices: add!, expand, expand!, reset!, str, update, update!
import ..QuantumOperators: scalartype
import ..Spatials: dlmsave

export Action, Algorithm, Assignment, CategorizedGenerator, Data, Eager, ExpansionStyle, Formula, Frontend, Generator, Lazy, OperatorGenerator, Parameters
export eager, lazy, checkoptions, datatype, hasoption, options, optionsinfo, qldload, qldsave, run!

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

Convert a set of `Parameters` to a string with each number hosting at most `ndecimal` decimal places. Here, the `select` function can select the key-value pairs to be contained by the keys.
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
    Parameters(ops::OperatorSet) -> NamedTuple{(), Tuple{}}

Get the parameters of an `OperatorSet`, which is defined to be an empty `NamedTuple`.
"""
@inline Parameters(ops::OperatorSet) = Parameters()

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
    function Formula{V}(expression::Function, parameters::Parameters) where V
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
@inline Base.show(io::IO, m::MIME"text/latex", gen::Generator) = show(io, m, latexify(expand(gen)))
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
    *(cat::CategorizedGenerator, factor::Number) -> CategorizedGenerator
    *(factor::Number, cat::CategorizedGenerator) -> CategorizedGenerator

Multiply a categorized generator of (representations of) quantum operators with a factor.
"""
@inline Base.:*(cat::CategorizedGenerator, factor::Number) = factor * cat
@inline function Base.:*(factor::Number, cat::CategorizedGenerator)
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
@inline function opsmatch(ops₁::NamedTuple, ops₂::NamedTuple, ::Val{key}) where key
    content = (get(ops₁, key, nothing), get(ops₂, key, nothing))
    return any(isnothing, content) || opsmatch(content...)
end
@inline opsmatch(ops₁, ops₂) = ops₁==ops₂ || zero(ops₁)==ops₁ || zero(ops₂)==ops₂
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
    OperatorGenerator{V<:Operators, CG<:CategorizedGenerator{V}, B<:Bond, H<:Hilbert, TS<:ZeroAtLeast{Term}} <: Generator{V}

A generator of operators based on the terms, bonds and Hilbert space of a quantum lattice system.
"""
struct OperatorGenerator{V<:Operators, CG<:CategorizedGenerator{V}, B<:Bond, H<:Hilbert, TS<:ZeroAtLeast{Term}} <: Generator{V}
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
@inline ExpansionStyle(::Type{<:OperatorGenerator{<:Operators, CG}}) where {CG<:CategorizedGenerator} = ExpansionStyle(CG)
@inline Base.valtype(::Type{<:OperatorGenerator{V}}) where {V<:Operators} = V
@inline expand(gen::OperatorGenerator, ::Lazy) = expand(gen.operators, lazy)

"""
    OperatorGenerator(bonds::Vector{<:Bond}, hilbert::Hilbert, terms::OneOrMore{Term}, boundary::Boundary=plain, style::ExpansionStyle=eager; half::Bool=false)

Construct a generator of quantum operators based on the input bonds, Hilbert space, terms and (twisted) boundary condition.

When the boundary condition is [`plain`](@ref), the boundary operators will be set to be empty for simplicity and efficiency.
"""
function OperatorGenerator(bonds::Vector{<:Bond}, hilbert::Hilbert, terms::OneOrMore{Term}, boundary::Boundary=plain, style::ExpansionStyle=eager; half::Bool=false)
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
    return OperatorGenerator(CategorizedGenerator(constops, alterops, boundops, parameters, boundary, style), bonds, hilbert, terms, half)
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
    Generator(operators::CategorizedGenerator{<:Operators}, bonds::Vector{<:Bond}, hilbert::Hilbert, terms::OneOrMore{Term}, half::Bool) -> OperatorGenerator
    Generator(bonds::Vector{<:Bond}, hilbert::Hilbert, terms::OneOrMore{Term}, boundary::Boundary=plain, style::ExpansionStyle=eager; half::Bool=false) -> OperatorGenerator

Construct an `OperatorGenerator`.
"""
@inline function Generator(operators::CategorizedGenerator{<:Operators}, bonds::Vector{<:Bond}, hilbert::Hilbert, terms::OneOrMore{Term}, half::Bool)
    return OperatorGenerator(operators, bonds, hilbert, terms, half)
end
@inline function Generator(bonds::Vector{<:Bond}, hilbert::Hilbert, terms::OneOrMore{Term}, boundary::Boundary=plain, style::ExpansionStyle=eager; half::Bool=false)
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
    Frontend

Frontend of algorithms applied to a quantum lattice system.
"""
abstract type Frontend end
@inline Base.:(==)(frontend₁::Frontend, frontend₂::Frontend) = ==(efficientoperations, frontend₁, frontend₂)
@inline Base.isequal(frontend₁::Frontend, frontend₂::Frontend) = isequal(efficientoperations, frontend₁, frontend₂)

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
    @recipe plot(assignment::Assignment)

Define the recipe for the visualization of an assignment of an algorithm.
"""
@recipe function plot(assignment::Assignment)
    title --> str(assignment)
    titlefontsize --> 10
    attr = seriestype(assignment.data)
    isnothing(attr) || begin
        seriestype --> attr
        attr==:path && begin
            legend --> false
            minorgrid --> true
            xminorticks --> 10
            yminorticks --> 10
        end
    end
    Tuple(assignment.data)
end
@inline seriestype(_...) = nothing
@inline seriestype(data::Data) = seriestype(Tuple(data)...)
@inline seriestype(::AbstractVector{<:Number}, ::Union{AbstractVector{<:Number}, AbstractMatrix{<:Number}}, _...) = :path
@inline seriestype(::AbstractVector{<:Number}, ::AbstractVector{<:Number}, ::Union{AbstractMatrix{<:Number}, AbstractArray{<:Number, 3}}, _...) = :heatmap

"""
    Algorithm{F<:Frontend, P<:Parameters, M<:Function} <: Function

An algorithm associated with an frontend.
"""
mutable struct Algorithm{F<:Frontend, P<:Parameters, M<:Function} <: Function
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
    @timeit alg.timer string(assign) begin
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
    @timeit alg.timer string(assign) begin
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

Add an assignment on a algorithm by providing the contents of the assignment, and run this assignment.
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
        @info "Assignment $name: time consumed $(time(alg.timer[string(assign)])/10^9)s."
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
    return string(append(prefix, "-"), string(obj), prepend(suffix, "-"), prepend(extension, "."))
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

Save arbitrary data to a qld file.
"""
function qldsave(filename::String, args...; mode::String="a+")
    @assert mode∈("a+", "w") "qldsave error: mode ($(repr(mode))) is not \"a+\" or \"w\"."
    @assert iseven(length(args)) "qldsave error: wrong formed input data."
    jldopen(filename, mode) do file
        io = IOBuffer()
        for i in 1:2:length(args)
            serialize(io, args[i+1])
            file[args[i]] = take!(io)
        end
    end
end

"""
    qldload(filename::String) -> Dict{String, Any}
    qldload(filename::String, name::String) -> Any
    qldload(filename::String, name₁::String, name₂::String, names::String...) -> Tuple

Load data from a qld file.
"""
function qldload(filename::String)
    result = jldopen(filename, "r") do file
        loadtodict!(Dict{String, Any}(), file)
    end
    for (key, bytes) in result
        result[key] = deserialize(IOBuffer(bytes))
    end
    return result
end
function qldload(filename::String, name::String)
    return jldopen(filename, "r") do file
        deserialize(IOBuffer(file[name]))
    end
end
function qldload(filename::String, name₁::String, name₂::String, names::String...)
    return jldopen(filename, "r") do file
        map((name₁, name₁, names...)) do name
            deserialize(IOBuffer(file[name]))
        end
    end
end

"""
    dlmsave(assignment::Assignment, delim='\t'; prefix::String="", suffix::String="", ndecimal::Int=10, select::Function=name::Symbol->true, front::String="", rear::String="")

Save the data of an assignment to a delimited file.
"""
@inline function dlmsave(assignment::Assignment, delim='\t'; prefix::String="", suffix::String="", ndecimal::Int=10, select::Function=name::Symbol->true, front::String="", rear::String="")
    dlmsave(joinpath(dirname(assignment), string(str(assignment; prefix=prefix, suffix=suffix, ndecimal=ndecimal, select=select, front=front, rear=rear), ".dlm")), Tuple(assignment.data)..., delim)
end

end  # module
