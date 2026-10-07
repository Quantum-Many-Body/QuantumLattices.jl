module Frameworks

using Base: @propagate_inbounds
using Base.Iterators: flatten, repeated
using HDF5: attrs, delete_object, h5open
using Latexify: latexify
using LinearAlgebra: norm
using Serialization: deserialize, serialize
using SHA: sha512
using StaticArrays: SVector
using TimerOutputs: @timeit, TimerOutput, time
using ..DegreesOfFreedom: CoordinatedIndex, Hilbert, Index, Term
using ..QuantumLattices: OneOrMore, ZeroAtLeast, ZeroOrMore, dimension, id, value
using ..QuantumOperators: LinearTransformation, Operator, OperatorPack, OperatorProd, OperatorSet, OperatorSum, Operators, identity, idtype, operatortype
using ..Spatials: AbstractLattice, Bond, Neighbors, Point, bonds, isintracell, isparallel, nneighbor, rcoordinate
using ..Toolkit: Float, atol, efficientoperations, parametertype, reparameter, rtol

import ..QuantumLattices: add!, expand, expand!, reset!, str, update, update!
import ..QuantumOperators: scalartype
import ..Spatials: Bond, dlmsave
import ..Toolkit: contenttoshow, showasleaf, showcontent

export Algorithm, Assignment, Boundary, CategorizedGenerator, Data, Eager, Embedding, ExpansionStyle, Formula, FrameworkElement, Generator, LatticeModel, Lazy, OperatorGenerator, Parameters, ParametricGenerator, StaticGenerator
export checkoptions, config, contenttocache, contenttoconfig, datatype, dependencytypes, eager, hasoption, lazy, options, optionsinfo, plain, qlcclean, qlclean, qlcsave, qldclean, qldsave, qlload, qlsave, run!, stamp, @delegate

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
    str(params::Parameters; ndecimal::Int=14, select::Function=name::Symbol->true, front::String="", rear::String="") -> String

Convert a set of `Parameters` to a string of the form `"front-name₁(value₁)name₂(value₂)-rear"`.

Each parameter `nameᵢ(valueᵢ)` is included only when `select(nameᵢ)` returns `true`. Numeric values are rounded to `ndecimal` decimal places. If `front` or `rear` is empty, the leading/trailing `"-"` is omitted.
"""
function str(params::Parameters; ndecimal::Int=14, select::Function=name::Symbol->true, front::String="", rear::String="")
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
    Boundary{Names}(values::AbstractVector{<:Number}, vectors::AbstractVector{<:AbstractVector{<:Number}}) where Names

Boundary twist of operators.
"""
struct Boundary{Names, D<:Number, V<:AbstractVector} <: LinearTransformation
    values::Vector{D}
    vectors::Vector{V}
    function Boundary{Names}(values::AbstractVector{<:Number}, vectors::AbstractVector{<:AbstractVector{<:Number}}) where Names
        @assert length(Names)==length(values)==length(vectors) "Boundary error: mismatched names, values and vectors."
        datatype = promote_type(eltype(values), Float)
        new{Names, datatype, eltype(vectors)}(convert(Vector{datatype}, values), vectors)
    end
end
@inline Base.:(==)(bound₁::Boundary, bound₂::Boundary) = keys(bound₁)==keys(bound₂) && ==(efficientoperations, bound₁, bound₂)
@inline Base.isequal(bound₁::Boundary, bound₂::Boundary) = isequal(keys(bound₁), keys(bound₂)) && isequal(efficientoperations, bound₁, bound₂)
@inline Base.valtype(::Type{<:Boundary}, M::Type{<:Operator}) = reparameter(M, :value, promote_type(Complex{Int}, scalartype(M)))
@inline Base.valtype(B::Type{<:Boundary}, MS::Type{<:Operators}) = (M = valtype(B, eltype(MS)); Operators{M, idtype(M)})
@inline contenttoshow(bound::Boundary) = (keys=keys(bound), values=bound.values, vectors=bound.vectors)

"""
    keys(bound::Boundary) -> ZeroAtLeast{Symbol}
    keys(::Type{<:Boundary{Names}}) where Names -> Names

Get the names of the boundary parameters.
"""
@inline Base.keys(bound::Boundary) = keys(typeof(bound))
@inline Base.keys(::Type{<:Boundary{Names}}) where Names = Names

"""
    (bound::Boundary)(operator::Operator; origin::Union{AbstractVector, Nothing}=nothing) -> Operator

Get the boundary twisted operator.
"""
@inline function (bound::Boundary)(operator::Operator; origin::Union{AbstractVector, Nothing}=nothing)
    values = isnothing(origin) ? bound.values : bound.values-origin
    return replace(operator, operator.value*exp(1im*mapreduce(u->angle(u, bound.vectors, values), +, id(operator))))
end

"""
    Parameters(bound::Boundary)

Get the parameters of the twisted boundary condition.
"""
@inline Parameters(bound::Boundary) = NamedTuple{keys(bound)}(ntuple(i->bound.values[i], Val(fieldcount(typeof(keys(bound))))))

"""
    update!(bound::Boundary; parameters...) -> Boundary

Update the values of the boundary twisted phase.
"""
@inline @generated function update!(bound::Boundary; parameters...)
    exprs = []
    for (i, name) in enumerate(QuoteNode.(keys(bound)))
        push!(exprs, :(bound.values[$i] = get(parameters, $name, bound.values[$i])))
    end
    return Expr(:block, exprs..., :(return bound))
end

"""
    merge!(bound::Boundary, another::Boundary) -> typeof(bound)

Merge the values and vectors of the twisted boundary condition from another one.
"""
@inline function Base.merge!(bound::Boundary, another::Boundary)
    @assert keys(bound)==keys(another) "merge! error: mismatched names of boundary parameters."
    bound.values .= another.values
    bound.vectors .= another.vectors
    return bound
end

"""
    reset!(bound::Boundary, values::AbstractVector{<:Number}) -> Boundary
    reset!(bound::Boundary, vectors::AbstractVector{<:AbstractVector{<:Number}}) -> Boundary
    reset!(bound::Boundary, values::AbstractVector{<:Number}, vectors::AbstractVector{<:AbstractVector{<:Number}}) -> Boundary

Reset the values or vectors of a twisted boundary condition in-place.

!!! note
    The plain boundary condition keeps plain even when reset with new values or new vectors.
"""
@inline function reset!(bound::Boundary, values::AbstractVector{<:Number})
    isempty(keys(bound)) && return bound
    bound.values .= values
    return bound
end
@inline function reset!(bound::Boundary, vectors::AbstractVector{<:AbstractVector{<:Number}})
    isempty(keys(bound)) && return bound
    bound.vectors .= vectors
    return bound
end
@inline function reset!(bound::Boundary, values::AbstractVector{<:Number}, vectors::AbstractVector{<:AbstractVector{<:Number}})
    isempty(keys(bound)) && return bound
    bound.values .= values
    bound.vectors .= vectors
    return bound
end

"""
    plain

Plain boundary condition without any twist.
"""
const plain = Boundary{()}(Float[], SVector{0, Float}[])
@inline Base.valtype(::Type{typeof(plain)}, M::Type{<:Operator}) = M
@inline Base.valtype(::Type{typeof(plain)}, M::Type{<:Operators}) = M
@inline (::typeof(plain))(operator::Operator; kwargs...) = operator
@inline showcontent(io::IO, ::Boundary{()}) = print(io, "plain")

# Embedding
"""
    Bond(m::Operator{<:Number, <:NTuple{2, CoordinatedIndex}}, neighbors::Neighbors; atol::Real=atol, rtol::Real=rtol) -> Bond

Reconstruct the unitcell bond from a rank-2 operator and a set of neighbors.

The bond length `‖r₂ − r₁‖` selects the neighbor order via [`isapprox`](@ref) against the neighbor lengths, and the result is a 1-point bond (on-site, order 0) or 2-point bond as appropriate.
"""
function Bond(m::Operator{<:Number, <:NTuple{2, CoordinatedIndex}}, neighbors::Neighbors; atol::Real=atol, rtol::Real=rtol)
    len = norm(rcoordinate(m))
    order = nothing
    for (i, neighbor) in neighbors
        if isapprox(len, neighbor; atol=atol, rtol=rtol)
            order = i
            break
        end
    end
    isnothing(order) && error("Embedding error: bond length $(len) does not match any neighbor order in neighbors $(neighbors).")
    if iszero(order)
        return Bond(Point(m[1].index.site, m[1].rcoordinate, m[1].icoordinate))
    else
        return Bond(
            Int(order),
            Point(m[1].index.site, m[1].rcoordinate, m[1].icoordinate),
            Point(m[2].index.site, m[2].rcoordinate, m[2].icoordinate)
        )
    end
end

"""
    Embedding{U<:AbstractLattice, N<:Neighbors, B<:Bond} <: LinearTransformation

Replicate a unitcell operator to every translation-equivalent copy inside a set of bonds.

# Fields
- `unitcell::U` — the unitcell, which must carry non-empty translation vectors.
- `neighbors::N` — neighbor order vs. bond length map, used for [`Bond`](@ref) reconstruction.
- `refs::Dict{Int, Vector{B}}` — unitcell reference bonds indexed by `bond.kind`, used for [`isparallel`](@ref)-based matching to find the correct mapping key for operators generated by external codes.
- `mapping::Dict{String, Vector{Tuple{B, Int}}}` — precomputed lookup. Keys are [`str`](@ref) string representations of unitcell bonds. Values are lists of `(lattice_bond, direction)` where `direction = +1` (forward) or `-1` (reverse) as returned by [`isparallel`](@ref).
"""
struct Embedding{U<:AbstractLattice, N<:Neighbors, B<:Bond} <: LinearTransformation
    unitcell::U
    neighbors::N
    refs::Dict{Int, Vector{B}}
    mapping::Dict{String, Vector{Tuple{B, Int}}}
end

"""
    Embedding(unitcell::AbstractLattice, bonds::AbstractVector{<:Bond}, neighbors::Neighbors; ndecimal::Int=14)

Construct an `Embedding` by precomputing the mapping from unitcell bonds to the target bonds.

`ndecimal` controls the number of decimal places for rounding coordinates in [`str`](@ref) mapping keys.

# Validation
- `!isempty(unitcell.vectors)`.
- Every target bond must match a unitcell bond; unmatched bonds raise an `AssertionError`.
"""
function Embedding(unitcell::AbstractLattice, targets::AbstractVector{<:Bond}, neighbors::Neighbors; ndecimal::Int=14)
    @assert !isempty(unitcell.vectors) "Embedding error: unitcell must have non-empty vectors."
    B = eltype(targets)
    refs = Dict{Int, Vector{B}}()
    mapping = Dict{String, Vector{Tuple{B, Int}}}()
    for ref in bonds(unitcell, neighbors)
        mapping[str(ref; ndecimal=ndecimal)] = Vector{Tuple{B, Int}}()
        push!(get!(()->Vector{B}(), refs, ref.kind), ref)
    end
    for bond in targets
        matched = false
        haskey(refs, bond.kind) && for ref in refs[bond.kind]
            dir = isparallel(ref, bond, unitcell.vectors, length(unitcell))
            if dir != 0
                push!(mapping[str(ref; ndecimal=ndecimal)], (bond, dir))
                matched = true
                break
            end
        end
        @assert matched "Embedding error: bond $bond matches no unitcell bond."
    end
    return Embedding(unitcell, neighbors, refs, mapping)
end

"""
    Embedding(unitcell::AbstractLattice, bonds::AbstractVector{<:Bond})

Construct an `Embedding` without precomputed neighbors.

The neighbors are derived from the bond lengths: on-site gets length 0, and for other bond the length is computed via `norm`([`rcoordinate`](@ref)(bond)).
"""
function Embedding(unitcell::AbstractLattice, bonds::AbstractVector{<:Bond})
    dict = Dict{Int, Float64}()
    for bond in bonds
        if !haskey(dict, bond.kind)
            dict[bond.kind] = length(bond)==2 ? norm(rcoordinate(bond)) : length(bond)==1 ? 0.0 : error("Embedding error: $bond not supported.")
        end
    end
    return Embedding(unitcell, bonds, Neighbors(dict))
end

"""
    Embedding(unitcell::AbstractLattice, lattice::AbstractLattice, neighbors::Neighbors)

Construct an `Embedding` from a lattice and precomputed neighbors, generating target bonds via [`bonds`](@ref)(lattice, neighbors).
"""
@inline function Embedding(unitcell::AbstractLattice, lattice::AbstractLattice, neighbors::Neighbors)
    @assert dimension(unitcell)==dimension(lattice) "Embedding error: mismatched space dimension."
    return Embedding(unitcell, bonds(lattice, neighbors), neighbors)
end

"""
    Embedding(unitcell::AbstractLattice, lattice::AbstractLattice, order::Int)

Construct an `Embedding` from a lattice and a neighbor order, computing neighbors via [`Neighbors`](@ref)(unitcell, order).
"""
function Embedding(unitcell::AbstractLattice, lattice::AbstractLattice, order::Int)
    @assert order ≥ 0 "Embedding error: order must be non-negative."
    return Embedding(unitcell, lattice, Neighbors(unitcell, order))
end

"""
    (em::Embedding)(m::Operator{<:Number, <:NTuple{2, CoordinatedIndex}}; atol::Real=atol, rtol::Real=rtol, ndecimal::Int=14) -> OperatorSum

Apply the embedding to a rank-2 operator.

The bond that generated the operator is reconstructed via [`Bond`](@ref)(m, neighbors), and [`isparallel`](@ref) looks up the matching unitcell reference bond in `refs`. The cached `mapping` then expands the operator to every translation-equivalent copy. Here, the combined direction determines whether the bond is emitted as-is (>0) or reversed (<0).
"""
function (em::Embedding)(m::Operator{<:Number, <:NTuple{2, CoordinatedIndex}}; atol::Real=atol, rtol::Real=rtol, ndecimal::Int=14)
    bond = Bond(m, em.neighbors; atol=atol, rtol=rtol)
    matched = nothing
    dir′ = 0
    haskey(em.refs, bond.kind) && for ref in em.refs[bond.kind]
        dir′ = isparallel(ref, bond, em.unitcell.vectors, length(em.unitcell))
        if dir′ != 0
            matched = ref
            break
        end
    end
    @assert !isnothing(matched) "Embedding error: operator was not generated from a unitcell bond."
    result = zero(em, m)
    for (bond, dir′′) in em.mapping[str(matched; ndecimal=ndecimal)]
        dir = dir′ * dir′′
        dir < 0 && (bond = reverse(bond))
        point₁ = bond[1]
        point₂ = length(bond) ≥ 2 ? bond[2] : bond[1]
        add!(result, Operator(
            m.value,
            CoordinatedIndex(Index(point₁.site, m[1].index.internal), point₁.rcoordinate, point₁.icoordinate),
            CoordinatedIndex(Index(point₂.site, m[2].index.internal), point₂.rcoordinate, point₂.icoordinate)
        ))
    end
    return result
end

"""
    valtype(::Type{<:Embedding}, M::Type{<:OperatorProd}) -> Type{<:OperatorSum}
    valtype(P::Type{<:Embedding}, M::Type{<:OperatorSet}) -> Type{<:OperatorSum}

Return the concrete `OperatorSum` type that `Embedding` produces when applied to an `OperatorProd` or `OperatorSet`.
"""
@inline Base.valtype(::Type{<:Embedding}, M::Type{<:OperatorProd}) = OperatorSum{M, idtype(M)}
@inline Base.valtype(P::Type{<:Embedding}, M::Type{<:OperatorSet}) = valtype(P, eltype(M))

"""
    FrameworkElement

Abstract supertype for all elements managed by the `Frameworks` submodule: lattice models, algorithms and assignments.

It provides a unified protocol: parameters ([`Parameters`](@ref)/[`update!`](@ref)), naming ([`str`](@ref)/[`basename`](@ref)/[`pathof`](@ref)), persistence ([`config`](@ref)/[`stamp`](@ref)/[`qlsave`](@ref)) and display, plus an optional `valtype` protocol (a type-level `valtype` automatically grants `scalartype`/`eltype`).
"""
abstract type FrameworkElement end
@inline Base.show(io::IO, element::FrameworkElement) = print(io, nameof(typeof(element)))
@inline showasleaf(::Type{<:FrameworkElement}) = false
@inline Base.show(io::IO, ::MIME"text/plain", element::FrameworkElement) = showcontent(io, element)

"""
    valtype(element::FrameworkElement)
    valtype(::Type{<:FrameworkElement})

Get the valtype of a framework element. The instance-level method forwards to the type-level one, which subtypes may implement.
"""
@inline Base.valtype(element::FrameworkElement) = valtype(typeof(element))

"""
    scalartype(element::FrameworkElement)
    scalartype(::Type{T}) where {T<:FrameworkElement}

Get the scalar type of a framework element, derived from its valtype.
"""
@inline scalartype(element::FrameworkElement) = scalartype(typeof(element))
@inline scalartype(::Type{T}) where {T<:FrameworkElement} = scalartype(valtype(T))

"""
    eltype(element::FrameworkElement)
    eltype(::Type{T}) where {T<:FrameworkElement}

Get the eltype of a framework element, derived from its valtype.
"""
@inline Base.eltype(element::FrameworkElement) = eltype(typeof(element))
@inline Base.eltype(::Type{T}) where {T<:FrameworkElement} = eltype(valtype(T))

"""
    Parameters(element::FrameworkElement) -> NamedTuple

Get the parameters of a framework element.

Returns `element.parameters` if the type has a `:parameters` field, otherwise an empty `NamedTuple`.
"""
@inline @generated Parameters(element::FrameworkElement) = :parameters in fieldnames(element) ? :(element.parameters) : Parameters()

"""
    contenttoconfig(element::FrameworkElement) -> Tuple

Return the structural components that characterize a framework element.

Together with [`Parameters`](@ref), these components determine the element's identity.
Defaults to an empty `Tuple`. Subtypes should override this to specify which structural components to include.
"""
@inline contenttoconfig(element::FrameworkElement) = ()

"""
    config(element::FrameworkElement) -> String

Get the configuration fingerprint: the SHA-512 digest of [`contenttoconfig`](@ref).
"""
function config(element::FrameworkElement)
    io = IOBuffer()
    serialize(io, contenttoconfig(element))
    return bytes2hex(sha512(take!(io)))
end

"""
    stamp(element::FrameworkElement; ndecimal::Int=14, select::Function=name::Symbol->true, front::String="", rear::String="") -> String

Generate a deterministic stamp for a framework element, formed as `"config[-parameters]"`.

The stamp serves as the internal key in data/cache files.
`ndecimal`, `select`, `front`, and `rear` are passed to `str(::Parameters)`.
When `str(::Parameters)` returns an empty string, the intermediate `"-"` is omitted.
"""
function stamp(element::FrameworkElement; ndecimal::Int=14, select::Function=name::Symbol->true, front::String="", rear::String="")
    configuration = config(element)
    parameters = str(Parameters(element); ndecimal=ndecimal, select=select, front=front, rear=rear)
    return string(configuration, prepend(parameters, "-"))
end

"""
    contenttocache(element::FrameworkElement) -> NamedTuple

Return the cacheable content of a framework element as a `NamedTuple`.

Returns an empty `NamedTuple` by default.
Subtypes should override this to specify what should be cached.
"""
@inline contenttocache(element::FrameworkElement) = NamedTuple()

"""
    dirname(element::FrameworkElement) -> String

Get the dirname of the data/cache file of a framework element.

Defaults to `"."` if the type has no `dir` field.
"""
@inline @generated Base.dirname(element::FrameworkElement) = :dir in fieldnames(element) ? :(element.dir) : "."

"""
    basename(element::FrameworkElement, target::Symbol=:void; prefix::String="", suffix::String="") -> String

Get the basename of a framework element file, formed as `"prefix-string(element)-suffix.ext"`.

- `target::Symbol`: `:data` → `.qld`, `:cache` → `.qlc`, `:void` (default) → no extension.
- `prefix`, `suffix`: prepended/appended to the base name. The separator `"-"` is omitted when the corresponding part is empty.
"""
@inline function Base.basename(element::FrameworkElement, target::Symbol=:void; prefix::String="", suffix::String="")
    ext = target==:data ? "qld" : target==:cache ? "qlc" : ""
    return string(append(prefix, "-"), string(element), prepend(suffix, "-"), isempty(ext) ? "" : prepend(ext, "."))
end

"""
    pathof(element::FrameworkElement, target::Symbol; prefix::String="", suffix::String="") -> String

Get the full path of a framework element file.

Equivalent to `joinpath(dirname(element), basename(element, target; prefix, suffix))`.
"""
@inline function Base.pathof(element::FrameworkElement, target::Symbol; prefix::String="", suffix::String="")
    return joinpath(dirname(element), basename(element, target; prefix=prefix, suffix=suffix))
end

"""
    str(element::FrameworkElement; prefix::String="", suffix::String="", ndecimal::Int=14, select::Function=name::Symbol->true, front::String="", rear::String="") -> String

Get the string representation of a framework element, formed as `basename[-parameters]`.

- `prefix` and `suffix`: passed to `basename(::FrameworkElement)`.
- `ndecimal`, `select`, `front`, `rear`: passed to `str(::Parameters)`.
"""
@inline function str(element::FrameworkElement; prefix::String="", suffix::String="", ndecimal::Int=14, select::Function=name::Symbol->true, front::String="", rear::String="")
    base = basename(element; prefix=prefix, suffix=suffix)
    return string(base, prepend(str(Parameters(element); ndecimal=ndecimal, select=select, front=front, rear=rear), "-"))
end

"""
    qlsave(target::Symbol, element::FrameworkElement, elements::FrameworkElement...; prefix::String="", suffix::String="", ndecimal::Int=14, select::Function=name::Symbol->true, front::String="", rear::String="")

Save framework elements to data (`.qld`) or cache (`.qlc`) files.

- `target::Symbol`: `:data` or `:cache`.
- `prefix`, `suffix`: passed to `pathof` → `basename`.
- `ndecimal`, `select`, `front`, `rear`: passed to `stamp` → `str(::Parameters)`.
"""
function qlsave(target::Symbol, element::FrameworkElement, elements::FrameworkElement...; prefix::String="", suffix::String="", ndecimal::Int=14, select::Function=name::Symbol->true, front::String="", rear::String="")
    @assert target∈(:data, :cache) "qlsave error: target must be :data or :cache."
    for m in (element, elements...)
        content = target==:data ? m : contenttocache(m)
        path = pathof(m, target; prefix=prefix, suffix=suffix)
        qlsave(path, stamp(m; ndecimal=ndecimal, select=select, front=front, rear=rear), content)
    end
end

"""
    qldsave(element::FrameworkElement, elements::FrameworkElement...; kwargs...) -> String

Shortcut for `qlsave(:data, element, elements...; kwargs...)`.
"""
@inline qldsave(element::FrameworkElement, elements::FrameworkElement...; kwargs...) = qlsave(:data, element, elements...; kwargs...)

"""
    qlcsave(element::FrameworkElement, elements::FrameworkElement...; kwargs...) -> String

Shortcut for `qlsave(:cache, element, elements...; kwargs...)`.
"""
@inline qlcsave(element::FrameworkElement, elements::FrameworkElement...; kwargs...) = qlsave(:cache, element, elements...; kwargs...)

"""
    LatticeModel <: FrameworkElement

Abstract supertype for all representations of a quantum lattice system.

Subtypes must implement the type-level `valtype`. `Parameters` and `update!` should also be implemented as applicable.
"""
abstract type LatticeModel <: FrameworkElement end
@inline Base.:(==)(model₁::LatticeModel, model₂::LatticeModel) = ==(efficientoperations, model₁, model₂)
@inline Base.isequal(model₁::LatticeModel, model₂::LatticeModel) = isequal(efficientoperations, model₁, model₂)

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
@inline contenttoconfig(formula::Formula) = (formula.expression,)

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
@inline contenttoconfig(gen::StaticGenerator) = (gen.operators,)
@inline expand(gen::StaticGenerator, ::Lazy) = gen.operators

"""
    (transformation::LinearTransformation)(gen::StaticGenerator; kwargs...) -> StaticGenerator

Apply a linear transformation to a static generator of (representations of) quantum operators.
"""
@inline (transformation::LinearTransformation)(gen::StaticGenerator; kwargs...) = StaticGenerator(transformation(gen.operators; kwargs...))

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
@inline expand(cat::CategorizedGenerator{<:Any, <:Any, @NamedTuple{}, @NamedTuple{}, @NamedTuple{}}, ::Lazy) = cat.constops
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
@propagate_inbounds function Base.iterate(ee::CategorizedGeneratorExpand)
    v = iterate(ee.values)
    isnothing(v) && return nothing
    op = iterate(ee.ops)
    return op[1]*v[1], (v[2], op[2])
end
@propagate_inbounds function Base.iterate(ee::CategorizedGeneratorExpand, state)
    v = iterate(ee.values, state[1])
    isnothing(v) && return nothing
    op = iterate(ee.ops, state[2])
    return op[1]*v[1], (v[2], op[2])
end
@inline contenttoconfig(cat::CategorizedGenerator) = (cat.constops, cat.alterops, cat.boundops, cat.boundary)

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
    function OperatorGenerator(operators::CategorizedGenerator, bonds::Vector{<:Bond}, hilbert::Hilbert, terms::ZeroOrMore{Term}, half::Bool)
        terms = ZeroOrMore(terms)
        new{valtype(operators), typeof(operators), eltype(bonds), typeof(hilbert), typeof(terms)}(operators, bonds, hilbert, terms, half)
    end
end
@inline contenttoconfig(gen::OperatorGenerator) = (contenttoconfig(gen.operators), gen.half)
@inline expand(gen::OperatorGenerator, ::Lazy) = expand(gen.operators, lazy)
@inline contenttoshow(gen::OperatorGenerator) = (;
    bonds = gen.bonds,
    hilbert = gen.hilbert,
    terms = map(id, gen.terms),
    half = gen.half,
    operators = gen.operators,
)

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
    OperatorGenerator(lattice::AbstractLattice, hilbert::Hilbert, terms::OneOrMore{Term}, boundary::Boundary=plain; half::Bool=false, neighbors::Union{Int, Neighbors}=nneighbor(terms))

Convenience constructor that takes a lattice instead of a pre-computed bond list.

The required bonds are determined automatically from the terms via [`nneighbor`](@ref).
"""
@inline function OperatorGenerator(lattice::AbstractLattice, hilbert::Hilbert, terms::OneOrMore{Term}, boundary::Boundary=plain; half::Bool=false, neighbors::Union{Int, Neighbors}=nneighbor(terms))
    return OperatorGenerator(bonds(lattice, neighbors), hilbert, terms, boundary; half=half)
end

"""
    OperatorGenerator(operators::OperatorSet, bonds::Vector{<:Bond}, hilbert::Hilbert, terms::ZeroOrMore{Term}; half::Bool=false)

Construct an operator generator combining pre-computed operators with term-based operators.

Here, `operators` are treated as static operators (e.g., from [`Embedding`](@ref)) while `terms` provide parameterized operators.

!!! note
    Only `boundary=plain` is supported because pre-computed operators are added wholesale to `constops` without distinguishing boundary-crossing operators — a distinction that would require a key (like `id(term)`) which isn't available.
"""
@inline function OperatorGenerator(operators::OperatorSet, bonds::Vector{<:Bond}, hilbert::Hilbert, ::Tuple{}; half::Bool=false)
    return OperatorGenerator(CategorizedGenerator(operators, NamedTuple(), NamedTuple(), NamedTuple(), plain), bonds, hilbert, (), half)
end
function OperatorGenerator(operators::OperatorSet, bonds::Vector{<:Bond}, hilbert::Hilbert, terms::OneOrMore{Term}; half::Bool=false)
    terms = ZeroOrMore(terms)
    constops = Operators{promote_type(eltype(operators), mapreduce(term->operatortype(eltype(bonds), typeof(hilbert), typeof(term)), promote_type, terms))}()
    add!(constops, operators)
    emptybonds = eltype(bonds)[]
    map(term->expand!(constops, term, term.ismodulatable ? emptybonds : bonds, hilbert; half=half), terms)
    alterops = NamedTuple{map(id, terms)}(expansion(terms, emptybonds, bonds, hilbert, scalartype(constops); half=half))
    boundops = NamedTuple{map(id, terms)}(expansion(terms, emptybonds, hilbert, plain, scalartype(constops); half=half))
    parameters = NamedTuple{map(id, terms)}(map(value, terms))
    cat = CategorizedGenerator(constops, alterops, boundops, parameters, plain)
    return OperatorGenerator(cat, bonds, hilbert, terms, half)
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
    expand(gen::OperatorGenerator, name::Symbol) -> Operators

Expand an operator generator to get the operators of a specific term.
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
    Generator(lattice::AbstractLattice, hilbert::Hilbert, terms::OneOrMore{Term}, boundary::Boundary=plain; half::Bool=false, neighbors::Union{Int, Neighbors}=nneighbor(terms)) -> OperatorGenerator
    Generator(operators::OperatorSet, bonds::Vector{<:Bond}, hilbert::Hilbert, terms::ZeroOrMore{Term}; half::Bool=false) -> OperatorGenerator

Factory constructor for `Generator` subtypes.

Dispatches on the argument types to construct the appropriate concrete representation.
Unlike [`LatticeModel`](@ref), this factory only constructs `Generator` subtypes (not `Formula`).
"""
@inline Generator(ops::OperatorSet) = StaticGenerator(ops)
@inline Generator(constops, alterops::NamedTuple, boundops::NamedTuple, parameters::Parameters, boundary::Boundary) = CategorizedGenerator(constops, alterops, boundops, parameters, boundary)
@inline Generator(ops::CategorizedGenerator{<:Operators}, bonds::Vector{<:Bond}, hilbert::Hilbert, terms::OneOrMore{Term}, half::Bool) = OperatorGenerator(ops, bonds, hilbert, terms, half)
@inline Generator(bonds::Vector{<:Bond}, hilbert::Hilbert, terms::OneOrMore{Term}, boundary::Boundary=plain; half::Bool=false) = OperatorGenerator(bonds, hilbert, terms, boundary; half=half)
@inline Generator(lattice::AbstractLattice, hilbert::Hilbert, terms::OneOrMore{Term}, boundary::Boundary=plain; half::Bool=false, neighbors::Union{Int, Neighbors}=nneighbor(terms)) = OperatorGenerator(lattice, hilbert, terms, boundary; half=half, neighbors=neighbors)
@inline Generator(operators::OperatorSet, bonds::Vector{<:Bond}, hilbert::Hilbert, terms::ZeroOrMore{Term}; half::Bool=false) = OperatorGenerator(operators, bonds, hilbert, terms; half=half)

"""
    LatticeModel(expression::Function, parameters::Parameters) -> Formula
    LatticeModel(ops::OperatorSet) -> StaticGenerator
    LatticeModel(constops, alterops::NamedTuple, boundops::NamedTuple, parameters::Parameters, boundary::Boundary) -> CategorizedGenerator
    LatticeModel(ops::CategorizedGenerator{<:Operators}, bonds::Vector{<:Bond}, hilbert::Hilbert, terms::OneOrMore{Term}, half::Bool) -> OperatorGenerator
    LatticeModel(bonds::Vector{<:Bond}, hilbert::Hilbert, terms::OneOrMore{Term}, boundary::Boundary=plain; half::Bool=false) -> OperatorGenerator
    LatticeModel(lattice::AbstractLattice, hilbert::Hilbert, terms::OneOrMore{Term}, boundary::Boundary=plain; half::Bool=false, neighbors::Union{Int, Neighbors}=nneighbor(terms)) -> OperatorGenerator
    LatticeModel(operators::OperatorSet, bonds::Vector{<:Bond}, hilbert::Hilbert, terms::ZeroOrMore{Term}; half::Bool=false) -> OperatorGenerator

Unified factory constructor for `LatticeModel` subtypes.

Dispatches on the argument types to construct the appropriate concrete representation.
`Algorithm` and `Assignment` are not constructed via this factory.
"""
@inline LatticeModel(expression::Function, parameters::Parameters) = Formula(expression, parameters)
@inline LatticeModel(ops::OperatorSet) = StaticGenerator(ops)
@inline LatticeModel(constops, alterops::NamedTuple, boundops::NamedTuple, parameters::Parameters, boundary::Boundary) = CategorizedGenerator(constops, alterops, boundops, parameters, boundary)
@inline LatticeModel(ops::CategorizedGenerator{<:Operators}, bonds::Vector{<:Bond}, hilbert::Hilbert, terms::OneOrMore{Term}, half::Bool) = OperatorGenerator(ops, bonds, hilbert, terms, half)
@inline LatticeModel(bonds::Vector{<:Bond}, hilbert::Hilbert, terms::OneOrMore{Term}, boundary::Boundary=plain; half::Bool=false) = OperatorGenerator(bonds, hilbert, terms, boundary; half=half)
@inline LatticeModel(lattice::AbstractLattice, hilbert::Hilbert, terms::OneOrMore{Term}, boundary::Boundary=plain; half::Bool=false, neighbors::Union{Int, Neighbors}=nneighbor(terms)) = OperatorGenerator(lattice, hilbert, terms, boundary; half=half, neighbors=neighbors)
@inline LatticeModel(operators::OperatorSet, bonds::Vector{<:Bond}, hilbert::Hilbert, terms::ZeroOrMore{Term}; half::Bool=false) = OperatorGenerator(operators, bonds, hilbert, terms; half=half)

"""
    Data

Abstract type for the data of a task.
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
    Assignment{A, P<:Parameters, T<:Tuple, D<:Data} <: FrameworkElement

An assignment associated with a computation task (task for short).

The fields fall into four groups: identity (`name`), storage (`dir`, only used by `pathof`/persistence), computation (`task`, `parameters`, `dependencies`) and result (`data`).
"""
mutable struct Assignment{A, P<:Parameters, T<:Tuple, D<:Data} <: FrameworkElement
    const dir::String
    const name::Symbol
    const task::A
    parameters::P
    const dependencies::T
    data::D
    function Assignment(::Type{D}, dir::String, name::Symbol, task, parameters::Parameters, dependencies::ZeroAtLeast{Assignment}) where {D<:Data}
        new{typeof(task), typeof(parameters), typeof(dependencies), D}(dir, name, task, parameters, dependencies)
    end
end
@inline contenttoshow(assign::Assignment) = (; name=assign.name, task=assign.task, parameters=assign.parameters)

"""
    ==(assign₁::Assignment, assign₂::Assignment) -> Bool
    isequal(assign₁::Assignment, assign₂::Assignment) -> Bool

Judge whether two assignments are equivalent.

The compared fields are `(name, task, parameters, dependencies)` and `data`; `dir` is excluded. The `task` fields are compared fieldwise via `efficientoperations`, and an undefined `data` only matches another undefined `data`.
"""
@inline function Base.:(==)(assign₁::Assignment, assign₂::Assignment)
    ==((assign₁.name, assign₁.parameters, assign₁.dependencies), (assign₂.name, assign₂.parameters, assign₂.dependencies)) || return false
    ==(efficientoperations, assign₁.task, assign₂.task) || return false
    d₁, d₂ = isdefined(assign₁, :data), isdefined(assign₂, :data)
    d₁ == d₂ || return false
    !d₁ && return true
    return assign₁.data == assign₂.data
end
@inline function Base.isequal(assign₁::Assignment, assign₂::Assignment)
    isequal((assign₁.name, assign₁.parameters, assign₁.dependencies), (assign₂.name, assign₂.parameters, assign₂.dependencies)) || return false
    isequal(efficientoperations, assign₁.task, assign₂.task) || return false
    d₁, d₂ = isdefined(assign₁, :data), isdefined(assign₂, :data)
    d₁ == d₂ || return false
    !d₁ && return true
    return isequal(assign₁.data, assign₂.data)
end

"""
    update!(assign::Assignment; parameters...) -> Assignment

Update the parameters of an assignment.
"""
function update!(assign::Assignment; parameters...)
    length(parameters)>0 && (assign.parameters = update(assign.parameters; parameters...))
    return assign
end

"""
    show(io::IO, assign::Assignment)

Show an assignment.
"""
@inline Base.show(io::IO, assign::Assignment) = print(io, assign.name)

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
    print(io, "Assignment{<:$(nameof(fieldtype(A, :task)))} options:\n")
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
    dlmsave(assignment::Assignment, delim='\t'; prefix::String="", suffix::String="", ndecimal::Int=14, select::Function=name::Symbol->true, front::String="", rear::String="")

Save the data of an assignment to a delimited file.
"""
@inline function dlmsave(assignment::Assignment, delim='\t'; prefix::String="", suffix::String="", ndecimal::Int=14, select::Function=name::Symbol->true, front::String="", rear::String="")
    dlmsave(joinpath(dirname(assignment), string(str(assignment; prefix=prefix, suffix=suffix, ndecimal=ndecimal, select=select, front=front, rear=rear), ".dlm")), Tuple(assignment.data)..., delim)
end

"""
    Algorithm{L<:LatticeModel, P<:Parameters, M<:Function} <: FrameworkElement

An algorithm associated with a frontend.

An algorithm is a transparent proxy of its `frontend`: any pure query interface `f` of the frontend satisfies `f(alg, ...) ≡ f(alg.frontend, ...)` (see [`@delegate`](@ref)), with only the execution context (e.g. the `timer`) injected.

The fields fall into four groups: identity (`name`), storage (`dir`, only used by `pathof`/persistence), computation (`frontend`, `parameters`, `map`) and observation (`timer`).
"""
mutable struct Algorithm{L<:LatticeModel, P<:Parameters, M<:Function} <: FrameworkElement
    const dir::String
    const name::Symbol
    const frontend::L
    parameters::P
    const map::M
    const timer::TimerOutput
end

"""
    ==(algorithm₁::Algorithm, algorithm₂::Algorithm) -> Bool
    isequal(algorithm₁::Algorithm, algorithm₂::Algorithm) -> Bool

Judge whether two algorithms are equivalent.

The compared fields are `(name, frontend, parameters, map)`; `dir` (storage) and `timer` (observation) are excluded.
"""
@inline function Base.:(==)(algorithm₁::Algorithm, algorithm₂::Algorithm)
    return ==(
        (algorithm₁.name, algorithm₁.frontend, algorithm₁.parameters, algorithm₁.map),
        (algorithm₂.name, algorithm₂.frontend, algorithm₂.parameters, algorithm₂.map),
    )
end
@inline function Base.isequal(algorithm₁::Algorithm, algorithm₂::Algorithm)
    return isequal(
        (algorithm₁.name, algorithm₁.frontend, algorithm₁.parameters, algorithm₁.map),
        (algorithm₂.name, algorithm₂.frontend, algorithm₂.parameters, algorithm₂.map),
    )
end
@inline Base.valtype(::Type{<:Algorithm{L}}) where {L<:LatticeModel} = valtype(L)
@inline contenttoshow(alg::Algorithm) = (; name=alg.name, frontend=alg.frontend, parameters=alg.parameters)

"""
    Algorithm(name::Symbol, frontend::LatticeModel, parameters::Parameters=Parameters(frontend), map::Function=identity; dir::String=".", timer::TimerOutput=TimerOutput())

Construct an algorithm.
"""
@inline function Algorithm(name::Symbol, frontend::LatticeModel, parameters::Parameters=Parameters(frontend), map::Function=identity; dir::String=".", timer::TimerOutput=TimerOutput())
    return Algorithm(dir, name, frontend, parameters, map, timer)
end

"""
    @delegate function f(x::X, args...; kwargs...) ... end
    @delegate inject=(s₁, ...) function f(x::X, args...; kwargs...) ... end
    @delegate function f(::Type{<:X}, args...; kwargs...) ... end
    @delegate @wrapper f(x::X, args...; kwargs...) = ...

Annotate a pure query interface of a `LatticeModel` subtype `X` so that it is transparently delegated to [`Algorithm`](@ref).

The expansion is the original method plus a forwarding method:
- instance level (`x::X`): `f(x::Algorithm{<:X}, args...; kwargs...) = f(x.frontend, args...; kwargs...)`;
- type level (`::Type{<:X}`): `f(::Type{<:Algorithm{L}}, args...; kwargs...) where {L<:X} = f(L, args...; kwargs...)`.

The forwarding method keeps the original argument names (with type annotations), default values and definition form (long/short). Wrapping macros (e.g. `@inline`) are copied verbatim to the forwarding method; `@generated` is not supported.

With `inject=(s₁, ...)`, each listed symbol must be a keyword argument of the original signature (checked at macro expansion time; type annotations such as `timer::TimerOutput=tbatimer` are recognized): in the forwarding signature its default is replaced by the same-named field of the algorithm (`x.s₁`), and it is passed explicitly in the forwarding call. Injection is only available at the instance level.

# Boundaries
- Only annotate pure query interfaces; never annotate `update!`/`Parameters`.
- `@delegate` must be the outermost macro of the definition (wrapping macros such as `@inline` go inside it). A docstring can be attached directly; it binds to the original method.
"""
macro delegate(args...)
    isempty(args) && error("@delegate error: nothing to delegate.")
    inject = Symbol[]
    for arg in args[1:end-1]
        (arg isa Expr && arg.head == :(=) && arg.args[1] == :inject) || error("@delegate error: unsupported option `$arg`; only `inject=(...)` is allowed.")
        entries = arg.args[2]
        entries = entries isa Expr && entries.head == :tuple ? entries.args : Any[entries]
        for entry in entries
            entry isa QuoteNode && (entry = entry.value)
            entry isa Symbol || error("@delegate error: inject entries must be symbols, got `$entry`.")
            push!(inject, entry)
        end
    end
    ex = args[end]

    wrappers = Expr[]
    while ex isa Expr && ex.head == :macrocall
        name = ex.args[1]
        macname = name isa Symbol ? name : name isa GlobalRef ? name.name : nothing
        macname == Symbol("@generated") && error("@delegate error: `@generated` functions are not supported.")
        push!(wrappers, ex)
        ex = ex.args[end]
    end

    whereparams = Any[]
    if ex isa Expr && ex.head == :where
        append!(whereparams, ex.args[2:end])
        ex = ex.args[1]
    end
    (ex isa Expr && (ex.head == :function || ex.head == :(=))) || error("@delegate error: expected a function definition, got `$ex`.")
    long = ex.head == :function
    sig = ex.args[1]
    if sig isa Expr && sig.head == :where
        append!(whereparams, sig.args[2:end])
        sig = sig.args[1]
    end
    (sig isa Expr && sig.head == :call) || error("@delegate error: expected a call signature, got `$sig`.")

    fname = sig.args[1]
    posargs = sig.args[2:end]
    kwargs = Any[]
    if !isempty(posargs) && posargs[1] isa Expr && posargs[1].head == :parameters
        kwargs = posargs[1].args
        posargs = posargs[2:end]
    end
    isempty(posargs) && error("@delegate error: the signature must have at least one positional argument.")

    firstarg = posargs[1]
    (firstarg isa Expr && firstarg.head == :(::)) || error("@delegate error: the first positional argument must be annotated, e.g. `x::X` or `::Type{<:X}`.")
    if length(firstarg.args) == 1
        ann = firstarg.args[1]
        (ann isa Expr && ann.head == :curly && ann.args[1] == :Type) || error("@delegate error: an unnamed first positional argument is only supported in the type-level form `::Type{<:X}`.")
        inner = ann.args[2]
        X = inner isa Expr && inner.head == :<: ? inner.args[end] : inner
        typelevel = true
        xname = nothing
    else
        xname = firstarg.args[1]
        xname isa Symbol || error("@delegate error: the first positional argument must be annotated, e.g. `x::X` or `::Type{<:X}`.")
        ann = firstarg.args[2]
        X = ann isa Expr && ann.head == :<: ? ann.args[end] : ann
        typelevel = false
    end
    kwargname(kw::Symbol) = kw
    kwargname(kw::Expr) = kw.head == :kw ? kwargname(kw.args[1]) : kw.head == :(::) && !isempty(kw.args) && kw.args[1] isa Symbol ? kw.args[1] : nothing
    kwargname(kw) = nothing
    (typelevel && !isempty(inject)) && error("@delegate error: `inject` is only available at the instance level.")
    for entry in inject
        found = any(kw -> kwargname(kw) === entry, kwargs)
        found || error("@delegate error: `:$entry` declared in `inject` is not a keyword argument of the original signature.")
    end

    algtype = GlobalRef(Frameworks, :Algorithm)
    fsigpos, callpos = Any[], Any[]
    if typelevel
        L = :L in whereparams ? gensym(:L) : :L
        push!(whereparams, Expr(:<:, L, X))
        push!(fsigpos, Expr(:(::), Expr(:curly, :Type, Expr(:<:, Expr(:curly, algtype, Expr(:<:, L))))))
        push!(callpos, L)
    else
        push!(fsigpos, Expr(:(::), xname, Expr(:curly, algtype, Expr(:<:, X))))
        push!(callpos, Expr(:., xname, QuoteNode(:frontend)))
    end
    for (i, arg) in enumerate(posargs[2:end])
        if arg isa Symbol
            push!(fsigpos, arg)
            push!(callpos, arg)
        elseif arg isa Expr && arg.head == :(::) && length(arg.args) == 2
            push!(fsigpos, arg)
            push!(callpos, arg.args[1])
        elseif arg isa Expr && arg.head == :(::)
            name = Symbol("_arg", i)
            push!(fsigpos, Expr(:(::), name, arg.args[1]))
            push!(callpos, name)
        elseif arg isa Expr && arg.head == :kw
            inner = arg.args[1]
            if inner isa Symbol
                push!(fsigpos, arg)
                push!(callpos, inner)
            elseif inner isa Expr && inner.head == :(::) && length(inner.args) == 2
                push!(fsigpos, arg)
                push!(callpos, inner.args[1])
            elseif inner isa Expr && inner.head == :(::)
                name = Symbol("_arg", i)
                push!(fsigpos, Expr(:kw, Expr(:(::), name, inner.args[1]), arg.args[2]))
                push!(callpos, name)
            else
                error("@delegate error: unsupported positional argument `$arg`.")
            end
        elseif arg isa Expr && arg.head == :...
            inner = arg.args[1]
            if inner isa Symbol
                push!(fsigpos, arg)
                push!(callpos, Expr(:..., inner))
            elseif inner isa Expr && inner.head == :(::) && length(inner.args) == 2
                push!(fsigpos, arg)
                push!(callpos, Expr(:..., inner.args[1]))
            elseif inner isa Expr && inner.head == :(::)
                name = Symbol("_arg", i)
                push!(fsigpos, Expr(:..., Expr(:(::), name, inner.args[1])))
                push!(callpos, Expr(:..., name))
            else
                error("@delegate error: unsupported positional argument `$arg`.")
            end
        else
            error("@delegate error: unsupported positional argument `$arg`.")
        end
    end
    fsigkw, callkw = Any[], Any[]
    for kw in kwargs
        if kw isa Symbol
            push!(fsigkw, kw)
            push!(callkw, Expr(:kw, kw, kw))
        elseif kw isa Expr && kw.head == :kw
            inner = kw.args[1]
            if inner isa Symbol
                name, siginner = inner, inner
            elseif inner isa Expr && inner.head == :(::) && length(inner.args) == 2 && inner.args[1] isa Symbol
                name, siginner = inner.args[1], inner
            else
                error("@delegate error: unsupported keyword argument `$kw`.")
            end
            default = name in inject ? Expr(:., xname, QuoteNode(name)) : kw.args[2]
            push!(fsigkw, Expr(:kw, siginner, default))
            push!(callkw, Expr(:kw, name, name))
        elseif kw isa Expr && kw.head == :(::) && length(kw.args) == 2 && kw.args[1] isa Symbol
            push!(fsigkw, kw)
            push!(callkw, Expr(:kw, kw.args[1], kw.args[1]))
        elseif kw isa Expr && kw.head == :...
            push!(fsigkw, kw)
            push!(callkw, kw)
        else
            error("@delegate error: unsupported keyword argument `$kw`.")
        end
    end

    fsig = Expr(:call, fname)
    isempty(fsigkw) || push!(fsig.args, Expr(:parameters, fsigkw...))
    append!(fsig.args, fsigpos)
    fcall = Expr(:call, fname)
    isempty(callkw) || push!(fcall.args, Expr(:parameters, callkw...))
    append!(fcall.args, callpos)
    fbody = Expr(:block, __source__, fcall)

    function hasname(s::Symbol, ex)::Bool
        ex === s && return true
        ex isa QuoteNode && ex.value === s && return true
        ex isa Expr && any(arg->hasname(s, arg), ex.args) && return true
        return false
    end
    kept = Any[]
    while true
        newkept = Any[]
        for p in whereparams
            s = p isa Symbol ? p : p isa Expr && p.head == :<: ? p.args[1] : nothing
            (s === nothing || hasname(s, fsig) || any(q->hasname(s, q), kept)) && push!(newkept, p)
        end
        length(newkept) == length(kept) && break
        kept = newkept
    end
    whereparams = kept

    isempty(whereparams) || (fsig = Expr(:where, fsig, whereparams...))
    fdef = long ? Expr(:function, fsig, fbody) : Expr(:(=), fsig, fbody)
    for wrapper in reverse(wrappers)
        fdef = Expr(:macrocall, wrapper.args[1:end-1]..., fdef)
    end
    marked = Expr(:macrocall, GlobalRef(Core, Symbol("@__doc__")), __source__, args[end])
    return Expr(:block, __source__, esc(marked), esc(fdef))
end

"""
    datatype(::Type{A}, ::Type{F}) where {A, F<:LatticeModel}

Get the concrete subtype of `Data` produced by running a task of type `A` on a frontend of type `F`.

Explicit registration by defining a more specific method is optional; by default the data type is inferred from the return type of the corresponding `run!`.
"""
@inline function datatype(::Type{A}, ::Type{F}) where {A, F<:LatticeModel}
    D = Core.Compiler.return_type(run!, Tuple{Algorithm{F}, Assignment{A}})
    @assert isconcretetype(D) && D<:Data "datatype error: inference failed ($D) for task type $A on frontend type $F. Please define `datatype(::Type{$(nameof(A))}, ::Type{$(nameof(F))}) = YourData` explicitly, or check the type stability of the corresponding `run!`."
    return D
end

"""
    dependencytypes(::Type) -> Union{Nothing, Tuple}

Declare the expected dependency types of a task, checked when an assignment is registered on an algorithm.

Returns `nothing` by default, which means no validation. Opt in by defining e.g. `dependencytypes(::Type{<:MyTask}) = (SomeOtherTask,)`; the dependencies are then required to match in count and to satisfy `dep isa Assignment{<:Tᵢ}` elementwise. Only fixed-length dependency lists are supported.
"""
@inline dependencytypes(::Type) = nothing

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
    config(alg::Algorithm) -> String

Get the configuration label of an algorithm, inherited from its frontend.
"""
@inline config(alg::Algorithm) = config(alg.frontend)

"""
    summary(alg::Algorithm)

Provide a summary of an algorithm.
"""
function Base.summary(alg::Algorithm)
    @info "Summary of $(alg.name) with $(nameof(typeof(alg.frontend))) frontend:"
    @info string(alg.timer)
end

"""
    (alg::Algorithm)(assign::Assignment; checkoptions::Bool=true, options...) -> Assignment
    (assign::Assignment)(alg::Algorithm; checkoptions::Bool=true, options...) -> Assignment

Run an assignment based on an algorithm.

The object in parentheses is the authority on parameters: the first form uses the parameters of `assign` as the current parameters (the algorithm is updated to them), while the second uses those of `alg` (the assignment is updated to them). Dependencies are always recursed in the reverse form (`dependency(alg; ...)`), so that they inherit the parameter context of the call site.

Each assignment is a single-slot cache: the cached data is hit only when `assign.data` is defined and the stored parameters match the current context parameters within `(atol, rtol)` as judged by `match`; `task`/`dependencies` are const and never compared, options do not trigger recomputation, and dependency identities are not compared recursively (a hit `data` is a self-consistent snapshot regardless of the current state of the dependencies, while a miss is corrected by the local check of each node). As corollaries: ① parameters within the tolerance are regarded as identical; ② the only way to force recomputation is to construct a new assignment; ③ a shared dependency switching between different parameter contexts will thrash — duplicate the dependency instance if both must coexist.
"""
function (alg::Algorithm)(assign::Assignment; checkoptions::Bool=true, options...)
    @timeit alg.timer string(assign) begin
        checkoptions && Frameworks.checkoptions(typeof(assign); options...)
        ismatched = match(assign.parameters, alg.parameters)
        if !(isdefined(assign, :data) && ismatched)
            ismatched || update!(alg; assign.parameters...)
            map(dependency->dependency(alg; checkoptions=false, options...), assign.dependencies)
            @timeit alg.timer "run!" (assign.data = run!(alg, assign; options...))
        end
    end
    return assign
end
function (assign::Assignment)(alg::Algorithm; checkoptions::Bool=true, options...)
    @timeit alg.timer string(assign) begin
        checkoptions && Frameworks.checkoptions(typeof(assign); options...)
        ismatched = match(alg.parameters, assign.parameters)
        if !(isdefined(assign, :data) && ismatched)
            ismatched || update!(assign; alg.parameters...)
            map(dependency->dependency(alg; checkoptions=false, options...), assign.dependencies)
            @timeit alg.timer "run!" (assign.data = run!(alg, assign; options...))
        end
    end
    return assign
end

"""
    (alg::Algorithm)(name::Symbol, task, dependencies::ZeroOrMore{Assignment}; dir::String=alg.dir, delay::Bool=false, options...) -> Assignment
    (alg::Algorithm)(name::Symbol, task, parameters::Parameters=Parameters(), dependencies::ZeroOrMore{Assignment}=(); dir::String=alg.dir, delay::Bool=false, options...) -> Assignment

Construct an assignment based on an algorithm by providing the contents of the assignment, and run this assignment.

At construction time, the keys of the local `parameters` must be a subset of those of `alg.parameters`, and the dependencies must satisfy the [`dependencytypes`](@ref) declaration of the task (if any).
"""
@inline function (alg::Algorithm)(name::Symbol, task, dependencies::ZeroOrMore{Assignment}; dir::String=alg.dir, delay::Bool=false, options...)
    return alg(name, task, Parameters(), dependencies; delay=delay, dir=dir, options...)
end
@inline function (alg::Algorithm)(name::Symbol, task, parameters::Parameters=Parameters(), dependencies::ZeroOrMore{Assignment}=(); dir::String=alg.dir, delay::Bool=false, options...)
    unknown = setdiff(keys(parameters), keys(alg.parameters))
    @assert isempty(unknown) "registration error: unknown local parameter keys $unknown; the algorithm $(alg.name) only has parameters $(keys(alg.parameters))."
    dtypes = dependencytypes(typeof(task))
    if !isnothing(dtypes)
        deps = ZeroOrMore(dependencies)
        @assert length(deps)==length(dtypes) "registration error: task $(nameof(typeof(task))) expects $(length(dtypes)) dependencies of types $dtypes, but got $(length(deps)) dependencies of types $(map(typeof, deps))."
        for (i, (dep, T)) in enumerate(zip(deps, dtypes))
            @assert dep isa Assignment{<:T} "registration error: dependency $i of task $(nameof(typeof(task))) is expected to be an Assignment{<:$T}, but got $(typeof(dep))."
        end
    end
    assign = Assignment(datatype(typeof(task), typeof(alg.frontend)), dir, name, task, merge(alg.parameters, parameters), ZeroOrMore(dependencies))
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
    qlsave(filename::String, args...)

Save arbitrary data as key-value pairs to a file.

If the file does not exist, it is created. If a key already exists, it is replaced.
"""
function qlsave(filename::String, args...)
    @assert iseven(length(args)) "qlsave error: wrong formed input data."
    h5open(filename, "cw") do file
        for i in 1:2:length(args)
            key = string(args[i])
            haskey(file, key) && delete_object(file, key)
            io = IOBuffer()
            serialize(io, args[i+1])
            write(file, key, take!(io))
            attrs(file[key])["touchtime"] = Base.time()
        end
    end
end

"""
    qlload(filename::String) -> Dict{String, Any}
    qlload(filename::String, name::String) -> Any
    qlload(filename::String, name₁::String, name₂::String, names::String...) -> Tuple

Load data from a qld/qlc file.

With a single filename argument, returns all entries as a flat `Dict`.
With a key, returns that single entry.
With multiple keys, returns a `Tuple`.
"""
function qlload(filename::String)
    isfile(filename) || error("qlload error: file '$filename' not found.")
    result = Dict{String, Any}()
    h5open(filename, "r+") do file
        for name in keys(file)
            result[name] = load(file, name)
        end
    end
    return result
end
function qlload(filename::String, name::String)
    isfile(filename) || error("qlload error: file '$filename' not found.")
    return h5open(filename, "r+") do file
        load(file, name)
    end
end
function qlload(filename::String, name₁::String, name₂::String, names::String...)
    isfile(filename) || error("qlload error: file '$filename' not found.")
    return h5open(filename, "r+") do file
        map((name₁, name₂, names...)) do name
            load(file, name)
        end
    end
end
@inline function load(file, name)
    haskey(file, name) || error("qlload error: key '$name' not found.")
    dataset = file[name]
    attrs(dataset)["touchtime"] = Base.time()
    bytes = read(dataset)
    bytes isa Vector{UInt8} || error("qlload error: key '$name' is not a data entry.")
    return deserialize(IOBuffer(bytes))
end

"""
    qlclean(filename::String; age::Real=Inf, maxcount::Int=typemax(Int)) -> Int
    qlclean(target::Symbol, dir::String="."; age::Real=Inf, maxcount::Int=typemax(Int)) -> Int

Clean up entries in qld/qlc files.

- `filename`: clean a single file.
- `target=:data`/`:cache`: clean all `.qld`/`.qlc` files in `dir`.
- `age`: delete entries whose `touchtime` is older than `age` seconds (default `Inf`).
- `maxcount`: keep at most `maxcount` most recent entries per file (default `typemax(Int)`).

Files that become empty after cleaning are removed. Returns the number of entries deleted.
"""
function qlclean(filename::String; age::Real=Inf, maxcount::Int=typemax(Int))
    isfile(filename) || return 0
    deleted = 0
    nowtime = Base.time()
    entries = Pair{String, Float64}[]
    removable = false
    h5open(filename, "r+") do file
        for name in keys(file)
            touchtime = try attrs(file[name])["touchtime"] catch; 0.0 end
            push!(entries, name=>touchtime)
        end
        sort!(entries; by=pair->pair.second, rev=true)
        kept = 0
        for (name, touchtime) in entries
            if (nowtime-touchtime > age) || (kept >= maxcount)
                delete_object(file, name)
                deleted += 1
            else
                kept += 1
            end
        end
        removable = isempty(keys(file))
    end
    removable && rm(filename; force=true)
    return deleted
end
function qlclean(target::Symbol, dir::String="."; age::Real=Inf, maxcount::Int=typemax(Int))
    @assert target∈(:data, :cache) "qlclean error: target must be :data or :cache."
    extension = target==:data ? ".qld" : ".qlc"
    files = filter(file->endswith(file, extension), readdir(dir; join=true))
    return sum(qlclean(file; age=age, maxcount=maxcount) for file in files)
end

"""
    qldclean(dir::String="."; kwargs...) -> Int

Shortcut for `qlclean(:data, dir; kwargs...)`.
"""
@inline qldclean(dir::String="."; kwargs...) = qlclean(:data, dir; kwargs...)

"""
    qlcclean(dir::String="."; kwargs...) -> Int

Shortcut for `qlclean(:cache, dir; kwargs...)`.
"""
@inline qlcclean(dir::String="."; kwargs...) = qlclean(:cache, dir; kwargs...)

end  # module
