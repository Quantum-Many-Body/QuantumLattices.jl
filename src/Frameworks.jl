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
using ..Spatials: AbstractLattice, Bond, Neighbors, Point, bonds, isintracell, isparallel, minimumlengths, rcoordinate
using ..Toolkit: Float, atol, efficientoperations, indent, parametertype, reparameter, rtol

import ..QuantumLattices: add!, expand, expand!, reset!, str, update, update!
import ..QuantumOperators: scalartype
import ..Spatials: dlmsave
import ..Toolkit: contenttoshow, showasleaf, showcontent

export Action, Algorithm, Assignment, Boundary, CategorizedGenerator, Data, Eager, Embedding, ExpansionStyle, Formula, Frontend, Generator, LatticeModel, Lazy, OperatorGenerator, Parameters, ParametricGenerator, StaticGenerator
export checkoptions, config, contenttocache, contenttoconfig, datatype, eager, hasoption, lazy, options, optionsinfo, plain, qlcclean, qlclean, qlcsave, qldclean, qldsave, qlload, qlsave, run!, stamp

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

Convert a set of `Parameters` to a string of the form `"front-name₁(value₁)name₂(value₂)-rear"`.

Each parameter `nameᵢ(valueᵢ)` is included only when `select(nameᵢ)` returns `true`. Numeric values are rounded to `ndecimal` decimal places. If `front` or `rear` is empty, the leading/trailing `"-"` is omitted.
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
    Embedding{U<:AbstractLattice, L<:AbstractLattice, B<:Bond} <: LinearTransformation

Replicate a unitcell operator to every translation-equivalent copy inside a lattice.

# Fields
- `unitcell::U` — the unitcell lattice, which must carry non-empty translation vectors.
- `lattice::L` — the target lattice; every lattice bond must be translatable from a unitcell bond.
- `order::Int` — neighbor-shell scope (`0` = on-site only, `1` = NN, `2` = NNN, …).
- `lengths::Vector{Float64}` — cached neighbor lengths; `lengths[kind+1]` is the length of a `kind`th nearest neighbor of bond, computed by [`minimumlengths`](@ref).
- `mapping::Dict{B, Vector{Tuple{B, Int}}}` — precomputed lookup. Keys are unitcell bonds in forward orientation only (for 1-point bonds) or both forward and reverse orientations (for 2-point bonds). Values are lists of `(lattice_bond, direction)` where `direction = +1` (forward match) or `-1` (reverse match) as returned by [`isparallel`](@ref).
"""
struct Embedding{U<:AbstractLattice, L<:AbstractLattice, B<:Bond} <: LinearTransformation
    unitcell::U
    lattice::L
    order::Int
    lengths::Vector{Float64}
    mapping::Dict{B, Vector{Tuple{B, Int}}}
end

"""
    Embedding(unitcell::AbstractLattice, lattice::AbstractLattice, order::Int)

Construct an `Embedding` by precomputing the mapping from unitcell bonds to lattice bonds.

The construction computes neighbor lengths via [`minimumlengths`](@ref), enumerates all bonds up to `order` on the unitcell and lattice, and builds a mapping using [`isparallel`](@ref). Each lattice bond is iterated and matched to exactly one unitcell bond (or its reverse). For 1-point bonds the reverse entry is skipped since `reverse(bond) == bond`.

# Validation (all checked eagerly)
- `order ≥ 0`.
- `!isempty(unitcell.vectors)`.
- `dimension(unitcell) == dimension(lattice)`.
- Every lattice bond must match a unitcell bond; unmatched bonds raise an `AssertionError`.

# See also
- [`Embedding`](@ref) — the struct documentation.
- [`isparallel`](@ref) — the bond equivalence check reused from `Spatials`.
"""
function Embedding(unitcell::AbstractLattice, lattice::AbstractLattice, order::Int)
    @assert order ≥ 0 "Embedding error: order must be non-negative."
    @assert !isempty(unitcell.vectors) "Embedding error: unitcell must have non-empty vectors."
    @assert dimension(unitcell)==dimension(lattice) "Embedding error: mismatched space dimension."
    lengths = minimumlengths(unitcell.coordinates, unitcell.vectors, order)
    neighbors = Neighbors(lengths)
    refs = bonds(unitcell, neighbors)
    B = eltype(refs)
    mapping = Dict{B, Vector{Tuple{B, Int}}}()
    for bond in bonds(lattice, neighbors)
        matched = false
        for ref in refs
            dir = isparallel(ref, bond, unitcell.vectors, length(unitcell))
            if dir != 0
                push!(get!(()->Vector{Tuple{B, Int}}(), mapping, ref), (bond, dir))
                length(ref)==2 && push!(get!(()->Vector{Tuple{B, Int}}(), mapping, reverse(ref)), (bond, -dir))
                matched = true
                break
            end
        end
        @assert matched "Embedding error: lattice bond $bond matches no unitcell bond."
    end
    return Embedding(unitcell, lattice, order, lengths, mapping)
end

"""
    (em::Embedding)(m::Operator{<:Number, <:NTuple{2, CoordinatedIndex}}) -> OperatorSum

Apply the embedding to a rank-2 operator.

The bond that generated the operator is reconstructed from its `CoordinatedIndex` ids and the cached `lengths`: the spatial separation `‖r₂ − r₁‖` determines the neighbor order, and the appropriate 1-point or 2-point [`Bond`](@ref) is formed. On a cache hit in `mapping`, a copy is emitted at each matching lattice bond position:

- If `direction = +1`, the operator is placed at the lattice bond as-is.
- If `direction = -1`, the operator is placed at `reverse(bond)`.

The internal index and operator value are preserved from the unitcell operator.

# Errors
- `AssertionError` if the reconstructed bond is not found in `mapping` (operator was not generated from a unitcell bond), or if the bond length does not match any neighbor order in `lengths`.
"""
function (em::Embedding)(m::Operator{<:Number, <:NTuple{2, CoordinatedIndex}})
    len = norm(rcoordinate(m))
    neighbor = findfirst(l -> isapprox(l, len; atol=atol, rtol=rtol), em.lengths)
    isnothing(neighbor) && error("Embedding error: bond length $(len) does not match any neighbor order in lengths $(em.lengths).")
    ref = if neighbor == 1
        Bond(Point(m[1].index.site, m[1].rcoordinate, m[1].icoordinate))
    else
        Bond(
            neighbor - 1,
            Point(m[1].index.site, m[1].rcoordinate, m[1].icoordinate),
            Point(m[2].index.site, m[2].rcoordinate, m[2].icoordinate)
        )
    end
    @assert haskey(em.mapping, ref) "Embedding error: operator was not generated from a unitcell bond."
    result = zero(em, m)
    for (bond, dir) in em.mapping[ref]
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
    LatticeModel

Abstract supertype for all representations of a quantum lattice system.

Subtypes must implement `valtype`. `Parameters` and `update!` should also be implemented as applicable.
"""
abstract type LatticeModel end
@inline Base.:(==)(model₁::LatticeModel, model₂::LatticeModel) = ==(efficientoperations, model₁, model₂)
@inline Base.isequal(model₁::LatticeModel, model₂::LatticeModel) = isequal(efficientoperations, model₁, model₂)
@inline Base.show(io::IO, model::LatticeModel) = print(io, nameof(typeof(model)))
@inline showasleaf(::Type{<:LatticeModel}) = false
@inline Base.show(io::IO, ::MIME"text/plain", model::LatticeModel) = showcontent(io, model)

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
    Parameters(model::LatticeModel) -> NamedTuple

Get the parameters of a lattice model.

Returns `model.parameters` if the type has a `:parameters` field, otherwise an empty `NamedTuple`.
"""
@inline @generated Parameters(model::LatticeModel) = :parameters in fieldnames(model) ? :(model.parameters) : Parameters()

"""
    contenttoconfig(model::LatticeModel) -> Tuple

Return the structural components that characterize a lattice model.

Together with [`Parameters`](@ref), these components determine the model's identity.
Defaults to an empty `Tuple`. Subtypes should override this to specify which structural components to include.
"""
@inline contenttoconfig(model::LatticeModel) = ()

"""
    config(model::LatticeModel) -> String

Get the configuration fingerprint: the SHA-512 digest of [`contenttoconfig`](@ref).
"""
function config(model::LatticeModel)
    io = IOBuffer()
    serialize(io, contenttoconfig(model))
    return bytes2hex(sha512(take!(io)))
end

"""
    stamp(model::LatticeModel; ndecimal::Int=10, select::Function=name::Symbol->true, front::String="", rear::String="") -> String

Generate a deterministic stamp for a lattice model, formed as `"config[-parameters]"`.

The stamp serves as the internal key in data/cache files.
`ndecimal`, `select`, `front`, and `rear` are passed to `str(::Parameters)`.
When `str(::Parameters)` returns an empty string, the intermediate `"-"` is omitted.
"""
function stamp(model::LatticeModel; ndecimal::Int=10, select::Function=name::Symbol->true, front::String="", rear::String="")
    configuration = config(model)
    parameters = str(Parameters(model); ndecimal=ndecimal, select=select, front=front, rear=rear)
    return string(configuration, prepend(parameters, "-"))
end

"""
    contenttocache(model::LatticeModel) -> NamedTuple

Return the cacheable content of a lattice model as a `NamedTuple`.

Returns an empty `NamedTuple` by default.
Subtypes should override this to specify what should be cached.
"""
@inline contenttocache(model::LatticeModel) = NamedTuple()

"""
    dirname(model::LatticeModel) -> String

Get the dirname of the data/cache file of a lattice model.

Defaults to `"."` if the type has no `dir` field.
"""
@inline @generated Base.dirname(model::LatticeModel) = :dir in fieldnames(model) ? :(model.dir) : "."

"""
    basename(model::LatticeModel, target::Symbol=:void; prefix::String="", suffix::String="") -> String

Get the basename of a lattice model file, formed as `"prefix-string(model)-suffix.ext"`.

- `target::Symbol`: `:data` → `.qld`, `:cache` → `.qlc`, `:void` (default) → no extension.
- `prefix`, `suffix`: prepended/appended to the base name. The separator `"-"` is omitted when the corresponding part is empty.
"""
@inline function Base.basename(model::LatticeModel, target::Symbol=:void; prefix::String="", suffix::String="")
    ext = target==:data ? "qld" : target==:cache ? "qlc" : ""
    return string(append(prefix, "-"), string(model), prepend(suffix, "-"), isempty(ext) ? "" : prepend(ext, "."))
end

"""
    pathof(model::LatticeModel, target::Symbol; prefix::String="", suffix::String="") -> String

Get the full path of a lattice model file.

Equivalent to `joinpath(dirname(model), basename(model, target; prefix, suffix))`.
"""
@inline function Base.pathof(model::LatticeModel, target::Symbol; prefix::String="", suffix::String="")
    return joinpath(dirname(model), basename(model, target; prefix=prefix, suffix=suffix))
end

"""
    str(model::LatticeModel; prefix::String="", suffix::String="", ndecimal::Int=10, select::Function=name::Symbol->true, front::String="", rear::String="") -> String

Get the string representation of a lattice model, formed as `basename[-parameters]`.

- `prefix` and `suffix`: passed to `basename(::LatticeModel)`.
- `ndecimal`, `select`, `front`, `rear`: passed to `str(::Parameters)`.
"""
@inline function str(model::LatticeModel; prefix::String="", suffix::String="", ndecimal::Int=10, select::Function=name::Symbol->true, front::String="", rear::String="")
    base = basename(model; prefix=prefix, suffix=suffix)
    return string(base, prepend(str(Parameters(model); ndecimal=ndecimal, select=select, front=front, rear=rear), "-"))
end

"""
    qlsave(target::Symbol, model::LatticeModel, models::LatticeModel...; prefix::String="", suffix::String="", ndecimal::Int=10, select::Function=name::Symbol->true, front::String="", rear::String="")

Save lattice models to data (`.qld`) or cache (`.qlc`) files.

- `target::Symbol`: `:data` or `:cache`.
- `prefix`, `suffix`: passed to `pathof` → `basename`.
- `ndecimal`, `select`, `front`, `rear`: passed to `stamp` → `str(::Parameters)`.
"""
function qlsave(target::Symbol, model::LatticeModel, models::LatticeModel...; prefix::String="", suffix::String="", ndecimal::Int=10, select::Function=name::Symbol->true, front::String="", rear::String="")
    @assert target∈(:data, :cache) "qlsave error: target must be :data or :cache."
    for m in (model, models...)
        content = target==:data ? m : contenttocache(m)
        path = pathof(m, target; prefix=prefix, suffix=suffix)
        qlsave(path, stamp(m; ndecimal=ndecimal, select=select, front=front, rear=rear), content)
    end
end

"""
    qldsave(model::LatticeModel, models::LatticeModel...; kwargs...) -> String

Shortcut for `qlsave(:data, model, models...; kwargs...)`.
"""
@inline qldsave(model::LatticeModel, models::LatticeModel...; kwargs...) = qlsave(:data, model, models...; kwargs...)

"""
    qlcsave(model::LatticeModel, models::LatticeModel...; kwargs...) -> String

Shortcut for `qlsave(:cache, model, models...; kwargs...)`.
"""
@inline qlcsave(model::LatticeModel, models::LatticeModel...; kwargs...) = qlsave(:cache, model, models...; kwargs...)

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
    function OperatorGenerator(operators::CategorizedGenerator, bonds::Vector{<:Bond}, hilbert::Hilbert, terms::OneOrMore{Term}, half::Bool)
        terms = OneOrMore(terms)
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
    Assignment{A<:Action, P<:Parameters, M<:Function, T<:Tuple, D<:Data} <: LatticeModel

An assignment associated with an action.
"""
mutable struct Assignment{A<:Action, P<:Parameters, M<:Function, T<:Tuple, D<:Data} <: LatticeModel
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
@inline Base.valtype(::Type{<:Assignment{<:Action, <:Parameters, <:Function, <:Tuple, D}}) where {D<:Data} = D
@inline contenttoshow(assign::Assignment) = (; name=assign.name, action=assign.action, parameters=assign.parameters)

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
    dlmsave(assignment::Assignment, delim='\t'; prefix::String="", suffix::String="", ndecimal::Int=10, select::Function=name::Symbol->true, front::String="", rear::String="")

Save the data of an assignment to a delimited file.
"""
@inline function dlmsave(assignment::Assignment, delim='\t'; prefix::String="", suffix::String="", ndecimal::Int=10, select::Function=name::Symbol->true, front::String="", rear::String="")
    dlmsave(joinpath(dirname(assignment), string(str(assignment; prefix=prefix, suffix=suffix, ndecimal=ndecimal, select=select, front=front, rear=rear), ".dlm")), Tuple(assignment.data)..., delim)
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
@inline contenttoshow(alg::Algorithm) = (; name=alg.name, frontend=alg.frontend, parameters=alg.parameters)

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
