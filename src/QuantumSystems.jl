module QuantumSystems

using LinearAlgebra: dot, ishermitian, norm
using Printf: @printf
using SparseArrays: SparseMatrixCSC
using StaticArrays: SMatrix, SVector
using ..DegreesOfFreedom: CompositeIndex, CoordinatedIndex, Coupling, Hilbert, Index, InternalIndex, MatrixCouplingComponent, Ordinal, Pattern, SimpleInternal, Term, TermAmplitude, TermCoupling, indextype, @pattern
using ..QuantumLattices: OneAtLeast, ZeroAtLeast, decompose, rank, str
using ..QuantumOperators: LaTeX, Operator, OperatorIndex, OperatorProd, Operators, latexformat
using ..Spatials: Bond, Point, direction, isparallel, rcoordinate
using ..Toolkit: atol, efficientoperations, rtol, Float, VectorSpace, VectorSpaceDirectProducted, VectorSpaceStyle, delta, rawtype

import ..DegreesOfFreedom: MatrixCoupling, diagonalfields, internalindextype, isdefinite, patternrule, showablefields, statistics
import ..QuantumLattices: expand, expand!, kind, permute, shape
import ..QuantumOperators: latexname, matrix, script

# Canonical complex fermionic/bosonic systems
export σ⁰, σˣ, σʸ, σᶻ, σ⁺, σ⁻, σ¹¹, σ¹², σ²¹, σ²², annihilation, creation, latexofbosons, latexoffermions, latexofparticles, Lˣ, Lʸ, Lᶻ
export 𝕒, 𝕒⁺, 𝕒𝕒, 𝕒𝕒⁺, 𝕒⁺𝕒, 𝕒⁺𝕒⁺, 𝕔, 𝕔⁺, 𝕔𝕔, 𝕔𝕔⁺, 𝕔⁺𝕔, 𝕔⁺𝕔⁺, 𝕕, 𝕕⁺, 𝕕𝕕, 𝕕𝕕⁺, 𝕕⁺𝕕, 𝕕⁺𝕕⁺, Fock, FockIndex, isannihilation, iscreation, isnormalordered
export Coulomb, FockTerm, Hopping, Hubbard, InterOrbitalInterSpin, InterOrbitalIntraSpin, Onsite, PairHopping, Pairing, SpinFlip

# SU(2) spin systems
export Γˣ, Γʸ, Γᶻ, Γ′ˣ, Γ′ʸ, Γ′ᶻ, DMˣ, DMʸ, DMᶻ, Isingˣ, Isingʸ, Isingᶻ, latexofspins
export 𝕊, 𝕊ᵀ𝕊, SpinIndex, Spin, totalspin
export Γ, Γ′, DM, Heisenberg, Ising, Kitaev, SingleIonAnisotropy, SpinTerm, Zeeman

# Phononic systems
export latexofphonons
export 𝕦, 𝕦ᵀ𝕦, 𝕡, Phonon, PhononIndex
export Elastic, Kinetic, Hooke, PhononTerm

# Canonical complex fermionic/bosonic systems and hardcore bosonic systems
## FockIndex
"""
    annihilation

Indicate that the nambu index is annihilation.
"""
const annihilation = 1

"""
    creation

Indicate that the nambu index is creation.
"""
const creation = 2

"""
    FockIndex{T, O<:Union{Int, Symbol, Colon}, S<:Union{Rational{Int}, Symbol, Colon}} <: InternalIndex

Fock index, i.e., the internal index to specify the generators of a Fock space.
"""
struct FockIndex{T, O<:Union{Int, Symbol, Colon}, S<:Union{Rational{Int}, Symbol, Colon}} <: InternalIndex
    orbital::O
    spin::S
    nambu::Int
    function FockIndex{T}(orbital::Union{Int, Symbol, Colon}, spin::Union{Rational{Int}, Int, Symbol, Colon}, nambu::Int) where T
        @assert T∈(:f, :b, :) "FockIndex error: wrong statistics."
        isa(spin, Rational{Int}) && @assert spin.den∈(1, 2) "FockIndex error: wrong spin."
        @assert nambu∈(annihilation, creation) "FockIndex error: wrong input nambu($nambu)."
        isa(spin, Int) && (spin = convert(Rational{Int}, spin))
        new{T, typeof(orbital), typeof(spin)}(orbital, spin, nambu)
    end
end
### basic methods of concrete InternalIndex
@inline statistics(::Type{<:FockIndex}) = Symbol(":")
@inline statistics(::Type{<:FockIndex{T}}) where T = T
@inline isdefinite(::Type{<:FockIndex{T, Int, Rational{Int}} where T}) = true
@inline Base.:(==)(index₁::FockIndex, index₂::FockIndex) = statistics(index₁)==statistics(index₂) && ==(efficientoperations, index₁, index₂)
@inline Base.isequal(index₁::FockIndex, index₂::FockIndex) = isequal(statistics(index₁), statistics(index₂)) && isequal(efficientoperations, index₁, index₂)
@inline Base.hash(index::FockIndex, h::UInt) = hash((statistics(index), index.orbital, index.spin, index.nambu), h)
@inline Base.adjoint(index::FockIndex) = FockIndex{statistics(index)}(index.orbital, index.spin, 3-index.nambu)
@inline @generated function Base.replace(index::FockIndex; kwargs...)
    exprs = [:(get(kwargs, $name, getfield(index, $name))) for name in QuoteNode.(fieldnames(index))]
    return :(rawtype(typeof(index)){statistics(index)}($(exprs...)))
end
### requested by show
@inline showablefields(::Type{<:FockIndex}) = (:orbital, :spin)
### requested by Pattern
@inline diagonalfields(::Type{<:FockIndex}) = (:orbital, :spin)
### requested by MatrixCoupling
@inline internalindextype(::Type{FockIndex}, ::Type{O}, ::Type{S}, ::Type{Int}) where {O<:Union{Int, Symbol, Colon}, S<:Union{Rational{Int}, Symbol, Colon}} = FockIndex{:, O, S}
@inline internalindextype(::Type{FockIndex{T}}, ::Type{O}, ::Type{S}, ::Type{Int}) where {T, O<:Union{Int, Symbol, Colon}, S<:Union{Rational{Int}, Symbol, Colon}} = FockIndex{T, O, S}

"""
    FockIndex(orbital::Union{Int, Symbol, Colon}, spin::Union{Rational{Int}, Int, Symbol, Colon}, nambu::Int)
    FockIndex{T}(orbital::Union{Int, Symbol, Colon}, spin::Union{Rational{Int}, Int, Symbol, Colon}, nambu::Int) where T
    FockIndex{T, O, S}(orbital::Union{Int, Symbol, Colon}, spin::Union{Rational{Int}, Symbol, Colon}, nambu::Int) where {T, O, S}

Construct a Fock index.
"""
@inline FockIndex(orbital::Union{Int, Symbol, Colon}, spin::Union{Rational{Int}, Int, Symbol, Colon}, nambu::Int) = FockIndex{:}(orbital, spin, nambu)
@inline FockIndex{T, O, S}(orbital::Union{Int, Symbol, Colon}, spin::Union{Rational{Int}, Symbol, Colon}, nambu::Int) where {T, O, S} = FockIndex{T}(orbital, spin, nambu)

"""
    isannihilation(index::FockIndex) -> Bool
    isannihilation(index::Index) -> Bool
    isannihilation(index::CompositeIndex) -> Bool

Judge whether the nambu index is `annihilation`.
"""
@inline isannihilation(index::FockIndex) = index.nambu==annihilation
@inline isannihilation(index::Index) = isannihilation(InternalIndex(index))
@inline isannihilation(index::CompositeIndex) = isannihilation(Index(index))

"""
    iscreation(index::FockIndex) -> Bool
    iscreation(index::Index) -> Bool
    iscreation(index::CompositeIndex) -> Bool

Judge whether the nambu index is `creation`.
"""
@inline iscreation(index::FockIndex) = index.nambu==creation
@inline iscreation(index::Index) = iscreation(InternalIndex(index))
@inline iscreation(index::CompositeIndex) = iscreation(Index(index))

### convenient construction and string representation 
"""
    𝕔(orbital, spin) -> FockIndex{:f}
    𝕔(site, orbital, spin) -> Index{<:FockIndex{:f}}
    𝕔(site, orbital, spin, rcoordinate, icoordinate) -> CoordinatedIndex{<:Index{<:FockIndex{:f}}}

Convenient construction of `FockIndex{:f}`, `Index{<:FockIndex{:f}}`, `CoordinatedIndex{<:Index{<:FockIndex{f}}}` with the `nambu` attribute being `annihilation`.
"""
function 𝕔 end

"""
    𝕔⁺(orbital, spin) -> FockIndex{:f}
    𝕔⁺(site, orbital, spin) -> Index{<:FockIndex{:f}}
    𝕔⁺(site, orbital, spin, rcoordinate, icoordinate) -> CoordinatedIndex{<:Index{<:FockIndex{:f}}}

Convenient construction of `FockIndex{:f}`, `Index{<:FockIndex{:f}}`, `CoordinatedIndex{<:Index{<:FockIndex{f}}}` with the `nambu` attribute being `creation`.
"""
function 𝕔⁺ end

"""
    𝕒(orbital, spin) -> FockIndex{:b}
    𝕒(site, orbital, spin) -> Index{<:FockIndex{:b}}
    𝕒(site, orbital, spin, rcoordinate, icoordinate) -> CoordinatedIndex{<:Index{<:FockIndex{:b}}}

Convenient construction of `FockIndex{:b}`, `Index{<:FockIndex{:b}}`, `CoordinatedIndex{<:Index{<:FockIndex{:b}}}` with the `nambu` attribute being `annihilation`.
"""
function 𝕒 end

"""
    𝕒⁺(orbital, spin) -> FockIndex{:b}
    𝕒⁺(site, orbital, spin) -> Index{<:FockIndex{:b}}
    𝕒⁺(site, orbital, spin, rcoordinate, icoordinate) -> CoordinatedIndex{<:Index{<:FockIndex{:b}}}

Convenient construction of `FockIndex{:b}`, `Index{<:FockIndex{:b}}`, `CoordinatedIndex{<:Index{<:FockIndex{:b}}}` with the `nambu` attribute being `creation`.
"""
function 𝕒⁺ end

"""
    𝕕(orbital, spin) -> FockIndex{:}
    𝕕(site, orbital, spin) -> Index{<:FockIndex{:}}
    𝕕(site, orbital, spin, rcoordinate, icoordinate) -> CoordinatedIndex{<:Index{<:FockIndex{:}}}

Convenient construction of `FockIndex{:}`, `Index{<:FockIndex{:}}`, `CoordinatedIndex{<:Index{<:FockIndex{:}}}` with the `nambu` attribute being `annihilation`.
"""
function 𝕕 end

"""
    𝕕⁺(orbital, spin) -> FockIndex{:}
    𝕕⁺(site, orbital, spin) -> Index{<:FockIndex{:}}
    𝕕⁺(site, orbital, spin, rcoordinate, icoordinate) -> CoordinatedIndex{<:Index{<:FockIndex{:}}}

Convenient construction of `FockIndex{:}`, `Index{<:FockIndex{:}}`, `CoordinatedIndex{<:Index{<:FockIndex{:}}}` with the `nambu` attribute being `creation`.
"""
function 𝕕⁺ end

const _annihilation_ = (:𝕔, :𝕒, :𝕕), (QuoteNode(:f), QuoteNode(:b), :)
for (name, statistics) in zip(_annihilation_...)
    @eval @inline $name(orbital, spin) = FockIndex{$statistics}(orbital, spin, annihilation)
    @eval @inline $name(site, orbital, spin) = Index(site, FockIndex{$statistics}(orbital, spin, annihilation))
    @eval @inline $name(site, orbital, spin, rcoordinate, icoordinate) = CoordinatedIndex(Index(site, FockIndex{$statistics}(orbital, spin, annihilation)), rcoordinate, icoordinate)
end
const _creation_ = (:𝕔⁺, :𝕒⁺, :𝕕⁺), (QuoteNode(:f), QuoteNode(:b), :)
for (name, statistics) in zip(_creation_...)
    @eval @inline $name(orbital, spin) = FockIndex{$statistics}(orbital, spin, creation)
    @eval @inline $name(site, orbital, spin) = Index(site, FockIndex{$statistics}(orbital, spin, creation))
    @eval @inline $name(site, orbital, spin, rcoordinate, icoordinate) = CoordinatedIndex(Index(site, FockIndex{$statistics}(orbital, spin, creation)), rcoordinate, icoordinate)
end
@inline Base.getindex(::Type{OperatorIndex}, index::FockIndex{:f}) = isannihilation(index) ? "𝕔" : iscreation(index) ? "𝕔⁺" : error("wrong index.")
@inline Base.getindex(::Type{OperatorIndex}, index::FockIndex{:b}) = isannihilation(index) ? "𝕒" : iscreation(index) ? "𝕒⁺" : error("wrong index.")
@inline Base.getindex(::Type{OperatorIndex}, index::FockIndex{:}) = isannihilation(index) ? "𝕕" : iscreation(index) ? "𝕕⁺" : error("wrong index.")

### LaTeX format output
"""
    script(index::FockIndex, ::Val{:orbital}; kwargs...) -> String
    script(index::FockIndex, ::Val{:spin}; kwargs...) -> String
    script(index::FockIndex, ::Val{:spinsym}; kwargs...) -> String
    script(index::FockIndex, ::Val{:nambu}; kwargs...) -> String

Get the requested script of a Fock index.
"""
@inline script(index::FockIndex, ::Val{:orbital}; kwargs...) = str(index.orbital)
@inline script(index::FockIndex, ::Val{:spin}; kwargs...) = str(index.spin)
@inline function script(index::FockIndex, ::Val{:spinsym}; kwargs...)
    nspin = get(kwargs, :nspin, nothing)
    return (isnothing(nspin) || nspin==2) && index.spin==-1//2 ? "↓" : (isnothing(nspin) || nspin==2) && index.spin==1//2 ? "↑" : (isnothing(nspin) || nspin==1) && index.spin==0 ? "" : str(index.spin)
end
@inline script(index::FockIndex, ::Val{:nambu}; kwargs...) = iscreation(index) ? "\\dagger" : isannihilation(index) ? "" : str(index.nambu)

"""
    latexoffermions

Default LaTeX format of a fermionic index.
"""
const latexoffermions = LaTeX{(:nambu,), (:site, :orbital, :spinsym)}('c')
@inline latexname(::Type{<:FockIndex{:f}}) = Symbol("FockIndex{:f}")
@inline latexname(::Type{<:Index{<:FockIndex{:f}}}) = Symbol("Index{FockIndex{:f}}")
@inline latexname(::Type{<:CompositeIndex{<:Index{<:FockIndex{:f}}}}) = Symbol("CompositeIndex{Index{FockIndex{:f}}}")
latexformat(FockIndex{:f}, latexoffermions)
latexformat(Index{<:FockIndex{:f}}, latexoffermions)
latexformat(CompositeIndex{<:Index{<:FockIndex{:f}}}, latexoffermions)

"""
    latexofbosons

Default LaTeX format of a bosonic index.
"""
const latexofbosons = LaTeX{(:nambu,), (:site, :orbital, :spinsym)}('b')
@inline latexname(::Type{<:Index{<:FockIndex{:b}}}) = Symbol("Index{FockIndex{:b}}")
@inline latexname(::Type{<:CompositeIndex{<:Index{<:FockIndex{:b}}}}) = Symbol("CompositeIndex{Index{FockIndex{:b}}}")
@inline latexname(::Type{<:FockIndex{:b}}) = Symbol("FockIndex{:b}")
latexformat(FockIndex{:b}, latexofbosons)
latexformat(Index{<:FockIndex{:b}}, latexofbosons)
latexformat(CompositeIndex{<:Index{<:FockIndex{:b}}}, latexofbosons)

"""
    latexofparticles

Default LaTeX format of a wildcard Fock index.
"""
const latexofparticles = LaTeX{(:nambu,), (:site, :orbital, :spinsym)}('d')
@inline latexname(::Type{<:Index{<:FockIndex{:}}}) = Symbol("Index{FockIndex}")
@inline latexname(::Type{<:CompositeIndex{<:Index{<:FockIndex{:}}}}) = Symbol("CompositeIndex{Index{FockIndex}}")
@inline latexname(::Type{<:FockIndex{:}}) = Symbol("FockIndex")
latexformat(FockIndex{:}, latexofparticles)
latexformat(Index{<:FockIndex{:}}, latexofparticles)
latexformat(CompositeIndex{<:Index{<:FockIndex{:}}}, latexofparticles)

## Fock
"""
    Fock{T} <: SimpleInternal{FockIndex{T, Int, Rational{Int}}}

Fock space of Fock generators at a single point.
"""
struct Fock{T} <: SimpleInternal{FockIndex{T, Int, Rational{Int}}}
    norbital::Int
    nspin::Int
    function Fock{T}(norbital::Int, nspin::Int) where T
        @assert T∈(:f, :b) "Fock error: wrong statistics."
        new{T}(norbital, nspin)
    end
end
@inline shape(fock::Fock) = (1:fock.norbital, 1:fock.nspin, 1:2)
@inline Base.eltype(::Type{Fock}) = (FockIndex{T, Int, Rational{Int}} where T)
@inline Base.convert(::Type{<:CartesianIndex}, index::FockIndex{T}, fock::Fock{T}) where T = CartesianIndex(index.orbital, Int(index.spin+(fock.nspin-1)//2)+1, index.nambu)
@inline Base.convert(::Type{<:FockIndex}, index::CartesianIndex{3}, fock::Fock{T}) where T = FockIndex{T}(index[1], index[2]-1-(fock.nspin-1)//2, index[3])
@inline Base.summary(io::IO, fock::Fock) = @printf io "%s-element Fock{%s}" length(fock) repr(statistics(fock))
@inline Base.show(io::IO, fock::Fock) = @printf io "%s{%s}(%s)" fock|>typeof|>nameof repr(statistics(fock)) join(("$name=$(getfield(fock, name))" for name in fock|>typeof|>fieldnames), ", ")
@inline Base.match(::Type{<:FockIndex{:}}, ::Type{<:Fock{T}}) where T = true
@inline Base.match(::Type{<:FockIndex{T}}, ::Type{<:Fock{T}}) where T = true
@inline Base.match(::Type{<:FockIndex{T₁}}, ::Type{<:Fock{T₂}}) where {T₁, T₂} = false
### requested by ConstrainedInternal
@inline function shape(internal::Fock, index::FockIndex{T, <:Union{Int, Symbol, Colon}, <:Union{Rational{Int}, Symbol, Colon}}) where T
    return (fockshape(index.orbital, internal.norbital), fockshape(index.spin, internal.nspin), index.nambu:index.nambu)
end
@inline fockshape(::Union{Symbol, Colon}, n::Int) = 1:n
@inline fockshape(v::Int, n::Int) = ((@assert 0<v<n+1 "shape error: out of range."); v:v)
@inline function fockshape(v::Rational{Int}, n::Int)
    @assert abs(v)<=(n-1)//2 "shape error: out of range."
    index = Int(v+(n-1)//2)+1
    return index:index
end

## Boundary
"""
    angle(id::CoordinatedIndex{<:Index{<:FockIndex}}, vectors::AbstractVector{<:AbstractVector{<:Number}}, values::AbstractVector{<:Number}) -> Complex{<:Number}

Get the twist phase.
"""
function Base.angle(id::CoordinatedIndex{<:Index{<:FockIndex}}, vectors::AbstractVector{<:AbstractVector{<:Number}}, values::AbstractVector{<:Number})
    datatype = promote_type(eltype(values), Float)
    phase = length(vectors)==1 ? 2*convert(datatype, pi)*dot(decompose(id.icoordinate, vectors[1]), values) :
            length(vectors)==2 ? 2*convert(datatype, pi)*dot(decompose(id.icoordinate, vectors[1], vectors[2]), values) :
            length(vectors)==3 ? 2*convert(datatype, pi)*dot(decompose(id.icoordinate, vectors[1], vectors[2], vectors[3]), values) :
            error("angle error: not supported number of input basis vectors.")
    return isannihilation(id) ? phase : iscreation(id) ? -phase : error("angle error: not supported.")
end

## New methods
"""
    isnormalordered(opt::Operator{<:Number, <:ZeroAtLeast{Union{CompositeIndex{<:Index{<:FockIndex}}, Index{<:FockIndex}, FockIndex}}}) -> Bool

Judge whether an operator is normal ordered.
"""
function isnormalordered(opt::Operator{<:Number, <:ZeroAtLeast{Union{CompositeIndex{<:Index{<:FockIndex}}, Index{<:FockIndex}, FockIndex}}})
    flag = true
    for i = 1:rank(opt)
        flag && isannihilation(opt[i]) && (flag = false)
        flag || iscreation(opt[i]) && return false
    end
    return true
end

"""
    *(f₁::Operator{<:Number, <:OneAtLeast{FockIndex{:f}}}, f₂::Operator{<:Number, <:OneAtLeast{FockIndex{:f}}}) -> Operator
    *(f₁::Operator{<:Number, <:OneAtLeast{Index{<:FockIndex{:f}, Int}}}, f₂::Operator{<:Number, <:OneAtLeast{Index{<:FockIndex{:f}, Int}}}) -> Operator
    *(f₁::Operator{<:Number, <:OneAtLeast{CompositeIndex{<:Index{<:FockIndex{:f}, Int}}}}, f₂::Operator{<:Number, <:OneAtLeast{CompositeIndex{<:Index{<:FockIndex{:f}, Int}}}}) -> Operator

Get the multiplication of two fermionic Fock operators.
"""
const block = quote
    result = invoke(*, Tuple{OperatorProd, OperatorProd}, f₁, f₂)
    rank(f₁)>0 && rank(f₂)>0 && f₁[end]==f₂[1] && return replace(result, zero(valtype(result)))
    return result
end
@eval @inline Base.:*(f₁::Operator{<:Number, <:OneAtLeast{FockIndex{:f}}}, f₂::Operator{<:Number, <:OneAtLeast{FockIndex{:f}}}) = $block
@eval @inline Base.:*(f₁::Operator{<:Number, <:OneAtLeast{Index{<:FockIndex{:f}, Int}}}, f₂::Operator{<:Number, <:OneAtLeast{Index{<:FockIndex{:f}, Int}}}) = $block
@eval @inline Base.:*(f₁::Operator{<:Number, <:OneAtLeast{CompositeIndex{<:Index{<:FockIndex{:f}, Int}}}}, f₂::Operator{<:Number, <:OneAtLeast{CompositeIndex{<:Index{<:FockIndex{:f}, Int}}}}) = $block

## Permutation
"""
    permute(id₁::FockIndex, id₂::FockIndex) -> ZeroAtLeast{Operator}

Permute two Fock indexes and get the result.
"""
function permute(id₁::FockIndex{:f}, id₂::FockIndex{:f})
    @assert id₁ ≠ id₂ "permute error: two identical fermionic indexes should vanish due to the fermionic statistics."
    if id₁' == id₂
        return (Operator(1), Operator(-1, id₂, id₁))
    else
        return (Operator(-1, id₂, id₁),)
    end
end
function permute(id₁::FockIndex{:b}, id₂::FockIndex{:b})
    if id₁' == id₂
        return (Operator(iscreation(id₁) ? -1 : +1), Operator(1, id₂, id₁))
    else
        return (Operator(1, id₂, id₁),)
    end
end
@inline permute(id₁::FockIndex{:b}, id₂::FockIndex{:f}) = (Operator(1, id₂, id₁),)
@inline permute(id₁::FockIndex{:f}, id₂::FockIndex{:b}) = (Operator(1, id₂, id₁),)

## Coupling
### Pauli matrices
"""
    const σ⁰ = SparseMatrixCSC([1 0; 0 1])
    const σˣ = SparseMatrixCSC([0 1; 1 0])
    const σʸ = SparseMatrixCSC([0 -1im; 1im 0])
    const σᶻ = SparseMatrixCSC([1 0; 0 -1])
    const σ⁺ = SparseMatrixCSC([0 1; 0 0])
    const σ⁻ = SparseMatrixCSC([0 0; 1 0])
    const σ¹¹ = SparseMatrixCSC([1 0; 0 0])
    const σ¹² = SparseMatrixCSC([0 1; 0 0])
    const σ²¹ = SparseMatrixCSC([0 0; 1 0])
    const σ²² = SparseMatrixCSC([0 0; 0 1])

Pauli matrices σ⁰, σˣ, σʸ, σᶻ, σ⁺, σ⁻, σ¹¹, σ¹², σ²¹ and σ²².
"""
const σ⁰ = SparseMatrixCSC([1 0; 0 1])
const σˣ = SparseMatrixCSC([0 1; 1 0])
const σʸ = SparseMatrixCSC([0 -1im; 1im 0])
const σᶻ = SparseMatrixCSC([1 0; 0 -1])
const σ⁺ = SparseMatrixCSC([0 1; 0 0])
const σ⁻ = SparseMatrixCSC([0 0; 1 0])
const σ¹¹ = SparseMatrixCSC([1 0; 0 0])
const σ¹² = SparseMatrixCSC([0 1; 0 0])
const σ²¹ = SparseMatrixCSC([0 0; 1 0])
const σ²² = SparseMatrixCSC([0 0; 0 1])

### Rotation matrices
"""
    const Lˣ = SparseMatrixCSC([0 0 0; 0 0 1im; 0 -1im 0])
    const Lʸ = SparseMatrixCSC([0 0 -1im; 0 0 0; 1im 0 0])
    const Lᶻ = SparseMatrixCSC([0 1im 0; -1im 0 0; 0 0 0])

Three-dimensional rotation generators Lˣ, Lʸ and Lᶻ.
"""
const Lˣ = SparseMatrixCSC([0 0 0; 0 0 1im; 0 -1im 0])
const Lʸ = SparseMatrixCSC([0 0 -1im; 0 0 0; 1im 0 0])
const Lᶻ = SparseMatrixCSC([0 1im 0; -1im 0 0; 0 0 0])

### MatrixCoupling
const default_matrix = SparseMatrixCSC(hcat(1))
"""
    MatrixCoupling{F}(sites::Union{NTuple{2, Ordinal}, Colon}, orbital::Union{AbstractMatrix, Colon}, spin::Union{AbstractMatrix, Colon}, nambu::AbstractMatrix) where {F<:FockIndex}

Construct a matrix coupling for Fock systems.
"""
@inline function MatrixCoupling{F}(sites::Union{NTuple{2, Ordinal}, Colon}, orbital::Union{AbstractMatrix, Colon}, spin::Union{AbstractMatrix, Colon}, nambu::AbstractMatrix) where {F<:FockIndex}
    return MatrixCoupling{F}(sites, MatrixCouplingComponent(F, Val(:orbital), orbital), MatrixCouplingComponent(F, Val(:spin), spin), MatrixCouplingComponent(F, Val(:nambu), nambu))
end
@inline function MatrixCouplingComponent(::Type{<:FockIndex}, ::Val, ::Colon)
    return MatrixCouplingComponent(SVector(:), SVector(:), default_matrix)
end
@inline function MatrixCouplingComponent(::Type{<:FockIndex}, ::Val{:orbital}, matrix::AbstractMatrix)
    return MatrixCouplingComponent(1:size(matrix)[1], 1:size(matrix)[2], matrix)
end
@inline function MatrixCouplingComponent(::Type{<:FockIndex}, ::Val{:spin}, matrix::AbstractMatrix)
    return MatrixCouplingComponent((size(matrix)[1]-1)//2:-1:(1-size(matrix)[1])//2, (size(matrix)[2]-1)//2:-1:(1-size(matrix)[2])//2, matrix)
end
@inline function MatrixCouplingComponent(::Type{<:FockIndex}, ::Val{:nambu}, matrix::AbstractMatrix)
    @assert size(matrix)==(2, 2) "MatrixCouplingComponent error: for nambu subspace, the input matrix must be 2×2."
    return MatrixCouplingComponent(SVector(creation, annihilation), SVector(annihilation, creation), matrix)
end

"""
    𝕔⁺𝕔(sites::Union{NTuple{2, Ordinal}, Colon}, orbital::Union{AbstractMatrix, Colon}, spin::Union{AbstractMatrix, Colon}) -> MatrixCoupling
    𝕔⁺𝕔⁺(sites::Union{NTuple{2, Ordinal}, Colon}, orbital::Union{AbstractMatrix, Colon}, spin::Union{AbstractMatrix, Colon}) -> MatrixCoupling
    𝕔𝕔(sites::Union{NTuple{2, Ordinal}, Colon}, orbital::Union{AbstractMatrix, Colon}, spin::Union{AbstractMatrix, Colon}) -> MatrixCoupling
    𝕔𝕔⁺(sites::Union{NTuple{2, Ordinal}, Colon}, orbital::Union{AbstractMatrix, Colon}, spin::Union{AbstractMatrix, Colon}) -> MatrixCoupling

    𝕒⁺𝕒(sites::Union{NTuple{2, Ordinal}, Colon}, orbital::Union{AbstractMatrix, Colon}, spin::Union{AbstractMatrix, Colon}) -> MatrixCoupling
    𝕒⁺𝕒⁺(sites::Union{NTuple{2, Ordinal}, Colon}, orbital::Union{AbstractMatrix, Colon}, spin::Union{AbstractMatrix, Colon}) -> MatrixCoupling
    𝕒𝕒(sites::Union{NTuple{2, Ordinal}, Colon}, orbital::Union{AbstractMatrix, Colon}, spin::Union{AbstractMatrix, Colon}) -> MatrixCoupling
    𝕒𝕒⁺(sites::Union{NTuple{2, Ordinal}, Colon}, orbital::Union{AbstractMatrix, Colon}, spin::Union{AbstractMatrix, Colon}) -> MatrixCoupling

    𝕕⁺𝕕(sites::Union{NTuple{2, Ordinal}, Colon}, orbital::Union{AbstractMatrix, Colon}, spin::Union{AbstractMatrix, Colon}) -> MatrixCoupling
    𝕕⁺𝕕⁺(sites::Union{NTuple{2, Ordinal}, Colon}, orbital::Union{AbstractMatrix, Colon}, spin::Union{AbstractMatrix, Colon}) -> MatrixCoupling
    𝕕𝕕(sites::Union{NTuple{2, Ordinal}, Colon}, orbital::Union{AbstractMatrix, Colon}, spin::Union{AbstractMatrix, Colon}) -> MatrixCoupling
    𝕕𝕕⁺(sites::Union{NTuple{2, Ordinal}, Colon}, orbital::Union{AbstractMatrix, Colon}, spin::Union{AbstractMatrix, Colon}) -> MatrixCoupling

Construct a matrix coupling for Fock systems.
"""
@inline 𝕔⁺𝕔(sites::Union{NTuple{2, Ordinal}, Colon}, orbital::Union{AbstractMatrix, Colon}, spin::Union{AbstractMatrix, Colon}) = MatrixCoupling{FockIndex{:f}}(sites, orbital, spin, σ¹¹)
@inline 𝕔⁺𝕔⁺(sites::Union{NTuple{2, Ordinal}, Colon}, orbital::Union{AbstractMatrix, Colon}, spin::Union{AbstractMatrix, Colon}) = MatrixCoupling{FockIndex{:f}}(sites, orbital, spin, σ¹²)
@inline 𝕔𝕔(sites::Union{NTuple{2, Ordinal}, Colon}, orbital::Union{AbstractMatrix, Colon}, spin::Union{AbstractMatrix, Colon}) = MatrixCoupling{FockIndex{:f}}(sites, orbital, spin, σ²¹)
@inline 𝕔𝕔⁺(sites::Union{NTuple{2, Ordinal}, Colon}, orbital::Union{AbstractMatrix, Colon}, spin::Union{AbstractMatrix, Colon}) = MatrixCoupling{FockIndex{:f}}(sites, orbital, spin, σ²²)
@inline 𝕒⁺𝕒(sites::Union{NTuple{2, Ordinal}, Colon}, orbital::Union{AbstractMatrix, Colon}, spin::Union{AbstractMatrix, Colon}) = MatrixCoupling{FockIndex{:b}}(sites, orbital, spin, σ¹¹)
@inline 𝕒⁺𝕒⁺(sites::Union{NTuple{2, Ordinal}, Colon}, orbital::Union{AbstractMatrix, Colon}, spin::Union{AbstractMatrix, Colon}) = MatrixCoupling{FockIndex{:b}}(sites, orbital, spin, σ¹²)
@inline 𝕒𝕒(sites::Union{NTuple{2, Ordinal}, Colon}, orbital::Union{AbstractMatrix, Colon}, spin::Union{AbstractMatrix, Colon}) = MatrixCoupling{FockIndex{:b}}(sites, orbital, spin, σ²¹)
@inline 𝕒𝕒⁺(sites::Union{NTuple{2, Ordinal}, Colon}, orbital::Union{AbstractMatrix, Colon}, spin::Union{AbstractMatrix, Colon}) = MatrixCoupling{FockIndex{:b}}(sites, orbital, spin, σ²²)
@inline 𝕕⁺𝕕(sites::Union{NTuple{2, Ordinal}, Colon}, orbital::Union{AbstractMatrix, Colon}, spin::Union{AbstractMatrix, Colon}) = MatrixCoupling{FockIndex}(sites, orbital, spin, σ¹¹)
@inline 𝕕⁺𝕕⁺(sites::Union{NTuple{2, Ordinal}, Colon}, orbital::Union{AbstractMatrix, Colon}, spin::Union{AbstractMatrix, Colon}) = MatrixCoupling{FockIndex}(sites, orbital, spin, σ¹²)
@inline 𝕕𝕕(sites::Union{NTuple{2, Ordinal}, Colon}, orbital::Union{AbstractMatrix, Colon}, spin::Union{AbstractMatrix, Colon}) = MatrixCoupling{FockIndex}(sites, orbital, spin, σ²¹)
@inline 𝕕𝕕⁺(sites::Union{NTuple{2, Ordinal}, Colon}, orbital::Union{AbstractMatrix, Colon}, spin::Union{AbstractMatrix, Colon}) = MatrixCoupling{FockIndex}(sites, orbital, spin, σ²²)

## Term
"""
    Onsite(id::Symbol, value, coupling=Coupling(𝕕⁺(:, :, :), 𝕕(:, :, :)); ishermitian::Bool=true, amplitude::Union{Function, Nothing}=nothing, ismodulatable::Bool=true)

Onsite term.

Type alias for `Term{:Onsite, id, V, Int, C<:TermCoupling, A<:TermAmplitude}`.
"""
const Onsite{id, V, C<:TermCoupling, A<:TermAmplitude} = Term{:Onsite, id, V, Int, C, A}
@inline function Onsite(id::Symbol, value, coupling=Coupling(𝕕⁺(:, :, :), 𝕕(:, :, :)); ishermitian::Bool=true, amplitude::Union{Function, Nothing}=nothing, ismodulatable::Bool=true)
    return Term{:Onsite}(id, value, 0, coupling, ishermitian; amplitude=amplitude, ismodulatable=ismodulatable)
end

"""
    Hopping(id::Symbol, value, bondkind, coupling=Coupling(𝕕⁺(:, :, :), 𝕕(:, :, :)); amplitude::Union{Function, Nothing}=nothing, ismodulatable::Bool=true)

Hopping term.

Type alias for `Term{:Hopping, id, V, B, C<:TermCoupling, A<:TermAmplitude}`.
"""
const Hopping{id, V, B, C<:TermCoupling, A<:TermAmplitude} = Term{:Hopping, id, V, B, C, A}
@inline function Hopping(id::Symbol, value, bondkind, coupling=Coupling(𝕕⁺(:, :, :), 𝕕(:, :, :)); amplitude::Union{Function, Nothing}=nothing, ismodulatable::Bool=true)
    @assert bondkind≠0 "Hopping error: input bondkind (neighbor) cannot be 0. Use `Onsite` instead."
    return Term{:Hopping}(id, value, bondkind, coupling, false; amplitude=amplitude, ismodulatable=ismodulatable)
end

"""
    Pairing(id::Symbol, value, bondkind, coupling; amplitude::Union{Function, Nothing}=nothing, ismodulatable::Bool=true)

Pairing term.

Type alias for `Term{:Pairing, id, V, B, C<:TermCoupling, A<:TermAmplitude}`.
"""
const Pairing{id, V, B, C<:TermCoupling, A<:TermAmplitude} = Term{:Pairing, id, V, B, C, A}
@inline function Pairing(id::Symbol, value, bondkind, coupling; amplitude::Union{Function, Nothing}=nothing, ismodulatable::Bool=true)
    return Term{:Pairing}(id, value, bondkind, coupling, false; amplitude=amplitude, ismodulatable=ismodulatable)
end
function expand!(operators::Operators, term::Pairing, bond::Bond, hilbert::Hilbert; half::Bool=false)
    argtypes = Tuple{Operators, Term, Bond, Hilbert}
    invoke(expand!, argtypes, operators, term, bond, hilbert; half=half)
    length(bond)==2 && invoke(expand!, argtypes, operators, term, reverse(bond), hilbert; half=half)
    return operators
end

"""
    Hubbard(id::Symbol, value; amplitude::Union{Function, Nothing}=nothing, ismodulatable::Bool=true)

Hubbard term.

Type alias for `Term{:Hubbard, id, V, Int, C<:TermCoupling, A<:TermAmplitude}`.
"""
const Hubbard{id, V, C<:TermCoupling, A<:TermAmplitude} = Term{:Hubbard, id, V, Int, C, A}
@inline function Hubbard(id::Symbol, value; amplitude::Union{Function, Nothing}=nothing, ismodulatable::Bool=true)
    return Term{:Hubbard}(id, value, 0, Coupling(𝕕⁺(:, :, 1//2), 𝕕(:, :, 1//2), 𝕕⁺(:, :, -1//2), 𝕕(:, :, -1//2)), true; amplitude=amplitude, ismodulatable=ismodulatable)
end

"""
    InterOrbitalInterSpin(id::Symbol, value; amplitude::Union{Function, Nothing}=nothing, ismodulatable::Bool=true)

Interorbital-interspin term.

Type alias for `Term{:InterOrbitalInterSpin, id, V, Int, C<:TermCoupling, A<:TermAmplitude}`.
"""
const InterOrbitalInterSpin{id, V, C<:TermCoupling, A<:TermAmplitude} = Term{:InterOrbitalInterSpin, id, V, Int, C, A}
@inline function InterOrbitalInterSpin(id::Symbol, value; amplitude::Union{Function, Nothing}=nothing, ismodulatable::Bool=true)
    return Term{:InterOrbitalInterSpin}(id, value, 0, Coupling(@pattern(𝕕⁺(:, α, σ), 𝕕(:, α, σ), 𝕕⁺(:, β, σ′), 𝕕(:, β, σ′); constraint=α<β && σ≠σ′)), true; amplitude=amplitude, ismodulatable=ismodulatable)
end

"""
    InterOrbitalIntraSpin(id::Symbol, value; amplitude::Union{Function, Nothing}=nothing, ismodulatable::Bool=true)

Interorbital-intraspin term.

Type alias for `Term{:InterOrbitalIntraSpin, id, V, Int, C<:TermCoupling, A<:TermAmplitude}`.
"""
const InterOrbitalIntraSpin{id, V, C<:TermCoupling, A<:TermAmplitude} = Term{:InterOrbitalIntraSpin, id, V, Int, C, A}
@inline function InterOrbitalIntraSpin(id::Symbol, value; amplitude::Union{Function, Nothing}=nothing, ismodulatable::Bool=true)
    return Term{:InterOrbitalIntraSpin}(id, value, 0, Coupling(@pattern(𝕕⁺(:, α, σ), 𝕕(:, α, σ), 𝕕⁺(:, β, σ), 𝕕(:, β, σ); constraint=α<β)), true; amplitude=amplitude, ismodulatable=ismodulatable)
end

"""
    SpinFlip(id::Symbol, value; amplitude::Union{Function, Nothing}=nothing, ismodulatable::Bool=true)

Spin-flip term.

Type alias for `Term{:SpinFlip, id, V, Int, C<:TermCoupling, A<:TermAmplitude}`.
"""
const SpinFlip{id, V, C<:TermCoupling, A<:TermAmplitude} = Term{:SpinFlip, id, V, Int, C, A}
@inline function SpinFlip(id::Symbol, value; amplitude::Union{Function, Nothing}=nothing, ismodulatable::Bool=true)
    return Term{:SpinFlip}(id, value, 0, Coupling(@pattern(𝕕⁺(:, α, 1//2), 𝕕⁺(:, β, -1//2), 𝕕(:, α, -1//2), 𝕕(:, β, 1//2); constraint=α<β)), false; amplitude=amplitude, ismodulatable=ismodulatable)
end

"""
    PairHopping(id::Symbol, value; amplitude::Union{Function, Nothing}=nothing, ismodulatable::Bool=true)

Pair-hopping term.

Type alias for `Term{:PairHopping, id, V, Int, C<:TermCoupling, A<:TermAmplitude}`.
"""
const PairHopping{id, V, C<:TermCoupling, A<:TermAmplitude} = Term{:PairHopping, id, V, Int, C, A}
@inline function PairHopping(id::Symbol, value; amplitude::Union{Function, Nothing}=nothing, ismodulatable::Bool=true)
    return Term{:PairHopping}(id, value, 0, Coupling(@pattern(𝕕⁺(:, α, 1//2), 𝕕⁺(:, α, -1//2), 𝕕(:, β, -1//2), 𝕕(:, β, 1//2); constraint=α<β)), false; amplitude=amplitude, ismodulatable=ismodulatable)
end

"""
    Coulomb(id::Symbol, value, bondkind, coupling=Coupling(𝕕⁺(:, :, :), 𝕕(:, :, :))^2; ishermitian::Bool=true, amplitude::Union{Function, Nothing}=nothing, ismodulatable::Bool=true)

Coulomb term.

Type alias for `Term{:Coulomb, id, V, B, C<:TermCoupling, A<:TermAmplitude}`.
"""
const Coulomb{id, V, B, C<:TermCoupling, A<:TermAmplitude} = Term{:Coulomb, id, V, B, C, A}
@inline function Coulomb(id::Symbol, value, bondkind, coupling=Coupling(𝕕⁺(:, :, :), 𝕕(:, :, :))^2; ishermitian::Bool=true, amplitude::Union{Function, Nothing}=nothing, ismodulatable::Bool=true)
    return Term{:Coulomb}(id, value, bondkind, coupling, ishermitian; amplitude=amplitude, ismodulatable=ismodulatable)
end

"""
    FockTerm

Type alias for `Union{Onsite, Hopping, Pairing, Hubbard, InterOrbitalInterSpin, InterOrbitalIntraSpin, SpinFlip, PairHopping, Coulomb}`.
"""
const FockTerm = Union{Onsite, Hopping, Pairing, Hubbard, InterOrbitalInterSpin, InterOrbitalIntraSpin, SpinFlip, PairHopping, Coulomb}

# SU(2) spin systems
const spinajointmap = Dict('x'=>'x', 'y'=>'y', 'z'=>'z', '+'=>'-', '-'=>'+')

## SpinIndex
"""
    SpinIndex{S, T<:Union{Char, Symbol, Colon}} <: InternalIndex

Spin index, i.e., the internal index to specify the generators of a spin space.
"""
struct SpinIndex{S, T<:Union{Char, Symbol, Colon}} <: InternalIndex
    tag::T
    function SpinIndex{S}(tag::Union{Char, Symbol, Colon}) where S
        @assert isa(S, Rational{Int}) && S.den==2 || isa(S, Int) || S==Colon() "SpinIndex error: not supported spin($S)."
        isa(tag, Char) && @assert tag in ('x', 'y', 'z', '+', '-') "SpinIndex error: not supported tag($tag)."
        new{S, typeof(tag)}(tag)
    end
end
### basic methods of concrete InternalIndex
@inline statistics(::Type{<:SpinIndex}) = :b
@inline isdefinite(::Type{<:SpinIndex{T, Char} where T}) = true
@inline Base.:(==)(index₁::SpinIndex, index₂::SpinIndex) = totalspin(index₁)==totalspin(index₂) && ==(efficientoperations, index₁, index₂)
@inline Base.isequal(index₁::SpinIndex, index₂::SpinIndex) = isequal(totalspin(index₁), totalspin(index₂)) && isequal(efficientoperations, index₁, index₂)
@inline Base.hash(index::SpinIndex, h::UInt) = hash((totalspin(index), index.tag), h)
@inline Base.adjoint(index::SpinIndex) = SpinIndex{totalspin(index)}(spinajointmap[index.tag])
@inline @generated function Base.replace(index::SpinIndex; kwargs...)
    exprs = [:(get(kwargs, $name, getfield(index, $name))) for name in QuoteNode.(fieldnames(index))]
    return :(rawtype(typeof(index)){totalspin(index)}($(exprs...)))
end
### requested by MatrixCoupling
@inline internalindextype(::Type{SpinIndex}, ::Type{T}) where {T<:Union{Char, Symbol, Colon}} = SpinIndex{:, T}
@inline internalindextype(::Type{SpinIndex{S}}, ::Type{T}) where {S, T<:Union{Char, Symbol, Colon}} = SpinIndex{S, T}

"""
    SpinIndex(tag::Union{Char, Symbol, Colon})
    SpinIndex{S}(tag::Union{Char, Symbol, Colon}) where S
    SpinIndex{S, T}(tag::Union{Char, Symbol, Colon}) where {S, T}

Construct a spin index.
"""
@inline SpinIndex(tag::Union{Char, Symbol, Colon}) = SpinIndex{:}(tag)
@inline SpinIndex{S, T}(tag::Union{Char, Symbol, Colon}) where {S, T} = SpinIndex{S}(tag)

"""
    totalspin(::SpinIndex) -> Rational{Int}/Int/Colon
    totalspin(::Type{<:SpinIndex}) -> Rational{Int}/Int/Colon/Float64

    totalspin(::Index{<:SpinIndex}) -> Rational{Int}/Int/Colon
    totalspin(::Type{<:Index{<:SpinIndex}}) -> Rational{Int}/Int/Colon/Float64

    totalspin(::CompositeIndex{<:Index{<:SpinIndex}}) -> Rational{Int}/Int/Colon
    totalspin(::Type{<:CompositeIndex{<:Index{<:SpinIndex}}}) -> Rational{Int}/Int/Colon/Float64

Get the total spin.
"""
@inline totalspin(index::SpinIndex) = totalspin(typeof(index))
@inline totalspin(index::Index{<:SpinIndex}) = totalspin(typeof(index))
@inline totalspin(index::CompositeIndex{<:Index{<:SpinIndex}}) = totalspin(typeof(index))
@inline totalspin(::Type{<:SpinIndex}) = NaN
@inline totalspin(::Type{<:SpinIndex{S}}) where S = S
@inline totalspin(::Type{I}) where {I<:Index{<:SpinIndex}} = totalspin(internalindextype(I))
@inline totalspin(::Type{I}) where {I<:CompositeIndex{<:Index{<:SpinIndex}}} = totalspin(indextype(I))

### convenient construction and string representation 
"""
    𝕊(tag) -> SpinIndex
    𝕊(site, tag) -> Index{<:SpinIndex}
    𝕊(site, tag, rcoordinate, icoordinate) -> CoordinatedIndex{<:Index{<:SpinIndex}}

    𝕊{S}(tag) where S -> SpinIndex{S}
    𝕊{S}(site, tag) where S -> Index{<:SpinIndex{S}}
    𝕊{S}(site, tag, rcoordinate, icoordinate) where S -> CoordinatedIndex{<:Index{<:SpinIndex{S}}}

Convenient construction of `SpinIndex`, `Index{<:SpinIndex}`, `CoordinatedIndex{<:Index{<:SpinIndex}}`.
"""
struct 𝕊{S} <: Function end
@inline 𝕊(tag) = SpinIndex(tag)
@inline 𝕊(site, tag) = Index(site, SpinIndex(tag))
@inline 𝕊(site, tag, rcoordinate, icoordinate) = CoordinatedIndex(Index(site, SpinIndex(tag)), rcoordinate, icoordinate)
@inline 𝕊{S}(tag) where S = SpinIndex{S}(tag)
@inline 𝕊{S}(site, tag) where S = Index(site, SpinIndex{S}(tag))
@inline 𝕊{S}(site, tag, rcoordinate, icoordinate) where S = CoordinatedIndex(Index(site, SpinIndex{S}(tag)), rcoordinate, icoordinate)
@inline Base.getindex(::Type{OperatorIndex}, index::SpinIndex{:}) = "𝕊"
@inline Base.getindex(::Type{OperatorIndex}, index::SpinIndex) = "𝕊{$(totalspin(index))}"

### matrix
"""
    matrix(index::Union{SpinIndex{S, Char}, Index{SpinIndex{S, Char}}, CompositeIndex{<:Index{SpinIndex{S, Char}}}}, dtype::Type{<:Number}=ComplexF64) where S -> Matrix{dtype}

Get the matrix representation of an index acting on the local spin space.
"""
function matrix(index::Union{SpinIndex{S, Char}, Index{SpinIndex{S, Char}}, CompositeIndex{<:Index{SpinIndex{S, Char}}}}, dtype::Type{<:Number}=ComplexF64) where S
    N = Int(2*S+1)
    result = zeros(dtype, (N, N))
    spin = convert(dtype, S)
    index = InternalIndex(index)
    for i = 1:N, j = 1:N
        # row, col = N+1-i, N+1-j # Sᶻ in ascending order
        row, col = i, j # Sᶻ in descending order
        m, n = spin+1-i, spin+1-j
        result[row, col] = (index.tag == 'x') ? (delta(i+1, j)+delta(i, j+1))*sqrt(spin*(spin+1)-m*n)/2 :
            (index.tag == 'y') ? (delta(i+1, j)-delta(i, j+1))*sqrt(spin*(spin+1)-m*n)/2im :
            (index.tag == 'z') ? delta(i, j)*m :
            (index.tag == '+') ? delta(i+1, j)*sqrt(spin*(spin+1)-m*n) :
            delta(i, j+1)*sqrt(spin*(spin+1)-m*n)
    end
    return result
end

## LaTeX format output
"""
    script(index::SpinIndex, ::Val{:tag}; kwargs...) -> String

Get the requested script of a spin index.
"""
@inline script(index::SpinIndex{S, Colon}, ::Val{:tag}; kwargs...) where S = ":"
@inline script(index::SpinIndex, ::Val{:tag}; kwargs...) = string(index.tag)

"""
    latexofspins

Default LaTeX format of a spin index.
"""
const latexofspins = LaTeX{(:tag,), (:site,)}('S')
@inline latexname(::Type{<:SpinIndex}) = Symbol("SpinIndex")
@inline latexname(::Type{<:Index{<:SpinIndex}}) = Symbol("Index{SpinIndex}")
@inline latexname(::Type{<:CompositeIndex{<:Index{<:SpinIndex}}}) = Symbol("CompositeIndex{Index{SpinIndex}}")
latexformat(SpinIndex, latexofspins)
latexformat(Index{<:SpinIndex}, latexofspins)
latexformat(CompositeIndex{<:Index{<:SpinIndex}}, latexofspins)

## Spin
"""
    Spin{S} <: SimpleInternal{SpinIndex{S, Char}}

Spin space of spin generators at a single point.
"""
struct Spin{S} <: SimpleInternal{SpinIndex{S, Char}}
    function Spin{S}() where S
        @assert isa(S, Rational{Int}) && S.den==2 || isa(S, Int) "Spin error: not supported spin($S)."
        new{S}()
    end
end
@inline shape(::Spin) = (1:3,)
@inline Base.eltype(::Type{Spin}) = (SpinIndex{S, Char} where S)
@inline Base.convert(::Type{<:CartesianIndex}, index::SpinIndex, ::Spin) = CartesianIndex(Int(index.tag)-Int('x')+1)
@inline Base.convert(::Type{<:SpinIndex}, index::CartesianIndex{1}, sp::Spin) = SpinIndex{totalspin(sp)}(Char(index[1]+Int('x')-1))
@inline Base.summary(io::IO, spin::Spin) = @printf io "%s-element Spin{%s}" length(spin) totalspin(spin)
@inline Base.show(io::IO, spin::Spin) = @printf io "%s{%s}()" spin|>typeof|>nameof totalspin(spin)
@inline Base.match(::Type{<:SpinIndex{:}}, ::Type{<:Spin{S}}) where S = true
@inline Base.match(::Type{<:SpinIndex{S}}, ::Type{<:Spin{S}}) where S = true
@inline Base.match(::Type{<:SpinIndex{S₁}}, ::Type{<:Spin{S₂}}) where {S₁, S₂} = false
### requested by ConstrainedInternal
@inline function shape(::Spin, index::SpinIndex)
    @assert isa(index.tag, Char) "shape error: input index ($index) not definite."
    pos = Int(index.tag) - Int('x') + 1
    return (pos:pos,)
end

"""
    totalspin(::Spin) -> Rational{Int}/Int/Colon
    totalspin(::Type{<:Spin}) -> Rational{Int}/Int/Colon

Get the total spin.
"""
@inline totalspin(spin::Spin) = totalspin(typeof(spin))
@inline totalspin(::Type{<:Spin{S}}) where S = S

## Permutation
"""
    permute(id₁::SpinIndex, id₂::SpinIndex) -> ZeroAtLeast{Operator}

Permute two spin indexes and get the result.
"""
function permute(id₁::SpinIndex, id₂::SpinIndex)
    if id₁ ≠ id₂
        S = totalspin(id₁)
        if id₁.tag == 'x'
            id₂.tag=='y' && return (Operator(+1im, 𝕊{S}('z')), Operator(1, id₂, id₁))
            id₂.tag=='z' && return (Operator(-1im, 𝕊{S}('y')), Operator(1, id₂, id₁))
            id₂.tag=='+' && return (Operator(-1, 𝕊{S}('z')), Operator(1, id₂, id₁))
            id₂.tag=='-' && return (Operator(+1, 𝕊{S}('z')), Operator(1, id₂, id₁))
        elseif id₁.tag == 'y'
            id₂.tag=='x' && return (Operator(-1im, 𝕊{S}('z')), Operator(1, id₂, id₁))
            id₂.tag=='z' && return (Operator(+1im, 𝕊{S}('x')), Operator(1, id₂, id₁))
            id₂.tag=='+' && return (Operator(-1im, 𝕊{S}('z')), Operator(1, id₂, id₁))
            id₂.tag=='-' && return (Operator(-1im, 𝕊{S}('z')), Operator(1, id₂, id₁))
        elseif id₁.tag == 'z'
            id₂.tag=='x' && return (Operator(+1im, 𝕊{S}('y')), Operator(1, id₂, id₁))
            id₂.tag=='y' && return (Operator(-1im, 𝕊{S}('x')), Operator(1, id₂, id₁))
            id₂.tag=='+' && return (Operator(+1, 𝕊{S}('+')), Operator(1, id₂, id₁))
            id₂.tag=='-' && return (Operator(-1, 𝕊{S}('-')), Operator(1, id₂, id₁))
        elseif id₁.tag == '+'
            id₂.tag=='x' && return (Operator(+1, 𝕊{S}('z')), Operator(1, id₂, id₁))
            id₂.tag=='y' && return (Operator(+1im, 𝕊{S}('z')), Operator(1, id₂, id₁))
            id₂.tag=='z' && return (Operator(-1, 𝕊{S}('+')), Operator(1, id₂, id₁))
            id₂.tag=='-' && return (Operator(+2, 𝕊{S}('z')), Operator(1, id₂, id₁))
        elseif id₁.tag == '-'
            id₂.tag=='x' && return (Operator(-1, 𝕊{S}('z')), Operator(1, id₂, id₁))
            id₂.tag=='y' && return (Operator(1im, 𝕊{S}('z')), Operator(1, id₂, id₁))
            id₂.tag=='z' && return (Operator(+1, 𝕊{S}('-')), Operator(1, id₂, id₁))
            id₂.tag=='+' && return (Operator(-2, 𝕊{S}('z')), Operator(1, id₂, id₁))
        end
        error("permute error: not supported spin indexes.")
    else
        return (Operator(1, id₂, id₁),)
    end
end

## Coupling
### Spin coupling matrix
"""
    const Isingˣ = SparseMatrixCSC([1 0 0; 0 0 0; 0 0 0])
    const Isingʸ = SparseMatrixCSC([0 0 0; 0 1 0; 0 0 0])
    const Isingᶻ = SparseMatrixCSC([0 0 0; 0 0 0; 0 0 1])

Ising coupling matrices Isingˣ, Isingʸ and Isingᶻ.
"""
const Isingˣ = SparseMatrixCSC([1 0 0; 0 0 0; 0 0 0])
const Isingʸ = SparseMatrixCSC([0 0 0; 0 1 0; 0 0 0])
const Isingᶻ = SparseMatrixCSC([0 0 0; 0 0 0; 0 0 1])

"""
    const Γˣ = SparseMatrixCSC([0 0 0; 0 0 1; 0 1 0])
    const Γʸ = SparseMatrixCSC([0 0 1; 0 0 0; 1 0 0])
    const Γᶻ = SparseMatrixCSC([0 1 0; 1 0 0; 0 0 0])

Γ coupling matrices Γˣ, Γʸ and Γᶻ.
"""
const Γˣ = SparseMatrixCSC([0 0 0; 0 0 1; 0 1 0])
const Γʸ = SparseMatrixCSC([0 0 1; 0 0 0; 1 0 0])
const Γᶻ = SparseMatrixCSC([0 1 0; 1 0 0; 0 0 0])

"""
    const Γ′ˣ = SparseMatrixCSC([0 1 1; 1 0 0; 1 0 0])
    const Γ′ʸ = SparseMatrixCSC([0 1 0; 1 0 1; 0 1 0])
    const Γ′ᶻ = SparseMatrixCSC([0 0 1; 0 0 1; 1 1 0])

Γ′ coupling matrices Γ′ˣ, Γ′ʸ and Γ′ᶻ.
"""
const Γ′ˣ = SparseMatrixCSC([0 1 1; 1 0 0; 1 0 0])
const Γ′ʸ = SparseMatrixCSC([0 1 0; 1 0 1; 0 1 0])
const Γ′ᶻ = SparseMatrixCSC([0 0 1; 0 0 1; 1 1 0])

"""
    const DMˣ = SparseMatrixCSC([0 0 0; 0 0 1; 0 -1 0])
    const DMʸ = SparseMatrixCSC([0 0 -1; 0 0 0; 1 0 0])
    const DMᶻ = SparseMatrixCSC([0 1 0; -1 0 0; 0 0 0])

DM coupling matrices DMˣ, DMʸ and DMᶻ.
"""
const DMˣ = SparseMatrixCSC([0 0 0; 0 0 1; 0 -1 0])
const DMʸ = SparseMatrixCSC([0 0 -1; 0 0 0; 1 0 0])
const DMᶻ = SparseMatrixCSC([0 1 0; -1 0 0; 0 0 0])

### MatrixCoupling
"""
    MatrixCoupling{S}(sites::Union{NTuple{2, Ordinal}, Colon}, matrix::AbstractMatrix; rows::AbstractVector=SVector('x', 'y', 'z'), cols::AbstractVector=SVector('x', 'y', 'z')) where {S<:SpinIndex}

Construct a matrix coupling for spin systems.
"""
@inline function MatrixCoupling{S}(sites::Union{NTuple{2, Ordinal}, Colon}, matrix::AbstractMatrix; rows::AbstractVector=SVector('x', 'y', 'z'), cols::AbstractVector=SVector('x', 'y', 'z')) where {S<:SpinIndex}
    @assert size(matrix)==(length(rows), length(cols)) "MatrixCoupling error: mismatched input matrix and rows/cols."
    return MatrixCoupling{S}(sites, MatrixCouplingComponent(rows, cols, matrix))
end

"""
    𝕊ᵀ𝕊(sites::Union{NTuple{2, Ordinal}, Colon}, matrix::AbstractMatrix; rows::AbstractVector=SVector('x', 'y', 'z'), cols::AbstractVector=SVector('x', 'y', 'z')) -> MatrixCoupling

Construct a matrix coupling for spin system.
"""
@inline 𝕊ᵀ𝕊(sites::Union{NTuple{2, Ordinal}, Colon}, matrix::AbstractMatrix; rows::AbstractVector=SVector('x', 'y', 'z'), cols::AbstractVector=SVector('x', 'y', 'z')) = MatrixCoupling{SpinIndex}(sites, matrix; rows=rows, cols=cols)

## Term
"""
    SpinTerm(id::Symbol, value, bondkind, coupling; amplitude::Union{Function, Nothing}=nothing, ismodulatable::Bool=true)

Generic spin term.

Type alias for `Term{:SpinTerm, id, V, B, C<:TermCoupling, A<:TermAmplitude}`.
"""
const SpinTerm{id, V, B, C<:TermCoupling, A<:TermAmplitude} = Term{:SpinTerm, id, V, B, C, A}
@inline function SpinTerm(id::Symbol, value, bondkind, coupling; amplitude::Union{Function, Nothing}=nothing, ismodulatable::Bool=true)
    return Term{:SpinTerm}(id, value, bondkind, coupling, true; amplitude=amplitude, ismodulatable=ismodulatable)
end
@inline function patternrule(sites::OneAtLeast{Colon}, ::Val{:SpinTerm}, bondlength::Integer)
    bondlength==1 && return ntuple(i->Ordinal(1), Val(fieldcount(typeof(sites))))
    bondlength==2 && return ntuple(i->Ordinal(2-i%2), Val(fieldcount(typeof(sites))))
    error("patternrule error: not supported for a generic bond containing $bondlength points.")
end

"""
    Zeeman(
        id::Symbol, value, direction::Char, g::Number=1;
        amplitude::Union{Function, Nothing}=nothing, ismodulatable::Bool=true
    )
    Zeeman(
        id::Symbol, value, direction::Union{AbstractVector{<:Number}, Tuple{Number, Number}}, g::Union{Number, AbstractMatrix{<:Number}}=1;
        unit::Symbol=:degree, amplitude::Union{Function, Nothing}=nothing, ismodulatable::Bool=true
    )

Zeeman term.

Type alias for `Term{:Zeeman, id, V, Int, C<:TermCoupling, A<:TermAmplitude}`.
"""
const Zeeman{id, V, C<:TermCoupling, A<:TermAmplitude} = Term{:Zeeman, id, V, Int, C, A}
@inline function Zeeman(
    id::Symbol, value, direction::Char, g::Number=1;
    amplitude::Union{Function, Nothing}=nothing, ismodulatable::Bool=true
)
    @assert lowercase(direction)∈('x', 'y', 'z') "Zeeman error: not supported direction."
    coupling = Coupling(g, 𝕊(:, lowercase(direction)))
    return Term{:Zeeman}(id, value, 0, coupling, true; amplitude=amplitude, ismodulatable=ismodulatable)
end
@inline function Zeeman(
    id::Symbol, value, dir::Union{AbstractVector{<:Number}, Tuple{Number, Number}}, g::Union{Number, AbstractMatrix{<:Number}}=1;
    unit::Symbol=:degree, amplitude::Union{Function, Nothing}=nothing, ismodulatable::Bool=true
)
    couplings = dot(direction(dir, unit), Lande(g), SVector(Coupling(𝕊(:, 'x')), Coupling(𝕊(:, 'y')), Coupling(𝕊(:, 'z'))))
    return Term{:Zeeman}(id, value, 0, couplings, true; amplitude=amplitude, ismodulatable=ismodulatable)
end
@inline Lande(g::Number) = SMatrix{3, 3}(g, 0, 0, 0, g, 0, 0, 0, g)
@inline Lande(g::AbstractMatrix{<:Number}) = (@assert(size(g)==(3, 3), "Lande error: the g-tensor must be 3×3."); g)

"""
    SingleIonAnisotropy(id::Symbol, value, direction::Char; amplitude::Union{Function, Nothing}=nothing, ismodulatable::Bool=true)
    SingleIonAnisotropy(id::Symbol, value, matrix::AbstractMatrix{<:Number}; amplitude::Union{Function, Nothing}=nothing, ismodulatable::Bool=true)

Single ion anisotropy term.

Type alias for `Term{:SingleIonAnisotropy, id, V, Int, C<:TermCoupling, A<:TermAmplitude}`.
"""
const SingleIonAnisotropy{id, V, C<:TermCoupling, A<:TermAmplitude} = Term{:SingleIonAnisotropy, id, V, Int, C, A}
@inline function SingleIonAnisotropy(id::Symbol, value, direction::Char; amplitude::Union{Function, Nothing}=nothing, ismodulatable::Bool=true)
    @assert lowercase(direction)∈('x', 'y', 'z') "SingleIonAnisotropy error: not supported direction."
    coupling = Coupling{𝕊}(:, (lowercase(direction), lowercase(direction)))
    return Term{:SingleIonAnisotropy}(id, value, 0, coupling, true; amplitude=amplitude, ismodulatable=ismodulatable)
end
@inline function SingleIonAnisotropy(id::Symbol, value, matrix::AbstractMatrix{<:Number}; amplitude::Union{Function, Nothing}=nothing, ismodulatable::Bool=true)
    @assert ishermitian(matrix) "SingleIonAnisotropy error: the anisotropy matrix must be Hermitian."
    @assert size(matrix)==(3, 3) "SingleIonAnisotropy error: the anisotropy matrix must be 3×3."
    return Term{:SingleIonAnisotropy}(id, value, 0, 𝕊ᵀ𝕊(:, matrix), true; amplitude=amplitude, ismodulatable=ismodulatable)
end

"""
    Ising(id::Symbol, value, bondkind, direction::Char; amplitude::Union{Function, Nothing}=nothing, ismodulatable::Bool=true)

Ising term.

Type alias for `Term{:Ising, id, V, B, C<:TermCoupling, A<:TermAmplitude}`.
"""
const Ising{id, V, B, C<:TermCoupling, A<:TermAmplitude} = Term{:Ising, id, V, B, C, A}
@inline function Ising(id::Symbol, value, bondkind, direction::Char; amplitude::Union{Function, Nothing}=nothing, ismodulatable::Bool=true)
    @assert lowercase(direction)∈('x', 'y', 'z') "Ising error: not supported direction."
    coupling = Coupling{𝕊}(:, (lowercase(direction), lowercase(direction)))
    return Term{:Ising}(id, value, bondkind, coupling, true; amplitude=amplitude, ismodulatable=ismodulatable)
end

"""
    Heisenberg(id::Symbol, value, bondkind; form::Symbol=Symbol("+-z"), amplitude::Union{Function, Nothing}=nothing, ismodulatable::Bool=true)

Heisenberg term.

Type alias for `Term{:Heisenberg, id, V, B, C<:TermCoupling, A<:TermAmplitude}`.
"""
const Heisenberg{id, V, B, C<:TermCoupling, A<:TermAmplitude} = Term{:Heisenberg, id, V, B, C, A}
@inline function Heisenberg(id::Symbol, value, bondkind; form::Symbol=Symbol("+-z"), amplitude::Union{Function, Nothing}=nothing, ismodulatable::Bool=true)
    @assert form∈(:xyz, Symbol("+-z")) "Heisenberg error: form should :xyz or Symbol(\"+-z\")."
    couplings = if form==:xyz
        Coupling{𝕊}(1//1, :, ('x', 'x')) + Coupling{𝕊}(1//1, :, ('y', 'y')) + Coupling{𝕊}(1//1, :, ('z', 'z'))
    else
        Coupling{𝕊}(1//2, :, ('+', '-')) + Coupling{𝕊}(1//2, :, ('-', '+')) + Coupling{𝕊}(1//1, :, ('z', 'z'))
    end
    return Term{:Heisenberg}(id, value, bondkind, couplings, true; amplitude=amplitude, ismodulatable=ismodulatable)
end

"""
    Kitaev(
        id::Symbol, value, bondkind;
        x::AbstractVector{<:Union{Number, Tuple{Number, Number}, AbstractVector{<:Number}}},
        y::AbstractVector{<:Union{Number, Tuple{Number, Number}, AbstractVector{<:Number}}},
        z::AbstractVector{<:Union{Number, Tuple{Number, Number}, AbstractVector{<:Number}}},
        unit::Symbol=:degree,
        amplitude::Union{Function, Nothing}=nothing,
        ismodulatable::Bool=true
    )

Kitaev term. Since Kitaev term is symmetric on every bond, only one direction of a bond is needed. The inverse direction of a bond can be handled automatically by this function.

Here, `x`, `y` and `z` assign the x-bonds, y-bonds, and z-bonds, respectively, with each kind of bond can be
1) a `Number` specifying the azimuth angle of a bond in the 2-dimensional case, or
2) a `Tuple{Number, Number}` specifying the polar and azimuth angle pairs of a bond in the 3-dimensional case, or
3) an `AbstractVector{<:Number}` specifying the direction of a bond.

Type alias for `Term{:Kitaev, id, V, B, C<:TermCoupling, A<:TermAmplitude}`.
"""
const Kitaev{id, V, B, C<:TermCoupling, A<:TermAmplitude} = Term{:Kitaev, id, V, B, C, A}
function Kitaev(
    id::Symbol, value, bondkind;
    x::AbstractVector{<:Union{Number, Tuple{Number, Number}, AbstractVector{<:Number}}},
    y::AbstractVector{<:Union{Number, Tuple{Number, Number}, AbstractVector{<:Number}}},
    z::AbstractVector{<:Union{Number, Tuple{Number, Number}, AbstractVector{<:Number}}},
    unit::Symbol=:degree,
    amplitude::Union{Function, Nothing}=nothing,
    ismodulatable::Bool=true
)
    dirs = (x=direction.(x, unit), y=direction.(y, unit), z=direction.(z, unit))
    function kitaev(bond::Bond)
        coordinate = rcoordinate(bond)
        any(v->abs(isparallel(v, coordinate; atol=atol, rtol=rtol))==1, dirs.x) && return 𝕊ᵀ𝕊(:, Isingˣ)
        any(v->abs(isparallel(v, coordinate; atol=atol, rtol=rtol))==1, dirs.y) && return 𝕊ᵀ𝕊(:, Isingʸ)
        any(v->abs(isparallel(v, coordinate; atol=atol, rtol=rtol))==1, dirs.z) && return 𝕊ᵀ𝕊(:, Isingᶻ)
        error("Kitaev error: wrong bond.")
    end
    return Term{:Kitaev}(id, value, bondkind, kitaev, true; amplitude=amplitude, ismodulatable=ismodulatable)
end

"""
    Γ(
        id::Symbol, value, bondkind;
        x::AbstractVector{<:Union{Number, Tuple{Number, Number}, AbstractVector{<:Number}}},
        y::AbstractVector{<:Union{Number, Tuple{Number, Number}, AbstractVector{<:Number}}},
        z::AbstractVector{<:Union{Number, Tuple{Number, Number}, AbstractVector{<:Number}}},
        unit::Symbol=:degree,
        amplitude::Union{Function, Nothing}=nothing,
        ismodulatable::Bool=true
    )

Γ Term. Since Γ term is symmetric on every bond, only one direction of a bond is needed. The inverse direction of a bond can be handled automatically by this function.

Here, `x`, `y` and `z` assign the x-bonds, y-bonds, and z-bonds, respectively, with each kind of bond can be
1) a `Number` specifying the azimuth angle of a bond in the 2-dimensional case, or
2) a `Tuple{Number, Number}` specifying the polar and azimuth angle pairs of a bond in the 3-dimensional case, or
3) an `AbstractVector{<:Number}` specifying the direction of a bond.

Type alias for `Term{:Γ, id, V, B, C<:TermCoupling, A<:TermAmplitude}`.
"""
const Γ{id, V, B, C<:TermCoupling, A<:TermAmplitude} = Term{:Γ, id, V, B, C, A}
function Γ(
    id::Symbol, value, bondkind;
    x::AbstractVector{<:Union{Number, Tuple{Number, Number}, AbstractVector{<:Number}}},
    y::AbstractVector{<:Union{Number, Tuple{Number, Number}, AbstractVector{<:Number}}},
    z::AbstractVector{<:Union{Number, Tuple{Number, Number}, AbstractVector{<:Number}}},
    unit::Symbol=:degree,
    amplitude::Union{Function, Nothing}=nothing,
    ismodulatable::Bool=true
)
    dirs = (x=direction.(x, unit), y=direction.(y, unit), z=direction.(z, unit))
    function γ(bond::Bond)
        coordinate = rcoordinate(bond)
        any(v->abs(isparallel(v, coordinate; atol=atol, rtol=rtol))==1, dirs.x) && return 𝕊ᵀ𝕊(:, Γˣ)
        any(v->abs(isparallel(v, coordinate; atol=atol, rtol=rtol))==1, dirs.y) && return 𝕊ᵀ𝕊(:, Γʸ)
        any(v->abs(isparallel(v, coordinate; atol=atol, rtol=rtol))==1, dirs.z) && return 𝕊ᵀ𝕊(:, Γᶻ)
        error("Γ error: wrong bond.")
    end
    return Term{:Γ}(id, value, bondkind, γ, true; amplitude=amplitude, ismodulatable=ismodulatable)
end

"""
    Γ′(
        id::Symbol, value, bondkind;
        x::AbstractVector{<:Union{Number, Tuple{Number, Number}, AbstractVector{<:Number}}},
        y::AbstractVector{<:Union{Number, Tuple{Number, Number}, AbstractVector{<:Number}}},
        z::AbstractVector{<:Union{Number, Tuple{Number, Number}, AbstractVector{<:Number}}},
        unit::Symbol=:degree,
        amplitude::Union{Function, Nothing}=nothing,
        ismodulatable::Bool=true
    )

Γ′ Term. Since Γ′ term is symmetric on every bond, only one direction of a bond is needed. The inverse direction of a bond can be handled automatically by this function.

Here, `x`, `y` and `z` assign the x-bonds, y-bonds, and z-bonds, respectively, with each bond can be
1) a `Number` specifying the azimuth angle of a bond in the 2-dimensional case, or
2) a `Tuple{Number, Number}` specifying the polar and azimuth angle pairs of a bond in the 3-dimensional case, or
3) an `AbstractVector{<:Number}` specifying the direction of a bond.

Type alias for `Term{:Γ′, id, V, B, C<:TermCoupling, A<:TermAmplitude}`.
"""
const Γ′{id, V, B, C<:TermCoupling, A<:TermAmplitude} = Term{:Γ′, id, V, B, C, A}
function Γ′(
    id::Symbol, value, bondkind;
    x::AbstractVector{<:Union{Number, Tuple{Number, Number}, AbstractVector{<:Number}}},
    y::AbstractVector{<:Union{Number, Tuple{Number, Number}, AbstractVector{<:Number}}},
    z::AbstractVector{<:Union{Number, Tuple{Number, Number}, AbstractVector{<:Number}}},
    unit::Symbol=:degree,
    amplitude::Union{Function, Nothing}=nothing,
    ismodulatable::Bool=true
)
    dirs = (x=direction.(x, unit), y=direction.(y, unit), z=direction.(z, unit))
    function γ′(bond::Bond)
        coordinate = rcoordinate(bond)
        any(v->abs(isparallel(v, coordinate; atol=atol, rtol=rtol))==1, dirs.x) && return 𝕊ᵀ𝕊(:, Γ′ˣ)
        any(v->abs(isparallel(v, coordinate; atol=atol, rtol=rtol))==1, dirs.y) && return 𝕊ᵀ𝕊(:, Γ′ʸ)
        any(v->abs(isparallel(v, coordinate; atol=atol, rtol=rtol))==1, dirs.z) && return 𝕊ᵀ𝕊(:, Γ′ᶻ)
        error("Γ′ error: wrong bond.")
    end
    return Term{:Γ′}(id, value, bondkind, γ′, true; amplitude=amplitude, ismodulatable=ismodulatable)
end

"""
    DM(
        id::Symbol,
        value,
        bondkind,
        vectors::Pair{<:AbstractVector{<:Union{Number, Tuple{Number, Number}, AbstractVector{<:Number}}}, <:Union{Char, AbstractVector{<:Number}}}...;
        unit::Symbol=:degree,
        amplitude::Union{Function, Nothing}=nothing,
        ismodulatable::Bool=true
    )

DM term. Since DM term is antisymmetric on every bond, only the positive direction of a bond is needed. The negative direction of a bond can be handled automatically by this function.

Here, `vectors` specify the unit DM vector on every bond in the form `[bond₁, bond₂, ...]=>v`, where `bondᵢ` can be
1) a `Number` specifying the azimuth angle of a bond in the 2-dimensional case, or
2) a `Tuple{Number, Number}` specifying the polar and azimuth angle pairs of a bond in the 3-dimensional case, or
3) an `AbstractVector{<:Number}` specifying the direction of a bond;
and `v` can be
1) a `Char` of 'x', 'y' or 'z', indicating the unit DM vector on the set of bonds is along the x, y or z direction, or
2) an `AbstractVector{<:Number}`, specifying the direction of the DM vector on the set of bonds.

Type alias for `Term{:DM, id, V, B, C<:TermCoupling, A<:TermAmplitude}`.
"""
const DM{id, V, B, C<:TermCoupling, A<:TermAmplitude} = Term{:DM, id, V, B, C, A}
function DM(
    id::Symbol,
    value,
    bondkind,
    vectors::Pair{<:AbstractVector{<:Union{Number, Tuple{Number, Number}, AbstractVector{<:Number}}}, <:Union{Char, AbstractVector{<:Number}}}...;
    unit::Symbol=:degree,
    amplitude::Union{Function, Nothing}=nothing,
    ismodulatable::Bool=true
)
    dirs = [direction.(pair.first, unit)=>direction(pair.second) for pair in vectors]
    function dm(bond::Bond)
        coordinate = rcoordinate(bond)
        for pair in dirs
            for v in pair.first
                parallel = isparallel(v, coordinate; atol=atol, rtol=rtol)
                abs(parallel)==1 && return 𝕊ᵀ𝕊(:, parallel*(pair.second[1]*DMˣ + pair.second[2]*DMʸ + pair.second[3]*DMᶻ))
            end
        end
        error("dm error: wrong bond.")
    end
    return Term{:DM}(id, value, bondkind, dm, true; amplitude=amplitude, ismodulatable=ismodulatable)
end

# Phononic systems
## PhononIndex
"""
    PhononIndex{K, D<:Union{Char, Symbol, Colon}} <: InternalIndex

Phonon index, i.e., the internal index to specify the generators of the vibration of a lattice point.
"""
struct PhononIndex{K, D<:Union{Char, Symbol, Colon}} <: InternalIndex
    direction::D
    function PhononIndex{K}(direction::Union{Char, Symbol, Colon}) where K
        @assert K∈(:u, :p) "PhononIndex error: wrong kind($kind)."
        isa(direction, Char) && @assert direction∈('x', 'y', 'z') "PhononIndex error: wrong direction($direction)."
        new{K, typeof(direction)}(direction)
    end
end
### basic methods of concrete InternalIndex
@inline statistics(::Type{<:PhononIndex}) = :b
@inline isdefinite(::Type{<:PhononIndex{K, Char} where K}) = true
@inline Base.:(==)(index₁::PhononIndex, index₂::PhononIndex) = kind(index₁)==kind(index₂) && ==(efficientoperations, index₁, index₂)
@inline Base.isequal(index₁::PhononIndex, index₂::PhononIndex) = isequal(kind(index₁), kind(index₂)) && isequal(efficientoperations, index₁, index₂)
@inline Base.hash(index::PhononIndex, h::UInt) = hash((kind(index), index.direction), h)
@inline Base.adjoint(index::PhononIndex) = index
@inline @generated function Base.replace(index::PhononIndex; kwargs...)
    exprs = [:(get(kwargs, $name, getfield(index, $name))) for name in QuoteNode.(fieldnames(index))]
    return :(rawtype(typeof(index)){kind(index)}($(exprs...)))
end
### requested by MatrixCoupling
@inline internalindextype(::Type{PhononIndex{K}}, ::Type{D}) where {K, D<:Union{Char, Symbol, Colon}} = PhononIndex{K, D}

"""
    PhononIndex{K}(direction::Union{Char, Symbol, Colon}) where K
    PhononIndex{K, D}(direction::Union{Char, Symbol, Colon}) where {K, D}

Construct a phonon index.
"""
@inline PhononIndex{K, D}(direction::Union{Char, Symbol, Colon}) where {K, D} = PhononIndex{K}(direction)

### convenient construction and string representation 
"""
    𝕦(direction) -> PhononIndex{:u}
    𝕦(site, direction) -> Index{<:PhononIndex{:u}}
    𝕦(site, direction, rcoordinate, icoordinate) -> CoordinatedIndex{<:Index{<:PhononIndex{:u}}}

Convenient construction of `SpinIndex{:u}`, `Index{<:SpinIndex{:u}}`, `CoordinatedIndex{<:Index{<:SpinIndex{:u}}}`.
"""
@inline 𝕦(direction) = PhononIndex{:u}(direction)
@inline 𝕦(site, direction) = Index(site, PhononIndex{:u}(direction))
@inline 𝕦(site, direction, rcoordinate, icoordinate) = CoordinatedIndex(Index(site, PhononIndex{:u}(direction)), rcoordinate, icoordinate)

"""
    𝕡(direction) -> PhononIndex{:p}
    𝕡(site, direction) -> Index{<:PhononIndex{:p}}
    𝕡(site, direction, rcoordinate, icoordinate) -> CoordinatedIndex{<:Index{<:PhononIndex{:p}}}

Convenient construction of `SpinIndex{:p}`, `Index{<:SpinIndex{:p}}`, `CoordinatedIndex{<:Index{<:SpinIndex{:p}}}`.
"""
@inline 𝕡(direction) = PhononIndex{:p}(direction)
@inline 𝕡(site, direction) = Index(site, PhononIndex{:p}(direction))
@inline 𝕡(site, direction, rcoordinate, icoordinate) = CoordinatedIndex(Index(site, PhononIndex{:p}(direction)), rcoordinate, icoordinate)
@inline Base.getindex(::Type{OperatorIndex}, ::PhononIndex{:u}) = "𝕦"
@inline Base.getindex(::Type{OperatorIndex}, ::PhononIndex{:p}) = "𝕡"

"""
    kind(index::PhononIndex) -> Symbol
    kind(::Type{<:PhononIndex{K}}) where K -> Symbol

    kind(index::Index{<:PhononIndex}) -> Symbol
    kind(::Type{<:Index{<:PhononIndex{K}}}) where K -> Symbol

    kind(index::CoordinatedIndex{<:Index{<:PhononIndex}}) -> Symbol
    kind(::Type{<:CoordinatedIndex{<:Index{<:PhononIndex{K}}}}) where K -> Symbol

Get the kind of a phonon index.
"""
@inline kind(index::PhononIndex) = kind(typeof(index))
@inline kind(index::Index{<:PhononIndex}) = kind(typeof(index))
@inline kind(index::CoordinatedIndex{<:Index{<:PhononIndex}}) = kind(typeof(index))
@inline kind(::Type{<:PhononIndex}) = Symbol(":")
@inline kind(::Type{<:PhononIndex{K}}) where K = K
@inline kind(::Type{I}) where {I<:Index{<:PhononIndex}} = kind(internalindextype(I))
@inline kind(::Type{I}) where {I<:CoordinatedIndex{<:Index{<:PhononIndex}}} = kind(indextype(I))

## LaTeX format output
"""
    script(index::PhononIndex, ::Val{:direction}; kwargs...) -> String

Get the requested script of a phonon index.
"""
@inline script(index::PhononIndex{K, Colon}, ::Val{:direction}; kwargs...) where K = ":"
@inline script(index::PhononIndex, ::Val{:direction}; kwargs...) = string(index.direction)

@inline body(::PhononIndex{:u}) = "u"
@inline body(::PhononIndex{:p}) = "p"
@inline body(index::Index{<:PhononIndex}) = body(InternalIndex(index))
@inline body(index::CompositeIndex{<:Index{<:PhononIndex}}) = body(Index(index))
"""
    latexofphonons

Default LaTeX format of a phonon index.
"""
const latexofphonons = LaTeX{(:direction,), (:site,)}(body, "", "")
@inline latexname(::Type{<:PhononIndex}) = Symbol("PhononIndex")
@inline latexname(::Type{<:Index{<:PhononIndex}}) = Symbol("Index{PhononIndex}")
@inline latexname(::Type{<:CompositeIndex{<:Index{<:PhononIndex}}}) = Symbol("CompositeIndex{Index{PhononIndex}}")
latexformat(PhononIndex, latexofphonons)
latexformat(Index{<:PhononIndex}, latexofphonons)
latexformat(CompositeIndex{<:Index{<:PhononIndex}}, latexofphonons)

## Phonon
"""
    Phonon{K} <: SimpleInternal{PhononIndex{K, Char}}

Phonon space of lattice vibration generators at a single point.
"""
struct Phonon{K} <: SimpleInternal{PhononIndex{K, Char}}
    ndirection::Int
    function Phonon{K}(ndirection::Integer) where K
        @assert K∈(:u, :p, :) "Phonon error: wrong kind($K)."
        @assert ndirection∈(1, 2, 3) "Phonon error: wrong number of directions."
        new{K}(ndirection)
    end
end
@inline shape(pn::Phonon) = (1:pn.ndirection,)
@inline Base.convert(::Type{<:CartesianIndex}, index::PhononIndex{K, Char}, ::Phonon{K}) where K = CartesianIndex(Int(index.direction)-Int('x')+1)
@inline Base.convert(::Type{<:PhononIndex}, index::CartesianIndex{1}, pn::Phonon) = PhononIndex{kind(pn)}(Char(Int('x')+index[1]-1))
@inline Base.summary(io::IO, pn::Phonon) = @printf io "%s-element Phonon{%s}" length(pn) repr(kind(pn))
@inline Base.show(io::IO, pn::Phonon{K}) where K = @printf io "%s%s(%s)" pn|>typeof|>nameof (K==Colon() ? "" : "{$(repr(K))}") join(("$name=$(getfield(pn, name))" for name in pn|>typeof|>fieldnames), ", ")
@inline Base.show(io::IO, ::Type{Phonon{:}}) = @printf io "%s" "Phonon{:}"
@inline Base.match(::Type{<:PhononIndex{K}}, ::Type{<:Phonon{:}}) where K = true
@inline Base.match(::Type{<:PhononIndex{K}}, ::Type{<:Phonon{K}}) where K = true
@inline Base.match(::Type{<:PhononIndex{K₁}}, ::Type{<:Phonon{K₂}}) where {K₁, K₂} = false
@inline Base.filter(::Type{<:PhononIndex{K}}, ph::Phonon{:}) where K = Phonon{K}(ph.ndirection)
@inline Base.filter(::Type{<:PhononIndex{K}}, ::Type{Phonon{:}}) where K = Phonon{K}
### requested by ConstrainedInternal
@inline shape(pn::Phonon, index::PhononIndex) = (phononshape(index.direction, pn.ndirection),)
@inline phononshape(::Union{Symbol, Colon}, n::Int) = 1:n
@inline function phononshape(v::Char, n::Int)
    index = Int(v) - Int('x') + 1
    @assert 0<index<n+1 "shape error: out of range."
    return index:index
end

"""
    Phonon(ndirection::Integer)
    Phonon{K}(ndirection::Integer) where K

Construct a phonon space.
"""
@inline Phonon(ndirection::Integer) = Phonon{:}(ndirection)

"""
    kind(pn::Phonon)
    kind(::Type{<:Phonon{K}}) where K

Get the kind of a phonon space.
"""
@inline kind(pn::Phonon) = kind(typeof(pn))
@inline kind(::Type{<:Phonon{K}}) where K = K

## Permutation
"""
    permute(id₁::PhononIndex, id₂::PhononIndex) -> ZeroAtLeast{Operator}

Permute two phonon indexes and get the result.
"""
function permute(id₁::PhononIndex, id₂::PhononIndex)
    k₁, k₂ = kind(id₁), kind(id₂)
    if id₁.direction==id₂.direction && k₁≠k₂
        return (Operator(k₁==:u ? 1im : -1im), Operator(1, id₂, id₁))
    else
        return (Operator(1, id₂, id₁),)
    end
end

## Coupling
### MatrixCoupling
"""
    MatrixCoupling{P}(sites::Union{NTuple{2, Ordinal}, Colon}, matrix::AbstractMatrix; rows::Union{AbstractVector, Nothing}=nothing, cols::Union{AbstractVector, Nothing}=nothing) where {P<:PhononIndex{:u}}

Construct a set of `Coupling`s corresponding to the dynamical matrix of phonons.
"""
function MatrixCoupling{P}(sites::Union{NTuple{2, Ordinal}, Colon}, matrix::AbstractMatrix; rows::Union{AbstractVector, Nothing}=nothing, cols::Union{AbstractVector, Nothing}=nothing) where {P<:PhononIndex{:u}}
    @assert size(matrix)[1]∈(1, 2, 3) && size(matrix)[2]∈(1, 2, 3) "MatrixCoupling error: mismatched dimension of input matrix."
    isnothing(rows) && (rows = size(matrix)[1]==1 ? SVector('x') : size(matrix)[1]==2 ? SVector('x', 'y') : SVector('x', 'y', 'z'))
    isnothing(cols) && (cols = size(matrix)[2]==1 ? SVector('x') : size(matrix)[2]==2 ? SVector('x', 'y') : SVector('x', 'y', 'z'))
    @assert size(matrix)==(length(rows), length(cols)) "MatrixCoupling error: mismatched input matrix and rows/cols."
    return MatrixCoupling{P}(sites, MatrixCouplingComponent(rows, cols, matrix))
end

"""
    𝕦ᵀ𝕦(sites::Union{NTuple{2, Ordinal}, Colon}, matrix::AbstractMatrix; rows::Union{AbstractVector, Nothing}=nothing, cols::Union{AbstractVector, Nothing}=nothing) -> MatrixCoupling

Construct a set of `Coupling`s corresponding to the dynamical matrix of phonons.
"""
@inline function 𝕦ᵀ𝕦(sites::Union{NTuple{2, Ordinal}, Colon}, matrix::AbstractMatrix; rows::Union{AbstractVector, Nothing}=nothing, cols::Union{AbstractVector, Nothing}=nothing)
    return MatrixCoupling{PhononIndex{:u}}(sites, matrix; rows=rows, cols=cols)
end

### expand
"""
    expand(pnc::Coupling{<:Number, <:Pattern{<:NTuple{2, Index{<:PhononIndex{:u}}}}}, ::Val{:Hooke}, bond::Bond, hilbert::Hilbert) -> VectorSpace

Expand the default phonon potential coupling on a given bond.
"""
function expand(pnc::Coupling{<:Number, <:Pattern{<:NTuple{2, Index{<:PhononIndex{:u}}}}}, ::Val{:Hooke}, bond::Bond, hilbert::Hilbert)
    R̂ = rcoordinate(bond)/norm(rcoordinate(bond))
    @assert isapprox(pnc.value, 1, atol=atol, rtol=rtol) "expand error: wrong coefficient of Hooke coupling."
    @assert isa(pnc.pattern.indexes.sites, NTuple{2, Colon}) "expand error: the `:sites` attributes of the Hooke coupling pattern must be a 2-tuple of colons."
    pn₁ = filter(InternalIndex(pnc.pattern.indexes[1]), hilbert[bond[1].site])
    pn₂ = filter(InternalIndex(pnc.pattern.indexes[2]), hilbert[bond[2].site])
    @assert pn₁.ndirection==pn₂.ndirection==length(R̂) "expand error: mismatched number of directions."
    return PPExpand(R̂, (bond[1], bond[2]))
end
struct PPExpand{N, D<:Number} <: VectorSpace{Operator{D, ZeroAtLeast{CoordinatedIndex{Index{PhononIndex{:u, Char}, Int}, SVector{N, D}}, 2}}}
    direction::SVector{N, D}
    points::NTuple{2, Point{N, D}}
end
@inline VectorSpaceStyle(::Type{<:PPExpand}) = VectorSpaceDirectProducted(:forward)
@inline shape(pnce::PPExpand) = (1:length(pnce.direction), 1:length(pnce.direction), 1:4)
function Base.convert(::Type{<:Operator}, index::CartesianIndex{3}, pnce::PPExpand)
    direction₁ = Char(Int('x')+index[1]-1)
    direction₂ = Char(Int('x')+index[2]-1)
    coeff = index[3]∈(1, 4) ? 1 : -1
    pos₁, pos₂ = index[3]==1 ? (1, 1) : index[3]==2 ? (1, 2) : index[3]==3 ? (2, 1) : (2, 2)
    index₁ = 𝕦(pnce.points[pos₁].site, direction₁, pnce.points[pos₁].rcoordinate, pnce.points[pos₁].icoordinate)
    index₂ = 𝕦(pnce.points[pos₂].site, direction₂, pnce.points[pos₂].rcoordinate, pnce.points[pos₂].icoordinate)
    return Operator(pnce.direction[index[1]]*pnce.direction[index[2]]*coeff, index₁, index₂)
end

## Term
"""
    Kinetic(id::Symbol, value; amplitude::Union{Function, Nothing}=nothing, ismodulatable::Bool=true)

Kinetic energy of phonons.

Type alias for `Term{:Kinetic, id, V, Int, C<:TermCoupling, A<:TermAmplitude}`.
"""
const Kinetic{id, V, C<:TermCoupling, A<:TermAmplitude} = Term{:Kinetic, id, V, Int, C, A}
@inline function Kinetic(id::Symbol, value; amplitude::Union{Function, Nothing}=nothing, ismodulatable::Bool=true)
    return Term{:Kinetic}(id, value, 0, Coupling(𝕡(:, :), 𝕡(:, :)), true; amplitude=amplitude, ismodulatable=ismodulatable)
end

"""
    Hooke(id::Symbol, value, bondkind; amplitude::Union{Function, Nothing}=nothing, ismodulatable::Bool=true)

Potential energy of phonons by the Hooke's law.

Type alias for `Term{:Hooke, id, V, B, C<:TermCoupling, A<:TermAmplitude}`
"""
const Hooke{id, V, B, C<:TermCoupling, A<:TermAmplitude} = Term{:Hooke, id, V, B, C, A}
@inline function Hooke(id::Symbol, value, bondkind; amplitude::Union{Function, Nothing}=nothing, ismodulatable::Bool=true)
    return Term{:Hooke}(id, value, bondkind, Coupling(𝕦(:, :), 𝕦(:, :)), true; amplitude=amplitude, ismodulatable=ismodulatable)
end

"""
    Elastic(id::Symbol, value, bondkind, coupling; amplitude::Union{Function, Nothing}=nothing, ismodulatable::Bool=true)

Generic elastic energy of phonons.

Type alias for `Term{:Elastic, id, V, B, C<:TermCoupling, A<:TermAmplitude}`
"""
const Elastic{id, V, B, C<:TermCoupling, A<:TermAmplitude} = Term{:Elastic, id, V, B, C, A}
@inline function Elastic(id::Symbol, value, bondkind, coupling; amplitude::Union{Function, Nothing}=nothing, ismodulatable::Bool=true)
    value, factor = promote(value, 1//2)
    return Term{:Elastic, id}(value, bondkind, TermCoupling(coupling), TermAmplitude(amplitude), true, ismodulatable, factor)
end
function expand!(operators::Operators, term::Elastic, bond::Bond, hilbert::Hilbert; half::Bool=false)
    argtypes = Tuple{Operators, Term, Bond, Hilbert}
    invoke(expand!, argtypes, operators, term, bond, hilbert; half=half)
    invoke(expand!, argtypes, operators, term, reverse(bond), hilbert; half=half)
    return operators
end

"""
    PhononTerm

Type alias for `Union{Kinetic, Hooke, Elastic}`.
"""
const PhononTerm = Union{Kinetic, Hooke, Elastic}

end # module
