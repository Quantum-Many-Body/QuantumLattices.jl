module QuantumSystems

using LinearAlgebra: dot, ishermitian, norm
using Printf: @printf, @sprintf
using SparseArrays: SparseMatrixCSC
using StaticArrays: SMatrix, SVector
using ..DegreesOfFreedom: AbstractIndex, CompositeIndex, Component, InternalIndexProd, CoordinatedIndex, CompositeInternal, Coupling, Hilbert, ConstrainedInternal, Index, Ordinal, SimpleInternalIndex, SimpleInternal, Term, TermAmplitude, TermCoupling, @pattern
using ..QuantumLattices: decompose, dtype, kind
using ..QuantumOperators: ID, LaTeX, Operator, OperatorProd, Operators, latexformat
using ..Spatials: Bond, Point, direction, isparallel, rcoordinate
using ..Toolkit: atol, efficientoperations, rtol, Float, VectorSpace, VectorSpaceCartesian, VectorSpaceStyle, delta, getcontent, rawtype, tostr

import ..DegreesOfFreedom: MatrixCoupling, allequalfields, indextype, isdefinite, patternrule, statistics
import ..QuantumLattices: ⊗, expand, expand!, permute, rank
import ..QuantumOperators: latexname, matrix, optype, script
import ..Toolkit: shape

# Canonical complex fermionic/bosonic systems
export 𝕓, 𝕗, 𝔽, annihilation, creation, latexofbosons, latexoffermions, latexofparticles, @σ_str, @L_str
export Coulomb, FockIndex, Fock, FockTerm, Hopping, Hubbard, InterOrbitalInterSpin, InterOrbitalIntraSpin, Onsite, PairHopping, Pairing, SpinFlip, isannihilation, iscreation, isnormalordered

# SU(2) spin systems
export latexofspins, SID, Spin, totalspin, @Γ_str, @Γ′_str, @DM_str, @Heisenberg_str, @Ising_str
export DM, Heisenberg, Ising, Kitaev, SingleIonAnisotropy, SpinTerm, Zeeman, Γ, Γ′

# Phononic systems
export latexofphonons, Elastic, PID, Phonon, Kinetic, Hooke, PhononTerm

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
    FockIndex{T, O<:Union{Int, Symbol, Colon}, S<:Union{Rational{Int}, Symbol, Colon}, N<:Union{Int, Symbol, Colon}} <: SimpleInternalIndex

Fock index, i.e., the internal index of a Fock space.
"""
struct FockIndex{T, O<:Union{Int, Symbol, Colon}, S<:Union{Rational{Int}, Symbol, Colon}, N<:Union{Int, Symbol, Colon}} <: SimpleInternalIndex
    orbital::O
    spin::S
    nambu::N
    function FockIndex{T}(orbital::Union{Int, Symbol, Colon}, spin::Union{Rational{Int}, Int, Symbol, Colon}, nambu::Union{Int, Symbol, Colon}) where T
        @assert T∈(:f, :b, :) "FockIndex error: wrong statistics."
        isa(spin, Rational{Int}) && @assert spin.den∈(1, 2) "FockIndex error: wrong spin."
        isa(nambu, Int) && @assert nambu∈(annihilation, creation) "FockIndex error: wrong input nambu($nambu)."
        isa(spin, Int) && (spin = convert(Rational{Int}, spin))
        new{T, typeof(orbital), typeof(spin), typeof(nambu)}(orbital, spin, nambu)
    end
end
### basic methods of concrete SimpleInternalIndex
@inline statistics(::Type{<:FockIndex{T}}) where T = T
@inline isdefinite(::Type{<:FockIndex{T, Int, Rational{Int}} where T}) = true
@inline Base.:(==)(index₁::FockIndex, index₂::FockIndex) = statistics(index₁)==statistics(index₂) && ==(efficientoperations, index₁, index₂)
@inline Base.isequal(index₁::FockIndex, index₂::FockIndex) = isequal(statistics(index₁), statistics(index₂)) && isequal(efficientoperations, index₁, index₂)
@inline Base.hash(index::FockIndex, h::UInt) = hash((statistics(index), index.orbital, index.spin, index.nambu), h)
@inline Base.adjoint(index::FockIndex{T, <:Union{Int, Symbol, Colon}, <:Union{Rational{Int}, Symbol, Colon}, Int}) where T = FockIndex{T}(index.orbital, index.spin, 3-index.nambu)
@inline @generated function Base.replace(index::FockIndex; kwargs...)
    exprs = [:(get(kwargs, $name, getfield(index, $name))) for name in QuoteNode.(fieldnames(index))]
    return :(rawtype(typeof(index)){statistics(index)}($(exprs...)))
end
### requested by InternalPattern
@inline allequalfields(::Type{<:FockIndex}) = (:orbital, :spin)
### requested by MatrixCouplingProd
@inline indextype(::Type{FockIndex}, ::Type{O}, ::Type{S}, ::Type{N}) where {O<:Union{Int, Symbol, Colon}, S<:Union{Rational{Int}, Symbol, Colon}, N<:Union{Int, Symbol, Colon}} = FockIndex{:, O, S, N}
@inline indextype(::Type{FockIndex{T}}, ::Type{O}, ::Type{S}, ::Type{N}) where {T, O<:Union{Int, Symbol, Colon}, S<:Union{Rational{Int}, Symbol, Colon}, N<:Union{Int, Symbol, Colon}} = FockIndex{T, O, S, N}

"""
    FockIndex(orbital::Union{Int, Symbol, Colon}, spin::Union{Rational{Int}, Int, Symbol, Colon}, nambu::Union{Int, Symbol, Colon})
    FockIndex{T}(orbital::Union{Int, Symbol, Colon}, spin::Union{Rational{Int}, Int, Symbol, Colon}, nambu::Union{Int, Symbol, Colon}) where T
    FockIndex{T, O, S, N}(orbital::Union{Int, Symbol, Colon}, spin::Union{Rational{Int}, Symbol, Colon}, nambu::Union{Int, Symbol, Colon}) where {T, O, S, N}

Construct a Fock index.
"""
@inline FockIndex(orbital::Union{Int, Symbol, Colon}, spin::Union{Rational{Int}, Int, Symbol, Colon}, nambu::Union{Int, Symbol, Colon}) = FockIndex{:}(orbital, spin, nambu)
@inline FockIndex{T, O, S, N}(orbital::Union{Int, Symbol, Colon}, spin::Union{Rational{Int}, Symbol, Colon}, nambu::Union{Int, Symbol, Colon}) where {T, O, S, N} = FockIndex{T}(orbital, spin, nambu)

### convenient construction and string representation 
"""
    𝕗(orbital, spin, nambu) -> FockIndex{:f}
    𝕗(site, orbital, spin, nambu) -> Index{<:FockIndex{:f}}
    𝕗(site, orbital, spin, nambu, rcoordinate, icoordinate) -> CoordinatedIndex{<:Index{<:FockIndex{:f}}}

    𝕓(orbital, spin, nambu) -> FockIndex{:b}
    𝕓(site, orbital, spin, nambu) -> Index{<:FockIndex{:b}}
    𝕓(site, orbital, spin, nambu, rcoordinate, icoordinate) -> CoordinatedIndex{<:Index{<:FockIndex{:b}}}

    𝔽(orbital, spin, nambu) -> FockIndex{:}
    𝔽(site, orbital, spin, nambu) -> Index{<:FockIndex{:}}
    𝔽(site, orbital, spin, nambu, rcoordinate, icoordinate) -> CoordinatedIndex{<:Index{<:FockIndex{:}}}

Convenient construction of `FockIndex`, `Index{<:FockIndex}`, `CoordinatedIndex{<:Index{<:FockIndex}}`.
"""
const _fock_ = (:𝕗, :𝕓, :𝔽), (QuoteNode(:f), QuoteNode(:b), :)
for (name, statistics) in zip(_fock_...)
    @eval @inline $name(orbital, spin, nambu) = FockIndex{$statistics}(orbital, spin, nambu)
    @eval @inline $name(site, orbital, spin, nambu) = Index(site, FockIndex{$statistics}(orbital, spin, nambu))
    @eval @inline $name(site, orbital, spin, nambu, rcoordinate, icoordinate) = CoordinatedIndex(Index(site, FockIndex{$statistics}(orbital, spin, nambu)), rcoordinate, icoordinate)
end
@inline Base.getindex(::Type{AbstractIndex}, ::Type{I}) where {I<:Union{FockIndex{:f}, Index{<:FockIndex{:f}}, CoordinatedIndex{<:Index{<:FockIndex{:f}}}}} = 𝕗
@inline Base.getindex(::Type{AbstractIndex}, ::Type{I}) where {I<:Union{FockIndex{:b}, Index{<:FockIndex{:b}}, CoordinatedIndex{<:Index{<:FockIndex{:b}}}}} = 𝕓
@inline Base.getindex(::Type{AbstractIndex}, ::Type{I}) where {I<:Union{FockIndex{:}, Index{<:FockIndex{:}}, CoordinatedIndex{<:Index{<:FockIndex{:}}}}} = 𝔽
@inline Base.getindex(::Type{AbstractIndex}, ::typeof(𝕗)) = FockIndex{:f}
@inline Base.getindex(::Type{AbstractIndex}, ::typeof(𝕓)) = FockIndex{:b}
@inline Base.getindex(::Type{AbstractIndex}, ::typeof(𝔽)) = FockIndex{:}

"""
    isannihilation(index::FockIndex) -> Bool
    isannihilation(index::Index) -> Bool
    isannihilation(index::CompositeIndex) -> Bool

Judge whether the nambu index is `annihilation`.
"""
@inline isannihilation(index::FockIndex) = index.nambu==annihilation
@inline isannihilation(index::Index) = isannihilation(index.internal)
@inline isannihilation(index::CompositeIndex) = isannihilation(getcontent(index, :index))

"""
    iscreation(index::FockIndex) -> Bool
    iscreation(index::Index) -> Bool
    iscreation(index::CompositeIndex) -> Bool

Judge whether the nambu index is `creation`.
"""
@inline iscreation(index::FockIndex) = index.nambu==creation
@inline iscreation(index::Index) = iscreation(index.internal)
@inline iscreation(index::CompositeIndex) = iscreation(getcontent(index, :index))

### patternrule
"""
    patternrule(nambus::NTuple{N, Colon}, ::Val{}, ::Type{<:FockIndex}, ::Val{:nambu}) where N -> NTuple{N, Int}

Default pattern rule for the `:nambu` attribute of Fock indexes.
"""
@inline patternrule(::NTuple{N, Colon}, ::Val{}, ::Type{<:FockIndex}, ::Val{:nambu}) where N = ntuple(i->isodd(i) ? creation : annihilation, Val(N))

## LaTeX format output
"""
    script(index::FockIndex, ::Val{:orbital}; kwargs...) -> String
    script(index::FockIndex, ::Val{:spin}; kwargs...) -> String
    script(index::FockIndex, ::Val{:spinsym}; kwargs...) -> String
    script(index::FockIndex, ::Val{:nambu}; kwargs...) -> String

Get the requested script of a Fock index.
"""
@inline script(index::FockIndex, ::Val{:orbital}; kwargs...) = tostr(index.orbital)
@inline script(index::FockIndex, ::Val{:spin}; kwargs...) = tostr(index.spin)
@inline function script(index::FockIndex, ::Val{:spinsym}; kwargs...)
    nspin = get(kwargs, :nspin, nothing)
    return (isnothing(nspin) || nspin==2) && index.spin==-1//2 ? "↓" : (isnothing(nspin) || nspin==2) && index.spin==1//2 ? "↑" : (isnothing(nspin) || nspin==1) && index.spin==0 ? "" : tostr(index.spin)
end
@inline script(index::FockIndex, ::Val{:nambu}; kwargs...) = iscreation(index) ? "\\dagger" : isannihilation(index) ? "" : tostr(index.nambu)

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
const latexofparticles = LaTeX{(:nambu,), (:site, :orbital, :spinsym)}('f')
@inline latexname(::Type{<:Index{<:FockIndex{:}}}) = Symbol("Index{FockIndex}")
@inline latexname(::Type{<:CompositeIndex{<:Index{<:FockIndex{:}}}}) = Symbol("CompositeIndex{Index{FockIndex}}")
@inline latexname(::Type{<:FockIndex{:}}) = Symbol("FockIndex")
latexformat(FockIndex{:}, latexofparticles)
latexformat(Index{<:FockIndex{:}}, latexofparticles)
latexformat(CompositeIndex{<:Index{<:FockIndex{:}}}, latexofparticles)

## Fock
"""
    Fock{T} <: SimpleInternal{FockIndex{T, Int, Rational{Int}, Int}}

Fock space at a single point.
"""
struct Fock{T} <: SimpleInternal{FockIndex{T, Int, Rational{Int}, Int}}
    norbital::Int
    nspin::Int
    function Fock{T}(norbital::Int, nspin::Int) where T
        @assert T∈(:f, :b) "Fock error: wrong statistics."
        new{T}(norbital, nspin)
    end
end
@inline shape(fock::Fock) = (1:fock.norbital, 1:fock.nspin, 1:2)
@inline Base.eltype(::Type{Fock}) = (FockIndex{T, Int, Rational{Int}, Int} where T)
@inline Base.convert(::Type{<:CartesianIndex}, index::FockIndex{T}, fock::Fock{T}) where T = CartesianIndex(index.orbital, Int(index.spin+(fock.nspin-1)//2)+1, index.nambu)
@inline Base.convert(::Type{<:FockIndex}, index::CartesianIndex{3}, fock::Fock{T}) where T = FockIndex{T}(index[1], index[2]-1-(fock.nspin-1)//2, index[3])
@inline Base.summary(io::IO, fock::Fock) = @printf io "%s-element Fock{%s}" length(fock) repr(statistics(fock))
@inline Base.show(io::IO, fock::Fock) = @printf io "%s{%s}(%s)" fock|>typeof|>nameof repr(statistics(fock)) join(("$name=$(getfield(fock, name))" for name in fock|>typeof|>fieldnames), ", ")
@inline Base.match(::Type{<:FockIndex{:}}, ::Type{<:Fock{T}}) where T = true
@inline Base.match(::Type{<:FockIndex{T}}, ::Type{<:Fock{T}}) where T = true
@inline Base.match(::Type{<:FockIndex{T₁}}, ::Type{<:Fock{T₂}}) where {T₁, T₂} = false
### requested by ConstrainedInternal
@inline function shape(internal::Fock, index::FockIndex{T, <:Union{Int, Symbol, Colon}, <:Union{Rational{Int}, Symbol, Colon}, Int}) where T
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
    isnormalordered(opt::Operator{<:Number, <:ID{CompositeIndex{<:Index{<:FockIndex}}}}) -> Bool

Judge whether an operator is normal ordered.
"""
function isnormalordered(opt::Operator{<:Number, <:ID{CompositeIndex{<:Index{<:FockIndex}}}})
    flag = true
    for i = 1:rank(opt)
        flag && (opt[i].index.internal.nambu == annihilation) && (flag = false)
        flag || (opt[i].index.internal.nambu == creation) && return false
    end
    return true
end

"""
    *(f₁::Operator{<:Number, <:ID{FockIndex{:f}}}, f₂::Operator{<:Number, <:ID{FockIndex{:f}}}) -> Operator
    *(f₁::Operator{<:Number, <:ID{Index{<:FockIndex{:f}, Int}}}, f₂::Operator{<:Number, <:ID{Index{<:FockIndex{:f}, Int}}}) -> Operator
    *(f₁::Operator{<:Number, <:ID{CompositeIndex{<:Index{<:FockIndex{:f}, Int}}}}, f₂::Operator{<:Number, <:ID{CompositeIndex{<:Index{<:FockIndex{:f}, Int}}}}) -> Operator

Get the multiplication of two fermionic Fock operators.
"""
const block = quote
    result = invoke(*, Tuple{OperatorProd, OperatorProd}, f₁, f₂)
    rank(f₁)>0 && rank(f₂)>0 && f₁[end]==f₂[1] && return replace(result, value=zero(valtype(result)))
    return result
end
@eval @inline Base.:*(f₁::Operator{<:Number, <:ID{FockIndex{:f}}}, f₂::Operator{<:Number, <:ID{FockIndex{:f}}}) = $block
@eval @inline Base.:*(f₁::Operator{<:Number, <:ID{Index{<:FockIndex{:f}, Int}}}, f₂::Operator{<:Number, <:ID{Index{<:FockIndex{:f}, Int}}}) = $block
@eval @inline Base.:*(f₁::Operator{<:Number, <:ID{CompositeIndex{<:Index{<:FockIndex{:f}, Int}}}}, f₂::Operator{<:Number, <:ID{CompositeIndex{<:Index{<:FockIndex{:f}, Int}}}}) = $block
@inline Base.:*(m₁::Operator{<:Number, Tuple{}}, m₂::Operator{<:Number, Tuple{}}) = Operator(m₁.value*m₂.value)

## Permutation
"""
    permute(id₁::FockIndex, id₂::FockIndex) -> Tuple{Vararg{Operator}}

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
    if id₁'==id₂
        return (Operator(iscreation(id₁) ? -1 : +1), Operator(1, id₂, id₁))
    else
        return (Operator(1, id₂, id₁),)
    end
end
@inline permute(id₁::FockIndex{:b}, id₂::FockIndex{:f}) = (Operator(1, id₂, id₁),)
@inline permute(id₁::FockIndex{:f}, id₂::FockIndex{:b}) = (Operator(1, id₂, id₁),)

## Coupling
### MatrixCoupling
const default_matrix = SparseMatrixCSC(hcat(1))
"""
    MatrixCoupling(
        sites::Union{NTuple{2, Ordinal}, Colon},
        F::Union{Type{<:FockIndex}, typeof(𝕗), typeof(𝕓), typeof(𝔽)},
        orbital::Union{AbstractMatrix, Colon},
        spin::Union{AbstractMatrix, Colon},
        nambu::Union{AbstractMatrix, Colon}
    )

Construct a matrix coupling for Fock systems.
"""
@inline function MatrixCoupling(sites::Union{NTuple{2, Ordinal}, Colon}, F::Union{typeof(𝕗), typeof(𝕓), typeof(𝔽)}, orbital::Union{AbstractMatrix, Colon}, spin::Union{AbstractMatrix, Colon}, nambu::Union{AbstractMatrix, Colon})
    return MatrixCoupling(sites, AbstractIndex[F], orbital, spin, nambu)
end
@inline function MatrixCoupling(sites::Union{NTuple{2, Ordinal}, Colon}, ::Type{F}, orbital::Union{AbstractMatrix, Colon}, spin::Union{AbstractMatrix, Colon}, nambu::Union{AbstractMatrix, Colon}) where {F<:FockIndex}
    return MatrixCoupling(sites, F, Component(F, Val(:orbital), orbital), Component(F, Val(:spin), spin), Component(F, Val(:nambu), nambu))
end
@inline Component(::Type{<:FockIndex}, ::Val, ::Colon) = Component(SVector(:), SVector(:), default_matrix)
@inline Component(::Type{<:FockIndex}, ::Val{:orbital}, matrix::AbstractMatrix) = (check(matrix); Component(1:size(matrix)[1], 1:size(matrix)[2], matrix))
@inline Component(::Type{<:FockIndex}, ::Val{:spin}, matrix::AbstractMatrix) = (check(matrix); Component((size(matrix)[1]-1)//2:-1:(1-size(matrix)[1])//2, (size(matrix)[2]-1)//2:-1:(1-size(matrix)[2])//2, matrix))
@inline Component(::Type{<:FockIndex}, ::Val{:nambu}, matrix::AbstractMatrix) = (check(matrix); @assert size(matrix)==(2, 2) "Component error: for nambu subspace, the input matrix must be 2×2."; Component(1:1:2, 2:-1:1, matrix))
@inline check(matrix::AbstractMatrix) = @assert map(firstindex, axes(matrix))==(1, 1) "Component error: matrix check fails."

### Pauli matrices
"""
    σ"0" => SparseMatrixCSC([1 0; 0 1])
    σ"x" => SparseMatrixCSC([0 1; 1 0])
    σ"y" => SparseMatrixCSC([0 -1im; 1im 0])
    σ"z" => SparseMatrixCSC([1 0; 0 -1])
    σ"+" => SparseMatrixCSC([0 1; 0 0])
    σ"-" => SparseMatrixCSC([0 0; 1 0])
    σ"11" => SparseMatrixCSC([1 0; 0 0])
    σ"22" => SparseMatrixCSC([0 0; 0 1])

The Pauli matrix σ⁰, σˣ, σʸ, σᶻ, σ⁺, σ⁻, σ¹¹, σ²².
"""
macro σ_str(str::String)
    str=="0" && return SparseMatrixCSC([1 0; 0 1])
    str=="x" && return SparseMatrixCSC([0 1; 1 0])
    str=="y" && return SparseMatrixCSC([0 -1im; 1im 0])
    str=="z" && return SparseMatrixCSC([1 0; 0 -1])
    str=="+" && return SparseMatrixCSC([0 1; 0 0])
    str=="-" && return SparseMatrixCSC([0 0; 1 0])
    str=="11" && return SparseMatrixCSC([1 0; 0 0])
    str=="22" && return SparseMatrixCSC([0 0; 0 1])
    error("@σ_str error: wrong input string.")
end

### Rotation matrices
"""
    L"x" => SparseMatrixCSC([0 0 0; 0 0 1im; 0 -1im 0])
    L"y" => SparseMatrixCSC([0 0 -1im; 0 0 0; 1im 0 0])
    L"z" => SparseMatrixCSC([0 1im 0; -1im 0 0; 0 0 0])

The three-dimensional rotation generators.
"""
macro L_str(str::String)
    str=="x" && return SparseMatrixCSC([0 0 0; 0 0 1im; 0 -1im 0])
    str=="y" && return SparseMatrixCSC([0 0 -1im; 0 0 0; 1im 0 0])
    str=="z" && return SparseMatrixCSC([0 1im 0; -1im 0 0; 0 0 0])
    error("@L_str error: wrong input string.")
end

## Term
"""
    Onsite(id::Symbol, value, coupling=Coupling(𝔽(:, :, :, :), 𝔽(:, :, :, :)); ishermitian::Bool=true, amplitude::Union{Function, Nothing}=nothing, ismodulatable::Bool=true)

Onsite term.

Type alias for `Term{:Onsite, id, V, Int, C<:TermCoupling, A<:TermAmplitude}`.
"""
const Onsite{id, V, C<:TermCoupling, A<:TermAmplitude} = Term{:Onsite, id, V, Int, C, A}
@inline function Onsite(id::Symbol, value, coupling=Coupling(𝔽(:, :, :, :), 𝔽(:, :, :, :)); ishermitian::Bool=true, amplitude::Union{Function, Nothing}=nothing, ismodulatable::Bool=true)
    return Term{:Onsite}(id, value, 0, coupling, ishermitian; amplitude=amplitude, ismodulatable=ismodulatable)
end

"""
    Hopping(id::Symbol, value, bondkind, coupling=Coupling(𝔽(:, :, :, :), 𝔽(:, :, :, :)); amplitude::Union{Function, Nothing}=nothing, ismodulatable::Bool=true)

Hopping term.

Type alias for `Term{:Hopping, id, V, B, C<:TermCoupling, A<:TermAmplitude}`.
"""
const Hopping{id, V, B, C<:TermCoupling, A<:TermAmplitude} = Term{:Hopping, id, V, B, C, A}
@inline function Hopping(id::Symbol, value, bondkind, coupling=Coupling(𝔽(:, :, :, :), 𝔽(:, :, :, :)); amplitude::Union{Function, Nothing}=nothing, ismodulatable::Bool=true)
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
@inline patternrule(::NTuple{2, Colon}, ::Val{:Pairing}, ::Type{<:FockIndex}, ::Val{:nambu}) = (annihilation, annihilation)
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
    return Term{:Hubbard}(id, value, 0, Coupling(:, 𝔽, :, (1//2, 1//2, -1//2, -1//2), (2, 1, 2, 1)), true; amplitude=amplitude, ismodulatable=ismodulatable)
end

"""
    InterOrbitalInterSpin(id::Symbol, value; amplitude::Union{Function, Nothing}=nothing, ismodulatable::Bool=true)

Interorbital-interspin term.

Type alias for `Term{:InterOrbitalInterSpin, id, V, Int, C<:TermCoupling, A<:TermAmplitude}`.
"""
const InterOrbitalInterSpin{id, V, C<:TermCoupling, A<:TermAmplitude} = Term{:InterOrbitalInterSpin, id, V, Int, C, A}
@inline function InterOrbitalInterSpin(id::Symbol, value; amplitude::Union{Function, Nothing}=nothing, ismodulatable::Bool=true)
    return Term{:InterOrbitalInterSpin}(
        id, value, 0, Coupling(@pattern(𝔽(:, α, σ, 2), 𝔽(:, α, σ, 1), 𝔽(:, β, σ′, 2), 𝔽(:, β, σ′, 1); constraint=α<β && σ≠σ′)), true;
        amplitude=amplitude, ismodulatable=ismodulatable
    )
end

"""
    InterOrbitalIntraSpin(id::Symbol, value; amplitude::Union{Function, Nothing}=nothing, ismodulatable::Bool=true)

Interorbital-intraspin term.

Type alias for `Term{:InterOrbitalIntraSpin, id, V, Int, C<:TermCoupling, A<:TermAmplitude}`.
"""
const InterOrbitalIntraSpin{id, V, C<:TermCoupling, A<:TermAmplitude} = Term{:InterOrbitalIntraSpin, id, V, Int, C, A}
@inline function InterOrbitalIntraSpin(id::Symbol, value; amplitude::Union{Function, Nothing}=nothing, ismodulatable::Bool=true)
    return Term{:InterOrbitalIntraSpin}(
        id, value, 0, Coupling(@pattern(𝔽(:, α, σ, 2), 𝔽(:, α, σ, 1), 𝔽(:, β, σ, 2), 𝔽(:, β, σ, 1); constraint=α<β)), true;
        amplitude=amplitude, ismodulatable=ismodulatable
    )
end

"""
    SpinFlip(id::Symbol, value; amplitude::Union{Function, Nothing}=nothing, ismodulatable::Bool=true)

Spin-flip term.

Type alias for `Term{:SpinFlip, id, V, Int, C<:TermCoupling, A<:TermAmplitude}`.
"""
const SpinFlip{id, V, C<:TermCoupling, A<:TermAmplitude} = Term{:SpinFlip, id, V, Int, C, A}
@inline function SpinFlip(id::Symbol, value; amplitude::Union{Function, Nothing}=nothing, ismodulatable::Bool=true)
    return Term{:SpinFlip}(
        id, value, 0, Coupling(@pattern(𝔽(:, α, 1//2, 2), 𝔽(:, β, -1//2, 2), 𝔽(:, α, -1//2, 1), 𝔽(:, β, 1//2, 1); constraint=α<β)), false;
        amplitude=amplitude, ismodulatable=ismodulatable
    )
end

"""
    PairHopping(id::Symbol, value; amplitude::Union{Function, Nothing}=nothing, ismodulatable::Bool=true)

Pair-hopping term.

Type alias for `Term{:PairHopping, id, V, Int, C<:TermCoupling, A<:TermAmplitude}`.
"""
const PairHopping{id, V, C<:TermCoupling, A<:TermAmplitude} = Term{:PairHopping, id, V, Int, C, A}
@inline function PairHopping(id::Symbol, value; amplitude::Union{Function, Nothing}=nothing, ismodulatable::Bool=true)
    return Term{:PairHopping}(
        id, value, 0, Coupling(@pattern(𝔽(:, α, 1//2, 2), 𝔽(:, α, -1//2, 2), 𝔽(:, β, -1//2, 1), 𝔽(:, β, 1//2, 1); constraint=α<β)), false;
        amplitude=amplitude, ismodulatable=ismodulatable
    )
end

"""
    Coulomb(
        id::Symbol, value, bondkind, coupling=Coupling(𝔽(:, :, :, :), 𝔽(:, :, :, :))^2;
        ishermitian::Bool=true,
        amplitude::Union{Function, Nothing}=nothing,
        ismodulatable::Bool=true
    )

Coulomb term.

Type alias for `Term{:Coulomb, id, V, B, C<:TermCoupling, A<:TermAmplitude}`.
"""
const Coulomb{id, V, B, C<:TermCoupling, A<:TermAmplitude} = Term{:Coulomb, id, V, B, C, A}
@inline function Coulomb(
    id::Symbol, value, bondkind, coupling=Coupling(𝔽(:, :, :, :), 𝔽(:, :, :, :))^2;
    ishermitian::Bool=true,
    amplitude::Union{Function, Nothing}=nothing,
    ismodulatable::Bool=true
)
    return Term{:Coulomb}(id, value, bondkind, coupling, ishermitian; amplitude=amplitude, ismodulatable=ismodulatable)
end

"""
    FockTerm

Type alias for `Union{Onsite, Hopping, Pairing, Hubbard, InterOrbitalInterSpin, InterOrbitalIntraSpin, SpinFlip, PairHopping, Coulomb}`.
"""
const FockTerm = Union{Onsite, Hopping, Pairing, Hubbard, InterOrbitalInterSpin, InterOrbitalIntraSpin, SpinFlip, PairHopping, Coulomb}

# # SU(2) spin systems
# const sidtagmap = Dict(1=>'x', 2=>'y', 3=>'z', 4=>'+', 5=>'-')
# const sidseqmap = Dict(v=>k for (k, v) in sidtagmap)
# const sidajointmap = Dict('x'=>'x', 'y'=>'y', 'z'=>'z', '+'=>'-', '-'=>'+')
# const sidrepmap = Dict('x'=>'ˣ', 'y'=>'ʸ', 'z'=>'ᶻ', '+'=>'⁺', '-'=>'⁻', '0'=>'⁰')
# const sidreprevmap = Dict(v=>k for (k, v) in sidrepmap)

# ## SID
# """
#     SID{S, T<:Union{Char, Symbol, Colon}} <: SimpleInternalIndex

# The spin id.
# """
# struct SID{S, T<:Union{Char, Symbol, Colon}} <: SimpleInternalIndex
#     tag::T
#     function SID{S}(tag::Union{Char, Symbol, Colon}) where S
#         @assert isa(S, Rational{Int}) && S.den==2 || isa(S, Integer) || S==: "SID error: not supported spin($S)."
#         isa(tag, Char) && @assert tag in ('x', 'y', 'z', '+', '-') "SID error: not supported tag($tag)."
#         new{S, typeof(tag)}(tag)
#     end
# end
# @inline SID(tag::Union{Char, Symbol, Colon}) = SID{:}(tag)
# @inline Base.:(==)(sid₁::SID, sid₂::SID) = totalspin(sid₁)==totalspin(sid₂) && ==(efficientoperations, sid₁, sid₂)
# @inline Base.isequal(sid₁::SID, sid₂::SID) = isequal(totalspin(sid₁), totalspin(sid₂)) && isequal(efficientoperations, sid₁, sid₂)
# @inline Base.adjoint(sid::SID) = SID{totalspin(sid)}(sidajointmap[sid.tag])
# @inline Base.hash(sid::SID, h::UInt) = hash((totalspin(sid), sid.tag), h)
# @inline Base.show(io::IO, sid::SID) = @printf io "SID{%s}(%s)" totalspin(sid) default(sid.tag)
# @inline Base.show(io::IO, sid::SID{:}) = @printf io "SID(%s)" default(sid.tag)
# @inline @generated function Base.replace(sid::SID; kwargs...)
#     exprs = [:(get(kwargs, $name, getfield(sid, $name))) for name in QuoteNode.(fieldnames(sid))]
#     return :(rawtype(typeof(sid)){totalspin(sid)}($(exprs...)))
# end
# @inline statistics(::Type{<:SID}) = :b
# @inline SID(iid::SID, ::CompositeInternal) = iid

# """
#     totalspin(::SID) -> Rational{Int}/Int/Symbol
#     totalspin(::Type{<:SID}) -> Rational{Int}/Int/Symbol

#     totalspin(::Index{<:Union{Int, Colon}, <:SID}) -> Rational{Int}/Int/Symbol
#     totalspin(::Type{<:Index{<:Union{Int, Colon}, <:SID}}) -> Rational{Int}/Int/Symbol

#     totalspin(::CompositeIndex{<:Index{<:Union{Int, Colon}, <:SID}}) -> Rational{Int}/Int/Symbol
#     totalspin(::Type{<:CompositeIndex{<:Index{<:Union{Int, Colon}, <:SID}}}) -> Rational{Int}/Int/Symbol

# Get the total spin.
# """
# @inline totalspin(sid::SID) = totalspin(typeof(sid))
# @inline totalspin(index::Index{<:Union{Int, Colon}, <:SID}) = totalspin(typeof(index))
# @inline totalspin(index::CompositeIndex{<:Index{<:Union{Int, Colon}, <:SID}}) = totalspin(typeof(index))
# @inline totalspin(::Type{<:SID{S}}) where S = S
# @inline totalspin(::Type{<:Index{<:Union{Int, Colon}, <:SID{S}}}) where S = S
# @inline totalspin(::Type{<:CompositeIndex{<:Index{<:Union{Int, Colon}, <:SID{S}}}}) where S = S

# ### requested by Constraint
# @inline isdefinite(::Type{<:SID{T, Char} where T}) = true
# @inline indextype(::Type{SID}, ::Type{T}) where {T<:Union{Char, Symbol, Colon}} = SID{:, T}
# @inline indextype(::Type{SID{S}}, ::Type{T}) where {S, T<:Union{Char, Symbol, Colon}} = SID{S, T}

# ### matrix
# """
#     matrix(sid::SID{S, Char}, dtype::Type{<:Number}=Complex{Float}) where S -> Matrix{dtype}
#     matrix(index::Index{<:Union{Int, Colon}, <:SID}, dtype::Type{<:Number}=Complex{Float}) -> Matrix{dtype}
#     matrix(index::CompositeIndex{<:Index{<:Union{Int, Colon}, <:SID}}, dtype::Type{<:Number}=Complex{Float}) -> Matrix{dtype}

# Get the matrix representation of a sid.
# """
# function matrix(sid::SID{S, Char}, dtype::Type{<:Number}=Complex{Float}) where S
#     N = Int(2*S+1)
#     result = zeros(dtype, (N, N))
#     spin = convert(dtype, S)
#     for i = 1:N, j = 1:N
#         row, col = N+1-i, N+1-j
#         m, n = spin+1-i, spin+1-j
#         result[row, col] = (sid.tag == 'x') ? (delta(i+1, j)+delta(i, j+1))*sqrt(spin*(spin+1)-m*n)/2 :
#             (sid.tag == 'y') ? (delta(i+1, j)-delta(i, j+1))*sqrt(spin*(spin+1)-m*n)/2im :
#             (sid.tag == 'z') ? delta(i, j)*m :
#             (sid.tag == '+') ? delta(i+1, j)*sqrt(spin*(spin+1)-m*n) :
#             delta(i, j+1)*sqrt(spin*(spin+1)-m*n)
#     end
#     return result
# end
# @inline matrix(index::Index{<:Union{Int, Colon}, <:SID}, dtype::Type{<:Number}=Complex{Float}) = matrix(index.iid, dtype)
# @inline matrix(index::CompositeIndex{<:Index{<:Union{Int, Colon}, <:SID}}, dtype::Type{<:Number}=Complex{Float}) = matrix(getcontent(index, :index), dtype)

# ## Spin
# """
#     Spin{S} <: SimpleInternal{SID{S, Char}}

# The spin internal degrees of freedom.
# """
# struct Spin{S} <: SimpleInternal{SID{S, Char}}
#     function Spin{S}() where S
#         @assert isa(S, Rational{Int}) && S.den==2 || isa(S, Integer) "Spin error: not supported spin($S)."
#         new{S}()
#     end
# end
# @inline Base.eltype(::Type{Spin}) = (SID{S, Char} where S)
# @inline shape(sp::Spin) = (1:length(sidtagmap),)
# @inline Base.CartesianIndex(sid::SID, ::Spin) = CartesianIndex(sidseqmap[sid.tag])
# @inline SID(index::CartesianIndex{1}, sp::Spin) = SID{totalspin(sp)}(sidtagmap[index[1]])
# @inline Base.summary(io::IO, spin::Spin) = @printf io "%s-element Spin{%s}" length(spin) totalspin(spin)
# @inline Base.show(io::IO, spin::Spin) = @printf io "%s{%s}()" spin|>typeof|>nameof totalspin(spin)
# @inline Base.match(::Type{<:SID{:}}, ::Type{<:Spin{S}}) where S = true
# @inline Base.match(::Type{<:SID{S}}, ::Type{<:Spin{S}}) where S = true
# @inline Base.match(::Type{<:SID{S₁}}, ::Type{<:Spin{S₂}}) where {S₁, S₂} = false

# """
#     totalspin(::Spin) -> Rational{Int}/Int/Symbol
#     totalspin(::Type{<:Spin}) -> Rational{Int}/Int/Symbol

# Get the total spin.
# """
# @inline totalspin(spin::Spin) = totalspin(typeof(spin))
# @inline totalspin(::Type{<:Spin{S}}) where S = S

# ## LaTeX format output
# """
#     script(::Val{:tag}, sid::SID; kwargs...) -> String

# Get the requested script of an sid.
# """
# @inline script(::Val{:tag}, sid::SID; kwargs...) = sid.tag==(:) ? ":" : string(sid.tag)

# """
#     latexofspins

# The default LaTeX format of a spin index.
# """
# const latexofspins = LaTeX{(:tag,), (:site,)}('S')
# @inline latexname(::Type{<:Index{<:Union{Int, Colon}, <:SID}}) = Symbol("Index{Union{Int, Colon}, SID}")
# @inline latexname(::Type{<:CompositeIndex{<:Index{<:Union{Int, Colon}, <:SID}}}) = Symbol("CompositeIndex{Index{Union{Int, Colon}, SID}}")
# latexformat(Index{<:Union{Int, Colon}, <:SID}, latexofspins)
# latexformat(CompositeIndex{<:Index{<:Union{Int, Colon}, <:SID}}, latexofspins)
# @inline latexname(::Type{<:SID}) = Symbol("SID")
# latexformat(SID, LaTeX{(:tag,), ()}('S'))

# ## Permutation
# """
#     permute(id₁::SID, id₂::SID) -> Tuple{Vararg{Operator}}

# Permute two spin indexes and get the result.
# """
# function permute(id₁::SID, id₂::SID)
#     if id₁ ≠ id₂
#         S = totalspin(id₁)
#         if id₁.tag == 'x'
#             id₂.tag=='y' && return (Operator(+1im, SID{S}('z')), Operator(1, id₂, id₁))
#             id₂.tag=='z' && return (Operator(-1im, SID{S}('y')), Operator(1, id₂, id₁))
#             id₂.tag=='+' && return (Operator(-1, SID{S}('z')), Operator(1, id₂, id₁))
#             id₂.tag=='-' && return (Operator(+1, SID{S}('z')), Operator(1, id₂, id₁))
#         elseif id₁.tag == 'y'
#             id₂.tag=='x' && return (Operator(-1im, SID{S}('z')), Operator(1, id₂, id₁))
#             id₂.tag=='z' && return (Operator(+1im, SID{S}('x')), Operator(1, id₂, id₁))
#             id₂.tag=='+' && return (Operator(-1im, SID{S}('z')), Operator(1, id₂, id₁))
#             id₂.tag=='-' && return (Operator(-1im, SID{S}('z')), Operator(1, id₂, id₁))
#         elseif id₁.tag == 'z'
#             id₂.tag=='x' && return (Operator(+1im, SID{S}('y')), Operator(1, id₂, id₁))
#             id₂.tag=='y' && return (Operator(-1im, SID{S}('x')), Operator(1, id₂, id₁))
#             id₂.tag=='+' && return (Operator(+1, SID{S}('+')), Operator(1, id₂, id₁))
#             id₂.tag=='-' && return (Operator(-1, SID{S}('-')), Operator(1, id₂, id₁))
#         elseif id₁.tag == '+'
#             id₂.tag=='x' && return (Operator(+1, SID{S}('z')), Operator(1, id₂, id₁))
#             id₂.tag=='y' && return (Operator(+1im, SID{S}('z')), Operator(1, id₂, id₁))
#             id₂.tag=='z' && return (Operator(-1, SID{S}('+')), Operator(1, id₂, id₁))
#             id₂.tag=='-' && return (Operator(+2, SID{S}('z')), Operator(1, id₂, id₁))
#         elseif id₁.tag == '-'
#             id₂.tag=='x' && return (Operator(-1, SID{S}('z')), Operator(1, id₂, id₁))
#             id₂.tag=='y' && return (Operator(1im, SID{S}('z')), Operator(1, id₂, id₁))
#             id₂.tag=='z' && return (Operator(+1, SID{S}('-')), Operator(1, id₂, id₁))
#             id₂.tag=='+' && return (Operator(-2, SID{S}('z')), Operator(1, id₂, id₁))
#         end
#         error("permute error: not supported spin indexes.")
#     else
#         return (Operator(1, id₂, id₁),)
#     end
# end

# ## Coupling
# ### requested by ConstrainedInternal
# @inline shape(iidspace::ConstrainedInternal{<:SID, <:Spin}) = (sidrange(iidspace.iid.tag),)
# @inline sidrange(::Union{Symbol, Colon}) = 1:5
# @inline sidrange(tag::Char) = (pos=sidseqmap[tag]; pos:pos)

# ### MatrixCoupling
# """
#     MatrixCoupling(sites::Union{NTuple{2, Int}, Colon}, ::Type{S}, matrix::AbstractMatrix; rows::AbstractVector=SVector('x', 'y', 'z'), cols::AbstractVector=SVector('x', 'y', 'z')) where {S<:SID}

# Construct a set of `Coupling`s between two `Index{<:Union{Int, Colon}, <:SID}`s with the coefficients specified by a matrix.
# """
# function MatrixCoupling(sites::Union{NTuple{2, Int}, Colon}, ::Type{S}, matrix::AbstractMatrix; rows::AbstractVector=SVector('x', 'y', 'z'), cols::AbstractVector=SVector('x', 'y', 'z')) where {S<:SID}
#     @assert size(matrix)==(length(rows), length(cols)) "MatrixCoupling error: mismatched input matrix and rows/cols."
#     return MatrixCoupling(sites, S, Component(rows, cols, matrix))
# end

# ### Spin coupling matrix
# """
#     Heisenberg"" => SparseMatrixCSC([1 0 0; 0 1 0; 0 0 1])

# The Heisenberg coupling matrix.
# """
# macro Heisenberg_str(str::String)
#     str=="" && return SparseMatrixCSC([1 0 0; 0 1 0; 0 0 1])
#     error("@Heisenberg_str error: wrong input string.")
# end

# """
#     Ising"x" => SparseMatrixCSC([1 0 0; 0 0 0; 0 0 0])
#     Ising"y" => SparseMatrixCSC([0 0 0; 0 1 0; 0 0 0])
#     Ising"z" => SparseMatrixCSC([0 0 0; 0 0 0; 0 0 1])

# The Ising coupling matrix.
# """
# macro Ising_str(str::String)
#     str=="x" && return SparseMatrixCSC([1 0 0; 0 0 0; 0 0 0])
#     str=="y" && return SparseMatrixCSC([0 0 0; 0 1 0; 0 0 0])
#     str=="z" && return SparseMatrixCSC([0 0 0; 0 0 0; 0 0 1])
#     error("@Ising_str error: wrong input string.")
# end

# """
#     Γ"x" => SparseMatrixCSC([0 0 0; 0 0 1; 0 1 0])
#     Γ"y" => SparseMatrixCSC([0 0 1; 0 0 0; 1 0 0])
#     Γ"z" => SparseMatrixCSC([0 1 0; 1 0 0; 0 0 0])

# The Γ coupling matrix.
# """
# macro Γ_str(str::String)
#     str=="x" && return SparseMatrixCSC([0 0 0; 0 0 1; 0 1 0])
#     str=="y" && return SparseMatrixCSC([0 0 1; 0 0 0; 1 0 0])
#     str=="z" && return SparseMatrixCSC([0 1 0; 1 0 0; 0 0 0])
#     error("@Γ_str error: wrong input string.")
# end

# """
#     Γ′"x" => SparseMatrixCSC([0 1 1; 1 0 0; 1 0 0])
#     Γ′"y" => SparseMatrixCSC([0 1 0; 1 0 1; 0 1 0])
#     Γ′"z" => SparseMatrixCSC([0 0 1; 0 0 1; 1 1 0])

# The Γ′ coupling matrix.
# """
# macro Γ′_str(str::String)
#     str=="x" && return SparseMatrixCSC([0 1 1; 1 0 0; 1 0 0])
#     str=="y" && return SparseMatrixCSC([0 1 0; 1 0 1; 0 1 0])
#     str=="z" && return SparseMatrixCSC([0 0 1; 0 0 1; 1 1 0])
#     error("@Γ′_str error: wrong input string.")
# end

# """
#     DM"x" => SparseMatrixCSC([0 0 0; 0 0 1; 0 -1 0])
#     DM"y" => SparseMatrixCSC([0 0 -1; 0 0 0; 1 0 0])
#     DM"z" => SparseMatrixCSC([0 1 0; -1 0 0; 0 0 0])

# The DM coupling matrix.
# """
# macro DM_str(str::String)
#     str=="x" && return SparseMatrixCSC([0 0 0; 0 0 1; 0 -1 0])
#     str=="y" && return SparseMatrixCSC([0 0 -1; 0 0 0; 1 0 0])
#     str=="z" && return SparseMatrixCSC([0 1 0; -1 0 0; 0 0 0])
#     error("@DM_str error: wrong input string.")
# end

# ## Term
# """
#     SpinTerm(id::Symbol, value, bondkind, coupling; amplitude::Union{Function, Nothing}=nothing, ismodulatable::Bool=true)

# Generic spin term.

# Type alias for `Term{:SpinTerm, id, V, B, C<:TermCoupling, A<:TermAmplitude}`.
# """
# const SpinTerm{id, V, B, C<:TermCoupling, A<:TermAmplitude} = Term{:SpinTerm, id, V, B, C, A}
# @inline function SpinTerm(id::Symbol, value, bondkind, coupling; amplitude::Union{Function, Nothing}=nothing, ismodulatable::Bool=true)
#     return Term{:SpinTerm}(id, value, bondkind, coupling, true; amplitude=amplitude, ismodulatable=ismodulatable)
# end
# @inline function patternrule(::Val{:SpinTerm}, ::Val{termrank}, bondlength::Integer) where termrank
#     bondlength==1 && return ntuple(i->1, Val(termrank))
#     bondlength==2 && return ntuple(i->2-i%2, Val(termrank))
#     error("patternrule error: not supported for a generic bond containing $bondlength points.")
# end

# """
#     Zeeman(id::Symbol, value, direction::Char, g::Number=1; amplitude::Union{Function, Nothing}=nothing, ismodulatable::Bool=true)
#     Zeeman(id::Symbol, value, direction::Union{AbstractVector{<:Number}, Tuple{Number, Number}}, g::Union{Number, AbstractMatrix{<:Number}}=1; unit::Symbol=:degree, amplitude::Union{Function, Nothing}=nothing, ismodulatable::Bool=true)

# Zeeman term.

# Type alias for `Term{:Zeeman, id, V, Int, C<:TermCoupling, A<:TermAmplitude}`.
# """
# const Zeeman{id, V, C<:TermCoupling, A<:TermAmplitude} = Term{:Zeeman, id, V, Int, C, A}
# @inline function Zeeman(id::Symbol, value, direction::Char, g::Number=1; amplitude::Union{Function, Nothing}=nothing, ismodulatable::Bool=true)
#     @assert lowercase(direction)∈('x', 'y', 'z') "Zeeman error: not supported direction."
#     coupling = Coupling(g, :, SID, (lowercase(direction),))
#     return Term{:Zeeman}(id, value, 0, coupling, true; amplitude=amplitude, ismodulatable=ismodulatable)
# end
# @inline function Zeeman(id::Symbol, value, dir::Union{AbstractVector{<:Number}, Tuple{Number, Number}}, g::Union{Number, AbstractMatrix{<:Number}}=1; unit::Symbol=:degree, amplitude::Union{Function, Nothing}=nothing, ismodulatable::Bool=true)
#     couplings = dot(direction(dir, unit), Lande(g), SVector(Coupling(:, SID, ('x',)), Coupling(:, SID, ('y',)), Coupling(:, SID, ('z',))))
#     return Term{:Zeeman}(id, value, 0, couplings, true; amplitude=amplitude, ismodulatable=ismodulatable)
# end
# @inline Lande(g::Number) = SMatrix{3, 3}(g, 0, 0, 0, g, 0, 0, 0, g)
# @inline Lande(g::AbstractMatrix{<:Number}) = (@assert(size(g)==(3, 3), "Lande error: the g-tensor must be 3×3."); g)

# """
#     SingleIonAnisotropy(id::Symbol, value, direction::Char; amplitude::Union{Function, Nothing}=nothing, ismodulatable::Bool=true)
#     SingleIonAnisotropy(id::Symbol, value, matrix::AbstractMatrix{<:Number}; amplitude::Union{Function, Nothing}=nothing, ismodulatable::Bool=true)

# Single ion anisotropy term.

# Type alias for `Term{:SingleIonAnisotropy, id, V, Int, C<:TermCoupling, A<:TermAmplitude}`.
# """
# const SingleIonAnisotropy{id, V, C<:TermCoupling, A<:TermAmplitude} = Term{:SingleIonAnisotropy, id, V, Int, C, A}
# @inline function SingleIonAnisotropy(id::Symbol, value, direction::Char; amplitude::Union{Function, Nothing}=nothing, ismodulatable::Bool=true)
#     @assert lowercase(direction)∈('x', 'y', 'z') "SingleIonAnisotropy error: not supported direction."
#     coupling = Coupling(:, SID, (lowercase(direction), lowercase(direction)))
#     return Term{:SingleIonAnisotropy}(id, value, 0, coupling, true; amplitude=amplitude, ismodulatable=ismodulatable)
# end
# @inline function SingleIonAnisotropy(id::Symbol, value, matrix::AbstractMatrix{<:Number}; amplitude::Union{Function, Nothing}=nothing, ismodulatable::Bool=true)
#     @assert ishermitian(matrix) "SingleIonAnisotropy error: the anisotropy matrix must be Hermitian."
#     @assert size(matrix)==(3, 3) "SingleIonAnisotropy error: the anisotropy matrix must be 3×3."
#     couplings = dot(SVector(Coupling(:, SID, ('x',)), Coupling(:, SID, ('y',)), Coupling(:, SID, ('z',))), matrix, SVector(Coupling(:, SID, ('x',)), Coupling(:, SID, ('y',)), Coupling(:, SID, ('z',))))
#     return Term{:SingleIonAnisotropy}(id, value, 0, couplings, true; amplitude=amplitude, ismodulatable=ismodulatable)
# end

# """
#     Ising(id::Symbol, value, bondkind, direction::Char; amplitude::Union{Function, Nothing}=nothing, ismodulatable::Bool=true)

# Ising term.

# Type alias for `Term{:Ising, id, V, B, C<:TermCoupling, A<:TermAmplitude}`.
# """
# const Ising{id, V, B, C<:TermCoupling, A<:TermAmplitude} = Term{:Ising, id, V, B, C, A}
# @inline function Ising(id::Symbol, value, bondkind, direction::Char; amplitude::Union{Function, Nothing}=nothing, ismodulatable::Bool=true)
#     @assert lowercase(direction)∈('x', 'y', 'z') "Ising error: not supported direction."
#     coupling = Coupling(:, SID, (lowercase(direction), lowercase(direction)))
#     return Term{:Ising}(id, value, bondkind, coupling, true; amplitude=amplitude, ismodulatable=ismodulatable)
# end

# """
#     Heisenberg(id::Symbol, value, bondkind; form::Symbol=Symbol("+-z"), amplitude::Union{Function, Nothing}=nothing, ismodulatable::Bool=true)

# Heisenberg term.

# Type alias for `Term{:Heisenberg, id, V, B, C<:TermCoupling, A<:TermAmplitude}`.
# """
# const Heisenberg{id, V, B, C<:TermCoupling, A<:TermAmplitude} = Term{:Heisenberg, id, V, B, C, A}
# @inline function Heisenberg(id::Symbol, value, bondkind; form::Symbol=Symbol("+-z"), amplitude::Union{Function, Nothing}=nothing, ismodulatable::Bool=true)
#     @assert form∈(:xyz, Symbol("+-z")) "Heisenberg error: form should :xyz or Symbol(\"+-z\")."
#     couplings = if form==:xyz
#         Coupling(1//1, :, SID, ('x', 'x')) + Coupling(1//1, :, SID, ('y', 'y')) + Coupling(1//1, :, SID, ('z', 'z'))
#     else
#         Coupling(1//2, :, SID, ('+', '-')) + Coupling(1//2, :, SID, ('-', '+')) + Coupling(1//1, :, SID, ('z', 'z'))
#     end
#     return Term{:Heisenberg}(id, value, bondkind, couplings, true; amplitude=amplitude, ismodulatable=ismodulatable)
# end

# """
#     Kitaev(
#         id::Symbol, value, bondkind;
#         x::AbstractVector{<:Union{Number, Tuple{Number, Number}, AbstractVector{<:Number}}},
#         y::AbstractVector{<:Union{Number, Tuple{Number, Number}, AbstractVector{<:Number}}},
#         z::AbstractVector{<:Union{Number, Tuple{Number, Number}, AbstractVector{<:Number}}},
#         unit::Symbol=:degree,
#         amplitude::Union{Function, Nothing}=nothing,
#         ismodulatable::Bool=true
#     )

# Kitaev term. Since Kitaev term is symmetric on every bond, only one direction of a bond is needed. The inverse direction of a bond can be handled automatically by this function.

# Here, `x`, `y` and `z` assign the x-bonds, y-bonds, and z-bonds, respectively, with each kind of bond can be
# 1) a `Number` specifying the azimuth angle of a bond in the 2-dimensional case, or
# 2) a `Tuple{Number, Number}` specifying the polar and azimuth angle pairs of a bond in the 3-dimensional case, or
# 3) an `AbstractVector{<:Number}` specifying the direction of a bond.

# Type alias for `Term{:Kitaev, id, V, B, C<:TermCoupling, A<:TermAmplitude}`.
# """
# const Kitaev{id, V, B, C<:TermCoupling, A<:TermAmplitude} = Term{:Kitaev, id, V, B, C, A}
# function Kitaev(
#     id::Symbol, value, bondkind;
#     x::AbstractVector{<:Union{Number, Tuple{Number, Number}, AbstractVector{<:Number}}},
#     y::AbstractVector{<:Union{Number, Tuple{Number, Number}, AbstractVector{<:Number}}},
#     z::AbstractVector{<:Union{Number, Tuple{Number, Number}, AbstractVector{<:Number}}},
#     unit::Symbol=:degree,
#     amplitude::Union{Function, Nothing}=nothing,
#     ismodulatable::Bool=true
# )
#     dirs = (x=direction.(x, unit), y=direction.(y, unit), z=direction.(z, unit))
#     function kitaev(bond::Bond)
#         coordinate = rcoordinate(bond)
#         any(v->abs(isparallel(v, coordinate; atol=atol, rtol=rtol))==1, dirs.x) && return MatrixCoupling(: , SID, Ising"x")
#         any(v->abs(isparallel(v, coordinate; atol=atol, rtol=rtol))==1, dirs.y) && return MatrixCoupling(: , SID, Ising"y")
#         any(v->abs(isparallel(v, coordinate; atol=atol, rtol=rtol))==1, dirs.z) && return MatrixCoupling(: , SID, Ising"z")
#         error("Kitaev error: wrong bond.")
#     end
#     return Term{:Kitaev}(id, value, bondkind, kitaev, true; amplitude=amplitude, ismodulatable=ismodulatable)
# end

# """
#     Γ(
#         id::Symbol, value, bondkind;
#         x::AbstractVector{<:Union{Number, Tuple{Number, Number}, AbstractVector{<:Number}}},
#         y::AbstractVector{<:Union{Number, Tuple{Number, Number}, AbstractVector{<:Number}}},
#         z::AbstractVector{<:Union{Number, Tuple{Number, Number}, AbstractVector{<:Number}}},
#         unit::Symbol=:degree,
#         amplitude::Union{Function, Nothing}=nothing,
#         ismodulatable::Bool=true
#     )

# Γ Term. Since Γ term is symmetric on every bond, only one direction of a bond is needed. The inverse direction of a bond can be handled automatically by this function.

# Here, `x`, `y` and `z` assign the x-bonds, y-bonds, and z-bonds, respectively, with each kind of bond can be
# 1) a `Number` specifying the azimuth angle of a bond in the 2-dimensional case, or
# 2) a `Tuple{Number, Number}` specifying the polar and azimuth angle pairs of a bond in the 3-dimensional case, or
# 3) an `AbstractVector{<:Number}` specifying the direction of a bond.

# Type alias for `Term{:Γ, id, V, B, C<:TermCoupling, A<:TermAmplitude}`.
# """
# const Γ{id, V, B, C<:TermCoupling, A<:TermAmplitude} = Term{:Γ, id, V, B, C, A}
# function Γ(
#     id::Symbol, value, bondkind;
#     x::AbstractVector{<:Union{Number, Tuple{Number, Number}, AbstractVector{<:Number}}},
#     y::AbstractVector{<:Union{Number, Tuple{Number, Number}, AbstractVector{<:Number}}},
#     z::AbstractVector{<:Union{Number, Tuple{Number, Number}, AbstractVector{<:Number}}},
#     unit::Symbol=:degree,
#     amplitude::Union{Function, Nothing}=nothing,
#     ismodulatable::Bool=true
# )
#     dirs = (x=direction.(x, unit), y=direction.(y, unit), z=direction.(z, unit))
#     function γ(bond::Bond)
#         coordinate = rcoordinate(bond)
#         any(v->abs(isparallel(v, coordinate; atol=atol, rtol=rtol))==1, dirs.x) && return MatrixCoupling(: , SID, Γ"x")
#         any(v->abs(isparallel(v, coordinate; atol=atol, rtol=rtol))==1, dirs.y) && return MatrixCoupling(: , SID, Γ"y")
#         any(v->abs(isparallel(v, coordinate; atol=atol, rtol=rtol))==1, dirs.z) && return MatrixCoupling(: , SID, Γ"z")
#         error("Γ error: wrong bond.")
#     end
#     return Term{:Γ}(id, value, bondkind, γ, true; amplitude=amplitude, ismodulatable=ismodulatable)
# end

# """
#     Γ′(
#         id::Symbol, value, bondkind;
#         x::AbstractVector{<:Union{Number, Tuple{Number, Number}, AbstractVector{<:Number}}},
#         y::AbstractVector{<:Union{Number, Tuple{Number, Number}, AbstractVector{<:Number}}},
#         z::AbstractVector{<:Union{Number, Tuple{Number, Number}, AbstractVector{<:Number}}},
#         unit::Symbol=:degree,
#         amplitude::Union{Function, Nothing}=nothing,
#         ismodulatable::Bool=true
#     )

# Γ′ Term. Since Γ′ term is symmetric on every bond, only one direction of a bond is needed. The inverse direction of a bond can be handled automatically by this function.

# Here, `x`, `y` and `z` assign the x-bonds, y-bonds, and z-bonds, respectively, with each bond can be
# 1) a `Number` specifying the azimuth angle of a bond in the 2-dimensional case, or
# 2) a `Tuple{Number, Number}` specifying the polar and azimuth angle pairs of a bond in the 3-dimensional case, or
# 3) an `AbstractVector{<:Number}` specifying the direction of a bond.

# Type alias for `Term{:Γ′, id, V, B, C<:TermCoupling, A<:TermAmplitude}`.
# """
# const Γ′{id, V, B, C<:TermCoupling, A<:TermAmplitude} = Term{:Γ′, id, V, B, C, A}
# function Γ′(
#     id::Symbol, value, bondkind;
#     x::AbstractVector{<:Union{Number, Tuple{Number, Number}, AbstractVector{<:Number}}},
#     y::AbstractVector{<:Union{Number, Tuple{Number, Number}, AbstractVector{<:Number}}},
#     z::AbstractVector{<:Union{Number, Tuple{Number, Number}, AbstractVector{<:Number}}},
#     unit::Symbol=:degree,
#     amplitude::Union{Function, Nothing}=nothing,
#     ismodulatable::Bool=true
# )
#     dirs = (x=direction.(x, unit), y=direction.(y, unit), z=direction.(z, unit))
#     function γ′(bond::Bond)
#         coordinate = rcoordinate(bond)
#         any(v->abs(isparallel(v, coordinate; atol=atol, rtol=rtol))==1, dirs.x) && return MatrixCoupling(: , SID, Γ′"x")
#         any(v->abs(isparallel(v, coordinate; atol=atol, rtol=rtol))==1, dirs.y) && return MatrixCoupling(: , SID, Γ′"y")
#         any(v->abs(isparallel(v, coordinate; atol=atol, rtol=rtol))==1, dirs.z) && return MatrixCoupling(: , SID, Γ′"z")
#         error("Γ′ error: wrong bond.")
#     end
#     return Term{:Γ′}(id, value, bondkind, γ′, true; amplitude=amplitude, ismodulatable=ismodulatable)
# end

# """
#     DM(
#         id::Symbol,
#         value,
#         bondkind,
#         vectors::Pair{<:AbstractVector{<:Union{Number, Tuple{Number, Number}, AbstractVector{<:Number}}}, <:Union{Char, AbstractVector{<:Number}}}...;
#         unit::Symbol=:degree,
#         amplitude::Union{Function, Nothing}=nothing,
#         ismodulatable::Bool=true
#     )

# DM term. Since DM term is antisymmetric on every bond, only the positive direction of a bond is needed. The negative direction of a bond can be handled automatically by this function.

# Here, `vectors` specify the unit DM vector on every bond in the form `[bond₁, bond₂, ...]=>v`, where `bondᵢ` can be
# 1) a `Number` specifying the azimuth angle of a bond in the 2-dimensional case, or
# 2) a `Tuple{Number, Number}` specifying the polar and azimuth angle pairs of a bond in the 3-dimensional case, or
# 3) an `AbstractVector{<:Number}` specifying the direction of a bond;
# and `v` can be
# 1) a `Char` of 'x', 'y' or 'z', indicating the unit DM vector on the set of bonds is along the x, y or z direction, or
# 2) an `AbstractVector{<:Number}`, specifying the direction of the DM vector on the set of bonds.

# Type alias for `Term{:DM, id, V, B, C<:TermCoupling, A<:TermAmplitude}`.
# """
# const DM{id, V, B, C<:TermCoupling, A<:TermAmplitude} = Term{:DM, id, V, B, C, A}
# function DM(
#     id::Symbol,
#     value,
#     bondkind,
#     vectors::Pair{<:AbstractVector{<:Union{Number, Tuple{Number, Number}, AbstractVector{<:Number}}}, <:Union{Char, AbstractVector{<:Number}}}...;
#     unit::Symbol=:degree,
#     amplitude::Union{Function, Nothing}=nothing,
#     ismodulatable::Bool=true
# )
#     dirs = [direction.(pair.first, unit)=>direction(pair.second) for pair in vectors]
#     function dm(bond::Bond)
#         coordinate = rcoordinate(bond)
#         for pair in dirs
#             for v in pair.first
#                 parallel = isparallel(v, coordinate; atol=atol, rtol=rtol)
#                 abs(parallel)==1 && return MatrixCoupling(:, SID, parallel*(pair.second[1]*DM"x"+pair.second[2]*DM"y"+pair.second[3]*DM"z"))
#             end
#         end
#         error("dm error: wrong bond.")
#     end
#     return Term{:DM}(id, value, bondkind, dm, true; amplitude=amplitude, ismodulatable=ismodulatable)
# end

# # Phononic systems
# ## PID
# """
#     PID{D<:Union{Char, Symbol, Colon}} <: SimpleInternalIndex

# The phonon id.
# """
# struct PID{D<:Union{Char, Symbol, Colon}} <: SimpleInternalIndex
#     tag::Char
#     direction::D
#     function PID(tag::Char, direction::Union{Char, Symbol, Colon})
#         @assert tag∈('p', 'u') "PID error: wrong tag($tag)."
#         isa(direction, Char) && @assert direction∈('x', 'y', 'z') "PID error: wrong direction($direction)."
#         new{typeof(direction)}(tag, direction)
#     end
# end
# @inline Base.adjoint(pid::PID) = pid
# @inline statistics(::Type{<:PID}) = :b
# @inline Base.show(io::IO, pid::PID) = @printf io "PID(%s, %s)" repr(pid.tag) default(pid.direction)
# @inline PID(iid::PID, ::CompositeInternal) = iid

# ### requested by Constraint
# @inline isdefinite(::Type{PID{Char}}) = true
# @inline indextype(::Type{PID}, ::Type{Char}, ::Type{D}) where {D<:Union{Char, Symbol, Colon}} = PID{D}

# ## Phonon
# """
#     Phonon <: SimpleInternal{PID{Char}}

# The phonon internal degrees of freedom.
# """
# struct Phonon <: SimpleInternal{PID{Char}}
#     ndirection::Int
#     function Phonon(ndirection::Integer)
#         @assert ndirection∈(1, 2, 3) "Phonon error: wrong number of directions."
#         new(ndirection)
#     end
# end
# @inline shape(pn::Phonon) = (1:2, 1:pn.ndirection)
# @inline Base.CartesianIndex(pid::PID{Char}, ::Phonon) = CartesianIndex(pid.tag=='u' ? 1 : 2, Int(pid.direction)-Int('x')+1)
# @inline PID(index::CartesianIndex{2}, ::Phonon) = PID(index[1]==1 ? 'u' : 'p', Char(Int('x')+index[2]-1))

# ## LaTeX format output
# """
#     script(::Val{:direction}, pid::PID; kwargs...) -> String

# Get the requested script of an pid.
# """
# @inline script(::Val{:direction}, pid::PID; kwargs...) = pid.direction==(:) ? ":" : string(pid.direction)

# @inline body(pid::PID) = string(pid.tag)
# @inline body(index::Index{<:Union{Int, Colon}, <:PID}) = string(index.iid.tag)
# @inline body(index::CompositeIndex{<:Index{<:Union{Int, Colon}, <:PID}}) = string(getcontent(index, :index).iid.tag)
# """
#     latexofphonons

# The default LaTeX format of a phonon index.
# """
# const latexofphonons = LaTeX{(:direction,), (:site,)}(body, "", "")
# @inline latexname(::Type{<:Index{<:Union{Int, Colon}, <:PID}}) = Symbol("Index{Union{Int, Colon}, PID}")
# @inline latexname(::Type{<:CompositeIndex{<:Index{<:Union{Int, Colon}, <:PID}}}) = Symbol("CompositeIndex{Index{Union{Int, Colon}, PID}}")
# latexformat(Index{<:Union{Int, Colon}, <:PID}, latexofphonons)
# latexformat(CompositeIndex{<:Index{<:Union{Int, Colon}, <:PID}}, latexofphonons)
# @inline latexname(::Type{<:PID}) = Symbol("PID")
# latexformat(PID, LaTeX{(:direction,), ()}(body, "", ""))

# ## Permutation
# """
#     permute(id₁::PID, id₂::PID) -> Tuple{Vararg{Operator}}

# Permute two phonon indexes and get the result.
# """
# function permute(id₁::PID, id₂::PID)
#     if id₁.direction==id₂.direction && id₁.tag≠id₂.tag
#         return (Operator(id₁.tag=='u' ? 1im : -1im), Operator(1, id₂, id₁))
#     else
#         return (Operator(1, id₂, id₁),)
#     end
# end

# ## Coupling
# ### requested by ConstrainedInternal
# @inline function shape(iidspace::ConstrainedInternal{<:PID, Phonon})
#     return (iidspace.iid.tag=='u' ? (1:1) : (2:2), pidrange(iidspace.iid.direction, iidspace.internal.ndirection))
# end
# @inline pidrange(::Union{Symbol, Colon}, ndirection::Int) = 1:ndirection
# @inline function pidrange(direction::Char, ndirection::Int)
#     direction = Int(direction) - Int('x') + 1
#     @assert 0<direction<ndirection+1 "pidrange error: direction out of range."
#     return direction:direction
# end

# ### MatrixCoupling
# """
#     MatrixCoupling(sites::Union{NTuple{2, Int}, Colon}, ::Type{PID}, matrix::AbstractMatrix; rows::Union{AbstractVector, Nothing}=nothing, cols::Union{AbstractVector, Nothing}=nothing)

# Construct a set of `Coupling`s corresponding to the dynamical matrix of phonons.
# """
# function MatrixCoupling(sites::Union{NTuple{2, Int}, Colon}, ::Type{PID}, matrix::AbstractMatrix; rows::Union{AbstractVector, Nothing}=nothing, cols::Union{AbstractVector, Nothing}=nothing)
#     @assert size(matrix)[1]∈(1, 2, 3) && size(matrix)[2]∈(1, 2, 3) "MatrixCoupling error: mismatched dimension of input matrix."
#     isnothing(rows) && (rows = size(matrix)[1]==1 ? SVector('x') : size(matrix)[1]==2 ? SVector('x', 'y') : SVector('x', 'y', 'z'))
#     isnothing(cols) && (cols = size(matrix)[2]==1 ? SVector('x') : size(matrix)[2]==2 ? SVector('x', 'y') : SVector('x', 'y', 'z'))
#     @assert size(matrix)==(length(rows), length(cols)) "MatrixCoupling error: mismatched input matrix and rows/cols."
#     return MatrixCoupling(sites, PID, Component(SVector('u'), SVector('u'), default_matrix), Component(rows, cols, matrix))
# end

# ### expand
# """
#     expand(::Val{:Hooke}, pnc::Coupling{<:Number, <:NTuple{2, Index{<:Union{Int, Colon}, PID{Colon}}}}, bond::Bond, hilbert::Hilbert) -> PPExpand

# Expand the default phonon potential coupling on a given bond.
# """
# function expand(::Val{:Hooke}, pnc::Coupling{<:Number, <:NTuple{2, Index{<:Union{Int, Colon}, PID{Colon}}}}, bond::Bond, hilbert::Hilbert)
#     R̂ = rcoordinate(bond)/norm(rcoordinate(bond))
#     @assert pnc.indexes.iids.tags==('u', 'u') "expand error: wrong tags of Hooke coupling."
#     @assert isapprox(pnc.value, 1, atol=atol, rtol=rtol) "expand error: wrong coefficient of Hooke coupling."
#     @assert pnc.indexes.sites∈((1, 2), (:, :)) "expand error: wrong sites of Hooke coupling."
#     pn₁ = filter(pnc.indexes[1].iid, hilbert[bond[1].site])
#     pn₂ = filter(pnc.indexes[2].iid, hilbert[bond[2].site])
#     @assert pn₁.ndirection==pn₂.ndirection==length(R̂) "expand error: mismatched number of directions."
#     return PPExpand(R̂, (bond[1], bond[2]))
# end
# struct PPExpand{N, D<:Number} <: VectorSpace{Operator{D, ID{CoordinatedIndex{Index{Int, PID{Char}}, SVector{N, D}}, 2}}}
#     direction::SVector{N, D}
#     points::NTuple{2, Point{N, D}}
# end
# @inline VectorSpaceStyle(::Type{<:PPExpand}) = VectorSpaceCartesian()
# @inline shape(pnce::PPExpand) = (1:length(pnce.direction), 1:length(pnce.direction), 1:4)
# function Operator(index::CartesianIndex{3}, pnce::PPExpand)
#     direction₁ = Char(Int('x')+index[1]-1)
#     direction₂ = Char(Int('x')+index[2]-1)
#     coeff = index[3]∈(1, 4) ? 1 : -1
#     pos₁, pos₂ = index[3]==1 ? (1, 1) : index[3]==2 ? (1, 2) : index[3]==3 ? (2, 1) : (2, 2)
#     index₁ = CoordinatedIndex(Index(pnce.points[pos₁].site, PID('u', direction₁)), pnce.points[pos₁].rcoordinate, pnce.points[pos₁].icoordinate)
#     index₂ = CoordinatedIndex(Index(pnce.points[pos₂].site, PID('u', direction₂)), pnce.points[pos₂].rcoordinate, pnce.points[pos₂].icoordinate)
#     return Operator(pnce.direction[index[1]]*pnce.direction[index[2]]*coeff, index₁, index₂)
# end

# ## Term
# """
#     Kinetic(id::Symbol, value; amplitude::Union{Function, Nothing}=nothing, ismodulatable::Bool=true)

# Kinetic energy of phonons.

# Type alias for `Term{:Kinetic, id, V, Int, C<:TermCoupling, A<:TermAmplitude}`.
# """
# const Kinetic{id, V, C<:TermCoupling, A<:TermAmplitude} = Term{:Kinetic, id, V, Int, C, A}
# @inline function Kinetic(id::Symbol, value; amplitude::Union{Function, Nothing}=nothing, ismodulatable::Bool=true)
#     return Term{:Kinetic}(id, value, 0, Coupling(:, PID, ('p', 'p'), :), true; amplitude=amplitude, ismodulatable=ismodulatable)
# end

# """
#     Hooke(id::Symbol, value, bondkind; amplitude::Union{Function, Nothing}=nothing, ismodulatable::Bool=true)

# Potential energy of phonons by the Hooke's law.

# Type alias for `Term{:Hooke, id, V, B, C<:TermCoupling, A<:TermAmplitude}`
# """
# const Hooke{id, V, B, C<:TermCoupling, A<:TermAmplitude} = Term{:Hooke, id, V, B, C, A}
# @inline function Hooke(id::Symbol, value, bondkind; amplitude::Union{Function, Nothing}=nothing, ismodulatable::Bool=true)
#     return Term{:Hooke}(id, value, bondkind, Coupling(:, PID, ('u', 'u'), :), true; amplitude=amplitude, ismodulatable=ismodulatable)
# end

# """
#     Elastic(id::Symbol, value, bondkind, coupling; amplitude::Union{Function, Nothing}=nothing, ismodulatable::Bool=true)

# Generic elastic energy of phonons.

# Type alias for `Term{:Elastic, id, V, B, C<:TermCoupling, A<:TermAmplitude}`
# """
# const Elastic{id, V, B, C<:TermCoupling, A<:TermAmplitude} = Term{:Elastic, id, V, B, C, A}
# @inline function Elastic(id::Symbol, value, bondkind, coupling; amplitude::Union{Function, Nothing}=nothing, ismodulatable::Bool=true)
#     value, factor = promote(value, 1//2)
#     return Term{:Elastic, id}(value, bondkind, TermCoupling(coupling), TermAmplitude(amplitude), true, ismodulatable, factor)
# end
# function expand!(operators::Operators, term::Elastic, bond::Bond, hilbert::Hilbert; half::Bool=false)
#     argtypes = Tuple{Operators, Term, Bond, Hilbert}
#     invoke(expand!, argtypes, operators, term, bond, hilbert; half=half)
#     invoke(expand!, argtypes, operators, term, reverse(bond), hilbert; half=half)
#     return operators
# end
# """
#     PhononTerm

# Type alias for `Union{Kinetic, Hooke, Elastic}`.
# """
# const PhononTerm = Union{Kinetic, Hooke, Elastic}

end # module
