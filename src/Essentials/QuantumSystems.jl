module QuantumSystems

using LinearAlgebra: dot, norm
using Printf: @printf, @sprintf
using SparseArrays: SparseMatrixCSC
using StaticArrays: SVector
using ..DegreesOfFreedom: wildcard, AbstractCompositeIndex, Component, CompositeIID, CompositeIndex, CompositeInternal, Diagonal, Hilbert, IIDSpace, Index, SimpleIID, SimpleInternal, Term, TermAmplitude, TermCoupling, TermModulate, @iids
using ..QuantumOperators: ID, LaTeX, Operator, OperatorProd, Operators, latexformat
using ..Spatials: Bond, Point, rcoordinate
using ...Essentials: dtype, kind
using ...Interfaces: decompose, dimension
using ...Prerequisites: atol, rtol, Float, decimaltostr, delta
using ...Prerequisites.Traits: efficientoperations, getcontent, rawtype
using ...Prerequisites.VectorSpaces: VectorSpace, VectorSpaceCartesian, VectorSpaceStyle

import ..QuantumOperators: latexname, matrix, optype, script
import ..DegreesOfFreedom: Constraint, Coupling, MatrixCoupling, constrainttype, couplingcenters, iidtype, isconcreteiid, statistics
import ...Interfaces: ⊗, ⋅, expand, expand!, permute, rank
import ...Prerequisites.VectorSpaces: shape

# Canonical complex fermionic/bosonic systems
export σ⁰, σˣ, σʸ, σᶻ, σ⁺, σ⁻, annihilation, creation, latexofbosons, latexoffermions, latexofparticles
export Coulomb, FID, Fock, FockTerm, Hopping, Hubbard, InterOrbitalInterSpin, InterOrbitalIntraSpin, Onsite, PairHopping, Pairing, SpinFlip, isnormalordered

# SU(2) spin systems
export latexofspins, SID, Spin, SpinTerm, totalspin, @dm_str, @gamma_str, @heisenberg_str, @ising_str

# Phononic systems
export latexofphonons, Elastic, PID, Phonon, Kinetic, Hooke, PhononTerm

# Canonical complex fermionic/bosonic systems and hardcore bosonic systems
## FID
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

@inline default(::typeof(:)) = ":"
@inline default(value) = value
"""
    FID{T, O<:Union{Int, Symbol, typeof(:)}, S<:Union{Int, Symbol, typeof(:)}, N<:Union{Int, Symbol, typeof(:)}} <: SimpleIID

The Fock id.
"""
struct FID{T, O<:Union{Int, Symbol, typeof(:)}, S<:Union{Int, Symbol, typeof(:)}, N<:Union{Int, Symbol, typeof(:)}} <: SimpleIID
    orbital::O
    spin::S
    nambu::N
    function FID{T}(orbital::Union{Int, Symbol, typeof(:)}, spin::Union{Int, Symbol, typeof(:)}, nambu::Union{Int, Symbol, typeof(:)}) where T
        @assert T∈(:f, :b, wildcard) "FID error: wrong statistics."
        isa(nambu, Int) && @assert nambu∈(1, 2) "FID error: wrong input nambu($nambu)."
        new{T, typeof(orbital), typeof(spin), typeof(nambu)}(orbital, spin, nambu)
    end
end
@inline FID(orbital::Union{Int, Symbol, typeof(:)}, spin::Union{Int, Symbol, typeof(:)}, nambu::Union{Int, Symbol, typeof(:)}) = FID{wildcard}(orbital, spin, nambu)
@inline Base.:(==)(fid₁::FID, fid₂::FID) = statistics(fid₁)==statistics(fid₂) && ==(efficientoperations, fid₁, fid₂)
@inline Base.isequal(fid₁::FID, fid₂::FID) = isequal(statistics(fid₁), statistics(fid₂)) && isequal(efficientoperations, fid₁, fid₂)
@inline Base.hash(fid::FID, h::UInt) = hash((statistics(fid), fid.orbital, fid.spin, fid.nambu), h)
@inline Base.show(io::IO, fid::FID) = @printf io "FID{%s}(%s)" repr(statistics(fid)) join((fid.orbital|>default, fid.spin|>default, fid.nambu|>default), ", ")
@inline Base.show(io::IO, fid::FID{wildcard}) = @printf io "FID(%s)" join((fid.orbital|>default, fid.spin|>default, fid.nambu|>default), ", ")
@inline Base.adjoint(fid::FID{T, <:Union{Int, Symbol, typeof(:)}, <:Union{Int, Symbol, typeof(:)}, Int}) where T = FID{T}(fid.orbital, fid.spin, 3-fid.nambu)
@inline @generated function Base.replace(fid::FID; kwargs...)
    exprs = [:(get(kwargs, $name, getfield(fid, $name))) for name in QuoteNode.(fieldnames(fid))]
    return :(rawtype(typeof(fid)){statistics(fid)}($(exprs...)))
end
@inline statistics(::Type{<:FID{T}}) where T = T
@inline FID(iid::FID, ::CompositeInternal) = iid

### required by Constraint
@inline isconcreteiid(::Type{<:FID{T, Int, Int, Int} where T}) = true
@inline iidtype(::Type{FID}, ::Type{O}, ::Type{S}, ::Type{N}) where {O<:Union{Int, Symbol, typeof(:)}, S<:Union{Int, Symbol, typeof(:)}, N<:Union{Int, Symbol, typeof(:)}} = FID{wildcard, O, S, N}
@inline iidtype(::Type{FID{T}}, ::Type{O}, ::Type{S}, ::Type{N}) where {T, O<:Union{Int, Symbol, typeof(:)}, S<:Union{Int, Symbol, typeof(:)}, N<:Union{Int, Symbol, typeof(:)}} = FID{T, O, S, N}

## Fock
"""
    Fock{T} <: SimpleInternal{FID{T, Int, Int, Int}}

The Fock internal degrees of freedom.
"""
struct Fock{T} <: SimpleInternal{FID{T, Int, Int, Int}}
    norbital::Int
    nspin::Int
    function Fock{T}(norbital::Int, nspin::Int) where T
        @assert T∈(:f, :b) "Fock error: wrong statistics."
        new{T}(norbital, nspin)
    end
end
@inline Base.eltype(::Type{Fock}) = (FID{T, Int, Int, Int} where T)
@inline shape(fock::Fock) = (1:fock.norbital, 1:fock.nspin, 1:2)
@inline Base.CartesianIndex(fid::FID{T}, fock::Fock{T}) where T = CartesianIndex(fid.orbital, fid.spin, fid.nambu)
@inline FID(index::CartesianIndex{3}, fock::Fock) = FID{statistics(fock)}(index[1], index[2], index[3])
@inline Base.summary(io::IO, fock::Fock) = @printf io "%s-element Fock{%s}" length(fock) repr(statistics(fock))
@inline Base.show(io::IO, fock::Fock) = @printf io "%s{%s}(%s)" fock|>typeof|>nameof repr(statistics(fock)) join(("$name=$(getfield(fock, name))" for name in fock|>typeof|>fieldnames), ", ")
@inline Base.match(::Type{<:FID{wildcard}}, ::Type{<:Fock{T}}) where T = true
@inline Base.match(::Type{<:FID{T}}, ::Type{<:Fock{T}}) where T = true
@inline Base.match(::Type{<:FID{T₁}}, ::Type{<:Fock{T₂}}) where {T₁, T₂} = false

## LaTeX format output
"""
    script(::Val{:orbital}, fid::FID; kwargs...) -> Int
    script(::Val{:spint}, fid::FID; kwargs...) -> Int
    script(::Val{:spsym}, fid::FID; kwargs...) -> String
    script(::Val{:nambu}, fid::FID; kwargs...) -> String

Get the required script of an fid.
"""
@inline script(::Val{:orbital}, fid::FID; kwargs...) = default(fid.orbital)
@inline script(::Val{:spint}, fid::FID; kwargs...) = default(fid.spin)
@inline script(::Val{:spsym}, fid::FID; kwargs...) = fid.spin==1 ? "↓" : fid.spin==2 ? "↑" : string(default(fid.spin))
@inline script(::Val{:nambu}, fid::FID; kwargs...) = fid.nambu==creation ? "\\dagger" : fid.nambu==annihilation ? "" : string(default(fid.nambu))

"""
    script(::Val{:site}, index::Index{<:FID}; kwargs...) -> Int
    script(attr::Val, index::Index{<:FID}; kwargs...) -> Union{Int, String}

Get the required script of a Fock index.
"""
@inline script(::Val{:site}, index::Index{<:FID}; kwargs...) = index.site
@inline script(attr::Val, index::Index{<:FID}; kwargs...) = script(attr, index.iid; kwargs...)

"""
    latexoffermions

The default LaTeX format of a fermionic index.
"""
const latexoffermions = LaTeX{(:nambu,), (:site, :orbital, :spsym)}('c')
@inline latexname(::Type{<:Index{<:FID{:f}}}) = Symbol("Index{FID{:f}}")
@inline latexname(::Type{<:AbstractCompositeIndex{<:Index{<:FID{:f}}}}) = Symbol("AbstractCompositeIndex{Index{FID{:f}}}")
latexformat(Index{<:FID{:f}}, latexoffermions)
latexformat(AbstractCompositeIndex{<:Index{<:FID{:f}}}, latexoffermions)
@inline latexname(::Type{<:FID{:f}}) = Symbol("FID{:f}")
latexformat(FID{:f}, LaTeX{(:nambu,), (:orbital, :spsym)}('c'))

"""
    latexofbosons

The default LaTeX format of a bosonic index.
"""
const latexofbosons = LaTeX{(:nambu,), (:site, :orbital, :spsym)}('b')
@inline latexname(::Type{<:Index{<:FID{:b}}}) = Symbol("Index{FID{:b}}")
@inline latexname(::Type{<:AbstractCompositeIndex{<:Index{<:FID{:b}}}}) = Symbol("AbstractCompositeIndex{Index{FID{:b}}}")
latexformat(Index{<:FID{:b}}, latexofbosons)
latexformat(AbstractCompositeIndex{<:Index{<:FID{:b}}}, latexofbosons)
@inline latexname(::Type{<:FID{:b}}) = Symbol("FID{:b}")
latexformat(FID{:b}, LaTeX{(:nambu,), (:orbital, :spsym)}('b'))

"""
    latexofparticles

The default LaTeX format of a wildcard Fock index.
"""
const latexofparticles = LaTeX{(:nambu,), (:site, :orbital, :spsym)}('f')
@inline latexname(::Type{<:Index{<:FID{wildcard}}}) = Symbol("Index{FID}")
@inline latexname(::Type{<:AbstractCompositeIndex{<:Index{<:FID{wildcard}}}}) = Symbol("AbstractCompositeIndex{Index{FID}}")
latexformat(Index{<:FID{wildcard}}, latexofparticles)
latexformat(AbstractCompositeIndex{<:Index{<:FID{wildcard}}}, latexofparticles)
@inline latexname(::Type{<:FID{wildcard}}) = Symbol("FID")
latexformat(FID{wildcard}, LaTeX{(:nambu,), (:orbital, :spsym)}('f'))

## Boundary
"""
    angle(id::CompositeIndex{<:Index{<:FID}}, vectors::AbstractVector{<:AbstractVector{<:Number}}, values::AbstractVector{<:Number}) -> Complex{<:Number}

Get the twist phase corresponding to a Fock index.
"""
function Base.angle(id::CompositeIndex{<:Index{<:FID}}, vectors::AbstractVector{<:AbstractVector{<:Number}}, values::AbstractVector{<:Number})
    datatype = promote_type(eltype(values), Float)
    phase = length(vectors)==1 ? 2*convert(datatype, pi)*dot(decompose(id.icoordinate, vectors[1]), values) :
            length(vectors)==2 ? 2*convert(datatype, pi)*dot(decompose(id.icoordinate, vectors[1], vectors[2]), values) :
            length(vectors)==3 ? 2*convert(datatype, pi)*dot(decompose(id.icoordinate, vectors[1], vectors[2], vectors[3]), values) :
            error("angle error: not supported number of input basis vectors.")
    return id.index.iid.nambu==annihilation ? phase : id.index.iid.nambu==creation ? -phase : error("angle error: not supported Fock index.")
end

## New methods
"""
    isnormalordered(opt::Operator{<:Number, <:ID{AbstractCompositeIndex{<:Index{<:FID}}}}) -> Bool

Judge whether an operator is normal ordered.
"""
function isnormalordered(opt::Operator{<:Number, <:ID{AbstractCompositeIndex{<:Index{<:FID}}}})
    flag = true
    for i = 1:rank(opt)
        flag && (opt.id[i].index.iid.nambu == annihilation) && (flag = false)
        flag || (opt.id[i].index.iid.nambu == creation) && return false
    end
    return true
end

"""
    *(
        f1::Operator{<:Number, <:ID{AbstractCompositeIndex{<:Index{<:FID{:f}}}}},
        f2::Operator{<:Number, <:ID{AbstractCompositeIndex{<:Index{<:FID{:f}}}}}
    ) -> Union{Nothing, Operator}

Get the multiplication of two fermionic Fock operators.
"""
@inline function Base.:*(
        f1::Operator{<:Number, <:ID{AbstractCompositeIndex{<:Index{<:FID{:f}}}}},
        f2::Operator{<:Number, <:ID{AbstractCompositeIndex{<:Index{<:FID{:f}}}}}
        )
    rank(f1)>0 && rank(f2)>0 && f1.id[end]==f2.id[1] && return nothing
    return invoke(*, Tuple{OperatorProd, OperatorProd}, f1, f2)
end

## Permutation
"""
    permute(id₁::CompositeIndex{<:Index{<:FID{:f}}}, id₂::CompositeIndex{<:Index{<:FID{:f}}}) -> Tuple{Vararg{Operator}}

Permute two fermionic indexes and get the result.
"""
function permute(id₁::CompositeIndex{<:Index{<:FID{:f}}}, id₂::CompositeIndex{<:Index{<:FID{:f}}})
    @assert id₁.index≠id₂.index || id₁.rcoordinate≠id₂.rcoordinate || id₁.icoordinate≠id₂.icoordinate "permute error: permuted ids should not be equal to each other."
    if id₁.index'==id₂.index && id₁.rcoordinate==id₂.rcoordinate && id₁.icoordinate==id₂.icoordinate
        return (Operator(1), Operator(-1, id₂, id₁))
    else
        return (Operator(-1, id₂, id₁),)
    end
end

"""
    permute(id₁::CompositeIndex{<:Index{<:FID{:b}}}, id₂::CompositeIndex{<:Index{<:FID{:b}}}) -> Tuple{Vararg{Operator}}

Permute two bosonic indexes and get the result.
"""
function permute(id₁::CompositeIndex{<:Index{<:FID{:b}}}, id₂::CompositeIndex{<:Index{<:FID{:b}}})
    @assert id₁.index≠id₂.index || id₁.rcoordinate≠id₂.rcoordinate || id₁.icoordinate≠id₂.icoordinate "permute error: permuted ids should not be equal to each other."
    if id₁.index'==id₂.index && id₁.rcoordinate==id₂.rcoordinate && id₁.icoordinate==id₂.icoordinate
        if id₁.index.iid.nambu == creation
            return (Operator(1), Operator(1, id₂, id₁))
        else
            return (Operator(-1), Operator(1, id₂, id₁))
        end
    else
        return (Operator(1, id₂, id₁),)
    end
end

## Coupling
### requested by IIDSpace
@inline function shape(iidspace::IIDSpace{<:FID, <:Fock})
    @assert isa(iidspace.iid.nambu, Int) "shape error: nambu not determined."
    obrange = fidshape(iidspace.iid.orbital, iidspace.internal.norbital)
    sprange = fidshape(iidspace.iid.spin, iidspace.internal.nspin)
    phrange = iidspace.iid.nambu:iidspace.iid.nambu
    return (obrange, sprange, phrange)
end
@inline fidshape(::Union{Symbol, typeof(:)}, n::Int) = 1:n
@inline fidshape(v::Int, n::Int) where field = ((@assert 0<v<n+1 "shape error: $field out of range."); v:v)

### Coupling
"""
    Coupling{F}(orbital::Union{NTuple{N, Int}, typeof(:)}, spin::Union{NTuple{N, Int}, typeof(:)}, nambu::Union{NTuple{N, Int}, typeof(:)}) where {F<:FID, N}
    Coupling{F}(value, orbital::Union{NTuple{N, Int}, typeof(:)}, spin::Union{NTuple{N, Int}, typeof(:)}, nambu::Union{NTuple{N, Int}, typeof(:)}) where {F<:FID, N}

Construct a `Coupling` among different `FID`s with the given `orbital`, `spin` and `nambu`.
"""
@inline Coupling{F}(orbital::Union{NTuple{N, Int}, typeof(:)}, spin::Union{NTuple{N, Int}, typeof(:)}, nambu::Union{NTuple{N, Int}, typeof(:)}) where {F<:FID, N} = Coupling{F}(1, orbital, spin, nambu)
@inline function Coupling{F}(value, orbital::Union{NTuple{N, Int}, typeof(:)}, spin::Union{NTuple{N, Int}, typeof(:)}, nambu::Union{NTuple{N, Int}, typeof(:)}) where {F<:FID, N}
    return Coupling(value, ID(F, default(orbital, N|>Val), default(spin, N|>Val), default(nambu, N|>Val)))
end
@inline default(fields, ::Val) = fields
@inline default(::typeof(:), N::Val) = ntuple(i->:, N)

### requested by expand of Coupling
@inline function CompositeIID(coupling::Coupling{V, I}, info::Val=Val(:term)) where {V, I<:ID{FID}}
    return CompositeIID(map((fid::FID, order::Int)->FID{statistics(fid)}(fid.orbital, fid.spin, nambu(info, fid.nambu, order)), coupling.iids, ntuple(i->i, Val(rank(I)))))
end
@inline nambu(::Val, ::typeof(:), order::Int) = order%2==1 ? creation : annihilation
@inline nambu(::Val, nambu::Int, ::Int) = nambu

### Constraint
@inline Constraint(::NTuple{N, FID{S, typeof(:), Int} where S}) where N = Constraint{N}(Diagonal(:orbital))
@inline Constraint(::NTuple{N, FID{S, Int, typeof(:)} where S}) where N = Constraint{N}(Diagonal(:spin))
@inline Constraint(::NTuple{N, FID{S, typeof(:), typeof(:)} where S}) where N = Constraint{N}(Diagonal(:orbital, :spin))

### MatrixCoupling
const default_element = [:]
const default_matrix = SparseMatrixCSC(hcat(1))
"""
    MatrixCoupling{F}(orbital::Union{AbstractMatrix, typeof(:)}, spin::Union{AbstractMatrix, typeof(:)}, nambu::Union{AbstractMatrix, typeof(:)}) where {F<:FID}

Construct a set of `Coupling`s between two `FID`s with the coefficients specified by matrices acting on separated internal spaces.
"""
@inline function MatrixCoupling{F}(orbital::Union{AbstractMatrix, typeof(:)}, spin::Union{AbstractMatrix, typeof(:)}, nambu::Union{AbstractMatrix, typeof(:)}) where {F<:FID}
    return MatrixCoupling{F}(Component(F, Val(:orbital), orbital), Component(F, Val(:spin), spin), Component(F, Val(:nambu), nambu))
end
@inline Component(::Type{<:FID}, ::Val, matrix::AbstractMatrix) = Component(1:size(matrix)[1], 1:size(matrix)[2], matrix)
@inline Component(::Type{<:FID}, ::Val, ::typeof(:)) = Component(default_element, default_element, default_matrix)
@inline Component(::Type{<:FID}, ::Val{:nambu}, matrix::AbstractMatrix) = (@assert size(matrix)==(2, 2) "Component error: for nambu subspace, the input matrix must be 2×2."; Component(1:1:2, 2:-1:1, matrix))
@inline constrainttype(::Type{<:MatrixCoupling{<:FID, <:Tuple{Component{typeof(:)}, Component{Int}, Component}}}) = Constraint{(2,), 1, Tuple{Diagonal{(:orbital,)}}}
@inline constrainttype(::Type{<:MatrixCoupling{<:FID, <:Tuple{Component{Int}, Component{typeof(:)}, Component}}}) = Constraint{(2,), 1, Tuple{Diagonal{(:spin,)}}}
@inline constrainttype(::Type{<:MatrixCoupling{<:FID, <:Tuple{Component{typeof(:)}, Component{typeof(:)}, Component}}}) = Constraint{(2,), 1, Tuple{Diagonal{(:orbital, :spin)}}}

### Pauli matrices
"""
    const σ⁰ = SparseMatrixCSC([1 0; 0 1])
    const σˣ = SparseMatrixCSC([0 1; 1 0])
    const σʸ = SparseMatrixCSC([0 1im; -1im 0])
    const σᶻ = SparseMatrixCSC([-1 0; 0 1])
    const σ⁺ = SparseMatrixCSC([0 0; 1 0])
    const σ⁻ = SparseMatrixCSC([0 1; 0 0])

The Pauli matrix σ⁰, σˣ, σʸ, σᶻ, σ⁺, σ⁻.
"""
const σ⁰ = SparseMatrixCSC([1 0; 0 1])
const σˣ = SparseMatrixCSC([0 1; 1 0])
const σʸ = SparseMatrixCSC([0 1im; -1im 0])
const σᶻ = SparseMatrixCSC([-1 0; 0 1])
const σ⁺ = SparseMatrixCSC([0 0; 1 0])
const σ⁻ = SparseMatrixCSC([0 1; 0 0])

## Term
"""
    Onsite(id::Symbol, value;  coupling=Coupling{FID}(:, :, (2, 1)), amplitude::Union{Function, Nothing}=nothing, ishermitian::Bool=true, modulate::Union{Function, Bool}=false)

Onsite term.

Type alias for `Term{:Onsite, id, V, Int, C<:TermCoupling, A<:TermAmplitude, M<:TermModulate}`.
"""
const Onsite{id, V, C<:TermCoupling, A<:TermAmplitude, M<:TermModulate} = Term{:Onsite, id, V, Int, C, A, M}
@inline function Onsite(id::Symbol, value; coupling=Coupling{FID}(:, :, (2, 1)), amplitude::Union{Function, Nothing}=nothing, ishermitian::Bool=true, modulate::Union{Function, Bool}=false)
    return Term{:Onsite}(id, value, 0, coupling, ishermitian; amplitude=amplitude, modulate=modulate)
end

"""
    Hopping(id::Symbol, value, bondkind; coupling=Coupling{FID}(:, :, (2, 1)), amplitude::Union{Function, Nothing}=nothing, modulate::Union{Function, Bool}=false)

Hopping term.

Type alias for `Term{:Hopping, id, V, B, C<:TermCoupling, A<:TermAmplitude, M<:TermModulate}`.
"""
const Hopping{id, V, B, C<:TermCoupling, A<:TermAmplitude, M<:TermModulate} = Term{:Hopping, id, V, B, C, A, M}
@inline function Hopping(id::Symbol, value, bondkind; coupling=Coupling{FID}(:, :, (2, 1)), amplitude::Union{Function, Nothing}=nothing, modulate::Union{Function, Bool}=false)
    @assert bondkind≠0 "Hopping error: input bondkind (neighbor) cannot be 0. Use `Onsite` instead."
    return Term{:Hopping}(id, value, bondkind, coupling, false; amplitude=amplitude, modulate=modulate)
end

"""
    Pairing(id::Symbol, value, bondkind; coupling=Coupling{FID}(:, :, (1, 1)), amplitude::Union{Function, Nothing}=nothing, modulate::Union{Function, Bool}=false)

Pairing term.

Type alias for `Term{:Pairing, id, V, B, C<:TermCoupling, A<:TermAmplitude, M<:TermModulate}`.
"""
const Pairing{id, V, B, C<:TermCoupling, A<:TermAmplitude, M<:TermModulate} = Term{:Pairing, id, V, B, C, A, M}
@inline function Pairing(id::Symbol, value, bondkind; coupling=Coupling{FID}(:, :, (1, 1)), amplitude::Union{Function, Nothing}=nothing, modulate::Union{Function, Bool}=false)
    return Term{:Pairing}(id, value, bondkind, coupling, false; amplitude=amplitude, modulate=modulate)
end
@inline nambu(::Val{:Pairing}, ::typeof(:), ::Int) = annihilation
function expand!(operators::Operators, term::Pairing, bond::Bond, hilbert::Hilbert; half::Bool=false)
    argtypes = Tuple{Operators, Term, Bond, Hilbert}
    invoke(expand!, argtypes, operators, term, bond, hilbert; half=half)
    length(bond)==2 && invoke(expand!, argtypes, operators, term, reverse(bond), hilbert; half=half)
    return operators
end

"""
    Hubbard(id::Symbol, value; amplitude::Union{Function, Nothing}=nothing, modulate::Union{Function, Bool}=false)

Hubbard term.

Type alias for `Term{:Hubbard, id, V, Int, C<:TermCoupling, A<:TermAmplitude, M<:TermModulate}`.
"""
const Hubbard{id, V, C<:TermCoupling, A<:TermAmplitude, M<:TermModulate} = Term{:Hubbard, id, V, Int, C, A, M}
@inline function Hubbard(id::Symbol, value; amplitude::Union{Function, Nothing}=nothing, modulate::Union{Function, Bool}=false)
    return Term{:Hubbard}(id, value, 0, Coupling{FID}(:, (2, 2, 1, 1), (2, 1, 2, 1)), true; amplitude=amplitude, modulate=modulate)
end

"""
    InterOrbitalInterSpin(id::Symbol, value; amplitude::Union{Function, Nothing}=nothing, modulate::Union{Function, Bool}=false)

Interorbital-interspin term.

Type alias for `Term{:InterOrbitalInterSpin, id, V, Int, C<:TermCoupling, A<:TermAmplitude, M<:TermModulate}`.
"""
const InterOrbitalInterSpin{id, V, C<:TermCoupling, A<:TermAmplitude, M<:TermModulate} = Term{:InterOrbitalInterSpin, id, V, Int, C, A, M}
@inline function InterOrbitalInterSpin(id::Symbol, value; amplitude::Union{Function, Nothing}=nothing, modulate::Union{Function, Bool}=false)
    return Term{:InterOrbitalInterSpin}(id, value, 0, Coupling(@iids(FID(α, σ₁, 2), FID(α, σ₁, 1), FID(β, σ₂, 2), FID(β, σ₂, 1); constraint=α<β && σ₁≠σ₂)), true; amplitude=amplitude, modulate=modulate)
end

"""
    InterOrbitalIntraSpin(id::Symbol, value; amplitude::Union{Function, Nothing}=nothing, modulate::Union{Function, Bool}=false)

Interorbital-intraspin term.

Type alias for `Term{:InterOrbitalIntraSpin, id, V, Int, C<:TermCoupling, A<:TermAmplitude, M<:TermModulate}`.
"""
const InterOrbitalIntraSpin{id, V, C<:TermCoupling, A<:TermAmplitude, M<:TermModulate} = Term{:InterOrbitalIntraSpin, id, V, Int, C, A, M}
@inline function InterOrbitalIntraSpin(id::Symbol, value; amplitude::Union{Function, Nothing}=nothing, modulate::Union{Function, Bool}=false)
    return Term{:InterOrbitalIntraSpin}(id, value, 0, Coupling(@iids(FID(α, σ, 2), FID(α, σ, 1), FID(β, σ, 2), FID(β, σ, 1); constraint=α<β)), true; amplitude=amplitude, modulate=modulate)
end

"""
    SpinFlip(id::Symbol, value; amplitude::Union{Function, Nothing}=nothing, modulate::Union{Function, Bool}=false)

Spin-flip term.

Type alias for `Term{:SpinFlip, id, V, Int, C<:TermCoupling, A<:TermAmplitude, M<:TermModulate}`.
"""
const SpinFlip{id, V, C<:TermCoupling, A<:TermAmplitude, M<:TermModulate} = Term{:SpinFlip, id, V, Int, C, A, M}
@inline function SpinFlip(id::Symbol, value; amplitude::Union{Function, Nothing}=nothing, modulate::Union{Function, Bool}=false)
    return Term{:SpinFlip}(id, value, 0, Coupling(@iids(FID(α, 2, 2), FID(β, 1, 2), FID(α, 1, 1), FID(β, 2, 1); constraint=α<β)), false; amplitude=amplitude, modulate=modulate)
end

"""
    PairHopping(id::Symbol, value; amplitude::Union{Function, Nothing}=nothing, modulate::Union{Function, Bool}=false)

Pair-hopping term.

Type alias for `Term{:PairHopping, id, V, Int, C<:TermCoupling, A<:TermAmplitude, M<:TermModulate}`.
"""
const PairHopping{id, V, C<:TermCoupling, A<:TermAmplitude, M<:TermModulate} = Term{:PairHopping, id, V, Int, C, A, M}
@inline function PairHopping(id::Symbol, value; amplitude::Union{Function, Nothing}=nothing, modulate::Union{Function, Bool}=false)
    return Term{:PairHopping}(id, value, 0, Coupling(@iids(FID(α, 2, 2), FID(α, 1, 2), FID(β, 1, 1), FID(β, 2, 1); constraint=α<β)), false; amplitude=amplitude, modulate=modulate)
end

"""
    Coulomb(id::Symbol, value, bondkind; coupling=Coupling{FID}(:, :, (2, 1))*Coupling{FID}(:, :, (2, 1)), amplitude::Union{Function, Nothing}=nothing, ishermitian::Bool=true, modulate::Union{Function, Bool}=false)

Coulomb term.

Type alias for `Term{:Coulomb, id, V, B, C<:TermCoupling, A<:TermAmplitude, M<:TermModulate}`.
"""
const Coulomb{id, V, B, C<:TermCoupling, A<:TermAmplitude, M<:TermModulate} = Term{:Coulomb, id, V, B, C, A, M}
@inline function Coulomb(id::Symbol, value, bondkind; coupling=Coupling{FID}(:, :, (2, 1))*Coupling{FID}(:, :, (2, 1)), amplitude::Union{Function, Nothing}=nothing, ishermitian::Bool=true, modulate::Union{Function, Bool}=false)
    return Term{:Coulomb}(id, value, bondkind, coupling, ishermitian; amplitude=amplitude, modulate=modulate)
end

"""
    FockTerm

Type alias for `Union{Onsite, Hopping, Pairing, Hubbard, InterOrbitalInterSpin, InterOrbitalIntraSpin, SpinFlip, PairHopping, Coulomb}`.
"""
const FockTerm = Union{Onsite, Hopping, Pairing, Hubbard, InterOrbitalInterSpin, InterOrbitalIntraSpin, SpinFlip, PairHopping, Coulomb}

# SU(2) spin systems
const sidtagmap = Dict(1=>'x', 2=>'y', 3=>'z', 4=>'+', 5=>'-')
const sidseqmap = Dict(v=>k for (k, v) in sidtagmap)
const sidajointmap = Dict('x'=>'x', 'y'=>'y', 'z'=>'z', '+'=>'-', '-'=>'+')
const sidrepmap = Dict('x'=>'ˣ', 'y'=>'ʸ', 'z'=>'ᶻ', '+'=>'⁺', '-'=>'⁻', '0'=>'⁰')
const sidreprevmap = Dict(v=>k for (k, v) in sidrepmap)

## SID
"""
    SID{S, T<:Union{Char, Symbol, typeof(:)}} <: SimpleIID

The spin id.
"""
struct SID{S, T<:Union{Char, Symbol, typeof(:)}} <: SimpleIID
    tag::T
    function SID{S}(tag::Union{Char, Symbol, typeof(:)}) where S
        @assert isa(S, Rational{Int}) && S.den==2 || isa(S, Integer) || S==wildcard "SID error: not supported spin($S)."
        isa(tag, Char) && @assert tag in ('x', 'y', 'z', '+', '-') "SID error: not supported tag($tag)."
        new{S, typeof(tag)}(tag)
    end
end
@inline SID(tag::Union{Char, Symbol, typeof(:)}) = SID{wildcard}(tag)
@inline Base.:(==)(sid₁::SID, sid₂::SID) = totalspin(sid₁)==totalspin(sid₂) && ==(efficientoperations, sid₁, sid₂)
@inline Base.isequal(sid₁::SID, sid₂::SID) = isequal(totalspin(sid₁), totalspin(sid₂)) && isequal(efficientoperations, sid₁, sid₂)
@inline Base.adjoint(sid::SID) = SID{totalspin(sid)}(sidajointmap[sid.tag])
@inline Base.hash(sid::SID, h::UInt) = hash((totalspin(sid), sid.tag), h)
@inline Base.show(io::IO, sid::SID) = @printf io "SID{%s}(%s)" totalspin(sid) repr(sid.tag)
@inline Base.show(io::IO, sid::SID{wildcard}) = @printf io "SID(%s)" repr(sid.tag)
@inline @generated function Base.replace(sid::SID; kwargs...)
    exprs = [:(get(kwargs, $name, getfield(sid, $name))) for name in QuoteNode.(fieldnames(sid))]
    return :(rawtype(typeof(sid)){totalspin(sid)}($(exprs...)))
end
@inline totalspin(sid::SID) = totalspin(typeof(sid))
@inline totalspin(::Type{<:SID{S}}) where S = S
@inline statistics(::Type{<:SID}) = :b
@inline SID(iid::SID, ::CompositeInternal) = iid

### required by Constraint
@inline isconcreteiid(::Type{<:SID{T, Char} where T}) = true
@inline iidtype(::Type{SID}, ::Type{T}) where {T<:Union{Char, Symbol, typeof(:)}} = SID{wildcard, T}
@inline iidtype(::Type{SID{S}}, ::Type{T}) where {S, T<:Union{Char, Symbol, typeof(:)}} = SID{S, T}

### matrix
"""
    matrix(sid::SID{S, Char}, dtype::Type{<:Number}=Complex{Float}) where S -> Matrix{dtype}

Get the matrix representation of a sid.
"""
function matrix(sid::SID{S, Char}, dtype::Type{<:Number}=Complex{Float}) where S
    N = Int(2*S+1)
    result = zeros(dtype, (N, N))
    spin = convert(dtype, S)
    for i = 1:N, j = 1:N
        row, col = N+1-i, N+1-j
        m, n = spin+1-i, spin+1-j
        result[row, col] = (sid.tag == 'x') ? (delta(i+1, j)+delta(i, j+1))*sqrt(spin*(spin+1)-m*n)/2 :
            (sid.tag == 'y') ? (delta(i+1, j)-delta(i, j+1))*sqrt(spin*(spin+1)-m*n)/2im :
            (sid.tag == 'z') ? delta(i, j)*m :
            (sid.tag == '+') ? delta(i+1, j)*sqrt(spin*(spin+1)-m*n) :
            delta(i, j+1)*sqrt(spin*(spin+1)-m*n)
    end
    return result
end

## Spin
"""
    Spin{S} <: SimpleInternal{SID{S, Char}}

The spin internal degrees of freedom.
"""
struct Spin{S} <: SimpleInternal{SID{S, Char}}
    function Spin{S}() where S
        @assert isa(S, Rational{Int}) && S.den==2 || isa(S, Integer) "Spin error: not supported spin($S)."
        new{S}()
    end
end
@inline Base.eltype(::Type{Spin}) = (SID{S, Char} where S)
@inline shape(sp::Spin) = (1:length(sidtagmap),)
@inline Base.CartesianIndex(sid::SID, ::Spin) = CartesianIndex(sidseqmap[sid.tag])
@inline SID(index::CartesianIndex{1}, sp::Spin) = SID{totalspin(sp)}(sidtagmap[index[1]])
@inline Base.summary(io::IO, spin::Spin) = @printf io "%s-element Spin{%s}" length(spin) totalspin(spin)
@inline totalspin(spin::Spin) = totalspin(typeof(spin))
@inline totalspin(::Type{<:Spin{S}}) where S = S
@inline Base.show(io::IO, spin::Spin) = @printf io "%s{%s}()" spin|>typeof|>nameof totalspin(spin)
@inline Base.match(::Type{<:SID{wildcard}}, ::Type{<:Spin{S}}) where S = true
@inline Base.match(::Type{<:SID{S}}, ::Type{<:Spin{S}}) where S = true
@inline Base.match(::Type{<:SID{S₁}}, ::Type{<:Spin{S₂}}) where {S₁, S₂} = false

## LaTeX format output
"""
    script(::Val{:tag}, sid::SID; kwargs...) -> Char

Get the required script of an sid.
"""
@inline script(::Val{:tag}, sid::SID; kwargs...) = sid.tag

"""
    script(::Val{:site}, index::Index{<:SID}; kwargs...) -> Int
    script(attr::Val, index::Index{<:SID}; kwargs...) -> Union{Int, Char}

Get the required script of a spin index.
"""
@inline script(::Val{:site}, index::Index{<:SID}; kwargs...) = index.site
@inline script(attr::Val, index::Index{<:SID}; kwargs...) = script(attr, index.iid; kwargs...)

"""
    latexofspins

The default LaTeX format of a spin index.
"""
const latexofspins = LaTeX{(:tag,), (:site,)}('S')
@inline latexname(::Type{<:Index{<:SID}}) = Symbol("Index{SID}")
@inline latexname(::Type{<:AbstractCompositeIndex{<:Index{<:SID}}}) = Symbol("AbstractCompositeIndex{Index{SID}}")
latexformat(Index{<:SID}, latexofspins)
latexformat(AbstractCompositeIndex{<:Index{<:SID}}, latexofspins)
@inline latexname(::Type{<:SID}) = Symbol("SID")
latexformat(SID, LaTeX{(:tag,), ()}('S'))

## Permutation
"""
    permute(id₁::CompositeIndex{<:Index{SID{S, Char}}}, id₂::CompositeIndex{<:Index{SID{S, Char}}}) where S -> Tuple{Vararg{Operator}}

Permute two spin indexes and get the result.
"""
function permute(id₁::CompositeIndex{<:Index{SID{S, Char}}}, id₂::CompositeIndex{<:Index{SID{S, Char}}}) where S
    @assert id₁.index≠id₂.index || id₁.rcoordinate≠id₂.rcoordinate || id₁.icoordinate≠id₂.icoordinate "permute error: permuted ids should not be equal to each other."
    if id₁.index.site==id₂.index.site && id₁.rcoordinate==id₂.rcoordinate && id₁.icoordinate==id₂.icoordinate
        if id₁.index.iid.tag == 'x'
            id₂.index.iid.tag=='y' && return (Operator(+1im, permutespinindex(id₁, 'z')), Operator(1, id₂, id₁))
            id₂.index.iid.tag=='z' && return (Operator(-1im, permutespinindex(id₁, 'y')), Operator(1, id₂, id₁))
            id₂.index.iid.tag=='+' && return (Operator(-1, permutespinindex(id₁, 'z')), Operator(1, id₂, id₁))
            id₂.index.iid.tag=='-' && return (Operator(+1, permutespinindex(id₁, 'z')), Operator(1, id₂, id₁))
        elseif id₁.index.iid.tag == 'y'
            id₂.index.iid.tag=='x' && return (Operator(-1im, permutespinindex(id₁, 'z')), Operator(1, id₂, id₁))
            id₂.index.iid.tag=='z' && return (Operator(+1im, permutespinindex(id₁, 'x')), Operator(1, id₂, id₁))
            id₂.index.iid.tag=='+' && return (Operator(-1im, permutespinindex(id₁, 'z')), Operator(1, id₂, id₁))
            id₂.index.iid.tag=='-' && return (Operator(-1im, permutespinindex(id₁, 'z')), Operator(1, id₂, id₁))
        elseif id₁.index.iid.tag == 'z'
            id₂.index.iid.tag=='x' && return (Operator(+1im, permutespinindex(id₁, 'y')), Operator(1, id₂, id₁))
            id₂.index.iid.tag=='y' && return (Operator(-1im, permutespinindex(id₁, 'x')), Operator(1, id₂, id₁))
            id₂.index.iid.tag=='+' && return (Operator(+1, id₂), Operator(1, id₂, id₁))
            id₂.index.iid.tag=='-' && return (Operator(-1, id₂), Operator(1, id₂, id₁))
        elseif id₁.index.iid.tag == '+'
            id₂.index.iid.tag=='x' && return (Operator(+1, permutespinindex(id₁, 'z')), Operator(1, id₂, id₁))
            id₂.index.iid.tag=='y' && return (Operator(+1im, permutespinindex(id₁, 'z')), Operator(1, id₂, id₁))
            id₂.index.iid.tag=='z' && return (Operator(-1, id₁), Operator(1, id₂, id₁))
            id₂.index.iid.tag=='-' && return (Operator(+2, permutespinindex(id₁, 'z')), Operator(1, id₂, id₁))
        elseif id₁.index.iid.tag == '-'
            id₂.index.iid.tag=='x' && return (Operator(-1, permutespinindex(id₁, 'z')), Operator(1, id₂, id₁))
            id₂.index.iid.tag=='y' && return (Operator(1im, permutespinindex(id₁, 'z')), Operator(1, id₂, id₁))
            id₂.index.iid.tag=='z' && return (Operator(+1, id₁), Operator(1, id₂, id₁))
            id₂.index.iid.tag=='+' && return (Operator(-2, permutespinindex(id₁, 'z')), Operator(1, id₂, id₁))
        end
    else
        return (Operator(1, id₂, id₁),)
    end
end
@inline permutespinindex(id::CompositeIndex{<:Index{<:SID}}, tag::Char) = replace(id, index=replace(id.index, iid=replace(id.index.iid, tag=tag)))

## Coupling
### required by IIDSpace
@inline shape(iidspace::IIDSpace{<:SID, <:Spin}) = (sidrange(iidspace.iid.tag),)
@inline sidrange(::Union{Symbol, typeof(:)}) = 1:5
@inline sidrange(tag::Char) = (pos=sidseqmap[tag]; pos:pos)

### Coupling
"""
    Coupling{S}(tag::NTuple{N, Char}) where {S<:SID, N}
    Coupling{S}(value, tag::NTuple{N, Char}) where {S<:SID, N}

Construct a `Coupling` among different `SID`s with the given `tag`.
"""
@inline Coupling{S}(tag::NTuple{N, Char}) where {S<:SID, N} = Coupling{S}(1, tag)
@inline Coupling{S}(value, tag::NTuple{N, Char}) where {S<:SID, N} = Coupling(value, ID(S, tag))

### MatrixCoupling
const default_tag = ['x', 'y', 'z']
"""
    MatrixCoupling{S}(tag::AbstractMatrix) where {S<:SID}

Construct a set of `Coupling`s between two `SID`s with the coefficients specified by a matrix.
"""
@inline MatrixCoupling{S}(matrix::AbstractMatrix) where {S<:SID} = MatrixCoupling{S}(Component(default_tag, default_tag, matrix))

### Spin coupling matrix
"""
    heisenberg"" -> SparseMatrixCSC

The Heisenberg coupling matrix.
"""
macro heisenberg_str(str::String)
    str=="" && return SparseMatrixCSC([1 0 0; 0 1 0; 0 0 1])
    error("@heisenberg_str error: wrong input string.")
end

"""
    ising"x" -> SparseMatrixCSC
    ising"y" -> SparseMatrixCSC
    ising"z" -> SparseMatrixCSC

The Ising coupling matrix.
"""
macro ising_str(str::String)
    str=="x" && return SparseMatrixCSC([1 0 0; 0 0 0; 0 0 0])
    str=="y" && return SparseMatrixCSC([0 0 0; 0 1 0; 0 0 0])
    str=="z" && return SparseMatrixCSC([0 0 0; 0 0 0; 0 0 1])
    error("@ising_str error: wrong input string.")
end

"""
    gamma"x" -> SparseMatrixCSC
    gamma"y" -> SparseMatrixCSC
    gamma"z" -> SparseMatrixCSC

The Gamma coupling matrix.
"""
macro gamma_str(str::String)
    str=="x" && return SparseMatrixCSC([0 0 0; 0 0 1; 0 1 0])
    str=="y" && return SparseMatrixCSC([0 0 1; 0 0 0; 1 0 0])
    str=="z" && return SparseMatrixCSC([0 1 0; 0 1 0; 0 0 0])
    error("@gamma_str error: wrong input string.")
end

"""
    dm"x" -> SparseMatrixCSC
    dm"y" -> SparseMatrixCSC
    dm"z" -> SparseMatrixCSC

The DM coupling matrix.
"""
macro dm_str(str::String)
    str=="x" && return SparseMatrixCSC([0 0 0; 0 0 1; 0 -1 0])
    str=="y" && return SparseMatrixCSC([0 0 -1; 0 0 0; 1 0 0])
    str=="z" && return SparseMatrixCSC([0 1 0; 0 -1 0; 0 0 0])
    error("@dm_str error: wrong input string.")
end

## Term
"""
    SpinTerm(id::Symbol, value, bondkind, coupling; amplitude::Union{Function, Nothing}=nothing, modulate::Union{Function, Bool}=false)

Spin term.

Type alias for `Term{:SpinTerm, id, V, B, C<:TermCoupling, A<:TermAmplitude, M<:TermModulate}`.
"""
const SpinTerm{id, V, B, C<:TermCoupling, A<:TermAmplitude, M<:TermModulate} = Term{:SpinTerm, id, V, B, C, A, M}
@inline function SpinTerm(id::Symbol, value, bondkind, coupling; amplitude::Union{Function, Nothing}=nothing, modulate::Union{Function, Bool}=false)
    return Term{:SpinTerm}(id, value, bondkind, coupling, true; amplitude=amplitude, modulate=modulate)
end
@inline function couplingcenters(sc::Coupling, bond::Bond, ::Val{:SpinTerm})
    length(bond)==1 && return ntuple(i->1, Val(rank(sc)))
    length(bond)==2 && return ntuple(i->2-i%2, Val(rank(sc)))
    error("couplingcenters error: not supported for a generic bond containing $(length(bond)) points.")
end

# Phononic systems
## PID
"""
    PID{D<:Union{Char, Symbol, typeof(:)}} <: SimpleIID

The phonon id.
"""
struct PID{D<:Union{Char, Symbol, typeof(:)}} <: SimpleIID
    tag::Char
    direction::D
    function PID(tag::Char, direction::Union{Char, Symbol, typeof(:)})
        @assert tag∈('p', 'u') "PID error: wrong tag($tag)."
        isa(direction, Char) && @assert direction∈('x', 'y', 'z') "PID error: wrong direction($direction)."
        new{typeof(direction)}(tag, direction)
    end
end
@inline Base.adjoint(pid::PID) = pid
@inline statistics(::Type{<:PID}) = :b
@inline PID(iid::PID, ::CompositeInternal) = iid

### required by Constraint
@inline isconcreteiid(::Type{PID{Char}}) = true
@inline iidtype(::Type{PID}, ::Type{Char}, ::Type{D}) where {D<:Union{Char, Symbol, typeof(:)}} = PID{D}

## Phonon
"""
    Phonon <: SimpleInternal{PID{Char}}

The phonon internal degrees of freedom.
"""
struct Phonon <: SimpleInternal{PID{Char}}
    ndirection::Int
    function Phonon(ndirection::Integer)
        @assert ndirection∈(1, 2, 3) "Phonon error: wrong number of directions."
        new(ndirection)
    end
end
@inline shape(pn::Phonon) = (1:2, 1:pn.ndirection)
@inline Base.CartesianIndex(pid::PID{Char}, ::Phonon) = CartesianIndex(pid.tag=='u' ? 1 : 2, Int(pid.direction)-Int('x')+1)
@inline PID(index::CartesianIndex{2}, ::Phonon) = PID(index[1]==1 ? 'u' : 'p', Char(Int('x')+index[2]-1))

## LaTeX format output
"""
    script(::Val{:BD}, pid::PID, l::LaTeX) -> Char
    script(::Val{:direction}, pid::PID; kwargs...) -> Int

Get the required script of an pid.
"""
@inline script(::Val{:BD}, pid::PID, l::LaTeX) = l.body[pid.tag]
@inline script(::Val{:direction}, pid::PID; kwargs...) = pid.direction

"""
    script(::Val{:BD}, index::Index{<:PID}, l::LaTeX) -> Char
    script(::Val{:BD}, index::AbstractCompositeIndex{<:Index{<:PID}}, l::LaTeX) -> Char
    script(::Val{:site}, index::Index{<:PID}; kwargs...) -> Int
    script(::Val{:direction}, index::Index{<:PID}; kwargs...) -> Char

Get the required script of a phonon index.
"""
@inline script(::Val{:BD}, index::Index{<:PID}, l::LaTeX) = l.body[index.iid.tag]
@inline script(::Val{:BD}, index::AbstractCompositeIndex{<:Index{<:PID}}, l::LaTeX) = l.body[getcontent(index, :index).iid.tag]
@inline script(::Val{:site}, index::Index{<:PID}; kwargs...) = index.site
@inline script(::Val{:direction}, index::Index{<:PID}; kwargs...) = index.iid.direction

"""
    latexofphonons

The default LaTeX format of a phonon index.
"""
const latexofphonons = LaTeX{(:direction,), (:site,)}(Dict('p'=>'p', 'u'=>'u'), "", "")
@inline latexname(::Type{<:Index{<:PID}}) = Symbol("Index{PID}")
@inline latexname(::Type{<:AbstractCompositeIndex{<:Index{<:PID}}}) = Symbol("AbstractCompositeIndex{Index{PID}}")
latexformat(Index{<:PID}, latexofphonons)
latexformat(AbstractCompositeIndex{<:Index{<:PID}}, latexofphonons)
@inline latexname(::Type{<:PID}) = Symbol("PID")
latexformat(PID, LaTeX{(:direction,), ()}(Dict('p'=>'p', 'u'=>'u'), "", ""))

## Permutation
"""
    permute(id₁::CompositeIndex{<:Index{PID{Char}}}, id₂::CompositeIndex{<:Index{PID{Char}}}) -> Tuple{Vararg{Operator}}

Permute two phonon indexes and get the result.
"""
function permute(id₁::CompositeIndex{<:Index{PID{Char}}}, id₂::CompositeIndex{<:Index{PID{Char}}})
    if id₁.index.iid.direction==id₂.index.iid.direction && id₁.index.iid.tag≠id₂.index.iid.tag && id₁.rcoordinate==id₂.rcoordinate && id₁.icoordinate==id₂.icoordinate
        if id₁.index.iid.tag=='u'
            return (Operator(1im), Operator(1, id₂, id₁))
        else
            return (Operator(-1im), Operator(1, id₂, id₁))
        end
    else
        return (Operator(1, id₂, id₁),)
    end
end

## Coupling
### required by IIDSpace
@inline function shape(iidspace::IIDSpace{<:PID, Phonon})
    return (iidspace.iid.tag=='u' ? (1:1) : (2:2), pidrange(iidspace.iid.direction, iidspace.internal.ndirection))
end
@inline pidrange(::Union{Symbol, typeof(:)}, ndirection::Int) = 1:ndirection
@inline function pidrange(direction::Char, ndirection::Int)
    direction = Int(direction) - Int('x') + 1
    @assert 0<direction<ndirection+1 "pidrange error: direction out of range."
    return direction:direction
end

### Coupling
"""
    Coupling{PID}(tag::NTuple{N, Char}, direction::Union{NTuple{N, Char}, typeof(:)}) where N
    Coupling{PID}(value, tag::NTuple{N, Char}, direction::Union{NTuple{N, Char}, typeof(:)}) where N

Construct a `Coupling` among different `PID`s with the given `tag` and `direction`.
"""
@inline Coupling{PID}(tag::NTuple{N, Char}, direction::Union{NTuple{N, Char}, typeof(:)}) where N = Coupling{PID}(1, tag, direction)
@inline Coupling{PID}(value, tag::NTuple{N, Char}, direction::Union{NTuple{N, Char}, typeof(:)}) where N = Coupling(value, ID(PID, tag, default(direction, N|>Val)))

### Constraint
@inline Constraint(::NTuple{N, PID{typeof(:)}}) where N = Constraint{N}(Diagonal(:direction))

### MatrixCoupling
const default_u = ['u']
const default_1du = ['x']
const default_2du = ['x', 'y']
const default_3du = ['x', 'y', 'z']
"""
    MatrixCoupling{PID}(matrix::AbstractMatrix)

Construct a set of `Coupling`s corresponding to the dynamical matrix of phonons.
"""
function MatrixCoupling{PID}(matrix::AbstractMatrix)
    @assert size(matrix)[1]∈(1, 2, 3) && size(matrix)[2]∈(1, 2, 3) "MatrixCoupling error: mismatched dimension of input matrix."
    left = size(matrix)[1]==1 ? default_1du : size(matrix)[1]==2 ? default_2du : default_3du
    right = size(matrix)[2]==1 ? default_1du : size(matrix)[2]==2 ? default_2du : default_3du
    return MatrixCoupling{PID}(Component(default_u, default_u, default_matrix), Component(left, right, matrix))
end

### expand
"""
    expand(::Val{:Hooke}, pnc::Coupling{<:Number, NTuple{2, PID{typeof(:)}}}, bond::Bond, hilbert::Hilbert) -> PPExpand

Expand the default phonon potential coupling on a given bond.
"""
function expand(::Val{:Hooke}, pnc::Coupling{<:Number, NTuple{2, PID{typeof(:)}}}, bond::Bond, hilbert::Hilbert)
    R̂ = rcoordinate(bond)/norm(rcoordinate(bond))
    @assert pnc.iids.tags==('u', 'u') "expand error: wrong tags of phonon coupling."
    @assert isapprox(pnc.value, 1, atol=atol, rtol=rtol) "expand error: wrong coefficient of phonon coupling."
    pn₁ = filter(pnc.iids[1], hilbert[bond[1].site])
    pn₂ = filter(pnc.iids[2], hilbert[bond[2].site])
    @assert pn₁.ndirection==pn₂.ndirection==length(R̂) "expand error: mismatched number of directions."
    return PPExpand(R̂, (bond[1], bond[2]))
end
struct PPExpand{N, D<:Number} <: VectorSpace{Operator{D, ID{CompositeIndex{Index{PID{Char}}, SVector{N, D}}, 2}}}
    direction::SVector{N, D}
    points::NTuple{2, Point{N, D}}
end
@inline VectorSpaceStyle(::Type{<:PPExpand}) = VectorSpaceCartesian()
@inline shape(pnce::PPExpand) = (1:length(pnce.direction), 1:length(pnce.direction), 1:4)
function Operator(index::CartesianIndex{3}, pnce::PPExpand)
    direction₁ = Char(Int('x')+index[1]-1)
    direction₂ = Char(Int('x')+index[2]-1)
    coeff = index[3]∈(1, 4) ? 1 : -1
    pos₁, pos₂ = index[3]==1 ? (1, 1) : index[3]==2 ? (1, 2) : index[3]==3 ? (2, 1) : (2, 2)
    index₁ = CompositeIndex(Index(pnce.points[pos₁].site, PID('u', direction₁)), pnce.points[pos₁].rcoordinate, pnce.points[pos₁].icoordinate)
    index₂ = CompositeIndex(Index(pnce.points[pos₂].site, PID('u', direction₂)), pnce.points[pos₂].rcoordinate, pnce.points[pos₂].icoordinate)
    return Operator(pnce.direction[index[1]]*pnce.direction[index[2]]*coeff, index₁, index₂)
end

## Term
"""
    Kinetic(id::Symbol, value; amplitude::Union{Function, Nothing}=nothing, modulate::Union{Function, Bool}=false)

Kinetic energy of phonons.

Type alias for `Term{:Kinetic, id, V, Int, C<:TermCoupling, A<:TermAmplitude, M<:TermModulate}`.
"""
const Kinetic{id, V, C<:TermCoupling, A<:TermAmplitude, M<:TermModulate} = Term{:Kinetic, id, V, Int, C, A, M}
@inline function Kinetic(id::Symbol, value; amplitude::Union{Function, Nothing}=nothing, modulate::Union{Function, Bool}=false)
    return Term{:Kinetic}(id, value, 0, Coupling{PID}(('p', 'p'), :), true; amplitude=amplitude, modulate=modulate)
end

"""
    Hooke(id::Symbol, value, bondkind; amplitude::Union{Function, Nothing}=nothing, modulate::Union{Function, Bool}=false)

Potential energy of phonons by the Hooke's law.

Type alias for `Term{:Hooke, id, V, B, C<:TermCoupling, A<:TermAmplitude, M<:TermModulate}`
"""
const Hooke{id, V, B, C<:TermCoupling, A<:TermAmplitude, M<:TermModulate} = Term{:Hooke, id, V, B, C, A, M}
@inline function Hooke(id::Symbol, value, bondkind; amplitude::Union{Function, Nothing}=nothing, modulate::Union{Function, Bool}=false)
    return Term{:Hooke}(id, value, bondkind, Coupling{PID}(('u', 'u'), :), true; amplitude=amplitude, modulate=modulate)
end

"""
    Elastic(id::Symbol, value, bondkind; coupling, amplitude::Union{Function, Nothing}=nothing, modulate::Union{Function, Bool}=false)

Generic elastic energy of phonons.

Type alias for `Term{:Elastic, id, V, B, C<:TermCoupling, A<:TermAmplitude, M<:TermModulate}`
"""
const Elastic{id, V, B, C<:TermCoupling, A<:TermAmplitude, M<:TermModulate} = Term{:Elastic, id, V, B, C, A, M}
@inline function Elastic(id::Symbol, value, bondkind; coupling, amplitude::Union{Function, Nothing}=nothing, modulate::Union{Function, Bool}=false)
    value, factor = promote(value, 1//2)
    return Term{:Elastic, id}(value, bondkind, TermCoupling(coupling), TermAmplitude(amplitude), true, TermModulate(id, modulate), factor)
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
