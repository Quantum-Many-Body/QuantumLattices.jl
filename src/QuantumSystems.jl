module QuantumSystems

using LinearAlgebra: dot, norm
using Printf: @printf, @sprintf
using SparseArrays: SparseMatrixCSC
using StaticArrays: SVector
using ..DegreesOfFreedom: wildcard, AbstractCompositeIndex, Component, CompositeIID, CompositeIndex, CompositeInternal, Coupling, Hilbert, IIDSpace, Index, SimpleIID, SimpleInternal, Term, TermAmplitude, TermCoupling, TermModulate, @indexes
using ..QuantumLattices: decompose, dtype, kind
using ..QuantumOperators: ID, LaTeX, Operator, OperatorProd, Operators, latexformat
using ..Spatials: Bond, Point, rcoordinate
using ..Toolkit: atol, efficientoperations, rtol, Float, VectorSpace, VectorSpaceCartesian, VectorSpaceStyle, decimaltostr, delta, getcontent, rawtype

import ..DegreesOfFreedom: MatrixCoupling, diagonalizablefields, iidtype, isdefinite, sitestructure, statistics
import ..QuantumLattices: ⊗, ⋅, expand, expand!, permute, rank
import ..QuantumOperators: latexname, matrix, optype, script
import ..Toolkit: shape

# Canonical complex fermionic/bosonic systems
export annihilation, creation, latexofbosons, latexoffermions, latexofparticles, @σ_str, @L_str
export Coulomb, FID, Fock, FockTerm, Hopping, Hubbard, InterOrbitalInterSpin, InterOrbitalIntraSpin, Onsite, PairHopping, Pairing, SpinFlip, isannihilation, iscreation, isnormalordered

# SU(2) spin systems
export latexofspins, SID, Spin, SpinTerm, totalspin, @Γ_str, @DM_str, @Heisenberg_str, @Ising_str

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

@inline default(::Colon) = ":"
@inline default(value::Char) = repr(value)
@inline default(value) = string(value)
@inline default(value::Rational{Int}) = value.den==1 ? repr(value.num) : repr(value)
"""
    FID{T, O<:Union{Int, Symbol, Colon}, S<:Union{Rational{Int}, Symbol, Colon}, N<:Union{Int, Symbol, Colon}} <: SimpleIID

The Fock id.
"""
struct FID{T, O<:Union{Int, Symbol, Colon}, S<:Union{Rational{Int}, Symbol, Colon}, N<:Union{Int, Symbol, Colon}} <: SimpleIID
    orbital::O
    spin::S
    nambu::N
    function FID{T}(orbital::Union{Int, Symbol, Colon}, spin::Union{Rational{Int}, Symbol, Colon}, nambu::Union{Int, Symbol, Colon}) where T
        @assert T∈(:f, :b, wildcard) "FID error: wrong statistics."
        isa(spin, Rational{Int}) && @assert spin.den∈(1, 2) "FID error: wrong spin."
        isa(nambu, Int) && @assert nambu∈(1, 2) "FID error: wrong input nambu($nambu)."
        new{T, typeof(orbital), typeof(spin), typeof(nambu)}(orbital, spin, nambu)
    end
end
@inline FID{T}(orbital::Union{Int, Symbol, Colon}, spin::Int, nambu::Union{Int, Symbol, Colon}) where T = FID{T}(orbital, spin//1, nambu)
@inline FID(orbital::Union{Int, Symbol, Colon}, spin::Union{Rational{Int}, Int, Symbol, Colon}, nambu::Union{Int, Symbol, Colon}) = FID{wildcard}(orbital, spin, nambu)
@inline Base.:(==)(fid₁::FID, fid₂::FID) = statistics(fid₁)==statistics(fid₂) && ==(efficientoperations, fid₁, fid₂)
@inline Base.isequal(fid₁::FID, fid₂::FID) = isequal(statistics(fid₁), statistics(fid₂)) && isequal(efficientoperations, fid₁, fid₂)
@inline Base.hash(fid::FID, h::UInt) = hash((statistics(fid), fid.orbital, fid.spin, fid.nambu), h)
@inline Base.show(io::IO, fid::FID) = @printf io "FID{%s}(%s)" repr(statistics(fid)) join((fid.orbital|>default, fid.spin|>default, fid.nambu|>default), ", ")
@inline Base.show(io::IO, fid::FID{wildcard}) = @printf io "FID(%s)" join((fid.orbital|>default, fid.spin|>default, fid.nambu|>default), ", ")
@inline Base.adjoint(fid::FID{T, <:Union{Int, Symbol, Colon}, <:Union{Rational{Int}, Symbol, Colon}, Int}) where T = FID{T}(fid.orbital, fid.spin, 3-fid.nambu)
@inline @generated function Base.replace(fid::FID; kwargs...)
    exprs = [:(get(kwargs, $name, getfield(fid, $name))) for name in QuoteNode.(fieldnames(fid))]
    return :(rawtype(typeof(fid)){statistics(fid)}($(exprs...)))
end
@inline statistics(::Type{<:FID{T}}) where T = T
@inline FID(iid::FID, ::CompositeInternal) = iid

### requested by Constraint
@inline isdefinite(::Type{<:FID{T, Int, Rational{Int}} where T}) = true
@inline iidtype(::Type{FID}, ::Type{O}, ::Type{S}, ::Type{N}) where {O<:Union{Int, Symbol, Colon}, S<:Union{Rational{Int}, Symbol, Colon}, N<:Union{Int, Symbol, Colon}} = FID{wildcard, O, S, N}
@inline iidtype(::Type{FID{T}}, ::Type{O}, ::Type{S}, ::Type{N}) where {T, O<:Union{Int, Symbol, Colon}, S<:Union{Rational{Int}, Symbol, Colon}, N<:Union{Int, Symbol, Colon}} = FID{T, O, S, N}

"""
    isannihilation(fid::FID) -> Bool
    isannihilation(index::Index) -> Bool
    isannihilation(index::AbstractCompositeIndex) -> Bool

Judge whether the nambu index is `annihilation`.
"""
@inline isannihilation(fid::FID) = fid.nambu==annihilation
@inline isannihilation(index::Index) = isannihilation(index.iid)
@inline isannihilation(index::AbstractCompositeIndex) = isannihilation(getcontent(index, :index))

"""
    iscreation(fid::FID) -> Bool
    iscreation(index::Index) -> Bool
    iscreation(index::AbstractCompositeIndex) -> Bool

Judge whether the nambu index is `creation`.
"""
@inline iscreation(fid::FID) = fid.nambu==creation
@inline iscreation(index::Index) = iscreation(index.iid)
@inline iscreation(index::AbstractCompositeIndex) = iscreation(getcontent(index, :index))

## Fock
"""
    Fock{T} <: SimpleInternal{FID{T, Int, Rational{Int}, Int}}

The Fock internal degrees of freedom.
"""
struct Fock{T} <: SimpleInternal{FID{T, Int, Rational{Int}, Int}}
    norbital::Int
    nspin::Int
    function Fock{T}(norbital::Int, nspin::Int) where T
        @assert T∈(:f, :b) "Fock error: wrong statistics."
        new{T}(norbital, nspin)
    end
end
@inline Base.eltype(::Type{Fock}) = (FID{T, Int, Rational{Int}, Int} where T)
@inline shape(fock::Fock) = (1:fock.norbital, 1:fock.nspin, 1:2)
@inline Base.CartesianIndex(fid::FID{T}, fock::Fock{T}) where T = CartesianIndex(fid.orbital, Int(fid.spin+(fock.nspin-1)//2)+1, fid.nambu)
@inline FID(index::CartesianIndex{3}, fock::Fock) = FID{statistics(fock)}(index[1], index[2]-1-(fock.nspin-1)//2, index[3])
@inline Base.summary(io::IO, fock::Fock) = @printf io "%s-element Fock{%s}" length(fock) repr(statistics(fock))
@inline Base.show(io::IO, fock::Fock) = @printf io "%s{%s}(%s)" fock|>typeof|>nameof repr(statistics(fock)) join(("$name=$(getfield(fock, name))" for name in fock|>typeof|>fieldnames), ", ")
@inline Base.match(::Type{<:FID{wildcard}}, ::Type{<:Fock{T}}) where T = true
@inline Base.match(::Type{<:FID{T}}, ::Type{<:Fock{T}}) where T = true
@inline Base.match(::Type{<:FID{T₁}}, ::Type{<:Fock{T₂}}) where {T₁, T₂} = false

## LaTeX format output
"""
    script(::Val{:orbital}, fid::FID; kwargs...) -> String
    script(::Val{:spint}, fid::FID; kwargs...) -> String
    script(::Val{:spinsym}, fid::FID; kwargs...) -> String
    script(::Val{:nambu}, fid::FID; kwargs...) -> String

Get the requested script of an fid.
"""
@inline script(::Val{:orbital}, fid::FID; kwargs...) = default(fid.orbital)
@inline script(::Val{:spin}, fid::FID; kwargs...) = default(fid.spin)
@inline script(::Val{:spinsym}, fid::FID; kwargs...) = fid.spin==-1//2 ? "↓" : fid.spin==1//2 ? "↑" : fid.spin==0 ? "" : default(fid.spin)
@inline script(::Val{:nambu}, fid::FID; kwargs...) = iscreation(fid) ? "\\dagger" : isannihilation(fid) ? "" : default(fid.nambu)

"""
    latexoffermions

The default LaTeX format of a fermionic index.
"""
const latexoffermions = LaTeX{(:nambu,), (:site, :orbital, :spinsym)}('c')
@inline latexname(::Type{<:Index{<:Union{Int, Colon}, <:FID{:f}}}) = Symbol("Index{Union{Int, Colon}, FID{:f}}")
@inline latexname(::Type{<:AbstractCompositeIndex{<:Index{<:Union{Int, Colon}, <:FID{:f}}}}) = Symbol("AbstractCompositeIndex{Index{Union{Int, Colon}, FID{:f}}}")
latexformat(Index{<:Union{Int, Colon}, <:FID{:f}}, latexoffermions)
latexformat(AbstractCompositeIndex{<:Index{<:Union{Int, Colon}, <:FID{:f}}}, latexoffermions)
@inline latexname(::Type{<:FID{:f}}) = Symbol("FID{:f}")
latexformat(FID{:f}, LaTeX{(:nambu,), (:orbital, :spinsym)}('c'))

"""
    latexofbosons

The default LaTeX format of a bosonic index.
"""
const latexofbosons = LaTeX{(:nambu,), (:site, :orbital, :spinsym)}('b')
@inline latexname(::Type{<:Index{<:Union{Int, Colon}, <:FID{:b}}}) = Symbol("Index{Union{Int, Colon}, FID{:b}}")
@inline latexname(::Type{<:AbstractCompositeIndex{<:Index{<:Union{Int, Colon}, <:FID{:b}}}}) = Symbol("AbstractCompositeIndex{Index{Union{Int, Colon}, FID{:b}}}")
latexformat(Index{<:Union{Int, Colon}, <:FID{:b}}, latexofbosons)
latexformat(AbstractCompositeIndex{<:Index{<:Union{Int, Colon}, <:FID{:b}}}, latexofbosons)
@inline latexname(::Type{<:FID{:b}}) = Symbol("FID{:b}")
latexformat(FID{:b}, LaTeX{(:nambu,), (:orbital, :spinsym)}('b'))

"""
    latexofparticles

The default LaTeX format of a wildcard Fock index.
"""
const latexofparticles = LaTeX{(:nambu,), (:site, :orbital, :spinsym)}('f')
@inline latexname(::Type{<:Index{<:Union{Int, Colon}, <:FID{wildcard}}}) = Symbol("Index{Union{Int, Colon}, FID}")
@inline latexname(::Type{<:AbstractCompositeIndex{<:Index{<:Union{Int, Colon}, <:FID{wildcard}}}}) = Symbol("AbstractCompositeIndex{Index{Union{Int, Colon}, FID}}")
latexformat(Index{<:Union{Int, Colon}, <:FID{wildcard}}, latexofparticles)
latexformat(AbstractCompositeIndex{<:Index{<:Union{Int, Colon}, <:FID{wildcard}}}, latexofparticles)
@inline latexname(::Type{<:FID{wildcard}}) = Symbol("FID")
latexformat(FID{wildcard}, LaTeX{(:nambu,), (:orbital, :spinsym)}('f'))

## Boundary
"""
    angle(id::CompositeIndex{<:Index{Int, <:FID}}, vectors::AbstractVector{<:AbstractVector{<:Number}}, values::AbstractVector{<:Number}) -> Complex{<:Number}

Get the twist phase corresponding to a Fock index.
"""
function Base.angle(id::CompositeIndex{<:Index{Int, <:FID}}, vectors::AbstractVector{<:AbstractVector{<:Number}}, values::AbstractVector{<:Number})
    datatype = promote_type(eltype(values), Float)
    phase = length(vectors)==1 ? 2*convert(datatype, pi)*dot(decompose(id.icoordinate, vectors[1]), values) :
            length(vectors)==2 ? 2*convert(datatype, pi)*dot(decompose(id.icoordinate, vectors[1], vectors[2]), values) :
            length(vectors)==3 ? 2*convert(datatype, pi)*dot(decompose(id.icoordinate, vectors[1], vectors[2], vectors[3]), values) :
            error("angle error: not supported number of input basis vectors.")
    return id.index.iid.nambu==annihilation ? phase : id.index.iid.nambu==creation ? -phase : error("angle error: not supported Fock index.")
end

## New methods
"""
    isnormalordered(opt::Operator{<:Number, <:ID{AbstractCompositeIndex{<:Index{Int, <:FID}}}}) -> Bool

Judge whether an operator is normal ordered.
"""
function isnormalordered(opt::Operator{<:Number, <:ID{AbstractCompositeIndex{<:Index{Int, <:FID}}}})
    flag = true
    for i = 1:rank(opt)
        flag && (opt.id[i].index.iid.nambu == annihilation) && (flag = false)
        flag || (opt.id[i].index.iid.nambu == creation) && return false
    end
    return true
end

"""
    *(
        f1::Operator{<:Number, <:ID{AbstractCompositeIndex{<:Index{Int, <:FID{:f}}}}},
        f2::Operator{<:Number, <:ID{AbstractCompositeIndex{<:Index{Int, <:FID{:f}}}}}
    ) -> Union{Nothing, Operator}

Get the multiplication of two fermionic Fock operators.
"""
@inline function Base.:*(
        f1::Operator{<:Number, <:ID{AbstractCompositeIndex{<:Index{Int, <:FID{:f}}}}},
        f2::Operator{<:Number, <:ID{AbstractCompositeIndex{<:Index{Int, <:FID{:f}}}}}
        )
    rank(f1)>0 && rank(f2)>0 && f1.id[end]==f2.id[1] && return nothing
    return invoke(*, Tuple{OperatorProd, OperatorProd}, f1, f2)
end

## Permutation
"""
    permute(id₁::CompositeIndex{<:Index{Int, <:FID{:f}}}, id₂::CompositeIndex{<:Index{Int, <:FID{:f}}}) -> Tuple{Vararg{Operator}}

Permute two fermionic indexes and get the result.
"""
function permute(id₁::CompositeIndex{<:Index{Int, <:FID{:f}}}, id₂::CompositeIndex{<:Index{Int, <:FID{:f}}})
    @assert id₁.index≠id₂.index || id₁.rcoordinate≠id₂.rcoordinate || id₁.icoordinate≠id₂.icoordinate "permute error: permuted ids should not be equal to each other."
    if id₁.index'==id₂.index && id₁.rcoordinate==id₂.rcoordinate && id₁.icoordinate==id₂.icoordinate
        return (Operator(1), Operator(-1, id₂, id₁))
    else
        return (Operator(-1, id₂, id₁),)
    end
end

"""
    permute(id₁::CompositeIndex{<:Index{Int, <:FID{:b}}}, id₂::CompositeIndex{<:Index{Int, <:FID{:b}}}) -> Tuple{Vararg{Operator}}

Permute two bosonic indexes and get the result.
"""
function permute(id₁::CompositeIndex{<:Index{Int, <:FID{:b}}}, id₂::CompositeIndex{<:Index{Int, <:FID{:b}}})
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
@inline fidshape(::Union{Symbol, Colon}, n::Int) = 1:n
@inline fidshape(v::Int, n::Int) = ((@assert 0<v<n+1 "shape error: out of range."); v:v)
@inline function fidshape(v::Rational{Int}, n::Int)
    @assert abs(v)<=(n-1)//2 "shape error: out of range."
    index = Int(v+(n-1)//2)+1
    return index:index
end

### requested by expand of Coupling
@inline function CompositeIID(coupling::Coupling{V, I}, info::Val=Val(:term)) where {V, I<:ID{<:Index{<:Union{Int, Colon}, <:FID}}}
    return CompositeIID(map((fid::FID, order::Int)->FID{statistics(fid)}(fid.orbital, fid.spin, nambu(info, fid.nambu, order)), coupling.indexes.iids, ntuple(i->i, Val(rank(I)))))
end
@inline nambu(::Val, ::Colon, order::Int) = order%2==1 ? creation : annihilation
@inline nambu(::Val, nambu::Int, ::Int) = nambu

### diagonalizablefields
@inline diagonalizablefields(::Type{<:Index{<:Union{Int, Colon}, <:FID}}) = (:orbital, :spin)

### MatrixCoupling
const default_matrix = SparseMatrixCSC(hcat(1))
"""
    MatrixCoupling(sites::Union{NTuple{2, Int}, Colon}, ::Type{F}, orbital::Union{AbstractMatrix, Colon}, spin::Union{AbstractMatrix, Colon}, nambu::Union{AbstractMatrix, Colon}) where {F<:FID}

Construct a set of `Coupling`s between two `Index{<:Union{Int, Colon}, <:FID}`s with the coefficients specified by matrices acting on separated internal spaces.
"""
@inline function MatrixCoupling(sites::Union{NTuple{2, Int}, Colon}, ::Type{F}, orbital::Union{AbstractMatrix, Colon}, spin::Union{AbstractMatrix, Colon}, nambu::Union{AbstractMatrix, Colon}) where {F<:FID}
    return MatrixCoupling(sites, F, Component(F, Val(:orbital), orbital), Component(F, Val(:spin), spin), Component(F, Val(:nambu), nambu))
end
@inline Component(::Type{<:FID}, ::Val, ::Colon) = Component(SVector(:), SVector(:), default_matrix)
@inline Component(::Type{<:FID}, ::Val{:orbital}, matrix::AbstractMatrix) = Component(1:size(matrix)[1], 1:size(matrix)[2], matrix)
@inline Component(::Type{<:FID}, ::Val{:spin}, matrix::AbstractMatrix) = Component((size(matrix)[1]-1)//2:-1:(1-size(matrix)[1])//2, (size(matrix)[2]-1)//2:-1:(1-size(matrix)[2])//2, matrix)
@inline Component(::Type{<:FID}, ::Val{:nambu}, matrix::AbstractMatrix) = (@assert size(matrix)==(2, 2) "Component error: for nambu subspace, the input matrix must be 2×2."; Component(1:1:2, 2:-1:1, matrix))

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
    Onsite(id::Symbol, value, coupling=Coupling(Index(:, FID(:, :, :)), Index(:, FID(:, :, :))); ishermitian::Bool=true, amplitude::Union{Function, Nothing}=nothing, modulate::Union{Function, Bool}=true)

Onsite term.

Type alias for `Term{:Onsite, id, V, Int, C<:TermCoupling, A<:TermAmplitude, M<:TermModulate}`.
"""
const Onsite{id, V, C<:TermCoupling, A<:TermAmplitude, M<:TermModulate} = Term{:Onsite, id, V, Int, C, A, M}
@inline function Onsite(id::Symbol, value, coupling=Coupling(Index(:, FID(:, :, :)), Index(:, FID(:, :, :))); ishermitian::Bool=true, amplitude::Union{Function, Nothing}=nothing, modulate::Union{Function, Bool}=true)
    return Term{:Onsite}(id, value, 0, coupling, ishermitian; amplitude=amplitude, modulate=modulate)
end

"""
    Hopping(id::Symbol, value, bondkind, coupling=Coupling(Index(:, FID(:, :, :)), Index(:, FID(:, :, :))); amplitude::Union{Function, Nothing}=nothing, modulate::Union{Function, Bool}=true)

Hopping term.

Type alias for `Term{:Hopping, id, V, B, C<:TermCoupling, A<:TermAmplitude, M<:TermModulate}`.
"""
const Hopping{id, V, B, C<:TermCoupling, A<:TermAmplitude, M<:TermModulate} = Term{:Hopping, id, V, B, C, A, M}
@inline function Hopping(id::Symbol, value, bondkind, coupling=Coupling(Index(:, FID(:, :, :)), Index(:, FID(:, :, :))); amplitude::Union{Function, Nothing}=nothing, modulate::Union{Function, Bool}=true)
    @assert bondkind≠0 "Hopping error: input bondkind (neighbor) cannot be 0. Use `Onsite` instead."
    return Term{:Hopping}(id, value, bondkind, coupling, false; amplitude=amplitude, modulate=modulate)
end

"""
    Pairing(id::Symbol, value, bondkind, coupling; amplitude::Union{Function, Nothing}=nothing, modulate::Union{Function, Bool}=true)

Pairing term.

Type alias for `Term{:Pairing, id, V, B, C<:TermCoupling, A<:TermAmplitude, M<:TermModulate}`.
"""
const Pairing{id, V, B, C<:TermCoupling, A<:TermAmplitude, M<:TermModulate} = Term{:Pairing, id, V, B, C, A, M}
@inline function Pairing(id::Symbol, value, bondkind, coupling; amplitude::Union{Function, Nothing}=nothing, modulate::Union{Function, Bool}=true)
    return Term{:Pairing}(id, value, bondkind, coupling, false; amplitude=amplitude, modulate=modulate)
end
@inline nambu(::Val{:Pairing}, ::Colon, ::Int) = annihilation
function expand!(operators::Operators, term::Pairing, bond::Bond, hilbert::Hilbert; half::Bool=false)
    argtypes = Tuple{Operators, Term, Bond, Hilbert}
    invoke(expand!, argtypes, operators, term, bond, hilbert; half=half)
    length(bond)==2 && invoke(expand!, argtypes, operators, term, reverse(bond), hilbert; half=half)
    return operators
end

"""
    Hubbard(id::Symbol, value; amplitude::Union{Function, Nothing}=nothing, modulate::Union{Function, Bool}=true)

Hubbard term.

Type alias for `Term{:Hubbard, id, V, Int, C<:TermCoupling, A<:TermAmplitude, M<:TermModulate}`.
"""
const Hubbard{id, V, C<:TermCoupling, A<:TermAmplitude, M<:TermModulate} = Term{:Hubbard, id, V, Int, C, A, M}
@inline function Hubbard(id::Symbol, value; amplitude::Union{Function, Nothing}=nothing, modulate::Union{Function, Bool}=true)
    return Term{:Hubbard}(id, value, 0, Coupling(:, FID, :, (1//2, 1//2, -1//2, -1//2), (2, 1, 2, 1)), true; amplitude=amplitude, modulate=modulate)
end

"""
    InterOrbitalInterSpin(id::Symbol, value; amplitude::Union{Function, Nothing}=nothing, modulate::Union{Function, Bool}=true)

Interorbital-interspin term.

Type alias for `Term{:InterOrbitalInterSpin, id, V, Int, C<:TermCoupling, A<:TermAmplitude, M<:TermModulate}`.
"""
const InterOrbitalInterSpin{id, V, C<:TermCoupling, A<:TermAmplitude, M<:TermModulate} = Term{:InterOrbitalInterSpin, id, V, Int, C, A, M}
@inline function InterOrbitalInterSpin(id::Symbol, value; amplitude::Union{Function, Nothing}=nothing, modulate::Union{Function, Bool}=true)
    return Term{:InterOrbitalInterSpin}(
        id, value, 0, Coupling(@indexes(Index(:, FID(α, σ₁, 2)), Index(:, FID(α, σ₁, 1)), Index(:, FID(β, σ₂, 2)), Index(:, FID(β, σ₂, 1)); constraint=α<β && σ₁≠σ₂)), true;
        amplitude=amplitude, modulate=modulate
    )
end

"""
    InterOrbitalIntraSpin(id::Symbol, value; amplitude::Union{Function, Nothing}=nothing, modulate::Union{Function, Bool}=true)

Interorbital-intraspin term.

Type alias for `Term{:InterOrbitalIntraSpin, id, V, Int, C<:TermCoupling, A<:TermAmplitude, M<:TermModulate}`.
"""
const InterOrbitalIntraSpin{id, V, C<:TermCoupling, A<:TermAmplitude, M<:TermModulate} = Term{:InterOrbitalIntraSpin, id, V, Int, C, A, M}
@inline function InterOrbitalIntraSpin(id::Symbol, value; amplitude::Union{Function, Nothing}=nothing, modulate::Union{Function, Bool}=true)
    return Term{:InterOrbitalIntraSpin}(
        id, value, 0, Coupling(@indexes(Index(:, FID(α, σ, 2)), Index(:, FID(α, σ, 1)), Index(:, FID(β, σ, 2)), Index(:, FID(β, σ, 1)); constraint=α<β)), true;
        amplitude=amplitude, modulate=modulate
    )
end

"""
    SpinFlip(id::Symbol, value; amplitude::Union{Function, Nothing}=nothing, modulate::Union{Function, Bool}=true)

Spin-flip term.

Type alias for `Term{:SpinFlip, id, V, Int, C<:TermCoupling, A<:TermAmplitude, M<:TermModulate}`.
"""
const SpinFlip{id, V, C<:TermCoupling, A<:TermAmplitude, M<:TermModulate} = Term{:SpinFlip, id, V, Int, C, A, M}
@inline function SpinFlip(id::Symbol, value; amplitude::Union{Function, Nothing}=nothing, modulate::Union{Function, Bool}=true)
    return Term{:SpinFlip}(
        id, value, 0, Coupling(@indexes(Index(:, FID(α, 1//2, 2)), Index(:, FID(β, -1//2, 2)), Index(:, FID(α, -1//2, 1)), Index(:, FID(β, 1//2, 1)); constraint=α<β)), false;
        amplitude=amplitude, modulate=modulate
    )
end

"""
    PairHopping(id::Symbol, value; amplitude::Union{Function, Nothing}=nothing, modulate::Union{Function, Bool}=true)

Pair-hopping term.

Type alias for `Term{:PairHopping, id, V, Int, C<:TermCoupling, A<:TermAmplitude, M<:TermModulate}`.
"""
const PairHopping{id, V, C<:TermCoupling, A<:TermAmplitude, M<:TermModulate} = Term{:PairHopping, id, V, Int, C, A, M}
@inline function PairHopping(id::Symbol, value; amplitude::Union{Function, Nothing}=nothing, modulate::Union{Function, Bool}=true)
    return Term{:PairHopping}(
        id, value, 0, Coupling(@indexes(Index(:, FID(α, 1//2, 2)), Index(:, FID(α, -1//2, 2)), Index(:, FID(β, -1//2, 1)), Index(:, FID(β, 1//2, 1)); constraint=α<β)), false;
        amplitude=amplitude, modulate=modulate
    )
end

"""
    Coulomb(
        id::Symbol, value, bondkind, coupling=Coupling(Index(:, FID(:, :, :)), Index(:, FID(:, :, :)))^2;
        ishermitian::Bool=true,
        amplitude::Union{Function, Nothing}=nothing,
        modulate::Union{Function, Bool}=true
    )

Coulomb term.

Type alias for `Term{:Coulomb, id, V, B, C<:TermCoupling, A<:TermAmplitude, M<:TermModulate}`.
"""
const Coulomb{id, V, B, C<:TermCoupling, A<:TermAmplitude, M<:TermModulate} = Term{:Coulomb, id, V, B, C, A, M}
@inline function Coulomb(
    id::Symbol, value, bondkind, coupling=Coupling(Index(:, FID(:, :, :)), Index(:, FID(:, :, :)))^2;
    ishermitian::Bool=true,
    amplitude::Union{Function, Nothing}=nothing,
    modulate::Union{Function, Bool}=true
)
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
    SID{S, T<:Union{Char, Symbol, Colon}} <: SimpleIID

The spin id.
"""
struct SID{S, T<:Union{Char, Symbol, Colon}} <: SimpleIID
    tag::T
    function SID{S}(tag::Union{Char, Symbol, Colon}) where S
        @assert isa(S, Rational{Int}) && S.den==2 || isa(S, Integer) || S==wildcard "SID error: not supported spin($S)."
        isa(tag, Char) && @assert tag in ('x', 'y', 'z', '+', '-') "SID error: not supported tag($tag)."
        new{S, typeof(tag)}(tag)
    end
end
@inline SID(tag::Union{Char, Symbol, Colon}) = SID{wildcard}(tag)
@inline Base.:(==)(sid₁::SID, sid₂::SID) = totalspin(sid₁)==totalspin(sid₂) && ==(efficientoperations, sid₁, sid₂)
@inline Base.isequal(sid₁::SID, sid₂::SID) = isequal(totalspin(sid₁), totalspin(sid₂)) && isequal(efficientoperations, sid₁, sid₂)
@inline Base.adjoint(sid::SID) = SID{totalspin(sid)}(sidajointmap[sid.tag])
@inline Base.hash(sid::SID, h::UInt) = hash((totalspin(sid), sid.tag), h)
@inline Base.show(io::IO, sid::SID) = @printf io "SID{%s}(%s)" totalspin(sid) default(sid.tag)
@inline Base.show(io::IO, sid::SID{wildcard}) = @printf io "SID(%s)" default(sid.tag)
@inline @generated function Base.replace(sid::SID; kwargs...)
    exprs = [:(get(kwargs, $name, getfield(sid, $name))) for name in QuoteNode.(fieldnames(sid))]
    return :(rawtype(typeof(sid)){totalspin(sid)}($(exprs...)))
end
@inline totalspin(sid::SID) = totalspin(typeof(sid))
@inline totalspin(::Type{<:SID{S}}) where S = S
@inline statistics(::Type{<:SID}) = :b
@inline SID(iid::SID, ::CompositeInternal) = iid

### requested by Constraint
@inline isdefinite(::Type{<:SID{T, Char} where T}) = true
@inline iidtype(::Type{SID}, ::Type{T}) where {T<:Union{Char, Symbol, Colon}} = SID{wildcard, T}
@inline iidtype(::Type{SID{S}}, ::Type{T}) where {S, T<:Union{Char, Symbol, Colon}} = SID{S, T}

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
    script(::Val{:tag}, sid::SID; kwargs...) -> String

Get the requested script of an sid.
"""
@inline script(::Val{:tag}, sid::SID; kwargs...) = sid.tag==(:) ? ":" : string(sid.tag)

"""
    latexofspins

The default LaTeX format of a spin index.
"""
const latexofspins = LaTeX{(:tag,), (:site,)}('S')
@inline latexname(::Type{<:Index{<:Union{Int, Colon}, <:SID}}) = Symbol("Index{Union{Int, Colon}, SID}")
@inline latexname(::Type{<:AbstractCompositeIndex{<:Index{<:Union{Int, Colon}, <:SID}}}) = Symbol("AbstractCompositeIndex{Index{Union{Int, Colon}, SID}}")
latexformat(Index{<:Union{Int, Colon}, <:SID}, latexofspins)
latexformat(AbstractCompositeIndex{<:Index{<:Union{Int, Colon}, <:SID}}, latexofspins)
@inline latexname(::Type{<:SID}) = Symbol("SID")
latexformat(SID, LaTeX{(:tag,), ()}('S'))

## Permutation
"""
    permute(id₁::CompositeIndex{<:Index{Int, SID{S, Char}}}, id₂::CompositeIndex{<:Index{Int, SID{S, Char}}}) where S -> Tuple{Vararg{Operator}}

Permute two spin indexes and get the result.
"""
function permute(id₁::CompositeIndex{<:Index{Int, SID{S, Char}}}, id₂::CompositeIndex{<:Index{Int, SID{S, Char}}}) where S
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
@inline permutespinindex(id::CompositeIndex{<:Index{Int, <:SID}}, tag::Char) = replace(id, index=replace(id.index, iid=replace(id.index.iid, tag=tag)))

## Coupling
### requested by IIDSpace
@inline shape(iidspace::IIDSpace{<:SID, <:Spin}) = (sidrange(iidspace.iid.tag),)
@inline sidrange(::Union{Symbol, Colon}) = 1:5
@inline sidrange(tag::Char) = (pos=sidseqmap[tag]; pos:pos)

### MatrixCoupling
"""
    MatrixCoupling(sites::Union{NTuple{2, Int}, Colon}, ::Type{S}, matrix::AbstractMatrix) where {S<:SID}

Construct a set of `Coupling`s between two `Index{<:Union{Int, Colon}, <:SID}`s with the coefficients specified by a matrix.
"""
@inline MatrixCoupling(sites::Union{NTuple{2, Int}, Colon}, ::Type{S}, matrix::AbstractMatrix) where {S<:SID} = MatrixCoupling(sites, S, Component(SVector('x', 'y', 'z'), SVector('x', 'y', 'z'), matrix))

### Spin coupling matrix
"""
    Heisenberg"" => SparseMatrixCSC([1 0 0; 0 1 0; 0 0 1])

The Heisenberg coupling matrix.
"""
macro Heisenberg_str(str::String)
    str=="" && return SparseMatrixCSC([1 0 0; 0 1 0; 0 0 1])
    error("@Heisenberg_str error: wrong input string.")
end

"""
    Ising"x" => SparseMatrixCSC([1 0 0; 0 0 0; 0 0 0])
    Ising"y" => SparseMatrixCSC([0 0 0; 0 1 0; 0 0 0])
    Ising"z" => SparseMatrixCSC([0 0 0; 0 0 0; 0 0 1])

The Ising coupling matrix.
"""
macro Ising_str(str::String)
    str=="x" && return SparseMatrixCSC([1 0 0; 0 0 0; 0 0 0])
    str=="y" && return SparseMatrixCSC([0 0 0; 0 1 0; 0 0 0])
    str=="z" && return SparseMatrixCSC([0 0 0; 0 0 0; 0 0 1])
    error("@Ising_str error: wrong input string.")
end

"""
    Γ"x" => SparseMatrixCSC([0 0 0; 0 0 1; 0 1 0])
    Γ"y" => SparseMatrixCSC([0 0 1; 0 0 0; 1 0 0])
    Γ"z" => SparseMatrixCSC([0 1 0; 1 0 0; 0 0 0])

The Γ coupling matrix.
"""
macro Γ_str(str::String)
    str=="x" && return SparseMatrixCSC([0 0 0; 0 0 1; 0 1 0])
    str=="y" && return SparseMatrixCSC([0 0 1; 0 0 0; 1 0 0])
    str=="z" && return SparseMatrixCSC([0 1 0; 1 0 0; 0 0 0])
    error("@Γ_str error: wrong input string.")
end

"""
    DM"x" => SparseMatrixCSC([0 0 0; 0 0 1; 0 -1 0])
    DM"y" => SparseMatrixCSC([0 0 -1; 0 0 0; 1 0 0])
    DM"z" => SparseMatrixCSC([0 1 0; -1 0 0; 0 0 0])

The DM coupling matrix.
"""
macro DM_str(str::String)
    str=="x" && return SparseMatrixCSC([0 0 0; 0 0 1; 0 -1 0])
    str=="y" && return SparseMatrixCSC([0 0 -1; 0 0 0; 1 0 0])
    str=="z" && return SparseMatrixCSC([0 1 0; -1 0 0; 0 0 0])
    error("@DM_str error: wrong input string.")
end

## Term
"""
    SpinTerm(id::Symbol, value, bondkind, coupling; amplitude::Union{Function, Nothing}=nothing, modulate::Union{Function, Bool}=true)

Spin term.

Type alias for `Term{:SpinTerm, id, V, B, C<:TermCoupling, A<:TermAmplitude, M<:TermModulate}`.
"""
const SpinTerm{id, V, B, C<:TermCoupling, A<:TermAmplitude, M<:TermModulate} = Term{:SpinTerm, id, V, B, C, A, M}
@inline function SpinTerm(id::Symbol, value, bondkind, coupling; amplitude::Union{Function, Nothing}=nothing, modulate::Union{Function, Bool}=true)
    return Term{:SpinTerm}(id, value, bondkind, coupling, true; amplitude=amplitude, modulate=modulate)
end
@inline function sitestructure(::Val{:SpinTerm}, ::Val{termrank}, bondlength::Integer) where termrank
    bondlength==1 && return ntuple(i->1, Val(termrank))
    bondlength==2 && return ntuple(i->2-i%2, Val(termrank))
    error("sitestructure error: not supported for a generic bond containing $bondlength points.")
end

# Phononic systems
## PID
"""
    PID{D<:Union{Char, Symbol, Colon}} <: SimpleIID

The phonon id.
"""
struct PID{D<:Union{Char, Symbol, Colon}} <: SimpleIID
    tag::Char
    direction::D
    function PID(tag::Char, direction::Union{Char, Symbol, Colon})
        @assert tag∈('p', 'u') "PID error: wrong tag($tag)."
        isa(direction, Char) && @assert direction∈('x', 'y', 'z') "PID error: wrong direction($direction)."
        new{typeof(direction)}(tag, direction)
    end
end
@inline Base.adjoint(pid::PID) = pid
@inline statistics(::Type{<:PID}) = :b
@inline Base.show(io::IO, pid::PID) = @printf io "PID(%s, %s)" repr(pid.tag) default(pid.direction)
@inline PID(iid::PID, ::CompositeInternal) = iid

### requested by Constraint
@inline isdefinite(::Type{PID{Char}}) = true
@inline iidtype(::Type{PID}, ::Type{Char}, ::Type{D}) where {D<:Union{Char, Symbol, Colon}} = PID{D}

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
    script(::Val{:direction}, pid::PID; kwargs...) -> String

Get the requested script of an pid.
"""
@inline script(::Val{:direction}, pid::PID; kwargs...) = pid.direction==(:) ? ":" : string(pid.direction)

"""
    script(::Val{:BD}, pid::PID, l::LaTeX) -> String
    script(::Val{:BD}, index::Index{<:Union{Int, Colon}, <:PID}, l::LaTeX) -> String
    script(::Val{:BD}, index::AbstractCompositeIndex{<:Index{<:Union{Int, Colon}, <:PID}}, l::LaTeX) -> String

Get the requested script of a phonon index.
"""
@inline script(::Val{:BD}, pid::PID, l::LaTeX) = l.body[pid.tag]
@inline script(::Val{:BD}, index::Index{<:Union{Int, Colon}, <:PID}, l::LaTeX) = l.body[index.iid.tag]
@inline script(::Val{:BD}, index::AbstractCompositeIndex{<:Index{<:Union{Int, Colon}, <:PID}}, l::LaTeX) = l.body[getcontent(index, :index).iid.tag]

"""
    latexofphonons

The default LaTeX format of a phonon index.
"""
const latexofphonons = LaTeX{(:direction,), (:site,)}(Dict('p'=>"p", 'u'=>"u"), "", "")
@inline latexname(::Type{<:Index{<:Union{Int, Colon}, <:PID}}) = Symbol("Index{Union{Int, Colon}, PID}")
@inline latexname(::Type{<:AbstractCompositeIndex{<:Index{<:Union{Int, Colon}, <:PID}}}) = Symbol("AbstractCompositeIndex{Index{Union{Int, Colon}, PID}}")
latexformat(Index{<:Union{Int, Colon}, <:PID}, latexofphonons)
latexformat(AbstractCompositeIndex{<:Index{<:Union{Int, Colon}, <:PID}}, latexofphonons)
@inline latexname(::Type{<:PID}) = Symbol("PID")
latexformat(PID, LaTeX{(:direction,), ()}(Dict('p'=>"p", 'u'=>"u"), "", ""))

## Permutation
"""
    permute(id₁::CompositeIndex{<:Index{Int, PID{Char}}}, id₂::CompositeIndex{<:Index{Int, PID{Char}}}) -> Tuple{Vararg{Operator}}

Permute two phonon indexes and get the result.
"""
function permute(id₁::CompositeIndex{<:Index{Int, PID{Char}}}, id₂::CompositeIndex{<:Index{Int, PID{Char}}})
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
### requested by IIDSpace
@inline function shape(iidspace::IIDSpace{<:PID, Phonon})
    return (iidspace.iid.tag=='u' ? (1:1) : (2:2), pidrange(iidspace.iid.direction, iidspace.internal.ndirection))
end
@inline pidrange(::Union{Symbol, Colon}, ndirection::Int) = 1:ndirection
@inline function pidrange(direction::Char, ndirection::Int)
    direction = Int(direction) - Int('x') + 1
    @assert 0<direction<ndirection+1 "pidrange error: direction out of range."
    return direction:direction
end

### MatrixCoupling
"""
    MatrixCoupling(sites::Union{NTuple{2, Int}, Colon}, ::Type{PID}, matrix::AbstractMatrix)

Construct a set of `Coupling`s corresponding to the dynamical matrix of phonons.
"""
function MatrixCoupling(sites::Union{NTuple{2, Int}, Colon}, ::Type{PID}, matrix::AbstractMatrix)
    @assert size(matrix)[1]∈(1, 2, 3) && size(matrix)[2]∈(1, 2, 3) "MatrixCoupling error: mismatched dimension of input matrix."
    left = size(matrix)[1]==1 ? SVector('x') : size(matrix)[1]==2 ? SVector('x', 'y') : SVector('x', 'y', 'z')
    right = size(matrix)[2]==1 ? SVector('x') : size(matrix)[2]==2 ? SVector('x', 'y') : SVector('x', 'y', 'z')
    return MatrixCoupling(sites, PID, Component(SVector('u'), SVector('u'), default_matrix), Component(left, right, matrix))
end

### expand
"""
    expand(::Val{:Hooke}, pnc::Coupling{<:Number, <:NTuple{2, Index{<:Union{Int, Colon}, PID{Colon}}}}, bond::Bond, hilbert::Hilbert) -> PPExpand

Expand the default phonon potential coupling on a given bond.
"""
function expand(::Val{:Hooke}, pnc::Coupling{<:Number, <:NTuple{2, Index{<:Union{Int, Colon}, PID{Colon}}}}, bond::Bond, hilbert::Hilbert)
    R̂ = rcoordinate(bond)/norm(rcoordinate(bond))
    @assert pnc.indexes.iids.tags==('u', 'u') "expand error: wrong tags of Hooke coupling."
    @assert isapprox(pnc.value, 1, atol=atol, rtol=rtol) "expand error: wrong coefficient of Hooke coupling."
    @assert pnc.indexes.sites∈((1, 2), (:, :)) "expand error: wrong sites of Hooke coupling."
    pn₁ = filter(pnc.indexes[1].iid, hilbert[bond[1].site])
    pn₂ = filter(pnc.indexes[2].iid, hilbert[bond[2].site])
    @assert pn₁.ndirection==pn₂.ndirection==length(R̂) "expand error: mismatched number of directions."
    return PPExpand(R̂, (bond[1], bond[2]))
end
struct PPExpand{N, D<:Number} <: VectorSpace{Operator{D, ID{CompositeIndex{Index{Int, PID{Char}}, SVector{N, D}}, 2}}}
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
    Kinetic(id::Symbol, value; amplitude::Union{Function, Nothing}=nothing, modulate::Union{Function, Bool}=true)

Kinetic energy of phonons.

Type alias for `Term{:Kinetic, id, V, Int, C<:TermCoupling, A<:TermAmplitude, M<:TermModulate}`.
"""
const Kinetic{id, V, C<:TermCoupling, A<:TermAmplitude, M<:TermModulate} = Term{:Kinetic, id, V, Int, C, A, M}
@inline function Kinetic(id::Symbol, value; amplitude::Union{Function, Nothing}=nothing, modulate::Union{Function, Bool}=true)
    return Term{:Kinetic}(id, value, 0, Coupling(:, PID, ('p', 'p'), :), true; amplitude=amplitude, modulate=modulate)
end

"""
    Hooke(id::Symbol, value, bondkind; amplitude::Union{Function, Nothing}=nothing, modulate::Union{Function, Bool}=true)

Potential energy of phonons by the Hooke's law.

Type alias for `Term{:Hooke, id, V, B, C<:TermCoupling, A<:TermAmplitude, M<:TermModulate}`
"""
const Hooke{id, V, B, C<:TermCoupling, A<:TermAmplitude, M<:TermModulate} = Term{:Hooke, id, V, B, C, A, M}
@inline function Hooke(id::Symbol, value, bondkind; amplitude::Union{Function, Nothing}=nothing, modulate::Union{Function, Bool}=true)
    return Term{:Hooke}(id, value, bondkind, Coupling(:, PID, ('u', 'u'), :), true; amplitude=amplitude, modulate=modulate)
end

"""
    Elastic(id::Symbol, value, bondkind, coupling; amplitude::Union{Function, Nothing}=nothing, modulate::Union{Function, Bool}=true)

Generic elastic energy of phonons.

Type alias for `Term{:Elastic, id, V, B, C<:TermCoupling, A<:TermAmplitude, M<:TermModulate}`
"""
const Elastic{id, V, B, C<:TermCoupling, A<:TermAmplitude, M<:TermModulate} = Term{:Elastic, id, V, B, C, A, M}
@inline function Elastic(id::Symbol, value, bondkind, coupling; amplitude::Union{Function, Nothing}=nothing, modulate::Union{Function, Bool}=true)
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
